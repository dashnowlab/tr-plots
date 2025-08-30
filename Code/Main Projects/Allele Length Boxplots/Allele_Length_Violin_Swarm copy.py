"""
---------------------------------------------
 Script: Allele Length Violin+Swarm Generator
 Purpose:
   - Reads a CSV of allele lengths with ancestry and metadata
   - Flags pathogenic individuals using per-locus thresholds
   - Aggregates counts (same as boxplot script for consistency)
   - Builds horizontal violin plots with overlaid points ("swarm")
   - Saves plots as PNG and HTML
   - In test mode, previews a subset for quick checks
---------------------------------------------
"""

import pandas as pd
import re
import plotly.express as px
import os
import ast

# --- TEST MODE ---
TEST_MODE = True               # Toggle this flag for quick testing (preview only)
TEST_LIMIT = 3                 # How many (gene, disease) plots to generate in test mode
SAVE_TEST_OUTPUTS = False      # Toggle saving plots when in test mode

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130               # fixed header space (annotations), keeps plot area aligned
PNG_SCALE = 2

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots/"

DATA_PATH = f"{BASE_DIR}Data/Other Data/83_loci_503_samples_withancestrycolumns.csv"
OUTPUT_DIR = f"{BASE_DIR}Results/Plots/Allele_Length_Violin_Swarm"

# Normal mode: save into PNG and HTML subfolders
OUTPUT_DIR_PNG = os.path.join(OUTPUT_DIR, "PNG")
OUTPUT_DIR_HTML = os.path.join(OUTPUT_DIR, "HTML")
os.makedirs(OUTPUT_DIR_PNG, exist_ok=True)
os.makedirs(OUTPUT_DIR_HTML, exist_ok=True)

# If test mode: override to a single test_outputs folder
if TEST_MODE:
    OUTPUT_DIR = os.path.join(OUTPUT_DIR, "test_outputs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    OUTPUT_DIR_PNG = OUTPUT_DIR
    OUTPUT_DIR_HTML = OUTPUT_DIR

# --- POPULATION PALETTE (1kG-style) ---
POP_COLOR = {
    'EUR': '#1f77b4',  # blue
    'EAS': '#2ca02c',  # green
    'SAS': '#9467bd',  # purple
    'AMR': '#d62728',  # red
    'AFR': '#ff7f0e',  # orange/yellow
    'All': '#7f7f7f',  # gray for aggregates
    'Unknown': '#ff99cc',  # fallback for unknowns (optional)
}
SUPERPOP_ORDER = ['All', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'Unknown']

# --- Load data ---
df = pd.read_csv(DATA_PATH)

# --- Add 'All' population ---
df_all = df.copy()
df_all["SuperPop"] = "All"
df = pd.concat([df, df_all], ignore_index=True)

# --- Flag pathogenic individuals using thresholds ---
df['Is pathogenic'] = (
    (df['Allele length'] >= df['Pathogenic min']) &
    (df['Allele length'] <= df['Pathogenic max'])
)

# --- Aggregate setup (to match boxplot script; not strictly required for violin) ---
total_count = (
    df.groupby(['Disease', 'Gene', 'SuperPop', 'Allele length'])['Sample ID']
      .nunique()
      .reset_index(name='total_count')
)

def count_affected(df_in, inh_mode, min_alleles):
    df_inh = df_in[df_in['Inheritance'].astype(str) == f"['{inh_mode}']"]
    df_path = df_inh[df_inh['Is pathogenic'] == True]
    grouped = (
        df_path.groupby(['Disease', 'Gene', 'SuperPop', 'Sample ID'])
               .size()
               .reset_index(name='path_count')
    )
    affected = grouped[grouped['path_count'] >= min_alleles]
    return (
        affected.groupby(['Disease', 'Gene', 'SuperPop'])['Sample ID']
                .nunique()
                .reset_index(name=f'{inh_mode.lower()}_affected_counts')
    )

ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

df_agg = total_count.copy()
for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
    df_agg = df_agg.merge(
        affected_df,
        on=['Disease', 'Gene', 'SuperPop'],
        how='left'
    )

for inh in ['ad', 'ar', 'xd', 'xr']:
    df_agg[f'{inh}_affected_counts'] = df_agg[f'{inh}_affected_counts'].fillna(0).astype(int)
    df_agg[f'{inh}_percentage'] = df_agg[f'{inh}_affected_counts'] / df_agg['total_count'] * 100

def _wrap_to_lines(s: str, max_len: int = 110):
    parts = [p.strip() for p in s.split(",")]
    lines, line = [], ""
    for seg in parts:
        seg = seg.strip()
        if not seg:
            continue
        add = (", " if line else "") + seg
        if len(line) + len(seg) + (2 if line else 0) > max_len:
            if line:
                lines.append(line)
            line = seg
        else:
            line += add if line else seg
    if line:
        lines.append(line)
    return lines

def create_violin_swarm(filtered_df, gene, disease, original_df):
    """
    Build a horizontal violin plot with overlaid points ("swarm"):
      - X: Allele length
      - Y: SuperPop
    Fixed-size header via annotations to keep plot area aligned across figures.
    """
    if filtered_df.empty:
        return None

    raw_inh = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]['Inheritance']
    try:
        inh_list = ast.literal_eval(raw_inh.iloc[0]) if not raw_inh.empty else []
        inheritance = inh_list[0] if inh_list else ''
    except Exception:
        inheritance = str(raw_inh.iloc[0]) if not raw_inh.empty else ''
    inheritance = inheritance if inheritance in ['AD', 'AR', 'XD', 'XR'] else 'Unknown'

    present = filtered_df['SuperPop'].dropna().unique().tolist()
    ordered_categories = [c for c in SUPERPOP_ORDER if c in present]

    fdf = filtered_df.copy()
    fdf['SuperPop'] = pd.Categorical(fdf['SuperPop'], categories=ordered_categories, ordered=True)

    relevant_samples = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)].copy()
    relevant_samples['SuperPop'] = pd.Categorical(relevant_samples['SuperPop'], categories=ordered_categories, ordered=True)

    counts = (
        relevant_samples
        .groupby('SuperPop', observed=True)['Sample ID']
        .nunique()
        .reindex(ordered_categories)
    ).dropna().astype(int)
    pop_desc = ', '.join(f"{pop}: {counts.loc[pop]}" for pop in counts.index)
    pop_lines = _wrap_to_lines(pop_desc, max_len=110)

    # ---- Build the violin + swarm (no built-in title) ----
    fig = px.violin(
        fdf,
        x='Allele length',
        y='SuperPop',
        color='SuperPop',
        category_orders={"SuperPop": ordered_categories},
        color_discrete_map=POP_COLOR,
        box=False,                 # violin only; box off (can set True if desired)
        points='all',              # show points (swarm-like)
        title=None,
    )

    # Make the point cloud look denser (swarm-ish)
    fig.update_traces(
        jitter=0.35,                # horizontal jitter of points
        marker=dict(size=4, line=dict(width=0)),
        meanline_visible=True
    )

    # Consistent sizing & style
    fig.update_layout(
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        margin=dict(t=TOP_MARGIN, r=40, b=60, l=80),
        autosize=False,
        plot_bgcolor='white',
        showlegend=False,  # keep consistent
        font=dict(color='black'),
        xaxis=dict(title="Allele Length", ticks='outside', showline=True, linecolor='black'),
        yaxis=dict(title="", ticks='outside', showline=True, linecolor='black'),
    )

    # ---- Fixed-position header via annotations ----
    main_title = "<b>Allele Lengths per Population</b>"
    subtitle_lines = [f"{gene} - {disease}"] + \
                     [f"Total Individuals per Population: {pop_lines[0]}"] + \
                     pop_lines[1:] + \
                     [f"Inheritance: {inheritance}"]

    annos = [
        dict(
            text=main_title, x=0, xref="paper", xanchor="left",
            y=1.4, yref="paper", yanchor="top",
            showarrow=False, align="left", font=dict(size=18)
        )
    ]
    y0 = 1.32
    for i, line in enumerate(subtitle_lines):
        annos.append(
            dict(
                text=f"<span style='font-size:12px'>{line}</span>",
                x=0, xref="paper", xanchor="left",
                y=y0 - 0.04*i, yref="paper", yanchor="top",
                showarrow=False, align="left"
            )
        )
    fig.update_layout(annotations=annos)

    return fig

# --- Generate plots ---
printed = 0
for gene in df_agg['Gene'].unique():
    diseases = df_agg[df_agg['Gene'] == gene]['Disease'].unique()
    for disease in diseases:
        sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]
        fig = create_violin_swarm(sub_df.copy(), gene, disease, df)
        if not fig:
            continue

        safe_gene = re.sub(r'[\\/]', '_', gene)
        safe_disease = re.sub(r'[\\/]', '_', disease)
        png_path = os.path.join(OUTPUT_DIR_PNG, f"{safe_gene}_{safe_disease}_allele_length_violin_swarm.png")
        html_path = os.path.join(OUTPUT_DIR_HTML, f"{safe_gene}_{safe_disease}_allele_length_violin_swarm.html")

        if TEST_MODE:
            print(f"Previewing: {gene} / {disease}")
            fig.show()
            if SAVE_TEST_OUTPUTS:
                fig.write_html(html_path)
                fig.write_image(png_path, width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
            printed += 1
            if printed >= TEST_LIMIT:
                break
        else:
            fig.write_html(html_path)
            fig.write_image(png_path, width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
    if TEST_MODE and printed >= TEST_LIMIT:
        break

# --- Finished ---
print("--- Test mode ON: Test completed ---" if TEST_MODE else "--- Done ---")
