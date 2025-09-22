import os
import re
import ast
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio

# --- TEST MODE / RENDERER ---
TEST_MODE = False               # set True if you want a quick on-screen preview
TEST_LIMIT = 3
SAVE_TEST_OUTPUTS = False
pio.renderers.default = "browser"  # ensures fig.show() opens from terminal

# --- VISUAL OPTIONS ---
POINT_SIZE = 6
POINT_OPACITY = 0.75
VIOLIN_SHOW_BOX = True
JITTER = 0.35
INCLUDE_AGGREGATE_ALL = True   # <-- add an "All" super-pop row per (Gene,Disease)

# --- POPULATION PALETTE (1kG-style) ---
POP_COLOR = {
    'EUR': '#1f77b4',  # blue
    'EAS': '#2ca02c',  # green
    'SAS': '#9467bd',  # purple
    'AMR': '#d62728',  # red
    'AFR': '#ff7f0e',  # orange/yellow
    'All': '#7f7f7f',  # gray for aggregates
}
# Put "All" at the top, then a readable order
SUPERPOP_ORDER = ['All', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'Unknown']

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots/"
DATA_PATH = f"{BASE_DIR}Data/Other Data/83_loci_503_samples_withancestrycolumns.csv"
OUTPUT_DIR = f"{BASE_DIR}Results/Plots/Allele_Length_Violin_Swarm"
OUTPUT_DIR_PNG = os.path.join(OUTPUT_DIR, "PNG")
OUTPUT_DIR_HTML = os.path.join(OUTPUT_DIR, "HTML")
os.makedirs(OUTPUT_DIR_PNG, exist_ok=True)
os.makedirs(OUTPUT_DIR_HTML, exist_ok=True)

if TEST_MODE:
    OUTPUT_DIR = os.path.join(OUTPUT_DIR, "test_outputs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    OUTPUT_DIR_PNG = OUTPUT_DIR
    OUTPUT_DIR_HTML = OUTPUT_DIR

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 125
PNG_SCALE = 2

# --- Load data ---
df = pd.read_csv(DATA_PATH)

# --- Optional: create an aggregated "All" super-pop row per (Gene,Disease) ---
if INCLUDE_AGGREGATE_ALL:
    # Keep columns needed in plots/labels; duplicate rows with SuperPop='All'
    keep_cols = df.columns.tolist()
    agg = df.copy()
    agg['SuperPop'] = 'All'
    df = pd.concat([df, agg], ignore_index=True)

# --- Ensure SuperPop has a stable categorical order ---
if 'SuperPop' in df.columns:
    # include only known categories + any unexpected ones appended
    uniques = df['SuperPop'].astype(str).unique().tolist()
    cats = [c for c in SUPERPOP_ORDER if c in uniques]
    extras = [c for c in uniques if c not in cats]
    df['SuperPop'] = pd.Categorical(df['SuperPop'], categories=cats + extras, ordered=True)

# --- Aggregate setup (kept for compatibility) ---
total_count = (
    df.groupby(['Disease', 'Gene', 'SuperPop', 'Allele length'], observed=False)['Sample ID']
      .nunique()
      .reset_index(name='total_count')
)

def count_affected(df_in, inh_mode, min_alleles):
    # Pathogenic bands may be absent for some rows; guard with non-null mask
    m = (
        df_in['Inheritance'].astype(str).eq(f"['{inh_mode}']")
        & df_in['Pathogenic min'].notna()
        & df_in['Pathogenic max'].notna()
    )
    df_inh = df_in[m].copy()
    if df_inh.empty:
        return pd.DataFrame(columns=['Disease','Gene','SuperPop',f'{inh_mode.lower()}_affected_counts'])

    df_path = df_inh[
        (df_inh['Allele length'] >= df_inh['Pathogenic min'])
        & (df_inh['Allele length'] <= df_inh['Pathogenic max'])
    ].copy()

    if df_path.empty:
        return pd.DataFrame(columns=['Disease','Gene','SuperPop',f'{inh_mode.lower()}_affected_counts'])

    grouped = (
        df_path.groupby(['Disease', 'Gene', 'SuperPop', 'Sample ID'], observed=False)
               .size()
               .reset_index(name='path_count')
    )
    affected = grouped[grouped['path_count'] >= min_alleles]
    return (
        affected.groupby(['Disease', 'Gene', 'SuperPop'], observed=False)['Sample ID']
                .nunique()
                .reset_index(name=f'{inh_mode.lower()}_affected_counts')
    )

# Compute affected counts (safe even if empty)
ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

# Merge into summary and compute percentages
df_agg = total_count.copy()
for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
    if not affected_df.empty:
        df_agg = df_agg.merge(
            affected_df,
            on=['Disease', 'Gene', 'SuperPop'],
            how='left'
        )

for inh in ['ad', 'ar', 'xd', 'xr']:
    if f'{inh}_affected_counts' not in df_agg.columns:
        df_agg[f'{inh}_affected_counts'] = 0
    df_agg[f'{inh}_affected_counts'] = df_agg[f'{inh}_affected_counts'].fillna(0).astype(int)
    df_agg[f'{inh}_percentage'] = np.where(
        df_agg['total_count'] > 0,
        df_agg[f'{inh}_affected_counts'] / df_agg['total_count'] * 100,
        np.nan
    )

# --- Plotting helper ---
def create_violin_beeswarm(original_df, gene, disease):
    subset = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)].copy()
    if subset.empty:
        return None

    # Inheritance parsing (optional; only used in title)
    raw_inh = subset['Inheritance']
    try:
        inh_list = ast.literal_eval(raw_inh.iloc[0]) if not raw_inh.empty else []
        inheritance = inh_list[0] if inh_list else ''
    except Exception:
        inheritance = str(raw_inh.iloc[0]) if not raw_inh.empty else ''
    if not inheritance:
        inheritance = "Unknown"
    if inheritance not in ['AD', 'AR', 'XD', 'XR', 'Unknown']:
        inheritance = "Unknown"

    # Totals per population (for subtitle)
    relevant = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]
    unique_pops = (
        relevant.groupby('SuperPop', observed=False)['Sample ID']
                .nunique()
                .reset_index(name='pop_total')
                .sort_values(by='SuperPop')
    )
    pop_desc = ', '.join(f"{row['SuperPop']}: {row['pop_total']}" for _, row in unique_pops.iterrows())

    # Wrap subtitle nicely
    max_line_len = 107
    pop_lines, line = [], ""
    for segment in pop_desc.split(", "):
        if len(line) + len(segment) > max_line_len:
            pop_lines.append(line.strip(", "))
            line = ""
        line += segment + ", "
    if line:
        pop_lines.append(line.strip(", "))
    pop_size_html = "<br>".join([
        f"<span style='font-size:12px'>Total Individuals per Population: {ln}</span>" if i == 0
        else f"<span style='font-size:12px'>{ln}</span>"
        for i, ln in enumerate(pop_lines)
    ])

    disease_html = f"<span style='font-size:12px'>{gene} - {disease}</span>"
    title_html = (
        f"<span style='font-size:18px; font-weight:bold'>Allele Lengths per Population</span><br>"
        f"{disease_html}<br>{pop_size_html}<br>"
        f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
    )

    # Violin colored by SuperPop (uses our palette)
    fig = px.violin(
        subset,
        x='Allele length',
        y='SuperPop',
        orientation='h',
        box=VIOLIN_SHOW_BOX,
        points=False,
        color='SuperPop',
        category_orders={'SuperPop': [c for c in SUPERPOP_ORDER if c in subset['SuperPop'].cat.categories]},
        color_discrete_map=POP_COLOR,
        title=title_html,
    )

    # Beeswarm points, also colored by SuperPop
    strip_fig = px.strip(
        subset,
        x='Allele length',
        y='SuperPop',
        orientation='h',
        color='SuperPop',
        category_orders={'SuperPop': [c for c in SUPERPOP_ORDER if c in subset['SuperPop'].cat.categories]},
        color_discrete_map=POP_COLOR,
        hover_data=['Sample ID', 'Allele length', 'SuperPop', 'Gene', 'Disease']
    )

    # Manual jitter to avoid vertical stacking
    for tr in strip_fig.data:
        x_vals = np.array(tr.x, dtype=float)
        if JITTER and JITTER > 0:
            noise = np.random.uniform(-JITTER, JITTER, size=x_vals.shape)
            x_vals = x_vals + noise
        tr.x = x_vals
        tr.update(
            marker=dict(size=POINT_SIZE, opacity=POINT_OPACITY, line=dict(width=0)),
            showlegend=False
        )
        fig.add_trace(tr)

    fig.update_layout(
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        margin=dict(t=TOP_MARGIN),
        autosize=False,
        xaxis=dict(title="Allele Length", ticks='outside', showline=True, linecolor='black', zeroline=False),
        yaxis=dict(title="", ticks='outside', showline=True, linecolor='black'),
        plot_bgcolor='white',
        font=dict(color='black'),
        title=dict(y=0.95),
        showlegend=False
    )

    return fig

# --- Generate plots ---
printed = 0
pairs = df[['Gene', 'Disease']].drop_duplicates().itertuples(index=False, name=None)
for gene, disease in pairs:
    fig = create_violin_beeswarm(df, gene, disease)
    if fig is None:
        continue

    safe_gene = re.sub(r'[\\/]', '_', gene)
    safe_disease = re.sub(r'[\\/]', '_', disease)
    png_path = os.path.join(OUTPUT_DIR_PNG, f"{safe_gene}_{safe_disease}_allele_length_violin_beeswarm.png")
    html_path = os.path.join(OUTPUT_DIR_HTML, f"{safe_gene}_{safe_disease}_allele_length_violin_beeswarm.html")

    if TEST_MODE:
        print(f"Previewing: {gene} / {disease}")
        fig.show()  # should open in your default browser
        if SAVE_TEST_OUTPUTS:
            fig.write_html(html_path)
            try:
                fig.write_image(png_path, width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
            except Exception as e:
                print(f"[PNG export skipped] {e}")
        printed += 1
        if printed >= TEST_LIMIT:
            break
    else:
        fig.write_html(html_path)
        try:
            fig.write_image(png_path, width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
        except Exception as e:
            print(f"[PNG export skipped] {e}")

print("--- Test mode ON: Test completed ---" if TEST_MODE else "--- Done ---")
