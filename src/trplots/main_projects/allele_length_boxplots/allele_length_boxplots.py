"""
---------------------------------------------
 Script: Allele Length Boxplot Generator
 Purpose:
   - Reads a CSV of allele lengths with ancestry and metadata
   - Flags pathogenic individuals using per-locus thresholds
   - Aggregates counts and computes % pathogenic by population
   - Builds horizontal boxplots of allele length by population
   - Saves plots as PNG and HTML
   - In test mode, previews a subset for quick checks
---------------------------------------------
"""

import re
import ast
import pandas as pd
import plotly.express as px
import plotly.io as pio

from trplots.config import (
    OTHER_DATA,                 # data/other_data
    ENSURE_DIR                  # results helper
)

# make Plotly open figures in the browser by default
pio.renderers.default = "browser"

# --- TEST MODE ---
TEST_MODE = True               # Toggle this flag for quick testing (preview only)
TEST_LIMIT = 3                 # How many (gene, disease) plots to generate in test mode
SAVE_TEST_OUTPUTS = False       # Toggle saving plots when in test mode

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130               # fixed header space (annotations), keeps plot area aligned
PNG_SCALE = 2

# --- Input paths ---
DATA_PATH = OTHER_DATA / "83_loci_503_samples_withancestrycolumns.csv"

# --- Output paths ---
OUTPUT_DIR_PNG  = ENSURE_DIR("plots", "allele_length_boxplots", "allele_length_boxplots", "png")
OUTPUT_DIR_HTML = ENSURE_DIR("plots", "allele_length_boxplots", "allele_length_boxplots", "html")

# If test mode: override to a single test_outputs folder under results/plots/allele_length_boxplots
if TEST_MODE:
    TEST_OUTPUT   = ENSURE_DIR("plots", "allele_length_boxplots", "allele_length_boxplots", "test_outputs")
    OUTPUT_DIR_PNG = TEST_OUTPUT
    OUTPUT_DIR_HTML = TEST_OUTPUT

# --- POPULATION PALETTE (1kG-style) ---
POP_COLOR = {
    'EUR': '#1f77b4',  # blue
    'EAS': '#2ca02c',  # green
    'SAS': '#9467bd',  # purple
    'AMR': '#d62728',  # red
    'AFR': '#ff7f0e',  # orange/yellow
    'Unknown': '#7f7f7f',  # gray for unknowns
    'All': '#ff99cc',  # fallback for all
}
SUPERPOP_ORDER = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'Unknown', 'All']

# --- Helper: count affected individuals per inheritance model ---
def count_affected(df_in, inh_mode, min_alleles):
    """
    For a given inheritance model (e.g., 'AD'), count individuals
    who meet the 'Is pathogenic' criteria.
    - AD/XD: need >=1 pathogenic allele (min_alleles=1)
    - AR/XR: need >=2 pathogenic alleles (min_alleles=2)
    NOTE: This function matches stringified lists like "['AD']" as in source data.
    """
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

def _wrap_to_lines(s: str, max_len: int = 110):
    """Split a comma-separated string into ~max_len lines, returning a list of lines."""
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

def create_boxplot(filtered_df, gene, disease, original_df):
    """
    Build a horizontal boxplot:
      - X: Allele length
      - Y: SuperPop
    Uses 1kG-style colors and a readable SuperPop order with 'All' first.
    Fixed-size header via annotations to keep plot area aligned across figures.
    """
    if filtered_df.empty:
        return None

    # Inheritance for subtitle; data stores as stringified list
    raw_inh = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]['Inheritance']
    try:
        inh_list = ast.literal_eval(raw_inh.iloc[0]) if not raw_inh.empty else []
        inheritance = inh_list[0] if inh_list else ''
    except Exception:
        inheritance = str(raw_inh.iloc[0]) if not raw_inh.empty else ''
    inheritance = inheritance if inheritance in ['AD', 'AR', 'XD', 'XR'] else 'Unknown'

    # Enforce order (only present categories)
    present = filtered_df['SuperPop'].dropna().unique().tolist()
    ordered_categories = [c for c in SUPERPOP_ORDER if c in present]

    # Ordered categoricals
    fdf = filtered_df.copy()
    fdf['SuperPop'] = pd.Categorical(fdf['SuperPop'], categories=ordered_categories, ordered=True)

    relevant_samples = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)].copy()
    relevant_samples['SuperPop'] = pd.Categorical(relevant_samples['SuperPop'], categories=ordered_categories, ordered=True)

    # Totals per population (subtitle), same order as plot
    counts = (
        relevant_samples
        .groupby('SuperPop', observed=True)['Sample ID']
        .nunique()
        .reindex(ordered_categories)
    ).dropna().astype(int)
    pop_desc = ', '.join(f"{pop}: {counts.loc[pop]}" for pop in counts.index)
    pop_lines = _wrap_to_lines(pop_desc, max_len=110)

    # ---- Build the boxplot (no built-in title) ----
    fig = px.box(
        fdf,
        x='Allele length',
        y='SuperPop',
        color='SuperPop',
        category_orders={"SuperPop": ordered_categories},
        color_discrete_map=POP_COLOR,
        title=None,
        points=False,
    )

    # Consistent sizing & style
    fig.update_layout(
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        margin=dict(t=TOP_MARGIN, r=40, b=60, l=80),
        autosize=False,
        plot_bgcolor='white',
        showlegend=False,
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
                y=y0 - 0.06*i, yref="paper", yanchor="top",
                showarrow=False, align="left"
            )
        )
    fig.update_layout(annotations=annos)
    return fig

def main():
    # Load data
    df = pd.read_csv(DATA_PATH)

    # --- Add 'All' population ---
    df_all = df.copy()
    df_all["SuperPop"] = "All"
    df = pd.concat([df, df_all], ignore_index=True)

    # --- Flag pathogenic individuals using thresholds ---
    # True if allele length is within [Pathogenic min, Pathogenic max]
    if 'Allele length' not in df.columns or 'Pathogenic min' not in df.columns or 'Pathogenic max' not in df.columns:
        raise ValueError("DATA file missing one of required columns: 'Allele length', 'Pathogenic min', 'Pathogenic max'")

    df['Is pathogenic'] = (
        (df['Allele length'] >= df['Pathogenic min']) &
        (df['Allele length'] <= df['Pathogenic max'])
    )

    # --- Aggregate setup ---
    # Count unique samples per (Disease, Gene, SuperPop, Allele length)
    total_count = (
        df.groupby(['Disease', 'Gene', 'SuperPop', 'Allele length'])['Sample ID']
          .nunique()
          .reset_index(name='total_count')
    )

    # --- Compute affected counts for each model ---
    ad_affected = count_affected(df, 'AD', 1)
    ar_affected = count_affected(df, 'AR', 2)
    xd_affected = count_affected(df, 'XD', 1)
    xr_affected = count_affected(df, 'XR', 2)

    # --- Merge affected counts into main summary ---
    df_agg = total_count.copy()
    for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
        df_agg = df_agg.merge(
            affected_df,
            on=['Disease', 'Gene', 'SuperPop'],
            how='left'
        )

    # Fill missing counts with 0 and compute percentages
    for inh in ['ad', 'ar', 'xd', 'xr']:
        df_agg[f'{inh}_affected_counts'] = df_agg[f'{inh}_affected_counts'].fillna(0).astype(int)
        df_agg[f'{inh}_percentage'] = df_agg[f'{inh}_affected_counts'] / df_agg['total_count'] * 100

    # --- Generate plots ---
    printed = 0
    for gene in df_agg['Gene'].unique():
        diseases = df_agg[df_agg['Gene'] == gene]['Disease'].unique()
        for disease in diseases:
            sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]
            fig = create_boxplot(sub_df.copy(), gene, disease, df)
            if not fig:
                continue

            safe_gene = re.sub(r'[\\/]', '_', gene)
            safe_disease = re.sub(r'[\\/]', '_', disease)

            png_path  = OUTPUT_DIR_PNG  / f"{safe_gene}_{safe_disease}_allele_length_boxplot.png"
            html_path = OUTPUT_DIR_HTML / f"{safe_gene}_{safe_disease}_allele_length_boxplot.html"

            if TEST_MODE:
                print(f"Previewing: {gene} / {disease}")
                fig.show(renderer="browser")
                if SAVE_TEST_OUTPUTS:
                    fig.write_html(str(html_path))
                    try:
                        fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
                    except Exception as e:
                        print(f"Failed to write PNG {png_path}: {e} — ensure 'kaleido' is installed")
                printed += 1
                if printed >= TEST_LIMIT:
                    break
            else:
                fig.write_html(str(html_path))
                try:
                    fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
                except Exception as e:
                    print(f"Failed to write PNG {png_path}: {e} — ensure 'kaleido' is installed")
        if TEST_MODE and printed >= TEST_LIMIT:
            break

    print("--- Test mode ON: Test completed ---" if TEST_MODE else "--- Done ---")

if __name__ == "__main__":
    main()
