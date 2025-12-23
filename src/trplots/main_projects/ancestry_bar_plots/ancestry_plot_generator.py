"""
---------------------------------------------
 Script: Ancestry Bar Plot Generator
 Purpose:
   - Reads a dataset containing ancestry and genotype information 
     for multiple samples
   - Cleans and standardizes metadata (population labels, sex, inheritance)
   - Counts pathogenic individuals under AD/AR/XD/XR models
   - Computes percentages and 95% binomial confidence intervals
   - Builds horizontal bar plots with dropdowns for Pore type and Sex
   - Saves plots as interactive HTML (and shows a preview)
---------------------------------------------
"""

import pandas as pd
import re
import plotly.graph_objects as go
import plotly.io as pio
import os
import ast
import numpy as np
from statsmodels.stats import proportion
import argparse

pio.renderers.default = "browser"

# --- TEST MODE ---
TEST_MODE = True                # Toggle this flag for quick testing (only one plot generated)
TEST_LIMIT = 3                   # How many (gene,disease) plots in test mode
SAVE_TEST_OUTPUTS = False         # Toggle saving plots when in test mode

def parse_args():
    p = argparse.ArgumentParser(description="Ancestry bar plot generator")
    p.add_argument("--test", dest="test", action="store_true", help="Enable test mode")
    p.add_argument("--no-test", dest="test", action="store_false", help="Disable test mode")
    p.set_defaults(test=TEST_MODE)
    p.add_argument("--test-limit", dest="test_limit", type=int, default=TEST_LIMIT, help="How many plots to generate in test mode")
    p.add_argument("--save-test-outputs", dest="save_test_outputs", action="store_true", default=SAVE_TEST_OUTPUTS, help="Save outputs even when in test mode")
    p.add_argument("--output-html", dest="output_html", type=str, default=None, help="Override HTML output directory")
    p.add_argument("--output-png", dest="output_png", type=str, default=None, help="Override PNG output directory")
    return p.parse_args()

from trplots.config import OTHER_DATA, OUTPUT_BASE
from pathlib import Path

# --- File locations ---
ALLELE_SPREADSHEET_PATH = "/Users/annelisethorn/Documents/github/tr-plots/data/other_data/allele_spreadsheet.xlsx"
OUTPUT_BASE = Path(OUTPUT_BASE) / "plots" / "new_ancestry_plot_generator"

if TEST_MODE:
    # In test mode: just one folder, no subfolders
    OUTPUT_DIR = os.path.join(OUTPUT_BASE, "test_outputs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    OUTPUT_HTML_DIR = OUTPUT_DIR
    OUTPUT_PNG_DIR  = OUTPUT_DIR
else:
    # Normal mode: structured HTML/PNG subfolders
    OUTPUT_DIR = OUTPUT_BASE
    OUTPUT_HTML_DIR = os.path.join(OUTPUT_DIR, "HTML")
    OUTPUT_PNG_DIR  = os.path.join(OUTPUT_DIR, "PNG")
    os.makedirs(OUTPUT_HTML_DIR, exist_ok=True)
    os.makedirs(OUTPUT_PNG_DIR,  exist_ok=True)

# --- Load allele spreadsheet ---
df = pd.read_excel(ALLELE_SPREADSHEET_PATH, engine="openpyxl")

# --- Normalize column names (case-insensitive lookup) ---
cols_lower = {c.lower(): c for c in df.columns}
def col(*names):
    for n in names:
        if n.lower() in cols_lower:
            return cols_lower[n.lower()]
    return None

# Map required columns
sample_id_col = col("Sample ID", "Sample_ID", "sample_id")
gene_col = col("Gene")
disease_col = col("Disease")
repeat_col = col("Repeat Count", "RepeatCount", "repeat_count")
pathogenic_min_col = col("Pathogenic min", "Pathogenic_min", "pathogenic_min")
pathogenic_max_col = col("Pathogenic max", "Pathogenic_max", "pathogenic_max")
inheritance_col = col("Inheritance")
sex_col = col("Sex")
population_col = col("Population", "Cleaned population description", "Ancestry")
pore_col = col("Pore")

if not all([sample_id_col, gene_col, disease_col, repeat_col]):
    raise SystemExit(f"Spreadsheet missing required columns. Found: {list(df.columns)}")

# Rename columns for consistency with downstream code
df = df.rename(columns={
    sample_id_col: "Sample ID",
    gene_col: "Gene",
    disease_col: "Disease",
    repeat_col: "Repeat count",
    pathogenic_min_col: "Pathogenic min",
    pathogenic_max_col: "Pathogenic max",
    inheritance_col: "Inheritance",
    sex_col: "Sex",
    population_col: "Cleaned population description",
    pore_col: "Pore"
})

# Fill missing Pore/Sex/Population with defaults
if "Pore" not in df.columns or df["Pore"].isna().all():
    df["Pore"] = "All Types"
if "Sex" not in df.columns or df["Sex"].isna().all():
    df["Sex"] = "All Sexes"
if "Cleaned population description" not in df.columns or df["Cleaned population description"].isna().all():
    df["Cleaned population description"] = "All"

df["Pore"] = df["Pore"].fillna("All Types").astype(str)
df["Sex"] = df["Sex"].fillna("All Sexes").astype(str)
df["Cleaned population description"] = df["Cleaned population description"].fillna("All").astype(str)

# --- Flag pathogenic individuals (based on repeat counts) ---
df['Is pathogenic'] = df.apply(
    lambda row: (
        row.get('Pathogenic min') <= row.get('Repeat count') <= row.get('Pathogenic max')
        if (pd.notna(row.get('Pathogenic min')) and pd.notna(row.get('Pathogenic max')) and pd.notna(row.get('Repeat count')))
        else np.nan
    ),
    axis=1
)

# --- Clean inheritance column (remove stringified lists like "['AD']") ---
def clean_inheritance(x):
    if isinstance(x, str) and x.startswith('['):
        try:
            val = ast.literal_eval(x)
            return val[0] if val else np.nan
        except (ValueError, SyntaxError):
            return np.nan
    return x

if "Inheritance" in df.columns:
    mask = df['Inheritance'].astype(str).str.startswith('[')
    df.loc[mask, 'Inheritance'] = df.loc[mask, 'Inheritance'].apply(clean_inheritance)
else:
    df["Inheritance"] = "Unknown"

# --- Standardize sex labels ---
sex_map = {
    'm': 'Male', 'male': 'Male', 'M': 'Male', 'MALE': 'Male',
    'f': 'Female', 'female': 'Female', 'F': 'Female', 'FEMALE': 'Female',
    'u': 'Unknown', 'unknown': 'Unknown', 'U': 'Unknown', 'UNK': 'Unknown', 'Unk': 'Unknown'
}
df['Sex'] = (
    df['Sex']
    .astype(str)
    .str.strip()
    .map(lambda s: sex_map.get(s, s))
)
df.loc[df['Sex'].isin(['', 'nan', 'None']), 'Sex'] = 'Unknown'

# --- Aggregate counts per Disease / Gene / Population / Pore / Sex ---
total_count = df.groupby(
    ['Disease', 'Gene', 'Cleaned population description', 'Pore', 'Sex']
)['Sample ID'].nunique().reset_index(name='total_count')

# --- Count pathogenic individuals under each inheritance model ---
def count_affected(df_in, inh_mode, min_alleles):
    df_inh = df_in[df_in['Inheritance'] == inh_mode]
    df_path = df_inh[df_inh['Is pathogenic'] == True]
    grouped = df_path.groupby(
        ['Disease', 'Gene', 'Cleaned population description', 'Pore', 'Sex', 'Sample ID']
    ).size().reset_index(name='path_count')
    affected = grouped[grouped['path_count'] >= min_alleles]
    return affected.groupby(
        ['Disease', 'Gene', 'Cleaned population description', 'Pore', 'Sex']
    )['Sample ID'].nunique().reset_index(name=f'{inh_mode.lower()}_affected_counts')

ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

# --- Merge affected counts into main summary ---
df_agg = total_count.copy()
for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
    df_agg = df_agg.merge(
        affected_df,
        on=['Disease', 'Gene', 'Cleaned population description', 'Pore', 'Sex'],
        how='left'
    )

# --- Compute percentages and 95% CIs ---
for inh in ['ad', 'ar', 'xd', 'xr']:
    df_agg[f'{inh}_affected_counts'] = df_agg[f'{inh}_affected_counts'].fillna(0).astype(int)
    df_agg[f'{inh}_percentage'] = df_agg[f'{inh}_affected_counts'] / df_agg['total_count'] * 100

def binomial_ci(x, n, confidence=0.95):
    if n == 0:
        return (0.0, 0.0)
    lower, upper = proportion.proportion_confint(x, n, alpha=1-confidence, method='wilson')
    return lower * 100, upper * 100

for inh in ['ad', 'ar', 'xd', 'xr']:
    lowers, uppers = [], []
    for x, n in zip(df_agg[f'{inh}_affected_counts'], df_agg['total_count']):
        lo, up = binomial_ci(x, n)
        lowers.append(lo)
        uppers.append(up)
    df_agg[f'{inh}_ci_lower'] = lowers
    df_agg[f'{inh}_ci_upper'] = uppers

# --- Compute "All" population summary (collapsed across subgroups) ---
group_cols = ['Gene', 'Disease', 'Pore']
agg_cols = ['total_count'] + \
           [f'{inh}_affected_counts' for inh in ['ad', 'ar', 'xd', 'xr']] + \
           [f'{inh}_percentage' for inh in ['ad', 'ar', 'xd', 'xr']]

df_all_rows = (
    df_agg[df_agg['Cleaned population description'] != 'All']
    .groupby(group_cols, as_index=False)[agg_cols]
    .sum()
)
df_all_rows['Cleaned population description'] = 'All'  # acts across all sexes

# --- Recalculate confidence intervals for 'All' ---
for inh in ['ad', 'ar', 'xd', 'xr']:
    lowers, uppers = [], []
    for x, n in zip(df_all_rows[f'{inh}_affected_counts'], df_all_rows['total_count']):
        lo, up = binomial_ci(x, n)
        lowers.append(lo)
        uppers.append(up)
    df_all_rows[f'{inh}_ci_lower'] = lowers
    df_all_rows[f'{inh}_ci_upper'] = uppers

# --- Append 'All' back into main aggregated table ---
df_agg = pd.concat([df_agg, df_all_rows], ignore_index=True)

# --- (Optional) Helper to compute fully combined rows across pore types ---
def add_combined_row(group_df):
    """
    Creates combined rows for each (Gene, Disease) across pore types,
    adding 'All' population and recomputing percentages and CIs.
    (Kept here for completeness; not used downstream.)
    """
    result = []
    for (gene, disease), group in group_df.groupby(['Gene', 'Disease']):
        for pore in ['R9', 'R10', 'All Types']:
            if pore == 'All Types':
                df_sub = df[(df['Gene'] == gene) & (df['Disease'] == disease)]
            else:
                df_sub = df[(df['Gene'] == gene) & (df['Disease'] == disease) & (df['Pore'] == pore)]

            if df_sub.empty:
                continue

            unique_samples = df_sub['Sample ID'].nunique()

            combined = {
                'Gene': gene,
                'Disease': disease,
                'Pore': pore,
                'Cleaned population description': 'All',
                'total_count': unique_samples
            }

            for inh in ['ad', 'ar', 'xd', 'xr']:
                inh_mode = inh.upper()
                min_alleles = 1 if inh in ['ad', 'xd'] else 2

                df_inh = df_sub[df_sub['Inheritance'] == inh_mode]
                df_path = df_inh[df_inh['Is pathogenic'] == True]
                grouped = df_path.groupby('Sample ID').size().reset_index(name='path_count')
                affected = grouped[grouped['path_count'] >= min_alleles]
                affected_count = affected['Sample ID'].nunique()

                combined[f'{inh}_affected_counts'] = affected_count
                combined[f'{inh}_percentage'] = (affected_count / unique_samples * 100 if unique_samples > 0 else 0)
                lo, hi = binomial_ci(affected_count, unique_samples)
                combined[f'{inh}_ci_lower'] = lo
                combined[f'{inh}_ci_upper'] = hi

            result.append(combined)
    return pd.DataFrame(result)

# --- Define population order for plotting ---
all_pops = sorted(set(df['Cleaned population description'].dropna().astype(str).tolist()))
if 'All' in all_pops:
    all_pops.remove('All')
superpop_order = ['All', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

# --- Build plot for a (Gene, Disease) pair with Pore+Sex dropdowns ---
def create_horizontal_bar_plot(filtered_df, gene, disease, original_df, superpop_order):
    """
    Builds a horizontal bar plot showing % pathogenic individuals across populations,
    with dropdowns for Pore type and Sex. Ensures all populations appear on Y axis.
    """
    # Determine inheritance mode for labeling
    raw_inh = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]['Inheritance']
    inheritance = raw_inh.iloc[0] if not raw_inh.empty else ''
    if inheritance not in ['AD', 'AR', 'XD', 'XR']:
        print(f"Skipping unknown inheritance for {gene} / {disease}")
        return None

    inh = inheritance.lower()
    pore_order = ['All Types', 'R9', 'R10']
    sex_order = ['All Sexes', 'Male', 'Female', 'Unknown']

    # Build traces and dynamic subtitles for every (pore, sex) combo
    trace_map = {}
    subtitle_map = {}

    for pore in pore_order:
        for sex in sex_order:
            df_sel = filtered_df.copy()
            if pore != 'All Types':
                df_sel = df_sel[df_sel['Pore'] == pore]
            if sex != 'All Sexes':
                df_sel = df_sel[df_sel['Sex'] == sex]

            # Ignore pre-made 'All' rows; recompute per selection
            df_sel = df_sel[df_sel['Cleaned population description'] != 'All']

            # Aggregate counts for this selection (or empty frame if none)
            if not df_sel.empty:
                grp = df_sel.groupby('Cleaned population description', as_index=False).agg({
                    'total_count': 'sum',
                    f'{inh}_affected_counts': 'sum'
                })
            else:
                grp = pd.DataFrame(columns=['Cleaned population description', 'total_count', f'{inh}_affected_counts'])

            # Fresh 'All' population row for this (pore, sex)
            all_total = grp['total_count'].sum() if not grp.empty else 0
            all_affected = grp[f'{inh}_affected_counts'].sum() if not grp.empty else 0
            all_row = pd.DataFrame({
                'Cleaned population description': ['All'],
                'total_count': [all_total],
                f'{inh}_affected_counts': [all_affected]
            })
            df_grouped = pd.concat([all_row, grp], ignore_index=True)

            # Ensure all populations appear (use passed-in order)
            df_grouped = (
                df_grouped
                .set_index('Cleaned population description')
                .reindex(superpop_order)
                .reset_index()
            )

            # Fill missing counts with zero
            for col in ['total_count', f'{inh}_affected_counts']:
                if col in df_grouped.columns:
                    df_grouped[col] = df_grouped[col].fillna(0)

            # Percentages & CIs
            with np.errstate(divide='ignore', invalid='ignore'):
                df_grouped['percentage'] = np.where(
                    df_grouped['total_count'] > 0,
                    df_grouped[f'{inh}_affected_counts'] / df_grouped['total_count'] * 100,
                    0.0
                )
            df_grouped['percentage_display'] = df_grouped['percentage'].apply(
                lambda x: f"{int(x)}" if x == int(x) else f"{x:.2f}"
            )
            df_grouped['affected'] = df_grouped[f'{inh}_affected_counts']

            # CI per row
            ci_bounds = df_grouped.apply(
                lambda row: binomial_ci(int(row['affected']), int(row['total_count'])),
                axis=1
            )
            df_grouped['ci_lower'] = [lo for lo, hi in ci_bounds]
            df_grouped['ci_upper'] = [hi for lo, hi in ci_bounds]
            df_grouped['error_plus'] = df_grouped['ci_upper'] - df_grouped['percentage']
            df_grouped['error_minus'] = df_grouped['percentage'] - df_grouped['ci_lower']

            trace_map[(pore, sex)] = df_grouped

            # Dynamic subtitle with total counts per population (wrapped)
            pop_lines, line = [], ""
            for _, row in df_grouped.iterrows():
                part = f"{row['Cleaned population description']}: {int(row['total_count'])}"
                if len(line + part) > 100:
                    pop_lines.append(line.rstrip(', '))
                    line = ""
                line += part + ", "
            if line:
                pop_lines.append(line.rstrip(', '))
            subtitle_map[(pore, sex)] = "<br>".join(
                [f"<span style='font-size:12px'>{l}</span>" for l in pop_lines]
            )

    if not trace_map:
        return None

    # Build the figure (one trace per (pore, sex))
    fig = go.Figure()
    trace_keys = list(trace_map.keys())

    for key in trace_keys:
        df_grouped = trace_map[key]
        fig.add_bar(
            x=df_grouped['percentage'],
            y=df_grouped['Cleaned population description'],
            orientation='h',
            error_x=dict(
                array=df_grouped['error_plus'],
                arrayminus=df_grouped['error_minus'],
                thickness=1
            ),
            marker_color="lightblue",
            customdata=np.stack([
                df_grouped['affected'],
                df_grouped['total_count'],
                df_grouped['percentage_display'],
                df_grouped['ci_lower'],
                df_grouped['ci_upper']
            ], axis=-1),
            hovertemplate=(
                "Population: %{y}<br>"
                "# Pathogenic Individuals: %{customdata[0]}<br>"
                "# Total Individuals: %{customdata[1]}<br>"
                "Pathogenic Genotype: %{customdata[2]}%<br>"
                "95% CI: [%{customdata[3]:.2f}%, %{customdata[4]:.2f}%]<extra></extra>"
            ),
            visible=False
        )

    # Default selection: All Types • All Sexes (if present)
    default_key = ('All Types', 'All Sexes') if ('All Types', 'All Sexes') in trace_keys else trace_keys[0]
    default_idx = trace_keys.index(default_key)
    fig.data[default_idx].visible = True

    # Dropdown to switch (pore, sex); also updates dynamic subtitle
    buttons = []
    for i, (pore, sex) in enumerate(trace_keys):
        visible = [j == i for j in range(len(trace_keys))]
        buttons.append(dict(
            label=f"{pore} • {sex}",
            method="update",
            args=[
                {"visible": visible},
                {"annotations": [
                    dict(
                        text="Total Individuals per Population: " + subtitle_map[(pore, sex)],
                        x=-0.06, y=1.15, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=12, color="black"),
                        align="left", xanchor="left", name="pop_subtitle"
                    ),
                    dict(
                        text="Pore Type & Sex:",
                        x=1.15, y=1.23, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=12, color="black"),
                        xanchor="center"
                    )
                ]}
            ]
        ))

    # Layout & axes styling
    fig.update_layout(
        width=900, height=500, margin=dict(t=125), bargap=0.5,
        plot_bgcolor='white', showlegend=False, font=dict(color='black'),
        title=dict(
            y=0.95,
            text=(
                f"<span style='font-size:18px; font-weight:bold'>Pathogenic Genotype Distribution</span><br>"
                f"<span style='font-size:12px'>{gene} - {disease}</span><br>"
                f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
            )
        ),
        updatemenus=[dict(
            active=default_idx, buttons=buttons, direction="down",
            showactive=True, x=1.15, xanchor="center", y=1.15, yanchor="top"
        )],
        annotations=[
            dict(
                text="Total Individuals per Population: " + subtitle_map[default_key],
                x=-0.0555, y=1.15, xref="paper", yref="paper",
                showarrow=False, font=dict(size=12, color="black"),
                align="left", xanchor="left", name="pop_subtitle"
            ),
            dict(
                text="Pore & Sex:", x=1.15, y=1.23, xref="paper", yref="paper",
                showarrow=False, font=dict(size=12, color="black"),
                xanchor="center"
            )
        ],
        xaxis=dict(
            range=[0, 100], tickmode='linear', tick0=0, dtick=10,
            ticks='outside', showline=True, linecolor='black', linewidth=1,
            zeroline=False, tickfont=dict(color='black'),
            title=dict(text="Pathogenic Genotypes (%)", font=dict(color='black'))
        ),
        yaxis=dict(
            title="",
            showline=True,
            linecolor='black',
            categoryorder='array',
            categoryarray=superpop_order   # force order from caller
        )
    )
    return fig

# --- Generate plots ---
def main():
    # Determine dynamic population order from the loaded dataframe `df`
    all_pops = sorted(set(df['Cleaned population description'].dropna().astype(str).tolist()))
    # Preserve a canonical order for major superpops, append any extras found in data
    canonical = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    extras = [p for p in all_pops if p not in canonical and p != 'All']
    superpop_order = ['All'] + [p for p in canonical if p in all_pops] + extras

    # --- Generate plots (previously top-level code) ---
    made = 0
    for gene in df_agg['Gene'].unique():
        for disease in df_agg[df_agg['Gene'] == gene]['Disease'].unique():
            sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]

            # Skip if no meaningful data
            if sub_df['total_count'].sum() == 0:
                continue

            fig = create_horizontal_bar_plot(sub_df.copy(), gene, disease, df, superpop_order)
            if not fig:
                continue

            safe_gene = re.sub(r'[\\/]', '_', gene)
            safe_disease = re.sub(r'[\\/]', '_', disease)

            html_path = os.path.join(OUTPUT_HTML_DIR, f"{safe_gene}_{safe_disease}_ancestry_plot.html")
            png_path  = os.path.join(OUTPUT_PNG_DIR, f"{safe_gene}_{safe_disease}_ancestry_plot.png")

            # Save files
            fig.write_html(html_path)
            try:
                fig.write_image(png_path, scale=2)   # scale=2 makes higher resolution
            except Exception:
                pass

            if TEST_MODE:
                # Preview only in test mode
                print(f"Previewing: {gene} / {disease}")
                fig.show(renderer="browser")

                # Optional save in test mode
                if SAVE_TEST_OUTPUTS:
                    print(f"Saving: {gene} / {disease}")
                    fig.write_html(html_path)

                made += 1
                if made >= TEST_LIMIT:
                    break
        if TEST_MODE and made >= TEST_LIMIT:
            break

    # --- Finished ---
    if TEST_MODE:
        print("--- Test mode ON: Test completed ---")
    else:
        print("--- Done ---")

# Guard to avoid execution on import
if __name__ == "__main__":
    main()