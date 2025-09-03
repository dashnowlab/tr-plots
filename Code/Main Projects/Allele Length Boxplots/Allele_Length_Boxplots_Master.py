"""
---------------------------------------------
 Script: Allele Length Boxplot Generator (Locus-based)
 Purpose:
   - Reads an Excel master spreadsheet of allele lengths w/ ancestry & metadata
   - Defines a Locus per row (Chromosome:Position; fallback to Disease ID)
   - Flags pathogenic individuals using per-locus thresholds
   - Aggregates counts and computes % pathogenic by population
   - Builds horizontal boxplots of allele length by population
   - Saves plots as PNG and HTML (one per locus)
   - Test mode: previews a subset
---------------------------------------------
"""

import pandas as pd
import re
import plotly.express as px
import os
import ast
import numpy as np

# --- TEST MODE ---
TEST_MODE = False               # Toggle this flag for quick testing (preview only)
TEST_LIMIT = 3                  # How many locus plots to generate in test mode
SAVE_TEST_OUTPUTS = False       # Toggle saving plots when in test mode

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 150               # a little extra header space to fit the locus line
PNG_SCALE = 2

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots/"
DATA_PATH = f"{BASE_DIR}Data/Other Data/Allele_Master_Spreadsheet.xlsx"
SHEET_NAME = "Sheet1"

OUTPUT_DIR = f"{BASE_DIR}Results/Plots/Allele_Length_Boxplots_Locus2"
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

# --- Helpers ---
def _norm_inh(val):
    """
    Normalize Inheritance values to one of {AD, AR, XD, XR} or 'Unknown'.
    Accepts values like "['AD']", "AD", "['AR','something']".
    """
    if pd.isna(val):
        return 'Unknown'
    modes = {'AD', 'AR', 'XD', 'XR'}
    s = str(val).strip()
    try:
        parsed = ast.literal_eval(s)
        if isinstance(parsed, (list, tuple)):
            for m in parsed:
                m = str(m).strip().upper()
                if m in modes:
                    return m
        s = str(parsed)
    except Exception:
        pass
    s = s.replace('[', '').replace(']', '').replace("'", '').replace('"', '').strip().upper()
    first = s.split(',')[0].strip()
    return first if first in modes else 'Unknown'

def _wrap_to_lines(s: str, max_len: int = 110):
    """Split a comma-separated string into ~max_len lines, returning a list of lines."""
    parts = [p.strip() for p in s.split(",")]
    lines, line = [], ""
    for seg in parts:
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

def _make_locus(row):
    """
    Preferred locus ID: chr{Chromosome}:{Position}
    Fallbacks if missing: Disease ID, else Gene|Disease|Chromosome:Position
    """
    chrom = str(row.get('Chromosome', '')).strip()
    pos = str(row.get('Position', '')).strip()
    disease_id = str(row.get('Disease ID', '')).strip()

    chrom_valid = chrom not in ('', 'nan', 'None')
    pos_valid = pos not in ('', 'nan', 'None')

    if chrom_valid and pos_valid:
        # Ensure 'Chr' prefix and integer-like position formatting if possible
        chrom_str = f"Chr{chrom}" if not str(chrom).lower().startswith('chr') else str(chrom)
        try:
            pos_int = int(float(pos))
            pos_str = f"{pos_int}"
        except Exception:
            pos_str = pos
        return f"{chrom_str}:{pos_str}"

    if disease_id and disease_id.lower() != 'nan':
        return disease_id

    # final fallback (should be rare)
    gene = str(row.get('Gene', '')).strip()
    disease = str(row.get('Disease', '')).strip()
    return f"{gene}|{disease}|{chrom}:{pos}"

# --- Load data ---
df = pd.read_excel(DATA_PATH, sheet_name=SHEET_NAME)

# Ensure numeric columns are numeric
for col in ['Allele length', 'Pathogenic min', 'Pathogenic max']:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

# Drop rows with missing essentials
df = df.dropna(subset=['Gene', 'Disease', 'SuperPop', 'Sample ID',
                       'Allele length', 'Pathogenic min', 'Pathogenic max'])

# Normalize inheritance
df['Inheritance_norm'] = df['Inheritance'].apply(_norm_inh)

# Create Locus column
df['Locus'] = df.apply(_make_locus, axis=1)

# --- Add 'All' population ---
df_all = df.copy()
df_all["SuperPop"] = "All"
df = pd.concat([df, df_all], ignore_index=True)

# --- Flag pathogenic individuals using thresholds ---
df['Is pathogenic'] = (
    (df['Allele length'] >= df['Pathogenic min']) &
    (df['Allele length'] <= df['Pathogenic max'])
)

# --- Aggregate per Locus ---
group_keys_base = ['Disease', 'Gene', 'Locus', 'SuperPop']

# Total unique samples per allele length bin (per locus & pop)
total_count = (
    df.groupby(group_keys_base + ['Allele length'], observed=True)['Sample ID']
      .nunique()
      .reset_index(name='total_count')
)

def count_affected(df_in, inh_mode, min_alleles):
    """
    Count affected individuals per locus & pop given inheritance model:
      - AD/XD: >=1 pathogenic allele
      - AR/XR: >=2 pathogenic alleles
    """
    df_inh = df_in[df_in['Inheritance_norm'] == inh_mode]
    df_path = df_inh[df_inh['Is pathogenic'] == True]
    grouped = (
        df_path.groupby(group_keys_base + ['Sample ID'], observed=True)
               .size()
               .reset_index(name='path_count')
    )
    affected = grouped[grouped['path_count'] >= min_alleles]
    return (
        affected.groupby(group_keys_base, observed=True)['Sample ID']
                .nunique()
                .reset_index(name=f'{inh_mode.lower()}_affected_counts')
    )

# Affected counts
ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

# Merge into main summary (keys include Locus)
df_agg = total_count.copy()
for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
    df_agg = df_agg.merge(affected_df, on=group_keys_base, how='left')

# Fill missing counts with 0 and compute percentages
for inh in ['ad', 'ar', 'xd', 'xr']:
    colc = f'{inh}_affected_counts'
    df_agg[colc] = df_agg[colc].fillna(0).astype(int)
    df_agg[f'{inh}_percentage'] = df_agg[colc] / df_agg['total_count'] * 100

def create_boxplot(filtered_df, gene, disease, locus, original_df):
    """
    Build a horizontal boxplot for one Locus:
      - X: Allele length
      - Y: SuperPop
    """
    if filtered_df.empty:
        return None

    # Subset rows for this locus
    subset = original_df[
        (original_df['Gene'] == gene) &
        (original_df['Disease'] == disease) &
        (original_df['Locus'] == locus)
    ]

    # Inheritance (normalized)
    raw_inh = subset['Inheritance_norm']
    inheritance = raw_inh.iloc[0] if not raw_inh.empty else 'Unknown'
    if inheritance not in ['AD', 'AR', 'XD', 'XR']:
        inheritance = 'Unknown'

    # SuperPop presence & order
    present = filtered_df['SuperPop'].dropna().unique().tolist()
    ordered_categories = [c for c in SUPERPOP_ORDER if c in present]

    fdf = filtered_df.copy()
    fdf['SuperPop'] = pd.Categorical(fdf['SuperPop'], categories=ordered_categories, ordered=True)

    relevant_samples = subset.copy()
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

    # ---- Build the boxplot ----
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
    subtitle_lines = [
        f"{gene} - {disease}",
        f"{locus}",
    ]
    if pop_lines:
        subtitle_lines.append(f"Total Individuals per Population: {pop_lines[0]}")
        subtitle_lines.extend(pop_lines[1:])
    subtitle_lines.append(f"Inheritance: {inheritance}")

    annos = [
        dict(
            text=main_title, x=0, xref="paper", xanchor="left",
            y=1.42, yref="paper", yanchor="top",
            showarrow=False, align="left", font=dict(size=18)
        )
    ]
    y0 = 1.34
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

# --- Generate plots: one per Locus ---
printed = 0
# Group to unique triplets
for (gene, disease, locus), _grp in df_agg.groupby(['Gene', 'Disease', 'Locus'], sort=False):
    sub_df = df_agg[
        (df_agg['Gene'] == gene) &
        (df_agg['Disease'] == disease) &
        (df_agg['Locus'] == locus)
    ]
    fig = create_boxplot(sub_df.copy(), gene, disease, locus, df)
    if fig is None:
        continue

    safe_gene = re.sub(r'[\\/]', '_', str(gene))
    safe_disease = re.sub(r'[\\/]', '_', str(disease))
    safe_locus = re.sub(r'[\\/:*?"<>|]+', '_', str(locus))[:120]  # trim long names
    png_path = os.path.join(OUTPUT_DIR_PNG, f"{safe_gene}_{safe_disease}_{safe_locus}_allele_length_boxplot.png")
    html_path = os.path.join(OUTPUT_DIR_HTML, f"{safe_gene}_{safe_disease}_{safe_locus}_allele_length_boxplot.html")

    if TEST_MODE:
        print(f"Previewing: {gene} / {disease} / {locus}")
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

if TEST_MODE:
    print("--- Test mode ON: Test completed ---")
else:
    print("--- Done ---")
