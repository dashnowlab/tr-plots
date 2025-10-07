"""
---------------------------------------------
 Script: Allele Length Boxplot Generator (Locus-based, All Alleles)
 Purpose:
   - Reads an Excel master spreadsheet of allele lengths w/ ancestry & metadata
   - Defines a Locus per row (Chromosome:Position; fallback to Disease ID)
   - Builds horizontal boxplots of **all** allele lengths by population
   - Forces ALL populations to show even if a population has no data for a locus
   - Saves plots as PNG and HTML (one per locus)
   - Test mode: previews a subset
---------------------------------------------
"""

import re
import ast
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from pathlib import Path

# import shared paths/helpers from config
from trplots.config import (
    OTHER_DATA,                
    ENSURE_DIR,
    ALLELE_LENGTH_PLOTS_OUTPUT
)

pio.renderers.default = "browser"

# --- TEST MODE ---
TEST_MODE = False               # Toggle this flag for quick testing (preview only)
TEST_LIMIT = 3                 # How many locus plots to generate in test mode
SAVE_TEST_OUTPUTS = True       # Toggle saving plots when in test mode (PNG/HTML)

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130               # space for the header lines
PNG_SCALE = 2

# --- File locations (via config) ---
# Excel master lives under data/other_data
DATA_PATH = (OTHER_DATA / "allele_master_spreadsheet2.xlsx")
SHEET_NAME = "Sheet1"

# --- Output roots (under results/...) ---
# Normal (non-test) output -> results/plots/allele_length_boxplots/allele_length_boxplots_master/{png,html}
OUTPUT_DIR_PNG  = ENSURE_DIR("plots", "allele_length_boxplots", "allele_length_boxplots_master", "png")
OUTPUT_DIR_HTML = ENSURE_DIR("plots", "allele_length_boxplots", "allele_length_boxplots_master", "html")

# If test mode: override to a single test_outputs folder under results/plots/allele_length_boxplots
if TEST_MODE:
    TEST_OUTPUT = ENSURE_DIR("plots", "allele_length_boxplots", "allele_length_boxplots_master", "test_outputs")
    OUTPUT_DIR_PNG = TEST_OUTPUT
    OUTPUT_DIR_HTML = TEST_OUTPUT

# --- POPULATION PALETTE (1kG-style) ---
POP_COLOR = {
    'EUR': '#1f77b4',  # blue
    'EAS': '#2ca02c',  # green
    'SAS': '#9467bd',  # purple
    'AMR': '#d62728',  # red
    'AFR': '#ff7f0e',  # orange/yellow
    'Unknown': '#7f7f7f',  # gray for unknown
    'All': '#ff99cc',  # pink for all
}
SUPERPOP_ORDER = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'Unknown', 'All']

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
        chrom_str = f"Chr{chrom}" if not str(chrom).lower().startswith('chr') else str(chrom)
        try:
            pos_int = int(float(pos))
            pos_str = f"{pos_int}"
        except Exception:
            pos_str = pos
        return f"{chrom_str}:{pos_str}"

    if disease_id and disease_id.lower() != 'nan':
        return disease_id

    gene = str(row.get('Gene', '')).strip()
    disease = str(row.get('Disease', '')).strip()
    return f"{gene}|{disease}|{chrom}:{pos}"

# --- Load data ---
df = pd.read_excel(DATA_PATH, sheet_name=SHEET_NAME)

# Ensure numeric columns are numeric
if 'Allele length' in df.columns:
    df['Allele length'] = pd.to_numeric(df['Allele length'], errors='coerce')

# Normalize inheritance (used only for subtitle context)
if 'Inheritance' in df.columns:
    df['Inheritance_norm'] = df['Inheritance'].apply(_norm_inh)
else:
    df['Inheritance_norm'] = 'Unknown'

# Create Locus column
df['Locus'] = df.apply(_make_locus, axis=1)

# --- Add 'All' population (aggregate superpop row copied from the raw table) ---
df_all = df.copy()
df_all["SuperPop"] = "All"
df = pd.concat([df, df_all], ignore_index=True)

# --- Safety: ensure required columns exist ---
for required in ['Gene', 'Disease', 'Locus', 'SuperPop', 'Sample ID', 'Allele length']:
    if required not in df.columns:
        raise ValueError(f"Missing required column: '{required}'")

# -------------------- Plotting (ALL alleles) --------------------
def create_boxplot_all(gene, disease, locus, raw_df, show_no_data_note=True, show_points=False):
    """
    Build a horizontal boxplot for one Locus using **all allele lengths**,
    showing ALL populations even if they have no data.
    """
    subset = raw_df[
        (raw_df['Gene'] == gene) &
        (raw_df['Disease'] == disease) &
        (raw_df['Locus'] == locus)
    ].copy()

    ordered_categories = SUPERPOP_ORDER[:]
    subset['SuperPop'] = pd.Categorical(
        subset['SuperPop'],
        categories=ordered_categories,
        ordered=True
    )

    # If subset is empty (no rows at all), create a dummy frame to keep axis happy
    plot_df = subset if not subset.empty else pd.DataFrame({
        'Allele length': [np.nan],
        'SuperPop': pd.Categorical(['All'], categories=ordered_categories, ordered=True)
    })

    # Inheritance (normalized)
    base_inh = raw_df[
        (raw_df['Gene'] == gene) & (raw_df['Disease'] == disease) & (raw_df['Locus'] == locus)
    ]['Inheritance_norm']
    inheritance = base_inh.iloc[0] if not base_inh.empty else 'Unknown'
    if inheritance not in ['AD', 'AR', 'XD', 'XR']:
        inheritance = 'Unknown'

    # Counts per population (unique individuals)
    counts = (
        subset.groupby('SuperPop', observed=True)['Sample ID']
              .nunique()
              .reindex(ordered_categories)
              .fillna(0)
              .astype(int)
    )
    pop_desc = ', '.join(f"{pop}: {counts.loc[pop]}" for pop in ordered_categories)
    pop_lines = _wrap_to_lines(pop_desc, max_len=110)

    # --------- Build figure with explicit box stats + separate outlier points ----------
    fig = go.Figure()

    def five_number_summary(x):
        x = np.asarray(x, dtype=float)
        x = x[~np.isnan(x)]
        if x.size == 0:
            return None
        q1 = np.percentile(x, 25)
        med = np.percentile(x, 50)
        q3 = np.percentile(x, 75)
        iqr = q3 - q1
        lf = q1 - 1.5 * iqr
        uf = q3 + 1.5 * iqr
        lower_whisk = x[x >= lf].min() if np.any(x >= lf) else x.min()
        upper_whisk = x[x <= uf].max() if np.any(x <= uf) else x.max()
        outliers = x[(x < lf) | (x > uf)]
        return q1, med, q3, lower_whisk, upper_whisk, outliers

    for pop in ordered_categories:
        xs = subset.loc[subset["SuperPop"] == pop, "Allele length"].to_numpy(dtype=float)
        xs = xs[~np.isnan(xs)]
        if xs.size == 0:
            fig.add_trace(
                go.Box(
                    orientation="h",
                    y=[pop],
                    x=None,
                    name="",
                    marker_color=POP_COLOR.get(pop, "#7f7f7f"),
                    showlegend=False,
                    boxpoints=False,
                    hoverinfo="skip"
                )
            )
            continue

        stats = five_number_summary(xs)
        if stats is None:
            continue
        q1, med, q3, lw, uw, out = stats

        fig.add_trace(
            go.Box(
                orientation="h",
                name="",
                y=[pop],
                q1=[q1],
                median=[med],
                q3=[q3],
                lowerfence=[lw],
                upperfence=[uw],
                marker_color=POP_COLOR.get(pop, "#7f7f7f"),
                boxpoints=False,
                showlegend=False,
                hoveron="boxes",
                hovertemplate=(
                    "Population: %{y}<br>"
                    f"Min: {lw:.3g}<br>"
                    f"Q1:  {q1:.3g}<br>"
                    f"Median: {med:.3g}<br>"
                    f"Q3:  {q3:.3g}<br>"
                    f"Max: {uw:.3g}<extra></extra>"
                ),
            )
        )

        if out.size:
            fig.add_trace(
                go.Scatter(
                    x=out,
                    y=[pop] * out.size,
                    mode="markers",
                    name="",
                    marker=dict(size=6, color=POP_COLOR.get(pop, "#7f7f7f")),
                    showlegend=False,
                    hovertemplate="Allele length: %{x}<br>Population: %{y}<extra></extra>",
                )
            )

    fig.update_yaxes(
        categoryorder='array',
        categoryarray=ordered_categories,
        tickmode='array',
        tickvals=ordered_categories,
        autorange='reversed',
        title=""
    )

    missing = [c for c in ordered_categories if counts.get(c, 0) == 0]
    if show_no_data_note and missing:
        for cat in missing:
            fig.add_annotation(
                x=0.0, xref='paper', xanchor='left',
                y=cat, yref='y', yanchor='middle',
                text="<span style='font-size:11px;color:#888'>no data</span>",
                showarrow=False, align='left'
            )

    fig.update_layout(
        hovermode="closest",
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        margin=dict(t=TOP_MARGIN, r=40, b=60, l=90),
        autosize=False,
        plot_bgcolor='white',
        showlegend=False,
        font=dict(color='black'),
        xaxis=dict(title="Allele Length (bp)", ticks='outside', showline=True, linecolor='black'),
        yaxis=dict(ticks='outside', showline=True, linecolor='black'),
    )

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

# --- Generate plots: one per Locus (ALL alleles) ---
printed = 0

unique_triplets = (
    df[['Gene', 'Disease', 'Locus']]
    .dropna()
    .drop_duplicates()
)

for _, row in unique_triplets.iterrows():
    gene, disease, locus = row['Gene'], row['Disease'], row['Locus']
    fig = create_boxplot_all(gene, disease, locus, df, show_no_data_note=True, show_points=False)

    # Safe filenames
    safe_gene = re.sub(r'[\\/]', '_', str(gene))
    safe_disease = re.sub(r'[\\/]', '_', str(disease))
    safe_locus = re.sub(r'[\\/:*?"<>|]+', '_', str(locus))[:120]

    png_path  = OUTPUT_DIR_PNG  / f"{safe_gene}_{safe_disease}_{safe_locus}_allele_length_boxplot_ALL.png"
    html_path = OUTPUT_DIR_HTML / f"{safe_gene}_{safe_disease}_{safe_locus}_allele_length_boxplot_ALL.html"

    if TEST_MODE:
        print(f"Previewing: {gene} / {disease} / {locus}")
        fig.show()
        if SAVE_TEST_OUTPUTS:
            fig.write_html(str(html_path))
            fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
        printed += 1
        if printed >= TEST_LIMIT:
            break
    else:
        fig.write_html(str(html_path))
        fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)

if TEST_MODE:
    print("--- Test mode ON: Test completed ---")
else:
    print("--- Done ---")
