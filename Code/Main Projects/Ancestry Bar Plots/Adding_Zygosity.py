"""
---------------------------------------------
 Script: Ancestry Bar Plot Generator (with Zygosity)
 Purpose:
   - Reads a dataset containing ancestry and genotype information 
     for multiple samples
   - Cleans and standardizes metadata (population labels, sex, inheritance, zygosity)
   - Counts pathogenic individuals under AD/AR/XD/XR models (zygosity-aware)
   - Computes percentages and 95% binomial confidence intervals
   - Builds two-panel plots:
       Left: % affected with CIs (dropdown: Pore, Sex)
       Right: Zygosity composition among pathogenic rows (HOM/HET/HEMI)
   - Saves plots as interactive HTML and PNG (and shows a preview in TEST_MODE)
---------------------------------------------
"""

import pandas as pd
import re
import plotly.graph_objects as go
import os
import ast
import numpy as np
from statsmodels.stats import proportion
from plotly.subplots import make_subplots

# --- TEST MODE ---
TEST_MODE = True                # Toggle this flag for quick testing (only one plot generated)
TEST_LIMIT = 1                  # How many (gene,disease) plots in test mode
SAVE_TEST_OUTPUTS = True        # Toggle saving plots when in test mode

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots"

SEQ_DATA_PATH = f"{BASE_DIR}/Data/Sequencing Data/83 Loci 503 Samples/83_loci_503_samples_with_sex.xlsx"
PORE_PATH = f"{BASE_DIR}/Data/Other Data/1KGP_ONT_500_Summary_Sample_ID_Pore.csv"
OUTPUT_BASE = os.path.join(BASE_DIR, "Results/Plots/Ancestry_Bar_Plots")

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

# --- Load main dataset ---
df = pd.read_excel(SEQ_DATA_PATH)

# --- Merge pore information (R9 / R10) ---
pore_df = pd.read_csv(PORE_PATH)
pore_df['Base Sample ID'] = pore_df['Sample_ID'].apply(lambda x: x.split('-')[0] if isinstance(x, str) else x)
df['Base Sample ID'] = df['Sample ID'].apply(lambda x: x.split('-')[0])
df = df.merge(pore_df[['Base Sample ID', 'Pore']], on='Base Sample ID', how='left')

# --- Use ancestry information ---
if 'SuperPop' in df.columns:
    df['Cleaned population description'] = df['SuperPop']
elif 'Cleaned population description' not in df.columns:
    raise KeyError("Missing both 'SuperPop' and 'Cleaned population description' columns.")

# --- Flag pathogenic individuals (based on repeat counts) ---
df['Is pathogenic'] = df.apply(
    lambda row: (
        row['Pathogenic min'] <= row['Repeat count'] <= row['Pathogenic max']
        if pd.notna(row['Pathogenic min']) and pd.notna(row['Pathogenic max']) and pd.notna(row['Repeat count'])
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

mask = df['Inheritance'].astype(str).str.startswith('[')
df.loc[mask, 'Inheritance'] = df.loc[mask, 'Inheritance'].apply(clean_inheritance)

# --- Standardize sex labels ---
if 'Sex' not in df.columns:
    df['Sex'] = 'All Sexes'
else:
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

# --- ZYGOSITY: normalize and derive allele weights ---------------------------
# If missing, create column (will behave as unknown)
if 'Zygosity' not in df.columns:
    df['Zygosity'] = np.nan

zyg_map = {
    'het': 'HET', 'heterozygous': 'HET',
    'hom': 'HOM', 'homozygous': 'HOM',
    'hemi': 'HEMI', 'hemizygous': 'HEMI',
    'unknown': 'UNK', 'unk': 'UNK', 'na': 'UNK', 'nan': 'UNK', 'none': 'UNK', '': 'UNK'
}

def normalize_zygosity(x):
    if pd.isna(x):
        return 'UNK'
    s = str(x).strip().lower()
    return zyg_map.get(s, s.upper())

df['Zygosity'] = df['Zygosity'].apply(normalize_zygosity)

# Per-row allele weight (only counts if row is pathogenic)
allele_weight = {'HET': 1, 'HEMI': 1, 'HOM': 2}
df['pathogenic_alleles'] = np.where(
    df['Is pathogenic'] == True,
    df['Zygosity'].map(allele_weight).fillna(0).astype(int),
    0
)

# --- Aggregate counts per Disease / Gene / Population / Pore / Sex ---
total_count = df.groupby(
    ['Disease', 'Gene', 'Cleaned population description', 'Pore', 'Sex']
)['Sample ID'].nunique().reset_index(name='total_count')

# --- Count pathogenic individuals under each inheritance model (zygosity-aware) ---
def count_affected(df_in, inh_mode):
    """
    Zygosity-aware affected counting.
    AD/XD: affected if sum(pathogenic_alleles) >= 1
    AR:    affected if sum(pathogenic_alleles) >= 2
    XR:    affected if (Male: >=1) (Female: >=2) (Unknown: >=2, conservative)
    """
    df_inh = df_in[df_in['Inheritance'] == inh_mode].copy()
    if df_inh.empty:
        return pd.DataFrame(columns=[
            'Disease','Gene','Cleaned population description','Pore','Sex', f'{inh_mode.lower()}_affected_counts'
        ])

    # Sum allele counts per sample across rows
    per_sample = df_inh.groupby(
        ['Disease','Gene','Cleaned population description','Pore','Sex','Sample ID'],
        as_index=False
    )['pathogenic_alleles'].sum()

    if inh_mode in ('AD', 'XD'):
        per_sample['is_affected'] = per_sample['pathogenic_alleles'] >= 1
    elif inh_mode == 'AR':
        per_sample['is_affected'] = per_sample['pathogenic_alleles'] >= 2
    elif inh_mode == 'XR':
        def xr_rule(row):
            if row['Sex'] == 'Male':
                return row['pathogenic_alleles'] >= 1
            elif row['Sex'] == 'Female':
                return row['pathogenic_alleles'] >= 2
            else:
                return row['pathogenic_alleles'] >= 2
        per_sample['is_affected'] = per_sample.apply(xr_rule, axis=1)
    else:
        per_sample['is_affected'] = False

    affected = per_sample[per_sample['is_affected']]
    out = affected.groupby(
        ['Disease','Gene','Cleaned population description','Pore','Sex']
    )['Sample ID'].nunique().reset_index(name=f'{inh_mode.lower()}_affected_counts')
    return out

ad_affected = count_affected(df, 'AD')
ar_affected = count_affected(df, 'AR')
xd_affected = count_affected(df, 'XD')
xr_affected = count_affected(df, 'XR')

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

# --- Define population order for plotting ---
all_pops = sorted(set(df['Cleaned population description'].dropna().astype(str).tolist()))
if 'All' in all_pops:
    all_pops.remove('All')
superpop_order = ['All', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

# --- Helper: zygosity breakdowns for hover & composition (pathogenic rows only) ---
def zygosity_breakdown(original_df, gene, disease, pore, sex, inheritance, superpop_order):
    """
    Returns DF with columns:
      Cleaned population description, zyg_HOM, zyg_HET, zyg_HEMI
    Counts are number of PATHOGENIC rows per population by zygosity.
    Safe if 'Zygosity' is missing (returns zeros).
    """
    if 'Zygosity' not in original_df.columns:
        return pd.DataFrame({
            'Cleaned population description': superpop_order,
            'zyg_HOM': 0, 'zyg_HET': 0, 'zyg_HEMI': 0
        })

    df_sel = original_df[
        (original_df['Gene'] == gene) &
        (original_df['Disease'] == disease) &
        (original_df['Inheritance'] == inheritance) &
        (original_df['Is pathogenic'] == True)
    ].copy()

    if pore != 'All Types':
        df_sel = df_sel[df_sel['Pore'] == pore]
    if sex != 'All Sexes':
        df_sel = df_sel[df_sel['Sex'] == sex]

    if df_sel.empty:
        return pd.DataFrame({
            'Cleaned population description': superpop_order,
            'zyg_HOM': 0, 'zyg_HET': 0, 'zyg_HEMI': 0
        })

    # Normalize zyg labels
    zmap = {
        'hom':'HOM','homozygous':'HOM',
        'het':'HET','heterozygous':'HET',
        'hemi':'HEMI','hemizygous':'HEMI'
    }
    df_sel['Zygosity'] = df_sel['Zygosity'].astype(str).str.strip().str.lower().map(lambda s: zmap.get(s, s.upper()))
    ztab = (
        df_sel
        .groupby(['Cleaned population description','Zygosity'])
        .size()
        .unstack(fill_value=0)
    )

    for need in ['HOM','HET','HEMI']:
        if need not in ztab.columns:
            ztab[need] = 0

    out = (ztab.reset_index()
           .rename(columns={'HOM':'zyg_HOM','HET':'zyg_HET','HEMI':'zyg_HEMI'})
           .set_index('Cleaned population description')
           .reindex(superpop_order)
           .reset_index())
    return out

# --- Build plot for a (Gene, Disease) pair with Pore+Sex dropdowns (two panels) ---
def create_horizontal_bar_plot(filtered_df, gene, disease, original_df):
    """
    Builds a two-panel figure:
      Left: % pathogenic individuals across populations (with 95% CIs)
      Right: Zygosity composition (HOM/HET/HEMI) among pathogenic rows (100% stacked)
    Includes dropdowns for Pore type and Sex. Ensures all populations appear on Y axis.
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

            # Ensure all populations appear (fixed order)
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

            # Zygosity composition for hover & right panel
            ztab = zygosity_breakdown(original_df, gene, disease, pore, sex, inheritance, superpop_order)
            df_grouped = df_grouped.merge(ztab, on='Cleaned population description', how='left')
            for col in ['zyg_HOM','zyg_HET','zyg_HEMI']:
                if col not in df_grouped or df_grouped[col].isna().all():
                    df_grouped[col] = 0
                else:
                    df_grouped[col] = df_grouped[col].fillna(0).astype(int)

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

    # --- Build the figure with two panels ---
    fig = make_subplots(
        rows=1, cols=2,
        shared_yaxes=True,
        horizontal_spacing=0.12,
        column_widths=[0.6, 0.4],
        specs=[[{"type": "bar"}, {"type": "bar"}]],
        subplot_titles=(
            "Pathogenic Genotypes (%)",
            "Zygosity Composition (Pathogenic Rows)"
        )
    )

    trace_keys = list(trace_map.keys())

    # We'll add 4 traces per (pore, sex) view:
    #   - Left panel: 1 bar (affected % with CI)
    #   - Right panel: 3 stacked bars (HOM/HET/HEMI, normalized to 100%)
    for key in trace_keys:
        df_grouped = trace_map[key].copy()

        # Left panel (% affected)
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
                df_grouped['ci_upper'],
                df_grouped['zyg_HOM'],
                df_grouped['zyg_HET'],
                df_grouped['zyg_HEMI']
            ], axis=-1),
            hovertemplate=(
                "Population: %{y}<br>"
                "# Pathogenic Individuals: %{customdata[0]}<br>"
                "# Total Individuals: %{customdata[1]}<br>"
                "Pathogenic Genotype: %{customdata[2]}%<br>"
                "95% CI: [%{customdata[3]:.2f}%, %{customdata[4]:.2f}%]<br>"
                "<b>Pathogenic zygosity rows</b>: HOM %{customdata[5]}, "
                "HET %{customdata[6]}, HEMI %{customdata[7]}<extra></extra>"
            ),
            visible=False,
            row=1, col=1
        )

        # Right panel (composition among pathogenic rows)
        z_total = (df_grouped['zyg_HOM'] + df_grouped['zyg_HET'] + df_grouped['zyg_HEMI']).replace(0, np.nan)
        for comp_col, name in [('zyg_HOM','HOM'), ('zyg_HET','HET'), ('zyg_HEMI','HEMI')]:
            comp_pct = (df_grouped[comp_col] / z_total * 100).fillna(0.0)
            fig.add_bar(
                x=comp_pct,
                y=df_grouped['Cleaned population description'],
                orientation='h',
                name=name,
                hovertemplate=(
                    "Population: %{y}<br>"
                    f"{name} among pathogenic rows: %{x:.2f}%"
                    "<extra></extra>"
                ),
                visible=False,
                row=1, col=2
            )

    # Default selection: All Types • All Sexes (if present)
    default_key = ('All Types', 'All Sexes') if ('All Types', 'All Sexes') in trace_keys else trace_keys[0]
    default_block_start = trace_keys.index(default_key) * 4  # 4 traces per state
    for i in range(4):
        fig.data[default_block_start + i].visible = True

    # Dropdown to switch (pore, sex)
    buttons = []
    for s_idx, (pore, sex) in enumerate(trace_keys):
        vis = [False] * len(fig.data)
        block_start = s_idx * 4
        for k in range(4):
            vis[block_start + k] = True
        buttons.append(dict(
            label=f"{pore} • {sex}",
            method="update",
            args=[
                {"visible": vis},
                {"annotations": [
                    dict(
                        text="Total Individuals per Population: " + subtitle_map[(pore, sex)],
                        x=-0.06, y=1.15, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=12, color="black"),
                        align="left", xanchor="left"
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
        width=1100, height=520, margin=dict(t=125), showlegend=True,
        barmode='stack',  # right panel stacking
        title=dict(
            y=0.95,
            text=(
                f"<span style='font-size:18px; font-weight:bold'>Pathogenic Genotype Distribution</span><br>"
                f"<span style='font-size:12px'>{gene} - {disease}</span><br>"
                f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
            )
        ),
        updatemenus=[dict(
            active=trace_keys.index(default_key),
            buttons=buttons, direction="down",
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
        ]
    )

    # Axes
    fig.update_xaxes(title_text="Pathogenic Genotypes (%)", row=1, col=1, range=[0,100], dtick=10,
                     ticks='outside', showline=True, linecolor='black')
    fig.update_yaxes(title_text="", row=1, col=1, categoryorder='array', categoryarray=superpop_order,
                     showline=True, linecolor='black')

    fig.update_xaxes(title_text="Zygosity Composition (%, pathogenic rows)", row=1, col=2, range=[0,100], dtick=20,
                     ticks='outside', showline=True, linecolor='black')
    fig.update_yaxes(title_text="", row=1, col=2, categoryorder='array', categoryarray=superpop_order,
                     showline=True, linecolor='black')

    return fig

# --- Generate plots ---
made = 0
for gene in df_agg['Gene'].unique():
    for disease in df_agg[df_agg['Gene'] == gene]['Disease'].unique():
        sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]

        # Skip if no meaningful data
        if sub_df['total_count'].sum() == 0:
            continue

        fig = create_horizontal_bar_plot(sub_df.copy(), gene, disease, df)
        if not fig:
            continue

        safe_gene = re.sub(r'[\\/]', '_', gene)
        safe_disease = re.sub(r'[\\/]', '_', disease)

        html_path = os.path.join(OUTPUT_HTML_DIR, f"{safe_gene}_{safe_disease}_ancestry_plot.html")
        png_path  = os.path.join(OUTPUT_PNG_DIR, f"{safe_gene}_{safe_disease}_ancestry_plot.png")

        # Save files
        fig.write_html(html_path)
        fig.write_image(png_path, scale=2)   # scale=2 makes higher resolution

        if TEST_MODE:
            print(f"Previewing: {gene} / {disease}")
            fig.show()

            if SAVE_TEST_OUTPUTS:
                print(f"Saving: {gene} / {disease}")
                fig.write_html(html_path)

            made += 1
            if made >= TEST_LIMIT:
                break
        else:
            fig.write_html(html_path)

    if TEST_MODE and made >= TEST_LIMIT:
        break

# --- Finished ---
if TEST_MODE:
    print("--- Test mode ON: Test completed ---")
else:
    print("--- Done ---")