"""
---------------------------------------------
 Script: Allele Length Boxplot Generator
 Purpose:
   - Reads a CSV of allele lengths with ancestry and metadata
   - Flags pathogenic individuals using per-locus thresholds
   - Aggregates counts and computes % pathogenic by population
   - Builds horizontal boxplots of allele length by population
   - Saves plots as PNG and HTML 
   - In test mode, previews a single plot for quick checks
---------------------------------------------
"""

import pandas as pd
import re
import plotly.express as px
import os
import ast

# --- TEST MODE ---
TEST_MODE = False              # Toggle this flag for quick testing (preview only)
TEST_LIMIT = 1                # How many (gene, disease) plots to generate in test mode
SAVE_TEST_OUTPUTS = True      # Toggle saving plots when in test mode

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots/"

DATA_PATH = f"{BASE_DIR}Code/Matching Files/CSVs/83_loci_503_samples_withancestrycolumns.csv"
OUTPUT_DIR = f"{BASE_DIR}Results/Plots/Allele_Length_Boxplots"

# Normal mode: save into PNG and HTML subfolders
OUTPUT_DIR_PNG = os.path.join(OUTPUT_DIR, "PNG")
OUTPUT_DIR_HTML = os.path.join(OUTPUT_DIR, "HTML")
os.makedirs(OUTPUT_DIR_PNG, exist_ok=True)
os.makedirs(OUTPUT_DIR_HTML, exist_ok=True)

# If test mode: override to a single test_outputs folder
if TEST_MODE:
    OUTPUT_DIR = os.path.join(OUTPUT_DIR, "test_outputs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save both PNG + HTML into this one folder
    OUTPUT_DIR_PNG = OUTPUT_DIR
    OUTPUT_DIR_HTML = OUTPUT_DIR

# --- Load data ---
df = pd.read_csv(DATA_PATH)

# --- Flag pathogenic individuals using thresholds ---
# True if allele length is within [Pathogenic min, Pathogenic max]
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

# Helper: count affected individuals per inheritance model
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

# --- Plotting helper ---
def create_boxplot(filtered_df, gene, disease, original_df):
    """
    Build a horizontal boxplot:
      - X: Allele length
      - Y: SuperPop
    Title includes gene/disease, total counts per population, and inheritance mode.
    """
    if filtered_df.empty:
        return None

    # Get inheritance mode for title; data sometimes stores as stringified list
    raw_inh = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]['Inheritance']
    try:
        inh_list = ast.literal_eval(raw_inh.iloc[0]) if not raw_inh.empty else []
        inheritance = inh_list[0] if inh_list else ''
    except Exception:
        inheritance = str(raw_inh.iloc[0]) if not raw_inh.empty else ''
    if not inheritance:
        inheritance = "Unknown"
    if inheritance not in ['AD', 'AR', 'XD', 'XR', 'Unknown']:
        print(f"Skipping unknown inheritance: {gene}")
        return None

    # Totals per population (for subtitle)
    relevant_samples = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]
    unique_pops = (
        relevant_samples
        .groupby('SuperPop')['Sample ID']
        .nunique()
        .reset_index(name='pop_total')
        .sort_values(by='SuperPop')
    )
    pop_desc = ', '.join(f"{row['SuperPop']}: {row['pop_total']}" for _, row in unique_pops.iterrows())

    # Wrap long population summary
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

    # Wrap "GENE - DISEASE" if long
    max_line_len_disease = 80
    disease_line, disease_lines = "", []
    for word in f"{gene} - {disease}".split(" "):
        if len(disease_line) + len(word) + 1 > max_line_len_disease:
            disease_lines.append(disease_line.strip())
            disease_line = ""
        disease_line += word + " "
    if disease_line:
        disease_lines.append(disease_line.strip())
    disease_html = "<br>".join(f"<span style='font-size:12px'>{ln}</span>" for ln in disease_lines)

    # Compose title
    title_html = (
        f"<span style='font-size:18px; font-weight:bold'>Allele Lengths per Population</span><br>"
        f"{disease_html}<br>{pop_size_html}<br>"
        f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
    )

    # Build the boxplot
    fig = px.box(
        filtered_df,
        x='Allele length',
        y='SuperPop',
        color_discrete_sequence=["blue"],
        title=title_html,
    )

    # Style
    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=150),
        bargap=0.5,
        xaxis=dict(title="Allele Length", ticks='outside', showline=True, linecolor='black'),
        yaxis=dict(title="", ticks='outside', showline=True, linecolor='black'),
        plot_bgcolor='white',
        showlegend=False,
        font=dict(color='black'),
        title=dict(y=0.95)
    )
    return fig

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
        png_path = os.path.join(OUTPUT_DIR_PNG, f"{safe_gene}_{safe_disease}_allele_length_boxplot.png")
        html_path = os.path.join(OUTPUT_DIR_HTML, f"{safe_gene}_{safe_disease}_allele_length_boxplot.html")

        if TEST_MODE:
            # Preview in test mode
            print(f"Previewing: {gene} / {disease}")
            fig.show()

            # Save test outputs (controlled by SAVE_TEST_OUTPUTS at top of script)
            if SAVE_TEST_OUTPUTS:
                print(f"Saving: {gene} / {disease}")
                fig.write_html(html_path)
                fig.write_image(png_path, format="png", scale=2)

            printed += 1
            if printed >= TEST_LIMIT:
                break
        else:
            # Non-test mode: save quietly (no printing/showing)
            fig.write_html(html_path)
            fig.write_image(png_path, format="png", scale=2)
    if TEST_MODE and printed >= TEST_LIMIT:
        break

# --- Finished ---
if TEST_MODE:
    print("--- Test mode ON: Test completed ---")
else:
    print("--- Done ---")
