"""
---------------------------------------------
# Script: Allele Length Violin + Beeswarm Plot Generator
# Purpose:
#   - Reads a CSV of allele lengths with ancestry and metadata
#   - Flags pathogenic individuals using per-locus thresholds
#   - Aggregates counts and computes % pathogenic by population
#   - Builds horizontal violin + beeswarm plots of allele length by population
#   - Saves plots as PNG and HTML 
#   - In test mode, previews a single plot for quick checks
---------------------------------------------
"""

import os
import re
import ast
import numpy as np
import pandas as pd
import plotly.express as px

# --- TEST MODE ---
TEST_MODE = False              
TEST_LIMIT = 1               
SAVE_TEST_OUTPUTS = True      

# --- VISUAL OPTIONS ---
COLOR_BY_PATHOGENIC = True    
POINT_SIZE = 6                
POINT_OPACITY = 0.75          
VIOLIN_SHOW_BOX = True        
JITTER = 0.35                 

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

# --- Load data ---
df = pd.read_csv(DATA_PATH)

# --- Flag pathogenic individuals ---
df['Is pathogenic'] = (
    (df['Allele length'] >= df['Pathogenic min']) &
    (df['Allele length'] <= df['Pathogenic max'])
)

# --- Map to categorical labels (Benign, Pathogenic, Unknown) ---
def map_status(val):
    if pd.isna(val):
        return "Unknown"
    return "Pathogenic" if val else "Benign"

df['Pathogenic status'] = df['Is pathogenic'].apply(map_status)
df['Pathogenic status'] = pd.Categorical(
    df['Pathogenic status'],
    categories=["Benign", "Pathogenic", "Unknown"],
    ordered=True
)

# --- Aggregate setup (kept for completeness/compatibility) ---
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

# Compute affected counts
ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

# Merge into summary
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

# --- Plotting helper ---
def create_violin_beeswarm(original_df, gene, disease, color_by_pathogenic=True):
    subset = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)].copy()
    if subset.empty:
        return None

    raw_inh = subset['Inheritance']
    try:
        inh_list = ast.literal_eval(raw_inh.iloc[0]) if not raw_inh.empty else []
        inheritance = inh_list[0] if inh_list else ''
    except Exception:
        inheritance = str(raw_inh.iloc[0]) if not raw_inh.empty else ''
    if not inheritance:
        inheritance = "Unknown"
    if inheritance not in ['AD', 'AR', 'XD', 'XR', 'Unknown']:
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

    # Wrap long population summary like in the boxplot generator
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

    fig = px.violin(
        subset,
        x='Allele length',
        y='SuperPop',
        orientation='h',
        box=VIOLIN_SHOW_BOX,
        points=False,
        color_discrete_sequence=["lightblue"],
        title=title_html,
    )

    # Use Pathogenic status for coloring
    color_arg = 'Pathogenic status' if color_by_pathogenic else None

    strip_fig = px.strip(
        subset,
        x='Allele length',
        y='SuperPop',
        orientation='h',
        color=color_arg,
        category_orders={'Pathogenic status': ["Benign", "Pathogenic", "Unknown"]},
        hover_data=['Sample ID', 'Allele length', 'SuperPop', 'Gene', 'Disease']
    )

    for tr in strip_fig.data:
        x_vals = np.array(tr.x, dtype=float)
        if JITTER and JITTER > 0:
            noise = np.random.uniform(-JITTER, JITTER, size=x_vals.shape)
            x_vals = x_vals + noise
        tr.x = x_vals
        tr.update(
            marker=dict(size=POINT_SIZE, opacity=POINT_OPACITY, line=dict(width=0)),
            showlegend=color_by_pathogenic
        )
        fig.add_trace(tr)

    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=150, r=20, b=40, l=80),
        xaxis=dict(title="Allele Length", ticks='outside', showline=True, linecolor='black', zeroline=False),
        yaxis=dict(title="", ticks='outside', showline=True, linecolor='black'),
        plot_bgcolor='white',
        font=dict(color='black'),
        title=dict(y=0.95)
    )
    if not color_by_pathogenic:
        fig.update_layout(showlegend=False)

    return fig

# --- Generate plots ---
printed = 0
pairs = df[['Gene', 'Disease']].drop_duplicates().itertuples(index=False, name=None)
for gene, disease in pairs:
    fig = create_violin_beeswarm(df, gene, disease, color_by_pathogenic=COLOR_BY_PATHOGENIC)
    if not fig:
        continue

    safe_gene = re.sub(r'[\\/]', '_', gene)
    safe_disease = re.sub(r'[\\/]', '_', disease)
    png_path = os.path.join(OUTPUT_DIR_PNG, f"{safe_gene}_{safe_disease}_allele_length_violin_beeswarm.png")
    html_path = os.path.join(OUTPUT_DIR_HTML, f"{safe_gene}_{safe_disease}_allele_length_violin_beeswarm.html")

    if TEST_MODE:
        print(f"Previewing: {gene} / {disease}")
        fig.show()
        if SAVE_TEST_OUTPUTS:
            fig.write_html(html_path)
            try:
                fig.write_image(png_path, format="png", scale=2)
            except Exception as e:
                print(f"[PNG export skipped] {e}")
        printed += 1
        if printed >= TEST_LIMIT:
            break
    else:
        fig.write_html(html_path)
        try:
            fig.write_image(png_path, format="png", scale=2)
        except Exception as e:
            print(f"[PNG export skipped] {e}")

# --- Finished ---
if TEST_MODE:
    print("--- Test mode ON: Test completed ---")
else:
    print("--- Done ---")
