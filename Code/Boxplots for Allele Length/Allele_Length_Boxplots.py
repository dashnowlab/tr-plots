import pandas as pd
import re
import plotly.express as px
import os
import ast

# Load data
# df = pd.read_csv('/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_503_samples_withancestrycolumns.csv')
df = pd.read_excel('/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_503_samples_withancestrycolumns_2.xlsx')


# Compute 'Is pathogenic' from allele length and pathogenic thresholds
df['Is pathogenic'] = (
    (df['Allele length'] >= df['Pathogenic min']) &
    (df['Allele length'] <= df['Pathogenic max'])
)

# Output paths
OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Allele_Length_Boxplots/PNG/Switched_Axes/83_loci_503_samples_corrected"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_DIR2 = "/Users/annelisethorn/Documents/Anschutz/Plots/Allele_Length_Boxplots/HTML/Switched_Axes/83_loci_503_samples_corrected"
os.makedirs(OUTPUT_DIR2, exist_ok=True)

# Group by Disease, Gene, SuperPop, Allele length
agg_df = df.groupby(['Disease', 'Gene', 'SuperPop', 'Allele length'])

# Count total unique samples
total_count = df.groupby(['Disease', 'Gene', 'SuperPop', 'Allele length'])['Sample ID'].nunique().reset_index(name='total_count')

def count_affected(df, inh_mode, min_alleles):
    df_inh = df[df['Inheritance'].astype(str) == f"['{inh_mode}']"]
    df_path = df_inh[df_inh['Is pathogenic'] == True]
    grouped = df_path.groupby(['Disease', 'Gene', 'SuperPop', 'Sample ID']).size().reset_index(name='path_count')
    affected = grouped[grouped['path_count'] >= min_alleles]
    affected_counts = affected.groupby(['Disease', 'Gene', 'SuperPop'])['Sample ID'].nunique().reset_index(name=f'{inh_mode.lower()}_affected_counts')
    return affected_counts

# Affected counts
ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

# Merge and fill NaNs
df_agg = total_count.copy()
for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
    df_agg = df_agg.merge(affected_df, on=['Disease', 'Gene', 'SuperPop'], how='left')

for inh in ['ad', 'ar', 'xd', 'xr']:
    df_agg[f'{inh}_affected_counts'] = df_agg[f'{inh}_affected_counts'].fillna(0).astype(int)
    df_agg[f'{inh}_percentage'] = df_agg[f'{inh}_affected_counts'] / df_agg['total_count'] * 100

def create_boxplot(filtered_df, gene, disease, original_df):
    if filtered_df.empty:
        return None

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

    relevant_samples = df[(df['Gene'] == gene) & (df['Disease'] == disease)]

    unique_pops = (
        relevant_samples
        .groupby('SuperPop')['Sample ID']
        .nunique()
        .reset_index(name='pop_total')
        .sort_values(by='SuperPop')
    )

    pop_desc = ', '.join(f"{row['SuperPop']}: {row['pop_total']}" for _, row in unique_pops.iterrows())

    max_line_len = 107
    pop_lines = []
    line = ""
    for segment in pop_desc.split(", "):
        if len(line) + len(segment) > max_line_len:
            pop_lines.append(line.strip(", "))
            line = ""
        line += segment + ", "
    if line:
        pop_lines.append(line.strip(", "))

    pop_size_html = "<br>".join([
        f"<span style='font-size:12px'>Total Individuals per Population: {line}</span>" if i == 0
        else f"<span style='font-size:12px'>{line}</span>"
        for i, line in enumerate(pop_lines)
    ])

    # Wrap the disease name if too long
    max_line_len_disease = 80
    disease_line = ""
    disease_lines = []

    for word in f"{gene} - {disease}".split(" "):
        if len(disease_line) + len(word) + 1 > max_line_len_disease:
            disease_lines.append(disease_line.strip())
            disease_line = ""
        disease_line += word + " "
    if disease_line:
        disease_lines.append(disease_line.strip())

    disease_html = "<br>".join(f"<span style='font-size:12px'>{line}</span>" for line in disease_lines)

    title_html = (
        f"<span style='font-size:18px; font-weight:bold'>Allele Lengths per Population</span><br>"
        f"{disease_html}<br>{pop_size_html}<br>"
        f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
    )

    fig = px.box(
        filtered_df,
        x='Allele length',
        y='SuperPop',
        color_discrete_sequence=["blue"],
        title=title_html,
    )

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

# Save plots
for gene in df_agg['Gene'].unique():
    diseases = df_agg[df_agg['Gene'] == gene]['Disease'].unique()
    for disease in diseases:
        sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]
        fig = create_boxplot(sub_df.copy(), gene, disease, df)
        if fig:
            safe_gene = re.sub(r'[\\/]', '_', gene)
            safe_disease = re.sub(r'[\\/]', '_', disease)
            fig.write_image(os.path.join(OUTPUT_DIR, f"{safe_gene}_{safe_disease}_allele_length_boxplot.png"), format="png", scale=2)
            fig.write_html(os.path.join(OUTPUT_DIR2, f"{safe_gene}_{safe_disease}_allele_length_boxplot.html"))

print("--- Saved plots ---")
