import pandas as pd
import re
import plotly.graph_objects as go
import os
import ast
import numpy as np
from statsmodels.stats import proportion

TEST_MODE = True

# Load data
df = pd.read_excel('/Users/annelisethorn/Documents/GitHub/tr-plots/Code/Matching Files/Excels/83_loci_503_samples_withancestrycolumns.xlsx')

# Output directories
OUTPUT_DIR = '/Users/annelisethorn/Documents/GitHub/tr-plots/Plots/Ancestry_Plots/83_loci_503_samples_subtitletest'
OUTPUT_DIR2 = '/Users/annelisethorn/Documents/GitHub/tr-plots/Plots/Ancestry_Plots/83_loci_503_samples_subtitletest'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR2, exist_ok=True)

# Merge Pore info
pore_df = pd.read_csv('/Users/annelisethorn/Documents/GitHub/tr-plots/Datasets/Other/1KGP_ONT_500_Summary_Sample_ID_Pore.csv')
pore_df['Base Sample ID'] = pore_df['Sample_ID'].apply(lambda x: x.split('-')[0] if isinstance(x, str) else x)
df['Base Sample ID'] = df['Sample ID'].apply(lambda x: x.split('-')[0])
df = df.merge(pore_df[['Base Sample ID', 'Pore']], on='Base Sample ID', how='left')

# Use SuperPop
if 'SuperPop' in df.columns:
    df['Cleaned population description'] = df['SuperPop']
elif 'Cleaned population description' not in df.columns:
    raise KeyError("Missing both 'SuperPop' and 'Cleaned population description' columns.")

# Is pathogenic
df['Is pathogenic'] = df.apply(
    lambda row: (
        row['Pathogenic min'] <= row['Repeat count'] <= row['Pathogenic max']
        if pd.notna(row['Pathogenic min']) and pd.notna(row['Pathogenic max']) and pd.notna(row['Repeat count'])
        else np.nan
    ),
    axis=1
)

# Clean inheritance
def clean_inheritance(x):
    if isinstance(x, str) and x.startswith('['):
        try:
            val = ast.literal_eval(x)
            return val[0] if val else np.nan
        except Exception:
            return np.nan
    return x

mask = df['Inheritance'].astype(str).str.startswith('[')
df.loc[mask, 'Inheritance'] = df.loc[mask, 'Inheritance'].apply(clean_inheritance)

# Count total
total_count = df.groupby(['Disease', 'Gene', 'Cleaned population description', 'Pore'])['Sample ID'].nunique().reset_index(name='total_count')

# Count affected
def count_affected(df, inh_mode, min_alleles):
    df_inh = df[df['Inheritance'] == inh_mode]
    df_path = df_inh[df_inh['Is pathogenic'] == True]
    grouped = df_path.groupby(['Disease', 'Gene', 'Cleaned population description', 'Pore', 'Sample ID']).size().reset_index(name='path_count')
    affected = grouped[grouped['path_count'] >= min_alleles]
    return affected.groupby(['Disease', 'Gene', 'Cleaned population description', 'Pore'])['Sample ID'].nunique().reset_index(name=f'{inh_mode.lower()}_affected_counts')

ad_affected = count_affected(df, 'AD', 1)
ar_affected = count_affected(df, 'AR', 2)
xd_affected = count_affected(df, 'XD', 1)
xr_affected = count_affected(df, 'XR', 2)

df_agg = total_count.copy()
for affected_df in [ad_affected, ar_affected, xd_affected, xr_affected]:
    df_agg = df_agg.merge(affected_df, on=['Disease', 'Gene', 'Cleaned population description', 'Pore'], how='left')

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

# "All" population
group_cols = ['Gene', 'Disease', 'Pore']
agg_cols = ['total_count'] + \
           [f'{inh}_affected_counts' for inh in ['ad', 'ar', 'xd', 'xr']] + \
           [f'{inh}_percentage' for inh in ['ad', 'ar', 'xd', 'xr']]
df_all_rows = df_agg[df_agg['Cleaned population description'] != 'All'].groupby(group_cols, as_index=False)[agg_cols].sum()
df_all_rows['Cleaned population description'] = 'All'
for inh in ['ad', 'ar', 'xd', 'xr']:
    lowers, uppers = [], []
    for x, n in zip(df_all_rows[f'{inh}_affected_counts'], df_all_rows['total_count']):
        lo, up = binomial_ci(x, n)
        lowers.append(lo)
        uppers.append(up)
    df_all_rows[f'{inh}_ci_lower'] = lowers
    df_all_rows[f'{inh}_ci_upper'] = uppers

df_agg = pd.concat([df_agg, df_all_rows], ignore_index=True)

# Plot function
def create_horizontal_bar_plot(filtered_df, gene, disease, original_df):
    inheritance_raw = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]['Inheritance']
    inheritance = inheritance_raw.iloc[0] if not inheritance_raw.empty else ''
    if inheritance not in ['AD', 'AR', 'XD', 'XR']:
        return None
    inh = inheritance.lower()

    trace_order = ['All Types', 'R9', 'R10']
    trace_map = {}
    subtitles = []

    for pore_type in trace_order:
        df_sub = filtered_df.copy() if pore_type == 'All Types' else filtered_df[filtered_df['Pore'] == pore_type]
        if df_sub.empty:
            continue
        df_grouped = df_sub.groupby('Cleaned population description', as_index=False).agg({
            'total_count': 'sum',
            f'{inh}_affected_counts': 'sum'
        })
        df_grouped['percentage'] = df_grouped[f'{inh}_affected_counts'] / df_grouped['total_count'] * 100
        df_grouped['percentage_display'] = df_grouped['percentage'].apply(lambda x: f"{int(x)}" if x == int(x) else f"{x:.2f}")
        df_grouped['affected'] = df_grouped[f'{inh}_affected_counts']
        ci_bounds = df_grouped.apply(lambda row: binomial_ci(row['affected'], row['total_count']), axis=1)
        df_grouped['ci_lower'] = [lo for lo, hi in ci_bounds]
        df_grouped['ci_upper'] = [hi for lo, hi in ci_bounds]
        df_grouped['error_plus'] = df_grouped['ci_upper'] - df_grouped['percentage']
        df_grouped['error_minus'] = df_grouped['percentage'] - df_grouped['ci_lower']
        trace_map[pore_type] = df_grouped

        # Subtitle
        pop_lines = []
        line = ""
        for _, row in df_grouped.iterrows():
            part = f"{row['Cleaned population description']}: {int(row['total_count'])}"
            if len(line + part) > 100:
                pop_lines.append(line.rstrip(', '))
                line = ""
            line += part + ", "
        if line:
            pop_lines.append(line.rstrip(', '))
        subtitle = "<br>".join([f"<span style='font-size:12px'>{l}</span>" for l in pop_lines])
        subtitles.append(subtitle)

    fig = go.Figure()
    for i, pore_type in enumerate(trace_order):
        if pore_type not in trace_map:
            continue
        df_grouped = trace_map[pore_type]
        fig.add_bar(
            x=df_grouped['percentage'],
            y=df_grouped['Cleaned population description'],
            orientation='h',
            error_x=dict(array=df_grouped['error_plus'], arrayminus=df_grouped['error_minus']),
            customdata=np.stack([
                df_grouped['affected'],
                df_grouped['total_count'],
                df_grouped['percentage_display'],
                df_grouped['ci_lower'],
                df_grouped['ci_upper']
            ], axis=-1),
            hovertemplate="Population: %{y}<br># Pathogenic: %{customdata[0]}<br># Total: %{customdata[1]}<br>% Pathogenic: %{customdata[2]}%<br>95% CI: [%{customdata[3]:.2f}%, %{customdata[4]:.2f}%]<extra></extra>",
            visible=(i == 0)
        )

    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=125),
        bargap=0.5,
        plot_bgcolor='white',
        showlegend=False,
        title=dict(
            text=f"<span style='font-size:18px; font-weight:bold'>Pathogenic Genotype Distribution</span><br>"
                 f"<span style='font-size:12px'>{gene} - {disease}</span><br>"
                 f"<span style='font-size:12px'>Inheritance: {inheritance}</span>",
            y=0.95
        ),
        updatemenus=[
            dict(
                active=0,
                buttons=[
                    dict(
                        label=pore,
                        method="update",
                        args=[
                            {"visible": [j == i for j in range(len(trace_map))]},
                            {"annotations": [
                                dict(
                                    text=subtitles[i],
                                    x=-0.054,
                                    y=1.27,
                                    xref="paper",
                                    yref="paper",
                                    showarrow=False,
                                    font=dict(size=12),
                                    align="left"
                                ),
                                dict(
                                    text="Pore Type:",
                                    x=1.15,
                                    y=1.23,
                                    xref="paper",
                                    yref="paper",
                                    showarrow=False,
                                    font=dict(size=12, color="black"),
                                    xanchor="center"
                                )
                            ]}
                        ]
                    )
                    for i, pore in enumerate(trace_map)
                ],
                direction="down",
                showactive=True,
                x=1.1,
                xanchor="center",
                y=1.27,
                yanchor="top"
            )
        ],
        annotations=[
            dict(
                text=subtitles[0],
                x=-0.054,
                y=1.27,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=12),
                align="left"
            )
        ],
        xaxis=dict(
            range=[0, 100],
            tickmode='linear',
            tick0=0,
            dtick=10,
            ticks='outside',
            showline=True,
            linecolor='black',
            linewidth=1,
            zeroline=False,
            tickfont=dict(color='black'),
            titlefont=dict(color='black'),
            title="Pathogenic Genotypes (%)"
        ),
        yaxis=dict(
            title="",
            showline=True,
            linecolor='black'
        )
    )
    return fig

# --- Generate ---
printed = False
for gene in df_agg['Gene'].unique():
    for disease in df_agg[df_agg['Gene'] == gene]['Disease'].unique():
        sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]
        if sub_df['total_count'].sum() == 0:
            continue
        fig = create_horizontal_bar_plot(sub_df.copy(), gene, disease, df)
        if fig:
            safe_gene = re.sub(r'[\\/]', '_', gene)
            safe_disease = re.sub(r'[\\/]', '_', disease)
            if not TEST_MODE:
                fig.write_html(os.path.join(OUTPUT_DIR2, f"{safe_gene}_{safe_disease}_ancestry_plot.html"))
            else:
                fig.show()
                printed = True
                break
    if TEST_MODE and printed:
        break

print("--- Test completed ---" if TEST_MODE else "--- Saved all plots ---")
