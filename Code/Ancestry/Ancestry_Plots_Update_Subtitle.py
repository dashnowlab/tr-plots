import pandas as pd
import re
import plotly.graph_objects as go
import os
import ast
import numpy as np
from statsmodels.stats import proportion

# Enable test mode to generate only one plot
TEST_MODE = True

# Load data
df = pd.read_csv('.../Matching Files/CSVs/83_loci_503_samples_withancestrycolumns.csv')

# Output directories
OUTPUT_DIR = '.../.../Plots/Ancestry_Plots/HTML/83_loci_503_samples_Subtitle_Test_2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Merge Pore information
pore_df = pd.read_csv('.../..../Datasets/Other/1KGP_ONT_500_Summary_Sample_ID_Pore.csv')
pore_df['Base Sample ID'] = pore_df['Sample_ID'].apply(lambda x: x.split('-')[0] if isinstance(x, str) else x)
df['Base Sample ID'] = df['Sample ID'].apply(lambda x: x.split('-')[0])
df = df.merge(pore_df[['Base Sample ID', 'Pore']], on='Base Sample ID', how='left')

# Use ancestry column if present
if 'SuperPop' in df.columns:
    df['Cleaned population description'] = df['SuperPop']
elif 'Cleaned population description' not in df.columns:
    raise KeyError("Missing both 'SuperPop' and 'Cleaned population description' columns.")

# Add 'Is pathogenic'
df['Is pathogenic'] = df.apply(
    lambda row: (
        row['Pathogenic min'] <= row['Repeat count'] <= row['Pathogenic max']
        if pd.notna(row['Pathogenic min']) and pd.notna(row['Pathogenic max']) and pd.notna(row['Repeat count'])
        else np.nan
    ),
    axis=1
)

# Normalize inheritance field (remove stringified lists)
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

# Aggregation
total_count = df.groupby(['Disease', 'Gene', 'Cleaned population description', 'Pore'])['Sample ID'].nunique().reset_index(name='total_count')

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

# Compute "All" population as sum of subgroups, per Gene / Disease / Pore
group_cols = ['Gene', 'Disease', 'Pore']
agg_cols = ['total_count'] + \
           [f'{inh}_affected_counts' for inh in ['ad', 'ar', 'xd', 'xr']] + \
           [f'{inh}_percentage' for inh in ['ad', 'ar', 'xd', 'xr']]

# Sum values and label as 'All' population
df_all_rows = (
    df_agg[df_agg['Cleaned population description'] != 'All']
    .groupby(group_cols, as_index=False)[agg_cols]
    .sum()
)
df_all_rows['Cleaned population description'] = 'All'

# Recalculate confidence intervals for 'All'
for inh in ['ad', 'ar', 'xd', 'xr']:
    lowers, uppers = [], []
    for x, n in zip(df_all_rows[f'{inh}_affected_counts'], df_all_rows['total_count']):
        lo, up = binomial_ci(x, n)
        lowers.append(lo)
        uppers.append(up)
    df_all_rows[f'{inh}_ci_lower'] = lowers
    df_all_rows[f'{inh}_ci_upper'] = uppers

# Append to main df_agg
df_agg = pd.concat([df_agg, df_all_rows], ignore_index=True)

def add_combined_row(group_df):
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
                if inh in ['ad', 'xd']:
                    min_alleles = 1
                else:
                    min_alleles = 2

                df_inh = df_sub[df_sub['Inheritance'] == inh_mode]
                df_path = df_inh[df_inh['Is pathogenic'] == True]
                grouped = df_path.groupby('Sample ID').size().reset_index(name='path_count')
                affected = grouped[grouped['path_count'] >= min_alleles]
                affected_count = affected['Sample ID'].nunique()

                combined[f'{inh}_affected_counts'] = affected_count
                combined[f'{inh}_percentage'] = (
                    affected_count / unique_samples * 100 if unique_samples > 0 else 0
                )
                lo, hi = binomial_ci(affected_count, unique_samples)
                combined[f'{inh}_ci_lower'] = lo
                combined[f'{inh}_ci_upper'] = hi

            result.append(combined)
    return pd.DataFrame(result)

pop_counts_strs = []

def create_horizontal_bar_plot(filtered_df, gene, disease, original_df):
    raw_inh = original_df[(original_df['Gene'] == gene) & (original_df['Disease'] == disease)]['Inheritance']
    inheritance = raw_inh.iloc[0] if not raw_inh.empty else ''
    if inheritance not in ['AD', 'AR', 'XD', 'XR']:
        print(f"Skipping unknown inheritance for {gene} / {disease}")
        return None

    inh = inheritance.lower()
    trace_map = {}
    trace_order = ['All Types', 'R9', 'R10']

    for pore_type in trace_order:
        if pore_type == 'All Types':
            df_pore = filtered_df.copy()
        else:
            df_pore = filtered_df[filtered_df['Pore'] == pore_type]

        if df_pore.empty:
            continue

        df_grouped = df_pore.groupby('Cleaned population description', as_index=False).agg({
            'total_count': 'sum',
            f'{inh}_affected_counts': 'sum'
        })

        # Move 'All' to the bottom of the DataFrame (which becomes top in the plot)
        df_grouped['sort_order'] = df_grouped['Cleaned population description'].apply(lambda x: 0 if x == 'All' else 1)
        df_grouped = df_grouped.sort_values(by='sort_order').drop(columns='sort_order')


        df_grouped['percentage'] = df_grouped[f'{inh}_affected_counts'] / df_grouped['total_count'] * 100
        df_grouped['percentage_display'] = df_grouped['percentage'].apply(lambda x: f"{int(x)}" if x == int(x) else f"{x:.2f}")
        df_grouped['affected'] = df_grouped[f'{inh}_affected_counts']
        ci_bounds = df_grouped.apply(lambda row: binomial_ci(row['affected'], row['total_count']), axis=1)
        df_grouped['ci_lower'] = [lo for lo, hi in ci_bounds]
        df_grouped['ci_upper'] = [hi for lo, hi in ci_bounds]
        df_grouped['error_plus'] = df_grouped['ci_upper'] - df_grouped['percentage']
        df_grouped['error_minus'] = df_grouped['percentage'] - df_grouped['ci_lower']
        trace_map[pore_type] = df_grouped

        # Build population counts subtitle for All Types
        pop_counts_str = ""
        if 'All Types' in trace_map:
            pop_df = trace_map['All Types'][['Cleaned population description', 'total_count']].copy()
            pop_df = pop_df.sort_values(by='Cleaned population description', key=lambda x: x != 'All')  # Move 'All' to top

            # Combine into a readable string with line breaks
            pop_lines = []
            line = ""
            for idx, row in pop_df.iterrows():
                part = f"{row['Cleaned population description']}: {int(row['total_count'])}"
                if len(line + part) > 100:
                    pop_lines.append(line.rstrip(', '))
                    line = ""
                line += part + ", "
            if line:
                pop_lines.append(line.rstrip(', '))

            pop_counts_str = "<br>".join([
                f"<span style='font-size:12px'>Total Individuals per Population: {line}</span>" if i == 0 else
                f"<span style='font-size:12px'>{line}</span>"
                for i, line in enumerate(pop_lines)
            ])
        else:
            pop_counts_str = "<span style='font-size:12px'>No population count data.</span>"
        
        pop_counts_strs.append(pop_counts_str)

    fig = go.Figure()
    trace_labels = []

    for key in trace_order:
        if key in trace_map:
            df_grouped = trace_map[key]
            fig.add_bar(
                x=df_grouped['percentage'],
                y=df_grouped['Cleaned population description'],
                orientation='h',
                error_x=dict(array=df_grouped['error_plus'], arrayminus=df_grouped['error_minus'], thickness=1),
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
            trace_labels.append(key)

    if 'All Types' in trace_labels:
        fig.data[trace_labels.index('All Types')].visible = True

    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=125),
        bargap=0.5,
        plot_bgcolor='white',
        showlegend=False,
        font=dict(color='black'),
        title=dict(
            y=0.95,
            text=(
                f"<span style='font-size:18px; font-weight:bold'>Pathogenic Genotype Distribution</span><br>"
                f"<span style='font-size:12px'>{gene} - {disease}</span><br>"
                f"{pop_counts_str}<br>"
                f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
            )
        ),
        updatemenus=[
            dict(
                active=trace_labels.index('All Types'),
                buttons=[
                    dict(
                        label=label,
                        method="update",
                        args=[
                            {"visible": [i == j for j in range(len(trace_labels))]},
                            {"annotations": [
                                dict(
                                    text=pop_counts_strs[i],
                                    x=0,
                                    y=1.03,
                                    xref="paper",
                                    yref="paper",
                                    showarrow=False,
                                    font=dict(size=12, color="black"),
                                    align="left",
                                    xanchor="left",
                                    name="pop_subtitle"
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
                    for i, label in enumerate(trace_labels)
                ],
                direction="down",
                showactive=True,
                x=1.15,
                xanchor="center",
                y=1.15,
                yanchor="top"
            )
        ],

        annotations=[
            dict(
                text=pop_counts_strs[trace_labels.index('All Types')],
                x=0,
                y=1.03,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=12, color="black"),
                align="left",
                xanchor="left",
                name="pop_subtitle"
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

# --- GENERATE PLOTS ---
printed = False
for gene in df_agg['Gene'].unique():
    for disease in df_agg[df_agg['Gene'] == gene]['Disease'].unique():
        sub_df = df_agg[(df_agg['Gene'] == gene) & (df_agg['Disease'] == disease)]

        # Skip if no meaningful data
        if sub_df['total_count'].sum() == 0:
            continue

        fig = create_horizontal_bar_plot(sub_df.copy(), gene, disease, df)
        if fig:
            safe_gene = re.sub(r'[\\/]', '_', gene)
            safe_disease = re.sub(r'[\\/]', '_', disease)

            print(f"Plot generated for: {gene} / {disease}")
            fig.show()  # For visual testing (or comment this if you only want to print)

            if not TEST_MODE:
                fig.write_html(os.path.join(OUTPUT_DIR, f"{safe_gene}_{safe_disease}_ancestry_plot.html"))
            else:
                printed = True
                break  # Only do one plot in test mode
    if TEST_MODE and printed:
        break

print("--- Test completed ---" if TEST_MODE else "--- Saved all plots ---")