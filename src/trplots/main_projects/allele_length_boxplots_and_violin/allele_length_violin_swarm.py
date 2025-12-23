import re
import ast
import pandas as pd
import plotly.express as px
import plotly.io as pio
from pathlib import Path
import argparse

# Pull paths from the shared config
pio.renderers.default = "browser"

# Match boxplot script path logic and behavior
BASE_DIR = Path("/Users/annelisethorn/Documents/github/tr-plots")

# --- Defaults: input spreadsheet and sheet name ---
DATA_PATH = BASE_DIR / "data" / "other_data" / "allele_spreadsheet.xlsx"
SHEET_NAME = "Integrated Alleles"

# --- Output directories ---
OUTPUT_ROOT = BASE_DIR / "results" / "plots" / "allele_length_violin_swarm_plots"
OUTPUT_DIR_PNG  = OUTPUT_ROOT / "png"
OUTPUT_DIR_HTML = OUTPUT_ROOT / "html"
TEST_OUTPUT_DIR = OUTPUT_ROOT / "test_outputs"

# Ensure the standard output dirs exist up front
OUTPUT_DIR_PNG.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR_HTML.mkdir(parents=True, exist_ok=True)

# --- TEST MODE ---
TEST_MODE = True
TEST_LIMIT = 2
SAVE_TEST_OUTPUTS = True

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130
PNG_SCALE = 2


# --- POPULATION PALETTE (1kG-style) ---
POP_COLOR = {
    'EUR': '#1f77b4',   # blue
    'EAS': '#2ca02c',   # green
    'SAS': '#9467bd',   # purple
    'AMR': '#d62728',   # red
    'AFR': '#ff7f0e',   # orange/yellow
    'All': '#7f7f7f',   # gray for aggregates
    'Unknown': '#ff99cc',  # fallback for unknowns (optional)
}
SUPERPOP_ORDER = ['All', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'Unknown']

# --- Helper functions (Performance optimizations retained) ---

def _wrap_to_lines(s: str, max_len: int = 110):
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

def create_violin_swarm(locus_df, locus_key_data):
    """Build a horizontal violin plot with overlaid points (swarm) using pre-calculated data."""
    gene, disease = locus_key_data['key']
    counts = locus_key_data['counts']
    inheritance = locus_key_data['inheritance']

    if locus_df.empty:
        return None

    ordered_categories = counts.index.tolist()
    fdf = locus_df.copy()
    # Ensure SuperPop is a categorical type for correct ordering in Plotly
    fdf['SuperPop'] = pd.Categorical(fdf['SuperPop'], categories=ordered_categories, ordered=True)

    # Prepare population description from pre-calculated counts
    pop_desc = ', '.join(f"{pop}: {counts.loc[pop]}" for pop in ordered_categories if pop in counts.index)
    pop_lines = _wrap_to_lines(pop_desc, max_len=110)

    # ---- Build the violin + swarm (no built-in title) ----
    fig = px.violin(
        fdf,
        x='Allele length',
        y='SuperPop',
        color='SuperPop',
        category_orders={"SuperPop": ordered_categories},
        color_discrete_map=POP_COLOR,
        box=False,
        points='all',
        title=None,
        orientation='h',
    )

    # Make the point cloud look denser (swarm-ish)
    fig.update_traces(
        jitter=0.35,
        marker=dict(size=4, line=dict(width=0)),
        meanline_visible=True
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
                     ([f"Total Individuals per Population: {pop_lines[0]}"] if pop_lines else []) + \
                     (pop_lines[1:] if len(pop_lines) > 1 else []) + \
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
                y=y0 - 0.04*i, yref="paper", yanchor="top",
                showarrow=False, align="left"
            )
        )
    fig.update_layout(annotations=annos)

    return fig

def parse_args():
    """Parse CLI args to override test mode, limits, data path and outputs."""
    p = argparse.ArgumentParser(description="Allele length violin+swarm generator")
    p.add_argument("--test", dest="test", action="store_true", help="Enable test mode")
    p.add_argument("--no-test", dest="test", action="store_false",
                   help="Disable test mode", default=not TEST_MODE)
    p.add_argument("--limit", type=int, default=TEST_LIMIT,
                   help="Set the test limit for number of plots")
    p.add_argument("--data-path", type=str, default=str(DATA_PATH),
                   help="Path to allele_spreadsheet.xlsx (Excel)")
    p.add_argument("--sheet", type=str, default=SHEET_NAME,
                   help="Excel sheet name (default: 'Integrated Alleles')")
    p.add_argument("--output-dir", type=str, default=None,
                   help="If provided, write all outputs (PNG/HTML) here instead of the default png/html dirs")
    p.set_defaults(test=TEST_MODE)
    return p.parse_args()

# --- Main function: move top-level processing here ---
def main():
    args = parse_args()

    # --- Test mode overrides ---
    global TEST_MODE, TEST_LIMIT
    TEST_MODE = args.test
    TEST_LIMIT = args.limit

    # Default output destinations
    output_dir_png = OUTPUT_DIR_PNG
    output_dir_html = OUTPUT_DIR_HTML
    if args.output_dir:
        outdir = Path(args.output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        output_dir_png = outdir
        output_dir_html = outdir


    # Data Excel override (pathlib.Path compatible)
    data_xlsx = Path(args.data_path)
    sheet_name = args.sheet

    # Load data from the integrated spreadsheet
    print(f"Loading data from: {data_xlsx} [sheet: {sheet_name}]")
    needed_cols = [
        'Gene', 'Disease', 'SuperPop', 'Allele length', 'Inheritance',
        'Pathogenic min', 'Pathogenic max', 'Sample ID', 'Sample ID Cleaned'
    ]
    df = pd.read_excel(
        data_xlsx,
        sheet_name=sheet_name,
        usecols=lambda c: True if c in needed_cols else False,
        engine="openpyxl"
    )

    # Prefer cleaned sample IDs if available
    if 'Sample ID Cleaned' in df.columns:
        df.rename(columns={'Sample ID Cleaned': 'SampleIDForCounts'}, inplace=True)
    elif 'Sample ID' in df.columns:
        df.rename(columns={'Sample ID': 'SampleIDForCounts'}, inplace=True)
    else:
        df['SampleIDForCounts'] = pd.NA

    # Coerce numeric fields and clean lengths
    df['Allele length'] = pd.to_numeric(df.get('Allele length'), errors='coerce')
    df.loc[~pd.Series(pd.notna(df['Allele length'])), 'Allele length'] = pd.NA
    df.loc[df['Allele length'] <= 0, 'Allele length'] = pd.NA
    for col in ['Pathogenic min', 'Pathogenic max']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # --- Add 'All' population (Vectorized) ---
    df_all = df.assign(SuperPop="All")
    df = pd.concat([df, df_all], ignore_index=True)
    
    # Assign 'SuperPop' as a categorical type for consistent ordering/grouping
    df['SuperPop'] = pd.Categorical(df['SuperPop'], categories=SUPERPOP_ORDER, ordered=True)

    # --- Flag pathogenic individuals using thresholds (Vectorized) ---
    if {'Pathogenic min', 'Pathogenic max'}.issubset(df.columns):
        df['Is pathogenic'] = (
            (df['Allele length'] >= df['Pathogenic min']) &
            (df['Allele length'] <= df['Pathogenic max'])
        )
    else:
        df['Is pathogenic'] = pd.NA

    # --- Pre-calculate all data outside the plot loop (Performance optimization) ---

    # 1. Total unique individuals per Locus/SuperPop (Required for header)
    locus_pop_counts = (
        df.groupby(['Gene', 'Disease', 'SuperPop'], observed=True)['SampleIDForCounts']
          .nunique(dropna=True)
          .rename('total_individuals')
    )
    
    # 2. Inheritance Mode (Required for header)
    inheritance_map = {}
    for (gene, disease), group in df.groupby(['Gene', 'Disease']):
        # Find the first non-null 'Inheritance' value
        raw_inh = group['Inheritance'].dropna().iloc[0] if not group['Inheritance'].dropna().empty else ''
        try:
            # Safely evaluate string to list and get the first element
            inh_list = ast.literal_eval(raw_inh) if isinstance(raw_inh, str) and raw_inh.startswith('[') else [raw_inh]
            inheritance = str(inh_list[0]) if inh_list and str(inh_list[0]) else ''
        except Exception:
            inheritance = str(raw_inh) if raw_inh else ''
            
        inheritance = inheritance if inheritance in ['AD', 'AR', 'XD', 'XR'] else 'Unknown'
        inheritance_map[(gene, disease)] = inheritance

    # 3. Locus-specific dataframes (for plotting)
    locus_data_map = {}
    for (gene, disease), group in df.groupby(['Gene', 'Disease']):
        # Filter for the plotting columns
        plot_df = group[['Allele length', 'SuperPop']].copy()
        
        # Get the individual counts for the current locus
        counts = locus_pop_counts.loc[(gene, disease)]
        
        # Drop 'Unknown' and reindex to get the sorted list of present SuperPops
        present_pops = counts.drop('Unknown', errors='ignore').index.tolist()
        ordered_categories = [c for c in SUPERPOP_ORDER if c in present_pops]
        
        # Re-order counts based on SUPERPOP_ORDER for the final output
        counts_reordered = counts.reindex(ordered_categories).dropna().astype(int)

        locus_data_map[(gene, disease)] = {
            'df': plot_df,
            'counts': counts_reordered,
            'inheritance': inheritance_map[(gene, disease)],
            'key': (gene, disease)
        }

    # --- Generate plots using the pre-calculated map ---
    printed = 0
    locus_keys = list(locus_data_map.keys())
    
    for gene, disease in locus_keys:
        locus_key = (gene, disease)
        locus_data = locus_data_map[locus_key]

        fig = create_violin_swarm(locus_data['df'], locus_data)
        if not fig:
            continue

        safe_gene = re.sub(r'[\\/]', '_', gene)
        safe_disease = re.sub(r'[\\/]', '_', disease)

        png_path = (output_dir_png / f"{safe_gene}_{safe_disease}.png")
        html_path = (output_dir_html / f"{safe_gene}_{safe_disease}.html")

        if TEST_MODE:
            print(f"Previewing: {gene} / {disease}")
            fig.show(renderer="browser")
            if SAVE_TEST_OUTPUTS:
                TEST_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
                fig.write_html(str(TEST_OUTPUT_DIR / html_path.name))
                fig.write_image(str(TEST_OUTPUT_DIR / png_path.name), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
            printed += 1
            if printed >= TEST_LIMIT:
                break
        else:
            fig.write_html(str(html_path))
            fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)

    # --- Finished ---
    print("--- Test mode ON: Test completed ---" if TEST_MODE else "--- Done ---")


# Guard to prevent execution on import
if __name__ == "__main__":
    main()