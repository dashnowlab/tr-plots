import re
import ast
import pandas as pd
import plotly.express as px
import plotly.io as pio
from pathlib import Path
import argparse

# Pull paths from the shared config
# NOTE: Removed MASTER_ALLELE_SPREADSHEET.
from trplots.config import (
    OTHER_DATA, # <--- Used to construct the input path
    ENSURE_DIR,
    ALLELE_LENGTH_PLOTS_OUTPUT,
)

# make Plotly open figures in the browser by default
pio.renderers.default = "browser"

# --- TEST MODE ---
TEST_MODE = True              # Toggle this flag for quick testing (preview only)
TEST_LIMIT = 3                 # How many (gene, disease) plots to generate in test mode
SAVE_TEST_OUTPUTS = False      # Toggle saving plots when in test mode

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130               # fixed header space (annotations), keeps plot area aligned
PNG_SCALE = 2

# --- File locations (via config) ---
# UPDATED: Revert to previous file name, assuming it's in the OTHER_DATA directory.
if OTHER_DATA is None:
    raise FileNotFoundError("The 'OTHER_DATA' path is missing in trplots.config. Cannot locate input CSV.")
DATA_CSV = OTHER_DATA / "83_loci_503_samples_withancestrycolumns.csv"


# Default output roots under results/plots/allele_length_boxplots/violin_swarm/...
OUTPUT_ROOT = ALLELE_LENGTH_PLOTS_OUTPUT / "violin_swarm"
# ENSURE_DIR is correctly imported and used below.
OUTPUT_DIR_PNG  = ALLELE_LENGTH_PLOTS_OUTPUT / "violin_swarm" / "png"
OUTPUT_DIR_HTML = ALLELE_LENGTH_PLOTS_OUTPUT / "violin_swarm" / "html"


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
    """
    Build a horizontal violin plot with overlaid points ("swarm")
    using pre-calculated data.
    """
    gene, disease = locus_key_data['key']
    counts = locus_key_data['counts']
    inheritance = locus_key_data['inheritance']
    
    if locus_df.empty:
        return None

    # Get ordered categories from pre-calculated counts
    ordered_categories = counts.index.tolist()
    
    # Prepare dataframe for plotting
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
    p.add_argument("--data-csv", type=str, default=str(DATA_CSV),
                   help="Override the input CSV file")
    p.add_argument("--output-dir", type=str, default=str(OUTPUT_ROOT),
                   help="Override the output directory")
    p.set_defaults(test=TEST_MODE) 
    return p.parse_args()

# --- Main function: move top-level processing here ---
def main():
    args = parse_args()

    # --- Test mode overrides ---
    global TEST_MODE, TEST_LIMIT, OUTPUT_ROOT, OUTPUT_DIR_PNG, OUTPUT_DIR_HTML
    TEST_MODE = args.test
    TEST_LIMIT = args.limit
    
    # Update paths using ENSURE_DIR from the config
    OUTPUT_ROOT = Path(args.output_dir)
    
    # We must ensure the actual output directories exist using the config's helper
    if not TEST_MODE:
        # Use ENSURE_DIR to guarantee these paths exist if not in test mode
        base_parts = list(ALLELE_LENGTH_PLOTS_OUTPUT.parts[len(ALLELE_LENGTH_PLOTS_OUTPUT.parents[1].parts):])
        OUTPUT_DIR_PNG = ENSURE_DIR(*(base_parts + ["violin_swarm", "png"]))
        OUTPUT_DIR_HTML = ENSURE_DIR(*(base_parts + ["violin_swarm", "html"]))
    else:
        # If test mode is set, use the "test_outputs" folder and ensure it exists
        OUTPUT_ROOT = ENSURE_DIR("plots", "allele_length_boxplots", "violin_swarm", "test_outputs")
        OUTPUT_DIR_PNG = OUTPUT_ROOT
        OUTPUT_DIR_HTML = OUTPUT_ROOT


    # Data CSV override (pathlib.Path compatible)
    global DATA_CSV
    DATA_CSV = Path(args.data_csv)

    # Load data
    print(f"Loading data from: {DATA_CSV}")
    df = pd.read_csv(DATA_CSV)

    # --- Add 'All' population (Vectorized) ---
    df_all = df.assign(SuperPop="All")
    df = pd.concat([df, df_all], ignore_index=True)
    
    # Assign 'SuperPop' as a categorical type for consistent ordering/grouping
    df['SuperPop'] = pd.Categorical(df['SuperPop'], categories=SUPERPOP_ORDER, ordered=True)

    # --- Flag pathogenic individuals using thresholds (Vectorized) ---
    df['Is pathogenic'] = (
        (df['Allele length'] >= df['Pathogenic min']) &
        (df['Allele length'] <= df['Pathogenic max'])
    )

    # --- Pre-calculate all data outside the plot loop (Performance optimization) ---

    # 1. Total unique individuals per Locus/SuperPop (Required for header)
    locus_pop_counts = (
        df.groupby(['Gene', 'Disease', 'SuperPop'], observed=True)['Sample ID']
          .nunique()
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

        png_path = (OUTPUT_DIR_PNG / f"{safe_gene}_{safe_disease}_allele_length_violin_swarm.png")
        html_path = (OUTPUT_DIR_HTML / f"{safe_gene}_{safe_disease}_allele_length_violin_swarm.html")

        if TEST_MODE:
            print(f"Previewing: {gene} / {disease}")
            fig.show(renderer="browser")
            if SAVE_TEST_OUTPUTS:
                fig.write_html(str(html_path))
                fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
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