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
import argparse
from pathlib import Path

# import shared paths/helpers from config
from trplots.config import (
    OTHER_DATA,                
    ENSURE_DIR
)

pio.renderers.default = "browser"

# --- Input file locations (via config) ---
DATA_PATH = (OTHER_DATA / "master_allele_spreadsheet.xlsx")
SHEET_NAME = "Integrated Alleles"

# --- Output roots (Simplified/centralized definition) ---
# Ensuring base directories exist once outside of the plotting loop
OUTPUT_DIR_PNG_BASE  = ENSURE_DIR("plots", "allele_length_boxplots", "png")
OUTPUT_DIR_HTML_BASE = ENSURE_DIR("plots", "allele_length_boxplots", "html")

# --- TEST MODE ---
TEST_MODE_DEFAULT = True         # Default state
TEST_LIMIT_DEFAULT = 8            # Default limit
SAVE_TEST_OUTPUTS_DEFAULT = False # Default save state

def parse_args():
    """Parse command line arguments to override defaults."""
    p = argparse.ArgumentParser(description="Allele length boxplot generator")
    p.add_argument("--test", dest="test", action="store_true", help="Enable test mode")
    p.add_argument("--no-test", dest="test", action="store_false", help="Disable test mode")
    p.set_defaults(test=TEST_MODE_DEFAULT)
    p.add_argument("--test-limit", type=int, default=TEST_LIMIT_DEFAULT, help="Number of plots to generate in test mode")
    p.add_argument("--save-test-outputs", action="store_true", default=SAVE_TEST_OUTPUTS_DEFAULT,
                   help="Save outputs even when running in test mode")
    # Use the new DATA_PATH as the default for the CLI argument
    p.add_argument("--data-path", type=str, default=str(DATA_PATH),
                   help="Override the input data path (Excel)")
    p.add_argument("--output-dir", type=str, default=None,
                   help="Optional: override the output directory (single folder used for PNG/HTML)")
    return p.parse_args()

# --- Figure sizing (standardized) ---
FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130               # space for the header lines
PNG_SCALE = 2

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

# --- Helpers (Unchanged) ---
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
            # Use float conversion first to safely handle scientific notation/decimals
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

# -------------------- Plotting (ALL alleles) - Optimized Function --------------------
def create_boxplot_all(locus_data, locus_metadata, show_no_data_note=True, show_points=False):
    """
    Build a horizontal boxplot for one Locus using **pre-calculated data**.
    """
    
    # Extract metadata from the pre-calculated dictionary
    gene = locus_metadata['gene']
    disease = locus_metadata['disease']
    locus = locus_metadata['locus']
    inheritance = locus_metadata['inheritance']
    counts = locus_metadata['counts']
    ordered_categories = locus_metadata['ordered_categories']
    
    # locus_data is the pre-filtered DataFrame subset (Allele length, SuperPop)
    
    # --------- Helper function for boxplot statistics (Tukey method) ----------
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
        
        # Calculate whiskers (minimum/maximum values within 1.5 * IQR)
        lower_whisk = x[x >= lf].min() if np.any(x >= lf) else x.min()
        upper_whisk = x[x <= uf].max() if np.any(x <= uf) else x.max()
        outliers = x[(x < lf) | (x > uf)]
        
        return q1, med, q3, lower_whisk, upper_whisk, outliers

    # --------- Build figure with explicit box stats + separate outlier points ----------
    fig = go.Figure()

    # Group the pre-filtered data for fast retrieval inside the loop
    locus_groups = locus_data.groupby('SuperPop', observed=True)['Allele length']

    for pop in ordered_categories:
        
        # Get data for the population
        if pop in locus_groups.groups:
            xs = locus_groups.get_group(pop).to_numpy(dtype=float)
            xs = xs[~np.isnan(xs)]
        else:
            xs = np.array([]) 

        if xs.size == 0:
            # Add an empty box to ensure the y-axis is present for this population
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

        # Add the main box trace
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
                boxpoints=False, # We handle points/outliers separately
                showlegend=False,
                hoveron="boxes",
                hovertemplate=None # Use default info
            )
        )

        # Add outliers as scatter points
        if show_points and out.size:
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

    # Axis ordering (must be done outside the loop)
    fig.update_yaxes(
        categoryorder='array',
        categoryarray=ordered_categories,
        tickmode='array',
        tickvals=ordered_categories,
        autorange='reversed',
        title=""
    )

    # Counts per population (unique individuals) -> pre-calculated
    pop_desc = ', '.join(f"{pop}: {counts.get(pop, 0)}" for pop in ordered_categories)
    pop_lines = _wrap_to_lines(pop_desc, max_len=110)
    
    # Annotations for "no data"
    missing = [c for c in ordered_categories if counts.get(c, 0) == 0]
    if show_no_data_note and missing:
        for cat in missing:
            fig.add_annotation(
                x=0.0, xref='paper', xanchor='left',
                y=cat, yref='y', yanchor='middle',
                text="<span style='font-size:11px;color:#888'>no data</span>",
                showarrow=False, align='left'
            )

    # Layout (unchanged)
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

    # Fixed-position header via annotations
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

# --- Main entrypoint (make script import-safe) ---
def main(args=None):
    
    # Parse arguments first
    args = parse_args()
    
    # Establish final test mode and path settings
    test_mode = args.test
    test_limit = args.test_limit
    save_test_outputs = args.save_test_outputs
    
    data_path = Path(args.data_path)
    
    # Initialize output paths based on their base definitions
    output_dir_png = OUTPUT_DIR_PNG_BASE
    output_dir_html = OUTPUT_DIR_HTML_BASE

    # Apply output-dir override if provided
    if args.output_dir:
        # If output_dir is provided, all outputs go to this single folder
        outdir = Path(args.output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        output_dir_png = outdir
        output_dir_html = outdir

    # Apply test mode override
    if test_mode:
        # If test mode is ON, and no custom output_dir was given, use the 'test_outputs' subfolder
        if not args.output_dir:
            test_output = ENSURE_DIR("plots", "allele_length_boxplots", "test_outputs")
            output_dir_png = test_output
            output_dir_html = test_output
        print(f"--- Test mode ON (limit={test_limit}) ---")
    
    # 2. Load data
    try:
        df = pd.read_excel(data_path, sheet_name=SHEET_NAME)
    except FileNotFoundError:
        print(f"Error: Input data file not found at: {data_path}")
        return
    except ValueError:
        print(f"Error: Sheet '{SHEET_NAME}' not found in file: {data_path}. Check SHEET_NAME variable.")
        return

    # --- Data Cleaning and Preparation (Vectorized and done once) ---
    
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
    if 'SuperPop' not in df.columns:
         raise ValueError("Missing 'SuperPop' column in input data.")
         
    df = pd.concat([df, df_all], ignore_index=True)
    
    # Ensure SuperPop is categorical for consistent sorting in groups
    df['SuperPop'] = pd.Categorical(df['SuperPop'], categories=SUPERPOP_ORDER, ordered=True)

    # --- Safety: ensure required columns exist ---
    for required in ['Gene', 'Disease', 'Locus', 'SuperPop', 'Sample ID Cleaned', 'Allele length', 'Inheritance_norm']:
        if required not in df.columns:
            raise ValueError(f"Missing required column: '{required}'. Please check the input file structure.")

    # --- 3. Pre-calculate Locus Data and Metadata (Performance gain) ---
    
    # Identify unique loci (Gene, Disease, Locus triplet)
    # Group the main DataFrame once by the Locus identifier for fast subsetting
    locus_groups = df.groupby(['Gene', 'Disease', 'Locus'], observed=True)
    locus_data_map = {}
    
    for (gene, disease, locus), group_df in locus_groups:
        
        # Calculate unique individual counts for all populations present in this locus
        counts_series = (
            group_df.groupby('SuperPop', observed=True)['Sample ID Cleaned']
                    .nunique()
        ).reindex(SUPERPOP_ORDER).fillna(0).astype(int)
        
        # Extract inheritance mode (guaranteed to exist due to cleanup step)
        inheritance = group_df['Inheritance_norm'].iloc[0]

        # Store the pre-filtered subset and metadata
        locus_data_map[(gene, disease, locus)] = {
            # Only store the necessary columns for plotting (Allele length, SuperPop)
            'data': group_df[['Allele length', 'SuperPop']].copy(),
            'metadata': {
                'gene': gene,
                'disease': disease,
                'locus': locus,
                'inheritance': inheritance,
                'counts': counts_series.to_dict(),
                'ordered_categories': SUPERPOP_ORDER
            }
        }
        
    # --- 4. Generate plots: one per Locus (using pre-calculated data) ---
    printed = 0

    for _, data_entry in locus_data_map.items():
        gene = data_entry['metadata']['gene']
        disease = data_entry['metadata']['disease']
        locus = data_entry['metadata']['locus']
        
        # Skip if no data exists 
        if data_entry['data'].empty or data_entry['data']['Allele length'].dropna().empty:
            continue

        # Call the optimized plotting function
        fig = create_boxplot_all(
            data_entry['data'], 
            data_entry['metadata'], 
            show_no_data_note=True, 
            show_points=False 
        )

        # Safe filenames
        safe_gene = re.sub(r'[\\/]', '_', str(gene))
        safe_disease = re.sub(r'[\\/]', '_', str(disease))
        safe_locus = re.sub(r'[\\/:*?"<>|]+', '_', str(locus))[:120]

        # Use the finalized output directories
        png_path  = output_dir_png  / f"{safe_gene}_{safe_disease}_{safe_locus}_allele_length_boxplot.png"
        html_path = output_dir_html / f"{safe_gene}_{safe_disease}_{safe_locus}_allele_length_boxplot.html"

        if test_mode:
            print(f"Previewing: {gene} / {disease} / {locus}")
            fig.show()
            if save_test_outputs:
                try:
                    fig.write_html(str(html_path))
                except Exception as e:
                    print(f"Failed to write HTML {html_path}: {e}")
                try:
                    fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
                except Exception as e:
                    print(f"Failed to write PNG {png_path}: {e} — ensure 'kaleido' is installed")
            printed += 1
            if printed >= test_limit:
                break
        else:
            try:
                fig.write_html(str(html_path))
            except Exception as e:
                print(f"Failed to write HTML {html_path}: {e}")
            try:
                fig.write_image(str(png_path), width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
            except Exception as e:
                print(f"Failed to write PNG {png_path}: {e} — ensure 'kaleido' is installed")

    if test_mode:
        print(f"--- Test mode ON: Test completed (limit={test_limit}) ---")
    else:
        print("--- Done ---")
 
if __name__ == "__main__":
    main()