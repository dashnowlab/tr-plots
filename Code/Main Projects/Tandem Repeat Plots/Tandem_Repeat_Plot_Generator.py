"""
---------------------------------------------
 Script: Tandem Repeat Plot Generator
 Purpose:
   - Reads a VCF file that contains repeat lengths 
     for many genetic loci
   - Matches each variant in the VCF to metadata stored in a JSON file
     (info about genes, diseases, and thresholds)
   - For each locus, it collects repeat counts across samples
   - Plots a bar chart of allele repeat counts
   - Adds reference lines/ranges showing "benign", "intermediate", 
     or "pathogenic" thresholds from the metadata
   - Saves the plots as PNG (image) and HTML (interactive) files.
---------------------------------------------
"""

import json
import os
import pysam
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import re

# --- TEST MODE ---
TEST_MODE = True   # Toggle this flag for quick testing
TEST_LIMIT = 5     # how many VCF records to process when testing

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots"

VCF_PATH = f"{BASE_DIR}/Data/Sequencing Data/83 Loci 503 Samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"
METADATA_PATH = f"{BASE_DIR}/Data/Other Data/STRchive-loci.json"
OUTPUT_DIR = f"{BASE_DIR}/Results/Plots/Tandem_Repeats_Plots"

# Make sure the output folder exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# If test mode, use a subfolder for test outputs
if TEST_MODE:
    OUTPUT_DIR = os.path.join(OUTPUT_DIR, "test_outputs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Load the metadata file ---
with open(METADATA_PATH, "r") as file:
    loci_data = json.load(file)  # list of loci with info like gene name, disease, thresholds

# --- Helper function to draw reference lines or ranges on plots ---
def add_range_marker_or_line(fig, x0, x1, label_id, label_positions,
                              chart_width_px, x_span,
                              base_y_line=1.02, base_y_label=1.1,
                              bump_y=0.08, bump_px_threshold=30,
                              display_label=None, level=None):
    """
    Adds a marker or line to the plot:
    - If x0 == x1, it draws a vertical dashed line (for a single cutoff point).
    - If x0 != x1, it draws a horizontal line above the bars (for a range).
    It also places a text label so the user knows what the line/range means.
    """

    # Decide where to place the label (centered if a single line, left if a range)
    if x0 == x1:
        x_label = x0
        xanchor = "center"
    else:
        x_label = x0
        xanchor = "left"

    # Calculate how wide labels are in data space (so they don’t overlap too much)
    pixels_per_unit = chart_width_px / x_span
    label_pixel_width = len(label_id) * 7
    label_data_width = label_pixel_width / pixels_per_unit

    # Track if this label overlaps with others, and bump it higher if so
    bump_level = level if level is not None else 0
    if level is None:
        for other_x0, other_x1 in label_positions.values():
            overlap = (min(x1, other_x1) - max(x0, other_x0)) * pixels_per_unit
            if overlap >= -bump_px_threshold:
                bump_level += 1

    # Adjust vertical placement depending on bump level
    y_line = base_y_line + bump_y * bump_level
    y_label = base_y_label + bump_y * bump_level

    # Make a readable label
    label = display_label or label_id.capitalize()
    unique_key = f"{label_id}_{x0}_{x1}"
    label_positions[unique_key] = (x0, x1)

    # Add either a vertical line (if single point) or horizontal line (if range)
    if x0 == x1:
        fig.add_shape(
            type="line", x0=x0, x1=x0, y0=0, y1=1,
            xref="x", yref="paper",
            line=dict(color="black", width=1, dash="dot")
        )
        fig.add_annotation(
            x=x_label, y=y_label,
            text=label, showarrow=False,
            xref="x", yref="paper",
            font=dict(color="black", size=12),
            xanchor=xanchor
        )
    else:
        fig.add_shape(
            type="line", x0=x0, x1=x1, y0=y_line, y1=y_line,
            xref="x", yref="paper",
            line=dict(color="black", width=1)
        )
        fig.add_annotation(
            x=x_label, y=y_label,
            text=label, showarrow=False,
            xref="x", yref="paper",
            font=dict(color="black", size=12),
            xanchor=xanchor
        )

# --- Read the VCF file ---
vcf_in = pysam.VariantFile(VCF_PATH)

total_vcf_loci = 0       # keep track of how many sites are in the VCF
matched_vcf_loci = 0     # how many of those matched the JSON metadata

# --- Go through each record (variant) in the VCF ---
for i, record in enumerate(vcf_in.fetch()):              # enumerate so we can count
    if TEST_MODE and i >= TEST_LIMIT:                    # stop early in test mode
        break
    total_vcf_loci += 1
    chrom = record.chrom
    pos = record.pos   

    # Find the metadata entry that matches this chromosome + position
    matched_locus = next(
        (entry for entry in loci_data
         if chrom == entry["chrom"] and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
        None
    )
    if not matched_locus:
        print(f"{chrom}:{pos} not matched in loci_data")
        continue
    matched_vcf_loci += 1

    # Get the repeat motif (the repeating DNA sequence)
    motif = matched_locus.get("reference_motif_reference_orientation")
    if isinstance(motif, list):
        motif_str = "".join(motif)   # join list into a string if needed (metadata sometimes stores as list)
    else:
        motif_str = str(motif)
    motif_length = len(motif_str)

    # --- Collect repeat counts for each sample ---
    allele_lengths = []
    repeat_counts = []
    for sample in record.samples.values():
        al_lengths = sample.get("AL")  # AL = allele lengths in bases
        if al_lengths:
            allele_lengths.extend([length for length in al_lengths if length is not None])
            repeat_counts.extend([length // motif_length for length in al_lengths if length is not None])

    if not repeat_counts:
        continue  # skip if no counts found

    # --- Get extra metadata for this locus ---
    gene = matched_locus.get("gene", "Unknown")
    disease = matched_locus.get("disease", "Unknown")
    benign_min = matched_locus.get("benign_min")
    benign_max = matched_locus.get("benign_max")
    pathogenic_min = matched_locus.get("pathogenic_min")
    pathogenic_max = matched_locus.get("pathogenic_max")
    intermediate_min = matched_locus.get("intermediate_min")
    intermediate_max = matched_locus.get("intermediate_max")

    # --- Put repeat counts into a dataframe for plotting ---
    df = pd.DataFrame({"Repeat Count": repeat_counts})
    min_bin, max_bin = min(repeat_counts), max(repeat_counts)

    # Format a subtitle for the plot
    disease_clean = re.sub(r"\s+;", ";", disease) if isinstance(disease, str) else disease
    subtitle = f"{chrom.title()} - {pos} | {gene} - {disease_clean if isinstance(disease_clean, str) else disease_clean}"

    # Wrap long subtitles onto multiple lines
    max_line_len = 107
    subtitle_lines, line = [], ""
    for segment in subtitle.split(" "):
        if len(line) + len(segment) + 1 > max_line_len:
            subtitle_lines.append(line.rstrip())
            line = ""
        line += segment + " "
    if line:
        subtitle_lines.append(line.rstrip())
    subtitle_html = "<br>".join(f"<span style='font-size:12px; color:black'>{l}</span>" for l in subtitle_lines)

    # --- Make histogram plot of repeat counts ---
    fig = px.histogram(
        df, x="Repeat Count",
        nbins=(max_bin - min_bin + 1),
        title=f"<span style='font-size:18px; font-weight:bold; color:black'>Allele Size Distribution</span><br>{subtitle_html}",
    )

    # --- Style the plot (size, labels, axes, background, etc.) ---
    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=125),
        title_font_size=18,
        xaxis_title="Repeat Count",
        yaxis_title="Allele Count",
        bargap=0.7,
        plot_bgcolor='white',
        title=dict(y=0.95, font=dict(color='black'))
    )
    fig.update_xaxes(showline=True, linecolor='black', zeroline=False)
    fig.update_yaxes(showline=True, linecolor='black', zeroline=False)
    fig.update_traces(hovertemplate="Repeat Count=%{x}<br>Allele Count=%{y}<extra></extra>")
    fig.update_layout(xaxis_range=[None, max(repeat_counts) + 2.5])

    # --- Add reference ranges (benign, intermediate, pathogenic) ---
    x_range = fig.layout.xaxis.range
    x_span = (x_range[1] - x_range[0]) if x_range and None not in x_range else max(repeat_counts) - min(repeat_counts) + 2.5
    label_positions = {}

    if benign_min is not None and benign_max is not None:
        add_range_marker_or_line(fig, benign_min, benign_max, "Benign", label_positions,
                                 chart_width_px=fig.layout.width, x_span=x_span, level=0)
    if intermediate_min is not None and intermediate_max is not None:
        add_range_marker_or_line(fig, intermediate_min, intermediate_max, "Intermediate", label_positions,
                                 chart_width_px=fig.layout.width, x_span=x_span, level=1)

    # Special case for VWA1 gene (two pathogenic ranges)
    if gene == "VWA1":
        add_range_marker_or_line(fig, 1, 1.85, "Pathogenic", label_positions,
                                 chart_width_px=fig.layout.width, x_span=x_span, level=2)
        add_range_marker_or_line(fig, 2.15, 3, "Pathogenic", label_positions,
                                 chart_width_px=fig.layout.width, x_span=x_span, level=2)
    elif pathogenic_min is not None and pathogenic_max is not None:
        add_range_marker_or_line(fig, pathogenic_min, pathogenic_max, "Pathogenic", label_positions,
                                 chart_width_px=fig.layout.width, x_span=x_span, level=2)

    # --- Save plot as PNG and HTML ---
    fig.write_image(os.path.join(OUTPUT_DIR, f"{gene}_{chrom}_{pos}_allele_dist.png"),
                    format="png", width=900, height=500, scale=2)
    fig.write_html(os.path.join(OUTPUT_DIR, f"{gene}_{chrom}_{pos}_allele_dist.html"))

# --- Finished ---
if TEST_MODE:
    print(f"--- Test mode ON: processed only the first {TEST_LIMIT} records ---")
else:
    print("--- Done ---")