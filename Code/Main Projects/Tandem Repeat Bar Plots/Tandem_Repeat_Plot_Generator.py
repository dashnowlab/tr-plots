"""
---------------------------------------------
 Script: Tandem Repeat Plot Generator
 Purpose:
   - Reads a VCF file that contains repeat lengths for many loci
   - Matches each variant in the VCF to JSON metadata (gene, disease, thresholds)
   - For each locus, collects repeat counts across samples
   - Plots a histogram of allele repeat counts
   - Adds reference lines/ranges for benign/intermediate/pathogenic thresholds
   - Saves plots as PNG and HTML (non-test mode saves silently)
   - Test mode: preview limited plots; optionally save to test_outputs/
---------------------------------------------
"""

import json
import os
import re
import numpy as np
import pandas as pd
import pysam
import plotly.express as px
import plotly.graph_objects as go

# --- TEST MODE ---
TEST_MODE = True          # Quick testing: preview, limit work
TEST_LIMIT = 3            # Number of VCF records to process in test mode
SAVE_TEST_OUTPUTS = False  # If True, also save files when TEST_MODE is on

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots"

VCF_PATH = f"{BASE_DIR}/Data/Sequencing Data/83 Loci 503 Samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"
METADATA_PATH = f"{BASE_DIR}/Data/Other Data/STRchive-loci.json"
OUTPUT_BASE = os.path.join(BASE_DIR, "Results/Plots/Tandem_Repeats_Bar_Plots")

if TEST_MODE:
    OUTPUT_DIR = os.path.join(OUTPUT_BASE, "test_outputs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    OUTPUT_HTML_DIR = OUTPUT_DIR
    OUTPUT_PNG_DIR  = OUTPUT_DIR
else:
    OUTPUT_DIR = OUTPUT_BASE
    OUTPUT_HTML_DIR = os.path.join(OUTPUT_DIR, "HTML")
    OUTPUT_PNG_DIR  = os.path.join(OUTPUT_DIR, "PNG")
    os.makedirs(OUTPUT_HTML_DIR, exist_ok=True)
    os.makedirs(OUTPUT_PNG_DIR,  exist_ok=True)

# --- Load the metadata file ---
with open(METADATA_PATH, "r") as file:
    loci_data = json.load(file)  # list of loci with info like gene name, disease, thresholds

# --- Helpers ---
def norm_chrom(s: str) -> str:
    """Return chromosome name with 'chr' prefix (case-insensitive)."""
    s = str(s)
    return s if s.lower().startswith("chr") else f"chr{s}"

def title_outside_parens(s: str) -> str:
    """Title-case text outside (...) while preserving case inside parentheses."""
    if not isinstance(s, str):
        return s
    parts = re.split(r'(\([^)]*\))', s)  # keep paren groups
    return ''.join(p if (p.startswith('(') and p.endswith(')')) else p.title()
                   for p in parts)

def add_range_marker_or_line(fig, x0, x1, label_id, label_positions,
                              chart_width_px, x_span,
                              base_y_line=1.02, base_y_label=1.1,
                              bump_y=0.08, bump_px_threshold=30,
                              display_label=None, level=None):
    """
    If x0==x1 -> vertical dashed line (single cutoff)
    Else       -> horizontal line above bars (range)
    Also adds a text label, bumping up if overlapping.
    """
    if x0 == x1:
        x_label = x0
        xanchor = "center"
    else:
        x_label = x0
        xanchor = "left"

    pixels_per_unit = chart_width_px / max(x_span, 1e-9)

    bump_level = level if level is not None else 0
    if level is None:
        for other_x0, other_x1 in label_positions.values():
            overlap = (min(x1, other_x1) - max(x0, other_x0)) * pixels_per_unit
            if overlap >= -bump_px_threshold:
                bump_level += 1

    y_line = base_y_line + bump_y * bump_level
    y_label = base_y_label + bump_y * bump_level

    label = display_label or label_id.capitalize()
    unique_key = f"{label_id}_{x0}_{x1}"
    label_positions[unique_key] = (x0, x1)

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

total_vcf_loci = 0
matched_vcf_loci = 0

# Choose an iterator that works with or without an index
try:
    iterator = vcf_in.fetch()
except ValueError:
    iterator = vcf_in  # fall back: sequential

# --- Iterate through VCF records ---
for i, record in enumerate(iterator):
    if TEST_MODE and i >= TEST_LIMIT:
        break

    total_vcf_loci += 1
    pos = record.pos
    chrom = norm_chrom(record.chrom)

    # Find matching locus in metadata (normalize JSON chrom too)
    matched_locus = next(
        (entry for entry in loci_data
         if norm_chrom(entry["chrom"]) == chrom and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
        None
    )
    if not matched_locus:
        print(f"{chrom}:{pos} not matched in loci_data")
        continue
    matched_vcf_loci += 1

    # Motif (repeating DNA sequence)
    motif = matched_locus.get("reference_motif_reference_orientation")
    motif_str = "".join(motif) if isinstance(motif, list) else str(motif or "")
    if not motif_str:
        print(f"Skip {chrom}:{pos} — empty motif")
        continue
    motif_length = len(motif_str)
    if motif_length == 0:
        print(f"Skip {chrom}:{pos} — motif length 0")
        continue

    # Collect repeat counts across samples
    repeat_counts = []
    for sample in record.samples.values():
        al_lengths = sample.get("AL")  # allele lengths in bases
        if al_lengths:
            repeat_counts.extend([length // motif_length for length in al_lengths if length is not None])

    if not repeat_counts:
        continue

    # Locus metadata
    gene = matched_locus.get("gene", "Unknown")
    disease = matched_locus.get("disease", "Unknown")
    benign_min = matched_locus.get("benign_min")
    benign_max = matched_locus.get("benign_max")
    pathogenic_min = matched_locus.get("pathogenic_min")
    pathogenic_max = matched_locus.get("pathogenic_max")
    intermediate_min = matched_locus.get("intermediate_min")
    intermediate_max = matched_locus.get("intermediate_max")

    # Prepare dataframe
    df = pd.DataFrame({"Repeat Count": repeat_counts})
    min_bin, max_bin = int(min(repeat_counts)), int(max(repeat_counts))

    # Subtitle (wrapped, preserve case inside parentheses)
    disease_text = disease if isinstance(disease, str) else str(disease)
    disease_text = re.sub(r"\s+;", ";", disease_text)
    disease_for_title = title_outside_parens(disease_text)

    subtitle = f"{chrom} - {pos} | {gene} - {disease_for_title}"

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

    # --- Histogram ---
    min_bin, max_bin = int(min(repeat_counts)), int(max(repeat_counts))
    fig = px.histogram(
        df, x="Repeat Count",
        nbins=(max_bin - min_bin + 1),   # 1 bin per integer (Plotly sizes bars)
        title=(
            "<span style='font-size:18px; font-weight:bold; color:black'>"
            "Allele Size Distribution</span><br>" + subtitle_html
        ),
    )
    fig.update_traces(marker_color="lightblue")
    fig.update_layout(bargap=0.7)

    # --- Axis: keep median/p95 logic for the right edge ---
    arr = np.asarray(repeat_counts, dtype=float)
    median = float(np.median(arr))
    p95 = float(np.quantile(arr, 0.95))
    x_right = max(2 * median + 0.5, p95 + 0.5, 5.5)
    fig.update_xaxes(range=[0, x_right])

    # Styling
    fig.update_layout(
        width=900, height=500, margin=dict(t=125), title_font_size=18,
        xaxis_title="Repeat Count", yaxis_title="Allele Count",
        plot_bgcolor='white', title=dict(y=0.95, font=dict(color='black'))
    )
    fig.update_xaxes(showline=True, linecolor='black', zeroline=False, showgrid=False, ticks="outside")
    fig.update_yaxes(showline=True, linecolor='black', zeroline=False, showgrid=True, gridcolor='rgba(0,0,0,0.08)', ticks="outside")
    fig.update_traces(hovertemplate="Repeat Count=%{x}<br>Allele Count=%{y}<extra></extra>")

    # Reference ranges
    x_range = fig.layout.xaxis.range
    x_span = (x_range[1] - x_range[0]) if x_range and None not in x_range else ((max_bin + 2.5) - 0)
    label_positions = {}

    # Draw benign (always at level 0)
    if benign_min is not None and benign_max is not None:
        add_range_marker_or_line(
            fig, benign_min, benign_max, "Benign", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=0
        )

    # Draw intermediate (if present, use level 1)
    has_intermediate = (intermediate_min is not None and intermediate_max is not None)
    if has_intermediate:
        add_range_marker_or_line(
            fig, intermediate_min, intermediate_max, "Intermediate", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=1
        )

    # If there's no intermediate range for this locus, pathogenic moves to level 1; otherwise level 2
    pathogenic_level = 1 if not has_intermediate else 2

    # Draw pathogenic
    if gene == "VWA1":
        # VWA1 has two pathogenic windows in your metadata/spec
        add_range_marker_or_line(
            fig, 1, 1, "Pathogenic", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=pathogenic_level
        )
        add_range_marker_or_line(
            fig, 3, 3, "Pathogenic", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=pathogenic_level
        )
    elif pathogenic_min is not None and pathogenic_max is not None:
        add_range_marker_or_line(
            fig, pathogenic_min, pathogenic_max, "Pathogenic", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=pathogenic_level
        )

    # --- Save / Preview ---
    safe_gene = re.sub(r'[\\/]', '_', gene)
    html_path = os.path.join(OUTPUT_HTML_DIR, f"{safe_gene}_{chrom}_{pos}_allele_dist.html")
    png_path  = os.path.join(OUTPUT_PNG_DIR,  f"{safe_gene}_{chrom}_{pos}_allele_dist.png")

    if TEST_MODE:
        fig.show()
        if SAVE_TEST_OUTPUTS:
            fig.write_html(html_path)
            try:
                fig.write_image(png_path, format="png", width=900, height=500, scale=2)
            except Exception as e:
                print(f"PNG export failed for {png_path}: {e} (HTML saved)")
    else:
        fig.write_html(html_path)
        try:
            fig.write_image(png_path, format="png", width=900, height=500, scale=2)
        except Exception as e:
            print(f"PNG export failed for {png_path}: {e} (HTML saved)")

# --- Finished ---
if TEST_MODE:
    print(f"--- Test mode ON: processed only the first {TEST_LIMIT} records; matched {matched_vcf_loci} ---")
else:
    print(f"--- Done: processed {total_vcf_loci}, matched {matched_vcf_loci} ---")
