import json
import os
import pysam
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import re

# File paths
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

JSON_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/Other/STRchive-loci.json"

# OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Tandem_Repeats_Plots/XMax_Repeat_Count/68_loci_100_samples_tr_plots"
# OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Tandem_Repeats_Plots/XMax_Repeat_Count/83_loci_88_samples_tr_plots"
OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/Tandem_Repeats_Plots/XMax_Repeat_Count/83_loci_503_samples_tr_plots2"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load JSON metadata
with open(JSON_PATH, "r") as file:
    loci_data = json.load(file)

def add_range_marker_or_line(fig, x0, x1, label_id, label_positions,
                              chart_width_px, x_span,
                              base_y_line=1.02, base_y_label=1.1,
                              bump_y=0.08, bump_px_threshold=30,
                              display_label=None, level=None):
    if x0 == x1:
        x_label = x0
        xanchor = "center"
    else:
        x_label = x0
        xanchor = "left"

    pixels_per_unit = chart_width_px / x_span
    label_pixel_width = len(label_id) * 7
    label_data_width = label_pixel_width / pixels_per_unit

    x0_padded = x0 - label_data_width / 2
    x1_padded = x1 + label_data_width / 2

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

vcf_in = pysam.VariantFile(VCF_PATH)

total_vcf_loci = 0
matched_vcf_loci = 0

for record in vcf_in.fetch():
    total_vcf_loci += 1
    chrom = record.chrom
    pos = record.pos

    matched_locus = next(
        (entry for entry in loci_data
         if chrom == entry["chrom"] and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
        None
    )
    if not matched_locus:
        print(f"{chrom}:{pos} not matched in loci_data")
        continue
    matched_vcf_loci += 1

    motif = matched_locus.get("reference_motif_reference_orientation")
    if isinstance(motif, list):
        motif_str = "".join(motif)
    else:
        motif_str = str(motif)
    motif_length = len(motif_str)
    motif = motif_str

    allele_lengths = []
    repeat_counts = []
    for sample in record.samples.values():
        al_lengths = sample.get("AL")
        if al_lengths:
            allele_lengths.extend([length for length in al_lengths if length is not None])
            repeat_counts.extend([length // motif_length for length in al_lengths if length is not None])

    if not repeat_counts:
        continue

    gene = matched_locus.get("gene", "Unknown")
    disease = matched_locus.get("disease", "Unknown")
    benign_min = matched_locus.get("benign_min")
    benign_max = matched_locus.get("benign_max")
    pathogenic_min = matched_locus.get("pathogenic_min")
    pathogenic_max = matched_locus.get("pathogenic_max")
    intermediate_min = matched_locus.get("intermediate_min")
    intermediate_max = matched_locus.get("intermediate_max")

    df = pd.DataFrame({"Repeat Count": repeat_counts})
    min_bin, max_bin = min(repeat_counts), max(repeat_counts)

    disease_clean = re.sub(r"\s+;", ";", disease) if isinstance(disease, str) else disease
    subtitle = f"{chrom.title()} - {pos} | {gene} - {disease_clean if isinstance(disease_clean, str) else disease_clean}"
    max_line_len = 107
    subtitle_lines = []
    line = ""

    for segment in subtitle.split(" "):
        if len(line) + len(segment) + 1 > max_line_len:
            subtitle_lines.append(line.rstrip())
            line = ""
        line += segment + " "
    if line:
        subtitle_lines.append(line.rstrip())
    subtitle_html = "<br>".join(f"<span style='font-size:12px'>{l}</span>" for l in subtitle_lines)

    fig = px.histogram(
        df, x="Repeat Count",
        nbins=(max_bin - min_bin + 1),
        title=f"<span style='font-size:18px; font-weight:bold'>Allele Size Distribution</span><br>{subtitle_html}",
    )

    fig.update_layout(
        width=900,
        height=500,
        margin=dict(t=125),
        title_font_size=18,
        title_font=dict(weight="bold"),
        title_font_color="black",
        xaxis_title="Repeat Count",
        yaxis_title="Allele Count",
        bargap=0.7,
        plot_bgcolor='white',
        title=dict(y=0.95)
    )

    fig.update_xaxes(
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        tickcolor='black',
        showline=True,
        linecolor='black',
        linewidth=1,
        zeroline=False
    )

    fig.update_yaxes(
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        tickcolor='black',
        showline=True,
        linecolor='black',
        linewidth=1,
        zeroline=False
    )

    fig.update_traces(
        hovertemplate="Repeat Count=%{x}<br>Allele Count=%{y}<extra></extra>"
    )

    fig.update_layout(xaxis_range=[None, max(repeat_counts) + 2.5])
    x_range = fig.layout.xaxis.range
    x_span = (x_range[1] - x_range[0]) if x_range and None not in x_range else max(repeat_counts) - min(repeat_counts) + 2.5

    label_positions = {}
    if benign_min is not None and benign_max is not None:
        add_range_marker_or_line(
            fig, benign_min, benign_max, "Benign", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=0
        )

    if intermediate_min is not None and intermediate_max is not None:
        add_range_marker_or_line(
            fig, intermediate_min, intermediate_max, "Intermediate", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=1
        )

    if gene == "VWA1":
        add_range_marker_or_line(
            fig, 1, 1.85, "Pathogenic", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=2
        )
        add_range_marker_or_line(
            fig, 2.15, 3, "Pathogenic", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=2
        )
    elif pathogenic_min is not None and pathogenic_max is not None:
        add_range_marker_or_line(
            fig, pathogenic_min, pathogenic_max, "Pathogenic", label_positions,
            chart_width_px=fig.layout.width, x_span=x_span, level=2
        )

    fig.write_image(os.path.join(OUTPUT_DIR, f"{gene}_{chrom}_{pos}_allele_dist.png"),
                    format="png", width=900, height=500, scale=2)
    fig.write_html(os.path.join(OUTPUT_DIR, f"{gene}_{chrom}_{pos}_allele_dist.html"))

print("--- Done ---")