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
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import plotly.express as px
import plotly.io as pio
import argparse


# force Plotly to open figures in the system browser
pio.renderers.default = "browser"

# --- TEST MODE ---
TEST_MODE = False          # Quick testing: preview, limit work
TEST_LIMIT = 3            # Number of VCF records to process in test mode
SAVE_TEST_OUTPUTS = False  # If True, also save files when TEST_MODE is on

# --- File locations ---
from trplots.config import VCF_PATH, JSON_PATH, OUTPUT_BASE  

OUTPUT_BASE = Path(OUTPUT_BASE/"Plots/Path_Ref_Motif_Tandem_Repeats_Copy_Plots_v9")

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
with open(JSON_PATH, "r") as file:
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

def parse_args():
    p = argparse.ArgumentParser(description="Tandem repeat plot generator")
    p.add_argument("--test", dest="test", action="store_true", help="Enable test mode")
    p.add_argument("--no-test", dest="test", action="store_false", help="Disable test mode")
    p.set_defaults(test=TEST_MODE)
    p.add_argument("--test-limit", dest="test_limit", type=int, default=TEST_LIMIT)
    p.add_argument("--save-test-outputs", dest="save_test_outputs", action="store_true", default=SAVE_TEST_OUTPUTS)
    return p.parse_args()

def main(args=None):
    if args is None:
        args = parse_args()
    global TEST_MODE, TEST_LIMIT, SAVE_TEST_OUTPUTS
    TEST_MODE = bool(args.test)
    TEST_LIMIT = int(args.test_limit)
    SAVE_TEST_OUTPUTS = bool(args.save_test_outputs)

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

        # Collect all pathogenic/reference motif candidates (normalized). If none -> skip.
        raw_motif_field = (
            matched_locus.get("pathogenic_motif_reference_orientation")
            or matched_locus.get("pathogenic_motif_gene_orientation")
            or matched_locus.get("reference_motif_reference_orientation")
        )
        candidates = []
        if isinstance(raw_motif_field, list):
            for m in raw_motif_field:
                if isinstance(m, str) and m.strip():
                    candidates.append(re.sub(r"\s+", "", m).upper())
        elif isinstance(raw_motif_field, str) and raw_motif_field.strip():
            candidates.append(re.sub(r"\s+", "", raw_motif_field).upper())

        # dedupe while preserving order
        seen = set()
        candidates = [m for m in candidates if not (m in seen or seen.add(m))]
        if not candidates:
            print(f"Skip {chrom}:{pos} — empty pathogenic/reference motif")
            continue

        # For each motif candidate, compute repeat counts across samples and build a combined DataFrame.
        rows = []
        for motif in candidates:
            mlen = len(motif)
            if mlen <= 0:
                continue
            for sample_name, sample in record.samples.items():
                al_lengths = sample.get("AL")
                if al_lengths is None:
                    continue
                # normalize AL to iterable (some VCF callers store single int)
                if isinstance(al_lengths, (list, tuple, np.ndarray)):
                    al_iter = al_lengths
                else:
                    al_iter = [al_lengths]
                for length in al_iter:
                    if length is None:
                        continue
                    try:
                        allele_len = int(length)
                    except (ValueError, TypeError):
                        continue
                    # compute repeat count per motif for this allele (float, two decimals)
                    try:
                        rc = round(float(allele_len) / float(mlen), 2)
                    except Exception:
                        continue
                    rows.append({
                        "Repeat Count": rc,
                        "Motif": motif,
                        "Allele Length": allele_len,
                        "Sample": sample_name
                    })

        if not rows:
            # nothing to plot for any motif
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

        # Prepare dataframe (rows across motifs)
        df = pd.DataFrame(rows)
        df['Motif'] = df['Motif'].astype(str)


        # Subtitle (wrapped, preserve case inside parentheses)
        disease_text = disease if isinstance(disease, str) else str(disease)
        disease_text = re.sub(r"\s+;", ";", disease_text)
        disease_for_title = title_outside_parens(disease_text)

        subtitle = f"{chrom} - {pos} | {gene} - {disease_for_title}"

        # 107 is arbitrary to keep subtitle within plot width
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
        min_bin, max_bin = int(df['Repeat Count'].min()), int(df['Repeat Count'].max())
        # Use motif as color and create a stacked bar chart when multiple motifs are present.
        # (nbins chosen so bins align with integer repeat counts)
        fig = px.histogram(
            df, x="Repeat Count", color="Motif",
            nbins=(max_bin - min_bin + 1),
            title=(
                "<span style='font-size:18px; font-weight:bold; color:black'>"
                "Allele Size Distribution</span><br>" + subtitle_html
            ),
            color_discrete_sequence=px.colors.qualitative.Plotly
        )
        # Make bars stacked (motifs stack per bin) and remove transparency so stack segments are solid
        fig.update_layout(barmode='stack', bargap=0.1)
        fig.update_traces(opacity=1, marker_line_width=0)

        # --- Axis: keep median/p95 logic for the right edge ---
        arr = np.asarray(df['Repeat Count'].values, dtype=float)
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
        # safe short gene string for filenames
        safe_gene = re.sub(r'[^A-Za-z0-9]', '_', str(gene))[:10]
        safe_disease = re.sub(r'[^A-Za-z0-9]', '_', str(disease))[:10]
        safe_chrom = re.sub(r'[^A-Za-z0-9]', '_', str(chrom))[:10]
        safe_pos = re.sub(r'[^A-Za-z0-9]', '_', str(pos))[:10]

        if TEST_MODE:
            # Preview: small, no file save
            fig.update_layout(width=900, height=500)
            fig.show()
        else:
            # Full save: PNG and HTML
            out_prefix = os.path.join(OUTPUT_DIR, f"{safe_gene}_{safe_disease}_{safe_chrom}_{safe_pos}")
            html_path = out_prefix + ".html"
            png_path = out_prefix + ".png"

            # Write HTML first (full title preserved)
            try:
                fig.write_html(html_path, include_plotlyjs="cdn")
            except Exception as e:
                print(f"HTML save failed for {html_path}: {e}")

            # Shorten title before PNG export to avoid kaleido temp-file name limits,
            # then restore original title afterwards.
            orig_title = None
            short_title = f"{safe_gene}_{safe_chrom}_{safe_pos}"
            try:
                if getattr(fig.layout, "title", None) and getattr(fig.layout.title, "text", None):
                    orig_title = fig.layout.title.text
                    fig.update_layout(title_text=short_title)
                fig.write_image(png_path, width=900, height=500)
            except Exception as e:
                print(f"PNG export failed for {png_path}: {e} (HTML saved: {html_path})")
            finally:
                if orig_title is not None:
                    try:
                        fig.update_layout(title_text=orig_title)
                    except Exception:
                        pass

    print(f"Total VCF loci: {total_vcf_loci}, Matched loci: {matched_vcf_loci}")

    # --- Cleanup ---
    vcf_in.close()

if __name__ == "__main__":
    main()
