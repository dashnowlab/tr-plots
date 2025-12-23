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
import plotly.graph_objects as go
import plotly.io as pio
import argparse
from typing import Optional
import math
import typing

# force Plotly to open figures in the system browser
pio.renderers.default = "browser"

# --- TEST MODE ---
TEST_MODE = True          # Quick testing: preview, limit work
TEST_LIMIT = 2            # Number of VCF records to process in test mode
SAVE_TEST_OUTPUTS = True  # If True, also save files when TEST_MODE is on

# --- File locations ---
from trplots.config import VCF_PATH, JSON_PATH, OUTPUT_BASE  

OUTPUT_BASE = Path(OUTPUT_BASE/"Plots/path_ref_motif_tandem_repeats_plots")

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
    # defensive defaults: fig.layout.width or x_span may be None when figure not fully rendered.
    if x0 == x1:
        x_label = x0
        xanchor = "center"
    else:
        x_label = x0
        xanchor = "left"

    # ensure numeric chart width and x_span
    try:
        chart_width_px = float(chart_width_px) if chart_width_px is not None else 900.0
    except Exception:
        chart_width_px = 900.0
    try:
        x_span = float(x_span) if (x_span is not None and not (isinstance(x_span, (list, tuple)))) else x_span
    except Exception:
        x_span = None
    # if x_span is a list/tuple (e.g. axis.range) extract extent, else default to 10 units
    if isinstance(x_span, (list, tuple)) and len(x_span) >= 2:
        try:
            x_span = float(x_span[1]) - float(x_span[0])
        except Exception:
            x_span = None
    if x_span is None or x_span <= 0:
        x_span = max(10.0, abs(x1 - x0) + 1.0)

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

def _longest_pure_repeat_count(sequence: Optional[str], motif: Optional[str]) -> float:
    """
    Return the longest contiguous run as repeat units (float/int) of `motif` found in `sequence`.
    - If sequence is present and motif present but no pure runs found -> return 0.0
    - If sequence missing or motif missing/unknown -> return NaN (cannot compute)
    """
    if sequence is None or motif is None:
        return float('nan')
    motif_clean = re.sub(r"\s+", "", str(motif)).upper()
    if not motif_clean or motif_clean == 'UNKNOWN':
        return float('nan')
    seq_upper = str(sequence).upper()

    # Build motif regex treating 'N' as single-base wildcard '.'
    motif_regex = ''.join('.' if ch == 'N' else re.escape(ch) for ch in motif_clean)
    try:
        pat = re.compile(f"(?:{motif_regex})+", flags=re.IGNORECASE)
    except re.error:
        return float('nan')

    best_units = 0
    for m in pat.finditer(seq_upper):
        run_bp = len(m.group(0))
        units = run_bp // len(motif_clean)
        if units > best_units:
            best_units = units

    # If sequence was valid but no run found, return 0.0 (user requested 0 instead of NaN)
    return float(best_units) if best_units > 0 else 0.0

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

    # Print test mode status
    if TEST_MODE:
        print(f"--- Test mode ON (limit={TEST_LIMIT}, save_test_outputs={SAVE_TEST_OUTPUTS}) ---")

    # --- Read allele spreadsheet generated by allele_spreadsheet.py and use its repeat-count logic ---
    spreadsheet_default = "/Users/annelisethorn/Documents/github/tr-plots/data/other_data/allele_spreadsheet.xlsx"
    sheet_path = Path(spreadsheet_default)
    if not sheet_path.exists():
        raise SystemExit(f"Allele spreadsheet not found: {sheet_path}")

    master_df = pd.read_excel(sheet_path, engine="openpyxl")
    # normalize cols (case-insensitive lookup)
    cols_lower = {c.lower(): c for c in master_df.columns}
    def col(*names):
        for n in names:
            if n.lower() in cols_lower:
                return cols_lower[n.lower()]
        return None
    chrom_col = col("Chromosome", "Chr", "chrom")
    pos_col = col("Position (Start)", "Position", "Pos", "position")
    gene_col = col("Gene")
    disease_col = col("Disease")
    motif_col = col("Motif")
    repeat_col = col("Repeat Count", "RepeatCount", "repeat_count")
    allele_seq_col = col("Allele Sequence", "Allele_Sequence", "allele_sequence")

    if not chrom_col or not pos_col:
        raise SystemExit("Spreadsheet must contain Chromosome and Position columns")

    master_df[chrom_col] = master_df[chrom_col].astype(str).apply(norm_chrom)

    # build interval index from loci_data for fast lookup
    import bisect, time
    loci_index = {}
    for entry in loci_data:
        ch = norm_chrom(entry.get("chrom", ""))
        s = int(entry.get("start_hg38", 0))
        e = int(entry.get("stop_hg38", 0))
        loci_index.setdefault(ch, []).append((s, e, entry))
    for ch, lst in loci_index.items():
        lst.sort(key=lambda x: x[0])
        starts = [t[0] for t in lst]
        loci_index[ch] = (starts, lst)
    def find_locus(ch, pos):
        if ch not in loci_index:
            return None
        starts, lst = loci_index[ch]
        idx = bisect.bisect_right(starts, pos)
        for j in range(idx-1, -1, -1):
            s,e,entry = lst[j]
            if e >= pos:
                return entry
        return None

    total_loci = 0
    matched_vcf_loci = 0
    skipped_loci = []
    previewed_count = 0
    grouped = master_df.groupby([chrom_col, pos_col], sort=False)
    start_time = time.perf_counter()
    for (chrom, pos), group in grouped:
        total_loci += 1
        try:
            pos_int = int(pos)
        except Exception:
            print(f"Skipping invalid position: {chrom}:{pos}")
            continue

        matched_locus = find_locus(chrom, pos_int)
        if not matched_locus:
            skipped_loci.append(f"{chrom}:{pos_int}")
            continue
        matched_vcf_loci += 1

        gene = matched_locus.get("gene", group[gene_col].iloc[0] if gene_col and gene_col in group else "Unknown")
        disease = matched_locus.get("disease", group[disease_col].iloc[0] if disease_col and disease_col in group else "Unknown")

        # Print preview in test mode
        if TEST_MODE:
            print(f"Previewing: {gene} / {disease} / {chrom}:{pos_int}")


        # assemble base_df with repeat counts, motifs, and allele sequence where available
        base_rows = []
        for _, r in group.iterrows():
            allele_seq = None
            if allele_seq_col and allele_seq_col in r and pd.notna(r[allele_seq_col]):
                allele_seq = r[allele_seq_col]

            motif_val = r[motif_col] if motif_col and motif_col in r and pd.notna(r[motif_col]) else None
            motif_clean = re.sub(r"\s+", "", str(motif_val)).upper() if motif_val is not None else None
            if motif_clean == "UNKNOWN":
                motif_clean = None

            repeat_val = None
            if repeat_col and repeat_col in r and pd.notna(r[repeat_col]):
                repeat_val = r[repeat_col]
            elif allele_seq is not None and motif_clean:
                repeat_val = _longest_pure_repeat_count(allele_seq, motif_clean)

            if repeat_val is None and allele_seq is None:
                continue

            base_rows.append({
                "Repeat Count": repeat_val,
                "Motif": motif_clean,
                "Allele Sequence": allele_seq
            })

        if not base_rows:
            print(f"Skipping {chrom}:{pos_int} — unable to derive Repeat Count for any allele rows")
            continue
        base_df = pd.DataFrame(base_rows)

        # filter to pathogenic/reference motifs from JSON if present
        raw_motif_field = (
            matched_locus.get("pathogenic_motif_reference_orientation")
            or matched_locus.get("pathogenic_motif_gene_orientation")
            or matched_locus.get("reference_motif_reference_orientation")
        )
        motifs = []
        if isinstance(raw_motif_field, list):
            for m in raw_motif_field:
                if isinstance(m, str) and m.strip():
                    motif_clean = re.sub(r"\s+", "", m).upper()
                    if motif_clean and motif_clean.upper() != "UNKNOWN":
                        motifs.append(motif_clean)
        elif isinstance(raw_motif_field, str) and raw_motif_field.strip():
            motif_clean = re.sub(r"\s+", "", raw_motif_field).upper()
            if motif_clean and motif_clean.upper() != "UNKNOWN":
                motifs.append(motif_clean)
        seen = set()
        motifs = [m for m in motifs if not (m in seen or seen.add(m))]
        
        # KEEP ALL motifs (do NOT pick just one)
        if not motifs:
            print(f"Skipping {chrom}:{pos_int} ({gene}) — no valid pathogenic reference motifs")
            continue

        # Build plot_df so every pathogenic motif from JSON appears (color distinguishes motif).
        expanded_rows = []
        for _, r in base_df.iterrows():
            allele_seq = r.get("Allele Sequence")
            base_repeat = r.get("Repeat Count")
            base_motif = r.get("Motif")

            for motif in motifs:
                if allele_seq is not None and pd.notna(allele_seq):
                    rc = _longest_pure_repeat_count(allele_seq, motif)
                elif base_motif == motif and pd.notna(base_repeat):
                    rc = base_repeat
                else:
                    rc = float("nan")
                expanded_rows.append({"Repeat Count": rc, "Motif": motif})

        plot_df = pd.DataFrame(expanded_rows)
        plot_df = plot_df.dropna(subset=["Repeat Count"])
        if plot_df.empty:
            print(f"Skipping {chrom}:{pos_int} ({gene}) — no usable repeat counts after filtering motifs")
            continue

        # histogram and plotting (reuse existing layout logic)
        min_bin = int(np.floor(plot_df["Repeat Count"].min()))
        max_bin = int(np.ceil(plot_df["Repeat Count"].max()))
        nbins = max(1, max_bin - min_bin + 1)
        subtitle = f"{chrom} - {pos_int} | {gene} - {title_outside_parens(str(disease))}"
        subtitle_html = "<br>".join(f"<span style='font-size:12px; color:black'>{l}</span>" for l in [subtitle])

        fig = px.histogram(
            plot_df, x="Repeat Count", color="Motif",
            nbins=nbins,
            title=(
                "<span style='font-size:18px; font-weight:bold; color:black'>"
                "Allele Size Distribution</span><br>" + subtitle_html
            ),
            color_discrete_sequence=px.colors.qualitative.Plotly
        )
        fig.update_layout(barmode='stack', bargap=0.1)
        fig.update_traces(opacity=1, marker_line_width=0)

        arr = np.asarray(plot_df["Repeat Count"].values, dtype=float)
        median = float(np.median(arr))
        p95 = float(np.quantile(arr, 0.95))
        x_right = max(2 * median + 0.5, p95 + 0.5, 5.5)
        fig.update_xaxes(range=[0, x_right])

        # obtain ranges from matched_locus
        benign_min = matched_locus.get("benign_min")
        benign_max = matched_locus.get("benign_max")
        intermediate_min = matched_locus.get("intermediate_min")
        intermediate_max = matched_locus.get("intermediate_max")
        pathogenic_min = matched_locus.get("pathogenic_min")
        pathogenic_max = matched_locus.get("pathogenic_max")

        x_range = fig.layout.xaxis.range
        x_span = (x_range[1] - x_range[0]) if x_range and None not in x_range else ((max_bin + 2.5) - 0)
        label_positions = {}
        if benign_min is not None and benign_max is not None:
            add_range_marker_or_line(fig, benign_min, benign_max, "Benign", label_positions, chart_width_px=fig.layout.width, x_span=x_span, level=0)
        has_intermediate = (intermediate_min is not None and intermediate_max is not None)
        if has_intermediate:
            add_range_marker_or_line(fig, intermediate_min, intermediate_max, "Intermediate", label_positions, chart_width_px=fig.layout.width, x_span=x_span, level=1)
        pathogenic_level = 1 if not has_intermediate else 2
        if gene == "VWA1":
            add_range_marker_or_line(fig, 1, 1, "Pathogenic", label_positions, chart_width_px=fig.layout.width, x_span=x_span, level=pathogenic_level)
            add_range_marker_or_line(fig, 3, 3, "Pathogenic", label_positions, chart_width_px=fig.layout.width, x_span=x_span, level=pathogenic_level)
        elif pathogenic_min is not None and pathogenic_max is not None:
            add_range_marker_or_line(fig, pathogenic_min, pathogenic_max, "Pathogenic", label_positions, chart_width_px=fig.layout.width, x_span=x_span, level=pathogenic_level)

        fig.update_layout(
            width=900, height=500, margin=dict(t=125),
            xaxis_title="Repeat Count", yaxis_title="Allele Count",
            plot_bgcolor='white'
        )
        fig.update_xaxes(showline=True, linecolor='black', zeroline=False, showgrid=False, ticks="outside")
        fig.update_yaxes(showline=True, linecolor='black', zeroline=False, showgrid=True, gridcolor='rgba(0,0,0,0.08)', ticks="outside")

        # --- Save / Preview ---
        safe_gene = re.sub(r'[^A-Za-z0-9_-]+', '_', str(gene) or "gene")[:50]
        html_path = os.path.join(OUTPUT_HTML_DIR, f"{safe_gene}_{chrom}_{pos_int}_allele_dist.html")
        png_path  = os.path.join(OUTPUT_PNG_DIR,  f"{safe_gene}_{chrom}_{pos_int}_allele_dist.png")

        if TEST_MODE:
            try:
                fig.show(renderer="browser")
            except Exception:
                pass
            previewed_count += 1
            if SAVE_TEST_OUTPUTS:
                fig.write_html(html_path)
                orig_title = None
                try:
                    if getattr(fig.layout, "title", None) and getattr(fig.layout.title, "text", None):
                        orig_title = fig.layout.title.text
                        fig.update_layout(title_text=f"{safe_gene}_{chrom}_{pos_int}")
                    fig.write_image(png_path, width=900, height=500)
                except Exception as e:
                    print(f"PNG export failed for {png_path}: {e} (HTML saved)")
                finally:
                    if orig_title is not None:
                        fig.update_layout(title_text=orig_title)
            if previewed_count >= TEST_LIMIT:
                break
        else:
            fig.write_html(html_path)
            orig_title = None
            try:
                if getattr(fig.layout, "title", None) and getattr(fig.layout.title, "text", None):
                    orig_title = fig.layout.title.text
                    fig.update_layout(title_text=f"{safe_gene}_{chrom}_{pos_int}")
                fig.write_image(png_path, width=900, height=500)
            except Exception as e:
                print(f"PNG export failed for {png_path}: {e} (HTML saved)")
            finally:
                if orig_title is not None:
                    fig.update_layout(title_text=orig_title)

    # --- Finished ---
    elapsed = time.perf_counter() - start_time
    print(f"Total loci in spreadsheet: {total_loci}, matched loci in JSON: {matched_vcf_loci}, elapsed {elapsed:.1f}s")
    if skipped_loci:
        print(f"Skipped {len(skipped_loci)} loci not found in JSON (example first 200):")
        for s in sorted(skipped_loci)[:200]:
            print("  " + s)
    return

if __name__ == "__main__":
    args = parse_args()
    main(args)
