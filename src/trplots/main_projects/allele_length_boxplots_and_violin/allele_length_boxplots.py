# ---------------------------------------------
#  Script: Allele Length Boxplot Generator (Locus-based, All Alleles)
# ---------------------------------------------
import re
import ast
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import argparse
from pathlib import Path

pio.renderers.default = "browser"

BASE_DIR = Path("/Users/annelisethorn/Documents/github/tr-plots")

# Default to the integrated spreadsheet produced by allele_spreadsheet.py
DATA_PATH = BASE_DIR / "data" / "other_data" / "allele_spreadsheet.xlsx"
SHEET_NAME = "Integrated Alleles"

# --- Output directories (subfolders for PNG and HTML) ---
OUTPUT_ROOT = (
    Path("/Users/annelisethorn/Documents/github/tr-plots")
    / "results" / "plots" / "allele_length_boxplots"
)
OUTPUT_DIR_PNG  = OUTPUT_ROOT / "png"
OUTPUT_DIR_HTML = OUTPUT_ROOT / "html"
TEST_OUTPUT_DIR = OUTPUT_ROOT / "test_outputs"

# Ensure the standard output dirs exist up front
OUTPUT_DIR_PNG.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR_HTML.mkdir(parents=True, exist_ok=True)
# test_outputs will be created on demand (only if we actually save test outputs)

def ENSURE_DIR(*parts):
    p = BASE_DIR.joinpath(*parts)
    p.mkdir(parents=True, exist_ok=True)
    return p

TEST_MODE_DEFAULT = True
TEST_LIMIT_DEFAULT = 2
SAVE_TEST_OUTPUTS_DEFAULT = True

def parse_args():
    p = argparse.ArgumentParser(description="Allele length boxplot generator")
    p.add_argument("--test", dest="test", action="store_true", help="Enable test mode")
    p.add_argument("--no-test", dest="test", action="store_false", help="Disable test mode")
    p.set_defaults(test=TEST_MODE_DEFAULT)
    p.add_argument("--test-limit", type=int, default=TEST_LIMIT_DEFAULT)
    p.add_argument("--save-test-outputs", action="store_true", default=SAVE_TEST_OUTPUTS_DEFAULT,
                   help="In test mode, also save PNG/HTML to test_outputs/")
    p.add_argument("--data-path", type=str, default=str(DATA_PATH))
    p.add_argument("--sheet", type=str, default=SHEET_NAME, help="Excel sheet name")
    p.add_argument("--output-dir", type=str, default=None,
                   help="If provided, write all outputs (PNG/HTML) here instead of the default png/html dirs")
    return p.parse_args()

FIG_WIDTH = 900
FIG_HEIGHT = 500
TOP_MARGIN = 130
PNG_SCALE = 2

POP_COLOR = {
    'All': '#ff99cc',   # pink for aggregates
    'EUR': '#1f77b4',   # blue
    'EAS': '#2ca02c',   # green
    'SAS': '#9467bd',   # purple
    'AMR': '#d62728',   # red
    'AFR': '#ff7f0e',   # orange/yellow
    'Unknown': '#7f7f7f',   # gray for for unknowns
}
SUPERPOP_ORDER = ['All', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'Unknown']

def _norm_inh(val):
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

def _make_locus_vec(df: pd.DataFrame) -> pd.Series:
    chrom = df.get('Chromosome')
    pos = df.get('Position')
    disease_id = df.get('Disease ID')

    chrom_s = chrom.astype(str).str.strip().fillna("")
    pos_s = pos.astype(str).str.strip().fillna("")
    did_s = disease_id.astype(str).str.strip().fillna("")

    chrom_ok = (~chrom_s.isin(["", "nan", "None"]))
    pos_ok = (~pos_s.isin(["", "nan", "None"]))

    chrom_str = np.where(
        chrom_s.str.lower().str.startswith("chr"),
        chrom_s,
        "Chr" + chrom_s
    )
    def _coerce_pos(x):
        try:
            return str(int(float(x)))
        except Exception:
            return x
    pos_str = pos_s.map(_coerce_pos)

    locus_from_cp = chrom_str + ":" + pos_str
    locus = np.where(chrom_ok & pos_ok, locus_from_cp,
                     np.where(~did_s.isin(["", "nan"]), did_s,
                              (df.get('Gene', "").astype(str).str.strip().fillna("") + "|" +
                               df.get('Disease', "").astype(str).str.strip().fillna("") + "|" +
                               chrom_s + ":" + pos_s)))
    return pd.Series(locus, index=df.index)

def create_boxplot_all(locus_data, locus_metadata, show_no_data_note=True, show_points=False):
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
        lower_whisk = x[x >= lf].min() if np.any(x >= lf) else x.min()
        upper_whisk = x[x <= uf].max() if np.any(x <= uf) else x.max()
        outliers = x[(x < lf) | (x > uf)]
        return q1, med, q3, lower_whisk, upper_whisk, outliers

    gene = locus_metadata['gene']
    disease = locus_metadata['disease']
    locus = locus_metadata['locus']
    inheritance = locus_metadata['inheritance']
    counts = locus_metadata['counts']
    ordered_categories = locus_metadata['ordered_categories']

    fig = go.Figure()
    locus_groups = locus_data.groupby('SuperPop', observed=True)['Allele length']

    # Keep 'Unknown' category present to ensure a visible bar even when empty

    for pop in ordered_categories:
        if pop in locus_groups.groups:
            xs = locus_groups.get_group(pop).to_numpy(dtype=float)
            xs = xs[~np.isnan(xs)]
        else:
            xs = np.array([])

        if xs.size == 0:
            # No data for this category: do not add a box trace.
            # The y-axis already lists categories via categoryarray/tickvals,
            # so the label will remain visible without a bar.
            continue

        stats = five_number_summary(xs)
        if stats is None:
            continue
        q1, med, q3, lw, uw, out = stats

        fig.add_trace(go.Box(
            orientation="h", name="", y=[pop],
            q1=[q1], median=[med], q3=[q3],
            lowerfence=[lw], upperfence=[uw],
            marker_color=POP_COLOR.get(pop, "#7f7f7f"),
            boxpoints=False, showlegend=False, hoveron="boxes",
            hovertemplate=None
        ))

        if show_points and out.size:
            fig.add_trace(go.Scatter(
                x=out, y=[pop] * out.size, mode="markers", name="",
                marker=dict(size=6, color=POP_COLOR.get(pop, "#7f7f7f")),
                showlegend=False,
                hovertemplate="Allele length: %{x}<br>Population: %{y}<extra></extra>",
            ))

    fig.update_yaxes(
        categoryorder='array', categoryarray=ordered_categories,
        tickmode='array', tickvals=ordered_categories,
        autorange='reversed', title=""
    )

    pop_desc = ', '.join(f"{pop}: {counts.get(pop, 0)}" for pop in ordered_categories)
    pop_lines = _wrap_to_lines(pop_desc, max_len=110)

    # Do not add "no data" annotations; keep axis labels without bars for empty categories

    fig.update_layout(
        hovermode="closest", width=FIG_WIDTH, height=FIG_HEIGHT,
        margin=dict(t=TOP_MARGIN, r=40, b=60, l=90),
        autosize=False, plot_bgcolor='white',
        showlegend=False, font=dict(color='black'),
        xaxis=dict(title="Allele Length (bp)", ticks='outside', showline=True, linecolor='black'),
        yaxis=dict(ticks='outside', showline=True, linecolor='black'),
    )

    main_title = "<b>Allele Lengths per Population</b>"
    subtitle_lines = [f"{gene} - {disease}", f"{locus}"]
    if pop_lines:
        subtitle_lines.append(f"Total Individuals per Population: {pop_lines[0]}")
        subtitle_lines.extend(pop_lines[1:])
    subtitle_lines.append(f"Inheritance: {inheritance}")

    annos = [dict(text=main_title, x=0, xref="paper", xanchor="left",
                  y=1.42, yref="paper", yanchor="top",
                  showarrow=False, align="left", font=dict(size=18))]
    y0 = 1.34
    for i, line in enumerate(subtitle_lines):
        annos.append(dict(
            text=f"<span style='font-size:12px'>{line}</span>",
            x=0, xref="paper", xanchor="left",
            y=y0 - 0.06*i, yref="paper", yanchor="top",
            showarrow=False, align="left"
        ))
    fig.update_layout(annotations=annos)
    return fig

def main(args=None):
    args = parse_args()

    test_mode = args.test
    test_limit = args.test_limit
    save_test_outputs = args.save_test_outputs
    data_path = Path(args.data_path)
    sheet_name = args.sheet

    # Default output destinations
    output_dir_png = OUTPUT_DIR_PNG
    output_dir_html = OUTPUT_DIR_HTML

    # If user forces a single output dir, we honor it
    if args.output_dir:
        outdir = Path(args.output_dir); outdir.mkdir(parents=True, exist_ok=True)
        output_dir_png = outdir
        output_dir_html = outdir

    # In test mode: preview, and optionally save into test_outputs/
    if test_mode:
        print(f"--- Test mode ON (limit={test_limit}, save_test_outputs={save_test_outputs}) ---")
        test_save_dir = None
        if save_test_outputs:
            TEST_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            test_save_dir = TEST_OUTPUT_DIR

    needed_cols = [
        'Gene','Disease','Disease ID','Chromosome','Position',
        'SuperPop','Sample ID Cleaned','Allele length','Inheritance'
    ]
    try:
        df = pd.read_excel(
            data_path, sheet_name=sheet_name,
            usecols=lambda c: True if c in needed_cols else False,
            engine="openpyxl"
        )
    except FileNotFoundError:
        print(f"Error: Input data file not found at: {data_path}"); return
    except ValueError as e:
        print(f"Error reading sheet '{sheet_name}': {e}"); return

    if 'Allele length' not in df.columns:
        raise ValueError("Missing 'Allele length' column")
    df['Allele length'] = pd.to_numeric(df['Allele length'], errors='coerce')
    df.loc[~np.isfinite(df['Allele length']), 'Allele length'] = np.nan
    df.loc[df['Allele length'] <= 0, 'Allele length'] = np.nan

    df['Inheritance_norm'] = df['Inheritance'].map(_norm_inh) if 'Inheritance' in df.columns else 'Unknown'

    df['Locus'] = _make_locus_vec(df)

    if 'SuperPop' not in df.columns:
        raise ValueError("Missing 'SuperPop' column in input data.")
    df['SuperPop'] = pd.Categorical(df['SuperPop'], categories=SUPERPOP_ORDER, ordered=True)

    # EARLY ROW FILTERING: keep only rows that can produce a plot
    keep_mask = (
        df['Allele length'].notna() &
        df['SuperPop'].notna() &
        df['Locus'].astype(str).str.len().gt(0)
    )
    df = df.loc[keep_mask].copy()

    df_min = df[['Gene','Disease','Locus','SuperPop','Sample ID Cleaned','Allele length','Inheritance_norm']]

    df_all = df_min.copy()
    df_all['SuperPop'] = 'All'
    df_all['SuperPop'] = pd.Categorical(df_all['SuperPop'], categories=SUPERPOP_ORDER, ordered=True)
    df_min = pd.concat([df_min, df_all], ignore_index=True)

    for required in ['Gene','Disease','Locus','SuperPop','Sample ID Cleaned','Allele length','Inheritance_norm']:
        if required not in df_min.columns:
            raise ValueError(f"Missing required column: '{required}'")

    locus_data_map = {}
    for (gene, disease, locus), group_df in df_min.groupby(['Gene','Disease','Locus'], observed=True):
        counts_series = (
            group_df.groupby('SuperPop', observed=True)['Sample ID Cleaned']
                    .nunique(dropna=True)
        ).reindex(SUPERPOP_ORDER).fillna(0).astype(int)

        inheritance = group_df['Inheritance_norm'].iloc[0]

        locus_data_map[(gene, disease, locus)] = {
            'data': group_df[['Allele length','SuperPop']].copy(),
            'metadata': {
                'gene': gene,
                'disease': disease,
                'locus': locus,
                'inheritance': inheritance,
                'counts': counts_series.to_dict(),
                'ordered_categories': SUPERPOP_ORDER
            }
        }

    made = 0
    for (_, data_entry) in locus_data_map.items():
        meta = data_entry['metadata']
        gene, disease, locus = meta['gene'], meta['disease'], meta['locus']

        if data_entry['data']['Allele length'].dropna().empty:
            continue

        fig = create_boxplot_all(
            data_entry['data'], meta,
            show_no_data_note=True, show_points=False
        )

        safe = lambda s: re.sub(r'[\\/:*?"<>|]+', '_', str(s))
        png_name  = f"{safe(gene)}_{safe(disease)}_{safe(locus)[:120]}_allele_length_boxplot.png"
        html_name = f"{safe(gene)}_{safe(disease)}_{safe(locus)[:120]}_allele_length_boxplot.html"

        if test_mode:
            # Always preview in test mode
            print(f"Previewing: {gene} / {disease} / {locus}")
            fig.show()

            # Only save when --save-test-outputs is enabled
            if save_test_outputs:
                try: fig.write_html(str((TEST_OUTPUT_DIR / html_name)))
                except Exception as e: print(f"Failed to write HTML {TEST_OUTPUT_DIR / html_name}: {e}")
                try: fig.write_image(str((TEST_OUTPUT_DIR / png_name)),
                                     width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
                except Exception as e: print(f"Failed to write PNG {TEST_OUTPUT_DIR / png_name}: {e} — ensure 'kaleido' is installed")

            made += 1
            if made >= test_limit:
                break

        else:
            # Normal (non-test) run: write to png/ and html/ (or --output-dir if given)
            try: fig.write_html(str((output_dir_html / html_name)))
            except Exception as e: print(f"Failed to write HTML {output_dir_html / html_name}: {e}")
            try: fig.write_image(str((output_dir_png / png_name)),
                                 width=FIG_WIDTH, height=FIG_HEIGHT, scale=PNG_SCALE)
            except Exception as e: print(f"Failed to write PNG {output_dir_png / png_name}: {e} — ensure 'kaleido' is installed")

    print(f"--- {'Test mode ON: Test completed' if test_mode else 'Done'} ---")

if __name__ == "__main__":
    main()
