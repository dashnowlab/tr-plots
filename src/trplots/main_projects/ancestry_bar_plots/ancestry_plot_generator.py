"""
---------------------------------------------
 Script: Ancestry Bar Plot Generator
 Purpose:
   - Reads a dataset containing ancestry and genotype information for multiple samples
   - Cleans and standardizes metadata (population labels, sex, inheritance, pore)
   - Counts pathogenic individuals under AD/AR/XD/XR models
   - Computes percentages and 95% binomial confidence intervals (Wilson)
   - Builds horizontal bar plots with dropdowns for Pore type and Sex
   - Saves plots as interactive HTML (and shows a preview in test mode)
---------------------------------------------
"""

import os
import re
import ast
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from statsmodels.stats import proportion

pio.renderers.default = "browser"

# ----------------------------
# Defaults
# ----------------------------
DEFAULT_TEST_MODE = True
DEFAULT_TEST_LIMIT = 2
DEFAULT_SAVE_TEST_OUTPUTS = True

# Repository base directory (go up from this file to repo root)
BASE_DIR = Path(__file__).resolve().parents[4]

# Default to the integrated spreadsheet produced by allele_spreadsheet.py
DATA_PATH = BASE_DIR / "data" / "other_data" / "allele_spreadsheet.xlsx"
SHEET_NAME = "Integrated Alleles"

# --- Figure sizing (standardized across all main_projects)
FIG_WIDTH = 900
FIG_HEIGHT = 500

# --- Superpopulation palette/order (shared across ancestry-aware plots)
POP_COLOR = {
    "All": "#ff99cc",   # pink for aggregates
    "EUR": "#1f77b4",   # blue
    "EAS": "#2ca02c",   # green
    "SAS": "#9467bd",   # purple
    "AMR": "#d62728",   # red
    "AFR": "#ff7f0e",   # orange/yellow
    "Unknown": "#7f7f7f",   # gray for for unknowns
}
SUPERPOP_ORDER_BASE = ["All", "AFR", "AMR", "EAS", "EUR", "SAS", "Unknown"]

# If you have trplots.config available, keep this:
import sys
sys.path.append(str(BASE_DIR / "src"))
from trplots.config import OUTPUT_BASE
OUTPUT_BASE = Path(OUTPUT_BASE) / "plots" / "ancestry_plots"

# ----------------------------
# CLI args
# ----------------------------
def parse_args():
    """
    Define and parse CLI arguments.
    Flags:
    - --test / --no-test: enable/disable test mode
    - --test-limit: number of (gene,disease) plots to preview in test mode
    - --save-test-outputs: save files while in test mode
    - --allele-spreadsheet: path to the integrated allele spreadsheet
    - --output-dir: override base output directory
    """
    p = argparse.ArgumentParser(description="Ancestry bar plot generator")
    p.add_argument("--test", dest="test", action="store_true", help="Enable test mode")
    p.add_argument("--no-test", dest="test", action="store_false", help="Disable test mode")
    p.set_defaults(test=DEFAULT_TEST_MODE)

    p.add_argument("--test-limit", dest="test_limit", type=int, default=DEFAULT_TEST_LIMIT,
                   help="How many (gene,disease) plots to generate in test mode")
    p.add_argument("--save-test-outputs", dest="save_test_outputs", action="store_true",
                   default=DEFAULT_SAVE_TEST_OUTPUTS, help="Save outputs even in test mode")

    p.add_argument("--allele-spreadsheet", dest="allele_spreadsheet", type=str,
                   default=str(DATA_PATH), help="Path to allele spreadsheet")

    p.add_argument("--output-dir", dest="output_dir", type=str, default=None,
                   help="Override output directory (base).")
    return p.parse_args()


# ----------------------------
# Helpers: column lookup + cleaning
# ----------------------------
def build_col_lookup(df: pd.DataFrame):
    """
    Build a case-insensitive column resolver for flexible spreadsheets.
    Returns a function `col(*names)` that returns the exact column name
    from the DataFrame if any of the provided aliases are present.
    """
    cols_lower = {c.lower().strip(): c for c in df.columns}
    def col(*names):
        for n in names:
            key = n.lower().strip()
            if key in cols_lower:
                return cols_lower[key]
        return None
    return col

def clean_inheritance(x):
    """
    Normalize inheritance strings, handling list-like strings such as
    "['AD', 'AR']" by selecting the first entry. Returns a single label
    (e.g., 'AD') or NaN on failure.
    """
    # Handles "['AD']" or "['AD', 'AR']" etc
    if isinstance(x, str):
        s = x.strip()
        if s.startswith("[") and s.endswith("]"):
            try:
                val = ast.literal_eval(s)
                if isinstance(val, (list, tuple)) and len(val) > 0:
                    return str(val[0]).strip()
            except Exception:
                return np.nan
        return s
    return x

def normalize_sex(x):
    """
    Coerce free-text sex values into 'Male', 'Female', or 'Unknown'.
    Accepts common abbreviations (m/f) and missing values.
    """
    if pd.isna(x):
        return "Unknown"
    s = str(x).strip().lower()
    sex_map = {
        "m": "Male", "male": "Male",
        "f": "Female", "female": "Female",
        "u": "Unknown", "unk": "Unknown", "unknown": "Unknown", "": "Unknown", "none": "Unknown", "nan": "Unknown"
    }
    return sex_map.get(s, str(x).strip() if str(x).strip() else "Unknown")

def normalize_pore(x):
    """
    Normalize sequencing pore identifiers into 'R9', 'R10', or 'All Types'.
    Leaves other values as-is for future bucketing.
    """
    # Force into R9 / R10 / All Types
    if pd.isna(x):
        return "All Types"
    s = str(x).strip().upper()
    if "R10" in s:
        return "R10"
    if "R9" in s:
        return "R9"
    if s in {"ALL", "ALL TYPES", "ALLTYPES"}:
        return "All Types"
    # If it's some other value, keep it but you may want to bucket it:
    return s

def normalize_population(x):
    """
    Try to coerce messy ancestry text into superpop-like labels.
    Keeps exact codes (AFR/AMR/EAS/EUR/SAS/Unknown/All) if already present.
    """
    if pd.isna(x):
        return "Unknown"
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return "Unknown"

    upper = s.upper()

    # Already canonical
    canonical = {"AFR", "AMR", "EAS", "EUR", "SAS", "UNKNOWN", "ALL"}
    if upper in canonical:
        return "Unknown" if upper == "UNKNOWN" else ("All" if upper == "ALL" else upper)

    # Common patterns
    # (If your spreadsheet has full strings like "African", "European", etc.)
    if "AFR" in upper or "AFRIC" in upper:
        return "AFR"
    if "AMR" in upper or "AMERI" in upper or "LATIN" in upper:
        return "AMR"
    if "EAS" in upper or "EAST ASIA" in upper or "CHINESE" in upper or "JAPAN" in upper or "KORE" in upper:
        return "EAS"
    if "EUR" in upper or "EURO" in upper or "WHITE" in upper:
        return "EUR"
    if "SAS" in upper or "SOUTH ASIA" in upper or "INDIA" in upper or "PAKIST" in upper or "BANGLA" in upper:
        return "SAS"

    # If no pattern matches, treat as Unknown
    return "Unknown"


def binomial_ci(x, n, confidence=0.95):
    """
    Compute Wilson score interval for `x` successes out of `n` trials.
    Returns (lower%, upper%) as floats in percent space.
    """
    if n <= 0:
        return (0.0, 0.0)
    lower, upper = proportion.proportion_confint(x, n, alpha=1-confidence, method="wilson")
    return lower * 100, upper * 100


# ----------------------------
# Core summarization (per gene/disease, with pore+sex selection)
# ----------------------------
def compute_group_summary(df_gene: pd.DataFrame, inheritance: str, superpop_order: list[str]) -> pd.DataFrame:
    """
    Returns a table with rows = populations (including 'All'),
    columns = total_count, affected_count, percentage, CI bounds, error bars.
    """
    inh = inheritance.upper()
    min_alleles = 1 if inh in {"AD", "XD"} else 2

    # Total individuals per pop (ALL alleles, regardless of pathogenic)
    totals = (
        df_gene.groupby("Population")["Sample ID"]
        .nunique()
        .reset_index(name="total_count")
    )

    # Affected individuals per pop for this inheritance model
    df_inh = df_gene[df_gene["Inheritance"] == inh]
    df_path = df_inh[df_inh["Is pathogenic"] == True]

    if df_path.empty:
        affected = pd.DataFrame({"Population": [], "affected": []})
    else:
        per_sample = (
            df_path.groupby(["Population", "Sample ID"])
            .size()
            .reset_index(name="path_alleles")
        )
        affected_samples = per_sample[per_sample["path_alleles"] >= min_alleles]
        affected = (
            affected_samples.groupby("Population")["Sample ID"]
            .nunique()
            .reset_index(name="affected")
        )

    out = totals.merge(affected, on="Population", how="left")
    out["affected"] = out["affected"].fillna(0).astype(int)

    # Add 'All' row freshly (sum across pops)
    all_row = pd.DataFrame({
        "Population": ["All"],
        "total_count": [int(out["total_count"].sum())],
        "affected": [int(out["affected"].sum())],
    })
    out = pd.concat([all_row, out], ignore_index=True)

    # Ensure all populations appear, in fixed order
    out = (
        out.set_index("Population")
        .reindex(superpop_order)
        .reset_index()
    )
    out["total_count"] = out["total_count"].fillna(0).astype(int)
    out["affected"] = out["affected"].fillna(0).astype(int)

    # Percent + CI
    with np.errstate(divide="ignore", invalid="ignore"):
        out["percentage"] = np.where(
            out["total_count"] > 0,
            out["affected"] / out["total_count"] * 100,
            0.0
        )

    ci = out.apply(lambda r: binomial_ci(int(r["affected"]), int(r["total_count"])), axis=1)
    out["ci_lower"] = [lo for lo, hi in ci]
    out["ci_upper"] = [hi for lo, hi in ci]
    out["error_plus"] = out["ci_upper"] - out["percentage"]
    out["error_minus"] = out["percentage"] - out["ci_lower"]

    # Display formatting
    out["percentage_display"] = out["percentage"].map(
        lambda x: f"{int(x)}" if float(x).is_integer() else f"{x:.2f}"
    )
    return out


def build_subtitle_counts(summary_df: pd.DataFrame) -> str:
    """
    Build a compact HTML string showing total individuals per population.
    Line-wraps to avoid overly wide titles in the figure.
    """
    parts = [f"{row['Population']}: {int(row['total_count'])}" for _, row in summary_df.iterrows()]
    # Wrap roughly every ~100 chars
    lines, cur = [], ""
    for p in parts:
        if len(cur) + len(p) + 2 > 100:
            if cur:
                lines.append(cur.rstrip(", "))
            cur = ""
        cur += p + ", "
    if cur:
        lines.append(cur.rstrip(", "))
    return "<br>".join([f"<span style='font-size:12px'>{l}</span>" for l in lines])


def create_horizontal_bar_plot(df_raw: pd.DataFrame, gene: str, disease: str, superpop_order: list[str]):
        """
        Create a horizontal bar plot of pathogenic genotype percentages
        per population for a given `gene`/`disease` pair.
        - Uses dropdown to switch Pore type (All Types / R9 / R10)
        - Y-axis category order enforces presence of 'Unknown' label even
            when counts are zero (no visible bar is drawn for zero values).
        Returns a Plotly Figure or None if data is insufficient.
        """
    df_gd = df_raw[(df_raw["Gene"] == gene) & (df_raw["Disease"] == disease)].copy()
    if df_gd.empty:
        return None

    # Pick a single inheritance label for this plot (most common non-null)
    inh_counts = (
        df_gd["Inheritance"].dropna().astype(str).str.upper()
        .value_counts()
    )
    if inh_counts.empty:
        return None

    inheritance = inh_counts.index[0]
    if inheritance not in {"AD", "AR", "XD", "XR"}:
        # skip weird/unknown inheritance modes
        return None

    pore_order = ["All Types", "R9", "R10"]

    trace_map = {}
    subtitle_map = {}

    for pore in pore_order:
        df_sel = df_gd.copy()

        # Normalize pore field before filtering
        df_sel["Pore"] = df_sel["Pore"].apply(normalize_pore)

        if pore != "All Types":
            df_sel = df_sel[df_sel["Pore"] == pore]

        # If empty, still build a zero summary to keep axes consistent
        if df_sel.empty:
            summary = pd.DataFrame({
                "Population": superpop_order,
                "total_count": [0] * len(superpop_order),
                "affected": [0] * len(superpop_order),
                "percentage": [0.0] * len(superpop_order),
                "ci_lower": [0.0] * len(superpop_order),
                "ci_upper": [0.0] * len(superpop_order),
                "error_plus": [0.0] * len(superpop_order),
                "error_minus": [0.0] * len(superpop_order),
                "percentage_display": ["0"] * len(superpop_order),
            })
        else:
            summary = compute_group_summary(df_sel, inheritance, superpop_order)

        trace_map[pore] = summary
        subtitle_map[pore] = build_subtitle_counts(summary)

    fig = go.Figure()
    trace_keys = list(trace_map.keys())

    for pore in trace_keys:
        summary = trace_map[pore]
        fig.add_bar(
            x=summary["percentage"],
            y=summary["Population"],
            orientation="h",
            error_x=dict(
                array=summary["error_plus"],
                arrayminus=summary["error_minus"],
                thickness=1
            ),
            marker_color="lightblue",
            customdata=np.stack([
                summary["affected"],
                summary["total_count"],
                summary["percentage_display"],
                summary["ci_lower"],
                summary["ci_upper"],
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

    default_key = "All Types"
    if default_key not in trace_keys:
        default_key = trace_keys[0]
    default_idx = trace_keys.index(default_key)
    fig.data[default_idx].visible = True

    buttons = []
    for i, pore in enumerate(trace_keys):
        visible = [j == i for j in range(len(trace_keys))]
        buttons.append(dict(
            label=pore,
            method="update",
            args=[
                {"visible": visible},
                {"annotations": [
                    dict(
                        text="Total Individuals per Population: " + subtitle_map[pore],
                        x=-0.06, y=1.15, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=12, color="black"),
                        align="left", xanchor="left"
                    ),
                    dict(
                        text="Pore Type:",
                        x=1.15, y=1.23, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=12, color="black"),
                        xanchor="center"
                    )
                ]}
            ]
        ))

    fig.update_layout(
        width=FIG_WIDTH, height=FIG_HEIGHT, margin=dict(t=125),
        bargap=0.5, plot_bgcolor="white", showlegend=False,
        font=dict(color="black"),
        title=dict(
            y=0.95,
            text=(
                "<span style='font-size:18px; font-weight:bold'>Pathogenic Genotype Distribution</span><br>"
                f"<span style='font-size:12px'>{gene} - {disease}</span><br>"
                f"<span style='font-size:12px'>Inheritance: {inheritance}</span>"
            )
        ),
        updatemenus=[dict(
            active=default_idx, buttons=buttons, direction="down",
            showactive=True, x=1.15, xanchor="center", y=1.15, yanchor="top"
        )],
        annotations=[
            dict(
                text="Total Individuals per Population: " + subtitle_map[default_key],
                x=-0.0555, y=1.15, xref="paper", yref="paper",
                showarrow=False, font=dict(size=12, color="black"),
                align="left", xanchor="left"
            ),
            dict(
                text="Pore Type:",
                x=1.15, y=1.23, xref="paper", yref="paper",
                showarrow=False, font=dict(size=12, color="black"),
                xanchor="center"
            )
        ],
        xaxis=dict(
            range=[0, 100], tickmode="linear", tick0=0, dtick=10,
            ticks="outside", showline=True, linecolor="black", linewidth=1,
            zeroline=False, tickfont=dict(color="black"),
            title=dict(text="Pathogenic Genotypes (%)", font=dict(color="black"))
        ),
        yaxis=dict(
            title="",
            showline=True,
            linecolor="black",
            categoryorder="array",
            categoryarray=superpop_order
        )
    )
    return fig


# ----------------------------
# Main
# ----------------------------
def main():
    """
    Entry point: loads spreadsheet, normalizes fields, computes per-population
    summaries, and writes HTML/PNG outputs. In test mode, previews a limited
    number of plots and writes to 'test_outputs'.
    """
    args = parse_args()

    test_mode = bool(args.test)
    test_limit = int(args.test_limit)
    save_test_outputs = bool(args.save_test_outputs)

    allele_path = args.allele_spreadsheet
    if not os.path.exists(allele_path):
        raise SystemExit(f"Allele spreadsheet not found: {allele_path}")

    out_base = Path(args.output_dir) if args.output_dir else OUTPUT_BASE
    if test_mode:
        output_dir = out_base / "test_outputs"
        output_html_dir = output_dir
        output_png_dir = output_dir
    else:
        output_dir = out_base
        output_html_dir = output_dir / "html"
        output_png_dir = output_dir / "png"

    output_html_dir.mkdir(parents=True, exist_ok=True)
    output_png_dir.mkdir(parents=True, exist_ok=True)

    # --- Load allele spreadsheet ---
    df = pd.read_excel(allele_path, sheet_name=SHEET_NAME, engine="openpyxl")
    col = build_col_lookup(df)

    # More robust population & pore detection
    sample_id_col = col("Sample ID", "Sample_ID", "sample_id")
    gene_col = col("Gene")
    disease_col = col("Disease")
    repeat_col = col("Repeat Count", "RepeatCount", "repeat_count", "Repeat count")
    pathogenic_min_col = col("Pathogenic min", "Pathogenic_min", "pathogenic_min")
    pathogenic_max_col = col("Pathogenic max", "Pathogenic_max", "pathogenic_max")
    inheritance_col = col("Inheritance", "inheritance")
    sex_col = col("Sex", "sex", "Gender", "gender")
    population_col = col(
        "Cleaned population description", "Population", "Ancestry",
        "SuperPop", "Super Population", "Superpopulation", "Super-population",
        "SubPop", "Sub Population", "Subpopulation"
    )
    pore_col = col("Pore", "Pore Type", "PoreType", "Pore type")

    required = [sample_id_col, gene_col, disease_col, repeat_col]
    if not all(required):
        raise SystemExit(f"Spreadsheet missing required columns. Found: {list(df.columns)}")

    df = df.rename(columns={
        sample_id_col: "Sample ID",
        gene_col: "Gene",
        disease_col: "Disease",
        repeat_col: "Repeat count",
        pathogenic_min_col: "Pathogenic min",
        pathogenic_max_col: "Pathogenic max",
        inheritance_col: "Inheritance",
        sex_col: "Sex",
        population_col: "Population",
        pore_col: "Pore",
    })

    # If some optional cols missing, create them
    if "Population" not in df.columns:
        df["Population"] = "Unknown"
    if "Pore" not in df.columns:
        df["Pore"] = "All Types"
    if "Sex" not in df.columns:
        df["Sex"] = "Unknown"
    if "Inheritance" not in df.columns:
        df["Inheritance"] = np.nan

    # Clean fields
    df["Inheritance"] = df["Inheritance"].apply(clean_inheritance).astype(str).str.upper().replace({"NAN": np.nan})
    df["Sex"] = df["Sex"].apply(normalize_sex)
    df["Pore"] = df["Pore"].apply(normalize_pore)
    df["Population"] = df["Population"].apply(normalize_population)

    # --- Flag pathogenic alleles ---
    def is_pathogenic_row(row):
        mn = row.get("Pathogenic min")
        mx = row.get("Pathogenic max")
        rc = row.get("Repeat count")
        if pd.notna(mn) and pd.notna(mx) and pd.notna(rc):
            try:
                return float(mn) <= float(rc) <= float(mx)
            except Exception:
                return np.nan
        return np.nan

    df["Is pathogenic"] = df.apply(is_pathogenic_row, axis=1)

    # Build superpop_order (always include canonical populations in standard order, Unknown at bottom)
    all_pops = sorted(set(df["Population"].dropna().astype(str).tolist()))
    base_no_unknown = [p for p in SUPERPOP_ORDER_BASE if p != "Unknown"]
    extras = [p for p in all_pops if p not in base_no_unknown and p not in {"All", "Unknown"}]
    superpop_order = base_no_unknown + extras + ["Unknown"]

    # --- Generate plots ---
    made = 0
    for gene in df["Gene"].dropna().astype(str).unique():
        df_gene = df[df["Gene"].astype(str) == str(gene)]
        for disease in df_gene["Disease"].dropna().astype(str).unique():
            df_gd = df_gene[df_gene["Disease"].astype(str) == str(disease)]
            if df_gd["Sample ID"].nunique() == 0:
                continue

            fig = create_horizontal_bar_plot(df, gene, disease, superpop_order)
            if fig is None:
                continue

            safe_gene = re.sub(r"[\\/]", "_", str(gene))
            safe_disease = re.sub(r"[\\/]", "_", str(disease))

            html_path = output_html_dir / f"{safe_gene}_{safe_disease}_ancestry_plot.html"
            png_path = output_png_dir / f"{safe_gene}_{safe_disease}_ancestry_plot.png"

            fig.write_html(str(html_path))
            try:
                fig.write_image(str(png_path), scale=2)
            except Exception:
                pass

            if test_mode:
                print(f"Previewing: {gene} / {disease}")
                fig.show(renderer="browser")
                if save_test_outputs:
                    print(f"Saved: {html_path}")
                made += 1
                if made >= test_limit:
                    break

        if test_mode and made >= test_limit:
            break

    print("\n--- Done ---" + (" (test mode)" if test_mode else ""))


if __name__ == "__main__":
    main()
