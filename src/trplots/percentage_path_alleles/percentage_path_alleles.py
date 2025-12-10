"""
Percentage of Pathogenic Alleles per Gene + Disease

Purpose:
    Load an allele spreadsheet (Excel/CSV), determine pathogenic alleles,
    group by Gene and Disease, and report counts and percent pathogenic.

Key behaviors / features:
    - Reads .xlsx or .csv input (via --input/-i).
    - Pathogenic flagging prefers the "Is Pathogenic" column; otherwise uses
      Repeat count within Pathogenic min/max (inclusive).
    - Groups by Gene, Disease (keeps NaNs via dropna=False) and computes:
        total_alleles, pathogenic_alleles, percent_pathogenic.
    - Writes XLSX output to --output/-o or alongside the input.
    - Creates parent output directories if missing.
"""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

# ---------- Defaults (override via CLI) ----------
# Base everything off the repo root so paths stay portable when this file moves.
BASE_DIR = Path(__file__).resolve().parents[3]  # .../tr-plots
DATA_DIR = BASE_DIR / "data"
# Default input/output paths (can be overridden via CLI flags)
DEFAULT_INPUT = DATA_DIR / "other_data/allele_spreadsheet.xlsx"
DEFAULT_OUTPUT = DATA_DIR / "other_data/percentage_pathogenic_by_gene_disease.xlsx"

# ---------- CLI ----------
def parse_args():
    # Simple CLI: input spreadsheet and optional output path
    p = argparse.ArgumentParser(description="Compute percentage of pathogenic alleles per Gene+Disease")
    p.add_argument("--input", "-i", type=str, default=DEFAULT_INPUT,
                   help="Input spreadsheet (.xlsx or .csv)")
    p.add_argument("--output", "-o", type=str, default=None,
                   help="Output Excel path (defaults to data/other_data/percentage_pathogenic_by_gene_disease.xlsx)")
    return p.parse_args()

def is_pathogenic_series(df: pd.DataFrame) -> pd.Series:
    """
    Return a boolean Series marking rows as pathogenic.
    Priority:
      1) Use explicit 'Is Pathogenic' if present.
      2) Else, derive from Repeat count vs Pathogenic min/max (inclusive).
      3) Else, all False.
    """
    # Prefer explicit 'Is Pathogenic' column if present
    if "Is Pathogenic" in df.columns:
        s = df["Is Pathogenic"].fillna("").astype(str).str.strip().str.upper()
        return s.isin({"YES", "Y", "1", "TRUE", "T"})
    # Fallback: use Pathogenic min/max and Repeat count if available
    if {"Pathogenic min", "Pathogenic max", "Repeat count"}.issubset(df.columns):
        rc = pd.to_numeric(df["Repeat count"], errors="coerce")
        pmin = pd.to_numeric(df["Pathogenic min"], errors="coerce")
        pmax = pd.to_numeric(df["Pathogenic max"], errors="coerce")
        # True when repeat count is present and within available bounds
        return (~rc.isna()) & (
            (((~pmin.isna()) & (rc >= pmin)) | pmin.isna()) &
            (((~pmax.isna()) & (rc <= pmax)) | pmax.isna())
        )
    # No way to determine: all non-pathogenic
    return pd.Series([False] * len(df), index=df.index)

def main():
    args = parse_args()
    inp = Path(args.input)
    if not inp.exists():
        raise SystemExit(f"Input not found: {inp}")

    # Load input (Excel or CSV)
    if inp.suffix.lower() in (".xls", ".xlsx"):
        df = pd.read_excel(inp, engine="openpyxl")
    else:
        df = pd.read_csv(inp)

    # Basic sanity check
    if "Gene" not in df.columns or "Disease" not in df.columns:
        raise SystemExit("Input must contain 'Gene' and 'Disease' columns")

    # Compute pathogenic boolean per row
    pathogenic_mask = is_pathogenic_series(df)
    df["_is_pathogenic_bool"] = pathogenic_mask

    # Group by Gene/Disease (keep NaNs) and count totals/pathogenic
    grouped = df.groupby(["Gene", "Disease"], dropna=False).agg(
        total_alleles=("Gene", "count"),
        pathogenic_alleles=("_is_pathogenic_bool", "sum")
    ).reset_index()

    # Percent pathogenic (safe divide; fill NaNs with 0.0)
    grouped["percent_pathogenic"] = (
        np.round((grouped["pathogenic_alleles"] / grouped["total_alleles"]) * 100, 2)
        .fillna(0.0)
    )

    # Resolve output path (force .xlsx, default to data/other_data/percentage_pathogenic_by_gene_disease.xlsx)
    outp = Path(args.output) if args.output else DEFAULT_OUTPUT.with_suffix(".xlsx")
    # If a non-xlsx extension was provided, replace it with .xlsx
    if outp.suffix.lower() != ".xlsx":
        outp = outp.with_suffix(".xlsx")

    outp.parent.mkdir(parents=True, exist_ok=True)

    # Write Excel only
    grouped.to_excel(outp, index=False, engine="openpyxl")

    if outp:
        print(f"Created spreadsheet: {outp}")
    else:
        print("No spreadsheet created.")

if __name__ == "__main__":
    main()