import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def parse_args():
    p = argparse.ArgumentParser(description="Compute percentage of pathogenic alleles per Gene+Disease")
    p.add_argument("--input", "-i",
                   default="/Users/annelisethorn/Documents/github/tr-plots/data/other_data/all_path_ref_motif_allele_spreadsheet.xlsx",
                   help="Input spreadsheet (.xlsx or .csv)")
    p.add_argument("--output", "-o",
                   default=None,
                   help="Output CSV/Excel path (defaults to input folder /percentage_pathogenic_by_gene_disease.csv)")
    return p.parse_args()

def is_pathogenic_series(df: pd.DataFrame) -> pd.Series:
    # Prefer explicit 'Is Pathogenic' column if present
    if "Is Pathogenic" in df.columns:
        s = df["Is Pathogenic"].fillna("").astype(str).str.strip().str.upper()
        return s.isin({"YES", "Y", "1", "TRUE", "T"})
    # Fallback: use Pathogenic min/max and Repeat count if available
    if {"Pathogenic min", "Pathogenic max", "Repeat count"}.issubset(df.columns):
        rc = pd.to_numeric(df["Repeat count"], errors="coerce")
        pmin = pd.to_numeric(df["Pathogenic min"], errors="coerce")
        pmax = pd.to_numeric(df["Pathogenic max"], errors="coerce")
        # For each row: True if rc is not NaN and within available bounds (inclusive)
        return (~rc.isna()) & (
            ((~pmin.isna()) & (rc >= pmin) | pmin.isna()) &
            ((~pmax.isna()) & (rc <= pmax) | pmax.isna())
        )
    # No way to determine: return all False
    return pd.Series([False] * len(df), index=df.index)

def main():
    args = parse_args()
    inp = Path(args.input)
    if not inp.exists():
        raise SystemExit(f"Input not found: {inp}")

    if inp.suffix.lower() in (".xls", ".xlsx"):
        df = pd.read_excel(inp, engine="openpyxl")
    else:
        df = pd.read_csv(inp)

    # ensure Gene and Disease columns exist
    if "Gene" not in df.columns or "Disease" not in df.columns:
        raise SystemExit("Input must contain 'Gene' and 'Disease' columns")

    pathogenic_mask = is_pathogenic_series(df)
    df["_is_pathogenic_bool"] = pathogenic_mask

    # Group and compute totals; remove earlier placeholder that referenced a non-existent column.
    grouped = df.groupby(["Gene", "Disease"], dropna=False).agg(
        total_alleles=("Gene", "count"),
        pathogenic_alleles=("_is_pathogenic_bool", "sum")
    ).reset_index()

    # percent_pathogenic: handle any potential zero-division and produce 0.0 for empty groups
    grouped["percent_pathogenic"] = (
        np.round((grouped["pathogenic_alleles"] / grouped["total_alleles"].replace({0: np.nan})) * 100, 2)
        .fillna(0.0)
    )

    # output path
    if args.output:
        outp = Path(args.output)
    else:
        outp = inp.parent / "percentage_pathogenic_by_gene_disease.csv"

    # ensure directory exists
    outp.parent.mkdir(parents=True, exist_ok=True)

    # save both CSV and XLSX for convenience
    grouped.to_csv(outp, index=False)
    try:
        xls_path = outp.with_suffix(".xlsx")
        grouped.to_excel(xls_path, index=False)
    except Exception:
        pass

    print(f"Wrote {len(grouped)} rows to {outp} (and .xlsx where possible)")

if __name__ == "__main__":
    main()