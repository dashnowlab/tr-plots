#!/usr/bin/env python3
"""
---------------------------------------------
 Script: Motif Metrics Calculator
 Purpose:
   - Compute motif-aware metrics from STR VCF.
   - Outputs allele-level and per-sample motif repeat statistics.
---------------------------------------------
"""

import argparse
import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional

# ====== Base paths ======
BASE_DIR = Path("/Users/annelisethorn/Documents/github/tr-plots")

# Input VCF (default)
VCF_FILE = (
    BASE_DIR
    / "data"
    / "sequencing_data"
    / "81_loci_502_samples"
    / "1000g-ont-strchive-81_loci_502_samples_81224_alleles.vcf.gz"
)

# Output root directory
OUT_ROOT = BASE_DIR / "results" / "motif_counts"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

# ====== Optional: pandas for Excel output ======
try:
    import pandas as pd
except Exception:
    pd = None


# ---------- Utility functions ----------
def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_info_field(info_str: str) -> Dict[str, str]:
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
    return info


def full_pure_count(seq: str, motif: str) -> int:
    """Return n if seq == motif * n; else 0."""
    if not seq or not motif:
        return 0
    m = len(motif)
    if m == 0 or len(seq) % m != 0:
        return 0
    n = len(seq) // m
    return n if motif * n == seq else 0


def max_consecutive_motif_run(seq: str, motif: str) -> int:
    """Return the maximum uninterrupted motif repeat run within the sequence."""
    if not seq or not motif:
        return 0
    m, best = len(motif), 0
    i = 0
    while i <= len(seq) - m:
        run = 0
        while i + (run + 1) * m <= len(seq) and seq[i + run * m : i + (run + 1) * m] == motif:
            run += 1
        if run > 0:
            best = max(best, run)
            i += run * m
        else:
            i += 1
    return best


def safe_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None


def write_table(out_prefix: Path, name: str, rows: List[Dict[str, object]]):
    csv_path = out_prefix.with_suffix(f".{name}.csv")
    if not rows:
        return
    headers = list(rows[0].keys())
    with open(csv_path, "w") as f:
        f.write(",".join(headers) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(h, "")) for h in headers) + "\n")

    if pd is not None:
        try:
            pd.DataFrame(rows).to_excel(out_prefix.with_suffix(f".{name}.xlsx"), index=False)
        except Exception:
            pass


# ---------- Main processing ----------
def process_vcf(vcf_path: Path, out_prefix: Path, per_sample: bool = False, pathogenic_threshold: Optional[float] = None):
    allele_rows = []
    sample_rows = []

    with open_maybe_gzip(vcf_path) as f:
        samples = []
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue

            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue
            CHROM, POS, ID, REF, ALT, QUAL, FILT, INFO = fields[:8]
            format_str = fields[8] if len(fields) > 8 else ""
            sample_data = fields[9:] if len(fields) > 9 else []

            info = parse_info_field(INFO)
            motif = info.get("MOTIF", "")
            motif_len = len(motif)
            AN = info.get("AN", "")
            ACs = info.get("AC", "").split(",") if "AC" in info else []
            alts = ALT.split(",") if ALT else []

            # --- REF allele ---
            allele_rows.append({
                "CHROM": CHROM,
                "POS": POS,
                "AlleleType": "REF",
                "AlleleSeq": REF,
                "Motif": motif,
                "MotifLen": motif_len,
                "AlleleLen": len(REF),
                "ApproxRepeatCount": len(REF)/motif_len if motif_len else None,
                "PureCount": full_pure_count(REF, motif),
                "MaxMotifRun": max_consecutive_motif_run(REF, motif),
                "IsPure": full_pure_count(REF, motif) > 0,
                "AC": "",
                "AN": AN
            })

            # --- ALT alleles ---
            for i, alt in enumerate(alts):
                allele_rows.append({
                    "CHROM": CHROM,
                    "POS": POS,
                    "AlleleType": f"ALT{i+1}",
                    "AlleleSeq": alt,
                    "Motif": motif,
                    "MotifLen": motif_len,
                    "AlleleLen": len(alt),
                    "ApproxRepeatCount": len(alt)/motif_len if motif_len else None,
                    "PureCount": full_pure_count(alt, motif),
                    "MaxMotifRun": max_consecutive_motif_run(alt, motif),
                    "IsPure": full_pure_count(alt, motif) > 0,
                    "AC": ACs[i] if i < len(ACs) else "",
                    "AN": AN
                })

    # --- Write output ---
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_table(out_prefix, "alleles", allele_rows)
    print(f"✅ Saved allele-level motif data to: {out_prefix}.alleles.csv/.xlsx")


# ---------- CLI ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute motif-aware metrics from STR VCF.")
    parser.add_argument("--vcf", default=str(VCF_FILE), help="Path to VCF (.vcf or .vcf.gz)")
    parser.add_argument("--out", default=str(OUT_ROOT / "motif_metrics"), help="Output file prefix")
    args = parser.parse_args()

    process_vcf(Path(args.vcf), Path(args.out))
