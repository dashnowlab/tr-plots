"""
---------------------------------------------
 Script: VCF → Locus–Sample Table (CSV/Excel)
 Purpose:
   - Reads a VCF of tandem-repeat calls and a JSON of locus metadata
   - Joins population info (SubPop/SuperPop) from a CSV
   - Computes repeat counts per allele using motif length
   - Flags alleles within pathogenic thresholds when available
   - Expands per-sample arrays → one row per allele
   - Saves CSV and Excel
   - Test mode: limit variants and optionally save to test_outputs/
---------------------------------------------
"""

import json
import os
import pysam
import pandas as pd
import argparse
from pathlib import Path

# --- TEST MODE ---
TEST_MODE = True                 # Quick testing: process only some VCF records
TEST_LIMIT = 100                 # Max VCF variants to process in test mode
SAVE_TEST_OUTPUTS = False         # Save files in test mode if True

# --- File locations (use central config) ---
from trplots.config import SEQ_DATA, OTHER_DATA, OUTPUT_BASE

# Paths
VCF_PATH = SEQ_DATA / "83_loci_503_samples" / "1000g-ont-strchive-83_loci_503_samples.vcf.gz"
JSON_PATH = OTHER_DATA / "strchive-loci.json"
# the summary file in repo is named '1kgp_ont_500_summary_-_sheet1.csv'
CSV_PATH = OTHER_DATA / "1kgp_ont_500_summary_-_sheet1.csv"

# Output roots (normal mode)
OUTPUT_BASE = OUTPUT_BASE / "matching_files_outputs"
OUTPUT_DIR_CSV = OUTPUT_BASE / "csvs"
OUTPUT_DIR_XLSX = OUTPUT_BASE / "excels"
os.makedirs(OUTPUT_DIR_CSV, exist_ok=True)
os.makedirs(OUTPUT_DIR_XLSX, exist_ok=True)

# If test mode, write into a single test_outputs folder
if TEST_MODE:
    OUTPUT_DIR = OUTPUT_BASE / "test_outputs"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    CSV_OUT = OUTPUT_DIR / "vcf_locus_sample_table.csv"
    XLSX_OUT = OUTPUT_DIR / "vcf_locus_sample_table.xlsx"
else:
    CSV_OUT = OUTPUT_DIR_CSV / "vcf_locus_sample_table.csv"
    XLSX_OUT = OUTPUT_DIR_XLSX / "vcf_locus_sample_table.xlsx"

# --- Load population info (CSV) ---
sample_info = pd.read_csv(CSV_PATH)

# Keep Sample_ID, SubPop, SuperPop (first, third, fourth columns in your file)
sample_info = sample_info.iloc[:, [0, 2, 3]]
sample_info.columns = ['Sample_ID', 'SubPop', 'SuperPop']

# Map Sample_ID → {SubPop, SuperPop}
population_map = (
    sample_info
    .drop_duplicates(subset="Sample_ID")
    .set_index("Sample_ID")[["SubPop", "SuperPop"]]
    .to_dict(orient="index")
)

# --- Load locus metadata (JSON) ---
with open(JSON_PATH, "r") as f:
    loci_data = json.load(f)

# --- Read VCF header + iterate records ---
def parse_args():
    p = argparse.ArgumentParser(description="Process VCF and produce allele/ancestry CSV")
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

    vcf_in = pysam.VariantFile(VCF_PATH)
    vcf_sample_ids = list(vcf_in.header.samples)                  # e.g., ["HG00096-1", ...]
    base_sample_ids = [sid.split('-')[0] for sid in vcf_sample_ids]  # normalized IDs

    # Collections / counters
    records = []
    total_vcf_loci = 0
    matched_vcf_loci = 0

    # --- Iterate through VCF ---
    for i, record in enumerate(vcf_in.fetch()):
        if TEST_MODE and i >= TEST_LIMIT:
            break

        total_vcf_loci += 1
        chrom = record.chrom
        pos   = record.pos

        # Match current variant to a locus entry in JSON
        matched_locus = next(
            (entry for entry in loci_data
             if chrom == entry["chrom"] and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
            None
        )
        if not matched_locus:
            print(f"VCF locus {chrom}:{pos} not matched in loci_data")
            continue
        matched_vcf_loci += 1

        # Motif and length
        motif = matched_locus.get("reference_motif_reference_orientation")
        motif_str = "".join(motif) if isinstance(motif, list) else str(motif)
        motif_len = len(motif_str)

        # Thresholds (may be missing)
        pathogenic_min = matched_locus.get("pathogenic_min", None)
        pathogenic_max = matched_locus.get("pathogenic_max", None)

        # --- Per-sample fields (safe: drop None alleles) ---
        for sid, sample in record.samples.items():
            base_id = sid.split('-')[0]
            al_lengths = sample.get("AL")  # allele lengths in bases
            if not al_lengths:
                continue

            # Keep only non-null allele lengths and cast to int
            clean_lengths = [int(x) for x in al_lengths if x is not None]

            if not clean_lengths:
                continue

            # Repeat counts per allele (integer division by motif_len)
            repeat_counts = [length // motif_len for length in clean_lengths]

            # Population info
            pop_info = population_map.get(base_id, {"SubPop": "Unknown", "SuperPop": "Unknown"})

            # Pathogenic flag per allele, only if thresholds provided
            if (pathogenic_min is not None) and (pathogenic_max is not None):
                is_pathogenic = [pathogenic_min <= rc <= pathogenic_max for rc in repeat_counts]
            else:
                is_pathogenic = [False] * len(repeat_counts)

            # Store compact; NOTE: we only store cleaned values (no "None" strings)
            records.append({
                "locus_chrom": chrom,
                "locus_pos": pos,
                "sample_id": sid,
                "gene": matched_locus.get("gene", "Unknown"),
                "disease": matched_locus.get("disease", "Unknown"),
                "population": pop_info["SubPop"],
                "population_description": pop_info["SuperPop"],
                "motif": motif_str,
                "motif_length": motif_len,
                "allele_lengths": ";".join(map(str, clean_lengths)),
                "repeat_counts": ";".join(map(str, repeat_counts)),
                "is_pathogenic": ";".join("True" if b else "False" for b in is_pathogenic),
                "pathogenic_min": pathogenic_min,
                "pathogenic_max": pathogenic_max
            })

    # --- Expand compact entries to one row per allele (robust to stray values) ---
    expanded_records = []
    skipped = 0

    for rec in records:
        allele_lengths = rec["allele_lengths"].split(";") if rec["allele_lengths"] else []
        repeat_counts  = rec["repeat_counts"].split(";")  if rec["repeat_counts"]  else []
        is_pathogenic  = rec["is_pathogenic"].split(";")  if rec["is_pathogenic"]  else []

        # Arrays must align
        if not (len(allele_lengths) == len(repeat_counts) == len(is_pathogenic)):
            skipped += 1
            continue

        for a, rc, ip in zip(allele_lengths, repeat_counts, is_pathogenic):
            # Extra safety: skip any residual non-numeric entries
            if a == "" or rc == "":
                continue
            try:
                a_int = int(a)
                rc_int = int(rc)
            except ValueError:
                skipped += 1
                continue

            expanded_records.append({
                "gene": rec["gene"],
                "locus_chrom": rec["locus_chrom"],
                "locus_pos": rec["locus_pos"],
                "sample_id": rec["sample_id"],
                "disease": rec["disease"],
                "population": rec["population"],
                "population_description": rec["population_description"],
                "motif": rec["motif"],
                "motif_length": rec["motif_length"],
                "allele_length": a_int,
                "repeat_count": rc_int,
                "pathogenic_min": rec["pathogenic_min"],
                "pathogenic_max": rec["pathogenic_max"],
                "is_pathogenic": (ip == "True")
            })

    # --- Build DataFrame ---
    df = pd.DataFrame(expanded_records)

    # --- Save outputs ---
    if TEST_MODE:
        print(f"Preview rows: {len(df)}  |  Skipped misaligned: {skipped}")
        if SAVE_TEST_OUTPUTS:
            df.to_csv(CSV_OUT, index=False)
            df.to_excel(XLSX_OUT, index=False)
            print(f"[TEST] Saved CSV  → {CSV_OUT}")
            print(f"[TEST] Saved XLSX → {XLSX_OUT}")
    else:
        df.to_csv(CSV_OUT, index=False)
        df.to_excel(XLSX_OUT, index=False)
        print(f"Saved CSV  → {CSV_OUT}")
        print(f"Saved XLSX → {XLSX_OUT}")

    # --- Summary ---
    if {'locus_chrom', 'locus_pos'}.issubset(df.columns):
        n_unique = df[['locus_chrom', 'locus_pos']].drop_duplicates().shape[0]
        print(f"Unique loci in output: {n_unique}")
    else:
        print("Columns 'locus_chrom' and/or 'locus_pos' missing in output DataFrame.")

    print(f"Total loci visited in VCF: {total_vcf_loci}")
    print(f"Loci matched to metadata:  {matched_vcf_loci}")
    print(f"Misaligned compact rows skipped during expansion: {skipped}")
    print("--- Done ---")


if __name__ == "__main__":
    args = parse_args()
    main(args)
