"""
---------------------------------------------
 Script: VCF (gz) → CSV + Excel
 Purpose:
   - Reads a compressed VCF (.vcf.gz)
   - Writes a flat CSV of the VCF table (header + rows)
   - Exports the CSV to Excel (.xlsx)
---------------------------------------------
"""

import os
import gzip
import csv
import pandas as pd
from pathlib import Path

from trplots.config import SEQ_DATA, OUTPUT_BASE

# --- File locations ---
VCF_PATH = SEQ_DATA / "83_loci_503_samples" / "1000g-ont-strchive-83_loci_503_samples.vcf.gz"

CSV_DIR = Path(OUTPUT_BASE) / "reading_files_outputs" / "make_VCF_CSV_Excel_output"
EXCEL_DIR = CSV_DIR
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(EXCEL_DIR, exist_ok=True)

# Output filenames
CSV_PATH = CSV_DIR / "83_loci_503_samples_CSV.csv"
EXCEL_PATH = EXCEL_DIR / "83_loci_503_samples_Excel.xlsx"

# --- Convert VCF (gz) to CSV ---
with gzip.open(VCF_PATH, "rt") as vcf_in, open(CSV_PATH, "w", newline="") as csv_out:
    writer = csv.writer(csv_out)

    for line in vcf_in:
        line = line.strip()
        if line.startswith("##"):
            continue  # skip meta-information lines
        elif line.startswith("#CHROM"):
            header = line.lstrip("#").split("\t")
            writer.writerow(header)
        else:
            fields = line.split("\t")
            writer.writerow(fields)

# --- Save CSV to Excel ---
df = pd.read_csv(CSV_PATH)
df.to_excel(EXCEL_PATH, index=False)

print(f"--- Saved CSV ---")
print(f"--- Saved XLSX ---")
