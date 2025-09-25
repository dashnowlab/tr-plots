"""
---------------------------------------------
 Script: Extract Sample_ID + Pore
 Purpose:
   - Reads the 1KGP/ONT sample summary CSV
   - Keeps only Sample_ID and Pore columns
   - Saves a clean CSV for downstream merges
   - (VCF path is validated/opened to ensure consistency of dataset paths)
---------------------------------------------
"""


import os
import pandas as pd
import pysam
from pathlib import Path

from trplots.config import SEQ_DATA, OTHER_DATA, OUTPUT_BASE

VCF_PATH = SEQ_DATA / "83_loci_503_samples" / "1000g-ont-strchive-83_loci_503_samples.vcf.gz"
INPUT_CSV = OTHER_DATA / "1kgp_ont_500_summary_-_sheet1.csv"
OUTPUT_DIR = Path(OUTPUT_BASE) / "reading_files_outputs" / "check_sample_ids_output"
OUTPUT_CSV = OUTPUT_DIR / "1KGP_ONT_500_Summary_Sample_ID_Pore.csv"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Validate/open VCF (not used further; just to confirm path is correct) ---
_ = pysam.VariantFile(VCF_PATH)

# --- Load the CSV (header row is the 2nd row in the file) ---
df = pd.read_csv(INPUT_CSV, header=1)

# --- Keep only needed columns ---
required_cols = ["Sample_ID", "Pore"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise KeyError(f"Missing expected columns in input CSV: {missing}")

df = df[required_cols]

# --- Save cleaned CSV ---
df.to_csv(OUTPUT_CSV, index=False)
print(f"--- Done ---")
