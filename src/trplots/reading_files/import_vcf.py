"""
---------------------------------------------
 Script: VCF → Long-Format Table
 Purpose:
   - Reads an indexed VCF
   - Flattens VARIANT x SAMPLE data into long format
   - Optionally includes INFO and FORMAT fields
   - Saves the result as a CSV
---------------------------------------------
"""

import os
import pysam
import pandas as pd

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots"

# Dataset
VCF_PATH = f"{BASE_DIR}/Data/Sequencing Data/83 Loci 503 Samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# Output
OUTPUT_DIR = os.path.join(BASE_DIR, "Results/Reading Files Outputs/import_vcf_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "vcf_long_format.csv")

# --- Open VCF (requires .tbi/.csi index) ---
vcf_in = pysam.VariantFile(VCF_PATH)
samples = list(vcf_in.header.samples)

# --- Collect long-format rows ---
records = []
for record in vcf_in.fetch():
    # Core variant fields
    variant_info = {
        "CHROM": record.chrom,                                    # chromosome
        "POS": record.pos,                                        # 1-based position
        "ID": record.id,                                          # variant ID (can be None)
        "REF": record.ref,                                        # reference allele
        "ALT": ",".join(map(str, record.alts or [])),             # alt alleles joined
        "QUAL": record.qual,                                      # quality (float or None)
        "FILTER": ";".join(record.filter.keys()) if record.filter else "",  # filters
    }

    # Include INFO fields
    for key, value in record.info.items():
        variant_info[key] = value

    # One row per sample with FORMAT fields expanded
    for sample in samples:
        row = variant_info.copy()
        row["SAMPLE"] = sample

        sample_data = record.samples[sample]
        for fmt_key, fmt_value in sample_data.items():
            row[fmt_key] = fmt_value

        records.append(row)

# --- To DataFrame ---
df_long = pd.DataFrame(records)

# Optional: multi-index for fast slicing by locus + sample
# df_long.set_index(["CHROM", "POS", "SAMPLE"], inplace=True)

# --- Save ---
df_long.to_csv(OUTPUT_CSV, index=False)
print(f"--- Done ---")