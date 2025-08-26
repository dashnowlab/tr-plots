"""
---------------------------------------------
 Script: Check Specific Loci in VCF
 Purpose:
   - Defines a set of user-specified loci (chromosome, position)
   - Opens a target VCF file
   - Checks which of the specified loci are present in the file
   - Reports whether each locus is FOUND or MISSING
---------------------------------------------
"""

import pysam
import os

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/Anschutz"

# Choose the dataset of interest
VCF_PATH = f"{BASE_DIR}/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# --- User-defined loci of interest (chromosome, position) ---
USER_LOCI = {
    ("chr3", 63912684),
    ("chr13", 70139383),
    ("chr13", 70139429),
    ("chrX", 147912049),
    ("chr9", 69037286),
    ("chr9", 69037304),
    ("chr5", 10356343),
    ("chr1", 155188505),
    ("chr1", 156591765),
    ("chr15", 88569433),
    ("chr20", 4699397),
    ("chr20", 4699493),
    ("chr17", 17808358),
    ("chrX", 71453054),
}

# --- Open the VCF file ---
vcf = pysam.VariantFile(VCF_PATH)

# --- Check loci presence ---
found = set()
for record in vcf.fetch():
    if (record.chrom, record.pos) in USER_LOCI:
        found.add((record.chrom, record.pos))

# --- Report results ---
print("=== Loci Presence Report ===")
for chrom, pos in USER_LOCI:
    status = "FOUND" if (chrom, pos) in found else "MISSING"
    print(f"{chrom}:{pos}\t{status}")

print("--- Done ---")
