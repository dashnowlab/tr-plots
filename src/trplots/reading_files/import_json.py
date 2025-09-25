"""
---------------------------------------------
 Script: Inspect STRchive Loci JSON
 Purpose:
   - Loads loci metadata from STRchive JSON file
   - Extracts chromosome, start position, and 
     benign/pathogenic threshold ranges
   - Prints formatted output for review
---------------------------------------------
"""

import json
from trplots.config import OTHER_DATA

# --- File location ---
JSON_PATH = OTHER_DATA / "strchive-loci.json"

# --- Load JSON metadata ---
with open(JSON_PATH, "r") as file:
    loci_data = json.load(file)

# --- Helper for formatting values ---
def fmt(val):
    """Return int if available, else 'null'."""
    return int(val) if val is not None else "null"

# --- Extract and print relevant fields ---
for entry in loci_data:
    chrom = entry.get("chrom")
    start = entry.get("start_hg38")
    benign_min = entry.get("benign_min")
    benign_max = entry.get("benign_max")
    pathogenic_min = entry.get("pathogenic_min")
    pathogenic_max = entry.get("pathogenic_max")

    print(
        f"Chromosome: {chrom}, Start: {start}, "
        f"Benign: {fmt(benign_min)}-{fmt(benign_max)}, "
        f"Pathogenic: {fmt(pathogenic_min)}-{fmt(pathogenic_max)}"
    )

print("--- Done ---")
