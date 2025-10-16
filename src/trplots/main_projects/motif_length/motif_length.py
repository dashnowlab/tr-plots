"""
---------------------------------------------
 Script: Longest Motif Lengths
 Purpose:
   - Extract the longest motif lengths from a VCF file.
---------------------------------------------
"""

import gzip
import re
from collections import defaultdict
from pathlib import Path

from readchar import config

# import shared paths/helpers from config
from trplots.config import (
    SEQ_DATA,                
    ENSURE_DIR
)

vcf_path = SEQ_DATA / "83_loci_503_samples" / "1000g-ont-strchive-83_loci_503_samples.vcf.gz"

longest_motifs = defaultdict(lambda: {"motif": "", "length": 0})

with gzip.open(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        chrom = parts[0]
        pos = parts[1]
        info_field = parts[7]

        info_dict = dict(item.split("=", 1) for item in info_field.split(";") if "=" in item)
        
        # Prefer 'MOTIF' if present; fallback to 'RU'
        motif = info_dict.get("MOTIF") or info_dict.get("RU")
        gene = info_dict.get("GENE") or f"{chrom}:{pos}"

        if motif:
            motif_len = len(motif)
            if motif_len > longest_motifs[gene]["length"]:
                longest_motifs[gene] = {"motif": motif, "length": motif_len}

# Display results
# for gene, data in longest_motifs.items():
#     print(f"{gene}\tMotif: {data['motif']}\tLength: {data['length']}")

# Save to CSV
out_file = config.RESULT_PATH("longest_motifs.csv")
with out_file.open("w") as out:
    out.write("Gene_or_Locus,Motif,Length\n")
    for gene, data in longest_motifs.items():
        out.write(f"{gene},{data['motif']},{data['length']}\n")
