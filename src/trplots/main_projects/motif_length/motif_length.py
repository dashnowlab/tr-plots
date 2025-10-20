"""
---------------------------------------------
 Script: Longest Motif Lengths
 Purpose:
   - Extract the longest motif lengths from a VCF file.
---------------------------------------------
"""

import gzip
import re
import argparse
from collections import defaultdict
from pathlib import Path

# import shared paths/helpers from config
from trplots.config import SEQ_DATA, ENSURE_DIR

def parse_args():
    p = argparse.ArgumentParser(description="Extract longest motif lengths from a VCF")
    p.add_argument("--vcf", type=str, default=str(SEQ_DATA / "83_loci_503_samples" / "1000g-ont-strchive-83_loci_503_samples.vcf.gz"),
                   help="Path to input VCF (gzipped)")
    p.add_argument("--out", type=str, default=None, help="Output CSV path (default: next to VCF)")
    return p.parse_args()

def extract_longest_motifs(vcf_path: Path):
    longest_motifs = defaultdict(lambda: {"motif": "", "length": 0})

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

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

    return longest_motifs

def main(args=None):
    if args is None:
        args = parse_args()

    vcf_path = Path(args.vcf)
    if not vcf_path.exists():
        print(f"ERROR: VCF not found: {vcf_path}")
        return

    out_path = Path(args.out) if args.out else (vcf_path.parent / "longest_motifs.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    longest_motifs = extract_longest_motifs(vcf_path)

    with out_path.open("w", encoding="utf-8") as out:
        out.write("Gene_or_Locus,Motif,Length\n")
        for gene, data in sorted(longest_motifs.items()):
            out.write(f"{gene},{data['motif']},{data['length']}\n")

    print(f"Saved {len(longest_motifs)} motifs to {out_path}")

if __name__ == "__main__":
    main()
