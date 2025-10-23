#!/usr/bin/env python3
"""
Find longest tandem repeat motif in REF and ALT allele sequences from a gzipped VCF.

Outputs a CSV with one row per allele (REF and each ALT):
  locus, allele_label, allele_seq, motif_unit, unit_len, repeat_count, total_repeat_len, start0, end0

Options:
  --max-unit     max motif unit length to test (default 10)
  --min-repeats  minimum repeats to report (default 2; set to 1 to capture single occurrences)
  --top          print top N loci by total_repeat_len to stdout (default 20)
"""
import argparse
import gzip
import re
from pathlib import Path
from typing import Tuple, Optional

def find_longest_tandem_repeat(seq: str, max_unit_len: int = 10, min_repeats: int = 2) -> Tuple[Optional[str], int, int, int, Optional[int], Optional[int]]:
    """
    Search seq for the longest tandem repeat.

    Returns:
      (motif_unit, unit_len, repeat_count, total_repeat_len, start_index0, end_index0)
      If none found meeting min_repeats, returns (None, 0, 0, 0, None, None).
    """
    seq = (seq or "").upper()
    n = len(seq)
    if n == 0:
        return (None, 0, 0, 0, None, None)

    best = (None, 0, 0, 0, None, None)  # motif, unit_len, repeat_count, total_len, start, end
    max_unit_len = min(max_unit_len, n)

    valid_base_re = re.compile(r'^[ACGT]+$')

    for unit_len in range(1, max_unit_len + 1):
        # gather unique candidate units in this allele to avoid duplicate work
        units = set()
        for i in range(0, n - unit_len + 1):
            unit = seq[i:i+unit_len]
            if not valid_base_re.match(unit):
                continue
            units.add(unit)

        if not units:
            continue

        # pattern: one or more repeats; enforce min_repeats if >1
        for unit in units:
            if min_repeats > 1:
                pattern = re.compile(r'(?:' + re.escape(unit) + r')' + '{' + str(min_repeats) + r',}')
            else:
                pattern = re.compile(r'(?:' + re.escape(unit) + r')+')

            for m in pattern.finditer(seq):
                run_bases = len(m.group(0))
                repeat_count = run_bases // unit_len
                start = m.start()
                end = m.end()  # python slice end (exclusive)
                # choose best by total run length, tiebreaker larger unit_len then earlier start
                if (run_bases > best[3]) or (run_bases == best[3] and unit_len > best[1]) or (run_bases == best[3] and unit_len == best[1] and (best[4] is None or start < best[4])):
                    best = (unit, unit_len, repeat_count, run_bases, start, end)

    return best

def parse_info_field(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out

def process_vcf(vcf_path: Path, out_csv: Path, max_unit_len: int = 10, min_repeats: int = 2, top_n: int = 20):
    header = ["locus", "allele_label", "allele_seq", "motif_unit", "unit_len", "repeat_count", "total_repeat_len", "start0", "end0"]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    with gzip.open(vcf_path, "rt") as fh, out_csv.open("w", encoding="utf-8") as out:
        out.write(",".join(header) + "\n")
        for ln in fh:
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, vid, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]
            locus = vid if vid and vid != "." else f"{chrom}:{pos}"
            alts = [a for a in alt_field.split(",") if a and a != "."]

            alleles = [("REF", ref)] + [(f"ALT{i+1}", a) for i,a in enumerate(alts)]
            for allele_label, seq in alleles:
                motif, unit_len, repeat_count, total_len, start0, end0 = find_longest_tandem_repeat(seq, max_unit_len=max_unit_len, min_repeats=min_repeats)
                motif_str = motif if motif is not None else ""
                start_s = "" if start0 is None else str(start0)
                end_s = "" if end0 is None else str(end0)
                out_row = [locus, allele_label, seq, motif_str, str(unit_len), str(repeat_count), str(total_len), start_s, end_s]
                out.write(",".join(x.replace(",", ";") for x in out_row) + "\n")  # replace commas in allele seqs
                rows.append((locus, allele_label, motif_str, unit_len, repeat_count, total_len, start0, end0))

    # print brief top-N summary
    rows_sorted = sorted(rows, key=lambda x: x[5], reverse=True)  # by total_repeat_len
    print(f"Wrote detailed results to {out_csv}")
    print(f"Top {top_n} alleles by total_repeat_len (locus, allele_label, motif, unit_len, repeat_count, total_len, start0, end0):")
    for r in rows_sorted[:top_n]:
        print(",".join(map(str, r)))

def main():
    p = argparse.ArgumentParser(description="Find longest motif unit in REF/ALT allele sequences from VCF")
    p.add_argument("--vcf", type=str, required=True, help="Path to gzipped VCF")
    p.add_argument("--out", type=str, default=None, help="Output CSV path")
    p.add_argument("--max-unit", type=int, default=10, help="Maximum motif unit length to test")
    p.add_argument("--min-repeats", type=int, default=2, help="Minimum number of repeats to report (set 1 to include single occurrences)")
    p.add_argument("--top", type=int, default=20, help="How many top alleles to print as summary")
    args = p.parse_args()

    vcf_path = Path(args.vcf)
    if not vcf_path.exists():
        raise SystemExit(f"VCF not found: {vcf_path}")
    out_path = Path(args.out) if args.out else (vcf_path.parent / "allele_motif_units.csv")
    process_vcf(vcf_path, out_path, max_unit_len=args.max_unit, min_repeats=args.min_repeats, top_n=args.top)

if __name__ == "__main__":
    main()
