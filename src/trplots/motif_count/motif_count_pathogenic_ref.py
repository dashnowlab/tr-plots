#!/usr/bin/env python3
"""
---------------------------------------------
 Script: Motif Metrics Calculator (uses JSON pathogenic motif)
 Purpose:
   - Compute motif-aware metrics from STR VCF.
   - Uses locus JSON's pathogenic_motif_gene_orientation (gene orientation)
     as the motif when available (e.g., DAB1 -> "ATTT C" visual kerning ≈ "ATTT C" -> "ATTT C" -> "ATTT C" sans spaces).
   - Falls back to VCF INFO.MOTIF if JSON motif not found.
   - Outputs allele-level motif repeat statistics (+ MotifSource for auditing).
---------------------------------------------
"""

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Any

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

# Loci JSON (default)
LOCI_JSON = BASE_DIR / "data" / "other_data" / "strchive-loci.json"

# Output root directory
OUT_ROOT = BASE_DIR / "results" / "motif_counts_pathogenic_ref"
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


def write_table(out_prefix: Path, rows: List[Dict[str, object]], name: Optional[str] = None):
    """
    Write rows to CSV and (optionally) XLSX.
    If `name` is provided, files are written as out_prefix.{name}.csv/.xlsx
    Otherwise files are written as out_prefix.csv/.xlsx
    """
    if not rows:
        return

    headers = list(rows[0].keys())
    if name:
        csv_path = out_prefix.with_suffix(f".{name}.csv")
        xlsx_path = out_prefix.with_suffix(f".{name}.xlsx")
    else:
        csv_path = out_prefix.with_suffix(".csv")
        xlsx_path = out_prefix.with_suffix(".xlsx")

    with open(csv_path, "w") as f:
        f.write(",".join(headers) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(h, "")) for h in headers) + "\n")

    if pd is not None:
        try:
            pd.DataFrame(rows).to_excel(xlsx_path, index=False)
        except Exception:
            pass


# ---------- JSON motif loading ----------
def _clean_motif(m: str) -> str:
    """Uppercase, strip whitespace, and remove spaces from a motif."""
    if not isinstance(m, str):
        return ""
    return re.sub(r"\s+", "", m.upper())


def _extract_gene(obj: Dict[str, Any]) -> Optional[str]:
    """Try common keys to find a gene symbol in a JSON object."""
    for k in ("gene", "gene_symbol", "Gene", "GeneSymbol", "name", "locus", "Locus"):
        if k in obj and isinstance(obj[k], str) and obj[k].strip():
            return obj[k].strip().upper()
    return None


def load_pathogenic_motif_map(json_path: Path) -> Dict[str, str]:
    """
    Build {GENE -> motif} map from loci JSON.

    Looks for key 'pathogenic_motif_reference_orientation' (list[str] or str).
    If multiple motifs present, pick the longest (usually the clinically relevant expansion unit).
    Traverses lists/dicts recursively so it’s resilient to schema differences.
    """
    if not json_path.exists():
        return {}

    with open(json_path, "r") as f:
        data = json.load(f)

    gene_to_motif: Dict[str, str] = {}

    def visit(node: Any, current_gene: Optional[str] = None):
        if isinstance(node, dict):
            gene_here = _extract_gene(node) or current_gene

            # Primary key we care about
            if "pathogenic_motif_reference_orientation" in node:
                raw = node["pathogenic_motif_reference_orientation"]
                motifs: List[str] = []
                if isinstance(raw, list):
                    motifs = [_clean_motif(x) for x in raw if isinstance(x, str)]
                elif isinstance(raw, str):
                    motifs = [_clean_motif(raw)]
                motifs = [m for m in motifs if m]  # drop empties
                if motifs:
                    # choose the longest motif (more specific), else first
                    chosen = sorted(motifs, key=len, reverse=True)[0]
                    if gene_here:
                        gene_to_motif.setdefault(gene_here, chosen)

            # Recurse
            for v in node.values():
                visit(v, gene_here)

        elif isinstance(node, list):
            for v in node:
                visit(v, current_gene)

    visit(data, None)
    return gene_to_motif


# ---------- Loci JSON loading for gene/disease lookup ----------
def _norm_chrom(c: str) -> str:
    if not isinstance(c, str):
        return str(c)
    return c.lower().replace("chr", "").lstrip("0")  # normalize to e.g. '1', 'X'


def load_loci_list(json_path: Path) -> List[Dict[str, Any]]:
    """
    Load loci JSON and return list of entries with normalized chrom, numeric start/stop, gene, disease.
    """
    if not json_path.exists():
        return []
    with open(json_path, "r") as f:
        data = json.load(f)
    entries: List[Dict[str, Any]] = []

    def visit(node: Any):
        if isinstance(node, dict):
            chrom = node.get("chrom") or node.get("chromosome") or node.get("Chr") or None
            start = node.get("start_hg38") or node.get("start") or None
            stop = node.get("stop_hg38") or node.get("stop") or None
            gene = _extract_gene(node) or None
            disease = node.get("disease") or node.get("Disease") or None

            if chrom and start is not None and stop is not None:
                try:
                    entries.append({
                        "chrom": _norm_chrom(str(chrom)),
                        "start": int(start),
                        "stop": int(stop),
                        "gene": gene,
                        "disease": disease
                    })
                except Exception:
                    pass

            for v in node.values():
                visit(v)
        elif isinstance(node, list):
            for v in node:
                visit(v)

    visit(data)
    return entries


# ---------- Gene extraction from VCF ----------
GENE_KEYS_IN_INFO = ("GENE", "GENE_SYMBOL", "Gene", "GeneSymbol", "LOCUS", "Locus")

GENE_FROM_ID_RE = re.compile(r"^[A-Za-z0-9_.:-]*([A-Z0-9]{2,})")


def gene_from_vcf_record(info: Dict[str, str], record_id: str) -> Optional[str]:
    # 1) Try standard INFO keys
    for k in GENE_KEYS_IN_INFO:
        if k in info and info[k].strip():
            return info[k].strip().upper()

    # 2) Try to pull a gene-like token out of the ID (e.g., "DAB1:..." or "...|DAB1|...")
    if record_id and isinstance(record_id, str):
        # favor all-caps gene tokens separated by common delimiters
        tokens = re.split(r"[:|/_\-.,]", record_id)
        for t in tokens:
            tt = t.strip().upper()
            if re.fullmatch(r"[A-Z0-9]{2,}", tt):
                return tt

    return None


# ---------- Main processing ----------
def process_vcf(
    vcf_path: Path,
    out_prefix: Path,
    loci_json: Optional[Path] = None,
    per_sample: bool = True,   # default -> generate one row per sample allele (2 alleles/sample)
    pathogenic_threshold: Optional[float] = None,
):
    # Load JSON motif map
    motif_map = load_pathogenic_motif_map(loci_json) if loci_json else {}
    # Load loci list for gene/disease lookup
    loci_list = load_loci_list(loci_json) if loci_json else []

    allele_rows: List[Dict[str, object]] = []
    sample_rows: List[Dict[str, object]] = []

    with open_maybe_gzip(vcf_path) as f:
        samples: List[str] = []
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

            # Determine gene and motif preference
            # First try to find locus in JSON by chrom/pos
            chrom_norm = _norm_chrom(CHROM)
            try:
                pos_int = int(POS)
            except Exception:
                pos_int = None

            matched = None
            if pos_int is not None:
                for loc in loci_list:
                    if loc["chrom"] == chrom_norm and loc["start"] <= pos_int <= loc["stop"]:
                        matched = loc
                        break

            if matched:
                gene = (matched.get("gene") or "").upper() if matched.get("gene") else ""
                disease = matched.get("disease") or ""
            else:
                # fallback to existing heuristics
                gene = gene_from_vcf_record(info, ID) or ""
                disease = ""

            json_motif = _clean_motif(motif_map.get(gene, "")) if gene else ""

            if json_motif:
                motif = json_motif
                motif_source = "json_pathogenic"
            else:
                motif = _clean_motif(info.get("MOTIF", ""))
                motif_source = "vcf_INFO.MOTIF" if motif else "none"

            motif_len = len(motif)

            AN = info.get("AN", "")
            ACs = info.get("AC", "").split(",") if "AC" in info else []
            alts = ALT.split(",") if ALT else []

            # If per_sample requested, expand to one row per sample allele (diploid -> ~2 * n_samples rows)
            if per_sample and samples and sample_data:
                fmt_keys = format_str.split(':') if format_str else []
                for si, sd in enumerate(sample_data):
                    sample_name = samples[si] if si < len(samples) else f"sample_{si}"
                    vals = sd.split(':') if sd else []
                    samp_map = dict(zip(fmt_keys, vals))

                    gt = samp_map.get("GT", "") or samp_map.get("GQ", "")  # prefer GT; fallback harmless
                    # normalize GT delimiters
                    gt = gt.replace("|", "/")
                    alleles_idx = gt.split("/") if gt else []
                    # if no GT, try to fall back to AL or other allele-level field; otherwise skip
                    if not alleles_idx or alleles_idx == ['.']:
                        # try AL (allele lengths) if present as single values per sample (comma/semicolon separated)
                        al_field = samp_map.get("AL") or samp_map.get("ALLELE_LENGTH") or ""
                        al_values = re.split(r"[,;]", al_field) if al_field else []
                        if al_values:
                            # create one row per reported AL (treat as haplotypes)
                            for h_idx, alv in enumerate(al_values, start=1):
                                try:
                                    allele_len = int(alv)
                                except Exception:
                                    allele_len = None
                                repeat_count = (allele_len / motif_len) if (allele_len and motif_len) else None
                                allele_rows.append({
                                    "CHROM": CHROM,
                                    "POS": POS,
                                    "GENE": gene,
                                    "Disease": disease,
                                    "Sample": sample_name,
                                    "Haplotype": h_idx,
                                    "AlleleSeq": "",
                                    "AlleleLen": allele_len,
                                    "RepeatCount": repeat_count,
                                    "Motif": motif,
                                    "MotifLen": motif_len
                                })
                        # nothing to do for this sample
                        continue

                    # iterate allele indices from GT (one row per allele)
                    for h_idx, aidx in enumerate(alleles_idx, start=1):
                        if aidx in ('.', ''):
                            continue
                        try:
                            ai = int(aidx)
                        except Exception:
                            continue
                        if ai == 0:
                            seq = REF
                        else:
                            alt_index = ai - 1
                            seq = alts[alt_index] if alt_index < len(alts) else ""

                        allele_len = len(seq) if seq else None
                        # if sample AL field exists and provides length for this haplotype, prefer it
                        if "AL" in samp_map and samp_map["AL"]:
                            al_vals = re.split(r"[,;]", samp_map["AL"])
                            if len(al_vals) >= h_idx:
                                try:
                                    allele_len = int(al_vals[h_idx - 1])
                                except Exception:
                                    pass

                        repeat_count = (allele_len / motif_len) if (allele_len and motif_len) else None
                        allele_rows.append({
                            "CHROM": CHROM,
                            "POS": POS,
                            "GENE": gene,
                            "Disease": disease,
                            "Sample": sample_name,
                            "Haplotype": h_idx,
                            "AlleleSeq": seq,
                            "AlleleLen": allele_len,
                            "RepeatCount": repeat_count,
                            "Motif": motif,
                            "MotifLen": motif_len,
                            "MotifSource": motif_source,
                        })
            else:
                # fallback: keep locus-level REF+ALT summary rows (previous behavior)
                allele_rows.append({
                    "CHROM": CHROM,
                    "POS": POS,
                    "GENE": gene,
                    "Disease": disease,
                    "AlleleType": "REF",
                    "AlleleSeq": REF,
                    "Motif": motif,
                    "MotifLen": motif_len,
                    "MotifSource": motif_source,
                    "AlleleLen": len(REF),
                    "ApproxRepeatCount": (len(REF) / motif_len) if motif_len else None,
                    "PureCount": full_pure_count(REF, motif) if motif_len else 0,
                    "MaxMotifRun": max_consecutive_motif_run(REF, motif) if motif_len else 0,
                    "IsPure": (full_pure_count(REF, motif) > 0) if motif_len else False,
                    "AC": "",
                    "AN": AN,
                })

                for i, alt in enumerate(alts):
                    allele_rows.append({
                        "CHROM": CHROM,
                        "POS": POS,
                        "GENE": gene,
                        "Disease": disease,
                        "AlleleType": f"ALT{i+1}",
                        "AlleleSeq": alt,
                        "Path Ref Motif": motif,
                        "MotifLen": motif_len,
                        "MotifSource": motif_source,
                        "AlleleLen": len(alt),
                        "ApproxRepeatCount": (len(alt) / motif_len) if motif_len else None,
                        "PureCount": full_pure_count(alt, motif) if motif_len else 0,
                        "MaxMotifRun": max_consecutive_motif_run(alt, motif) if motif_len else 0,
                        "IsPure": (full_pure_count(alt, motif) > 0) if motif_len else False,
                        "AC": ACs[i] if i < len(ACs) else "",
                        "AN": AN,
                    })

    # --- Write output ---
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    # write without the ".alleles" suffix and use the CLI-provided prefix.
    write_table(out_prefix, allele_rows, name=None)
    print(f"--- Saved allele-level motif data to: {out_prefix}.csv/.xlsx ---")

# ---------- CLI ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute motif-aware metrics from STR VCF (preferring JSON pathogenic motifs).")
    parser.add_argument("--vcf", default=str(VCF_FILE), help="Path to VCF (.vcf or .vcf.gz)")
    parser.add_argument("--loci-json", default=str(LOCI_JSON), help="Path to loci JSON with pathogenic_motif_gene_orientation")
    parser.add_argument("--out", default=str(OUT_ROOT / "motif_metrics_pathogenic_ref"), help="Output file prefix")
    parser.add_argument("--per-sample", action="store_true", help="Expand output to one row per sample allele (default: enabled)")
    args = parser.parse_args()

    process_vcf(Path(args.vcf), Path(args.out), Path(args.loci_json), per_sample=True)
