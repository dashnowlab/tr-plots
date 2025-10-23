"""
---------------------------------------------
 Script: Creating Master Allele Spreadsheet (fixed AL handling)
 Purpose:
    - Integrate VCF and locus metadata into a comprehensive allele spreadsheet
    - Include demographic and technology data from sample summary CSV
    - Save as Excel file with clear column structure

 Key fixes vs prior version:
    - Interpret FORMAT/AL as *per-haplotype lengths* aligned to GT haplotypes,
      NOT as a global map keyed by ALT index.
    - Robust handling of phased/unphased GT, haploid calls, and missing alleles.
    - Optional: drop one or more samples by id via --drop-sample (repeatable).
---------------------------------------------
"""

import argparse
import gzip
import json
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd
from openpyxl.utils import get_column_letter  # noqa: F401  (kept for parity; not used directly here)


# ---------- Defaults (override via CLI) ----------
VCF_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/sequencing_data/81_loci_502_samples/1000g-ont-strchive-81_loci_502_samples_81224_alleles.vcf.gz'
LOCI_JSON = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive-loci.json'
SAMPLE_SUMMARY_CSV = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/1kgp_ont_500_summary_-_sheet1.csv'
OUTPUT_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/master_allele_spreadsheet_copy_81_loci_502_samples_copy.xlsx'


# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description="Create integrated master allele spreadsheet")
    p.add_argument("--vcf", type=str, default=VCF_FILE, help="Input VCF (gzipped or plain)")
    p.add_argument("--loci-json", type=str, default=LOCI_JSON, help="Loci JSON metadata")
    p.add_argument("--sample-csv", type=str, default=SAMPLE_SUMMARY_CSV, help="Sample summary CSV (header on row 2)")
    p.add_argument("--output", type=str, default=OUTPUT_FILE, help="Output Excel file path")
    p.add_argument("--drop-sample", action="append", default=[],
                   help="Sample ID(s) to drop (exact matches to VCF header sample names). Can be used multiple times.")
    return p.parse_args()


# ---------- Data Loading ----------
def load_loci_data(json_file: str) -> Dict[str, dict]:
    """
    Load STR loci info from JSON into a dict with lookups by 'id' and 'chrom:pos'.
    """
    print("Loading Loci Data...")
    with open(json_file, 'r') as f:
        data = json.load(f)

    loci_map: Dict[str, dict] = {}
    missing_flank = 0

    for entry in data:
        flank_motif_value = entry.get('flank_motif')
        if flank_motif_value is None:
            missing_flank += 1
        elif not isinstance(flank_motif_value, str):
            flank_motif_value = str(flank_motif_value)

        entry_data = {
            'Gene': entry.get('gene'),
            'Disease': entry.get('disease'),
            'Disease ID': entry.get('disease_id'),
            'Motif': (entry.get('reference_motif_reference_orientation') or [None])[0],
            'Motif length': entry.get('motif_len'),
            'Benign min': entry.get('benign_min'),
            'Benign max': entry.get('benign_max'),
            'Pathogenic min': entry.get('pathogenic_min'),
            'Pathogenic max': entry.get('pathogenic_max'),
            'Inheritance': ', '.join(entry.get('inheritance', [])),
            'Flank Motif Structure': flank_motif_value
        }

        id_key = entry['id']
        pos_key = f"{entry['chrom'].replace('chr', '')}:{entry['start_hg38']}"
        loci_map[id_key] = entry_data
        loci_map[pos_key] = entry_data

    print(f"Loaded {len(data)} unique loci (with {len(loci_map)} total lookup keys).")
    if missing_flank > 0:
        print(f"DIAGNOSTIC: {missing_flank} loci missing 'Flank Motif Structure' in JSON.")
    return loci_map


def load_all_sample_info(summary_csv_file: str) -> Dict[str, dict]:
    """
    Loads sample demographic and technology data directly from the summary CSV.
    The header for the required columns is on the second row (index 1).
    """
    print(f"Loading All Sample Data from {summary_csv_file}...")
    df_summary = pd.read_csv(summary_csv_file, header=1)

    df_summary.rename(columns={
        'Sample_ID': 'Sample ID Cleaned',
        'Sex': 'Sex',
        'SubPop': 'SubPop',
        'SuperPop': 'SuperPop',
        'Pore': 'Pore'
    }, inplace=True)

    columns_to_keep = ['Sample ID Cleaned', 'SubPop', 'SuperPop', 'Sex', 'Pore']
    df_sample_info = (
        df_summary[columns_to_keep]
        .drop_duplicates(subset=['Sample ID Cleaned'], keep='first')
        .set_index('Sample ID Cleaned')
    )

    for col in ['Sex', 'SubPop', 'SuperPop', 'Pore']:
        df_sample_info[col] = df_sample_info[col].fillna('Unknown')

    sample_map = df_sample_info.to_dict('index')
    print(f"Loaded {len(sample_map)} unique sample demographic entries.")
    return sample_map


# ---------- Helpers ----------
_GT_SPLIT_RE = re.compile(r'[\/|]')

def _safe_float_len(x) -> float:
    try:
        return float(x)
    except Exception:
        return 0.0


def _calc_motif_len(locus_data: dict, motif_from_info: Optional[str], ref: str) -> float:
    """
    Determine motif length with robust fallbacks.
    """
    motif_len = locus_data.get('Motif length')
    if isinstance(motif_len, (int, float)) and motif_len > 0:
        return float(motif_len)

    # try from locus motif string
    motif_to_use = locus_data.get('Motif') or motif_from_info
    if motif_to_use:
        return float(len(motif_to_use))

    # fallback to ref length (harmless default)
    return float(len(ref)) if ref else 0.0


def _hap_lengths_from_AL(AL_field: str) -> List[Optional[int]]:
    """
    Parse FORMAT/AL into per-haplotype integer lengths.
    Returns a list of up to 2 elements (diploid), padded with None if missing.
    """
    if not AL_field or AL_field == '.':
        return [None, None]
    toks = [t for t in AL_field.split(',') if t != '' and t != '.']
    out: List[Optional[int]] = []
    for t in toks[:2]:
        try:
            out.append(int(t))
        except Exception:
            out.append(None)
    while len(out) < 2:
        out.append(None)
    return out


def _clean_sample_id(raw: str) -> str:
    """
    Extract HGxxxx / GMxxxx / NAxxxx if present; else use leading token before - or _.
    """
    m = re.match(r'^(HG\d+|GM\d+|NA\d+)', raw)
    return m.group(0) if m else raw.split('-')[0].split('_')[0]


# ---------- Core ----------
def parse_vcf_and_merge(
    vcf_file: str,
    loci_map: Dict[str, dict],
    sample_map: Dict[str, dict],
    drop_samples: List[str]
) -> pd.DataFrame:
    """
    Parse VCF into an allele-level table and merge locus + demographic metadata.
    Uses per-haplotype AL lengths aligned with GT haplotypes.
    """
    print("Parsing VCF and building allele table...")
    allele_records: List[dict] = []
    vcf_header_sample_ids: List[str] = []

    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    DEFAULT_LOCUS_DATA = {
        'Gene': 'Unknown', 'Disease': 'Unknown', 'Disease ID': 'Unknown',
        'Motif': None, 'Motif length': 0,
        'Benign min': None, 'Benign max': None,
        'Pathogenic min': None, 'Pathogenic max': None,
        'Inheritance': 'Unknown', 'Flank Motif Structure': None
    }

    total_sites = 0

    try:
        with opener(vcf_file, mode, encoding='utf-8') as f:
            for line in f:
                if line.startswith('##'):
                    continue

                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    vcf_header_sample_ids = header[9:]
                    print(f"VCF Header Sample Count: {len(vcf_header_sample_ids)}")

                    if drop_samples:
                        missing = [s for s in drop_samples if s not in vcf_header_sample_ids]
                        if missing:
                            print(f"NOTE: --drop-sample not in header (will be ignored): {missing}")

                    continue

                if not vcf_header_sample_ids:
                    continue  # not past header yet

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 10:
                    continue

                chrom, pos, identifier, ref, alt_str, qual, filt, info_str, fmt_str = fields[:9]
                sample_fields = fields[9:]
                total_sites += 1

                # INFO -> dict
                info_dict = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info_dict[k] = v
                motif_from_info = info_dict.get('MOTIF')

                # locus lookup by id or chrom:pos (chrom without "chr")
                locus_key_id = identifier if identifier and identifier != '.' else None
                locus_key_pos = f"{chrom.replace('chr', '')}:{pos}"

                locus_data = loci_map.get(locus_key_id) or loci_map.get(locus_key_pos) or DEFAULT_LOCUS_DATA
                if locus_data is DEFAULT_LOCUS_DATA and motif_from_info and not locus_data['Motif']:
                    locus_data = dict(DEFAULT_LOCUS_DATA)
                    locus_data['Motif'] = motif_from_info

                motif_len = _calc_motif_len(locus_data, motif_from_info, ref)
                fmt_keys = fmt_str.split(':')

                # iterate samples
                for i, sample_field in enumerate(sample_fields):
                    sample_id_full = vcf_header_sample_ids[i]
                    if sample_id_full in drop_samples:
                        continue  # skip dropped samples entirely

                    sample_id_clean = _clean_sample_id(sample_id_full)

                    demo = sample_map.get(sample_id_clean, {
                        'SubPop': 'Unknown', 'SuperPop': 'Unknown', 'Sex': 'Unknown', 'Pore': 'Unknown'
                    })

                    fmt_vals = sample_field.split(':')
                    fmt_data = dict(zip(fmt_keys, fmt_vals))

                    gt = fmt_data.get('GT', './.')
                    if gt == '.' or gt == './.' or gt == '.|.':
                        continue  # fully missing

                    gt_toks = _GT_SPLIT_RE.split(gt)[:2]  # diploid max
                    while len(gt_toks) < 2:
                        gt_toks.append('.')  # pad if haploid

                    al_haps = _hap_lengths_from_AL(fmt_data.get('AL', '.'))

                    # per haplotype allele
                    for hap_idx, allele_tok in enumerate(gt_toks):
                        if allele_tok == '.':
                            continue  # this haplotype missing

                        if allele_tok == '0':
                            allele_len = len(ref) if ref else None
                        else:
                            allele_len = al_haps[hap_idx]  # per-haplotype length

                        if allele_len is None:
                            continue  # cannot resolve; skip just this haplotype

                        repeat_count = round(allele_len / motif_len, 2) if motif_len and motif_len > 0 else 'N/A'

                        allele_records.append({
                            'Chromosome': chrom.replace('chr', ''),
                            'Position (Start)': pos,
                            'Gene': locus_data['Gene'],
                            'Disease': locus_data['Disease'],
                            'Disease ID': locus_data['Disease ID'],
                            'Motif': locus_data['Motif'],
                            'Motif length': motif_len,
                            'Flank Motif Structure': locus_data['Flank Motif Structure'],
                            'Benign min': locus_data['Benign min'],
                            'Benign max': locus_data['Benign max'],
                            'Pathogenic min': locus_data['Pathogenic min'],
                            'Pathogenic max': locus_data['Pathogenic max'],
                            'Inheritance': locus_data['Inheritance'],
                            'Allele Length (bp)': allele_len,
                            'Repeat count (Calc)': repeat_count,
                            'Sample ID (Raw)': sample_id_full,
                            'Sample ID (Cleaned)': sample_id_clean,
                            'SubPop': demo.get('SubPop', 'Unknown'),
                            'SuperPop': demo.get('SuperPop', 'Unknown'),
                            'Sex': demo.get('Sex', 'Unknown'),
                            'Pore': demo.get('Pore', 'Unknown'),
                        })

    except FileNotFoundError:
        print(f"Error: VCF not found: {vcf_file}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Unexpected error while parsing VCF: {e}")
        return pd.DataFrame()

    if not allele_records:
        print("No allele records parsed.")
        return pd.DataFrame()

    df = pd.DataFrame(allele_records)
    n_unique_samples = df['Sample ID (Raw)'].nunique()
    print(f"Finished parsing VCF.")
    print(f"  Total variant records parsed: {total_sites}")
    print(f"  Generated allele rows: {len(df)}")
    print(f"  Unique samples represented: {n_unique_samples}")
    return df


# ---------- Orchestration ----------
def create_integrated_spreadsheet(loci_json: str, vcf_file: str, sample_csv: str,
                                  output_file: str, drop_samples: List[str]) -> Optional[str]:
    loci_map = load_loci_data(loci_json)
    sample_map = load_all_sample_info(sample_csv)

    df_integrated = parse_vcf_and_merge(vcf_file, loci_map, sample_map, drop_samples)

    if df_integrated.empty:
        print("No data parsed; skipping Excel write.")
        return None

    # Final column selection / rename
    column_rename_map = {
        'Chromosome': 'Chromosome',
        'Position (Start)': 'Position',
        'Gene': 'Gene',
        'Disease': 'Disease',
        'Disease ID': 'Disease ID',
        'Motif': 'Motif',
        'Motif length': 'Motif length',
        'Allele Length (bp)': 'Allele length',
        'Repeat count (Calc)': 'Repeat count',
        'Benign min': 'Benign min',
        'Benign max': 'Benign max',
        'Pathogenic min': 'Pathogenic min',
        'Pathogenic max': 'Pathogenic max',
        'Inheritance': 'Inheritance',
        'Sample ID (Raw)': 'Sample ID',
        'Sample ID (Cleaned)': 'Sample ID Cleaned',
        'SubPop': 'SubPop',
        'SuperPop': 'SuperPop',
        'Sex': 'Sex',
        'Flank Motif Structure': 'Flank Motif',
        'Pore': 'Pore'
    }

    present = [c for c in column_rename_map.keys() if c in df_integrated.columns]
    df_final = df_integrated[present].rename(columns=column_rename_map)
    desired_order = list(column_rename_map.values())
    df_final = df_final.reindex(columns=desired_order)

    # Write Excel
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_final.to_excel(writer, index=False, sheet_name='Integrated Alleles')

    print(f"\nSuccess! Integrated spreadsheet saved to '{output_file}'.")
    return output_file


# ---------- Main ----------
def main():
    args = parse_args()
    out = create_integrated_spreadsheet(
        loci_json=args.loci_json,
        vcf_file=args.vcf,
        sample_csv=args.sample_csv,
        output_file=args.output,
        drop_samples=args.drop_sample,
    )
    if out:
        print(f"Created spreadsheet: {out}")
    else:
        print("No spreadsheet created.")


if __name__ == "__main__":
    main()
