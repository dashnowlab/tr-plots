"""
Master Allele Spreadsheet Creator

Purpose:
    Parse a VCF of STR loci (gzipped or plain), merge locus metadata and sample
    demographics (CSV + 1kGP), and produce an integrated per-haplotype Excel
    spreadsheet with allele lengths, repeat counts and QC flags.

Key behaviors / features:
    - Reads loci JSON, sample summary CSV (header row 2), and 1kGP sample-info (tab-delim).
    - Parses VCF GT and per-haplotype AL fields; supports REF (0) and ALT alleles.
    - Converts missing or <=0 AL values to NaN and records QC reasons.
    - Derives motif length from loci metadata, info.MOTIF, or REF sequence.
    - Calculates repeat count = allele_length / motif_length (rounded to 2 decimals).
    - Emits one row per non-missing haplotype with Haplotype (1/2), Sample ID raw/cleaned,
      demographics (SubPop/SuperPop/Sex/Pore), locus annotations, QC Flag and QC Reasons.
    - Backfills SubPop/SuperPop/Sex from 1kGP data when CSV values are unknown/missing.
    - Flags loci whose Inheritance contains both AD and AR with QC reason
      'inheritance_AD_and_AR'.
    - Supports --drop-sample (repeatable) to skip exact VCF header sample names.
    - Writes an Excel file (sheet 'Integrated Alleles'). Numeric NaNs for
      Allele length and Repeat count are written as the string "NaN" for visibility.
    - Expects VCF sample names often beginning with HG*, GM*, NA*; a cleaning helper
      strips common suffixes to match sample CSV / 1kGP keys.
"""

import argparse
import gzip
import json
import os
import re
from typing import Dict, List, Optional

import pandas as pd
from openpyxl.utils import get_column_letter  # noqa: F401

# ---------- Defaults (override via CLI) ----------
VCF_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/sequencing_data/81_loci_502_samples/1000g-ont-strchive-81_loci_502_samples_81224_alleles.vcf.gz'
LOCI_JSON = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive-loci.json'
SAMPLE_SUMMARY_CSV = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/1kgp_ont_500_summary_-_sheet1.csv'
KGP_SAMPLE_INFO = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/1000_genomes_20130606_sample_info.txt'
OUTPUT_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/all_path_ref_motif_allele_spreadsheet.xlsx'

# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description="Create integrated master allele spreadsheet")
    p.add_argument("--vcf", type=str, default=VCF_FILE, help="Input VCF (gzipped or plain)")
    p.add_argument("--loci-json", type=str, default=LOCI_JSON, help="Loci JSON metadata")
    p.add_argument("--sample-csv", type=str, default=SAMPLE_SUMMARY_CSV, help="Sample summary CSV (header on row 2)")
    p.add_argument("--kgp-sample-info", type=str, default=KGP_SAMPLE_INFO,
                   help="1000 Genomes sample-info TXT (20130606)")
    p.add_argument("--output", type=str, default=OUTPUT_FILE, help="Output Excel file path")
    p.add_argument("--drop-sample", action="append", default=[],
                   help="Sample ID(s) to drop (exact match to VCF header sample names). Repeatable.")
    return p.parse_args()

# ---------- Data Loading ----------
def load_loci_data(json_file: str) -> Dict[str, dict]:
    with open(json_file, 'r') as f:
        data = json.load(f)

    loci_map: Dict[str, dict] = {}
    for entry in data:
        flank_motif_value = entry.get('flank_motif')
        if flank_motif_value is not None and not isinstance(flank_motif_value, str):
            flank_motif_value = str(flank_motif_value)

        # ONLY collect pathogenic_motif_reference_orientation (do not fall back to gene/reference keys).
        raw_pathogenic = entry.get('pathogenic_motif_reference_orientation')

        # Normalize and collect all candidate pathogenic motifs (may be multiple).
        pathogenic_motifs: List[str] = []
        if isinstance(raw_pathogenic, list):
            for m in raw_pathogenic:
                if isinstance(m, str) and m.strip():
                    pathogenic_motifs.append(re.sub(r"\s+", "", m).upper())
        elif isinstance(raw_pathogenic, str) and raw_pathogenic.strip():
            pathogenic_motifs.append(re.sub(r"\s+", "", raw_pathogenic).upper())
        # Deduplicate while preserving order
        seen_m = set()
        pathogenic_motifs = [m for m in pathogenic_motifs if not (m in seen_m or seen_m.add(m))]
        # If no pathogenic motifs are present, keep list empty and mark Motif as 'Unknown'
 
        entry_data = {
            'Gene': entry.get('gene'),
            'Disease': entry.get('disease'),
            'Disease ID': entry.get('disease_id'),
            # keep first motif for backward compatibility; if none present mark 'Unknown'
            'Motif': pathogenic_motifs[0] if pathogenic_motifs else 'Unknown',
            'Pathogenic Motifs': pathogenic_motifs,
            # original 'Motif length' kept as locus annotation (may be None)
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
    print(f"Loaded loci: {len(data)} (lookup keys: {len(loci_map)})")
    return loci_map

def load_all_sample_info(summary_csv_file: str) -> Dict[str, dict]:
    """
    Load demographics/tech from your CSV (header is on row 2).
    """
    print(f"Loading summary CSV: {summary_csv_file}")
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

    result = df_sample_info.to_dict('index')
    print(f"Summary CSV entries: {len(result)}")
    return result

# --- 1000G population → superpopulation map ---
_POP_TO_SUPERPOP = {
    # AFR
    'YRI': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR',
    # AMR
    'CLM': 'AMR', 'MXL': 'AMR', 'PEL': 'AMR', 'PUR': 'AMR',
    # EAS
    'CHB': 'EAS', 'CHS': 'EAS', 'JPT': 'EAS', 'CDX': 'EAS', 'KHV': 'EAS',
    # EUR
    'CEU': 'EUR', 'TSI': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR',
    # SAS
    'GIH': 'SAS', 'PJL': 'SAS', 'BEB': 'SAS', 'STU': 'SAS', 'ITU': 'SAS',
}

def _norm_sex(x: Optional[str]) -> str:
    """Normalize to 'XX' for female, 'XY' for male, else 'Unknown'."""
    if x is None:
        return 'Unknown'
    s = str(x).strip().lower()
    if s in ['f', 'female']:
        return 'XX'
    if s in ['m', 'male']:
        return 'XY'
    return 'Unknown'

def load_1kgp_sample_info(kgp_txt_path: str) -> Dict[str, dict]:
    """
    Read 1000 Genomes 20130606 sample-info (tab-delimited).
    Returns {Sample: {'SubPop_from_1kGP', 'SuperPop_from_1kGP', 'Sex_from_1kGP'}}
    """
    print(f"Loading 1kGP sample-info: {kgp_txt_path}")
    df = pd.read_csv(kgp_txt_path, sep='\t', dtype=str)

    # Required columns
    if 'Sample' not in df.columns or 'Population' not in df.columns:
        raise ValueError("1kGP sample-info file must contain 'Sample' and 'Population' columns.")

    # Optional 'Gender' (or sometimes 'Sex'); normalize to XX/XY if present
    if 'Gender' not in df.columns:
        if 'Sex' in df.columns:
            df['Gender'] = df['Sex']
        else:
            df['Gender'] = None

    df['Population'] = df['Population'].astype(str)
    df['SuperPop_from_1kGP'] = df['Population'].map(_POP_TO_SUPERPOP).fillna('Unknown')
    df['SubPop_from_1kGP'] = df['Population']
    df['Sex_from_1kGP'] = df['Gender'].apply(_norm_sex)

    df_out = df[['Sample', 'SubPop_from_1kGP', 'SuperPop_from_1kGP', 'Sex_from_1kGP']].drop_duplicates('Sample').set_index('Sample')
    result = df_out.to_dict('index')
    print(f"1kGP entries: {len(result)}")
    return result

# ---------- Helpers ----------
_GT_SPLIT_RE = re.compile(r'[\/|]')

def _calc_motif_len(locus_data: dict, motif_from_info: Optional[str], ref: str) -> float:
    """
    Determine motif length:
      - Prefer explicit locus 'Motif' only when it's a real motif (not the sentinel 'Unknown')
      - Otherwise use motif_from_info (INFO.MOTIF)
      - Fallback to REF length
    """
    motif_len = locus_data.get('Motif length')
    if isinstance(motif_len, (int, float)) and motif_len > 0:
        return float(motif_len)

    motif_to_use = locus_data.get('Motif')
    if motif_to_use and str(motif_to_use).strip().upper() != 'UNKNOWN':
        return float(len(motif_to_use))

    if motif_from_info:
        return float(len(re.sub(r"\s+", "", str(motif_from_info)).upper()))

    return float(len(ref)) if ref else 0.0

def _hap_lengths_from_AL(AL_field: str) -> List[Optional[int]]:
    if not AL_field or AL_field == '.':
        return [None, None]
    toks = [t for t in AL_field.split(',') if t not in ['', '.']]
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
    m = re.match(r'^(HG\d+|GM\d+|NA\d+)', raw)
    return m.group(0) if m else raw.split('-')[0].split('_')[0]

# NEW: helper to detect ambiguous inheritance (AD and AR both present)
_INHER_SPLIT_RE = re.compile(r'[,/;| ]+')
def _has_both_ad_ar(inheritance_text: Optional[str]) -> bool:
    if not inheritance_text:
        return False
    toks = [t.strip().lower() for t in _INHER_SPLIT_RE.split(str(inheritance_text)) if t.strip()]
    return ('ad' in toks) and ('ar' in toks)

# ---------- Core ----------
def parse_vcf_and_merge(
    vcf_file: str,
    loci_map: Dict[str, dict],
    sample_map_csv: Dict[str, dict],
    sample_map_1kgp: Dict[str, dict],
    drop_samples: List[str]
) -> pd.DataFrame:

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
                        print(f"NOTE: --drop-sample not in header (ignored): {missing}")
                continue

            if not vcf_header_sample_ids:
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 10:
                continue

            chrom, pos, identifier, ref, alt_str, qual, filt, info_str, fmt_str = fields[:9]
            sample_fields = fields[9:]
            total_sites += 1

            # INFO
            info_dict = {}
            for item in info_str.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
            motif_from_info = info_dict.get('MOTIF')

            # locus mapping by id or chrom:pos
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
                    continue

                sample_id_clean = _clean_sample_id(sample_id_full)

                # Start with CSV demographics
                demo = {'SubPop': 'Unknown', 'SuperPop': 'Unknown', 'Sex': 'Unknown', 'Pore': 'Unknown'}
                if sample_map_csv.get(sample_id_clean):
                    demo.update(sample_map_csv[sample_id_clean])

                # Backfill from 1kGP ONLY when Unknown/empty
                kgp = sample_map_1kgp.get(sample_id_clean, {})

                def _is_unknown(x):
                    return (x is None) or (str(x).strip().lower() in ['', 'unknown', 'nan'])

                # SubPop/SuperPop
                if _is_unknown(demo.get('SubPop')) and kgp.get('SubPop_from_1kGP'):
                    demo['SubPop'] = kgp['SubPop_from_1kGP']

                if _is_unknown(demo.get('SuperPop')):
                    sp = kgp.get('SuperPop_from_1kGP')
                    if sp:
                        demo['SuperPop'] = sp

                # Sex (XX/XY backfill only)
                if _is_unknown(demo.get('Sex')) and kgp.get('Sex_from_1kGP'):
                    demo['Sex'] = kgp['Sex_from_1kGP']

                fmt_vals = sample_field.split(':')
                fmt_data = dict(zip(fmt_keys, fmt_vals))

                gt = fmt_data.get('GT', './.')
                if gt in ('.', './.', '.|.'):
                    continue

                # Split into haplotypes (diploid max)
                gt_toks = re.split(r'[\/|]', gt)[:2]
                while len(gt_toks) < 2:
                    gt_toks.append('.')

                al_haps = _hap_lengths_from_AL(fmt_data.get('AL', '.'))

                # --- inheritance ambiguity check (once per row) ---
                inheritance_text = locus_data.get('Inheritance', '')
                inheritance_ambig = _has_both_ad_ar(inheritance_text)

                for hap_idx, allele_tok in enumerate(gt_toks):
                    if allele_tok == '.':
                        continue

                    # --- Compute allele length per haplotype with QC ---
                    qc_reasons = []

                    if allele_tok == '0':
                        # reference allele -> length = len(REF)
                        allele_len = len(ref) if ref is not None else float('nan')
                        if not (isinstance(allele_len, (int, float)) and allele_len > 0):
                            qc_reasons.append('ref_len_le0_or_missing')
                            allele_len = float('nan')
                    else:
                        # ALT allele -> from AL per-haplotype
                        allele_len = al_haps[hap_idx]
                        if (allele_len is None) or (allele_len <= 0):
                            qc_reasons.append('AL_missing_or_le0')
                            allele_len = float('nan')  # force NaN for <=0/missing

                    # Repeat count: only if both allele_len and motif_len are valid (>0)
                    if motif_len and isinstance(motif_len, (int, float)) and motif_len > 0:
                        repeat_count = (round(allele_len / motif_len, 2)
                                        if (isinstance(allele_len, (int, float)) and allele_len > 0)
                                        else float('nan'))
                    else:
                        qc_reasons.append('motif_len_le0_or_missing')
                        repeat_count = float('nan')

                    # NEW: append inheritance ambiguity reason
                    if inheritance_ambig:
                        qc_reasons.append('inheritance_AD_and_AR')

                    # NEW: flag if benign range falls entirely within pathogenic range
                    # (adds a QC reason so QC Flag becomes 'Yes' for that row)
                    try:
                        bmin = locus_data.get('Benign min')
                        bmax = locus_data.get('Benign max')
                        pmin = locus_data.get('Pathogenic min')
                        pmax = locus_data.get('Pathogenic max')
                        bmin_v = float(bmin) if bmin is not None else None
                        bmax_v = float(bmax) if bmax is not None else None
                        pmin_v = float(pmin) if pmin is not None else None
                        pmax_v = float(pmax) if pmax is not None else None
                        if (bmin_v is not None and bmax_v is not None and
                                pmin_v is not None and pmax_v is not None):
                            # check that benign interval lies entirely within pathogenic interval
                            if (bmin_v >= pmin_v) and (bmax_v <= pmax_v):
                                qc_reasons.append('benign_range_within_pathogenic_range')
                    except Exception:
                        pass

                    # Determine motif candidates to expand rows: use Pathogenic Motifs if present,
                    # otherwise fall back to locus 'Motif' (single) or motif from INFO.
                    motif_candidates = locus_data.get('Pathogenic Motifs') or []
                    if not motif_candidates:
                        if locus_data.get('Motif'):
                            motif_candidates = [locus_data.get('Motif')]
                        elif motif_from_info:
                            motif_candidates = [re.sub(r"\s+", "", str(motif_from_info)).upper()]
                        else:
                            motif_candidates = []

                    # If no motif candidates, still emit a row with motif empty (preserve behavior)
                    if not motif_candidates:
                        motif_candidates = [None]

                    for motif_val in motif_candidates:
                        # compute per-row motif length (prefer explicit motif string)
                        if motif_val:
                            motif_len_used = float(len(motif_val))
                        else:
                            motif_len_used = motif_len  # previous locus-level fallback

                        # recompute repeat_count based on motif_len_used if needed
                        if motif_len_used and isinstance(motif_len_used, (int, float)) and motif_len_used > 0:
                            repeat_count_used = (round(allele_len / motif_len_used, 2)
                                                 if (isinstance(allele_len, (int, float)) and allele_len > 0)
                                                 else float('nan'))
                        else:
                            repeat_count_used = float('nan')

                        # Determine pathogenicity: leave blank if not pathogenic, 'Yes' if within pathogenic range.
                        is_pathogenic = ''
                        if isinstance(repeat_count_used, (int, float)) and not pd.isna(repeat_count_used):
                            pmin = locus_data.get('Pathogenic min')
                            pmax = locus_data.get('Pathogenic max')
                            try:
                                pmin_val = float(pmin) if pmin is not None else None
                            except Exception:
                                pmin_val = None
                            try:
                                pmax_val = float(pmax) if pmax is not None else None
                            except Exception:
                                pmax_val = None

                            if pmin_val is not None:
                                if pmax_val is not None:
                                    if (repeat_count_used >= pmin_val) and (repeat_count_used <= pmax_val):
                                        is_pathogenic = 'Yes'
                                else:
                                    if repeat_count_used >= pmin_val:
                                        is_pathogenic = 'Yes'
                            else:
                                if pmax_val is not None and repeat_count_used <= pmax_val:
                                    is_pathogenic = 'Yes'

                        allele_records.append({
                            'Chromosome': chrom.replace('chr', ''),
                            'Position (Start)': pos,
                            'Gene': locus_data['Gene'],
                            'Disease': locus_data['Disease'],
                            'Disease ID': locus_data['Disease ID'],
                            'Motif': motif_val,
                            'Motif length': motif_len_used,
                            'Flank Motif Structure': locus_data['Flank Motif Structure'],
                            'Benign min': locus_data['Benign min'],
                            'Benign max': locus_data['Benign max'],
                            'Pathogenic min': locus_data['Pathogenic min'],
                            'Pathogenic max': locus_data['Pathogenic max'],
                            'Inheritance': inheritance_text,
                            'Allele Length (bp)': allele_len,        # NaN when invalid
                            'Repeat count (Calc)': repeat_count_used,     # NaN when invalid
                            'Is Pathogenic': is_pathogenic,
                            'Haplotype': hap_idx + 1,                # 1 or 2
                            'Sample ID (Raw)': sample_id_full,
                            'Sample ID (Cleaned)': sample_id_clean,
                            'SubPop': demo.get('SubPop', 'Unknown'),
                            'SuperPop': demo.get('SuperPop', 'Unknown'),
                            'Sex': demo.get('Sex', 'Unknown'),
                            'Pore': demo.get('Pore', 'Unknown'),
                            # --- QC flags ---
                            'QC Flag': 'Yes' if qc_reasons else '',
                            'QC Reasons': ';'.join(qc_reasons) if qc_reasons else '',
                        })

    if not allele_records:
        print("No allele records parsed.")
        return pd.DataFrame()

    df = pd.DataFrame(allele_records)
    print("Finished parsing VCF.")
    print(f"  Total variant records parsed: {total_sites}")
    print(f"  Generated allele rows: {len(df)}")
    print(f"  Unique samples represented: {df['Sample ID (Raw)'].nunique()}")
    return df

# ---------- Orchestration ----------
def create_integrated_spreadsheet(loci_json: str, vcf_file: str, sample_csv: str,
                                  kgp_txt: str, output_file: str, drop_samples: List[str]) -> Optional[str]:
    loci_map = load_loci_data(loci_json)
    sample_map_csv = load_all_sample_info(sample_csv)
    sample_map_1kgp = load_1kgp_sample_info(kgp_txt)

    df_integrated = parse_vcf_and_merge(vcf_file, loci_map, sample_map_csv, sample_map_1kgp, drop_samples)
    if df_integrated.empty:
        print("No data parsed; skipping Excel write.")
        return None

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
        'Is Pathogenic': 'Is Pathogenic',
        'Haplotype': 'Haplotype',
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
        'Pore': 'Pore',
        'QC Flag': 'QC Flag',
        'QC Reasons': 'QC Reasons',
    }

    present = [c for c in column_rename_map if c in df_integrated.columns]
    df_final = df_integrated[present].rename(columns=column_rename_map)

    desired_order = [
        'Chromosome', 'Position', 'Gene', 'Disease', 'Disease ID',
        'Motif', 'Motif length', 'Allele length', 'Repeat count',
        'Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max',
        'Is Pathogenic', 'Inheritance', 'Haplotype', 'Sample ID', 'Sample ID Cleaned',
        'SubPop', 'SuperPop', 'Sex', 'Flank Motif', 'Pore',
        'QC Flag', 'QC Reasons'
    ]
    df_final = df_final.reindex(columns=desired_order)

    # Write "NaN" strings to Excel for visibility (keep numeric NaN in memory)
    for col in ['Allele length', 'Repeat count']:
        if col in df_final.columns:
            df_final[col] = df_final[col].where(pd.notna(df_final[col]), 'NaN')

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
        kgp_txt=args.kgp_sample_info,
        output_file=args.output,
        drop_samples=args.drop_sample,
    )
    if out:
        print(f"Created spreadsheet: {out}")
    else:
        print("No spreadsheet created.")

if __name__ == "__main__":
    main()
