#!/usr/bin/env python3
"""
Master Allele Spreadsheet Creator
(with Pure Motif Repeat columns + pathogenic flags + robust Excel writer)

Purpose:
  Parse a VCF of STR loci (gzipped or plain), merge locus metadata and sample
  demographics (CSV + 1kGP), and produce an integrated per-haplotype Excel
  spreadsheet with:
    - Allele lengths and length-based repeat counts
    - Pure-motif metrics (pure exact count, max pure run, is-pure flag)
    - Pathogenic flags based on locus pathogenic ranges
      * IsPathogenic_ByRepeatCount
      * IsPathogenic_ByMaxPureRun
      * IsPathogenic_Agreement (Yes if the two match; No if both known but differ; Unknown otherwise)
    - QC flags and reasons

Robust Excel export:
  - Tries xlsxwriter first (fast), falls back to openpyxl
  - Chunking (default 50k rows per sheet)
  - Optional --nan-as-string to render NaN as "NaN" in Excel

Assumptions:
  - For pathogenic calls, repeat-based thresholds use:
      Pathogenic min (inclusive) and/or Pathogenic max (inclusive) if present.
    If only min is provided -> >= min is Yes.
    If only max is provided -> <= max is Yes.
    If both missing -> Unknown.
"""

import argparse
import gzip
import json
import os
import re
import shutil
import tempfile
from typing import Dict, List, Optional, Tuple

import pandas as pd
from openpyxl.utils import get_column_letter  # noqa: F401

# ---------- Defaults (override via CLI) ----------
VCF_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/sequencing_data/81_loci_502_samples/1000g-ont-strchive-81_loci_502_samples_81224_alleles.vcf.gz'
LOCI_JSON = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive-loci.json'
SAMPLE_SUMMARY_CSV = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/1kgp_ont_500_summary_-_sheet1.csv'
KGP_SAMPLE_INFO = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/1000_genomes_20130606_sample_info.txt'
OUTPUT_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/master_allele_spreadsheet_w_motif_info_v2.xlsx'

# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description="Create integrated master allele spreadsheet (pure motif + pathogenic flags)")
    p.add_argument("--vcf", type=str, default=VCF_FILE, help="Input VCF (gzipped or plain)")
    p.add_argument("--loci-json", type=str, default=LOCI_JSON, help="Loci JSON metadata")
    p.add_argument("--sample-csv", type=str, default=SAMPLE_SUMMARY_CSV, help="Sample summary CSV (header on row 2)")
    p.add_argument("--kgp-sample-info", type=str, default=KGP_SAMPLE_INFO,
                   help="1000 Genomes sample-info TXT (20130606)")
    p.add_argument("--output", type=str, default=OUTPUT_FILE, help="Output Excel file path")
    p.add_argument("--drop-sample", action="append", default=[],
                   help="Sample ID(s) to drop (exact match to VCF header sample names). Repeatable.")
    # Excel export options
    p.add_argument("--nan-as-string", action="store_true",
                   help="Render NaNs as the literal string 'NaN' in Excel (slower, larger files).")
    p.add_argument("--rows-per-sheet", type=int, default=50000,
                   help="Max rows per sheet for chunking large DataFrames.")
    p.add_argument("--excel-engine", type=str, default="auto",
                   choices=["auto", "xlsxwriter", "openpyxl"],
                   help="Preferred Excel writer engine. 'auto' tries xlsxwriter then openpyxl.")
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
    'YRI': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR', 'ACB': 'AFR', 'ASW': 'AFR',
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

    df_out = df[['Sample', 'SubPop_from_1kGP', 'SuperPop_from_1kGP', 'Sex_from_1kGP']]\
        .drop_duplicates('Sample').set_index('Sample')
    result = df_out.to_dict('index')
    print(f"1kGP entries: {len(result)}")
    return result

# ---------- Helpers ----------
_GT_SPLIT_RE = re.compile(r'[\/|]')
_INHER_SPLIT_RE = re.compile(r'[,/;| ]+')

def _has_both_ad_ar(inheritance_text: Optional[str]) -> bool:
    if not inheritance_text:
        return False
    toks = [t.strip().lower() for t in _INHER_SPLIT_RE.split(str(inheritance_text)) if t.strip()]
    return ('ad' in toks) and ('ar' in toks)

def _calc_motif_len(locus_data: dict, motif_from_info: Optional[str], ref: str) -> float:
    """
    Determine motif length (bp). If locus metadata gives motif_len, use it.
    Else if motif seq known (from locus or INFO), use its length.
    Else fall back to REF length (so length-based repeat count still computes).
    """
    motif_len = locus_data.get('Motif length')
    if isinstance(motif_len, (int, float)) and motif_len > 0:
        return float(motif_len)
    motif_to_use = locus_data.get('Motif') or motif_from_info
    if motif_to_use:
        return float(len(motif_to_use))
    return float(len(ref)) if ref else 0.0

def _get_motif_seq(locus_data: dict, motif_from_info: Optional[str]) -> Optional[str]:
    """Return the motif sequence string if known, else None."""
    motif_to_use = locus_data.get('Motif') or motif_from_info
    if motif_to_use and isinstance(motif_to_use, str) and len(motif_to_use) > 0:
        return motif_to_use
    return None

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

def _full_pure_count(seq: str, motif: str) -> Optional[int]:
    """Return n if seq == motif * n; else None."""
    if not seq or not motif:
        return None
    m = len(motif)
    if m == 0 or (len(seq) % m) != 0:
        return None
    n = len(seq) // m
    return n if (motif * n) == seq else None

def _max_consecutive_motif_run(seq: str, motif: str) -> Optional[int]:
    """Return the maximum uninterrupted motif repeat run within seq; None if motif unknown."""
    if not seq or not motif:
        return None
    m = len(motif)
    if m == 0:
        return None
    best = 0
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

def _pathogenic_call(value: Optional[float],
                     pathogenic_min: Optional[float],
                     pathogenic_max: Optional[float]) -> str:
    """
    Ternary pathogenic call for a numeric value against locus pathogenic bounds.
    - If both min and max are None -> 'Unknown'
    - If value is NaN/None -> 'Unknown'
    - If only min present -> 'Yes' if value >= min else 'No'
    - If only max present -> 'Yes' if value <= max else 'No'
    - If both present -> 'Yes' if min <= value <= max else 'No'
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "Unknown"

    has_min = pathogenic_min is not None and str(pathogenic_min) != 'nan'
    has_max = pathogenic_max is not None and str(pathogenic_max) != 'nan'

    if not has_min and not has_max:
        return "Unknown"

    if has_min and not has_max:
        try:
            return "Yes" if v >= float(pathogenic_min) else "No"
        except Exception:
            return "Unknown"

    if has_max and not has_min:
        try:
            return "Yes" if v <= float(pathogenic_max) else "No"
        except Exception:
            return "Unknown"

    try:
        return "Yes" if float(pathogenic_min) <= v <= float(pathogenic_max) else "No"
    except Exception:
        return "Unknown"

def _agreement3(a: Optional[str], b: Optional[str]) -> str:
    """
    Return 'Yes' if a and b are equal and both not Unknown.
    Return 'No'  if both known and differ.
    Else 'Unknown'.
    """
    a = (a or "Unknown").strip()
    b = (b or "Unknown").strip()
    if a == "Unknown" or b == "Unknown":
        return "Unknown"
    return "Yes" if a == b else "No"

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
            # If unknown locus but INFO carries MOTIF, record it so we can still compute purity
            if locus_data is DEFAULT_LOCUS_DATA and motif_from_info and not locus_data['Motif']:
                locus_data = dict(DEFAULT_LOCUS_DATA)
                locus_data['Motif'] = motif_from_info

            motif_len = _calc_motif_len(locus_data, motif_from_info, ref)
            motif_seq = _get_motif_seq(locus_data, motif_from_info)  # may be None
            fmt_keys = fmt_str.split(':')

            # Build allele index -> sequence map for this site
            alt_list = alt_str.split(',') if alt_str else []
            allele_index_to_seq = {0: ref}
            for idx, alt_seq in enumerate(alt_list, start=1):
                allele_index_to_seq[idx] = alt_seq

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

                if _is_unknown(demo.get('SubPop')) and kgp.get('SubPop_from_1kGP'):
                    demo['SubPop'] = kgp['SubPop_from_1kGP']
                if _is_unknown(demo.get('SuperPop')):
                    sp = kgp.get('SuperPop_from_1kGP')
                    if sp:
                        demo['SuperPop'] = sp
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

                # Per-haplotype allele lengths from AL (bp)
                al_haps = _hap_lengths_from_AL(fmt_data.get('AL', '.'))

                # --- inheritance ambiguity check (once per row) ---
                inheritance_text = locus_data.get('Inheritance', '')
                inheritance_ambig = _has_both_ad_ar(inheritance_text)

                # Iterate haplotypes
                for hap_idx, allele_tok in enumerate(gt_toks):
                    if allele_tok == '.':
                        continue

                    qc_reasons = []

                    # Determine allele sequence by GT index (0=REF, k=ALTk)
                    allele_seq = None
                    if allele_tok.isdigit():
                        aidx = int(allele_tok)
                        allele_seq = allele_index_to_seq.get(aidx)
                    else:
                        aidx = None  # unexpected

                    # Allele length per haplotype with QC
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

                    # Repeat count (length-based)
                    if motif_len and isinstance(motif_len, (int, float)) and motif_len > 0:
                        repeat_count = (round(allele_len / motif_len, 2)
                                        if (isinstance(allele_len, (int, float)) and allele_len > 0)
                                        else float('nan'))
                    else:
                        qc_reasons.append('motif_len_le0_or_missing')
                        repeat_count = float('nan')

                    # --- Pure motif metrics (sequence-based) ---
                    pure_count_exact: Optional[int] = None
                    max_pure_run: Optional[int] = None
                    is_pure_flag = 'No'

                    if motif_seq:
                        if allele_seq:
                            pure_count_exact = _full_pure_count(allele_seq, motif_seq)
                            max_pure_run = _max_consecutive_motif_run(allele_seq, motif_seq)
                            if pure_count_exact is not None and pure_count_exact > 0:
                                is_pure_flag = 'Yes'
                        else:
                            qc_reasons.append('allele_seq_missing_for_GT_index')
                            pure_count_exact = None
                            max_pure_run = None
                    else:
                        qc_reasons.append('motif_seq_missing')

                    # --- Pathogenic calls (repeat-based ranges compare to REPEAT COUNTS) ---
                    pmin = locus_data.get('Pathogenic min')
                    pmax = locus_data.get('Pathogenic max')

                    call_by_repeat = _pathogenic_call(repeat_count, pmin, pmax)
                    call_by_maxrun = _pathogenic_call(max_pure_run, pmin, pmax)
                    agree_flag = _agreement3(call_by_repeat, call_by_maxrun)

                    if inheritance_ambig:
                        qc_reasons.append('inheritance_AD_and_AR')

                    allele_records.append({
                        'Chromosome': chrom.replace('chr', ''),
                        'Position (Start)': pos,
                        'Gene': locus_data['Gene'],
                        'Disease': locus_data['Disease'],
                        'Disease ID': locus_data['Disease ID'],
                        'Motif': motif_seq if motif_seq else locus_data.get('Motif'),
                        'Motif length': motif_len,
                        'Flank Motif Structure': locus_data['Flank Motif Structure'],
                        'Benign min': locus_data['Benign min'],
                        'Benign max': locus_data['Benign max'],
                        'Pathogenic min': pmin,
                        'Pathogenic max': pmax,
                        'Inheritance': inheritance_text,

                        'Allele Length (bp)': allele_len,             # NaN when invalid
                        'Repeat count (Calc)': repeat_count,           # NaN when invalid

                        # Pure motif repeat columns
                        'Pure repeat count (Exact)': (float(pure_count_exact) if (pure_count_exact is not None) else float('nan')),
                        'Max pure motif run': (float(max_pure_run) if (max_pure_run is not None) else float('nan')),
                        'Is Pure Motif': is_pure_flag,

                        # Pathogenic flags
                        'IsPathogenic_ByRepeatCount': call_by_repeat,
                        'IsPathogenic_ByMaxPureRun': call_by_maxrun,
                        'IsPathogenic_Agreement': agree_flag,

                        'Haplotype': hap_idx + 1,                     # 1 or 2
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

# ---------- Robust Excel writer ----------
def _write_excel_chunked(df_integrated: pd.DataFrame,
                         df_summary: Optional[pd.DataFrame],
                         out_path: str,
                         engine_preference: str = "auto",
                         nan_as_string: bool = False,
                         rows_per_sheet: int = 50000) -> None:
    """
    Robust Excel writer:
      - tries xlsxwriter first (fast), falls back to openpyxl
      - splits big tables across multiple sheets
      - optional 'NaN' text rendering without mutating dtypes in memory
    """
    def _write_with_engine(engine_name: str, tmp_path: str):
        with pd.ExcelWriter(tmp_path, engine=engine_name) as writer:
            # Integrated Alleles (chunked)
            n = len(df_integrated)
            if n == 0:
                df_integrated.head(0).to_excel(writer, index=False, sheet_name="Integrated Alleles")
            else:
                start = 0
                part = 1
                while start < n:
                    end = min(start + rows_per_sheet, n)
                    chunk = df_integrated.iloc[start:end]
                    sheet = "Integrated Alleles" if part == 1 else f"Integrated Alleles ({part})"
                    na_rep = "NaN" if nan_as_string else None
                    chunk.to_excel(writer, index=False, sheet_name=sheet, na_rep=na_rep)
                    start = end
                    part += 1

            # Optional Locus Summary sheet
            if df_summary is not None and not df_summary.empty:
                na_rep = "NaN" if nan_as_string else None
                df_summary.to_excel(writer, index=False, sheet_name="Locus Summary", na_rep=na_rep)

    if engine_preference == "xlsxwriter":
        engines = ["xlsxwriter"]
    elif engine_preference == "openpyxl":
        engines = ["openpyxl"]
    else:
        engines = ["xlsxwriter", "openpyxl"]

    last_err = None
    with tempfile.TemporaryDirectory() as td:
        tmp_path = os.path.join(td, "tmp.xlsx")
        for eng in engines:
            try:
                _write_with_engine(eng, tmp_path)
                shutil.move(tmp_path, out_path)
                print(f"✅ Excel saved with engine='{eng}': {out_path}")
                return
            except Exception as e:
                print(f"Writer engine '{eng}' failed: {e}")
                last_err = e
    raise last_err or RuntimeError("Excel writer failed with all engines")

# ---------- Orchestration ----------
def create_integrated_spreadsheet(loci_json: str, vcf_file: str, sample_csv: str,
                                  kgp_txt: str, output_file: str, drop_samples: List[str],
                                  excel_engine: str, nan_as_string: bool, rows_per_sheet: int) -> Optional[str]:
    loci_map = load_loci_data(loci_json)
    sample_map_csv = load_all_sample_info(sample_csv)
    sample_map_1kgp = load_1kgp_sample_info(kgp_txt)

    df_integrated = parse_vcf_and_merge(vcf_file, loci_map, sample_map_csv, sample_map_1kgp, drop_samples)
    if df_integrated.empty:
        print("No data parsed; skipping Excel write.")
        return None

    # Rename/order for final Integrated sheet
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
        'Pure repeat count (Exact)': 'Pure repeat count (Exact)',
        'Max pure motif run': 'Max pure motif run',
        'Is Pure Motif': 'Is Pure Motif',
        'IsPathogenic_ByRepeatCount': 'IsPathogenic_ByRepeatCount',
        'IsPathogenic_ByMaxPureRun': 'IsPathogenic_ByMaxPureRun',
        'IsPathogenic_Agreement': 'IsPathogenic_Agreement',
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
        'Motif', 'Motif length',
        'Allele length', 'Repeat count',
        'Pure repeat count (Exact)', 'Max pure motif run', 'Is Pure Motif',
        'IsPathogenic_ByRepeatCount', 'IsPathogenic_ByMaxPureRun', 'IsPathogenic_Agreement',
        'Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max',
        'Inheritance', 'Haplotype', 'Sample ID', 'Sample ID Cleaned',
        'SubPop', 'SuperPop', 'Sex', 'Flank Motif', 'Pore',
        'QC Flag', 'QC Reasons'
    ]
    df_final = df_final.reindex(columns=desired_order)

    # Optional: quick Locus Summary
    try:
        df_summary = (
            df_final.groupby(['Chromosome', 'Position', 'Gene', 'Motif', 'Motif length'], dropna=False)
            .agg(
                n_alleles=('Allele length', 'size'),
                n_pure=('Is Pure Motif', lambda s: (s == 'Yes').sum()),
                pct_pure=('Is Pure Motif', lambda s: round((s == 'Yes').mean() * 100, 2)),
                median_repeat=('Repeat count', lambda s: pd.to_numeric(s, errors='coerce').median()),
                n_path_by_repeat=('IsPathogenic_ByRepeatCount', lambda s: (s == 'Yes').sum()),
                n_path_by_maxrun=('IsPathogenic_ByMaxPureRun', lambda s: (s == 'Yes').sum()),
                n_agree=('IsPathogenic_Agreement', lambda s: (s == 'Yes').sum()),
                n_disagree=('IsPathogenic_Agreement', lambda s: (s == 'No').sum()),
            )
            .reset_index()
            .sort_values(['Chromosome', 'Position'])
        )
    except Exception:
        df_summary = None

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    _write_excel_chunked(
        df_integrated=df_final,
        df_summary=df_summary,
        out_path=output_file,
        engine_preference=excel_engine,
        nan_as_string=nan_as_string,
        rows_per_sheet=rows_per_sheet
    )

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
        excel_engine=args.excel_engine,
        nan_as_string=args.nan_as_string,
        rows_per_sheet=args.rows_per_sheet,
    )
    if out:
        print(f"Created spreadsheet: {out}")
    else:
        print("No spreadsheet created.")

if __name__ == "__main__":
    main()
