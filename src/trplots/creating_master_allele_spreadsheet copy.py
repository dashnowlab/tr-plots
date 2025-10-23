# NOTE TO SELF: need to update file paths to use config.py style paths

"""
---------------------------------------------
 Script: Creating Master Allele Spreadsheet
 Purpose:
    - Integrate VCF and locus metadata into a comprehensive allele spreadsheet
    - Include demographic and technology data from sample summary CSV
    - Save as Excel file with clear column structure
---------------------------------------------
"""

import pandas as pd
import json
import csv
import io
import re
import os
import gzip
from openpyxl.utils import get_column_letter
import argparse

# --- File Paths ---
# UPDATED: Use the specific VCF file from the latest context
VCF_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/sequencing_data/83_loci_503_samples/1000g-ont-strchive-83_loci_503_samples.vcf.gz'
LOCI_JSON = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive-loci.json'

# Data now comes from the summary CSV, including Pop/Sex/Pore
SAMPLE_SUMMARY_CSV = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/1kgp_ont_500_summary_-_sheet1.csv'
OUTPUT_FILE = '/Users/annelisethorn/Documents/github/tr-plots/data/other_data/master_allele_spreadsheet.xlsx'

def parse_args():
    p = argparse.ArgumentParser(description="Create integrated master allele spreadsheet")
    p.add_argument("--vcf", type=str, default=VCF_FILE, help="Input VCF (gzipped or plain)")
    p.add_argument("--loci-json", type=str, default=LOCI_JSON, help="Loci JSON metadata")
    p.add_argument("--sample-csv", type=str, default=SAMPLE_SUMMARY_CSV, help="Sample summary CSV (header on row 2)")
    p.add_argument("--output", type=str, default=OUTPUT_FILE, help="Output Excel file path")
    return p.parse_args()

# --- Data Loading and Parsing Functions ---

def load_loci_data(json_file):
    """
    Loads STR Loci information from JSON into a dict. 
    UPDATED to pull 'flank_motif' for the Flank Motif Structure column.
    """
    print("Loading Loci Data...")
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    loci_map = {}
    missing_flank_count = 0
    
    for entry in data:
        # --- START: Updated to pull 'flank_motif' as requested ---
        flank_motif_value = entry.get('flank_motif')
        
        # Track missing value and ensure it's a string if present but complex
        if flank_motif_value is None:
            missing_flank_count += 1
        elif not isinstance(flank_motif_value, str):
            # Convert complex structures (like lists/dicts) into a printable string
            flank_motif_value = str(flank_motif_value)
        # --- END: Updated to pull 'flank_motif' ---

        entry_data = {
            'Gene': entry.get('gene'),
            'Disease': entry.get('disease'),
            'Disease ID': entry.get('disease_id'),
            'Motif': entry.get('reference_motif_reference_orientation', [None])[0],
            'Motif length': entry.get('motif_len'),
            'Benign min': entry.get('benign_min'),
            'Benign max': entry.get('benign_max'),
            'Pathogenic min': entry.get('pathogenic_min'),
            'Pathogenic max': entry.get('pathogenic_max'),
            'Inheritance': ', '.join(entry.get('inheritance', [])),
            'Flank Motif Structure': flank_motif_value # Use the 'flank_motif' value
        }
        
        id_key = entry['id']
        pos_key = f"{entry['chrom'].replace('chr', '')}:{entry['start_hg38']}"
        
        loci_map[id_key] = entry_data
        loci_map[pos_key] = entry_data
        
    print(f"Loaded {len(data)} unique loci (with {len(loci_map)} total lookup keys).")
    if missing_flank_count > 0:
        print(f"DIAGNOSTIC: Found {missing_flank_count} loci in the JSON with missing 'Flank Motif Structure' (Flank Motif).")
    return loci_map

def load_all_sample_info(summary_csv_file):
    """
    Loads sample demographic and technology data directly from the summary CSV.
    The header for the required columns is on the second row (index 1).
    """
    print(f"Loading All Sample Data from {summary_csv_file}...")
    
    # Load the CSV, specifying the header is on the second row (index 1)
    df_summary = pd.read_csv(summary_csv_file, header=1)
    
    # Select and rename the required columns
    df_summary.rename(columns={
        'Sample_ID': 'Sample ID Cleaned',
        'Sex': 'Sex',
        'SubPop': 'SubPop',
        'SuperPop': 'SuperPop',
        'Pore': 'Pore'
    }, inplace=True)
    
    # Filter down to the essential columns and set index
    columns_to_keep = ['Sample ID Cleaned', 'SubPop', 'SuperPop', 'Sex', 'Pore']
    df_sample_info = df_summary[columns_to_keep].drop_duplicates(
        subset=['Sample ID Cleaned'], keep='first'
    ).set_index('Sample ID Cleaned')

    # Ensure no missing values (fill with 'Unknown')
    df_sample_info['Sex'] = df_sample_info['Sex'].fillna('Unknown')
    df_sample_info['SubPop'] = df_sample_info['SubPop'].fillna('Unknown')
    df_sample_info['SuperPop'] = df_sample_info['SuperPop'].fillna('Unknown')
    df_sample_info['Pore'] = df_sample_info['Pore'].fillna('Unknown')

    # Convert to dictionary
    sample_map = df_sample_info.to_dict('index')
    
    print(f"Loaded {len(sample_map)} unique sample demographic entries.")
    return sample_map

def parse_vcf_and_merge(vcf_file, loci_map, sample_map):
    """
    Parses the VCF file, flattens it to an allele table, and merges metadata.
    """
    print("Parsing VCF and building allele table...")
    allele_records = []
    vcf_header_sample_ids = []
    
    try:
        opener = gzip.open if vcf_file.endswith('.gz') else open
        mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
        DEFAULT_LOCUS_DATA = {
            'Gene': 'Unknown', 'Disease': 'Unknown', 'Disease ID': 'Unknown',
            'Motif': None, 'Motif length': 0, 
            'Benign min': None, 'Benign max': None, 
            'Pathogenic min': None, 'Pathogenic max': None, 
            'Inheritance': 'Unknown', 'Flank Motif Structure': None
        }

        with opener(vcf_file, mode, encoding='utf-8') as f:
            for line in f:
                if line.startswith('##'): continue
                
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    vcf_header_sample_ids = header[9:] 
                    print(f"VCF Header Sample Count: {len(vcf_header_sample_ids)}")
                    continue
                
                if not vcf_header_sample_ids: continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 10: continue
                    
                chrom, pos, identifier, ref, alt_str, qual, filter_val, info_str, fmt_str = fields[:9]
                
                info_dict = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info_dict[k] = v
                motif_from_info = info_dict.get('MOTIF')
                
                # Dual Lookup Strategy for Locus Data
                locus_key_id = identifier if identifier and identifier != '.' else None
                locus_key_pos = f"{chrom.replace('chr', '')}:{pos}"

                locus_data = loci_map.get(locus_key_id)
                if locus_data is None:
                    locus_data = loci_map.get(locus_key_pos)
                    
                if locus_data is None:
                    locus_data = DEFAULT_LOCUS_DATA.copy()
                    locus_data['Motif'] = motif_from_info 
                
                # Motif Length Sanitization
                motif_len = locus_data.get('Motif length')
                if motif_len is None or not isinstance(motif_len, (int, float)) or motif_len <= 0:
                    try:
                        motif_len = float(motif_len) if motif_len else 0
                    except ValueError:
                        motif_len = 0
                    if motif_len <= 0:
                        motif_to_use = locus_data.get('Motif') or motif_from_info
                        if motif_to_use and len(motif_to_use) > 0:
                            motif_len = len(motif_to_use)
                        else:
                            motif_len = 3 

                # Process Genotypes per Sample
                genotype_data = fields[9:]
                
                for i, sample_genotype_str in enumerate(genotype_data):
                    sample_id_full = vcf_header_sample_ids[i]
                    
                    match = re.match(r'^(HG\d+|GM\d+|NA\d+)', sample_id_full)
                    sample_id_clean = match.group(0) if match else sample_id_full.split('-')[0].split('_')[0]

                    # Fetch demographic data from the consolidated map
                    demo_data = sample_map.get(sample_id_clean, {
                        'SubPop': 'Unknown', 'SuperPop': 'Unknown', 'Sex': 'Unknown', 'Pore': 'Unknown'
                    })

                    format_keys = fmt_str.split(':')
                    format_values = sample_genotype_str.split(':')
                    fmt_data = dict(zip(format_keys, format_values))
                    
                    # Split the GT field (e.g., '0/1' -> ['0', '1'], './.' -> ['.', '.'])
                    gt_field = fmt_data.get('GT', './.')
                    gt_indices = re.split(r'[/|]', gt_field)
                    
                    allele_lengths_str = fmt_data.get('AL', '.')
                    
                    # 1. Build length map: Index 0 is always Ref length
                    len_map = {0: len(ref)}
                    try:
                        allele_lengths = [int(l) for l in allele_lengths_str.split(',') if l and l != '.']
                        for k, length in enumerate(allele_lengths):
                            len_map[k + 1] = length
                    except ValueError: pass
                    
                    # 2. Iterate over the GT indices to extract each allele's length
                    for allele_type in gt_indices:
                        
                        # CRITICAL FILTERING STEP: If the individual allele index is missing ('.'), skip
                        if allele_type == '.': continue
                        
                        # Fixes: Ensure individual index is an integer
                        try:
                            allele_type_int = int(allele_type)
                        except ValueError: 
                            # If GT is malformed, skip the allele
                            continue 
                        
                        allele_length = None
                        
                        # FIX: Logic to prioritize reference length (index 0) if length map misses it
                        if allele_type_int == 0:
                            # Use the reference sequence length directly, ensuring it's not empty
                            allele_length = len(ref) if ref else None
                        else:
                            # For alternate alleles (index 1+), fetch from AL field map
                            allele_length = len_map.get(allele_type_int)

                        # IMPORTANT: Skip if allele length still couldn't be determined
                        if allele_length is None: continue
                            
                        repeat_count = round(allele_length / motif_len, 2) if motif_len and motif_len > 0 else 'N/A'
                        
                        record = {
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
                            'Allele Length (bp)': allele_length, 
                            'Repeat count (Calc)': repeat_count,
                            'Sample ID (Raw)': sample_id_full,
                            'Sample ID (Cleaned)': sample_id_clean,
                            'SubPop': demo_data.get('SubPop', 'Unknown'),
                            'SuperPop': demo_data.get('SuperPop', 'Unknown'), 
                            'Sex': demo_data.get('Sex', 'Unknown'),
                            'Pore': demo_data.get('Pore', 'Unknown'),
                        }
                        allele_records.append(record)
    
    except FileNotFoundError:
        print(f"Error: The VCF file was not found at the specified path: {vcf_file}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An unexpected error occurred during VCF parsing: {e}")
        return pd.DataFrame()

    if allele_records:
        df_alleles = pd.DataFrame(allele_records)
        processed_sample_ids = df_alleles['Sample ID (Raw)'].unique()
        
        print(f"Finished parsing VCF. Generated {len(allele_records)} individual allele calls.")
        print(f"Total Unique Samples in final spreadsheet: {len(processed_sample_ids)}")
    else:
        df_alleles = pd.DataFrame()
        
    return df_alleles

def create_integrated_spreadsheet(loci_json, vcf_file, output_file):
    """Main function to orchestrate data loading, parsing, and saving."""
    try:
        # 1. Load Data
        loci_map = load_loci_data(loci_json)
        
        # 2. Load Sample Demographics
        sample_map = load_all_sample_info(SAMPLE_SUMMARY_CSV)
        
        # 3. Parse VCF and Merge Data
        df_integrated = parse_vcf_and_merge(vcf_file, loci_map, sample_map)

        if df_integrated.empty:
            print("No data was processed from the VCF file. Skipping Excel creation.")
            return

        # --- Final Formatting and Saving ---
        column_rename_map = {
            'Chromosome': 'Chromosome', 'Position (Start)': 'Position', 'Gene': 'Gene', 
            'Disease': 'Disease', 'Disease ID': 'Disease ID', 'Motif': 'Motif', 
            'Motif length': 'Motif length', 'Allele Length (bp)': 'Allele length', 
            'Repeat count (Calc)': 'Repeat count', 'Benign min': 'Benign min', 
            'Benign max': 'Benign max', 'Pathogenic min': 'Pathogenic min', 
            'Pathogenic max': 'Pathogenic max', 'Inheritance': 'Inheritance',
            'Sample ID (Raw)': 'Sample ID', 'Sample ID (Cleaned)': 'Sample ID Cleaned', 
            'SubPop': 'SubPop', 'SuperPop': 'SuperPop', 'Sex': 'Sex',
            'Flank Motif Structure': 'Flank Motif',
            'Pore': 'Pore'
        }
        
        columns_to_select = [col for col in column_rename_map.keys() if col in df_integrated.columns]
        
        df_final = df_integrated[columns_to_select].rename(columns=column_rename_map)
        desired_order = list(column_rename_map.values())
        df_final = df_final[desired_order]
        
        # Save to Excel
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df_final.to_excel(writer, index=False, sheet_name='Integrated Alleles')
        
        print(f"\nSuccess! Integrated spreadsheet saved to '{output_file}'.")
        return output_file

    except FileNotFoundError as e:
        print(f"Error: One of the input files was not found: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# --- Execution ---
def main(args=None):
    if args is None:
        args = parse_args()

    # Allow CLI overrides of the module-level SAMPLE_SUMMARY_CSV
    global SAMPLE_SUMMARY_CSV
    SAMPLE_SUMMARY_CSV = args.sample_csv

    result = create_integrated_spreadsheet(args.loci_json, args.vcf, args.output)
    if result:
        print(f"Created spreadsheet: {result}")
    else:
        print("No spreadsheet created.")


if __name__ == "__main__":
    main()
