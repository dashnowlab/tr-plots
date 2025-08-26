"""
---------------------------------------------
 Script: CSV + JSON + Ancestry Merger
 Purpose:
   - Reads debug CSV with allele + locus info
   - Merges ancestry information from TSV
   - Cleans and standardizes population labels
   - Enriches rows with inheritance info from JSON
   - Computes motif length & repeat counts
   - Flags pathogenic individuals based on thresholds
   - Saves enriched dataset as CSV and Excel
   - Test mode: limit rows and optionally save to test_outputs/
---------------------------------------------
"""

import os
import json
import numpy as np
import pandas as pd

# --- TEST MODE ---
TEST_MODE = False          # Toggle for quick checks
TEST_LIMIT = 100        # Number of rows to process in test mode
SAVE_TEST_OUTPUTS = True  # If True, write files in test mode; otherwise just preview

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots"

CSV_PATH  = f"{BASE_DIR}/Results/Debug/debug_info_83_loci_503_samples.csv"
JSON_PATH = f"{BASE_DIR}/Data/Other Data/STRchive-loci.json"
TSV_PATH  = f"{BASE_DIR}/Data/Other Data/sample_information.tsv"

# Output roots (normal mode)
OUTPUT_BASE     = f"{BASE_DIR}/Results/Matching Files Outputs"
OUTPUT_DIR_CSV  = os.path.join(OUTPUT_BASE, "CSVs")
OUTPUT_DIR_XLSX = os.path.join(OUTPUT_BASE, "Excels")
os.makedirs(OUTPUT_DIR_CSV, exist_ok=True)
os.makedirs(OUTPUT_DIR_XLSX, exist_ok=True)

# In test mode, write everything to a single test_outputs folder
if TEST_MODE:
    TEST_OUT = os.path.join(OUTPUT_BASE, "test_outputs")
    os.makedirs(TEST_OUT, exist_ok=True)
    OUTPUT_DIR_CSV  = TEST_OUT
    OUTPUT_DIR_XLSX = TEST_OUT

OUTPUT_CSV   = os.path.join(OUTPUT_DIR_CSV,  "83_loci_503_samples_withancestrycolumns.csv")
OUTPUT_EXCEL = os.path.join(OUTPUT_DIR_XLSX, "83_loci_503_samples_withancestrycolumns.xlsx")

# --- Load data ---
csv_data = pd.read_csv(CSV_PATH, sep=",", header=0)
with open(JSON_PATH, "r") as file:
    json_data = json.load(file)
tsv_data = pd.read_csv(TSV_PATH, sep="\t", skiprows=1)

# If test mode, work on a head() slice to speed things up
work_df = csv_data.copy()
if TEST_MODE:
    work_df = work_df.head(TEST_LIMIT)

# --- Add ancestry population columns ---
def add_population_info(csv_df, tsv_df):
    id_col    = tsv_df.columns[0]
    sub_col   = tsv_df.columns[2]
    super_col = tsv_df.columns[3]

    population_map = (
        tsv_df.set_index(id_col)[[sub_col, super_col]]
              .to_dict(orient='index')
    )

    csv_df['Sample ID'] = csv_df['Sample ID'].astype(str).str.strip()
    csv_df['Population'] = csv_df['Sample ID'].map(lambda x: population_map.get(x, {}).get(sub_col, 'Unknown'))
    csv_df['Population description'] = csv_df['Sample ID'].map(lambda x: population_map.get(x, {}).get(super_col, 'Unknown'))
    return csv_df

# --- Normalize population names ---
def clean_population_names(name):
    mapping = {
        'Unknown': 'Unknown',
        'Finnish in Finland': 'Finnish',
        'Han Chinese South, China': 'East Asian',
        'Puerto Rican in Puerto Rico': 'Admixed American',
        'Colombian in Medellin, Colombia': 'Admixed American',
        'African Caribbean in Barbados': 'African/African American',
        'Peruvian in Lima, Peru': 'Admixed American',
        'Kinh in Ho Chi Minh City, Vietnam': 'East Asian',
        'Gambian in Western Division Ð Mandinka': 'African/African American',
        'Punjabi in Lahore, Pakistan': 'South Asian',
        'Esan in Nigeria': 'African/African American',
        'Mende in Sierra Leone': 'African/African American',
        'Sri Lankan Tamil in the UK': 'South Asian',
        'Bengali in Bangladesh': 'South Asian',
        'Han Chinese in Beijing, China': 'East Asian',
        'Japanese in Tokyo, Japan': 'East Asian',
        'Luhya in Webuye, Kenya': 'African/African American',
        'Toscani in Italia': 'European (non Finnish)'
    }
    return mapping.get(name, name)

def add_population_info_with_cleaned(csv_df):
    csv_df['Cleaned population description'] = csv_df['Population description'].apply(clean_population_names)
    return csv_df

# --- Add inheritance info from JSON ---
def add_json_info_by_chrom_pos(csv_df, json_data):
    json_df = pd.DataFrame(json_data)
    json_df['chrom']      = json_df['chrom'].astype(str).str.replace('chr', '', regex=False)
    json_df['start_hg38'] = pd.to_numeric(json_df['start_hg38'], errors='coerce')
    json_df['stop_hg38']  = pd.to_numeric(json_df['stop_hg38'], errors='coerce')

    csv_df['Chromosome'] = csv_df['Chromosome'].astype(str).str.replace('chr', '', regex=False)
    csv_df['Position']   = pd.to_numeric(csv_df['Position'], errors='coerce')

    results = []
    for _, row in csv_df.iterrows():
        chrom, pos = row['Chromosome'], row['Position']
        match = json_df[
            (json_df['chrom'] == chrom) &
            (json_df['start_hg38'] <= pos) &
            (json_df['stop_hg38'] >= pos)
        ]
        if not match.empty:
            row['Inheritance'] = match.iloc[0]['inheritance']
        results.append(row)

    return pd.DataFrame(results)

# --- Apply the enrichment ---
df = add_population_info(work_df, tsv_data)
df = add_population_info_with_cleaned(df)
df = add_json_info_by_chrom_pos(df, json_data)

# Drop extra column if present
df.drop(columns=['ALLR'], errors='ignore', inplace=True)

# --- Derived fields ---
df['Motif length']  = df['Motif'].astype(str).str.len()
df['Allele length'] = pd.to_numeric(df['Allele length'], errors='coerce')

df['Repeat count'] = df.apply(
    lambda row: row['Allele length'] / row['Motif length'] if row['Motif length'] > 0 else None,
    axis=1
)

def flag_pathogenic(row):
    if pd.isna(row['Pathogenic min']) or pd.isna(row['Pathogenic max']) or pd.isna(row['Repeat count']):
        return np.nan
    return row['Pathogenic min'] <= row['Repeat count'] <= row['Pathogenic max']

df['Is pathogenic']   = df.apply(flag_pathogenic, axis=1)
df['Pathogenic min']  = pd.to_numeric(df['Pathogenic min'], errors='coerce')
df['Pathogenic max']  = pd.to_numeric(df['Pathogenic max'], errors='coerce')

# --- Reorder columns for readability ---
desired_order = [
    'Chromosome', 'Position', 'Gene', 'Disease', 'Motif',
    'Allele length', 'Motif length', 'Repeat count',
    'Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max',
    'Is pathogenic', 'Sample ID'
]
present_cols = [c for c in desired_order if c in df.columns]
remaining    = [c for c in df.columns if c not in present_cols]
df = df[present_cols + remaining]

# --- Save outputs ---
if TEST_MODE:
    print(f"[TEST] Rows processed: {len(df)} (limit={TEST_LIMIT})")
    if SAVE_TEST_OUTPUTS:
        df.to_csv(OUTPUT_CSV, index=False)
        df.to_excel(OUTPUT_EXCEL, index=False)
        print(f"[TEST] Saved CSV  → {OUTPUT_CSV}")
        print(f"[TEST] Saved XLSX → {OUTPUT_EXCEL}")
else:
    df.to_csv(OUTPUT_CSV, index=False)
    df.to_excel(OUTPUT_EXCEL, index=False)
    print(f"Saved CSV  → {OUTPUT_CSV}")
    print(f"Saved XLSX → {OUTPUT_EXCEL}")

# --- Finished ---
print("--- Done ---")
