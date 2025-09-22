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
   - Fills all required columns, blank if missing, 'unknown' if unknown
   - Saves enriched dataset as CSV and Excel
   - Test mode: limit rows and optionally save to test_outputs/
---------------------------------------------
"""

import os
import json
import numpy as np
import pandas as pd

# --- TEST MODE ---
TEST_MODE = False
TEST_LIMIT = 100
SAVE_TEST_OUTPUTS = True

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots"

CSV_PATH  = f"{BASE_DIR}/Results/Debug/debug_info_83_loci_503_samples.csv"
JSON_PATH = f"{BASE_DIR}/Data/Other Data/STRchive-loci.json"
TSV_PATH  = f"{BASE_DIR}/Data/Other Data/sample_information.tsv"

OUTPUT_BASE     = f"{BASE_DIR}/Results/Matching Files Outputs"
OUTPUT_DIR_CSV  = os.path.join(OUTPUT_BASE, "CSVs")
OUTPUT_DIR_XLSX = os.path.join(OUTPUT_BASE, "Excels")
os.makedirs(OUTPUT_DIR_CSV, exist_ok=True)
os.makedirs(OUTPUT_DIR_XLSX, exist_ok=True)

if TEST_MODE:
    TEST_OUT = os.path.join(OUTPUT_BASE, "test_outputs")
    os.makedirs(TEST_OUT, exist_ok=True)
    OUTPUT_DIR_CSV  = TEST_OUT
    OUTPUT_DIR_XLSX = TEST_OUT

OUTPUT_CSV   = os.path.join(OUTPUT_DIR_CSV,  "83_loci_503_samples_withancestrycolumns_updated.csv")
OUTPUT_EXCEL = os.path.join(OUTPUT_DIR_XLSX, "83_loci_503_samples_withancestrycolumns_updated.xlsx")

# --- Load data ---
csv_data = pd.read_csv(CSV_PATH, sep=",", header=0)
with open(JSON_PATH, "r") as file:
    json_data = json.load(file)
tsv_data = pd.read_csv(TSV_PATH, sep="\t", skiprows=1)

work_df = csv_data.copy()
if TEST_MODE:
    work_df = work_df.head(TEST_LIMIT)

# --- Add ancestry population columns ---
def add_population_info(csv_df, tsv_df):
    id_col    = tsv_df.columns[0]
    sub_col   = tsv_df.columns[2]
    super_col = tsv_df.columns[3]
    sex_col   = tsv_df.columns[1] if len(tsv_df.columns) > 1 else None

    population_map = (
        tsv_df.set_index(id_col)[[sub_col, super_col] + ([sex_col] if sex_col else [])]
              .to_dict(orient='index')
    )

    csv_df['Sample ID'] = csv_df['Sample ID'].astype(str).str.strip()
    csv_df['SubPop'] = csv_df['Sample ID'].map(lambda x: population_map.get(x, {}).get(sub_col, 'unknown'))
    csv_df['SuperPop'] = csv_df['Sample ID'].map(lambda x: population_map.get(x, {}).get(super_col, 'unknown'))
    if sex_col:
        csv_df['Sex'] = csv_df['Sample ID'].map(lambda x: population_map.get(x, {}).get(sex_col, 'unknown'))
    else:
        csv_df['Sex'] = 'unknown'
    return csv_df

# --- Normalize population names ---
def clean_population_names(name):
    mapping = {
        'Unknown': 'unknown',
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
    return mapping.get(name, name if name else 'unknown')

def add_population_info_with_cleaned(csv_df):
    csv_df['SubPop'] = csv_df['SubPop'].apply(clean_population_names)
    csv_df['SuperPop'] = csv_df['SuperPop'].apply(clean_population_names)
    return csv_df

# --- Add inheritance info from JSON ---
def add_json_info_by_chrom_pos(csv_df, json_data):
    json_df = pd.DataFrame(json_data)
    json_df['chrom']      = json_df['chrom'].astype(str).str.replace('chr', '', regex=False)
    json_df['start_hg38'] = pd.to_numeric(json_df['start_hg38'], errors='coerce')
    json_df['stop_hg38']  = pd.to_numeric(json_df['stop_hg38'], errors='coerce')

    csv_df['Chromosome'] = csv_df['Chromosome'].astype(str).str.replace('chr', '', regex=False)
    csv_df['Position']   = pd.to_numeric(csv_df['Position'], errors='coerce')

    inheritance_list = []
    for _, row in csv_df.iterrows():
        chrom, pos = row['Chromosome'], row['Position']
        match = json_df[
            (json_df['chrom'] == chrom) &
            (json_df['start_hg38'] <= pos) &
            (json_df['stop_hg38'] >= pos)
        ]
        if not match.empty:
            inheritance_list.append(match.iloc[0].get('inheritance', 'unknown'))
        else:
            inheritance_list.append('unknown')
    csv_df['Inheritance'] = inheritance_list
    return csv_df

# --- Fill all required columns, blank if missing, 'unknown' if unknown ---
def fill_column(df, col, default=''):
    if col not in df.columns:
        df[col] = default
    df[col] = df[col].replace({np.nan: '', None: '', 'nan': '', 'NaN': '', pd.NA: ''})
    df[col] = df[col].replace({'Unknown': 'unknown', 'unknown': 'unknown'})
    return df

# --- Apply the enrichment ---
df = add_population_info(work_df, tsv_data)
df = add_population_info_with_cleaned(df)
df = add_json_info_by_chrom_pos(df, json_data)

# --- Derived fields ---
df['Motif length']  = df['Motif'].astype(str).str.len().replace(0, '')
df['Allele length'] = pd.to_numeric(df['Allele length'], errors='coerce')
df['Repeat count'] = df.apply(
    lambda row: row['Allele length'] / row['Motif length'] if row['Motif length'] not in ['', 0] else '',
    axis=1
)
df['Sample ID Cleaned'] = df['Sample ID'].astype(str).str.strip()

# --- Outlier columns ---
df = fill_column(df, 'Is Outlier', '')
df = fill_column(df, 'Is Extreme Outlier', '')

# --- Disease ID ---
df = fill_column(df, 'Disease ID', '')

# --- Flank Motif ---
df = fill_column(df, 'Flank Motif', '')

# --- Pathogenic/Benign columns ---
for col in ['Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max']:
    df[col] = pd.to_numeric(df.get(col, ''), errors='coerce').replace({np.nan: '', None: ''})

def flag_pathogenic(row):
    if row['Pathogenic min'] == '' or row['Pathogenic max'] == '' or row['Repeat count'] == '':
        return ''
    try:
        return row['Pathogenic min'] <= row['Repeat count'] <= row['Pathogenic max']
    except Exception:
        return ''
df['Is pathogenic'] = df.apply(flag_pathogenic, axis=1)

# --- Ensure all required columns exist and are filled ---
required_cols = [
    'Chromosome', 'Position', 'Gene', 'Disease', 'Disease ID', 'Motif', 'Motif length',
    'Allele length', 'Repeat count', 'Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max',
    'Inheritance', 'Sample ID', 'Sample ID Cleaned', 'SubPop', 'SuperPop', 'Sex',
    'Flank Motif', 'Is Outlier', 'Is Extreme Outlier'
]
for col in required_cols:
    df = fill_column(df, col, '')

# --- Reorder columns for readability ---
df = df[required_cols + [c for c in df.columns if c not in required_cols]]

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

print("--- Done ---")
