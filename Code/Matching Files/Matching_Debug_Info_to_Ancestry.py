import pandas as pd
import json
import numpy as np

# File paths
CSV_PATH = "/Users/annelisethorn/Documents/Anschutz/Debug/debug_info_68_loci_100_samples.csv"
# CSV_PATH = "/Users/annelisethorn/Documents/Anschutz/Debug/debug_info_83_loci_88_samples.csv"
# CSV_PATH = "/Users/annelisethorn/Documents/Anschutz/Debug/debug_info_83_loci_503_samples.csv"

JSON_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json"
TSV_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/sample_information.tsv"

# Load CSV (VCF data) and use the first row as column headers
csv_data = pd.read_csv(CSV_PATH, sep=",", header=0)

# Make a copy so we don't modify the original data
csv_data_copy = csv_data.copy()

# Load JSON metadata for loci
with open(JSON_PATH, "r") as file:
    json_data = json.load(file)

# Load Sample Population Information (ancestry data)
tsv_data = pd.read_csv(TSV_PATH, sep="\t", skiprows=1)

# Add ancestry population columns
def add_population_info(csv_df_copy, tsv_df):
    tsv_first_col = tsv_df.columns[0]
    tsv_third_col = tsv_df.columns[2]
    tsv_fourth_col = tsv_df.columns[3]

    population_map = tsv_df.set_index(tsv_first_col)[[tsv_third_col, tsv_fourth_col]].to_dict(orient='index')

    csv_df_copy['Sample ID'] = csv_df_copy['Sample ID'].astype(str).str.strip()
    csv_df_copy['Population'] = csv_df_copy['Sample ID'].map(lambda x: population_map.get(x, {}).get(tsv_third_col, 'Unknown'))
    csv_df_copy['Population description'] = csv_df_copy['Sample ID'].map(lambda x: population_map.get(x, {}).get(tsv_fourth_col, 'Unknown'))

    return csv_df_copy

# Normalize population descriptions
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

def add_population_info_with_cleaned(csv_df_copy, tsv_df):
    csv_df_copy['Cleaned population description'] = csv_df_copy['Population description'].apply(clean_population_names)
    return csv_df_copy

# Add inheritance info from JSON using CHROM and POS
def add_json_info_by_chrom_pos(csv_df_copy, json_data):
    json_df = pd.DataFrame(json_data)
    json_df['chrom'] = json_df['chrom'].astype(str).str.replace('chr', '', regex=False)
    json_df['start_hg38'] = pd.to_numeric(json_df['start_hg38'], errors='coerce')
    json_df['stop_hg38'] = pd.to_numeric(json_df['stop_hg38'], errors='coerce')

    csv_df_copy['Chromosome'] = csv_df_copy['Chromosome'].astype(str).str.replace('chr', '', regex=False)
    csv_df_copy['Position'] = pd.to_numeric(csv_df_copy['Position'], errors='coerce')

    field = 'inheritance'
    output_col = 'Inheritance'

    # Prepare a list to collect results
    results = []

    # For each row in csv_df, find the matching interval in json_df
    for idx, row in csv_df_copy.iterrows():
        chrom = row['Chromosome']
        pos = row['Position']
        match = json_df[
            (json_df['chrom'] == chrom) &
            (json_df['start_hg38'] <= pos) &
            (json_df['stop_hg38'] >= pos)
        ]
        if not match.empty:
            # If multiple matches, take the first
            match_row = match.iloc[0]
            row[output_col] = match_row[field]
        results.append(row)

    merged = pd.DataFrame(results)
    return merged

# Apply the enrichment steps
csv_data_withancestrycolumns = add_population_info(csv_data_copy, tsv_data)
csv_data_withancestrycolumns = add_population_info_with_cleaned(csv_data_withancestrycolumns, tsv_data)
csv_data_withancestrycolumns = add_json_info_by_chrom_pos(csv_data_withancestrycolumns, json_data)

# Drop the ALLR column if it exists
csv_data_withancestrycolumns.drop(columns=['ALLR'], errors='ignore', inplace=True)

# Add Motif length column
csv_data_withancestrycolumns['Motif length'] = csv_data_withancestrycolumns['Motif'].astype(str).str.len()

# Convert Allele length to numeric if it's not already
csv_data_withancestrycolumns['Allele length'] = pd.to_numeric(csv_data_withancestrycolumns['Allele length'], errors='coerce')

# Calculate REPEAT_COUNT = AL / MOTIF_LENGTH
csv_data_withancestrycolumns['Repeat count'] = csv_data_withancestrycolumns.apply(
    lambda row: row['Allele length'] / row['Motif length'] if row['Motif length'] > 0 else None,
    axis=1
)

def flag_pathogenic(row):
    """
    • If either PATHOGENIC_MIN or PATHOGENIC_MAX is missing → return np.nan
    • Otherwise return True/False depending on the repeat-count range
    """
    if (
        pd.isna(row['Pathogenic min']) or
        pd.isna(row['Pathogenic max']) or
        pd.isna(row['Repeat count'])
    ):
        return np.nan          # will appear as a blank cell in the CSV
    return row['Pathogenic min'] <= row['Repeat count'] <= row['Pathogenic max']

csv_data_withancestrycolumns['Is pathogenic'] = (
    csv_data_withancestrycolumns.apply(flag_pathogenic, axis=1)
)

# Ensure Pathogenic min and Pathogenic max are numeric
csv_data_withancestrycolumns['Pathogenic min'] = pd.to_numeric(csv_data_withancestrycolumns['Pathogenic min'], errors='coerce')
csv_data_withancestrycolumns['Pathogenic max'] = pd.to_numeric(csv_data_withancestrycolumns['Pathogenic max'], errors='coerce')

# Define the exact column order
desired_order = [
    'Chromosome', 'Position', 'Gene', 'Disease', 'Motif',
    'Allele length', 'Motif length', 'Repeat count',
    'Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max',
    'Is pathogenic', 'Sample ID'
]

present_desired_cols = [c for c in desired_order if c in csv_data_withancestrycolumns.columns]
remaining_cols       = [c for c in csv_data_withancestrycolumns.columns
                        if c not in present_desired_cols]

final_col_order = present_desired_cols + remaining_cols


csv_data_withancestrycolumns = csv_data_withancestrycolumns[present_desired_cols + remaining_cols]

# Export to CSVs
# output_csv_path = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/68_loci_100_samples_withancestrycolumns.csv"
# output_csv_path = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_88_samples_withancestrycolumns.csv"
output_csv_path = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_503_samples_withancestrycolumns.csv"

csv_data_withancestrycolumns.to_csv(
    output_csv_path,
    index=False
)

# Export to Excel
# output_excel_path = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/68_loci_100_samples_withancestrycolumns.xlsx"
# output_excel_path = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_88_samples_withancestrycolumns.xlsx"
output_excel_path = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_503_samples_withancestrycolumns.xlsx"

csv_data_withancestrycolumns.to_excel(
    output_excel_path,
    index=False
)

# Print progress
print("--- Done ---")