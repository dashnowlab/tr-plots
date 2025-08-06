import pandas as pd
import json
import numpy as np

# File paths
CSV_PATH = "Data/debug_info_83_loci_503_samples.csv"
JSON_PATH = "Data/STRchive-loci-v2.7.0.json"
SAMPLE_METADATA_PATH = "Data/igsr_samples.tsv" # Data is from here: https://www.internationalgenome.org/data-portal/sample

# Load CSV (VCF data) and use the first row as column headers
csv_data = pd.read_csv(CSV_PATH, sep=",", header=0)

# Make a copy so we don't modify the original data
csv_data_copy = csv_data.copy()

# Load JSON metadata for loci
with open(JSON_PATH, "r") as file:
    json_data = json.load(file)

# Load Sample Population Information (ancestry data)
sample_data = pd.read_csv(SAMPLE_METADATA_PATH, sep="\t", skiprows=1)

# Pre-create population mapping dictionary
def create_population_mapping(sample_metadata_df):
    """Create mapping dictionaries once instead of repeatedly accessing DataFrame"""
    tsv_first_col = sample_metadata_df.columns[0]
    tsv_third_col = sample_metadata_df.columns[2]
    tsv_fourth_col = sample_metadata_df.columns[3]

    population_dict = {}
    description_dict = {}

    for _, row in sample_metadata_df.iterrows():
        sample_id = row[tsv_first_col]
        population_dict[sample_id] = row[tsv_third_col]
        description_dict[sample_id] = row[tsv_fourth_col]

    return population_dict, description_dict

# Vectorized population mapping
def add_population_info_optimized(csv_df_copy, sample_metadata_df):
    population_dict, description_dict = create_population_mapping(sample_metadata_df)
    # remove everything in sample ID after the first -
    csv_df_copy['Sample ID'] = csv_df_copy['Sample ID'].astype(str).str.strip().str.split('-').str[0]
    # Replace GM with NA at the start of Sample ID
    csv_df_copy['Sample ID'] = csv_df_copy['Sample ID'].str.replace('^GM', 'NA', regex=True).str.strip()
    csv_df_copy['Population'] = csv_df_copy['Sample ID'].map(population_dict).fillna('Unknown')
    csv_df_copy['Population description'] = csv_df_copy['Sample ID'].map(description_dict).fillna('Unknown')

    return csv_df_copy

# Pre-create cleaning mapping dictionary
POPULATION_CLEANING_MAP = {
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

def add_population_info_with_cleaned_optimized(csv_df_copy):
    # Use vectorized map instead of apply
    csv_df_copy['Cleaned population description'] = csv_df_copy['Population description'].map(
        POPULATION_CLEANING_MAP
    ).fillna(csv_df_copy['Population description'])
    return csv_df_copy

# Faster interval matching
def add_json_info_by_chrom_pos_optimized(csv_df_copy, json_data):
    json_df = pd.DataFrame(json_data)
    json_df['chrom'] = json_df['chrom'].astype(str).str.replace('chr', '', regex=False)
    json_df['start_hg38'] = pd.to_numeric(json_df['start_hg38'], errors='coerce')
    json_df['stop_hg38'] = pd.to_numeric(json_df['stop_hg38'], errors='coerce')

    csv_df_copy['Chromosome'] = csv_df_copy['Chromosome'].astype(str).str.replace('chr', '', regex=False)
    csv_df_copy['Position'] = pd.to_numeric(csv_df_copy['Position'], errors='coerce')

    # Create a mapping dictionary for much faster lookup
    inheritance_map = {}

    for _, json_row in json_df.iterrows():
        chrom = json_row['chrom']
        start = json_row['start_hg38']
        stop = json_row['stop_hg38']
        inheritance = json_row.get('inheritance', None)

        if chrom not in inheritance_map:
            inheritance_map[chrom] = []
        inheritance_map[chrom].append((start, stop, inheritance))

    # Sort intervals by start position for each chromosome
    for chrom in inheritance_map:
        inheritance_map[chrom].sort(key=lambda x: x[0])
    
    def find_inheritance(row):
        chrom = row['Chromosome']
        pos = row['Position']

        if pd.isna(pos) or chrom not in inheritance_map:
            return None

        # Binary search would be even faster, but linear search through sorted intervals is good enough for this sized dataset
        for start, stop, inheritance in inheritance_map[chrom]:
            if start <= pos <= stop:
                return inheritance
        return None

    csv_df_copy['Inheritance'] = csv_df_copy.apply(find_inheritance, axis=1)
    return csv_df_copy

# Main processing
print("Adding population info...")
csv_data_withancestrycolumns = add_population_info_optimized(csv_data_copy, sample_data)

print("Cleaning population names...")
csv_data_withancestrycolumns = add_population_info_with_cleaned_optimized(csv_data_withancestrycolumns)

print("Adding inheritance info...")
csv_data_withancestrycolumns = add_json_info_by_chrom_pos_optimized(csv_data_withancestrycolumns, json_data)

print("Processing remaining columns...")

# Drop the ALLR column if it exists
csv_data_withancestrycolumns.drop(columns=['ALLR'], errors='ignore', inplace=True)

# Add Motif length column
csv_data_withancestrycolumns['Motif length'] = csv_data_withancestrycolumns['Motif'].astype(str).str.len()

# Convert Allele length to numeric if it's not already
csv_data_withancestrycolumns['Allele length'] = pd.to_numeric(csv_data_withancestrycolumns['Allele length'], errors='coerce')

# Remove any existing incorrect column
csv_data_withancestrycolumns.drop(columns=['Repeat count'], errors='ignore', inplace=True)

# Vectorized repeat count calculation
mask = csv_data_withancestrycolumns['Motif length'] > 0
csv_data_withancestrycolumns['Repeat count'] = np.where(
    mask,
    (csv_data_withancestrycolumns['Allele length'] / csv_data_withancestrycolumns['Motif length']).round(2),
    np.nan
)

# OPTIMIZED: Vectorized pathogenic flagging
csv_data_withancestrycolumns['Pathogenic min'] = pd.to_numeric(csv_data_withancestrycolumns['Pathogenic min'], errors='coerce')
csv_data_withancestrycolumns['Pathogenic max'] = pd.to_numeric(csv_data_withancestrycolumns['Pathogenic max'], errors='coerce')

# Vectorized pathogenic check
valid_ranges = (
    csv_data_withancestrycolumns['Pathogenic min'].notna() & 
    csv_data_withancestrycolumns['Pathogenic max'].notna() & 
    csv_data_withancestrycolumns['Repeat count'].notna()
)

csv_data_withancestrycolumns['Is pathogenic'] = np.where(
    valid_ranges,
    (csv_data_withancestrycolumns['Pathogenic min'] <= csv_data_withancestrycolumns['Repeat count']) & 
    (csv_data_withancestrycolumns['Repeat count'] <= csv_data_withancestrycolumns['Pathogenic max']),
    np.nan
)

# Define the exact column order
desired_order = [
    'Chromosome', 'Position', 'Gene', 'Disease', 'Motif',
    'Allele length', 'Motif length', 'Repeat count',
    'Benign min', 'Benign max', 'Pathogenic min', 'Pathogenic max',
    'Is pathogenic', 'Sample ID'
]

present_desired_cols = [c for c in desired_order if c in csv_data_withancestrycolumns.columns]
remaining_cols = [c for c in csv_data_withancestrycolumns.columns if c not in present_desired_cols]

csv_data_withancestrycolumns = csv_data_withancestrycolumns[present_desired_cols + remaining_cols]

# Export to CSV
output_path = "Output/83_loci_503_samples_withancestrycolumns_2.csv"
print(f"Exporting to output file: {output_path}")

csv_data_withancestrycolumns.to_csv(
    output_path,
    index=False
)

print("--- Done ---")
