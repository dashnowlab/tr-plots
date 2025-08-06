from re import match
import pandas as pd
import json
import numpy as np

# File paths
CSV_PATH = '/Users/annelisethorn/Documents/Anschutz/Debug/debug_info_83_loci_503_samples.csv'
CSV_PATH_2 = '/Users/annelisethorn/Documents/Anschutz/Datasets/Ancestry_1KGP_ONT_500_Summary - Sheet1.csv'
JSON_PATH = '/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json'

# Load data
csv_data = pd.read_csv(CSV_PATH)
csv_df_copy = csv_data.copy()

sample_info_csv = pd.read_csv(CSV_PATH_2)

# Add population info
def add_population_info(df, sample_info_csv):
    first_col = sample_info_csv.columns[0]
    third_col = sample_info_csv.columns[2]
    fourth_col = sample_info_csv.columns[3]

    population_map = sample_info_csv.set_index(first_col)[[third_col, fourth_col]].to_dict(orient='index')

    df['Sample ID'] = df['Sample ID'].astype(str).str.strip()
    df['Sample ID Cleaned'] = df['Sample ID'].apply(remove_characters_after_first_dash)

    df['SubPop'] = df['Sample ID Cleaned'].map(lambda x: population_map.get(x, {}).get(third_col, 'Unknown'))
    df['SuperPop'] = df['Sample ID Cleaned'].map(lambda x: population_map.get(x, {}).get(fourth_col, 'Unknown'))

    return df

# Remove characters in Sample IDs in csv_df_copy at the first -
def remove_characters_after_first_dash(sample_id):
    if '-' in sample_id:
        return sample_id.split('-')[0]
    return sample_id

# Match Sample IDs and add SubPop and SuperPop columns
df = add_population_info(csv_df_copy, sample_info_csv) 

# Match Chromosome and position with JSON data and add Inheritance and Disease columns
def match_chrom_pos_and_add_info(df, json_path):
    with open(json_path, 'r') as f:
        loci_data = json.load(f)

    # Create a dictionary for quick lookup
    loci_dict = {(locus['chrom'], locus['start_hg38'], locus['stop_hg38']): locus for locus in loci_data}

    # Initialize new columns
    df['Inheritance'] = pd.Series(dtype='object')
    df['Disease'] = pd.Series(dtype='object')


    for index, row in df.iterrows():
        chrom = row['Chromosome']
        pos = row['Position']

        # Find matching locus
        matched_locus = next((locus for (locus_chrom, start, stop), locus in loci_dict.items() 
                              if locus_chrom == chrom and start <= pos <= stop), None)

        if matched_locus:
            df.at[index, 'Inheritance'] = matched_locus.get('inheritance', 'Unknown')
            df.at[index, 'Disease'] = matched_locus.get('disease', 'Unknown')

    return df

# Match Chromosome and position with JSON data and add Inheritance and Disease columns
df = match_chrom_pos_and_add_info(df, JSON_PATH)

# Export
df.to_csv("/Users/annelisethorn/Documents/Anschutz/6.24/83_loci_503_samples_withancestrycolumns_5.csv", index=False)
df.to_excel("/Users/annelisethorn/Documents/Anschutz/6.24/83_loci_503_samples_withancestrycolumns_5.xlsx", index=False)

print("--- Done ---")
