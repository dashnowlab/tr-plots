"""
---------------------------------------------
 Script: Read TSV File
 Purpose:
   - Loads a tab-separated file (TSV)
   - Extracts and prints selected columns
---------------------------------------------
"""

import csv
import os

# --- File location ---
TSV_PATH = "/Users/annelisethorn/Documents/GitHub/tr-plots/Data/Other Data/sample_information.tsv"

# --- Helper: read TSV into list of rows ---
def read_tsv(file_path):
    data = []
    with open(file_path, "r") as file:
        tsv_reader = csv.reader(file, delimiter="\t")
        for row in tsv_reader:
            data.append(row)
    return data

# --- Load TSV ---
tsv_data = read_tsv(TSV_PATH)

# --- Extract columns (skip header row) ---
first_column = [row[0] for row in tsv_data[1:]]  
second_column = [row[1] for row in tsv_data[1:]]  
third_column = [row[2] for row in tsv_data[1:]]  
fourth_column = [row[3] for row in tsv_data[1:]]  

# --- Print results ---
print("First column:", first_column[:10], "...")   # show first 10 only
print("Second column:", second_column[:10], "...")
print("Third column:", third_column[:10], "...")
print("Fourth column:", fourth_column[:10], "...")

print(f"--- Done. Parsed {len(tsv_data)-1} rows ---")
