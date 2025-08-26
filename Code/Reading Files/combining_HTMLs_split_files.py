"""
---------------------------------------------
 Script: Split Merged HTML Plots into Chunks
 Purpose:
   - Reads multiple individual HTML plot files from a directory
   - Extracts and flattens their <body> contents
   - Splits combined content into evenly sized chunks
   - Writes each chunk into its own HTML file
   - Saves the chunked outputs into a specified subdirectory
---------------------------------------------
"""

import os
import math
from bs4 import BeautifulSoup

# --- File locations ---
BASE_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots/Plots/Ancestry_Plots/HTML"
INPUT_DIR = os.path.join(BASE_DIR, "83_loci_503_samples_People_AllColumn")
OUTPUT_DIR = os.path.join(INPUT_DIR, "SplitFiles")

# --- Prepare output directory ---
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Collect valid HTML files ---
html_files = sorted([f for f in os.listdir(INPUT_DIR) if f.lower().endswith(".html")])

# --- Load and store <body> contents from each file ---
parsed_bodies = []
for file in html_files:
    file_path = os.path.join(INPUT_DIR, file)
    with open(file_path, "r", encoding="utf-8") as html_file:
        file_soup = BeautifulSoup(html_file.read(), "html.parser")
        if file_soup.body:
            parsed_bodies.append(file_soup.body.contents)

# --- Flatten list of contents across all files ---
flat_contents = [item for sublist in parsed_bodies for item in sublist]

# --- Split into 5 evenly sized chunks ---
chunk_size = math.ceil(len(flat_contents) / 5)
chunks = [flat_contents[i:i + chunk_size] for i in range(0, len(flat_contents), chunk_size)]

# --- Write each chunk into its own HTML file ---
for i, chunk in enumerate(chunks, start=1):
    output_doc = BeautifulSoup("<html><body></body></html>", "html.parser")
    output_doc.body.extend(chunk)

    output_path = os.path.join(OUTPUT_DIR, f"{i}_83_loci_503_samples_ancestry_plots.html")
    with open(output_path, "w", encoding="utf-8") as out_file:
        out_file.write(output_doc.prettify())

print(f"--- Split {len(html_files)} input files into {len(chunks)} chunked HTML files ---")