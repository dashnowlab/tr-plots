"""
---------------------------------------------
 Script: Merge HTML Plot Files
 Purpose:
   - Reads multiple individual HTML plot files from a directory
   - Extracts their <body> contents
   - Appends them together into a single combined HTML document
   - Saves the merged file into the specified output directory
---------------------------------------------
"""

import os
from bs4 import BeautifulSoup
from pathlib import Path
from trplots.config import OUTPUT_BASE

# --- File locations ---
BASE_DIR = Path(OUTPUT_BASE) / "plots" / "allele_length_boxplots_master"
INPUT_DIR = BASE_DIR / "HTML"
OUTPUT_DIR = INPUT_DIR / "1_file"

# --- Initialize an empty HTML document ---
output_doc = BeautifulSoup("<html><body></body></html>", "html.parser")

# --- Collect and merge HTML files ---
for file in sorted(os.listdir(INPUT_DIR)):
    if not file.lower().endswith(".html"):
        continue

    file_path = os.path.join(INPUT_DIR, file)
    with open(file_path, "r", encoding="utf-8") as html_file:
        file_soup = BeautifulSoup(html_file.read(), "html.parser")
        if file_soup.body:
            output_doc.body.extend(file_soup.body.contents)

# --- Save combined HTML output ---
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_path = os.path.join(OUTPUT_DIR, "83_loci_503_samples_OneFile.html")

with open(output_path, "w", encoding="utf-8") as out_file:
    out_file.write(output_doc.prettify())

print(f"--- Done ---")