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
import argparse
from pathlib import Path
from trplots.config import OUTPUT_BASE

# --- File locations ---
BASE_DIR = Path(OUTPUT_BASE) / "plots" / "allele_length_boxplots_master"
INPUT_DIR = BASE_DIR / "HTML"
OUTPUT_DIR = INPUT_DIR / "1_file"

def parse_args():
    p = argparse.ArgumentParser(description="Merge individual HTML plots into single HTML")
    p.add_argument("--input-dir", type=str, default=str(INPUT_DIR), help="Directory containing HTML files")
    p.add_argument("--output-file", type=str, default=None, help="Output merged HTML path")
    return p.parse_args()

def main(args=None):
    if args is None:
        args = parse_args()
    input_dir = Path(args.input_dir)
    output_file = Path(args.output_file) if args.output_file else (OUTPUT_DIR / "83_loci_503_samples_OneFile.html")

    # --- Initialize an empty HTML document ---
    output_doc = BeautifulSoup("<html><body></body></html>", "html.parser")

    # --- Ensure input exists ---
    if not input_dir.exists():
        print(f"ERROR: input directory not found: {input_dir}")
        return

    # --- Collect and merge HTML files ---
    files = sorted([f for f in input_dir.iterdir() if f.suffix.lower() == ".html"])
    for file_path in files:
        with open(file_path, "r", encoding="utf-8") as html_file:
            file_soup = BeautifulSoup(html_file.read(), "html.parser")
            if file_soup.body:
                output_doc.body.extend(file_soup.body.contents)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as out_file:
        out_file.write(output_doc.prettify())

    print(f"Saved merged file to: {output_file}")
    print(f"--- Done ---")

if __name__ == "__main__":
    args = parse_args()
    main(args)