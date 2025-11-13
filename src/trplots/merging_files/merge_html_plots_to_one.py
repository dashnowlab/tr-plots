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
import math

# --- File locations ---
# /Users/annelisethorn/Documents/github/tr-plots/results/plots/Path_Ref_Motif_Tandem_Repeats_Plots_v2/HTML
#
# BASE_DIR = Path(OUTPUT_BASE) / "plots" / "Path_Ref_Motif_Tandem_Repeats_Plots_v2"
# INPUT_DIR = BASE_DIR / "HTML"
# Use the actual directory you provided (falls back to OUTPUT_BASE if needed)
INPUT_DIR = Path("/Users/annelisethorn/Documents/github/tr-plots/results/plots/Path_Ref_Motif_Tandem_Repeats_Plots_v2/HTML")
if not INPUT_DIR.exists():
    BASE_DIR = Path(OUTPUT_BASE) / "plots" / "Path_Ref_Motif_Tandem_Repeats_Plots_v2"
    INPUT_DIR = BASE_DIR / "HTML"

OUTPUT_DIR = INPUT_DIR / "1_file"

def parse_args():
    p = argparse.ArgumentParser(description="Merge individual HTML plots into single HTML")
    p.add_argument("--input-dir", type=str, default=str(INPUT_DIR), help="Directory containing HTML files")
    p.add_argument("--output-file", type=str, default=None, help="Output merged HTML path")
    p.add_argument("--parts", type=int, default=5, help="Split merged output into N files (default: 5)")
    return p.parse_args()

def main(args=None):
    if args is None:
        args = parse_args()
    input_dir = Path(args.input_dir)
    output_file_arg = args.output_file
    parts = max(1, int(args.parts or 1))

    # --- Initialize an empty HTML document ---
    # (we'll create separate docs per part)

    # --- Ensure input exists ---
    if not input_dir.exists():
        print(f"ERROR: input directory not found: {input_dir}")
        return

    # --- Collect and merge HTML files ---
    files = sorted([f for f in input_dir.iterdir() if f.suffix.lower() == ".html"])
    total = len(files)
    if total == 0:
        print(f"No HTML files found in {input_dir}")
        return

    # compute chunking
    chunk_size = math.ceil(total / parts)
    for part_idx in range(parts):
        start = part_idx * chunk_size
        if start >= total:
            break
        end = min(start + chunk_size, total)
        chunk_files = files[start:end]

        # build output path for this part
        if output_file_arg:
            # allow formatting with {part}
            if "{part}" in output_file_arg:
                out_path = Path(output_file_arg.format(part=part_idx + 1))
            else:
                base = Path(output_file_arg)
                out_path = base.parent / f"{base.stem}_part{part_idx + 1}.html"
        else:
            out_path = OUTPUT_DIR / f"Path_Ref_Motif_Allele_Dist_part{part_idx + 1}.html"

        # create a fresh document and append chunk contents
        output_doc = BeautifulSoup("<html><body></body></html>", "html.parser")
        for file_path in chunk_files:
            with open(file_path, "r", encoding="utf-8") as html_file:
                file_soup = BeautifulSoup(html_file.read(), "html.parser")
                if file_soup.body:
                    output_doc.body.extend(file_soup.body.contents)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as out_file:
            out_file.write(output_doc.prettify())
        print(f"Saved merged part {part_idx + 1}/{parts} to: {out_path}")
    print(f"--- Done: merged {total} files into {min(parts, math.ceil(total/chunk_size))} parts ---")

if __name__ == "__main__":
    args = parse_args()
    main(args)