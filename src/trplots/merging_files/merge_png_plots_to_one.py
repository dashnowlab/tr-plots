"""
---------------------------------------------
 Script: Combine PNG Plots into Single HTML
 Purpose:
   - Reads multiple PNG files from a directory
   - Embeds each image in an HTML <img> tag with caption
   - Writes all images into a single HTML file
---------------------------------------------
"""

from pathlib import Path
from trplots.config import OUTPUT_BASE
import argparse
import os


def parse_args():
    p = argparse.ArgumentParser(description="Combine PNG plots into single HTML")
    p.add_argument("--input-dir", type=str, default=str(INPUT_DIR), help="PNG input directory")
    p.add_argument("--output-html", type=str, default=str(OUTPUT_HTML), help="Output HTML path")
    return p.parse_args()

# --- File locations ---
BASE_DIR = Path(OUTPUT_BASE) / "plots" / "allele_length_boxplots_master"
INPUT_DIR = BASE_DIR / "PNG"
OUTPUT_HTML = INPUT_DIR / "83_loci_503_samples_OneFile.html"

def main(args=None):
    if args is None:
        args = parse_args()
    input_dir = Path(args.input_dir)
    output_html = Path(args.output_html)

    # Ensure input directory exists
    if not input_dir.exists():
        print(f"ERROR: input directory not found: {input_dir}")
        return

    # --- Collect PNG files ---
    png_files = sorted([f for f in os.listdir(input_dir) if f.lower().endswith(".png")])

    # --- Build HTML content ---
    html = ['<html><head><title>PNG Plots</title></head><body>']
    for png in png_files:
        caption = os.path.splitext(png)[0]
        img_path = Path(input_dir) / png
        html.append('<div style="margin-bottom:40px;text-align:center;">')
        html.append(f'<img src="file://{img_path.resolve()}" style="max-width:90vw;max-height:80vh;display:block;margin:auto;" alt="{caption}"/>')
        html.append(f'<div style="font-family:Helvetica,Arial,sans-serif;font-size:14px;margin-top:8px;">{caption}</div>')
        html.append('</div>')
    html.append('</body></html>')

    output_html.parent.mkdir(parents=True, exist_ok=True)
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write('\n'.join(html))

    print(f"--- Combined {len(png_files)} PNG files into {output_html} ---")

if __name__ == "__main__":
    args = parse_args()
    main(args)
