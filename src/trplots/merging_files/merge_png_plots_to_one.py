"""
---------------------------------------------
 Script: Combine PNG Plots into Single HTML
 Purpose:
   - Reads multiple PNG files from a directory
   - Embeds each image in an HTML <img> tag with caption
   - Writes all images into a single HTML file
---------------------------------------------
"""

import os
from pathlib import Path
from trplots.config import OUTPUT_BASE

# --- File locations ---
BASE_DIR = Path(OUTPUT_BASE) / "plots" / "allele_length_boxplots_master"
INPUT_DIR = BASE_DIR / "PNG"
OUTPUT_HTML = INPUT_DIR / "83_loci_503_samples_OneFile.html"

# --- Collect PNG files ---
png_files = sorted([f for f in os.listdir(INPUT_DIR) if f.lower().endswith(".png")])

# --- Build HTML content ---
html = ['<html><head><title>PNG Plots</title></head><body>']
for png in png_files:
    caption = os.path.splitext(png)[0]
    img_path = Path(INPUT_DIR) / png
    html.append('<div style="margin-bottom:40px;text-align:center;">')
    html.append(f'<img src="file://{img_path.resolve()}" style="max-width:90vw;max-height:80vh;display:block;margin:auto;" alt="{caption}"/>')
    html.append(f'<div style="font-family:Helvetica,Arial,sans-serif;font-size:14px;margin-top:8px;">{caption}</div>')
    html.append('</div>')
html.append('</body></html>')

# --- Write combined HTML file ---
with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
    f.write('\n'.join(html))

print(f"--- Combined {len(png_files)} PNG files into {OUTPUT_HTML} ---")
