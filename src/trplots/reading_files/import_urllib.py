"""
---------------------------------------------
 Script: Download 1000 Genomes Sample Info
 Purpose:
   - Retrieves the sample info file from the
     1000 Genomes Project FTP server
   - Saves it locally for downstream processing
---------------------------------------------
"""

import urllib.request
from pathlib import Path
from trplots.config import OTHER_DATA

# --- File locations ---
URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt"
OUTPUT_DIR = OTHER_DATA
OUTPUT_FILE = OUTPUT_DIR / "1000_genomes_20130606_sample_info.txt"

# --- Prepare output directory ---
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Download file ---
urllib.request.urlretrieve(URL, OUTPUT_FILE)

print(f"--- Done ---")
