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
import os

# --- File locations ---
URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt"
OUTPUT_DIR = "/Users/annelisethorn/Documents/GitHub/tr-plots/Data/Other Data"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "1000_Genomes_20130606_sample_info.txt")

# --- Prepare output directory ---
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Download file ---
urllib.request.urlretrieve(URL, OUTPUT_FILE)

print(f"--- Done ---")
