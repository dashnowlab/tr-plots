import os
from bs4 import BeautifulSoup
import math
import zipfile

# Initialize an empty list to hold parsed HTML file bodies
parsed_bodies = []

input_dir = '/Users/annelisethorn/Documents/Anschutz/Plots/Ancestry_Plots/HTML/83_loci_503_samples_People_AllColumn'
output_dir = '/Users/annelisethorn/Documents/Anschutz/Plots/Ancestry_Plots/HTML/83_loci_503_samples_People_AllColumn/SplitFiles'

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Sort and read all valid HTML files
html_files = sorted([f for f in os.listdir(input_dir) if f.lower().endswith('.html')])

# Load and store body contents of each HTML file
for file in html_files:
    file_path = os.path.join(input_dir, file)
    with open(file_path, 'r', encoding='utf-8') as html_file:
        file_soup = BeautifulSoup(html_file.read(), "html.parser")
        if file_soup.body:
            parsed_bodies.append(file_soup.body.contents)

# Flatten list of contents
flat_contents = [item for sublist in parsed_bodies for item in sublist]

# Split contents into parts
chunk_size = math.ceil(len(flat_contents) / 5)
chunks = [flat_contents[i:i + chunk_size] for i in range(0, len(flat_contents), chunk_size)]

# Write each chunk to a separate HTML file in the output directory
for i, chunk in enumerate(chunks, start=1):
    output_doc = BeautifulSoup("<html><body></body></html>", "html.parser")
    output_doc.body.extend(chunk)
    output_path = os.path.join(output_dir, f'{i}_83_loci_503_samples_ancestry_plots.html')
    with open(output_path, 'w', encoding='utf-8') as out_file:
        out_file.write(output_doc.prettify())
