import os
from bs4 import BeautifulSoup

# Initialize an empty HTML document
output_doc = BeautifulSoup("<html><body></body></html>", "html.parser")

input_dir = '/Users/annelisethorn/Documents/Anschutz/Plots/Ancestry_Plots/HTML/83_loci_503_samples_Pore_Dropdown_Correct'

output_dir = '/Users/annelisethorn/Documents/Anschutz/Plots/Ancestry_Plots/HTML/83_loci_503_samples_Pore_Dropdown_Correct'

# Sort the files for consistent order
for file in sorted(os.listdir(input_dir)):
    if not file.lower().endswith('.html'):
        continue

    file_path = os.path.join(input_dir, file)
    with open(file_path, 'r', encoding='utf-8') as html_file:
        file_soup = BeautifulSoup(html_file.read(), "html.parser")
        if file_soup.body:
            output_doc.body.extend(file_soup.body.contents)

# Output the combined result to a file in the output_dir
output_path = os.path.join(output_dir, '83_loci_503_samples_OneFile.html')
os.makedirs(output_dir, exist_ok=True)
with open(output_path, 'w', encoding='utf-8') as out_file:
    out_file.write(output_doc.prettify())
