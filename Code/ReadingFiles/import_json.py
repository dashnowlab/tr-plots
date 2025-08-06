import json

# Open and read the JSON file
with open('/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json', 'r') as file:
    data = json.load(file)

# Print the data
# print(data)

# Extract and print relevant fields
for entry in data:
    chrom = entry.get('chrom')
    start = entry.get('start_hg38')
    benign_min = entry.get('benign_min')
    benign_max = entry.get('benign_max')
    pathogenic_min = entry.get('pathogenic_min')
    pathogenic_max = entry.get('pathogenic_max')

    # Format values: use int if available, else 'null'
    def fmt(val):
        return int(val) if val is not None else 'null'
    
    print(f"Chromosome: {chrom}, Start: {start}, "
          f"Benign: {fmt(benign_min)}-{fmt(benign_max)}, "
          f"Pathogenic: {fmt(pathogenic_min)}-{fmt(pathogenic_max)}")