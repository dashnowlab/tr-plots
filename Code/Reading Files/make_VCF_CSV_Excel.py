import gzip
import csv
import pandas as pd

# vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# csv_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/CSVs/68_loci_100_samples_CSV.csv"
# csv_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/CSVs/83_loci_88_samples_CSV.csv"
csv_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/CSVs/83_loci_503_samples_CSV.csv"

# excel_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/Excels/68_loci_100_samples_Excel.xlsx"
# excel_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/Excels/83_loci_88_samples_Excel.xlsx"
excel_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/Excels/83_loci_503_samples_Excel.xlsx"

# Open the compressed VCF file
with gzip.open(vcf_path, 'rt') as vcf_in, open(csv_path, 'w', newline='') as csv_out:
    writer = csv.writer(csv_out)

    for line in vcf_in:
        line = line.strip()
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            header = line.lstrip("#").split("\t")
            writer.writerow(header)
        else:
            fields = line.split("\t")
            writer.writerow(fields)

# Save as Excel
df = pd.read_csv(csv_path)
df.to_excel(excel_path, index=False)
print(f'--- Done ---')
