import pysam
import pandas as pd

# Path to your uncompressed VCF file
vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
# vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# Path to sample ancestry, technology, etc. file
csv_path = '/Users/annelisethorn/Documents/Anschutz/Datasets/1KGP_ONT_500_Summary - Sheet1.csv'

# Open VCF
vcf = pysam.VariantFile(vcf_path)

# Open CSV the second row is the header
# Read the CSV file
df = pd.read_csv(csv_path, header=1)

# print(df.head())

# Drop all columns that aren't 'Sample_ID' or 'Pore'
df = df[['Sample_ID', 'Pore']]

# print(df.head())

# Export the DataFrame to a CSV file
df.to_csv('/Users/annelisethorn/Documents/Anschutz/Datasets/1KGP_ONT_500_Summary_Sample_ID_Pore.csv', index=False)

