import pysam

# Path to the VCF file
VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# Open the VCF file with pysam
vcf = pysam.VariantFile(VCF_PATH)

# Use a set to track unique loci as (chromosome, position)
unique_loci = set()

# Iterate through each record in the VCF
for record in vcf.fetch():
    unique_loci.add((record.contig, record.pos))

# Print the result
print(f"Number of unique loci: {len(unique_loci)}")
