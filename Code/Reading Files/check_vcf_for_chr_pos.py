import pysam

# Define loci of interest (chromosome, position)
user_loci = {
    ("chr3", 63912684),
    ("chr13", 70139383),
    ("chr13", 70139429),
    ("chrX", 147912049),
    ("chr9", 69037286),
    ("chr9", 69037304),
    ("chr5", 10356343),
    ("chr1", 155188505),
    ("chr1", 156591765),
    ("chr15", 88569433),
    ("chr20", 4699397),
    ("chr20", 4699493),
    ("chr17", 17808358),
    ("chrX", 71453054),
}

# Path to your uncompressed VCF file
# vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
vcf_path = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# Open VCF
vcf = pysam.VariantFile(vcf_path)

# Check which loci are present
found = set()
for record in vcf.fetch():
    if (record.chrom, record.pos) in user_loci:
        found.add((record.chrom, record.pos))

# Report results
for chrom, pos in user_loci:
    status = "FOUND" if (chrom, pos) in found else "MISSING"
    print(f"{chrom}:{pos}\t{status}")
