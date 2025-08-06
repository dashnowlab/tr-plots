import pysam
import pandas as pd

# Path to the VCF file
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68Samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/88_samples/Datasets/88_samples/1000g-ONT-88Samples.vcf.gz"
VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/503_samples/1000g-ONT-STRchive-503Samples.vcf.gz"

# Open the indexed VCF file
vcf_in = pysam.VariantFile(VCF_PATH)

samples = list(vcf_in.header.samples)

# Collect long-format records
records = []

for record in vcf_in.fetch():
    variant_info = {
        "CHROM": record.contig,
        "POS": record.pos,
        "ID": record.id,
        "REF": record.ref,
        "ALT": ','.join(str(alt) for alt in record.alts),
        "QUAL": record.qual,
        "FILTER": ';'.join(record.filter.keys()),
    }

    # Include INFO fields if desired
    for key, value in record.info.items():
        variant_info[key] = value

    # For each sample, create a separate row
    for sample in samples:
        sample_data = record.samples[sample]
        row = variant_info.copy()
        row["SAMPLE"] = sample

        for fmt_key, fmt_value in sample_data.items():
            row[fmt_key] = fmt_value

        records.append(row)

# Convert to DataFrame
df_long = pd.DataFrame(records)

# Optional: Set multi-index for fast access
df_long.set_index(["CHROM", "POS", "SAMPLE"], inplace=True)

# # Preview
# print(df_long.head(600))

# Save to CSV
df_long.to_csv("vcf_long_format2.csv", index=True)

print(f'--- Done ---')
