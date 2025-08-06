import pandas as pd
import pysam
import json

VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

# Load the JSON
with open("/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json") as f:
    loci_data = json.load(f)

# Parse the compressed and indexed VCF
vcf_in = pysam.VariantFile(VCF_PATH)

# Map matches
matches = []
for record in vcf_in:
    chrom = record.chrom
    pos = record.pos
    for locus in loci_data:
        if chrom == locus["chrom"] and locus["start_hg38"] <= pos <= locus["stop_hg38"]:
            matches.append({
                "chrom": record.chrom,
                "pos": record.pos,
                # "TRID": record.info.get("TRID"),
                # "matched_disease": locus["disease"],
                # "gene": locus["gene"],
                # "motif": locus["reference_motif_reference_orientation"]
            })

df_matches = pd.DataFrame(matches)
# print(df_matches)

# Save DataFrame to CSV
df_matches.to_csv("/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/68_loci_100_samples_Matching_JSON.csv", index=False)
# df_matches.to_csv("/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_88_samples_Matching_JSON.csv", index=False)
# df_matches.to_csv("/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_503_samples_Matching_JSON.csv", index=False)

# Save DataFrame to Excel
df_matches.to_excel("/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/68_loci_100_samples_Matching_JSON.xlsx", index=False)
# df_matches.to_excel("/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_88_samples_Matching_JSON.xlsx", index=False)
# df_matches.to_excel("/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_503_samples_Matching_JSON.xlsx", index=False)
print(f'--- Done ---')
