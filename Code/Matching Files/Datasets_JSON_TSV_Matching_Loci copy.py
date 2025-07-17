import json
import pysam
import pandas as pd

# File Paths
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/68_loci_100_samples/100HPRC.trgt-v0.8.0.STRchive.sorted-68_loci_100_samples.vcf.gz"
# VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_88_samples/1000g-ONT-83_loci_88_samples.vcf.gz"
VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/83_loci_503_samples/1000g-ONT-STRchive-83_loci_503_samples.vcf.gz"

JSON_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json"
CSV_PATH = '/Users/annelisethorn/Documents/Anschutz/Datasets/Ancestry_1KGP_ONT_500_Summary - Sheet1.csv'

# Load Sample Population Information from CSV
sample_info_CSV = pd.read_csv(CSV_PATH)

# Use sample_id 
sample_info_CSV = sample_info_CSV.iloc[:, [0, 2, 3]]
sample_info_CSV.columns = ['Sample_ID', 'SubPop', 'SuperPop']

# Create Mapping Dictionary (remove duplicates first)
population_map = (
    sample_info_CSV
    .drop_duplicates(subset="Sample_ID")
    .set_index("Sample_ID")[["SubPop", "SuperPop"]]
    .to_dict(orient="index")
)

# Load Loci Metadata
with open(JSON_PATH, "r") as f:
    loci_data = json.load(f)

# Read VCF
vcf_in = pysam.VariantFile(VCF_PATH)
vcf_sample_ids = list(vcf_in.header.samples)
base_sample_ids = [sid.split('-')[0] for sid in vcf_sample_ids]  # Normalize

# Locus-Sample Mapping Collection
records = []

total_vcf_loci = 0
matched_vcf_loci = 0

for record in vcf_in.fetch():
    total_vcf_loci += 1
    chrom = record.chrom
    pos = record.pos

    matched_locus = next(
        (entry for entry in loci_data
         if chrom == entry["chrom"] and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
        None
    )
    # if matched_locus:
    #     print(f"Matched VCF locus {chrom}:{pos}")

    if not matched_locus:
        print(f"VCF locus {chrom}:{pos} not matched in loci_data")
        continue
    matched_vcf_loci += 1

    # Calculate motif length as # of bases
    motif = matched_locus.get("reference_motif_reference_orientation")
    # If motif is a list (e.g., ['GCC']), join to string and strip quotes/brackets
    if isinstance(motif, list):
        motif_str = "".join(motif)
    else:
        motif_str = str(motif)
    motif_length = len(motif_str)
    motif = motif_str

    # print(f"Processing {chrom}:{pos} - {motif} ({motif_length})")

    for sid, sample in record.samples.items():
        base_id = sid.split('-')[0]
        al_lengths = sample.get("AL")
        if not al_lengths:
            continue

        repeat_counts = [(allele_len) // motif_length for allele_len in al_lengths if allele_len]
        pop_info = population_map.get(base_id, {"SubPop": "Unknown", "SuperPop": "Unknown"})

        # Example: Assume matched_locus has 'pathogenic_min' and 'pathogenic_max' keys
        pathogenic_min = matched_locus.get("pathogenic_min", None)
        pathogenic_max = matched_locus.get("pathogenic_max", None)

        # Ensure pathogenic_min is less than or equal to pathogenic_max
        # Determine for each repeat count if it falls within the pathogenic range
        if pathogenic_min is not None and pathogenic_max is not None:
            is_pathogenic = [
            pathogenic_min <= rc <= pathogenic_max for rc in repeat_counts
            ]
        else:
            is_pathogenic = [False] * len(repeat_counts)
        
        # # double check repeat_counts are correctly calculated
        # print(f"Sample {sid} - Repeat Counts: {repeat_counts}, Pathogenic Min: {pathogenic_min}, Pathogenic Max: {pathogenic_max}, Pathogenic Pathogenic: {is_pathogenic}")

        records.append({
            "locus_chrom": chrom,
            "locus_pos": pos,
            "sample_id": sid,
            "disease": matched_locus.get("disease", "Unknown"),
            "gene": matched_locus.get("gene", "Unknown"),
            "population": pop_info["SubPop"],
            "population_description": pop_info["SuperPop"],
            "motif": motif,
            "motif_length": motif_length,
            "allele_lengths": ";".join(map(str, al_lengths)),
            "repeat_counts": ";".join(map(str, repeat_counts)),
            "is_pathogenic": ";".join(map(str, is_pathogenic)),
            "pathogenic_min": pathogenic_min,
            "pathogenic_max": pathogenic_max
        })

# Ensure that each allele length and repeat count is represented in the final DataFrame
expanded_records = []
skipped = 0
for record in records:
    allele_lengths = record["allele_lengths"].split(";")
    repeat_counts = record["repeat_counts"].split(";")
    is_pathogenic = record["is_pathogenic"].split(";")
    if len(allele_lengths) != len(repeat_counts) or len(allele_lengths) != len(is_pathogenic):
        skipped += 1
        continue

    for i in range(len(allele_lengths)):
        expanded_records.append({
            "gene": record["gene"],
            "locus_chrom": record["locus_chrom"],
            "locus_pos": record["locus_pos"],
            "sample_id": record["sample_id"],
            "disease": record["disease"],
            "population": record["population"],
            "population_description": record["population_description"],
            "motif": record["motif"],
            "motif_length": record["motif_length"],
            "allele_length": allele_lengths[i],
            "repeat_count": repeat_counts[i],
            "pathogenic_min": record["pathogenic_min"],
            "pathogenic_max": record["pathogenic_max"],
            "is_pathogenic": is_pathogenic[i] == "True"  # Convert to boolean if needed
        })  

# # Save to CSV
df = pd.DataFrame(expanded_records)
# output_file = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_88_samples_Matching_JSON_TSV.csv"
output_file = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/CSVs/83_loci_503_samples_Matching_JSON_CSV.csv"
df.to_csv(output_file, index=False)
print(f"Saved CSV File")

# Save to Excel
# excel_output_file = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_88_samples_Matching_JSON_TSV.xlsx"
excel_output_file = "/Users/annelisethorn/Documents/Anschutz/Code/Matching Files/Excels/83_loci_503_samples_Matching_JSON_CSV.xlsx"
df.to_excel(excel_output_file, index=False)
print(f"Saved Excel File")

# After loading df
if 'locus_chrom' in df.columns and 'locus_pos' in df.columns:
    unique_loci = set(zip(df['locus_chrom'], df['locus_pos']))
    print(f"Number of unique loci in processed DataFrame: {len(unique_loci)}")
else:
    print("locus_chrom and locus_pos columns not found in DataFrame.")

print(f"Total loci in VCF: {total_vcf_loci}")
print(f"Loci matched to loci_data: {matched_vcf_loci}")

print("--- Done ---")
