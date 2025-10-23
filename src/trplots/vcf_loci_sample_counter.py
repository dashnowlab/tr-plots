import gzip
from collections import defaultdict
import os

# The VCF file provided by the user (it was inside a .vcf.gz)
VCF_FILE_PATH = "/Users/annelisethorn/Documents/github/tr-plots/data/sequencing_data/83_loci_503_samples/1000g-ont-strchive-83_loci_503_samples.vcf.gz"

def count_samples_per_locus(file_path):
    """
    Reads a VCF file, identifies sample IDs, and counts the number
    of samples that have a genotype call (i.e., not './.') for each locus.
    
    Locus is defined as CHROM:POS.
    """
    
    print(f"Reading VCF file: {os.path.basename(file_path)}...")
    
    # Check if the file path implies gzipped content
    is_gzipped = file_path.endswith('.gz')
    
    # Use appropriate file opening mechanism
    if is_gzipped:
        open_func = gzip.open
        mode = 'rt' # read text mode for gzip
    else:
        open_func = open
        mode = 'r'

    locus_sample_counts = defaultdict(int)
    sample_ids = []
    
    try:
        with open_func(file_path, mode) as f:
            for line in f:
                # 1. Skip metadata lines
                if line.startswith('##'):
                    continue
                
                # 2. Identify the header line and extract sample IDs
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    # Sample IDs start from the 10th column (index 9)
                    sample_ids = parts[9:]
                    num_samples = len(sample_ids)
                    print(f"Found {num_samples} total sample IDs in the header.")
                    continue
                
                # 3. Process data lines (one per locus)
                if sample_ids:
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    pos = parts[1]
                    locus_key = f"{chrom}:{pos}"
                    
                    # Genotype data starts at index 9
                    genotype_data = parts[9:]
                    
                    # Count samples with a non-missing genotype call
                    # The VCF standard format for missing genotype is "./."
                    samples_with_calls = 0
                    
                    # The genotype format is parts[8]. We need to find the index of 'GT'
                    format_keys = parts[8].split(':')
                    try:
                        gt_index = format_keys.index('GT')
                    except ValueError:
                        # If GT is not present in the FORMAT field, assume genotype is missing for all
                        gt_index = -1 
                        
                    for sample_genotype_string in genotype_data:
                        if gt_index != -1:
                            # Extract GT field from the sample string (e.g., '0/1:100:...')
                            sample_fields = sample_genotype_string.split(':')
                            if len(sample_fields) > gt_index:
                                genotype = sample_fields[gt_index]
                                # Check for missing genotype call. Genotypes like '0/0', '0/1', '1/1', etc. are calls.
                                if genotype not in ('./.', '.|.', '.'):
                                    samples_with_calls += 1
                                # NOTE: This simple check assumes any non-missing GT means the sample was genotyped.
                                # For complex analysis (e.g., quality filters), more detailed parsing is needed, 
                                # but this addresses the user's immediate question.
                            
                    locus_sample_counts[locus_key] = samples_with_calls

    except FileNotFoundError:
        print(f"Error: The VCF file was not found at the expected path: {file_path}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return None

    return locus_sample_counts

# --- Execution ---
counts = count_samples_per_locus(VCF_FILE_PATH)

if counts:
    print("\n--- Samples per Locus (CHROM:POS) ---")
    
    # Determine the width for alignment
    max_locus_len = max(len(locus) for locus in counts.keys())
    
    # Sort by chromosome and position
    sorted_loci = sorted(counts.items(), key=lambda item: (item[0].split(':')[0], int(item[0].split(':')[1])))
    
    # Print formatted output
    for locus, count in sorted_loci:
        print(f"{locus.ljust(max_locus_len)}: {count} samples")

    print(f"\nTotal unique loci processed: {len(counts)}")
