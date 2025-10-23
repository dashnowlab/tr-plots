import gzip
import os
import re

# Use the VCF file path from your context
VCF_FILE_PATH = '/Users/annelisethorn/Documents/github/tr-plots/data/sequencing_data/83_loci_503_samples/1000g-ont-strchive-83_loci_503_samples.vcf.gz'

def count_total_alleles(file_path):
    """
    Reads a VCF file and counts the total number of non-missing allele calls 
    across all samples and all loci.
    """
    
    print(f"Reading VCF file: {os.path.basename(file_path)}...")
    
    # Determine file opening function
    is_gzipped = file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    total_allele_calls = 0
    num_loci = 0
    vcf_header_found = False
    
    try:
        with open_func(file_path, mode) as f:
            for line in f:
                # 1. Skip metadata lines
                if line.startswith('##'):
                    continue
                
                # 2. Identify the header line
                if line.startswith('#CHROM'):
                    vcf_header_found = True
                    continue
                
                # 3. Process data lines (locus records)
                if vcf_header_found:
                    parts = line.strip().split('\t')
                    if len(parts) < 10: continue
                        
                    format_keys = parts[8].split(':')
                    
                    # Find the index of the Genotype (GT) field
                    try:
                        gt_index = format_keys.index('GT')
                    except ValueError:
                        # If GT is missing, we cannot count alleles for this locus
                        continue
                        
                    num_loci += 1
                    
                    # Genotype data starts at index 9
                    genotype_data = parts[9:]
                    
                    for sample_genotype_string in genotype_data:
                        sample_fields = sample_genotype_string.split(':')
                        
                        # Ensure we have enough fields to look up GT
                        if len(sample_fields) > gt_index:
                            genotype = sample_fields[gt_index]
                            
                            # Split genotype by separator ('/' or '|')
                            gt_indices = re.split(r'[/|]', genotype)
                            
                            # Count each called allele
                            for allele_index in gt_indices:
                                # '.' indicates a missing allele (e.g., in a diploid call like '0/.' or a haploid call like '.')
                                if allele_index != '.':
                                    total_allele_calls += 1
                                # NOTE: We do not check for '0' (reference) or '1+' (alternate) here, 
                                # we only check if the call itself is NOT missing.

        print(f"Total loci processed: {num_loci}")
        return total_allele_calls

    except FileNotFoundError:
        print(f"Error: The VCF file was not found at the expected path: {file_path}")
        return 0
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return 0

# --- Execution ---
total_count = count_total_alleles(VCF_FILE_PATH)

if total_count > 0:
    print(f"\nTotal number of called alleles in the VCF file: {total_count}")
