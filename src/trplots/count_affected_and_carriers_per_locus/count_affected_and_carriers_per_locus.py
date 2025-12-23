"""
---------------------------------------------
 Script: Calculate Affected & Carrier Counts per Locus
 Purpose:
   - Loads allele spreadsheet with repeat counts and pathogenic thresholds
   - Flags alleles as pathogenic based on gene-specific rules
   - Determines if samples are AFFECTED based on inheritance models
   - Counts carriers and affected samples per locus
   - Saves summary table as Excel
---------------------------------------------
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================================
# CONFIGURATION & FILE PATHS
# ============================================================================

# --- Load allele spreadsheet ---
# This spreadsheet contains per-allele data including repeat counts, genes,
# diseases, inheritance models, and pathogenic thresholds.
ALLELE_SPREADSHEET_PATH = "/Users/annelisethorn/Documents/github/tr-plots/data/other_data/allele_spreadsheet.xlsx"
OUTPUT_PATH = Path(ALLELE_SPREADSHEET_PATH).parent / "affected_carrier_count_per_loci.xlsx"

df = pd.read_excel(ALLELE_SPREADSHEET_PATH, engine="openpyxl")

# ============================================================================
# COLUMN MAPPING
# ============================================================================

# --- Normalize column names (case-insensitive lookup) ---
# Different spreadsheets may use varying column name conventions (e.g., "Gene" vs "gene").
# This helper function allows flexible column lookup regardless of case or minor naming differences.
cols_lower = {c.lower(): c for c in df.columns}

def col(*names):
    """
    Return the actual column name from the spreadsheet for the first matching name (case-insensitive).
    If none of the provided names match, return None.
    """
    for n in names:
        if n.lower() in cols_lower:
            return cols_lower[n.lower()]
    return None

# Map required columns to standardized names
chrom_col          = col("Chromosome", "Chr", "chrom")
pos_col            = col("Position (Start)", "Position", "Pos", "position")
gene_col           = col("Gene")
disease_col        = col("Disease")
sample_id_col      = col("Sample ID", "Sample_ID", "sample_id")
sample_id_cleaned_col = col("Sample ID Cleaned", "Sample ID (Cleaned)", "Sample_ID_Cleaned", "sample_id_cleaned")
repeat_col         = col("Repeat Count", "RepeatCount", "repeat_count")
pathogenic_min_col = col("Pathogenic min", "Pathogenic_min", "pathogenic_min")
pathogenic_max_col = col("Pathogenic max", "Pathogenic_max", "pathogenic_max")
inheritance_col    = col("Inheritance")
sex_col            = col("Sex")

# Ensure all critical columns are present
if not all([gene_col, disease_col, repeat_col, pathogenic_min_col, pathogenic_max_col, sample_id_col]):
    raise SystemExit(
        "Spreadsheet missing required columns. Found: "
        f"{list(df.columns)}"
    )

# ============================================================================
# COLUMN RENAMING & STANDARDIZATION
# ============================================================================

# --- Rename columns for consistency ---
# Standardize column names to a uniform format for downstream processing.
df = df.rename(columns={
    gene_col: "Gene",
    disease_col: "Disease",
    sample_id_col: "Sample ID",
    repeat_col: "Repeat count",
    pathogenic_min_col: "Pathogenic min",
    pathogenic_max_col: "Pathogenic max"
})

# Use cleaned sample ID if available, otherwise use raw sample ID
if sample_id_cleaned_col:
    df = df.rename(columns={sample_id_cleaned_col: "Sample ID (Cleaned)"})
    # Use cleaned ID for grouping
    df["Sample ID for Analysis"] = df["Sample ID (Cleaned)"].fillna(df["Sample ID"])
else:
    df["Sample ID for Analysis"] = df["Sample ID"]

# Handle inheritance column (default to "Unknown" if missing)
if inheritance_col:
    df = df.rename(columns={inheritance_col: "Inheritance"})
else:
    df["Inheritance"] = "Unknown"

# Handle sex column (default to "Unknown" if missing)
if sex_col:
    df = df.rename(columns={sex_col: "Sex"})
else:
    df["Sex"] = "Unknown"

# --- Define grouping columns for locus identification ---
# If chromosome and position are available, group by (Chromosome, Position, Gene, Disease).
# Otherwise, group by (Gene, Disease).
if chrom_col and pos_col:
    df = df.rename(columns={chrom_col: "Chromosome", pos_col: "Position"})
    group_cols = ["Chromosome", "Position", "Gene", "Disease"]
else:
    group_cols = ["Gene", "Disease"]

# ============================================================================
# SEX STANDARDIZATION
# ============================================================================

# --- Standardize sex labels ---
# Map various representations of sex to standardized labels: "Male", "Female", "Unknown".
# This ensures consistent handling of sex-linked inheritance models (e.g., XR).
sex_map = {
    "m": "Male", "male": "Male", "M": "Male", "MALE": "Male",
    "f": "Female", "female": "Female", "F": "Female", "FEMALE": "Female",
    "xx": "Female", "XX": "Female",
    "xy": "Male",   "XY": "Male",
    "u": "Unknown", "unknown": "Unknown", "U": "Unknown",
    "UNK": "Unknown", "Unk": "Unknown"
}

df["Sex"] = (
    df["Sex"]
    .astype(str)
    .str.strip()
    .map(lambda s: sex_map.get(s, s))
)
df.loc[df["Sex"].isin(["", "nan", "None"]), "Sex"] = "Unknown"

# ============================================================================
# ALLELE-LEVEL PATHOGENIC FLAGGING
# ============================================================================

# --- Flag pathogenic alleles (allele-level) ---
def is_allele_pathogenic(row):
    """
    Determine if a single allele is pathogenic based on repeat count and gene-specific rules.
    
    Special handling for VWA1:
        - VWA1 has discrete pathogenic repeat counts (e.g., 1 and 3).
        - Pathogenic if repeat count == pathogenic_min OR == pathogenic_max (not a range).
    
    For all other genes:
        - Pathogenic if repeat count is within the inclusive range [pathogenic_min, pathogenic_max].
    
    Returns:
        bool: True if allele is pathogenic, False otherwise.
    """
    gene = row.get("Gene", "")
    repeat_count = row.get("Repeat count")
    path_min = row.get("Pathogenic min")
    path_max = row.get("Pathogenic max")

    # Skip if any required data is missing
    if not all([pd.notna(repeat_count), pd.notna(path_min), pd.notna(path_max)]):
        return False

    if str(gene).strip().upper() == "VWA1":
        # VWA1: pathogenic if repeat count equals pathogenic_min OR pathogenic_max
        return repeat_count == path_min or repeat_count == path_max
    else:
        # Other genes: pathogenic if within range [pathogenic_min, pathogenic_max]
        return path_min <= repeat_count <= path_max

# Apply pathogenic flagging to each allele (row) in the dataframe
df["Is pathogenic"] = df.apply(is_allele_pathogenic, axis=1)

# ============================================================================
# SAMPLE-LEVEL AFFECTED STATUS (INHERITANCE MODEL LOGIC)
# ============================================================================

# --- Determine if sample is AFFECTED based on inheritance model ---
def is_affected_by_inheritance(sample_df: pd.DataFrame) -> bool:
    """
    Determine if a sample is AFFECTED (i.e., clinically pathogenic) based on the inheritance model
    and the number of pathogenic alleles the sample carries.
    
    Inheritance Models:
        - AD (Autosomal Dominant):        Affected if >= 1 pathogenic allele
        - AR (Autosomal Recessive):       Affected if >= 2 pathogenic alleles
        - XD (X-linked Dominant):         Affected if >= 1 pathogenic allele (any sex)
        - XR (X-linked Recessive):        Affected if:
                                            - Male:   >= 1 pathogenic allele
                                            - Female: >= 2 pathogenic alleles
                                            - Unknown sex: >= 2 pathogenic alleles (conservative)
        - Unknown / other:                Default to AD (>= 1 pathogenic allele)
    
    Args:
        sample_df: DataFrame containing all allele rows for a single sample at a single locus.
    
    Returns:
        bool: True if sample is affected, False otherwise.
    """
    if sample_df.empty:
        return False

    # Extract inheritance model, pathogenic allele count, and sex for this sample
    inheritance = sample_df["Inheritance"].iloc[0] if "Inheritance" in sample_df.columns else "Unknown"
    pathogenic_count = sample_df["Is pathogenic"].sum()
    sex = sample_df["Sex"].iloc[0] if "Sex" in sample_df.columns else "Unknown"

    # Normalize inheritance and sex strings
    inheritance = str(inheritance).strip().upper() if pd.notna(inheritance) else "Unknown"
    sex = str(sex).strip()

    # Apply inheritance model logic
    if inheritance in ["AD", "XD"]:
        # Autosomal Dominant or X-linked Dominant: >=1 pathogenic allele
        return pathogenic_count >= 1

    elif inheritance == "AR":
        # Autosomal Recessive: >=2 pathogenic alleles
        return pathogenic_count >= 2

    elif inheritance == "XR":
        # X-linked Recessive: sex-dependent threshold
        if sex == "Male":
            return pathogenic_count >= 1
        elif sex == "Female":
            return pathogenic_count >= 2
        else:
            # Unknown sex: conservative choice, require >=2
            return pathogenic_count >= 2

    else:
        # Default: treat as AD (>=1 pathogenic allele)
        return pathogenic_count >= 1

# ============================================================================
# PER-LOCUS ANALYSIS: CARRIERS & AFFECTED COUNTS
# ============================================================================

# --- Count carriers & affected samples per locus ---
# A "carrier" is any sample with >= 1 pathogenic allele.
# An "affected" sample meets the inheritance model's clinical threshold (e.g., >=2 for AR).
results = []

for locus_key, locus_df in df.groupby(group_cols):
    # Per-sample classification within this locus (using cleaned sample ID)
    sample_info = {}

    for sample_id, sample_group in locus_df.groupby("Sample ID for Analysis"):
        pathogenic_alleles_count = int(sample_group["Is pathogenic"].sum())
        is_carrier = pathogenic_alleles_count >= 1  # Carrier: has at least 1 pathogenic allele
        is_affected = is_affected_by_inheritance(sample_group)  # Affected: meets inheritance threshold
        sex = sample_group["Sex"].iloc[0] if "Sex" in sample_group.columns else "Unknown"

        # Store per-sample info for later printout
        sample_info[sample_id] = {
            "carrier": is_carrier,
            "affected": is_affected,
            "pathogenic_alleles": pathogenic_alleles_count,
            "sex": sex,
            "pathogenic_repeats": list(sample_group.loc[sample_group["Is pathogenic"], "Repeat count"].values)
        }

    # Aggregate counts for this locus
    total_samples = len(sample_info)
    # Carriers: only those with >=1 pathogenic allele who are NOT affected
    carrier_samples = sum(1 for v in sample_info.values() if v["carrier"] and not v["affected"])
    affected_samples = sum(1 for v in sample_info.values() if v["affected"])

    # Get lists of carrier (non-affected) and affected sample IDs (using cleaned IDs)
    carrier_sample_ids = [s for s, v in sample_info.items() if v["carrier"] and not v["affected"]]
    affected_sample_ids = [s for s, v in sample_info.items() if v["affected"]]
    
    # Convert lists to comma-separated strings for spreadsheet
    carrier_ids_str = ", ".join(sorted(carrier_sample_ids)) if carrier_sample_ids else ""
    affected_ids_str = ", ".join(sorted(affected_sample_ids)) if affected_sample_ids else ""

    # Allele-level counts
    total_alleles = len(locus_df)
    pathogenic_alleles = int(locus_df["Is pathogenic"].sum())

    # Inheritance model for locus (assume consistent within a locus)
    inheritance = locus_df["Inheritance"].iloc[0] if "Inheritance" in locus_df.columns else "Unknown"
    # Ensure inheritance is never NaN, None, or empty string
    if pd.isna(inheritance) or str(inheritance).strip() == "":
        inheritance = "Unknown"

    # Build result row for this locus
    if isinstance(locus_key, tuple):
        row = dict(zip(group_cols, locus_key))
    else:
        row = {group_cols[0]: locus_key}

    row.update({
        "Inheritance": inheritance,
        "Total Samples": total_samples,
        "Total Carrier Samples": carrier_samples,
        "Percent Carrier Samples": np.round((carrier_samples / total_samples) * 100, 4),
        "Carrier Sample IDs": carrier_ids_str,
        "Total Affected Samples": affected_samples,
        "Percent Affected Samples": np.round((affected_samples / total_samples) * 100, 4),
        "Affected Sample IDs": affected_ids_str,
        "Total Alleles": total_alleles,
        "Total Pathogenic Alleles": pathogenic_alleles,
        "Percent Pathogenic Alleles": np.round((pathogenic_alleles / total_alleles) * 100, 4)
    })

    results.append(row)

# ============================================================================
# OUTPUT: SAVE SUMMARY TABLE AS EXCEL
# ============================================================================

# Create summary dataframe
result = pd.DataFrame(results)

# Sort by affected sample count (descending), then carrier count
result = result.sort_values(
    ["Total Affected Samples", "Total Carrier Samples"],
    ascending=[False, False]
)

# Save results to Excel only
result.to_excel(OUTPUT_PATH, index=False, engine="openpyxl")
print(f"Pathogenic counts per loci written to: {OUTPUT_PATH}")
# print(f"\nTotal loci analyzed: {len(result)}")
# print(f"Loci with >=1 affected sample: {len(result[result['Total_Affected_Samples'] > 0])}")
# print(f"Loci with >=1 carrier sample: {len(result[result['Total_Carrier_Samples'] > 0])}")