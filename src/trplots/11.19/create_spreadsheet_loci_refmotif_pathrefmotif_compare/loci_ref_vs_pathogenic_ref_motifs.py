import json
import pandas as pd
from pathlib import Path
from trplots.config import OTHER_DATA

# Path to your JSON file
JSON_PATH = "/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive-loci.json"

# Helper to normalize motif fields (string or list) into a single string
def normalize_motifs(value):
    """
    Accepts a string, list of strings, or None.
    Returns a single string with all motifs (if multiple) joined by commas.
    """
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        # keep order, drop empty/None, strip whitespace
        cleaned = [str(m).strip() for m in value if m]
        return ", ".join(cleaned)
    # Fallback for weird types
    return str(value)

def main():
    # Load JSON
    with open(JSON_PATH, "r") as f:
        loci = json.load(f)

    rows = []
    for entry in loci:
        chrom = entry.get("chrom")
        pos_hg38 = entry.get("start_hg38")
        gene = entry.get("gene")
        disease = entry.get("disease")
        
        ref_motif = normalize_motifs(
            entry.get("reference_motif_reference_orientation")
        )
        path_ref_motif = normalize_motifs(
            entry.get("pathogenic_motif_reference_orientation")
        )

        rows.append({
            "Chrom": chrom,
            "Position_hg38": pos_hg38,
            "Gene": gene,
            "Disease": disease,
            "Ref Motif (ref orientation)": ref_motif,
            "Pathogenic Motif (ref orientation)": path_ref_motif,
        })

    df = pd.DataFrame(rows)

    # Show a preview
    print(df.head())

    # Save to CSV
    output_excel = Path(OTHER_DATA) / "loci_ref_vs_pathogenic_ref_motifs.xlsx"
    output_excel.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(output_excel, index=False)
    print(f"\nWrote {len(df)} rows to {output_excel}")

if __name__ == "__main__":
    main()
