import json
import pandas as pd

JSON_PATH = "/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive-loci.json"
OUTPUT_XLSX = "/Users/annelisethorn/Documents/github/tr-plots/data/other_data/strchive_loci_full_table.xlsx"

def flatten(entry, parent_key="", sep="."):
    """
    Recursively flattens a nested dict or list.
    Lists become comma-separated strings.
    """
    items = {}
    for key, value in entry.items():
        new_key = f"{parent_key}{sep}{key}" if parent_key else key

        if isinstance(value, dict):
            # recurse into dictionaries
            items.update(flatten(value, new_key, sep=sep))

        elif isinstance(value, list):
            # lists become comma-separated, converting inner values too
            cleaned_list = []
            for v in value:
                if isinstance(v, dict):
                    cleaned_list.append(
                        "; ".join(f"{kk}:{vv}" for kk, vv in v.items())
                    )
                else:
                    cleaned_list.append(str(v))
            items[new_key] = ", ".join(cleaned_list)

        else:
            # normal value
            items[new_key] = value

    return items

def main():
    # Load JSON
    with open(JSON_PATH, "r") as f:
        data = json.load(f)

    rows = []
    for entry in data:
        flat = flatten(entry)
        rows.append(flat)

    df = pd.DataFrame(rows)

    # Write Excel
    df.to_excel(OUTPUT_XLSX, index=False)

    print(f"✓ Wrote {len(df)} loci with ALL fields to {OUTPUT_XLSX}")
    print("✓ Columns:", len(df.columns))

if __name__ == "__main__":
    main()
