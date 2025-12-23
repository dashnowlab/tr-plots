# tr-plots

Population-scale analysis and plotting of tandem repeat (TR) loci from 1kGP ONT cohorts. This repo assembles per-haplotype allele spreadsheets from STR VCFs, merges metadata and demographics, and produces interactive HTML/PNG plots and summary tables.

**Contents**
- Data in [data](data): sequencing VCFs and metadata.
- Source in [src/trplots](src/trplots): ETL + plotting scripts.
- Outputs in [results](results): per-script subfolders for HTML/PNG/CSV/XLSX.

**Quick Start**
- macOS + Python 3.10–3.12 recommended.
- Create a virtual env and install dependencies.
- Generate the master allele spreadsheet, then run plot generators.

---

**Setup**
- **Create env**:
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```
- **Configure paths**: Edit [src/trplots/config.py](src/trplots/config.py) to point at your desired inputs/outputs.
	- **Key variables**: `VCF_PATH`, `JSON_PATH`, `SEQ_DATA`, `OTHER_DATA`, `OUTPUT_BASE`.
- **Data layout** (provided in this repo):
	- VCFs: [data/sequencing_data/81_loci_502_samples](data/sequencing_data/81_loci_502_samples) and [data/sequencing_data/83_loci_503_samples](data/sequencing_data/83_loci_503_samples)
	- Metadata: [data/other_data](data/other_data) including `strchive-loci.json`, sample summaries and 1kGP sample info.

---

**Workflow Overview**
- **1. Build master allele spreadsheet**: Parses STR VCF, merges loci metadata and sample info, computes allele lengths and repeat counts.
	- Script: [src/trplots/create_master_allele_spreadsheet/allele_spreadsheet.py](src/trplots/create_master_allele_spreadsheet/allele_spreadsheet.py)
	- Inputs: VCF (`--vcf`), loci JSON (`--loci-json`), sample summary CSV (`--sample-csv`), 1kGP sample-info TXT (`--kgp-sample-info`).
	- Output: Excel at `data/other_data/allele_spreadsheet.xlsx` (default, override with `--output`).

**2. Generate plots**: Creates interactive HTML and PNG plots under `results/plots/<project>`.
	- Ancestry bars: [src/trplots/main_projects/ancestry_bar_plots/ancestry_plot_generator.py](src/trplots/main_projects/ancestry_bar_plots/ancestry_plot_generator.py)
	- Tandem repeat distributions: [src/trplots/main_projects/tandem_repeat_bar_plots/path_ref_motif_tandem_repeat_plot_generator.py](src/trplots/main_projects/tandem_repeat_bar_plots/path_ref_motif_tandem_repeat_plot_generator.py)
	- Allele length box/violin plots:
		- [src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_boxplots.py](src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_boxplots.py)
		- [src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_violin_swarm.py](src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_violin_swarm.py)

- **3. Summaries and utilities**:
	- Percentage pathogenic by gene+disease: [src/trplots/percentage_path_alleles/percentage_path_alleles.py](src/trplots/percentage_path_alleles/percentage_path_alleles.py)
	- Count affected/carriers per locus: [src/trplots/count_affected_and_carriers_per_locus/count_affected_and_carriers_per_locus.py](src/trplots/count_affected_and_carriers_per_locus/count_affected_and_carriers_per_locus.py)
	- Match VCF loci to JSON + ancestry: [src/trplots/matching_files/vcf_loci_matching_json_metadata_csv_ancestry.py](src/trplots/matching_files/vcf_loci_matching_json_metadata_csv_ancestry.py)
	- Merge many plots:
		- HTML → one file: [src/trplots/merging_files/merge_html_plots_to_one.py](src/trplots/merging_files/merge_html_plots_to_one.py)
		- PNG → single HTML: [src/trplots/merging_files/merge_png_plots_to_one.py](src/trplots/merging_files/merge_png_plots_to_one.py)

---

**Run Instructions**
- From the repo root, run scripts directly:

1) Build master allele spreadsheet
```bash
python src/trplots/create_master_allele_spreadsheet/allele_spreadsheet.py \
	--vcf data/sequencing_data/81_loci_502_samples/1000g-ont-strchive-81_loci_502_samples_81224_alleles.vcf.gz \
	--loci-json data/other_data/strchive-loci.json \
	--sample-csv data/other_data/1kgp_ont_500_summary_-_sheet1.csv \
	--kgp-sample-info data/other_data/1000_genomes_20130606_sample_info.txt \
	--output data/other_data/allele_spreadsheet.xlsx
```

2) Ancestry bar plots (respects `--test`, `--test-limit`, `--save-test-outputs`, `--output-dir`)
```bash
python src/trplots/main_projects/ancestry_bar_plots/ancestry_plot_generator.py --test --test-limit 10
# Optional: override output base
python src/trplots/main_projects/ancestry_bar_plots/ancestry_plot_generator.py --output-dir results/plots/ancestry_plots
```

3) Tandem repeat plots (consumes allele spreadsheet; respects `--test`, `--test-limit`, `--output-dir`)
```bash
python src/trplots/main_projects/tandem_repeat_bar_plots/path_ref_motif_tandem_repeat_plot_generator.py --test --test-limit 10 \
	--output-dir results/plots/path_ref_motif_tandem_repeats_plots
```

4) Allele length box/violin plots
```bash
python src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_boxplots.py --test --test-limit 10
python src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_violin_swarm.py --test --test-limit 10
```

5) Matching VCF loci to JSON + ancestry (uses defaults from `trplots.config`)
```bash
python src/trplots/matching_files/vcf_loci_matching_json_metadata_csv_ancestry.py --test --test-limit 100
```

6) Percentage pathogenic alleles (by Gene+Disease)
```bash
python src/trplots/percentage_path_alleles/percentage_path_alleles.py \
	--input data/other_data/allele_spreadsheet.xlsx \
	--output data/other_data/percentage_pathogenic_by_gene_disease.xlsx
```

7) Merge plot outputs
```bash
# Merge HTML plots produced by ancestry/tandem-repeat generators
python src/trplots/merging_files/merge_html_plots_to_one.py --input-dir results/plots/ancestry_plot_generator/HTML --parts 5

# Combine PNGs into a single HTML
python src/trplots/merging_files/merge_png_plots_to_one.py --input-dir results/plots/ancestry_plot_generator/PNG
```

---

**Script Behaviors (Overview)**
- **allele_spreadsheet.py**: Builds `allele_spreadsheet.xlsx` by integrating VCF haplotypes, loci metadata, and demographics; computes repeat counts and QC flags.
- **ancestry_plot_generator.py**: Produces ancestry bar charts; outputs to `results/plots/ancestry_plots/{html,png}` (test mode → `test_outputs`).
- **tandem_repeat_bar_plots/** (`path_ref_motif_tandem_repeat_plot_generator.py`): Generates histograms and repeat-count distributions per locus/motif; outputs to `results/plots/path_ref_motif_tandem_repeats_plots/{html,png}` (test mode → `test_outputs`).
- **allele_length_boxplots_and_violin/**: Creates box/violin + swarm plots of allele lengths; outputs to `results/plots/{allele_length_boxplots,allele_length_violin_swarm_plots}/{html,png}` (test mode → `test_outputs`).
- **motif_count/**: Computes motif-aware metrics from STR VCF using pathogenic motifs from JSON; supports `--vcf` and related inputs (see script help).
- **matching_files/**: Aligns VCF loci to JSON metadata and ancestry, writing intermediate CSVs and debug info.
- **merging_files/**: Utilities to merge many plot files into single HTML documents.
- **percentage_path_alleles.py**: Summarizes percentage of pathogenic alleles per Gene+Disease from the master spreadsheet.
- **count_affected_and_carriers_per_locus.py**: Tallies affected/carrier counts per locus from integrated data.

---

**Tips & Notes**
- Many scripts support quick `--test` runs with `--test-limit N` to validate setup before full processing.
- Default input/output locations are configurable via [src/trplots/config.py](src/trplots/config.py).
- Outputs are organized under [results](results) at `results/plots/<project>/{html,png}`; test runs save under `test_outputs`.
- If `pysam` fails to install on macOS, ensure Xcode Command Line Tools are installed (`xcode-select --install`) and try reinstalling `pysam` in the virtual env.

---

**Development**
- Format/lint as needed; scripts are standalone and use `argparse` for CLI.
- To discover flags for any script: `python <script.py> --help`.
