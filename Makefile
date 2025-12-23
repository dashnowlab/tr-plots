.PHONY: help test-all test-ancestry test-boxplots test-violin test-tr full-all full-ancestry full-boxplots full-violin full-tr

# Configuration
PY ?= python
TEST_LIMIT ?= 2

# Script paths
ANCESTRY := src/trplots/main_projects/ancestry_bar_plots/ancestry_plot_generator.py
BOXPLOTS := src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_boxplots.py
VIOLIN   := src/trplots/main_projects/allele_length_boxplots_and_violin/allele_length_violin_swarm.py
TR_PATHREF := src/trplots/main_projects/tandem_repeat_bar_plots/path_ref_motif_tandem_repeat_plot_generator.py

help:
	@echo "Targets:"
	@echo "  test-all           Run all plot generators in test mode (TEST_LIMIT=$(TEST_LIMIT))"
	@echo "  test-ancestry      Run ancestry plots in test mode"
	@echo "  test-boxplots      Run boxplots in test mode"
	@echo "  test-violin        Run violin+swarm in test mode"
	@echo "  test-tr            Run tandem repeat (path_ref) in test mode"
	@echo "  full-all           Run all plot generators (full)"
	@echo "  full-ancestry      Run ancestry plots (full)"
	@echo "  full-boxplots      Run boxplots (full)"
	@echo "  full-violin        Run violin+swarm (full)"
	@echo "  full-tr            Run tandem repeat (path_ref) (full)"

test-all: test-ancestry test-boxplots test-violin test-tr

test-ancestry:
	$(PY) $(ANCESTRY) --test --test-limit $(TEST_LIMIT)

test-boxplots:
	$(PY) $(BOXPLOTS) --test --test-limit $(TEST_LIMIT)

test-violin:
	$(PY) $(VIOLIN) --test --test-limit $(TEST_LIMIT)

test-tr:
	$(PY) $(TR_PATHREF) --test --test-limit $(TEST_LIMIT)

full-all: full-ancestry full-boxplots full-violin full-tr

full-ancestry:
	$(PY) $(ANCESTRY)

full-boxplots:
	$(PY) $(BOXPLOTS)

full-violin:
	$(PY) $(VIOLIN)

full-tr:
	$(PY) $(TR_PATHREF)
