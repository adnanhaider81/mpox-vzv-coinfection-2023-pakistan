SHELL := /bin/bash

ENV_NAME := mpox-vzv-env
THREADS ?= 4

.PHONY: env run dry test clean help

help:
	@echo "make env   - create conda env"
	@echo "make run   - run full Snakemake workflow"
	@echo "make dry   - dry run to show planned steps"
	@echo "make test  - quick sanity check on plotting"
	@echo "make clean - remove work and results"

env:
	conda env create -f env/environment.yml || echo "Env may already exist"
	@echo "Activate with: conda activate $(ENV_NAME)"

run:
	@if [ -z "$$NCBI_EMAIL" ]; then echo "Set NCBI_EMAIL before running"; exit 1; fi
	snakemake -s workflow/Snakefile -c $(THREADS) --printshellcmds

dry:
	snakemake -s workflow/Snakefile -n -c $(THREADS)

test:
	python -m pip install -r env/requirements.txt
	python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
	@echo "Wrote results-example/example_plot.png"

clean:
	rm -rf work results results-example logs .snakemake
