SHELL := /bin/bash

.PHONY: setup download metadata de report all

setup:
	Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org'); if (file.exists('renv.lock')) renv::restore(prompt = FALSE) else renv::init(bare = TRUE)"

download:
	bash scripts/00_download_geo.sh

metadata:
	Rscript scripts/01_make_metadata.R

de:
	Rscript scripts/02_run_de_limma.R

report:
	Rscript scripts/03_make_report.R

pathways: results/tables/pathways_hallmark.csv results/figures/pathways_hallmark.png

results/tables/pathways_hallmark.csv results/figures/pathways_hallmark.png: scripts/04_pathways.R results/tables/de_limma_results.csv
	Rscript scripts/04_pathways.R

all: setup download metadata de pathways report
