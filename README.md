# GSE224867 DLBCL Doxorubicin RNA-seq

Reproducible analysis scaffold for processing and differential expression analysis of GEO series **GSE224867** in DLBCL under doxorubicin treatment.

## Question & design

This project identifies a core transcriptional response to **doxorubicin (ADR)** across **18 DLBCL cell lines**, using a paired design (**UT vs ADR** within each cell line) and linear modeling that blocks for cell line.

## Project layout

- `data/raw/`: downloaded raw GEO files (ignored by git)
- `data/metadata/samples.csv`: sample metadata generated from expression matrix column names
- `scripts/`: pipeline scripts
- `reports/report.Rmd`: analysis report source
- `results/`: output tables and figures

## Quick start

```bash
make setup
make all
open results/report.html
```

## Outputs

After `make all`, you get:

- `results/tables/de_limma_results.csv` — differential expression results (ADR vs UT)
- `results/figures/pca.png` — PCA of the top variable genes
- `results/figures/volcano.png` — volcano plot
- `results/report.html` — HTML report (open it in the browser)

### Snapshot

## Methods (high level)

Expression values are loaded from the GEO supplementary matrix (`norm_counts`). Differential expression is modeled with **limma** using the design `~ cell_line + treatment` (UT as reference). If values are integer-like, the pipeline can switch to **voom+limma**.

## Make targets

- `setup`: initialize/restore `renv`
- `download`: download GEO raw data
- `metadata`: generate sample metadata
- `de`: run differential expression step
- `report`: render report
- `all`: run full pipeline

## Limitations / next steps

- The GEO file used here is labeled `norm_counts` (not raw counts). This repo demonstrates a clean, reproducible workflow and paired modeling; a future extension would quantify from FASTQ/SRA to obtain raw counts and run DESeq2/edgeR end-to-end.
