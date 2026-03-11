# GSE224867 DLBCL Doxorubicin RNA-seq
### Paired bulk RNA-seq re-analysis of doxorubicin response across 18 DLBCL cell lines

A reproducible bulk RNA-seq analysis pipeline for **GSE224867**, built around the published processed expression matrix and a paired differential expression design.

---

## Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Pipeline Architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Run](#run)
- [Results](#results)
- [Repository Structure](#repository-structure)
- [Limitations](#limitations)
- [Potential Improvements](#potential-improvements)

---

## Overview

This repo identifies a core transcriptional response to **doxorubicin (ADR)** in diffuse large B-cell lymphoma by comparing:

- `UT` vs `ADR`
- within the same cell line
- across **18 matched DLBCL cell lines** (`36` samples total)

The analysis is based on the GEO supplementary matrix:

- `GSE224867_norm_counts_DLBCL_TIS.tsv.gz`

and uses:

- paired modeling with `~ cell_line + treatment`
- `limma` on `log2(x + 1)` when the input is already normalized
- optional `voom + limma` fallback if the matrix is integer-like
- Hallmark pathway enrichment with `fgsea`

Current outputs include:

- PCA
- volcano plot
- Hallmark pathway bar plot
- differential expression table
- rendered HTML report

---

## Dataset

| Field | Value |
|---|---|
| GEO accession | [GSE224867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224867) |
| Disease context | Diffuse large B-cell lymphoma |
| Perturbation | Doxorubicin (`ADR`) |
| Control | Untreated (`UT`) |
| Input data used here | GEO processed `norm_counts` matrix |
| Samples | `36` |
| Matched cell lines | `18` |
| Design | Paired / blocked by cell line |

Treatment balance:

- `18` `UT`
- `18` `ADR`

---

## Pipeline Architecture

```text
GEO processed expression matrix
        |
        |  scripts/00_download_geo.sh
        v
Downloaded norm_counts matrix
        |
        |  scripts/01_make_metadata.R
        v
Sample metadata parsed from matrix columns
        |
        |  scripts/02_run_de_limma.R
        |  model: ~ cell_line + treatment
        v
Differential expression results + PCA + volcano plot
        |
        |  scripts/04_pathways.R
        v
Hallmark pathway enrichment (fgsea)
        |
        |  scripts/03_make_report.R
        v
HTML report
```

---

## Requirements

This project is R-based and uses `renv` for reproducibility.

Core requirements:

- R
- `renv`
- packages restored from `renv.lock`

Setup:

```bash
make setup
```

---

## Run

Full pipeline:

```bash
make all
```

Main make targets:

- `make setup` - restore the R environment with `renv`
- `make download` - download the GEO processed matrix
- `make metadata` - generate `data/metadata/samples.csv`
- `make de` - run paired differential expression with `limma`
- `make pathways` - run Hallmark enrichment with `fgsea`
- `make report` - render the HTML report
- `make all` - run the full workflow end to end

Open the final report:

```bash
open results/report.html
```

---

## Results

### Summary metrics

| Metric | Value |
|---|---|
| Samples analyzed | `36` |
| Matched cell lines | `18` |
| Genes tested | `37,581` |
| Genes with `FDR < 0.05` | `5,747` |
| Genes with `FDR < 0.05` and `|logFC| >= 1` | `934` |

### Main outputs

| Output | File |
|---|---|
| Differential expression table | `results/tables/de_limma_results.csv` |
| PCA plot | `results/figures/pca.png` |
| Volcano plot | `results/figures/volcano.png` |
| Hallmark pathway table | `results/tables/pathways_hallmark.csv` |
| Hallmark pathway plot | `results/figures/pathways_hallmark.png` |
| Volcano labels table | `results/tables/volcano_labels.csv` |
| HTML report | `results/report.html` |

### Pathway signal snapshot

The current Hallmark enrichment output shows strong pathway structure rather than a weak diffuse response. Top pathways include:

- negative enrichment of proliferation / growth programs such as:
  - `HALLMARK_MYC_TARGETS_V1`
  - `HALLMARK_E2F_TARGETS`
  - `HALLMARK_OXIDATIVE_PHOSPHORYLATION`
- positive enrichment of inflammatory and interferon programs such as:
  - `HALLMARK_INTERFERON_GAMMA_RESPONSE`
  - `HALLMARK_INTERFERON_ALPHA_RESPONSE`
  - `HALLMARK_INFLAMMATORY_RESPONSE`

This gives the README a cleaner story: **ADR suppresses proliferative programs while inducing stress / immune-like response programs** in this paired DLBCL panel.

---

## Repository Structure

```text
gse224867-dlbcl-doxo-rnaseq/
├── data/
│   ├── metadata/
│   │   └── samples.csv
│   └── raw/
│       └── GSE224867_norm_counts_DLBCL_TIS.tsv.gz
├── reports/
│   └── report.Rmd
├── results/
│   ├── figures/
│   │   ├── pca.png
│   │   ├── volcano.png
│   │   └── pathways_hallmark.png
│   ├── tables/
│   │   ├── de_limma_results.csv
│   │   ├── pathways_hallmark.csv
│   │   └── volcano_labels.csv
│   └── report.html
├── scripts/
│   ├── 00_download_geo.sh
│   ├── 01_make_metadata.R
│   ├── 02_run_de_limma.R
│   ├── 03_make_report.R
│   └── 04_pathways.R
├── Makefile
├── renv.lock
└── README.md
```

---

## Limitations

- This is **not** a raw FASTQ-to-counts pipeline.
- The starting point is GEO `norm_counts`, so this repo does not re-run alignment or raw-count quantification.
- Because the input matrix is processed rather than raw counts, `DESeq2` / count-based modeling is not the primary analysis path here.
- Gene labels in the DE output are still Ensembl-oriented; a cleaner gene-symbol annotation layer would improve interpretability.

---

## Potential Improvements

If you want this repo to feel more complete and publication-ready, I would prioritize these additions:

### README / presentation

- Add figure thumbnails directly into the README for `pca.png`, `volcano.png`, and `pathways_hallmark.png`.
- Add a short biological interpretation paragraph under the results section.
- Add exact citation information for the original GSE224867 study.

### Analysis outputs

- Add an **MA plot** to complement the volcano plot.
- Add a **sample-to-sample correlation heatmap** or clustered distance heatmap.
- Add a **paired line plot** for a handful of top genes to show within-cell-line ADR shifts directly.
- Add a **heatmap of top DE genes** across all 36 samples with treatment annotation.
- Add a **ranked GSEA running-score plot** for 2-3 top Hallmark pathways rather than only the NES bar chart.
- Add a **direction summary plot** showing counts of up- vs down-regulated genes at different thresholds.

### Interpretation / downstream biology

- Map Ensembl IDs to symbols in the final DE table.
- Add over-representation analysis for:
  - DNA damage response
  - apoptosis
  - p53 signaling
  - cell-cycle checkpoints
- Add a comparison of effect sizes across cell lines to identify:
  - shared ADR response genes
  - heterogeneous / context-specific responses

### Reproducibility

- Save a compact `sessionInfo()` or package-version file into `results/`.
- Add a small `config/` layer for thresholds such as:
  - FDR cutoff
  - volcano label count
  - pathway collection choice

