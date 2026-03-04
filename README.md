# GSE224867 DLBCL Doxorubicin RNA-seq

Reproducible analysis scaffold for processing and differential expression analysis of GEO series **GSE224867** in DLBCL under doxorubicin treatment.

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
```

## Make targets

- `setup`: initialize/restore `renv`
- `download`: download GEO raw data
- `metadata`: generate sample metadata
- `de`: run differential expression step
- `report`: render report
- `all`: run full pipeline
