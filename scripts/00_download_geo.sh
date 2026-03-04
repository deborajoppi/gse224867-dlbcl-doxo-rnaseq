#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224867/suppl/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"
OUT="data/raw/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"

curl -L "$URL" -o "$OUT"
echo "Saved: $OUT"
