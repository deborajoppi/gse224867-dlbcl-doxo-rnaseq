cat > scripts/01_make_metadata.R <<'EOF'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

infile <- "data/raw/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"
outfile <- "data/metadata/samples.csv"

if (!file.exists(infile)) {
  stop("Input file not found: ", infile, "\nRun: make download", call. = FALSE)
}

# Read only header line (tab-separated)
con <- gzfile(infile, open = "rt")
header <- readLines(con, n = 1)
close(con)

cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
if (length(cols) < 3) stop("Unexpected file format: too few columns in header.", call. = FALSE)

sample_ids <- cols[-1]  # first column is gene id
m <- regexec("^(.*)_(UT|ADR)$", sample_ids)
parts <- regmatches(sample_ids, m)

ok <- lengths(parts) == 3
if (!all(ok)) {
  bad <- sample_ids[!ok]
  stop("Some sample names do not match 'CELL_UT'/'CELL_ADR' pattern:\n- ",
       paste(bad, collapse = "\n- "), call. = FALSE)
}

cell_line <- vapply(parts, `[[`, character(1), 2)
treatment <- vapply(parts, `[[`, character(1), 3)

df <- data.frame(
  sample_id = sample_ids,
  cell_line = factor(cell_line),
  treatment = factor(treatment, levels = c("UT", "ADR")),
  stringsAsFactors = FALSE
)

dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
write.csv(df, outfile, row.names = FALSE)
message("Wrote: ", outfile, " (n=", nrow(df), ")")
EOF

chmod +x scripts/01_make_metadata.R
