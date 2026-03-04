#!/usr/bin/env Rscript
infile <- "data/raw/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"
outfile <- "data/metadata/samples.csv"

if (!file.exists(infile)) stop("Missing input: ", infile, "\nRun: make download", call. = FALSE)

con <- gzfile(infile, "rt")
header <- readLines(con, n = 1)
close(con)

cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
sample_ids <- cols[-1]  # first column is gene id

m <- regexec("^(.*)_(UT|ADR)$", sample_ids)
parts <- regmatches(sample_ids, m)
ok <- lengths(parts) == 3
if (!all(ok)) {
  bad <- sample_ids[!ok]
  stop("Bad sample names (expected CELL_UT or CELL_ADR):\n- ", paste(bad, collapse = "\n- "), call. = FALSE)
}

cell_line <- vapply(parts, `[[`, character(1), 2)
treatment <- vapply(parts, `[[`, character(1), 3)

df <- data.frame(
  sample_id = sample_ids,
  cell_line = factor(cell_line),
  treatment = factor(treatment, levels = c("UT", "ADR")),
  stringsAsFactors = FALSE
)

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
write.csv(df, outfile, row.names = FALSE)
message("Wrote ", nrow(df), " samples to ", outfile)
