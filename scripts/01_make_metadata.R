#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

raw_dir <- "data/raw"
out_file <- "data/metadata/samples.csv"

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

candidates <- list.files(raw_dir, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
if (length(candidates) == 0) {
  stop("No TSV/TSV.GZ files found in data/raw/. Run `make download` or add a file manually.")
}

input_file <- candidates[[1]]
con <- if (grepl("\\.gz$", input_file)) gzfile(input_file, open = "rt") else file(input_file, open = "rt")
on.exit(close(con), add = TRUE)

header <- readLines(con, n = 1)
if (length(header) == 0) stop("Input file appears empty: ", input_file)

cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
if (length(cols) <= 1) stop("Expected tab-delimited matrix with sample columns in: ", input_file)

metadata <- data.frame(sample_id = cols[-1], stringsAsFactors = FALSE)
write.csv(metadata, out_file, row.names = FALSE)

message("Wrote ", nrow(metadata), " samples to ", out_file)
