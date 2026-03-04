#!/usr/bin/env Rscript

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required. Install it with install.packages('rmarkdown').")
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  input = "reports/report.Rmd",
  output_file = "report.html",
  output_dir = "results",
  quiet = TRUE
)

message("Rendered results/report.html")
