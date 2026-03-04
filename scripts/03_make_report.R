#!/usr/bin/env Rscript
stopifnot(requireNamespace("rmarkdown", quietly = TRUE))
rmarkdown::render(
  input = "reports/report.Rmd",
  output_file = "report.html",
  output_dir = "results",
  quiet = TRUE
)
message("Wrote results/report.html")
