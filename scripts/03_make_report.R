#!/usr/bin/env Rscript
stopifnot(requireNamespace("rmarkdown", quietly = TRUE))

rmarkdown::render(
  input = "reports/report.Rmd",
  output_file = "report.html",
  output_dir = "results",
  knit_root_dir = normalizePath(".", winslash = "/", mustWork = TRUE),
  quiet = FALSE
)

message("Wrote results/report.html")
