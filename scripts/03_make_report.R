#!/usr/bin/env Rscript
stopifnot(requireNamespace("rmarkdown", quietly = TRUE))

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
res  <- normalizePath(file.path(root, "results"), winslash = "/", mustWork = TRUE)
inp  <- normalizePath(file.path(root, "reports", "report.Rmd"), winslash = "/", mustWork = TRUE)

rmarkdown::render(
  input = inp,
  output_file = "report.html",
  output_dir = res,
  knit_root_dir = res,
  quiet = FALSE
)

message("Wrote results/report.html")
