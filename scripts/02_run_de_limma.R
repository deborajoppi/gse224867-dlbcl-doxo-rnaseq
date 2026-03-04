#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Placeholder DE table scaffold; replace with limma workflow.
de_table <- data.frame(
  gene = character(),
  logFC = numeric(),
  P.Value = numeric(),
  adj.P.Val = numeric(),
  stringsAsFactors = FALSE
)

write.csv(de_table, "results/tables/de_limma_results.csv", row.names = FALSE)
message("Created placeholder results/tables/de_limma_results.csv")
