#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_pkgs <- c("fgsea", "msigdbr", "AnnotationDbi", "org.Hs.eg.db")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

in_file <- file.path("results", "tables", "de_limma_results.csv")
out_table <- file.path("results", "tables", "pathways_hallmark.csv")
out_plot <- file.path("results", "figures", "pathways_hallmark.png")

dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("Input DE file not found: ", in_file)
}

de <- utils::read.csv(in_file, check.names = FALSE)
if (nrow(de) == 0) {
  stop("Input DE file has no rows: ", in_file)
}

detect_symbol_column <- function(df) {
  preferred <- c("SYMBOL", "symbol", "gene_symbol", "GeneSymbol", "hgnc_symbol")
  hit <- preferred[preferred %in% names(df)]
  if (length(hit) > 0) return(hit[[1]])
  return(NULL)
}

detect_ensembl_column <- function(df) {
  name_hits <- grep("ensembl|ensg", names(df), ignore.case = TRUE, value = TRUE)
  if (length(name_hits) > 0) return(name_hits[[1]])

  char_cols <- names(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
  for (col in char_cols) {
    vals <- as.character(df[[col]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) == 0) next
    clean <- sub("\\..*$", "", vals)
    frac_ensg <- mean(grepl("^ENSG[0-9]+$", clean, ignore.case = TRUE))
    if (is.finite(frac_ensg) && frac_ensg >= 0.5) return(col)
  }
  return(NULL)
}

symbol_col <- detect_symbol_column(de)
gene_symbols <- NULL

if (!is.null(symbol_col)) {
  gene_symbols <- as.character(de[[symbol_col]])
} else {
  ensembl_col <- detect_ensembl_column(de)
  if (is.null(ensembl_col)) {
    stop(
      "Could not detect a gene identifier column. Expected SYMBOL or ENSEMBL-like IDs.",
      " Available columns: ", paste(names(de), collapse = ", ")
    )
  }

  ensembl_ids <- as.character(de[[ensembl_col]])
  ensembl_ids <- sub("\\..*$", "", ensembl_ids)

  gene_symbols <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  gene_symbols <- unname(as.character(gene_symbols))
}

rank_col <- NULL
if ("t" %in% names(de)) {
  rank_col <- "t"
} else if ("logFC" %in% names(de)) {
  rank_col <- "logFC"
}

if (is.null(rank_col)) {
  stop(
    "Could not find a ranking column. Expected 't' or 'logFC'.",
    " Available columns: ", paste(names(de), collapse = ", ")
  )
}

rank_values <- suppressWarnings(as.numeric(de[[rank_col]]))

rank_df <- data.frame(
  symbol = gene_symbols,
  score = rank_values,
  stringsAsFactors = FALSE
)

rank_df <- rank_df[!is.na(rank_df$symbol) & nzchar(rank_df$symbol), , drop = FALSE]
rank_df <- rank_df[is.finite(rank_df$score), , drop = FALSE]

if (nrow(rank_df) == 0) {
  stop("No valid genes remained after identifier processing and rank filtering.")
}

# Keep one value per symbol for fgsea (largest absolute statistic).
rank_df <- rank_df[order(abs(rank_df$score), decreasing = TRUE), , drop = FALSE]
rank_df <- rank_df[!duplicated(rank_df$symbol), , drop = FALSE]

stats <- rank_df$score
names(stats) <- rank_df$symbol
stats <- sort(stats, decreasing = TRUE)

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
if (nrow(hallmark) == 0) {
  stop("msigdbr returned no Hallmark gene sets for Homo sapiens.")
}

pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

fgsea_res <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = stats
)

fgsea_res <- fgsea_res[order(fgsea_res$padj, fgsea_res$pval), , drop = FALSE]
# CSV can't store list-columns (e.g. leadingEdge). Convert any list columns to strings.
list_cols <- vapply(fgsea_res, is.list, logical(1))
if (any(list_cols)) {
  for (nm in names(fgsea_res)[list_cols]) {
    fgsea_res[[nm]] <- vapply(
      fgsea_res[[nm]],
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  }
}
utils::write.csv(fgsea_res, out_table, row.names = FALSE)

if (nrow(fgsea_res) == 0) {
  stop("fgsea completed but returned zero rows.")
}

plot_df <- fgsea_res[order(fgsea_res$NES, decreasing = TRUE), , drop = FALSE]
plot_df <- utils::head(plot_df, 15)

grDevices::png(filename = out_plot, width = 2000, height = 1400, res = 200)
op <- graphics::par(no.readonly = TRUE)
on.exit({
  graphics::par(op)
  grDevices::dev.off()
}, add = TRUE)

cols <- ifelse(plot_df$NES >= 0, "#D55E00", "#0072B2")
graphics::par(mar = c(5, 16, 4, 2))
graphics::barplot(
  rev(plot_df$NES),
  names.arg = rev(plot_df$pathway),
  horiz = TRUE,
  las = 1,
  col = rev(cols),
  xlab = "Normalized Enrichment Score (NES)",
  main = "Top 15 Hallmark Pathways (fgsea)"
)
graphics::abline(v = 0, lty = 2, col = "gray30")

message("Wrote: ", out_table)
message("Wrote: ", out_plot)
