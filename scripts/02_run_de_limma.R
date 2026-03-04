#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

infile <- "data/raw/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"
meta_file <- "data/metadata/samples.csv"

out_table <- "results/tables/de_limma_results.csv"
out_pca <- "results/figures/pca.png"
out_volcano <- "results/figures/volcano.png"

stopifnot(file.exists(infile))
stopifnot(file.exists(meta_file))

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

dt <- fread(cmd = sprintf("gzip -dc %s", shQuote(infile)))
gene_id <- dt[[1]]
expr <- as.matrix(dt[, -1, with = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- gene_id

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id","cell_line","treatment") %in% names(meta)))
meta$cell_line <- factor(meta$cell_line)
meta$treatment <- factor(meta$treatment, levels = c("UT","ADR"))

# align
meta <- meta[match(colnames(expr), meta$sample_id), ]
if (any(is.na(meta$sample_id))) stop("Some expression columns not found in metadata.", call. = FALSE)
stopifnot(identical(meta$sample_id, colnames(expr)))

# check integer-like on subset
xsub <- expr[seq_len(min(1000, nrow(expr))), , drop = FALSE]
integer_like <- all(abs(xsub - round(xsub)) < 1e-6, na.rm = TRUE)

design <- model.matrix(~ cell_line + treatment, data = meta)
coef_name <- "treatmentADR"

if (integer_like && requireNamespace("edgeR", quietly = TRUE)) {
  message("Integer-like -> edgeR + voom + limma")
  y <- edgeR::DGEList(counts = expr)
  y <- edgeR::calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  emat <- v$E
} else {
  message("Non-integer / norm_counts -> limma on log2(x+1)")
  emat <- log2(expr + 1)
  fit <- eBayes(lmFit(emat, design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
}

tt$gene_id <- rownames(tt)
tt <- tt[, c("gene_id", setdiff(colnames(tt), "gene_id"))]
write.csv(tt, out_table, row.names = FALSE)
message("Wrote: ", out_table)

# PCA
vars <- apply(emat, 1, var, na.rm = TRUE)
top <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
pca <- prcomp(t(emat[top, , drop = FALSE]), scale. = TRUE)
cols <- as.integer(meta$treatment)

png(out_pca, width = 900, height = 700)
plot(pca$x[,1], pca$x[,2], pch = 19, col = cols,
     xlab = "PC1", ylab = "PC2", main = "PCA (top variable genes)")
legend("topright", legend = levels(meta$treatment), pch = 19, col = seq_along(levels(meta$treatment)))
dev.off()

# Volcano (color significant up/down)
logFC <- tt$logFC
adjP  <- tt$adj.P.Val

sig_cut <- 0.05
lfc_cut <- 1

is_sig <- !is.na(adjP) & (adjP < sig_cut)
is_up  <- is_sig & !is.na(logFC) & (logFC >  lfc_cut)
is_down<- is_sig & !is.na(logFC) & (logFC < -lfc_cut)

# default color for non-significant
pt_col <- rep("grey70", length(logFC))
pt_col[is_up]   <- "red"
pt_col[is_down] <- "blue"

png(out_volcano, width = 900, height = 700)
plot(
  logFC, -log10(adjP),
  pch = 19,
  col = pt_col,
  xlab = "log2FC (ADR vs UT)",
  ylab = "-log10(adj.P.Val)",
  main = "Volcano"
)
abline(v = c(-lfc_cut, lfc_cut), lty = 2)
abline(h = -log10(sig_cut), lty = 2)
legend(
  "topright",
  legend = c("Up (sig)", "Down (sig)", "Not sig"),
  col = c("red", "blue", "grey70"),
  pch = 19,
  bty = "n"
)
dev.off()
message("Saved: ", out_volcano)
