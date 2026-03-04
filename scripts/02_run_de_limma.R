#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing package: data.table")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Missing package: limma")
})

infile <- "data/raw/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"
meta_file <- "data/metadata/samples.csv"

out_table <- "results/tables/de_limma_results.csv"
out_pca <- "results/figures/pca.png"
out_volcano <- "results/figures/volcano.png"

if (!file.exists(infile)) stop("Missing input: ", infile, "\nRun: make download", call. = FALSE)
if (!file.exists(meta_file)) stop("Missing metadata: ", meta_file, "\nRun: make metadata", call. = FALSE)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Read gz fast (works on mac/linux)
cmd <- sprintf("gzip -dc %s", shQuote(infile))
dt <- data.table::fread(cmd = cmd)

gene_id <- dt[[1]]
expr <- as.matrix(dt[, -1, with = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- gene_id

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta$treatment <- factor(meta$treatment, levels = c("UT", "ADR"))
meta$cell_line <- factor(meta$cell_line)

# Align metadata to expression columns
if (!all(colnames(expr) %in% meta$sample_id)) {
  missing <- setdiff(colnames(expr), meta$sample_id)
  stop("Metadata missing samples:\n- ", paste(missing, collapse = "\n- "), call. = FALSE)
}
meta <- meta[match(colnames(expr), meta$sample_id), ]
stopifnot(identical(meta$sample_id, colnames(expr)))

# Decide method: integer-like? -> voom+limma; else log-limma
nr <- nrow(expr)
sub_idx <- seq_len(min(1000, nr))
xsub <- expr[sub_idx, , drop = FALSE]
integer_like <- all(abs(xsub - round(xsub)) < 1e-6, na.rm = TRUE)

design <- model.matrix(~ cell_line + treatment, data = meta)
coef_name <- "treatmentADR"
if (!coef_name %in% colnames(design)) stop("Design missing ", coef_name, call. = FALSE)

if (integer_like && requireNamespace("edgeR", quietly = TRUE)) {
  message("Integer-like values detected -> edgeR + voom + limma")
  y <- edgeR::DGEList(counts = expr)
  y <- edgeR::calcNormFactors(y)
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  emat <- v$E
} else {
  message("Using limma on log2(x+1) (norm_counts / non-integers)")
  emat <- log2(expr + 1)
  fit <- limma::lmFit(emat, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
}

tt$gene_id <- rownames(tt)
tt <- tt[, c("gene_id", setdiff(colnames(tt), "gene_id"))]
write.csv(tt, out_table, row.names = FALSE)
message("Wrote: ", out_table, " (n=", nrow(tt), ")")

# PCA (top variable genes)
vars <- apply(emat, 1, var, na.rm = TRUE)
top <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
pca <- prcomp(t(emat[top, , drop = FALSE]), scale. = TRUE)

col_vec <- as.integer(meta$treatment)
png(out_pca, width = 900, height = 700)
plot(pca$x[,1], pca$x[,2],
     pch = 19, col = col_vec,
     xlab = sprintf("PC1 (%.1f%%)", 100 * (pca$sdev[1]^2 / sum(pca$sdev^2))),
     ylab = sprintf("PC2 (%.1f%%)", 100 * (pca$sdev[2]^2 / sum(pca$sdev^2))),
     main = "PCA (top variable genes)")
legend("topright", legend = levels(meta$treatment), pch = 19, col = seq_along(levels(meta$treatment)))
dev.off()
message("Saved: ", out_pca)

# Volcano
logFC <- tt$logFC
adjP <- tt$adj.P.Val
png(out_volcano, width = 900, height = 700)
plot(logFC, -log10(adjP),
     pch = 19,
     xlab = "log2 Fold Change (ADR vs UT)",
     ylab = "-log10(adj.P.Val)",
     main = "Volcano plot")
abline(v = c(-1, 1), lty = 2)
abline(h = -log10(0.05), lty = 2)
dev.off()
message("Saved: ", out_volcano)
