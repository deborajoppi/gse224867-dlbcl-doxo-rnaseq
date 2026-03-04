cat > scripts/02_run_de_limma.R <<'EOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing package: data.table")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Missing package: limma")
})

infile <- "data/raw/GSE224867_norm_counts_DLBCL_TIS.tsv.gz"
meta_file <- "data/metadata/samples.csv"
out_table <- "results/tables/DE_ADR_vs_UT.csv"
out_pca <- "results/figures/pca.png"
out_volcano <- "results/figures/volcano.png"

if (!file.exists(infile)) stop("Input file not found: ", infile, call. = FALSE)
if (!file.exists(meta_file)) stop("Metadata not found: ", meta_file, "\nRun: make metadata", call. = FALSE)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Fast read gz via system gzip (portable on mac/linux)
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
  stop("Metadata missing some samples:\n- ", paste(missing, collapse = "\n- "), call. = FALSE)
}
meta <- meta[match(colnames(expr), meta$sample_id), ]
stopifnot(identical(meta$sample_id, colnames(expr)))

# Check integer-like on a subset (helps decide voom vs log-limma)
nr <- nrow(expr)
sub_idx <- seq_len(min(1000, nr))
xsub <- expr[sub_idx, , drop = FALSE]
integer_like <- all(abs(xsub - round(xsub)) < 1e-6, na.rm = TRUE)

design <- model.matrix(~ cell_line + treatment, data = meta)
coef_name <- "treatmentADR"
if (!coef_name %in% colnames(design)) {
  stop("Could not find coefficient ", coef_name, " in design matrix.", call. = FALSE)
}

if (integer_like && requireNamespace("edgeR", quietly = TRUE)) {
  message("Data look integer-like -> using edgeR + voom + limma.")
  y <- edgeR::DGEList(counts = expr)
  y <- edgeR::calcNormFactors(y)
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  emat <- v$E
} else {
  message("Data not integer-like (or edgeR missing) -> using limma on log2(x+1).")
  emat <- log2(expr + 1)
  fit <- limma::lmFit(emat, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
}

tt$gene_id <- rownames(tt)
tt <- tt[, c("gene_id", setdiff(colnames(tt), "gene_id"))]
write.csv(tt, out_table, row.names = FALSE)
message("Wrote: ", out_table, " (n=", nrow(tt), ")")

# PCA on top variable genes
vars <- apply(emat, 1, var, na.rm = TRUE)
top <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
pca <- prcomp(t(emat[top, , drop = FALSE]), scale. = TRUE)

png(out_pca, width = 900, height = 700)
plot(pca$x[,1], pca$x[,2],
     pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", 100 * (pca$sdev[1]^2 / sum(pca$sdev^2))),
     ylab = sprintf("PC2 (%.1f%%)", 100 * (pca$sdev[2]^2 / sum(pca$sdev^2))),
     main = "PCA (top variable genes)")
legend("topright", legend = levels(meta$treatment), pch = 19, col = seq_along(levels(meta$treatment)))
# color points by treatment (after legend, so default plot shows points first)
points(pca$x[,1], pca$x[,2], pch = 19, col = as.integer(meta$treatment))
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
EOF

chmod +x scripts/02_run_de_limma.R
