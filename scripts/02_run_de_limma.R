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

# Volcano (layered plotting + labels + write label table)
logFC <- tt$logFC
adjP  <- tt$adj.P.Val
gene_ids <- tt$gene_id

sig_cut <- 0.05
lfc_cut <- 1

is_sig  <- !is.na(adjP) & adjP < sig_cut
is_up   <- is_sig & !is.na(logFC) & logFC >  lfc_cut
is_down <- is_sig & !is.na(logFC) & logFC < -lfc_cut

# Ensembl -> HUGO (handles ENSG..._at and version suffix)
sym <- gene_ids
clean <- gene_ids
clean <- sub("_at$", "", clean)
clean <- sub("\\..*$", "", clean)

if (any(grepl("^ENSG", clean))) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install mapping pkgs: BiocManager::install(c('AnnotationDbi','org.Hs.eg.db'))")
  }
  is_ens <- grepl("^ENSG", clean)
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = clean[is_ens],
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  sym[is_ens] <- ifelse(is.na(mapped) | mapped == "", gene_ids[is_ens], mapped)
}

# choose labels = top 10 up + top 10 down by adj.P.Val (among is_up / is_down)
pick_top <- function(idx, n=10) {
  if (!length(idx)) return(integer(0))
  idx <- idx[order(adjP[idx], na.last = NA)]
  idx[seq_len(min(n, length(idx)))]
}
lab_up   <- pick_top(which(is_up), 10)
lab_down <- pick_top(which(is_down), 10)
lab_idx  <- unique(c(lab_up, lab_down))

# write out which genes are being labeled (debug + transparency)
if (length(lab_idx) > 0) {
  write.csv(
    data.frame(
      gene_id = gene_ids[lab_idx],
      symbol  = sym[lab_idx],
      logFC   = logFC[lab_idx],
      adjP    = adjP[lab_idx]
    ),
    file = "results/tables/volcano_labels.csv",
    row.names = FALSE
  )
}

png(out_volcano, width = 900, height = 700)

# 1) background: all points in very light grey with transparency
plot(
  logFC, -log10(adjP),
  pch = 16,
  col = grDevices::rgb(0.6, 0.6, 0.6, 0.15),
  cex = 0.6,
  xlab = "log2FC (ADR vs UT)",
  ylab = "-log10(adj.P.Val)",
  main = "Volcano"
)

# 2) overlay significant up/down on top (bigger + stronger alpha)
if (any(is_up)) {
  points(logFC[is_up], -log10(adjP[is_up]),
         pch = 16, col = grDevices::rgb(1, 0, 0, 0.75), cex = 0.8)
}
if (any(is_down)) {
  points(logFC[is_down], -log10(adjP[is_down]),
         pch = 16, col = grDevices::rgb(0, 0, 1, 0.75), cex = 0.8)
}

abline(v = c(-lfc_cut, lfc_cut), lty = 2)
abline(h = -log10(sig_cut), lty = 2)

legend(
  "topright",
  legend = c("Up (adj.P<0.05 & logFC>1)", "Down (adj.P<0.05 & logFC<-1)", "Other"),
  col = c(grDevices::rgb(1,0,0,0.9), grDevices::rgb(0,0,1,0.9), grDevices::rgb(0.6,0.6,0.6,0.4)),
  pch = 16,
  bty = "n"
)

# labels (draw last so they appear on top)
if (length(lab_idx) > 0) {
  text(
    x = logFC[lab_idx],
    y = -log10(adjP[lab_idx]),
    labels = sym[lab_idx],
    pos = 3,
    cex = 0.85
  )
}

dev.off()
message("Saved: ", out_volcano)
