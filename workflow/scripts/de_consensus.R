#!/usr/bin/env Rscript
# Consensus differential expression: DESeq2 (primary) + edgeR-QLF + limma-voom.
# Operates on RAW integer counts, WITHIN a single assay, with optional covariates.
# Primary contrast = (first non-reference level of `factor`) vs reference level.
# Refs: Love et al. 2014 Genome Biol 15:550 (DESeq2); Robinson et al. 2010 Bioinformatics
#       26:139 + Chen et al. 2016 F1000Res (edgeR-QLF); Law et al. 2014 Genome Biol 15:R29
#       (limma-voom); Conesa et al. 2016 Genome Biol 17:13 (best practices).
# NOTE: validate on a real cohort before production. Handles the common 2-level case;
# extend for >2 levels / multiple contrasts as needed.
suppressMessages({library(DESeq2); library(edgeR); library(limma)})

a <- commandArgs(trailingOnly = TRUE)
matrix_f <- a[1]; samples_f <- a[2]; assay <- a[3]; outdir <- a[4]
factor_col <- a[5]; reference <- a[6]
covariates <- if (length(a) >= 7 && nzchar(a[7])) strsplit(a[7], ",")[[1]] else character(0)
min_count <- as.numeric(a[8]); fdr <- as.numeric(a[9]); lfc <- as.numeric(a[10])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- load harmonized RAW-count matrix (the '# assay' line is a comment) -------
counts <- read.table(matrix_f, header = TRUE, sep = "\t", row.names = 1,
                     comment.char = "#", check.names = FALSE)
meta <- read.table(samples_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta <- meta[meta$assay == assay, , drop = FALSE]
keep_samples <- intersect(meta$sample, colnames(counts))
stopifnot(length(keep_samples) >= 2)
counts <- round(as.matrix(counts[, keep_samples, drop = FALSE]))
meta <- meta[match(keep_samples, meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample

# factor + reference level; drop covariates that are constant within this assay subset
grp <- factor(meta[[factor_col]]); grp <- relevel(grp, ref = reference)
meta[[factor_col]] <- grp
covariates <- covariates[sapply(covariates, function(c) c %in% colnames(meta) &&
                                  length(unique(meta[[c]])) > 1)]
levs <- setdiff(levels(grp), reference)
stopifnot(length(levs) >= 1)
target <- levs[1]                                   # primary contrast: target vs reference
design_rhs <- paste(c(covariates, factor_col), collapse = " + ")
form <- as.formula(paste("~", design_rhs))
message("DE assay=", assay, " contrast=", target, "_vs_", reference, " design=~", design_rhs)

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = form)
dds <- dds[rowSums(counts(dds)) >= min_count, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c(factor_col, target, reference))
res <- as.data.frame(res)
deseq2 <- data.frame(gene_id = rownames(res), log2FC = res$log2FoldChange, padj = res$padj)
write.table(deseq2, file.path(outdir, "deseq2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# variance-stabilized matrix for QC/visualization
try(write.table(as.data.frame(assay(vst(dds, blind = FALSE))),
                file.path(outdir, "vst_matrix.tsv"), sep = "\t", quote = FALSE), silent = TRUE)

# ---- edgeR-QLF (TMM) --------------------------------------------------------
y <- DGEList(counts = counts, group = grp)
y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                              # TMM (Robinson & Oshlack 2010)
mm <- model.matrix(form, data = meta)
y <- estimateDisp(y, mm)
fit <- glmQLFit(y, mm)
coef_name <- paste0(factor_col, target)
qlf <- glmQLFTest(fit, coef = coef_name)
et <- topTags(qlf, n = Inf)$table
edger <- data.frame(gene_id = rownames(et), log2FC = et$logFC, padj = et$FDR)
write.table(edger, file.path(outdir, "edger.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- limma-voom -------------------------------------------------------------
v <- voom(y, mm)
lf <- eBayes(lmFit(v, mm))
tt <- topTable(lf, coef = coef_name, number = Inf, sort.by = "none")
limma_res <- data.frame(gene_id = rownames(tt), log2FC = tt$logFC, padj = tt$adj.P.Val)
write.table(limma_res, file.path(outdir, "limma.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- consensus (>= 2 of 3 methods significant) ------------------------------
sig <- function(d) d$gene_id[which(!is.na(d$padj) & d$padj < fdr & abs(d$log2FC) >= lfc)]
calls <- list(DESeq2 = sig(deseq2), edgeR = sig(edger), limma = sig(limma_res))
all_genes <- unique(c(deseq2$gene_id, edger$gene_id, limma_res$gene_id))
nmeth <- sapply(all_genes, function(g) sum(sapply(calls, function(s) g %in% s)))
m <- function(d) d[match(all_genes, d$gene_id), c("log2FC", "padj")]
cons <- data.frame(gene_id = all_genes,
                   deseq2_log2FC = m(deseq2)$log2FC, deseq2_padj = m(deseq2)$padj,
                   edger_log2FC = m(edger)$log2FC,   edger_padj = m(edger)$padj,
                   limma_log2FC = m(limma_res)$log2FC, limma_padj = m(limma_res)$padj,
                   n_methods_sig = nmeth, contrast = paste0(target, "_vs_", reference))
cons <- cons[order(-cons$n_methods_sig, cons$deseq2_padj), ]
cons <- cons[cons$n_methods_sig >= 2, ]
write.table(cons, file.path(outdir, "consensus_DE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("consensus DE genes (>=2 methods): ", nrow(cons))
