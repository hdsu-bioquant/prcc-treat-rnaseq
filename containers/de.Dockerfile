# Differential-expression image for RNA-Seq-pRCC (DESeq2 + edgeR + limma-voom).
# Build, then convert to a SIF for Apptainer use (Docker only as a build source):
#   docker build -f containers/de.Dockerfile -t rnaseq-prcc-de:1.0 .
#   apptainer build de.sif docker-daemon://rnaseq-prcc-de:1.0
# Then point IMG["de"] in workflow/rules/common.smk at the SIF (or keep the
# bioconductor base + this layer in your registry) and pin the digest for release.
FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN R -e "BiocManager::install(c('DESeq2','edgeR','limma','apeglm'), update=FALSE, ask=FALSE)" \
 && R -e "stopifnot(all(c('DESeq2','edgeR','limma','apeglm') %in% rownames(installed.packages())))"
