# de.smk — OPTIONAL differential expression (downstream, central, WITHIN assay)
# Consensus of DESeq2 (primary) + edgeR-QLF + limma-voom on RAW counts.
# Not part of the GDC standard (GDC stops at counts); run on the harmonized matrix.
# Refs: Love et al. 2014; Robinson et al. 2010 / Chen et al. 2016; Law et al. 2014; Conesa et al. 2016.

rule differential_expression:
    input:
        matrix  = join(RESULTS, "matrix/gene_counts_matrix.tsv"),
        samples = config["samples"]
    output:
        consensus = join(RESULTS, "de/{assay}/consensus_DE.tsv")
    params:
        outdir     = join(RESULTS, "de/{assay}"),
        factor     = config["de"]["factor"],
        reference  = config["de"]["reference_level"],
        covariates = config["de"]["covariates"],
        min_count  = config["de"]["min_count"],
        fdr        = config["de"]["fdr"],
        lfc        = config["de"]["lfc"]
    wildcard_constraints:
        assay = "|".join([re.escape(a) for a in DE_ASSAYS]) if DE_ASSAYS else "noassay"
    container: IMG["de"]
    resources:
        mem_mb = 16000, runtime = 180, slurm_partition = "single"
    shell:
        "Rscript workflow/scripts/de_consensus.R {input.matrix} {input.samples} {wildcards.assay} "
        "{params.outdir} {params.factor} {params.reference} '{params.covariates}' "
        "{params.min_count} {params.fdr} {params.lfc}"
