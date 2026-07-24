# quantify.smk — pRCC-TREAT primary + secondary quantification
# PRIMARY (both branches): STAR --quantMode GeneCounts -> raw counts (unstranded primary,
#          stranded columns kept for diagnostics). No dedup, no normalization here.
# SECONDARY: QuantSeq UMI-dedup -> HTSeq-count; (optional, off) full-length FPKM/FPKM-UQ/TPM.

# ---- PRIMARY: STAR gene counts, both branches -------------------------------#
rule star_gene_counts_fl:
    input:
        counts = join(RESULTS, "full_length/{sample}/{sample}.ReadsPerGene.out.tab")
    output:
        tsv = join(RESULTS, "full_length/{sample}/{sample}.star_gene_counts.tsv")
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30, slurm_partition = "single"
    shell:
        "python workflow/scripts/extract_star_counts.py {input.counts} {output.tsv}"

rule star_gene_counts_qs:
    # QuantSeq PRIMARY = STAR GeneCounts on UMI-extracted (if any) but NON-deduplicated reads.
    input:
        counts = join(RESULTS, "quantseq/{sample}/{sample}.ReadsPerGene.out.tab")
    output:
        tsv = join(RESULTS, "quantseq/{sample}/{sample}.star_gene_counts.tsv")
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30, slurm_partition = "single"
    shell:
        "python workflow/scripts/extract_star_counts.py {input.counts} {output.tsv}"

# ---- SECONDARY: QuantSeq UMI-deduplicated HTSeq counts ----------------------#
rule htseq_dedup_qs:
    input:
        bam = join(RESULTS, "quantseq/{sample}/{sample}.dedup.bam"),
        gtf = GTF
    output:
        counts = join(RESULTS, "quantseq/{sample}/{sample}.dedup_htseq_counts.tsv")
    params:
        strand = QS_HTSEQ_STRAND      # QuantSeq FWD -> "yes" (forward)
    threads: 2
    container: IMG["htseq"]
    resources:
        mem_mb = 8000, runtime = 240, slurm_partition = "single"
    shell:
        "htseq-count -f bam -r pos -s {params.strand} -t exon -i gene_id -m union "
        "--nonunique none {input.bam} {input.gtf} > {output.counts}"

# ---- OPTIONAL (off by default): GDC-style FPKM/FPKM-UQ/TPM ------------------#
rule fpkm_tpm_fl:
    input:
        counts  = join(RESULTS, "full_length/{sample}/{sample}.ReadsPerGene.out.tab"),
        lengths = GENE_LENGTHS
    output:
        tsv = join(RESULTS, "full_length/{sample}/{sample}.augmented_star_gene_counts.tsv")
    params:
        col = config["full_length"]["count_column"]
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30, slurm_partition = "single"
    shell:
        "python workflow/scripts/augment_star_counts.py "
        "{input.counts} {input.lengths} {params.col} {output.tsv}"

# ---- cohort matrices -------------------------------------------------------#
rule merge_count_matrix:
    # PRIMARY cohort matrix: STAR unstranded raw counts, both branches.
    input:
        fl = expand(join(RESULTS, "full_length/{s}/{s}.star_gene_counts.tsv"), s=FL_SAMPLES),
        qs = expand(join(RESULTS, "quantseq/{s}/{s}.star_gene_counts.tsv"), s=QS_SAMPLES)
    output:
        matrix = join(RESULTS, "matrix/gene_counts_matrix.tsv")
    params:
        samplesheet = config["samples"], results = RESULTS
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60, slurm_partition = "single"
    shell:
        "python workflow/scripts/merge_counts.py {params.samplesheet} {params.results} {output.matrix}"

rule merge_umidedup_matrix:
    # SECONDARY cohort matrix: QuantSeq UMI-deduplicated HTSeq counts.
    input:
        qs = expand(join(RESULTS, "quantseq/{s}/{s}.dedup_htseq_counts.tsv"), s=QS_SAMPLES)
    output:
        matrix = join(RESULTS, "matrix/quantseq_umidedup_matrix.tsv")
    params:
        samplesheet = config["samples"], results = RESULTS
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60, slurm_partition = "single"
    shell:
        "python workflow/scripts/merge_htseq.py {params.samplesheet} {params.results} dedup {output.matrix}"
