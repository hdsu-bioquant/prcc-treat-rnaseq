# quantify.smk — canonical per-library expression products and run-level matrices
# Primary expression basis: STAR unstranded raw GeneCounts for every library.
# UMI-bearing libraries additionally receive molecule-level HTSeq counts.

# ---- UMI molecule counts ---------------------------------------------------#
rule htseq_dedup_fl:
    input:
        bam = join(RESTRICTED, "libraries/{library}/alignments/umi_dedup.bam"),
        gtf = GTF
    output:
        counts = temp(join(INTERMEDIATE, "libraries/{library}/quantification/umi_molecule_counts.tsv"))
    params:
        strand = htseq_strand
    wildcard_constraints:
        library = FL_UMI_LIBRARY_PATTERN
    threads: 2
    container: IMG["htseq"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.counts})"
        htseq-count -f bam -r pos -s {params.strand} -t exon -i gene_id -m union \
          --nonunique none {input.bam} {input.gtf} > {output.counts}
        """

rule htseq_dedup_qs:
    input:
        bam = join(RESTRICTED, "libraries/{library}/alignments/umi_dedup.bam"),
        gtf = GTF
    output:
        counts = temp(join(INTERMEDIATE, "libraries/{library}/quantification/umi_molecule_counts.tsv"))
    params:
        strand = htseq_strand
    wildcard_constraints:
        library = QS_UMI_LIBRARY_PATTERN
    threads: 2
    container: IMG["htseq"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.counts})"
        htseq-count -f bam -r pos -s {params.strand} -t exon -i gene_id -m union \
          --nonunique none {input.bam} {input.gtf} > {output.counts}
        """

# ---- canonical per-library expression table -------------------------------#
rule gene_expression:
    input:
        unpack(gene_expression_inputs)
    output:
        tsv = join(RESULTS, "libraries/{library}/gene_expression.tsv")
    params:
        assay = lambda wc: library_assay(wc.library),
        umi_counts = lambda wc: (
            join(INTERMEDIATE, "libraries", wc.library, "quantification", "umi_molecule_counts.tsv")
            if library_has_umi(wc.library) else "-"
        )
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})"
        python workflow/scripts/build_gene_expression.py \
          {input.counts} {input.lengths} {params.assay} {params.umi_counts} {output.tsv}
        """

# ---- run-level matrices ----------------------------------------------------#
rule merge_raw_count_matrix:
    input:
        expand(join(RESULTS, "libraries/{lib}/gene_expression.tsv"), lib=LIBRARIES)
    output:
        matrix = join(RESULTS, "matrices/raw_gene_counts.tsv")
    params:
        samplesheet = config["samples"], results = RESULTS
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        "python workflow/scripts/merge_expression.py {params.samplesheet} {params.results} unstranded {output.matrix}"

rule merge_umi_molecule_matrix:
    input:
        expand(join(RESULTS, "libraries/{lib}/gene_expression.tsv"), lib=UMI_LIBRARIES)
    output:
        matrix = join(RESULTS, "matrices/umi_molecule_counts.tsv")
    params:
        samplesheet = config["samples"], results = RESULTS
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        "python workflow/scripts/merge_expression.py {params.samplesheet} {params.results} umi_molecule_count {output.matrix}"
