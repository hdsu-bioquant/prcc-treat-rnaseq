# quantify.smk — canonical per-library expression products and run-level matrices
# Primary expression basis: STAR unstranded raw GeneCounts for every library.
# UMI-bearing libraries additionally receive molecule-level HTSeq counts.

# ---- UMI molecule counts ---------------------------------------------------#
rule htseq_dedup_umi:
    input:
        bam = join(RESTRICTED, "libraries/{library}/alignments/umi_dedup.bam"),
        gtf = GTF
    output:
        counts = temp(join(INTERMEDIATE, "libraries/{library}/quantification/umi_molecule_counts.tsv"))
    params:
        strand = htseq_strand
    wildcard_constraints:
        library = UMI_LIBRARY_PATTERN
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
        ),
        script = join(SCRIPT_DIR, "build_gene_expression.py")
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})"
        python {params.script:q} \
          {input.counts:q} {input.lengths:q} {params.assay:q} {params.umi_counts:q} {output.tsv:q}
        """

# ---- run-level matrices ----------------------------------------------------#
rule merge_raw_count_matrix:
    input:
        expand(join(RESULTS, "libraries/{lib}/gene_expression.tsv"), lib=LIBRARIES)
    output:
        matrix = join(RESULTS, "matrices/raw_gene_counts.tsv")
    params:
        samplesheet = config["samples"],
        results = RESULTS,
        script = join(SCRIPT_DIR, "merge_expression.py")
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        "python {params.script:q} {params.samplesheet:q} {params.results:q} unstranded {output.matrix:q}"

rule merge_umi_molecule_matrix:
    input:
        expand(join(RESULTS, "libraries/{lib}/gene_expression.tsv"), lib=UMI_LIBRARIES)
    output:
        matrix = join(RESULTS, "matrices/umi_molecule_counts.tsv")
    params:
        samplesheet = config["samples"],
        results = RESULTS,
        script = join(SCRIPT_DIR, "merge_expression.py")
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        "python {params.script:q} {params.samplesheet:q} {params.results:q} umi_molecule_count {output.matrix:q}"
