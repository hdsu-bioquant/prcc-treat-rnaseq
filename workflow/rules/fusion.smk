# fusion.smk — OPTIONAL gene-fusion detection (full-length branch only)
# GDC fusion tools: STAR-Fusion v1.6 + Arriba v1.1.0. QuantSeq 3' SE data cannot
# yield split/spanning chimeric reads, so fusion is full-length-only by design.

rule arriba:
    input:
        bam   = join(RESULTS, "full_length/{library}/{library}.Aligned.sortedByCoord.bam"),
        gtf   = GTF,
        fasta = FASTA
    output:
        fusions = join(RESULTS, "full_length/{library}/fusion/{library}.arriba.fusions.tsv")
    params:
        discarded = join(RESULTS, "full_length/{library}/fusion/{library}.arriba.discarded.tsv")
    threads: 4
    container: IMG["arriba"]
    resources:
        mem_mb = 32000, runtime = 240
    shell:
        "arriba -x {input.bam} -g {input.gtf} -a {input.fasta} "
        "-o {output.fusions} -O {params.discarded}"

rule star_fusion:
    input:
        chim = join(RESULTS, "full_length/{library}/{library}.Chimeric.out.junction")
    output:
        pred = join(RESULTS, "full_length/{library}/fusion/{library}.starfusion.predictions.tsv")
    params:
        ctat   = config.get("ctat_genome_lib", ""),   # external CTAT lib (not a tracked input)
        outdir = join(RESULTS, "full_length/{library}/fusion/starfusion")
    threads: 4
    container: IMG["starfusion"]
    resources:
        mem_mb = 32000, runtime = 480
    shell:
        r"""
        STAR-Fusion --chimeric_junction {input.chim} --genome_lib_dir {params.ctat} \
                    --CPU {threads} --output_dir {params.outdir}
        cp {params.outdir}/star-fusion.fusion_predictions.tsv {output.pred}
        """
