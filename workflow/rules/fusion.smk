# fusion.smk — OPTIONAL gene-fusion detection (full-length branch only)
# Outputs remain restricted by default until the consortium finalizes the sharing policy.

rule arriba:
    input:
        bam   = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        gtf   = GTF,
        fasta = FASTA
    output:
        fusions = join(RESTRICTED, "libraries/{library}/features/fusion/arriba.fusions.tsv")
    params:
        discarded = join(RESTRICTED, "libraries/{library}/features/fusion/arriba.discarded.tsv")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    threads: 4
    container: IMG["arriba"]
    resources:
        mem_mb = 32000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.fusions})"
        arriba -x {input.bam} -g {input.gtf} -a {input.fasta} \
          -o {output.fusions} -O {params.discarded}
        """

rule star_fusion:
    input:
        chim = join(INTERMEDIATE, "libraries/{library}/star/Chimeric.out.junction")
    output:
        pred = join(RESTRICTED, "libraries/{library}/features/fusion/starfusion.predictions.tsv")
    params:
        ctat   = config.get("ctat_genome_lib", ""),
        outdir = lambda wc: join(INTERMEDIATE, "libraries", wc.library, "starfusion")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    threads: 4
    container: IMG["starfusion"]
    resources:
        mem_mb = 32000, runtime = 480
    shell:
        r"""
        mkdir -p {params.outdir} "$(dirname {output.pred})"
        STAR-Fusion --chimeric_junction {input.chim} --genome_lib_dir {params.ctat} \
                    --CPU {threads} --output_dir {params.outdir}
        cp {params.outdir}/star-fusion.fusion_predictions.tsv {output.pred}
        """
