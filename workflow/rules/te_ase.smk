# te_ase.smk — OPTIONAL transposable-element and allele-specific-expression
# modules (full-length branch only; optional toolkit).

rule tecount:
    # TEtranscripts/TEcount: TE-family quantification. Requires a TE annotation GTF
    # (config: te_gtf), supplied externally (TEtranscripts curated TE GTF).
    input:
        bam = join(RESULTS, "full_length/{sample}/{sample}.Aligned.sortedByCoord.bam"),
        gtf = GTF
    output:
        tab = join(RESULTS, "full_length/{sample}/te/{sample}.TEcount.cntTable")
    params:
        te_gtf  = config["te_gtf"],
        project = join(RESULTS, "full_length/{sample}/te/{sample}.TEcount")
    container: IMG["tetx"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        TEcount --sortByPos --BAM {input.bam} --GTF {input.gtf} --TE {params.te_gtf} \
                --project {params.project}
        test -s {output.tab}
        """

rule ase_readcounter:
    # GATK ASEReadCounter: allele-specific expression at heterozygous germline sites.
    # Requires an EXTERNAL germline VCF (config: ase_germline_vcf) — germline variant
    # calling is DNA-seq/WES and is OUT OF SCOPE for this RNA-seq pipeline.
    input:
        bam   = join(RESULTS, "full_length/{sample}/{sample}.Aligned.sortedByCoord.bam"),
        bai   = join(RESULTS, "full_length/{sample}/{sample}.Aligned.sortedByCoord.bam.bai"),
        fasta = FASTA
    output:
        tsv = join(RESULTS, "full_length/{sample}/ase/{sample}.ASEReadCounter.tsv")
    params:
        vcf = config["ase_germline_vcf"]
    container: IMG["gatk"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        "gatk ASEReadCounter -R {input.fasta} -I {input.bam} -V {params.vcf} -O {output.tsv}"
