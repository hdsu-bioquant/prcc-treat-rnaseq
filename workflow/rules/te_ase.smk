# te_ase.smk — OPTIONAL transposable-element and allele-specific-expression
# modules (full-length only). These outputs remain restricted by default pending
# a later scope/privacy review.

rule tecount:
    input:
        bam = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        gtf = GTF
    output:
        tab = join(RESTRICTED, "libraries/{library}/features/te/TEcount.cntTable")
    params:
        te_gtf  = config.get("te_gtf", ""),
        project = lambda wc: join(RESTRICTED, "libraries", wc.library, "features", "te", "TEcount")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    container: IMG["tetx"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.tab})"
        TEcount --sortByPos --BAM {input.bam} --GTF {input.gtf} --TE {params.te_gtf} \
                --project {params.project}
        test -s {output.tab}
        """

rule ase_readcounter:
    input:
        bam   = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        bai   = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam.bai"),
        fasta = FASTA
    output:
        tsv = join(RESTRICTED, "libraries/{library}/features/ase/ASEReadCounter.tsv")
    params:
        vcf = config.get("ase_germline_vcf", "")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    container: IMG["gatk"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})"
        gatk ASEReadCounter -R {input.fasta} -I {input.bam} -V {params.vcf} -O {output.tsv}
        """
