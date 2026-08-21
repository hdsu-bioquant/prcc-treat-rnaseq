# qc.smk — FastQC per sample + a single MultiQC report (GDC reports FastQC too)

def fastqc_inputs(wc):
    fqs = [samples.loc[wc.library, "fq1"]]
    if samples.loc[wc.library, "fq2"] != "-":
        fqs.append(samples.loc[wc.library, "fq2"])
    return fqs

rule fastqc:
    input:
        fastqc_inputs
    output:
        done = touch(join(RESULTS, "qc/fastqc/{library}.done"))
    params:
        outdir = join(RESULTS, "qc/fastqc")
    container: IMG["fastqc"]
    threads: 2
    resources:
        mem_mb = 4000, runtime = 60
    shell:
        "mkdir -p {params.outdir} && fastqc -o {params.outdir} -t {threads} {input}"

rule multiqc:
    input:
        fastqc = expand(join(RESULTS, "qc/fastqc/{s}.done"), s=LIBRARIES),
        fl     = expand(join(RESULTS, "full_length/{s}/{s}.star_gene_counts.tsv"), s=FL_LIBRARIES),
        qs     = expand(join(RESULTS, "quantseq/{s}/{s}.star_gene_counts.tsv"), s=QS_LIBRARIES),
    output:
        html = join(RESULTS, "qc/multiqc_report.html")
    params:
        scan   = RESULTS,
        outdir = join(RESULTS, "qc")
    container: IMG["multiqc"]
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        "multiqc {params.scan} -o {params.outdir} -n multiqc_report.html -f"

# ---- optional: RSeQC read_distribution (full-length PE only) ----------------#
GENE_BED12 = join(REFDIR, "gencode.v36.bed12")

rule gtf_to_bed12:
    input:
        gtf = GTF
    output:
        bed = GENE_BED12
    container: IMG["py"]
    shell:
        "python workflow/scripts/gtf_to_bed12.py {input.gtf} {output.bed}"

rule rseqc_read_distribution:
    input:
        bam = join(RESULTS, "full_length/{library}/{library}.Aligned.sortedByCoord.bam"),
        bed = GENE_BED12
    output:
        txt = join(RESULTS, "full_length/{library}/qc/{library}.rseqc.read_distribution.txt")
    container: IMG["rseqc"]
    resources:
        mem_mb = 8000, runtime = 120
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} > {output.txt}"
