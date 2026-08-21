# qc.smk — canonical QC metrics + MultiQC; detailed raw QC remains restricted


def fastqc_inputs(wc):
    fqs = [samples.loc[wc.library, "fq1"]]
    if samples.loc[wc.library, "fq2"] != "-":
        fqs.append(samples.loc[wc.library, "fq2"])
    return fqs


rule fastqc:
    input:
        fastqc_inputs
    output:
        outdir = directory(join(RESTRICTED, "libraries/{library}/qc/fastqc"))
    container: IMG["fastqc"]
    threads: 2
    resources:
        mem_mb = 4000, runtime = 60
    shell:
        r"""
        rm -rf {output.outdir}
        mkdir -p {output.outdir}
        fastqc -o {output.outdir} -t {threads} {input}
        """


rule qc_metrics:
    input:
        star = join(RESTRICTED, "libraries/{library}/logs/{library}.Log.final.out"),
        expr = join(RESULTS, "libraries/{library}/gene_expression.tsv")
    output:
        tsv = join(RESULTS, "libraries/{library}/qc_metrics.tsv")
    params:
        sample_id = lambda wc: biological_sample_id(wc.library),
        assay = lambda wc: samples.loc[wc.library, "assay"],
        layout = lambda wc: samples.loc[wc.library, "layout"],
        has_umi = lambda wc: samples.loc[wc.library, "has_umi"]
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})"
        python workflow/scripts/build_qc_metrics.py \
          {input.star} {input.expr} {wildcards.library} {params.sample_id} \
          {params.assay} {params.layout} {params.has_umi} {output.tsv}
        """


rule merge_qc_metrics:
    input:
        expand(join(RESULTS, "libraries/{lib}/qc_metrics.tsv"), lib=LIBRARIES)
    output:
        tsv = join(RESULTS, "qc/qc_metrics.tsv")
    params:
        samplesheet = config["samples"], results = RESULTS
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        "python workflow/scripts/merge_qc_metrics.py {params.samplesheet} {params.results} {output.tsv}"


rule qc_multiqc_custom:
    input:
        join(RESULTS, "qc/qc_metrics.tsv")
    output:
        tsv = temp(join(INTERMEDIATE, "qc/prcc_qc_metrics_mqc.tsv"))
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})"
        python workflow/scripts/qc_to_multiqc.py {input} {output.tsv}
        """


rule multiqc:
    input:
        fastqc = expand(join(RESTRICTED, "libraries/{lib}/qc/fastqc"), lib=LIBRARIES),
        star = expand(join(RESTRICTED, "libraries/{lib}/logs/{lib}.Log.final.out"), lib=LIBRARIES),
        custom = join(INTERMEDIATE, "qc/prcc_qc_metrics_mqc.tsv"),
        versions = join(INTERMEDIATE, "qc/prcc_pipeline_mqc_versions.yml"),
        fastp = (
            expand(join(RESTRICTED, "libraries/{lib}/qc/{lib}.fastp.json"), lib=FL_LIBRARIES)
            if config.get("full_length", {}).get("trim_adapters", False) else []
        ),
        rseqc = (
            expand(join(RESULTS, "libraries/{lib}/{lib}.rseqc_read_distribution.txt"), lib=FL_LIBRARIES)
            if config.get("modules", {}).get("rseqc", False) else []
        )
    output:
        html = join(RESULTS, "qc/multiqc_report.html"),
        data = directory(join(RESTRICTED, "qc/multiqc_data"))
    params:
        build = join(INTERMEDIATE, "qc/multiqc_build"),
        config = os.path.abspath(os.path.join(workflow.basedir, "config", "multiqc_portable.yaml"))
    container: IMG["multiqc"]
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        r"""
        rm -rf {params.build}
        mkdir -p {params.build} "$(dirname {output.html})" "$(dirname {output.data})"
        multiqc {input} -c {params.config} -o {params.build} -n multiqc_report.html -f
        cp {params.build}/multiqc_report.html {output.html}

        # MultiQC names its data directory from the report stem (for example
        # multiqc_report_data in v1.21). Keep the pRCC-RNA-Seq output contract
        # stable as restricted/qc/multiqc_data regardless of MultiQC version.
        data_src="{params.build}/multiqc_report_data"
        if [[ ! -d "$data_src" ]]; then
            data_src="$(find {params.build} -mindepth 1 -maxdepth 1 -type d -name '*_data' -print -quit)"
        fi
        if [[ -z "$data_src" || ! -d "$data_src" ]]; then
            echo "ERROR: MultiQC data directory not found under {params.build}" >&2
            exit 1
        fi

        rm -rf {output.data}
        mv "$data_src" {output.data}
        """


# ---- optional: RSeQC read_distribution (full-length PE only) ---------------#
GENE_BED12 = join(REFDIR, "annotation.bed12")

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
        bam = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        bed = GENE_BED12
    output:
        txt = join(RESULTS, "libraries/{library}/{library}.rseqc_read_distribution.txt")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    container: IMG["rseqc"]
    resources:
        mem_mb = 8000, runtime = 120
    shell:
        r"""
        mkdir -p "$(dirname {output.txt})"
        read_distribution.py -i {input.bam} -r {input.bed} > {output.txt}
        """
