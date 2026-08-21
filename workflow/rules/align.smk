# align.smk — STAR (GDC recipe) for both assay branches
# UMI extraction is library-level and may precede either branch.

# ---------------------------- full-length branch ----------------------------#
rule umi_extract_fl:
    input:
        fq1 = raw_fq1,
        fq2 = raw_fq2
    output:
        fq1 = join(RESULTS, "full_length/{library}/umi/{library}_R1.umi.fastq.gz"),
        fq2 = join(RESULTS, "full_length/{library}/umi/{library}_R2.umi.fastq.gz")
    params:
        pattern = lambda wc: library_umi_pattern(wc.library)
    log:
        join(RESULTS, "full_length/{library}/umi/{library}.umi_extract.log")
    container: IMG["umitools"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        "umi_tools extract --stdin {input.fq1} --bc-pattern={params.pattern} "
        "--stdout {output.fq1} --read2-in {input.fq2} --read2-out {output.fq2} "
        "--ignore-read-pair-suffixes --log {log}"

rule fastp_fl:
    # OFF by default (GDC does not trim). If enabled, UMI extraction happens first.
    input:
        unpack(fl_pretrim_input)
    output:
        fq1  = join(RESULTS, "full_length/{library}/trim/{library}_R1.trim.fastq.gz"),
        fq2  = join(RESULTS, "full_length/{library}/trim/{library}_R2.trim.fastq.gz"),
        html = join(RESULTS, "full_length/{library}/trim/{library}.fastp.html"),
        json = join(RESULTS, "full_length/{library}/trim/{library}.fastp.json")
    threads: 4
    container: IMG["fastp"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        "fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} "
        "--detect_adapter_for_pe -w {threads} -h {output.html} -j {output.json}"

rule star_fl:
    input:
        unpack(fl_star_input),
        idx = STAR_IDX_DONE
    output:
        bam    = join(RESULTS, "full_length/{library}/{library}.Aligned.out.bam"),
        counts = join(RESULTS, "full_length/{library}/{library}.ReadsPerGene.out.tab"),
        tx     = join(RESULTS, "full_length/{library}/{library}.Aligned.toTranscriptome.out.bam"),
        chim   = join(RESULTS, "full_length/{library}/{library}.Chimeric.out.junction"),
        log    = join(RESULTS, "full_length/{library}/{library}.Log.final.out")
    params:
        idxdir = STAR_IDX,
        prefix = join(RESULTS, "full_length/{library}/{library}."),
        extra  = config["star"]["gdc_params"],
        rg     = lambda wc: "ID:%s SM:%s PL:ILLUMINA" % (wc.library, biological_sample_id(wc.library)),
        tmp    = lambda wc: join(TMPDIR, "star_fl_" + wc.library)
    threads: config["star"]["threads"]
    container: IMG["star"]
    resources:
        mem_mb = 64000, runtime = 1440
    shell:
        r"""
        mkdir -p "$(dirname {params.tmp})"
        rm -rf {params.tmp}
        STAR --runThreadN {threads} --genomeDir {params.idxdir} \
             --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat \
             --outSAMattrRGline {params.rg} \
             --outFileNamePrefix {params.prefix} --outTmpDir {params.tmp} \
             {params.extra}
        """

rule sort_index_fl:
    input:
        bam = join(RESULTS, "full_length/{library}/{library}.Aligned.out.bam")
    output:
        bam = join(RESULTS, "full_length/{library}/{library}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "full_length/{library}/{library}.Aligned.sortedByCoord.bam.bai")
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && samtools index {output.bam}"

rule umi_dedup_fl:
    input:
        bam = join(RESULTS, "full_length/{library}/{library}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "full_length/{library}/{library}.Aligned.sortedByCoord.bam.bai")
    output:
        bam = join(RESULTS, "full_length/{library}/{library}.dedup.bam")
    log:
        join(RESULTS, "full_length/{library}/{library}.umi_dedup.log")
    container: IMG["umitools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        "umi_tools dedup --paired -I {input.bam} -S {output.bam} --log {log}"

# ------------------------------ QuantSeq branch -----------------------------#
rule umi_extract_qs:
    input:
        fq1 = raw_fq1
    output:
        fq = join(RESULTS, "quantseq/{library}/{library}.umi.fastq.gz")
    params:
        pattern = lambda wc: library_umi_pattern(wc.library)
    log:
        join(RESULTS, "quantseq/{library}/{library}.umi_extract.log")
    container: IMG["umitools"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        "umi_tools extract --stdin {input.fq1} --bc-pattern={params.pattern} "
        "--stdout {output.fq} --log {log}"

rule bbduk_qs:
    input:
        fq = qs_trim_input
    output:
        fq = join(RESULTS, "quantseq/{library}/{library}.trim.fastq.gz")
    threads: 4
    container: IMG["bbmap"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        # Lexogen QuantSeq FWD: right-trim 3' adapter read-through + poly(A).
        r"""
        RES=$(ls -d /usr/local/opt/bbmap*/resources 2>/dev/null | head -1)
        bbduk.sh in={input.fq} out={output.fq} \
          ref="$RES/truseq_rna.fa.gz" \
          literal=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
          k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 threads={threads}
        """

rule star_qs:
    input:
        fq  = join(RESULTS, "quantseq/{library}/{library}.trim.fastq.gz"),
        idx = STAR_IDX_DONE
    output:
        bam    = join(RESULTS, "quantseq/{library}/{library}.Aligned.out.bam"),
        counts = join(RESULTS, "quantseq/{library}/{library}.ReadsPerGene.out.tab"),
        log    = join(RESULTS, "quantseq/{library}/{library}.Log.final.out")
    params:
        idxdir = STAR_IDX,
        prefix = join(RESULTS, "quantseq/{library}/{library}."),
        extra  = config["star"]["gdc_params"],
        rg     = lambda wc: "ID:%s SM:%s PL:ILLUMINA" % (wc.library, biological_sample_id(wc.library)),
        tmp    = lambda wc: join(TMPDIR, "star_qs_" + wc.library)
    threads: config["star"]["threads"]
    container: IMG["star"]
    resources:
        mem_mb = 64000, runtime = 1440
    shell:
        r"""
        mkdir -p "$(dirname {params.tmp})"
        rm -rf {params.tmp}
        STAR --runThreadN {threads} --genomeDir {params.idxdir} \
             --readFilesIn {input.fq} --readFilesCommand zcat \
             --outSAMattrRGline {params.rg} \
             --outFileNamePrefix {params.prefix} --outTmpDir {params.tmp} \
             {params.extra}
        """

rule sort_index_qs:
    input:
        bam = join(RESULTS, "quantseq/{library}/{library}.Aligned.out.bam")
    output:
        bam = join(RESULTS, "quantseq/{library}/{library}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "quantseq/{library}/{library}.Aligned.sortedByCoord.bam.bai")
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && samtools index {output.bam}"

rule umi_dedup_qs:
    input:
        bam = join(RESULTS, "quantseq/{library}/{library}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "quantseq/{library}/{library}.Aligned.sortedByCoord.bam.bai")
    output:
        bam = join(RESULTS, "quantseq/{library}/{library}.dedup.bam")
    log:
        join(RESULTS, "quantseq/{library}/{library}.umi_dedup.log")
    container: IMG["umitools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        "umi_tools dedup -I {input.bam} -S {output.bam} --log {log}"
