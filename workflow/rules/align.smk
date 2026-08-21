# align.smk — assay-aware preprocessing + STAR using the common library output layout
# UMI extraction is library-level and may precede either assay branch.

# ---------------------------- full-length branch ----------------------------#
rule umi_extract_fl:
    input:
        fq1 = raw_fq1,
        fq2 = raw_fq2
    output:
        fq1 = temp(join(INTERMEDIATE, "libraries/{library}/preprocess/R1.umi.fastq.gz")),
        fq2 = temp(join(INTERMEDIATE, "libraries/{library}/preprocess/R2.umi.fastq.gz"))
    params:
        pattern = lambda wc: library_umi_pattern(wc.library)
    log:
        join(RESTRICTED, "libraries/{library}/logs/umi_extract.log")
    wildcard_constraints:
        library = FL_UMI_LIBRARY_PATTERN
    container: IMG["umitools"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.fq1})" "$(dirname {log})"
        umi_tools extract --stdin {input.fq1} --bc-pattern={params.pattern} \
          --stdout {output.fq1} --read2-in {input.fq2} --read2-out {output.fq2} \
          --ignore-read-pair-suffixes --log {log}
        """

rule fastp_fl:
    # OFF by default (GDC does not trim). If enabled, UMI extraction happens first.
    input:
        unpack(fl_pretrim_input)
    output:
        fq1  = temp(join(INTERMEDIATE, "libraries/{library}/preprocess/R1.trim.fastq.gz")),
        fq2  = temp(join(INTERMEDIATE, "libraries/{library}/preprocess/R2.trim.fastq.gz")),
        html = join(RESTRICTED, "libraries/{library}/qc/{library}.fastp.html"),
        json = join(RESTRICTED, "libraries/{library}/qc/{library}.fastp.json")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    threads: 4
    container: IMG["fastp"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.fq1})" "$(dirname {output.html})"
        fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} \
          --detect_adapter_for_pe -w {threads} -h {output.html} -j {output.json}
        """

rule star_fl:
    input:
        unpack(fl_star_input),
        idx = STAR_IDX_DONE
    output:
        bam    = temp(join(INTERMEDIATE, "libraries/{library}/star/Aligned.out.bam")),
        counts = temp(join(INTERMEDIATE, "libraries/{library}/star/ReadsPerGene.out.tab")),
        tx     = temp(join(INTERMEDIATE, "libraries/{library}/star/Aligned.toTranscriptome.out.bam")),
        chim   = temp(join(INTERMEDIATE, "libraries/{library}/star/Chimeric.out.junction")),
        log    = join(RESTRICTED, "libraries/{library}/logs/{library}.Log.final.out")
    params:
        idxdir = STAR_IDX,
        stardir = lambda wc: join(INTERMEDIATE, "libraries", wc.library, "star"),
        extra  = config["star"]["gdc_params"],
        rg     = lambda wc: "ID:%s SM:%s PL:ILLUMINA" % (wc.library, biological_sample_id(wc.library)),
        tmp    = lambda wc: join(TMPDIR, "star_fl_" + wc.library)
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    threads: config["star"]["threads"]
    container: IMG["star"]
    resources:
        mem_mb = 64000, runtime = 1440
    shell:
        r"""
        mkdir -p {params.stardir} "$(dirname {params.tmp})" "$(dirname {output.log})"
        rm -rf {params.tmp}
        STAR --runThreadN {threads} --genomeDir {params.idxdir} \
             --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat \
             --outSAMattrRGline {params.rg} \
             --outFileNamePrefix {params.stardir}/ --outTmpDir {params.tmp} \
             {params.extra}
        cp {params.stardir}/Log.final.out {output.log}
        """

rule sort_index_fl:
    input:
        bam = join(INTERMEDIATE, "libraries/{library}/star/Aligned.out.bam")
    output:
        bam = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        bai = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam.bai")
    wildcard_constraints:
        library = FL_LIBRARY_PATTERN
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})"
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

rule umi_dedup_fl:
    input:
        bam = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        bai = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam.bai")
    output:
        bam = join(RESTRICTED, "libraries/{library}/alignments/umi_dedup.bam")
    log:
        join(RESTRICTED, "libraries/{library}/logs/umi_dedup.log")
    wildcard_constraints:
        library = FL_UMI_LIBRARY_PATTERN
    container: IMG["umitools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
        umi_tools dedup --paired -I {input.bam} -S {output.bam} --log {log}
        """

# ------------------------------ QuantSeq branch -----------------------------#
rule umi_extract_qs:
    input:
        fq1 = raw_fq1
    output:
        fq = temp(join(INTERMEDIATE, "libraries/{library}/preprocess/R1.umi.fastq.gz"))
    params:
        pattern = lambda wc: library_umi_pattern(wc.library)
    log:
        join(RESTRICTED, "libraries/{library}/logs/umi_extract.log")
    wildcard_constraints:
        library = QS_UMI_LIBRARY_PATTERN
    container: IMG["umitools"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.fq})" "$(dirname {log})"
        umi_tools extract --stdin {input.fq1} --bc-pattern={params.pattern} \
          --stdout {output.fq} --log {log}
        """

rule bbduk_qs:
    input:
        fq = qs_trim_input
    output:
        fq = temp(join(INTERMEDIATE, "libraries/{library}/preprocess/R1.trim.fastq.gz"))
    wildcard_constraints:
        library = QS_LIBRARY_PATTERN
    threads: 4
    container: IMG["bbmap"]
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        # Lexogen QuantSeq FWD: right-trim 3' adapter read-through + poly(A).
        r"""
        mkdir -p "$(dirname {output.fq})"
        RES=$(ls -d /usr/local/opt/bbmap*/resources 2>/dev/null | head -1)
        bbduk.sh in={input.fq} out={output.fq} \
          ref="$RES/truseq_rna.fa.gz" \
          literal=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
          k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 threads={threads}
        """

rule star_qs:
    input:
        fq  = join(INTERMEDIATE, "libraries/{library}/preprocess/R1.trim.fastq.gz"),
        idx = STAR_IDX_DONE
    output:
        bam    = temp(join(INTERMEDIATE, "libraries/{library}/star/Aligned.out.bam")),
        counts = temp(join(INTERMEDIATE, "libraries/{library}/star/ReadsPerGene.out.tab")),
        tx     = temp(join(INTERMEDIATE, "libraries/{library}/star/Aligned.toTranscriptome.out.bam")),
        chim   = temp(join(INTERMEDIATE, "libraries/{library}/star/Chimeric.out.junction")),
        log    = join(RESTRICTED, "libraries/{library}/logs/{library}.Log.final.out")
    params:
        idxdir = STAR_IDX,
        stardir = lambda wc: join(INTERMEDIATE, "libraries", wc.library, "star"),
        extra  = config["star"]["gdc_params"],
        rg     = lambda wc: "ID:%s SM:%s PL:ILLUMINA" % (wc.library, biological_sample_id(wc.library)),
        tmp    = lambda wc: join(TMPDIR, "star_qs_" + wc.library)
    wildcard_constraints:
        library = QS_LIBRARY_PATTERN
    threads: config["star"]["threads"]
    container: IMG["star"]
    resources:
        mem_mb = 64000, runtime = 1440
    shell:
        r"""
        mkdir -p {params.stardir} "$(dirname {params.tmp})" "$(dirname {output.log})"
        rm -rf {params.tmp}
        STAR --runThreadN {threads} --genomeDir {params.idxdir} \
             --readFilesIn {input.fq} --readFilesCommand zcat \
             --outSAMattrRGline {params.rg} \
             --outFileNamePrefix {params.stardir}/ --outTmpDir {params.tmp} \
             {params.extra}
        cp {params.stardir}/Log.final.out {output.log}
        """

rule sort_index_qs:
    input:
        bam = join(INTERMEDIATE, "libraries/{library}/star/Aligned.out.bam")
    output:
        bam = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        bai = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam.bai")
    wildcard_constraints:
        library = QS_LIBRARY_PATTERN
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})"
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

rule umi_dedup_qs:
    input:
        bam = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam"),
        bai = join(RESTRICTED, "libraries/{library}/alignments/genomic.sorted.bam.bai")
    output:
        bam = join(RESTRICTED, "libraries/{library}/alignments/umi_dedup.bam")
    log:
        join(RESTRICTED, "libraries/{library}/logs/umi_dedup.log")
    wildcard_constraints:
        library = QS_UMI_LIBRARY_PATTERN
    container: IMG["umitools"]
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
        umi_tools dedup -I {input.bam} -S {output.bam} --log {log}
        """
