# align.smk — STAR (GDC recipe) for both assay branches
# Full-length poly-A (PE): optional fastp -> STAR 2-pass (GDC params) -> sort
# QuantSeq 3' (SE+UMI): umi_tools extract -> BBDuk polyA/adapter -> STAR SE -> umi_tools dedup

# ---------------------------- full-length branch ----------------------------#
rule fastp_fl:
    # OFF by default (GDC does not trim). Enabled via full_length.trim_adapters: true
    input:
        fq1 = raw_fq1,
        fq2 = raw_fq2
    output:
        fq1  = join(RESULTS, "full_length/{sample}/trim/{sample}_R1.trim.fastq.gz"),
        fq2  = join(RESULTS, "full_length/{sample}/trim/{sample}_R2.trim.fastq.gz"),
        html = join(RESULTS, "full_length/{sample}/trim/{sample}.fastp.html"),
        json = join(RESULTS, "full_length/{sample}/trim/{sample}.fastp.json")
    threads: 4
    container: IMG["fastp"]
    resources:
        mem_mb = 8000, runtime = 240, slurm_partition = "single"
    shell:
        "fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} "
        "--detect_adapter_for_pe -w {threads} -h {output.html} -j {output.json}"

rule star_fl:
    input:
        unpack(fl_star_input),
        idx = STAR_IDX_DONE
    output:
        bam    = join(RESULTS, "full_length/{sample}/{sample}.Aligned.out.bam"),
        counts = join(RESULTS, "full_length/{sample}/{sample}.ReadsPerGene.out.tab"),
        tx     = join(RESULTS, "full_length/{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
        chim   = join(RESULTS, "full_length/{sample}/{sample}.Chimeric.out.junction"),
        log    = join(RESULTS, "full_length/{sample}/{sample}.Log.final.out")
    params:
        idxdir = STAR_IDX,
        prefix = join(RESULTS, "full_length/{sample}/{sample}."),
        extra  = config["star"]["gdc_params"],
        rg     = lambda wc: "ID:%s SM:%s PL:ILLUMINA" % (wc.sample, wc.sample),
        tmp    = lambda wc: join(TMPDIR, "star_fl_" + wc.sample)
    threads: config["star"]["threads"]
    container: IMG["star"]
    resources:
        mem_mb = 64000, runtime = 1440, slurm_partition = "single"
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
        bam = join(RESULTS, "full_length/{sample}/{sample}.Aligned.out.bam")
    output:
        bam = join(RESULTS, "full_length/{sample}/{sample}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "full_length/{sample}/{sample}.Aligned.sortedByCoord.bam.bai")
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 240, slurm_partition = "single"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && samtools index {output.bam}"

# ------------------------------ QuantSeq branch -----------------------------#
rule umi_extract_qs:
    input:
        fq1 = raw_fq1
    output:
        fq = join(RESULTS, "quantseq/{sample}/{sample}.umi.fastq.gz")
    params:
        pattern = config["quantseq"]["umi_pattern"]
    log:
        join(RESULTS, "quantseq/{sample}/{sample}.umi_extract.log")
    container: IMG["umitools"]
    resources:
        mem_mb = 8000, runtime = 240, slurm_partition = "single"
    shell:
        "umi_tools extract --stdin {input.fq1} --bc-pattern={params.pattern} "
        "--stdout {output.fq} --log {log}"

rule bbduk_qs:
    input:
        fq = qs_trim_input          # umi-extracted R1 (has_umi) or raw R1 (no UMI)
    output:
        fq = join(RESULTS, "quantseq/{sample}/{sample}.trim.fastq.gz")
    threads: 4
    container: IMG["bbmap"]
    resources:
        mem_mb = 8000, runtime = 240, slurm_partition = "single"
    shell:
        # Lexogen QuantSeq FWD: right-trim the 3' adapter read-through + polyA.
        # BBDuk's bundled adapter file lives in the bbmap install's resources/ dir inside
        # the container (NOT the CWD), so resolve it; polyA is trimmed via a literal poly-A kmer.
        r"""
        RES=$(ls -d /usr/local/opt/bbmap*/resources 2>/dev/null | head -1)
        bbduk.sh in={input.fq} out={output.fq} \
          ref="$RES/truseq_rna.fa.gz" \
          literal=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
          k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 threads={threads}
        """

rule star_qs:
    input:
        fq  = join(RESULTS, "quantseq/{sample}/{sample}.trim.fastq.gz"),
        idx = STAR_IDX_DONE
    output:
        bam    = join(RESULTS, "quantseq/{sample}/{sample}.Aligned.out.bam"),
        counts = join(RESULTS, "quantseq/{sample}/{sample}.ReadsPerGene.out.tab"),
        log    = join(RESULTS, "quantseq/{sample}/{sample}.Log.final.out")
    params:
        idxdir = STAR_IDX,
        prefix = join(RESULTS, "quantseq/{sample}/{sample}."),
        extra  = config["star"]["gdc_params"],
        rg     = lambda wc: "ID:%s SM:%s PL:ILLUMINA" % (wc.sample, wc.sample),
        tmp    = lambda wc: join(TMPDIR, "star_qs_" + wc.sample)
    threads: config["star"]["threads"]
    container: IMG["star"]
    resources:
        mem_mb = 64000, runtime = 1440, slurm_partition = "single"
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
        bam = join(RESULTS, "quantseq/{sample}/{sample}.Aligned.out.bam")
    output:
        bam = join(RESULTS, "quantseq/{sample}/{sample}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "quantseq/{sample}/{sample}.Aligned.sortedByCoord.bam.bai")
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 240, slurm_partition = "single"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && samtools index {output.bam}"

rule umi_dedup_qs:
    # UMI-aware deduplication (coordinate-only dedup is WRONG for UMI data).
    input:
        bam = join(RESULTS, "quantseq/{sample}/{sample}.Aligned.sortedByCoord.bam"),
        bai = join(RESULTS, "quantseq/{sample}/{sample}.Aligned.sortedByCoord.bam.bai")
    output:
        bam = join(RESULTS, "quantseq/{sample}/{sample}.dedup.bam")
    log:
        join(RESULTS, "quantseq/{sample}/{sample}.umi_dedup.log")
    container: IMG["umitools"]
    resources:
        mem_mb = 16000, runtime = 240, slurm_partition = "single"
    shell:
        "umi_tools dedup -I {input.bam} -S {output.bam} --log {log}"
