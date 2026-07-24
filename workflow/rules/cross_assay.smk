# cross_assay.smk — Branch-A SECONDARY: 3'-restricted gene-count matrix (HTSeq)
# Emulates QuantSeq 3'-tag counting: count full-length reads only in each gene's 3'
# exonic window. To keep HTSeq (scheme tool) fast on full-depth PE BAMs, we FIRST
# subset the BAM to the 3' windows (samtools view -L), so HTSeq only processes the
# few % of reads that fall there (~minutes instead of ~hours). Off by default.

rule make_3prime_windows:
    input:
        gtf = GTF
    output:
        gtf = THREEP_GTF,
        bed = THREEP_BED
    params:
        w = config["cross_assay"]["window_bp"]
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 30, slurm_partition = "single"
    shell:
        "python workflow/scripts/make_3prime_windows.py {input.gtf} {output.gtf} {output.bed} {params.w}"

rule subset_3prime_bam:
    # Pre-filter the full-length BAM to reads overlapping the 3' windows, then NAME-sort.
    # Name-sorting is essential: -L subsetting orphans many mates (mate outside the
    # windows), and HTSeq -r pos would buffer them -> OOM. With a name-sorted BAM +
    # HTSeq -r name, mates are adjacent / orphans flush immediately -> flat memory, fast.
    input:
        bam = join(RESULTS, "full_length/{sample}/{sample}.Aligned.sortedByCoord.bam"),
        bed = THREEP_BED
    output:
        bam = temp(join(RESULTS, "full_length/{sample}/{sample}.3prime_subset_ns.bam"))
    threads: 4
    container: IMG["samtools"]
    resources:
        mem_mb = 16000, runtime = 120, slurm_partition = "single"
    shell:
        "samtools view -u -L {input.bed} {input.bam} | "
        "samtools sort -n -@ {threads} -m 2G -o {output.bam} -"

rule htseq_3prime_fl:
    # HTSeq over the 3' GTF, on the NAME-sorted pre-filtered BAM (-r name -> flat memory).
    # unstranded (-s no) to remove length AND strand-protocol differences vs QuantSeq.
    input:
        bam = join(RESULTS, "full_length/{sample}/{sample}.3prime_subset_ns.bam"),
        gtf = THREEP_GTF
    output:
        counts = join(RESULTS, "full_length/{sample}/{sample}.3prime_htseq_counts.tsv")
    threads: 2
    container: IMG["htseq"]
    resources:
        mem_mb = 8000, runtime = 120, slurm_partition = "single"
    shell:
        "htseq-count -f bam -r name -s no -t exon -i gene_id -m union "
        "--nonunique none {input.bam} {input.gtf} > {output.counts}"

rule merge_3prime_matrix:
    input:
        fl = expand(join(RESULTS, "full_length/{s}/{s}.3prime_htseq_counts.tsv"), s=FL_SAMPLES)
    output:
        matrix = join(RESULTS, "matrix/three_prime_restricted_matrix.tsv")
    params:
        samplesheet = config["samples"], results = RESULTS
    container: IMG["py"]
    resources:
        mem_mb = 8000, runtime = 60, slurm_partition = "single"
    shell:
        "python workflow/scripts/merge_htseq.py {params.samplesheet} {params.results} 3prime {output.matrix}"
