# references.smk — local development index generation + derived gene metadata
#
# Production GDC references are installed before analysis and validated during
# workflow initialization in common.smk. Normal Snakemake execution never
# downloads the GDC FASTA/GTF/STAR index.

if REFERENCE_MODE == "local":
    rule local_star_index:
        input:
            fasta = FASTA,
            gtf   = GTF
        output:
            done = STAR_IDX_DONE
        params:
            idx = STAR_IDX,
            sjdb = int(config["reference"].get("sjdb_overhang", 100)),
            sa_arg = ("--genomeSAindexNbases %s" % config["reference"]["genome_sa_index_nbases"]
                      if "genome_sa_index_nbases" in config["reference"] else ""),
            chrbin_arg = ("--genomeChrBinNbits %s" % config["reference"]["genome_chr_bin_nbits"]
                          if "genome_chr_bin_nbits" in config["reference"] else "")
        threads: config["star"]["threads"]
        container: IMG["star"]
        resources:
            mem_mb = 8000, runtime = 60
        shell:
            r"""
            mkdir -p {params.idx}
            STAR --runThreadN {threads} --runMode genomeGenerate \
                 --genomeDir {params.idx} --genomeFastaFiles {input.fasta} \
                 --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdb} \
                 {params.sa_arg} {params.chrbin_arg}
            test -s {output.done}
            """


rule gene_lengths:
    # Non-overlapping (union-exon) gene length per gene_id from the active GTF,
    # used as the length/gene-metadata input for full-length GDC-style FPKM / FPKM-UQ / TPM.
    input:
        gtf = GTF
    output:
        tsv = GENE_LENGTHS
    container: IMG["py"]
    shell:
        "python workflow/scripts/gene_lengths.py {input.gtf} {output.tsv}"
