# references.smk — GDC production references or a local test/reference bundle
# Production default: exact GDC GRCh38.d1.vd1 + GENCODE v36 resources.
# Local mode: use supplied FASTA/GTF and build a STAR index (used by synthetic tests).

REFERENCE_MODE = config["reference"].get("mode", "gdc")

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

elif REFERENCE_MODE == "gdc":
    rule gdc_genome_fasta:
        output:
            fasta = FASTA
        params:
            url = config["reference"]["genome_fasta_tar_url"],
            d   = REFDIR
        shell:
            r"""
            mkdir -p {params.d}
            wget -O {params.d}/GRCh38.d1.vd1.fa.tar.gz "{params.url}"
            tar -xzf {params.d}/GRCh38.d1.vd1.fa.tar.gz -C {params.d}
            test -s {output.fasta}
            """

    rule gdc_gtf:
        output:
            gtf = GTF
        params:
            url = config["reference"]["gtf_gz_url"],
            d   = REFDIR
        shell:
            r"""
            mkdir -p {params.d}
            wget -O {params.d}/gencode.v36.annotation.gtf.gz "{params.url}"
            gunzip -f {params.d}/gencode.v36.annotation.gtf.gz
            test -s {output.gtf}
            """

    rule gdc_star_index:
        # Use the GDC-distributed STAR 2.7.5c index (sjdbOverhang 100) for GDC adherence.
        output:
            done = STAR_IDX_DONE
        params:
            url = config["reference"]["star_index_tgz_url"],
            idx = STAR_IDX,
            d   = REFDIR
        shell:
            r"""
            mkdir -p {params.idx}
            wget -O {params.d}/star_index.tgz "{params.url}"
            tar -xzf {params.d}/star_index.tgz -C {params.idx} --strip-components=1
            test -s {output.done}
            """
else:
    raise ValueError("Unsupported reference.mode: %s (expected 'gdc' or 'local')" % REFERENCE_MODE)

rule gene_lengths:
    # Non-overlapping (union-exon) gene length per gene_id from the active GTF,
    # used as the length term in the optional GDC-style FPKM / FPKM-UQ / TPM formulas.
    input:
        gtf = GTF
    output:
        tsv = GENE_LENGTHS
    container: IMG["py"]
    shell:
        "python workflow/scripts/gene_lengths.py {input.gtf} {output.tsv}"
