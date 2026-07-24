# references.smk — obtain the EXACT GDC reference resources (GDC adherence)
# GRCh38.d1.vd1 genome, GENCODE v36 GTF, and the GDC-built STAR 2.7.5c index.
# Source: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files

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
    # Use the GDC-distributed STAR 2.7.5c index (sjdbOverhang 100) for byte-exact
    # GDC adherence. (To rebuild instead, swap this for a STAR --runMode genomeGenerate.)
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

rule gene_lengths:
    # Non-overlapping (union-exon) gene length per gene_id from GENCODE v36,
    # used as the "Gene Length (L)" term in the GDC FPKM / FPKM-UQ / TPM formulas.
    input:
        gtf = GTF
    output:
        tsv = GENE_LENGTHS
    container: IMG["py"]
    shell:
        "python workflow/scripts/gene_lengths.py {input.gtf} {output.tsv}"
