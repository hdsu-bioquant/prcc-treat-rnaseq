# pRCC-RNA-Seq

RNA-seq processing pipeline developed for harmonized papillary renal-cell carcinoma analyses, with a supported core covering full-length paired-end and QuantSeq single-end libraries, optional UMI handling, GDC-aligned expression processing, standardized QC, and portable run outputs.

The pipeline is designed for both independent use and harmonized multi-site operation. The documentation is therefore split by use case.

## Documentation

### General users

Use [`docs/users/general/`](docs/users/general/README.md) for a conventional pipeline guide covering installation, configuration, execution, outputs, and troubleshooting. General users normally run with `consortium_run: false` when deliberately departing from the consortium resource contract.

### Consortium users

Use [`docs/users/consortium/`](docs/users/consortium/README.md) for controlled, self-contained operating documentation. Consortium sites should follow the maintained SOPs and IO definitions for harmonized runs.

### Maintainers

Release policy, qualification-baseline maintenance, and technical debt are documented under [`docs/maintainers/`](docs/maintainers/README.md).

## Supported core workflow

- `full_length + paired`
- `quantseq + single`
- library-specific UMI handling with a fixed contiguous 5′ UMI on R1 or R2 and optional adjacent discard bases
- GDC-aligned STAR processing using GRCh38.d1.vd1, GENCODE v36, and the maintained GDC STAR 2.7.5c index for consortium operation
- per-library canonical expression tables, run-level count matrices, standardized QC, and provenance

Optional fusion, TE, ASE, and RSeQC modules are not part of the supported consortium core unless explicitly stated otherwise in the applicable release documentation.

## Quick start

Start from the maintained run templates:

```bash
mkdir -p /path/to/run
cp templates/config.yaml /path/to/run/config.yaml
cp templates/samples.tsv /path/to/run/samples.tsv
```

Validate the sample sheet:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
```

Perform a dry run:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  -n
```

For installation and execution details, continue with the appropriate documentation branch above.

## Repository-maintained operational components

- controller environment: [`environments/controller.yaml`](environments/controller.yaml)
- run templates: [`templates/`](templates/README.md)
- execution-profile templates: [`templates/profiles/`](templates/profiles/README.md)
- software/container manifest: `workflow/config/software_versions.yaml`
- GDC resource installation: [`resources/`](resources/README.md)
- read-only site installation/run preflight: `scripts/verify_installation.py`
- maintainer repository consistency check: `tests/release/check_release_consistency.py`
- qualification tests: [`tests/`](tests/README.md)

Pipeline release identity is maintained in `workflow/release.yaml`; copied run configurations do not control it.
