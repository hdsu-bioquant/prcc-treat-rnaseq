# General user guide

This guide is for users running pRCC-RNA-Seq outside a controlled consortium workflow. It describes the maintained pipeline interface without imposing consortium-specific harmonization requirements.

## Guide

1. [`installation.md`](installation.md) — controller environment, containers, references, and installation preflight.
2. [`configuration.md`](configuration.md) — sample sheet, run configuration, supported assay/layout combinations, and optional modules.
3. [`running-the-pipeline.md`](running-the-pipeline.md) — dry runs, local/SLURM execution, and monitoring.
4. [`outputs.md`](outputs.md) — output layout, canonical expression products, QC, and provenance.
5. [`troubleshooting.md`](troubleshooting.md) — common installation and execution failures.

The maintained templates remain the executable interface:

- `templates/config.yaml`
- `templates/samples.tsv`
- `templates/profiles/local/`
- `templates/profiles/slurm/`

For harmonized consortium operation, use the separate [`../consortium/`](../consortium/README.md) documentation rather than this guide.
