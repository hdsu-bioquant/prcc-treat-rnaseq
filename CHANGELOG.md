# Changelog

All notable changes to pRCC-TREAT-RNA-Seq will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and release versions follow [Semantic Versioning](https://semver.org/). Maintainer rules for
updating this file are defined in [`docs/maintainers/release-policy.md`](docs/maintainers/release-policy.md#changelog-maintenance).

## [Unreleased]

## [1.0.0-rc.1] - 2026-08-31

### Added

- Assay-aware Snakemake workflow for the supported consortium core:
  - full-length paired-end RNA-seq;
  - QuantSeq single-end RNA-seq;
  - optional library-specific UMI handling with the maintained fixed 5' UMI model.
- GDC-aligned STAR two-pass processing with maintained GRCh38.d1.vd1 / GENCODE v36 resources.
- Canonical per-library gene-expression tables and run-level count matrices.
- GDC-style FPKM, FPKM-UQ, and TPM outputs for full-length libraries.
- Deterministic UMI molecule counting with maintained UMI-tools behavior.
- Standardized QC and portable `results/`, site-retained `restricted/`, and disposable `intermediate/` output tiers.
- Maintained controller environment and local/SLURM execution-profile templates.
- Read-only site installation preflight and explicit container version probes.
- Deterministic synthetic qualification and realistic public-data qualification against maintained baselines.
- Controlled consortium documentation with numbered SOPs and input/output definitions.
- General user documentation for non-consortium operation.
- Maintainer release, qualification-baseline, release-check, and technical-debt documentation.
- MIT license and machine-readable citation metadata.

### Changed

- Renamed the pipeline from `pRCC-RNA-Seq` to `pRCC-TREAT-RNA-Seq`; the maintained repository slug is `prcc-treat-rnaseq`.
- Updated the maintained controller to Snakemake 9.20.0, which restores minute-based interpretation of bare `runtime` values supplied through profiles/CLI.
- Removed compatibility handling for the unpublished temporary `umitools-1.1.6.sif` development filename; the maintained local image name is `umitools.sif`.
- Pipeline release identity is maintained by the repository rather than copied run configuration.
- Qualification-baseline maintenance is explicit and separated from synthetic fixture integrity.
- Consortium onboarding is organized as a lightweight first-line software installation followed by reference/resource and realistic site qualification.
- Refreshed the maintained synthetic and realistic deterministic validation baselines for the release-identity rename; only `run/config.yaml` changed and maintained computational/QC checksums remained unchanged.

