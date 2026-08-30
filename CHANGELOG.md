# Changelog

All notable changes to pRCC-RNA-Seq will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and release versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

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

- Pipeline release identity is maintained by the repository rather than copied run configuration.
- Qualification-baseline maintenance is explicit and separated from synthetic fixture integrity.
- Consortium onboarding is organized as a lightweight first-line software installation followed by reference/resource and realistic site qualification.

### Notes

- Optional fusion, TE, ASE, and RSeQC modules are not part of the supported consortium core unless explicitly stated for a release.
- The first release-candidate version and release date will be assigned when the release candidate is cut.
