# Consortium IO definition — run configuration

| Field | Value |
|---|---|
| Document ID | IO-02 |
| Status | Pilot |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `1.0.0-rc.1` |
| Last revised | 2026-08-30 |

## Contract

Consortium production runs start from `templates/config.yaml` and keep `consortium_run: true`.

Run-owned fields normally edited for each run are:

- `samples` — path to the copied `samples.tsv`;
- `output` — run-wise output root;
- optional `tmpdir` — STAR scratch location;
- `reference.dir` — site location of the maintained GDC resource bundle when different from the repository-relative default;
- `star.threads` — execution/resource setting appropriate for the site.

Absolute paths are recommended for run-owned files and FASTQs.

## Maintained consortium settings

The following define the harmonized production contract and are not site-tuning parameters:

- `consortium_run: true`;
- `reference.mode: gdc`;
- genome `GRCh38.d1.vd1.fa`;
- annotation `gencode.v36.annotation.gtf`;
- STAR index `star-2.7.5c_GRCh38.d1.vd1_gencode.v36`;
- `reference.sjdb_overhang: 100`;
- maintained `star.gdc_params`;
- `full_length.trim_adapters: false`;
- `quantseq.bbduk_polyA: true`;
- optional modules disabled unless the applicable consortium release explicitly states otherwise.

The site reference installation must have a valid qualification stamp corresponding to the repository-maintained `resources/gdc_installed_reference.sha256` manifest.

## Pipeline release identity

Pipeline name, release identity, and output-contract version are maintainer-owned in `workflow/release.yaml`. They are not run-config fields and must not be added to copied configurations.

## Site execution settings

Scheduler/account/partition/resource limits/latency handling/Apptainer bind paths belong in the copied execution profile, not in run-owned biological configuration.
