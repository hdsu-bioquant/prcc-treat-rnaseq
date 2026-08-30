# SOP-04 — Results delivery

| Field | Value |
|---|---|
| SOP ID | SOP-04 |
| Status | Draft |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Define the standard package transferred from a consortium site after a successful production run.

## Portable package

The intended consortium transfer unit is the run's `results/` directory. Its maintained contract is defined in `../io-definitions/results-contract.md` and `../io-definitions/qc-and-run-metadata.md`.

Before transfer, confirm that `results/` contains the expected per-library expression tables, run-level matrices, QC summaries, and run/provenance metadata.

## Restricted and intermediate content

Do not transfer `restricted/` or `intermediate/` as part of the standard portable results package.

`restricted/` can contain original sample-sheet information, source-path information, and optional module outputs. It remains at the originating site unless separate data-sharing approval explicitly covers the relevant content.

`intermediate/` contains working files and is not part of the delivery contract.

## Delivery record

The receiving consortium process should retain enough information to associate the delivered `results/` package with the site, pipeline release/provenance, and run identity without requiring transfer of restricted source data.
