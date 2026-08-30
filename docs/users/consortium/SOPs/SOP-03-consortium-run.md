# SOP-03 — Consortium run

| Field | Value |
|---|---|
| SOP ID | SOP-03 |
| Status | Draft |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Run consortium RNA-seq libraries using the maintained harmonized processing contract.

## Preconditions

- site has passed SOP-02 for the applicable pipeline release/environment;
- run inputs satisfy the consortium sample-sheet contract;
- persistent site profile is available;
- GDC reference qualification stamp remains valid.

## Prepare the run

```bash
mkdir -p /path/to/run
cp templates/config.yaml /path/to/run/config.yaml
cp templates/samples.tsv /path/to/run/samples.tsv
```

Edit the copied files according to `../io-definitions/sample-sheet.md` and `../io-definitions/run-configuration.md`.

For consortium production, `consortium_run` must remain `true` and the maintained harmonization settings must not be changed.

Validate the sample sheet:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
```

Run preflight:

```bash
python scripts/verify_installation.py \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm
```

Perform a dry run:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm \
  -n
```

## Execute

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm \
  --keep-going
```

Site-specific execution resources may be adjusted in the copied profile. Scientific/harmonization settings belong to the maintained run contract and must remain unchanged unless a consortium-approved analysis deviation is explicitly recorded.

## Completion checks

- workflow exits successfully;
- expected library outputs exist for every input library;
- run-level matrices/QC/provenance are present;
- portable versus restricted data boundaries are reviewed before transfer.

Continue with SOP-04 for delivery.
