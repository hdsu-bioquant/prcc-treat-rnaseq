# SOP-03 — Consortium run

| Field | Value |
|---|---|
| SOP ID | SOP-03 |
| Status | Draft |
| Document version | 0.2 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Run consortium RNA-seq libraries using the maintained harmonized processing contract and the persistent site installation qualified under SOP-02.

## Preconditions

- site has passed SOP-02 for the applicable pipeline release/environment;
- the qualified GDC reference installation remains available;
- the persistent local or SLURM execution profile that passed qualification is available;
- run inputs satisfy the consortium sample-sheet contract.

## 1. Select the qualified execution profile

Reuse the profile that passed SOP-02.

For a SLURM site, for example:

```bash
PROFILE="$HOME/.config/snakemake/prcc-rnaseq-slurm"
```

For a local/workstation execution environment:

```bash
PROFILE="$HOME/.config/snakemake/prcc-rnaseq-local"
```

Do not switch production to a materially different site profile/environment without applying the upgrade/requalification guidance in SOP-05.

## 2. Prepare the run-owned configuration

```bash
RUN_DIR=/absolute/path/to/run
mkdir -p "$RUN_DIR"
cp templates/config.yaml "$RUN_DIR/config.yaml"
cp templates/samples.tsv "$RUN_DIR/samples.tsv"
```

Edit the copied files according to `../io-definitions/sample-sheet.md` and `../io-definitions/run-configuration.md`.

For consortium production, `consortium_run` must remain `true` and the maintained harmonization settings must not be changed.

Use canonical absolute paths where practical, especially for run configuration, FASTQs, references, output, and scratch locations that must be visible from compute nodes and containers.

## 3. Validate the run inputs

```bash
python workflow/scripts/sample_sheet.py "$RUN_DIR/samples.tsv"
```

## 4. Run preflight against the production run and qualified profile

```bash
python scripts/verify_installation.py \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE"
```

The preflight is read-only. Resolve installation/profile failures rather than compensating by changing harmonization settings.

## 5. Perform a dry run

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --dry-run
```

## 6. Execute

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --keep-going
```

Site-specific execution resources may be adjusted in the copied persistent profile. Scientific/harmonization settings belong to the maintained consortium run contract and must remain unchanged unless a consortium-approved analysis deviation is explicitly recorded.

## Completion checks

- workflow exits successfully;
- expected library outputs exist for every input library;
- run-level matrices, QC, and provenance are present;
- the run used the qualified pipeline release/reference/profile combination or an explicitly requalified replacement; and
- portable versus restricted data boundaries are reviewed before transfer.

Continue with SOP-04 for delivery.
