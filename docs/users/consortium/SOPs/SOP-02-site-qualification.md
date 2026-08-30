# SOP-02 — Site qualification

| Field | Value |
|---|---|
| SOP ID | SOP-02 |
| Status | Draft |
| Document version | 0.2 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Complete the consortium-specific site setup and demonstrate that the site's production-style installation reproduces the maintained realistic public-data baseline.

This SOP deliberately teaches the same ownership and execution model used later for consortium production: persistent GDC resources and execution profile, run-owned configuration and sample sheet, explicit preflight, direct Snakemake execution, and maintained deterministic validation.

## Preconditions

- SOP-01 completed successfully for the applicable pipeline release;
- enough local/HPC storage and compute capacity for human STAR alignment;
- the repository is visible from the machines that will execute jobs.

## Procedure

### 1. Install and qualify the maintained GDC reference bundle

Choose the persistent site location for the maintained GDC resources and install them from the repository root:

```bash
GDC_DIR=/absolute/path/to/gdc
bash resources/get_gdc_references.sh "$GDC_DIR"
```

Then qualify the extracted installation against the repository-maintained canonical manifest:

```bash
bash resources/verify_gdc_references.sh --qualify "$GDC_DIR"
```

Consortium operation uses the maintained extracted GDC bundle and pre-built STAR index. Do not rebuild the production STAR index locally.

### 2. Obtain the exact realistic qualification FASTQs

```bash
bash tests/real/get_test_data.sh
```

The qualification uses the exact pinned public FASTQ bytes defined under `tests/real/`. Do not substitute files obtained through another retrieval route.

At any later time the local copies can be rechecked without network access:

```bash
bash tests/real/get_test_data.sh --verify-only
```

### 3. Create the persistent site execution profile

Choose the profile matching the production execution environment.

For SLURM:

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/slurm ~/.config/snakemake/prcc-rnaseq-slurm
PROFILE="$HOME/.config/snakemake/prcc-rnaseq-slurm"
```

For local/workstation execution:

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/local ~/.config/snakemake/prcc-rnaseq-local
PROFILE="$HOME/.config/snakemake/prcc-rnaseq-local"
```

Edit the copied profile for site execution settings such as resource limits, scheduler/account settings where applicable, latency, and Apptainer/Singularity bind roots. Do not edit the maintained templates in place and do not put scheduler/site infrastructure settings into run-owned `config.yaml`.

The profile that passes this realistic qualification should be retained and reused for consortium production on the same execution environment.

### 4. Create a fresh realistic qualification run directory

Resolve the physical repository and GDC paths and create a run-owned directory outside the maintained repository:

```bash
REPO="$(realpath /absolute/path/to/pRCC-RNA-Seq)"
GDC_DIR="$(realpath /absolute/path/to/gdc)"
RUN_DIR=/absolute/path/to/pRCC-TREAT_realistic_qualification

mkdir -p "$RUN_DIR"
cp "$REPO/tests/real/config.yaml"  "$RUN_DIR/config.yaml"
cp "$REPO/tests/real/samples.tsv" "$RUN_DIR/samples.tsv"
```

Edit only the site/run paths marked for editing in the copied `config.yaml`:

```yaml
samples: /absolute/path/to/pRCC-TREAT_realistic_qualification/samples.tsv
output: /absolute/path/to/pRCC-TREAT_realistic_qualification/output

reference:
  dir: /absolute/path/to/gdc
```

An absolute `tmpdir:` may also be added when dedicated temporary storage is required. Keep `consortium_run: true` and do not alter the maintained GDC filenames, scientific STAR parameters, assay-processing settings, or optional-module settings for this qualification.

In the copied `samples.tsv`, change only the FASTQ path fields to the exact files downloaded by `tests/real/get_test_data.sh`. Use canonical absolute paths visible from execution nodes and inside the container runtime. Do not change library identifiers, accessions, assay/layout/strandedness, UMI fields, or row order.

This manual copy/edit process intentionally rehearses the same run-owned configuration model used for production data.

### 5. Run the installation/run preflight

From the physical repository root:

```bash
cd -P "$REPO"
```

Run the read-only preflight against the copied qualification config and the persistent site profile:

```bash
python scripts/verify_installation.py \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE"
```

The preflight checks controller/runtime conformity, required core containers and tool-version probes, GDC reference structure/qualification-stamp identity, and basic execution-profile availability. It does not install or repair anything.

Use `--full-reference-check` only when a complete canonical installed-reference SHA256 verification is specifically required.

### 6. Validate inputs and perform a profile-enabled dry run

Recheck the pinned public inputs:

```bash
bash tests/real/get_test_data.sh --verify-only
```

Validate the copied sample sheet:

```bash
python workflow/scripts/sample_sheet.py "$RUN_DIR/samples.tsv"
```

Then dry-run with the same persistent profile intended for production:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --dry-run
```

Resolve profile, bind-path, scheduler, filesystem, or resource problems here rather than changing maintained harmonization settings.

### 7. Execute the realistic qualification

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --keep-going
```

The expected output root is:

```text
$RUN_DIR/output/
├── results/
├── restricted/
└── intermediate/
```

### 8. Validate against the maintained realistic baseline

```bash
python tests/real/validate_results.py --run-dir "$RUN_DIR"
```

Acceptance criterion: the validator completes successfully and the generated deterministic validation manifest agrees with `tests/real/expected/validation_checksums.sha256`.

Consortium sites must not update maintained expected baselines. A mismatch must be investigated and reported to the maintainers.

## Acceptance criteria

SOP-02 is complete when:

- the GDC reference installation has a valid qualification stamp for the maintained canonical manifest;
- the exact realistic qualification inputs pass integrity checks;
- a persistent local or SLURM profile has been configured for the actual production environment;
- the copied qualification run passes installation preflight and sample-sheet validation;
- the profile-enabled Snakemake execution completes successfully; and
- `tests/real/validate_results.py` passes against the maintained realistic baseline.

Only after SOP-02 passes should the site process restricted consortium samples using SOP-03.

## Technical companion documentation

`tests/real/README.md` documents the realistic fixture and validator in additional implementation detail. The normative consortium procedure is this SOP.
