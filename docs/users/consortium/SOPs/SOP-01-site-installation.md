# SOP-01 — Site software installation

| Field | Value |
|---|---|
| SOP ID | SOP-01 |
| Status | Draft |
| Document version | 0.2 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Establish and verify the lightweight software layer required before consortium-specific reference installation and realistic site qualification.

Completion of this SOP demonstrates that the designated pipeline release, controller environment, container runtime, and maintained core containers can execute the synthetic workflow successfully. The larger GDC reference bundle and production-style execution profile are introduced in SOP-02.

## Procedure

### 1. Obtain the designated pipeline release

Use the exact pipeline release designated for the consortium exercise. For a published release or release candidate, obtain either:

- the corresponding release archive supplied by the maintainers; or
- a Git checkout explicitly switched to the designated release tag.

Do not perform consortium qualification or production from a moving development branch unless the maintainers have explicitly designated that exact revision for a pilot.

After obtaining the repository, enter its root directory and confirm the maintained release identity:

```bash
cat workflow/release.yaml
```

The `pipeline_release` value must match the release or release candidate communicated by the maintainers.

### 2. Create the maintained controller environment

```bash
conda env create -f environments/controller.yaml
conda activate prcc-rnaseq-controller
```

Confirm the controller tools are available:

```bash
python --version
snakemake --version
```

The maintained environment currently uses Python 3.13.x and Snakemake 9.19.0. The exact Python patch version may vary within the maintained minor line.

### 3. Make Apptainer or Singularity available

Confirm that either runtime is on `PATH` and functional:

```bash
apptainer --version
```

or:

```bash
singularity --version
```

Only one compatible runtime is required.

### 4. Obtain and verify the maintained core containers

From the repository root:

```bash
bash containers/pull_images.sh
```

The helper uses `workflow/config/software_versions.yaml` and runs the declared version probes for maintained default/core containers. Exact SIF byte identity is not a consortium requirement; successful tool-version probes and subsequent qualification are the relevant conformity checks.

### 5. Run the synthetic qualification

```bash
bash tests/synthetic/run_test.sh
```

The synthetic test is intentionally small and does not require the production GDC reference bundle or a site production profile. It exercises the maintained core workflow routes and compares deterministic outputs with the repository-maintained synthetic baseline.

## Acceptance criteria

SOP-01 is complete when:

- the designated repository release is available at the site;
- the maintained controller environment can be activated;
- Python and Snakemake report the expected maintained versions;
- Apptainer or Singularity is available;
- required core containers are present and pass their declared version probes; and
- `bash tests/synthetic/run_test.sh` completes successfully against the maintained baseline.

This establishes the first-line software installation only. Continue with **SOP-02 — Site qualification** to install and qualify the GDC resources, create the persistent local or SLURM execution profile, rehearse run configuration, perform installation preflight, and complete realistic production-style qualification.
