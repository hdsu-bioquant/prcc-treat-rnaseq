# SOP-01 — Site installation

| Field | Value |
|---|---|
| SOP ID | SOP-01 |
| Status | Draft |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Establish the persistent software, reference, container, and execution-profile installation required before consortium qualification or production runs.

## Required persistent components

- repository checkout for the applicable pipeline release;
- Conda controller environment from `environments/controller.yaml`;
- Apptainer or Singularity;
- maintained GDC reference installation;
- maintained default/core containers;
- copied local or SLURM execution profile appropriate for the site.

## Procedure

### 1. Create the controller environment

```bash
conda env create -f environments/controller.yaml
conda activate prcc-rnaseq-controller
python --version
snakemake --version
```

### 2. Make Apptainer/Singularity available

Confirm that either `apptainer` or `singularity` is on `PATH` and reports its version.

### 3. Install and qualify the GDC reference bundle

```bash
bash resources/get_gdc_references.sh /path/to/gdc
bash resources/verify_gdc_references.sh --qualify /path/to/gdc
```

The consortium reference contract is defined in `../io-definitions/run-configuration.md`. Sites must use the maintained extracted GDC bundle and pre-built STAR index; do not rebuild the production STAR index locally.

### 4. Obtain the maintained default containers

```bash
bash containers/pull_images.sh
```

The maintained software manifest is `workflow/config/software_versions.yaml`. Required default containers must pass their declared version probes. Exact SIF byte identity is not a consortium requirement.

### 5. Create a persistent execution profile

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/slurm ~/.config/snakemake/prcc-rnaseq-slurm
# or, for a workstation:
cp -r templates/profiles/local ~/.config/snakemake/prcc-rnaseq-local
```

Edit the copied profile for site scheduler/resources/bind paths. Do not place site scheduler settings in run-owned `config.yaml`.

### 6. Run installation preflight

After preparing a run config for qualification, execute:

```bash
python scripts/verify_installation.py \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm
```

The preflight is read-only and must not be used as an installer or repair mechanism.

## Acceptance criteria

- maintained controller versions satisfy the preflight;
- Apptainer/Singularity is available;
- GDC reference structure and qualification stamp are valid;
- all required core containers are present and pass expected version probes;
- copied execution profile is readable and appropriate for the site.

Proceed to SOP-02 only after these criteria are met.
