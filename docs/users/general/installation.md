# General installation

## Obtain the pipeline

Use a published release/archive when reproducibility matters, or clone the repository and explicitly check out the desired release/tag. Avoid treating a moving development branch as a stable analysis release.

After obtaining the repository, run subsequent commands from its root unless stated otherwise.

## Controller environment

The maintained Conda specification is `environments/controller.yaml`:

```bash
conda env create -f environments/controller.yaml
conda activate prcc-rnaseq-controller
```

The maintained environment uses Python 3.13.x, Snakemake 9.19.0, and the SLURM executor plugin. Sites using another scheduler may maintain an equivalent controller environment.

## Apptainer or Singularity

Make either `apptainer` or `singularity` available on `PATH`. Workflow software is executed in containers; STAR, samtools, HTSeq, UMI-tools, fastp, BBDuk, FastQC, and MultiQC do not need separate host installations.

## Containers

Pre-pull the maintained default images with:

```bash
bash containers/pull_images.sh
```

The maintained software manifest is `workflow/config/software_versions.yaml`. Default containers have version probes used by the pull helper and installation preflight. The pipeline does not require byte-identical SIF files between installations.

## References

For GDC-aligned processing, install the maintained bundle with:

```bash
bash resources/get_gdc_references.sh resources/gdc
```

General users may use `reference.mode: local` for deliberate custom-reference analyses. Such runs should set `consortium_run: false` and are outside the consortium reference contract.

## Installation preflight

After a run config and any persistent execution profile have been prepared, the read-only preflight can check controller versions, reference installation, required local containers/tool versions, and the selected profile:

```bash
python scripts/verify_installation.py \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm
```

Use `--full-reference-check` when a complete GDC installed-reference SHA256 verification is desired.
