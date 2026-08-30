# Running the pipeline

## Dry run

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  -n
```

## Execution profiles

Copy and edit a maintained profile once for the compute environment:

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/local ~/.config/snakemake/prcc-rnaseq-local
# or
cp -r templates/profiles/slurm ~/.config/snakemake/prcc-rnaseq-slurm
```

Run with the copied profile:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm \
  --keep-going
```

Profile settings are site/execution settings and should not be placed in the biological run configuration.

Before a production-style run, the read-only installation preflight can validate the copied config/profile against the local installation.

See `templates/profiles/README.md` for profile-specific details.
