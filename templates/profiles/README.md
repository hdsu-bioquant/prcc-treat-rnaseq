# pRCC-RNA-Seq execution-profile templates

These profiles configure **where and with what infrastructure resources Snakemake runs**.
They are intentionally separate from the biological run interface in `templates/config.yaml`
and `templates/samples.tsv`.

Choose one starting point:

- `local/config.yaml` — one Linux workstation using local CPU/RAM plus Apptainer/Singularity;
- `slurm/config.yaml` — a shared-filesystem SLURM cluster using the Snakemake SLURM executor plugin plus Apptainer/Singularity.

A profile is normally copied **once per user or compute environment** and reused for many
runs. It is not copied into the portable `results/` package and is not part of the scientific
harmonization contract.

## 1. Copy a profile; do not edit the maintained template

A convenient Snakemake user-profile location on Linux is `~/.config/snakemake/`.
For a workstation:

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/local ~/.config/snakemake/prcc-rnaseq-local
```

For SLURM:

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/slurm ~/.config/snakemake/prcc-rnaseq-slurm
```

You may instead keep the copied profile in a site or project directory and pass its path to
`--profile`.

## 2. Settings shared by both profiles

Both templates enable:

- Apptainer/Singularity software deployment;
- `rerun-incomplete: true`;
- printed shell commands;
- `rerun-triggers: mtime`;
- conservative fallback `mem_mb` / `runtime` values for lightweight rules.

The workflow itself declares resources for compute-heavy rules. In particular, current human
STAR alignment jobs request **64,000 MB** and use `star.threads` from the run configuration.
Changing execution limits does not change the GDC STAR parameter recipe.

### Apptainer bind paths

Apptainer normally exposes the current working directory and the user's home directory. If
FASTQs, references, output, scratch, or the repository live on other filesystem roots, add
those roots to `apptainer-args` in the copied profile, for example:

```yaml
apptainer-args: "--bind /shared/data,/scratch,/reference"
```

Use paths that are valid on the machine(s) where jobs actually execute. On SLURM that means
the paths must be visible from compute nodes, not only the login node.

## 3. Local workstation profile

Edit `cores` to the number of CPU cores Snakemake may use concurrently. Edit
`resources.mem_mb` to the amount of RAM, in MB, that you are willing to make available to the
workflow scheduler.

Example for a 32-core / 96-GB workstation:

```yaml
cores: 16
resources:
  mem_mb: 90000
```

`resources.mem_mb` is important locally: it prevents Snakemake from starting combinations of
jobs whose declared memory reservations exceed the workstation-wide limit. Because STAR jobs
currently request 64,000 MB, the local limit must be at least 64,000 MB for a production human
run. In practice, leaving RAM for the OS and Snakemake controller is advisable; an 80–96 GB or
larger workstation is the practical target for this workflow.

Also ensure that `cores` is not lower than `star.threads` in the run-specific `config.yaml`.
The production template currently uses `star.threads: 8`.

Launch from the repository root:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-local \
  --keep-going
```

Command-line options override the copied profile, so a one-off test can temporarily use fewer
cores without editing it, for example `--cores 8`.

## 4. SLURM profile

The SLURM profile requires Snakemake >=8.6 and the executor plugin in the environment that
runs the Snakemake controller:

```bash
conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

The template's `jobs` value limits how many scheduler jobs Snakemake may have active/submitted
at once. Tune it to local policy.

Most rules already declare `threads`, `mem_mb`, and `runtime`; the SLURM executor uses these as
per-job scheduler requests. The profile's `default-resources` are only fallbacks for rules that
do not declare a value.

If your cluster requires a partition and/or account, uncomment them in the copied profile:

```yaml
default-resources:
  mem_mb: 4000
  runtime: 60
  slurm_partition: "compute"
  slurm_account: "project_account"
```

Do not invent these names: use the values supplied by your HPC administrators. A site with
special high-memory or QoS requirements can add standard Snakemake `set-resources` overrides
to its copied profile without modifying workflow rules.

Launch from the repository root:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm \
  --keep-going
```

The Snakemake controller should remain running for the duration of the workflow. Use `tmux`, a
site-supported persistent login/session mechanism, or submit the controller itself as an HPC
job if required by local policy.

## 5. Validate a copied profile before production

First validate the biological run files independently:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
snakemake --snakefile workflow/Snakefile --configfile /path/to/run/config.yaml -n
```

Then check the profile with a profile-enabled dry-run:

```bash
# local
snakemake --snakefile workflow/Snakefile --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-local -n

# SLURM
snakemake --snakefile workflow/Snakefile --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-slurm -n
```

Before analyzing consortium data on a new installation, run the frozen synthetic qualification
test described in `tests/synthetic/README.md`. That test is intentionally independent of these
production profiles: every partner runs the same fixed two-core synthetic workflow. The site's
actual copied local/SLURM profile is exercised later by the realistic `tests/real/`
production-style qualification.

## 6. Resource expectations

Execution-profile resources describe machine/scheduler capacity; rule resources describe the
requirements of individual jobs. Do not lower a rule's production memory request merely to fit
a smaller machine.

Current supported planning targets are:

| Use case | CPU | Memory / scheduling expectation |
|---|---:|---|
| Frozen synthetic qualification | fixed at 2 cores | tiny synthetic reference; a 32-GB workstation is comfortably suitable for testing |
| Human GDC run on one workstation | at least 4 usable cores | must be able to schedule a 64,000-MB STAR job |
| Recommended human workstation | 8+ usable cores | approximately 80–96 GB+ physical RAM to leave operating-system headroom |
| Human GDC run on SLURM | site-dependent | selected worker nodes/partition must satisfy the workflow's per-job resource requests |

The synthetic row is an installation/regression target, **not** a statement that 32 GB is
sufficient for production human alignment. Precise production peak-RSS/runtime values should
be refined from the realistic integration test and real consortium runs rather than inferred
from the miniature fixture.

## 7. What belongs where

Keep these concepts separate:

| File/location | Purpose | Typical lifetime |
|---|---|---|
| `templates/config.yaml` → copied run config | biological/scientific run settings + run paths | one per run |
| `templates/samples.tsv` → copied sample sheet | sequencing-library metadata + FASTQ paths | one per run |
| `templates/profiles/local/` | official workstation profile template | copied once per workstation/user |
| `templates/profiles/slurm/` | official SLURM profile template | copied once per site/user |
| `~/.config/snakemake/prcc-rnaseq-*` (or another external path) | actual site/user profile copy | reused across runs; not tracked in this repository |

Do not put patient/library metadata, FASTQ lists, or scientific STAR parameters into an
execution profile. Conversely, scheduler partition/account names and site-specific Apptainer
bind paths do not belong in the biological run configuration.
