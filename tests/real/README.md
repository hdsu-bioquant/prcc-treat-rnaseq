# Realistic public-data qualification

This is the production-style qualification exercise for pRCC-RNA-Seq. It complements the
fully synthetic regression test rather than replacing it.

The qualification has three purposes at the same time:

1. run the workflow on small **real human RNA-seq data** using the maintained GDC reference bundle;
2. test cross-run/cross-site harmonization of the deliberate deterministic output subset; and
3. train/check the operational workflow partners will later use for consortium datasets.

For that reason, setup is intentionally **documentation-led**. There is no script that creates a
Snakemake profile or silently prepares a run directory. A partner manually creates/reuses the
site execution profile, manually copies the maintained qualification run files into a fresh run
directory, edits only the documented path fields, and invokes Snakemake directly.

Automation is retained where it improves reproducibility rather than hiding the operating model:

- `get_test_data.sh` obtains/verifies the exact public FASTQ bytes;
- `workflow/scripts/sample_sheet.py` validates the copied library sheet;
- `validate_results.py` checks qualification inputs, outputs, provenance, checksums and the frozen
  realistic baseline.

## Qualification libraries

| library | public run | route | role |
|---|---|---|---|
| `REAL_FL` | `SRR493372` | full-length, paired-end, no UMI | realistic human full-length/GDC path |
| `REAL_QS_UMI` | `SRR16932032` | QuantSeq FWD, single-end, UMI | real 6-nt UMI + 4-nt post-UMI spacer path |

The exact compressed ENA FASTQs are part of the test definition. `data_manifest.tsv` pins the
URL, filename, compressed byte size and MD5 for all three files. Qualification does not use SRA
Toolkit and does not dynamically discover accessions.

`SRR16932032` is 101 nt/read. In the pinned FASTQ the first six bases are the UMI and bases 7-10
are the observed `TATA` spacer. The maintained `samples.tsv` describes this semantically as:

```text
umi_pattern       NNNNNN
umi_location      read1_start
umi_discard_bases 4
```

The pipeline therefore removes four adjacent non-UMI bases without hard-coding the spacer
sequence. Exact raw-to-extracted UMI transformation remains regression-tested by the synthetic
fixture; this realistic test checks the exact public input architecture plus production-like
end-to-end behavior.

---

## Operational ownership model

Keep three lifetimes separate:

| ownership/lifetime | examples | expected handling |
|---|---|---|
| Maintainer-owned | repository, `templates/`, `tests/real/`, maintained baselines | versioned and changed through normal repository maintenance |
| Site/user-owned, persistent | installed GDC bundle, copied Snakemake profile | configure once and reuse |
| Run-owned | copied `config.yaml`, copied `samples.tsv`, output, run logs | create fresh for each dataset/run |

The realistic qualification exercises this operating model. Consortium SOPs can use the same
ownership boundaries with `templates/config.yaml` and `templates/samples.tsv` for restricted datasets.

---

## 0. Prerequisites

Before starting, the site should have:

- the maintained Snakemake controller environment from `environments/controller.yaml` (or a site-equivalent environment that has been explicitly validated);
- Apptainer/Singularity available on the machines that execute jobs;
- the repository available on a filesystem visible to those jobs;
- the pinned default container images pre-pulled for the current checkout (`bash containers/pull_images.sh`);
- the exact GDC resource bundle installed and passing the fast structural verification;
- enough local/HPC resources for human STAR alignment.

From the repository root, a site installation can be checked structurally with:

```bash
bash resources/verify_gdc_references.sh /absolute/path/to/gdc
```

The canonical installed-reference SHA256 manifest is frozen in the repository. A consortium site
must qualify its installation against that manifest once before running:

```bash
bash resources/verify_gdc_references.sh --qualify /absolute/path/to/gdc
```

Normal consortium runs then check the small qualification stamp plus reference structure rather
than hashing the full human reference bundle every time.

---

## 1. Download the exact public test FASTQs once

From the repository root:

```bash
bash tests/real/get_test_data.sh
```

The script downloads to:

```text
tests/real/data/
├── SRR493372_1.fastq.gz
├── SRR493372_2.fastq.gz
└── SRR16932032.fastq.gz
```

For every file it requires all of the following:

- pinned ENA URL;
- exact compressed filename;
- exact compressed byte size;
- published MD5;
- successful `gzip -t`.

It resumes partial downloads, retries transient failures, skips already-valid files, and performs
one clean retry if a resumed file reaches completion but fails integrity verification.

At any later time, verify locally without network access:

```bash
bash tests/real/get_test_data.sh --verify-only
```

Do not substitute files obtained through SRA Toolkit or accession re-discovery. The qualification
is defined by these exact compressed bytes.

---

## 2. Configure the site execution profile once

This step is deliberately manual because the profile is site infrastructure, not run biology.
Choose the profile template that matches the environment and copy the **directory** once.

### SLURM example

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/slurm ~/.config/snakemake/prcc-rnaseq-slurm
```

Then edit:

```text
~/.config/snakemake/prcc-rnaseq-slurm/config.yaml
```

for the site: scheduler job limits, partition/account if required, latency, and any Apptainer bind
roots needed to expose the repository, test FASTQs, references, scratch and run output on compute
nodes. Prefer binding the canonical shared-filesystem root to the **same absolute path** inside the
container. Do not encode the location of one particular pipeline checkout with an Apptainer
`--pwd`; the profile should remain reusable across pipeline checkouts and runs.

### Local workstation example

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/local ~/.config/snakemake/prcc-rnaseq-local
```

Edit the copied profile's workstation-wide CPU/RAM limits and bind paths as needed.

Do **not** edit the maintained files under `templates/profiles/`. Once this profile successfully
runs the realistic qualification, retain it for later production runs on the same environment.
See `templates/profiles/README.md` for the full profile contract.

For the commands below, define the chosen profile path, for example:

```bash
PROFILE="$HOME/.config/snakemake/prcc-rnaseq-slurm"
```

---

## 3. Create a fresh qualification run directory

Choose a run-owned directory outside the maintained repository. For example:

```bash
# Resolve the repository through its physical/canonical filesystem path.
# This matters on HPC systems where a convenient login-node symlink may not
# exist under the same absolute name inside Apptainer on compute nodes.
REPO="$(realpath /absolute/path/to/pRCC-RNA-Seq)"
RUN_DIR=/absolute/path/to/pRCC-TREAT_realistic_qualification
GDC_DIR="$(realpath /absolute/path/to/gdc)"

mkdir -p "$RUN_DIR"
cp "$REPO/tests/real/config.yaml"  "$RUN_DIR/config.yaml"
cp "$REPO/tests/real/samples.tsv" "$RUN_DIR/samples.tsv"
```

These two files are the **maintained qualification run definition**. They are intentionally
separate from the generic `templates/config.yaml` and `templates/samples.tsv` used to start a new
consortium dataset, but follow the same interface and copy-edit-run workflow.

### Edit the copied `config.yaml`

Change only the explicitly marked site/run paths:

```yaml
samples: /absolute/path/to/pRCC-TREAT_realistic_qualification/samples.tsv
output: /absolute/path/to/pRCC-TREAT_realistic_qualification/output

reference:
  dir: /absolute/path/to/gdc
```

If the site needs dedicated temporary storage, it may also add an absolute `tmpdir:`. Keep
`consortium_run: true`, GDC filenames, STAR scientific parameters, trimming settings and optional
module settings unchanged for this qualification.

### Edit the copied `samples.tsv`

Change only the FASTQ path fields so they point to the exact files downloaded in Step 1. With the
recommended repository-local test data location these are:

```text
$REPO/tests/real/data/SRR493372_1.fastq.gz
$REPO/tests/real/data/SRR493372_2.fastq.gz
$REPO/tests/real/data/SRR16932032.fastq.gz
```

Use canonical absolute paths (for example, the output of `realpath`), not login-node-only
symlink aliases. The same absolute path must be visible on the execution node and inside
Apptainer. Do not change library IDs, accessions, assay/layout/strandedness, UMI fields or row
order.

This manual copy/edit step is intentional. It rehearses the same run-owned configuration workflow
used later for restricted consortium datasets.

---

## 4. Preflight the copied run

Run these commands from the repository's **physical** root:

```bash
cd -P "$REPO"
printf 'logical PWD: %s\nphysical PWD: ' "$PWD"
pwd -P
```

The displayed paths should be the same canonical shared-filesystem path. If they differ, use the
physical path reported by `pwd -P` for the repository, references, FASTQs and run configuration.

First re-check the exact public inputs without network access:

```bash
bash tests/real/get_test_data.sh --verify-only
```

Run the read-only installation preflight against the copied config and the persistent site profile:

```bash
python scripts/verify_installation.py \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE"
```

This includes the routine GDC structure/qualification-stamp check, controller-version checks, default
container version probes, and basic profile checks. Use `--full-reference-check` only when a full
canonical installed-reference SHA256 verification is required.

Validate the copied sample sheet:

```bash
python workflow/scripts/sample_sheet.py "$RUN_DIR/samples.tsv"
```

Then perform a profile-enabled dry-run using the actual site profile that will be retained for
production:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --dry-run
```

Resolve profile/bind/resource problems here. Do not compensate for infrastructure problems by
changing the scientific harmonization settings in the run configuration.

---

## 5. Run the realistic qualification directly with Snakemake

Still from the repository root:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --keep-going
```

This is intentionally the same invocation style used for a normal production run. There is no
qualification wrapper that creates or modifies the run/profile on the user's behalf.

The run should create:

```text
$RUN_DIR/output/
├── results/
├── restricted/
└── intermediate/
```

---

## 6. Validate the completed run

The realistic baseline is already maintainer-frozen. Ordinary partner/release qualification uses:

```bash
python tests/real/validate_results.py --run-dir "$RUN_DIR"
```

The validator checks, among other things:

- that the run-owned config/sample sheet conform to the qualification procedure;
- `consortium_run: true` and the exact GDC reference identity;
- exact byte size + MD5 + gzip integrity of the pinned qualification FASTQs;
- the observed QuantSeq UMI architecture and deterministic UMI-extraction QC;
- realistic non-empty raw gene counts and UMI molecule counts;
- canonical per-library expression semantics;
- stable QC tables and restricted BAM/FastQC products;
- portable run/reference/software provenance;
- MultiQC content;
- package checksums and the generated deterministic validation manifest; and
- byte-for-byte agreement with `tests/real/expected/validation_checksums.sha256`.

A successful run therefore demonstrates the pinned public qualification inputs, production-style
run ownership, the site's persistent execution profile, qualified GDC reference identity, realistic
full-length and QuantSeq+UMI processing, the portable/restricted output contract, and deterministic
agreement with the maintained baseline.

### Maintainer-only baseline changes

`--skip-frozen-baseline` is not part of normal partner qualification. Maintainers may use it only
during a deliberate, reviewed baseline-change exercise after an ordinary run has exposed an
understood expected difference. Preview or accept the generated realistic manifest with:

```bash
bash tests/maintainers/update_validation_baseline.sh realistic "$RUN_DIR"
bash tests/maintainers/update_validation_baseline.sh realistic "$RUN_DIR" --apply
```

An unchanged generated manifest needs no repeat run. A metadata/provenance-only difference may be
accepted after review. If deterministic computational products changed, establish byte-for-byte
reproducibility with two new clean runs of the finalized implementation before updating the
baseline. The full procedure is in `docs/maintainers/qualification-baselines.md`.

---

## Partner qualification summary

A new partner/site follows the operational exercise above:

1. install the maintained GDC bundle and run `verify_gdc_references.sh --qualify` once;
2. copy/configure a local or SLURM execution profile once and retain it;
3. download/verify the exact qualification FASTQs;
4. create a fresh qualification run directory;
5. copy `tests/real/config.yaml` + `samples.tsv` there and edit only documented paths;
6. validate/dry-run/execute directly with the persistent profile;
7. run `validate_results.py` **without** `--skip-frozen-baseline`.

Passing therefore checks both computational harmonization and compliance with the intended run
setup procedure.

---

## Tracked versus local files

Maintained in Git:

```text
tests/real/
├── README.md
├── config.yaml
├── samples.tsv
├── data_manifest.tsv
├── get_test_data.sh
├── validate_results.py
└── expected/
    ├── README.md
    └── validation_checksums.sha256   # maintainer-frozen realistic baseline
```

Local/generated and Git-ignored:

```text
tests/real/data/       # exact public FASTQs downloaded once
```

The actual qualification run directory should normally live outside the maintained repository.
The persistent site execution profile also normally lives outside the repository (for example under
`~/.config/snakemake/`).

The old draft `tests/real/manifest.tsv` is superseded by `data_manifest.tsv` and should be removed.
