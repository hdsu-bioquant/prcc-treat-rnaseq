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

## Operational model being qualified

Keep three lifetimes separate:

| ownership/lifetime | examples | expected handling |
|---|---|---|
| Maintainer-owned | repository, `templates/`, `tests/real/`, frozen baselines | do not edit in place |
| Site/user-owned, persistent | installed GDC bundle, copied Snakemake profile | configure once and reuse |
| Run-owned | copied `config.yaml`, copied `samples.tsv`, output, run logs | create fresh for each dataset/run |

The realistic qualification is intended to be the first guided exercise in that operating model.
Later consortium SOPs can use the same steps with `templates/config.yaml` and
`templates/samples.tsv` for real restricted datasets.

---

## 0. Prerequisites

Before starting, the site should have:

- a working Snakemake controller environment (consortium-tested version: 9.19.0);
- Apptainer/Singularity available on the machines that execute jobs;
- the repository available on a filesystem visible to those jobs;
- the exact GDC resource bundle installed and passing the fast structural verification;
- enough local/HPC resources for human STAR alignment.

From the repository root, a current development installation can be checked with:

```bash
bash resources/verify_gdc_references.sh /absolute/path/to/gdc
```

During initial maintainer development, the canonical installed-reference SHA256 manifest is
intentionally not frozen yet. Once it exists in a release, a consortium site must qualify its
installation once before running:

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
nodes.

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
REPO=/absolute/path/to/pRCC-RNA-Seq
RUN_DIR=/absolute/path/to/pRCC-TREAT_realistic_qualification
GDC_DIR=/absolute/path/to/gdc

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

Use absolute paths. Do not change library IDs, accessions, assay/layout/strandedness, UMI fields or
row order.

This manual copy/edit step is intentional. It rehearses the same run-owned configuration workflow
used later for restricted consortium datasets.

---

## 4. Preflight the copied run

Run these commands from the repository root:

```bash
cd "$REPO"
```

First re-check the exact public inputs without network access:

```bash
bash tests/real/get_test_data.sh --verify-only
```

Check the installed GDC structure:

```bash
bash resources/verify_gdc_references.sh "$GDC_DIR"
```

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

### During initial maintainer qualification

Before the realistic frozen baseline has been created:

```bash
python tests/real/validate_results.py \
  --run-dir "$RUN_DIR" \
  --skip-frozen-baseline
```

This skips **only** comparison with the not-yet-frozen realistic validation hashes. It still checks:

- that the run-owned config/sample sheet conform to the qualification procedure;
- `consortium_run: true` and the exact GDC reference identity;
- exact byte size + MD5 + gzip integrity of the FASTQs referenced by the copied sample sheet;
- the observed 101-nt QuantSeq reads and 6-nt UMI + 4-nt `TATA` architecture;
- realistic non-empty raw gene counts and UMI molecule counts;
- canonical per-library expression semantics;
- stable QC tables and restricted BAM/FastQC products;
- portable run/reference/software provenance;
- MultiQC content;
- package checksums and the generated deterministic validation manifest.

No arbitrary mapping-rate pass/fail threshold is imposed at this stage. The exact qualified outputs
will be captured by the deterministic baseline after reproducibility is demonstrated.

### After the frozen baseline exists

Partner/release qualification uses the same command **without** the development option:

```bash
python tests/real/validate_results.py --run-dir "$RUN_DIR"
```

---

## 7. Maintainer-only: demonstrate two clean deterministic runs

This section is for freezing the first realistic baseline, not for ordinary partner operation.

After a successful first run with `--skip-frozen-baseline`, preserve its generated deterministic
manifest outside `output/`:

```bash
cp "$RUN_DIR/output/results/run/validation_checksums.sha256" \
   "$RUN_DIR/validation_checksums.run1.sha256"
```

Delete the first run's output so the second execution is clean while keeping the same run-owned
config/sample sheet and the same persistent site profile:

```bash
rm -rf "$RUN_DIR/output"
```

Repeat Steps 4-6, again using `--skip-frozen-baseline` for validation.

Then require byte-for-byte equality of the deterministic validation manifests:

```bash
cmp "$RUN_DIR/validation_checksums.run1.sha256" \
    "$RUN_DIR/output/results/run/validation_checksums.sha256" \
  && echo "PASS: realistic deterministic validation manifest reproduced"
```

If the manifests differ, investigate. Do not immediately replace a baseline.

---

## 8. Maintainer-only: freeze the realistic output and GDC reference baselines

Only after the two clean runs above reproduce exactly:

```bash
cp "$RUN_DIR/output/results/run/validation_checksums.sha256" \
   "$REPO/tests/real/expected/validation_checksums.sha256"
```

The exact GDC installation used for those successful realistic runs has now earned reference
qualification. Generate the **maintainer-owned** installed-reference manifest once:

```bash
cd "$REPO"
bash resources/verify_gdc_references.sh \
  --write-canonical-manifest "$GDC_DIR"
```

Then qualify the current site installation against the newly frozen canonical manifest:

```bash
bash resources/verify_gdc_references.sh --qualify "$GDC_DIR"
```

Commit/freeze `resources/gdc_installed_reference.sha256` only after this sequence. Partners should
never generate their own canonical baseline; they qualify their installation against the shared
maintainer file.

---

## 9. Final normal qualification after freezing

Remove the previous output and run once more using the same copied run files and persistent profile:

```bash
rm -rf "$RUN_DIR/output"

cd "$REPO"
snakemake \
  --snakefile workflow/Snakefile \
  --configfile "$RUN_DIR/config.yaml" \
  --profile "$PROFILE" \
  --keep-going

python tests/real/validate_results.py --run-dir "$RUN_DIR"
```

At this point a successful run proves together:

- exact public qualification input bytes;
- production-style run-owned configuration handling;
- the site's persistent execution profile;
- the exact consortium GDC reference identity and qualification stamp;
- realistic FL and QuantSeq+UMI processing;
- the portable/restricted output contract; and
- byte-for-byte agreement of the deterministic canonical result subset.

---

## What partners do after the baselines are frozen

A new partner/site follows the same operational exercise, except the maintainer-only freeze steps
are omitted:

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
    └── validation_checksums.sha256   # added only after maintainer qualification
```

Local/generated and Git-ignored:

```text
tests/real/data/       # exact public FASTQs downloaded once
```

The actual qualification run directory should normally live outside the maintained repository.
The persistent site execution profile also normally lives outside the repository (for example under
`~/.config/snakemake/`).

The old draft `tests/real/manifest.tsv` is superseded by `data_manifest.tsv` and should be removed.
