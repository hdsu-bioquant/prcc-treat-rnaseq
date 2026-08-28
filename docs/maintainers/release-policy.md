# Release, qualification, and requalification policy

This document defines repository-level release metadata and the validation expected when maintained
code, references, dependencies, provenance, or portable outputs change.

## Version convention

- Pipeline releases use Semantic Versioning (`MAJOR.MINOR.PATCH`).
- Git/repository release tags use the `v` prefix, for example `v1.0.0`.
- Pre-release candidates use the SemVer pre-release form, for example `v1.0.0-rc.1`.
- Machine-readable pipeline identity is maintainer-owned in `workflow/release.yaml` and is not a
  run configuration field.

## Release qualification

A release may be described as **qualified for the core expression workflow** only when all of the
following are true for the exact source revision and release metadata being distributed:

1. the supported core scope is documented;
2. the synthetic qualification passes;
3. the realistic public-data qualification passes against the maintained GDC reference bundle;
4. the canonical installed-reference manifest is present and the site used for realistic
   qualification is qualified against it;
5. the maintained deterministic synthetic and realistic baselines correspond to that release; and
6. release documentation and templates describe the same supported scope and operational model.

The supported core scope is `full_length + paired` and `quantseq + single`, with the maintained UMI
architecture. Optional fusion, TE, ASE, and other modules require their own validation before they
are included in the same qualification claim.

## Pipeline version versus output-contract version

The pipeline version identifies the released implementation and may change for code,
documentation, operational, dependency, or provenance changes.

The output-contract version changes only when the portable contract itself becomes incompatible or
materially changes, for example:

- adding, removing, or renaming required portable files;
- changing required columns or their meaning;
- changing the semantics of existing portable values; or
- changing the `results/`, `restricted/`, or `intermediate/` boundary in a way that affects the
  portable deliverable.

A pipeline version change does **not** automatically require an output-contract version change.
The output-contract version is stored in `workflow/release.yaml`.

## Validation levels for repository changes

The appropriate validation level is determined by what a change can affect. These levels are
default maintenance guidance; maintainers may run additional validation when warranted.

For Levels A and B, begin with one ordinary qualification against the maintained expected
baselines. An unchanged deterministic manifest is already evidence that the maintained expected
result reproduced; repeated clean runs are reserved for establishing a **new** computational
baseline when an intentional change alters deterministic computational outputs.

### Level A — scientific/output-affecting changes

Run synthetic and realistic qualification when changes can affect scientific results or their
deterministic generation, including changes to:

- scientific algorithms, command-line parameters, filtering, trimming, alignment, counting,
  normalization, UMI extraction/deduplication, or representative-selection policy;
- supported assay, layout, or UMI semantics;
- genome, annotation, STAR index, or their scientific interpretation;
- analysis container/tool versions that can affect scientific outputs;
- deterministic ordering/canonicalization logic; or
- portable expression/QC computation semantics.

If both qualifications reproduce their maintained deterministic baselines, no additional clean run
is required merely because the change was potentially output-affecting.

If an intentional change produces a new deterministic computational result, first investigate and
finalize the intended implementation. Then establish reproducibility of each affected qualification
baseline with two fresh clean runs whose generated validation manifests agree byte-for-byte before
updating the maintained expected checksums. See `qualification-baselines.md`.

### Level B — deterministic metadata/provenance changes

Run synthetic and realistic qualification when a change affects files included in
`results/run/validation_checksums.sha256` even if scientific values are not expected to change.
Confirm that observed differences are limited to the intended metadata, provenance, or contract
change before updating an affected baseline.

Examples include:

- changing the maintained controller environment (Python/Snakemake/executor-plugin constraints) when
  scientific computation is not expected to change;
- changing `pipeline_release` in `workflow/release.yaml` because release identity is recorded in the
  deterministic portable run metadata;
- changing stable portable run metadata included in the deterministic validation set; or
- intentional output-contract metadata changes that do not alter scientific computation.

A reviewed baseline difference limited to directly generated deterministic metadata/provenance does
not by itself require two computational runs. If the observed differences extend into computational
products or reveal a change in how those products are generated, treat the change as Level A for
baseline-update purposes.

### Level C — repository/static changes

For changes confined to prose documentation, comments, maintainer guides, static checks, citation
or licensing metadata, or code that is not on an enabled supported-core execution path, static and
syntax checks are normally sufficient.

If such a change alters generated portable files or execution behavior, use Level B or Level A as
appropriate.

## Baseline update policy

Maintained baselines are versioned expected outputs used to detect unintended drift. They may be
updated when an intentional repository change changes the expected deterministic result.

The normal sequence is:

1. run the applicable qualification once against the existing maintained baseline;
2. if the generated validation manifest is unchanged, accept the qualification result and stop;
3. if it differs, identify every changed deterministic file and establish whether the difference is
   intended;
4. retain the existing baseline for unintended or unexplained differences;
5. for an intended metadata/provenance-only difference, review and update the affected baseline;
6. for an intended computational-output difference, finalize the implementation and obtain
   byte-for-byte agreement from two fresh clean runs of each affected qualification before updating
   the affected baseline; and
7. rerun ordinary validation against the updated baseline.

The first failing/mismatching run is diagnostic. When two-run reproducibility is required, perform
the two clean runs after the intended implementation/configuration is finalized rather than counting
a pre-fix exploratory run as one of them.

Use `tests/maintainers/update_validation_baseline.sh` and the procedure in
`qualification-baselines.md` for deterministic output baselines. A failing qualification without an
understood intended change is a signal to investigate rather than a reason to update a baseline.

## Release-candidate procedure

For a release candidate intended to carry the core qualification claim:

1. run repository/static consistency checks;
2. set the intended release identity in `workflow/release.yaml`;
3. run ordinary synthetic qualification against the maintained baseline;
4. run ordinary realistic qualification against the maintained baseline;
5. investigate any deterministic differences and apply the baseline-update policy above only where
   the differences are understood and intended;
6. rerun ordinary qualification against any updated baselines; and
7. tag the exact reviewed source revision using the corresponding `v...` release tag.

A second clean workflow run is therefore not automatic for a release candidate: it is required when
an intentional computational change creates a new deterministic expected result that must be shown
to reproduce before its baseline is accepted.
