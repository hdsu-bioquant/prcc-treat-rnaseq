# SOP-02 — Site qualification

| Field | Value |
|---|---|
| SOP ID | SOP-02 |
| Status | Draft |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Demonstrate that a site installation reproduces the maintained deterministic synthetic and realistic qualification baselines before restricted consortium data are processed.

## Preconditions

- SOP-01 completed;
- maintained GDC references qualified at the site;
- persistent execution profile prepared;
- default/core containers available and version-probed successfully.

## Procedure

### 1. Synthetic qualification

From the repository root:

```bash
bash tests/synthetic/run_test.sh
```

Acceptance criterion: the synthetic workflow and validator complete successfully against the maintained fixture and expected deterministic-output baseline.

### 2. Realistic public-data qualification

Follow `tests/real/README.md` to obtain the pinned public FASTQs, create a fresh realistic run directory, use the site's real persistent execution profile, run the workflow, and validate:

```bash
python tests/real/validate_results.py --run-dir "$RUN_DIR"
```

Acceptance criterion: the completed run matches the maintained realistic deterministic baseline.

## Qualification rule

A normal qualification requires one successful run against each maintained baseline. A second clean run is not required when the authoritative expected checksums are reproduced unchanged.

If a deterministic mismatch occurs, stop and investigate. Consortium sites must not update maintained expected baselines. Baseline changes are maintainer operations governed by `docs/maintainers/qualification-baselines.md`.

## Completion

Only after both synthetic and realistic qualification pass should the site process restricted consortium samples.
