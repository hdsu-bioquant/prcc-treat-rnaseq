# SOP-05 — Upgrade and requalification

| Field | Value |
|---|---|
| SOP ID | SOP-05 |
| Status | Draft |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Define the site procedure when adopting a new maintained pipeline release or when a persistent site installation changes in a way that can affect qualification.

## Procedure

1. Obtain the target repository release.
2. Review release notes and consortium documentation applicable to that release.
3. Update/recreate the controller environment if required by the release.
4. Rerun `bash containers/pull_images.sh` to obtain/verify maintained default containers.
5. Apply any maintained reference-installation changes and requalify references if required.
6. Review the copied persistent execution profile for compatibility with the target controller/Snakemake/plugin versions.
7. Run the read-only installation preflight.
8. Repeat SOP-02 synthetic and realistic qualification against the target release's maintained baselines.
9. Process new restricted consortium data only after qualification passes.

A site does not update repository-maintained expected validation checksums. If a maintained release intentionally changes deterministic outputs, the corresponding new baselines are supplied by maintainers with that release.
