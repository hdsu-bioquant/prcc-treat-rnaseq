# SOP-05 — Upgrade and requalification

| Field | Value |
|---|---|
| SOP ID | SOP-05 |
| Status | Draft |
| Document version | 0.2 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Purpose

Define the site procedure when adopting a new maintained pipeline release or when persistent site infrastructure changes in a way that can affect the previously qualified execution environment.

## Procedure

1. Obtain the exact target pipeline release/tag or release archive designated by the maintainers.
2. Review release notes and the consortium documentation applicable to that release.
3. Repeat the software-layer steps from SOP-01 as required: controller environment, container runtime, core containers, and synthetic qualification.
4. Apply any maintained GDC reference changes and requalify the installed reference bundle when required by the release.
5. Review or recreate the persistent local/SLURM execution profile when controller, scheduler, filesystem, bind-path, or execution-policy changes require it.
6. Repeat the relevant SOP-02 realistic qualification steps using a fresh run directory, including full installation/run preflight against the persistent profile.
7. Process new restricted consortium data only after required qualification passes.

A site does not update repository-maintained expected validation checksums. If a maintained release intentionally changes deterministic outputs, the corresponding reviewed baselines are supplied by the maintainers with that release.
