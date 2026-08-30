# Maintainer documentation

- [`release-policy.md`](release-policy.md) — versioning, qualification levels, and release/requalification policy.
- [`qualification-baselines.md`](qualification-baselines.md) — maintenance of deterministic validation and reference baselines.
- [`controlled-documentation.md`](controlled-documentation.md) — versioning/status/applicability rules for consortium SOPs and IO contracts.
- [`technical-debt.md`](technical-debt.md) — known maintenance leftovers that should be addressed deliberately rather than hidden in user documentation.

User-facing documentation lives under `docs/users/`. Consortium SOPs and IO definitions are controlled separately from the conventional general-user guide. Changelog maintenance is defined as part of the release process in `release-policy.md`.

## Root release metadata

- `CHANGELOG.md` — notable release-level changes and release history;
- `CITATION.cff` — machine-readable software citation metadata;
- `LICENSE` — repository software license.

These files are reviewed as release-facing metadata and are checked for basic consistency by the repository consistency check.

## Repository consistency check

Before release preparation and after repository-maintenance changes, run:

```bash
python tests/release/check_release_consistency.py
```

This is a **maintainer-facing static repository check**. It verifies maintained release metadata, controller/software declarations, frozen qualification artifacts, supported-core template agreement, and controlled-documentation structure. It does not inspect a site's installed references, containers, or execution environment and does not execute the scientific workflow.

By contrast, `scripts/verify_installation.py` is a **site/user-facing runtime preflight**. It checks an actual installation (controller/runtime, local containers/tool versions, references, and optionally an execution profile) before realistic qualification or production execution.
