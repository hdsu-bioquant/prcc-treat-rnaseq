# Release-level static checks

Run from the repository root:

```bash
python tests/release/check_release_consistency.py
```

The check is intentionally lightweight. It verifies release metadata ownership, the maintained
controller/preflight specification, required default-container version probes, agreement of the
supported-core template and realistic qualification settings, presence of maintained qualification
artifacts, optional-module defaults, and selected stale-documentation regressions. It does not
replace synthetic or realistic computational qualification.
