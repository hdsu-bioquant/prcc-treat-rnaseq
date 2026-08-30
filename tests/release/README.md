# Release/static consistency checks

Run from the repository root:

```bash
python tests/release/check_release_consistency.py
```

`check_release_consistency.py` is a maintainer-facing **static repository check**. It verifies that repository-maintained declarations and controlled artifacts agree, including release metadata, controller constraints, core container version probes, production/realistic config invariants, frozen qualification artifacts, and documentation structure.

It does not inspect a site's installed GDC bundle, execute containers, validate a copied execution profile, or run the scientific workflow. Those tasks belong to site qualification and workflow tests.

Do not confuse it with:

```bash
python scripts/verify_installation.py ...
```

`verify_installation.py` is a user/site-facing **runtime preflight** for an actual installation. It checks the active controller/runtime, local core containers and version probes, GDC reference installation when supplied by config/path, and an optional local/SLURM profile.

The repository consistency check is part of release preparation and should also be run after maintenance changes that can introduce declaration/documentation drift.
