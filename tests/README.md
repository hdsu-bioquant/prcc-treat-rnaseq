# pRCC-RNA-Seq qualification tests

> Consortium onboarding uses the synthetic fixture in SOP-01 and the realistic fixture in SOP-02. The controlled procedure is defined in `docs/users/consortium/SOPs/`; the files under `tests/` document and implement the qualification fixtures.


The repository contains two complementary qualification layers with deliberately different jobs.

## `tests/synthetic/` — deterministic workflow regression

The synthetic fixture is tiny, version-controlled, independent of the production GDC bundle, and
runs with a fixed local execution model. It exercises all currently implemented assay/layout/UMI
routes, exact UMI preprocessing, the portable/restricted output structure, biological expectations,
and the frozen deterministic output subset.

Consortium sites run this as the final acceptance step of SOP-01; maintainers also use it after relevant workflow changes:

```bash
bash tests/synthetic/run_test.sh
```

See `tests/synthetic/README.md` for details.

## `tests/real/` — production-style public-data qualification

The realistic qualification uses two exact MD5-pinned public ENA runs and the pre-installed
maintained GDC reference bundle. It deliberately exercises the site's **actual persistent
local/SLURM execution profile** and the same manual copy-edit-run configuration workflow expected
for future consortium datasets.

It is documentation-led rather than setup-script-led: partners copy/configure a profile once,
copy the maintained realistic `config.yaml` and `samples.tsv` into a fresh run directory, edit only
documented paths, invoke Snakemake directly, and run the maintained validator.

Start by downloading the exact public inputs:

```bash
bash tests/real/get_test_data.sh
```

Consortium sites follow `docs/users/consortium/SOPs/SOP-02-site-qualification.md`; `tests/real/README.md` provides the fixture/validator implementation details and general technical companion procedure.

The synthetic test answers "does the workflow implementation behave correctly?"; the realistic
test additionally answers "does this site run the consortium production interface against the
qualified GDC bundle and reproduce the harmonized realistic outputs?"

## Maintainer baseline updates

Ordinary qualification always compares against the maintained expected checksums. If an intentional
change requires a deterministic baseline update, use the review/apply helper in
`tests/maintainers/update_validation_baseline.sh` and follow
`docs/maintainers/qualification-baselines.md`. The helper does not update synthetic input-fixture or
GDC reference integrity manifests.
