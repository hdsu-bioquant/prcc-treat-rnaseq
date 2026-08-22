# Pipeline validation suites

pRCC-RNA-Seq uses two complementary validation levels.

## 1. `synthetic/` — installation and regression smoke test

The synthetic suite is tiny, committed to Git, uses its own miniature reference, and covers
every supported core route (FL/QS, paired/single, UMI/no-UMI). Run it with:

```bash
bash tests/synthetic/run_test.sh
```

It validates exact known counts **and the same `results/` / `restricted/` /
`intermediate/` output contract used in normal analyses**. The deterministic canonical output
subset is additionally compared byte-for-byte against the frozen reference hashes in
`synthetic/expected/validation_checksums.sha256`. See
[`synthetic/README.md`](synthetic/README.md) for details.

For consortium qualification, the synthetic test answers:

> Can this installation execute every core pRCC-RNA-Seq pathway and produce the canonical
> portable result package correctly?

## 2. `real/` — realistic GDC-reference harmonization test

The real-data suite is intentionally separate. Its public FASTQs are downloaded rather
than stored in Git and are processed using the fixed production GDC reference resources and
the site's normal execution profile.

Unlike the synthetic smoke test, this should ultimately be launched in the **same way as a
real production run**, not through a special test wrapper. It answers:

> Can this site's actual HPC/workstation, profile, containers, storage, and fixed GDC
> resources reproduce the consortium-standard outputs on realistic human RNA-seq data?

The real fixture/download procedure and cross-site checksum baseline will be finalized after
the production input/output templates are frozen.
