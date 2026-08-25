# Production GDC reference installation

The maintained GDC human reference bundle is installed **before** analysis. Normal Snakemake execution does not download GDC references and does not rebuild the GDC STAR index.

## Two integrity questions

The installation uses two deliberately separate checks:

1. **Official GDC archive identity (MD5).** `gdc_resources.tsv` contains the pinned GDC URLs, archive names and official GDC MD5 values. `get_gdc_references.sh` verifies those MD5s before extraction. This establishes that the downloaded archives are exactly the maintained GDC resources.
2. **Consortium installed-bundle identity (SHA256).** After the realistic production qualification has passed, maintainers will generate one canonical `resources/gdc_installed_reference.sha256` file from that qualified extracted installation. Partners verify against that maintainer-owned manifest; they do not generate their own baseline.

The second layer is intentionally **not frozen yet during development**. Until `gdc_installed_reference.sha256` exists in the repository, consortium-mode runs enforce the exact maintained filenames and fast structural checks but do not require a qualification stamp.

## Install once

From the repository root:

```bash
bash resources/get_gdc_references.sh resources/gdc
```

A different shared/site reference directory may be supplied. The installer retains the proven direct `wget` transfer behavior, verifies every archive against the official GDC MD5 before extraction, retains the archives, and reuses already valid archives/installed targets.

If a maintainer-frozen canonical installed-reference manifest is present in the repository, the installer additionally verifies the extracted bundle against it and writes a small qualification stamp inside the installed reference directory. It does **not** re-hash the extracted bundle at the start of every pipeline run.

## Routine checks

Fast structural verification:

```bash
bash resources/verify_gdc_references.sh /path/to/gdc
```

Recheck retained downloaded archives against the official GDC MD5s:

```bash
bash resources/verify_gdc_references.sh --archives /path/to/gdc
```

After a canonical manifest has been frozen, qualify a newly installed or copied consortium bundle once:

```bash
bash resources/verify_gdc_references.sh --qualify /path/to/gdc
```

This performs the expensive installed-file SHA256 verification and writes `.prcc_treat_reference_qualification.tsv`. Normal `consortium_run: true` workflow startup then checks only the tiny qualification stamp against the canonical manifest identity, plus fast reference-structure checks.

A full `--qualify` check should be repeated after copying/migrating the reference installation, after suspected storage corruption/modification, or when the maintained reference bundle changes.

## Maintainer-only canonical manifest generation

Do **not** generate a canonical installed manifest merely because the files have downloaded successfully. The intended sequence is:

1. install the exact GDC archives and verify official MD5s;
2. pass the realistic production qualification against that installation;
3. generate the canonical installed-reference manifest once:

```bash
bash resources/verify_gdc_references.sh --write-canonical-manifest /path/to/qualified/gdc
```

4. review/freeze `resources/gdc_installed_reference.sha256` in the maintained repository;
5. qualify the site installation against the now-frozen manifest:

```bash
bash resources/verify_gdc_references.sh --qualify /path/to/qualified/gdc
```

Partners should receive the frozen canonical manifest with the pipeline and verify against it. They should not replace it with a locally generated manifest.
