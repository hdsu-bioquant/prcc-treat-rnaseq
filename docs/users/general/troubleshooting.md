# General troubleshooting

## Sample-sheet validation fails

Run the validator directly to obtain the focused error message:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
```

Common causes include duplicate/unsafe identifiers, inconsistent assay/layout combinations, invalid FASTQ paths, and incomplete UMI metadata.

## References fail validation

For GDC resources, first run the fast structural check:

```bash
bash resources/verify_gdc_references.sh /path/to/gdc
```

For a complete installed-reference integrity check:

```bash
python scripts/verify_installation.py --configfile /path/to/run/config.yaml --full-reference-check
```

## Container/version check fails

Refresh maintained default images:

```bash
bash containers/pull_images.sh
```

Then rerun `scripts/verify_installation.py`. Existing images are checked using their declared tool-version probes.

## Profile or scheduler failures

Use `templates/profiles/README.md` to separate site scheduler/account/partition/bind-path settings from biological run configuration.

## Qualification failures

Synthetic and realistic qualification procedures are documented under `tests/`. Maintainers should not replace expected validation baselines until a mismatch has been understood; see `docs/maintainers/qualification-baselines.md`.
