# General configuration

Start every run from copied templates:

```bash
mkdir -p /path/to/run
cp templates/config.yaml /path/to/run/config.yaml
cp templates/samples.tsv /path/to/run/samples.tsv
```

Relative paths are resolved from the process working directory, not from the location of the copied YAML/TSV files. Absolute paths are recommended for run-owned data.

## Supported core library combinations

- `full_length` with `paired`
- `quantseq` with `single`

One sample-sheet row represents one sequencing library. `library_id` is the workflow/output key; `sample_id` identifies the biological sample and may repeat across libraries.

Validate the completed sheet before launching:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
```

See `templates/README.md` for the complete field descriptions and validation rules.

## Consortium switch

`consortium_run` is an explicit boolean.

- `true` enables the maintained consortium reference/harmonization checks.
- `false` is appropriate for deliberate non-consortium or custom-reference analyses.

## Optional modules

Fusion, TE, ASE, and RSeQC modules are disabled by default. Their availability does not imply consortium qualification or support.
