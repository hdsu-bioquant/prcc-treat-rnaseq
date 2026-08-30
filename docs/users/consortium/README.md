# Consortium user documentation

This documentation defines the controlled operating interface for harmonized multi-site pRCC-RNA-Seq analyses. It is written for pRCC-TREAT but is intentionally usable by other consortia that adopt the same maintained processing contract.

The consortium documentation is self-contained. Normative instructions do not depend on the general-user guide.

## Required onboarding sequence

Consortium sites should proceed in this order:

1. [`SOP-01 — Site software installation`](SOPs/SOP-01-site-installation.md): obtain the designated release, establish the controller/runtime/container software layer, and pass the synthetic qualification.
2. [`SOP-02 — Site qualification`](SOPs/SOP-02-site-qualification.md): install/qualify GDC resources, establish the persistent local or SLURM profile, rehearse run configuration, run preflight, and pass realistic public-data qualification.
3. [`SOP-03 — Consortium run`](SOPs/SOP-03-consortium-run.md): process restricted consortium libraries using the qualified installation/profile.
4. [`SOP-04 — Results delivery`](SOPs/SOP-04-results-delivery.md): prepare the portable results package for transfer.
5. [`SOP-05 — Upgrade and requalification`](SOPs/SOP-05-upgrade-requalification.md): adopt later releases or material site-infrastructure changes.

The separation between SOP-01 and SOP-02 is deliberate: sites first prove the lightweight software stack with the synthetic fixture before investing in the larger GDC reference installation and production-style realistic qualification.

## IO definitions

- [`sample-sheet.md`](io-definitions/sample-sheet.md)
- [`run-configuration.md`](io-definitions/run-configuration.md)
- [`results-contract.md`](io-definitions/results-contract.md)
- [`qc-and-run-metadata.md`](io-definitions/qc-and-run-metadata.md)

The SOPs define required procedures. The IO definitions define the maintained input, configuration, output, QC, and run-metadata contracts used by those procedures.
