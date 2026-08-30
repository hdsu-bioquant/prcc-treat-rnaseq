# User documentation

pRCC-TREAT-RNA-Seq has two user-documentation tracks because independent analyses and harmonized consortium analyses have different operational requirements.

## General users

[`general/`](general/README.md) is a conventional pipeline guide for users applying pRCC-TREAT-RNA-Seq outside a controlled consortium workflow. It explains supported inputs, installation, configuration, execution, outputs, and troubleshooting. General users may deliberately adapt supported resources and settings and normally set `consortium_run: false` when they do so.

## Consortium users

[`consortium/`](consortium/README.md) is the controlled documentation set for multi-site harmonized operation. It is intentionally self-contained and does not rely on general-user prose for normative instructions.

Consortium documentation consists of:

- numbered, versioned [`SOPs/`](consortium/SOPs/) describing required procedures;
- [`io-definitions/`](consortium/io-definitions/) defining the maintained input, configuration, output, QC, and run-metadata contracts.
