# Controlled consortium documentation

This policy defines lightweight maintenance rules for the consortium-facing documents under
`docs/users/consortium/`. It applies to both numbered SOPs and IO definitions. The aim is to make
site instructions and interface contracts traceable to a pipeline release without introducing a
separate document-management system outside the repository.

## Controlled document sets

The controlled consortium documentation consists of:

- `docs/users/consortium/SOPs/` — required operating procedures;
- `docs/users/consortium/io-definitions/` — maintained input, configuration, output, QC, and
  run-metadata contracts.

The general-user documentation under `docs/users/general/` is intentionally not controlled by this
policy. Consortium documents remain self-contained and must not depend on general-user prose for
normative instructions.

## Required metadata

Each controlled document records, near the top of the file:

- a stable document identifier (`SOP-..` or `IO-..`);
- status;
- document version;
- owner;
- applicable pipeline release; and
- last-revised date.

The document version is independent of the pipeline version. The pipeline release identifies the
software implementation; the document version identifies the controlled text/contract.

## Document status

Use the following status values:

- `Draft` — under active maintainer development and not yet issued for partner use;
- `Pilot` — issued for a release-candidate or controlled partner pilot;
- `Approved` — maintained procedure/contract for routine consortium use;
- `Superseded` — no longer current for new work.

Historical copies do not need to be duplicated in the live documentation tree. Git tags/releases
retain the exact document versions issued with earlier pipeline releases.

## Document versioning

Use lightweight `MAJOR.MINOR` document versions.

- Increment **MINOR** for compatible clarifications, added operational detail, wording corrections,
  or other changes that do not materially alter the required procedure or IO contract.
- Increment **MAJOR** when a document changes required site actions, responsibilities, data-sharing
  boundaries, input semantics, or an IO contract in a materially incompatible way.
- Pre-release documents may remain in the `0.x` series. Promotion of status or update of the
  applicable pipeline-release field does not by itself require a document-version increment when
  the normative content is unchanged.

Update `Last revised` whenever the controlled file changes.

## Release applicability

`Applicable pipeline release` should name the exact release or release family for which the
controlled document is issued. During repository development it may be `development`. Before a
partner pilot or production release, maintainers update the field and status to match the intended
release.

A controlled document change does not automatically imply a pipeline-version or qualification
change. Validation level is determined by the underlying repository change as defined in
`release-policy.md`. Pure documentation/metadata edits are normally Level C; documentation that
accompanies a behavioral or output change follows the level of that underlying change.

## Maintainer review

Before issuing a release or release candidate:

1. confirm that each consortium SOP and IO definition contains complete controlled metadata;
2. confirm that status and applicable release are appropriate for the intended distribution;
3. check that normative values agree with maintained templates, manifests, and workflow behavior;
4. run `python tests/release/check_release_consistency.py`; and
5. commit document changes before tagging the release so the tag retains the exact issued text.

The repository consistency check validates document presence and metadata structure, but maintainer
review remains responsible for scientific/operational correctness of the prose.
