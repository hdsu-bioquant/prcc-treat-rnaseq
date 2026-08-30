# Technical debt and deferred cleanup

This file records small infrastructure/maintenance leftovers that are intentionally not changed as part of documentation-only work. It is not a release roadmap and should contain concrete repository issues rather than conversational development history.

## UMI-tools compatibility remnants

`containers/pull_images.sh` still contains compatibility handling for the temporary historical filename `containers/sif/umitools-1.1.6.sif`. Normal maintained operation uses `containers/sif/umitools.sif`. The compatibility branch can be removed in a future behavior-changing cleanup once support for older checkouts is no longer useful.

The synthetic and realistic validators also contain explicit UMI-tools version/source expectations in addition to the maintained `workflow/config/software_versions.yaml` manifest. Consider consolidating these assertions so the manifest remains the single version source while preserving validation coverage.

