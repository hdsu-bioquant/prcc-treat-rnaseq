#!/usr/bin/env python3
"""Pipeline-native UMI metadata and UMI-tools extraction compilation.

The public sample-sheet interface intentionally describes biological/technical
semantics rather than UMI-tools command syntax.  The currently supported class is:

    one fixed-length contiguous UMI at the 5' start of read 1 or read 2,
    optionally followed by a fixed number of additional bases to discard.

Future UMI architectures should extend this semantic layer rather than exposing
UMI-tools regexes in partner-facing metadata.
"""

from __future__ import annotations

from dataclasses import dataclass
import re

UMI_COLUMNS = ("umi_pattern", "umi_location", "umi_discard_bases")
ALLOWED_UMI_LOCATIONS = {"read1_start", "read2_start"}
UMI_PATTERN_RE = re.compile(r"^N+$")


@dataclass(frozen=True)
class UmiSpec:
    """Validated pipeline-native description of one supported UMI architecture."""

    pattern: str
    location: str
    discard_bases: int

    @property
    def length(self) -> int:
        return len(self.pattern)


@dataclass(frozen=True)
class UmiToolsExtractSpec:
    """Internal UMI-tools extraction representation compiled from :class:`UmiSpec`."""

    regex: str
    umi_read: int


def parse_umi_spec(pattern: str, location: str, discard_bases: str | int) -> UmiSpec:
    """Validate semantic UMI metadata and return a normalized :class:`UmiSpec`.

    This intentionally supports only the class claimed by the current pipeline:
    fixed-length ``N+`` UMIs at the 5' start of read 1 or read 2, with an optional
    fixed number of immediately adjacent bases removed after the UMI.
    """

    pattern = str(pattern).strip().upper()
    location = str(location).strip().lower()
    discard_text = str(discard_bases).strip()

    if not UMI_PATTERN_RE.fullmatch(pattern):
        raise ValueError(
            f"umi_pattern='{pattern}' is unsupported; currently use one or more 'N' "
            "characters to describe a fixed-length contiguous UMI (for example NNNNNN)"
        )

    if location not in ALLOWED_UMI_LOCATIONS:
        raise ValueError(
            f"umi_location='{location}' is unsupported; currently supported: "
            + ", ".join(sorted(ALLOWED_UMI_LOCATIONS))
        )

    if not discard_text.isdigit():
        raise ValueError(
            f"umi_discard_bases='{discard_text}' is invalid; use a non-negative integer"
        )

    return UmiSpec(
        pattern=pattern,
        location=location,
        discard_bases=int(discard_text),
    )


def compile_umitools_extract(spec: UmiSpec) -> UmiToolsExtractSpec:
    """Compile one pipeline-native UMI specification for ``umi_tools extract``.

    Regex mode lets UMI-tools distinguish molecular-identifier bases (``umi_1``)
    from additional bases that must simply be removed (``discard_1``).  The regex
    is an internal implementation detail and is never required in samples.tsv.
    """

    regex = rf"^(?P<umi_1>.{{{spec.length}}})"
    if spec.discard_bases:
        regex += rf"(?P<discard_1>.{{{spec.discard_bases}}})"

    umi_read = 1 if spec.location == "read1_start" else 2
    return UmiToolsExtractSpec(regex=regex, umi_read=umi_read)
