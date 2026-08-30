#!/usr/bin/env python3
"""Static release-consistency checks for maintained release metadata and core defaults."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_yaml(path: Path):
    with path.open() as fh:
        return yaml.safe_load(fh)


release_path = ROOT / "workflow" / "release.yaml"
release = load_yaml(release_path)
if not isinstance(release, dict):
    fail("workflow/release.yaml must be a YAML mapping")
for key in ("pipeline_name", "repository_slug", "pipeline_release", "output_contract"):
    if key not in release:
        fail(f"workflow/release.yaml missing {key}")
pipeline_name = release["pipeline_name"]
if not isinstance(pipeline_name, str) or not pipeline_name.strip():
    fail("workflow/release.yaml pipeline_name must be a non-empty string")
if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", pipeline_name):
    fail("workflow/release.yaml pipeline_name must be a portable display identifier")
repository_slug = release["repository_slug"]
if not isinstance(repository_slug, str) or not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", repository_slug):
    fail("workflow/release.yaml repository_slug must be a lowercase hyphenated repository identifier")
pipeline_release = release["pipeline_release"]
if not isinstance(pipeline_release, str) or not pipeline_release.strip():
    fail("workflow/release.yaml pipeline_release must be a non-empty string")
if pipeline_release != "development" and not re.fullmatch(r"\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?", pipeline_release):
    fail("workflow/release.yaml pipeline_release must be development or a SemVer-like release")
if int(release["output_contract"]) != 1:
    fail("output contract changed from 1; update policy/docs deliberately before release")


controller_env_path = ROOT / "environments" / "controller.yaml"
controller_env = load_yaml(controller_env_path)
if not isinstance(controller_env, dict):
    fail("environments/controller.yaml must be a YAML mapping")
if controller_env.get("name") != f"{repository_slug}-controller":
    fail("controller environment name must agree with workflow/release.yaml repository_slug")
controller_deps = controller_env.get("dependencies", [])
if not isinstance(controller_deps, list):
    fail("environments/controller.yaml dependencies must be a list")
controller_specs = {dep for dep in controller_deps if isinstance(dep, str)}
for required in ("python=3.13", "snakemake=9.20.0"):
    if required not in controller_specs:
        fail(f"maintained controller environment missing required constraint: {required}")
if not any(dep.startswith("snakemake-executor-plugin-slurm") for dep in controller_specs):
    fail("maintained controller environment must include the SLURM executor plugin")

preflight = ROOT / "scripts" / "verify_installation.py"
if not preflight.is_file():
    fail("read-only installation preflight helper is missing")

software_manifest = load_yaml(ROOT / "workflow" / "config" / "software_versions.yaml")
tools = software_manifest.get("tools", {}) if isinstance(software_manifest, dict) else {}
if not tools:
    fail("software/container manifest has no tools")
for name, spec in tools.items():
    if not isinstance(spec, dict):
        fail(f"software manifest entry {name} must be a mapping")
    if spec.get("pull_default", False) and not str(spec.get("version_probe", "")).strip():
        fail(f"default container {name} must declare a version_probe for preflight verification")

configs = {
    "template": load_yaml(ROOT / "templates" / "config.yaml"),
    "realistic": load_yaml(ROOT / "tests" / "real" / "config.yaml"),
}
for label, cfg in configs.items():
    if "pipeline_release" in cfg:
        fail(f"{label} run config must not contain maintainer-owned pipeline_release")
    if cfg.get("consortium_run") is not True:
        fail(f"{label} config must retain consortium_run: true")
    if cfg.get("full_length", {}).get("trim_adapters") is not False:
        fail(f"{label} config full_length.trim_adapters drifted")
    if cfg.get("quantseq", {}).get("bbduk_polyA") is not True:
        fail(f"{label} config quantseq.bbduk_polyA drifted")
    if any(bool(v) for v in cfg.get("modules", {}).values()):
        fail(f"{label} config enables an optional unqualified module")

ref_fields = ("mode", "genome_fasta", "gtf", "star_index", "sjdb_overhang")
for field in ref_fields:
    if configs["template"]["reference"].get(field) != configs["realistic"]["reference"].get(field):
        fail(f"template/realistic reference.{field} disagree")
if configs["template"]["star"].get("gdc_params") != configs["realistic"]["star"].get("gdc_params"):
    fail("template/realistic star.gdc_params disagree")

frozen = (
    ROOT / "resources" / "gdc_installed_reference.sha256",
    ROOT / "tests" / "synthetic" / "expected" / "validation_checksums.sha256",
    ROOT / "tests" / "real" / "expected" / "validation_checksums.sha256",
)
for path in frozen:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"required frozen qualification artifact missing/empty: {path.relative_to(ROOT)}")


synthetic_fixture_manifest = ROOT / "tests" / "synthetic" / "checksums.sha256"
if not synthetic_fixture_manifest.is_file() or synthetic_fixture_manifest.stat().st_size == 0:
    fail("synthetic fixture integrity manifest is missing or empty")
fixture_manifest_text = synthetic_fixture_manifest.read_text()
if "tests/synthetic/expected/validation_checksums.sha256" in fixture_manifest_text:
    fail(
        "synthetic fixture checksums must not include the separately maintained validation-output baseline"
    )

baseline_helper = ROOT / "tests" / "maintainers" / "update_validation_baseline.sh"
if not baseline_helper.is_file():
    fail("maintainer validation-baseline update helper is missing")

# Root release-facing metadata must exist and agree on stable project identity/license.
license_path = ROOT / "LICENSE"
citation_path = ROOT / "CITATION.cff"
changelog_path = ROOT / "CHANGELOG.md"
for path in (license_path, citation_path, changelog_path):
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"required release metadata missing/empty: {path.relative_to(ROOT)}")

license_text = license_path.read_text()
if not license_text.startswith("MIT License\n"):
    fail("LICENSE must contain the maintained MIT license text")
if not re.search(r"Copyright \(c\) [^\n]*Health Data Science Unit", license_text):
    fail("LICENSE must identify Health Data Science Unit as the maintained copyright holder")

citation = load_yaml(citation_path)
if not isinstance(citation, dict):
    fail("CITATION.cff must be a YAML mapping")
if str(citation.get("cff-version", "")) != "1.2.0":
    fail("CITATION.cff must use CFF 1.2.0")
if citation.get("title") != pipeline_name:
    fail("CITATION.cff title must agree with workflow/release.yaml pipeline_name")
if citation.get("license") != "MIT":
    fail("CITATION.cff license must agree with the repository MIT license")
if citation.get("repository-code") != f"https://github.com/hdsu-bioquant/{repository_slug}":
    fail("CITATION.cff repository-code must agree with workflow/release.yaml repository_slug")
if pipeline_release != "development":
    if str(citation.get("version", "")) != pipeline_release:
        fail("CITATION.cff version must agree with workflow/release.yaml pipeline_release")
    if not citation.get("date-released"):
        fail("CITATION.cff must contain date-released for a release candidate/release")
authors = citation.get("authors", [])
if not isinstance(authors, list) or not authors:
    fail("CITATION.cff must contain software authors")
for family in ("Bökenkamp", "Schwarz", "Herrmann"):
    if not any(isinstance(author, dict) and author.get("family-names") == family for author in authors):
        fail(f"CITATION.cff missing maintained software author: {family}")

changelog_text = changelog_path.read_text()
if not changelog_text.startswith("# Changelog\n"):
    fail("CHANGELOG.md must start with the Changelog heading")
if f"All notable changes to {pipeline_name}" not in changelog_text:
    fail("CHANGELOG.md project name must agree with workflow/release.yaml pipeline_name")
if "## [Unreleased]" not in changelog_text:
    fail("CHANGELOG.md must retain an Unreleased section")
if "### Notes" in changelog_text:
    fail("CHANGELOG.md should use release-change categories rather than a free-form Notes section")
if pipeline_release != "development":
    release_date = str(citation.get("date-released"))
    if f"## [{pipeline_release}] - {release_date}" not in changelog_text:
        fail("CHANGELOG.md release heading/date must agree with workflow/release.yaml and CITATION.cff")

required_docs = (
    ROOT / "docs" / "users" / "README.md",
    ROOT / "docs" / "users" / "general" / "README.md",
    ROOT / "docs" / "users" / "general" / "installation.md",
    ROOT / "docs" / "users" / "general" / "configuration.md",
    ROOT / "docs" / "users" / "general" / "running-the-pipeline.md",
    ROOT / "docs" / "users" / "general" / "outputs.md",
    ROOT / "docs" / "users" / "general" / "troubleshooting.md",
    ROOT / "docs" / "users" / "consortium" / "README.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-01-site-installation.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-02-site-qualification.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-03-consortium-run.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-04-results-delivery.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-05-upgrade-requalification.md",
    ROOT / "docs" / "users" / "consortium" / "io-definitions" / "sample-sheet.md",
    ROOT / "docs" / "users" / "consortium" / "io-definitions" / "run-configuration.md",
    ROOT / "docs" / "users" / "consortium" / "io-definitions" / "results-contract.md",
    ROOT / "docs" / "users" / "consortium" / "io-definitions" / "qc-and-run-metadata.md",
    ROOT / "docs" / "maintainers" / "README.md",
    ROOT / "docs" / "maintainers" / "controlled-documentation.md",
    ROOT / "docs" / "maintainers" / "technical-debt.md",
)
for path in required_docs:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"required documentation missing/empty: {path.relative_to(ROOT)}")

def parse_controlled_table(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != 2 or cells[0] in {"Field", "---"}:
            continue
        values[cells[0]] = cells[1].strip("`")
    return values


controlled_docs: list[tuple[Path, str, str]] = []
for path in sorted((ROOT / "docs" / "users" / "consortium" / "SOPs").glob("SOP-*.md")):
    match = re.match(r"(SOP-\d+)-", path.name)
    if not match:
        fail(f"controlled SOP filename lacks stable SOP ID: {path.relative_to(ROOT)}")
    controlled_docs.append((path, "SOP ID", match.group(1)))

io_ids = {
    "sample-sheet.md": "IO-01",
    "run-configuration.md": "IO-02",
    "results-contract.md": "IO-03",
    "qc-and-run-metadata.md": "IO-04",
}
for name, doc_id in io_ids.items():
    controlled_docs.append((ROOT / "docs" / "users" / "consortium" / "io-definitions" / name, "Document ID", doc_id))

seen_ids: set[str] = set()
for path, id_field, expected_id in controlled_docs:
    metadata = parse_controlled_table(path)
    for field in (id_field, "Status", "Document version", "Owner", "Applicable pipeline release", "Last revised"):
        if not metadata.get(field):
            fail(f"controlled document missing metadata field {field}: {path.relative_to(ROOT)}")
    if metadata[id_field] != expected_id:
        fail(f"controlled document ID disagrees with maintained identity: {path.relative_to(ROOT)}")
    if expected_id in seen_ids:
        fail(f"duplicate controlled document ID: {expected_id}")
    seen_ids.add(expected_id)
    if metadata["Status"] not in {"Draft", "Pilot", "Approved", "Superseded"}:
        fail(f"invalid controlled document status in {path.relative_to(ROOT)}: {metadata['Status']}")
    if not re.fullmatch(r"\d+\.\d+", metadata["Document version"]):
        fail(f"invalid controlled document version in {path.relative_to(ROOT)}: {metadata['Document version']}")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", metadata["Last revised"]):
        fail(f"invalid Last revised date in {path.relative_to(ROOT)}: {metadata['Last revised']}")
    if pipeline_release != "development":
        if metadata["Applicable pipeline release"] != pipeline_release:
            fail(f"controlled document applicability must match {pipeline_release}: {path.relative_to(ROOT)}")
        if metadata["Status"] not in {"Pilot", "Approved"}:
            fail(f"release-facing controlled document must be Pilot or Approved: {path.relative_to(ROOT)}")

# Controlled consortium documentation is self-contained and must not depend on
# the unversioned/convenience general-user guide for normative instructions.
for path in (ROOT / "docs" / "users" / "consortium").rglob("*.md"):
    text = path.read_text()
    if "docs/users/general" in text or "../general" in text:
        fail(f"consortium documentation must not depend on general-user docs: {path.relative_to(ROOT)}")

# The pre-RC rename/controller upgrade should be complete outside explicit changelog history.
legacy_markers = ("pRCC-RNA-Seq", "prcc-rnaseq", "prcc_rnaseq", "snakemake=9.19.0", "Snakemake 9.19.0")
for path in ROOT.rglob("*"):
    if not path.is_file() or path == changelog_path or path == Path(__file__).resolve():
        continue
    if path.suffix not in {".md", ".py", ".sh", ".yaml", ".yml", ".cff"} and path.name not in {".gitignore", "README.md"}:
        continue
    try:
        text = path.read_text()
    except UnicodeDecodeError:
        continue
    for marker in legacy_markers:
        if marker in text:
            fail(f"legacy pipeline/controller identity remains in {path.relative_to(ROOT)}: {marker}")

# Keep development archaeology out of normal user-facing documentation.
user_docs = [ROOT / "README.md", *list((ROOT / "docs" / "users").rglob("*.md"))]
stale_user_phrases = (
    "umitools-1.1.6.sif",
    "UMI-tools 1.1.4",
    "expected during development",
    "not frozen yet",
    "will be frozen",
)
for path in user_docs:
    text = path.read_text()
    for phrase in stale_user_phrases:
        if phrase in text:
            fail(f"stale development wording remains in user documentation {path.relative_to(ROOT)}: {phrase}")

root_readme = (ROOT / "README.md").read_text()
if pipeline_name not in root_readme:
    fail("top-level README must contain the maintained pipeline name")
for required_link in (
    "docs/users/general/README.md",
    "docs/users/consortium/README.md",
    "docs/maintainers/README.md",
    "CITATION.cff",
    "LICENSE",
    "CHANGELOG.md",
):
    if required_link not in root_readme:
        fail(f"top-level README missing documentation entry point: {required_link}")

for path in (
    ROOT / "README.md",
    ROOT / "templates" / "profiles" / "README.md",
    ROOT / "docs" / "users" / "general" / "installation.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-02-site-qualification.md",
    ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-03-consortium-run.md",
):
    if "scripts/verify_installation.py" not in path.read_text():
        fail(f"installation preflight is not documented in {path.relative_to(ROOT)}")

sop1 = (ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-01-site-installation.md").read_text()
for required in ("workflow/release.yaml", "containers/pull_images.sh", "tests/synthetic/run_test.sh"):
    if required not in sop1:
        fail(f"SOP-01 missing first-line installation step: {required}")
if "scripts/verify_installation.py" in sop1:
    fail("SOP-01 must end at lightweight synthetic qualification; run-specific preflight belongs in SOP-02")
if "resources/get_gdc_references.sh" in sop1:
    fail("SOP-01 must not install the heavyweight GDC reference bundle; that belongs in SOP-02")

sop2 = (ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-02-site-qualification.md").read_text()
for required in (
    "resources/get_gdc_references.sh",
    "resources/verify_gdc_references.sh --qualify",
    "tests/real/get_test_data.sh",
    "templates/profiles/slurm",
    "templates/profiles/local",
    "scripts/verify_installation.py",
    "tests/real/validate_results.py",
):
    if required not in sop2:
        fail(f"SOP-02 missing consortium site-qualification step: {required}")

sop3 = (ROOT / "docs" / "users" / "consortium" / "SOPs" / "SOP-03-consortium-run.md").read_text()
for required in (f"{repository_slug}-slurm", f"{repository_slug}-local", 'PROFILE='):
    if required not in sop3:
        fail(f"SOP-03 must document both qualified SLURM/local profile paths: {required}")

for path in (
    ROOT / "docs" / "maintainers" / "README.md",
    ROOT / "tests" / "release" / "README.md",
):
    if "tests/release/check_release_consistency.py" not in path.read_text():
        fail(f"maintainer repository consistency check is not documented in {path.relative_to(ROOT)}")

maintainer_readme = (ROOT / "docs" / "maintainers" / "README.md").read_text()
if "controlled-documentation.md" not in maintainer_readme:
    fail("maintainer README must link the controlled-documentation policy")
release_policy_text = (ROOT / "docs" / "maintainers" / "release-policy.md").read_text()
if "## Changelog maintenance" not in release_policy_text:
    fail("release policy must define changelog maintenance")

print(
    "PASS: release/citation/license metadata, controller/preflight policy, supported-core configs, maintained qualification artifacts, and controlled documentation are consistent"
)
