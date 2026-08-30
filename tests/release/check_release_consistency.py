#!/usr/bin/env python3
"""Static release-consistency checks for maintained release metadata and core defaults."""

from __future__ import annotations

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
for key in ("pipeline_name", "pipeline_release", "output_contract"):
    if key not in release:
        fail(f"workflow/release.yaml missing {key}")
if release["pipeline_name"] != "pRCC-RNA-Seq":
    fail("workflow/release.yaml pipeline_name must be pRCC-RNA-Seq")
if not isinstance(release["pipeline_release"], str) or not release["pipeline_release"].strip():
    fail("workflow/release.yaml pipeline_release must be a non-empty string")
if int(release["output_contract"]) != 1:
    fail("output contract changed from 1; update policy/docs deliberately before release")


controller_env_path = ROOT / "environments" / "controller.yaml"
controller_env = load_yaml(controller_env_path)
if not isinstance(controller_env, dict):
    fail("environments/controller.yaml must be a YAML mapping")
controller_deps = controller_env.get("dependencies", [])
if not isinstance(controller_deps, list):
    fail("environments/controller.yaml dependencies must be a list")
controller_specs = {dep for dep in controller_deps if isinstance(dep, str)}
for required in ("python=3.13", "snakemake=9.19.0"):
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
    ROOT / "docs" / "maintainers" / "technical-debt.md",
)
for path in required_docs:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"required documentation missing/empty: {path.relative_to(ROOT)}")

# Controlled consortium documentation is self-contained and must not depend on
# the unversioned/convenience general-user guide for normative instructions.
for path in (ROOT / "docs" / "users" / "consortium").rglob("*.md"):
    text = path.read_text()
    if "docs/users/general" in text or "../general" in text:
        fail(f"consortium documentation must not depend on general-user docs: {path.relative_to(ROOT)}")

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
for required_link in (
    "docs/users/general/README.md",
    "docs/users/consortium/README.md",
    "docs/maintainers/README.md",
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
for required in ("prcc-rnaseq-slurm", "prcc-rnaseq-local", 'PROFILE='):
    if required not in sop3:
        fail(f"SOP-03 must document both qualified SLURM/local profile paths: {required}")

for path in (
    ROOT / "docs" / "maintainers" / "README.md",
    ROOT / "tests" / "release" / "README.md",
):
    if "tests/release/check_release_consistency.py" not in path.read_text():
        fail(f"maintainer repository consistency check is not documented in {path.relative_to(ROOT)}")

print(
    "PASS: release metadata, controller/preflight policy, supported-core configs, maintained qualification artifacts, and controlled documentation structure are consistent"
)
