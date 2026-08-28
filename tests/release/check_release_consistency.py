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
if preflight.stat().st_mode & 0o111 == 0:
    fail("installation preflight helper is not executable")

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
if baseline_helper.stat().st_mode & 0o111 == 0:
    fail("maintainer validation-baseline update helper is not executable")

stale_phrases = {
    ROOT / "README.md": (
        "do not generate/freeze that canonical SHA256 manifest yet",
        "will be frozen only after the realistic production qualification succeeds",
    ),
    ROOT / "resources" / "README.md": (
        "intentionally **not frozen yet during development**",
    ),
    ROOT / "tests" / "real" / "expected" / "README.md": (
        "validation_checksums.sha256` is intentionally absent",
    ),
}
for path, phrases in stale_phrases.items():
    text = path.read_text()
    for phrase in phrases:
        if phrase in text:
            fail(f"stale qualification wording remains in {path.relative_to(ROOT)}: {phrase}")

for path in (ROOT / "README.md", ROOT / "templates" / "profiles" / "README.md"):
    if "scripts/verify_installation.py" not in path.read_text():
        fail(f"installation preflight is not documented in {path.relative_to(ROOT)}")

print(
    "PASS: release metadata, controller/preflight policy, supported-core configs, maintained qualification artifacts, and key documentation are consistent"
)
