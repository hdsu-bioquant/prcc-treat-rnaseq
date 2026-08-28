#!/usr/bin/env python3
"""Read-only installation preflight for pRCC-RNA-Seq.

The preflight validates a controller environment, a site GDC reference
installation, required local containers and their expected tool versions, and
(optionally) an execution profile. It does not install, pull, qualify, repair, or rewrite any of them.
An explicit --report path writes a copy of the report.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import yaml
from packaging.specifiers import SpecifierSet
from packaging.version import Version

ROOT = Path(__file__).resolve().parents[1]
ENV_SPEC = ROOT / "environments" / "controller.yaml"
RELEASE_SPEC = ROOT / "workflow" / "release.yaml"
SOFTWARE_SPEC = ROOT / "workflow" / "config" / "software_versions.yaml"
REF_VERIFY = ROOT / "resources" / "verify_gdc_references.sh"
REF_MANIFEST = ROOT / "resources" / "gdc_installed_reference.sha256"
REF_STAMP = ".prcc_treat_reference_qualification.tsv"


@dataclass
class Check:
    status: str
    component: str
    observed: str
    expected: str = ""
    detail: str = ""


def load_yaml(path: Path):
    with path.open() as handle:
        return yaml.safe_load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def conda_dependency_specs() -> dict[str, str]:
    env = load_yaml(ENV_SPEC)
    deps = env.get("dependencies", []) if isinstance(env, dict) else []
    result: dict[str, str] = {}
    for raw in deps:
        if not isinstance(raw, str):
            continue
        for index, char in enumerate(raw):
            if char in "=<>!~":
                result[raw[:index]] = raw[index:]
                break
        else:
            result[raw] = ""
    return result


def version_matches(installed: str, spec: str) -> bool:
    if not spec:
        return True
    if spec.startswith("=") and not spec.startswith("=="):
        wanted = spec[1:]
        # Conda's single '=' commonly fixes a version or a version prefix.
        return installed == wanted or installed.startswith(wanted + ".")
    return Version(installed) in SpecifierSet(spec)


def package_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def run_capture(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def read_stamp(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t", 1)
        if len(fields) == 2:
            values[fields[0]] = fields[1]
    return values


def find_profile_file(raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw).expanduser().resolve()
    if path.is_file():
        return path
    if path.is_dir():
        for name in ("profile.yaml", "config.yaml"):
            candidate = path / name
            if candidate.is_file():
                return candidate
    return path


def check_controller(checks: list[Check]) -> None:
    specs = conda_dependency_specs()
    py_spec = specs.get("python", "")
    py_version = platform.python_version()
    checks.append(Check("PASS" if version_matches(py_version, py_spec) else "FAIL",
                        "controller.python", py_version, py_spec))

    for dist, component in (
        ("snakemake", "controller.snakemake"),
        ("snakemake-executor-plugin-slurm", "controller.slurm_plugin"),
    ):
        installed = package_version(dist)
        expected = specs.get(dist, "")
        if installed is None:
            checks.append(Check("FAIL", component, "not installed", expected))
        else:
            checks.append(Check("PASS" if version_matches(installed, expected) else "FAIL",
                                component, installed, expected))

    runtime = shutil.which("apptainer") or shutil.which("singularity")
    if runtime is None:
        checks.append(Check("FAIL", "container.runtime", "not found", "apptainer or singularity on PATH"))
    else:
        result = run_capture([runtime, "--version"])
        observed = result.stdout.strip().replace("\t", " ") if result.returncode == 0 else "version command failed"
        checks.append(Check("PASS" if result.returncode == 0 else "FAIL",
                            "container.runtime", observed, Path(runtime).name))


def reference_from_config(configfile: Path) -> tuple[Path | None, bool, str]:
    cfg = load_yaml(configfile)
    if not isinstance(cfg, dict):
        return None, True, "invalid config mapping"
    ref = cfg.get("reference", {})
    if not isinstance(ref, dict):
        return None, bool(cfg.get("consortium_run", True)), "missing reference mapping"
    mode = str(ref.get("mode", "gdc"))
    refdir = ref.get("dir")
    if not refdir:
        return None, bool(cfg.get("consortium_run", True)), f"reference.dir missing (mode={mode})"
    return Path(str(refdir)).expanduser().resolve(), bool(cfg.get("consortium_run", True)), mode


def check_reference(
    checks: list[Check],
    refdir: Path | None,
    require_stamp: bool,
    full_reference_check: bool,
) -> None:
    if refdir is None:
        checks.append(Check("WARN", "reference", "not checked", "supply --reference-dir or --configfile"))
        return

    if not refdir.is_dir():
        checks.append(Check("FAIL", "reference.directory", str(refdir), "existing GDC reference directory"))
        return

    command = ["bash", str(REF_VERIFY)]
    if full_reference_check:
        command.append("--canonical")
    command.append(str(refdir))
    result = run_capture(command)
    checks.append(Check("PASS" if result.returncode == 0 else "FAIL",
                        "reference.structure" if not full_reference_check else "reference.canonical_sha256",
                        str(refdir),
                        "fast structure" if not full_reference_check else "canonical manifest",
                        result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""))

    if not require_stamp:
        checks.append(Check("SKIP", "reference.qualification_stamp", "not required by config"))
        return

    if not REF_MANIFEST.is_file():
        checks.append(Check("FAIL", "reference.canonical_manifest", "missing", str(REF_MANIFEST.relative_to(ROOT))))
        return
    manifest_hash = sha256_file(REF_MANIFEST)
    checks.append(Check("PASS", "reference.canonical_manifest", manifest_hash,
                        str(REF_MANIFEST.relative_to(ROOT))))

    stamp_path = refdir / REF_STAMP
    if not stamp_path.is_file():
        checks.append(Check("FAIL", "reference.qualification_stamp", "missing", REF_STAMP,
                            "run resources/verify_gdc_references.sh --qualify once for this installation"))
        return
    stamp = read_stamp(stamp_path)
    observed = stamp.get("canonical_manifest_sha256", "missing")
    checks.append(Check("PASS" if observed == manifest_hash else "FAIL",
                        "reference.qualification_stamp", observed, manifest_hash))


def check_profile(checks: list[Check], raw_profile: str | None) -> None:
    profile = find_profile_file(raw_profile)
    if profile is None:
        checks.append(Check("WARN", "execution.profile", "not checked", "supply --profile for a production preflight"))
        return
    if not profile.is_file():
        checks.append(Check("FAIL", "execution.profile", str(profile), "readable profile YAML"))
        return
    try:
        cfg = load_yaml(profile)
    except Exception as exc:  # pragma: no cover - diagnostic path
        checks.append(Check("FAIL", "execution.profile", str(profile), "valid YAML", str(exc)))
        return
    if not isinstance(cfg, dict):
        checks.append(Check("FAIL", "execution.profile", str(profile), "YAML mapping"))
        return
    executor = str(cfg.get("executor", "local"))
    checks.append(Check("PASS", "execution.profile", str(profile), "readable YAML", f"executor={executor}"))
    if executor == "slurm":
        for command in ("sbatch", "squeue"):
            resolved = shutil.which(command)
            checks.append(Check("PASS" if resolved else "FAIL", f"slurm.{command}", resolved or "not found", "on PATH"))


def tool_rows(all_containers: bool) -> Iterable[tuple[str, dict]]:
    manifest = load_yaml(SOFTWARE_SPEC)
    tools = manifest.get("tools", {}) if isinstance(manifest, dict) else {}
    for name, spec in tools.items():
        if not isinstance(spec, dict):
            continue
        if all_containers or bool(spec.get("pull_default", False)):
            yield str(name), spec


def check_containers(checks: list[Check], all_containers: bool) -> None:
    runtime = shutil.which("apptainer") or shutil.which("singularity")
    if runtime is None:
        checks.append(Check("FAIL", "containers", "not checked", "container runtime unavailable"))
        return

    sif_dir = ROOT / "containers" / "sif"
    for name, spec in tool_rows(all_containers):
        expected_version = str(spec.get("version", ""))
        source = str(spec.get("uri", ""))
        local_name = str(spec.get("local_sif", f"{name}.sif"))
        image = sif_dir / local_name
        label = f"container.{name}"
        if not image.is_file() or image.stat().st_size == 0:
            checks.append(Check("FAIL", label, "missing", str(image.relative_to(ROOT)), source))
            continue

        details = [f"source={source}", f"image={image.relative_to(ROOT)}"]
        probe = str(spec.get("version_probe", "")).strip()
        if probe:
            result = run_capture([runtime, "exec", str(image), "bash", "-lc", probe])
            observed = result.stdout.strip().replace("\n", " | ")
            ok = result.returncode == 0 and expected_version in observed
            checks.append(Check("PASS" if ok else "FAIL", label, observed or "probe produced no output",
                                expected_version, "; ".join(details + [f"probe={probe}"])))
        else:
            checks.append(Check("WARN", label, str(image.relative_to(ROOT)), expected_version,
                                "; ".join(details + ["no version_probe declared"])))


def source_revision() -> str:
    git = shutil.which("git")
    if not git or not (ROOT / ".git").exists():
        return "unavailable"
    result = run_capture([git, "-C", str(ROOT), "rev-parse", "HEAD"])
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def render_report(checks: list[Check]) -> str:
    release = load_yaml(RELEASE_SPEC)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "pRCC-RNA-Seq installation preflight",
        f"timestamp_utc\t{timestamp}",
        f"pipeline_release\t{release.get('pipeline_release', 'unknown')}",
        f"output_contract\t{release.get('output_contract', 'unknown')}",
        f"source_revision\t{source_revision()}",
        "",
        "status\tcomponent\tobserved\texpected\tdetail",
    ]
    for check in checks:
        clean = [
            check.status,
            check.component,
            check.observed.replace("\t", " ").replace("\n", " | "),
            check.expected.replace("\t", " ").replace("\n", " | "),
            check.detail.replace("\t", " ").replace("\n", " | "),
        ]
        lines.append("\t".join(clean))
    failed = sum(c.status == "FAIL" for c in checks)
    warned = sum(c.status == "WARN" for c in checks)
    lines.extend(["", f"summary\t{'FAIL' if failed else 'PASS'}\tfailures={failed}\twarnings={warned}"])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--reference-dir", help="installed GDC reference directory")
    source.add_argument("--configfile", help="run config from which to read reference.dir and consortium_run")
    parser.add_argument("--profile", help="copied Snakemake profile directory or YAML file")
    parser.add_argument("--full-reference-check", action="store_true",
                        help="perform the expensive full canonical reference SHA256 verification")
    parser.add_argument("--all-containers", action="store_true",
                        help="also verify containers for optional modules")
    parser.add_argument("--report", help="write the same preflight report to this path")
    args = parser.parse_args()

    refdir: Path | None = None
    require_stamp = True
    if args.configfile:
        configfile = Path(args.configfile).expanduser().resolve()
        if not configfile.is_file():
            parser.error(f"configfile does not exist: {configfile}")
        refdir, require_stamp, mode = reference_from_config(configfile)
        if mode != "gdc":
            refdir = None
            require_stamp = False
    elif args.reference_dir:
        refdir = Path(args.reference_dir).expanduser().resolve()
    else:
        default = ROOT / "resources" / "gdc"
        refdir = default if default.is_dir() else None

    checks: list[Check] = []
    check_controller(checks)
    check_reference(checks, refdir, require_stamp, args.full_reference_check)
    check_containers(checks, args.all_containers)
    check_profile(checks, args.profile)

    report = render_report(checks)
    sys.stdout.write(report)
    if args.report:
        report_path = Path(args.report).expanduser()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(report)
        print(f"Report written to: {report_path}", file=sys.stderr)

    return 1 if any(c.status == "FAIL" for c in checks) else 0


if __name__ == "__main__":
    raise SystemExit(main())
