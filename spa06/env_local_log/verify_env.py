#!/usr/bin/env python3
# %%
import json
import os
from pathlib import Path
import subprocess
import sys

# %%
SOURCE_ROOT = Path("/software/miniconda3/envs")
TARGET_ROOT = Path("/home/lcpyspa/miniconda3/envs")
MODULES = {
    "commot_env": [
        "commot", "scanpy", "anndata", "numpy", "pandas", "scipy",
        "sklearn", "matplotlib", "seaborn", "numba",
    ],
    "pyscenic": [
        "pyscenic", "scanpy", "anndata", "numpy", "pandas", "scipy",
        "sklearn", "matplotlib", "arboreto", "ctxcore",
    ],
    "spatrack_env": [
        "spaTrack", "scanpy", "anndata", "numpy", "pandas", "scipy",
        "sklearn", "matplotlib", "torch",
    ],
    "steropy_env": [
        "stereo", "scanpy", "anndata", "numpy", "pandas", "scipy",
        "sklearn", "matplotlib", "seaborn",
    ],
    "stlearn": [
        "stlearn", "scanpy", "anndata", "numpy", "pandas", "scipy",
        "sklearn", "matplotlib", "seaborn",
    ],
}


# %%
def load_records(prefix):
    records = {}
    channels = {}
    for path in sorted((prefix / "conda-meta").glob("*.json")):
        with path.open(encoding="utf-8") as handle:
            record = json.load(handle)
        key = (
            record.get("name"),
            record.get("version"),
            record.get("build"),
        )
        records[key] = path.name
        channels[key] = record.get("channel")
    return records, channels


def run(command):
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=600,
        env={**os.environ, "PYTHONNOUSERSITE": "1"},
    )
    return result.returncode, result.stdout.strip()


def python_report(python, modules):
    code = """
import importlib
import json
import sys

reports = {}
for module_name in sys.argv[1:]:
    try:
        module = importlib.import_module(module_name)
        reports[module_name] = {
            "ok": True,
            "version": str(getattr(module, "__version__", "")),
        }
    except Exception as exc:
        reports[module_name] = {
            "ok": False,
            "error": "{}: {}".format(type(exc).__name__, exc),
        }
print(json.dumps({
    "python": sys.version,
    "executable": sys.executable,
    "imports": reports,
}, indent=2, sort_keys=True))
"""
    return run([str(python), "-c", code, *modules])


def pip_report(python):
    return run([str(python), "-m", "pip", "list", "--format=json"])


# %%
def main():
    env_name = sys.argv[1]
    source = SOURCE_ROOT / env_name
    target = TARGET_ROOT / env_name
    source_records, source_channels = load_records(source)
    target_records, target_channels = load_records(target)
    missing = sorted(source_records.keys() - target_records.keys())
    extra = sorted(target_records.keys() - source_records.keys())
    channel_differences = sorted(
        key
        for key in source_records.keys() & target_records.keys()
        if source_channels[key] != target_channels[key]
    )
    source_pip_list_rc, source_pip_list_output = pip_report(source / "bin" / "python")
    target_pip_list_rc, target_pip_list_output = pip_report(target / "bin" / "python")
    source_pip_packages = {
        item["name"].lower(): item["version"]
        for item in json.loads(source_pip_list_output)
    }
    target_pip_packages = {
        item["name"].lower(): item["version"]
        for item in json.loads(target_pip_list_output)
    }
    pip_missing = sorted(source_pip_packages.keys() - target_pip_packages.keys())
    pip_extra = sorted(target_pip_packages.keys() - source_pip_packages.keys())
    pip_version_differences = sorted(
        name
        for name in source_pip_packages.keys() & target_pip_packages.keys()
        if source_pip_packages[name] != target_pip_packages[name]
    )

    python_rc, python_output = python_report(
        target / "bin" / "python",
        MODULES[env_name],
    )
    source_pip_rc, source_pip = run(
        [str(source / "bin" / "python"), "-m", "pip", "check"]
    )
    target_pip_rc, target_pip = run(
        [str(target / "bin" / "python"), "-m", "pip", "check"]
    )

    print("ENV", env_name)
    print("source:", source)
    print("target:", target)
    print("source_package_count:", len(source_records))
    print("target_package_count:", len(target_records))
    print("missing_package_count:", len(missing))
    print("extra_package_count:", len(extra))
    print("channel_metadata_mismatch_count:", len(channel_differences))
    print("source_pip_package_count:", len(source_pip_packages))
    print("target_pip_package_count:", len(target_pip_packages))
    print("pip_missing_packages:", pip_missing)
    print("pip_extra_packages:", pip_extra)
    print("pip_version_differences:", [
        {
            "name": name,
            "source": source_pip_packages[name],
            "target": target_pip_packages[name],
        }
        for name in pip_version_differences
    ])
    if missing:
        print("missing_packages:", json.dumps(missing, indent=2))
    if extra:
        print("extra_packages:", json.dumps(extra, indent=2))
    print("package_set_match:", not missing and not extra)
    print("python_check_rc:", python_rc)
    print(python_output)
    print("source_pip_check_rc:", source_pip_rc)
    print(source_pip)
    print("target_pip_check_rc:", target_pip_rc)
    print(target_pip)
    print("pip_check_matches_source:", source_pip_rc == target_pip_rc and source_pip == target_pip)

    package_ok = not missing and not extra
    pip_versions_ok = (
        source_pip_list_rc == 0
        and target_pip_list_rc == 0
        and not pip_missing
        and not pip_extra
        and not pip_version_differences
    )
    imports_ok = (
        python_rc == 0
        and '"ok": false' not in python_output.lower()
    )
    pip_ok = target_pip_rc == 0 or source_pip_rc == target_pip_rc
    overall_ok = package_ok and pip_versions_ok and imports_ok and pip_ok
    print("migration_check:", "PASS" if overall_ok else "FAIL")
    return 0 if overall_ok else 1


# %%
if __name__ == "__main__":
    raise SystemExit(main())
