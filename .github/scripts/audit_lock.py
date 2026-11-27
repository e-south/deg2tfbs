"""deg2tfbs
deg2tfbs/.github/scripts/audit_lock.py

Audit all external releases in the universal lock or installed project environment.

Module Author(s): Eric J. South
"""

import argparse
import json
import re
import subprocess
import sys
import sysconfig
import tomllib
from collections import defaultdict
from collections.abc import Iterator
from importlib.metadata import distributions
from itertools import zip_longest
from pathlib import Path
from tempfile import TemporaryDirectory

CPU_PACKAGES = {"torch", "torchaudio", "torchvision"}


def locked_releases(project: Path) -> dict[tuple[str, str], str]:
    lock = tomllib.loads((project / "uv.lock").read_text())
    releases = {}
    for package in lock["package"]:
        name, version, source = package["name"], package["version"], package["source"]
        if name == "deg2tfbs" and source == {"editable": "."}:
            continue
        public_version = version
        if source == {"registry": "https://download.pytorch.org/whl/cpu"}:
            # PyPI advisories track the public release of each official CPU build.
            if name not in CPU_PACKAGES or not version.endswith("+cpu"):
                raise ValueError(f"Unrecognized CPU build: {name} {version}")
            public_version = version.removesuffix("+cpu")
        elif source != {"registry": "https://pypi.org/simple"} or "+" in version:
            raise ValueError(f"Unrecognized release: {name} {version} from {source}")
        releases[name, version] = public_version
    return releases


def installed_releases(project: Path, locked: dict[tuple[str, str], str]) -> Iterator[tuple[str, str]]:
    if Path(sys.prefix).resolve() != (project / ".venv").resolve():
        raise ValueError("Run the installed audit with the project's .venv interpreter")

    paths = {sysconfig.get_path(key) for key in ("purelib", "platlib")}
    found_project = False
    for dist in distributions(path=sorted(paths)):
        name = re.sub(r"[-_.]+", "-", dist.metadata["Name"]).lower()
        version = dist.version
        origin = json.loads(dist.read_text("direct_url.json") or "{}")
        if name == "deg2tfbs":
            if origin.get("url") != project.as_uri() or "dir_info" not in origin:
                raise ValueError(f"Expected the project built from {project}, got {origin}")
            found_project = True
            continue
        if origin:
            raise ValueError(f"Unrecognized direct install: {name} {version} from {origin}")
        if "+" in version or (name in CPU_PACKAGES and sys.platform == "linux"):
            if (
                sys.platform != "linux"
                or name not in CPU_PACKAGES
                or not version.endswith("+cpu")
                or (name, version) not in locked
            ):
                raise ValueError(f"Unrecognized installed build: {name} {version}")
            version = locked[name, version]
        yield name, version
    if not found_project:
        raise ValueError(f"Expected an installed copy of {project}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installed", action="store_true", help="audit this project's .venv")
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[2]
    locked = locked_releases(project)
    if args.installed:
        releases = installed_releases(project, locked)
    else:
        releases = ((name, public) for (name, _), public in locked.items())

    versions = defaultdict(set)
    for name, version in releases:
        versions[name].add(f"{name}=={version}")
    count = sum(map(len, versions.values()))
    source = "the installation" if args.installed else "the lock"
    print(f"Auditing {count} external releases from {source}", flush=True)

    # pip-audit requires at most one version per name in each invocation.
    failed = False
    with TemporaryDirectory() as directory:
        requirements = Path(directory) / "requirements.txt"
        for batch in zip_longest(*(sorted(pins) for pins in versions.values())):
            requirements.write_text("\n".join(pin for pin in batch if pin is not None) + "\n")
            result = subprocess.run([
                "uvx", "pip-audit@2.10.1", "--strict", "--no-deps", "--disable-pip",
                "-r", str(requirements),
            ])
            if result.returncode:
                failed = True
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
