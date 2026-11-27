"""deg2tfbs
deg2tfbs/.github/scripts/audit_lock.py

Audit every external release in the universal uv lock against PyPI.

Module Author(s): Eric J. South
"""

import subprocess
import tomllib
from collections import defaultdict
from itertools import zip_longest
from pathlib import Path
from tempfile import TemporaryDirectory


def main() -> None:
    lock = Path(__file__).resolve().parents[2] / "uv.lock"
    packages = tomllib.loads(lock.read_text())["package"]
    versions = defaultdict(set)
    for package in packages:
        name, version, source = package["name"], package["version"], package["source"]
        if name == "deg2tfbs" and source == {"editable": "."}:
            continue
        if source == {"registry": "https://download.pytorch.org/whl/cpu"}:
            # PyPI advisories track the public release of each official CPU build.
            if name not in {"torch", "torchaudio", "torchvision"} or not version.endswith("+cpu"):
                raise ValueError(f"Unrecognized CPU build: {name} {version}")
            version = version.removesuffix("+cpu")
        elif source != {"registry": "https://pypi.org/simple"}:
            raise ValueError(f"Unrecognized source for {name}: {source}")
        versions[name].add(f"{name}=={version}")

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
