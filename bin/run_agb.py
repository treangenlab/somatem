#!/usr/bin/env python3
"""Run legacy AGB reproducibly from a Nextflow task."""

import argparse
import gzip
import os
from pathlib import Path
import shutil
import subprocess

import agb_src
import pkg_resources


def copy_input(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".gz":
        with gzip.open(str(source), "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(str(source), str(destination))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--info", required=True, type=Path)
    parser.add_argument("--assembler", default="flye")
    parser.add_argument("--threads", required=True, type=int)
    args = parser.parse_args()

    input_dir = Path("agb_input")
    copy_input(args.gfa, input_dir / "assembly_graph.gfa")
    copy_input(args.graph, input_dir / "assembly_graph.gv")
    copy_input(args.info, input_dir / "assembly_info.txt")

    package_copy = Path("agb_src")
    if package_copy.exists():
        shutil.rmtree(str(package_copy))
    shutil.copytree(str(Path(agb_src.__file__).resolve().parent), str(package_copy))

    utils = package_copy / "scripts" / "utils.py"
    utils.write_text(
        utils.read_text(encoding="utf-8").replace("os.io.open", "open"),
        encoding="utf-8",
    )

    executable = shutil.which("agb.py")
    if not executable:
        raise RuntimeError("agb.py was not found in PATH")

    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(Path.cwd())
    subprocess.run(
        [
            executable,
            "-a",
            args.assembler,
            "-i",
            str(input_dir),
            "-t",
            str(args.threads),
        ],
        check=True,
        env=environment,
    )

    viewer = Path("agb_output/viewer.html")
    if not viewer.is_file() or viewer.stat().st_size == 0:
        raise RuntimeError("AGB did not create agb_output/viewer.html")

    version = pkg_resources.get_distribution("agb").version
    Path("versions.yml").write_text(
        'AGB:\n    agb: "{}"\n'.format(version), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
