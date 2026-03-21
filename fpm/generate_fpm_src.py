#!/usr/bin/env python3
"""Generate an fpm-buildable source tree in fpm-src/."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys


FPM_DIR = Path(__file__).resolve().parent
REPO_ROOT = FPM_DIR.parent
OUT = FPM_DIR / "fpm-src"
SRC_IN = REPO_ROOT / "src"
SRC_OUT = OUT / "src"
EXAMPLE_IN = REPO_ROOT / "example"
EXAMPLE_OUT = OUT / "example"
LIB_INFO_TEMPLATE = SRC_IN / "fvtkhdf_lib_info.F90.in"

SRC_COPY_FILES = (
    "vtkhdf_assert.F90",
    "vtkhdf_assert.inc",
    "vtkhdf_vtk_cell_types.F90",
    "vtkhdf_h5_c_binding.F90",
    "vtkhdf_version_param.F90",
    "vtkhdf_h5e_context_type.F90",
    "vtkhdf_h5_c_shim.c",
)

FYPP_FILES = (
    "vtkhdf_h5.F90.fypp",
    "vtkhdf_ug_type.F90.fypp",
    "vtkhdf_ug_file_type.F90.fypp",
    "vtkhdf_mb_file_type.F90.fypp",
    "vtkhdf_ctx_type.F90.fypp",
)

EXAMPLE_GLOBS = (
    "*.F90",
    "*.f90",
)

EXAMPLE_COPY_FILES = (
    "README.md",
)


def copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def run_fypp_with_args(src: Path, dst: Path, mpi: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["fypp"]
    if mpi:
        cmd.append("-DUSE_MPI")
    cmd.extend((str(src), str(dst)))
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def load_version_fields() -> dict[str, str]:
    manifest_text = (FPM_DIR / "fpm.toml").read_text()
    match = re.search(r'^\s*version\s*=\s*"([^"]+)"\s*$', manifest_text, re.MULTILINE)
    if match is None:
        raise ValueError("could not find version in fpm.toml")

    version = match.group(1)
    parts = version.split(".")
    if len(parts) != 3 or not all(part.isdigit() for part in parts):
        raise ValueError(f"unsupported version format in fpm.toml: {version!r}")

    major, minor, patch = parts
    return {
        "@PROJECT_VERSION@": version,
        "@PROJECT_VERSION_MAJOR@": major,
        "@PROJECT_VERSION_MINOR@": minor,
        "@PROJECT_VERSION_PATCH@": patch,
    }


def render_lib_info(dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    text = LIB_INFO_TEMPLATE.read_text()
    for old, new in load_version_fields().items():
        text = text.replace(old, new)
    dst.write_text(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an fpm-buildable source tree in fpm-src/."
    )
    parser.add_argument(
        "--features",
        default="",
        help="comma-separated fpm feature names to specialize generated sources",
    )
    parser.add_argument(
        "--profile",
        default="",
        help="fpm profile name to mirror when specializing generated sources",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="deprecated alias for --features mpi",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    features = {item.strip() for item in args.features.split(",") if item.strip()}
    if "mpi" in args.profile:
        features.add("mpi")
    if args.mpi:
        features.add("mpi")
    mpi_enabled = any("mpi" in feature for feature in features)

    if shutil.which("fypp") is None:
        print("error: fypp not found in PATH", file=sys.stderr)
        return 1

    if not (FPM_DIR / "fpm.toml").is_file():
        print("error: fpm.toml not found in fpm/", file=sys.stderr)
        return 1

    if OUT.exists():
        shutil.rmtree(OUT)

    SRC_OUT.mkdir(parents=True)

    copy_file(FPM_DIR / "fpm.toml", OUT / "fpm.toml")
    copy_file(FPM_DIR / "README.md", OUT / "README.md")
    copy_file(REPO_ROOT / "LICENSE.md", OUT / "LICENSE.md")

    for name in SRC_COPY_FILES:
        copy_file(SRC_IN / name, SRC_OUT / name)

    for name in FYPP_FILES:
        src = SRC_IN / name
        dst = SRC_OUT / name.removesuffix(".fypp")
        run_fypp_with_args(src, dst, mpi=mpi_enabled)

    render_lib_info(SRC_OUT / "fvtkhdf_lib_info.F90")

    for pattern in EXAMPLE_GLOBS:
        for src in sorted(EXAMPLE_IN.glob(pattern)):
            copy_file(src, EXAMPLE_OUT / src.name)

    for name in EXAMPLE_COPY_FILES:
        copy_file(EXAMPLE_IN / name, EXAMPLE_OUT / name)

    feature_note = f" ({','.join(sorted(features))})" if features else ""
    mode = " (MPI)" if mpi_enabled and not feature_note else feature_note
    print(f"Generated {OUT.relative_to(FPM_DIR)}/{mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
