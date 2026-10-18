"""Discovery of external programs and Python packages (used by the pipeline and `mdbuddy doctor`)."""

from __future__ import annotations

import importlib
import importlib.metadata
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

GMX_CANDIDATES = ("gmx", "gmx_mpi", "gmx_d", "gmx_mpi_d")
DSSP_CANDIDATES = ("mkdssp", "dssp")
# vermouth parses classic DSSP output and is tested with these versions
DSSP_TESTED = ("2.2.1", "3.0.0")


class ToolError(RuntimeError):
    pass


def _run(cmd: list[str], timeout: float = 30.0) -> str:
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, timeout=timeout)
        return proc.stdout or ""
    except (OSError, subprocess.SubprocessError):
        return ""


def find_program(name: str) -> Optional[str]:
    """Find an executable on PATH or next to the running Python (an un-activated venv/conda env)."""
    found = shutil.which(name)
    if found:
        return found
    candidate = Path(sys.executable).parent / name
    return str(candidate) if candidate.exists() and os.access(candidate, os.X_OK) else None


def find_gmx() -> Optional[str]:
    env = os.environ.get("MDBUDDY_GMX") or os.environ.get("GMX")
    if env and shutil.which(env):
        return shutil.which(env)
    for name in GMX_CANDIDATES:
        path = shutil.which(name)
        if path:
            return path
    return None


def gmx_version(gmx: str) -> Optional[str]:
    out = _run([gmx, "--version"])
    m = re.search(r"GROMACS version:\s*(\S+)", out)
    return m.group(1) if m else None


def find_dssp() -> Optional[str]:
    env = os.environ.get("MDBUDDY_DSSP")
    if env and shutil.which(env):
        return shutil.which(env)
    for name in DSSP_CANDIDATES:
        path = shutil.which(name)
        if path:
            return path
    return None


def dssp_version(exe: str) -> Optional[str]:
    out = _run([exe, "--version"])
    m = re.search(r"(\d+\.\d+\.\d+)", out)
    return m.group(1) if m else None


def module_version(name: str, dist: Optional[str] = None) -> Optional[str]:
    try:
        return importlib.metadata.version(dist or name)
    except importlib.metadata.PackageNotFoundError:
        pass
    try:
        mod = importlib.import_module(name)
    except Exception:
        return None
    return getattr(mod, "__version__", "unknown")


def has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


def resolve_dssp(setting: str) -> list[str]:
    """Translate ``coarse_grain.dssp`` into martinize2 arguments.

    auto   -> a DSSP 2.x/3.x executable if found, else mdtraj, else a found DSSP 4.x
    mdtraj -> ``-dssp`` without executable (vermouth uses mdtraj)
    none   -> no secondary structure assignment (not recommended for Martini)
    ss:<string> -> ``-ss <string>``
    anything else -> treated as a DSSP executable path/name
    """
    if setting == "none":
        return []
    if setting.startswith("ss:"):
        return ["-ss", setting[3:]]
    if setting == "mdtraj":
        if not has_module("mdtraj"):
            raise ToolError("coarse_grain.dssp is 'mdtraj' but mdtraj is not installed (pip install mdtraj)")
        return ["-dssp"]
    if setting == "auto":
        exe = find_dssp()
        if exe and dssp_version(exe) in DSSP_TESTED:
            return ["-dssp", exe]
        if has_module("mdtraj"):
            return ["-dssp"]
        if exe:
            return ["-dssp", exe]
        raise ToolError("no secondary-structure tool found: install mdtraj (`pip install mdtraj`) or DSSP 3.0.0 "
                        "(`conda install -c salilab dssp`), or set coarse_grain.dssp")
    exe = shutil.which(setting)
    if not exe:
        raise ToolError(f"coarse_grain.dssp: executable '{setting}' not found")
    return ["-dssp", exe]
