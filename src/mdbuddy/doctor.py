"""`mdbuddy doctor`: check that everything needed for a build is installed.

Run this first. It never installs anything; it reports what is missing and
how to install it so the user can decide.
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import urllib.request
from dataclasses import asdict, dataclass
from typing import Optional

from . import tools
from .fetch import USER_AGENT


@dataclass
class Check:
    name: str
    status: str  # ok | missing | warn
    required: bool
    version: Optional[str] = None
    detail: str = ""
    fix: str = ""


def _ver_tuple(v: Optional[str]) -> tuple:
    out = []
    for part in (v or "").replace("-", ".").split("."):
        if part.isdigit():
            out.append(int(part))
        else:
            break
    return tuple(out)


def check_python() -> Check:
    v = sys.version.split()[0]
    ok = sys.version_info >= (3, 9)
    return Check("python", "ok" if ok else "missing", True, v, sys.executable,
                 "" if ok else "Python >= 3.9 is required")


def check_module(name: str, dist: Optional[str], required: bool, fix: str, purpose: str) -> Check:
    v = tools.module_version(name, dist)
    if v is None:
        return Check(name, "missing", required, None, purpose, fix)
    return Check(name, "ok", required, v, purpose)


def check_martinize2() -> Check:
    exe = shutil.which("martinize2")
    if not exe:
        return Check("martinize2", "missing", True, None, "coarse-grains the protein",
                     "pip install vermouth  (provides the martinize2 command)")
    return Check("martinize2", "ok", True, tools.module_version("vermouth"), exe)


def check_gmx() -> Check:
    exe = tools.find_gmx()
    if not exe:
        return Check("gromacs", "missing", False, None,
                     "only needed to run/check the simulation, not to build the inputs",
                     "conda install -c conda-forge gromacs  (or `module load gromacs` on a cluster)")
    v = tools.gmx_version(exe)
    if _ver_tuple(v) and _ver_tuple(v) < (2021,):
        return Check("gromacs", "warn", False, v, exe,
                     "GROMACS >= 2021 is needed for the C-rescale barostat used in equilibration")
    return Check("gromacs", "ok", False, v, exe)


def check_dssp() -> Check:
    exe = tools.find_dssp()
    has_mdtraj = tools.has_module("mdtraj")
    if not exe and not has_mdtraj:
        return Check("dssp", "missing", True, None, "secondary structure for martinize2",
                     "pip install mdtraj   or   conda install -c salilab dssp  (DSSP 3.0.0)")
    if exe:
        v = tools.dssp_version(exe)
        if v in tools.DSSP_TESTED:
            return Check("dssp", "ok", True, v, exe)
        if has_mdtraj:
            return Check("dssp", "ok", True, v, f"{exe} is DSSP {v} (untested by vermouth); mdtraj will be used")
        return Check("dssp", "warn", True, v, exe,
                     f"vermouth is tested with DSSP {', '.join(tools.DSSP_TESTED)}; DSSP {v} may fail. "
                     "Install mdtraj (pip install mdtraj) as a fallback, or DSSP 3.0.0 (conda install -c salilab dssp)")
    return Check("dssp", "ok", True, tools.module_version("mdtraj"), "using mdtraj for DSSP")


def check_cache() -> Check:
    from .forcefields import cache_dir
    d = cache_dir()
    try:
        d.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=d):
            pass
        return Check("cache", "ok", True, None, str(d))
    except OSError as exc:
        return Check("cache", "missing", True, None, f"{d}: {exc}", "set MDBUDDY_CACHE to a writable directory")


NETWORK_HOSTS = {
    "OPM": "https://opm-assets.storage.googleapis.com/pdb/1bl8.pdb",
    "RCSB": "https://files.rcsb.org/header/1BL8.pdb",
    "GitHub raw (force fields)": "https://raw.githubusercontent.com/marrink-lab/martini-forcefields/main/README.md",
}


def check_network() -> list[Check]:
    out = []
    for label, url in NETWORK_HOSTS.items():
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=10):
                out.append(Check(f"net:{label}", "ok", False, None, url))
        except Exception as exc:  # noqa: BLE001 - report any failure
            out.append(Check(f"net:{label}", "warn", False, None, f"{url}: {exc}",
                             "download the file elsewhere and use protein.source: file / forcefield_files.ff_dir"))
    return out


def run_checks(network: bool = False) -> list[Check]:
    checks = [
        check_python(),
        check_module("yaml", "PyYAML", True, "pip install pyyaml", "reads spec files"),
        check_martinize2(),
        check_module("vermouth", None, True, "pip install vermouth", "martinize2 library"),
        check_module("COBY", "COBY", True, "pip install COBY", "builds membrane, inserts protein, solvates"),
        check_dssp(),
        check_gmx(),
        check_module("mdtraj", None, False, "pip install mdtraj", "optional DSSP backend"),
        check_cache(),
    ]
    if network:
        checks += check_network()
    return checks


def ready(checks: list[Check]) -> bool:
    return all(c.status != "missing" for c in checks if c.required)


def report(checks: list[Check], as_json: bool = False) -> str:
    if as_json:
        return json.dumps({"ready": ready(checks), "checks": [asdict(c) for c in checks]}, indent=2)
    symbol = {"ok": "OK  ", "missing": "MISS", "warn": "WARN"}
    lines = ["MDBuddy environment check", ""]
    for c in checks:
        req = "required" if c.required else "optional"
        lines.append(f"[{symbol[c.status]}] {c.name:<28} {c.version or '':<14} ({req}) {c.detail}")
        if c.fix and c.status != "ok":
            lines.append(f"       -> {c.fix}")
    lines.append("")
    lines.append("Ready to build." if ready(checks) else
                 "Not ready: install the missing required items above, then run `mdbuddy doctor` again.")
    return "\n".join(lines)
