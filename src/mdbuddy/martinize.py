"""Coarse-grain the protein with martinize2 (vermouth)."""

from __future__ import annotations

import shutil
from pathlib import Path

from . import topology
from .tools import find_program, resolve_dssp


class MartinizeError(RuntimeError):
    pass


def merge_args(spec, chains: list[str]) -> list[str]:
    mc = spec.cg["merge_chains"]
    if mc == "auto":
        # An elastic network / Go model only acts within a molecule, so an oligomer
        # must be one molecule or its subunits can drift apart.
        if len(chains) > 1 and spec.cg["model"] in ("elastic", "go"):
            return ["-merge", ",".join(chains)]
        return []
    if mc == "all":
        return ["-merge", ",".join(chains)] if len(chains) > 1 else []
    if mc == "none":
        return []
    args = []
    for group in mc:
        args += ["-merge", ",".join(group)]
    return args


def build_command(spec, in_pdb: str, chains: list[str], dssp: list[str]) -> list[str]:
    cg = spec.cg
    cmd = [
        find_program("martinize2") or "martinize2",
        "-f", in_pdb,
        "-x", "cg.pdb",
        "-o", "cg.top",
        "-ff", spec.ff_info["martinize_ff"],
        "-name", cg["molecule_name"],
        "-ignh",
        "-cys", str(cg["cys"]),
        "-p", cg["position_restraints"],
        "-pf", f"{cg['posres_force']:g}",
        "-maxwarn", str(cg["maxwarn"]),
    ]
    cmd += dssp
    if cg["model"] == "elastic":
        e = cg["elastic"]
        cmd += ["-elastic", "-ef", f"{e['force']:g}", "-el", f"{e['lower']:g}", "-eu", f"{e['upper']:g}"]
    elif cg["model"] == "go":
        g = cg["go"]
        cmd += ["-go", "-go-eps", f"{g['epsilon']:g}", "-go-low", f"{g['low']:g}",
                "-go-up", f"{g['up']:g}", "-go-res-dist", str(g["res_dist"])]
    cmd += merge_args(spec, chains)
    cmd += [str(a) for a in cg["extra_args"]]
    return cmd


def run(spec, lock, entry, clean_pdb: Path, chains: list[str], work: Path, toppar: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    local_in = work / "protein_clean.pdb"
    shutil.copyfile(clean_pdb, local_in)
    dssp = resolve_dssp(spec.cg["dssp"])
    cmd = build_command(spec, local_in.name, chains, dssp)
    (work / "martinize2_cmd.sh").write_text("#!/bin/bash\n# Re-run exactly what MDBuddy ran\n"
                                            + " ".join(_q(c) for c in cmd) + "\n")
    lock.run(entry, cmd, cwd=work, log_name="martinize2.log")

    cg_pdb, cg_top = work / "cg.pdb", work / "cg.top"
    if not cg_pdb.exists() or not cg_top.exists():
        raise MartinizeError("martinize2 finished but cg.pdb / cg.top are missing; see martinize2.log")

    nan_beads = [l[12:26].strip() for l in cg_pdb.read_text().splitlines()
                 if l.startswith(("ATOM", "HETATM")) and "nan" in l[30:54].lower()]
    if nan_beads:
        raise MartinizeError(f"martinize2 wrote nan coordinates for {len(nan_beads)} bead(s) (e.g. {nan_beads[:3]}): "
                             "their atoms are missing from the input. Keep protein.repair: auto so they are rebuilt.")
    molecules = topology.read_molecules(cg_top)
    itps = [inc for inc in topology.read_includes(cg_top) if not inc.startswith("martini")]
    if not molecules:
        raise MartinizeError("martinize2 wrote no [ molecules ]; see martinize2.log")
    copied = []
    for itp in itps:
        src = work / itp
        if not src.exists():
            raise MartinizeError(f"martinize2 topology includes {itp} but the file is missing")
        shutil.copyfile(src, toppar / src.name)
        copied.append(toppar / src.name)
    go_files = []
    if spec.cg["model"] == "go":
        for name in ("go_atomtypes.itp", "go_nbparams.itp"):
            src = work / name
            if not src.exists():
                raise MartinizeError(f"Go model requested but martinize2 did not write {name}")
            shutil.copyfile(src, toppar / name)
            go_files.append(toppar / name)

    warnings = [l for l in (work / "martinize2.log").read_text().splitlines() if "WARNING" in l.upper()]
    return {
        "cg_pdb": cg_pdb,
        "molecules": molecules,
        "itps": [p.name for p in copied],
        "go_files": [p.name for p in go_files],
        "outputs": [cg_pdb, cg_top, *copied, *go_files],
        "warnings": warnings,
        "dssp": dssp,
    }


def _q(s: str) -> str:
    import shlex
    return shlex.quote(s)
