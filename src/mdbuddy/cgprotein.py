"""Pre-coarse-grained proteins (``cg_proteins:``): structures that are already Martini models.

Typical use: protein islands or clusters taken from an earlier simulation, which should be
inserted into a new (larger, or differently composed) membrane exactly as arranged. Nothing
is re-coarse-grained or re-oriented. MDBuddy copies the structure and its itps into the build,
checks that the bead count and bead names match the topology, and hands the structure to COBY:

* ``placement: keep``   COBY's ``center_protein:False``: coordinates are kept relative to the
                        structure's own box centre (gro box line or PDB CRYST1), so with
                        ``membrane.box: from_structure`` they are kept exactly.
* ``placement: center`` the structure's centre of geometry goes to the box centre (+ cx/cy/cz).
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Optional

from . import topology
from .provenance import sha256
from .structure import read_gro, read_pdb


class CGProteinError(RuntimeError):
    pass


def read_beads(path: Path) -> tuple[list[tuple[int, str, str, float, float, float]], Optional[list[float]]]:
    """[(resid, resname, name, x, y, z)] in nm and the box in nm (or None) of a CG .gro or .pdb."""
    path = Path(path)
    if path.suffix.lower() == ".gro":
        _, atoms, box = read_gro(path)  # all nine values of a triclinic box are kept
        return [(a.resid, a.resname, a.name, a.x, a.y, a.z) for a in atoms], box or None
    atoms, header = read_pdb(path)
    box = None
    for line in header:
        if line.startswith("CRYST1"):
            try:
                box = [float(line[6:15]) / 10, float(line[15:24]) / 10, float(line[24:33]) / 10]
            except ValueError:
                box = None
            if box and min(box) <= 0.11:  # placeholder CRYST1 1 1 1
                box = None
    return [(a.resid, a.resname, a.name, a.x / 10, a.y / 10, a.z / 10) for a in atoms], box


def read_structure(path: Path) -> tuple[list[str], list[tuple[float, float, float]], Optional[list[float]]]:
    """(bead names, coordinates in nm, box in nm or None) of a CG .gro or .pdb."""
    beads, box = read_beads(path)
    return [b[2] for b in beads], [b[3:] for b in beads], box


def write_coby_pdb(path: Path, beads, box: Optional[list[float]]) -> None:
    """The structure COBY reads: PDB with CRYST1. COBY 1.0.15's .gro reader skips the box line, so
    center_protein:False would ignore a .gro box; its PDB reader uses CRYST1 correctly."""
    lines = []
    if box:
        lines.append("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1" % (box[0] * 10, box[1] * 10, box[2] * 10,
                                                                               90.0, 90.0, 90.0))
    for i, (resid, resname, name, x, y, z) in enumerate(beads, start=1):
        if len(name) > 4 or len(resname) > 4:
            raise CGProteinError(f"bead {i} ({resname} {name}): names longer than 4 characters do not fit a PDB file")
        field = f" {name:<3}" if len(name) < 4 else name
        lines.append("ATOM  %5d %4s %-4s%5d    %8.3f%8.3f%8.3f  1.00  0.00" % (
            i % 100000, field, resname, resid % 10000, x * 10, y * 10, z * 10))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def suggest_timestep(max_fs: float, output_ps: float) -> Optional[float]:
    """Largest sensible time step at or below max_fs whose energy interval still divides output_ps.

    simulation.output_ps has to be a whole number of nstcalcenergy (100) steps, so suggesting a
    time step that does not divide it would just trade one rejected spec for another.
    """
    for ts in (20.0, 15.0, 10.0, 8.0, 5.0, 4.0, 2.0, 1.0):
        if ts > max_fs:
            continue
        steps = output_ps / (100 * ts / 1000.0)
        if abs(steps - round(steps)) < 1e-9:
            return ts
    return None


def stiff_bond_warning(info: topology.TopologyInfo, molecules, sim,
                       ff_masses: Optional[dict] = None) -> Optional[str]:
    """Warn when a supplied topology's bonds are too stiff for the chosen time step.

    martinize2 writes [ constraints ] for the Martini protein backbone, but hand-made itps
    sometimes use very stiff bonds instead (e.g. 150000 kJ/mol/nm^2 in an OMP island built for
    an earlier study). grompp then emits a *warning* about the oscillational period, which fails
    the build's own -maxwarn 0 check, and its message does not say what to do. Say it here, with
    the safe maximum and a time step that actually validates.
    """
    masses = dict(ff_masses or {})
    masses.update(info.atomtype_masses)  # the entry's own itps win
    worst = None
    for name, _count in molecules:
        mt = info.moleculetypes.get(name)
        if mt is None:
            continue
        found = topology.fastest_bond(mt, masses)
        if found and (worst is None or found[3] < worst[1][3]):
            worst = (name, found)
    if worst is None:
        return None
    name, (i, j, k, period_ps) = worst
    timestep_fs = float(sim["timestep_fs"])
    safe_fs = period_ps / 5.0 * 1000.0
    if timestep_fs <= safe_fs + 1e-9:
        return None
    pick = suggest_timestep(safe_fs, float(sim["output_ps"]))
    fix = (f"set simulation.timestep_fs: {pick:g}" if pick else
           f"lower simulation.timestep_fs below {safe_fs:.1f} (and adjust simulation.output_ps "
           "to a whole number of 100-step intervals)")
    return (f"{name} bond {i}-{j} has k = {k:g} kJ/mol/nm2, an oscillational period of "
            f"{period_ps:.3f} ps, so GROMACS wants a time step <= {safe_fs:.1f} fs but "
            f"simulation.timestep_fs is {timestep_fs:g}. grompp will warn and the -maxwarn 0 "
            f"check will fail. Either {fix}, or rebuild the itps with [ constraints ] for the "
            "backbone as martinize2 does.")


def box_from_structure(spec) -> list[float]:
    """Box (nm) of the first cg_proteins structure: membrane.box: from_structure."""
    first = spec.resolve_path(spec.cg_proteins[0]["structure"])
    _, _, box = read_structure(first)
    if not box:
        raise CGProteinError(f"membrane.box: from_structure, but {first.name} has no box "
                             "(use a .gro, or a PDB with a CRYST1 record)")
    if len(box) > 3 and any(abs(v) > 1e-6 for v in box[3:]):
        raise CGProteinError(f"{first.name} has a triclinic box; MDBuddy builds rectangular boxes only. "
                             "Give membrane.box: [x, y, z] instead.")
    return [round(float(v), 3) for v in box[:3]]


def _itps_from_top(top: Path) -> tuple[list[Path], list[tuple[str, int]]]:
    """Protein itps (those defining a molecule listed in [ molecules ]) and the molecules of a .top."""
    molecules = topology.read_molecules(top)
    names = {n for n, _ in molecules}
    itps = []
    for inc in topology.read_includes(top):
        path = (top.parent / inc).resolve()
        if not path.exists():
            continue  # force-field files the .top was built against; MDBuddy brings its own
        info = topology.parse(path)
        if set(info.moleculetypes) & names:
            itps.append(path)
    return itps, molecules


def prepare(spec, work: Path, toppar: Path, existing_itps: list[str]) -> dict:
    """Copy every cg_proteins entry into ``work`` (02_cg) and its itps into ``toppar``; check them."""
    work.mkdir(parents=True, exist_ok=True)
    toppar.mkdir(parents=True, exist_ok=True)
    taken = {name: None for name in existing_itps}  # itp file name -> sha256 (None: not ours)
    # bead masses from the force-field files already in toppar, for the stiff-bond check
    ff_info = topology.TopologyInfo()
    for name in existing_itps:
        if (toppar / name).exists():
            topology.parse(toppar / name, info=ff_info)
    ff_masses = ff_info.atomtype_masses
    moltype_source: dict[str, str] = {}
    entries, all_itps, records = [], [], []
    for i, e in enumerate(spec.cg_proteins):
        where = f"cg_proteins[{i}]"
        src = spec.resolve_path(e["structure"])
        if not src.exists():
            raise CGProteinError(f"{where}: structure not found: {src}")
        if e.get("topology"):
            top = spec.resolve_path(e["topology"])
            if not top.exists():
                raise CGProteinError(f"{where}: topology not found: {top}")
            itps, molecules = _itps_from_top(top)
            if not molecules:
                raise CGProteinError(f"{where}: {top.name} has no [ molecules ] section")
            if e.get("itps"):
                itps = [spec.resolve_path(p) for p in e["itps"]]
            if e.get("moleculetypes"):
                molecules = [(n, int(c)) for n, c in e["moleculetypes"]]
        else:
            itps = [spec.resolve_path(p) for p in e["itps"]]
            molecules = [(n, int(c)) for n, c in e["moleculetypes"]]
        missing = [str(p) for p in itps if not p.exists()]
        if missing:
            raise CGProteinError(f"{where}: itp not found: {', '.join(missing)}")

        info = topology.TopologyInfo()
        for p in itps:
            topology.parse(p, info=info)
            nested = topology.read_includes(p)
            if nested:
                raise CGProteinError(f"{where}: {p.name} #includes {nested}; list the included itps "
                                     "directly under itps instead")
        undefined = [n for n, _ in molecules if n not in info.moleculetypes]
        if undefined:
            raise CGProteinError(f"{where}: moleculetypes {undefined} are not defined in {[p.name for p in itps]}")

        names, xyz, box = read_structure(src)
        expected = [a for n, c in molecules for _ in range(c) for a in info.moleculetypes[n].atom_names]
        if len(expected) != len(names):
            raise CGProteinError(f"{where}: {src.name} has {len(names)} beads but the moleculetypes "
                                 f"{molecules} have {len(expected)}")
        bad = next((k for k, (a, b) in enumerate(zip(expected, names)) if a[:5] != b[:5]), None)
        if bad is not None:
            raise CGProteinError(f"{where}: bead {bad + 1} is '{names[bad]}' in {src.name} but '{expected[bad]}' "
                                 "in the topology; the structure and moleculetypes must be in the same order")
        if e["placement"] == "keep" and not box:
            raise CGProteinError(f"{where}: placement 'keep' needs the structure's box (a .gro, or a PDB "
                                 "with CRYST1); COBY keeps coordinates relative to that box centre")

        shutil.copyfile(src, work / f"cg_protein_{i}_input{src.suffix.lower()}")
        dest = work / f"cg_protein_{i}.pdb"
        write_coby_pdb(dest, read_beads(src)[0], box[:3] if box else None)
        copied = []
        for p in itps:
            digest = sha256(p)
            name = p.name
            if name in taken and taken[name] != digest:
                name = f"cg{i}_{p.name}"
            if name not in taken:
                shutil.copyfile(p, toppar / name)
                taken[name] = digest
                all_itps.append(name)
            copied.append(name)
            for mt in topology.parse(p).moleculetypes:
                if mt in moltype_source and moltype_source[mt] != digest:
                    raise CGProteinError(f"moleculetype {mt} is defined differently in two cg_proteins itps; "
                                         "rename one of them (moleculetype names must be unique in a topology)")
                moltype_source[mt] = digest
        has_posres = any("POSRES" in p.read_text() for p in itps)
        stiff = stiff_bond_warning(info, molecules, spec.sim, ff_masses)
        records.append({"structure": str(src), "sha256": sha256(src), "copied_to": dest.name,
                        "posres_in_itps": has_posres, "stiff_bond_warning": stiff,
                        "itps": [{"file": n, "sha256": taken[n]} for n in copied],
                        "molecules": molecules, "beads": len(names), "box_nm": box[:3] if box else None,
                        "placement": e["placement"]})
        entries.append({"structure": dest, "molecules": molecules, "xyz": xyz, "box": box, **{
            k: e[k] for k in ("placement", "cx", "cy", "cz", "rx", "ry", "rz", "coby_args")}})
    return {"entries": entries, "itps": all_itps,
            "molecules": [m for e in entries for m in e["molecules"]], "records": records}


def entries_from_state(spec, work: Path, state_entries: list[dict]) -> list[dict]:
    """Re-read the copied structures (for --from membrane)."""
    out = []
    for e in state_entries:
        path = work / e["structure"]
        names, xyz, box = read_structure(path)
        out.append({**e, "structure": path, "xyz": xyz, "box": box,
                    "molecules": [tuple(m) for m in e["molecules"]]})
    return out


def coby_arg(entry: dict, rel: str) -> str:
    """COBY protein argument for one pre-coarse-grained structure."""
    parts = [f"file:{rel}", *(f"moleculetypes:{n}:{c}" for n, c in entry["molecules"])]
    if entry["placement"] == "keep":
        parts.append("center_protein:False")
    for k in ("cx", "cy", "cz", "rx", "ry", "rz"):
        if float(entry[k]):
            parts.append(f"{k}:{float(entry[k]):g}")
    if entry["coby_args"]:
        parts.append(entry["coby_args"])
    return " ".join(parts)
