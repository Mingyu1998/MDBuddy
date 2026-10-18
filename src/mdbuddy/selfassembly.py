"""Did the bilayer self-assemble? Where is the protein? (`mdbuddy selfassembly check`)

Reads the last frame of a self-assembly run, tests each box axis for a complete periodic
bilayer (tail ends in the core, head groups in two leaflets, full in-plane coverage), and maps
the bilayer normal and midplane back to the frame of the original input structure, so the system
can be rebuilt with COBY (any composition, including asymmetric) using ``orientation: manual``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional

import yaml

from . import topology
from .structure import read_gro, read_pdb

HEAD_NAMES = ("PO4", "PO41", "ROH")  # phosphate (or cholesterol hydroxyl) marks the head group
CORE_HALF_NM = 1.2  # tail ends within this distance of the midplane count as "in the core"
HEAD_RANGE_NM = (1.2, 3.2)  # head groups this far from the midplane count as "in a leaflet"


class SelfAssemblyError(RuntimeError):
    pass


def _np():
    import numpy as np
    return np


def _kabsch(P, Q):
    """R, t minimising |R P + t - Q| (rows are points)."""
    np = _np()
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def _unwrap_chain(xyz, box):
    """Make consecutive points whole under PBC (each relative to the previous one)."""
    np = _np()
    out = xyz.copy()
    for i in range(1, len(out)):
        delta = out[i] - out[i - 1]
        out[i] = out[i - 1] + delta - box * np.round(delta / box)
    return out


def _molecule_slices(molecules, info):
    start = 0
    for name, count in molecules:
        n = info.moleculetypes[name].natoms
        for _ in range(count):
            yield name, start, start + n
            start += n


def analyse(gro: Path, top: Path, lipid_names: list[str], protein_moltypes: list[str]) -> dict:
    np = _np()
    info = topology.parse(top)
    molecules = info.molecules
    _, atoms, box = read_gro(gro)
    L = np.array(box[:3])
    xyz = np.array([[a.x, a.y, a.z] for a in atoms])

    heads, ends, lipid_beads, lip_names = [], [], [], []
    bb = []
    for name, s, e in _molecule_slices(molecules, info):
        names = info.moleculetypes[name].atom_names
        if name in lipid_names:
            mol = xyz[s:e]
            h = next((i for i, n in enumerate(names) if n in HEAD_NAMES), 0)
            rel = mol - mol[h]
            rel -= L * np.round(rel / L)
            far = int(np.argmax(np.linalg.norm(rel, axis=1)))
            heads.append(mol[h])
            ends.append(mol[h] + rel[far])
            lipid_beads.append(mol)
            lip_names.append(name)
        elif name in protein_moltypes:
            bb += [s + i for i, n in enumerate(names) if n == "BB"]
    if not heads:
        raise SelfAssemblyError("no lipids found in the topology")
    heads, ends = np.array(heads), np.array(ends)
    allbeads = np.vstack(lipid_beads)
    prot = xyz[bb] if bb else np.zeros((0, 3))

    axes = []
    for k in range(3):
        ang = 2 * math.pi * allbeads[:, k] / L[k]
        c = (math.atan2(np.sin(ang).mean(), np.cos(ang).mean()) / (2 * math.pi)) * L[k] % L[k]
        dh = heads[:, k] - c
        dh -= L[k] * np.round(dh / L[k])
        de = ends[:, k] - c
        de -= L[k] * np.round(de / L[k])
        core = float((np.abs(de) < CORE_HALF_NM).mean())
        leaf = float(((np.abs(dh) > HEAD_RANGE_NM[0]) & (np.abs(dh) < HEAD_RANGE_NM[1])).mean())
        o = [i for i in range(3) if i != k]
        da = allbeads[:, k] - c
        da -= L[k] * np.round(da / L[k])
        pts = [allbeads[np.abs(da) < 2.0][:, o]]
        if len(prot):
            dp = prot[:, k] - c
            dp -= L[k] * np.round(dp / L[k])
            pts.append(prot[np.abs(dp) < 2.0][:, o])
        pts = np.vstack(pts)
        nb = [max(1, int(L[i] // 1.0)) for i in o]
        hist, _, _ = np.histogram2d(pts[:, 0] % L[o[0]], pts[:, 1] % L[o[1]], bins=nb,
                                    range=[[0, L[o[0]]], [0, L[o[1]]]])
        coverage = float((hist > 0).mean())
        up, lo = dh[dh > 0], dh[dh < 0]
        axes.append({"axis": "xyz"[k], "centre_nm": float(c), "tail_ends_in_core": round(core, 3),
                     "heads_in_leaflets": round(leaf, 3), "in_plane_coverage": round(coverage, 3),
                     "score": core * leaf * coverage,
                     "thickness_PP_nm": round(float(up.mean() - lo.mean()), 2) if len(up) and len(lo) else None,
                     "lipids_per_leaflet": [int((dh > 0).sum()), int((dh < 0).sum())]})
    best = max(axes, key=lambda a: a["score"])
    formed = best["tail_ends_in_core"] > 0.8 and best["heads_in_leaflets"] > 0.8 and best["in_plane_coverage"] > 0.9
    result = {"bilayer_formed": formed, "normal_axis": best["axis"], "axes": axes,
              "thickness_PP_nm": best["thickness_PP_nm"], "lipids_per_leaflet": best["lipids_per_leaflet"]}
    if not formed:
        result["advice"] = ("no complete bilayer yet: extend the run (gmx convert-tpr -extend), or check the "
                            "composition/temperature; the bilayer may also have formed as a vesicle or stack")
    return result | {"_protein_bb": bb, "_box": L.tolist()}


def check_build(build: Path, gro: Optional[Path] = None) -> dict:
    """Analyse a self-assembly build directory and map the result to the input frame."""
    np = _np()
    build = Path(build)
    spec_data = yaml.safe_load((build / "spec.resolved.yaml").read_text())
    state = json.loads((build / ".mdbuddy_state.json").read_text())
    gmx = build / "gromacs"
    gro = Path(gro) if gro else gmx / "step7_self_assembly.gro"
    if not gro.exists():
        raise SelfAssemblyError(f"{gro} not found: run the self-assembly first (gromacs/run.sh) or pass --gro")
    lipid_names = list(dict.fromkeys(list((spec_data["membrane"]["upper"] or spec_data["membrane"]["lipids"]))
                                     + list((spec_data["membrane"]["lower"] or spec_data["membrane"]["lipids"]))))
    protein_moltypes = [m[0] for m in state["protein_molecules"]]
    res = analyse(gro, gmx / "topol.top", lipid_names, protein_moltypes)
    bb, L = res.pop("_protein_bb"), np.array(res.pop("_box"))
    if not res["bilayer_formed"] or not bb:
        return res

    k = "xyz".index(res["normal_axis"])
    c = next(a["centre_nm"] for a in res["axes"] if a["axis"] == res["normal_axis"])
    _, fin_atoms, _ = read_gro(gro)
    _, ini_atoms, ini_box = read_gro(gmx / "step5_input.gro")
    fin = _unwrap_chain(np.array([[fin_atoms[i].x, fin_atoms[i].y, fin_atoms[i].z] for i in bb]), L)
    ini = _unwrap_chain(np.array([[ini_atoms[i].x, ini_atoms[i].y, ini_atoms[i].z] for i in bb]),
                        np.array(ini_box[:3]))
    R, t = _kabsch(fin, ini)  # final -> initial box frame
    rmsd = float(np.sqrt((((fin @ R.T + t) - ini) ** 2).sum(axis=1).mean()))

    e = np.zeros(3)
    e[k] = 1.0
    cent = fin.mean(axis=0)
    plane_k = c + L[k] * round((cent[k] - c) / L[k])
    point = cent.copy()
    point[k] = plane_k
    depth = (fin - point) @ e  # BB distance from the midplane (nm)

    # positive-inside: K/R backbone beads just outside the core on each side
    cg_atoms, _ = read_pdb(build / "02_cg" / "cg.pdb")
    bb_res = [a.resname for a in cg_atoms if a.name == "BB"]
    if len(bb_res) == len(bb):
        band = [(d, r) for d, r in zip(depth, bb_res) if 1.5 < abs(d) < 3.0 and r in ("LYS", "ARG")]
        kr_pos, kr_neg = sum(1 for d, _ in band if d > 0), sum(1 for d, _ in band if d < 0)
    else:
        kr_pos = kr_neg = 0
    sign = -1.0 if kr_pos > kr_neg else 1.0  # normal points away from the K/R-rich (cytoplasmic) side

    n_ini = R @ (sign * e)
    p_ini = R @ point + t
    cg_bb = np.array([[a.x, a.y, a.z] for a in cg_atoms if a.name == "BB"]) / 10.0
    t0 = ini.mean(axis=0) - cg_bb.mean(axis=0)  # cg frame (nm) -> initial box frame
    p_or = (p_ini - t0) * 10.0  # oriented (01_orient) frame, Angstrom
    rec = state.get("orientation") or {}
    Ro = np.array(rec.get("rotation", np.eye(3)))
    origin = np.array(rec.get("origin_A", [0.0, 0.0, 0.0])) if "rotation" in rec else np.zeros(3)
    shift = float(rec.get("z_shift_A", 0.0))
    n_in = Ro.T @ n_ini
    p_in = Ro.T @ (p_or - np.array([0.0, 0.0, shift])) + origin

    tm = int((np.abs(depth) < 1.5).sum())
    spans = bool(depth.min() < -1.5 and depth.max() > 1.5)
    res.update({
        "protein_backbone_rmsd_to_start_nm": round(rmsd, 3),
        "protein_backbone_beads_in_core": tm,
        "protein_spans_bilayer": spans,
        "kr_beads_per_side": {"+normal": kr_pos if sign > 0 else kr_neg, "-normal": kr_neg if sign > 0 else kr_pos},
        "membrane_normal_input_frame": [round(float(v), 4) for v in n_in],
        "membrane_centre_input_frame_A": [round(float(v), 3) for v in p_in],
        "rebuild_spec": {"protein": {"orientation": "manual", "manual": {
            "normal": [round(float(v), 4) for v in n_in], "center_point": [round(float(v), 3) for v in p_in]}}},
    })
    if not spans:
        res["advice"] = ("the protein does not span the bilayer: it may be peripheral, or the bilayer formed "
                         "next to it; inspect the trajectory before using this orientation")
    return res
