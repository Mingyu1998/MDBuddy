"""Place a protein in a membrane when no OPM orientation exists.

Two ways, both producing an "OPM-like" frame (membrane normal = +z, membrane centre at z = 0):

* ``manual``: the user gives rotations / an axis to align / a z shift (CHARMM-GUI "Orient" style),
  optionally chosen interactively with the HTML viewer (``mdbuddy orient view``).
* ``hydrophobic``: an educated guess. Every membrane normal, depth and hydrophobic thickness is
  scored with a simple implicit-membrane model: the water->octanol transfer free energy of each
  residue (Wimley-White scale) weighted by its relative solvent accessibility, switched on inside
  the hydrophobic slab. This is the idea behind PPM/OPM, but much simpler, so treat the result as a
  starting point and check it (``mdbuddy orient view``) - especially for multi-domain proteins.

The "inside" (cytoplasmic) side is chosen with the positive-inside rule (more Lys/Arg next to the
cytoplasmic interface) and put at -z, so the spec's ``lower`` leaflet is the cytoplasmic one.

numpy is required (installed with COBY/vermouth); scipy is used when available.
"""

from __future__ import annotations

import math
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Optional

from .structure import Atom

# Wimley & White water->octanol transfer free energies, kcal/mol (Wimley, Creamer & White,
# Biochemistry 35, 5109 (1996)); charged forms for D, E, K, R; neutral His (pH 7).
WW_OCTANOL = {
    "ALA": 0.50, "ARG": 1.81, "ASN": 0.85, "ASP": 3.64, "CYS": -0.02, "GLN": 0.77, "GLU": 3.63,
    "GLY": 1.15, "HIS": 0.11, "ILE": -1.12, "LEU": -1.25, "LYS": 2.80, "MET": -0.67, "PHE": -1.71,
    "PRO": 0.14, "SER": 0.46, "THR": 0.25, "TRP": -2.09, "TYR": -0.71, "VAL": -0.46,
}
# Theoretical maximum residue accessibility, A^2 (Tien et al., PLoS ONE 8, e80635 (2013))
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104,
    "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
    "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
}
ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
           "CYX": "CYS", "MSE": "MET"}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
PROBE = 1.4
INTERFACE_WIDTH = 1.5  # A, softness of the slab boundary
POSITIVE_BAND = 15.0  # A beyond the hydrophobic boundary counted for the positive-inside rule


class OrientError(RuntimeError):
    pass


def _np():
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - numpy ships with COBY/vermouth
        raise OrientError("membrane orientation needs numpy (pip install numpy)") from exc
    return np


# --------------------------------------------------------------------------- residues ----------

@dataclass
class Residue:
    chain: str
    resid: int
    icode: str
    resname: str
    atoms: list[int] = field(default_factory=list)  # indices into the atom list


def residues_of(atoms: list[Atom]) -> list[Residue]:
    out: "OrderedDict[tuple, Residue]" = OrderedDict()
    for i, a in enumerate(atoms):
        if a.record != "ATOM":
            continue
        name = ALIASES.get(a.resname, a.resname)
        if name not in WW_OCTANOL:
            continue
        key = (a.chain, a.resid, a.icode)
        if key not in out:
            out[key] = Residue(a.chain, a.resid, a.icode, name)
        out[key].atoms.append(i)
    return list(out.values())


def _element(a: Atom) -> str:
    if a.element:
        return a.element.upper()[:1]
    return a.name.lstrip("0123456789")[:1].upper()


def relative_sasa(atoms: list[Atom], residues: list[Residue], n_points: int = 96):
    """Shrake-Rupley relative accessibility per residue (0 buried .. 1 fully exposed)."""
    np = _np()
    xyz = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
    radii = np.array([VDW.get(_element(a), 1.8) + PROBE for a in atoms])
    # golden-spiral unit sphere
    k = np.arange(n_points) + 0.5
    phi = np.arccos(1 - 2 * k / n_points)
    theta = math.pi * (1 + 5 ** 0.5) * k
    sphere = np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(xyz)
        neighbours = tree.query_ball_point(xyz, r=radii + radii.max())
    except ImportError:
        neighbours = None
    atom_area = np.zeros(len(atoms))
    for i in range(len(atoms)):
        if neighbours is not None:
            nb = np.array([j for j in neighbours[i] if j != i], dtype=int)
        else:
            d = np.linalg.norm(xyz - xyz[i], axis=1)
            nb = np.where((d < radii + radii[i]) & (d > 0))[0]
        pts = xyz[i] + radii[i] * sphere
        if len(nb):
            d2 = ((pts[:, None, :] - xyz[nb][None, :, :]) ** 2).sum(-1)
            free = (d2 >= radii[nb][None, :] ** 2).all(axis=1)
        else:
            free = np.ones(n_points, dtype=bool)
        atom_area[i] = 4 * math.pi * radii[i] ** 2 * free.mean()
    rsa = np.array([min(1.0, atom_area[r.atoms].sum() / MAX_ASA[r.resname]) for r in residues])
    return rsa


def _side_chain_centres(atoms: list[Atom], residues: list[Residue]):
    np = _np()
    centres = []
    for r in residues:
        sc = [atoms[i] for i in r.atoms if atoms[i].name not in ("N", "C", "O", "OXT")]
        if r.resname == "GLY" or not sc:
            sc = [atoms[i] for i in r.atoms]
        centres.append([sum(a.x for a in sc) / len(sc), sum(a.y for a in sc) / len(sc),
                        sum(a.z for a in sc) / len(sc)])
    return np.array(centres)


# ------------------------------------------------------------------------- transforms ----------

def rotation_to_z(normal):
    """Rotation matrix taking unit vector ``normal`` onto +z (Rodrigues)."""
    np = _np()
    n = np.asarray(normal, float) / np.linalg.norm(normal)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(n, z)
    c = float(np.dot(n, z))
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (np.linalg.norm(v) ** 2))


def euler_xyz(rx: float, ry: float, rz: float):
    """R = Rz . Ry . Rx (degrees): rotate about x, then y, then z (fixed axes)."""
    np = _np()
    a, b, g = (math.radians(v) for v in (rx, ry, rz))
    Rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    Ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    Rz = np.array([[math.cos(g), -math.sin(g), 0], [math.sin(g), math.cos(g), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def apply(atoms: list[Atom], R, origin) -> list[Atom]:
    """New atoms with x' = R (x - origin)."""
    np = _np()
    out = []
    for a in atoms:
        x = R @ (np.array([a.x, a.y, a.z]) - origin)
        b = Atom(**{**a.__dict__})
        b.x, b.y, b.z = (float(v) for v in x)
        out.append(b)
    return out


_SEL = re.compile(r"^(?:(?P<chain>[A-Za-z0-9]):)?(?P<ranges>[0-9,\-]+)$")


def select(atoms: list[Atom], selection: str) -> list[Atom]:
    """'A:10-30,45' or '10-30' (all chains); several selections separated by spaces."""
    picked = []
    for part in selection.split():
        m = _SEL.match(part)
        if not m:
            raise OrientError(f"bad selection '{part}': use e.g. 'A:10-30,45' or '10-30'")
        ids: set[int] = set()
        for rng in m.group("ranges").split(","):
            lo, _, hi = rng.partition("-")
            ids.update(range(int(lo), int(hi or lo) + 1))
        picked += [a for a in atoms if a.resid in ids and (m.group("chain") in (None, a.chain))]
    if not picked:
        raise OrientError(f"selection '{selection}' matches no atoms")
    return picked


def _centroid(atoms: list[Atom]):
    np = _np()
    return np.array([[a.x, a.y, a.z] for a in atoms]).mean(axis=0)


def manual(atoms: list[Atom], settings: dict) -> tuple[list[Atom], dict]:
    """CHARMM-GUI-style manual placement. Applied in this order (all optional):

      normal: [x, y, z]        membrane normal in the input frame (towards the non-cytoplasmic side)
      align: [selA, selB]      or: make the vector centroid(selA) -> centroid(selB) the normal
      center_point: [x, y, z]  membrane centre in the input frame, Angstrom (default: protein centre)
      rotate: [rx, ry, rz]     extra rotation in degrees (x, then y, then z) about the membrane centre
      center: sel              put the centroid of this selection at the membrane midplane (z = 0)
      shift_z: dz              finally move the protein by dz Angstrom along the normal

    The result is in the OPM-like frame: normal = +z, membrane centre at z = 0.
    """
    np = _np()
    settings = {k: v for k, v in settings.items() if v is not None}
    origin = np.array(settings["center_point"], float) if settings.get("center_point") else _centroid(atoms)
    R = np.eye(3)
    if settings.get("normal"):
        vec = np.array(settings["normal"], float)
        if np.linalg.norm(vec) < 1e-9:
            raise OrientError("manual.normal must not be the zero vector")
        R = rotation_to_z(vec)
    elif settings.get("align"):
        sel_a, sel_b = settings["align"]
        vec = _centroid(select(atoms, sel_b)) - _centroid(select(atoms, sel_a))
        if np.linalg.norm(vec) < 1e-6:
            raise OrientError("manual.align: the two selections have the same centroid")
        R = rotation_to_z(vec)
    if settings.get("rotate"):
        R = euler_xyz(*settings["rotate"]) @ R
    placed = apply(atoms, R, origin)
    dz = 0.0
    if settings.get("center"):
        dz -= _centroid(select(placed, settings["center"]))[2]
    dz += float(settings.get("shift_z") or 0.0)
    for a in placed:
        a.z += dz
    record = {"method": "manual", "settings": settings, "rotation": R.tolist(),
              "origin_A": [float(v) for v in origin], "z_shift_A": dz}
    return placed, record


def membrane_in_input_frame(record: dict) -> tuple[list[float], list[float]]:
    """(normal, centre point) of the membrane expressed in the input frame of a placement record."""
    np = _np()
    R = np.array(record["rotation"])
    origin = np.array(record["origin_A"])
    shift = float(record.get("z_shift_A", 0.0))
    normal = R.T @ np.array([0.0, 0.0, 1.0])
    centre = origin + R.T @ np.array([0.0, 0.0, -shift])
    return [float(v) for v in normal], [float(v) for v in centre]


# ----------------------------------------------------------------------- hydrophobic ----------

def _sphere_points(n: int, hemisphere: bool = True):
    np = _np()
    k = np.arange(n) + 0.5
    zc = 1 - k / n if hemisphere else 1 - 2 * k / n
    r = np.sqrt(1 - zc ** 2)
    theta = math.pi * (1 + 5 ** 0.5) * k
    return np.stack([r * np.cos(theta), r * np.sin(theta), zc], axis=1)


def _cone(axis, max_deg: float, n: int):
    np = _np()
    R = rotation_to_z(axis).T  # maps +z onto axis
    cosmax = math.cos(math.radians(max_deg))
    k = np.arange(n) + 0.5
    zc = 1 - (1 - cosmax) * k / n
    r = np.sqrt(1 - zc ** 2)
    theta = math.pi * (1 + 5 ** 0.5) * k
    pts = np.stack([r * np.cos(theta), r * np.sin(theta), zc], axis=1)
    return np.vstack([axis / np.linalg.norm(axis), pts @ R.T])


def _scores(proj, weights, shifts, thicknesses):
    """E[d, D] = sum_i w_i * slab(proj_i - d, D)."""
    np = _np()
    z = np.abs(proj[None, None, :] - shifts[:, None, None])
    half = thicknesses[None, :, None] / 2
    slab = 1.0 / (1.0 + np.exp((z - half) / INTERFACE_WIDTH))
    return slab @ weights


@dataclass
class Placement:
    normal: list[float]  # input frame, unit vector (towards the non-cytoplasmic side)
    center_A: list[float]  # membrane centre point in the input frame
    thickness_A: float  # hydrophobic thickness
    energy_kcal: float  # implicit-membrane transfer energy of the best placement
    uncertainty_deg: float  # max angle between near-optimal normals (within 1 kcal/mol)
    core_residues: int  # exposed residues inside the hydrophobic core
    tm_segments: list[str]
    kr_inside: int
    kr_outside: int
    flipped_by_positive_inside: bool
    verdict: str
    notes: list[str]

    def as_dict(self) -> dict:
        return dict(self.__dict__)


HULL_DEPTH = 10.0  # A: residues deeper than this inside the convex hull cannot face lipids


def hull_depth(points):
    """Distance of each point inside the convex hull to the hull surface (A); None without scipy."""
    np = _np()
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return None
    hull = ConvexHull(points)
    # facet equations: n.x + b <= 0 inside, n unit length
    return np.min(-(points @ hull.equations[:, :3].T + hull.equations[:, 3]), axis=1)


def hydrophobic(atoms: list[Atom], topology: str = "positive_inside",
                thickness_range: tuple[float, float] = (20.0, 40.0),
                pore_correction: str = "auto") -> tuple[list[Atom], dict]:
    """Educated guess of the membrane placement from exposed hydrophobicity.

    pore_correction: "auto" (retry with lipid-facing weighting when the plain scan finds no
    clear belt, e.g. beta-barrels whose charged pores cancel it), "on" or "off".
    """
    if pore_correction not in ("auto", "on", "off"):
        raise OrientError("pore_correction must be auto, on or off")
    if pore_correction == "auto":
        placed, record = _hydrophobic(atoms, topology, thickness_range, False)
        if not record["placement"]["verdict"].startswith("not a"):
            return placed, record
        placed2, record2 = _hydrophobic(atoms, topology, thickness_range, True)
        if record2["placement"]["verdict"].startswith("not a"):
            return placed, record
        record2["placement"]["notes"].insert(0, "plain scan found no hydrophobic belt; used pore correction "
                                                "(only residues near the outer surface count), as for beta-barrels")
        return placed2, record2
    return _hydrophobic(atoms, topology, thickness_range, pore_correction == "on")


def _hydrophobic(atoms, topology, thickness_range, pore: bool):
    np = _np()
    residues = residues_of(atoms)
    if len(residues) < 10:
        raise OrientError("too few standard residues to estimate a membrane placement")
    rsa = relative_sasa(atoms, residues)
    dg = np.array([WW_OCTANOL[r.resname] for r in residues])
    pos = _side_chain_centres(atoms, residues)
    if pore:
        depth = hull_depth(pos)
        if depth is None:
            raise OrientError("pore correction needs scipy (pip install scipy)")
        rsa = rsa / (1.0 + np.exp((depth - HULL_DEPTH) / 2.0))
    weights = rsa * dg
    centre = pos.mean(axis=0)
    rel = pos - centre
    extent = float(np.abs(rel).max())

    # coarse: 1000 directions, fixed thickness
    dirs = _sphere_points(1000)
    shifts = np.arange(-min(extent, 40.0), min(extent, 40.0) + 1e-9, 2.0)
    coarse = np.array([[30.0]])
    best = []
    for n in dirs:
        e = _scores(rel @ n, weights, shifts, coarse[0])
        i = int(np.argmin(e[:, 0]))
        best.append((float(e[i, 0]), n, float(shifts[i])))
    best.sort(key=lambda t: t[0])

    # refine around the 6 best distinct directions: finer normals, shifts and thickness
    thick = np.arange(thickness_range[0], thickness_range[1] + 1e-9, 1.0)
    seeds, solutions = [], []
    for e0, n0, d0 in best:
        if all(abs(float(np.dot(n0, s[0]))) < math.cos(math.radians(15)) for s in seeds):
            seeds.append((n0, d0))
        if len(seeds) == 6:
            break
    for s, d0 in seeds:
        fine_shifts = d0 + np.arange(-6.0, 6.0 + 1e-9, 0.5)
        for n in _cone(s, 10.0, 80):
            e = _scores(rel @ n, weights, fine_shifts, thick)
            i, j = np.unravel_index(int(np.argmin(e)), e.shape)
            solutions.append((float(e[i, j]), n / np.linalg.norm(n), float(fine_shifts[i]), float(thick[j])))
    solutions.sort(key=lambda t: t[0])
    e_best, n_best, d_best, D_best = solutions[0]
    near = [s for s in solutions if s[0] <= e_best + 1.0]
    spread = max(math.degrees(math.acos(min(1.0, abs(float(np.dot(n_best, s[1])))))) for s in near)

    m_center = centre + d_best * n_best
    R = rotation_to_z(n_best)
    z = (pos - m_center) @ n_best
    half = D_best / 2
    up = [r for r, zz in zip(residues, z) if half < zz < half + POSITIVE_BAND and r.resname in ("LYS", "ARG")]
    lo = [r for r, zz in zip(residues, z) if -half - POSITIVE_BAND < zz < -half and r.resname in ("LYS", "ARG")]
    flipped = False
    notes = []
    if topology == "positive_inside":
        if len(up) > len(lo):  # cytoplasm (K/R rich) must be at -z
            R = np.diag([1.0, -1.0, -1.0]) @ R
            n_best = -n_best
            up, lo = lo, up
            flipped = True
        if abs(len(up) - len(lo)) <= 2:
            notes.append(f"positive-inside rule is weak here ({len(lo)} vs {len(up)} Lys/Arg); "
                         "the cytoplasmic side may be the other one")
    elif topology == "flip":
        R = np.diag([1.0, -1.0, -1.0]) @ R
        n_best = -n_best
        up, lo = lo, up
        flipped = True

    placed = apply(atoms, R, m_center)
    # residues in the core and TM segments (by CA)
    zc = (pos - m_center) @ n_best
    core = int(((np.abs(zc) < half) & (rsa > 0.15)).sum())
    segments = _tm_segments(residues, atoms, R, m_center, half)

    if e_best > -5 or not segments:
        verdict = "not a transmembrane protein (or no clear hydrophobic belt): consider peripheral/soluble"
    elif spread > 20:
        verdict = "membrane protein, but the orientation is poorly defined - check it visually"
    else:
        verdict = "transmembrane placement found"
    if extent > 2.5 * D_best:
        notes.append("large protein relative to the membrane: soluble domains can bias this simple model")

    result = Placement(
        normal=[float(v) for v in n_best], center_A=[float(v) for v in m_center], thickness_A=float(D_best),
        energy_kcal=round(e_best, 2), uncertainty_deg=round(spread, 1), core_residues=core,
        tm_segments=segments, kr_inside=len(lo), kr_outside=len(up), flipped_by_positive_inside=flipped,
        verdict=verdict, notes=notes)
    record = {"method": "hydrophobic", "pore_correction": pore, "rotation": R.tolist(), "origin_A": m_center.tolist(),
              "placement": result.as_dict()}
    return placed, record


def _tm_segments(residues, atoms, R, origin, half, min_len: int = 8) -> list[str]:
    np = _np()
    segs, cur = [], []
    prev = None
    for r in residues:
        ca = [atoms[i] for i in r.atoms if atoms[i].name == "CA"]
        if not ca:
            continue
        z = float((R @ (np.array([ca[0].x, ca[0].y, ca[0].z]) - origin))[2])
        inside = abs(z) < half
        contiguous = prev is not None and prev.chain == r.chain and r.resid == prev.resid + 1
        if inside and (contiguous or not cur):
            cur.append(r)
        else:
            if len(cur) >= min_len:
                segs.append(f"{cur[0].chain}:{cur[0].resid}-{cur[-1].resid}")
            cur = [r] if inside else []
        prev = r
    if len(cur) >= min_len:
        segs.append(f"{cur[0].chain}:{cur[0].resid}-{cur[-1].resid}")
    return segs
