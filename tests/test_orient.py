"""Membrane placement without OPM data: hydrophobic guess, manual placement, viewer."""

import json
import math
import re

import pytest

np = pytest.importorskip("numpy")

from mdbuddy import orient, viewer  # noqa: E402
from mdbuddy.structure import Atom  # noqa: E402


def toy_bundle(rotation, shift):
    """Four straight helices along +z: Leu belt |z| < 14 A, Glu above (+z), Lys/Arg below (-z)."""
    atoms, resid = [], 0
    for hx, hy in ((7, 7), (-7, 7), (-7, -7), (7, -7)):
        for i in range(36):
            resid += 1
            z = 1.5 * i - 26.0
            name = "LEU" if abs(z) < 14 else ("GLU" if z > 0 else ("LYS" if i % 2 else "ARG"))
            ang = math.radians(100 * i)
            ca = np.array([hx + 2.3 * math.cos(ang), hy + 2.3 * math.sin(ang), z])
            out = np.array([math.cos(ang), math.sin(ang), 0.0])
            pts = {"N": ca + [0, 0, -0.9], "CA": ca, "C": ca + [0, 0, 0.9], "O": ca + [0.5, 0, 1.2],
                   "CB": ca + 1.5 * out, "CG": ca + 2.8 * out, "CD": ca + 4.0 * out}
            for nm, p in pts.items():
                q = rotation @ p + shift
                atoms.append(Atom("ATOM", 0, nm, "", name, "A", resid, "", *[float(v) for v in q],
                                  element=nm[0]))
    return atoms


def test_hydrophobic_recovers_known_placement():
    R = orient.euler_xyz(35, -50, 70)
    t = np.array([12.0, -4.0, 30.0])
    placed, rec = orient.hydrophobic(toy_bundle(R, t))
    p = rec["placement"]
    n = np.array(p["normal"])
    true_n = R @ [0, 0, 1]
    assert math.degrees(math.acos(min(1, abs(n @ true_n)))) < 5
    assert abs((np.array(p["center_A"]) - t) @ true_n) < 2.0
    assert n @ true_n > 0  # Lys/Arg side (-z of the toy) is cytoplasmic -> normal points to the Glu side
    assert p["verdict"] == "transmembrane placement found"
    assert 20 <= p["thickness_A"] <= 36
    # placed frame: Leu belt around z = 0, Lys below, Glu above
    leu = [a.z for a in placed if a.resname == "LEU" and a.name == "CA"]
    lys = [a.z for a in placed if a.resname == "LYS" and a.name == "CA"]
    assert abs(sum(leu) / len(leu)) < 2 and max(lys) < 0


def test_hydrophobic_topology_flip_and_input_frame_roundtrip():
    atoms = toy_bundle(np.eye(3), np.zeros(3))
    _, rec = orient.hydrophobic(atoms, topology="flip")
    assert rec["placement"]["normal"][2] < -0.9
    p1, rec = orient.hydrophobic(atoms)
    n, c = orient.membrane_in_input_frame(rec)
    p2, _ = orient.manual(atoms, {"normal": n, "center_point": c})
    assert max(abs(a.z - b.z) for a, b in zip(p1, p2)) < 1e-6


def test_no_belt_is_not_transmembrane():
    R = np.eye(3)
    atoms = [a for a in toy_bundle(R, np.zeros(3))]
    for a in atoms:
        if a.resname == "LEU":
            a.resname = "GLU"
    _, rec = orient.hydrophobic(atoms, pore_correction="off")
    assert rec["placement"]["verdict"].startswith("not a")


def test_manual_align_rotate_center_shift():
    atoms = toy_bundle(orient.euler_xyz(90, 0, 0), np.zeros(3))  # helices now along -y
    placed, rec = orient.manual(atoms, {"align": ["1-2", "35-36"], "center": "10-25", "shift_z": 2.0})
    ca = [a for a in placed if a.name == "CA" and a.resid <= 36]
    assert ca[-1].z - ca[0].z > 45  # helix 1 now runs along +z
    mid = [a.z for a in placed if 10 <= a.resid <= 25]
    assert abs(sum(mid) / len(mid) - 2.0) < 1e-6
    rotated, _ = orient.manual(atoms, {"rotate": [-90, 0, 0]})  # undo the +90: x' = x0 - centroid(x0)
    x0 = np.array([[a.x, a.y, a.z] for a in toy_bundle(np.eye(3), np.zeros(3))])
    assert np.allclose(np.array([[a.x, a.y, a.z] for a in rotated]), x0 - x0.mean(axis=0), atol=1e-6)


def test_selection_parsing():
    atoms = toy_bundle(np.eye(3), np.zeros(3))
    assert {a.resid for a in orient.select(atoms, "A:1-3,5")} == {1, 2, 3, 5}
    assert {a.resid for a in orient.select(atoms, "7 9-10")} == {7, 9, 10}
    with pytest.raises(orient.OrientError):
        orient.select(atoms, "A:x-y")
    with pytest.raises(orient.OrientError):
        orient.select(atoms, "B:1-5")


def test_viewer_embeds_structure_and_placement(tmp_path):
    atoms = toy_bundle(np.eye(3), np.zeros(3))
    _, rec = orient.hydrophobic(atoms)
    out = viewer.write_html(tmp_path / "v.html", atoms, rec, thickness_A=28, title="toy")
    html = out.read_text()
    assert viewer.THREEDMOL_URL in html and "<title>Membrane Placement Viewer</title>" in html
    data = json.loads(re.search(r"window\.MDBUDDY = (\{.*?\});</script>", html).group(1))
    assert len(data["lines"]) == len(atoms) == len(data["xyz"])
    assert np.allclose(data["rotation"], rec["rotation"]) and data["thickness"] == 28
