"""Shared fixtures. The pipeline test uses tiny stand-ins for martinize2, COBY and the
force-field downloads so the plumbing can be tested without the real tools."""

from __future__ import annotations

import hashlib
import os
import stat
import textwrap
from pathlib import Path

import pytest

from mdbuddy import forcefields
from mdbuddy.structure import Atom, write_pdb


def make_atom(record, name, resname, chain, resid, x, y, z, altloc="", element=""):
    return Atom(record=record, serial=0, name=name, altloc=altloc, resname=resname, chain=chain,
                resid=resid, icode="", x=x, y=y, z=z, element=element)


@pytest.fixture
def opm_pdb(tmp_path: Path) -> Path:
    """Two-chain toy 'membrane protein' in OPM style (membrane centre at z=0), complete heavy atoms."""
    from mdbuddy.prepare import HEAVY_ATOMS

    def residue(chain, resid, resname, x0, z, skip=()):
        out = []
        for k, name in enumerate(("N", "CA", "C", "O") + HEAVY_ATOMS[resname]):
            if name not in skip:
                out.append(make_atom("ATOM", name, resname, chain, resid, x0 + 0.3 * k, 0.0, z,
                                     element=name[0]))
        return out

    atoms = []
    for chain, x0 in (("A", 0.0), ("B", 10.0)):
        for i, resname in enumerate(("ALA", "LEU", "LYS"), start=1):
            z = -5.0 + 3.8 * (i - 1)
            atoms += residue(chain, i, resname, x0, z, skip=("CB",) if resname == "LYS" else ())
            atoms.append(make_atom("ATOM", "H", resname, chain, i, x0, 1.0, z, element="H"))
        # altloc pair on the last CB
        atoms.append(make_atom("ATOM", "CB", "LYS", chain, 3, x0 + 1, 0, 3.0, altloc="A", element="C"))
        atoms.append(make_atom("ATOM", "CB", "LYS", chain, 3, x0 + 1.2, 0, 3.0, altloc="B", element="C"))
    # a gap in chain B: residue 10 far away
    for k, name in enumerate(("N", "CA", "C", "O")):
        atoms.append(make_atom("ATOM", name, "GLY", "B", 10, 30.0 + 0.3 * k, 0.0, 0.0, element=name[0]))
    atoms.append(make_atom("HETATM", "O", "HOH", "W", 1, 5.0, 5.0, 20.0, element="O"))
    atoms.append(make_atom("HETATM", "C1", "RET", "A", 300, 1.0, 1.0, 1.0, element="C"))
    atoms.append(make_atom("HETATM", "O", "DUM", "", 400, 0.0, 0.0, 16.2, element="O"))
    path = tmp_path / "toy_opm.pdb"
    write_pdb(path, atoms, header=["REMARK      1/2 of bilayer thickness:   16.2"])
    return path


FAKE_MARTINIZE = r'''#!{python}
"""Stand-in for martinize2: one BB bead per residue, one molecule (chains merged)."""
import sys
args = sys.argv[1:]
def opt(flag, default=None):
    return args[args.index(flag) + 1] if flag in args else default
inp, xout, top, name = opt("-f"), opt("-x"), opt("-o"), opt("-name", "molecule")
open("args.txt", "w").write(" ".join(args))
res = []
for line in open(inp):
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        res.append((line[17:20], line[21], int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54])))
with open(xout, "w") as fh:
    for i, (rn, ch, ri, x, y, z) in enumerate(res, 1):
        fh.write("ATOM  %5d  BB  %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n" % (i, rn, ch, ri, x, y, z))
    fh.write("END\n")
mol = name + "_0"
with open(mol + ".itp", "w") as fh:
    fh.write("#ifndef POSRES_FC\n#define POSRES_FC 1000\n#endif\n\n[ moleculetype ]\n%s 1\n\n[ atoms ]\n" % mol)
    for i, (rn, ch, ri, x, y, z) in enumerate(res, 1):
        q = 1.0 if rn == "LYS" else 0.0
        fh.write("%d P5 %d %s BB %d %.1f\n" % (i, ri, rn, i, q))
with open(top, "w") as fh:
    fh.write('#include "martini.itp"\n#include "%s.itp"\n\n[ system ]\nx\n\n[ molecules ]\n%s 1\n' % (mol, mol))
print("INFO done")
'''

FAKE_COBY = r'''
"""Stand-in for COBY: protein + POPC (membrane or flooded) + water + ions, neutral."""
def _write(prefix, atoms, box, mols):
    with open(prefix + ".gro", "w") as fh:
        fh.write("fake\n%d\n" % len(atoms))
        for i, (ri, rn, an) in enumerate(atoms, 1):
            fh.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (ri, rn, an, i, 1.0, 1.0, 1.0))
        fh.write("%10.5f%10.5f%10.5f\n" % tuple(box))


def COBY(**kw):
    import json
    if "protein" not in kw:  # lipid template build (self-assembly)
        _write(kw["out_sys"], [(1, "POPC", "NC3"), (1, "POPC", "PO4"), (2, "POPC", "NC3"), (2, "POPC", "PO4")],
               kw["box"], [])
        return
    json.dump(kw, open("coby_kwargs.json", "w"))
    prot = kw["protein"].split()[0][len("file:"):]
    beads = [l for l in open(prot) if l.startswith("ATOM")]
    moltypes = [t.split(":")[1] for t in kw["protein"].split() if t.startswith("moleculetypes:")]
    nlys = sum(1 for l in beads if l[17:20] == "LYS")
    atoms = [(1, l[17:20], "BB") for l in beads]
    r = 1
    mols = [(m, 1) for m in moltypes]
    nlip = 4
    if "flooding" in kw:
        nlip = sum(int(t.split(":")[2]) for t in kw["flooding"].split() if t.startswith("solute:"))
    for _ in range(nlip):
        r += 1
        atoms += [(r, "POPC", "NC3"), (r, "POPC", "PO4")]
    mols.append(("POPC", nlip))
    for _ in range(10):
        r += 1
        atoms.append((r, "W", "W"))
    mols.append(("W", 10))
    for _ in range(2):
        r += 1
        atoms.append((r, "NA", "NA"))
    for _ in range(2 + nlys):
        r += 1
        atoms.append((r, "CL", "CL"))
    mols += [("NA", 2), ("CL", 2 + nlys)]
    _write(kw["out_sys"], atoms, kw["box"], mols)
    with open(kw["out_top"], "w") as fh:
        fh.write("[ system ]\nfake\n\n[ molecules ]\n")
        for n, c in mols:
            fh.write("%s %d\n" % (n, c))
    open(kw["out_log"], "w").write("fake COBY\n")
'''

FF_FILES = {
    "martini_v3.0.0.itp": """
        [ defaults ]
        1 1

        [ atomtypes ]
        P5 72.0 0.000 A 0.0 0.0

        [ nonbond_params ]
        P5 P5 1 0.0 0.0
    """,
    "martini_v3.0.0_ions_v1.itp": """
        [ moleculetype ]
        NA 1
        [ atoms ]
        1 TQ5 1 NA NA 1 1.0
        [ moleculetype ]
        CL 1
        [ atoms ]
        1 TQ5 1 CL CL 1 -1.0
    """,
    "martini_v3.0.0_solvents_v1.itp": """
        [ moleculetype ]
        W 1
        [ atoms ]
        1 W 1 W W 1 0
    """,
    "martini_v3.0.0_ffbonded_v2.itp": """
        #define b_fake 1 0.47 1250
    """,
    "martini_v3.0.0_phospholipids_PC_v2.itp": """
        [ moleculetype ]
        POPC 1
        [ atoms ]
        1 Q1 1 POPC NC3 1 1.0
        2 Q5 1 POPC PO4 2 -1.0
    """,
}


@pytest.fixture
def fake_tools(tmp_path: Path, monkeypatch):
    """Put fake martinize2 and COBY on the path and serve fake force-field files."""
    import sys

    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "martinize2"
    exe.write_text(FAKE_MARTINIZE.replace("{python}", sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")

    pydir = tmp_path / "py"
    pydir.mkdir()
    (pydir / "COBY.py").write_text(FAKE_COBY)
    monkeypatch.setenv("PYTHONPATH", str(pydir))

    ffdir = tmp_path / "ffsrc"
    ffdir.mkdir()
    files = {}
    for name, content in FF_FILES.items():
        p = ffdir / name
        p.write_text(textwrap.dedent(content).lstrip())
        files[name] = forcefields.FFFile(name, p.as_uri(), hashlib.sha256(p.read_bytes()).hexdigest(), "core")
    lipid_files = ("martini_v3.0.0_ffbonded_v2.itp", "martini_v3.0.0_phospholipids_PC_v2.itp")
    monkeypatch.setitem(forcefields.REGISTRY, "martini3", [f for n, f in files.items() if n not in lipid_files])
    pc = files["martini_v3.0.0_phospholipids_PC_v2.itp"]
    monkeypatch.setitem(forcefields.LIPID_SETS, "lipidome", [
        files["martini_v3.0.0_ffbonded_v2.itp"],
        forcefields.FFFile(pc.name, pc.url, pc.sha256, "lipid"),
    ])
    monkeypatch.setenv("MDBUDDY_CACHE", str(tmp_path / "cache"))
    monkeypatch.delenv("MDBUDDY_GMX", raising=False)
    monkeypatch.delenv("GMX", raising=False)
    return {"bin": bindir, "ffdir": ffdir}
