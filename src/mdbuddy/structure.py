"""Minimal, dependency-free PDB/GRO handling.

Only what the pipeline needs: read ATOM/HETATM records, filter them, write
them back, and compute geometry (centre of geometry, extents).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

WATER_RESNAMES = {"HOH", "WAT", "TIP3", "SOL", "DOD", "H2O"}
DUMMY_RESNAMES = {"DUM"}  # OPM membrane boundary pseudo-atoms


@dataclass
class Atom:
    record: str  # ATOM or HETATM
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resid: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""

    @classmethod
    def from_pdb_line(cls, line: str) -> "Atom":
        def num(s: str, default: float = 0.0) -> float:
            s = s.strip()
            return float(s) if s else default

        serial_txt = line[6:11].strip()
        return cls(
            record=line[0:6].strip(),
            serial=int(serial_txt) if serial_txt.isdigit() else 0,
            name=line[12:16].strip(),
            altloc=line[16:17].strip(),
            resname=line[17:21].strip(),
            chain=line[21:22].strip(),
            resid=int(line[22:26]),
            icode=line[26:27].strip(),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            occupancy=num(line[54:60], 1.0),
            bfactor=num(line[60:66], 0.0),
            element=line[76:78].strip() if len(line) >= 78 else "",
        )

    def to_pdb_line(self, serial: Optional[int] = None) -> str:
        serial = self.serial if serial is None else serial
        name = self.name
        # PDB convention: names shorter than 4 chars start in column 14
        name_field = f" {name:<3}" if len(name) < 4 else name[:4]
        return (
            f"{self.record:<6}{serial % 100000:>5} {name_field}{self.altloc or ' '}"
            f"{self.resname:>3} {self.chain or ' '}{self.resid % 10000:>4}{self.icode or ' '}   "
            f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}{self.occupancy:6.2f}{self.bfactor:6.2f}"
            f"          {self.element:>2}"
        )

    @property
    def is_hydrogen(self) -> bool:
        if self.element:
            return self.element.upper() in ("H", "D")
        return self.name.lstrip("0123456789").upper().startswith("H")


def read_pdb(path: str | Path) -> tuple[list[Atom], list[str]]:
    """Return (atoms, remark/header lines). Only the first MODEL is read."""
    atoms: list[Atom] = []
    header: list[str] = []
    with open(path) as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec in ("ATOM", "HETATM"):
                atoms.append(Atom.from_pdb_line(line.rstrip("\n")))
            elif rec == "ENDMDL":
                break
            elif rec in ("REMARK", "HEADER", "TITLE", "COMPND", "CRYST1"):
                header.append(line.rstrip("\n"))
    return atoms, header


def write_pdb(path: str | Path, atoms: Iterable[Atom], header: Iterable[str] = ()) -> None:
    lines = list(header)
    prev_chain = None
    for i, atom in enumerate(atoms, start=1):
        if prev_chain is not None and atom.chain != prev_chain and atom.record == "ATOM":
            lines.append("TER")
        prev_chain = atom.chain
        lines.append(atom.to_pdb_line(serial=i))
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def centre_of_geometry(atoms: Iterable[Atom]) -> tuple[float, float, float]:
    xs, ys, zs, n = 0.0, 0.0, 0.0, 0
    for a in atoms:
        xs += a.x
        ys += a.y
        zs += a.z
        n += 1
    if n == 0:
        raise ValueError("no atoms")
    return xs / n, ys / n, zs / n


def extents(atoms: Iterable[Atom]) -> tuple[float, float, float]:
    atoms = list(atoms)
    if not atoms:
        raise ValueError("no atoms")
    return (
        max(a.x for a in atoms) - min(a.x for a in atoms),
        max(a.y for a in atoms) - min(a.y for a in atoms),
        max(a.z for a in atoms) - min(a.z for a in atoms),
    )


def opm_half_thickness(header: Iterable[str]) -> Optional[float]:
    """Parse 'REMARK 1/2 of bilayer thickness: 16.2' written by OPM (Angstrom)."""
    for line in header:
        if "1/2 of bilayer thickness" in line:
            try:
                return float(line.split(":")[1].split()[0])
            except (IndexError, ValueError):
                return None
    return None


@dataclass
class GroAtom:
    resid: int
    resname: str
    name: str
    index: int
    x: float
    y: float
    z: float


def read_gro(path: str | Path) -> tuple[str, list[GroAtom], list[float]]:
    """Return (title, atoms, box). Coordinates in nm."""
    with open(path) as fh:
        title = fh.readline().rstrip("\n")
        natoms = int(fh.readline().strip())
        atoms = []
        for _ in range(natoms):
            line = fh.readline()
            atoms.append(GroAtom(
                resid=int(line[0:5]),
                resname=line[5:10].strip(),
                name=line[10:15].strip(),
                index=int(line[15:20]),
                x=float(line[20:28]),
                y=float(line[28:36]),
                z=float(line[36:44]),
            ))
        box = [float(v) for v in fh.readline().split()]
    return title, atoms, box
