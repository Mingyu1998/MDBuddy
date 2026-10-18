"""Build the membrane, insert the protein and solvate with COBY.

MDBuddy writes a standalone ``coby_build.py`` and runs it, so the exact COBY
call is kept with the outputs and can be re-run by hand.
"""

from __future__ import annotations

import math
import shutil
import sys
from pathlib import Path
from typing import Optional

from . import cgprotein, lipids
from .provenance import sha256
from .structure import centre_of_geometry, extents, read_pdb


class MembraneError(RuntimeError):
    pass


def _fmt(x: float) -> str:
    return f"{x:.3f}".rstrip("0").rstrip(".")


def box_size(spec, cg_atoms=None, extent_nm: Optional[tuple[float, float, float]] = None) -> list[float]:
    """Box in nm. 'auto' = protein extent + padding (square in x/y); 'from_structure' = first cg_proteins box."""
    m = spec.membrane
    if m["box"] == "from_structure":
        return cgprotein.box_from_structure(spec)
    if m["box"] != "auto":
        return [float(v) for v in m["box"]]
    ex, ey, ez = extent_nm if extent_nm else (e / 10.0 for e in extents(cg_atoms))
    xy = max(ex, ey) + 2 * m["lipid_padding"]
    z = max(ez + 2 * m["water_padding"], m["min_box_z"])
    # round up to 0.5 nm so boxes are tidy and stable across tiny coordinate changes
    return [math.ceil(xy * 2) / 2, math.ceil(xy * 2) / 2, math.ceil(z * 2) / 2]


def protein_offset_nm(spec, cg_atoms, membrane_center_z_A: float = 0.0) -> float:
    """z offset (nm) that puts the membrane midplane of the input frame at the box centre.

    COBY places the protein's centre of geometry at box centre + (cx, cy, cz) and
    the membrane midplane at the box centre. In the input frame the midplane is
    at ``membrane_center_z`` (0 for OPM), so cz = COG_z - membrane_center_z.
    """
    _, _, cz = centre_of_geometry(cg_atoms)
    return (cz - float(membrane_center_z_A)) / 10.0


def _composition(comp: dict, lib: str, lipid_apl: Optional[dict] = None) -> str:
    """COBY lipid sub-arguments, each with the COBY parameter library that matches its itp.

    ``lipid_apl`` gives single lipids their own area per lipid (nm^2), which changes how many of
    them fit in the leaflet (their density) while the ratios stay as given.
    """
    parts = []
    for name, ratio in comp.items():
        arg = f"lipid:{name}:{_fmt(float(ratio))}"
        if lipids.lookup(name, lib):
            arg += f":params:{lipids.coby_params(name, lib)}"
        if lipid_apl and name in lipid_apl:
            arg += f":apl:{_fmt(float(lipid_apl[name]))}"
        parts.append(arg)
    return " ".join(parts)


def membrane_arg(spec, extra: Optional[str] = None) -> str:
    """The generated membrane: membrane-wide options first (COBY applies later ones to a leaflet)."""
    m = spec.membrane
    lib = lipids.library_for(spec)
    leaflets = spec.leaflets()
    upper_apl = m["upper_apl"] or m["apl"]
    lower_apl = m["lower_apl"] or m["apl"]
    head = []
    if float(m["cz"]):
        head.append(f"cz:{_fmt(float(m['cz']))}")
    extra = spec.data["coby"]["membrane_args"] if extra is None else extra
    if extra:
        head.append(extra)
    return " ".join(head + [
        f"leaflet:upper {_composition(leaflets['upper'], lib, m['lipid_apl'])} apl:{_fmt(upper_apl)}",
        f"leaflet:lower {_composition(leaflets['lower'], lib, m['lipid_apl'])} apl:{_fmt(lower_apl)}",
    ])


def protein_arg(spec, cg_pdb_rel: str, molecules: list[tuple[str, int]], cz: float,
                extra: Optional[str] = None) -> str:
    moltypes = " ".join(f"moleculetypes:{name}:{count}" for name, count in molecules)
    extra = spec.data["coby"]["protein_args"] if extra is None else extra
    return f"file:{cg_pdb_rel} {moltypes} cz:{cz:.3f}" + (f" {extra}" if extra else "")


def solvation_arg(spec, extra: Optional[str] = None) -> str:
    s = spec.solvent
    extra = spec.data["coby"]["solvation_args"] if extra is None else extra
    return (f"solv:{s['water']} pos:{s['cation']} neg:{s['anion']} "
            f"salt_molarity:{_fmt(float(s['salt_molarity']))}" + (f" {extra}" if extra else ""))


class _UserFiles:
    """Copies files named in raw COBY strings (file:X, import_library) into 03_membrane/user_files."""

    def __init__(self, spec, work: Path):
        self.spec, self.dir, self.records = spec, work / "user_files", []

    def _copy(self, path: str) -> str:
        src = self.spec.resolve_path(path)
        if not src.exists():
            raise MembraneError(f"file referenced in the coby section not found: {src}")
        self.dir.mkdir(parents=True, exist_ok=True)
        digest = sha256(src)
        for r in self.records:
            if r["sha256"] == digest:
                return r["copied_to"]
        name = src.name
        if (self.dir / name).exists():
            name = f"{len(self.records)}_{src.name}"
        shutil.copyfile(src, self.dir / name)
        rel = f"user_files/{name}"
        self.records.append({"source": str(src), "copied_to": rel, "sha256": digest})
        return rel

    def string(self, text: str) -> str:
        return " ".join(f"file:{self._copy(t[5:])}" if t.startswith("file:") and len(t) > 5 else t
                        for t in str(text).split())

    def value(self, key: str, value):
        if isinstance(value, (list, tuple)):
            return [self.value(key, v) for v in value]
        if isinstance(value, str):
            if key == "import_library" and value.endswith(".py") and not value.startswith("file:"):
                return self._copy(value)
            return self.string(value)
        return value


def fill_params(raw: str, lib: str) -> str:
    """Add the matching COBY parameter library to lipid tokens of a raw membrane string that lack one."""
    out = []
    for tok in raw.split():
        parts = tok.split(":")
        kind = parts[0].lower()
        if kind in ("lipid", "lipid_extra") and len(parts) > 1 and "params" not in parts:
            name = (parts[2 + parts[1:].index("name")] if "name" in parts[1:-1] else
                    parts[1] if kind == "lipid" else None)
            if name and lipids.lookup(name, lib):
                tok += f":params:{lipids.coby_params(name, lib)}"
        out.append(tok)
    return " ".join(out)


def _extent_nm(xyz) -> tuple[float, float, float]:
    return tuple(max(p[k] for p in xyz) - min(p[k] for p in xyz) for k in range(3))


def write_script(spec, work: Path, cg_pdb: Optional[Path], molecules, itp_files: list[Path],
                 membrane_center_z_A: float = 0.0, cg_entries: list[dict] = ()) -> tuple[Path, dict]:
    """Write coby_build.py. ``cg_pdb``/``molecules``: the martinized protein (None if there is none);
    ``cg_entries``: pre-coarse-grained proteins from cgprotein.prepare."""
    work.mkdir(parents=True, exist_ok=True)
    rel = lambda p: Path(_relpath(p, work)).as_posix()
    files = _UserFiles(spec, work)
    proteins, cz = [], None
    if cg_pdb is not None:
        atoms, _ = read_pdb(cg_pdb)
        if not atoms:
            raise MembraneError(f"{cg_pdb} contains no beads")
        box = box_size(spec, atoms)
        cz = protein_offset_nm(spec, atoms, membrane_center_z_A)
        if abs(cz) > box[2] / 2:
            raise MembraneError(
                f"protein would be shifted by {cz:.1f} nm, outside the box; the input does not look "
                "membrane-oriented. Use an OPM structure or set protein.membrane_center_z.")
        proteins.append(protein_arg(spec, rel(cg_pdb), molecules, cz,
                                    files.string(spec.data["coby"]["protein_args"])))
    else:
        box = box_size(spec, extent_nm=_extent_nm(cg_entries[0]["xyz"]) if cg_entries else None)
    notes = []
    for e in cg_entries:
        proteins.append(cgprotein.coby_arg({**e, "coby_args": files.string(e["coby_args"])}, rel(e["structure"])))
        if e["placement"] == "keep" and any(abs(x - y) > 1e-3 for x, y in zip(e["box"][:3], box)):
            notes.append(f"{e['structure'].name}: box {e['box'][:3]} nm differs from the new box; its coordinates "
                         "are kept relative to the box centre")

    coby = spec.data["coby"]
    lib = lipids.library_for(spec)
    membranes = [fill_params(files.string(m), lib) for m in coby["membranes"]] or [membrane_arg(spec, files.string(coby["membrane_args"]))]
    kwargs = {
        "randseed": spec.seed,
        "box": box,
        "box_type": spec.membrane["box_type"],
        "params": spec.ff_info["coby_params"],
        "membrane": membranes[0] if len(membranes) == 1 else membranes,
    }
    if proteins:
        kwargs["protein"] = proteins[0] if len(proteins) == 1 else proteins
    kwargs["solvation"] = solvation_arg(spec, files.string(coby["solvation_args"]))
    for key, value in coby["kwargs"].items():
        kwargs[key] = files.value(key, value)
    if "molecule_import" in kwargs:
        # an imported lipid is looked up under the parameter library of the membrane it goes into
        # (the system's, as MDBuddy writes no params: for lipids it does not know)
        imports = kwargs["molecule_import"]
        imports = [imports] if isinstance(imports, str) else imports
        kwargs["molecule_import"] = [
            m + f" params:{spec.ff_info['coby_params']}"
            if "library_type" in m and "lipid" in m and "params:" not in m else m for m in imports]
    kwargs.update({
        "itp_input": [f"include:{rel(p)}" for p in itp_files],
        "out_sys": "system",
        "out_top": "coby_topol.top",
        "out_log": "coby.log",
        "sn": spec.name,
    })
    body = ",\n".join(f"    {k}={v!r}" for k, v in kwargs.items())
    script = (
        "# Generated by MDBuddy. Run from this directory: python coby_build.py\n"
        "import COBY\n\n"
        f"COBY.COBY(\n{body},\n)\n"
    )
    path = work / "coby_build.py"
    path.write_text(script)
    info = {"box_nm": box, "protein_cz_nm": round(cz, 3) if cz is not None else None,
            "proteins": len(proteins), "membranes": len(membranes), "notes": notes, "coby_kwargs": kwargs}
    if files.records:
        info["user_files"] = files.records
    return path, info


def _relpath(path: Path, start: Path) -> str:
    import os
    return os.path.relpath(path, start)


def run(spec, lock, entry, work: Path, script: Path) -> dict:
    lock.run(entry, [sys.executable, script.name], cwd=work, log_name="coby_stdout.log")
    gro, top = work / "system.gro", work / "coby_topol.top"
    if not gro.exists() or not top.exists():
        raise MembraneError("COBY finished but system.gro / coby_topol.top are missing; see coby.log")
    outputs = [gro, top]
    if (work / "system.pdb").exists():
        outputs.append(work / "system.pdb")
    return {"gro": gro, "top": top, "outputs": outputs}


# --- self-assembly ------------------------------------------------------------------------------

W_PER_NM3 = 55.56 / 4 * 6.02214076e23 * 1e-24  # Martini W beads per nm^3 at COBY's default water molarity
BEAD_VOLUME_NM3 = 0.1  # rough CG bead volume for box sizing


def _cross_section_nm2(cg_atoms) -> float:
    """Protein area in the membrane plane (beads within 1.5 nm of z = 0, oriented frame)."""
    pts = [(a.x / 10, a.y / 10) for a in cg_atoms if abs(a.z / 10) < 1.5]
    if len(pts) < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(pts).volume)  # 2D "volume" is the area
    except Exception:
        xs, ys = zip(*pts)
        return (max(xs) - min(xs)) * (max(ys) - min(ys)) * 0.7


def self_assembly_layout(spec, cg_atoms) -> dict:
    """Box, lipid numbers and water target for a protein + lipid 'soup'."""
    m, sa = spec.membrane, spec.data["self_assembly"]
    ex, ey, ez = (e / 10.0 for e in extents(cg_atoms))
    if m["box"] != "auto":
        box = [float(v) for v in m["box"]]
        area = box[0] * box[1]
    else:
        xy = math.ceil((max(ex, ey) + 2 * m["lipid_padding"]) * 2) / 2
        area = xy * xy
        box = [xy, xy, None]
    a_prot = _cross_section_nm2(cg_atoms)
    comp = spec.leaflets()["upper"]
    total_ratio = sum(float(v) for v in comp.values())
    n_total = (int(round(2 * max(area - a_prot, 0.25 * area) / float(m["apl"])))
               if sa["lipid_count"] == "auto" else int(sa["lipid_count"]))
    counts = {name: max(1, int(round(n_total * float(r) / total_ratio))) for name, r in comp.items()}
    lib = lipids.library_for(spec)
    lipid_volume = sum(n * (lipids.lookup(name, lib) or {}).get("beads", 12) * BEAD_VOLUME_NM3
                       for name, n in counts.items())
    water = int(sa["water_per_lipid"] * sum(counts.values()))
    volume = len(cg_atoms) * BEAD_VOLUME_NM3 + lipid_volume + water / W_PER_NM3
    if box[2] is None:
        z = max(volume / area, ez + 3.0, 10.0)
        box[2] = math.ceil(z * 2) / 2
    warnings = []
    if box[2] < 4.0 + 2.0:
        warnings.append("box is too thin for a bilayer plus water; increase membrane.box z")
    return {"box_nm": box, "lipid_counts": counts, "water_target": water,
            "protein_cross_section_nm2": round(a_prot, 1), "warnings": warnings}


SA_SCRIPT = '''\
# Generated by MDBuddy: protein + randomly placed lipids + water for self-assembly.
# Run from this directory: python coby_build.py
import COBY

SEED = {seed}
ITPS = {itps!r}
LIPIDS = {lipids!r}  # (name, COBY parameter library, number)


def single_lipid(name, params):
    """One lipid with COBY's own bead names and geometry, cut from a tiny bilayer."""
    COBY.COBY(box=[4, 4, 9], membrane=f"lipid:{{name}}:params:{{params}}", itp_input=ITPS,
              out_sys=f"template_{{name}}", out_log=f"template_{{name}}.log", sn=name, randseed=SEED, verbose=0)
    lines = open(f"template_{{name}}.gro").read().splitlines()
    first = lines[2][0:5]
    mol = [l for l in lines[2:-1] if l[0:5] == first]
    with open(f"{{name}}_single.gro", "w") as fh:
        fh.write(f"{{name}}\\n{{len(mol)}}\\n" + "\\n".join(mol) + "\\n   5.00000   5.00000   5.00000\\n")
    return f"file:{{name}}_single.gro moleculetype:{{name}} params:MDBUDDY_SA library_types:solute"


imports = [single_lipid(name, params) for name, params, _ in LIPIDS]
COBY.COBY(
    randseed=SEED,
    box={box!r},
    box_type="rectangular",
    params={params!r},
    protein={protein!r},
    molecule_import=imports,
    flooding=" ".join(f"solute:{{name}}:{{n}}" for name, _, n in LIPIDS) + " count:True params:MDBUDDY_SA",
    solvation={solvation!r},
    itp_input=ITPS,
    out_sys="system",
    out_top="coby_topol.top",
    out_log="coby.log",
    sn={name!r},
)
'''


def write_selfassembly_script(spec, work: Path, cg_pdb: Path, molecules, itp_files: list[Path]) -> tuple[Path, dict]:
    work.mkdir(parents=True, exist_ok=True)
    atoms, _ = read_pdb(cg_pdb)
    if not atoms:
        raise MembraneError(f"{cg_pdb} contains no beads")
    layout = self_assembly_layout(spec, atoms)
    lib = lipids.library_for(spec)
    rel = lambda p: Path(_relpath(p, work)).as_posix()
    moltypes = " ".join(f"moleculetypes:{name}:{count}" for name, count in molecules)
    script = SA_SCRIPT.format(
        seed=spec.seed, itps=[f"include:{rel(p)}" for p in itp_files],
        lipids=[(n, lipids.coby_params(n, lib) if lipids.lookup(n, lib) else spec.ff_info["coby_params"], c)
                for n, c in layout["lipid_counts"].items()],
        box=layout["box_nm"], params=spec.ff_info["coby_params"],
        protein=f"file:{rel(cg_pdb)} {moltypes}", solvation=solvation_arg(spec), name=spec.name)
    path = work / "coby_build.py"
    path.write_text(script)
    return path, layout
