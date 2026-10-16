"""System specification: the single YAML file that fully describes a build.

Everything MDBuddy does is driven by a spec. Defaults are filled in here, so the
resolved spec written next to the outputs (``spec.resolved.yaml``) is enough to
rebuild the same system.
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

SPEC_VERSION = 1

FORCEFIELDS = {
    # name in spec -> martinize2 -ff, COBY parameter library, elastic defaults
    "martini3": {"martinize_ff": "martini3001", "coby_params": "default",
                 "elastic": {"force": 700.0, "lower": 0.5, "upper": 0.9}},
    "martini22": {"martinize_ff": "martini22", "coby_params": "M2",
                  "elastic": {"force": 500.0, "lower": 0.5, "upper": 0.9}},
    "elnedyn22": {"martinize_ff": "elnedyn22", "coby_params": "M2",
                  "elastic": {"force": 500.0, "lower": 0.0, "upper": 0.9}},
}
PROTEIN_MODELS = ("elastic", "go", "none")
SOURCES = ("opm", "rcsb", "file")
ORIENTATIONS = ("opm", "as_is")
POSRES = ("none", "backbone", "all")

DEFAULTS: dict[str, Any] = {
    "mdbuddy": SPEC_VERSION,
    "name": None,
    "seed": 2026,
    "output_dir": None,
    "protein": {
        "source": "opm",
        "id": None,
        "path": None,
        "chains": None,
        "keep_hetatm": [],
        "orientation": None,  # filled from source: opm -> "opm", else "as_is"
        "membrane_center_z": 0.0,  # Angstrom, in the input frame
    },
    "coarse_grain": {
        "forcefield": "martini3",
        "model": "elastic",
        "elastic": None,  # filled from forcefield defaults
        "go": {"epsilon": 9.414, "low": 0.3, "up": 1.1, "res_dist": 3},
        "merge_chains": "auto",
        "dssp": "auto",
        "position_restraints": "backbone",
        "posres_force": 1000.0,
        "cys": "auto",
        "molecule_name": "protein",
        "maxwarn": 0,
        "extra_args": [],
    },
    "membrane": {
        "lipids": {"POPC": 1},
        "upper": None,
        "lower": None,
        "apl": 0.6,
        "upper_apl": None,
        "lower_apl": None,
        "box": "auto",
        "box_type": "rectangular",
        "lipid_padding": 3.0,
        "water_padding": 2.5,
        "min_box_z": 10.0,
    },
    "solvent": {
        "water": "W",
        "cation": "NA",
        "anion": "CL",
        "salt_molarity": 0.15,
    },
    "forcefield_files": {
        "optional": [],
        "extra_itps": [],
        "ff_dir": None,
    },
    "simulation": {
        "temperature": 310.0,
        "pressure": 1.0,
        "production_ns": 1000.0,
        "equilibration": "standard",
        "coupling_groups": ["SOLU_MEMB", "SOLV"],
        "gromacs_maxwarn": 0,
    },
}


class SpecError(ValueError):
    """Raised when a spec is invalid. ``problems`` lists every issue found."""

    def __init__(self, problems: list[str]):
        self.problems = problems
        super().__init__("Invalid MDBuddy spec:\n  - " + "\n  - ".join(problems))


def _merge(defaults: dict, user: dict, path: str, problems: list[str]) -> dict:
    out = copy.deepcopy(defaults)
    for key, value in user.items():
        if key not in defaults:
            problems.append(f"unknown key '{path}{key}'")
            continue
        if isinstance(defaults[key], dict) and key not in ("lipids",) and isinstance(value, dict):
            out[key] = _merge(defaults[key], value, f"{path}{key}.", problems)
        else:
            out[key] = value
    return out


@dataclass
class Spec:
    data: dict
    source_path: Optional[Path] = None
    problems: list[str] = field(default_factory=list)

    # Convenience accessors -------------------------------------------------
    @property
    def name(self) -> str:
        return self.data["name"]

    @property
    def protein(self) -> dict:
        return self.data["protein"]

    @property
    def cg(self) -> dict:
        return self.data["coarse_grain"]

    @property
    def membrane(self) -> dict:
        return self.data["membrane"]

    @property
    def solvent(self) -> dict:
        return self.data["solvent"]

    @property
    def ff_files(self) -> dict:
        return self.data["forcefield_files"]

    @property
    def sim(self) -> dict:
        return self.data["simulation"]

    @property
    def ff_info(self) -> dict:
        return FORCEFIELDS[self.cg["forcefield"]]

    @property
    def seed(self) -> int:
        return int(self.data["seed"])

    def leaflets(self) -> dict[str, dict[str, float]]:
        m = self.membrane
        return {
            "upper": m["upper"] or m["lipids"],
            "lower": m["lower"] or m["lipids"],
        }

    def lipid_names(self) -> list[str]:
        names: list[str] = []
        for comp in self.leaflets().values():
            for lipid in comp:
                if lipid not in names:
                    names.append(lipid)
        return names

    def output_dir(self, base: Optional[Path] = None) -> Path:
        out = self.data["output_dir"] or self.name
        out = Path(out)
        if not out.is_absolute():
            root = base or (self.source_path.parent if self.source_path else Path.cwd())
            out = root / out
        return out

    def resolve_path(self, p: str) -> Path:
        path = Path(p).expanduser()
        if not path.is_absolute() and self.source_path is not None:
            path = self.source_path.parent / path
        return path

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.data, sort_keys=False)


_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def _validate(d: dict, problems: list[str]) -> None:
    p, cg, m, s, sim = (d["protein"], d["coarse_grain"], d["membrane"],
                        d["solvent"], d["simulation"])

    if not d["name"] or not _NAME_RE.match(str(d["name"])):
        problems.append("'name' is required and may only contain letters, digits, '_', '-', '.'")

    if p["source"] not in SOURCES:
        problems.append(f"protein.source must be one of {SOURCES}")
    elif p["source"] in ("opm", "rcsb") and not p["id"]:
        problems.append(f"protein.id is required for source '{p['source']}'")
    elif p["source"] == "file" and not p["path"]:
        problems.append("protein.path is required for source 'file'")
    if p["orientation"] not in ORIENTATIONS:
        problems.append(f"protein.orientation must be one of {ORIENTATIONS}")
    if p["chains"] is not None and not isinstance(p["chains"], list):
        problems.append("protein.chains must be a list of chain IDs or null")

    if cg["forcefield"] not in FORCEFIELDS:
        problems.append(f"coarse_grain.forcefield must be one of {tuple(FORCEFIELDS)}")
    if cg["model"] not in PROTEIN_MODELS:
        problems.append(f"coarse_grain.model must be one of {PROTEIN_MODELS}")
    if cg["model"] == "go" and cg["forcefield"] != "martini3":
        problems.append("the Go model (coarse_grain.model: go) is only available with martini3")
    if cg["position_restraints"] not in POSRES:
        problems.append(f"coarse_grain.position_restraints must be one of {POSRES}")
    if not _NAME_RE.match(str(cg["molecule_name"])):
        problems.append("coarse_grain.molecule_name may only contain letters, digits, '_', '-', '.'")
    mc = cg["merge_chains"]
    if not (mc in ("auto", "all", "none") or isinstance(mc, list)):
        problems.append("coarse_grain.merge_chains must be auto, all, none, or a list like [[A, B], [C, D]]")

    for leaflet, comp in (("lipids", m["lipids"]), ("upper", m["upper"]), ("lower", m["lower"])):
        if comp is None:
            continue
        if not isinstance(comp, dict) or not comp:
            problems.append(f"membrane.{leaflet} must be a mapping of lipid name -> ratio")
            continue
        for lipid, ratio in comp.items():
            if not isinstance(ratio, (int, float)) or ratio <= 0:
                problems.append(f"membrane.{leaflet}.{lipid}: ratio must be a positive number")
    box = m["box"]
    if box != "auto":
        if (not isinstance(box, list) or len(box) != 3
                or not all(isinstance(v, (int, float)) and v > 0 for v in box)):
            problems.append("membrane.box must be 'auto' or [x, y, z] in nm")
    if m["box_type"] not in ("rectangular",):
        problems.append("membrane.box_type: only 'rectangular' is supported for membranes right now")
    for key in ("apl", "lipid_padding", "water_padding"):
        if not isinstance(m[key], (int, float)) or m[key] <= 0:
            problems.append(f"membrane.{key} must be a positive number")

    if not isinstance(s["salt_molarity"], (int, float)) or s["salt_molarity"] < 0:
        problems.append("solvent.salt_molarity must be >= 0")

    if sim["equilibration"] not in ("standard",):
        problems.append("simulation.equilibration: only 'standard' is available")
    if not isinstance(sim["production_ns"], (int, float)) or sim["production_ns"] <= 0:
        problems.append("simulation.production_ns must be > 0")
    groups = sim["coupling_groups"]
    if groups not in (["SOLU_MEMB", "SOLV"], ["PROT", "MEMB", "SOLV"], ["System"]):
        problems.append("simulation.coupling_groups must be [SOLU_MEMB, SOLV], [PROT, MEMB, SOLV] or [System]")


def from_dict(user: dict, source_path: Optional[Path] = None) -> Spec:
    problems: list[str] = []
    if not isinstance(user, dict):
        raise SpecError(["spec must be a YAML mapping"])
    version = user.get("mdbuddy", SPEC_VERSION)
    if version != SPEC_VERSION:
        problems.append(f"unsupported spec version mdbuddy: {version} (expected {SPEC_VERSION})")
    data = _merge(DEFAULTS, user, "", problems)

    if data["protein"]["orientation"] is None:
        data["protein"]["orientation"] = "opm" if data["protein"]["source"] == "opm" else "as_is"
    ff = data["coarse_grain"]["forcefield"]
    if ff in FORCEFIELDS:
        elastic = dict(FORCEFIELDS[ff]["elastic"])
        user_elastic = data["coarse_grain"]["elastic"] or {}
        unknown = set(user_elastic) - set(elastic)
        if unknown:
            problems.append(f"unknown key(s) in coarse_grain.elastic: {sorted(unknown)}")
        elastic.update(user_elastic)
        data["coarse_grain"]["elastic"] = elastic
    if isinstance(data["protein"]["id"], str):
        data["protein"]["id"] = data["protein"]["id"].lower()

    _validate(data, problems)
    if problems:
        raise SpecError(problems)
    return Spec(data=data, source_path=source_path)


def load(path: str | Path) -> Spec:
    path = Path(path).resolve()
    with open(path) as fh:
        user = yaml.safe_load(fh) or {}
    return from_dict(user, source_path=path)


TEMPLATE = """\
# MDBuddy system spec. Run `mdbuddy build {name}.yaml` to build it.
# Every option is documented in skills/mdbuddy/references/spec_reference.md
mdbuddy: 1
name: {name}
seed: 2026

protein:
  source: opm            # opm (oriented, recommended) | rcsb | file
  id: {pdb_id}           # PDB ID for opm/rcsb
  # path: my_protein.pdb # for source: file
  chains: null           # e.g. [A, B]; null keeps all chains
  # orientation: as_is   # required for non-OPM input; then set membrane_center_z (Angstrom)

coarse_grain:
  forcefield: martini3   # martini3 | martini22 | elnedyn22
  model: elastic         # elastic | go (martini3 only) | none
  # elastic: {{force: 700, lower: 0.5, upper: 0.9}}
  position_restraints: backbone

membrane:
  lipids: {{POPC: 1}}      # symmetric; or give upper: / lower: separately
  # upper: {{POPC: 7, CHOL: 3}}
  # lower: {{POPE: 5, POPS: 2, CHOL: 3}}
  apl: 0.6               # nm^2 per lipid
  box: auto              # or [x, y, z] in nm

solvent:
  salt_molarity: 0.15    # NaCl, mol/L, added on top of neutralising ions

simulation:
  temperature: 310       # K
  production_ns: 1000
"""
