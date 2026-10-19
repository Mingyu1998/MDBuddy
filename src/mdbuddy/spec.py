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
    # scfix: whether the force field defines the side-chain angle/torsion corrections that
    # martinize2 >= 0.11 applies by default. Only Martini 3 does; the Martini 2 force fields
    # need -noscfix or martinize2 warns and (with -maxwarn 0) writes nothing.
    "martini3": {"martinize_ff": "martini3001", "coby_params": "default", "scfix": True,
                 "elastic": {"force": 700.0, "lower": 0.5, "upper": 0.9}},
    "martini22": {"martinize_ff": "martini22", "coby_params": "M2", "scfix": False,
                  "elastic": {"force": 500.0, "lower": 0.5, "upper": 0.9}},
    "elnedyn22": {"martinize_ff": "elnedyn22", "coby_params": "M2", "scfix": False,
                  "elastic": {"force": 500.0, "lower": 0.0, "upper": 0.9}},
}
PROTEIN_MODELS = ("elastic", "go", "none")
SOURCES = ("opm", "rcsb", "file", "none")
ORIENTATIONS = ("opm", "as_is", "manual", "hydrophobic", "self_assembly")
TOPOLOGY_RULES = ("positive_inside", "as_found", "flip")
POSRES = ("none", "backbone", "all")
TC_GROUP_CHOICES = (["Protein", "MEMBRANE", "SOLVENT"], ["SOLU_MEMB", "SOLVENT"], ["MEMBRANE", "SOLVENT"], ["System"])
CG_PROTEIN_KEYS = {"structure", "topology", "itps", "moleculetypes", "placement", "cx", "cy", "cz", "rx", "ry", "rz",
                   "coby_args"}
# COBY arguments MDBuddy generates itself; coby.kwargs may not override them
COBY_RESERVED = {"box", "box_type", "protein", "membrane", "solvation", "itp_input", "out_sys", "out_top", "out_log",
                 "sn", "randseed", "params"}
NON_RECURSIVE = ("lipids", "lipid_apl", "kwargs")
COMM_GROUP_CHOICES = (["SOLU_MEMB", "SOLVENT"], ["System"])

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
        "repair": "auto",  # rebuild missing heavy atoms with PDBFixer: auto | off
        "orientation": None,  # filled from source: opm -> "opm", else "as_is"
        "membrane_center_z": 0.0,  # Angstrom, in the input frame (orientation: as_is)
        # orientation: manual (applied in this order; all optional)
        "manual": {"normal": None, "center_point": None, "align": None, "rotate": None,
                   "center": None, "shift_z": None},
        "topology": "positive_inside",  # orientation: hydrophobic -> which side is cytoplasmic (-z)
        "pore_correction": "auto",  # orientation: hydrophobic -> lipid-facing weighting (beta-barrels)
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
        "box": "auto",  # auto | [x, y, z] nm | from_structure (box of the first cg_proteins structure)
        "cz": 0.0,  # nm: membrane midplane offset from the box centre along z
        "lipid_apl": {},  # per-lipid area per lipid (nm^2), e.g. {POPE: 0.62}; tunes local lipid density
        "box_type": "rectangular",
        "lipid_padding": 3.0,
        "water_padding": 2.5,
        "min_box_z": 10.0,
    },
    # Pre-coarse-grained proteins (e.g. protein islands from earlier simulations), inserted as given:
    # [{structure, topology | itps + moleculetypes, placement: keep | center, cx..rz, coby_args}]
    "cg_proteins": [],
    # Direct COBY control, all recorded: extra sub-arguments, raw membrane strings, other COBY keywords
    "coby": {
        "membrane_args": "",  # appended to the generated membrane argument (e.g. "hole:circle:...")
        "membranes": [],  # full raw COBY membrane strings; replace the generated membrane (several allowed)
        "protein_args": "",  # appended to the atomistic protein's COBY argument
        "solvation_args": "",  # appended to the solvation argument
        "kwargs": {},  # any other COBY keyword argument (e.g. import_library, molecule_import, flooding)
    },
    "self_assembly": {  # protein.orientation: self_assembly
        "lipid_count": "auto",  # or a number; auto fills the box area at membrane.apl
        "water_per_lipid": 25,  # CG water beads per lipid (Martini W = 4 waters)
        "run_ns": 1000,  # length of the self-assembly run
        "restrain_protein": True,  # backbone restraints during self-assembly (protein stays put)
        "align": "hydrophobic",  # initial protein alignment: hydrophobic | as_is
    },
    "solvent": {
        "water": "W",
        "cation": "NA",
        "anion": "CL",
        "salt_molarity": 0.15,
    },
    "forcefield_files": {
        "martini3_lipids": "lipidome",
        "optional": [],
        "extra_itps": [],
        "ff_dir": None,
    },
    "simulation": {
        "temperature": 310.0,
        "pressure": 1.0,
        "production_ns": 1000.0,
        "equilibration": "standard",
        "output_ps": 500.0,
        "coupling_groups": ["Protein", "MEMBRANE", "SOLVENT"],
        "comm_groups": ["SOLU_MEMB", "SOLVENT"],
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
        if isinstance(defaults[key], dict) and key not in NON_RECURSIVE and isinstance(value, dict):
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
        """Lipids of the generated membrane plus those named in raw coby.membranes strings."""
        names: list[str] = []
        if not self.data["coby"]["membranes"]:
            for comp in self.leaflets().values():
                for lipid in comp:
                    if lipid not in names:
                        names.append(lipid)
        for raw in self.data["coby"]["membranes"]:
            for tok in str(raw).split():
                parts = tok.split(":")
                name = None
                if parts[0].lower() == "lipid" and len(parts) > 1:
                    name = parts[2 + parts[1:].index("name")] if "name" in parts[1:-1] else parts[1]
                elif parts[0].lower() == "lipid_extra" and "name" in parts[1:-1]:
                    name = parts[2 + parts[1:].index("name")]
                if name and name not in names:
                    names.append(name)
        return names

    @property
    def has_atomistic_protein(self) -> bool:
        return self.protein["source"] != "none"

    @property
    def cg_proteins(self) -> list[dict]:
        return self.data["cg_proteins"]

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
    cgp = d["cg_proteins"]
    if not isinstance(cgp, list):
        problems.append("cg_proteins must be a list")
        cgp = []
    for i, e in enumerate(cgp):
        where = f"cg_proteins[{i}]"
        if not isinstance(e, dict) or not e.get("structure"):
            problems.append(f"{where}: needs 'structure' (a CG .pdb or .gro)")
            continue
        unknown = set(e) - CG_PROTEIN_KEYS
        if unknown:
            problems.append(f"{where}: unknown key(s) {sorted(unknown)}")
        if not e.get("topology") and not (e.get("itps") and e.get("moleculetypes")):
            problems.append(f"{where}: give 'topology' (a .top with [ molecules ]) or both 'itps' and 'moleculetypes'")
        if e.get("moleculetypes") is not None and not all(
                isinstance(m, list) and len(m) == 2 and isinstance(m[1], int) for m in e["moleculetypes"]):
            problems.append(f"{where}.moleculetypes must be a list of [name, count] in structure order")
        if e.get("placement", "keep") not in ("keep", "center"):
            problems.append(f"{where}.placement must be 'keep' (file coordinates) or 'center'")
    if p["source"] == "none":
        if p["orientation"] == "self_assembly":
            problems.append("self-assembly needs an atomistic protein (protein.source opm/rcsb/file)")
    if cgp and p["orientation"] == "self_assembly":
        problems.append("cg_proteins cannot be combined with self-assembly yet")
    if d["membrane"]["box"] == "auto" and p["source"] == "none" and not cgp:
        problems.append("a membrane-only system needs an explicit membrane.box: [x, y, z]")
    if d["membrane"]["box"] == "auto" and (p["source"] != "none") + len(cgp) > 1:
        problems.append("membrane.box: auto sizes the box around one protein; with more than one protein "
                        "set membrane.box to [x, y, z] or from_structure")
    if d["membrane"]["box"] == "from_structure" and not cgp:
        problems.append("membrane.box: from_structure needs cg_proteins (it uses the first structure's box)")
    la = d["membrane"]["lipid_apl"]
    if not isinstance(la, dict) or not all(isinstance(v, (int, float)) and v > 0 for v in la.values()):
        problems.append("membrane.lipid_apl must map lipid names to areas per lipid in nm^2")
    c = d["coby"]
    for key in ("membrane_args", "protein_args", "solvation_args"):
        if not isinstance(c[key], str):
            problems.append(f"coby.{key} must be a string of COBY sub-arguments")
    if not isinstance(c["membranes"], list) or not all(isinstance(m, str) for m in c["membranes"]):
        problems.append("coby.membranes must be a list of COBY membrane strings")
    if not isinstance(c["kwargs"], dict):
        problems.append("coby.kwargs must be a mapping of COBY keyword arguments")
    else:
        # COBY also treats any keyword starting with prot/memb/solv as protein/membrane/solvation
        clash = sorted(k for k in c["kwargs"] if k in COBY_RESERVED or k.startswith(("prot", "memb", "solv")))
        if clash:
            problems.append(f"coby.kwargs may not set {clash}: MDBuddy generates them (use membrane.*, "
                            "coby.membranes, coby.*_args or the spec instead)")
    if p["repair"] not in ("auto", "off"):
        problems.append("protein.repair must be 'auto' (rebuild missing heavy atoms with PDBFixer) or 'off'")
    if p["orientation"] not in ORIENTATIONS:
        problems.append(f"protein.orientation must be one of {ORIENTATIONS}")
    if p["topology"] not in TOPOLOGY_RULES:
        problems.append(f"protein.topology must be one of {TOPOLOGY_RULES}")
    if p["pore_correction"] not in ("auto", "on", "off"):
        problems.append("protein.pore_correction must be auto, on or off")
    man = p["manual"]
    if p["orientation"] == "manual":
        if not any(v is not None for v in man.values()):
            problems.append("protein.orientation: manual needs at least one of protein.manual.normal, "
                            "center_point, align, rotate, center, shift_z (or use `mdbuddy orient view`)")
        for key, n in (("normal", 3), ("center_point", 3), ("rotate", 3), ("align", 2)):
            v = man[key]
            if v is not None and (not isinstance(v, list) or len(v) != n):
                problems.append(f"protein.manual.{key} must be a list of {n} values")
        if man["normal"] is not None and man["align"] is not None:
            problems.append("protein.manual: use either normal or align, not both")
    sa = d["self_assembly"]
    if sa["lipid_count"] != "auto" and (not isinstance(sa["lipid_count"], int) or sa["lipid_count"] < 20):
        problems.append("self_assembly.lipid_count must be 'auto' or an integer >= 20")
    if not isinstance(sa["water_per_lipid"], (int, float)) or sa["water_per_lipid"] < 5:
        problems.append("self_assembly.water_per_lipid must be >= 5")
    if sa["align"] not in ("hydrophobic", "as_is"):
        problems.append("self_assembly.align must be 'hydrophobic' or 'as_is'")
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
    if box not in ("auto", "from_structure"):
        if (not isinstance(box, list) or len(box) != 3
                or not all(isinstance(v, (int, float)) and v > 0 for v in box)):
            problems.append("membrane.box must be 'auto', 'from_structure' or [x, y, z] in nm")
    if m["box_type"] not in ("rectangular",):
        problems.append("membrane.box_type: only 'rectangular' is supported for membranes right now")
    for key in ("apl", "lipid_padding", "water_padding"):
        if not isinstance(m[key], (int, float)) or m[key] <= 0:
            problems.append(f"membrane.{key} must be a positive number")

    if d["forcefield_files"]["martini3_lipids"] not in ("lipidome", "release2021"):
        problems.append("forcefield_files.martini3_lipids must be 'lipidome' (Martini 3 lipidome, 2025) "
                        "or 'release2021' (lipids of the original Martini 3 release)")
    if not isinstance(s["salt_molarity"], (int, float)) or s["salt_molarity"] < 0:
        problems.append("solvent.salt_molarity must be >= 0")

    if sim["equilibration"] not in ("standard",):
        problems.append("simulation.equilibration: only 'standard' is available")
    if not isinstance(sim["production_ns"], (int, float)) or sim["production_ns"] <= 0:
        problems.append("simulation.production_ns must be > 0")
    if sim["coupling_groups"] not in TC_GROUP_CHOICES:
        problems.append("simulation.coupling_groups must be one of "
                        + ", ".join("[" + ", ".join(g) + "]" for g in TC_GROUP_CHOICES))
    if sim["comm_groups"] not in COMM_GROUP_CHOICES:
        problems.append("simulation.comm_groups must be one of "
                        + ", ".join("[" + ", ".join(g) + "]" for g in COMM_GROUP_CHOICES))
    out_ps = sim["output_ps"]
    # production dt is 20 fs and energies are computed every 100 steps -> multiples of 2 ps
    if (not isinstance(out_ps, (int, float)) or out_ps <= 0
            or abs(out_ps / 2.0 - round(out_ps / 2.0)) > 1e-9):
        problems.append("simulation.output_ps must be a positive multiple of 2 ps")


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
    for e in data["cg_proteins"] if isinstance(data["cg_proteins"], list) else []:
        if isinstance(e, dict):
            e.setdefault("placement", "keep")
            for k in ("cx", "cy", "cz", "rx", "ry", "rz"):
                e.setdefault(k, 0.0)
            e.setdefault("coby_args", "")
    no_protein = data["protein"]["source"] == "none" and not data["cg_proteins"]
    if no_protein and "coupling_groups" not in (user.get("simulation") or {}):
        data["simulation"]["coupling_groups"] = ["MEMBRANE", "SOLVENT"]

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
  # orientation: opm     # opm | as_is | manual | hydrophobic | self_assembly
  #   no OPM entry?  hydrophobic = educated guess from exposed hydrophobicity,
  #   manual = you place it (`mdbuddy orient view protein.pdb` helps), self_assembly = let lipids find it

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
  output_ps: 500         # trajectory/energy output interval in production
"""
