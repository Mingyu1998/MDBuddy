#!/usr/bin/env python3
"""Regenerate src/mdbuddy/data/lipids.json from upstream sources.

Run this when a pinned source changes. It needs openpyxl, plus local checkouts of:
  * Martini-Force-Field-Initiative/M3-Lipid-Parameters (naming table + lipidome itps)
  * the COBY source package directory (the one containing molecule_definitions/)
  * the 2021 Martini 3 release phospholipid + sterol itps and the pinned Martini 2 lipid itps

    python tools/build_lipid_db.py --lipidome M3-Lipid-Parameters --coby COBY/src/COBY \
        --v1 martini_v3.0.0_phospholipids_v1.itp --sterol martini_v3.0_sterols_v1.0.itp \
        --m2 martini_v2.0_POPE_02.itp martini_v2.0_POPG_02.itp martini_v2.0_CDL0_01.itp \
        --lipidome-commit <sha> --coby-version 1.0.15

Nothing here is guessed: chemistry comes from the official naming spreadsheet or the itp
headers, bead lists from the itps, and COBY support from COBY's own definition files.
Hand-curated corrections are in CURATED below, each with its reason.
"""

from __future__ import annotations

import argparse
import collections
import importlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path

NON_LIPID_FILES = ("ions", "solvents")
NON_LIPID_PREFIX = ("INO",)  # free inositol sugars in the PI file

# Corrections to the spreadsheet's "Corresponding Martini 2 name" column and entries the
# spreadsheet does not cover. Keys are lipidome names.
CURATED_LIPIDOME = {
    # PIPs: the spreadsheet copies the M3 name into the M2 column (it lists M2 POP4..POP7,
    # which do not exist). Martini 2 has one generic model per phosphorylation level; ring
    # positions are not encoded in the M2 itps. Mapping below is by phosphorylation level.
    "POP1": {"m2": [], "note": "PI(3)P. Martini 2 has a single generic PIP (POP1) without ring position."},
    "POP4": {"m2": ["POP1"], "note": "PI(4)P. Closest M2: generic PIP POP1 (ring position not encoded in M2)."},
    "POP5": {"m2": [], "note": "PI(5)P. Martini 2 has a single generic PIP (POP1) without ring position."},
    "POP2": {"m2": [], "note": "PI(3,4)P2 - NOT the common PI(4,5)P2, which is POP6 in the lipidome."},
    "POP6": {"m2": ["POP2"], "note": "PI(4,5)P2. Closest M2: generic PIP2 POP2 (ring positions not encoded in M2)."},
    "POP7": {"m2": [], "note": "PI(3,5)P2."},
    "POP3": {"m2": ["POP3"], "note": "PI(3,4,5)P3."},
    "SAP1": {"m2": [], "note": "PI(3)P."},
    "SAP4": {"m2": ["PAP1"], "note": "PI(4)P. Closest M2: generic PIP PAP1."},
    "SAP5": {"m2": [], "note": "PI(5)P."},
    "SAP2": {"m2": [], "note": "PI(3,4)P2 - NOT PI(4,5)P2 (that is SAP6)."},
    "SAP6": {"m2": ["PAP2"], "note": "PI(4,5)P2. Closest M2: generic PIP2 PAP2."},
    "SAP7": {"m2": [], "note": "PI(3,5)P2."},
    "SAP3": {"m2": ["PAP3"], "note": "PI(3,4,5)P3."},
    # Cardiolipin is not in the spreadsheet.
    "TOCL": {"class": "CL", "chemistry": "tetra-C18:1", "common": "tetraoleoyl cardiolipin (TOCL)",
             "m2": ["CDL2", "CDL1", "CDL0"],
             "note": "Always 2- in Martini 3. M2 CDL0/CDL1/CDL2 are the 0/1-/2- protonation states with four "
                     "5-bead tails; only CDL2 has the same charge. Closest M3 match by tail type."},
    "TMCL": {"class": "CL", "chemistry": "tetra-C14:0", "common": "tetramyristoyl cardiolipin (TMCL)",
             "m2": [], "note": "Always 2- in Martini 3."},
    "CHOL": {"class": "sterol", "chemistry": "cholesterol", "common": "cholesterol", "m2": ["CHOL"],
             "note": "Martini 3 cholesterol v1.0 (Borges-Araujo 2023), a different bead topology from Martini 2 CHOL."},
}

CURATED_M2 = {
    "CDL0": {"class": "CL", "chemistry": "cardiolipin, four 5-bead mono-unsaturated tails", "charge": 0.0,
             "common": "cardiolipin 0 (both phosphates neutral)"},
    "CDL1": {"class": "CL", "chemistry": "cardiolipin, four 5-bead mono-unsaturated tails", "charge": -1.0,
             "common": "cardiolipin 1-"},
    "CDL2": {"class": "CL", "chemistry": "cardiolipin, four 5-bead mono-unsaturated tails", "charge": -2.0,
             "common": "cardiolipin 2-"},
    "CHOL": {"class": "sterol", "chemistry": "cholesterol", "common": "cholesterol"},
    "POP1": {"note": "Generic PIP (single phosphate); ring position not encoded. Often used as PI(4)P."},
    "POP2": {"note": "Generic PIP2; ring positions not encoded. Often used as PI(4,5)P2 -> lipidome POP6."},
    "POP3": {"note": "Generic PIP3, PI(3,4,5)P3."},
}


def read_itp(path: Path) -> dict:
    """moleculetype -> {atoms, charge, description} for one itp."""
    out, sec, cur, comments = {}, None, None, []
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if line.startswith(";"):
            text = line.lstrip(";").strip()
            if text:
                comments.append(text)
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m:
            sec = m.group(1)
            continue
        if not line or line.startswith("#"):
            continue
        fields = line.split(";")[0].split()
        if sec == "moleculetype":
            cur = fields[0]
            desc = ""
            for c in reversed(comments):
                low = c.lower()
                if ("corresponding to" in c or low.startswith(("a general model", "a martini"))
                        or "lipid topology for" in low):
                    desc = c
                    break
            out[cur] = {"atoms": [], "charge": 0.0, "description": desc}
            comments = []
        elif sec == "atoms" and cur:
            out[cur]["atoms"].append(fields[4])
            if len(fields) > 6:
                out[cur]["charge"] += float(fields[6])
    for v in out.values():
        v["charge"] = round(v["charge"], 3)
    return out


def load_coby(coby_dir: Path) -> dict:
    """params -> lipid -> bead list, read from COBY's definition files (COBY itself is not run)."""
    stub = Path(tempfile.mkdtemp())
    (stub / "COBY").mkdir()
    (stub / "COBY" / "__init__.py").write_text("")
    for sub in ("molecule_definitions", "fragment_definitions"):
        os.symlink(coby_dir / sub, stub / "COBY" / sub)
    sys.path.insert(0, str(stub))
    md = importlib.import_module("COBY.molecule_definitions")
    lib: dict = collections.defaultdict(dict)
    for (_, params), d in md.lipid_scaffolds.items():
        for key, val in d.get("lipids", {}).items():
            if isinstance(key, tuple) and key[1] == "beads":
                lib[params][key[0]] = [t for t in val.split() if t != "-"]
    return lib


def compare(coby_beads, itp_atoms) -> dict:
    """COBY bead list vs itp atom list -> status (+ positional renames when safe)."""
    if coby_beads is None:
        return {"status": "missing"}
    if coby_beads == itp_atoms:
        return {"status": "exact"}
    if len(coby_beads) == len(itp_atoms) and sorted(coby_beads) != sorted(itp_atoms):
        # Same count, same order of bead roles, only some names differ: safe positional rename
        renames = {i: [a, b] for i, (a, b) in enumerate(zip(coby_beads, itp_atoms)) if a != b}
        if all(a not in itp_atoms for a, _ in renames.values()):
            return {"status": "rename", "renames": {str(k): v for k, v in renames.items()}}
    reason = ("bead count differs" if len(coby_beads) != len(itp_atoms) else "bead order differs")
    return {"status": "unsupported", "reason": f"COBY {reason} from the itp"}


def table_rows(xlsx: Path, sheet: str, start: int, cols: dict) -> list[dict]:
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, data_only=True, read_only=True)
    rows, cur_file = [], None
    for r in list(wb[sheet].iter_rows(values_only=True))[start:]:
        if "file" in cols and isinstance(r[cols["file"]], str) and r[cols["file"]].endswith(".itp"):
            cur_file = r[cols["file"]]
            continue
        if r[cols["name"]]:
            row = {k: r[i] for k, i in cols.items() if k != "file"}
            row["file"] = cur_file
            rows.append(row)
    return rows


def clean(s):
    return re.sub(r"\s+", " ", str(s).replace("\\n", " ")).strip() if s else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lipidome", type=Path, required=True)
    ap.add_argument("--coby", type=Path, required=True)
    ap.add_argument("--v1", type=Path, required=True)
    ap.add_argument("--sterol", type=Path, required=True)
    ap.add_argument("--m2", type=Path, nargs="+", required=True)
    ap.add_argument("--lipidome-commit", required=True)
    ap.add_argument("--coby-version", required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "src/mdbuddy/data/lipids.json")
    a = ap.parse_args()

    coby = load_coby(a.coby)
    xlsx = a.lipidome / "tools" / "CG-Martini-lipids-naming.xlsx"

    # ---- Martini 3 lipidome ---------------------------------------------------------
    m3_rows = {r["name"]: r for r in table_rows(xlsx, "Martini 3 - All lipids", 5, {
        "file": 2, "cls": 2, "name": 5, "sn1": 6, "sn2": 7, "chemistry": 9, "common": 10,
        "m2": 26, "c36": 27})}
    lipidome = {}
    for itp in sorted((a.lipidome / "ITPs").glob("*.itp")):
        if any(t in itp.name for t in NON_LIPID_FILES) or itp.name in ("martini_v3.0.0.itp",
                                                                      "martini_v3.0.0_ffbonded_v2.itp"):
            continue
        for name, mol in read_itp(itp).items():
            if name.startswith(NON_LIPID_PREFIX):
                continue
            row = m3_rows.get(name, {})
            cls = row.get("cls") or re.sub(r"^martini_v3\.0\.0_|_v\d+\.itp$", "", itp.name)
            entry = {
                # True when the Martini 2 correspondence is known (naming table or curated below)
                "m2_mapping_known": bool(row) or "m2" in CURATED_LIPIDOME.get(name, {}),
                "class": cls,
                "chemistry": clean(row.get("chemistry")),
                "common": clean(row.get("common")),
                "sn1": row.get("sn1"), "sn2": row.get("sn2"),
                "charge": mol["charge"], "beads": len(mol["atoms"]),
                "itp": itp.name, "charmm36": clean(row.get("c36")),
                "m2": [clean(row["m2"])] if row.get("m2") else [],
                "description": mol["description"] or None,
                "coby_params": "LTF",
                "coby": compare(coby["LTF"].get(name), mol["atoms"]),
            }
            if name in CURATED_LIPIDOME:
                entry.update(CURATED_LIPIDOME[name])
            if name == "CHOL":
                entry["coby_params"] = "default"
                entry["coby"] = compare(coby["default"].get(name), mol["atoms"])
                entry["itp"] = a.sterol.name  # identical parameters to the lipidome copy
            lipidome[name] = entry

    # ---- Martini 2 (naming sheet) ----------------------------------------------------
    m2_rows = table_rows(xlsx, "Martini 2 - All lipids", 5, {
        "cls": 1, "name": 4, "sn1": 5, "sn2": 6, "chemistry": 8, "common": 9, "description": 10})
    pinned_m2 = {}
    for p in a.m2:
        for name, mol in read_itp(p).items():
            pinned_m2[name] = (p.name, mol)
    martini2 = {}
    for r in m2_rows:
        martini2[r["name"]] = {"class": r["cls"], "chemistry": clean(r["chemistry"]), "common": clean(r["common"]),
                               "sn1": r["sn1"], "sn2": r["sn2"], "description": clean(r["description"])}
    for name, extra in CURATED_M2.items():
        martini2.setdefault(name, {}).update(extra)
    for name, entry in martini2.items():
        if name in pinned_m2:
            itp_name, mol = pinned_m2[name]
            entry["pinned_itp"] = itp_name
            entry["charge"] = mol["charge"]
            entry["coby"] = compare(coby["M2"].get(name), mol["atoms"])
        else:
            entry["coby"] = ({"status": "unverified", "reason": "COBY has a definition; no pinned itp to compare"}
                             if name in coby["M2"] else {"status": "missing"})
        entry["coby_params"] = "M2"

    # ---- Martini 3, 2021 release (uses Martini 2 naming) -------------------------------
    v1 = read_itp(a.v1)
    v1.update(read_itp(a.sterol))
    release2021 = {}
    for name, mol in v1.items():
        same = martini2.get(name, {})
        release2021[name] = {
            "class": same.get("class"), "chemistry": same.get("chemistry"), "common": same.get("common"),
            "charge": mol["charge"], "beads": len(mol["atoms"]), "itp": (a.sterol if name == "CHOL" else a.v1).name,
            "description": mol["description"] or None,
            "coby_params": "default", "coby": compare(coby["default"].get(name), mol["atoms"]),
        }

    db = {
        "_about": "Generated by tools/build_lipid_db.py - do not edit by hand.",
        "sources": {
            "lipidome": f"Martini-Force-Field-Initiative/M3-Lipid-Parameters@{a.lipidome_commit}",
            "naming_table": "tools/CG-Martini-lipids-naming.xlsx in the lipidome repository",
            "coby": f"COBY {a.coby_version}",
        },
        "libraries": {"martini3_lipidome": lipidome, "martini3_2021": release2021, "martini2": martini2},
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(db, indent=1, sort_keys=True) + "\n")
    stats = {lib: collections.Counter(v["coby"]["status"] for v in d.values()) for lib, d in db["libraries"].items()}
    print(f"wrote {a.out}")
    for lib, c in stats.items():
        print(f"  {lib}: {sum(c.values())} lipids, COBY {dict(c)}")


if __name__ == "__main__":
    main()
