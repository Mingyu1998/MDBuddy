"""Clean the input structure before martinize2.

Removes OPM dummy atoms, waters, hydrogens and (by default) all hetero groups,
selects chains, resolves alternate locations, and reports gaps so the agent
can tell the user about them.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

from .structure import DUMMY_RESNAMES, WATER_RESNAMES, Atom, opm_half_thickness, read_pdb, write_pdb

# Gap detection: CA-CA distance above this (Angstrom) between consecutive residues
CA_GAP_CUTOFF = 4.2


class PrepareError(RuntimeError):
    pass


def _pick_altlocs(atoms: list[Atom]) -> list[Atom]:
    """Keep the first alternate location of each atom (usually 'A')."""
    seen: set = set()
    kept = []
    for a in atoms:
        key = (a.chain, a.resid, a.icode, a.name)
        if a.altloc:
            if key in seen:
                continue
            a.altloc = ""
        seen.add(key)
        kept.append(a)
    return kept


def _find_gaps(atoms: list[Atom]) -> list[dict]:
    gaps = []
    by_chain: "OrderedDict[str, list[Atom]]" = OrderedDict()
    for a in atoms:
        if a.record == "ATOM" and a.name == "CA":
            by_chain.setdefault(a.chain, []).append(a)
    for chain, cas in by_chain.items():
        for prev, cur in zip(cas, cas[1:]):
            d = ((cur.x - prev.x) ** 2 + (cur.y - prev.y) ** 2 + (cur.z - prev.z) ** 2) ** 0.5
            if d > CA_GAP_CUTOFF:
                gaps.append({"chain": chain, "after_resid": prev.resid, "before_resid": cur.resid,
                             "missing_residues": max(cur.resid - prev.resid - 1, 0),
                             "ca_distance": round(d, 2)})
    return gaps


def prepare(spec, in_pdb: Path, out_pdb: Path) -> dict:
    atoms, header = read_pdb(in_pdb)
    p = spec.protein
    report: dict = {"input_atoms": len(atoms)}

    report["opm_half_thickness_A"] = opm_half_thickness(header)

    dropped: dict[str, int] = {}

    def drop(reason: str) -> None:
        dropped[reason] = dropped.get(reason, 0) + 1

    keep_het = {r.upper() for r in p["keep_hetatm"]}
    chains = set(p["chains"]) if p["chains"] else None
    kept: list[Atom] = []
    for a in atoms:
        if a.resname in DUMMY_RESNAMES:
            drop("opm_dummy")
        elif a.resname in WATER_RESNAMES:
            drop("water")
        elif chains is not None and a.chain not in chains:
            drop("chain_not_selected")
        elif a.is_hydrogen:
            drop("hydrogen")
        elif a.record == "HETATM" and a.resname.upper() not in keep_het:
            drop(f"hetatm:{a.resname}")
        else:
            kept.append(a)
    kept = _pick_altlocs(kept)

    if not kept:
        raise PrepareError("no atoms left after cleaning; check protein.chains and the input file")
    present_chains = list(OrderedDict.fromkeys(a.chain for a in kept))
    if chains is not None:
        missing = sorted(chains - set(present_chains))
        if missing:
            raise PrepareError(f"chains {missing} not found in the input (found {present_chains})")

    # martinize2 needs chain IDs to separate molecules
    for a in kept:
        if not a.chain:
            a.chain = "A"

    gaps = _find_gaps(kept)
    write_pdb(out_pdb, kept, header=[h for h in header if h.startswith("CRYST1")])
    report.update({
        "output_atoms": len(kept),
        "chains": present_chains,
        "dropped": dropped,
        "gaps": gaps,
        "kept_hetatm": sorted({a.resname for a in kept if a.record == "HETATM"}),
    })
    return report
