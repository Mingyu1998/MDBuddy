"""The build pipeline: fetch -> clean -> forcefield -> martinize -> membrane -> gromacs -> check.

Each stage writes into its own folder under the output directory and records
what it did in ``mdbuddy.lock.json``. Stages can be re-run from any point with
``--from``; earlier results are read back from disk.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Callable, Optional

from . import assemble, cgprotein, citations, doctor, fetch, forcefields, lipids, martinize, membrane, prepare, topology
from .provenance import Lock, sha256
from .structure import read_pdb
from . import tools
from .tools import find_gmx

STAGES = ("fetch", "clean", "orient", "forcefield", "martinize", "membrane", "gromacs", "check")

DIRS = {
    "fetch": "00_input",
    "clean": "01_clean",
    "orient": "01_orient",
    "martinize": "02_cg",
    "membrane": "03_membrane",
    "gromacs": "gromacs",
}


class BuildError(RuntimeError):
    pass


def _state_path(out: Path) -> Path:
    return out / ".mdbuddy_state.json"


def _load_state(out: Path) -> dict:
    p = _state_path(out)
    return json.loads(p.read_text()) if p.exists() else {}


def _save_state(out: Path, state: dict) -> None:
    _state_path(out).write_text(json.dumps(state, indent=2, default=str) + "\n")


def user_lipid_definitions(spec) -> dict[str, str]:
    """Moleculetypes defined in forcefield_files.extra_itps (name -> file), for custom lipids."""
    found: dict[str, str] = {}
    for extra in spec.ff_files["extra_itps"]:
        path = spec.resolve_path(extra)
        if path.exists():
            for name in topology.parse(path).moleculetypes:
                found.setdefault(name, path.name)
    return found


def coby_definitions(spec):
    """Molecules the coby section gives COBY a structure for: names from molecule_import, or True
    when an import_library file is given (its contents are only known to COBY)."""
    kw = spec.data["coby"]["kwargs"]
    if kw.get("import_library"):
        return True
    imports = kw.get("molecule_import") or []
    names = set()
    for m in [imports] if isinstance(imports, str) else imports:
        for tok in str(m).split():
            key, _, value = tok.partition(":")
            if key in ("moleculetype", "moleculetypes", "name", "refname") and value:
                names.add(value)
    return names


def build(spec, from_stage: str = "fetch", until_stage: str = "check", grompp_check: bool = True,
          log: Callable[[str], None] = print, force: bool = False) -> dict:
    if from_stage not in STAGES or until_stage not in STAGES:
        raise BuildError(f"stages are: {', '.join(STAGES)}")
    start, stop = STAGES.index(from_stage), STAGES.index(until_stage)
    out = spec.output_dir()

    if start == 0 and out.exists() and any(out.iterdir()):
        if not force:
            raise BuildError(f"output directory {out} already exists and is not empty; "
                             "use --force to overwrite it or --from <stage> to resume")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "spec.resolved.yaml").write_text(spec.to_yaml())

    lock_path = out / "mdbuddy.lock.json"
    lock = Lock(lock_path)
    if start > 0 and lock_path.exists():
        lock.data = json.loads(lock_path.read_text())
        lock.data["stages"] = [s for s in lock.data["stages"] if STAGES.index(s["stage"]) < start]
    state = _load_state(out) if start > 0 else {}
    toppar = out / DIRS["gromacs"] / "toppar"

    lock.data["tools"] = {c.name: c.version for c in doctor.run_checks() if c.version}
    lock.data["spec"] = spec.data
    lock.save()

    # Lipid names are checked before anything runs: unknown names, cross-version name clashes,
    # and lipids COBY cannot build (see lipids.py)
    lib = lipids.library_for(spec)
    lipid_check = lipids.check(spec.lipid_names(), lib, custom=user_lipid_definitions(spec),
                               coby_definitions=coby_definitions(spec))
    lock.data["lipids"] = {"library": lib, "warnings": lipid_check.warnings,
                           "coby_renames": {k: {str(i): list(v) for i, v in t.items()}
                                            for k, t in lipid_check.renames.items()}}
    lock.save()
    if spec.protein["orientation"] == "self_assembly":
        sa_check = lipids.self_assembly_check(spec.leaflets()["upper"], spec.leaflets()["lower"], lib,
                                              float(spec.sim["temperature"]))
        lipid_check.errors += sa_check.errors
        lipid_check.warnings += sa_check.warnings
        lock.data["lipids"]["warnings"] = lipid_check.warnings
        lock.save()
    for w in lipid_check.warnings:
        log(f"    LIPID NOTE {w}")
    if lipid_check.errors:
        raise BuildError("membrane lipids:\n  - " + "\n  - ".join(lipid_check.errors))

    for i, stage in enumerate(STAGES):
        if i < start or i > stop:
            continue
        log(f"==> [{stage}]")
        entry = lock.stage(stage)

        if stage in ("fetch", "clean", "orient") and not spec.has_atomistic_protein:
            log("    skipped: no atomistic protein (protein.source: none)")
            lock.finish(entry)
            _save_state(out, state)
            continue

        if stage == "fetch":
            origin = fetch.fetch_structure(spec, out / DIRS["fetch"])
            lock.data["inputs"] = [origin]
            lock.finish(entry, [out / DIRS["fetch"] / "input.pdb"], base=out)
            log(f"    {origin.get('url') or origin.get('path')}")

        elif stage == "clean":
            d = out / DIRS["clean"]
            d.mkdir(parents=True, exist_ok=True)
            report = prepare.prepare(spec, out / DIRS["fetch"] / "input.pdb", d / "protein_clean.pdb")
            incomplete = report["incomplete_residues"]
            if incomplete:
                summary = ", ".join(f"{r['chain']}:{r['resname']}{r['resid']}" for r in incomplete[:12])
                summary += f" (+{len(incomplete) - 12} more)" if len(incomplete) > 12 else ""
                log(f"    {len(incomplete)} residue(s) with missing heavy atoms: {summary}")
                mode = spec.protein["repair"]
                if mode == "off":
                    raise BuildError("residues with missing heavy atoms cannot be coarse-grained (martinize2 "
                                     "writes nan coordinates for beads without atoms). Set protein.repair: auto "
                                     "(needs `pip install pdbfixer`) or complete the structure yourself.")
                if not tools.has_module("pdbfixer"):
                    raise BuildError("missing heavy atoms must be rebuilt before martinize2, and PDBFixer is not "
                                     "installed: `pip install pdbfixer` (or conda install -c conda-forge pdbfixer), "
                                     "then rerun with --from clean")
                raw = d / "protein_clean_unrepaired.pdb"
                (d / "protein_clean.pdb").rename(raw)
                prepare.repair_with_pdbfixer(lock, entry, raw, d / "protein_clean.pdb", spec.seed)
                still = prepare.find_incomplete(read_pdb(d / "protein_clean.pdb")[0])
                if still:
                    raise BuildError(f"PDBFixer left {len(still)} residue(s) incomplete; see 01_clean/repair_pdbfixer.log")
                report["repaired_with"] = "pdbfixer"
                state["repaired_with"] = "pdbfixer"
                report["output_atoms"] = len(read_pdb(d / "protein_clean.pdb")[0])
                log(f"    rebuilt missing atoms with PDBFixer -> 01_clean/protein_clean.pdb")
            (d / "clean_report.json").write_text(json.dumps(report, indent=2) + "\n")
            state["chains"] = report["chains"]
            state["opm_half_thickness_A"] = report["opm_half_thickness_A"]
            entry["notes"].append(report)
            lock.finish(entry, [d / "protein_clean.pdb"], base=out)
            log(f"    chains {report['chains']}, {report['output_atoms']} atoms, "
                f"{len(report['gaps'])} gap(s)")
            for gap in report["gaps"]:
                log(f"    WARNING gap in chain {gap['chain']} between residues {gap['after_resid']} and "
                    f"{gap['before_resid']} (CA-CA {gap['ca_distance']} A)")
            if spec.protein["orientation"] == "opm" and report["opm_half_thickness_A"] is None:
                log("    WARNING orientation is 'opm' but the file has no OPM bilayer REMARK; "
                    "check that the structure is membrane-oriented")

        elif stage == "orient":
            info = place_protein(spec, out / DIRS["clean"] / "protein_clean.pdb", out / DIRS["orient"], log)
            state["membrane_center_z_A"] = info["membrane_center_z_A"]
            state["orientation"] = info["record"]
            entry["notes"].append(info["record"])
            lock.finish(entry, [out / DIRS["orient"] / "protein_oriented.pdb"], base=out)

        elif stage == "forcefield":
            records = forcefields.install_forcefield(spec, toppar)
            if spec.cg["model"] == "go":
                main = toppar / records[0]["file"]
                forcefields.patch_for_go(main)
                records[0]["patched_for_go"] = True
                records[0]["sha256_after_patch"] = sha256(main)
            lock.data["forcefield_files"] = records
            forcefields.check_molecules(spec, toppar, [r["file"] for r in records])
            state["ff_includes"] = [r["file"] for r in records]
            lock.finish(entry, [toppar / r["file"] for r in records], base=out)
            log(f"    {len(records)} parameter files -> {toppar.relative_to(out)}")

        elif stage == "martinize":
            outputs = []
            state["protein_molecules"], state["protein_itps"] = [], []
            if spec.has_atomistic_protein:
                res = martinize.run(spec, lock, entry, out / DIRS["orient"] / "protein_oriented.pdb",
                                    state["chains"], out / DIRS["martinize"], toppar)
                state["protein_molecules"] = res["molecules"]
                state["protein_itps"] = res["itps"]
                entry["notes"].append({"dssp": res["dssp"], "warnings": res["warnings"]})
                outputs += res["outputs"]
                log(f"    molecules {res['molecules']}")
                for w in res["warnings"][:10]:
                    log(f"    {w.strip()}")
            state["cg_molecules"], state["cg_itps"], state["cg_entries"] = [], [], []
            if spec.cg_proteins:
                try:
                    cg = cgprotein.prepare(spec, out / DIRS["martinize"], toppar,
                                           state["ff_includes"] + state["protein_itps"])
                except cgprotein.CGProteinError as exc:
                    raise BuildError(str(exc)) from exc
                state["cg_molecules"], state["cg_itps"] = cg["molecules"], cg["itps"]
                state["cg_entries"] = [{**{k: v for k, v in e.items() if k not in ("xyz", "box")},
                                        "structure": e["structure"].name} for e in cg["entries"]]
                entry["notes"].append({"cg_proteins": cg["records"]})
                outputs += [e["structure"] for e in cg["entries"]] + [toppar / f for f in cg["itps"]]
                for r in cg["records"]:
                    log(f"    pre-coarse-grained {r['copied_to']}: {r['beads']} beads, {r['molecules']}, "
                        f"placement {r['placement']}")
                    if not r["posres_in_itps"] and spec.cg["position_restraints"] != "none":
                        log(f"    NOTE its itps have no #ifdef POSRES block, so it is not restrained during "
                            "equilibration (add one if the arrangement must be held)")
            if not state["protein_molecules"] and not state["cg_molecules"]:
                log("    no protein: membrane-only system")
            lock.finish(entry, outputs, base=out)

        elif stage == "membrane":
            d = out / DIRS["membrane"]
            itp_files = [toppar / f for f in state["ff_includes"] + state["protein_itps"] + state.get("cg_itps", [])]
            molecules = [tuple(m) for m in state["protein_molecules"]]
            cg_pdb = out / DIRS["martinize"] / "cg.pdb" if spec.has_atomistic_protein else None
            if spec.protein["orientation"] == "self_assembly":
                script, info = membrane.write_selfassembly_script(spec, d, cg_pdb, molecules, itp_files)
                log(f"    self-assembly soup: {info['lipid_counts']} in a {info['box_nm']} nm box "
                    f"(~{info['water_target']} water beads requested)")
                for w in info["warnings"]:
                    log(f"    WARNING {w}")
            else:
                cg_entries = cgprotein.entries_from_state(spec, out / DIRS["martinize"], state.get("cg_entries", []))
                try:
                    script, info = membrane.write_script(spec, d, cg_pdb, molecules, itp_files,
                                                         membrane_center_z_A=state.get("membrane_center_z_A", 0.0),
                                                         cg_entries=cg_entries)
                except cgprotein.CGProteinError as exc:
                    raise BuildError(str(exc)) from exc
                log(f"    box {info['box_nm']} nm, {info['membranes']} membrane(s), {info['proteins']} protein "
                    "structure(s)" + (f", protein z-offset {info['protein_cz_nm']} nm"
                                      if info["protein_cz_nm"] is not None else ""))
                for n in info["notes"]:
                    log(f"    NOTE {n}")
                for r in info.get("user_files", []):
                    log(f"    user file {r['source']} -> 03_membrane/{r['copied_to']}")
            entry["notes"].append({k: v for k, v in info.items() if k != "coby_kwargs"})
            res = membrane.run(spec, lock, entry, d, script)
            lock.finish(entry, [script, *res["outputs"]], base=out)

        elif stage == "gromacs":
            protein_moltypes = [m[0] for m in state["protein_molecules"] + state.get("cg_molecules", [])]
            res = assemble.assemble(spec, out / DIRS["membrane"], out / DIRS["gromacs"],
                                    state["ff_includes"], state["protein_itps"] + state.get("cg_itps", []),
                                    protein_moltypes)
            state["system_molecules"] = res["molecules"]
            entry["notes"].append({"validation": res["validation"], "index": res["index"]})
            used = frozenset({state["repaired_with"]} if state.get("repaired_with") else ())
            (out / "CITATIONS.md").write_text(citations.markdown(spec, used))
            lock.finish(entry, res["outputs"], base=out)
            v = res["validation"]
            log(f"    {v['atoms']} beads, net charge {v['net_charge']:+g} e")
            for lipid, n in (v.get("coby_renames") or {}).items():
                log(f"    renamed {n} atoms of {lipid} to match its itp (known COBY naming issue)")
            for w in v["warnings"]:
                log(f"    WARNING {w}")

        elif stage == "check":
            result = grompp_dry_run(out / DIRS["gromacs"], lock, entry) if grompp_check else None
            entry["notes"].append({"grompp": result})
            lock.finish(entry)
            if result is None:
                log("    skipped grompp check")
            else:
                log(f"    grompp {result}")

        _save_state(out, state)

    return {"output_dir": out, "lock": lock_path, "state": state}


def place_protein(spec, clean_pdb: Path, work: Path, log: Callable[[str], None]) -> dict:
    """The orient stage: write 01_orient/protein_oriented.pdb in the frame the membrane builder expects."""
    from . import orient, viewer
    from .structure import write_pdb

    work.mkdir(parents=True, exist_ok=True)
    out_pdb = work / "protein_oriented.pdb"
    atoms, header = read_pdb(clean_pdb)
    p = spec.protein
    mode = p["orientation"]
    centre_z = 0.0
    if mode in ("opm", "as_is"):
        shutil.copyfile(clean_pdb, out_pdb)
        centre_z = float(p["membrane_center_z"])
        record = {"method": mode, "membrane_center_z_A": centre_z}
        log(f"    {mode}: using the input frame (membrane centre at z = {centre_z:g} A)")
    elif mode == "manual":
        placed, record = orient.manual(atoms, p["manual"])
        write_pdb(out_pdb, placed)
        log("    manual placement applied")
    else:  # hydrophobic, or the initial alignment for self_assembly
        if mode == "self_assembly" and spec.data["self_assembly"]["align"] == "as_is":
            shutil.copyfile(clean_pdb, out_pdb)
            record = {"method": "as_is (self-assembly start)"}
            log("    self-assembly: protein kept in the input frame")
        else:
            placed, record = orient.hydrophobic(atoms, topology=p["topology"], pore_correction=p["pore_correction"])
            pl = record["placement"]
            log(f"    hydrophobic placement: {pl['verdict']}")
            log(f"    thickness {pl['thickness_A']:g} A, transfer energy {pl['energy_kcal']:g} kcal/mol, "
                f"orientation spread {pl['uncertainty_deg']:g} deg, TM segments {', '.join(pl['tm_segments']) or '-'}")
            log(f"    cytoplasmic side (-z) by positive-inside rule: {pl['kr_inside']} vs {pl['kr_outside']} Lys/Arg")
            for n in pl["notes"]:
                log(f"    NOTE {n}")
            if pl["verdict"].startswith("not a"):
                if mode == "self_assembly":
                    shutil.copyfile(clean_pdb, out_pdb)
                    record = {"method": "as_is (self-assembly start; no hydrophobic belt found)"}
                else:
                    raise BuildError("no clear transmembrane placement was found. Check the structure, place it "
                                     "yourself (orientation: manual; `mdbuddy orient view` helps) or let the "
                                     "lipids find it (orientation: self_assembly).")
            else:
                write_pdb(out_pdb, placed)
    normal, centre = (orient.membrane_in_input_frame(record) if "rotation" in record else (None, None))
    if normal is not None:
        record["membrane_normal_input_frame"] = normal
        record["membrane_centre_input_frame_A"] = centre
    (work / "orient_report.json").write_text(json.dumps(record, indent=2) + "\n")
    if mode in ("manual", "hydrophobic", "self_assembly"):
        viewer.write_html(work / "orientation_view.html", read_pdb(clean_pdb)[0], record,
                          thickness_A=(record.get("placement") or {}).get("thickness_A", 30.0), title=spec.name)
        log("    check it: open 01_orient/orientation_view.html in a browser (adjust and copy the spec lines)")
    return {"membrane_center_z_A": centre_z, "record": record}


def grompp_check_stages(gmx_dir: Path) -> list[str]:
    """One mdp of each kind: the first EM, the first MD stage (gen-vel), and the final run."""
    names = sorted(p.stem for p in gmx_dir.glob("step*.mdp") if not p.stem.endswith(".mdout"))
    md = [n for n in names if not n.endswith("minimization")]
    return list(dict.fromkeys([names[0], md[0], md[-1]])) if md else names[:1]


def grompp_dry_run(gmx_dir: Path, lock: Lock, entry: dict) -> Optional[str]:
    """Run grompp (-maxwarn 0) on representative stages to catch topology/mdp errors."""
    gmx = find_gmx()
    if not gmx:
        return "not run: GROMACS not found (the inputs are still complete)"
    check_dir = gmx_dir / ".check"
    check_dir.mkdir(exist_ok=True)
    stages = grompp_check_stages(gmx_dir)
    for stage in stages:
        cmd = [gmx, "grompp", "-f", f"../{stage}.mdp", "-c", "../step5_input.gro",
               "-r", "../step5_input.gro", "-p", "../topol.top", "-n", "../index.ndx",
               "-o", "check.tpr", "-po", "check.mdout.mdp", "-maxwarn", "0"]
        lock.run(entry, cmd, cwd=check_dir, log_name=f"grompp_{stage}.log")
        for f in ("check.tpr", "check.mdout.mdp"):
            (check_dir / f).unlink(missing_ok=True)
    return f"passed with 0 warnings for {', '.join(stages)}"
