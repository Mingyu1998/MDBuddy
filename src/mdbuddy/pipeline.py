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

from . import assemble, citations, doctor, fetch, forcefields, martinize, membrane, prepare
from .provenance import Lock, sha256
from .structure import read_pdb
from .tools import find_gmx

STAGES = ("fetch", "clean", "forcefield", "martinize", "membrane", "gromacs", "check")

DIRS = {
    "fetch": "00_input",
    "clean": "01_clean",
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

    for i, stage in enumerate(STAGES):
        if i < start or i > stop:
            continue
        log(f"==> [{stage}]")
        entry = lock.stage(stage)

        if stage == "fetch":
            origin = fetch.fetch_structure(spec, out / DIRS["fetch"])
            lock.data["inputs"] = [origin]
            lock.finish(entry, [out / DIRS["fetch"] / "input.pdb"], base=out)
            log(f"    {origin.get('url') or origin.get('path')}")

        elif stage == "clean":
            d = out / DIRS["clean"]
            d.mkdir(parents=True, exist_ok=True)
            report = prepare.prepare(spec, out / DIRS["fetch"] / "input.pdb", d / "protein_clean.pdb")
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
            res = martinize.run(spec, lock, entry, out / DIRS["clean"] / "protein_clean.pdb",
                                state["chains"], out / DIRS["martinize"], toppar)
            state["protein_molecules"] = res["molecules"]
            state["protein_itps"] = res["itps"]
            entry["notes"].append({"dssp": res["dssp"], "warnings": res["warnings"]})
            lock.finish(entry, res["outputs"], base=out)
            log(f"    molecules {res['molecules']}")
            for w in res["warnings"][:10]:
                log(f"    {w.strip()}")

        elif stage == "membrane":
            d = out / DIRS["membrane"]
            itp_files = [toppar / f for f in state["ff_includes"] + state["protein_itps"]]
            script, info = membrane.write_script(spec, d, out / DIRS["martinize"] / "cg.pdb",
                                                 [tuple(m) for m in state["protein_molecules"]], itp_files)
            entry["notes"].append({"box_nm": info["box_nm"], "protein_cz_nm": info["protein_cz_nm"]})
            res = membrane.run(spec, lock, entry, d, script)
            lock.finish(entry, [script, *res["outputs"]], base=out)
            log(f"    box {info['box_nm']} nm, protein z-offset {info['protein_cz_nm']} nm")

        elif stage == "gromacs":
            protein_moltypes = [m[0] for m in state["protein_molecules"]]
            res = assemble.assemble(spec, out / DIRS["membrane"], out / DIRS["gromacs"],
                                    state["ff_includes"], state["protein_itps"], protein_moltypes)
            state["system_molecules"] = res["molecules"]
            entry["notes"].append({"validation": res["validation"], "index": res["index"]})
            (out / "CITATIONS.md").write_text(citations.markdown(spec))
            lock.finish(entry, res["outputs"], base=out)
            v = res["validation"]
            log(f"    {v['atoms']} beads, net charge {v['net_charge']:+g} e")
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


# One mdp of each kind: soft-core EM, first NPT stage (gen-vel, C-rescale), production
GROMPP_CHECK_STAGES = ("step6.0_minimization", "step6.2_equilibration", "step7_production")


def grompp_dry_run(gmx_dir: Path, lock: Lock, entry: dict) -> Optional[str]:
    """Run grompp (-maxwarn 0) on representative stages to catch topology/mdp errors."""
    gmx = find_gmx()
    if not gmx:
        return "not run: GROMACS not found (the inputs are still complete)"
    check_dir = gmx_dir / ".check"
    check_dir.mkdir(exist_ok=True)
    for stage in GROMPP_CHECK_STAGES:
        cmd = [gmx, "grompp", "-f", f"../{stage}.mdp", "-c", "../step5_input.gro",
               "-r", "../step5_input.gro", "-p", "../topol.top", "-n", "../index.ndx",
               "-o", "check.tpr", "-po", "check.mdout.mdp", "-maxwarn", "0"]
        lock.run(entry, cmd, cwd=check_dir, log_name=f"grompp_{stage}.log")
        for f in ("check.tpr", "check.mdout.mdp"):
            (check_dir / f).unlink(missing_ok=True)
    return f"passed with 0 warnings for {', '.join(GROMPP_CHECK_STAGES)}"
