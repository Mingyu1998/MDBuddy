"""Turn COBY's output into a clean GROMACS input folder (CHARMM-GUI-style naming).

gromacs/
  step5_input.gro / .pdb     assembled system
  topol.top, toppar/         topology with pinned force-field files
  index.ndx                  Protein, MEMBRANE, SOLVENT, SOLU_MEMB, BB, System
  step6.*_*.mdp, step7_production.mdp
  run.sh                     grompp + mdrun for every stage
"""

from __future__ import annotations

import shutil
from pathlib import Path

from . import lipids, mdp, topology
from .structure import read_gro

# Names COBY may write that differ from the moleculetype names in Martini 2 files
MOLECULE_ALIASES = {"martini2": {"NA+": "NA", "CL-": "CL"}}


class AssembleError(RuntimeError):
    pass


def alias_renames(aliases: dict, names, gmx_dir: Path, includes: list[str]) -> dict[str, str]:
    """Which of COBY's molecule names have to be renamed to match the parameter files.

    The two Martini 2 ion files in circulation disagree: COBY 1.0.15 and the copy MDBuddy
    pins (COBY's M2_BacterialMembrane tutorial) call the ions NA/CL, while the cgmartini.nl
    martini_v2.0_ions.itp calls them NA+/CL-. So only rename a name the parameter files do
    not define into one they do; with either file the other spelling is left alone.
    """
    candidates = {n: aliases[n] for n in set(names) if n in aliases}
    if not candidates:
        return {}
    info = topology.TopologyInfo()
    for inc in includes:
        if (gmx_dir / inc).exists():
            topology.parse(gmx_dir / inc, info=info)
    return {n: a for n, a in candidates.items()
            if n not in info.moleculetypes and a in info.moleculetypes}


def classify(spec, protein_moltypes: list[str], name: str) -> str:
    if name in protein_moltypes:
        return "Protein"
    if name in spec.lipid_names():
        return "MEMBRANE"
    return "SOLVENT"


def write_index(path: Path, spec, molecules, info: topology.TopologyInfo,
                protein_moltypes: list[str]) -> dict[str, int]:
    groups: dict[str, list[int]] = {"System": [], "Protein": [], "MEMBRANE": [], "SOLVENT": [], "BB": []}
    atom = 0
    for name, count in molecules:
        mt = info.moleculetypes[name]
        cls = classify(spec, protein_moltypes, name)
        for _ in range(count):
            for aname in mt.atom_names:
                atom += 1
                groups["System"].append(atom)
                groups[cls].append(atom)
                if cls == "Protein" and aname == "BB":
                    groups["BB"].append(atom)
    groups["SOLU_MEMB"] = sorted(groups["Protein"] + groups["MEMBRANE"])
    order = ["System", "Protein", "MEMBRANE", "SOLU_MEMB", "SOLVENT", "BB"]
    lines = []
    for g in order:
        lines.append(f"[ {g} ]")
        ids = groups[g]
        for i in range(0, len(ids), 15):
            lines.append(" ".join(f"{n:>6}" for n in ids[i:i + 15]))
        lines.append("")
    path.write_text("\n".join(lines))
    return {g: len(groups[g]) for g in order}


def fix_coby_names(gro_path: Path, pdb_path: Path, molecules, info: topology.TopologyInfo,
                   renames: dict[str, dict[int, tuple[str, str]]]) -> dict[str, int]:
    """Rename atoms COBY names differently from the itp (known COBY definition typos).

    Only exact, pre-verified renames from data/lipids.json are applied: same bead count, same
    bead order, an atom is renamed only if it carries the known wrong name. Returns counts.
    """
    if not renames:
        return {}
    lines = gro_path.read_text().splitlines()
    natoms = int(lines[1])
    new_names: dict[int, str] = {}
    counts: dict[str, int] = {}
    i = 0
    for name, count in molecules:
        n = info.moleculetypes[name].natoms if name in info.moleculetypes else 0
        table = renames.get(name)
        for _ in range(count):
            if table:
                for pos, (wrong, right) in table.items():
                    if i + pos < natoms and lines[2 + i + pos][10:15].strip() == wrong:
                        new_names[i + pos] = right
                        counts[name] = counts.get(name, 0) + 1
            i += n
    if not new_names:
        return {}
    for idx, right in new_names.items():
        line = lines[2 + idx]
        lines[2 + idx] = line[:10] + f"{right:>5}" + line[15:]
    gro_path.write_text("\n".join(lines) + "\n")
    if pdb_path.exists():
        out, k = [], -1
        for line in pdb_path.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                k += 1
                if k in new_names:
                    right = new_names[k]
                    field = f" {right:<3}" if len(right) < 4 else right[:4]
                    line = line[:12] + field + line[16:]
            out.append(line)
        pdb_path.write_text("\n".join(out) + "\n")
    return counts


def validate(gro_path: Path, molecules, info: topology.TopologyInfo) -> dict:
    """Pure-Python consistency checks (no GROMACS needed). Returns {'errors', 'warnings', ...}."""
    errors, warnings = [], []
    if info.missing_includes:
        errors.append(f"missing #include files: {info.missing_includes}")
    if info.defaults_count != 1:
        errors.append(f"expected exactly one [ defaults ] section, found {info.defaults_count}")
    for name, first, second in info.duplicates:
        errors.append(f"moleculetype {name} defined twice ({first} and {second})")
    undefined = [n for n, _ in molecules if n not in info.moleculetypes]
    if undefined:
        errors.append(f"molecules not defined in any included itp: {sorted(set(undefined))}")
        return {"errors": errors, "warnings": warnings}

    _, atoms, box = read_gro(gro_path)
    expected = sum(info.moleculetypes[n].natoms * c for n, c in molecules)
    if expected != len(atoms):
        errors.append(f"topology has {expected} atoms but {gro_path.name} has {len(atoms)}")
        return {"errors": errors, "warnings": warnings}

    mismatches = 0
    first_mismatch = None
    i = 0
    for name, count in molecules:
        mt = info.moleculetypes[name]
        for _ in range(count):
            for aname in mt.atom_names:
                if atoms[i].name != aname[:5]:
                    mismatches += 1
                    if first_mismatch is None:
                        first_mismatch = (i + 1, name, aname, atoms[i].name)
                i += 1
    if mismatches:
        idx, name, want, got = first_mismatch
        errors.append(f"{mismatches} atom-name mismatches between topology and coordinates; first at atom "
                      f"{idx} ({name}: topology '{want}', coordinates '{got}')")

    charge = sum(info.moleculetypes[n].charge * c for n, c in molecules)
    if abs(charge) > 1e-3:
        warnings.append(f"system net charge is {charge:+.3f} e (expected 0)")
    return {"errors": errors, "warnings": warnings, "atoms": len(atoms), "net_charge": round(charge, 4),
            "box_nm": box[:3]}


RUN_SH = """#!/bin/bash
# Generated by MDBuddy for '{name}'. Runs minimisation, equilibration and production.
# Usage: ./run.sh                     (all stages)
#        ./run.sh step6.3             (from a stage to the end)
#        ./run.sh step6.0 step6.7     (from one stage to another, both included)
#        ./run.sh '' step6.7          (up to a stage: everything but production)
# Stage names may be abbreviated to a unique prefix, e.g. step7.
# Environment: GMX (default: gmx), MDRUN_FLAGS (e.g. "-nt 16 -gpu_id 0"), MAXWARN (default: {maxwarn})
set -euo pipefail
cd "$(dirname "$0")"
GMX=${{GMX:-gmx}}
MDRUN_FLAGS=${{MDRUN_FLAGS:-}}
MAXWARN=${{MAXWARN:-{maxwarn}}}
START=${{1:-}}
END=${{2:-}}

stages=({stages})

index_of() {{  # first stage whose name starts with $1, or an error
    local pat=$1 i
    for i in "${{!stages[@]}}"; do
        if [[ "${{stages[$i]}}" == "$pat"* ]]; then printf '%s\n' "$i"; return 0; fi
    done
    echo "run.sh: no stage matches '$pat'. Stages: ${{stages[*]}}" >&2
    return 1
}}

from=0
to=$(( ${{#stages[@]}} - 1 ))
if [ -n "$START" ]; then from=$(index_of "$START"); fi
if [ -n "$END" ]; then to=$(index_of "$END"); fi
if [ "$from" -gt "$to" ]; then
    echo "run.sh: ${{stages[$from]}} comes after ${{stages[$to]}}, nothing to run" >&2
    exit 1
fi

prev=step5_input
if [ "$from" -gt 0 ]; then prev=${{stages[$(( from - 1 ))]}}; fi
for (( i = from; i <= to; i++ )); do
    s=${{stages[$i]}}
    echo ">>> $s"
    $GMX grompp -f "$s.mdp" -c "$prev.gro" -r step5_input.gro -p topol.top -n index.ndx \\
        -o "$s.tpr" -po "$s.mdout.mdp" -maxwarn "$MAXWARN"
    $GMX mdrun -deffnm "$s" $MDRUN_FLAGS
    prev=$s
done
"""

README = """# {name}: GROMACS input generated by MDBuddy

| Item | Value |
|---|---|
| Force field | {ff} ({model} protein model) |
| Box | {box} nm |
| Beads | {atoms} |
| Net charge | {charge} e |

Molecules:

{molecules}

## Run

```bash
./run.sh                      # everything, in order
MDRUN_FLAGS="-nt 8" ./run.sh  # pass flags to mdrun
./run.sh step7                # only production (after equilibration)
./run.sh step6.0 step6.7      # minimisation + equilibration, stop before production
./run.sh '' step6.7           # same, written as "up to step6.7"
```

| File | Purpose |
|---|---|
| step5_input.gro | assembled system (start coordinates and restraint reference) |
| topol.top, toppar/ | topology; every parameter file is listed with its SHA-256 in ../mdbuddy.lock.json |
| index.ndx | groups Protein, MEMBRANE, SOLVENT, SOLU_MEMB, BB, System |
{stage_rows}

Rebuild from scratch: `mdbuddy build ../spec.resolved.yaml`
"""


def assemble(spec, coby_dir: Path, gmx_dir: Path, ff_includes: list[str], protein_itps: list[str],
             protein_moltypes: list[str]) -> dict:
    gmx_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(coby_dir / "system.gro", gmx_dir / "step5_input.gro")
    if (coby_dir / "system.pdb").exists():
        shutil.copyfile(coby_dir / "system.pdb", gmx_dir / "step5_input.pdb")

    family = "martini3" if spec.cg["forcefield"] == "martini3" else "martini2"
    includes = [f"toppar/{f}" for f in ff_includes] + [f"toppar/{f}" for f in protein_itps]
    raw = topology.read_molecules(coby_dir / "coby_topol.top")
    ion_renames = alias_renames(MOLECULE_ALIASES.get(family, {}), [n for n, _ in raw],
                                gmx_dir, includes)
    molecules = topology.merge_consecutive((ion_renames.get(n, n), c) for n, c in raw)

    defines = ["GO_VIRT"] if spec.cg["model"] == "go" else []
    top = gmx_dir / "topol.top"
    topology.write_top(top, includes, molecules, spec.name, defines=defines,
                       header=[f"Generated by MDBuddy for '{spec.name}'",
                               "Do not edit by hand: change the spec and rebuild."])

    info = topology.parse(top)
    renames = lipids.check(spec.lipid_names(), lipids.library_for(spec)).renames
    renamed = fix_coby_names(gmx_dir / "step5_input.gro", gmx_dir / "step5_input.pdb", molecules, info, renames)
    report = validate(gmx_dir / "step5_input.gro", molecules, info)
    report["coby_renames"] = renamed
    if report["errors"]:
        raise AssembleError("system failed validation:\n  - " + "\n  - ".join(report["errors"]))

    index_counts = write_index(gmx_dir / "index.ndx", spec, molecules, info, protein_moltypes)
    for key in ("coupling_groups", "comm_groups"):
        for g in spec.sim[key]:
            if g != "System" and index_counts.get(g, 0) == 0:
                raise AssembleError(f"index group {g} is empty; change simulation.{key}")

    posres = spec.cg["position_restraints"] != "none" and bool(protein_moltypes)
    # -DFLEXIBLE only helps where an itp defines an #ifdef FLEXIBLE block (Martini 3 proteins
    # and sterols do, Martini 2 does not). grompp warns about an unused -D macro, so writing it
    # blindly fails the -maxwarn 0 check.
    flexible = any("FLEXIBLE" in (gmx_dir / inc).read_text() for inc in includes
                   if (gmx_dir / inc).exists())
    stage_list = mdp.stages(spec)
    for st in stage_list:
        (gmx_dir / f"{st.name}.mdp").write_text(mdp.render(st, spec, posres, flexible))

    run_sh = gmx_dir / "run.sh"
    run_sh.write_text(RUN_SH.format(name=spec.name, maxwarn=spec.sim["gromacs_maxwarn"],
                                    stages=" ".join(st.name for st in stage_list)))
    run_sh.chmod(0o755)

    mol_lines = "\n".join(f"- {n}: {c}" for n, c in molecules)
    stage_rows = "\n".join(f"| {st.name}.mdp | {mdp.title(st, flexible)}"
                           + (f", backbone restraints {st.posres_fc:g} kJ/mol/nm²" if posres and st.posres_fc else "")
                           + f", {st.nsteps} steps |" for st in stage_list)
    (gmx_dir / "README.md").write_text(README.format(
        name=spec.name, ff=spec.cg["forcefield"], model=spec.cg["model"],
        box=" x ".join(f"{v:g}" for v in report["box_nm"][:3]), atoms=report["atoms"],
        charge=f"{report['net_charge']:+g}", molecules=mol_lines, stage_rows=stage_rows))

    outputs = [gmx_dir / "step5_input.gro", top, gmx_dir / "index.ndx", run_sh,
               *[gmx_dir / f"{st.name}.mdp" for st in stage_list]]
    return {"molecules": molecules, "validation": report, "index": index_counts,
            "stages": [st.name for st in stage_list], "outputs": outputs}
