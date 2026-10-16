# MDBuddy

**An agent skill and command-line tool for building Martini coarse-grained
membrane-protein systems for GROMACS. A reproducible alternative to
CHARMM-GUI Martini Maker.**

Give it a PDB ID and a description of the membrane. It returns a GROMACS folder
ready to run: coordinates, topology, force-field files, index groups, a
minimisation → equilibration → production `.mdp` ladder and a `run.sh`. Every
download, parameter file and command is recorded with checksums, so the same spec
rebuilds the same system.

```
PDB / OPM ──► clean ──► martinize2 ──► COBY (bilayer + insertion + solvent + ions) ──► GROMACS inputs ──► grompp check
```

Scope of v0.1: **membrane proteins from existing structures (OPM / RCSB / file),
Martini 3 and Martini 2.2, GROMACS.** Structure prediction (co-folding) input,
soluble proteins and other engines are planned.

## Two ways to use it

### 1. As a Claude Code skill (recommended)

The skill (`skills/mdbuddy/SKILL.md`) has the agent check your environment,
ask for the decisions that matter (force field, protein model, lipids per leaflet,
salt, temperature) with sensible defaults, write the spec, run the build stage by
stage, explain warnings such as missing loops, and report what was built.

Install as a plugin from this repository:

```
/plugin marketplace add Mingyu1998/MDBuddy
/plugin install mdbuddy@mdbuddy
```

or copy `skills/mdbuddy` into `~/.claude/skills/`. Then ask, for example:

> Build KcsA (1BL8) in a POPC:CHOL 7:3 bilayer with Martini 3 and 150 mM NaCl for GROMACS.

### 2. As a plain command-line tool

```bash
pip install -e ".[build]"         # mdbuddy + vermouth (martinize2) + COBY
mdbuddy doctor                    # what is installed / missing, with install commands
mdbuddy init kcsa --pdb-id 1bl8   # commented spec template
mdbuddy validate kcsa.yaml        # every setting, defaults filled in
mdbuddy build kcsa.yaml           # -> kcsa/gromacs/
cd kcsa/gromacs && ./run.sh
```

## Requirements

`mdbuddy doctor` checks all of these and prints the install command for anything missing:

| | Needed for | Install |
|---|---|---|
| Python ≥ 3.9, PyYAML | everything | |
| vermouth (martinize2) | coarse-graining | `pip install vermouth` |
| COBY (tested 1.0.15) | membrane, insertion, solvation | `pip install COBY` |
| DSSP 2.2.1/3.0.0 **or** mdtraj | secondary structure | `pip install mdtraj` or `conda install -c salilab dssp` |
| GROMACS ≥ 2021 | final `grompp` check and running | `conda install -c conda-forge gromacs` / `module load gromacs` |

MDBuddy never installs anything itself.

## A spec

```yaml
mdbuddy: 1
name: kcsa_popc_m3
protein:
  source: opm            # pre-oriented in the membrane
  id: 1bl8
coarse_grain:
  forcefield: martini3   # martini3 | martini22 | elnedyn22
  model: elastic         # elastic | go | none
membrane:
  upper: {POPC: 7, CHOL: 3}
  lower: {POPE: 5, POPS: 2, CHOL: 3}
solvent:
  salt_molarity: 0.15
simulation:
  temperature: 310
  production_ns: 1000
```

All options: [`skills/mdbuddy/references/spec_reference.md`](skills/mdbuddy/references/spec_reference.md).

Lipid names differ between Martini 2, the 2021 Martini 3 release and the 2025 Martini 3 lipidome. For example,
`DLPC`, `PIPC` and `POP2` are different molecules, and cardiolipin is `CDL2` vs `TOCL`. MDBuddy checks every
lipid against a table generated from the official lipidome naming spreadsheet and converts between them:

```bash
mdbuddy lipids show DLPC CDL2 POP2
mdbuddy lipids convert CDL2 POP2 PIPC --from m2 --to m3
```

Details and the full list of pitfalls: [`skills/mdbuddy/references/lipids.md`](skills/mdbuddy/references/lipids.md).
Background on the choices: [`skills/mdbuddy/references/decisions.md`](skills/mdbuddy/references/decisions.md).
Every mdp parameter, line by line: [`skills/mdbuddy/references/mdp_parameters.md`](skills/mdbuddy/references/mdp_parameters.md).
More examples: [`examples/`](examples).

## What you get

```
kcsa_popc_m3/
  spec.resolved.yaml    mdbuddy.lock.json    CITATIONS.md
  00_input/  01_clean/  02_cg/  03_membrane/
  gromacs/
    step5_input.gro  topol.top  toppar/  index.ndx
    step6.0_minimization.mdp (soft-core)  step6.1_minimization.mdp
    step6.2 … step6.7_equilibration.mdp  step7_production.mdp
    run.sh  README.md
```

File names follow CHARMM-GUI conventions, so existing analysis scripts keep working.

## Reproducibility

- Force-field files are downloaded from **pinned commits** (marrink-lab/martini-forcefields,
  MikkelDA/COBY) and **verified by SHA-256** before use. A mismatch stops the build.
- `mdbuddy.lock.json` records tool versions, the input structure checksum, every
  parameter file with its source and checksum, every command line, and output checksums.
- Random seeds (COBY placement, initial velocities) come from the spec.
- `03_membrane/coby_build.py` and `02_cg/martinize2_cmd.sh` let you re-run any step by hand.

## Development

```bash
pip install -e ".[dev]"
pytest
```

The tests use small stand-ins for martinize2, COBY and the force-field downloads, so
they run anywhere without the real tools.

## Status

v0.1.0. The pipeline has been tested end-to-end only with stand-in tools so far. The first
real builds (KcsA with Martini 3 and Martini 2.2) are next. Please open an issue with
the `mdbuddy.lock.json` of any build that fails.

## Citing

Each build writes a `CITATIONS.md` (Martini, martinize2, COBY, OPM, GROMACS). Please
cite the underlying tools. MDBuddy only puts them together.

## License

Apache-2.0
