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
PDB / OPM / file ──► clean (+ rebuild missing atoms) ──► orient ──► martinize2 ──► COBY ──► GROMACS inputs ──► grompp check
                                                          │
                     OPM frame · hydrophobic guess · manual (viewer) · lipid self-assembly
```

Scope: **membrane proteins from existing structures, including engineered or predicted ones without
OPM data. Martini 3 (2025 lipidome) and Martini 2.2, GROMACS.** Co-folding input, soluble proteins
and other engines are planned.

## Works with any AI agent

MDBuddy is a skill (open [Agent Skills](https://agentskills.io) `SKILL.md` format), an MCP server and a
CLI, so it works in Claude Code, OpenAI Codex, Gemini CLI, Cursor, VS Code/Copilot and other agents:

```bash
pip install -e ".[build]"                 # from a clone of this repository
mdbuddy skill install --agent all         # Claude Code, Codex, Gemini CLI skill folders
mdbuddy mcp-config codex                  # MCP snippet for codex | claude | gemini | cursor | vscode | claude-desktop
```

Claude Code users can also install the plugin: `/plugin marketplace add Mingyu1998/MDBuddy`, then
`/plugin install mdbuddy@mdbuddy`. Agents that read `AGENTS.md` find their instructions there. Details
for each agent are in [`docs/agents.md`](docs/agents.md). Then ask, for example:

> Build KcsA (1BL8) in a POPC:CHOL 7:3 bilayer with Martini 3 and 150 mM NaCl for GROMACS.

> My designed channel (design.pdb) has no OPM entry. Place it in a POPE:POPG 3:1 membrane.

The skill has the agent check your environment, agree the decisions that matter with you (sensible
defaults), write the spec, run the build stage by stage, explain warnings such as missing loops or
uncertain orientations, and report what was built.

## Or as a plain command-line tool

```bash
pip install -e ".[build]"         # mdbuddy + martinize2, COBY, mdtraj, PDBFixer
mdbuddy doctor                    # what is installed / missing, with install commands
mdbuddy init kcsa --pdb-id 1bl8   # commented spec template
mdbuddy validate kcsa.yaml        # every setting, defaults filled in
mdbuddy build kcsa.yaml           # -> kcsa/gromacs/
cd kcsa/gromacs && ./run.sh         # or ./run.sh step6.0 step6.7 for part of the ladder
```

## Proteins without OPM data

| `protein.orientation` | How |
|---|---|
| `hydrophobic` | implicit-membrane scan of exposed hydrophobicity plus the positive-inside rule. Checked on 14 OPM structures: 13 placed, normal within 6° for 11 |
| `manual` | `mdbuddy orient view protein.pdb` opens an interactive viewer (sliders against the hydrophobic slab) that prints the spec lines; or align residues/rotate/shift in the spec |
| `self_assembly` | protein + randomly placed lipids + water, the bilayer forms in a ~1 µs run. Then `mdbuddy selfassembly check` gives the orientation for a clean (also asymmetric) COBY rebuild. MDBuddy refuses lipids that cannot self-assemble (LPS, lipid A, gangliosides, asymmetric leaflets) |

```bash
mdbuddy orient guess design.pdb --html view.html
```

Method, benchmark and limitations: [`skills/mdbuddy/references/orientation.md`](skills/mdbuddy/references/orientation.md).

## Pre-coarse-grained proteins, large patches, custom lipids

For protein islands or clusters that are already Martini models (e.g. taken from an earlier
simulation), and for membranes built to a precise design:

| Spec | What it does |
|---|---|
| `protein.source: none` | no atomistic protein; martinize2 is not needed. Alone it builds a membrane-only patch |
| `cg_proteins: [{structure, topology, placement: keep}]` | inserts CG structures as arranged (COBY `center_protein:False`), checking bead count and names against the topology |
| `membrane.box: from_structure` | uses the box of the first CG structure, so kept coordinates stay exact |
| `membrane.lipid_apl: {LIPID: nm²}` | per-lipid area per lipid, to tune lipid density |
| `coby.membranes: [...]` | raw COBY membrane strings: several membranes (inner + outer, with `cz`), domains via hole + patch |
| `coby.kwargs` | other COBY keywords (`import_library`, `molecule_import`, `flooding`, ...). Referenced files are copied into the build and checksummed |
| custom lipids (e.g. LPS) | accepted when your itp (`forcefield_files.extra_itps`) defines them and `coby.kwargs` gives COBY their structure |

Examples: [`examples/protein_islands_m3.yaml`](examples/protein_islands_m3.yaml),
[`examples/membrane_patch_m3.yaml`](examples/membrane_patch_m3.yaml),
[`examples/cell_envelope_m22.yaml`](examples/cell_envelope_m22.yaml).

## Requirements

`mdbuddy doctor` checks all of these and prints the install command for anything missing:

| | Needed for | Install |
|---|---|---|
| Python ≥ 3.9, PyYAML | everything | |
| vermouth (martinize2) | coarse-graining | `pip install vermouth` |
| COBY (tested 1.0.15) | membrane, insertion, solvation | `pip install COBY` |
| DSSP 2.2.1/3.0.0 **or** mdtraj | secondary structure | `pip install mdtraj` or `conda install -c salilab dssp` |
| PDBFixer (OpenMM) | rebuilding missing side-chain atoms (most crystal structures) | `pip install pdbfixer` |
| matplotlib | imported by COBY 1.0.15 but not installed with it | `pip install matplotlib` |
| numpy (scipy optional) | orientation and self-assembly analysis | installed with COBY/vermouth |
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

Every changed or clashing name (generated from the official table):
[`docs/lipid_name_changes.md`](docs/lipid_name_changes.md). How MDBuddy handles them:
[`skills/mdbuddy/references/lipids.md`](skills/mdbuddy/references/lipids.md).
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
    step6.0_minimization.mdp (-DFLEXIBLE)  step6.1_minimization.mdp
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

v0.3.0. Tested end-to-end with the real tools (martinize2 0.15.0, COBY 1.0.15, PDBFixer/OpenMM 8):
- **KcsA (OPM):** Martini 3 / POPC, 20 incomplete residues rebuilt.
- **Bacteriorhodopsin:** randomly rotated, placed with `hydrophobic` (normal recovered within 0.7°).
- **KcsA self-assembly soup:** POPC, analysed on a synthetic assembled frame.
- **CG protein islands** (10 pre-coarse-grained rods kept in place, 0.0 nm displacement) in POPE/POPG/TOCL;
  a membrane-only 20 nm patch with a DPPC/CHOL domain (hole + patch); a custom lipid imported with
  `molecule_import`. COBY only, no martinize2.

GROMACS itself (grompp/mdrun) has not been run in the development environment yet. Please open an
issue with the `mdbuddy.lock.json` of any build that fails.

## Citing

Each build writes a `CITATIONS.md` (Martini, martinize2, COBY, OPM, GROMACS). Please
cite the underlying tools. MDBuddy only puts them together.

## License

Apache-2.0
