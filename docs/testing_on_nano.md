# Testing MDBuddy on nano: prompts for a coding agent

Paste these into a coding agent on nano (Claude Code, Codex, ...), **one at a time, in order**.
Each prompt is self-contained: the agent there has not seen the development conversation, but it
reads `AGENTS.md` and the skill from the repository.

Never tested with the real tools so far: GROMACS (grompp/mdrun), Martini 2.2 builds, the Go model,
a real self-assembly run, and the CG islands with real data. Prompts 1 and 2 matter most.

---

## 0. Setup and context (paste first)

```
I'm continuing development of MDBuddy (github.com/Mingyu1998/mdbuddy), a CLI + agent skill that
builds Martini CG membrane-protein systems for GROMACS (PDB -> martinize2 -> COBY -> GROMACS files).
Read AGENTS.md, README.md and skills/mdbuddy/SKILL.md first.

It was developed in a container without GROMACS, so the generated mdp/top files have never been
through grompp or mdrun. This machine (nano) has GROMACS, so the goal of the next tasks is to test
and fix that.

Setup:
- Clone the repo (or git pull if it's already here).
- Make a fresh venv and run: pip install -e ".[build,dev]"
- Run pytest, then: mdbuddy doctor --network

Report the doctor table (GROMACS version, COBY 1.0.15, martinize2/vermouth, DSSP/mdtraj, PDBFixer).

Rules for all following tasks:
- Work on a new branch per task, never on main. Don't push or merge without asking me.
- Ask before any run longer than ~30 min, and check nvidia-smi/top so you don't take resources
  other people are using.
- Never hand-edit build outputs. Fix the code or spec and rebuild.
- Keep the unit tests tool-free (fake tools in tests/conftest.py).
```

## 1. grompp + short MD on the basic system

```
Build examples/kcsa_popc_m3.yaml with `mdbuddy build` (the check stage runs grompp -maxwarn 0 on
three mdps). Then in the gromacs/ folder run the whole ladder with run.sh: EM steps 6.0 and 6.1,
equilibration 6.2-6.7 as written, and production cut to ~50,000 steps
(gmx convert-tpr -nsteps, or a temporary copy of the mdp).

Report:
- every grompp warning or note, and whether each one is our fault
- EM: final potential energy and Fmax
- equilibration: box, density and temperature stability (gmx energy)
- whether the protein stays in the membrane

Fix any mdp/top problem in src/mdbuddy (mdp.py, assemble.py, topology.py), add a unit test for
it, and update skills/mdbuddy/references/mdp_parameters.md if a parameter changes
(then run python tools/sync_skill.py).
```

## 2. Martini 2.2 and Go model

```
Test the other force-field paths end to end, with grompp and EM plus ~10k steps of the first NPT
stage each:
(a) examples/kcsa_asym_go_m3.yaml (GoMartini 3, virtual sites, GO_VIRT define)
(b) examples/kcsa_bacterial_m22.yaml (Martini 2.2 with elastic network). Martini 2 lipid itps are
    not pinned. Find the lab's copies (or the cgmartini.nl files) and pass them via
    forcefield_files.extra_itps. Note where they came from.
(c) the same protein with elnedyn22.

Check:
- the Go model: go_atomtypes/nbparams included, virtual-site exclusions right, no LINCS warnings
- Martini 2: ion names (NA+/CL- vs NA/CL aliases in assemble.py), the 20 fs timestep, the
  reaction-field settings in the M2 mdps

Fix, test and document as before.
```

## 3. Self-assembly (needs a longer run, so ask me first)

```
Test the self-assembly path for real. Build KcsA with protein.orientation: self_assembly (POPC,
323 K is fine) and run the step7_self_assembly stage. Start with ~200-500 ns and ask me before
going longer.

Then run `mdbuddy selfassembly check <build_dir>` and compare its membrane normal and centre with
the OPM orientation of 1BL8:
- angle between the normals
- centre offset in Å
- whether the in/out side matches (positive-inside)

Also check the timing (did the bilayer form, how long it took) and whether HEAD_NAMES /
CORE_HALF_NM / HEAD_RANGE_NM in selfassembly.py need tuning on real data. Finally, rebuild with
the printed `orientation: manual` lines and confirm it builds and passes grompp.
```

## 4. CG protein islands with real data

```
Test the cg_proteins feature (examples/protein_islands_m3.yaml, spec_reference.md section
"Pre-coarse-grained proteins") on real data. Take a martinize2-built CG protein, or a frame with
several proteins from an earlier simulation: a .gro with its box plus the .top/itps.

Build it with placement: keep and membrane.box: from_structure, and set membrane.cz to the old
membrane midplane minus half the box height.

Check:
- the islands keep their coordinates in step5_input.gro (max displacement should be 0)
- no lipid or water overlaps in EM
- POSRES holds the islands during equilibration (the itps need an #ifdef POSRES block)
- grompp passes with -maxwarn 0

Then test a large patch with membrane.lipid_apl tuning, and examples/membrane_patch_m3.yaml (a
domain made with hole + patch), through EM and a short NPT. Report how the domain boundary relaxes.
```

## 5. Cell envelope with LPS (once the postdoc's files are on nano)

```
Test examples/cell_envelope_m22.yaml with the postdoc's real LPS model:
- the itp goes in forcefield_files.extra_itps
- a single-molecule structure goes in coby.kwargs.molecule_import, with the right
  upbead/downbead bead numbers

Note that COBY silently keeps its own definition if the lipid name already exists in its library;
see troubleshooting.md. Check:
- both membranes are placed at their cz without overlap
- the counter-ions (LPS needs Ca2+/Mg2+; check that solvent/ion settings can express this, and add
  spec support if not)
- grompp/EM/short NPT
- the semi-isotropic coupling still makes sense with two membranes

If LPS charges or ions need spec changes, propose them to me before implementing.
```

## 6. Agent/MCP check (optional)

```
Register the MDBuddy MCP server with this agent (`mdbuddy mcp-config <agent>`, see docs/agents.md).
Then, using only the MCP tools, build examples/kcsa_popc_m3.yaml into a new folder, as a new user
would. Report anything confusing in the tool descriptions or the skill.
```
