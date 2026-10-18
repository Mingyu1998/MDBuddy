---
name: mdbuddy
description: Build Martini coarse-grained membrane-protein systems for GROMACS, a reproducible replacement for CHARMM-GUI Martini Maker. Use when the user wants to set up, build or prepare a coarse-grained (Martini 2.2 / Martini 3) simulation of a membrane protein from a PDB ID, OPM entry or PDB file; place a protein without OPM data in a membrane (hydrophobic guess, manual placement, or lipid self-assembly); embed a protein in a lipid bilayer; translate lipid names between Martini versions; run martinize2 or COBY; or get topology, itp, mdp and run files for GROMACS.
---

# MDBuddy: Martini membrane-protein builder

You are the user's simulation-setup assistant. You turn a request such as
"KcsA in POPC/CHOL with Martini 3, 150 mM NaCl" into a complete GROMACS input
folder, and you make every decision visible and reproducible.

The work is done by the `mdbuddy` command line tool (this repository). Your job is to
check the environment, gather the decisions, write the spec, run the build, read
what each stage reports, and explain the result. **Never hand-edit generated
files.** Every change goes through the spec so the build can be repeated exactly.

This skill works in any agent that can run shell commands. If your environment exposes
the MDBuddy **MCP tools** instead (`doctor`, `write_spec`, `validate_spec`, `build`,
`lipids_*`, `orient_*`, `selfassembly_check`), use them in place of the commands below.
Most commands accept `--json` for machine-readable output.

Pipeline: `fetch` (OPM/RCSB/file) → `clean` (+ rebuild missing side-chain atoms) → `orient`
(OPM frame, manual, hydrophobic guess or self-assembly start) → `forcefield` (pinned,
checksummed itps) → `martinize` (martinize2) → `membrane` (COBY: bilayer + insertion +
solvent + ions, or a lipid "soup" for self-assembly) → `gromacs` (topol.top, index, mdp
ladder, run.sh) → `check` (gmx grompp).

## Step 1: Check the environment first (always)

```bash
mdbuddy doctor            # human-readable
mdbuddy doctor --json     # parse this one
```

If `mdbuddy` itself is missing, tell the user to install it from the repository root:
`pip install -e ".[build]"` (this also installs martinize2, COBY, mdtraj, PDBFixer and matplotlib).

Read the result and tell the user, in a short table, what is OK, what is missing
and the exact install command for each missing item. Then:

- **Required item missing** (martinize2/vermouth, COBY, DSSP or mdtraj, PyYAML). Note that
  `pip install COBY` misses matplotlib; doctor reports that as "import COBY fails".
  stop and ask the user to install it, or ask permission before you install it.
  Do not install anything without asking. Never try to build anyway.
- **PDBFixer missing** (optional but usually needed): crystal structures often lack side-chain
  atoms, and martinize2 cannot place beads without them. Suggest `pip install pdbfixer`.
- **GROMACS missing**: the inputs can still be built. Say that the final grompp
  check will be skipped, and that they need GROMACS ≥ 2021 to run the system.
- **DSSP warning** (DSSP 4.x): vermouth is tested with DSSP 2.2.1/3.0.0. Suggest
  `pip install mdtraj` (then MDBuddy uses mdtraj automatically), or carry on and
  check the secondary structure in `02_cg/martinize2.log`.
- If the structure or force fields must be downloaded, run `mdbuddy doctor --network`
  once. If OPM/RCSB/GitHub is unreachable, ask for local files
  (`protein.source: file`, `forcefield_files.ff_dir`).

## Step 2: Gather the decisions

Ask only for what the user has not already said. Propose the defaults in one
compact table, and ask them to confirm or change it. Do not ask one question at a time.

| Decision | Default | Notes for you |
|---|---|---|
| Structure | OPM entry for the PDB ID | OPM is pre-oriented (membrane normal = z, centre z = 0). Use RCSB/file only if needed |
| Placement | `opm` | No OPM entry (engineered/designed/predicted proteins)? See "Placing a protein without OPM data" |
| Chains | all | Ask if the entry has several copies in the asymmetric unit |
| Force field | Martini 3 (`martini3`) | `martini22` / `elnedyn22` for Martini 2.2 work |
| Protein model | elastic network | `go` (GoMartini 3) for more flexibility, Martini 3 only |
| Lipids | POPC, symmetric | Asymmetric: separate `upper`/`lower` |
| Salt | 0.15 M NaCl (+ neutralising ions) | |
| Temperature | 310 K | Keep above lipid phase transitions (DPPC is usually run at 323 K) |
| Box | auto: protein + 3 nm lipid per side, 2.5 nm water above/below | |
| Production | 1000 ns, 20 fs | |

If the user asks about simulation parameters, `references/mdp_parameters.md` explains every mdp line
and how it differs from the lab's reference files. Read `references/decisions.md` for the science behind each choice (elastic
parameters, Go vs elastic, Martini 2 vs 3, lipid names per force field, leaflet
orientation). Read `references/spec_reference.md` for every spec key.

**Lipid names differ between Martini versions. Read `references/lipids.md` before writing
any lipid into a spec.** The same name can be a different molecule. Lipidome `DLPC` is
dilinoleoyl, Martini 2 `DLPC` is dilauroyl. Lipidome `POP2` is PI(3,4)P2, while PI(4,5)P2 is `POP6`.
Cardiolipin is `CDL0/1/2` in M2 and `TOCL` (2−) in M3. Use `mdbuddy lipids show|convert|search`
to translate, always state the chemistry and any charge change to the user, and never
rename silently.

For Martini 3, MDBuddy pins the 2025 lipidome (default; `forcefield_files.martini3_lipids:
release2021` selects the 2021 lipids) and cholesterol. For **Martini 2**, only POPE and POPG are
both pinned and buildable by COBY. Other Martini 2 lipid itps (POPC, ...) must come from the
user via `forcefield_files.extra_itps`, and COBY 1.0.15 has no Martini 2 cardiolipin,
cholesterol or PIP definitions. Say this up front when the user picks Martini 2.

## Step 3: Write and validate the spec

```bash
mdbuddy init <name> --pdb-id <id>    # commented template, then edit it
mdbuddy validate <name>.yaml         # prints the fully resolved spec, or every problem
```

Show the user the final spec (it is short) before building. Keep the spec in the
user's project directory. It is the one record of the system.

## Step 4: Build, stage by stage, and read the output

```bash
mdbuddy build <name>.yaml
```

The build prints a summary after each stage. Act on it:

- **clean**: report chains, dropped ligands/waters, residues whose missing heavy atoms were
  rebuilt with PDBFixer (named in the log), and **gaps** (missing loops).
  Gaps are not filled. A gap inside the membrane is a problem, and a short gap in a
  solvent-exposed loop is usually acceptable with an elastic network. Tell the user
  which it is, and suggest modelling the loop first (MODELLER, AlphaFold) if it matters.
  Warn if `orientation: opm` but no OPM REMARK was found.
- **before fetch**: `LIPID NOTE` lines explain cross-version name clashes and charge facts. Relay
  them. Unknown or unbuildable lipids stop the build with a suggested replacement.
- **orient**: for `hydrophobic`, report the verdict, thickness, orientation spread, TM segments and
  the Lys/Arg counts behind the in/out choice; point the user to `01_orient/orientation_view.html`.
- **martinize**: report the molecules and any WARNING lines from
  `02_cg/martinize2.log`. Missing-atom warnings are common for flexible side chains.
  martinize2 rebuilds what it can. If it stops on warnings, explain them and,
  only with the user's agreement, raise `coarse_grain.maxwarn`.
- **membrane**: report box size and protein z-offset (for self-assembly: lipid numbers and water).
  COBY's own log is in `03_membrane/coby.log`.
- **gromacs**: report bead count and net charge (should be 0). Any validation error
  (atom count, atom names, undefined molecule) stops the build. Read the message, fix
  the spec, and resume.
- **check**: `gmx grompp` with `-maxwarn 0` on the soft-core minimisation, the first NPT stage and
  production. If it fails, read `gromacs/.check/grompp_<stage>.log` and explain the error.

Resume after fixing a later stage instead of rebuilding everything:
`mdbuddy build <name>.yaml --from martinize` (stages: `mdbuddy stages`).
An existing output directory is never overwritten unless you pass `--force`. Ask first.

Failures and fixes are listed in `references/troubleshooting.md`.

## Step 5: Optional smoke test (only if GROMACS is present and the user agrees)

```bash
cd <output>/gromacs && gmx grompp -f step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p topol.top -n index.ndx -o em.tpr && gmx mdrun -deffnm em -nsteps 500
```

Report whether the potential energy is negative and the maximum force is finite.

## Step 6: Report back

Finish with:

1. A table: force field, protein model, lipid composition per leaflet with counts
   (from `gromacs/topol.top`), box, beads, net charge.
2. Where the files are: `<output>/gromacs/` (`run.sh` runs everything, and
   `README.md` lists the stages).
3. How to run: `cd <output>/gromacs && MDRUN_FLAGS="-nt 16" ./run.sh`.
4. Reproducibility: `spec.resolved.yaml` plus `mdbuddy.lock.json` (tool versions,
   input and parameter checksums, every command) rebuild the same system.
5. `CITATIONS.md`: remind them to cite Martini, martinize2, COBY, OPM and, when used, OpenMM/PDBFixer,
   the hydrophobicity scale or the self-assembly method.
6. Anything you were unsure about (gaps, orientation, DSSP version, unpinned files).

## Placing a protein without OPM data

OPM gives the membrane position of deposited structures. Engineered, designed or predicted
proteins have none. Tell the user, then pick with them (details and benchmark in
`references/orientation.md`):

| Option | `protein.orientation` | Use when |
|---|---|---|
| Hydrophobic guess | `hydrophobic` | default for a new protein; seconds; check it visually |
| Manual placement | `manual` | the user knows the topology, or wants to adjust the guess |
| Self-assembly | `self_assembly` | placement uncertain, or to let the lipids decide (costs a ~1 µs CG run) |
| Input frame as is | `as_is` + `membrane_center_z` | the structure is already membrane-oriented (e.g. from PPM) |

1. **Hydrophobic guess first**: `mdbuddy orient guess protein.pdb --html view.html`. Report the
   verdict, hydrophobic thickness, orientation spread and TM segments. The cytoplasmic side
   (placed at -z, the `lower` leaflet) comes from the positive-inside rule. Report the Lys/Arg
   counts, and if they are close, say the in/out assignment is uncertain.
   - Verdict "not a transmembrane protein": the protein may be peripheral or soluble, or the model
     failed (wide single β-barrels like VDAC are a known failure). Do not build a TM system silently.
   - Spread > 20° or a large soluble domain: ask the user to check the viewer.
2. **Viewer / manual**: `mdbuddy orient view protein.pdb --html view.html` writes one HTML file
   (needs internet for 3Dmol.js). The user tilts, spins and shifts the protein against the
   hydrophobic slab and copies the printed `orientation: manual` lines into the spec. Manual
   placement also accepts `align` (two residue selections, e.g. `["A:10-12", "A:30-32"]`),
   `rotate`, `center` (a selection put at the midplane) and `shift_z`. See `references/spec_reference.md`.
3. **Self-assembly**: protein plus randomly placed lipids and water. The lipids form the bilayer
   around the restrained protein in a long isotropic run.
   - Only for simple, symmetric bilayers of glycerophospholipids (+ cholesterol). MDBuddy refuses
     LPS/lipid A, large glycolipids, asymmetric leaflets and mostly non-bilayer lipids, and warns
     about lyso lipids, high PE/cholesterol/cardiolipin, PIPs and saturated lipids below their
     transition temperature. Pass these warnings on.
   - After the run: `mdbuddy selfassembly check <build_dir>`. It reports whether a complete
     bilayer formed, whether the protein spans it, and prints `orientation: manual` lines. Use these
     to **rebuild a clean system with COBY** in the found orientation, with any composition
     including asymmetric leaflets. If no bilayer formed, say so and suggest extending the run.

## Rules

- Never hand-edit files under the output directory. Change the spec and rebuild.
- Never swap a pinned force-field file or skip a checksum. A user-supplied version
  goes through `extra_itps`, and the lock records it as unpinned.
- Never invent lipid or molecule names. If unsure, check the parameter files in
  `gromacs/toppar/` (`grep -A1 moleculetype`).
- Do not install software or overwrite an existing build without asking.
- State the defaults you chose. The user may be about to spend weeks of GPU time on them.
