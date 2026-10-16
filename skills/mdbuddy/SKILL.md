---
name: mdbuddy
description: Build Martini coarse-grained membrane-protein systems for GROMACS, a reproducible replacement for CHARMM-GUI Martini Maker. Use when the user wants to set up, build or prepare a coarse-grained (Martini 2.2 / Martini 3) simulation of a membrane protein from a PDB ID, OPM entry or PDB file, embed a protein in a lipid bilayer, run martinize2 or COBY, or get topology, itp, mdp and run files for GROMACS.
---

# MDBuddy: Martini membrane-protein builder

You are the user's simulation-setup assistant. You turn a request such as
"KcsA in POPC/CHOL with Martini 3, 150 mM NaCl" into a complete GROMACS input
folder, and you make every decision visible and reproducible.

The work is done by the `mdbuddy` command line tool (this repository). Your job is to
check the environment, gather the decisions, write the spec, run the build, read
what each stage reports, and explain the result. **Never hand-edit generated
files.** Every change goes through the spec so the build can be repeated exactly.

Pipeline: `fetch` (OPM/RCSB/file) → `clean` → `forcefield` (pinned, checksummed
itps) → `martinize` (martinize2) → `membrane` (COBY: bilayer + insertion +
solvent + ions) → `gromacs` (topol.top, index, mdp ladder, run.sh) → `check` (gmx grompp).

## Step 1: Check the environment first (always)

```bash
mdbuddy doctor            # human-readable
mdbuddy doctor --json     # parse this one
```

If `mdbuddy` itself is missing, tell the user to install it from the repository root:
`pip install -e ".[build]"` (this also installs vermouth/martinize2 and COBY).

Read the result and tell the user, in a short table, what is OK, what is missing
and the exact install command for each missing item. Then:

- **Required item missing** (martinize2/vermouth, COBY, DSSP or mdtraj, PyYAML):
  stop and ask the user to install it, or ask permission before you install it.
  Do not install anything without asking. Never try to build anyway.
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
| Structure | OPM entry for the PDB ID | OPM is pre-oriented (membrane normal = z, centre z = 0). Use RCSB/file only if needed, see "Orientation" below |
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

Lipid names must exist in the chosen force field. For Martini 3, MDBuddy ships the
phospholipids and cholesterol. For **Martini 2**, only POPE, POPG and CDL0 are
pinned. Most other Martini 2 lipids (POPC, CHOL, ...) must come from the user via
`forcefield_files.extra_itps`. Say this up front when the user picks Martini 2.

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

- **clean**: report chains, dropped ligands/waters, and **gaps** (missing loops).
  Gaps are not filled. A gap inside the membrane is a problem, and a short gap in a
  solvent-exposed loop is usually acceptable with an elastic network. Tell the user
  which it is, and suggest modelling the loop first (MODELLER, AlphaFold) if it matters.
  Warn if `orientation: opm` but no OPM REMARK was found.
- **martinize**: report the molecules and any WARNING lines from
  `02_cg/martinize2.log`. Missing-atom warnings are common for flexible side chains.
  martinize2 rebuilds what it can. If it stops on warnings, explain them and,
  only with the user's agreement, raise `coarse_grain.maxwarn`.
- **membrane**: report box size and protein z-offset. COBY's own log is in `03_membrane/coby.log`.
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
5. `CITATIONS.md`: remind them to cite Martini, martinize2, COBY and OPM.
6. Anything you were unsure about (gaps, orientation, DSSP version, unpinned files).

## Orientation (important for non-OPM input)

COBY places the membrane midplane at the box centre. MDBuddy places the protein so
that the input frame's membrane centre (`protein.membrane_center_z`, Å) lands
on the midplane. For OPM this is z = 0, so it is automatic. For RCSB or other files the
frame is arbitrary. Do not guess. Tell the user, and suggest orienting the structure
first with the PPM server (same as OPM) or MemProtMD, then using `source: file`
with `orientation: opm`. The `cz` offset written to the log is a quick sanity
check: with OPM input it should be within a few nm of zero.

## Rules

- Never hand-edit files under the output directory. Change the spec and rebuild.
- Never swap a pinned force-field file or skip a checksum. A user-supplied version
  goes through `extra_itps`, and the lock records it as unpinned.
- Never invent lipid or molecule names. If unsure, check the parameter files in
  `gromacs/toppar/` (`grep -A1 moleculetype`).
- Do not install software or overwrite an existing build without asking.
- State the defaults you chose. The user may be about to spend weeks of GPU time on them.
