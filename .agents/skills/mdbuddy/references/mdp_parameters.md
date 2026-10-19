# MDP parameters, line by line

MDBuddy's `.mdp` files start from the lab's Martini inputs,
[DJ004/MolDynamics_gromacs `MDP_files/CGMD_*.mdp`](https://github.com/DJ004/MolDynamics_gromacs/tree/main/MDP_files)
(reviewed at commit `3db57dc`). Each parameter is listed below with the value MDBuddy writes, the value in those
files, and why. Rows marked **changed** are deliberate deviations. **fix** marks an error in the source files.

Lab files reviewed: `CGMD_minimization.mdp` (EM), `CGMD_NVT_equilibration.mdp` (NVT),
`CGMD_NPT_equilibration.mdp` (NPT-eq), `CGMD_NPT_production.mdp` (prod). The `ATMD_*` files are all-atom and do not
apply to Martini.

## Protocol

| MDBuddy stage | Lab equivalent | Notes |
|---|---|---|
| step6.0_minimization: EM, `-DFLEXIBLE` | EM | Constraints replaced by stiff harmonic bonds, so a freshly packed system relaxes without LINCS trouble. The lab file's soft-core block is **not** used, see below |
| step6.1_minimization: EM, constraints | none | Repeats minimisation with the constraints the MD stages actually use. CHARMM-GUI does the same |
| step6.2–6.6_equilibration: NPT, dt 2→20 fs, restraints 1000→50 | NVT + NPT-eq | Short time-step ramp, so a freshly built system does not blow up at 20 fs |
| step6.7_equilibration: NPT, 20 fs, 5 ns, restraints 10 | NPT-eq (5 ns, 20 fs) | Same length and time step as the lab NPT equilibration |
| step7_production: NPT, Parrinello-Rahman | prod | Same as lab file |

There is no separate NVT stage. With the C-rescale barostat the system can go straight from minimisation to
NPT (CHARMM-GUI Martini Maker does this too). The lab NVT file also has the electrostatics bug below, and it does
not generate velocities.

## Run control

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| title | (not written) | set | **changed**: obsolete in current GROMACS; grompp ignores it with a note |
| define (step6.0) | `-DPOSRES -DPOSRES_FC=1000 -DFLEXIBLE` | `;-DFLEXIBLE` (commented) | **fix**: the Martini 3 cholesterol itp says *"Add define=-DFLEXIBLE to select harmonic bonds for minimization purposes"*; martinize2 writes the same `#ifdef FLEXIBLE` block for the protein. Backbone restraints keep the protein in its OPM position |
| define (step6.1) | `-DPOSRES -DPOSRES_FC=1000` | — | **changed**: no `-DFLEXIBLE`, so the second EM runs with the constraints the MD stages use |
| define (EM, Martini 2) | `-DPOSRES -DPOSRES_FC=1000` | — | **changed**: no `-DFLEXIBLE` at all. No Martini 2 itp has an `#ifdef FLEXIBLE` block, and grompp counts an unused `-D` macro as a *warning*, which fails MDBuddy's own `-maxwarn 0` check. MDBuddy writes `-DFLEXIBLE` only when an included itp defines it, so step6.0 and step6.1 become the same minimisation run twice for Martini 2 (harmless; the second converges in one step) |
| comm-mode (EM) | `none` | (not set) | **changed**: with the GROMACS default (`Linear`) grompp emits a note about removing COM motion while position restraints are on. A minimisation has no COM motion to remove, so the note is pure noise |
| define (eq) | `-DPOSRES -DPOSRES_FC=<1000…10>` | `;-DPOSRE_BB` (commented) | martinize2 writes `#ifdef POSRES` with a `POSRES_FC` force constant, not `POSRE_BB` |
| define (prod) | empty | empty | same |
| integrator | steep / md | steep / md | same |
| tinit | default 0 | 0 (EM) | same |
| nsteps (EM) | 5000 | 5000 | same |
| emtol | 1000 | 1000 | same |
| emstep | 0.01 | default 0.01 | same, written explicitly |
| dt | 0.002→0.02 (eq), 0.02 (prod) | 0.02 | ramp only during equilibration |
| nstcomm | 100 | 100 | same |
| comm-grps | `SOLU_MEMB SOLVENT` | empty (= System) | **changed**: separate COM removal for membrane + protein and for solvent keeps the bilayer from drifting against the water (CHARMM-GUI default). Set `simulation.comm_groups: [System]` for the lab behaviour. Never put Protein in its own COM group |

## Output

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| nstxout / nstvout / nstfout | 0 | 0 | no .trr, same |
| nstlog, nstenergy (prod) | 25000 (500 ps) | 25000 | same, set by `simulation.output_ps` |
| nstxout-compressed (prod) | 25000 (500 ps) | 25000 | same, set by `simulation.output_ps` |
| compressed-x-precision | 100 (eq/prod), 1000 (EM) | 25000 (NPT), 100 (NVT), 1000 (EM) | **fix**: 25000 looks like a copy of the output frequency. Precision is *1/resolution in nm*. 100 (0.01 nm) is the Martini standard and gives smaller files |
| compressed-x-grps | default (all) | empty | same |
| energygrps | (not written) | System | not needed. The default is already one group |

## Neighbour searching

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| cutoff-scheme | Verlet | Verlet | same |
| nstlist | 20 | 20 | same; mdrun increases it automatically on GPUs |
| ns-type / ns_type | (not written) | grid | **changed**: obsolete since the group scheme was removed (GROMACS 2020); grompp ignores it with a note |
| pbc | xyz | xyz | same |
| verlet-buffer-tolerance | 0.005 | 0.005 (EM, NVT), default (NPT) | same |

## Electrostatics and van der Waals

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| coulombtype | reaction-field | Reaction-Field | same |
| rcoulomb | 1.1 | 1.1 | same |
| epsilon_r | 15 | 15 (EM, NPT, prod); **commented out in NVT** | **fix**: with `; epsilon_r` commented out the NVT file runs with the default epsilon_r = 1. Electrostatics are then 15× too strong |
| epsilon_rf | 0 | default 0 (EM, NPT, prod); **15 in NVT** | **fix**: Martini uses epsilon_rf = 0 (infinite). The NVT file puts 15 here instead of in epsilon_r. Written explicitly in every MDBuddy file |
| vdw-type | cutoff | Cut-off | same |
| vdw-modifier | Potential-shift | Potential-shift / Potential-shift-verlet | same (identical with the Verlet scheme) |
| rvdw | 1.1 | 1.1 | same |
| coulomb-modifier | (not written) | commented out | same; must not be used with reaction-field |

## Why there is no soft-core minimisation

The lab EM file minimises with a soft-core free-energy setup
(`free-energy = yes`, `init-lambda = 0.01`, `sc-alpha = 4`, `sc-power = 2`, `sc-coul = yes`,
`couple-moltype = system`, `couple-lambda0 = vdw-q`, `couple-lambda1 = none`, `couple-intramol = yes`)
to push apart overlapping beads. **MDBuddy does not write those lines**, because GROMACS cannot run them
for a Martini 3 protein system:

```
Fatal error:
Function type Restricted Angles does not support currentely free energy calculations
```

Martini 3 uses restricted bending angles (angle function type 10) for the protein backbone — martinize2
writes ~840 of them for KcsA — and for several lipids. `grompp` accepts `free-energy = yes` anyway and
writes a `.tpr`; `mdrun` then refuses it. Checked on GROMACS 2021.4: the error appears whatever is
coupled (`couple-moltype = system` or a single molecule type) and with `couple-intramol = no`, so there
is no setting that makes it work. The lab file predates Martini 3 proteins.

Instead, step6.0 minimises with `-DFLEXIBLE`, which turns the constraints in the protein and sterol itps
into stiff harmonic bonds, and step6.1 repeats the minimisation with the constraints restored. On the
KcsA/POPC test system step6.0 converged to Fmax < 1000 kJ/mol/nm in 2709 steps
(E_pot = -2.68e5 kJ/mol) and step6.1 was already converged (Fmax = 172 kJ/mol/nm), so the pair is enough
to clean up COBY's packing without soft-core.

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| refcoord-scaling (EM) | (not written) | all | only matters with pressure coupling, so dropped from EM |

## Martini 2 (martini22 / elnedyn22)

The mdp files are the same for Martini 2 and Martini 3, apart from `-DFLEXIBLE` above. That is
deliberate: the cut-offs and reaction-field settings MDBuddy writes
(`coulombtype = reaction-field`, `rcoulomb = rvdw = 1.1`, `epsilon_r = 15`, `epsilon_rf = 0`,
Verlet lists, `vdw-modifier = Potential-shift`) are the "new-rf" set that de Jong et al.
(Comput. Phys. Commun. 2016) derived **for Martini 2**, as the validated replacement for the old
group-cut-off/shifted-potential scheme that Martini 2 was originally published with. They carry
over to Martini 3 unchanged. `epsilon_r = 15` assumes standard (non-polarisable) water; with
polarisable water it must be 2.5.

The 20 fs production time step is the Martini 2 standard and is also what the lab's
`CGMD_NPT_production.mdp` uses. Checked on GROMACS 2021.4 for both `martini22` and `elnedyn22`:
EM converges and the first NPT stage runs with no LINCS warnings (constraint RMS deviation
~1e-8), so the step-size ramp in step6.2-6.6 is enough to reach it safely.

Ion names: COBY 1.0.15 and the Martini 2 ion file MDBuddy pins (from COBY's
`M2_BacterialMembrane` tutorial) both call the ions `NA` and `CL`, so no renaming happens. The
cgmartini.nl `martini_v2.0_ions.itp` calls them `NA+` and `CL-` instead; `assemble.MOLECULE_ALIASES`
only renames a molecule the parameter files do not define into one they do, so either file works.

## Temperature coupling

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| tcoupl | v-rescale | V-rescale | same |
| tc-grps | `Protein MEMBRANE SOLVENT` | `Protein MEMBRANE SOLVENT` (NVT: SYSTEM) | same; MDBuddy's index.ndx uses these group names |
| tau-t | 1.0 | 1.0 | same |
| ref-t | `simulation.temperature` (default 310 K) | 313 (NPT, prod), 310 (NVT) | the lab files disagree with each other (313 vs 310). MDBuddy uses one value from the spec |

## Pressure coupling

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| pcoupl (eq) | C-rescale | berendsen | **changed**: GROMACS ≥ 2021 issues a grompp *warning* for Berendsen (it does not give a correct ensemble), so a `-maxwarn 0` workflow fails. C-rescale (Bernetti & Bussi 2020) relaxes the box just as robustly and gives the correct ensemble |
| pcoupl (prod) | Parrinello-Rahman | parrinello-rahman | same |
| pcoupltype | semiisotropic | semiisotropic | same |
| tau-p | 4 (eq), 12 (prod) | 12 | eq shortened so the box relaxes within the short ramp stages. Production is the same |
| compressibility | 3e-4 3e-4 | 3e-4 3e-4 | same |
| ref-p | 1 1 | 1.0 1.0 | same |
| refcoord-scaling | com (restrained stages only) | com | same; omitted when there are no restraints |

## Velocities

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| gen-vel | yes in step6.2, else no | yes (NPT-eq), **no (NVT)**, no (prod) | **fix**: the lab NVT file starts from minimised coordinates without generating velocities |
| gen-temp | = ref-t | 313 / 310 | from the spec |
| gen-seed | `seed` from the spec (default 2026) | -1 (random) | **changed**: -1 makes every build different. A fixed seed makes it reproducible. Change `seed` for replicas |
| continuation | no / yes, paired with gen-vel | no / yes | same |

## Bonds

| Parameter | MDBuddy | Lab | Why |
|---|---|---|---|
| constraints | none | none | same (only `[ constraints ]` from the itps) |
| constraint-algorithm | Lincs | Lincs | same |
| lincs-order | 8 | 4 | **changed**: GROMACS recommends a higher order for large time steps with virtual sites and for accurate EM. Martini 3 cholesterol (and Go-model proteins) have virtual sites and constraint triangles. Cost is negligible in CG |
| lincs-iter | 2 | default 1 | **changed**: same reason |
| lincs-warnangle | 30 | 30 | same |

## Spec keys that control the mdp files

`simulation.temperature`, `pressure`, `production_ns`, `output_ps`, `coupling_groups`, `comm_groups`, `seed`,
`coarse_grain.position_restraints` (none → no POSRES defines). Everything else is fixed by design, and changing it
means editing `src/mdbuddy/mdp.py` and this document together.
