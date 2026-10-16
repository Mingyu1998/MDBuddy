# Spec reference

A spec is one YAML file. Anything left out takes the default shown. `mdbuddy validate
spec.yaml` prints the fully resolved spec. The same content is saved as
`spec.resolved.yaml` in the output folder.

```yaml
mdbuddy: 1                 # spec format version
name: kcsa_popc_m3         # letters, digits, _ - . ; also the default output folder
seed: 2026                 # COBY placement and gen-seed for velocities
output_dir: null           # default: ./<name> next to the spec

protein:
  source: opm              # opm | rcsb | file
  id: 1bl8                 # for opm / rcsb
  path: null               # for file (PDB format)
  chains: null             # [A, B] or null for all
  keep_hetatm: []          # HETATM residue names to keep (martinize2 must know them!)
  orientation: opm         # opm (input is membrane-oriented) | as_is; default from source
  membrane_center_z: 0.0   # Angstrom, membrane midplane in the input frame (0 for OPM)

coarse_grain:
  forcefield: martini3     # martini3 | martini22 | elnedyn22
  model: elastic           # elastic | go (martini3 only) | none
  elastic:                 # defaults depend on the force field
    force: 700             # kJ/mol/nm^2   (-ef)
    lower: 0.5             # nm            (-el)
    upper: 0.9             # nm            (-eu)
  go: {epsilon: 9.414, low: 0.3, up: 1.1, res_dist: 3}
  merge_chains: auto       # auto (merge oligomers when elastic/go) | all | none | [[A, B], [C, D]]
  dssp: auto               # auto | mdtraj | none | <dssp executable> | "ss:<DSSP string>"
  position_restraints: backbone   # none | backbone | all (used in minimisation/equilibration)
  posres_force: 1000       # default force in the itp; mdp files override per stage
  cys: auto                # disulfide detection (martinize2 -cys)
  molecule_name: protein   # moleculetype prefix -> protein_0, protein_1, ...
  maxwarn: 0               # martinize2 -maxwarn; raise only after reading the warnings
  extra_args: []           # extra martinize2 arguments, recorded in the lock

membrane:
  lipids: {POPC: 1}        # symmetric composition (ratios)
  upper: null              # per-leaflet composition; overrides lipids
  lower: null
  apl: 0.6                 # nm^2 per lipid
  upper_apl: null          # per-leaflet APL
  lower_apl: null
  box: auto                # or [x, y, z] in nm
  box_type: rectangular
  lipid_padding: 3.0       # nm of membrane around the protein (auto box)
  water_padding: 2.5       # nm of water above/below the protein (auto box)
  min_box_z: 10.0          # nm

solvent:
  water: W
  cation: NA
  anion: CL
  salt_molarity: 0.15      # mol/L, added on top of neutralising counter-ions

forcefield_files:
  optional: []             # martini3: sugars, small_molecules ; martini2: POPE, POPG, CDL0 (auto when used)
  extra_itps: []           # your own itp files (e.g. Martini 2 lipids); copied to toppar/ and included
  ff_dir: null             # folder with copies of the pinned files (offline use; checksums still verified)

simulation:
  temperature: 310         # K
  pressure: 1.0            # bar
  production_ns: 1000
  equilibration: standard  # see decisions.md
  output_ps: 500           # production output interval (xtc, energies, log); multiple of 2 ps
  coupling_groups: [Protein, MEMBRANE, SOLVENT]   # tc-grps; or [SOLU_MEMB, SOLVENT] or [System]
  comm_groups: [SOLU_MEMB, SOLVENT]  # comm-grps; or [System]
  gromacs_maxwarn: 0       # -maxwarn written into run.sh
```

## Output folder

```
<name>/
  spec.resolved.yaml       every setting, defaults filled in
  mdbuddy.lock.json        tool versions, input + parameter SHA-256s, every command run
  CITATIONS.md
  00_input/input.pdb       downloaded / copied structure
  01_clean/                cleaned structure + clean_report.json (dropped atoms, gaps)
  02_cg/                   martinize2 inputs/outputs, martinize2_cmd.sh, martinize2.log
  03_membrane/             coby_build.py (re-runnable), coby.log, system.gro/pdb
  gromacs/                 the folder to copy to the cluster
    step5_input.gro/.pdb, topol.top, toppar/, index.ndx (System, Protein, MEMBRANE, SOLU_MEMB, SOLVENT, BB),
    step6.0 (soft-core) and step6.1 _minimization.mdp, step6.2-6.7_equilibration.mdp, step7_production.mdp,
    run.sh, README.md
```
