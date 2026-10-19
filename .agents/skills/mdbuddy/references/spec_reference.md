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
  source: opm              # opm | rcsb | file | none (no atomistic protein: cg_proteins only, or membrane only)
  id: 1bl8                 # for opm / rcsb
  path: null               # for file (PDB format)
  chains: null             # [A, B] or null for all
  keep_hetatm: []          # HETATM residue names to keep (martinize2 must know them!)
  repair: auto             # rebuild missing heavy atoms with PDBFixer (+ local OpenMM relaxation) | off
  orientation: opm         # opm | as_is | manual | hydrophobic | self_assembly (default: opm for OPM, else as_is)
  membrane_center_z: 0.0   # as_is: Angstrom, membrane midplane in the input frame
  manual:                  # orientation: manual (see orientation.md); applied in this order
    normal: null           # [x, y, z] membrane normal in the input frame (towards the non-cytoplasmic side)
    align: null            # or ["A:10-12", "A:30-32"]: vector between selection centroids becomes the normal
    center_point: null     # [x, y, z] Angstrom, a point on the membrane midplane (default: protein centre)
    rotate: null           # [rx, ry, rz] degrees, extra rotation (x, then y, then z)
    center: null           # selection put at the midplane, e.g. "A:15-35"
    shift_z: null          # Angstrom, final shift along the normal
  topology: positive_inside  # hydrophobic: which side is cytoplasmic (-z): positive_inside | as_found | flip
  pore_correction: auto    # hydrophobic: lipid-facing weighting for pores/beta-barrels: auto | on | off

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
  box: auto                # auto (one protein) | [x, y, z] nm | from_structure (box of the first cg_proteins file)
  cz: 0.0                  # nm, membrane midplane relative to the box centre (generated membrane)
  lipid_apl: {}            # per-lipid area per lipid, nm^2, e.g. {CHOL: 0.4}: changes how many fit (density)
  box_type: rectangular
  lipid_padding: 3.0       # nm of membrane around the protein (auto box)
  water_padding: 2.5       # nm of water above/below the protein (auto box)
  min_box_z: 10.0          # nm

cg_proteins: []            # pre-coarse-grained proteins, inserted as given (no martinize2, no orientation)
# - structure: islands.gro # CG .gro or .pdb (paths relative to the spec)
#   topology: islands.top  # its .top: [ molecules ] in structure order, #includes of the protein itps
#   itps: null             # or list the itps yourself ...
#   moleculetypes: null    # ... with [[NAME, count], ...] in structure order (overrides the .top's)
#   placement: keep        # keep: COBY center_protein:False, coordinates kept relative to the file's box
#                          #   centre (gro box / PDB CRYST1 required) | center: centre of geometry at the box centre
#   cx: 0.0                # nm shifts and degree rotations applied by COBY (cx cy cz rx ry rz)
#   coby_args: ""          # any other COBY protein sub-arguments for this structure

coby:                      # direct COBY control; everything ends up in 03_membrane/coby_build.py and the lock
  membrane_args: ""        # added to the generated membrane, before the leaflets (membrane-wide), e.g.
                           #   "hole:circle:cx:0:cy:0:radius:5" or "patch:rectangle:cx:0:cy:0:xlength:20:ylength:20"
  membranes: []            # raw COBY membrane strings that REPLACE the generated one (several allowed, e.g.
                           #   inner + outer membrane with cz offsets, or a domain via hole + patch). MDBuddy
                           #   adds params: for known lipids that have none; membrane.lipids/upper/lower are ignored.
  protein_args: ""         # added to the atomistic protein's COBY argument
  solvation_args: ""       # added to the solvation argument
  kwargs: {}               # other COBY keywords, e.g. import_library: [defs.py], molecule_import: [...],
                           #   flooding, stacked_membranes. Not allowed: box, params, itp_input, out_*, sn,
                           #   randseed and keywords starting with prot/memb/solv (use the keys above).
                           # Files named as file:<path> (and import_library .py files) are copied into
                           #   03_membrane/user_files and recorded with their SHA-256.

self_assembly:             # only for protein.orientation: self_assembly (symmetric membranes only)
  lipid_count: auto        # or a number; auto covers the box area at membrane.apl
  water_per_lipid: 25      # CG water beads per lipid (1 bead = 4 waters)
  run_ns: 1000             # length of step7_self_assembly
  restrain_protein: true   # backbone restraints during assembly
  align: hydrophobic       # initial alignment: hydrophobic (TM axis along z) | as_is

solvent:
  water: W
  cation: NA
  anion: CL
  salt_molarity: 0.15      # mol/L, added on top of neutralising counter-ions

forcefield_files:
  martini3_lipids: lipidome  # lipidome (Martini 3 lipidome 2025, default) | release2021 (original M3 lipids)
  optional: []             # martini3: sugars, small_molecules ; martini2: POPE, POPG, CDL0 (auto when used)
  extra_itps: []           # your own itp files (e.g. Martini 2 lipids); copied to toppar/ and included.
                           # A lipid MDBuddy does not know (e.g. LPS) is accepted when one of these defines it
                           # AND coby.kwargs gives COBY its structure (molecule_import ... library_types:lipid,
                           # or import_library); MDBuddy then cannot check its naming beyond the atom-name check.
  ff_dir: null             # folder with copies of the pinned files (offline use; checksums still verified)

simulation:
  temperature: 310         # K
  pressure: 1.0            # bar
  production_ns: 1000
  equilibration: standard  # see decisions.md
  output_ps: 500           # production output interval (xtc, energies, log); multiple of 2 ps
  coupling_groups: [Protein, MEMBRANE, SOLVENT]   # tc-grps; or [SOLU_MEMB, SOLVENT], [MEMBRANE, SOLVENT]
                           #   (default when there is no protein) or [System]
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
  01_clean/                cleaned structure + clean_report.json (dropped atoms, gaps, rebuilt residues)
  01_orient/               protein_oriented.pdb, orient_report.json, orientation_view.html
  02_cg/                   martinize2 inputs/outputs, martinize2_cmd.sh, martinize2.log; cg_protein_<i>.gro/pdb
  03_membrane/             coby_build.py (re-runnable), coby.log, system.gro/pdb, user_files/ (files named in coby:)
                           (+ lipid templates for self-assembly)
  gromacs/                 the folder to copy to the cluster
    step5_input.gro/.pdb, topol.top, toppar/, index.ndx (System, Protein, MEMBRANE, SOLU_MEMB, SOLVENT, BB),
    step6.0 (-DFLEXIBLE) and step6.1 _minimization.mdp, step6.2-6.7_equilibration.mdp, step7_production.mdp,
    run.sh, README.md
```

## Pre-coarse-grained proteins and membrane-only systems

- `protein.source: none` skips fetch, clean, orient and martinize2 (martinize2/vermouth/DSSP are then
  not required by `mdbuddy doctor`). With no `cg_proteins` the system is membrane only and needs an
  explicit `membrane.box`.
- `cg_proteins` can be combined with an atomistic protein; then set `membrane.box` explicitly
  (`auto` sizes the box around one protein). Self-assembly does not take `cg_proteins`.
- Each entry is checked before COBY runs: bead count and bead names must match the moleculetypes in
  order, moleculetype names must be unique, and itps may not `#include` other files.
- `placement: keep` keeps coordinates relative to the box centre of the file; use
  `membrane.box: from_structure` to keep them exactly, and `membrane.cz` to put the membrane at the
  height the proteins expect (e.g. the old membrane midplane minus half the box height).
- Position restraints come from the itps (`#ifdef POSRES`). martinize2 writes them; for your own
  itps MDBuddy warns if there are none (the proteins are then free during equilibration).
- Examples: `examples/protein_islands_m3.yaml`, `examples/membrane_patch_m3.yaml`,
  `examples/cell_envelope_m22.yaml`.
