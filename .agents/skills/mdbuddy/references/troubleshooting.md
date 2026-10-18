# Troubleshooting

Exit codes: `1` = spec invalid or environment not ready, `2` = a stage failed (the
message says which, and points to the log).

## fetch
- **HTTP 404 from OPM**: the entry is not in OPM (or uses a different PDB ID, e.g. a
  representative structure). Ask the user, or use `source: rcsb` and read the orientation section.
- **Network errors**: download elsewhere, then `source: file`.

## clean
- **chains not found**: OPM files can re-letter chains. Look at `00_input/input.pdb`.
- **gaps**: reported, not fixed. Model missing loops before MDBuddy if they matter.
- **"no OPM bilayer REMARK"** with `orientation: opm`: the file is probably not oriented.

## clean (missing atoms)
- **"residues with missing heavy atoms"**: crystal structures often lack side-chain atoms. martinize2
  cannot place a bead without atoms (it writes `nan`), so MDBuddy rebuilds them with PDBFixer and
  relaxes only the rebuilt atoms with OpenMM (`01_clean/repair_pdbfixer.py`). Install with
  `pip install pdbfixer`. Missing *residues* (gaps) are never built.

## orient
- **hydrophobic: "not a transmembrane protein"**: no hydrophobic belt was found. The protein may be
  peripheral or soluble, or a wide single β-barrel (VDAC-like, a known failure). Use `manual` (with
  the viewer) or `self_assembly`.
- **"orientation is poorly defined"** (spread > 20°): check `01_orient/orientation_view.html`. Large
  soluble domains can bias the simple model.
- **in/out uncertain** (similar Lys/Arg counts): confirm the topology from biology, then use
  `topology: flip` if needed.
- **viewer shows nothing**: it loads 3Dmol.js from cdn.jsdelivr.net, so it needs internet.

## self-assembly
- **refused**: asymmetric leaflets, LPS/lipid A, gangliosides and mostly non-bilayer lipids cannot
  self-assemble. Self-assemble a simple symmetric bilayer to find the position, then rebuild with COBY.
- **`selfassembly check`: no complete bilayer**: extend the run (`gmx convert-tpr -s step7_self_assembly.tpr
  -extend 500000 -o next.tpr`), raise the temperature for saturated lipids, or reduce PE/cholesterol.
  A vesicle or bilayer stack is not usable. Rebuild with a different seed.
- **protein does not span the bilayer**: inspect the trajectory. It may be peripheral, or the bilayer
  formed beside it.

## forcefield
- **checksum mismatch**: the download or `ff_dir` copy differs from the pinned file.
  Do not bypass it. If the user really wants a different version, use `extra_itps`.
- **"not defined in any ... parameter file"**: a lipid/ion name has no itp. Fix the
  name (check `toppar/*.itp`) or add the file via `extra_itps`.

## martinize
- Read `02_cg/martinize2.log`.
- **Unknown residue / no mapping**: ligands or modified residues kept with
  `keep_hetatm` need martinize2 mappings. Usually drop them (default) or mutate
  (`extra_args: ["-mutate", "A-MSE12:MET"]`).
- **Missing atoms / "more warnings than allowed"**: inspect which residues. Raising
  `coarse_grain.maxwarn` is acceptable for missing side-chain atoms that martinize2 repairs.
  Say so in the report.
- **DSSP errors with DSSP 4.x**: set `coarse_grain.dssp: mdtraj` (pip install mdtraj),
  or use DSSP 3.0.0.

## membrane (COBY)
- Read `03_membrane/coby.log` and `coby_stdout.log`.
- **protein shifted outside the box**: the input is not membrane-oriented (see
  "Orientation" in SKILL.md), or `membrane_center_z` is wrong.
- **lipid not found in COBY library**: COBY's lipid names follow Martini naming. The
  library is selected by force field (`default` = Martini 3, `M2` = Martini 2).
- You can re-run COBY by hand: `cd 03_membrane && python coby_build.py`.
- **cg_proteins shifted by half a box**: should not happen with MDBuddy. COBY 1.0.15 ignores the
  box line of a `.gro` for `center_protein:False`, so MDBuddy always gives COBY a PDB with CRYST1
  (`02_cg/cg_protein_<i>.pdb`; your file is kept as `cg_protein_<i>_input.*`).
- **"has N beads but the moleculetypes have M" / bead name mismatch** (cg_proteins): the structure
  and the `[ molecules ]` (or `moleculetypes`) are not in the same order or count. Fix the input.
- **custom lipid: bead-count mismatch in COBY, or "detected multiple times"**: COBY already has a
  lipid of that name in its library and silently keeps its own (it generates many 4-letter
  phospholipid names, e.g. MYPC). Rename your lipid in the itp and the structure file.
- **membranes overlap**: with several `coby.membranes`, give each a `cz` (nm from the box centre)
  and make the box tall enough for both bilayers plus water between them.

## gromacs (validation)
- **atom-name mismatches**: the protein in `cg.pdb` does not match its itp (often a
  Go-model / virtual-site issue or a re-ordered file). Rebuild from `--from martinize`.
- **net charge != 0**: COBY could not neutralise (e.g. ion name not in the itps). Check
  `solvent.cation/anion`.

## check (grompp)
- Read `gromacs/.check/grompp_<stage>.log`. Common causes: GROMACS < 2021 (no
  C-rescale), a duplicated moleculetype from an `extra_itps` file, or `#include` of a
  file that also defines `[ defaults ]`.
