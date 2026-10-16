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

## gromacs (validation)
- **atom-name mismatches**: the protein in `cg.pdb` does not match its itp (often a
  Go-model / virtual-site issue or a re-ordered file). Rebuild from `--from martinize`.
- **net charge != 0**: COBY could not neutralise (e.g. ion name not in the itps). Check
  `solvent.cation/anion`.

## check (grompp)
- Read `gromacs/.check/grompp_check.log`. Common causes: GROMACS < 2021 (no
  C-rescale), a duplicated moleculetype from an `extra_itps` file, or `#include` of a
  file that also defines `[ defaults ]`.
