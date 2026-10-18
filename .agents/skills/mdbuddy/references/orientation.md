# Placing a protein in the membrane without OPM data

OPM provides the membrane position of deposited structures. Engineered, designed or predicted
(AlphaFold/co-folding) proteins have none. MDBuddy offers four ways to place them
(`protein.orientation`). All of them end in the same frame: membrane normal = +z,
membrane centre at z = 0, cytoplasmic side at −z (the `lower` leaflet).

| Mode | What it does | Cost | When |
|---|---|---|---|
| `hydrophobic` | implicit-membrane scan of exposed hydrophobicity | seconds | default first guess |
| `manual` | user-defined normal/centre, axis alignment, rotations, shifts | none | known topology, or adjusting the guess |
| `self_assembly` | protein + random lipids + water, the bilayer forms in MD | ~1 µs CG run | uncertain placement, let the lipids decide |
| `as_is` | trust the input frame (`membrane_center_z` in Å) | none | already oriented, e.g. by the PPM server |

The `orient` stage writes `01_orient/protein_oriented.pdb`, `orient_report.json` (method, rotation,
origin, normal and centre in the input frame) and `orientation_view.html` (interactive check).

## Hydrophobic guess (`orientation: hydrophobic`, `mdbuddy orient guess`)

The model is a simple implicit membrane, the same idea as PPM/OPM but much simpler:

- Each standard residue gets the water→octanol transfer free energy of the Wimley–White scale
  (Wimley, Creamer & White, *Biochemistry* 1996; charged D/E/K/R, neutral His). This is weighted by
  its relative solvent accessibility, from Shrake–Rupley SASA divided by the theoretical maximum
  of Tien et al., *PLoS ONE* 2013.
- A residue contributes when its side-chain centre is inside a hydrophobic slab of thickness D.
  The slab edge is a sigmoid of width 1.5 Å.
- All membrane normals (1,000 directions, then refinement within 10°), depths and thicknesses
  (20–40 Å) are scanned. The lowest total transfer energy wins.
- **Pore correction:** in β-barrels and wide channels, charged residues lining the water-filled pore
  are "exposed" and cancel the hydrophobic belt. If the plain scan finds no belt, MDBuddy repeats it
  counting only residues within 10 Å of the protein's convex hull, which are the lipid-facing ones.
  `pore_correction: on|off|auto` (default auto).
- **In/out:** the positive-inside rule. Lys/Arg within 15 Å outside each hydrophobic boundary are
  counted, and the richer side is placed at −z (cytoplasm). `topology: as_found | flip` overrides it.
- **Output:** verdict, hydrophobic thickness, transfer energy, orientation spread (largest angle
  between near-optimal solutions within 1 kcal/mol), TM segments and Lys/Arg counts.

### Benchmark against OPM (`tools/benchmark_orientation.py`)

14 OPM structures were randomly rotated and translated, then placed with the default settings.
"Tilt error" is the angle between the recovered and the OPM membrane normal. "Centre error" is
the offset along the normal. "In/out" says whether the cytoplasmic side matches OPM's assignment.

| PDB | Protein | Fold | Tilt error | Centre error | D found / OPM (Å) | In/out | Spread |
|---|---|---|---|---|---|---|---|
| 1bl8 | KcsA K⁺ channel | α | 0.7° | 0.0 Å | 36 / 32.4 | no | 7° |
| 1c3w | bacteriorhodopsin | α | 0.8° | 0.4 Å | 32 / 30.0 | yes | 6° |
| 1j4n | aquaporin-1 | α | 0.2° | 0.2 Å | 27 / 31.8 | yes | 3° |
| 1okc | ADP/ATP carrier | α | 3.6° | 1.6 Å | 25 / 29.6 | yes | 11° |
| 1ots | ClC transporter (with Fab) | α | 1.1° | 1.0 Å | 21 / 29.8 | yes | 8° |
| 1py6 | – | α | 25.6° | 0.3 Å | 33 / 29.6 | yes | 15° |
| 1u19 | rhodopsin | α | 1.5° | 0.3 Å | 33 / 31.2 | yes | 12° |
| 2bg9 | nicotinic acetylcholine receptor | α | 2.2° | 1.2 Å | 30 / 30.2 | no | 4° |
| 2nwl | GltPh transporter | α | 0.5° | 0.1 Å | 25 / 29.8 | yes | 2° |
| 2omf | OmpF porin (trimer) | β | 0.5° | 3.0 Å | 21 / 25.2 | no | 4° |
| 2rh1 | β2-adrenergic receptor | α | 0.2° | 0.2 Å | 31 / 31.8 | yes | 5° |
| 3emn | VDAC | β | **failed** | – | – | – | – |
| 3sn6 | β2AR–Gs complex | α | 5.9° | 0.3 Å | 30 / 31.0 | yes | 9° |
| 4hyj | – | α | 12.3° | 0.3 Å | 27 / 30.0 | yes | 22° |

Summary: 13 of 14 placed. The normal is within 6° for 11 of them, and the membrane centre is within
3 Å for all 13. The in/out side agrees with OPM for 10 of 13. KcsA and the nAChR have only a few,
evenly split Lys/Arg near the interface. VDAC, a wide single β-barrel with a charged pore and a
thin membrane, is not found. The thickness is a parameter of the crude model, not a prediction of
bilayer thickness. COBY builds the bilayer from the lipids themselves.

**What this means for you:** treat the result as a good starting point, never as ground truth.
Always look at it in the viewer. Report the in/out confidence (the Lys/Arg counts) to the user.
For proteins like VDAC, or when the verdict is "not a transmembrane protein", use manual placement
or self-assembly.

## Manual placement (`orientation: manual`, `mdbuddy orient view`)

`mdbuddy orient view protein.pdb --html view.html` writes one self-contained HTML page.
3Dmol.js is loaded from jsDelivr, so internet is needed. The page shows the protein coloured by
residue type (hydrophobic, Lys/Arg, Asp/Glu, Trp/Tyr) and the hydrophobic slab. It has sliders for
tilt, spin, depth and thickness, and starts from the hydrophobic guess (`--start input` for the raw
frame). It prints the spec lines, which MDBuddy applies to exactly the same geometry:

```yaml
protein:
  orientation: manual
  manual:
    normal: [0.0029, 0.5096, 0.8604]          # membrane normal, input frame, towards the outside
    center_point: [0.031, -2.348, -4.094]      # a point on the membrane midplane, Angstrom
```

Alternatives (applied in order): `align: ["A:10-12", "A:30-32"]` makes the vector from the first
selection's centroid to the second's the normal (CHARMM-GUI "align a vector along z").
`rotate: [rx, ry, rz]` adds rotations in degrees (x, then y, then z). `center: "A:15-35"` puts that
selection at the midplane, and `shift_z: 2.0` finally moves the protein along the normal.

## Self-assembly (`orientation: self_assembly`)

This is the coarse-grained self-assembly approach (Scott et al., *Structure* 2008, Sansom lab; used
at scale by MemProtMD, Stansfeld et al., *Structure* 2015).

**Build.** The protein is first aligned with the hydrophobic guess when one is found
(`self_assembly.align: as_is` keeps the input frame). It is placed in a box with lipids at random
positions and orientations, then filled with water and ions. COBY builds the soup: one correctly
named lipid of each type is cut from a tiny COBY bilayer, imported as a "solute" and flooded into
the box. The box is sized as follows:
- **x/y:** protein width + 2 × `membrane.lipid_padding`.
- **Lipid number:** enough to cover the x/y area minus the protein's cross-section at `membrane.apl`.
- **z:** room for the protein, the lipids and `self_assembly.water_per_lipid` CG waters per lipid
  (default 25 beads ≈ 100 waters).

**Protocol.** Soft-core EM → EM → 0.5 ns isotropic NPT (10 fs) → `step7_self_assembly` (20 fs,
isotropic C-rescale, `self_assembly.run_ns`, default 1 µs).
- The protein backbone is restrained throughout (`self_assembly.restrain_protein`), so the lipids
  assemble around a fixed protein.
- Pressure coupling is isotropic because the bilayer's orientation is only decided by the assembly.

**Check.** After the run, `mdbuddy selfassembly check <build_dir>` reads the last frame and:
1. **Tests each box axis for a complete periodic bilayer:** at least 80% of tail ends within
   1.2 nm of the midplane, at least 80% of head groups 1.2–3.2 nm from it, and at least 90% in-plane
   coverage.
2. **Reports the bilayer:** thickness (P–P), lipids per leaflet, whether the protein spans the
   bilayer, and its backbone RMSD to the start.
3. **Maps the orientation back:** the restrained protein is fitted back to its starting frame, and
   the bilayer normal and midplane are mapped to the frame of the original input. The in/out side
   uses the positive-inside rule.
4. **Prints `orientation: manual` lines.** Rebuild with them using COBY for a clean, equilibrated
   starting structure with **any** composition, including asymmetric leaflets.

**Which lipids can self-assemble.** MDBuddy checks the composition before building:

| | Lipids | Why |
|---|---|---|
| refused | LPS, lipid A, Re-LPS (RAMP/REMP), LOS, O-antigen, peptidoglycan, gangliosides (GM1/GM3) | huge, cation-bridged sugar head groups and very slow dynamics; they do not form a bilayer on CG timescales. Build that leaflet with COBY |
| refused | asymmetric `upper`/`lower` | lipids end up randomly in both leaflets |
| refused | > 50% glycerides, fatty acids, hydrocarbons | non-bilayer formers |
| warning | > 10% lysolipids | single tails favour micelles |
| warning | > 50% PE | cone shape favours non-lamellar phases |
| warning | > 40% cholesterol, > 25% cardiolipin | may stall or not form a bilayer |
| warning | PIPs | end up in both leaflets |
| warning | > 50% saturated lipids below 320 K | assemble above the phase transition (DPPC: 323 K) |

**Limitations.**
- The bilayer must form across the periodic box. Sometimes a vesicle, stack or ribbon forms
  instead, and `check` then reports "no complete bilayer".
- A strongly tilted restrained protein can end up in a bilayer that the box shape forces into a
  different plane. Aligning first with the hydrophobic guess (the default) avoids most of this.
- Self-assembled membranes are symmetric and well mixed, and contain no specific lipid–protein
  interactions beyond what forms during the run.
