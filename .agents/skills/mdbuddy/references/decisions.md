# Decision guide

Background for the choices in a spec. Use it to explain defaults to the user, not to
override what they ask for.

## Force field

| | Martini 3 (`martini3`) | Martini 2.2 (`martini22`) | ElNeDyn 2.2 (`elnedyn22`) |
|---|---|---|---|
| martinize2 `-ff` | `martini3001` | `martini22` | `elnedyn22` |
| martinize2 side-chain fix | default (on) | `-noscfix` | `-noscfix` |
| COBY parameter library | `default` | `M2` | `M2` |
| Elastic defaults (ef / el / eu) | 700 / 0.5 / 0.9 | 500 / 0.5 / 0.9 | 500 / 0.0 / 0.9 |
| Go model | yes | no | no |
| Lipid files shipped (pinned) | Martini 3 lipidome 2025 (PC/PE/PS/PG/PA/PI/PIPs/CL/SM/CER/glycerides…), CHOL; or the 2021 release lipids | core FF, ions, POPE, POPG, CDL0 | same as Martini 2.2 |

Martini 3 is the default for new work. Use Martini 2.2 when the user is extending or
comparing with earlier Martini 2 studies. Never mix the two in one system.

Martini 2 lipids beyond POPE/POPG/CDL0 (POPC, DPPC, CHOL, PIP2, ...) are distributed by
the Martini developers (cgmartini.nl). Ask the user for the itp file(s) and add them to
`forcefield_files.extra_itps`. MDBuddy records them with checksums, marked unpinned.
COBY 1.0.15 cannot place Martini 2 cardiolipin, cholesterol or PIPs at all.

The lipidome (Pedersen et al. 2025) is the default Martini 3 lipid set: it refines the 2021
lipids and adds hundreds more. Its names differ from Martini 2 / the 2021 release in places
(see `lipids.md`). Use `martini3_lipids: release2021` only to reproduce earlier Martini 3 work.

## Protein model

- **Elastic network** (default): harmonic bonds between backbone beads within
  `lower`–`upper` nm. It keeps the fold, but it also suppresses large conformational
  changes. With the defaults, `merge_chains: auto` makes an oligomer one molecule so the
  network spans the subunit interfaces. Otherwise subunits of a channel can drift apart.
  Set `merge_chains: none` only if the subunits should be independent.
- **GoMartini 3** (`model: go`): Lennard-Jones contacts via virtual sites, softer
  than an elastic network, and better for flexibility and conformational change. MDBuddy
  adds `#define GO_VIRT` to `topol.top` and hooks `go_atomtypes.itp` /
  `go_nbparams.itp` into `martini_v3.0.0.itp` (recorded in the lock file).
  The default `epsilon` is 9.414 kJ/mol.
- **none**: secondary structure only. Rarely what you want for a folded membrane protein.

ElNeDyn (`elnedyn22`) is itself an elastic-network force field, so `model: elastic` is the only
sensible choice for it. Checked with martinize2 0.15.0 on KcsA: `martini22 -elastic` writes 2161
rubber bands as bond function type 6 (no exclusions generated), while `elnedyn22 -elastic` writes
2138 of them as ordinary type-1 harmonic bonds, together with ElNeDyn's own heavier use of
constraints (744 vs 560) and fewer angles. Both are correct for their force field; the bands are
there either way, so do not look only for type 6 when checking an ElNeDyn topology.

Both Martini 2 force fields need martinize2's `-noscfix`. Since martinize2 0.11 the side-chain
angle/torsion corrections are applied by default, `martini22` and `elnedyn22` do not define them,
and the resulting warning stops the run at `-maxwarn 0`. MDBuddy passes `-noscfix` for them
automatically (`spec.FORCEFIELDS[...]["scfix"]`).

Secondary structure comes from DSSP, applied to the input structure. It is fixed for
the whole simulation in Martini, so the input conformation matters.

## Lipids and leaflets

- Leaflet `upper` is +z, `lower` is -z. Which biological side +z corresponds to depends
  on the OPM entry: check the topology ("in"/"out" side) on the entry page and tell the
  user. Do not assume. For asymmetric membranes this decides which leaflet gets which lipids.
- Ratios are relative, so `{POPC: 7, CHOL: 3}` means 70:30. COBY works out lipid
  numbers from the area per lipid (`apl`) and removes lipids where the protein is.
- A simple asymmetric plasma-membrane-like start (Martini 3 names):
  upper `{POPC: 7, CHOL: 3}`, lower `{POPE: 5, POPS: 2, CHOL: 3}`. Lipids without a pinned
  file (e.g. sphingomyelin, PIP2) need `extra_itps`.
- Bacterial inner-membrane-like (Martini 2 pinned files): `{POPE: 7, POPG: 2, CDL0: 1}`.
- Run above the phase transition of the saturated lipids. DPPC-rich membranes are
  usually simulated at 323 K.

## Box

`box: auto` gives a square x/y box of protein width + 2 × `lipid_padding` (3 nm), and
a height of protein height + 2 × `water_padding` (2.5 nm), at least `min_box_z`
(10 nm), rounded up to 0.5 nm. For oligomerisation or lipid-protein studies, increase
`lipid_padding`. Periodic images of the protein should be at least ~4 nm apart.

## Equilibration protocol (fixed in this version)

| Step | Integrator | dt | Backbone restraint (kJ/mol/nm²) | Steps |
|---|---|---|---|---|
| 6.0 | steepest descent, -DFLEXIBLE (harmonic bonds; Martini 3 only, see mdp_parameters.md) | – | 1000 | 5000 |
| 6.1 | steepest descent, constraints | – | 1000 | 5000 |
| 6.2 | md, C-rescale, gen-vel | 2 fs | 1000 | 25 000 (50 ps) |
| 6.3 | md, C-rescale | 5 fs | 500 | 25 000 |
| 6.4 | md, C-rescale | 10 fs | 250 | 25 000 |
| 6.5 | md, C-rescale | 15 fs | 100 | 25 000 |
| 6.6 | md, C-rescale | 20 fs | 50 | 25 000 |
| 6.7 | md, C-rescale | 20 fs | 10 | 250 000 (5 ns) |
| 7 | md, Parrinello-Rahman (tau 12 ps) | 20 fs | none | production_ns, output every output_ps |

The parameters follow the lab's Martini mdp files. Every line, and every deviation from those files, is explained in
`mdp_parameters.md`. Restraints use martinize2's `POSRES` / `POSRES_FC` blocks on the backbone (`-p backbone`).
Lipids are not restrained.
