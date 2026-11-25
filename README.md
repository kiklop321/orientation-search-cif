# Orientation search from a single 2D X-ray diffraction frame – CIF version

MATLAB tool that finds the Bragg spots on an image-plate frame of a rotating
crystal and determines the orientation of a crystal, **defined by a CIF
file**, that reproduces them. Any phase of any of the 230 space groups can
be used: the reciprocal lattice and the allowed reflections are built from
the cell, the symmetry operations and the atoms in the CIF.

Method: K. Salata et al., *The effect of Nb and O on the martensitic
transformation in Ti–Nb–O alloys*, J. Appl. Cryst. 59 (2026),
doi:[10.1107/S1600576726007739](https://doi.org/10.1107/S1600576726007739)
(Sections 3.1 and 4.3). Please cite it when you use the tool (`CITATION.cff`).

**Data.** The diffraction frames (`.img`) used by the examples and tests are
not part of this repository; see [`data/README.md`](data/README.md) for the
files and where they go. The code, the CIF files and the saved variants are
complete without them; `tests/test_cif.m`, `tests/test_cif_any.m`,
`tests/test_variants_mtex.m` and `show_variants_3d.m` run without any frame.

## Getting started in five steps

1. **Get the code:** "Code → Download ZIP" on GitHub (or `git clone`), unzip,
   open the folder `orientation-search-cif` in MATLAB (R2022b or newer with the
   Signal Processing, Statistics and Machine Learning and Aerospace toolboxes).
2. **Check it works:** `run tests/test_cif` – the last line must be
   `ALL TESTS true` (no data needed).
3. **Frames:** put the example frames into `data/` as listed in
   [`data/README.md`](data/README.md), or use your own `.img` frame.
4. **Settings:** open `run_orientation_search.m`; everything you set is in
   **section 0** – the frame (`image_file`), the phase (`cif_file`), the
   instrument (`S.wavelength`, `beam_center`), the rotation during the
   exposure (`S.rotation_range`) and how many grains to look for
   (`S.n_grains`). Every line has a comment saying what it means.
5. **Run** it (F5). The result is in `results/<name>/`: `orientation_overlay.png`,
   `summary.txt`, `spots_matched.csv`. If it does not look right, go through
   `step_by_step.m` section by section (Ctrl+Enter) – every figure says what
   to look for.

## Quick start

**New to the tool? Start with [GUIDE.md](GUIDE.md)** – a step-by-step guide from the
first run to your own frame, how to read the results and what to do when it
does not work. Variants of martensite / ω in a parent grain (MTEX) and the
pictures of an orientation: [GUIDE_VARIANTS.md](GUIDE_VARIANTS.md).

```matlab
cd orientation-search-cif
run tests/test_cif          % checks the CIF reader on the files in cif/ (seconds)
run_orientation_search      % orientation search, settings in section 0
```

To use another phase, put its CIF into `cif/` and set `cif_file` in section 0
of `run_orientation_search.m`.

## Folder structure

```
GUIDE.md                     step-by-step user guide
GUIDE_VARIANTS.md            guide to the variants with MTEX and to the orientation pictures
run_orientation_search.m     main script: settings (section 0) and the five steps
step_by_step.m               the same, one step at a time with a check figure after each (for finding problems)
make_variants_mtex.m         variants of a child phase in a parent crystal from two CIFs and an OR, with MTEX
show_variants_3d.m           open the 3D pictures and the viewer of saved variants (no MTEX needed)
visualization/               pictures of results
  animate_orientation.m      how n -> z (quaternion) and phi0 turn the unit cell: 6 stages + animation
  visualize_orientation.m    what q and phi0 mean: 3D picture of the cell, beam, rotation axis and image plate
  visualize_two_grains.m     orientation relationship of two grains (cells, pole figures, IPF, CSL)
examples/
  example_two_grains.m       search for several grains in one frame (two-grain Ti-24Nb)
tests/                       checks of the code (run after changes)
  test_cif.m                 CIF reader on the Ti phases, variants
  test_cif_any.m             CIF files of all crystal systems (cif/test_set/)
  test_beta_Ti.m             beta-Ti CIF on four single-grain Ti-24Nb / Ti-28Nb frames and the Ti-20Nb control
  test_variants_mtex.m       MTEX variants against variants_from_OR (needs MTEX)
variants/                    variants saved by make_variants_mtex.m (.mat for the search, _variants.csv,
                               _polefigures.png, _3d_overview/_3d_variants/_3d_variant1.png):
                               alpha2_in_beta (alpha'' in beta), alpha_hcp_in_beta (alpha' in beta, Burgers),
                               omega_in_beta (omega in beta)
cif/                         CIF files (from COD, plus alpha'' Ti-20Nb written by hand)
src/
  read_cif.m                 CIF reader: cell, symmetry operations, atoms
  spacegroup_ops.m           symmetry operations from a Hall / H-M symbol or number
  symmetry_data/syminfo.lib  CCP4 table of all 230 space groups in 541 settings
  expand_positions.m         all atoms of the unit cell from the asymmetric unit
  point_group_rotations.m    proper rotations of the crystal (Cartesian)
  phase_from_cif.m           phase for the search from a CIF file
  phase_reflections.m        allowed reflections and their reciprocal vectors
  orientation_pipeline.m     steps 1-5 for one frame and one phase
  element_Z.m                atomic numbers
  find_bragg_peaks.m         spot search on the frame
  search_orientation.m       two-stage orientation search
  orientation_grid.m         trial orientations in the fundamental region
  closest_quaternions.m      finer trial orientations around a solution
  match_spots.m              predicted spots of an orientation and their measured partners
  refine_cell.m              lattice parameters fitted to the 2theta of the matched spots
  refine_beam_center.m       optional: direct-beam position fitted to the matched spots
  plot_overlay.m             frame with measured and predicted spots
  describe_figure.m          explanatory text below a figure (step_by_step.m)
  ewald_angles.m             diffraction geometry of the rotating crystal
  calculate_reciprocal_lattice.m, short_quat.m
  variants_from_OR.m         optional: variants of a product phase in a parent crystal
  variants_mtex.m            optional: the same with MTEX (used by make_variants_mtex.m)
  plot_variants_3d.m         the variants in 3D in the parent cell, with the parallel planes / directions; viewer
external/                    Rigaku image reader (Z. Matěj), icosahedral mesh (A. Semechko)
data/                        example frames
  Ti20Nb0O.img               Ti-20Nb: alpha'' martensite, 12 variants of one beta grain
  beta_single_grain/         Ti-24Nb and Ti-28Nb, beta, one grain each:
    Ti28Nb0O.img             frame of the article (2023-03-22)
    Ti28Nb_1.img, Ti24Nb_1.img   2022-10-24
    Ti24Nb0O-50-40.img       2023-04-05, tube at 50 kV/40 mA
  two_grains/Ti24Nb0O.img    Ti-24Nb frame of the article: two beta grains
results/                     created by the script
```

## From the CIF to the reciprocal lattice

1. **Read** (`read_cif.m`): cell (a, b, c, α, β, γ), the symmetry operations
   (`_space_group_symop_operation_xyz` or `_symmetry_equiv_pos_as_xyz`, e.g.
   `-x+1/2, y, -z`), and the atom sites (label, element, fractional x y z,
   occupancy). Numbers with uncertainties (`3.1599(2)`), quoted values,
   text fields and element symbols with charges (`O2-`) are handled; several
   atoms on one site (e.g. Ti 0.8 / Nb 0.2) are kept with their occupancies.
   **CIF files without symmetry operations** (only a space-group symbol) are
   completed from the CCP4 table `syminfo.lib` (`spacegroup_ops.m`): by the
   Hall symbol, else by the Hermann–Mauguin symbol (`F d -3 m :1`,
   `R -3 c :H`, `P 1 21/c 1`, ...), else by the number. For rhombohedral
   groups the setting (:H / :R) follows from the cell; for groups with two
   origins and no choice given, origin 2 is used with a warning; a symbol not
   in the table (non-standard setting) falls back to the standard setting of
   the number with a warning. `cif.symmetry_source` records what was used.
2. **Expand** (`expand_positions.m`): every atom is transformed by every
   symmetry operation and mapped into the unit cell; duplicates are removed.
   Coordinates written with few digits (0.3333, 0.6667, ...) are set to the
   nearest multiple of 1/24 when closer than 5·10⁻⁴, otherwise forbidden
   reflections would not cancel exactly.
3. **Reflections** (`phase_reflections.m`): for all (h k l) up to `hkl_max`
   that can reach the Ewald sphere, F = Σ occ·Z·exp(2πi(hx + ky + lz)) is
   computed with the atomic number Z as scattering factor; reflections with
   |F|² ≈ 0 (systematic extinctions) are removed. The reciprocal vectors
   g = (h k l)·B, with B the reciprocal basis of the cell, are what is
   matched to the measured spots. Only positions are used, so no tables of
   atomic form factors are needed.
4. **Symmetry of the search** (`point_group_rotations.m`, `orientation_grid.m`):
   orientations related by a rotation of the crystal give the same pattern,
   so only the fundamental region of the crystal's rotation group is
   searched (1/24 of all directions for cubic, 1/12 hexagonal, 1/8
   tetragonal, 1/4 orthorhombic, 1/2 monoclinic, all for triclinic).

5. **Lattice parameters** (`refine_cell.m`): the cell of a CIF rarely fits
   the sample exactly (e.g. pure β-Ti vs. a Ti–Nb alloy), which shifts the
   predicted 2θ and loses matches at high angles. 2θ depends only on the
   d-spacing, i.e. only on the cell, so after the search the cell is fitted
   by least squares to the measured 2θ of the matched spots (1 parameter for
   a = b = c, 2 for a = b ≠ c, 3 otherwise; angles kept). The reflections are
   recomputed and the orientation refined around the previous one; this is
   repeated until the cell no longer changes (`S.refine_cell` = maximum number
   of rounds, 0 = keep the CIF cell). The refined cell and its standard
   deviation are written to `summary.txt`.

6. **Beam centre** (optional, `S.refine_beam_center = true`,
   `refine_beam_center.m`): the direct beam is hidden by the beam stop, and
   the value in the .img header ("Direct beam position") is the nominal one
   of the instrument, not a calibration. Without a calibrant the matched
   spots of the crystal itself fix the centre: their pixels are converted
   to (2θ, β) with a trial centre and compared with the predicted positions
   (free: row, column, and a small tilt of the crystal about x and y as
   nuisance parameters; a rotation about z, the rotation axis, changes only
   the angle φ at which a reflection diffracts, not its 2θ or β, so it cannot
   be fitted; residuals Δ2θ and Δβ·sin 2θ). A shift of the centre moves the
   spots in 2θ differently for different β, a change of the cell moves all
   of them alike, so both can be refined together; they are refined in
   turn in the loop of step 5. In the first round a coarse scan of
   ±`S.beam_center_scan` px (default 20) picks the start with the most
   matched spots. The sample–detector distance (127.26 mm, the radius of
   the cylindrical image plate, sample in its axis) is read from the header
   and not refined: it cannot be separated from the lattice parameters.

   Test on the four β frames: started from 1274, 3247 (found by making the
   2θ peaks as narrow as possible) all converge to
   [row, column] = [1275.0 ± 0.1, 3247.1 ± 0.2] (esd of one frame
   0.2–0.3 px), lattice parameter unchanged within 0.0001 Å; the α″ frame
   Ti20Nb0O (other sample, other phase, Cmcm CIF) gives [1275.1, 3247.2].
   Started from
   the header value [1287, 3246] (13 px off; with it only 15–31 spots match
   instead of 57–93), three of the four frames reach the same point
   (± 0.01 px); for Ti24Nb_1 the header centre already leads the orientation
   search to an orientation 3° off, and the fit stops at [1283.0, 3246.4]
   with 16 matched spots. The start should
   therefore be within about 10 px; the centre is a property of the
   instrument, so it is best refined once on a good single-crystal frame
   and then kept fixed (`S.refine_beam_center = false`).

## Settings (section 0 of `run_orientation_search.m`)

| Setting | Meaning |
|---|---|
| `image_file`, `name` | input frame, name of the output folder |
| `S.wavelength`, `beam_center` | X-ray wavelength (Å), direct-beam pixel [row, column] |
| `S.rotation_range` | rotation of the sample during the exposure relative to its centre φ0, e.g. `[-30 30]` (deg) |
| `S.phi0_range`, `S.phi0_step` | where and in which steps the centre φ0 is searched |
| `S.tth_range` | only predicted spots with 2θ in this range are matched |
| `plot_max_tth` | predicted reflections are drawn in the figure up to this 2θ (default 180 = all); does not affect the search |
| `P.*` | spot search: minimum height and prominence, merging distance, refinement window |
| `S.mesh_degree`, `S.tolerance` | fineness of the trial orientations, matching tolerance in 2θ and β (stage 1, 2) |
| `cif_file` | the phase; options of `phase_from_cif`: `'cell'` (lattice parameters of the sample), `'hkl_max'`: reflections −n…n (default 6) or `'auto'` = all that can diffract at the wavelength (d ≥ λ/2; for Mo Kα up to index ≈ 2a/λ, e.g. 10 for β-Ti, 14 for α″), more high-angle reflections but a slower search |
| `S.search_region` | `phase.parent_rotations` (fundamental region) or `'octant'` |
| `S.refine_cell` | maximum number of lattice-parameter refinements (stops when converged; 0 = none) |
| `S.refine_beam_center`, `S.beam_center_scan` | refine the direct-beam position to the matched spots (default false); half-width of the coarse scan before it (px, default 20, 0 = none) |
| `S.n_grains` | maximum number of grains: 1 = one crystal; > 1 = after each grain the search is repeated on the spots not yet explained (several grains of the same phase in the beam) |
| `S.min_spots_per_grain` | a further grain is accepted only if it explains at least this many spots (default 15) |

**Special case – variants.** If the crystal is a product of a phase
transformation in a parent crystal (as α″ in β in Ti–Nb), all its variants
can be predicted together and the orientation of the parent searched:
`phase_from_cif(cif, 'parent', parent_cif, 'OR', OR)` with the orientation
relationship as two pairs of parallel directions (see `variants_from_OR.m`).

**Variants with MTEX.** The variants can also be made with MTEX
(https://mtex-toolbox.github.io, tested with 6.1.1): `make_variants_mtex.m`
takes the two CIF files and the OR written as in the literature – a pair of
parallel planes and a pair of parallel directions, e.g. (110)β ∥ (001)α″,
[1−11]β ∥ [110]α″ – or an OR predefined in MTEX (`orientation.Burgers`,
`orientation.KurdjumovSachs`, ...). It computes all variants
(`variants(p2c, parent)` in MTEX), checks the OR on every variant, shows
them in pole figures of the parent crystal and saves them to
`variants/<name>.mat`; the search then uses them by

    phase = phase_from_cif(child_cif, 'variants', 'variants/alpha2_in_beta.mat');

MTEX supplies the symmetry and the OR only; atoms and structure factors of
the child still come from its CIF. The MTEX crystal frames are built from
the space group and cell of `read_cif` with X ∥ a, Z ∥ c* (the convention of
this tool; MTEX's own CIF reader does not know the newer CIF keys).
The script is generic (any parent, child and OR, 3- or 4-index notation,
or an OR predefined in MTEX such as `'Burgers'`, `'KurdjumovSachs'`). Besides
the `.mat` file for the search it writes `variants/<name>_variants.csv`: for
every variant the parent directions along the child axes, its Bunge Euler
angles in the parent frame, its misorientation to variant 1 and a check of
the OR; with `parent_orientation` (a search result or Euler angles of the
parent grain) also its Euler angles in the sample frame (X = beam, Z =
rotation axis). The misorientations between the variants are printed.

Checks (`tests/test_variants_mtex.m` and runs on Ti20Nb0O):
* β → α″ (OR of the article) and β → α hcp (Burgers, from MTEX and in
  4-index notation): the 12 MTEX variants are exactly the 12 of
  `variants_from_OR.m`.
* The full search on Ti20Nb0O (with lattice refinement) gives with the MTEX
  variants and with the own ones the same parent orientation (0.0000° apart,
  φ0 = 160.50°), 403 of 495 spots, the same cell, and the same 12 variant
  orientations in the sample (0.0000°). The independent single-crystal α″
  search on that frame lies 0.30° from variant 5.
* Misorientations between the variants: α″ 4.36°, 86.96°, 90°, 93.12°
  (the small one: two variants on the same {110}β, because b/a of α″ is not
  √2); hexagonal α (Burgers) 10.53°, 60°, 60.83°, 63.26°, 90° – the five
  classes known for Burgers variants.
* ω in β, (111)β ∥ (0001)ω, [1−10]β ∥ [11−20]ω: 4 variants, the c axis of
  each along one ⟨111⟩β, all 70.53° apart (the angle between ⟨111⟩).
* With the hexagonal variants the Ti20Nb0O frame (α″) gives only 42 of 495
  spots (α″ variants: 392 with the CIF cell) – the search tells the phases apart.

## Several grains in one frame

With `S.n_grains > 1` (see `examples/example_two_grains.m`) the spots explained by the
grains found so far are removed and the orientation search is repeated on
the remaining spots, with the refined cell. A new grain is accepted if it
explains at least `S.min_spots_per_grain` spots; the search stops at the
first grain that does not, or after `S.n_grains` grains. In
`orientation_overlay.png` every grain has its own colour (thick circles:
predicted and found; thin diamonds: predicted but not found), `spots_matched.csv`
has a column `grain`, and `summary.txt` lists the orientation, the matched
spots and the misorientation to grain 1 of every grain.

Example – the Ti-24Nb frame of the article (`data/two_grains/Ti24Nb0O.img`):

| | Grain 1 | Grain 2 | Grain 3 |
|---|---|---|---|
| matched spots | 76 | 38 | 7 (< 15, rejected) |
| φ0 | 326.0° | 32.6° | – |
| misorientation to grain 1 | – | 43.7° | – |


### Relationship of two grains (`visualization/visualize_two_grains.m`)

After `examples/example_two_grains.m`, `visualization/visualize_two_grains.m` shows how the two
grains are related, as in EBSD (`results/Ti24Nb0O_two_grains/grain_relationship.png`
and `relationship.txt`):

* the two unit cells in the sample frame (x = beam, z = rotation axis),
  faces in the IPF colour of the beam direction, with the best common
  plane and common direction drawn in yellow in both cells;
* pole figures {100}, {110}, {111} of both grains (projection along the
  beam); poles of the two grains closer than `tol` are ringed;
* inverse pole figure of the beam and of the rotation axis (cubic phases);
* misorientation angle and axis (smallest over the crystal symmetry),
  the nearest CSL relation with the Brandon criterion (cubic), and all
  pairs of parallel planes / directions with indices ≤ `max_index` within
  `tol`, together with the number of such pairs expected for two grains
  with a random misorientation (Monte Carlo, 500 orientations) - only a
  clearly larger number means a real relationship.

Ti24Nb0O: 43.7° about ≈ ⟨652⟩, no CSL relation (nearest Σ29, 46.4° ⟨221⟩, is
7.0° off, Brandon limit 2.8°), 3 parallel plane pairs within 3° against
3.4 ± 2.4 for random grains, i.e.
an ordinary high-angle boundary with no special relation. Diffraction gives
only the orientations; a common plane is a candidate for the boundary /
habit plane, but whether the boundary lies on it must be checked in the
microstructure.

## Output (`results/<name>/`)

| File | Content |
|---|---|
| `summary.txt` | phase, settings, best orientation, score, residuals |
| `orientation_overlay.png` | frame with measured spots (red +) and predicted reflections. One grain: black o = predicted and found, green o = predicted (recorded during the rotation, on the detector) but no spot within the tolerance. Several grains: one colour per grain, thick o = found, thin ◇ = not found |
| `spots_measured.csv` | 2θ, β of all spots found |
| `spots_matched.csv` | every matched prediction: grain, h k l, φ − φ0, predicted and measured 2θ and β, differences, index of the measured spot |
| `orientation.mat` | everything above plus the search result (`result`: best orientation and, of the last search, the tested orientations and their scores - after a lattice refinement only those of the fine stage) |

## Tests

`tests/test_cif.m` reads the CIF files in `cif/`:

| File | Phase | Space group | Source |
|---|---|---|---|
| `beta_Ti_COD9012924.cif` | β-Ti | Im-3m | COD 9012924 |
| `beta_TiNb_COD1538066.cif` | β Ti–Nb (Ti/Nb mixed site) | Im-3m | COD 1538066 |
| `alpha_Ti_hcp_COD9008517.cif` | α-Ti | P6₃/mmc | COD 9008517 |
| `omega_Ti_COD9012925.cif` | ω-Ti | P6/mmm | COD 9012925 |
| `alpha2_Ti20Nb_Cmcm.cif` | α″ Ti–20Nb | Cmcm | written for this tool (lattice parameters from Table 2 of the article) |

and checks the number of atoms in the cell, the number of rotations and the
systematic extinctions of each space group, and that the 12 α″ variants
generated from the β and α″ CIFs form a set that every rotation of β maps
onto itself (a rotation of the parent only permutes the variants).

`tests/test_beta_Ti.m` runs the β-Ti CIF (a = 3.282 Å, refined to every frame) on
the four single-grain frames of the β alloys and on the α″ control:

| Frame | Spots | Matched to one β crystal | Refined a (Å) | φ0 (deg) |
|---|---|---|---|---|
| Ti28Nb0O | 101 | 70 (69 %) | 3.2953 ± 0.0002 | 79.3 |
| Ti28Nb_1 | 74 | 65 (88 %) | 3.3004 ± 0.0002 | 230.8 |
| Ti24Nb_1 | 62 | 57 (92 %) | 3.2992 ± 0.0003 | 82.4 |
| Ti24Nb0O-50-40 | 162 | 92 (57 %) | 3.2965 ± 0.0002 | 237.1 |
| Ti20Nb0O (α″, control) | 495 | 31 (6 %) | – | 159.1 |

* **One grain per frame.** After removing the spots of the first crystal, a
  second orientation search on the remaining spots finds at most 0–6 spots
  (the level of random coincidences, cf. 6 % in the α″ control). The frame
  Ti24Nb0O of the article contains two grains (see *Several grains in one
  frame*) and is kept separately in `data/two_grains/`.
* A frame of the same position taken at 30 kV/20 mA (not included) gives
  the same orientation as Ti24Nb0O-50-40 within 0.67°. At 50 kV more spots
  remain unexplained (70 vs. 16 at 30 kV), likely λ/2 harmonics passed by
  the graphite (004) reflection.
* The lattice parameter of Ti28Nb0O agrees with that from the β peaks of its
  integrated 1D profile (3.2930 ± 0.0020 Å) within 0.0023 Å. Between frames
  of the same alloy, a varies by up to 0.4 %; the absolute lattice parameter
  from one frame therefore has a systematic uncertainty of about ±0.006 Å
  (sample position, detector distance). This does not affect the orientation.
* In the α″ control the few coinciding reflections still give the
  orientation of the parent β grain found with the α″ variants (within 2°).

`tests/test_cif_any.m` checks ten CIF files of all crystal systems from COD
(`cif/test_set/`): NaCl (Fm-3m), α-Fe (Im-3m, no symmetry operations in the
file), cristobalite (Fd-3m, origin 1), Mg (P6₃/mmc), corundum (R-3c,
hexagonal axes), quartz (P3₁21), rutile (P4₂/mnm), anatase (I4₁/amd,
origin 1), baddeleyite (P2₁/c) and albite (C-1, non-standard setting). For
every file the number of proper rotations of the point group and the
absence of all reflections forbidden by the space group (centring, glide
planes, screw axes) are checked, and, where the file lists its operations,
the operations taken from the symbol are compared with them.

COD files: Crystallography Open Database, https://www.crystallography.net (CC0).
`src/symmetry_data/syminfo.lib`: CCP4 symmetry library (CCP4, LGPL), copy
from https://www.mrc-lmb.cam.ac.uk/harry/imosflm/ver721/downloads/syminfo.lib.

## Limitations

* The lattice refinement (`refine_cell.m`) fits the lengths a, b, c; the
  angles of monoclinic and triclinic cells are kept as in the CIF.
* The allowed reflections are decided with f = Z; reflections that vanish
  only accidentally for the real form factors are kept (harmless for the
  search).
* One phase per run; several grains of that phase can be found with
  `S.n_grains`, spots of other phases stay unmatched.
* The crystal symmetry used (search region, misorientation) is the proper
  rotations of the CIF point group. For non-centrosymmetric groups (e.g.
  −43m, 4mm, 3m) the diffraction pattern has the higher symmetry of the
  Laue group, so two orientations that give the same pattern may be
  reported as different (e.g. 90° apart), and `visualization/visualize_two_grains.m`
  treats 23, m-3 and −43m crystals as non-cubic (no IPF, no CSL).
* `visualization/visualize_two_grains.m` assumes a single-crystal phase; with variants
  ('parent', 'OR') the printed indices would use the child cell for the
  parent orientation.
* `S.mesh_degree`: stage 2 takes 4^(d2+2−d1) orientations; the second
  value must not exceed the first by more than about 4 (default [7 10]).
* Spot merging: of a chain of three close peaks A > B > C (A–B and B–C
  closer than `P.min_distance_px`, A–C not) both B and C are removed.
* The handedness of the laboratory frame (sign of y and of φ) follows the
  Rigaku reader and has not been checked on a crystal of known orientation.
  Misorientations do not depend on it; absolute orientations (pole figures,
  IPF, "beam along [uvw]") would be mirrored if it were wrong.

## Requirements

MATLAB R2022b (tested) with the Signal Processing, Statistics and Machine
Learning and Aerospace toolboxes. Optional: MTEX 6.1.1 for
`make_variants_mtex.m` / `tests/test_variants_mtex.m`.

## Licence and third-party code

The code of this tool is released under the MIT licence (`LICENSE`).
It includes third-party code and data under their own terms:

| Folder / file | Author | Licence |
|---|---|---|
| `external/RigakuRapidUtils/rigaku.cpp`, `rigaku.h` | Rigaku Americas Corp., adapted by Z. Matěj | Rigaku free-software licence, `external/RigakuRapidUtils/license-rigaku.txt`. This software is based in part on software developed at Rigaku. |
| `external/RigakuRapidUtils/RAPIDread.*`, `RAPIDutils.m` | Z. Matěj (Charles University) | reader of Rigaku R-AXIS RAPID frames, used with the author's courtesy |
| `external/sphere/` | A. Semechko | MIT, `external/sphere/LICENSE.md` |
| `src/symmetry_data/syminfo.lib` | CCP4 | LGPL |
| `cif/*_COD*.cif`, `cif/test_set/COD_*.cif` | Crystallography Open Database | CC0 |
| MTEX (optional, not included) | R. Hielscher et al. | GPL, https://mtex-toolbox.github.io |
