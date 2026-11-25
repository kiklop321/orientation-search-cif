# User guide – step by step

This guide walks you through the tool from the first run to your own data.
For the method and every setting in detail see `README.md`.

**What the tool does.** You give it one 2D diffraction frame (Rigaku
R-AXIS RAPID `.img`) of a crystal that rotated during the exposure and a CIF
file of the phase. It finds the spots on the frame and the orientation of
the crystal that produces them, tells you which reflection (hkl) each spot
is, refines the lattice parameters, and optionally finds further grains and
the beam centre.

---

## 1. Before you start

* MATLAB R2022b or newer with the Signal Processing, Statistics and Machine
  Learning and Aerospace toolboxes.
* Open MATLAB and go to the folder:
  ```matlab
  cd 'C:\...\orientation-search-cif'
  ```
* Check that everything works:
  ```matlab
  test_cif          % a few seconds; the last line must be  ALL TESTS true
  ```

## 2. First run: the example

This example analyses `data/Ti20Nb0O.img` (α″ martensite) with
`cif/alpha2_Ti20Nb_Cmcm.cif`. Check that section 0 of
`run_orientation_search.m` has these settings (the file may have been
changed for another frame since):

```matlab
name        = 'Ti20Nb0O_alpha2_single';
image_file  = fullfile(here, 'data', 'Ti20Nb0O.img');
cif_file    = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');
phase       = phase_from_cif(cif_file, 'hkl_max', 6);      % no 'cell'
S.n_grains  = 1;
```

then run

```matlab
run_orientation_search
```

It takes about a minute. The command window shows the steps:

```
Phase: alpha'' Ti-20Nb (C m c m), single crystal      <- what was read from the CIF
Ti20Nb0O_alpha2_single: 495 spots found               <- spots on the frame
Best orientation q = [...], phi0 = 333.12 deg, 49 spots matched (39 s)
cell refinement 1: a b c = 3.1642 4.8596 4.6542 A ...  <- lattice parameters fitted to the spots
...
Saved to ...\results\Ti20Nb0O_alpha2_single
```

The figure `results/Ti20Nb0O_alpha2_single/orientation_overlay.png` shows:

| Mark | Meaning |
|---|---|
| red + | spot found on the frame |
| black o | reflection predicted by the orientation **and** a spot is there |
| green o | reflection predicted, but no spot within the tolerance |
| cyan + | beam centre |

A good solution: most black circles sit on red crosses, and the green ones
are mostly weak / high-angle reflections. (This frame contains 12 α″
variants, so one orientation explains only part of the spots – that is
expected here.)

## 3. Second example: β grains, one per frame

```matlab
test_beta_Ti
```

Runs the β-Ti CIF on four frames of β Ti–Nb (`data/beta_single_grain/`) and
on the α″ frame as a negative control, and prints a table:

```
frame        spots  matched_spots  fraction  phi0   ...  a_refined
Ti28Nb0O      101       70          0.69     79.5         3.2953
...
Ti20Nb0O      495       31          0.06    159.1         3.2815   <- wrong phase: only 6 %
```

The fraction of matched spots is the main quality check: 60–90 % for the
right phase, a few % for the wrong one.

## 4. Your own frame – step by step

All settings are in **section 0** of `run_orientation_search.m`. Open it and
go through it from top to bottom.

**4.1 The frame.** Copy the `.img` into `data/` and set
```matlab
name       = 'MySample';                              % output goes to results/MySample/
image_file = fullfile(here, 'data', 'MySample.img');
```

**4.2 The phase.** Get a CIF, e.g. from the Crystallography Open Database
(https://www.crystallography.net, search by elements), save it to `cif/`
and set
```matlab
cif_file = fullfile(here, 'cif', 'my_phase.cif');
phase = phase_from_cif(cif_file, 'hkl_max', 6);
```
If you know the lattice parameters of your sample, add them – the search is
then more reliable:
```matlab
phase = phase_from_cif(cif_file, 'hkl_max', 6, 'cell', [3.30 3.30 3.30 90 90 90]);
```

**4.3 The instrument.**
```matlab
S.wavelength = 0.7093;          % Mo K-alpha; Cu K-alpha = 1.5406
beam_center  = [1274, 3247];    % [row, column] of the direct beam
```
The sample–detector distance is read from the `.img` header. The beam
centre in the header is only nominal; [1274, 3247] is calibrated for our
RAPID. For another instrument or setup either determine it once (the
setting that makes the 2θ peaks narrowest) or let the tool refine it on a
good single-crystal frame:
```matlab
S.refine_beam_center = true;    % start must be within ~10 px; see README, step 6
```
then write the refined value (in `summary.txt`) into `beam_center` and
switch the refinement off again.

**4.4 The rotation.**
```matlab
S.rotation_range = [-30 30];    % the sample rotated by 60 deg during the exposure
S.phi0_range     = [0 360];     % where the centre of the rotation is searched
```
Set `rotation_range` to the oscillation of your measurement. Keep
`phi0_range = [0 360]`. φ0 is not the goniometer angle of your measurement:
the search looks only at one symmetry-equivalent copy of each orientation,
and that copy has its own φ0, so a range narrowed around the goniometer value
can exclude the solution. Narrowing is safe only around a φ0 found before
for the same crystal (e.g. another frame of the same grain).

**4.5 Run** `run_orientation_search` and look at the overlay and at
`summary.txt`.

**4.6 Is the result good?**

| Check | Good | Suspicious |
|---|---|---|
| matched / found spots | > 50 % for one crystal | < 15 % |
| `d2theta` residual (summary) | ± 0.03–0.1° | > 0.2° |
| refined cell | close to the CIF (< 2 %) | runs away |
| overlay | black circles on red crosses all over the frame | matches only at low angles |

## 5. Worked example: a Ti-24Nb frame (β phase)

This is section 4 done on a frame that is already in `data/`, so you can
compare your output with the one below. Open `run_orientation_search.m`
and set these lines in section 0 – everything else stays as in section 2:

```matlab
name        = 'Ti24Nb_1';
image_file  = fullfile(here, 'data', 'beta_single_grain', 'Ti24Nb_1.img');
cif_file    = fullfile(here, 'cif', 'beta_Ti_COD9012924.cif');                % beta-Ti (Im-3m)
phase       = phase_from_cif(cif_file, 'hkl_max', 6);
S.n_grains  = 1;
```

Why these: Ti-24Nb is β (bcc), so the phase is β-Ti. The CIF is pure
β-Ti (a = 3.282 Å); the alloy has a slightly larger cell, which the
lattice refinement (`S.refine_cell = 10`, already on) finds by itself.
The wavelength, beam centre and rotation range are the same as for all
frames in `data/`.

Run `run_orientation_search` (about 1 minute). **What to expect:**

```
Phase: Titanium (I m -3 m), single crystal
Ti24Nb_1: 62 spots found
Best orientation q = [0.97394 0.11049 -0.19809 0.00000], phi0 = 82.12 deg, 43 spots matched
cell refinement 1: a b c = 3.2979 ... rms d2theta 0.212 -> 0.047 deg
cell refinement 2: a b c = 3.2992 ... rms d2theta 0.035 -> 0.029 deg
cell refinement 3: a b c = 3.2992 ...                            <- converged
```

and in `results/Ti24Nb_1/summary.txt`:

| | Value | What it tells you |
|---|---|---|
| spots found | 62 | |
| matched | 57 of 62 (92 %); 43 before the cell refinement | one β crystal explains almost everything; the refined cell added 14 matches |
| refined a | 3.2992 ± 0.0003 Å | the lattice parameter of the alloy (CIF: 3.282 Å) |
| φ0 | 82.4° | centre of the rotation at which the crystal has this orientation |
| residuals | Δ2θ −0.002 ± 0.029°, Δβ 0.05 ± 0.10° | predicted and measured spots agree to a few hundredths of a degree |
| unexplained spots | 5 | weak spots, noise |

In the overlay almost every red cross carries a black circle; the green
circles are reflections at the edges of the rotation range or too weak to
be found. **This is what a good result looks like.**

**The same with the other Ti-24Nb frame** – the one with two grains:
```matlab
name       = 'Ti24Nb0O';
image_file = fullfile(here, 'data', 'two_grains', 'Ti24Nb0O.img');
```
With `S.n_grains = 1` the result is still a clean β orientation (76 spots,
a = 3.3084 Å), but **103 of 179 spots stay unexplained** – much more than
the 5 above. That is the sign of another crystal in the beam. Set
```matlab
S.n_grains = 3;     % the MAXIMUM number of grains, not the number you will get
```
The search finds grain 1, then searches again on the spots that grain 1
does not explain, and so on; it stops at the first grain that explains fewer
than `S.min_spots_per_grain` (15) spots, or after `S.n_grains` grains. So 3
only means "up to three" – with 10 the result here would be the same.
Run again (about 2 minutes). Expect:
```
Grain 2: search in the 103 unexplained spots
Grain 2: q = [...], phi0 = 32.62 deg, 38 spots, misorientation to grain 1 = 43.7 deg
Grain 3: best orientation explains only 7 spots (< 15) - no further grain
```
i.e. two β grains, 43.7° apart; a third grain is not supported by the data
(7 spots is the level of chance coincidences). The overlay shows grain 1 in
black and grain 2 in magenta. Section 6 explains how to look at the
relationship of the two grains.

**The other β frames** work the same way – change only `name` and
`image_file`:

| Frame (`data/beta_single_grain/`) | Spots | Matched | a (Å) |
|---|---|---|---|
| `Ti24Nb_1.img` | 62 | 57 | 3.2992 |
| `Ti24Nb0O-50-40.img` | 162 | 92 | 3.2965 |
| `Ti28Nb_1.img` | 74 | 65 | 3.3004 |
| `Ti28Nb0O.img` | 101 | 70 | 3.2953 |

(`test_beta_Ti.m` runs all of them at once.) If instead you keep the α″
CIF on a β frame, or the β CIF on the α″ frame `Ti20Nb0O.img`, only a few
per cent of the spots match – the quickest way to see that the phase is
wrong.

## 6. Several grains in the beam

If many spots stay unexplained, there may be more grains. Set
```matlab
S.n_grains = 3;                 % look for AT MOST three grains (stops earlier if no more are found)
S.min_spots_per_grain = 15;     % a grain must explain at least this many spots
```
Each grain gets its own colour in the overlay and its number in the
`grain` column of `spots_matched.csv`. The ready-made example is
```matlab
example_two_grains              % data/two_grains/Ti24Nb0O.img, finds 2 beta grains
```

**How are the two grains related?** Afterwards run
```matlab
visualize_two_grains
```
(set `result_dir` at its top to your results folder). It produces
`grain_relationship.png` – both unit cells, pole figures {100}, {110},
{111} and an inverse pole figure, as in EBSD – and `relationship.txt` with
the misorientation, the nearest CSL boundary (Σ3 twin etc.) and the planes /
directions that are parallel in both grains, compared with what two random
grains would give. Remember: diffraction gives the orientations, not the
position of the boundary; a common plane is only a candidate for the
boundary / habit plane.

### Martensite variants in a parent grain (MTEX)

(Everything in detail – all settings, the table, the 3D pictures and the
viewer, ω and hexagonal α′: **[GUIDE_VARIANTS.md](GUIDE_VARIANTS.md)**.)

If the frame shows martensite (e.g. α″) formed in one parent grain (β), all
its variants can be searched together; the result is the orientation of the
parent grain. The variants are made once with MTEX:

1. In `make_variants_mtex.m` set `mtex_path`, the two CIF files and the
   orientation relationship as a pair of parallel planes and a pair of
   parallel directions:
   ```matlab
   OR.parent_plane = [1 1 0];    OR.child_plane = [0 0 1];     % (110)beta || (001)alpha''
   OR.parent_dir   = [1 -1 1];   OR.child_dir   = [1 1 0];     % [1-11]beta || [110]alpha''
   ```
2. Run it. It prints a table of the variants – which parent direction each
   child axis lies along, Euler angles in the parent frame, misorientation to
   variant 1, a check of the OR – and the misorientations between the
   variants, shows them in pole figures of the parent and saves
   `variants/alpha2_in_beta.mat` (+ `_variants.csv`, `_polefigures.png`).
   With `parent_orientation` set (e.g. the `orientation.mat` of a search, or
   Euler angles of the parent grain) the table also gives every variant's
   Euler angles in the sample frame.
   For hexagonal martensite α′ use the α-Ti CIF and the Burgers OR in
   4-index notation (commented in the script):
   `OR.child_plane = [0 0 0 1]; OR.child_dir = [2 -1 -1 0];`
   For ω: the ω-Ti CIF, `OR.parent_plane = [1 1 1]; OR.child_plane = [0 0 0 1];
   OR.parent_dir = [1 -1 0]; OR.child_dir = [1 1 -2 0];` (4 variants).
   The script also draws the variants in 3D inside the parent cell: an
   overview (c axes of all variants), every variant with its pair of parallel
   planes (parent grey, child coloured) and directions (parent black, child
   dashed) and their indices, and a viewer window – choose the variant in the
   menu, turn it with the mouse.
3. In `run_orientation_search.m` use
   ```matlab
   cif_file = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');
   phase = phase_from_cif(cif_file, 'variants', fullfile(here, 'variants', 'alpha2_in_beta.mat'));
   ```
   Ti20Nb0O then gives one β orientation that explains 392 of the 495 spots
   with all 12 variants (one α″ orientation alone: 49).

## 7. What is in the results folder

| File | Content |
|---|---|
| `orientation_overlay.png` | frame with found and predicted spots |
| `summary.txt` | orientation (q, φ0), refined cell with errors, beam centre, residuals, unexplained spots |
| `spots_matched.csv` | one row per matched reflection: grain, h k l, predicted / measured 2θ and β, their differences, spot number |
| `spots_measured.csv` | all spots found, [2θ, β] in degrees |
| `orientation.mat` | everything for further work in MATLAB |

To use the orientation further in MATLAB:
```matlab
r  = load('results/MySample/orientation.mat');
g  = r.grains_saved(1);
Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
L  = Rz(g.phi0) * quat2dcm(g.q)';     % crystal axes -> instrument (x = beam, z = rotation axis)
L * [1; 0; 0]                         % where the crystal [100] points
```

## 8. When it does not work

**Find the step where it goes wrong: `step_by_step.m`.** It does the same as
`run_orientation_search.m`, but in sections that you run one by one with
Ctrl+Enter ("Run Section"); each draws its intermediate result and its
comment says what to look for:

| Section | Figure | Check |
|---|---|---|
| 1 frame | frame with the Debye rings of the phase | spots lie on the rings? (cell, beam centre) |
| 2 spots | found spots, threshold vs. background histogram | every visible spot has a +, no + on noise? |
| 2b one spot | counts around one spot, threshold line | why a weak spot is missed |
| 2c threshold | number of spots vs. `P.min_height`, on / off the rings | which threshold keeps real spots and drops noise |
| 3 2θ | 2θ of the spots against the allowed 2θ (no orientation yet) | phase, cell and geometry right? fraction on rings vs. chance |
| 4 search | predicted and measured spots in (2θ, β); best scores of the search | clearly one best orientation? |
| 5 residuals | Δ2θ, Δβ against 2θ and β | trend with 2θ → cell, sine in β → beam centre |
| 6 cell | residuals before / after the lattice refinement | |
| 7 final | the overlay of `run_orientation_search` | |

Example (Ti24Nb_1): the threshold 1.7 = 237 counts finds 62 spots, 58 of
them on an allowed 2θ (94 %, by chance 16 %); figure 2c shows that a lower
threshold adds mainly off-ring spots (1.5: 118 spots, 31 off the rings).

| Problem | What to try |
|---|---|
| few spots found / many false spots | change `P.min_height` (lower = more spots) and `P.min_prominence_*` |
| few matches, overlay shifted in 2θ | wrong cell: give `'cell'` in `phase_from_cif`; keep `S.refine_cell = 10` |
| matches only in one part of the frame | beam centre off: `S.refine_beam_center = true` |
| few matches everywhere | wrong phase, or wrong `S.rotation_range`; try a larger `S.tolerance` |
| several crystals in the beam | `S.n_grains > 1` |
| too slow | narrow `S.phi0_range`; `S.mesh_degree = [6 9]`; keep `'hkl_max', 6` |
| CIF error / warning | read the message: missing symmetry operations are taken from the space-group symbol; unknown elements are ignored |

## 9. Which file does what

| File | You run it to ... |
|---|---|
| `run_orientation_search.m` | analyse one frame (**the main script**) |
| `step_by_step.m` | the same step by step, to find where a problem comes from |
| `example_two_grains.m` | see the multi-grain search |
| `animate_orientation.m` | see how the vector n, the quaternion q and φ0 turn the unit cell (figure + GIF) |
| `visualize_orientation.m` | see in 3D what the found orientation (q, φ0) means in the experiment |
| `visualize_two_grains.m` | show the relationship of two grains |
| `show_variants_3d.m` | look at saved variants in 3D: choose a variant, turn it with the mouse (no MTEX needed) |
| `make_variants_mtex.m` | make the variants of a martensite in its parent (two CIFs + OR) with MTEX, for the search |
| `test_beta_Ti.m` | see the tool on four β frames and a wrong-phase control |
| `test_cif.m`, `test_cif_any.m` | check the CIF reader (after changes to the code) |
| `src/` | the functions – no need to open them for normal use |
