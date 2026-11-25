# Guide – variants of a phase in a parent crystal (MTEX) and pictures of orientations

This guide covers two related tools:

* **Part A – variants:** `make_variants_mtex.m` computes, with MTEX, all
  variants of a product phase (child, e.g. martensite α″, α′ or ω) formed in a
  parent crystal (e.g. β) from two CIF files and an orientation relationship
  (OR). It gives their orientations, shows them in pole figures and in 3D, and
  saves them for the orientation search. `show_variants_3d.m` opens the 3D
  pictures of saved variants without MTEX.
* **Part B – what an orientation means:** `visualization/animate_orientation.m`,
  `visualization/visualize_orientation.m` and section 8 of `step_by_step.m` show how the
  orientation found by the search (quaternion q and φ0) turns the unit cell
  in the experiment, and give its Euler angles.

For the orientation search itself see `GUIDE.md`; for the method `README.md`.

---

# Part A – variants with MTEX

## A1. What it is for

When a parent grain transforms (e.g. β → α″ martensite on quenching), the
product forms in several crystallographically equivalent orientations –
**variants**. The orientation relationship fixes how one variant lies in the
parent; the symmetry of the parent gives the others. For α″ in β there are
12, for α′ (Burgers) 12, for ω 4.

A diffraction frame of such a grain contains the spots of all variants
together. Instead of searching each variant separately, the orientation
search can predict all variants of one parent orientation at once and find
the orientation of the **parent** grain (`README.md`, "Special case –
variants"). `make_variants_mtex.m` makes these variants from two CIF files and
the OR written as in the literature.

## A2. What you need

* **MTEX** (free MATLAB toolbox for texture analysis):
  https://mtex-toolbox.github.io – download, unzip (tested: MTEX 6.1.1).
  You do not have to install it permanently: the script adds it to the path
  and starts it (`startup_mtex`) when it is not running yet.
* **Two CIF files:** the parent (e.g. β-Ti `cif/beta_Ti_COD9012924.cif`) and the
  child (e.g. α″ `cif/alpha2_Ti20Nb_Cmcm.cif`, α-Ti `cif/alpha_Ti_hcp_COD9008517.cif`, ω-Ti
  `cif/omega_Ti_COD9012925.cif`). MTEX uses only their space group and cell; the
  atoms of the child are used later by the search.
* **The orientation relationship:** a pair of parallel planes and a pair of
  parallel directions, e.g. {110}β ∥ (001)α″ and ⟨111⟩β ∥ ⟨110⟩α″.

## A3. Running `make_variants_mtex.m`

1. Open `make_variants_mtex.m` and go to **Settings**.
2. Set the path to MTEX (the folder that contains `startup_mtex.m`):
   ```matlab
   mtex_path = 'C:\...\mtex-6.1.1\mtex-6.1.1';
   ```
3. Set the two phases:
   ```matlab
   parent_cif = fullfile(here, 'cif', 'beta_Ti_COD9012924.cif');          % beta-Ti
   child_cif  = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');   % alpha''
   ```
4. Set the OR – **one member of each family**:
   ```matlab
   OR.parent_plane = [1 1 0];    OR.child_plane = [0 0 1];    % (110)_beta  || (001)_alpha''
   OR.parent_dir   = [1 -1 1];   OR.child_dir   = [1 1 0];    % [1-11]_beta || [110]_alpha''
   ```
   The direction must lie in the plane, in both phases: [1−11] lies in (110)
   because 1·1 + 1·(−1) + 0·1 = 0; [110] lies in (001). The other members of
   the families ({110}, ⟨111⟩) need not be written – they are exactly what
   gives the other variants.
5. Names for the tables and pictures, and the name of the output:
   ```matlab
   parent_label = 'beta';
   child_label  = 'alpha''''';        % alpha'' (two quotes inside a MATLAB string)
   name = 'alpha2_in_beta';           % -> variants/alpha2_in_beta.*
   ```
   `''` as a label takes the name and space group from the CIF.
6. Optional – the orientation of the parent grain (see A6):
   ```matlab
   parent_orientation = [];
   ```
7. Run the script (F5). It takes about half a minute (MTEX starts on the
   first run).

### Ways to give the OR

| Form | Example | When |
|---|---|---|
| plane + direction (3 indices) | `OR.parent_plane = [1 1 0]; OR.child_plane = [0 0 1]; OR.parent_dir = [1 -1 1]; OR.child_dir = [1 1 0];` | the usual way |
| plane + direction (4 indices, hexagonal) | `OR.child_plane = [0 0 0 1]; OR.child_dir = [2 -1 -1 0];` | hexagonal / trigonal child or parent; (hkil) and [UVTW], the third index must be −(first + second) |
| two direction pairs | `OR = struct('parent_dirs', [1 1 0; 1 -1 1], 'child_dirs', [0 0 1; 1 1 0]);` | as in `variants_from_OR.m` |
| predefined in MTEX | `OR = 'Burgers';` | bcc → hcp (Burgers); fcc → bcc: `'KurdjumovSachs'`, `'NishiyamaWassermann'`, `'Pitsch'`, `'Bain'`, `'GreningerTrojano'`; also `'PitschSchrader'`, `'ShojiNishiyama'` |

With a predefined OR the table has no "OR check" and the 3D pictures show no
planes/directions (the OR is not given as indices).

### Ready examples (commented in the script)

| Child | CIF | OR | Variants |
|---|---|---|---|
| α″ (Cmcm) | `alpha2_Ti20Nb_Cmcm.cif` | (110)β ∥ (001)α″, [1−11]β ∥ [110]α″ | 12 |
| α′ (P6₃/mmc), Burgers | `alpha_Ti_hcp_COD9008517.cif` | (110)β ∥ (0001)α′, [1−11]β ∥ [2−1−10]α′ | 12 |
| ω (P6/mmm) | `omega_Ti_COD9012925.cif` | (111)β ∥ (0001)ω, [1−10]β ∥ [11−20]ω | 4 |

The saved results of all three are in `variants/` (`alpha2_in_beta`,
`alpha_hcp_in_beta`, `omega_in_beta`).

## A4. What the script gives

### In the command window

```
beta -> alpha'': 12 variants
parent -> child misorientation: 45.05 deg (Euler 92.18 90.00 45.00 deg)

variant  child_001_along_parent  child_100_along_parent  phi1_parent_frame ... misorientation_to_variant_1  OR_check_deg
   1        [1 1 0]                 [0 0 1]                 135             ...        0                       0
   2        [0 -1 -1]               [-1 0 0]                0               ...       86.96                    0
   ...
Misorientations between variants (deg): 4.36, 86.96, 90.00, 93.12
```

| Column | Meaning |
|---|---|
| `variant` | number of the variant |
| `child_001_along_parent` | the parent direction along which the child [001] (c axis) lies |
| `child_100_along_parent` | the parent direction along which the child [100] (a axis) lies |
| `phi1/Phi/phi2_parent_frame` | Bunge Euler angles of the variant, with the parent crystal as the reference frame |
| `misorientation_to_variant_1` | the smallest rotation (over the child symmetry) from variant 1 to this variant (deg) |
| `OR_check_deg` | how well the variant satisfies the OR (largest angle between the given parallel pairs, allowing for the symmetry of both phases); should be 0 |
| `phi1/Phi/phi2_sample` | only with `parent_orientation`: Euler angles of the variant in the sample frame (A6) |

"Misorientations between variants" lists the different angles between all
pairs. They are characteristic of the OR and can be compared with the
literature or with EBSD: α″ in β 4.36°, 86.96°, 90°, 93.12° (4.36° = two
variants on the same {110}β; not 0 because b/a of α″ is not √2); α′ Burgers
10.53°, 60°, 60.83°, 63.26°, 90°; ω 70.53°.

### Files (folder `variants/`)

| File | Content |
|---|---|
| `<name>.mat` | the variants `V` (3×3×N, child → parent frame), the CIF paths, the OR, the table, the misorientations; **read by the search** |
| `<name>_variants.csv` | the table above (opens in Excel) |
| `<name>_polefigures.png` | pole figures (MTEX) of the child [001], [100], [010] of all variants in the parent frame; squares = parent ⟨110⟩ |
| `<name>_3d_overview.png` | parent cell with the c axes of all variants; variant 1 with its parallel planes / directions |
| `<name>_3d_variants.png` | every variant in the parent cell with its parallel planes / directions and their indices |
| `<name>_3d_variant1.png` | variant 1, large |

### Figures (stay open after the run)

1. **Pole figures (MTEX):** where the child axes of the variants point in the
   parent crystal; each variant in its colour with its number.
2. **Overview (3D):** left, the parent unit cell (grey) and, from its centre,
   the c axis of every variant (colour, with the parent direction). Right,
   variant 1 to scale inside the parent cell with the OR.
3. **All variants (3D):** a panel per variant; the title gives the members of
   the families that are parallel in this variant, e.g.
   `(0 -1 -1) beta || (0 0 1) alpha''  planes`,
   `[-1 1 -1] beta || [1 1 0] alpha''  directions`.
4. **Viewer:** one variant, large – choose it in the menu **variant:** at the
   top left, **turn it with the mouse** (rotation is on). Switching the
   variant keeps the view, so variants can be compared from the same side.

How to read the 3D pictures:

| Element | Meaning |
|---|---|
| grey cell | parent unit cell (axes a, b, c) |
| coloured cell | the child variant, **to scale** (e.g. ω is larger than β) |
| red / green / blue arrows | a / b / c axes (child: a′ b′ c′ in the large views) |
| grey section | the parent plane of the OR, cut through the cell centre |
| coloured section | the child plane of the OR, cut through the child cell |
| black line | the parent direction of the OR |
| coloured dashed line | the child direction of the OR – on the black line |

## A5. Using the variants in the orientation search

In `run_orientation_search.m` (section 0):

```matlab
cif_file = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');      % the CHILD
phase = phase_from_cif(cif_file, 'variants', fullfile(here, 'variants', 'alpha2_in_beta.mat'));
```

The parent CIF is taken from the `.mat` file. The search then predicts the
reflections of all variants and finds the orientation of the **parent**
grain; in the output `spots_matched.csv` the column `variant` tells which
variant explains each spot.

Example (Ti20Nb0O, `results/Ti20Nb0O_alpha2_variants/`): one β orientation
explains 403 of 495 spots with all 12 α″ variants; one α″ crystal alone
explains 49. The same result is obtained with the variants written by hand
in the tool (`variants_from_OR.m`) – checked to 0.0000°.

The variants can be used in `step_by_step.m` and `tests/test_beta_Ti.m` in the same
way (the `phase = ...` line).

## A6. The variants in the sample (orientation of the parent known)

If the orientation of the parent grain is known, the table also gives every
variant's Euler angles in the sample frame (X = beam, Z = rotation axis of
the sample, as in the search). In the settings:

```matlab
parent_orientation = fullfile(here, 'results', 'Ti20Nb0O_alpha2_variants', 'orientation.mat');
% or Bunge Euler angles of the parent:
parent_orientation = [111.01 16.08 49.49];
```

The `orientation.mat` can be the result of a search with these variants (it
then contains the parent orientation) or of a search of the parent phase
itself (e.g. β with the β CIF); grain 1 is used.

Check on Ti20Nb0O: the single-crystal α″ search (without variants) on the
same frame lies 0.30° from variant 5 computed this way.

## A7. Looking at saved variants: `show_variants_3d.m`

Opens the three 3D figures and the viewer of a saved `.mat` – no MTEX, no
recomputation:

```matlab
variants_file = fullfile(here, 'variants', 'alpha2_in_beta.mat');   % or omega_in_beta.mat, alpha_hcp_in_beta.mat
```

Run it (F5) and use the viewer: menu **variant:**, mouse to turn.

## A8. Checks

`tests/test_variants_mtex.m` (needs MTEX) compares the MTEX variants with those of
`variants_from_OR.m` for β → α″ (article OR), β → α′ (MTEX Burgers) and β → α′
(Burgers in 4-index notation): every MTEX variant must be exactly one of the
tool's own. Last line: `ALL TESTS true`.

## A9. When it does not work

| Message / problem | Cause and remedy |
|---|---|
| `MTEX is not on the MATLAB path` | `mtex_path` does not point to the folder with `startup_mtex.m` |
| `MTEX does not know the space group ...` | unusual space-group symbol in the CIF; check `_space_group_name_H-M_alt` / `_symmetry_space_group_name_H-M` |
| `Four-index notation needs the 3rd index = -(1st + 2nd)` | e.g. `[2 -1 -1 0]` is right, `[2 -1 0 0]` is not |
| error in `rotation.map` / OR check ≠ 0 | the direction does not lie in the plane (in one of the phases), or the angle between the two given vectors differs between parent and child |
| `MTEX has no predefined OR "..."` | misspelt name; the list is in the message |
| fewer / more variants than expected | the number follows from the symmetries: \|parent rotations\| / \|rotations common to both\|; check the space groups of the CIFs |
| pole figures empty when run with `matlab -batch` | MTEX lays them out for the window size; the script calls `drawNow`, in the normal MATLAB desktop they show |

Note: MTEX brings its own `normalize` function on the path; the tool does not
use MATLAB's `normalize`, so the search works with MTEX started.

---

# Part B – what an orientation means: pictures and Euler angles

The search describes an orientation by a quaternion **q** and an angle **φ0**:

* **n** – a direction fixed in the crystal (a unit vector: 2 independent
  parameters). The search tries directions n on the unit sphere.
* **q** – the shortest rotation that brings n to the rotation axis z of the
  sample: a tilt of the whole cell about the horizontal axis n × z (so the
  4th component of q is always 0).
* **φ0** – then a turn of the whole cell about z (1 parameter). n stays on z.

Together 3 parameters = one orientation (the orientation in the middle of the
exposure; during the exposure the sample turns φ0 ± 30°). Splitting it this
way is fast because a turn about z does not change where a reflection lands
on the detector (2θ, β), only the angle at which it diffracts.

## B1. `visualization/animate_orientation.m` – the cell turned by q and φ0

Settings: `result_dir` (a results folder with `orientation.mat`), `grain`,
`make_gif`. Shows six stages – start (a ∥ x = beam, b ∥ y, c ∥ z), the tilt by
q (n goes to z, orange arrow, its path dotted), the turn by φ0 about z (path
of a dotted), and the final cell seen along z. With `make_gif = true` the
movement is also played and saved as `.gif` and `.mp4` in `result_dir`.

## B2. `visualization/visualize_orientation.m` – the cell in the experiment

Settings: `result_dir`, `grain`, `n_rays`. Four panels: the cell inside the
cylindrical image plate with the beam, the beam centre (red dot) and the
middle of the plate (black ×); after q; after φ0; and the exposure with the
diffracted beams of the strongest matched reflections ending at their spots.
The command window gives the axis and angle of q, the crystal directions
along x, y, z and the Euler angles; the drawing is checked against the
predicted spots (difference in px).

## B3. `step_by_step.m`, section 8

The same two pictures for the orientation just found in sections 4–6. The
animation is off by default (`play_animation = false` in section 0).

## B4. Euler angles

`src/orientation_to_euler.m` converts (q, φ0) to Bunge Euler angles
(φ1, Φ, φ2) as used in EBSD:

1. q → rotation matrix `quat2dcm(q)`;
2. with φ0: `g = (Rz(φ0) * quat2dcm(q)')'` (sample → crystal);
3. Φ = acos g₃₃, φ1 = atan2(g₃₁, −g₃₂), φ2 = atan2(g₁₃, g₂₃) (`src/euler_bunge.m`).

Sample frame: X = beam, Z = rotation axis. In this parametrisation Φ and φ2
come from q alone and φ0 only adds to φ1. EBSD software usually uses another
sample frame (Z = normal of the sample surface); relate the frames before
comparing. The crystal symmetry gives equivalent sets (24 for cubic).

---

## Which file does what

| File | Use |
|---|---|
| `make_variants_mtex.m` | variants of a child phase in a parent from two CIFs and an OR (MTEX): table, pole figures, 3D, `.mat` for the search |
| `show_variants_3d.m` | 3D pictures and viewer of saved variants (no MTEX) |
| `tests/test_variants_mtex.m` | MTEX variants against the tool's own (needs MTEX) |
| `visualization/animate_orientation.m` | how q and φ0 turn the cell (figure, GIF/MP4) |
| `visualization/visualize_orientation.m` | the cell, beam and detector in 3D; Euler angles |
| `src/variants_mtex.m` | the MTEX part: crystal symmetries from the CIFs, OR, variants |
| `src/plot_variants_3d.m` | the 3D pictures and the viewer |
| `src/plot_orientation_steps.m`, `src/plot_orientation_3d.m` | the orientation pictures (used by the scripts above and `step_by_step.m`) |
| `src/orientation_to_euler.m`, `src/euler_bunge.m` | Euler angles |
| `variants/` | saved variants: α″, α′ (Burgers), ω in β |
