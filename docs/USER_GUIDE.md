# User guide

Start with the five steps in the [README](../README.md). This guide adds
what you need when you go further. Details of the method: [METHOD.md](METHOD.md).

## 1. Your own frame

Everything is in **section 0** of `run_orientation_search.m`:

| Setting | What to put | Example |
|---|---|---|
| `name` | name of the output folder `results/<name>/` | `'MySample'` |
| `image_file` | your `.img` frame | `fullfile(here, 'data', 'MySample.img')` |
| `cif_file` | CIF of the phase (e.g. from https://www.crystallography.net) | `fullfile(here, 'cif', 'beta_Ti_COD9012924.cif')` |
| `phase = phase_from_cif(...)` | add `'cell', [a b c α β γ]` if you know the lattice parameters of your sample | `phase_from_cif(cif_file, 'hkl_max', 6, 'cell', [3.30 3.30 3.30 90 90 90])` |
| `S.wavelength` | X-ray wavelength (Å) | `0.7093` (Mo Kα) |
| `beam_center` | pixel [row, column] of the direct beam; calibrated for our RAPID | `[1274, 3247]` |
| `S.rotation_range` | rotation of the sample during the exposure | `[-30 30]` |
| `S.n_grains` | **maximum** number of grains to look for | `1` |

Leave the rest as it is. Run (F5) and check:

| Check (`summary.txt`, overlay) | Good | Suspicious |
|---|---|---|
| matched / found spots | > 50 % for one crystal | < 15 % (wrong phase?) |
| `d2theta` residual | ± 0.03–0.1° | > 0.2° |
| overlay | black circles on red crosses over the whole frame | matches only at low angles |

## 2. Example to compare with: Ti-24Nb (β)

In section 0 set
```matlab
name       = 'Ti24Nb_1';
image_file = fullfile(here, 'data', 'beta_single_grain', 'Ti24Nb_1.img');
cif_file   = fullfile(here, 'cif', 'beta_Ti_COD9012924.cif');      % beta-Ti
```
and run. Expected: 62 spots found, **57 matched**, a = 3.2992 Å (the CIF of
pure β-Ti has 3.282 Å – the tool refines it), residuals Δ2θ ≈ 0.03°.

## 3. Several grains

If many spots stay unexplained, there may be another grain: set
`S.n_grains = 3` (up to three; the search stops when a further grain explains
fewer than 15 spots). Example: `data/two_grains/Ti24Nb0O.img` gives two β
grains (76 and 38 spots, 43.7° apart):
```matlab
run examples/example_two_grains
run visualization/visualize_two_grains     % pole figures, misorientation, CSL
```
Diffraction gives the orientations, not the boundary between the grains.

## 4. Martensite / ω in a parent grain (variants, MTEX)

A transformed grain contains several **variants** of the product (α″ in β:
12, α′: 12, ω: 4). The tool can predict all of them for one parent
orientation and so find the orientation of the **parent** grain.

1. Install MTEX (https://mtex-toolbox.github.io, unzip) and set `mtex_path`
   in `make_variants_mtex.m`.
2. Set the two CIF files and the orientation relationship – one parallel
   plane pair and one parallel direction pair (the direction must lie in the
   plane):

   | Child | `child_cif` | `OR.parent_plane / child_plane` | `OR.parent_dir / child_dir` |
   |---|---|---|---|
   | α″ (orthorhombic) | `alpha2_Ti20Nb_Cmcm.cif` | `[1 1 0]` / `[0 0 1]` | `[1 -1 1]` / `[1 1 0]` |
   | α′ (hexagonal, Burgers) | `alpha_Ti_hcp_COD9008517.cif` | `[1 1 0]` / `[0 0 0 1]` | `[1 -1 1]` / `[2 -1 -1 0]` |
   | ω (hexagonal) | `omega_Ti_COD9012925.cif` | `[1 1 1]` / `[0 0 0 1]` | `[1 -1 0]` / `[1 1 -2 0]` |

   (all three are ready in the script; or `OR = 'Burgers'`, `'KurdjumovSachs'`, … predefined in MTEX)
3. Run it. You get a table of the variants (`variants/<name>_variants.csv`),
   pole figures and 3D pictures, and a window where you choose a variant and
   turn it with the mouse. Saved variants can be opened again without MTEX:
   `show_variants_3d.m`.
4. Search with the variants – in `run_orientation_search.m`:
   ```matlab
   cif_file = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');
   phase = phase_from_cif(cif_file, 'variants', fullfile(here, 'variants', 'alpha2_in_beta.mat'));
   ```
   Example: `data/Ti20Nb0O.img` – one β orientation explains about 400 of the
   495 spots with all 12 α″ variants (one α″ crystal alone: 49).

## 5. What the orientation means

The orientation is a quaternion **q** (it tilts the crystal so that one
crystal direction n lies on the rotation axis z) and an angle **φ0** (then a
turn about z). To see it:

* `visualization/animate_orientation.m` – the cell turned by q, then by φ0 (figure + animation)
* `visualization/visualize_orientation.m` – the cell with the beam and the detector; prints the Euler angles (Bunge)
* `step_by_step.m`, section 8 – the same for the orientation just found

## 6. When it does not work

Run **`step_by_step.m`** section by section (Ctrl+Enter). Every figure has
a text under it saying what to look for.

| Problem | What to try |
|---|---|
| visible spots not found / many false spots | `P.min_height` (lower = more spots); `step_by_step` section 2c shows the best value |
| spots beside the Debye rings (`step_by_step` 1, 3) | wrong lattice parameter → `'cell'` in `phase_from_cif` |
| residuals go up and down with β (`step_by_step` 5) | beam centre → `S.refine_beam_center = true` once, then copy the value into `beam_center` |
| few matches everywhere | wrong phase or `S.rotation_range` |
| many spots unexplained | another grain (`S.n_grains`) or another phase |
| `Frame not found` | the frame is not in `data/` – see `data/README.md` |
| `Set mtex_path` | `make_variants_mtex.m` needs the folder of MTEX |
