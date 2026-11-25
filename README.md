# Orientation search from a single 2D X-ray diffraction frame

MATLAB tool: from one frame of a rotating crystal (Rigaku R-AXIS RAPID `.img`)
and a **CIF file** of the phase it finds the orientation of the crystal,
tells which reflection (hkl) each spot is and refines the lattice
parameters. It can also find several grains, and the variants of a
martensite / ω phase in its parent grain (with MTEX).

Method: K. Salata et al., J. Appl. Cryst. 59 (2026),
doi:[10.1107/S1600576726007739](https://doi.org/10.1107/S1600576726007739) – please cite it.

## Getting started

1. **Download** ("Code → Download ZIP" or `git clone`) and open the folder
   in MATLAB R2022b+ with the Signal Processing, Statistics and Machine
   Learning and Aerospace toolboxes.
2. **Check:** `run tests/test_cif` → last line `ALL TESTS true`.
3. **Frames:** two example frames are included (`data/`), so you can run
   the scripts at once; more in [`data/README.md`](data/README.md), or use
   your own.
4. **Settings:** open `run_orientation_search.m`; all you set is in
   **section 0** (frame, CIF, wavelength, beam centre, rotation range,
   number of grains) – every line has a comment.
5. **Run** (F5). Results: `results/<name>/`.

## What do you want to do?

| I want to … | Run |
|---|---|
| find the orientation of a crystal on my frame | **`run_orientation_search.m`** |
| see step by step why the result is wrong (spots, threshold, cell, beam centre) | `step_by_step.m` (section by section, Ctrl+Enter) |
| make the variants of a martensite / ω in a parent grain (two CIFs + OR) | `make_variants_mtex.m` (needs MTEX) |
| look at saved variants in 3D (turn them with the mouse) | `show_variants_3d.m` |
| see what the found orientation means (cell, beam, detector) | `visualization/visualize_orientation.m`, `visualization/animate_orientation.m` |
| compare two grains (pole figures, misorientation, CSL) | `visualization/visualize_two_grains.m` |

## Results (`results/<name>/`)

| File | Content |
|---|---|
| `orientation_overlay.png` | the frame with found spots (red +) and predicted reflections (black o = matched, green o = predicted but not found) |
| `summary.txt` | orientation, refined cell, residuals, unexplained spots |
| `spots_matched.csv` | every matched reflection: h k l, predicted and measured 2θ, β |
| `orientation.mat` | everything, for further work in MATLAB |

## More

* **[docs/USER_GUIDE.md](docs/USER_GUIDE.md)** – how to use it: your own
  frame, a worked example, several grains, variants with MTEX, when it does
  not work.
* **[docs/METHOD.md](docs/METHOD.md)** – how it works, all settings, tests,
  limitations, third-party code.

## Licence

MIT (`LICENSE`). Third-party code in `external/` and `src/symmetry_data/` is
under its own licence (see `docs/METHOD.md`). This software is based in part
on software developed at Rigaku.
