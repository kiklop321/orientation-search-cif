# Diffraction frames

Two example frames are included, so the tool can be tried at once:
`Ti20Nb0O.img` (used by `run_orientation_search.m` as delivered) and
`beta_single_grain/Ti24Nb_1.img` (used by `step_by_step.m` and the example
in `docs/USER_GUIDE.md`). The other frames (Rigaku R-AXIS RAPID `.img`,
2560 x 4600 px, 23.6 MB each) are not stored here; put them here with these
names – the scripts expect exactly this layout:

```
data/
  Ti20Nb0O.img                     Ti-20Nb: alpha'' martensite (12 variants of one beta grain)
                                   - run_orientation_search.m (as delivered),
                                     tests/test_beta_Ti.m (negative control), variant search
  beta_single_grain/               beta Ti-24Nb / Ti-28Nb, one grain each
    Ti24Nb_1.img                   - step_by_step.m, docs/USER_GUIDE.md section 2, tests/test_beta_Ti.m
    Ti24Nb0O-50-40.img             - tests/test_beta_Ti.m (tube at 50 kV / 40 mA)
    Ti28Nb_1.img                   - tests/test_beta_Ti.m
    Ti28Nb0O.img                   - tests/test_beta_Ti.m
  two_grains/
    Ti24Nb0O.img                   two beta grains - examples/example_two_grains.m, visualization/visualize_two_grains.m
```

Instrument: Rigaku R-AXIS RAPID, Mo K-alpha (0.7093 A), cylindrical image
plate of radius 127.26 mm (read from the header), direct beam at pixel
[row, column] = [1274, 3247] (calibrated; see README.md, step 6).

Where to get them: the frames of the article are published with its data
on Zenodo, doi:10.5281/zenodo.19087995; the others are available from the
author on request.

Your own frames can be anywhere: set `image_file` in section 0 of
`run_orientation_search.m` / `step_by_step.m`.
