# Issue #161 comment: true error against the requested tolerance on smooth PREM chords

A comment on #161 reported that on smooth PREM chords, the true error of `osc_prob_3nu_earth`
under `strategy='magnus'` exceeded the requested tolerance by up to 1.7× in 3 of 24 combinations.
That was measured against the tester's own layer-aligned product.  `rtol_gap.py` repeats the
tester's grid against the oracle `diagnostics.rst` prescribes: DOP853 at `rtol=1e-12`, run layer
by layer, and checked by tightening to `1e-13`.

The grid:
- 3 GeV, `electron_fraction=0.5`;
- `integration_method` `'gl'`, `'trapezoid'` and `'simpson'`;
- `magnus_exp_order` 2, 4, 6 and 8;
- `rtol = atol` 1e-4 and 1e-8;
- `costhz` -0.8 and -1.

That is 48 calls.

Run from the repository root, on `main` at `106ad0c`:

    python docs/dev/measurements/issue161_rtol_gap/rtol_gap.py

`output.txt` holds the run:
- The oracle moved by 3.3e-12 and 5.8e-13 when tightened.
- None of the 48 calls is outside the tolerance.
- True error ÷ tolerance: max 0.85 (GL, order 2, 1e-8), median 0.014, min 8e-8.

The tester's three exceedances do not reproduce against this oracle.  Their sizes, 1e-8 to
2e-8 at the 1e-8 request, match the offset between the package and the tester's Earth reference
that #167 records at any tolerance.  That points to the reference rather than the stopping rule.

This backs the paragraph "On smooth profiles..." under "What `rtol` and `atol` actually
control" in `docs/source/diagnostics.rst`.
