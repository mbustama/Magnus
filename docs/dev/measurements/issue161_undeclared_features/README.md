# Issue #161: undeclared jumps and narrow spikes on the refinement ladder

Branch `fix-161-ladder-undeclared-features`, against `main` at `22dfebf`.  The scripts are kept
as they were run; they read the scratch paths of the session, so they record the protocol.

## `population.py`: what fails, and how often

Every call at the default tolerance, 1e-3, scored against DOP853 at `rtol=1e-12`, split at the
features.
- **Castle wall** (the issue's Hamiltonian, written by hand, no `t_breakpoints`), through
  `osc_prob_energy_baseline`, 40 baselines from 1000 to 12000 km: 24 of 40 were off by more than
  1e-3 on `main`, up to 7.0e-2, all with only `MagnusConvergenceWarning`.  On the branch: 0.
- **Gaussian spikes** through `osc_prob_matter_std_potential` under `'auto'`: widths 2, 5, 20 and
  50 km at 500, 1000 and 2500 km along 4000 km, at 0.5, 1 and 3 GeV.  17 of 36 were off by more
  than 1e-3 on `main`, up to 1.4e-2, with no warning.  On the branch: 3, the 2 km spike at 2500 km
  at each energy, which falls between two probes of the 200-point grid.
- The issue's first proposal (an agreement to rounding counts as unproven) caught 19 of the 41
  wrong answers on `main` and fired on 6 of 27 correct calls, so it was not adopted.

## `benchmark.py`: nothing else changes

140 calls across the public API: vacuum, constant, exponential and Earth (2-5 flavors, standard,
NSI, LIV), the Sun and its models, smooth and spiked custom profiles, raw Hamiltonians,
`average=True`, explicit strategies, tight tolerances, the evolution operator, energy and
baseline scans.  `values` mode: 126 of 140 are bit-identical to `main`, values and warning
categories alike; the other 14 are the castle wall, the spikes and a density step, and
`changed_accuracy.py` shows every one of them moved toward the reference (castle wall 4.1e-2 to
4.5e-13, 5 km spike 1.1e-2 to 4.7e-8, step 2.9e-4 to 5.0e-13, 50 km spike 5.4e-6 to 1.4e-8).

`timing` mode, interleaved rounds with the order alternated, and `main` against `main` as the
noise reference: no change on any path with declared breakpoints, on the Sun under `'auto'`, or
in vacuum and constant density (0.997 to 1.013).  A single point that reaches the per-point ladder
without breakpoints -- a raw Hamiltonian, `strategy='magnus'`, a smooth custom profile under
`'auto'` -- is 2 % to 5 % slower, about 0.05 ms (0.7 ms on a 20 ms solar ladder), the cost of
keeping the last two levels' samples to compare; the author accepted it.

This backs `UNDECLARED_JUMP_STEP_TOLERANCE`, `AUTO_SHARP_FEATURE_SLABS_PER_PROBE`, the
`UnmarkedDiscontinuityWarning` docstring and row in `diagnostics.rst`, and the #161 entry of
`CHANGELOG.md`.
