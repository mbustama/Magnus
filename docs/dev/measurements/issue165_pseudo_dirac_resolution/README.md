# Issue #165: `average=True` on a pseudo-Dirac pair against a 50-digit reference

`average=True` on a pseudo-Dirac Hamiltonian, through the generic engine, takes each pair's
phase from the difference of two eigenvalues of the Hamiltonian.  `eigh` resolves that
difference only to about 1e-16 of the largest eigenvalue, so a splitting far below the standard
ones is partly lost.  `resolution.py` measures by how much, and checks that the new route,
`oscprob.osc_prob_pseudo_dirac_vacuum` with `average=True`
(`avgprob.phase_averaged_probabilities_pseudo_dirac`), does not lose it.

The setup:
- three flavors at the NuFIT defaults, E = 1 TeV, default spread;
- one pair at a time, on each of the three mass states;
- pair phases δL/2E of 0.5, 1, 3, 10, 20 and 30 rad;
- ratios δ/max m² from 1e-4 down to 4e-16.

The reference is the phase average of the exact eigenbasis, with every phase formed in 50-digit
arithmetic (`mpmath`).  Reported: the largest error over the 3x3 active block, the states and the
phases.

Run from the repository root, on branch `fix-165-pseudo-dirac-average` off `main` at `90f4893`:

    python docs/dev/measurements/issue165_pseudo_dirac_resolution/resolution.py

`output.txt` holds the run:
- The new route is within 1.1e-16 of the reference at every ratio.
- The generic route is off by 1.9e-8 down to a ratio of 1e-8.  That floor is not round-off in
  the eigenvalues: it is the central difference in ln E (step `_PHASE_SLOPE_STEP` = 1e-3) that
  estimates the phase slopes.  The new route uses the exact slope, dφ/dln E = −φ.
- Below 1e-8 the eigenvalue round-off takes over: 2.2e-7 at 1e-10, 6.9e-5 at 1e-12, 5.9e-3 at
  1e-14 and 0.25 at 4e-16.

This backs the numbers in the docstrings of `osc_prob_pseudo_dirac_vacuum` and
`phase_averaged_probabilities_pseudo_dirac`, the #165 entry of `CHANGELOG.md`, and the
pseudo-Dirac section of `docs/source/averaged_probability.rst`.
