# Changelog

All notable changes to Magνs are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.10.0] - 2026-07-30

This is the first version with a maintained changelog, and serves as the
project's baseline "first fully documented and tested" release: the entries
below summarize a single, large audit-and-modernization pass rather than
day-to-day incremental changes.

### Added

- Command-line calculator (`magnus prob`, also runnable as `python -m magnus`)
  for computing a single oscillation probability from the shell, covering
  vacuum, matter (constant/exponential density), Earth, and Sun, with
  standard, NSI, and LIV scenarios, for 2-5 flavors. See `docs/source/cli.rst`.
- Full `pytest` suite (`tests/`) and GitHub Actions CI: `tests.yml` (matrix
  across Python 3.10-3.12), `lint.yml` (Ruff), `pages.yml` (Sphinx ->
  GitHub Pages), `publish.yml` (PyPI on release).
- Gauss-Legendre commutator-free integrators (`integration_method='gl'`),
  silent vectorization of Hamiltonian/density-profile evaluation, an
  energy-batched scan engine for separable Hamiltonians, adaptive slab
  refinement with warm-starting across scan points, and slab edges aligned
  with PREM layer boundaries.
- `osc_prob_earth` and `osc_prob_sun`: generic entry points that accept an
  arbitrary user-supplied Hamiltonian in the Earth/Sun environments.
- Sphinx documentation (published to GitHub Pages) covering installation,
  quick start, the CLI, the full function listing, code architecture,
  methodology, and this changelog.
- Root-cause regression tests for every bug listed under Fixed below.

### Changed

- Collapsed roughly 1,150 lines of duplicated refinement/logging
  keyword-argument declarations across ~60 wrapper functions into a single
  source of truth (internally called the "G1" refactor), with a permanent
  test (`test_no_wrapper_redeclares_standard_refinement_kwargs`) guarding
  against the pattern recurring.
- Rewrote the Magnus-expansion numerical core: corrected higher-order term
  coefficients, added order-6 support, and restructured the term recursion
  and matrix exponentiation for batched/vectorized evaluation.
- Package layout consolidated under `src/magnus/` (src-layout) with proper
  `__init__.py` files; version metadata unified to a single source of truth.

### Fixed

- Two mixing-matrix formula bugs (`mixing_matrix_4x4` and `mixing_matrix_5x5`)
  that invalidated every sterile-neutrino (3+1, 3+2) calculation.
- `hamiltonian_2nu_nsi`'s `eps_aa` parameter was a silent no-op: it sat on
  both diagonal entries, making it a pure multiple of the identity (an
  unobservable global phase) with zero effect on any probability.
- `osc_prob_5nu_matter_nsi_exp_density` called the non-NSI inner function
  instead of the NSI one.
- A sign error in `hamiltonian_2nu_liv_energy_independent`'s off-diagonal
  term.
- Several `_nsi_td`/`_liv_td` position-dependent convenience functions
  crashed with `TypeError` on any call.
- `unpack_nsi_params_from_dict`/`unpack_liv_params_from_dict` silently
  returned `None` instead of raising for unsupported flavor counts.
- Missing `nubar` parameter on several matter/LIV wrapper functions.
- A boundary bug incorrectly rejecting `rho_central == 0.0` in the NSI
  exponential-density wrappers.
- `hamiltonian_3nu_liv` crashed with `TypeError` on every call (it forwarded
  an incomplete argument list to its own energy-independent helper).
- Dead/unreachable code and stale, copy-pasted docstrings (wrong flavor
  count, wrong matrix dimensions, description of the wrong scenario) across
  `hamiltonians{2,3,4,5}nu.py` and `oscprob.py`, found while writing complete
  docstrings for every function.
- NumPy 2.0 compatibility (removed deprecated type aliases).
