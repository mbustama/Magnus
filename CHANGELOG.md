# Changelog

All notable changes to Magνs are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `magnus.globaldefs.NUFIT_GLOBAL_FITS`: best-fit standard three-flavor
  oscillation parameters from every NuFit global-analysis release, v1.0
  (2012) through v6.1 (2025), by mass ordering and by each release's
  secondary category (`with_SK`/`without_SK` for v4.0+, `LEM`/`LID` for
  v2.1, `free_fluxes_rsbl`/`huber_fluxes_no_rsbl` for v1.0-v1.3).
  Transcribed directly from the official parameter tables at
  [nu-fit.org](http://www.nu-fit.org/?q=node/12).
- `magnus.globaldefs.load_nufit_params(version, ordering, category)`: loads
  a specific release/ordering/category as a plain `{s12, s23, s13, dCP,
  D21, D31}` dict, directly usable as keyword arguments to any
  `osc_prob_3nu_*` function.
- A "When is Magνs not the right tool?" section in the README and docs,
  covering quantum decoherence, open-system/bath coupling, neutrino decay,
  and self-consistent collective oscillations.
- `magnus.matter.exp_density_profile(density_matter_central, l_scale)`: builds
  an exponential density-profile callable tagged for the new fast
  interaction-picture integrator (see below); `osc_prob_2nu_matter_exp_density`
  and `osc_prob_2nu_sun` (and their NSI counterparts) now build their profile
  through it instead of an untagged lambda, so they pick up the speed-up
  automatically. A fast, closed-form interaction-picture Magnus integrator
  (`magnus.oscprob._osc_prob_ip_exp_dispatch`/`_osc_prob_ip_exp_core`) for
  two-flavor oscillations in a genuine exponential matter profile (the Sun):
  it factors the (possibly huge, at low energy) constant vacuum phase out of
  the Magnus expansion analytically, in closed form, instead of resolving it
  slab by slab, leaving only the matter-potential envelope -- whose exact
  exponential integral is also closed-form -- to be integrated. This directly
  fixes the bottleneck described under "Fixed" below for the two-flavor Sun
  wrappers: `osc_prob_2nu_sun`/`osc_prob_2nu_matter_exp_density` (and their
  NSI counterparts) now return correct, warning-free probabilities in a
  fraction of a second across the realistic solar-neutrino energy range
  (~0.1-18 MeV), where the general slab-refinement method previously either
  took many seconds to fail or returned a silently inaccurate result. The
  fast path is unconditionally exact in the accumulated vacuum phase and in
  the density profile's shape (no local-constant or local-linear
  approximation); its only approximation is truncating the resulting
  interaction-picture Magnus series at first order, which is validated
  against `solve_ivp` across the realistic energy range and is checked at
  runtime (both via the size of the neglected term and via successive-
  refinement agreement) before being trusted, with a transparent, lossless
  fallback to the general method whenever it is not (e.g., near an MSW
  resonance, where the matter term is no longer a small perturbation on the
  vacuum splitting). Three- to five-flavor scenarios are not yet covered:
  the neglected term's coefficient grows by roughly three orders of
  magnitude going from two to three flavors (one off-diagonal mass-pair vs.
  three, each with an O(1) diagonal mixing contribution), which pushes the
  slab count needed for a certified answer far past what stays fast; those
  wrappers therefore still use the general method unconditionally, exactly
  as before -- a natural target for a follow-up (e.g., a genuine adiabatic/
  WKB treatment, which does not have this scaling problem).
- `magnus.adiabatic`: a new module implementing exactly the adiabatic/WKB
  follow-up flagged above, generalized to *any* number of flavors and *any*
  Hamiltonian (not just the 2-flavor case the interaction-picture fast path
  above is restricted to). `adiabatic_propagator(H_func, l0, l1)` computes
  the evolution operator via pure adiabatic (instantaneous-eigenbasis)
  transport -- a dynamical phase (Simpson-integrated) plus a geometric
  (Berry) phase, captured implicitly via discrete parallel transport of the
  eigenvectors, with no restriction on the accumulated phase or on the
  Hamiltonian being real. `find_resonance_candidates`/
  `find_nonadiabatic_windows(H_func, l0, l1, threshold=0.1)` locate every
  position where adiabatic transport breaks down, for any pair of levels,
  via an *exact* Hellmann-Feynman diagnostic (no eigenvector finite
  differencing, which is gauge-ambiguous): candidates are exact critical
  points of each pairwise eigenvalue gap, and the adiabaticity verdict is
  the Landau-Zener-like `gamma_jk = |<v_j|dH/dl|v_k>| / gap_jk^2`; windows
  are grown to their physical width (independent of the search-grid
  spacing) and merged when they overlap, so any number of simultaneous or
  sequential resonances are handled uniformly.
  `hybrid_propagator(H_func, l0, l1, rtol=1e-3, atol=1e-3)` composes
  adiabatic transport with an exact local Magnus patch
  (`magnus.magnus.magnus_expansion_multislab`) at each non-adiabatic
  window, stitched via the exact composition law of quantum evolution
  (so the result is exactly unitary regardless of approximation accuracy),
  and self-certifies by tightening the adiabaticity threshold and two
  internal grid densities together until two successive results agree.
  Validated against `solve_ivp` on standard and BSM (NSI) 3-, 4-, and
  5-flavor Hamiltonians, and on synthetic cases with two independent and
  two merging resonances: 0-window cases match to ~2e-4 at 3,600x-25,800x
  the speed of `solve_ivp`; 1-2-window (patched) cases match to
  9e-4-2.9e-3 at 30x-91x the speed. See `docs/source/adiabatic_strategy.rst`
  for the full derivation and validation.
- `strategy` parameter (`'auto'` (default), `'hybrid'`, or `'magnus'`) on
  `osc_prob_matter_std_potential`, `osc_prob_matter_nsi`, `osc_prob_liv`,
  and every wrapper built on them (every `osc_prob_{2,3,4,5}nu_sun`,
  `osc_prob_{2,3,4,5}nu_sun_nsi`, and `osc_prob_{2,3,4,5}nu_sun_liv`
  function explicitly; every other wrapper via `**kwargs`), and also on
  `osc_prob_sun`/`osc_prob_earth` (the fully generic, arbitrary-Hamiltonian
  entry points, via a new `_osc_prob_hybrid_dispatch_generic`/
  `_osc_prob_with_potential` code path, since these bypass
  `osc_prob_matter_std_potential`/`_nsi`/`osc_prob_liv` entirely).
  `'magnus'` reproduces the exact pre-0.11.0 behavior; `'hybrid'`
  additionally tries `magnus.adiabatic.hybrid_propagator` for any
  position-dependent, breakpoint-free potential with a requested
  tolerance, returning a best-effort result plus the new
  `HybridCertificationWarning` if it fails to self-certify; `'auto'` tries
  `'hybrid'` first and falls back silently to the `'magnus'` strategies
  otherwise. Unlike the interaction-picture fast path, this applies to any
  number of flavors and does not require a tagged exponential profile. For
  `osc_prob_earth`, `'hybrid'`/`'auto'` almost always fall back to
  `'magnus'` in practice, since the PREM density profile's layer-boundary
  breakpoints are essentially always non-empty for a real trajectory --
  the general breakpoint-free requirement above, not a special case.
- `magnus.oscprob.HybridCertificationWarning` (subclasses
  `ToleranceNotAchievedWarning`): raised only when `strategy='hybrid'` is
  explicitly requested and the hybrid propagator fails to self-certify for
  at least one requested point.
- `tests/test_adiabatic.py`: unit tests for `magnus.adiabatic` (unitarity,
  resonance-candidate detection, window growth/merging in both directions,
  and `solve_ivp` cross-checks for 3-5 flavor Hamiltonians), plus new
  regression tests in `tests/test_oscprob.py`:
  `test_sun_2nu_default_strategy_avoids_tolerance_cap` demonstrates the fix
  directly through `osc_prob_2nu_sun` (the same (energy, baseline) point
  that `test_tolerance_cap_warns` shows still hits the refinement caps
  under `strategy='magnus'` is resolved, warning-free and matching
  `solve_ivp`, under the new default);
  `test_generic_osc_prob_sun_hybrid_strategy_resolves_hard_case` confirms
  the same fix through the fully generic `osc_prob_sun` entry point (a
  separate code path, `_osc_prob_hybrid_dispatch_generic`, not exercised
  by the other test); `test_generic_osc_prob_earth_strategy_falls_back_to_magnus`
  confirms `osc_prob_earth` is unaffected (PREM breakpoints disable the
  hybrid dispatch).
- `notebooks/11_magnus_adiabatic_hybrid_strategy.ipynb`: live comparison of
  `strategy='auto'`/`'hybrid'`/`'magnus'` for 2- through 5-flavor
  Hamiltonians (standard oscillations and an engineered BSM/NSI
  resonance), each cross-checked against a tight-tolerance `solve_ivp`
  ground truth in both runtime and accuracy, plus a real-data plot of the
  instantaneous eigenvalues and detected non-adiabatic window for the 3ν
  BSM case. Reproduces, live, the validation described in
  `docs/source/adiabatic_strategy.rst`.

### Changed

- Restructured the package to be flatter: `earth`, `globaldefs`, `magnus`
  (the numerical core), `matter`, `oscprob` (the main wrapper API), and
  `oscprobstd` are now flat sibling modules directly under `src/magnus/`,
  instead of each living inside its own single-file subpackage directory.
  `hamiltonians` remains a genuine subpackage (four distinct
  flavor-count-specific modules). Every `__init__.py` (top-level and
  `hamiltonians/`) now uses explicit, named imports and a hand-written
  `__all__` instead of `from .module import *` plus a `dir()`-computed
  `__all__`, which had been silently leaking implementation-detail names
  (`np`, `Optional`, `Callable`, ...) into the public namespace.
- `magnus.authors` and `magnus.version` are no longer part of the public
  API surface (excluded from `__all__` and from the Sphinx autoapi-generated
  docs), though they remain importable internally for the CLI's
  `--version` flag and `oscprob`'s banner-printing.

### Fixed

- ~47 docstring/type-annotation mismatches across the codebase, found by
  a systematic audit comparing every function's actual return/parameter
  type against what its docstring documented (not just the
  `numpy.ndarray`-documented-as-`list` case originally flagged in
  `hamiltonians{2,3,4,5}nu.py` and `oscprobstd.py`, but every such
  disagreement anywhere in `src/magnus`). Also fixes two unrelated
  copy-paste errors caught along the way in `oscprobstd.J()`'s docstring
  (a mislabeled parameter description and a wrong worked example).

### Removed

- `docs/source/sandbox/`, an untracked, unused pydata-theme experiment
  directory.

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
