# Changelog

All notable changes to Magνs are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [1.0.0rc1] - 2026-07-31

First public release candidate.  Magνs was developed privately up to this
point, so everything below is new to anyone outside the project.  The entries
are still grouped as Added/Changed/Fixed/Removed, describing the development
history that produced this release -- "Changed" and "Fixed" are relative to
earlier private states of the code, not to any published version -- because
that history is the most useful record of *why* the code looks the way it does.

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
- A short pronunciation note at the top of the README and of the docs
  landing page: Magνs is said just like "Magnus", with the neutrino symbol
  **ν** standing in for the "nu" syllable (with a nod to the Danish
  pronunciation, given where most of it was written).
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
  two merging resonances: 0-window cases match to ~2e-4 at 3,600x-4,800x
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
  `'magnus'` reproduces the exact behavior from before the adiabatic
  strategy was added; `'hybrid'`
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

- `max_n_slabs` is now method-aware: it defaults to None, meaning "use the cap
  that suits `integration_method`" -- 20000 for `'gl'`, 2000 for the
  cumulative-quadrature methods (`magnus.oscprob.MAX_N_SLABS_DEFAULT`). An
  explicit value is always used as given, so this changes nothing for callers
  who set it. A single cap could not serve both families: `'gl'` costs 1-3
  Hamiltonian evaluations per slab against the quadrature methods'
  `n_tpts_per_slab`, so at a shared cap of 2000 it hit the ceiling on problems
  it could resolve comfortably -- eV-scale sterile splittings over an
  Earth-crossing baseline need about 8,600 slabs -- and emitted
  `ToleranceNotAchievedWarning` on answers that were roughly 1,600x *more*
  accurate than the quadrature methods achieved within that same cap. Even at
  20000, `'gl'` remains the cheaper worst case (40,000-60,000 evaluations
  against the ~200,000 that 2000 quadrature slabs at 100 points already
  permit).
- **Breaking:** `integration_method` now defaults to `'gl'`
  (Gauss-Legendre commutator-free collocation) instead of `'trapezoid'`,
  everywhere it appears -- the 13 signature defaults across `magnus` and
  `oscprob`, and `magnus prob --integration-method`. For a Hamiltonian that
  is smooth within each slab, which layer-aligned slab edges make the common
  case, `'gl'` is simultaneously the faster and the more accurate choice: it
  needs only 1, 2, or 3 Hamiltonian evaluations per slab for orders <= 2,
  <= 4, <= 6, with its quadrature order matched exactly to the truncation
  order, where the cumulative-quadrature methods sample `n_tpts_per_slab`
  points and can let quadrature error dominate the truncation error. Numbers
  computed without passing `integration_method` explicitly will therefore
  change slightly, and two further behaviors switch on with it, both of which
  were already implemented and gated on `'gl'`: `n_tpts_per_slab` no longer
  participates in the adaptive refinement (accuracy is set by the slab count
  alone), and the physics-informed starting slab count from
  `magnus.suggest_n_slabs` is now applied by default. `'trapezoid'` and
  `'simpson'` remain fully supported and are the better choice when the
  Hamiltonian has a kink or a discontinuity *inside* a slab, where
  Gauss-Legendre loses its order advantage.

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
- Docstring equations that were still written as indented plain text now
  render as real LaTeX. The commit that converted ~600 instances of
  plain-text *symbol* notation to `:math:` roles left the bare *equations*
  alone; `magnus.py`'s module docstring (the Magnus expansion and its
  Bernoulli-number recursion), `MagnusConvergenceWarning`'s convergence
  criterion, `earth.prem_layer_edges_along_chord`'s chord-crossing
  quadratic, `matter.density_matter_func_exp`'s profile, and
  `oscprob._osc_prob_ip_exp_core`'s two interaction-picture equations are
  now `.. math::` blocks whose LaTeX matches `docs/source/methodology.rst`
  for the identical mathematics.
- Physical units carrying exponents are now math-mode throughout the
  docstrings, so they render as eV², g cm⁻³ rather than as the literal
  text `eV^2`, `g cm^{-3}`. The convention adopted, and applied uniformly
  across `oscprob`, `globaldefs`, `matter`, `earth`, and `cli`, is to
  math-wrap a whole unit group when — and only when — it carries an
  exponent (``[:math:`\text{g cm}^{-3}`]``), matching the pre-existing
  `\text{cm}^{-3}` usage in `oscprob`'s module docstring; exponent-free
  units (`[eV]`, `[radian]`) stay plain text, since they already render
  correctly. `docs/source/functions.rst`'s lone `:sup:`-based unit was
  converted to match. This also normalizes the `eV^-1`/`eV^{-1}` spelling
  inconsistency.
- All 70 remaining `>>>` doctest-style prompts (in `oscprob`'s module
  docstring and in `globaldefs.load_nufit_params`) are now
  `.. jupyter-execute::` blocks, the convention already used by the ~60
  `osc_prob_*` wrapper docstrings. No doctest runner has ever existed in
  this project — no `--doctest-modules`, no doctest CI step, no
  `doctest_namespace` fixture — so these examples were never executed or
  checked by anything, and had silently drifted out of sync with the code
  (see "Fixed" below). As executed cells they now show genuine, always-current
  output and will fail the docs build if they ever break again.
- **Breaking:** invalid input now raises `ValueError` instead of printing a
  message and calling `sys.exit(1)`. There were 62 such aborts across
  `oscprob` and `earth`; a library that terminates the interpreter cannot be
  recovered from in a notebook, a scan loop, or a caller that wants to fall
  back, and it made the failures impossible to assert on in tests. Most were
  already `raise ValueError(...)` inside a `try` block whose `except`
  immediately swallowed it and exited, so the messages are unchanged — they
  now propagate instead of being printed. `validate_input_battery` follows
  suit: it raises rather than returning `1`, so its return type is now `None`
  (it previously returned `0`/`1`, and every caller compared against `1`).
  Error messages raised as exceptions use the plain-text
  `gd.ERROR_MSG_NO_COLOR` prefix, since ANSI codes are meant for a terminal
  and end up in tracebacks and logs. The `magnus` CLI catches these and
  reports them as ordinary argument errors, so its behavior is unchanged.
- A single source of truth for the version number: the `version` field of
  `pyproject.toml`. `magnus/version.py` now resolves it via
  `importlib.metadata` when installed, falling back to parsing
  `pyproject.toml` when running off `src/` on the path, and
  `docs/source/conf.py` imports it rather than repeating it. The eleven
  decorative per-module `__version__` strings (which disagreed with each
  other — `"2.0"` in `magnus.py`, `"0.10.0"` in `oscprob.py`, `"1.0"`
  elsewhere — and which nothing ever read) are gone.
- The `strategy` keyword is now reachable from the command line as
  `magnus prob --strategy {auto,hybrid,magnus}`, for the environments where
  the Hamiltonian actually depends on position (`sun`, `earth`, and `matter`
  with `--density-profile exp`). It was previously Python-API-only, so CLI
  users silently got `'auto'` with no way to opt out or to force it.
- `globaldefs.set_color_output(enabled)` turns the ANSI color in the
  warning/error/tolerance prefixes on or off. `WARNING_MSG_NO_COLOR` and its
  siblings had existed and been exported since the beginning but were never
  used by anything: every call site hardcoded the colored variant, so there
  was no supported way to get clean output into a log file or a rendered
  notebook.
- `ruff check` is now blocking in CI rather than `continue-on-error`, with
  the rule configuration in `[tool.ruff.lint]` in `pyproject.toml`. Two
  codebase-wide conventions are exempted explicitly (`E741`, since `l` is the
  standard symbol for position here, and `E701` for the one-line cleanup
  guards); everything else was fixed, so the tree is clean and a new finding
  fails the build instead of being reported into a green checkmark.

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
- ~47 docstring/type-annotation mismatches across the codebase, found by
  a systematic audit comparing every function's actual return/parameter
  type against what its docstring documented (not just the
  `numpy.ndarray`-documented-as-`list` case originally flagged in
  `hamiltonians{2,3,4,5}nu.py` and `oscprobstd.py`, but every such
  disagreement anywhere in `src/magnus`). Also fixes two unrelated
  copy-paste errors caught along the way in `oscprobstd.J()`'s docstring
  (a mislabeled parameter description and a wrong worked example).
- Three latent errors in `oscprob`'s module-docstring examples, exposed by
  converting them to executed cells: `np.array([1.0, 10.0 100.0])` was
  missing a comma (a `SyntaxError`), `osc_prob_3nu_matter_constant_density`
  was called without its `oscprob.` prefix (a `NameError`), and the flavor
  indices were listed as "`NUE`, `NUMU`, and `NUMU`" instead of `NUTAU`.

- `adiabatic.hybrid_propagator` could report `certified=True` without having
  certified anything. Its three refinement knobs saturate at different
  iterations (`n_probe` at 5, `n_points` at 6, `threshold` at 11), so by
  iteration 12 all three were pinned: that iteration recomputed bit-identical
  inputs and the agreement test compared a result with itself. The loop now
  stops as soon as every knob has saturated and reports `certified=False`,
  which is the honest answer. Covered by a regression test that fails against
  the old code.
- `adiabatic.find_resonance_candidates` (and everything built on it) no longer
  evaluates the user's `H_func` outside the requested `[l0, l1]`. The
  finite-difference stencil reached to `l0 - h` and `l1 + h` at the endpoints,
  which can raise or return nonsense for a Hamiltonian defined only on its
  physical domain -- `earth.density_matter_func_prem`, for one, raises beyond
  `EARTH_RADIUS`. The stencil is now one-sided at the boundaries.
- The speedup chart in `docs/source/adiabatic_strategy.rst` disagreed with the
  validation table directly above it: it showed a 25,800x bar for the case the
  table reports as ~30x (25,800x was a different, unlisted measurement). Both
  now come from `VALIDATION_GRID` in the new `docs/make_figures.py`, so they
  cannot drift apart again.
- `docs/source/installation.rst` claimed `src/magnus/` had to be on the Python
  path "for a few modules that resolve sibling imports directly". No such
  import has existed since the package was flattened; only `src/` is needed.
  The README and `quickstart.rst` also still carried
  `sys.path.extend(['src', 'src/magnus'])  # until pip packaging lands`, long
  after packaging landed.
- `HybridCertificationWarning` is exported from `magnus.oscprob.__all__`. It
  was the only public class or function in the package missing from its
  module's `__all__`, despite being raised, documented, and cross-referenced.
- Four `Returns` sections said `np.narray` instead of `np.ndarray`, and
  `earth.coordinates_of_named_location`'s message said "the given name of the
  the location". Several other error messages had typos (`wil`, `lengh`, "only
  of the two", and a mangled function name in
  `values_to_unspecified_osc_params`).
- Every `:func:`/`:class:` cross-reference in the docs now resolves. Bare names
  in `quickstart.rst` and `adiabatic_strategy.rst` are qualified;
  sibling-module references in `hamiltonians{2,3}nu.py` are qualified; and
  references to private helpers (for which autoapi never emits targets, so they
  rendered as plain text) are now inline literals, which is what they should
  have been.
- Undocumented parameters: `tol` in `earth.density_matter_func_prem` and
  `earth.earth_radial_distance_from_depth`, and `A_eval_mode` in
  `magnus.magnus_expansion` and `magnus.magnus_expansion_multislab`.
- `osc_prob` and the four `osc_prob_{2,3,4,5}nu_vacuum` wrappers documented
  their parameters without types, unlike the other 167 functions in the
  package, so the rendered docs dropped the type column for the single most
  important function in the API. All 99 entries now carry the type from the
  signature.
- Stale "Routine listings": `globaldefs` was missing `load_nufit_params` and
  `matter` was missing `exp_density_profile`. `cli.build_parser` and `cli.main`
  were the only public functions with no `.. versionadded::` tag.
- `.github/workflows/tests.yml` only ran on pushes to `main` and `dev`, so
  every push to `dev-plotting` -- the branch this work happens on -- got no CI
  at all. It now covers `dev-*` too, and installs the package with
  `pip install -e .` rather than only its requirements, so the console script
  and the version lookup are exercised the way a user's install is.
- `docs/source/cli.rst` quoted `magnus prob --help` from a hand-copied paste
  that had gone stale (it never mentioned `--version`). It is now generated by
  `docs/regen_cli_help.py`, with a CI job that fails if the page and the parser
  disagree.

### Removed

- `docs/source/sandbox/`, an untracked, unused pydata-theme experiment
  directory.
- 1,257 lines of ad-hoc scratch code in the `if __name__ == "__main__":`
  blocks of `oscprob` (1,218 lines, ~8% of the file), `earth`, `matter`, and
  `magnus`. The first three could not run at all: executing those modules as
  scripts fails at import under the flattened package layout. `magnus`'s block
  did run, but it printed an unasserted comparison of the three integration
  methods that `tests/test_magnus_expansion.py` already parametrizes over with
  real assertions. The only `__main__` blocks left are the two genuine console
  entry points, `cli.py` and `__main__.py`.
- All `.. versionchanged::` directives. This being the first public release,
  there is no earlier published version for behavior to have changed *from*,
  so the six of them described how the code works rather than what changed.
  Their content is kept, reworded as `.. note::` blocks on the same functions.
- The eleven decorative per-module `__version__` strings, and the separate
  `[Unreleased]`/`[0.10.0]` changelog sections, which are consolidated here:
  nothing was published before this release, so the split served no reader.
