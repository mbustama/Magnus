# Changelog

All notable changes to Magνs are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **A compiled Cayley–Hamilton backend for the matrix exponential, selected by
  `magnus.magnus.EXPM_BACKEND`.**  `np.linalg.eigh` costs ~1.27 µs per 3×3
  *whatever the stack size* (measured 1.268 µs at N=108, 1.279 µs at N=4096 —
  flat, because it loops over LAPACK internally instead of vectorising).  The
  new `magnus.expmkernels` applies to `K` the polynomial interpolating
  `exp(-iλ)` on its spectrum instead: no eigenvectors, and the eigenvalues in
  closed form.  **6.8× on the exponential** at N=108 (162.6 → 23.8 µs), 7.3× at
  d=2, and **2.11× end to end** on a 60-energy PREM scan (9291 → 4409 µs, i.e.
  73.5 µs per energy).

  The gap between 6.8× and 2.11× is Amdahl's law: the exponential is about a
  third of a slab pass.  Quoting the 6.8× as a package speed-up would be quoting
  the wrong number.

  `'auto'` (the default) uses the kernel for 2×2 and 3×3 when numba is
  installed and `eigh` otherwise, and cannot fail; `'numba'` makes a missing
  numba an error rather than a silent downgrade; `'eigh'` is the reference route.
  **Dimensions 4 and 5 keep `eigh`** — there is no practical closed form for a
  4×4 or 5×5 Hermitian eigenproblem, so 4ν and 5ν stay correct and are not
  accelerated.  numba is an optional dependency (`pip install 'magnuspy[fast]'`),
  costing ~90 ms of `import magnus` when present plus a one-off ~0.7 s compile
  per kernel, cached to disk thereafter.

  Switching backend moves probabilities by at most 4.6e-15 across PREM chords,
  energy scans, NSI resonances, constant density and vacuum.

  Degeneracy is the whole risk in such a scheme, and two facts remove it.  A
  Hermitian matrix is never defective, so matching `exp` on the *distinct*
  eigenvalues is already exact and the confluent (Hermite) form is not needed.
  And with eigenvalues sorted and the spectrum shifted to put the median at
  zero, the one ill-conditioned coefficient multiplies a matrix whose norm
  shrinks with the same gap, so its contribution is bounded by `ε·gap` and
  *vanishes* as the gap closes.  There is therefore no tolerance, no crossover,
  and no near-degenerate branch to place: the error is 1e-16 at splittings of
  1e-2, 1e-6, 1e-10, 1e-14 and exactly zero alike.

  Two things this cost, both now pinned by tests.  The closed-form eigenvalues
  **degrade to ~4e-9 at an exact degeneracy** (`arccos` has infinite derivative
  where a repeated root sits) and the exponential stays at 2.5e-16 anyway,
  because interpolation error is *second* order in the displacement of a
  coalescing node; both halves are asserted, the sloppy one included.  And the
  kernel must read the **lower** triangle, because that is the one `eigh` reads
  (`UPLO` defaults to `'L'`) and `_expm_stack` admits input anti-Hermitian only
  to 1e-12 — a kernel reading the upper triangle exponentiates a different
  matrix on such input and the two backends diverge by ~2e-12, large enough to
  matter and small enough to look like rounding.

### Fixed

- **`_expm_stack`'s docstring claimed the `eigh` route was exactly unitary, and
  it is not.**  `U†U - I` measures 4e-16 for a single 3×3 and 4e-15 for a stack
  of 4096 — growing with stack size, never zero, because reconstruction from
  eigenvectors rounds like any other floating-point product.  The claim that
  probabilities "sum to 1 by construction" was the part worth correcting: they
  sum to 1 to about 1e-15, which is worth relying on, by rounding rather than by
  construction.  No behaviour changed by this entry.

- **The refinement ladder could stop while the answer was still outside the
  requested tolerance, and report success.**  ``t_breakpoints`` (the ~14 PREM
  layer crossings) are re-inserted into every level's grid, so at small counts
  the nominal refinement and the real one are different things: a nominal
  2 -> 3 slab step is a **16 -> 17 edge** step, a 6% refinement rather than a
  50% one.  Two grids differing by 6% agree for reasons unrelated to having
  converged, and ``np.allclose`` read that as success.

  An agreement now only counts when the two levels compared were genuinely
  different grids -- ``len(t_slab_edges)`` must grow by at least
  :data:`oscprob.MIN_EFFECTIVE_REFINEMENT` (1.25, measured).  Over 120 Earth
  configurations (costhz -0.15 to -0.99, 0.5-8 GeV, three tolerances, scored
  against a reference verified converged to 1e-13): one silent violation, 2.1x
  outside the tolerance asked for, and five more that missed but warned -- all
  six now zero.

  **This changes Earth results, by making them more accurate**: one of the
  twelve ``bitident`` rows moves, the Earth single point, by 2.3e-04.  It is
  inert wherever there are no breakpoints (solar is bit-identical), and does
  not reach the separable engine that answers Earth *scans*, which was checked
  separately and does not have the defect.  Cost is about 11% of wall clock
  (a 60-energy Earth scan goes from 9 ms to 10 ms) despite the median slab
  count rising from 9 to 21, because the added levels are the cheap small ones.

### Changed

- **``rtol``/``atol`` are documented for what they are: a stopping criterion,
  not an accuracy guarantee.**  The ladder halts when two successive levels
  agree; it never estimates the error of the answer it returns, which is a
  weaker promise than a stepping ODE integrator's ``rtol`` makes.  Corrected
  in ``osc_prob``, ``adiabatic.hybrid_propagator``, the CLI's ``--rtol``/
  ``--atol`` help, ``README.md``, ``architecture.rst`` (which said "until
  rtol/atol is met"), and a new section of ``implementation_details.rst`` that
  the others link to.  No behaviour changed by this entry.

- **``convergence_info`` reports what the ladder did.**  Alongside the existing
  ``n_slabs``/``n_tpts_per_slab`` it now carries ``n_slab_edges`` and
  ``n_slab_edges_previous`` (which make the real refinement step visible),
  ``n_slabs_previous``, ``n_tpts_per_slab_previous``, ``last_gap`` (None when
  only one level was ever computed), ``n_agreements``, and
  ``tolerance_achieved`` -- the programmatic form of
  ``ToleranceNotAchievedWarning``.

  It deliberately carries **no error estimate**.  Converting the gap into one by
  Richardson extrapolation, as the sibling NuOscProbExact does, was measured and
  rejected: Magnus has no stable convergence order (fitted on Earth chords it
  scatters from 1.4 to 7.2 against nominal orders of 2 and 4), and because
  breakpoints make the effective refinement ratio as low as 1.06 rather than
  1.5, dividing by ``r^p - 1`` under-reports the true error by 6-20x *even
  where the power law holds exactly*.  Under-reporting is the dangerous
  direction.

### Added

- **Palindromic density profiles are exploited on Earth chords.**  A chord through a
  spherically symmetric Earth meets every radius twice, so its density profile reads
  the same from either end.  The Magnus core now evaluates the Hamiltonian on the first
  half of such a slab chain and derives the rest by reversal, halving the calls to the
  caller's `H_func`.  Worth **1.4x-1.67x** on a single point and **1.56x-1.64x** on an
  energy scan when that Hamiltonian is expensive; plain PREM, whose density lookup is
  cheap, pays about 10% for it.  New public `magnus.magnus.USE_PALINDROME` (module
  switch, `True`) turns it off, and `magnus.magnus.palindromic()` is the predicate.

  The saving is halved Hamiltonian evaluations and nothing else, so it is worth what
  that Hamiltonian costs.  Standard PREM *scans* are unaffected: they are answered by
  the separable engine, which already evaluates the profile once and shares it across
  energies -- the same saving, taken earlier.

  Symmetry is **declared** by the Earth entry points, where it is a fact of chord
  geometry, not detected: detecting it would need the very evaluations the optimisation
  skips.  There is deliberately no user-facing way to declare it of an arbitrary
  profile.

  **This moves Earth single-point results by up to 8.6e-15 relative.**  The mirrored
  slab's nodes are reached by a different floating-point expression for the same real
  number, so the change is inherent rather than incidental.  `USE_PALINDROME = False`
  reproduces the previous numbers exactly.

### Changed

- **`oscprob.BATCH_WORKING_ENTRIES` lowered from 4,194,304 to 65,536** (about 67 MB to
  about 1 MB).  The batched scan engines are memory-bound, and the previous value was
  large enough that their working set spilled cache.  Measured across fifteen workloads
  on three engines, the new value is **1.19x-1.38x** quicker on Earth energy scans,
  1.06x-1.16x on cumulative baseline scans, flat on short scans and on the
  interaction-picture engine, and never slower anywhere.  **Bit-identical** at every
  budget tested -- tiles are independent and only concatenated -- so this changes no
  result.  Peak memory of a long scan drops accordingly.

### Fixed

- **The position-profile cache no longer hands out writable arrays.**  Values are
  returned by reference to every later caller asking for the same position grid, so a
  write through any one of them would have silently changed what the others received --
  and the cache sits under the matter term of the Hamiltonian, so the symptom would
  have been a wrong probability with nothing raised.  Cached arrays are now marked
  read-only, turning that into an exception at the point of the write.


## [1.0.0rc1] - 2026-07-31

First public release candidate.  Magνs was developed privately up to this
point, so everything below is new to anyone outside the project.  The entries
are still grouped as Added/Changed/Fixed/Removed, describing the development
history that produced this release -- "Changed" and "Fixed" are relative to
earlier private states of the code, not to any published version -- because
that history is the most useful record of *why* the code looks the way it does.

### Added

- **Notebook 09 is now a Lorentz-invariance-violation notebook.** It had been a
  scratch pad: no markdown, no figures, and no LIV content at all -- just
  `print(module.__all__)` and a `print(sys.path)` that had stopped working. It
  is now built around the one thing that makes LIV findable, the energy
  scaling: the vacuum term falls as `1/E` and the matter term is flat, so both
  switch *off* at high energy, while the LIV term grows as `E^n` and switches
  *on*. Four figures follow from that -- where standard oscillations stop and
  LIV keeps going, how the operator dimension `n_liv` sets that crossover, that
  matter does not rescue the standard prediction (`V_CC` does not grow with
  energy either), and how a null result at high energy becomes a limit on the
  LIV eigenvalue, with the `b*L*E^n = 1` estimate landing where the curve
  visibly departs.
- A `magnus.plotting` module of pre-packaged figures, so that a plot in the
  notebooks costs one call rather than the twenty-five to forty lines of
  `gridspec_kw`, tick locators, legend keywords and `savefig` that each figure
  used to carry. Taking stock of the fifty-odd notebook figures first showed
  that most are the same figure with different data: curves against baseline,
  against energy, against a sterile mixing angle, and the matrix-exponential
  convergence studies all reduce to one shape -- curves against a swept
  variable over an optional relative-error subpanel -- which is `plot_curves`,
  with `plot_probability_vs_baseline` and `plot_probability_vs_energy` as
  presets over it. Only three layouts are genuinely distinct and get their own
  functions: `plot_probability_with_profile`, `plot_biprobability` and
  `plot_oscillogram`, plus `plot_probability_with_average` for the decohered
  overlay. Defaults reproduce the existing house style exactly, and every
  function returns `(fig, ax)` so a packaged figure is a starting point rather
  than a dead end. `prob_label` absorbs the helper that had been copied into
  several notebooks, extended to cover the sterile states those notebooks
  needed but it did not.
- **`plotting.plot_curves_stacked`, for small multiples.** Auditing which
  notebook figures still built their own axes turned up one shape the module
  had missed: the same plot repeated once per case down a shared abscissa --
  notebook 07's four panels, one per detector. It is worth packaging for the
  same reason as the rest: the reader compares *between* panels, so every panel
  must carry identical limits, scales and tick spacings, while the abscissa
  labels, title, legend and ordinate label belong to exactly one panel each.
  Hand-built that is four formatting loops plus a frameless full-figure subplot
  added purely to hang a shared label on, and it is where a stack quietly stops
  being readable once one panel drifts. `legend_proxies` also retires a trick
  the notebook used: plotting dummy points outside the axis limits to
  manufacture legend handles for entries that describe a line *style* rather
  than any one curve. Notebook 07's two remaining hand-built figures now use
  this and `plot_curves`; the second turned out to need no new function at all,
  only a `gridspec_kw` that had never done anything.
- **`osc_prob_energy_baseline(..., cumulative=True)`: a whole baseline scan from one
  traversal of the profile.** The evolution operator is a time-ordered product, so
  `U(0->L2) = U(L1->L2) U(0->L1)`: every requested baseline is a *prefix* of the next,
  and recording the running product yields the whole scan at once instead of re-walking
  the profile N times. It is the `reduce` already in `osc_prob` with its intermediates
  kept rather than discarded.

  The grid is the union of the requested baselines (so each answer lands on a slab edge
  and is read off, never interpolated), a uniform accuracy grid, and any
  `t_breakpoints`. Sizing that accuracy grid is the one way a cumulative scan goes
  *silently* wrong — the traversal has nothing to compare itself against — so it is not
  guessed: one ordinary adaptive `osc_prob` call at the longest baseline reports the
  slab count it needed, which is the definition of the accuracy grid, and brings the
  existing safeguards and warnings with it. On a solar profile a plausible-looking guess
  of 2000 slabs is wrong by 1.6e-2 where the inherited number is right.

  Chunked traversal and conversion to probabilities at each snapshot are requirements
  rather than optimisations: they keep peak memory at `O(block) + O(result)` instead of
  holding N complex unitaries beside the answer.

  Measured against `solve_ivp`: a 1000-point solar scan takes 12.0 s per-point for an
  error of 5.6e-5, and 0.10 s cumulative for 5.1e-6 — **124x faster and 11x more
  accurate**. Opt-in rather than automatic: the two paths use different grids, so results
  differ within the requested tolerance.
- **Matplotlib is an optional dependency**, declared as the `plot` extra
  (`pip install 'magnuspy[plot]'`). The engine still needs only NumPy, SciPy
  and joblib: someone computing probabilities inside their own analysis code
  should not have to install a plotting stack. `magnus.plotting` imports
  cleanly without Matplotlib -- it defers the import into the calls that draw
  -- so `import magnus` works on a core-only install and only a plotting call
  raises, as `MatplotlibNotFoundError`, naming the command to fix it.

- The PyPI distribution is named **magnuspy**; the import package remains
  `magnus`. Plain `magnus` was already taken on PyPI by an unrelated project, so
  `pip install magnus` would have fetched someone else's package and the release
  workflow would have failed on upload. The two names are independent in Python
  packaging, so `pip install magnuspy` then `import magnus` is all that changes,
  and the console script is still `magnus`. `version.py` looks the distribution
  up under the new name: querying the import name would have raised
  PackageNotFoundError and fallen through to parsing `pyproject.toml`, which is
  absent from an installed wheel, so every installed user would have reported
  `0.0.0+unknown`.
- A license. Magνs is released under the GNU General Public License v3.0 only;
  `LICENSE` carries the full text, and it is declared in `pyproject.toml` as the
  PEP 639 SPDX expression `license = "GPL-3.0-only"`, which is what appears in
  the built distribution's metadata as `License-Expression`. It previously read
  `TBD`, so a `pip install` would have shipped the license text while declaring
  no license at all. The build requirement moves to `setuptools>=77`, the first
  version that understands the SPDX expression form; note that PEP 639 forbids
  pairing it with a `License ::` classifier, so there is deliberately none.
  Referenced from the README, the docs landing page, and both file trees.
- Status badges on the README and the docs landing page: CI tests, code quality,
  documentation, license, the supported Python version (3.10+, matching
  `requires-python` and the CI matrix rather than a copied-in default), the code
  style (ruff, which `lint.yml` actually enforces), and PyPI downloads via
  pepy.tech. The downloads badge reports on the distribution name, `magnuspy`,
  and stays blank until the first release is published -- pepy serves a 404 for a
  project PyPI does not know, so it renders as a broken image rather than as a
  count of zero.
- `magnus.expansionterms`: derives the terms of the Magnus expansion from the
  Bernoulli-number recursion symbolically, at any order, in exact rational
  arithmetic (`bernoulli`, `bernoulli_factor`, `omega_terms`, `magnus_terms`,
  `count_terms`, `format_term`, `print_magnus_terms`). The numerical core's
  coefficients were typed in and nothing checked them against the recursion they
  come from; the test suite now regenerates them and compares, order by order, at
  machine precision. See `docs/source/expansion_terms.rst`.
- Magnus orders 7 through 10. `MAGNUS_EXP_ORDER_MAX` rises from 6 to 10, and is
  now defined once in `magnus.magnus` and re-exported by `globaldefs`, which used
  to carry its own copy of the number. The default expansion order is unchanged.
  Orders 1-6 keep their hand-written expressions, which are hot and worth reading;
  7-10 are generated from the closed form of the recursion (every term is a
  right-nested chain of lower-order `Omega_m` around `A`, indexed by the
  compositions of `n-1`), with shared suffixes memoized so each distinct nested
  commutator is built once. Verified two ways: exact agreement with the symbolic
  generator, and measured convergence rates that keep improving with order (order
  8 reaches ~h^10 against an ODE ground truth where order 6 reaches ~h^8).
- `magnus.magnus.MagnusHighOrderCostWarning`, raised when an order above 6 is
  requested with a quadrature method. The number of terms roughly doubles per
  order, and the measured cost per slab is 2.7x order 6 at order 7, rising to
  ~17x at order 10. Higher order does converge faster in the slab width, so it is
  a trade rather than a mistake -- but narrowing the slabs at order 4 or 6 often
  reaches a given accuracy for less total work.
- Command-line calculator (`magnus prob`, also runnable as `python -m magnus`)
  for computing a single oscillation probability from the shell, covering
  vacuum, matter (constant/exponential density), Earth, and Sun, with
  standard, NSI, and LIV scenarios, for 2-5 flavors. See `docs/source/cli.rst`.
- Full `pytest` suite (`tests/`) and GitHub Actions CI: `tests.yml` (matrix
  across Python 3.10-3.12), `lint.yml` (Ruff), `pages.yml` (Sphinx ->
  GitHub Pages), `publish.yml` (PyPI on release).
- Test-coverage measurement. `pytest-cov` joins the `test` extra, the settings
  live in `[tool.coverage.run]` in `pyproject.toml`, and a separate Coverage job
  in `tests.yml` reports the figure on each run's summary page and uploads
  `coverage.xml` as an artifact. Branch coverage is enabled: line coverage alone
  overstates how well this package is tested, because `oscprob.py` is largely
  thin wrappers that one parametrized test sweeps in a single pass, whereas what
  matters is whether the dispatch chain, the refinement caps and the warning
  paths are each taken in both directions. The build fails below 90%, a floor set
  two points under the measured 92% so that it catches a regression without
  tripping on the fraction of a percent that moves between interpreters; it was
  left unset until the figure had settled, since a threshold invented before the
  first measurement either sits below the real figure and never fires, or above
  it and blocks unrelated work. The job also carries a Codecov upload step
  that stays dormant until a `CODECOV_TOKEN` secret exists, so nothing leaves the
  repository until public coverage reporting is deliberately switched on.
- Two structural sweeps closing what that first coverage run found: 23 of the 36
  NSI/LIV `osc_prob_*` wrappers, and 25 of the `hamiltonian_*` builders the
  package exports, were executed by nothing at all -- not by the tests, and not
  by the library either, which reaches for the `*_energy_independent` variants
  instead. Several of the wrappers even appeared in a parametrize list, but in
  tests that inspect a signature or the source text without ever calling the
  function, so a mistyped keyword in any of them would have shipped unnoticed.
  `test_every_bsm_wrapper_runs_and_is_unitary` now calls every one of them with
  non-zero NSI/LIV parameters and checks the result is a valid probability
  matrix, and `test_every_exported_hamiltonian_builder_is_hermitian` builds every
  exported Hamiltonian once, with complex off-diagonal couplings, and checks
  Hermiticity. Both discover their own subjects -- from the module and from
  `__all__` respectively -- so a name added later is swept without anyone
  remembering to extend a list. Both were verified to fail against a deliberately
  broken library before being kept.
- Tests for the paths that report failure rather than a result: the refusals of
  the closed-form interaction-picture integrator, both ways
  `adiabatic.hybrid_propagator` can decline to certify a non-converged patch,
  and what the dispatch layer then does with an uncertified point --
  `strategy='auto'` abandoning the batch for the general Magnus path,
  `strategy='hybrid'` keeping the answer and raising
  `HybridCertificationWarning`. None of these had ever been executed, which is a
  poor place for a blind spot in a package whose central claim is that it knows
  when it cannot verify its own answer. Also `magnus.version`'s two resolution
  routes, which must agree; the CLI's argument-error paths; the unknown-location
  error in `earth.coordinates_of_named_location`; the `verbose=2` banner and
  run-parameter dump; and `globaldefs.set_color_output`.
- `magnus.avgprob`, and the `average` keyword that reaches it from every
  oscillation-probability function through the shared `**kwargs` chain: the
  phase-averaged (fully decohered) probability, which is the exact
  `L/E -> infinity` limit astrophysical neutrinos arrive in. For a
  position-independent Hamiltonian it is closed-form, one eigendecomposition
  rather than the resolution of some 10^15 radians of phase, and for vacuum it
  depends on neither energy nor baseline, so one matrix serves an entire flux
  calculation. A position-dependent profile decoheres in the eigenbasis at
  production and is carried along the levels of the instantaneous Hamiltonian,
  with level-crossing probabilities taken from the convergence-checked Magnus
  patches of `magnus.adiabatic` rather than from a Landau-Zener formula -- they
  reproduce that formula to a few parts in a thousand, which is the check rather
  than the method. A profile with discontinuities (PREM) has no closed form and
  is averaged over an energy window instead, which is a different quantity and
  says so.
  Whether the averaged limit applies at all is decided per pair of eigenvalues
  from the phase they accumulate, not assumed: pairs that have neither decohered
  nor stayed coherent raise `magnus.oscprob.PhaseAveragingWarning`, and
  near-degenerate ones are summed coherently within blocks, since the naive
  incoherent sum returns a spurious mixture where the correct answer is that
  nothing oscillates. Documented in `docs/source/averaged_probability.rst`, with
  worked examples in `notebooks/10_magnus_averaged_probability.ipynb`.
- `magnus.oscprob.IP_EXP_N_SLABS_CAP` and `magnus.oscprob.IP_EXP_LOOP_CAP`: the
  interaction-picture integrator's slab and loop ceilings, previously written as
  bare numbers inside the function. Naming them changes no behaviour, and makes
  the method's conduct at the ceiling testable at a small cap -- reaching two
  million slabs to observe what happens at the boundary costs gigabytes and
  minutes, so with the values inlined those paths could not be tested at all.
- `docs/check_doc_snippets.py`, which executes every `.. jupyter-execute::`
  block in the documentation -- in the RST pages and in the docstrings autoapi
  renders -- and reports the page, line and traceback of any that fails. The
  documentation's examples are already run at build time, so a broken one fails
  CI; what this adds is the ability to find out in about a second, rather than
  from a full Sphinx build. The fast build used while writing docs stubs those
  directives out, so it validates the prose and the cross-references while
  saying nothing at all about the code, and a page can build cleanly while being
  broken.
- `tests/test_validation.py`, covering the guards that reject bad input: the
  per-flavour parameter dictionaries, flavour indices and expansion orders out
  of range, mismatched energy/baseline arrays, negative densities and
  composition ratios, and the Earth entry point's location, zenith-angle and
  baseline combinations. Each asserts that the input is refused as a
  `ValueError` naming the parameter at fault, rather than as whatever the
  interpreter happened to raise further downstream; writing them is what
  uncovered the `except KeyError` bug listed under Fixed.
- Three branches that cannot be reached are now marked `# pragma: no cover`,
  each with the argument for why written beside it: the `is_repeat` return in
  the interaction-picture integrator (below the ceiling the slab count strictly
  increases, at the ceiling every branch returns, so a repeat can never be
  compared), its no-progress guard (`round(2n) == n` has no solution for
  `n >= 1`, and it exists only to make a smaller growth factor safe), and the
  CLI's solar baseline check (the general one rejects that input first). The
  Windows-only ANSI setup line in `globaldefs` is marked too. They are kept
  rather than deleted: each becomes live again if a nearby constant changes.
- A sweep asserting that both construction routes of every builder offering the
  choice agree: the hardcoded expression for each matrix entry, and the explicit
  product of the mixing matrix, the mass matrix and the conjugate transpose that
  it was derived from. Nothing had been keeping the two from drifting apart for
  the Hamiltonian builders -- the existing checks covered only the 4x4 and 5x5
  mixing matrices -- and the product route was unexecuted for both vacuum
  energy-independent Hamiltonians, which `oscprob` calls on every run.
  Antineutrinos are included, since conjugating the mixing matrix is a separate
  line in several builders.
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
- `notebooks/12_magnus_adiabatic_hybrid_strategy.ipynb`: live comparison of
  `strategy='auto'`/`'hybrid'`/`'magnus'` for 2- through 5-flavor
  Hamiltonians (standard oscillations and an engineered BSM/NSI
  resonance), each cross-checked against a tight-tolerance `solve_ivp`
  ground truth in both runtime and accuracy, plus a real-data plot of the
  instantaneous eigenvalues and detected non-adiabatic window for the 3ν
  BSM case. Reproduces, live, the validation described in
  `docs/source/adiabatic_strategy.rst`.

### Changed

- `integration_method='gl'` now raises for orders above 6 rather than silently
  computing an order-6 result. The Gauss-Legendre commutator-free schemes are
  separately derived integrators, not products of the Magnus recursion, so they
  do not extend along with it. The check sits in `_gl_nodes` as well as in the
  input validator, since the validator is skipped when `validate_input=False` and
  that flag would otherwise reopen the silent-degradation path.
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
- The documentation build now runs on pull requests, not only on pushes to
  `main`. It was previously built only by `pages.yml`, so a page that failed to
  build, or silently dropped out of the sidebar, was not caught until after it
  had been merged. The pull-request build is deliberately stricter than the
  deploy build -- `-W` turns Sphinx warnings into errors, which is what catches
  a broken cross-reference or a page missing from the toctree -- while the
  deploy build stays permissive, so a late warning can never be what stops a
  release from publishing.
- `ruff check` is now blocking in CI rather than `continue-on-error`, with
  the rule configuration in `[tool.ruff.lint]` in `pyproject.toml`, including an
  explicit `select`. Leaving the selection to ruff's default is not safe for a
  blocking check: the default widened in ruff 0.16 and turned CI red on a
  codebase that had not changed. Two
  codebase-wide conventions are exempted explicitly (`E741`, since `l` is the
  standard symbol for position here, and `E701` for the one-line cleanup
  guards); everything else was fixed, so the tree is clean and a new finding
  fails the build instead of being reported into a green checkmark.

### Fixed

- **A batched solar call could exhaust the machine's memory.**
  `_osc_prob_ip_exp_core` — the closed-form interaction-picture integrator behind
  `osc_prob_2nu_sun` and its NSI/LIV siblings — built temporaries of shape
  `(n_energies, n_slabs, d, d)` while its ladder doubled `n_slabs` toward
  `IP_EXP_N_SLABS_CAP = 2_000_000`. The reasoning recorded beside that ceiling was
  about time only ("each slab costs one 2x2 eigendecomposition"), and never
  accounted for the energy count multiplying the working set. Measured, on the
  documented solar use case: **~1.3 GB per energy** — 1.56 GB at one energy, 5.34 GB
  at four, and a `MemoryError` of shape `(8, 2000000, 2, 2)` at eight. Nothing about
  the call is pathological; it is the advertised batched form of a public wrapper, and
  notebook 03 scans 1000 energies over exactly that range. The notebooks never hit it
  only because they use raw `osc_prob` loops. The working set is now tiled over both
  axes against a fixed budget (`BATCH_WORKING_ENTRIES`), giving a peak that is **flat
  in the energy count** — 79 MiB at 4 energies and 79 MiB at 256. The tiling is exact:
  the slab product is folded in the same order with the same parenthesis nesting, and
  a test pins the output at *bit equality* against an untiled run. Found while
  measuring the energy axis, and pre-existing — reproduced unchanged at `155e01e`.
- **The same integrator burned twenty-one refinement levels to reach refusals it could
  have predicted.** Certification requires `max|Omega_t|` below a trust threshold, and
  that maximum is bounded below by the diagonal entries, which have a closed form. When
  even that bound exceeds the threshold at the slab ceiling, no reachable slab count can
  certify. Two evaluations of the potential and no allocation now detect it, and the
  method refuses immediately instead of doubling its way to the same answer. It is a
  bound rather than an estimate, so it can only report "impossible" and cannot abandon a
  case that would have converged. At 10 MeV over a solar radius it fires for any
  tolerance of 1e-4 or tighter, or any baseline beyond two solar radii.
- **A scan whose result could not fit reported it as an out-of-memory kill.** Tiling
  bounds the engines' working set, but nothing shrinks the answer: N points over d
  flavors is `N*d*d` floats either way. `osc_prob_energy_baseline` now checks that
  against the operating system's free-memory figure and raises a `MemoryError` naming
  the size, rather than letting an overcommitting kernel take the machine down instead
  of the process. The check costs one multiply below a 64 MiB floor, and never blocks
  where free memory cannot be read.
- **A requested `n_slabs` was silently discarded whenever a tolerance was on,
  and that certified wrong answers.** `osc_prob` documented and implemented
  "if `rtol` or `atol` is given, `n_slabs` is ignored": the adaptive ladder
  started at `min_n_slabs = 1` regardless of what the caller asked for. The
  seed that replaced it, `magnus.suggest_n_slabs`, measures the *integral* of
  the Hamiltonian along the path, which is blind to structure that averages
  out. On notebook 03's castle-wall profile -- 50 square density walls -- the
  whole trajectory accumulates only ~9 radians, so a call asking for
  `n_slabs=150` was seeded with 2 slabs and stopped at 4. Four slabs cannot see
  fifty walls, and the ladder they sit on does not converge, it thrashes:
  0.43, 0.13, 0.13, 0.64, 0.12 at 2, 3, 4, 5, 6 slabs. The successive-iterate
  test fired on the accidental 3-vs-4 agreement and returned a probability
  wrong by **0.855**, with no warning. Tightening `rtol` does not help -- the
  comparison is between two answers that both failed to see the profile.
  `n_slabs` is now a *floor* on the ladder, in `osc_prob` and in both batched
  scan engines: refinement starts at `max(min_n_slabs, n_slabs)` and only ever
  climbs, clipped at `max_n_slabs` so a floor above the cap raises the existing
  not-achieved warning instead of stepping the ladder back down. With the
  default `n_slabs=1` the floor is inactive and nothing changes. The regression
  test is built on a `solve_ivp` oracle, not on agreement between two `osc_prob`
  grids -- that kind of agreement is what let this through.
- **Notebooks 02 and 03 shipped castle-wall figures drawn from those wrong
  probabilities.** Across a 6000-point baseline scan, 7.2% of points were off
  by more than 1e-2. The scans now get their `n_slabs=150` honoured, and they
  also pass the wall positions as `t_breakpoints`: the profile is a step
  function, and high-order quadrature reaches its nominal order only when the
  Hamiltonian is smooth inside each slab. Against a converged reference the
  worst point improves from 0.855 to 1.0e-3 -- and the scans run faster than
  the wrong version did, because slab edges placed on the discontinuities buy
  more accuracy per slab than piling on uniform slabs.
- **Eight cells in notebooks 02 and 03 taught a performance lesson that is not
  true.** They asserted that an `osc_prob` call costs "~4.5 ms, almost entirely
  fixed entry-path cost", evidenced by "n_slabs=1 and n_slabs=150 both measure
  4.5 ms per call, which is how you can tell the physics is not what costs",
  and predicted that looping instead of passing an array would turn a 2 s cell
  into 45 s. Measured: 0.26-0.46 ms per call, and on a profile that actually
  varies the same call costs 0.34 ms at `n_slabs=1` against 10.4 ms at
  `n_slabs=2000` -- the integration is nearly all of it. The loop-vs-array
  ratios are 1.4x-2.3x, not the ~22x claimed. Passing the array is still the
  right advice and the cells still do it; the numbers and the reasoning behind
  them are now the measured ones.
- **Every notebook that used matter effects was computing vacuum.** Notebooks
  01-10 built the coherent forward potential by calling
  `matter.num_density_e_func` with a density in g cm^-3 but without
  `density_matter_is_in_g_per_cm3=True` -- for constant densities and for
  `earth.density_matter_func_prem` alike. That yields `VCC = 8.8e-32 eV`
  instead of `3.8e-13 eV`, a factor of 4.3e18, so the matter term was ~20
  orders of magnitude below the vacuum one. In notebook 02 the "matter"
  probability came out bit-identical to the vacuum one. 19 call sites.
- **The matter potential was subtracted where it should be added.** For
  neutrinos the library computes `H_vac + h_matt(VCC)`; the antineutrino sign
  flip lives inside `VCC` (`matter.vcc_func_from_rho_func` applies
  `s = 1 if not nubar else -1`), and `matter.VCC_func` -- which the notebooks
  used -- always returns a positive potential. 45 sites in notebooks 02, 03,
  04, 07 and 08 carried a leading minus on a neutrino Hamiltonian. Checked
  against the closed-form `oscprobstd` result: the corrected sign agrees to
  2.4e-14, the old one is wrong by up to 135%. These two bugs masked each
  other -- with the matter term 20 orders down the sign was invisible, and the
  notebooks' own relative-error subpanels read a healthy 1e-12 because the
  standard formula was being fed the same wrong potential.
- **PREM was sampled at the centre of the Earth.** Notebooks 04, 05 and 07
  wrote `VCC_func_prem(r/gd.CONV_KM_TO_INV_EV)` where `r` was already in km,
  evaluating the profile at r ~ 1e-6 km and so using a constant 13.09 g cm^-3
  everywhere instead of the layered profile. The same notebooks used the
  correct `VCC_func_prem(r)` for their density *panels*, so the plotted profile
  disagreed with the physics behind it. Cross-checked against
  `osc_prob_3nu_earth`: the corrected Hamiltonian agrees to 3.2e-6, while the
  old one gave P(nu_mu -> nu_e) = 0.0001 against a true 0.0047. 7 sites.
- **A matter term discarded as dead code.** Eight Hamiltonians across notebooks
  05 and 07 read `return H_vac(...)` followed by a bare `+ H_matt` on the next
  line, which Python evaluates and throws away, so those Hamiltonians were pure
  vacuum. Joined onto the `return`.
- Notebook 07 had never run top to bottom: cell 14 referenced `baseline`, which
  is first assigned three cells later. Five figure titles in notebooks 07 and
  08 named the swept variable instead of the fixed one.

- A matter density that has already been converted to natural units, but is then
  declared to be in g cm^-3, is now flagged with `matter.DensityUnitWarning`
  instead of being converted a second time in silence. The two scales do not
  overlap -- the densest matter anyone models is some 1e15 g cm^-3, while any
  density from water upwards becomes 4.3e18 or more once converted -- so the
  check has three orders of magnitude of margin. It is worth having because the
  consequences do not look like a unit error: the matter term swamps every other
  scale, nu_e becomes an exact eigenstate, and the calculation returns a
  perfectly self-consistent `P_ee = 1`, which reads as a broken formula rather
  than as bad input.
- Passing a whole entry of `globaldefs.OSC_PARAMS_PREDEFINED` to a probability
  function, as `**OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']`, now raises a
  `ValueError` naming the two offending keys and pointing at
  `globaldefs.load_nufit_params`. Those entries carry `name` and `description`
  strings alongside the six mixing parameters; unchecked, they travelled the
  shared `**kwargs` chain until `magnus_expansion_multislab` rejected them,
  naming the one function in the chain with nothing to do with the mistake. The
  check sits in the four middle-layer functions rather than further down,
  because the averaged and ordinary paths diverge before the Magnus core is
  reached -- a guard placed later caught the strings on one path and ignored
  them on the other.
- `matter.vcc_func_from_rho_func` raised `TypeError: 'float' object is not
  callable` when given a constant electron number density together with
  `density_is_of_number_of_electrons=True`. `VCC_func` evaluates whatever it is
  handed at a position, and the sibling branch for a matter density wraps its
  constant in a function for exactly that reason; this one passed the bare number
  through. The documented way to supply a fixed electron number density therefore
  did not work at all. Found by the first test to reach the line.
- The two guards in `validate_input_osc_prob_earth` that reject a malformed
  `loc_ini`/`loc_fin` caught `KeyError`, which unpacking a sequence never raises,
  so neither could fire. A three-entry coordinate escaped as `too many values to
  unpack (expected 2)` and a non-iterable as `TypeError` -- the latter breaking
  the convention, established across the rest of the package, that bad input
  raises `ValueError` with a message naming the parameter at fault. They now
  catch `TypeError` and `ValueError`. Found by writing the first test that ever
  reached them.
- `osc_prob` raised `UnboundLocalError` instead of returning a probability when
  `max_num_loops < 1` was passed together with `validate_input=False` (the
  validator rejects that combination otherwise). The refinement-limit checks at
  the top of the loop could `return P` before the first iteration had produced
  one; they are refinement limits and only mean anything after a loop has run,
  so they are now guarded on that. The dead `iterate_over_magnus_exp_order`
  dispatch had been assigning `P` early, which hid the problem from static
  analysis until it was removed.
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
