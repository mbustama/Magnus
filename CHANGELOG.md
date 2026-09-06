# Changelog

All notable changes to Magνs are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [1.1.0] - 2026-09-06

### Changed

- **The sterile block of the matter projector now follows the Earth's layered
  composition.**  The matter Hamiltonian factorizes as `V_CC(l)` times a
  projector, and the density has always resolved the electron fraction per PREM
  layer -- but the projector was one constant matrix taking a single scalar
  `r = n_n/n_p` for the whole trajectory.  The two therefore described different
  media on any chord whose composition varies, which is every chord that reaches
  the core.

  Only the sterile states are affected: the active flavours all feel the same
  neutral-current potential, so it drops out as a phase, while a sterile state
  feels neither current and is left carrying `-V_NC = (r/2) V_CC`.  Three flavors
  are therefore untouched, exactly -- the projector's sterile block is empty
  there and the constant matrix is still returned.

  The cost of getting this wrong was larger than the shipped warning said.  Off
  resonance it is the ~2e-02 the warning quotes, but near the sterile matter
  resonance on a core-crossing chord it reaches **0.4 in probability**, flat
  under refinement, invariant in the mass splitting, surviving small mixing
  angles, and present at 3+2 as well.  No choice of scalar fixes it: the best
  value obtainable by scanning still leaves 7e-03 in that band.

  It costs nothing to fix.  Every Earth entry point already declares the PREM
  layer boundaries as slab edges, automatically, so each segment is homogeneous
  in composition and a per-segment projector is exact on the grid that exists --
  no regridding and no extra Hamiltonian evaluations.  Measured at parity with
  the constant-projector path.

  **This changes what an Earth call returns by default at four and five
  flavors.**  `ratio_number_neutrons_to_protons` defaults to `None` on the twelve
  Earth wrappers, resolving to the layered composition; passing a scalar
  explicitly reproduces the old behaviour exactly.  No signature changed: the
  parameter accepts a callable as well, and the general (non-Earth) entry points
  keep their scalar default, having no layered composition to resolve.

  `SterileMatterCompositionWarning` is retargeted rather than removed.  The
  default and a caller-supplied callable are silent; a scalar over a layered
  profile warns unconditionally, the previous 2% threshold having been measured
  passing chords whose error equalled the figure the warning quoted.

  `notebooks/sterile_projector_check.py` reproduces the defect and the fix in one
  command, three arms differing only in the projector.

## [1.0.13] - 2026-09-06

### Fixed

- The `HiddenFeatureWarning` told callers to pass `t_breakpoints` and then
  printed ones that did not work.  The suggestion was a pair-scale bracket, ten
  to a hundred times wider than the feature it was meant to straddle: a single
  point stayed at 3.0e-02, and a 60-point scan moved to 5.8e-02 **with no
  warning at all**, because passing `t_breakpoints` switches the scan the
  detector runs in off.  A caller who did what the message said therefore landed
  back where the warning exists to rescue them from, with a worse error and
  nothing on screen.

  The suggested edges are now localized: the flagged interval is re-sampled and
  the sub-interval carrying most of the variation is bracketed instead.
  Following the message takes the same cases from 3.0e-02 to 1.0e-04 and
  1.0e-02 to 9.9e-05, at a point and over a 60-point scan alike.  The
  localization runs only where the warning fires, so the detector's
  false-positive record is untouched by construction.

  The message also now says what to expect afterwards: breakpoints route the
  call to the general slab ladder, so accuracy lands near the requested
  tolerance rather than below it.

  Detection itself is unchanged.  The undetectable remainder is structural --
  a feature narrower than the reference grid cannot be found by refining a grid
  -- and no attempt was made to change that.

## [1.0.12] - 2026-09-05

### Changed

- The baseline scan folds its running product in a compiled kernel rather than a
  Python loop.  :func:`osc_prob_energy_baseline` answers every requested baseline
  from one traversal of the profile, because each answer is a prefix of the next;
  the loop that walked it snapshotted the running product at each requested
  distance, one Python iteration per slab.  Marginal cost per slab falls from
  about 1.2 us to 0.011 us at two flavours.

  End to end the gain grows with the grid, since the fold's share of the scan
  does: about **2.7x at 512 accumulation steps and 7.5x at 32 768** at two
  flavours, 1.5x and 3.9x at three.  It exceeds the survey's own 2-4x estimate
  for a reason worth recording -- the kernels of 1.0.6 through 1.0.11 had made
  building the operators so cheap that the Python fold had become the dominant
  cost of this scan.

  This needed its own kernel rather than reusing 1.0.11's.  That one folds
  ``acc <- acc @ U[k]`` with k descending; this one needs ``U[k] @ acc`` with k
  ascending, which is the mirror image -- a different product, not a different
  parenthesization -- and it must snapshot mid-fold at requested baselines, which
  a final-answer kernel cannot do.

  **Not bit-identical**, for the same reason as 1.0.11: NumPy routes these
  products through MKL's `zgemm`, which uses fused multiply-add, and a kernel
  compiled with `fastmath` off cannot reproduce that ordering.  Worst shift
  1.4e-14, which is twenty times smaller than 1.0.11's.  Unlike that release,
  neither side is systematically closer to exact -- against 40-digit mpmath the
  new fold wins at two flavours and the old at three -- so this is rounding
  exchanged rather than accuracy gained.  A numba-less install is exactly
  bit-identical to the old loop.  The scan-against-per-point invariant was
  checked on identical slab edges at all four flavour counts, 4.9e-15 to
  2.2e-14, against the suite's 1e-12 bar.

## [1.0.11] - 2026-09-05

### Changed

- The interaction-picture engine folds its slab operators in a compiled
  accumulator rather than a Python loop.  That loop ran one iteration per slab,
  each doing a single small matrix product and returning to Python, and measured
  56% of the engine's pass; the fold itself falls 18.3x, the engine 2.28x.  What
  a caller sees depends strongly on how many slabs are folded -- about 2.3x at
  32 768 slabs, and nothing measurable at 512 -- so the figure is a property of
  the request, not of the engine.

  This is the fourth loop of its kind to be compiled, after the separable scan's
  in 1.0.6 and the Gauss-Legendre expressions in 1.0.8, and it reaches a narrower
  set of calls than any of them: the engine serves only two-flavour requests on a
  profile tagged as exponential, with no supplied slab edges.

  **Not bit-identical**, and established as unachievable before it was written
  rather than discovered afterwards: NumPy routes these stacked 2x2 products
  through MKL's `zgemm`, which uses fused multiply-add, and a kernel compiled with
  `fastmath` off cannot reproduce that ordering.  Worst engine-output shift
  2.8e-13 over 34 configurations spanning slab counts 8 to 32 768 and tolerance
  ladders 1e-3 to 1e-9, with **every certification decision unchanged** and
  returned probabilities identical in end-to-end testing.  Scored against
  40-digit mpmath folds of the engine's own operators, the compiled fold errs at
  worst 2.4e-14 where the BLAS chain it replaces errs at 8.8e-14: the shift is
  dominated by the rounding of the code being removed, and the new fold is three
  to seven times closer to exact.  A numba-less install is exactly 0.0 from the
  old code.

## [1.0.10] - 2026-09-05

### Changed

- `_expm_stack` asks whether a stack is anti-Hermitian in one compiled pass
  rather than five full-stack temporaries and about seven traversals, and builds
  `K = 1j*Om` only after the branch is taken.  The framing stage falls 4.3-5.3x;
  a whole call falls about 1.1x, the difference being that the framing is only
  6-17% of a call.  Probabilities are unchanged: exactly 0.0 difference across
  128 configurations spanning both profiles, two to five flavours, four slab
  counts and all four Magnus orders.

  The kernel compares *squared* magnitudes and takes one `sqrt` at the end, so
  `scale` and `dev` can differ from NumPy's by up to 2 ulp.  Computing `abs()`
  per element instead is exact in intent but measured 0.47-0.78x -- slower than
  the NumPy expression it would replace -- and is not value-identical either, so
  there is no exact compiled alternative to prefer.  Those two scalars are used
  only to choose between the identity branch, the anti-Hermitian branch and the
  scipy fallback; the branch could change only for an input whose deviation-to-
  scale ratio sits within 2 ulp of 1e-12, four decades from where real input
  lands.  Neither value reaches a warning or a returned quantity.

  NaN and inf route exactly as before.  A plain maximum loop would not have done
  that -- it skips NaN where `numpy.max` propagates it, which would have sent a
  NaN stack to the eigensolver instead of the scipy fallback -- so the kernel
  carries a running sum of squared magnitudes, which vectorizes and cannot
  cancel, and checks it once at the end.

## [1.0.9] - 2026-09-05

### Changed

- Four- and five-flavour exponentials go to a batched Jacobi eigensolver rather
  than `numpy.linalg.eigh`.  `supports_dim` now answers yes for 2 through 5.
  The reasoning it used to carry -- that a 4x4 or 5x5 Hermitian eigenproblem has
  no practical closed form, so those dimensions keep `eigh` -- had a true premise
  and a conclusion that did not follow: what made them slow was never the missing
  closed form but `eigh`'s fixed per-matrix LAPACK overhead, about 2.3 us on a
  4x4, two thirds of a whole d = 4 pass.  The kernel is a cyclic complex-Hermitian
  Jacobi sweep that warm-starts each matrix from its predecessor's eigenvectors --
  consecutive matrices are consecutive slabs of one energy, so they arrive nearly
  diagonal -- and re-orthonormalizes that basis at every step, which is
  load-bearing rather than tidy: without it non-unitarity compounds along the
  chain, 2.3e-11 after 13 000 matrices against 3.9e-14 with it.

  Marginal cost per slab falls **1.81-1.95x at four flavours and 1.46-1.59x at
  five**, end to end on both profiles, arms interleaved in one process.

  **Unlike every other kernel in this series, this one is not bit-identical**, and
  it is not meant to be: it replaces an iterative eigensolver with a different
  iterative eigensolver.  Worst full-matrix probability shift measured 1.557e-12.
  It is held instead to `eigh`'s own accuracy class, and to the same 10x bar the
  suite already applies to the closed-form kernels: across 364 cells spanning the
  norm sweep from 1e-150 to 1e4, crossed with clustered spectra down to exactly
  degenerate, its error stays within 6.4x of `eigh`'s, both scored against
  `scipy.linalg.expm`.  Two and three flavours are untouched -- `max|dP|` is
  exactly 0.0 there, and their timings do not move.

  Being backward stable with no conditioning cliff, it needs no `SEV_TOL` gate of
  the kind the closed forms require; it returns sev = 0.0, escalating to `eigh`
  only if a matrix fails to converge in 30 sweeps, which no census has observed.
  It terminates when a sweep performs no rotation, at a per-element threshold of
  1e-30 relative.  That threshold is not the prototype's: the original squared
  off-norm floor sat inside the rounding equilibrium, so long chains declined to
  `eigh` almost always and the speed-up vanished silently.

## [1.0.8] - 2026-09-05

### Changed

- The order-4 Gauss-Legendre `Omega` is built by a compiled kernel that makes one
  pass over each slab's two node samples, fusing three steps the NumPy form pays
  separately: the constant-sample equality test, the commutator, and the linear
  combination.  The two samples are strided views, so the NumPy route copies both
  before the commutator can read them; the kernel reads them where they lie, and
  the equality scan stops at the first differing element rather than comparing
  every one.  Marginal cost per slab falls 1.52x at two flavours and 1.23-1.36x at
  three, less at higher flavour counts where the eigendecomposition dominates.
  Orders 2, 6 and 8 are untouched and keep the commutator kernel of 1.0.7.

  Unlike 1.0.6 and 1.0.7, this one introduces **no** new divergence between a
  numba-less install and a numba one: the kernel reproduces the compiled
  commutator's accumulation order and the same association of the scalar factors,
  so its output is bit-identical to the expression it replaces -- exactly 0.0
  difference across every shape, flavour, width form and edge case tested.
  Anything the kernel cannot take -- another dtype, another node count -- falls
  through to that expression unchanged.

- Orders 6 and 8 get the same treatment.  Order 6 collapses about nineteen
  full-stack temporaries into eight `d x d` scratch buffers hoisted out of the
  slab loop; order 8, about sixty-five into nineteen.  NumPy evaluates each
  operation across every slab before starting the next, so each intermediate is
  a full array written to memory and read back, while the kernel finishes one
  slab before moving on.  Marginal cost per slab falls 1.29x at order 6 and
  1.75x at order 8 end to end, median over eight profile/flavour cells.  Both
  are bit-identical to the expressions they replace, and the six chained
  commutators of order 8 are each stored where NumPy stores one rather than
  algebraically re-folded.

  Two rounding differences had to be worked around to reach that, both worth
  knowing if these kernels are ever extended.  NumPy divides a complex array by
  a real scalar through Smith's algorithm, which *multiplies by the reciprocal*;
  numba divides componentwise, and about two thirds of two million random values
  differ in the last bit.  Every `X/c` in these kernels is therefore written
  `X*(1.0/c)`.  Separately, numba's `**` does not round as
  `numpy.float64.__pow__` does -- 26% of random doubles differ in the last bit
  at `**3` -- so the order-8 weight-times-power scalars are precomputed as module
  constants in Python rather than raised inside the kernel.

  Order 2 was prototyped, measured and deliberately left alone: its branch is a
  single ufunc call with no boundaries to fuse, and at two flavours a kernel
  loses outright (0.82-1.04x per call), the dispatch's own fixed cost exceeding
  what streaming saves on a branch costing 0.01-0.05 us/slab.

- The set of accepted pass-through keyword names is computed once and cached
  rather than rebuilt by `inspect.signature` on every public call.  It was being
  rebuilt twice per call, once by the entry point and once by the `osc_prob` it
  delegates to, at about 221 us each; a single-point call is that much cheaper.
  This is fixed overhead, so it is invisible on large grids and worth most to
  callers making many small calls.  `cache_clear` remains available for tests
  that add or remove parameters at runtime.

## [1.0.7] - 2026-09-05

### Changed

- The Magnus term recursion and the Gauss-Legendre schemes compute their
  commutators in a compiled kernel that fuses both matrix products into one
  accumulation.  `magnus.commutator` itself is unchanged and stays the general
  pure-NumPy form; the kernel is a private path, and anything it cannot take
  falls through to the same expression.  Marginal cost per slab falls 3.1-3.4x
  at two flavours and 2.0-2.2x at three at order 4, less at higher flavour
  counts, and not at all on the cumulative path.

  As with 1.0.6, a numba-less install is no longer bit-identical to a numba one
  on these paths: worst observed shift 6.7e-14 across 36 configurations, every
  refinement decision unchanged.

## [1.0.6] - 2026-09-05

### Changed

- The separable energy scan composes its slab operators in a compiled kernel
  rather than a Python loop.  The association is unchanged -- the same left
  fold, earliest slab rightmost -- so this is not a reassociation; the kernel
  accumulates each matrix element as a compiled scalar sum where BLAS orders
  the same arithmetic its own way.  Marginal cost per slab falls about 3.5x at
  three flavours and 2.3x at two; the gain at four and five flavours is smaller
  and was not resolved on the machine used.  Without numba, or on a dtype the
  kernel was not built for, the original loop runs.

  One consequence is worth knowing: a numba-less install was bit-identical to a
  numba install on this path and may now differ at the 1e-14 level -- worst
  observed 1.28e-14 across 16 scan configurations, with every refinement
  decision, warning and slab count unchanged.  See `docs/source/performance.rst`.

## [Unreleased]

### Added

- Order-8 Gauss-Legendre collocation integrator, on four nodes.  `'gl'` now
  reaches orders 2, 4, 6 and 8 from 1, 2, 3 and 4 Hamiltonian evaluations per
  slab; `MAGNUS_EXP_ORDER_MAX_GL` rises from 6 to 8, and a request above 8
  raises as before.  Orders 6 and 8 use the commutator-optimal forms of Blanes,
  Casas & Ros, *BIT* **42**, 262 (2002), needing three and six commutators --
  the fewest possible at each order.  Verified against a `DOP853` reference:
  local error converges as `h^9` (measured slopes 8.3-9.7 over three
  Hamiltonians), and on a PREM chord the row sums hold to 1.8e-15.

### Changed

- The Gauss-Legendre schemes are described as *collocation* integrators
  throughout.  Several places had called them "commutator-free", which is the
  distinct family of Blanes & Moan (2006) and Alvermann & Fehske (2011) that
  replaces commutators with products of exponentials; these schemes are built
  from commutators, three at order 6 and six at order 8.

## [1.0.0] - 2026-08-11

### Fixed

- **The sterile states felt no matter at all on the NSI route.**
  `osc_prob_matter_nsi` built the standard part of its matter Hamiltonian as a
  literal `diag([1, 0, 0, 0])`, so the sterile states got zero where they carry
  `-V_NC = (r/2) V_CC`.  With every NSI coupling set to zero the route has to
  reproduce `osc_prob_matter_std_potential` exactly, and did not: it differed by
  **5.2e-02 at four flavors and 5.1e-02 at five**.  Three flavors agreed all
  along, which is why nothing caught it -- the only test that can see it needs a
  fourth state.  Anyone who has published 3+1 or 3+2 NSI-in-matter numbers from
  an earlier build should re-run them.

- **The solar LIV routes read the Sun's electron density as a mass density.**
  `osc_prob_Nnu_sun_liv` builds its profile from `NUM_DENSITY_E_SUN_CENTRAL`, an
  electron *number* density, but forwarded a `density_is_of_number_of_electrons`
  flag defaulting to False.  Against `osc_prob_Nnu_sun` with the couplings zeroed
  they differed by **0.69, 0.69, 0.45 and 0.43** at 2, 3, 4 and 5 flavors.
  `osc_prob_Nnu_sun` and `..._sun_nsi` never exposed those flags, which is why
  only the LIV family was affected.

- **A negative energy returned the antineutrino probability.**  `E < 0` flips the
  sign of the whole Hamiltonian, which is CP conjugation, so the call returned a
  unitary, entirely plausible answer to a question the caller did not ask --
  matching `nubar=True` to 1e-15 and differing from the intended answer by
  2.3e-02 in vacuum and 4.4e-02 in matter.  Non-positive energies are now
  rejected.

- **A NaN density was reported as a units mistake**, because `nan == 0.0` and
  `nan >= threshold` are both False and the guard fell through to its warning.
  Non-finite densities are now rejected, naming the real problem.

- **Notebook 12's four- and five-flavor ground truth** built its comparison
  Hamiltonian with the same zero-sterile literal, so `solve_ivp` was integrating
  the wrong problem and the error column blamed the strategy for the reference.
  The four-flavor standard case goes from `err_magnus = 1.73e-04` to 2.19e-06.

- **The quickstart's and README's code examples could not be run** -- one used
  names the page defines fifteen lines later, others passed `...` as a
  Hamiltonian body.  Both are now executed by the suite.

- **The documented solar averaging reduction outlived the script it came from.**
  The averaged-probability page, notebooks 13, 14 and 23 and the tutorials index
  all asserted "instantaneous 1.4e-03 -> averaged 2.6e-05, a **53x** reduction",
  citing `adversarial_batteries/avg_check.py`.  That script had gone stale at
  `0bf3a40`, which is on `main` and predates this audit; re-run, it gives
  6.000e-04 -> 7.110e-04 on the log-linear BS05 ray notebook 13 actually uses
  (8.889e-04 -> 6.051e-04 on the cubic-spline variant).  All four rows of the
  table are now measured values.

  The **rule** stated alongside them was the worse problem and is withdrawn:
  the docs told the reader to read the ratio of the two columns as a diagnostic,
  errors shrinking "more than twentyfold" being phase and the rest envelope.  A
  finite-window mean is an estimator with a bias of its own, and on a profile
  whose density varies across the window that bias does not shrink as the window
  widens -- on the solar ray the window mean drifts 0.5924 -> 0.6023 between six
  and forty-eight oscillation lengths, *away* from a limit rather than towards
  one.  So the ratio is meaningful only at fixed matter conditions, as in
  notebook 23, and the page now says so.  To obtain the averaged probability,
  ask for it: `average=True` evaluates the decohered limit in closed form and
  reproduces the adiabatic MSW expression to **3.33e-16** across 1--20 MeV.

- **The front page's solar timing comparison outlived its measurement too.**  The
  README and the docs index both read "40 averaged energies in 0.66 s, against
  131 s for 12 *instantaneous* ones from nuSQuIDS".  That was measured, but on a
  superseded design of notebook 25's section 10; the frozen
  `external_solar_nusquids.json` now sweeps four solver tolerances over a
  200-point grid targeting 40 energies, and contains no 12-energy run and nothing
  at 131 s.  Both pages now say what the notebook says and what its data
  supports: nuSQuIDS needs about **ten minutes** merely to reach the tolerance at
  which its output is a probability at all (568 s, the cheapest setting the
  generator marks physical), and a further factor of *N* to average the phase
  away.  Magnus's own 0.66 s is a live timing and re-measures at 0.68 s, which is
  run-to-run noise on one machine and is left as it stands.

- **`osc_prob_vacuum` documented twelve refinement knobs as forwarded that it
  deliberately drops.**  `t_slab_edges`, `magnus_exp_order`, `integration_method`,
  `rtol`, `atol`, both `growth_factor_*`, `max_num_loops` and the four
  `min`/`max_n_*` bounds each said "Forwarded to `osc_prob_energy_baseline`
  /`osc_prob`", while the body forwards only `n_jobs`, `validate_input`,
  `verbose`, `save_log` and `file_log`.  Not forwarding them is correct and the
  code says why -- a vacuum Hamiltonian is constant in position, so every point
  is exact with a single slab and there is nothing to refine -- so the
  **docstring** was the wrong text, and it is the one a caller reads before
  passing `rtol=1e-9` and getting silence.  No behavior changes.  `filename_log`
  and `close_file_log_upon_exit` were documented the same way and are likewise
  never read.

  Found by generalizing `tests/test_angles.py`'s "declares the parameter, never
  reads it" AST check from `angles` to every parameter in `src/magnus`.  The
  other 51 hits are the legitimate interface-signature pattern -- `_td` builders
  that take `l` and ignore it, constant-density closures -- or already documented
  as inert, such as `hamiltonian_2nu_liv`'s `nubar`, which has nothing to
  conjugate because the 2-flavor LIV rotation carries no phase.

- **`SEV_TOL`'s documented accuracy guarantee was not the measured one.**  The
  docstring claimed 2e-13 absolute across everything the gate admits; re-measured
  on two spectrum families with 40 random bases per rung, the admitted range
  reaches **2.0e-12**, and even the `m <= 1.1e3` corner the claim was calibrated
  on exceeds 2e-13 in about 1% of bases.  No value of the constant could have
  rescued it, because the test that pins the window requires the gate to sit
  above the very cell where the claim already fails, so the number was corrected
  instead: **5e-12** across the admitted range and **5e-13** in that corner, each
  about twice the worst measured.  The gate itself is unchanged at 1e4 and no
  result moves.

  Two earlier calibrations of this constant appeared to contradict each other --
  one putting the first unsafe cell at `m = 1.1e5`, the other at `4.4e3`.  They
  do not: they used different spectrum families under the same "scale 1e2" label,
  `[-s, -s(1-d), s]` spanning `2s` against `[0, d, S]` spanning `S`, a factor of
  four in `m`.  Compared at equal `m` they agree.  `m = tr(X^2)/6` is a spectral
  invariant and is the only fair label; "scale" is not.

- **The sterile matter term and the Earth's density describe different media.**  L37
  made the Earth's `Y_e` a function of radius and derived `r = (1 - Y_e)/Y_e` from it
  layer by layer, because the two are one statement about composition.  The sterile
  states' entry in the matter projector, `r/2`, could not follow: it is a single matrix
  for the whole chord, so it kept taking the caller's scalar, which defaults to 1.0 --
  isoscalar matter, i.e. `Y_e = 0.5`, precisely the uniform composition L37 replaced.

  So on every Earth chord with a sterile state the density and the projector disagree by
  construction.  Measured at `costhz = -0.95` with `s14 = 0.15`, `s24 = 0.10` and
  `D41 = 1 eV^2`, the isoscalar projector differs from one built with the core's own
  `r = 1.1478` by **2.1e-02** in `P(numu -> numu)` -- twenty times the default tolerance,
  and silent.  Three flavors never saw it, because the projector's sterile block is
  empty; that is the same reason `A2b` survived a max-effort review.

  No single `r` is right for a chord that crosses iron and rock, so this is **reported
  rather than resolved**: the Earth wrappers now raise
  `globaldefs.SterileMatterCompositionWarning` when the two disagree by more than 2%, and
  the message names the path-averaged ratio for the chord that was asked for, which
  silences it.  `electron_fraction=0.5` with the default `1.0` makes them agree exactly
  and reproduces the older uniform composition.  **No number changes unless you act on
  the warning.**

  The default was deliberately left at `1.0` rather than moved to the path-averaged value,
  because no single `r` is right for a chord crossing iron and rock and a self-consistent
  *looking* default would hide that.  The real fix is a position-dependent projector,
  `H_matt(l) = V_CC(l) x P(l)`, which is a structural change to the scalar-times-constant-
  matrix factorisation the matter path is built on; it is tracked for a future release as
  [issue #47](https://github.com/mbustama/Magnus/issues/47).

  The docstrings were also wrong about this.  Twelve Earth wrappers said the ratio "must
  match the value given to `vcc_func_from_rho_func`", which on the Earth path is
  unsatisfiable -- that function is handed a per-layer ratio the caller never supplies.

- **The sterile states' matter entry was unreachable on the Sun, and locked isoscalar.**
  `osc_prob_{4,5}nu_sun` and their `_nsi`/`_liv` variants never exposed
  `ratio_number_neutrons_to_protons`, and delegated without forwarding it -- two of them
  passing a hardcoded `1.0`.  So `r/2`, the sterile entry in the matter projector, was
  fixed at isoscalar for a medium that is nothing of the kind: the Sun is hydrogen-rich,
  `Y_e = (1 + X)/2` runs from about 0.68 at the center to 0.88 near the surface, and
  `r = (1 - Y_e)/Y_e` from about 0.47 down to 0.14.  Unlike the Earth, where 1.0 at
  least sits among the layer values, for the Sun it is outside the physical range
  entirely.

  The six wrappers now take it and forward it.  Left at the default the averaged solar
  survival probability moves by about **4e-03** at `s14 = 0.4`, above the 1e-3 default
  tolerance; at the nominal `s14 = 0.15` it is 7.9e-04, just under.  The solar profile is
  a fit to the electron *number* density, so `Y_e` is already inside it and there is
  nothing for the library to derive `r` from -- it has to be stated, which is why this is
  exposed rather than defaulted to a solar value.  **No number changes unless you set it.**

- **`t_breakpoints` was unusable on every Earth wrapper.**  Those wrappers place slab
  edges on the PREM shell crossings themselves, so an argument of the same name arrived
  in `**kwargs` and collided: the caller got `got multiple values for keyword argument
  't_breakpoints'` raised two layers down.  The keyword is listed as forwardable in this
  package's own unrecognized-keyword message, so it was reachable and broken.  A
  caller's breakpoints are now merged with the PREM crossings rather than either
  replacing the other -- dropping the crossings silently would be the very defect
  `t_breakpoints` exists to prevent.  Pass `t_slab_edges` to place every edge yourself.

- **`n_jobs` is documented as not being a pure performance knob.**  Splitting slabs
  across workers changes the order the arithmetic is done in, and the refinement
  ladder's stopping test compares successive levels, so it can stop one level earlier or
  later than the serial run.  Measured on a 3nu PREM chord over eight energies, serial
  against two workers differs by 1.2e-03 at the default `rtol = 1e-3`, 6.6e-08 at
  `rtol = 1e-6` and 5.6e-11 at `rtol = 1e-9` -- within the tolerance asked for, but not
  bitwise, which the docstring previously gave no hint of.

- **The command line could not reach half its own parameter sets.**
  `--osc-params-set` listed its choices by hand and offered only the NuFit 6.0
  entries, so after the default moved to 6.1 anyone asking for inverted ordering
  on the command line silently dropped a release behind the default.  The choices
  now come from `globaldefs.OSC_PARAMS_PREDEFINED` itself.

### Added

- **An `angles` keyword on every function that takes a mixing angle.**  It selects
  the convention the angles are stated in: `'sin'` (the default, and the only
  behavior before this) their sines, `'sin2'` their sines *squared* -- the form
  global fits are published in -- `'rad'` the angles in radians, or `'deg'` in
  degrees.  Under `'deg'` the CP phases are read as degrees too; under the other
  three they stay in radians, a sine being no way to state a phase.  Ninety-five
  functions take it, across `hamiltonians`, `oscprob`, `oscprobstd`, the
  `magnus prob` command line and `globaldefs.load_nufit_params`, and the four
  routes agree to 2.2e-13 end to end.  The default is a pass-through, so nothing
  that does not ask for it changes.

  It exists so a published parameter set can be typed in as published:
  `s12=0.308, ..., angles='sin2'`, or `s12=33.76, ..., dCP=212.0, angles='deg'`,
  rather than square-rooted and degree-converted by hand at the call site.  It
  also matches the keyword NuOscProbExact uses, so the two codes can be driven
  from one parameter set when they are compared.

  Four guards come with it, because a mis-stated convention is otherwise silent:
  an unrecognized value, a sine outside `[-1, 1]`, a negative `'sin2'`, and an
  angle above `2*pi` under `'rad'` all raise; a whole parameter set under one
  degree with `'deg'` raises `globaldefs.MixingAngleConventionWarning`, since
  `theta_13` at about 8.5 degrees is the smallest angle anyone measures and
  values that small are sines.  `'rad'` and `'sin'` are the one pairing no bound
  can separate -- `theta_12 = 0.589` against `sin theta_12 = 0.556` -- which is
  why `load_nufit_params` takes the keyword too: convert once, at the source,
  rather than by hand between the two.

- **`matter.matter_potential_projector` is exported.**  It is the one definition
  of the matter term's structure, and every place that rebuilt that structure by
  hand instead got the sterile entry wrong.  While it was unexported, autoapi
  did not document it and references to it from public docstrings resolved to
  nothing.

- **Input validation for the engine's own knobs**: non-positive `n_slabs`,
  `min_n_slabs`, `n_tpts_per_slab`, and a floor set above its own ceiling, were
  all accepted and then quietly ignored, so a typo looked like a setting that had
  been honored.  `rtol` and `atol` were already guarded; these now match.

- **A PEP 561 `py.typed` marker**, so the annotations the public API already
  carries are visible to mypy and pyright instead of being discarded.

- **`CITATION.cff`**, so GitHub offers "Cite this repository" and reference
  managers can read the metadata.

### Changed

- **The Earth's electron fraction is now resolved per PREM layer, and results
  change.**  Every Earth entry point assumed `Y_e = 0.5` -- exactly isoscalar
  matter, which nothing in the Earth is.  PREM is a density model and carries no
  composition, so `Y_e = <Z/A>` has to be supplied; it is now taken per layer, at
  radii that are already PREM boundaries:

  | layer | radii [km] | `Y_e` | material |
  |---|---|---|---|
  | core | `r <= 3480` | **0.4656** | iron |
  | mantle | 3480 - 6346.6 | **0.4957** | peridotite |
  | crust | 6346.6 - 6368 | **0.4952** | granitic |
  | ocean | `r > 6368` | **0.5551** | seawater (H has `Z/A = 1`) |

  The correction tracks how much core a chord crosses: at 1 GeV,
  `P(nu_mu -> nu_e)` falls to **23%** of its previous value at
  `cos(theta_z) = -1`, and moves by about 1% at -0.4.  **Anyone with published
  Earth numbers from an earlier build should re-run them**; passing
  `electron_fraction=0.5` reproduces the old uniform composition exactly.

  Each layer is settable (`electron_fraction_core` and friends).  Combining a
  per-layer value with the uniform `electron_fraction` is refused rather than
  silently resolved, and `Y_e` is validated as a fraction in `(0, 1]` -- 0.0 and
  5.0 used to be accepted, returning answers 0.51 and 0.74 away from the default.
  The neutron-to-proton ratio is derived from `Y_e` in the density conversion,
  `r = (1 - Y_e)/Y_e`, so a caller can no longer describe an iron core with
  isoscalar neutrons.

  Two caveats are documented rather than guessed at: the crust value differs from
  the mantle by 0.1%, so it exists for explicitness rather than effect; and PREM's
  ocean is a global average that a land-based baseline does not cross, for which
  `electron_fraction_ocean=Y_E_CRUST_PREM` is the right setting.

- **The solar LIV routes no longer accept four parameters they ignored.**
  `osc_prob_Nnu_sun_liv` took `electron_fraction`,
  `ratio_number_neutrons_to_protons` and two density flags and used none of them:
  passing `electron_fraction=0.25` changed the answer by exactly zero.  The solar
  profile is the standard exponential fit to the electron *number* density, so the
  mass-density conversion those describe never runs and `Y_e` is already inside
  the fit.  `osc_prob_Nnu_sun` and `..._sun_nsi` never exposed them.

- **The default oscillation parameters are NuFit 6.1, and there is now only one
  set of them.**  Omitting oscillation parameters used to fall back to a second
  copy of the numbers built from NuFit 6.0 constants, while
  `load_nufit_params()` with no arguments returned 6.1 -- so the same script got
  different answers depending on which door it came through, by **4.0e-03** in
  probability at 1 GeV over 1300 km.  `OSC_PARAMS_DEFAULT` is now derived from
  `load_nufit_params`, so the two cannot drift apart again.
  `OSC_PARAMS_NU_FIT_6_0_SK_NO` and `..._SK_IO` remain available by name for
  anyone reproducing an earlier number.  **Results that relied on the implicit
  default will change.**

- **The publish workflow gates on the tests, on `twine check --strict`, and on
  the release tag matching the packaged version.**  It previously built and
  uploaded on a release with no check of any kind, and PyPI does not allow
  re-uploading a version.

- **CI tests Python 3.13** as well as 3.10-3.12, since `requires-python` has no
  upper bound and pip will install on it either way.


- **Notebook 27, nine animated scenes, with the clips committed.**  Four are the
  scenes NuOscProbExact's notebook 19 draws, computed here so the two can be read
  side by side; the other five need something a closed-form slab code does not
  have -- a refinement ladder deciding it has converged, a front that travels, an
  observable that is an average rather than a value, and a Hamiltonian that varies
  along the path.  The notebook draws stills by default and hides the rendering
  behind `RENDER = True`, so CI never pays the hour it costs.  Seven shrunk GIFs
  are tracked in `img/` (14.7 MB); the raw renders go to `img/raw/`, which is
  gitignored, so a later render cannot silently replace 14.7 MB of committed files
  with the 224 MB it produces.  `tools/make_demo_video.py` owns the encoding.

- **Notebook 25 became an arbiter of when each code wins**, rather than a list of
  timings: reach (a slab product's error floors at 2.5e-11 on a smooth profile and
  then rises, while the Magnus expansion continues to 2.9e-13), generality (five
  flavors, where there is no comparison to draw), and pre-packaged observables
  (a solar average in 0.66 s against 130 s).  With it, supernova-shock comparisons
  against other codes, 3+1 and NSI in both the solar and shock settings, and
  probability-vs-energy panels for the shock.

- **A README image gallery**, ten figures lifted from the executed notebooks by
  `extract_gallery()`, so the front page shows the answers rather than describing
  them and cannot drift from what the notebooks produce.

- **An expansion-order section in notebook 24**: what the truncation order buys in
  accuracy, and what it costs in correctness, which is nothing -- every truncation
  lives in the Lie algebra, so the operator is unitary exactly rather than to the
  accuracy of the truncation.  With it, which engine the dispatcher picks and why,
  measured across four decades of tolerance.

- **`t_breakpoints`, `n_slabs` and `cumulative` reach the BSM wrappers.**  The
  keywords that decide a hard profile were in no signature that a caller of
  `osc_prob_matter_nsi` or `osc_prob_liv` could see, so the comparison those
  wrappers exist for could not be made on equal terms.

- **A `'constant'` engine: a position-independent Hamiltonian is answered in one
  batched exponential instead of one `osc_prob` call per point.**  When the
  matter potential does not vary with position, the Magnus series *terminates at
  its first term* — Ω₁ = −iHΔ and every higher Ω is a nested commutator of H
  with itself, hence zero — so `U = exp(-iHΔ)` is the exact answer and a whole
  energy scan is one stacked exponential.

  This case was previously turned away on purpose: `_osc_prob_scan_separable_dispatch`
  bailed on `not isinstance(VCC_func, Callable)`, its docstring saying "a constant
  potential falls back to the generic path".  So the easiest Hamiltonian there is
  took the slowest route available — a 60-energy scan made **18,000 `osc_prob`
  calls per 300 repetitions**, each rediscovering the same constancy and paying
  the full wrapper and refinement-ladder overhead.

  Measured against the route it replaces, interleaved with a control that came
  back at 1.00×:

  | flavors | matter scan | vacuum scan | single point |
  |---|---|---|---|
  | 2ν | **17.3×** | **24.7×** | 2.0× |
  | 3ν | **15.5×** | **18.9×** | 2.1× |
  | 4ν | 7.2× | 7.4× | 1.4× |
  | 5ν | 6.0× | 6.2× | 1.4× |

  4ν and 5ν gain less because they are on `eigh` rather than the
  Cayley–Hamilton kernel, which covers dimensions 2 and 3 only.  A 3ν
  constant-density scan is now **1.10 µs per energy against NuOscProbExact's
  1.44 µs batched and 13.25 µs looped**; a single point is 33.8 µs against its
  19.9 µs, the remainder being wrapper parameter resolution rather than
  arithmetic.  Results are bit-identical to the per-point route on every
  flavor count and both neutrino signs, and `n_slabs`, `n_tpts_per_slab`,
  `t_breakpoints` and `rtol`/`atol` are accepted and ignored because they can
  only ask for refinement of something already exact.

  **PREM and exponential profiles are untouched** and keep `separable`/`magnus`
  /`hybrid`: their potential varies with position, and a constant-H engine that
  captured one would propagate the whole trajectory with a single exponential —
  wrong by O(1) while still perfectly unitary.  A test asserts the engine
  identity, not merely the numbers.


- **A compiled Cayley–Hamilton backend for the matrix exponential, selected by
  `magnus.magnus.EXPM_BACKEND`.**  `np.linalg.eigh` costs ~1.27 µs per 3×3
  *whatever the stack size* (measured 1.268 µs at N=108, 1.279 µs at N=4096 —
  flat, because it loops over LAPACK internally instead of vectorizing).  The
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
  geometry, not detected: detecting it would need the very evaluations the optimization
  skips.  There is deliberately no user-facing way to declare it of an arbitrary
  profile.

  **This moves Earth single-point results by up to 8.6e-15 relative.**  The mirrored
  slab's nodes are reached by a different floating-point expression for the same real
  number, so the change is inherent rather than incidental.  `USE_PALINDROME = False`
  reproduces the previous numbers exactly.

### Changed

- **Per-call overhead cut across the wrappers, by caching what is pure and
  cheapening what is common.**  The largest single item in a single-point profile
  was `hamiltonian_3nu_vacuum_energy_independent` at ~15 µs, rebuilding the same
  PMNS matrix on every call; it is a pure function of eight scalars and is now
  memoized, handing back a copy so a caller writing into the result cannot
  poison the cache.  Likewise the constant-density branches of
  `matter.vcc_func_from_rho_func` (the callable branches are deliberately *not*
  cached: they return a closure over `rho_func` that callers tag with
  `is_exp_density_profile`, so caching those would trade microseconds for an
  aliasing bug).  `_n_required_params` cached `inspect.signature` weakly against
  the function — 42 µs of cumulative time per call, and once per point on the
  routes that legitimately loop.  `isinstance(x, typing.Callable)`, which routes
  through `ABCMeta.__instancecheck__`, replaced by the `callable()` builtin at 23
  sites (~9 `typing.__subclasscheck__` calls per invocation), and scalar fast
  paths added to `_normalize_energy_L` and the density-units guard.  No
  behavior changed by this entry.

- **A verbose run (`verbose >= 1`) takes the per-point route.**  The banner and
  run-parameter dump describe quantities the batched engines do not have —
  `magnus_exp_order`, slab counts, tolerances — so emitting them from a batched
  path would report a refinement ladder that never ran.


- **``rtol``/``atol`` are documented for what they are: a stopping criterion,
  not an accuracy guarantee.**  The ladder halts when two successive levels
  agree; it never estimates the error of the answer it returns, which is a
  weaker promise than a stepping ODE integrator's ``rtol`` makes.  Corrected
  in ``osc_prob``, ``adiabatic.hybrid_propagator``, the CLI's ``--rtol``/
  ``--atol`` help, ``README.md``, ``architecture.rst`` (which said "until
  rtol/atol is met"), and a new section of ``implementation_details.rst`` that
  the others link to.  No behavior changed by this entry.

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


- **`oscprob.BATCH_WORKING_ENTRIES` lowered from 4,194,304 to 65,536** (about 67 MB to
  about 1 MB).  The batched scan engines are memory-bound, and the previous value was
  large enough that their working set spilled cache.  Measured across fifteen workloads
  on three engines, the new value is **1.19x-1.38x** quicker on Earth energy scans,
  1.06x-1.16x on cumulative baseline scans, flat on short scans and on the
  interaction-picture engine, and never slower anywhere.  **Bit-identical** at every
  budget tested -- tiles are independent and only concatenated -- so this changes no
  result.  Peak memory of a long scan drops accordingly.

### Fixed

- **The memory guard read the host's free memory, not the cgroup's.**  Inside a
  container the two are unrelated, so a request that would be killed by the cgroup
  limit was waved through by a guard reading a number that did not apply to it.
  `_available_memory_bytes` now takes the minimum of the two, walking the cgroup
  ancestor chain and handling the v1 "no limit" sentinel.

- **The sterile state felt no medium, in six places.**  Four inline copies of the
  matter projector, plus the two flavor-specific builders, wrote the 3+1 matter
  term by hand and gave the sterile state a zero where it carries `-V_NC`.  All
  six now come from `matter.matter_potential_projector`, and a test fails if a
  seventh copy appears.  On a PREM chord this was worth 0.29 in probability, and
  it was flat in the requested tolerance, so no amount of refinement revealed it.

- **The evaluation-mode cache answered a question it was never asked, and the
  guarded version of it was dead code.**  `probe_eval_mode`'s `'constant'`
  verdict means "sampling A across [t0, t1] gave the same matrix every time" —
  a property of the function *on an interval*, not of the function.  Keyed on
  the function alone, a mode learned on a short baseline was served for a long
  one, and `_evaluate_A`'s `'constant'` branch then broadcast a single sample
  over a profile that genuinely varies, with no spot-check (unlike the
  `'vector'`/`'scalar'` hints, which self-heal).

  A two-layer profile written with the natural short-circuit — if every
  requested position falls in one layer, return that layer's matrix — probes as
  `'constant'` on a short interval and `'vector'` on a long one.  Looping
  baselines shortest-first then gave `P_ee = 0.906249` against a correct
  `0.903424`, row sums exactly 1.0, no warning, and **the right answer if the
  loop ran the other way round**.  Now keyed on `(function, t0, t1)`: measured
  loop-order dependence **5.8e-02 → 0.00e+00**.

  Three further defects were the same defect.  `oscprob._eval_mode_for` was a
  second, *unguarded* copy of this cache that reached across a module boundary
  into `magnus._EVAL_MODE_CACHE`; it guarded its store but not its lookup, so a
  callable Hamiltonian that cannot be weakly referenced (`__slots__`) or hashed
  (a `@dataclass`, or anything defining `__eq__`) raised `TypeError` from inside
  `osc_prob` where it had worked before.  The public `magnus.cached_eval_mode`,
  which *was* guarded, had zero live callers: its only call site sat in the
  `not callable(H_func)` arm of a ternary, and `H_func` is unconditionally
  rebound to a closure above it.  The copy is deleted, the original does the
  work via a new `key` parameter, and the comment claiming the cache keys by
  identity is corrected — `WeakKeyDictionary` keys by the referent's *equality*,
  which matters for a Hamiltonian defining `__eq__`.

  The interval key costs nothing measurable: a refinement ladder calling
  repeatedly at one interval still probes once, which is what the cache is for.


- **A second max-effort review, run from a fresh session, found four more; all are
  fixed here.**  The first review below was written *and* verified by the agent
  that wrote the code, which is the reason to look again.  None of the four
  changes a computed probability: they are a lost warning, an inconsistent error
  contract, an ignored backend request, and an unbounded cache.

  *The `DensityUnitWarning` repair in the review below was half a fix, and its
  other half was a new defect.*  `vcc_func_from_rho_func`'s constant-density cache
  skips the conversion that the two unit guards live inside, and the earlier
  repair re-emitted only the `density_matter_is_in_g_per_cm3=True` arm.  The arm
  it dropped is the dangerous one: an undeclared g cm^-3 density returns *exactly*
  the vacuum probability, so the warning is the only thing separating it from an
  answer.  Re-emitting from the cache site cannot work either — `warnings.warn`'s
  `stacklevel` attributes the call to a different frame, and the frame is part of
  the interpreter's registry key, so the imitation printed a *second* warning
  under the default filter where an uncached call printed one.  A density that
  would trip either guard is now not cached at all, so both keep firing from where
  they always did.  Measured against `main`: 3 of 3 identical calls warn under
  `simplefilter('always')`, 1 of 3 under the default filter, matching in both
  directions.

  *The constant engine accepted three refinement parameters the ladder rejects.*
  `max_n_slabs=0`, `rtol`/`atol` ≤ 0 and `max_num_loops=0` were all answered here
  while `osc_prob` raised `ValueError` for each, so whether a bad parameter was
  reported depended on whether the density happened to be constant.  The answers
  were never wrong — one exponential per point is exact either way — but the error
  contract was.  `_refinement_params_rejected` mirrors `osc_prob`'s validation and
  *declines*, so the caller sees `osc_prob`'s own message rather than a second
  wording of it.

  *`EXPM_BACKEND` did not cross a process boundary.*  loky re-imports magnus in
  every worker, where the switch is back at `'auto'`, so any call with
  `n_jobs != 1` silently ignored it — and it is the only backend control the
  `oscprob` wrappers expose, the `expm_backend` parameter reaching no further than
  the Magnus layer.  Worst for the one use the switch is documented for: a backend
  comparison run in parallel compared `'auto'` against itself.  Now carried by
  value into the worker and re-applied there.

  *The evaluation-mode cache's per-interval dictionary was unbounded.*  Weak keys
  bound the outer map, not the inner one, and a Hamiltonian defined at module
  scope never dies: 1000 distinct baselines through `osc_prob` retained 1000
  entries, about 184 KB, for the life of the process.  Now bounded at 256 and
  cleared wholesale, matching `_VACUUM_H_CACHE` and `_VCC_CONST_CACHE` — the
  ladder holds one entry no eviction can reach, and a scan uses each entry once,
  so neither population rewards a smarter policy.

  *A refinement bound named the wrong parameter, and had since long before this
  branch.*  Found while mirroring `osc_prob`'s validation for the fix above: one
  of its two ceiling checks tested `max_n_slabs` while its message named
  `max_n_tpts_per_slab`.  So `max_n_tpts_per_slab` was never validated at all —
  `0` and `1` were accepted — `max_n_slabs` was bounded at `> 2` while the message
  three lines above promised `> 1`, and a caller who passed `max_n_slabs=2` was
  refused in the name of a parameter they had not touched.  The rule both messages
  encode is that each ceiling clears its own floor (`min_n_slabs` defaults to 1,
  `min_n_tpts_per_slab` to 2), so the condition was the wrong half: it now tests
  `max_n_tpts_per_slab`.  `max_n_slabs=2` starts working and
  `max_n_tpts_per_slab <= 2` starts raising; nothing in the tests, the notebooks or
  the docs passed a value that newly raises.

  *The docs gate had been red since `2debd51`, and the branch did not know it.*  A
  `:func:` role in `implementation_details.rst` pointed at `magnus.magnus._expm_stack`,
  a private name autoapi does not document, which fails `-n -W`.  It went unnoticed
  because every check of that gate on this branch was an **incremental** build, and
  Sphinx does not re-read a file it believes unchanged — so the warning could not
  reappear once its file had been cached.  CI builds from a clean tree and would have
  caught it on the first push.  Now double backticks, which is what the same file does
  for `_osc_prob_scan_constant_h` and for `_expm_stack` itself nine lines later.  The
  gate is re-verified from `make clean`: 0 warnings.

  Examined and left alone: the Cayley-Hamilton algebra (`det X`, the divided
  differences, the root ordering and `Z²` re-derived independently against the
  code); the d = 2 kernel, which has no `SEV_TOL` gate but beats `eigh` two orders
  of magnitude beyond its documented range (1.9e-09 against 3.3e-09 at
  ‖K‖ = 1e7); `float(VCC_func)` on an array-valued density, unreachable because
  `validate_input_battery` rejects it first; and `verbose=None`, which already
  raised on `main`.  The constant engine was re-checked against `scipy.linalg.expm`
  rather than against another magnus engine — the two in-package routes share
  `_expm_stack` and agree to exactly 0.0, which is no evidence at all — and matches
  to 3.4e-15 across ν/ν̄, three densities, per-point baselines, L0 ≠ 0, a single
  scalar point and vacuum.

- **A max-effort code review of this branch found fifteen defects; all are fixed
  here or in the two commits below.**  Eight independent finder angles, every
  finding confirmed by execution rather than reading.  The two that mattered most
  were invisible to the tests that were supposed to catch them.

  *The compiled exponential was up to 7440x less accurate than `eigh` where a
  clustered spectrum meets a large norm* — 2.7e-07 against 3.0e-11 at
  ||K|| = 1e5, because `arccos` has infinite derivative at u = ±1.  It had been
  verified against random spectra at many norms, and separately at many
  eigenvalue separations at norm ~1; the damage needs both at once, which no
  single-axis sweep visits, and this file previously claimed the kernel was
  "the same order or slightly better at every norm" on that evidence.
  `expmkernels.SEV_TOL` now hands such matrices to `eigh`: worst absolute error
  over the whole separation-by-scale grid 8.7e-14, matrices declined on real work
  0.00%, speed unchanged.  The grid is now a test.

  *Two tests were vacuous.*  The antineutrino-sign test — written specifically to
  catch the `h_matt` bug in the commit below — passed identically with the engine
  under test disabled.  The agreement test compared the new engine against
  `osc_prob`'s *other* constant shortcut rather than the refinement ladder, so its
  tolerance could never fire; instrumented, the slab machinery ran zero times on
  either side.  Both now assert the route they claim to compare, and both were
  verified to fail when their target bug is reintroduced.

  Also: `expm_herm_stack` ignored `supports_dim` and handed 4x4 input to the 3x3
  kernel (error 2.4, unitarity 11.3, uninitialized eigenvalues) and segfaulted at
  d=1; the constant engine answered `L < L0` with the *transpose* of the right
  answer (29% off, row sums exactly 1) where `main` raised, accepted
  `magnus_exp_order=0`, and bypassed the output-size memory guard; `verbose=1`
  lost all its output; array-valued and 0-d-array parameters began raising in two
  new caches; and a cache hit silenced `DensityUnitWarning` after the first call.

  Four further findings inherited from earlier commits on this branch are recorded
  in `docs/dev/HANDOVER_OVERHEAD.md` and deliberately left for separate work.


- **`h_matt` meant two different things depending on the potential, and the new
  engine walked into it.**  `osc_prob_matter_nsi` and `osc_prob_liv` rebound
  `h_matt` to `VCC_func*h_matt` on their constant-potential branch, then passed
  that name to dispatchers documented to take `h_matt` as "the constant matrix
  multiplying `VCC_func(l)`" — which multiply by `VCC` themselves.  The
  separable engine never noticed, because the rebinding only happens on the
  constant branch it used to decline outright.

  The failure was invisible in the two ways that matter: `VCC²` is ~1e-25 rather
  than ~1e-13, so the matter term all but vanished and the answer stayed a
  plausible, unitary, nearly-vacuum probability; and `VCC²` has no sign, so the
  neutrino and antineutrino results came back **bit-identical**.  Standard
  constant-density matter was unaffected (that call site passes the bare
  projector), so testing the headline case alone would have missed it.  The
  scaled matrix now has its own name, and `tests/test_engines.py` compares the
  two routes across every scenario wrapper and both signs — verified to fail
  when the bug is reintroduced.

- **A new engine was invisible to the cross-check.**  `_CROSS_CHECK_FORCING`'s
  forbid lists and `ENGINE_FAMILIES` did not know about `'constant'`, and since
  it answers before `osc_prob_energy_baseline` — which is what records the
  payload the independent `expm` reference is built from — enabling it silently
  removed the only non-Magnus oracle in the table.  It is now listed in every
  other row's forbid set, and shares the `'exact'` family with `expm` rather
  than standing alone: the two use different exponential implementations but
  share the *assumption* that H is position-independent, and that assumption is
  the thing that could be wrong.


- **`_expm_stack`'s docstring claimed the `eigh` route was exactly unitary, and
  it is not.**  `U†U - I` measures 4e-16 for a single 3×3 and 4e-15 for a stack
  of 4096 — growing with stack size, never zero, because reconstruction from
  eigenvectors rounds like any other floating-point product.  The claim that
  probabilities "sum to 1 by construction" was the part worth correcting: they
  sum to 1 to about 1e-15, which is worth relying on, by rounding rather than by
  construction.  No behavior changed by this entry.

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
  rather than optimizations: they keep peak memory at `O(block) + O(result)` instead of
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
  bare numbers inside the function. Naming them changes no behavior, and makes
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
  per-flavor parameter dictionaries, flavor indices and expansion orders out
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
  by more than 1e-2. The scans now get their `n_slabs=150` honored, and they
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
- **PREM was sampled at the center of the Earth.** Notebooks 04, 05 and 07
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
