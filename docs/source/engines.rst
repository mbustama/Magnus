Engines and dispatch
======================

.. contents::
   :local:
   :depth: 2

Which engine answers a call, how the choice is made, and why each one
exists.  See :doc:`performance` for what they cost and :doc:`diagnostics`
for what to do when one warns.

.. _the-engines:

The engines
-------------

Seven engines can answer a request, and each declines the ones it cannot serve
honestly. Several of them share machinery, which is what makes the section after the
table necessary.

.. list-table::
   :header-rows: 1
   :widths: 18 30 28 24

   * - Engine
     - What it assumes
     - When it applies
     - When it declines
   * - **General Magnus ladder**
       (:func:`magnus.oscprob.osc_prob`)
     - Nothing beyond a Hermitian ``H(l)``.
     - Always. Every other engine falls back to it.
     - Never -- it is the terminal path.
   * - **Two-flavor interaction picture**
       (``_osc_prob_ip_exp_dispatch``)
     - A genuine exponential profile, built by
       :func:`magnus.matter.exp_density_profile`; exactly two flavors.
     - Single points and multi-energy scans at one baseline.
     - Non-exponential profiles, :math:`d > 2`, LIV, breakpoints, and whenever its own
       iteration fails to converge (typically near an MSW resonance).
   * - **Constant Hamiltonian**
       (``_osc_prob_scan_constant_h``)
     - ``V_CC`` does not depend on position, so neither does ``H``.
     - Vacuum and constant density, at any flavor count, for a single point or a scan,
       with per-point baselines allowed.
     - A position-dependent potential; user slab edges; parallel, logged or verbose runs.
   * - **Energy-batched separable scan**
       (``_osc_prob_scan_separable``)
     - ``H`` separates into an energy-dependent part and ``V_CC(l)`` times a constant
       matrix.
     - Many energies sharing one baseline.
     - Per-point baselines, user slab edges, parallel or logged runs, a constant
       potential (which the constant engine takes instead).
   * - **Cumulative baseline scan**
       (``_osc_prob_cumulative_scan``)
     - Baselines nest: :math:`U(0\to L_2) = U(L_1 \to L_2)\,U(0 \to L_1)`.
     - A baseline scan at a **single** energy, with a position-dependent ``H``.
     - Differing energies, ``t_slab_edges``, a baseline behind ``L0``, a constant ``H``.
   * - **Phase average**
       (:mod:`magnus.avgprob`)
     - The observable averages over the energy resolution; every interference term keeps
       its phase, weighted by the spread of that phase across ``average_spread``.
     - ``average=True``, on every entry point that takes the keyword.
     - Nothing -- but it warns where the result depends on the spread, and where the
       profile has a feature narrower than its 200-probe grid that could move
       probability, it takes the windows from the hybrid's refinement, or warns that
       none resolves it.
   * - **Adiabatic + Magnus hybrid**
       (:func:`magnus.adiabatic.hybrid_propagator`)
     - ``H`` is smooth at the scale of a 200-point probe grid.
     - Any dimension, any smooth position-dependent profile, with a requested tolerance.
     - Breakpoints or slab edges supplied, a constant potential, no requested tolerance,
       a profile that fails the resolution test, or failure to self-certify.

The eighth entry in the registry, ``scipy.linalg.expm``, never answers a request; it is
used as an oracle by
:func:`magnus.oscprob.cross_check_strategies` wherever it is *exact* -- a constant ``H``,
or a piecewise-constant one whose edges are declared.

Independence, and why it matters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The engines are **not** all independent of each other, and pretending otherwise would make
any cross-check between them worthless. :data:`magnus.oscprob.ENGINE_FAMILIES` records the
grouping the package will defend:

* ``'magnus-ladder'`` -- the general path, the cumulative scan and the separable scan. All
  three walk slabs with :func:`magnus.magnus.magnus_expansion_multislab`, and the cumulative
  scan additionally *sizes* its grid from an ordinary adaptive :func:`magnus.oscprob.osc_prob`
  probe, so it inherits that path's stopping rule as well.
* ``'interaction-picture'`` -- the two-flavor fast path. Same Magnus core, but the fast
  vacuum phase is factored out analytically first, so what it must resolve is a different
  function.
* ``'adiabatic'`` -- the hybrid strategy. A genuinely different method; its blind spots are
  the resonance detector's, not the quadrature's.
* ``'phase-average'`` -- the phase average. It propagates only across non-adiabatic
  windows, with the hybrid's Magnus patch, and carries every stretch between them
  analytically, so it shares no quadrature with the ladder; what it shares with the others
  is the eigendecomposition of the same ``H``.
* ``'exact'`` -- ``expm`` and the constant-Hamiltonian engine, independent of the rest.

Two engines in the same family can be wrong in the same way at the same time. Their
disagreement is informative; their agreement is not.


.. _dispatch-order:

Dispatch
----------

Every scenario wrapper (:func:`magnus.oscprob.osc_prob_matter_std_potential`,
:func:`magnus.oscprob.osc_prob_matter_nsi`, :func:`magnus.oscprob.osc_prob_liv`) tries the
engines in a fixed order, falling through on ``NotImplemented``:

.. list-table::
   :header-rows: 1
   :widths: 40 26 34

   * - Taken when
     - Engine
     - Why it is first
   * - ``average=True`` and the Hamiltonian is position-independent
       (on every entry point, ``osc_prob_energy_baseline`` and the
       Earth and Sun routes included)
     - closed-form phase average (``magnus.avgprob``)
     - No propagation at all; the phase average of a constant ``H`` is algebraic
   * - Smooth profile, a tolerance was requested, ``strategy != 'magnus'``,
       the scan is shorter than
       :data:`~magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`,
       and ``'auto'`` did not hand it to the ladder (below) --
       *and it certifies*
     - adiabatic + Magnus patch
     - Transports along the levels instead of resolving every oscillation
   * - Exponential profile, two flavors, not handed to the ladder --
       *and it converges*
     - interaction picture
     - An exact reference solution exists for this one case
   * - ``V_CC`` does not vary along the trajectory
     - constant Hamiltonian
     - The series terminates at its first term, so the answer is one exponential
   * - Many energies at a single baseline
     - energy-batched scan
     - One traversal serves every energy
   * - Baseline scan at one energy, at least
       :data:`~magnus.oscprob.CUMULATIVE_AUTO_MIN_POINTS` points
     - cumulative scan
     - Every baseline is a prefix of the longest one
   * - Otherwise, or whenever ``return_evolution_operator=True``
     - general Magnus ladder
     - Always applicable; the one that is never skipped, and the only one
       that forms the evolution operator

Each row falls through to the next on ``NotImplemented``, so the last row is
reached whenever nothing above it applies.


Two thresholds decide the seams, and both are constants with docstrings of their own:

* :data:`magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` = 8. Under
  ``strategy='auto'`` the hybrid strategy stands aside for a baseline scan of at least this
  many points, because the cumulative scan answers all of them from one traversal.  This was
  25; the constant's own docstring records why it moved, and why a later attempt to lower it
  to 1 was reverted.
* :data:`magnus.oscprob.CUMULATIVE_AUTO_MIN_POINTS` = 2. Below this there is no prefix to
  reuse.

**The ladder route of** ``'auto'`` (issue #70).  On a smooth profile the hybrid strategy's cost
is its window search, which does not follow the tolerance.  At a loose tolerance on a moderate
phase that makes it the slower route, so ``'auto'`` hands a request to the ladder ahead of the
hybrid when all three of these hold:

* ``min(rtol, atol)`` is at least :data:`magnus.oscprob.AUTO_LADDER_MIN_TOLERANCE` = 1e-6;
* the estimated accumulated phase (the integral of the spread of ``H``'s eigenvalues up to the
  longest baseline) is at most :data:`magnus.oscprob.AUTO_LADDER_MAX_PHASE` = 1e4 rad;
* the ladder's starting slab count is at most
  :data:`magnus.oscprob.AUTO_LADDER_MAX_FLOOR_FRACTION` = 1/4 of its cap, which keeps every
  solar path on the hybrid.

The ladder then runs at a tenth of the tolerance
(:data:`magnus.oscprob.AUTO_LADDER_TOLERANCE_MARGIN`), skips the interaction picture and starts
on slabs over which the Magnus series is guaranteed to converge.  The hybrid's test for an
undeclared density jump still runs, and still warns.  Over 20 smooth workloads with phases from
5 to 1.2e4 rad, the ladder was 2 to 60 times faster than the hybrid on a single point and 12 to
500 times faster per point of a 40-energy scan, within the tolerance on every one.  On the
profile of the paper's Fig. 1, its four scans of 140 energies take 40 ms of computation at the
default tolerance of 1e-3, where the hybrid took 8 s.

**The accuracy steps at the seam rather than varying smoothly, and that is by design.**
Adding one baseline to a scan just below it changes the answer, because it changes the engine.
Measured against ``solve_ivp``:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - Profile
     - err(N = 24)
     - err(N = 26)
     - Ratio
   * - solar exponential
     - 3.30e-05
     - 2.13e-08
     - 1 546×
   * - multi-resonance
     - 1.58e-03
     - 2.86e-09
     - 552 945×
   * - noisy
     - 6.27e-04
     - 1.04e-08
     - 60 418×
   * - castle wall + breakpoints
     - 2.80e-11
     - 2.80e-11
     - 1.0× (cumulative from N = 2)

The jump is always *toward* the truth, so this is a documentation matter rather than a
numerical one -- but a user scanning N and watching their answer move by five orders of
magnitude in accuracy will otherwise assume something is broken.

**Seeing which engine answered.** The fallbacks are silent by design: they happen on
ordinary calls and warning about them would be noise. Pass ``strategy_info`` to any of the
three scenario wrappers to see the route without changing it::

    info = {}
    P = oscprob.osc_prob_matter_std_potential(..., strategy_info=info)
    info['engine']      # 'hybrid', 'ip_exp', 'separable', 'constant',
                        # 'cumulative', 'magnus' or 'average'
    info['certified']   # for the hybrid strategy
    info['declined']    # [(engine, why it stood aside)]

This is the answer to "why did my result move?" and "why did this call get slow?", both of
which were previously unanswerable from outside the package.
