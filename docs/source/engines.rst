Engines and dispatch
====================

.. contents::
   :local:
   :depth: 2

This page explains which engine answers a call, how the choice is made, and why each one
exists.  See :doc:`performance` for what they cost and :doc:`diagnostics` for what to do
when one warns.

.. _the-engines:

The engines
-----------

.. figure:: ../../img/paper/strategies.svg
   :width: 100%
   :alt: The seven probability engines, in the order they are tried

   The seven engines, in the order they are tried.  A request that does not meet an
   engine's conditions (right) passes to the next, so the general Magnus ladder answers
   whatever no other engine takes.  Each row sketches what its engine does along the
   trajectory; shading is the matter density.  From the Magνs paper.

Seven engines can answer a request, from the most specialized to the most general, and the
first whose conditions the request meets answers it.  Several engines share machinery, which
is why `Independence, and why it matters`_ follows the table.

.. list-table::
   :header-rows: 1
   :widths: 4 18 30 24 24

   * -
     - Engine (name in ``strategy_info``)
     - What it does
     - When it applies
     - When it declines
   * - 1
     - **Averaged probability** (``'average'``; :mod:`magnus.avgprob`)
     - The phase average, by one of three routes: in closed form from one eigenbasis when
       ``H`` does not depend on position; transported along the instantaneous eigenstates,
       with a Magnus patch across each non-adiabatic crossing, when ``H`` varies smoothly;
       the mean over an energy window across declared discontinuities
       (:doc:`averaged_probability`).
     - ``average=True``, on every entry point that takes the keyword.
     - Never; it warns where the result depends on the spread, and where the profile has a
       feature narrower than its 200-probe grid that no refinement resolves.
   * - 2
     - **Adiabatic + Magnus** (``'hybrid'``; :func:`magnus.adiabatic.hybrid_propagator`)
     - Transport along the instantaneous eigenstates, with a Magnus patch across each
       non-adiabatic window (:doc:`adiabatic_strategy`).
     - Any smooth position-dependent ``H``, any number of flavors, with a requested tolerance.
     - Breakpoints or slab edges supplied, no requested tolerance, a profile it does not
       resolve on its probe grid; under ``'auto'``, a failed certification and the requests
       the rules below send elsewhere.
   * - 3
     - **Interaction picture** (``'ip_exp'``)
     - Factors out the vacuum phase (and the LIV term, if present) analytically.  Over
       each slab, it integrates the exponential matter envelope in closed form and keeps
       the first-order Magnus term.
     - Two flavors, a profile built by :func:`magnus.matter.exp_density_profile`, one
       baseline (a single point or an energy scan).
     - More than two flavors, any other profile (a tabulated solar model included),
       breakpoints or slab edges, and a failure to converge, as near an MSW resonance.
   * - 4
     - **Constant Hamiltonian** (``'constant'``)
     - One exponential, exact: the Magnus series terminates at its first term.
     - ``H`` does not vary along the path; any number of flavors, a single point or a scan,
       with per-point baselines allowed.
     - A position-dependent potential; user slab edges; parallel, logged or verbose runs; a
       refinement value ``osc_prob`` would reject (other refinement keywords are ignored).
   * - 5
     - **Energy-batched scan** (``'separable'``)
     - One set of slabs shared by every energy: the potential is sampled once per refinement
       level, and the energy is a batch dimension.
     - Many energies at one baseline, with ``H`` = energy-dependent part +
       :math:`V_{\rm CC}(l)` times a fixed matrix.
     - Per-point baselines, user slab edges, parallel or logged runs, a constant potential
       (engine 4 takes it).
   * - 6
     - **Cumulative scan** (``'cumulative'``)
     - One pass along the longest baseline, recording the running product at every requested
       baseline: :math:`\mathbb{U}(0\to L_2) = \mathbb{U}(L_1 \to L_2)\,\mathbb{U}(0 \to L_1)`.
     - Many baselines at one energy, with a position-dependent ``H``.
     - Differing energies, ``t_slab_edges``, a baseline behind ``L0``, a constant ``H``.
   * - 7
     - **General Magnus ladder** (``'magnus'``; :func:`magnus.oscprob.osc_prob`)
     - Slabs refined until two successive levels agree (:doc:`methodology`).
     - Always.
     - Never; it is the last engine.

Not every entry point tries every engine.  The matter scenario functions try the first five,
in order; ``osc_prob_energy_baseline`` tries the cumulative scan; ``osc_prob`` is the ladder
itself.  ``osc_prob_vacuum`` needs only the averaged probability and the constant engine.  A
Hamiltonian of your own, through ``osc_prob_earth`` or ``osc_prob_sun``, can reach the
averaged probability, the adiabatic engine, the cumulative scan and the ladder.  A request
for the evolution operator (``return_evolution_operator=True``) goes to the ladder, the only
engine that forms it.

The eighth entry in the registry, ``scipy.linalg.expm``, never answers a request; it is
used as an oracle by
:func:`magnus.oscprob.cross_check_strategies` wherever it is *exact* — a constant ``H``,
or a piecewise-constant one whose edges are declared.

Independence, and why it matters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The engines are **not** all independent of each other, and a cross-check between two
engines that share their machinery proves little.  :data:`magnus.oscprob.ENGINE_FAMILIES`
groups them into families:

* ``'magnus-ladder'`` — the general path, the cumulative scan and the separable scan. All
  three walk slabs with :func:`magnus.magnus.magnus_expansion_multislab`, and the cumulative
  scan additionally *sizes* its grid from an ordinary adaptive :func:`magnus.oscprob.osc_prob`
  probe, so it inherits that path's stopping rule as well.
* ``'interaction-picture'`` — the two-flavor fast path. It uses the same Magnus core, but the fast
  vacuum phase is factored out analytically first, so what it must resolve is a different
  function.
* ``'adiabatic'`` — the hybrid strategy. It is a different method, whose blind spots are
  the resonance detector's, not the quadrature's.
* ``'phase-average'`` — the phase average. It propagates only across non-adiabatic
  windows, with the hybrid's Magnus patch, and carries every stretch between them
  analytically.  The patch uses the ladder's Magnus slab kernel
  (:func:`magnus.magnus.magnus_expansion_multislab`), but only inside those windows; what the
  phase average shares with the others is the eigendecomposition of the same ``H``.
* ``'exact'`` — ``expm`` and the constant-Hamiltonian engine, independent of the rest.

Two engines in the same family can be wrong in the same way at the same time. Their
disagreement is informative; their agreement is not.


.. _dispatch-order:

Dispatch
--------

Each matter scenario function (:func:`magnus.oscprob.osc_prob_matter_std_potential`,
:func:`magnus.oscprob.osc_prob_matter_nsi`, :func:`magnus.oscprob.osc_prob_liv`) tries the
engines in a fixed order, falling through on ``NotImplemented``:

.. list-table::
   :header-rows: 1
   :widths: 40 26 34

   * - Taken when
     - Engine
     - Why it is first
   * - ``average=True`` (on every entry point that takes it, ``osc_prob_energy_baseline``
       and the Earth and Sun routes included)
     - averaged probability, by one of its three routes
     - It answers a different question from the other engines: the average, not the
       probability at one energy
   * - Smooth profile, a tolerance was requested, ``strategy != 'magnus'``,
       the scan is shorter than
       :data:`~magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` (than
       :data:`~magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS_TIGHT` below 1e-6),
       and ``'auto'`` did not hand it to the ladder (below) —
       *and it certifies*
     - adiabatic + Magnus patch
     - Transports along the levels instead of resolving every oscillation
   * - Exponential profile, two flavors, not handed to the ladder —
       *and it converges*
     - interaction picture
     - The vacuum phase is factored out, so the oscillation is never resolved slab by slab
   * - :math:`V_{\rm CC}` does not vary along the trajectory
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


Three thresholds decide where one engine hands a request to another.  Each is a constant
whose docstring records the measurements behind its value:

* :data:`magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` = 8. Under
  ``strategy='auto'`` the hybrid strategy stands aside for a baseline scan of at least this
  many points, because the cumulative scan answers all of them from one traversal.  The
  cumulative scan is the cheaper of the two at the median at every size measured.  A lower
  threshold would not help: below 8 points, the hybrid's request does not always pass to the
  cumulative scan, but to whichever engine applies next.
* :data:`magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS_TIGHT` = 2. The same threshold
  below a tolerance of 1e-6, where the hybrid strategy is the slower route for a baseline scan.
* :data:`magnus.oscprob.CUMULATIVE_AUTO_MIN_POINTS` = 2. Below this there is no prefix to
  reuse.

.. _auto-ladder-route:

**The ladder route of** ``'auto'``.  On a smooth profile, the hybrid strategy's cost is its
window search, which does not follow the tolerance, so at a loose tolerance on a moderate phase
it is the slower route, by one to two orders of magnitude.  Therefore, at a tolerance
``min(rtol, atol)`` of :data:`magnus.oscprob.AUTO_LADDER_MIN_TOLERANCE` = 1e-6 or looser,
``'auto'`` hands a request to the ladder ahead of the hybrid when both of these hold:

* the estimated accumulated phase (the integral of the spread of ``H``'s eigenvalues up to the
  longest baseline) is at most :data:`magnus.oscprob.AUTO_LADDER_MAX_PHASE` = 1e4 rad;
* the ladder's starting slab count is at most
  :data:`magnus.oscprob.AUTO_LADDER_MAX_FLOOR_FRACTION` = 1/4 of its cap, which keeps every
  solar path on the hybrid.

For an energy scan at one baseline that the energy-batched scan will take, the phase condition
is dropped: that engine shares its slabs across the energies, so a limit that prices the ladder
one point at a time does not apply.

The ladder then runs at a tenth of the tolerance
(:data:`magnus.oscprob.AUTO_LADDER_TOLERANCE_MARGIN`), skips the interaction picture and starts
on slabs over which the Magnus series is guaranteed to converge.  The hybrid's test for an
undeclared density jump still runs, and still warns.

**At a tighter tolerance**, the route stays open on ``integration_method='gl'`` at a single
baseline, with a phase limit that shrinks with the tolerance and the order.  The limit is
``AUTO_LADDER_MAX_PHASE*(tol/1e-6)**(1/p)``, with ``p`` the requested ``magnus_exp_order``,
capped at :data:`magnus.oscprob.AUTO_LADDER_TIGHT_MAX_PHASE` = 2000 rad so that partial solar
chords stay on the hybrid.  It shrinks because the ladder's slab count grows as
``tol**(-1/p)`` and the hybrid's cost does not.
The limit applies to an energy scan as well.  At these tolerances, the ladder runs at the
requested tolerance itself, because its levels are deep in the asymptotic regime, where the
difference between two of them already overestimates the error of the finer one.  For
example, at ``rtol = 1e-12``, ``atol = 1e-14`` and ``magnus_exp_order = 8``, the limit is
1000 rad, and the four curves of the validation example in :doc:`diagnostics` have phases of
10–78 rad, so they take this route.  A baseline scan goes to the cumulative scan at such
tolerances.

**Which engine answers a request.**  The table below combines these rules.  For each request
shape and tolerance ``min(rtol, atol)``, it gives the engine that answers a scenario function or
a wrapper under ``strategy='auto'``.  A Hamiltonian of your own
reaches fewer engines; see `The engines`_.  "Many energies" are at one
baseline, and "baselines" are at one energy.  "Too many slabs" means that the ladder would
start with more than a quarter of its slab cap, as across the Sun.  "The tightened limit" is
the phase limit of the paragraph above, on ``integration_method='gl'``; other quadratures keep
the adiabatic engine first below 1e-6.  ``strategy_info`` reports the engine that answered.

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Request
     - Tolerance ≥ 1e-6 (includes the default, 1e-3)
     - Tolerance < 1e-6
   * - ``average=True``
     - averaged probability
     - averaged probability
   * - Constant ``H``
     - constant Hamiltonian
     - constant Hamiltonian
   * - Smooth ``H``, one point
     - general ladder; adiabatic if the phase exceeds 1e4 or with too many slabs
     - general ladder if the phase is within the tightened limit; otherwise, or with too many
       slabs, adiabatic, and if that cannot certify, the ladder (the interaction picture for
       an exponential profile at two flavors)
   * - Smooth ``H``, many energies
     - energy-batched scan; adiabatic with too many slabs
     - energy-batched scan if the phase is within the tightened limit; otherwise, or with too
       many slabs, adiabatic, and if that cannot certify, the energy-batched scan
   * - Smooth ``H``, 2–7 baselines
     - cumulative scan; adiabatic if the phase exceeds 1e4 or with too many slabs
     - cumulative scan
   * - Smooth ``H``, 8 or more baselines
     - cumulative scan
     - cumulative scan
   * - Declared discontinuities (``t_breakpoints`` or ``t_slab_edges``)
     - ladder, energy-batched scan, or cumulative scan, by the shape of the request
     - the same

.. _strategy-magnus:

Passing ``strategy='magnus'`` uses only the Magnus engines: it keeps an energy scan on the
energy-batched scan whatever its phase, and it turns off the cumulative scan, so a baseline
scan is computed point by point.  This forces every result through the Magnus quadrature, but
it is not the fast choice for a baseline scan.

**The accuracy changes in a step at each threshold.**  Adding one baseline to a scan that sits
immediately below a threshold (from 7 to 8 baselines at 1e-6 and looser, from 1 to 2 below)
changes the answer, because it changes the engine.  In the cases measured, the step was toward
the more accurate answer, by up to six orders of magnitude; the docstring of
:data:`~magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` has the measurement.

**Seeing which engine answered.**  Falling through from one engine to the next raises no
warning, because it happens on ordinary calls.  Pass ``strategy_info`` to any scenario
function or wrapper, or to :func:`~magnus.oscprob.osc_prob_earth` or
:func:`~magnus.oscprob.osc_prob_sun`, to see the route without changing it::

    info = {}
    P = oscprob.osc_prob_matter_std_potential(..., strategy_info=info)
    info['engine']      # 'hybrid', 'ip_exp', 'separable', 'constant',
                        # 'cumulative', 'magnus' or 'average'
    info['certified']   # for the hybrid strategy
    info['declined']    # [(engine, why it gave up)], for engines that tried

A result that moves by more than the tolerance when a point is added to a scan, or a call
that suddenly costs more, usually reflects a change of engine, and the dictionary names it.
The one fallback that warns, once per session, is the adiabatic engine declining a profile
with a jump or with a feature too narrow for its grid.
