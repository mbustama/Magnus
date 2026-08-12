# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""adiabatic.py

Contains the adiabatic-transport-plus-Magnus-patch ("hybrid") propagator
used as an alternative strategy to compute the neutrino evolution
operator when the Hamiltonian is position-dependent and develops an
extreme accumulated phase (e.g., low-energy solar neutrinos crossing an
MSW resonance), the regime in which the plain slab-refinement Magnus
engine in :mod:`magnus.oscprob` needs a very large number of slabs (and
may raise ``ToleranceNotAchievedWarning``).

Physical idea: away from an eigenvalue crossing (or narrowly-avoided
crossing) of the instantaneous Hamiltonian :math:`H(l)`, the adiabatic
theorem says the evolution operator is well approximated by transport in
the *instantaneous eigenbasis* of :math:`H(l)` -- a dynamical phase (the
integral of the eigenvalues) plus a geometric (Berry) phase, both cheap
to compute on a coarse grid regardless of how large the accumulated
phase is. Near a genuine MSW resonance, the adiabatic approximation
breaks down over a narrow window, which is patched with an exact,
short-baseline Magnus computation (:func:`magnus.magnus.magnus_expansion_multislab`,
the package's own, already-unitary integrator). The two pieces are
stitched together with the exact composition law of quantum evolution,
:math:`U(l_2, l_0) = U(l_2, l_1) U(l_1, l_0)`, so the result is exactly
unitary regardless of the approximation's accuracy.

Where a patch is needed is decided by an *exact* Hellmann-Feynman
diagnostic (no finite-differenced eigenvectors, which are gauge-
ambiguous), so this applies to any Hermitian Hamiltonian of any
dimension, with any number of simultaneous or sequential resonances --
see :doc:`/adiabatic_strategy` for the full derivation, validation, and
worked examples.

This module is self-contained: it depends only on :mod:`magnus.magnus`
(the Magnus-expansion core), not on :mod:`magnus.oscprob`, so it can be
used directly on any Hamiltonian function, independently of the rest of
the oscillation-probability API. :mod:`magnus.oscprob` calls
:func:`hybrid_propagator` internally when ``strategy='hybrid'`` or
``strategy='auto'`` (the default) is passed to
:func:`magnus.oscprob.osc_prob_matter_std_potential`,
:func:`magnus.oscprob.osc_prob_matter_nsi`, or
:func:`magnus.oscprob.osc_prob_liv` (and, transitively, to every
``osc_prob_*_sun``/``osc_prob_*_sun_nsi``/``osc_prob_*_sun_liv`` wrapper),
and also when it is passed to the fully generic user-Hamiltonian entry
points :func:`magnus.oscprob.osc_prob_sun` and
:func:`magnus.oscprob.osc_prob_earth` (via
``magnus.oscprob._osc_prob_with_potential``).  For ``osc_prob_earth``
the hybrid path is normally declined, since a real Earth trajectory
supplies PREM layer breakpoints; see :doc:`/adiabatic_strategy`.

Routine listings
----------------

    * adiabatic_propagator - Evolution operator via pure adiabatic
           (instantaneous-eigenbasis) transport, no resonance patching
    * find_hidden_features - Detects structure too narrow for any grid
           this package lays down to sample
    * find_resonance_candidates - Locates every exact eigenvalue-gap
           critical point of H(l) via the Hellmann-Feynman theorem
    * find_nonadiabatic_windows - Filters/grows/merges candidates into
           position windows that need a Magnus patch
    * hybrid_propagator - Adiabatic transport with Magnus patches at
           non-adiabatic windows, self-certified against successive
           refinement of every internal tolerance knob
    * oscillation_sampling - Reports how coarsely a request samples the
           oscillation it is computing (cycles per step, Nyquist points)
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import simpson

import magnus.magnus as magnuscore


GAMMA_TO_ERROR = 0.85
r"""float: Module-level constant

Probability error to budget per unit of the adiabaticity parameter :math:`\gamma`, when
:func:`hybrid_propagator` decides whether a result with **no** non-adiabatic window may be
certified.  The rule is ``GAMMA_TO_ERROR * gamma_max <= atol + rtol``.

Write the pure adiabatic answer's error as :math:`|\Delta P| = k\,\gamma_\max`.  Certifying an
empty window list is safe exactly when :math:`k\,\gamma_\max \le` tolerance, so this constant
is an upper bound on :math:`k` -- and the only honest way to set it is to measure :math:`k`.

Measured over 149 configurations (``docs/dev/adversarial_batteries/gamma_slack_sweep.py``):
resonance width swept over a decade, d = 2-5, 5-80 MeV, 0.5-2 density scale heights, scoring the
pure adiabatic operator against ``solve_ivp``/DOP853.  **Only the small-**:math:`\gamma` **rows
matter**, because this rule governs the no-window case alone: once :math:`\gamma_\max` exceeds
the active threshold a window opens and the answer is patched exactly, whatever :math:`k` is
there.

========================== ======= ============ ============================
population                    n     max k        implied bound on gamma_max
========================== ======= ============ ============================
gamma_max < 1e-2              76    **0.812**    <= 1.23 x tolerance
gamma_max < 3e-3              35    0.679        <= 1.47 x tolerance
gamma_max < 1e-3              12    0.502        <= 1.99 x tolerance
all rows, incl. patched      149    1.136        <= 0.88 x tolerance
========================== ======= ============ ============================

:math:`k` falls towards ~0.5 as :math:`\gamma` shrinks, which is what the linear model predicts
asymptotically.  0.85 covers the worst case in the governed regime (0.812) with a little margin.

**This constant was wrong twice, in opposite directions.**  It began at 1.0 alongside a slack
factor of 2.0, derived from five points that all happened to sit at :math:`\gamma_\max <
10^{-3}` where :math:`k \approx 0.5`; that pair encoded :math:`k \le 0.5`, right for those five
and optimistic by up to 1.6x elsewhere in the governed regime.  Reading the *unrestricted*
maximum (1.136) then argued for a far stricter bound -- an over-correction, since those rows sit
at :math:`\gamma_\max \sim 0.2`, open a window immediately, and are never decided by this rule.
The slack factor is gone: it existed only to compensate for the mis-measured value, and with
:math:`k` measured in the regime that matters no fudge is needed.

.. versionadded:: 1.0.0
"""


RESOLUTION_RATIO = 0.70
r"""float: Module-level constant

Threshold of the probe-scale resolution test in ``_profile_is_resolved``, which decides
whether ``H_func`` is sampled finely enough for this module's finite-difference diagnostics
to mean anything.

Within each probe interval, the test asks **what fraction of the variation falls in one
half**.  For a :math:`C^1` Hamiltonian the two halves each carry about half, so the ratio tends
to **0.5**; a jump lands entirely inside one half, so the ratio tends to **1.0**.

The two limits are 0.5 and 1.0, but the honest threshold is set by measurement rather than by
their midpoint, because a jump *comparable in size to the local smooth variation* lands between
them: a jump :math:`J` on top of a smooth change :math:`S` within one interval gives
:math:`(J + S/2)/(J + S)`.  Measured:

=================================================== ==========
population                                          statistic
=================================================== ==========
192 smooth configurations (6 profile families plus
10 random Fourier sums, d = 2-5, 5-200 MeV)         **<= 0.602**
a jump 4.7x smaller than the steepest smooth step   0.773
15 random piecewise-constant profiles, d = 2 and 3  **1.000**
=================================================== ==========

0.70 sits in the gap with margin on both sides -- 16% above the smooth ceiling, 10% below the
weakest genuine discontinuity.  It also states what the test can and cannot catch: solving
:math:`(J + S/2)/(J + S) > t` gives :math:`J/S > (2t-1)/(1-t)`, so at 0.70 a jump must be at
least **1.33x the local smooth variation** to be seen.  A smaller one is genuinely
indistinguishable from steep smooth behaviour at that sampling density.

Two earlier formulations were wrong, both caught by measurement.  Comparing the *global* largest
adjacent change at two grid densities masks any jump smaller than the largest smooth variation
elsewhere on the path: a discontinuity 4.7x smaller than the steepest smooth step went
undetected, and the answer came back wrong by 2.0e-02, silently.  Comparing one half against the
*whole interval* fixes that but false-positives at a smooth turning point, where the interval's
net change is near zero while each half is not.  Comparing each half against the sum of the two
is immune to both: it measures concentration of variation, which is what a jump is.

.. versionadded:: 1.0.0
"""


LOCAL_JUMP_RATIO = 0.5
r"""float: Module-level constant

Threshold of the **local confirmation** in ``_profile_is_resolved``: having flagged a probe
interval whose variation is concentrated in one half, re-sample that interval alone on
``N_LOCAL_CONFIRM`` points and ask what fraction of its variation still falls in a single
adjacent step.

The two limits are far apart and are set by arithmetic, not by taste.  A jump is not diluted by
refinement -- one step still carries all of it -- so the fraction tends to **1.0**.  A smooth
feature is spread over the whole interval, so with 32 sub-steps the largest one carries roughly
:math:`1/32` of the variation.  Measured (``docs/dev/adversarial_batteries/resolution_fp.py``),
over the *flagged* intervals only, which is the population this constant decides:

=================================================== ======= ==========
population                                            n      statistic
=================================================== ======= ==========
10 smooth families x d = 2-5 x 3 energies x 12
sub-intervals (1440 configurations)                    79    **<= 0.087**
3 piecewise-constant families, same sweep              348   **>= 1.000**
=================================================== ======= ==========

0.5 sits between them with a factor of 5.7 of margin on the smooth side and 2.0 on the jump
side.  Over the same sweep the completed test reports **0 of 1440** smooth configurations as
unresolved, and every sub-interval that genuinely contains a jump as unresolved.

**Why this exists.**  Without it the concentration test is a maximum over intervals, and one
interval decides -- so the interval containing a smooth *turning point* decides.  There the two
halves are genuinely asymmetric: one nearly cancels while the other does not, and the ratio is a
draw in :math:`[0.5, 1]` that depends on where the extremum happens to fall inside its interval.
Refining the grid does not remove it; it re-draws it, which is why the two-stage
coarse-then-fine protocol could not separate the two either.  Measured on a Gaussian bump of
width :math:`10^{-2}` of the trajectory -- a profile the module answers to 3.5e-09 -- the
statistic hit 0.75 at the interval containing the peak, and 6 of 30 baselines of one ordinary
scan were declared discontinuous.  The local confirmation is the discriminator the concentration
ratio alone does not have: it asks whether the concentration *survives refinement*, which is the
one thing a jump does and a turning point does not.

Costs nothing on an ordinary call: it runs only on intervals the cheap test already flagged, and
on a profile that flags none it is never reached.

.. versionadded:: 1.0.0
"""


N_LOCAL_CONFIRM = 33
r"""int: Module-level constant

Points used to re-sample one flagged probe interval in ``_profile_is_resolved`` (see
:data:`LOCAL_JUMP_RATIO`).  32 sub-steps puts the smooth limit at :math:`1/32 \approx 0.03`,
comfortably below :data:`LOCAL_JUMP_RATIO`, and resolves a feature down to
:math:`1/(199 \times 32) \approx 1.6\times10^{-4}` of the trajectory -- the same scale as the
probe refinement ceiling ``max_n_probe = 6400``, so the local test does not claim to see
anything the caller's own refinement could not.

.. versionadded:: 1.0.0
"""


MAX_LOCAL_CONFIRMATIONS = 8
r"""int: Module-level constant

At most this many flagged intervals are confirmed locally, taken in decreasing order of how
much variation they carry.  A profile with hundreds of genuine jumps is found by the first one,
and a bound is needed because the flagged set is unbounded in principle.

.. versionadded:: 1.0.0
"""


HIDDEN_FEATURE_CONCENTRATION = 0.3
r"""float: Module-level constant

Threshold of :func:`find_hidden_features`, which looks for structure **below the scale any grid
in this package samples on**.  This is the one exposure the adversarial validation could not
close: a Gaussian narrower than the probe spacing is invisible to the hybrid detector, to the
general ladder's slab grid and to the cumulative scan alike, so all three agree and all three
are wrong -- silently, by up to 2.9e-02 against a requested 1e-3.

The statistic is **concentration**, not size, and that choice is the whole design.  Within each
interval of a reference grid, compare the total variation a much denser grid sees inside it with
the change its two endpoints show; the excess is variation hidden between reference samples.
Report the largest such excess as a fraction of the profile's total variation.

A first attempt used the cruder ``TV_dense/TV_reference`` and was wrong: a sinusoid at exactly
the probe spacing sends the denominator to zero and that ratio to :math:`10^{13}`, while the
package answers such a profile to ~1e-11.  What separates the two is *where* the hidden
variation sits.  An aliased sinusoid hides some in **every** interval, so its share of the total
is ~1/n_ref; one narrow bump hides all of it in **one**, so its share is ~1.  The excess is
summed over **adjacent pairs** of reference intervals, so that a feature landing on an interval
boundary is not halved by the split.

Measured over 67 profiles the package serves -- solar, multi-resonance, noisy, sinusoids at 1x,
2x and 1/2x the probe spacing, 400 crossings, a declared step, 30 random Fourier sums, and 30
Gaussian bumps of random width down to 1e-4 of the trajectory at random positions:

======================================================== ==============
population                                                concentration
======================================================== ==============
67 smooth/resolvable profiles                             max **0.060**
                                                          p99 0.047
======================================================== ==============

and, over 60 random positions each, the detection rate for features in the band no grid here
resolves:

============= ================ ================ ================
feature width detection at 0.2 detection at 0.3 detection at 0.5
============= ================ ================ ================
3e-5          0.70             **0.68**         0.55
1e-5          0.90             **0.90**         0.90
3e-6          0.82             **0.82**         0.82
1e-6          0.73             **0.73**         0.73
============= ================ ================ ================

**0.3 gives zero false positives over all 67 smooth profiles** -- five times the measured
ceiling -- at the best detection the margin allows; 0.2 buys two points of detection for half
the margin, and 0.5 costs thirteen.

The caller now varies the sampling density with the size of the request (see
:data:`N_HIDDEN_FEATURE_SUBDIVISION`), so the ceiling was re-measured at **every** density the
dispatcher can choose rather than only at the default: 0.0597 at 8 sub-steps, 0.0601 at 16,
0.0602 at 32 -- **0 of 67 in all three**.  The statistic is a *fraction* of the total variation,
which is why it barely moves: refining the dense grid adds the same variation to numerator and
denominator.

**This detects most of the class, not all of it, and the shortfall is structural.**  A feature
of width 3e-5 is right at the edge of what ``max_n_probe = 6400`` can partially resolve, so its
variation is partly visible to the reference grid and the statistic is diluted; a feature of
width 1e-6 is far below the *dense* spacing
(:math:`(l_1-l_0)/51192 \approx 2\times10^{-5}`), so whether a sample lands inside it is luck.
Between those two limits detection is ~0.8-0.9.  Against a prior state of **zero** detection and
a silent 2.9e-02 error, that is the improvement on offer; it is not a guarantee, and the
docstring says so rather than implying one.

**Why the reference grid is the ceiling and not the starting density.**  Against
``n_probe0 = 200`` the same statistic flags widths of 1e-3 and 1e-4 too -- correctly, in that
they *are* hidden at that density, and uselessly, in that the refinement resolves them.  Keying
it to ``max_n_probe`` asks the question that matters: is there structure left that no amount of
refinement will reach?

.. versionadded:: 1.0.0
"""


N_HIDDEN_FEATURE_SUBDIVISION = 8
r"""int: Module-level constant

Sub-steps per reference interval in :func:`find_hidden_features`, so 51 192 samples of the
profile in total.  Chosen on cost, because this runs on ordinary calls: the statistic is nearly
independent of it (the separation is 4.7x at 2 sub-steps and 4.6x at 32), while the cost is not.

============ ============== ===========
sub-steps     dense samples  scan cost
============ ============== ===========
4             25 597         0.22 ms
**8**         **51 193**     **0.37 ms**
16            102 385        1.30 ms
32            204 769        2.85 ms
============ ============== ===========

The jump past 8 is superlinear -- the arrays stop fitting in cache -- and 0.37 ms is about 3% of
an ordinary 13 ms single-point call, where 2.85 ms would be 20% and fail the package's own 10%
performance criterion.  What it costs in reach is the very narrowest features: the dense spacing
is :math:`(l_1-l_0)/51192`, and detection of anything below that is a matter of whether a sample
lands inside it (measured 0.73 at a width of 1e-6, against 0.90 at 1e-5).

:mod:`magnus.oscprob` scales this with the number of requested points -- 8 below four points,
16 below sixteen, 32 above -- because the scan runs **once per call** whatever the point count,
so its share of the work falls as the request grows.  That holds it under about 7 % of the call
at every size instead of spending 20 % of the cheapest one, and the false-positive rate was
re-measured at each of those three densities (0 of 67 every time; see
:data:`HIDDEN_FEATURE_CONCENTRATION`).

Raise it if you have reason to think the profile hides something finer; the scan is
:func:`find_hidden_features` and takes ``n_sub`` directly.

.. versionadded:: 1.0.0
"""


def _variation_steps(values: np.ndarray) -> np.ndarray:
    r"""Magnitude of the change between consecutive samples, for a scalar or matrix profile.

    .. versionadded:: 1.0.0
    """
    diffs = np.diff(np.asarray(values), axis=0)
    if diffs.ndim == 1:
        return np.abs(diffs)
    return np.max(np.abs(diffs), axis=tuple(range(1, diffs.ndim)))


def find_hidden_features(profile: Callable, l0: float, l1: float,
    n_ref: Optional[int] = 6400, n_sub: Optional[int] = None) -> Dict:
    r"""Looks for structure too narrow for any grid this package lays down to sample.

    See :data:`HIDDEN_FEATURE_CONCENTRATION` for the statistic, why it is a *concentration*
    rather than a size, and the measured separation.

    **Pass the scalar potential when there is one.**  For the separable Hamiltonians this package
    builds, :math:`H(l) = h_\mathrm{vac}/E + V_{CC}(l)\,P_{ee}` is affine in :math:`V_{CC}`, so
    every difference is :math:`|\Delta V_{CC}|` times a constant and the statistic is
    **identical** either way -- verified bit-for-bit at d = 2, 3 and 5.  It is also 18x cheaper
    (2.6 ms against 48 ms), because sampling a scalar avoids allocating a stack of
    :math:`2\times10^5` matrices.  This is a free choice, not an approximation.

    Cost is a single vectorized evaluation of ``profile`` and two array passes; no
    eigendecomposition.  It depends on the profile and the interval but **not** on energy, so a
    caller scanning many energies or baselines should run it once, not once per point.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    profile : Callable
        The position-dependent part of the problem: either the scalar potential (preferred) or
        ``H_func`` itself.  Must accept an array of positions.
    l0, l1 : float
        Interval to scan.
    n_ref : int, optional
        Reference grid, which should be the finest grid the caller's machinery can reach -- the
        hybrid strategy's ``max_n_probe``. Default: 6400.
    n_sub : int, optional
        Sub-steps per reference interval. Default: :data:`N_HIDDEN_FEATURE_SUBDIVISION`.

    Returns
    -------
    dict
        ``'concentration'`` (the statistic), ``'hidden'`` (whether it exceeds
        :data:`HIDDEN_FEATURE_CONCENTRATION`), ``'l_lo'``/``'l_hi'`` (the reference interval
        carrying the excess) and ``'l_centre'``.  On a constant profile, concentration 0.0.
    """
    n_sub = N_HIDDEN_FEATURE_SUBDIVISION if n_sub is None else n_sub
    quiet = {'concentration': 0.0, 'hidden': False,
             'l_lo': float(l0), 'l_hi': float(l1), 'l_centre': 0.5*(float(l0) + float(l1))}
    if (n_ref < 2) or (n_sub < 2) or (l1 == l0):
        return quiet

    dense = np.linspace(float(l0), float(l1), (n_ref - 1)*n_sub + 1)
    try:
        values = np.asarray(profile(dense))
    except Exception:                      # noqa: BLE001 -- any failure means "not vectorized"
        return quiet
    # This is a diagnostic; it must never break a call it was only meant to inspect.  A profile
    # that returns a bare scalar for array input indexed out of range here, and one containing
    # inf or nan produced a nan "concentration" that compared False by luck rather than by
    # design.  Both are now quiet refusals: a profile whose variation is not a finite number is
    # one this test has nothing to say about.
    if (values.ndim == 0) or (values.shape[0] != len(dense)):
        return quiet
    if not np.all(np.isfinite(values)):
        return quiet

    steps = _variation_steps(values)
    per_interval = steps.reshape(n_ref - 1, n_sub).sum(axis=1)
    endpoints = _variation_steps(values[::n_sub])
    total = float(per_interval.sum())
    scale = float(np.max(np.abs(values)))
    if total <= 1.0e-12*max(scale, 1.0e-300):
        return quiet

    # Summed over ADJACENT PAIRS of reference intervals.  A feature that happens to land on an
    # interval boundary splits its variation between the two, and the single-interval maximum
    # then halves -- measured as a drop in detection from 0.90 to 0.62 at a width of 1e-5.  A
    # spread-out profile gains only the same factor of two, from ~1/n_ref to ~2/n_ref, so the
    # separation is untouched.
    hidden = np.maximum(per_interval - endpoints, 0.0)
    ref = np.linspace(float(l0), float(l1), n_ref)
    if len(hidden) > 1:
        paired = hidden[:-1] + hidden[1:]
        i = int(np.argmax(paired))
        concentration, l_lo, l_hi = float(paired[i]/total), ref[i], ref[i + 2]
    else:
        i = int(np.argmax(hidden))
        concentration, l_lo, l_hi = float(hidden[i]/total), ref[i], ref[i + 1]
    return {'concentration': concentration,
            'hidden': bool(concentration > HIDDEN_FEATURE_CONCENTRATION),
            'l_lo': float(l_lo), 'l_hi': float(l_hi),
            'l_centre': float(0.5*(l_lo + l_hi))}


def oscillation_sampling(H_func: Callable, l0: float, l1: float,
    baselines: Optional[np.ndarray] = None, n_probe: Optional[int] = 8) -> Dict:
    r"""How finely does a scan sample the fastest oscillation on its trajectory?

    A long trajectory carries an enormous oscillation phase, and a scan over it usually samples
    that phase far too coarsely to represent it -- so the returned array is a set of correct
    values that must not be read as a curve.  For a solar or supernova problem the observable is
    normally the **phase-averaged** probability anyway; see :mod:`magnus.avgprob`, whose
    :func:`magnus.avgprob.coherence_report` decides which pairs of eigenvalues have decohered,
    and which are in neither limit.

    **This reports; it never warns.**  A Nyquist criterion is objectively correct and fires on
    essentially everything: measured over the physical profile families, a scan would need 4400
    points on a solar trajectory and **73 000** on a supernova ray to sample the fastest
    oscillation twice per cycle, so 44 of 45 realistic scan sizes are formally aliased
    (``docs/dev/adversarial_batteries/alias_fp.py``).  A warning that fires on 98 % of calls is
    noise however correct each firing is, so the information is offered here and in
    ``strategy_info`` instead, and the caller decides.

    Cost: ``n_probe`` evaluations of ``H_func`` and one ``eigvalsh`` each.  Measured against the
    scan it describes: **5.5 % of the cheapest scan** in the physical population (a 13 ms Earth
    chord), 1.8 % of a 21 ms one, and under 0.1 % of anything substantial
    (``docs/dev/adversarial_batteries/alias_cost.py``).  Callers who do not ask for
    ``strategy_info`` pay **nothing** -- ``oscprob`` only runs this when a report was requested.  Eight samples are enough: the fastest
    oscillation comes from the **largest** eigenvalue spread, which tracks the matter potential
    and varies smoothly, and 8 samples agree with 4096 to 2 % over every family measured.  (The
    *smallest* gap would need far denser sampling, because it has a sharp minimum at an MSW
    resonance -- but that sets the slowest oscillation, not the fastest, and is not what aliases.)

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        The Hamiltonian as a function of position, ``H(l)``, at the energy of interest [eV].
    l0, l1 : float
        Interval the scan traverses.
    baselines : np.ndarray, optional
        The requested baselines.  When given, the spacing and the aliasing verdict are computed
        from them; when omitted, only the trajectory-level quantities are returned.
    n_probe : int, optional
        Samples along the trajectory. Default: 8.

    Returns
    -------
    dict
        ``'oscillation_length'`` (shortest on the trajectory [eV^-1]),
        ``'cycles_over_trajectory'``, ``'nyquist_points'`` (points a scan would need to sample
        the fastest oscillation twice per cycle), and, when ``baselines`` is given,
        ``'spacing'``, ``'cycles_per_step'`` and ``'aliased'``.  Empty dict if the spectrum is
        degenerate or the interval has zero length.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus import adiabatic
        from magnus.hamiltonians import hamiltonians3nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        h_vac = np.asarray(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31']))
        e00 = np.diag([1.0, 0.0, 0.0])
        energy = 10.0e6

        def H(l):
            v = 1.0e-11*np.exp(-np.asarray(l, dtype=float)/gd.L_SCALE_SUN)
            return (1.0/energy)*h_vac + np.asarray(v)[..., None, None]*e00

        print(adiabatic.oscillation_sampling(H, 0.0, 3.0*gd.L_SCALE_SUN, 100))

    Reported, never warned about.  A trajectory of twenty thousand cycles
    sampled at a hundred points is not necessarily wrong -- an averaged
    observable may not care -- but it is worth knowing before trusting an
    instantaneous one.
"""
    l0, l1 = float(l0), float(l1)
    if (l1 == l0) or (n_probe < 1):
        return {}
    spread = 0.0
    for l in np.linspace(l0, l1, int(n_probe)):
        try:
            lam = np.linalg.eigvalsh(np.asarray(H_func(l)))
        except Exception:                  # noqa: BLE001 -- a diagnostic must never break a call
            return {}
        if not np.all(np.isfinite(lam)):
            return {}
        spread = max(spread, float(np.max(lam) - np.min(lam)))
    if spread <= 0.0:
        return {}

    osc = 2.0*np.pi/spread
    span = abs(l1 - l0)
    out = {'oscillation_length': osc,
           'cycles_over_trajectory': span/osc,
           'nyquist_points': int(np.ceil(2.0*span/osc)) + 1}
    if baselines is not None:
        Ls = np.atleast_1d(np.asarray(baselines, dtype=float))
        if Ls.size > 1:
            spacing = float(np.max(np.diff(np.sort(Ls))))
            out.update(spacing=spacing, cycles_per_step=spacing/osc,
                       aliased=bool(spacing > 0.5*osc))
    return out


THRESHOLD0_PROVENANCE = 0.1
r"""float: Module-level constant

The value :func:`hybrid_propagator` starts its adiabaticity-threshold ladder at, recorded here
with its measurement because the measurement is what stopped it being changed.

**What was measured.**  Over 3 profiles x d = 2, 3 x three requested tolerances
(``docs/dev/adversarial_batteries/constants_audit.py``), sweeping ``threshold0`` from 1 down to
1e-3 at a **fixed baseline**: accuracy is identical at every value in 16 of 18 rows, and at
``rtol <= 1e-3`` a lower start is up to **6.5x cheaper** (1.57 s to 0.24 s on a solar profile at
``rtol = 1e-5``).  That is because certifying an empty window list additionally requires
:math:`\gamma` to fit the tolerance (see :data:`GAMMA_TO_ERROR`), so the ladder reaches whatever
threshold the request needs regardless of where it starts; the start only decides how many
iterations that takes, and each one re-runs the detector at doubled ``n_probe`` and the transport
at doubled ``n_points``.  At ``rtol = 1e-2`` the sign flips -- there no window is needed at all,
so starting low opens one that is not, and 0.21 s becomes 0.95 s at d = 3.

**So the right value looked like a rule rather than a constant**, and one was built: start at
:math:`(\mathrm{atol} + \mathrm{rtol})/\texttt{GAMMA\_TO\_ERROR}`, the exact :math:`\gamma`
at which an empty window list stops being certifiable, clipped to ``[min_threshold, 0.1]``.  At
the default tolerance that is 2.4e-03 rather than 0.1.

**It was then rejected, on evidence the sweep above could not produce.**  Scored against
``solve_ivp`` on the package's bit-identity workloads, which include an **energy scan** the
fixed-baseline sweep did not:

=================================== ============ ============== ==============
workload                             ``t0=0.1``   the rule       verdict
=================================== ============ ============== ==============
single point, solar                  1.624e-06    1.184e-10      13711x better
sub-threshold scan, N = 8            3.220e-05    3.814e-05      1.2x worse
**energy scan, fixed baseline**      2.509e-05    **4.954e-04**  **20x worse**
=================================== ============ ============== ==============

All three stay inside the requested 1e-3, but 4.95e-04 spends half the budget where 2.5e-05
spent a fortieth.  The mechanism is visible once looked for: starting low opens a window on the
first iteration, and ``windows_next or windows_prev`` then short-circuits the :math:`\gamma`
check, so agreement can be accepted at a **coarser** transport grid than the old start would
have forced.  The saving and the loss are the same effect seen from two sides.

**The lesson is the one this package keeps re-learning**: the sweep that justified the rule ran
at a fixed baseline, and the row that refuted it was an energy scan.  A population that does not
contain the workload you are about to change is not evidence about it -- the same shape of
mistake that made :data:`GAMMA_TO_ERROR` wrong twice.  0.1 stays until a population that spans
scans as well as points says otherwise.

.. versionadded:: 1.0.0
"""


def _H_on_grid(H_func: Callable, ls: np.ndarray) -> np.ndarray:
    r"""``H_func`` at every position in ``ls``, in one vectorized call where possible.

    The Hamiltonians :mod:`magnus.oscprob` builds are written to accept an array of positions
    and return a stack of matrices (``vcc[..., None, None]*h_matt``), which is the same fast
    path :func:`magnus.magnus.magnus_expansion_multislab` relies on.  Calling such a function
    once per position instead costs the interpreter overhead of a Python loop over the whole
    grid -- measured at 1.2x on an ordinary single-point solar call, which is the difference
    between the resolution test being free and being noticeable.

    Falls back to the loop for a Hamiltonian that is only defined for scalar input, which is a
    supported (if slower) way to write one; :mod:`magnus.magnus` detects the same case and warns
    about it separately.

    .. versionadded:: 1.0.0
    """
    ls = np.asarray(ls, dtype=float)
    try:
        stacked = np.asarray(H_func(ls), dtype=complex)
        if stacked.ndim == 3 and stacked.shape[0] == len(ls):
            return stacked
    except Exception:                      # noqa: BLE001 -- any failure means "not vectorized"
        pass
    return np.array([np.asarray(H_func(l), dtype=complex) for l in ls])


def _profile_is_resolved(H_func: Callable, l0: float, l1: float, n_probe: int) -> bool:
    r"""Whether ``H_func`` is continuous at the scale this module samples it on.

    Everything in this module -- the Hellmann-Feynman derivative ``dH/dl``, the adiabaticity
    parameter built from it, and the parallel transport in :func:`adiabatic_propagator` -- assumes
    ``H_func`` varies smoothly between probe points.  On a piecewise-discontinuous profile that
    assumption fails silently rather than loudly: the finite differences return a finite number,
    the windows come back empty, and the pure adiabatic answer is returned with full confidence.
    Measured on an unmarked density step, that answer was wrong by **0.54** in probability while
    reporting ``certified=True``.

    :mod:`magnus.oscprob` guards this by declining the hybrid strategy when the caller passes
    ``t_breakpoints`` or ``t_slab_edges``.  That guard **fails open**: it declines when the caller
    *tells* it about the discontinuity and accepts when they do not, which is exactly backwards
    with respect to the risk.  This test replaces asking with measuring.

    Method: evaluate ``H`` on the probe grid and again at its midpoints, and compare the largest
    adjacent change at spacing ``h`` with the largest at spacing ``h/2``.  For a :math:`C^1`
    Hamiltonian the latter is about half the former; across a jump the two are equal, because a
    finer grid still straddles the jump.  See :data:`RESOLUTION_RATIO`.

    **What this cannot do.** A feature *narrower than the probe spacing* is generally not sampled
    by either grid, so neither sees it and this test reports "resolved".  That is a limit of any
    fixed grid, not of this test: the cure is a caller-supplied ``t_breakpoints`` at the feature,
    or a larger ``n_probe``.  The test catches discontinuities, which are always straddled by
    some interval, and sharp features comparable to the spacing.

    **Why the caller must not stop at the first failure.**  At one density this test cannot tell
    a genuine jump from a feature that is merely sharp *relative to that density*, and the two
    want opposite treatment: refinement resolves the second and never resolves the first.
    :func:`hybrid_propagator` therefore re-runs this at ``max_n_probe`` before concluding
    anything, which is what distinguishes them.  Skipping that step made a Gaussian of width
    :math:`10^{-3}(l_1-l_0)` -- which the refinement handles exactly, to 1.1e-11 -- get abandoned
    as if it were a step function.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square matrix.
    l0, l1 : float
        Interval over which ``H_func`` is used.
    n_probe : int
        Number of points on the coarse grid (the same grid
        :func:`find_resonance_candidates` uses).

    Returns
    -------
    bool
        False when ``H_func`` shows a jump at this probe scale, True otherwise (including for a
        constant Hamiltonian, where there is nothing to resolve).
    """
    if n_probe < 2:
        return True
    ls = np.linspace(l0, l1, n_probe)
    mids = 0.5*(ls[:-1] + ls[1:])
    Hc = _H_on_grid(H_func, ls)
    Hm = _H_on_grid(H_func, mids)

    # Per interval, how much of the variation falls in each half.
    first = np.max(np.abs(Hm - Hc[:-1]), axis=(1, 2))
    second = np.max(np.abs(Hc[1:] - Hm), axis=(1, 2))
    total = first + second

    # A constant (or numerically constant) Hamiltonian has nothing to resolve.  The floor is
    # scaled to the Hamiltonian itself, not absolute: these matrices carry physical magnitudes
    # spanning many orders.
    scale = np.max(np.abs(Hc))
    live = total > 1.0e-12*scale
    if not np.any(live):
        return True

    # Then drop intervals carrying far less variation than a typical one.  Near a smooth
    # turning point an interval's two halves can differ by orders of magnitude -- one of them
    # rounding to exactly zero -- while the interval as a whole moves by ~3% of typical, and
    # the ratio there is noise, not structure.  Measured on a sine at 28 samples per period,
    # that single interval drove the statistic to 1.0000 and would have declared every
    # oscillating profile discontinuous.  The median is used rather than the maximum precisely
    # so that a jump smaller than the steepest smooth step elsewhere still survives the cut --
    # masking it is the bug this whole test exists to avoid.
    live &= total > 0.25*np.median(total[live])
    if not np.any(live):
        return True

    ratio = np.where(live, np.maximum(first, second)/np.where(total > 0.0, total, 1.0), 0.0)
    flagged = np.where(ratio > RESOLUTION_RATIO)[0]
    if flagged.size == 0:
        return True

    # A concentrated interval is a candidate, not a verdict.  The one interval containing a
    # smooth turning point is genuinely asymmetric -- one half nearly cancels -- and since the
    # statistic above is a maximum over intervals, that single interval decides for the whole
    # profile.  Refining the probe grid re-draws where the extremum falls inside its interval
    # rather than removing the effect, so the caller's coarse-then-fine protocol cannot separate
    # the two either.  Re-sampling the flagged interval *alone* can: a jump is not diluted by
    # refinement and one step still carries all of it, while a smooth feature spreads over every
    # step.  See LOCAL_JUMP_RATIO for the measured separation.
    for i in flagged[np.argsort(-total[flagged])][:MAX_LOCAL_CONFIRMATIONS]:
        xs = np.linspace(ls[i], ls[i + 1], N_LOCAL_CONFIRM)
        steps = np.max(np.abs(np.diff(_H_on_grid(H_func, xs), axis=0)), axis=(1, 2))
        if steps.max() > LOCAL_JUMP_RATIO*total[i]:
            return False
    return True


def _eigs_along(H_func: Callable, ls: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    r"""Diagonalizes ``H_func`` on a grid, with eigenvectors phase-fixed by discrete parallel
    transport.

    Complex dtype is forced before every ``eigh`` call: a genuinely real-valued Hamiltonian (no
    CP violation) is a legitimate special case, not a reason to special-case the code path, and
    ``eigh`` on a real array returns real eigenvectors that later fail to hold a complex parallel-
    transport phase.

    At each step past the first, eigenvector ``k`` is multiplied by the complex phase that makes
    its overlap with the previous step's eigenvector ``k`` real and positive. This is the discrete
    analogue of parallel transport and implicitly captures the geometric (Berry) phase, exactly,
    with no separate formula: the dynamical phase (see :func:`adiabatic_propagator`) and this
    phase-fixing are jointly equivalent to solving the adiabatic evolution equation.

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square matrix.
    ls : np.ndarray
        Positions at which to diagonalize ``H_func``, shape ``(n,)``.

    Returns
    -------
    (np.ndarray, np.ndarray)
        Eigenvalues, shape ``(n, d)``, and phase-fixed eigenvectors, shape ``(n, d, d)`` (each
        ``W[i, :, k]`` is the ``k``-th eigenvector at ``ls[i]``).
    """
    n = len(ls)
    H0 = np.asarray(H_func(ls[0]), dtype=complex)
    d = H0.shape[-1]
    lam = np.empty((n, d))
    W = np.empty((n, d, d), dtype=complex)
    lam[0], W[0] = np.linalg.eigh(H0)
    for i in range(1, n):
        Hi = np.asarray(H_func(ls[i]), dtype=complex)
        li, Wi = np.linalg.eigh(Hi)
        for k in range(d):
            overlap = np.vdot(W[i - 1, :, k], Wi[:, k])
            phase = overlap / abs(overlap) if abs(overlap) > 1e-14 else 1.0
            Wi[:, k] *= np.conj(phase)
        lam[i], W[i] = li, Wi
    return lam, W


def adiabatic_propagator(H_func: Callable, l0: float, l1: float,
    n_points: Optional[int] = 201) -> np.ndarray:
    r"""Computes the evolution operator via pure adiabatic (instantaneous-eigenbasis) transport.

    Diagonalizes ``H_func`` on a grid of ``n_points`` positions between ``l0`` and ``l1``,
    integrates each eigenvalue's dynamical phase with Simpson's rule (trapezoidal quadrature
    leaves a spurious residual that can look like a physics limit but is pure quadrature error),
    and reassembles the evolution operator in the original (flavor) basis:

    .. math::

       U(l_1, l_0) \approx W(l_1)\, \mathrm{diag}\!\left(e^{-i\Phi_k}\right)\, W(l_0)^\dagger ,
       \qquad \Phi_k = \int_{l_0}^{l_1} \lambda_k(l)\, dl ,

    with :math:`W(l)` the (parallel-transported; see ``_eigs_along``) matrix of instantaneous
    eigenvectors of :math:`H(l)` and :math:`\lambda_k(l)` its eigenvalues. This is *exact* in the
    strict adiabatic limit (no eigenvalue crossing or narrowly-avoided crossing along the
    trajectory) and remains unitary by construction (a diagonal phase conjugated by unitary
    matrices) regardless of grid density -- the only thing ``n_points`` controls is how well the
    quadrature/parallel-transport approximate the continuum limit, not whether the result is
    unitary. See :func:`hybrid_propagator` for what to do when the trajectory does cross a
    resonance.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix. May be
        real- or complex-valued.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    n_points : int, optional
        Number of positions at which to diagonalize ``H_func`` between ``l0`` and ``l1``.
        Default: 201.

    Returns
    -------
    np.ndarray
        The evolution operator, exactly unitary.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus import adiabatic
        from magnus.hamiltonians import hamiltonians3nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        h_vac = np.asarray(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31']))
        e00 = np.diag([1.0, 0.0, 0.0])
        energy = 10.0e6

        def H(l):
            v = 1.0e-11*np.exp(-np.asarray(l, dtype=float)/gd.L_SCALE_SUN)
            return (1.0/energy)*h_vac + np.asarray(v)[..., None, None]*e00

        U = np.asarray(adiabatic.adiabatic_propagator(H, 0.0, 3.0*gd.L_SCALE_SUN))

        print('shape', U.shape)
        print('unitary to %.1e' % np.max(np.abs(U.conj().T @ U - np.eye(3))))

    Cheap however large the accumulated phase, because it transports in the
    instantaneous eigenbasis rather than resolving the oscillation.  That is
    exactly why it needs a patch wherever the transport stops being adiabatic.
"""
    if l1 == l0:
        d = np.asarray(H_func(l0)).shape[-1]
        return np.eye(d, dtype=complex)
    ls = np.linspace(l0, l1, n_points)
    lam, W = _eigs_along(H_func, ls)
    d = lam.shape[1]
    Phi = np.array([simpson(lam[:, k], x=ls) for k in range(d)])
    return W[-1] @ np.diag(np.exp(-1j * Phi)) @ W[0].conj().T


def _dH_dl(H_func: Callable, l: float, h: float,
    bounds: Optional[Tuple[float, float]] = None) -> np.ndarray:
    r"""Ordinary real central finite difference of ``H_func`` at ``l``, step ``h``.

    Deliberately **not** complex-step differentiation (``Im[H(l+ih)]/h``): that technique is
    valid only for functions that are real-valued at real input, and ``H_func`` here is routinely
    complex-valued at real ``l`` (e.g., any nonzero CP-violating phase). Applying it anyway
    divides the Hamiltonian's l-independent complex entries by the (tiny) step ``h``, which blows
    up to astronomically wrong derivatives -- a subtle, easy-to-miss trap, caught by comparing
    against this real, always-valid alternative. This real finite difference is robust because
    ``dH/dl`` is always smooth (tied to the smoothness of the underlying density/potential
    profile), independent of how sharp the resulting *eigenvalue* crossing is -- the sharpness
    lives entirely in the gap, in the denominator of ``_point_adiabaticity``, never in
    ``dH/dl`` itself.

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position.
    l : float
        Position at which to evaluate the derivative.
    h : float
        Absolute finite-difference step.
    bounds : (float, float), optional
        Interval ``(l0, l1)`` outside which ``H_func`` must never be evaluated. At a position
        within ``h`` of either end, the stencil is made one-sided so that it stays inside,
        rather than reaching past the boundary. Callers should pass this whenever ``H_func``
        is only defined on the physical domain -- a radial profile undefined for negative
        radius, or one that raises beyond a maximum radius, as
        :func:`magnus.earth.density_matter_func_prem` does. If None, an unclamped central
        difference is used.

    Returns
    -------
    np.ndarray
        Approximation to :math:`dH/dl` at ``l``.
    """
    lm, lp = l - h, l + h
    if bounds is not None:
        lo, hi = bounds
        if lm < lo:
            lm, lp = lo, min(lo + 2.0 * h, hi)
        elif lp > hi:
            lm, lp = max(hi - 2.0 * h, lo), hi
    span = lp - lm
    if span <= 0.0:
        return np.zeros_like(np.asarray(H_func(l), dtype=complex))
    Hp = np.asarray(H_func(lp), dtype=complex)
    Hm = np.asarray(H_func(lm), dtype=complex)
    return (Hp - Hm) / span


def find_resonance_candidates(H_func: Callable, l0: float, l1: float,
    n_probe: Optional[int] = 200, fd_step_frac: Optional[float] = 1e-6,
    info: Optional[Dict] = None) -> List[Dict]:
    r"""Locates every exact eigenvalue-gap critical point of ``H_func`` between ``l0`` and ``l1``.

    For every pair of levels :math:`(j, k)`, scans for sign changes of

    .. math::

       f_{jk}(l) = \langle v_j(l)|\, dH/dl\, |v_j(l)\rangle - \langle v_k(l)|\, dH/dl\, |v_k(l)\rangle ,

    refined by bisection to machine precision in position. By the Hellmann-Feynman theorem,
    :math:`d\lambda_k/dl = \langle v_k|\, dH/dl\, |v_k\rangle` *exactly* (no eigenvector finite
    difference, which would be gauge-ambiguous and fragile), so a sign change of :math:`f_{jk}`
    is an exact critical point of the gap :math:`\lambda_j - \lambda_k` -- a genuine crossing or
    near-crossing candidate, for *any* Hermitian ``H_func`` of *any* dimension, with no assumption
    of a separable or otherwise special structure. Every pair is scanned, so any number of
    simultaneous or sequential resonances (between any pair of levels) are all found.

    A returned candidate is a structural fact about ``H_func`` (an extremum of that pair's gap);
    whether it is actually non-adiabatic (whether it needs a Magnus patch) is a separate question,
    answered by :func:`find_nonadiabatic_windows`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    n_probe : int, optional
        Number of positions on the initial scan grid used to bracket sign changes. Default: 200.
    fd_step_frac : float, optional
        Finite-difference step for ``_dH_dl``, as a fraction of ``l1 - l0``. Default: 1e-6.

        **Provenance.**  A central difference trades truncation error (:math:`\sim h^2`) against
        subtractive cancellation (:math:`\sim \epsilon/h`), so there is an optimum, and this
        default had never been checked against it.  Measured
        (``docs/dev/adversarial_batteries/constants_audit.py``) against the **analytic**
        :math:`dH/dl` -- available because :math:`H = h_\mathrm{vac}/E + C\,n_e(l)\,P_{ee}` and
        the profiles used have closed-form :math:`n_e'`, so the reference carries no error of
        its own -- over d = 2, 3 and 10-50 MeV, as a fraction of the largest :math:`|dH/dl|` on
        the path:

        ========================= ============ ==================== =========================
        profile                    optimum      error at **1e-6**    error at the optimum
        ========================= ============ ==================== =========================
        solar exponential          1e-5         6.8e-11              1.9e-11
        sinusoid, period span/7    1e-6         3.4e-10              3.4e-10
        Gaussian bump, w = 1e-2    1e-7 - 1e-8  2.5e-09              2.6e-11
        ========================= ============ ==================== =========================

        **The optimum is not a single number** -- it moves with the profile's shortest length
        scale, which is what the theory predicts.  But the curve is flat enough that this does
        not matter: anywhere in :math:`10^{-8}` to :math:`10^{-5}` the error stays below
        :math:`3\times10^{-9}` relative on every profile measured, which is six orders below
        anything that could move a probability at the tolerances this package works to.  Outside
        that band it degrades fast in both directions -- 1.6e-04 at :math:`10^{-12}`
        (cancellation), 0.23 at :math:`10^{-2}` (truncation) -- so the band, not the value, is
        the thing to preserve.
    info : dict, optional
        If given, filled in place with the probe-grid quantities this function had to compute
        anyway -- ``'ls'`` (the grid), ``'lam'``, ``'W'`` (eigenvalues and eigenvectors, shapes
        ``(n, d)`` and ``(n, d, d)``) and ``'dH'``. :func:`find_nonadiabatic_windows` sweeps
        :math:`\gamma` on exactly this grid with exactly this finite-difference step, so
        without this it would recompute all of it: ~600 extra Hamiltonian evaluations and a
        second eigendecomposition, which measured as **1.4x** on an ordinary single-point solar
        call. Default: None.

    Returns
    -------
    list of dict
        One entry per candidate, with keys ``'l'`` (position), ``'j'``, ``'k'`` (the level
        indices, ``j < k``), and ``'gap'`` (:math:`\lambda_k - \lambda_j` at that position),
        sorted by position.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus import adiabatic
        from magnus.hamiltonians import hamiltonians3nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        h_vac = np.asarray(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31']))
        e00 = np.diag([1.0, 0.0, 0.0])
        energy = 10.0e6

        def H(l):
            v = 1.0e-11*np.exp(-np.asarray(l, dtype=float)/gd.L_SCALE_SUN)
            return (1.0/energy)*h_vac + np.asarray(v)[..., None, None]*e00

        for c in adiabatic.find_resonance_candidates(H, 0.0, 3.0*gd.L_SCALE_SUN):
            print('levels %d-%d cross at l/l_scale = %.3f, gap %.2e eV'
                  % (c['j'], c['k'], c['l']/gd.L_SCALE_SUN, c['gap']))

    A candidate is a critical point of a pairwise gap, found exactly through
    Hellmann-Feynman rather than by scanning for a minimum.  Whether it is
    actually non-adiabatic is a separate question --
    :func:`find_nonadiabatic_windows` answers it.
"""
    ls = np.linspace(l0, l1, n_probe)
    h = (l1 - l0) * fd_step_frac
    bounds = (l0, l1)
    # One vectorized call and one batched eigendecomposition, rather than a Python loop over
    # the grid: same values, and the difference is visible at the entry point.
    Hs = _H_on_grid(H_func, ls)
    d = Hs.shape[-1]
    lam, W = np.linalg.eigh(Hs)
    dH = np.array([_dH_dl(H_func, l, h, bounds) for l in ls])
    if info is not None:
        info.update(ls=ls, lam=lam, W=W, dH=dH)

    def f_pair(l: float, j: int, k: int) -> float:
        H = np.asarray(H_func(l), dtype=complex)
        _, Wi = np.linalg.eigh(H)
        dHl = _dH_dl(H_func, l, h, bounds)
        vj, vk = Wi[:, j], Wi[:, k]
        return float(np.real(np.vdot(vj, dHl @ vj) - np.vdot(vk, dHl @ vk)))

    candidates = []
    for j in range(d):
        for k in range(j + 1, d):
            fjk = np.real(np.einsum('ni,nij,nj->n', np.conj(W[:, :, j]), dH, W[:, :, j])
                - np.einsum('ni,nij,nj->n', np.conj(W[:, :, k]), dH, W[:, :, k]))
            sgn = np.sign(fjk)
            changes = np.where(np.diff(sgn) != 0)[0]
            for idx in changes:
                a, b = ls[idx], ls[idx + 1]
                fa, fb = f_pair(a, j, k), f_pair(b, j, k)
                if fa == 0.0:
                    l_star = a
                elif fb == 0.0:
                    l_star = b
                else:
                    for _ in range(60):
                        m = 0.5 * (a + b)
                        fm = f_pair(m, j, k)
                        if np.sign(fm) == np.sign(fa):
                            a, fa = m, fm
                        else:
                            b, fb = m, fm
                    l_star = 0.5 * (a + b)
                H_star = np.asarray(H_func(l_star), dtype=complex)
                lam_star = np.linalg.eigvalsh(H_star)
                gap = float(lam_star[k] - lam_star[j])
                candidates.append({'l': l_star, 'j': j, 'k': k, 'gap': gap})
    candidates.sort(key=lambda c: c['l'])
    return candidates


def _point_adiabaticity(H_func: Callable, l: float, j: int, k: int, fd_step: float,
    bounds: Optional[Tuple[float, float]] = None) -> float:
    r"""Adiabaticity parameter :math:`\gamma_{jk}(l) = |\langle v_j|\, dH/dl\, |v_k\rangle| / (\lambda_k - \lambda_j)^2`
    (Landau-Zener form), computed exactly from the Hellmann-Feynman off-diagonal matrix element --
    no eigenvector finite difference. Large :math:`\gamma` signals a narrowly-avoided (or exact)
    crossing where the adiabatic approximation breaks down; ``fd_step`` is an *absolute* step
    (unlike ``fd_step_frac`` elsewhere), since callers evaluate this at positions found by
    bisection, arbitrarily close together. ``bounds``, if given, keeps the finite-difference
    stencil inside the physical domain (see ``_dH_dl``).
    """
    H = np.asarray(H_func(l), dtype=complex)
    lam, W = np.linalg.eigh(H)
    dH = _dH_dl(H_func, l, fd_step, bounds)
    vj, vk = W[:, j], W[:, k]
    coupling = np.abs(np.vdot(vj, dH @ vk))
    gap = abs(lam[k] - lam[j])
    return coupling / gap**2 if gap > 0 else np.inf


def _estimate_window_bounds(H_func: Callable, l_star: float, j: int, k: int, l0: float, l1: float,
    threshold: float, fd_step: float, safety_factor: Optional[float] = 2.0,
    max_doublings: Optional[int] = 60) -> Tuple[float, float]:
    r"""Grows a window outward from a candidate position until the adiabaticity parameter drops
    below ``threshold``, then pads it by ``safety_factor``.

    Growing by doubling (rather than tying the window width to the search-grid spacing used by
    :func:`find_resonance_candidates`) makes the window a property of the physical transition
    width alone: the same physical case gives the same window regardless of ``n_probe``.
    """
    def grow(sign: float) -> float:
        step = fd_step
        l_edge = l_star
        for _ in range(max_doublings):
            l_try = l_star + sign * step
            if sign > 0 and l_try >= l1:
                return l1
            if sign < 0 and l_try <= l0:
                return l0
            if _point_adiabaticity(H_func, l_try, j, k, fd_step, (l0, l1)) < threshold:
                l_edge = l_try
                break
            step *= 2.0
        else:
            return l1 if sign > 0 else l0
        width = abs(l_edge - l_star)
        l_pad = l_star + sign * safety_factor * width
        return min(l_pad, l1) if sign > 0 else max(l_pad, l0)
    return grow(-1.0), grow(+1.0)


def find_nonadiabatic_windows(H_func: Callable, l0: float, l1: float,
    threshold: Optional[float] = 0.1, n_probe: Optional[int] = 200,
    fd_step_frac: Optional[float] = 1e-6,
    info: Optional[Dict] = None) -> Tuple[List[Tuple[float, float]], List[Dict]]:
    r"""Finds every position window along ``[l0, l1]`` where ``H_func`` needs a Magnus patch.

    Calls :func:`find_resonance_candidates`, evaluates the adiabaticity parameter
    :math:`\gamma_{jk}` (see ``_point_adiabaticity``) at each candidate, grows a window around
    every candidate with :math:`\gamma_{jk} > \text{threshold}` (see ``_estimate_window_bounds``),
    and merges any windows that overlap or touch -- so two (or more) resonances close enough
    together are correctly folded into a single patch, rather than either double-counted or
    (worse) silently dropped. This works for any number of simultaneous or sequential resonances,
    between any pair of levels.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    threshold : float, optional
        Adiabaticity parameter above which a candidate is treated as non-adiabatic. Default: 0.1.
    n_probe : int, optional
        Forwarded to :func:`find_resonance_candidates`. Default: 200.
    fd_step_frac : float, optional
        Forwarded to :func:`find_resonance_candidates`. Default: 1e-6.
    info : dict, optional
        If given, filled in place with diagnostics about this call, following the same
        out-parameter convention as ``convergence_info`` in :func:`magnus.oscprob.osc_prob`.
        Currently one key, ``'gamma_max'``: the largest adiabaticity parameter seen anywhere on
        the probe grid or at any candidate, over every level pair. It is ``inf`` if some pair's
        gap vanishes exactly. :func:`hybrid_propagator` uses it to decide whether an *empty*
        window list may be certified -- without it, "no window opened" is indistinguishable from
        "no window was looked for hard enough". Default: None.

    Returns
    -------
    (list of (float, float), list of dict)
        The merged, non-overlapping windows (each a ``(l_b, l_c)`` pair, sorted by position), and
        the candidate list from :func:`find_resonance_candidates`, each entry additionally
        carrying its evaluated ``'gamma'``.
    """
    probe = {}
    candidates = find_resonance_candidates(H_func, l0, l1, n_probe=n_probe,
        fd_step_frac=fd_step_frac, info=probe)
    fd_step = (l1 - l0) * fd_step_frac
    windows = []
    gamma_max = 0.0
    for c in candidates:
        gamma = _point_adiabaticity(H_func, c['l'], c['j'], c['k'], fd_step, (l0, l1))
        c['gamma'] = gamma
        gamma_max = max(gamma_max, gamma)
        if gamma > threshold:
            l_b, l_c = _estimate_window_bounds(H_func, c['l'], c['j'], c['k'], l0, l1, threshold,
                fd_step)
            windows.append([l_b, l_c])

    # Sweep the probe grid as well, not only the gap extrema above.  A gap extremum is where the
    # *gap* is stationary, which is not where gamma = |<v_j|dH/dl|v_k>| / gap^2 peaks: on a
    # rapidly varying profile the coupling can be large between the extrema, and evaluating only
    # at them understates the maximum badly -- measured at 196x on a 3nu NSI profile modulated by
    # a strong sine (3.6e-04 at the extrema against 7.0e-02 along the path).
    #
    # Without this the failure is silent rather than merely inaccurate.  No window ever opens, so
    # successive refinements differ only in the adiabatic-transport grid, converge to the same
    # wrong adiabatic limit, agree with each other, and hybrid_propagator certifies a result that
    # was off by 4.3e-02 against solve_ivp.  Lowering the threshold cannot rescue it, because the
    # threshold is only ever compared against values sampled where gamma happens to be small.
    #
    # Genuinely reuses what find_resonance_candidates already computed, rather than repeating
    # it.  The two agree exactly by construction -- same np.linspace(l0, l1, n_probe), same
    # (l1-l0)*fd_step_frac step, same bounds -- so this is an identity, not an approximation.
    # It was written as a fresh pass, which doubled the detector's Hamiltonian evaluations
    # (~600 of them) and added a second eigendecomposition for nothing: 1.4x at the entry
    # point, and the comment here used to claim the reuse that the code did not do.
    ls_probe = probe['ls']
    lam_p, W_p, dH_p = probe['lam'], probe['W'], probe['dH']
    d_p = lam_p.shape[-1]
    for j in range(d_p):
        for k in range(j + 1, d_p):
            vj, vk = W_p[:, :, j], W_p[:, :, k]
            coupling = np.abs(np.einsum('ni,nij,nj->n', np.conj(vj), dH_p, vk))
            gap = np.abs(lam_p[:, k] - lam_p[:, j])
            gamma_p = np.where(gap > 0.0, coupling/np.where(gap > 0.0, gap, 1.0)**2, np.inf)
            if gamma_p.size:
                gamma_max = max(gamma_max, float(np.max(gamma_p)))
            over = np.where(gamma_p > threshold)[0]
            if over.size == 0:
                continue
            # One window per *contiguous run* of exceedance, grown from that run's peak --
            # not one per exceeding point.  Growing from every point and merging pads each
            # window outward independently, and enough of them bridge the quiet stretch
            # between two genuinely separate crossings: on the two-crossing fixture in
            # tests/test_avgprob.py that collapsed 2 windows into 1 spanning almost the whole
            # trajectory, which destroys the crossing structure the averaged-probability
            # report is built on.  Runs also make this cheap: two growth searches there
            # rather than forty-two.
            for run in np.split(over, np.where(np.diff(over) != 1)[0] + 1):
                peak = int(run[np.argmax(gamma_p[run])])
                l_b, l_c = _estimate_window_bounds(H_func, float(ls_probe[peak]), j, k, l0, l1,
                    threshold, fd_step)
                windows.append([l_b, l_c])

    if info is not None:
        info['gamma_max'] = gamma_max

    if not windows:
        return [], candidates
    windows.sort()
    merged = [windows[0]]
    for w in windows[1:]:
        if w[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], w[1])
        else:
            merged.append(w)
    return [tuple(w) for w in merged], candidates


def _local_evolution_operator(H_func: Callable, l_b: float, l_c: float, magnus_exp_order: int,
    integration_method: str, n_slabs0: Optional[int] = 400, max_n_slabs: Optional[int] = 32_768,
    patch_atol: Optional[float] = 1e-7) -> Tuple[np.ndarray, bool]:
    r"""Computes the (exact, not adiabatic) evolution operator across a single non-adiabatic
    window, via the package's own Magnus kernel, doubling the slab count until convergence.

    ``patch_atol`` **provenance.**  Swept over 1e-5 to 1e-9 across 18 workloads
    (``constants_audit2.py``).  The worst error is 4.49e-04 at 1e-5, 1e-6 and the 1e-7 default,
    and 3.04e-04 at 1e-8 -- so the default has better than an order of magnitude of margin in
    the loose direction.  At **1e-9 it breaks**, and instructively: most rows improve sharply
    (a multi-resonance baseline scan goes 3.04e-04 to 1.15e-09) but one energy scan goes to
    **2.08e-02**, twenty times outside tolerance.  The mechanism is not this constant -- at 1e-9
    the patch cannot converge inside ``max_n_slabs``, so the hybrid strategy declines, and the
    energy-batched separable engine answers instead and is that much worse on that profile.  It
    warns (``MagnusConvergenceWarning``), so it is loud rather than silent.  What the row really
    measures is **fallback quality**, and it is recorded here because that is where it was
    found; the default is comfortably clear of it.

    ``n_slabs0`` **provenance.**  Swept over 100 to 1600 on the same workloads: worst error
    4.49e-04 at every value, since this only sets where the doubling starts.

    ``max_n_slabs`` is a **statement about when this method stops applying**, not a performance
    knob.  A patch is supposed to be a short, local repair of a narrow region where adiabatic
    transport fails; if it needs more slabs than a plain Magnus integration of the entire
    trajectory would, then the non-adiabatic region is not narrow and the hybrid strategy has no
    reason to exist for that request.  Returning ``False`` there is the honest answer: it makes
    :func:`hybrid_propagator` report ``certified=False``, and :func:`magnus.oscprob.osc_prob`
    then falls through to the general Magnus path -- which handles such cases correctly, and in
    the measured case 70x faster.

    Measured slab counts at convergence, which is where 32768 comes from:

    ======================================================= ===================
    patch                                                   slabs at convergence
    ======================================================= ===================
    multi-resonance profile, 8 patches (correct, 0.52 s)    800 - **12 800**
    sub-threshold bump, 2 patches (correct, 0.38 s)         800 - 3 200
    solar over 5 scale heights at rtol 1e-4, 1 patch
    covering 88% of the path                                **102 400**
    ======================================================= ===================

    The gap is a factor of eight, and 32768 sits inside it: 2.6x above the largest legitimate
    patch measured, 3.1x below the one that should decline.  The previous default of 500 000 was
    above everything, so nothing ever declined and that last row was patched at 45x the cost of
    simply handing the request to the general path.

    Uses :func:`magnus.magnus.magnus_expansion_multislab` directly (not
    :func:`magnus.oscprob.compute_evolution_operator_multiple_slabs`, to keep this module free of
    any dependency on :mod:`magnus.oscprob`).

    Returns
    -------
    (np.ndarray, bool)
        The evolution operator across ``[l_b, l_c]``, and whether it converged (agreed with the
        previous, half-as-fine slab count) within ``max_n_slabs``.
    """
    def U_at(n: int) -> np.ndarray:
        edges_lin = np.linspace(l_b, l_c, n + 1)
        edges = np.column_stack([edges_lin[:-1], edges_lin[1:]])

        def hh(t):
            return -1j * np.asarray(H_func(t))

        U_chain = magnuscore.magnus_expansion_multislab(hh, edges, n_tpts_per_slab=2,
            order=magnus_exp_order, integration_method=integration_method)
        return magnuscore.ordered_product(U_chain)

    n_slabs = n_slabs0
    U_prev = U_at(n_slabs)
    while n_slabs < max_n_slabs:
        n_slabs *= 2
        U_next = U_at(n_slabs)
        if np.max(np.abs(U_next - U_prev)) <= patch_atol:
            return U_next, True
        U_prev = U_next
    return U_prev, False


def _hybrid_propagator_once(H_func: Callable, l0: float, l1: float, threshold: float,
    n_probe: int, n_points: int, fd_step_frac: float, magnus_exp_order: int,
    integration_method: str) -> Tuple[np.ndarray, List[Tuple[float, float]], bool, float]:
    r"""One evaluation of the hybrid propagator at a fixed set of internal tolerance knobs (see
    :func:`hybrid_propagator` for the self-certifying refinement built on top of this).

    Returns the operator, the windows used, whether every local patch converged, and the largest
    adiabaticity parameter seen on the probe grid (which the caller needs to judge an *empty*
    window list -- see :data:`GAMMA_TO_ERROR`)."""
    info = {}
    windows, _ = find_nonadiabatic_windows(H_func, l0, l1, threshold=threshold, n_probe=n_probe,
        fd_step_frac=fd_step_frac, info=info)
    gamma_max = info.get('gamma_max', 0.0)
    if not windows:
        return (adiabatic_propagator(H_func, l0, l1, n_points=n_points), windows, True,
                gamma_max)
    d = np.asarray(H_func(l0), dtype=complex).shape[-1]
    U_total = np.eye(d, dtype=complex)
    cursor = l0
    all_patches_converged = True
    for (l_b, l_c) in windows:
        U_total = adiabatic_propagator(H_func, cursor, l_b, n_points=n_points) @ U_total
        U_patch, ok = _local_evolution_operator(H_func, l_b, l_c, magnus_exp_order,
            integration_method)
        all_patches_converged = all_patches_converged and ok
        U_total = U_patch @ U_total
        cursor = l_c
    U_total = adiabatic_propagator(H_func, cursor, l1, n_points=n_points) @ U_total
    return U_total, windows, all_patches_converged, gamma_max


def hybrid_propagator(H_func: Callable, l0: float, l1: float, rtol: Optional[float] = 1.e-3,
    atol: Optional[float] = 1.e-3, magnus_exp_order: Optional[int] = 6,
    integration_method: Optional[str] = 'gl', threshold0: Optional[float] = 0.1,
    min_threshold: Optional[float] = 1.e-6, n_probe0: Optional[int] = 200,
    max_n_probe: Optional[int] = 6400, n_points0: Optional[int] = 201,
    max_n_points: Optional[int] = 12864, fd_step_frac: Optional[float] = 1.e-6,
    max_iters: Optional[int] = 12,
    info: Optional[Dict] = None) -> Tuple[np.ndarray, List[Tuple[float, float]], bool]:
    r"""Computes the evolution operator via adiabatic transport, with a Magnus patch at every
    non-adiabatic window, self-certified against successive refinement of every internal
    tolerance knob.

    This is the main entry point of this module (see :doc:`/adiabatic_strategy` for the full
    derivation and validation). Given any Hermitian ``H_func`` of any dimension:

    #. Locates every non-adiabatic window along ``[l0, l1]`` (see
       :func:`find_nonadiabatic_windows`).
    #. If there are none, returns the pure adiabatic-transport operator (see
       :func:`adiabatic_propagator`).
    #. Otherwise, composes adiabatic transport between windows with an exact local Magnus patch
       *inside* each window (see ``_local_evolution_operator``), using the exact composition
       law of quantum evolution, :math:`U(l_2, l_0) = U(l_2, l_1)\, U(l_1, l_0)`, so the result is
       exactly unitary regardless of any approximation's accuracy.
    #. Self-certifies the result: a single fixed adiabaticity ``threshold`` (deciding which
       candidates count as non-adiabatic) is not safe in general -- too loose, and a genuine
       resonance is patched too narrowly or missed; too tight, and windows are patched
       needlessly, at some (still usually modest) extra cost. Rather than trust one fixed value,
       the whole computation (window threshold, adiabatic-transport grid density, and the probe
       grid used to *locate* candidates) is repeated with the knobs tightened together
       (threshold divided by 3, ``n_points``/``n_probe`` doubled) until two successive results
       agree within ``rtol``/``atol``, mirroring the successive-refinement discipline
       :func:`magnus.oscprob.osc_prob` already uses for the number of slabs.

       Each knob stops at its own ceiling (``min_threshold``, ``max_n_probe``,
       ``max_n_points``), which they reach at different iterations, so the later iterations
       tighten fewer knobs than the earlier ones. Once *all* of them have saturated, a further
       iteration would recompute bit-identical inputs and the agreement test would pass
       trivially, comparing a result with itself; the loop therefore stops at that point and
       reports ``certified=False`` rather than certifying on the strength of a comparison that
       carries no information.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix. May be
        real- or complex-valued, of any dimension.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    rtol : float, optional
        Relative tolerance on the *agreement* between successive refinement levels, and on the
        adiabaticity bound in ``_certified``.  Default: 1e-3.  Like every tolerance in this
        package it is a stopping rule rather than a guaranteed accuracy: the loop halts when
        two successive levels agree, and no error of the returned operator is ever estimated.
        See the ``rtol`` entry of :func:`magnus.oscprob.osc_prob` for what that does and does
        not promise.  ``certified`` is the flag that says whether the loop stopped because it
        agreed or because it ran out of room.
    atol : float, optional
        Absolute tolerance on the same agreement; see ``rtol``.  Default: 1e-3.
    magnus_exp_order : int, optional
        Magnus expansion order used for the local patch inside each non-adiabatic window.
        Default: 6.
    integration_method : str, optional
        Integration method used for the local patch ('gl', 'trapezoid', or 'simpson').
        Default: 'gl'.
    threshold0 : float, optional
        Starting adiabaticity threshold. Default: 0.1.  See :data:`THRESHOLD0_PROVENANCE` for
        what it was measured to do, and for why a tolerance-derived rule was built, tested and
        then **rejected**.

        **Provenance, and a reframing.**  This constant used to decide *correctness*: if
        :math:`\gamma` never crossed it, no window opened, successive refinements agreed with
        each other, and the answer was certified while wrong by up to 1.8e-02.  Since
        certification of an empty window list additionally requires :data:`GAMMA_TO_ERROR` *
        :math:`\gamma_\max` to fit the tolerance, that failure mode is closed and **this is now
        a cost knob**: it sets where the ``threshold /= 3`` ladder starts, not whether the
        result may be believed.

        Measured (``docs/dev/adversarial_batteries/constants_audit.py``) over 3 profiles x
        d = 2, 3 x three requested tolerances, sweeping ``threshold0`` from 1 down to 1e-3:

        * **Accuracy is identical at every value** in 16 of 18 rows.  The ladder reaches
          whatever threshold the tolerance requires regardless of where it starts.
        * **At ``rtol <= 1e-3``, lower is monotonically cheaper** -- up to **6.5x** (1.57 s to
          0.24 s on a solar profile at ``rtol = 1e-5``) -- because every step of the ladder
          re-runs the detector at doubled ``n_probe`` and the transport at doubled
          ``n_points``.  Starting low skips those iterations.
        * **At ``rtol = 1e-2`` the sign flips.**  There the tolerance does not require a window
          at all, so a low threshold opens one that is not needed: 0.21 s becomes 0.95 s at
          d = 3, buying an accuracy improvement (2.67e-03 to 1.44e-09) nobody asked for.

        So the brief's hypothesis is confirmed: **the right value is a rule, not a constant** --
        low when the requested tolerance is tight, high when it is loose.  0.1 is not the
        optimum at the default ``rtol = 1e-3``, where 0.01 is 2-3x cheaper at identical
        accuracy.  **It is deliberately left unchanged**: three profiles at one energy is
        precisely the size of population that made :data:`GAMMA_TO_ERROR` wrong twice, and
        retuning a default on it would repeat that mistake rather than learn from it.  The
        measurement is recorded here so the next person starts from evidence.
    min_threshold : float, optional
        Floor below which the threshold is not tightened further. Default: 1e-6.

        **Provenance.**  Swept over 1e-4 to 1e-8 across 18 ordinary workloads: the worst error
        is 4.49e-04 at every value, because the ladder stops long before reaching the floor.
        Even at ``rtol = atol = 1e-12`` on a multi-resonance profile it converges in 7
        iterations with a window open, and every value from 1e-4 to 1e-10 gives an identical
        answer in identical time.

        **The regime this constant governs**, found by construction rather than assumed: the
        floor is reached only when :math:`\gamma_\max` is *below* it -- so no window can ever
        open, however far the threshold falls -- **and** the requested tolerance is tighter than
        ``GAMMA_TO_ERROR`` :math:`\times \gamma_\max`, so the :math:`\gamma` rule cannot
        certify either.  An almost-flat profile
        (:math:`\gamma_\max = 3\times10^{-7}`) at ``rtol = atol = 1e-9`` satisfies both:

        ================= ========== ======== ========== ======
        ``min_threshold`` error      windows  iterations time
        ================= ========== ======== ========== ======
        1e-4              8.49e-13   0        9          3.3 s
        **1e-6**          8.49e-13   0        13         7.4 s
        1e-8              3.24e-12   1        13         7.7 s
        1e-10             3.24e-12   1        13         7.9 s
        ================= ========== ======== ========== ======

        So it does change behaviour there -- below :math:`\gamma_\max` a window opens -- but
        **not usefully**: the result is ``certified=False`` at every value, and the error is
        three orders inside the requested tolerance either way, with the window costing a
        factor of 2.4 in time and making the answer very slightly *worse*.  The floor decides
        how much work is done in a regime where the answer is already good and known to be
        uncertified; it does not decide correctness anywhere measured.
    n_probe0 : int, optional
        Starting number of positions used to locate resonance candidates. Default: 200.

        **Provenance.**  Swept over 50, 100, 200, 400, 800 across 18 workloads spanning single
        points, baseline scans **and energy scans**, 3 profile families, d = 2 and 3
        (``docs/dev/adversarial_batteries/constants_audit2.py``).  Worst error over all
        workloads: 4.98e-04, 5.22e-04, **4.49e-04**, 3.38e-04, 3.36e-04 -- flat within a factor
        of 1.6 and inside the requested 1e-3 everywhere.  Not load-bearing: the refinement
        doubles it, so the starting value only shifts which iteration finds a given feature.
    max_n_probe : int, optional
        Ceiling on the probe grid density. Default: 6400.  A cost ceiling rather than a
        calibration -- reaching it is reported rather than absorbed -- and it also sets what
        :func:`find_hidden_features` treats as resolvable.
    n_points0 : int, optional
        Starting number of positions used for adiabatic-transport quadrature. Default: 201.

        **Provenance.**  Swept over 51, 101, 201, 401, 801 on the same 18 workloads: the worst
        error is **4.49e-04 at every value, identical to three digits**.  The refinement doubles
        this too, so the starting value is invisible in the answer; it buys only iterations.
    max_n_points : int, optional
        Ceiling on the adiabatic-transport grid density. Default: 12864.  A cost ceiling, as
        ``max_n_probe``.
    fd_step_frac : float, optional
        Finite-difference step for the Hellmann-Feynman diagnostics, as a fraction of
        ``l1 - l0``. Default: 1e-6.
    max_iters : int, optional
        Maximum number of refinement iterations. Default: 12.
    info : dict, optional
        If given, filled in place with why this call ended as it did, following the same
        out-parameter convention as ``convergence_info`` in :func:`magnus.oscprob.osc_prob`.
        Keys: ``'resolved'`` (whether ``H_func`` passed the probe-scale resolution test -- see
        ``_profile_is_resolved``), ``'gamma_max'``, ``'n_windows'``, ``'iterations'``, and
        ``'patches_converged'``.  ``certified=False`` on its own does not say *which* of these
        failed, and the cures are different: an unresolved profile wants ``t_breakpoints``, an
        exhausted refinement wants a looser tolerance.  :mod:`magnus.oscprob` uses
        ``'resolved'`` to raise :class:`magnus.oscprob.UnmarkedDiscontinuityWarning` on the
        hybrid path instead of declining in silence. Default: None.

    Returns
    -------
    (np.ndarray, list of (float, float), bool)
        The evolution operator (exactly unitary regardless of ``certified``), the non-adiabatic
        windows used in the last iteration, and whether the result is certified (``True``).
        ``certified`` is ``False`` if the refinement exhausted ``max_iters``, if every knob
        reached its ceiling before two successive results agreed, if a local patch failed to
        converge within its own slab cap, or if ``H_func`` is not resolved at the probe scale
        (see ``_profile_is_resolved``) -- in all of these the returned operator is the best
        available estimate, still exactly unitary, but its accuracy is not certified to the
        requested tolerance.

    Notes
    -----
    Two successive results agreeing is **necessary but not sufficient**, and step 4 above states
    the reason narrowly. When no window opens at all, successive iterations differ only in the
    adiabatic-transport grid: they converge to the same adiabatic limit and agree with each
    other whether or not that limit is the right answer. Certifying an empty window list
    therefore additionally requires the adiabaticity parameter itself to be small enough for the
    requested tolerance (see :data:`GAMMA_TO_ERROR`); otherwise the loop keeps lowering the
    threshold until a window does open. Without that, a profile whose :math:`\gamma` stays just
    below ``threshold0`` everywhere is certified while wrong -- measured at 1.8e-02 against a
    requested 1e-3.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus import adiabatic
        from magnus.hamiltonians import hamiltonians3nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        h_vac = np.asarray(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31']))
        e00 = np.diag([1.0, 0.0, 0.0])
        energy = 10.0e6

        def H(l):
            v = 1.0e-11*np.exp(-np.asarray(l, dtype=float)/gd.L_SCALE_SUN)
            return (1.0/energy)*h_vac + np.asarray(v)[..., None, None]*e00

        U, windows, certified = adiabatic.hybrid_propagator(
            H, 0.0, 3.0*gd.L_SCALE_SUN)

        print('non-adiabatic windows :', len(windows))
        print('certified             :', certified)
        print('unitary to %.1e'
              % np.max(np.abs(np.asarray(U).conj().T @ U - np.eye(3))))

    ``certified`` is the value to check: it says the result agreed with itself
    under tightening, not that it is correct.  Zero windows means the whole
    trajectory was adiabatic, which for the real solar mixing angle is the
    usual answer.
"""
    threshold, n_probe, n_points = threshold0, n_probe0, n_points0

    # Everything below finite-differences H_func between probe points and assumes the result
    # means something.  On a piecewise-discontinuous profile it does not, and the failure is
    # silent: no window opens, the adiabatic answers at successive n_points agree with each
    # other, and this function would certify a result measured 0.54 wrong in probability.
    # Refusing to certify is enough to make the package safe, because osc_prob's
    # strategy='auto' treats an uncertified hybrid result as "this method does not fit" and
    # falls through to the general Magnus path, which handles such profiles correctly.
    #
    # Confirmed at max_n_probe before it is believed.  The cheap test at n_probe0 cannot
    # separate a genuine jump from a feature that is merely sharp at *that* density, and the
    # refinement below resolves the second while never resolving the first -- so failing at
    # one density is a reason to look harder, not a verdict.  Without the second stage a
    # Gaussian of width 1e-3 (l1-l0), which this module otherwise answers to 1.1e-11, was
    # abandoned as though it were a step.  The confirmation runs only on profiles that already
    # failed the cheap test, so ordinary calls never pay for it.
    resolved = (_profile_is_resolved(H_func, l0, l1, n_probe0)
                or _profile_is_resolved(H_func, l0, l1, max_n_probe))

    def report(n_windows: int, gamma_max: float, iterations: int, patches_ok: bool):
        if info is not None:
            info.update(resolved=bool(resolved), gamma_max=float(gamma_max),
                        n_windows=int(n_windows), iterations=int(iterations),
                        patches_converged=bool(patches_ok))

    U_prev, windows_prev, ok_prev, gamma_prev = _hybrid_propagator_once(H_func, l0, l1,
        threshold, n_probe, n_points, fd_step_frac, magnus_exp_order, integration_method)
    if not ok_prev or not resolved:
        report(len(windows_prev), gamma_prev, 1, ok_prev)
        return U_prev, windows_prev, False

    def adiabatic_is_good_enough(gamma_max: float) -> bool:
        r"""Whether a result with NO window may be certified on the strength of gamma alone.

        When no window opens, successive refinements differ only in the adiabatic-transport
        grid, so they converge to the same adiabatic limit and agree with each other whether or
        not that limit is right -- the agreement test carries no information about the thing
        that actually went wrong.  What does carry information is how non-adiabatic the path
        was: see :data:`GAMMA_TO_ERROR` for the measured relation between gamma_max and the
        error of the pure adiabatic answer.
        """
        return bool(GAMMA_TO_ERROR*gamma_max <= atol + rtol)

    iterations = 1
    for _ in range(max_iters):
        iterations += 1
        knobs_prev = (threshold, n_probe, n_points)
        threshold = max(threshold / 3.0, min_threshold)
        n_probe = min(n_probe * 2, max_n_probe)
        n_points = min(n_points * 2, max_n_points)
        if (threshold, n_probe, n_points) == knobs_prev:
            # Every knob has hit its ceiling. _hybrid_propagator_once is deterministic, so
            # rerunning it here would reproduce U_prev exactly and the agreement test below
            # would pass on a comparison of a result with itself -- which is no evidence of
            # convergence at all. Stop and report the result as uncertified instead.
            break
        U_next, windows_next, ok_next, gamma_next = _hybrid_propagator_once(H_func, l0, l1,
            threshold, n_probe, n_points, fd_step_frac, magnus_exp_order, integration_method)
        if not ok_next:
            report(len(windows_next), gamma_next, iterations, False)
            return U_next, windows_next, False
        if np.max(np.abs(U_next - U_prev)) <= atol + rtol * np.max(np.abs(U_prev)):
            # Agreement is necessary but not sufficient. If neither result patched anything,
            # both are pure adiabatic transport and their agreement is self-fulfilling; accept
            # it only when gamma says the adiabatic approximation was itself good enough. If
            # it does not, fall through and keep lowering the threshold, which is guaranteed to
            # open a window eventually since gamma_max is measured on the same grid the
            # threshold is compared against.
            if windows_next or windows_prev or adiabatic_is_good_enough(
                    max(gamma_next, gamma_prev)):
                report(len(windows_next), max(gamma_next, gamma_prev), iterations, True)
                return U_next, windows_next, True
        U_prev, windows_prev, ok_prev, gamma_prev = (U_next, windows_next, ok_next, gamma_next)

    report(len(windows_prev), gamma_prev, iterations, ok_prev)
    return U_prev, windows_prev, False


__all__ = [
    'adiabatic_propagator',
    'find_hidden_features',
    'find_resonance_candidates',
    'find_nonadiabatic_windows',
    'hybrid_propagator',
    'oscillation_sampling',
    # Documented as knobs -- each docstring carries the population it was
    # measured on -- and without this sphinx-autoapi does not document them,
    # which left every cross-reference to them rendering as dead text.
    'GAMMA_TO_ERROR',
    'RESOLUTION_RATIO',
    'LOCAL_JUMP_RATIO',
    'N_LOCAL_CONFIRM',
    'MAX_LOCAL_CONFIRMATIONS',
    'HIDDEN_FEATURE_CONCENTRATION',
    'N_HIDDEN_FEATURE_SUBDIVISION',
    'THRESHOLD0_PROVENANCE',
]
