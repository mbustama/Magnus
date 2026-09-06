Accuracy and diagnostics
==========================

.. contents::
   :local:
   :depth: 2

What you actually asked for when you passed ``rtol``, what each safeguard can
and cannot catch, and what every warning means.

Accuracy
----------

.. _what-rtol-atol-control:

What ``rtol`` and ``atol`` actually control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

They are a **stopping criterion, not an accuracy guarantee**, and the difference is worth
stating because the names invite the other reading.

The refinement ladder computes the probability matrix, grows ``n_slabs`` (and, for the
quadrature methods, ``n_tpts_per_slab``), recomputes, and stops when two successive levels
agree within ``atol + rtol*|P|``.  Nothing in that loop estimates the error of the answer
it returns.  A stepping ODE integrator's ``rtol`` is a different quantity: it bounds an
*estimated* local error per step, formed by comparing against an embedded lower-order
formula.  Magnus forms no such estimate; it infers convergence from agreement.

Usually that is conservative.  For a sequence converging as :math:`C n^{-p}` the
level-to-level gap overstates the error of the finer level, so an answer that stopped at
``rtol=1e-3`` is typically better than 1e-3.

**But agreement is evidence, not proof.**  On a sequence that is still jumping around, two
levels can agree by coincidence while both are far from the truth: measured on a sawtooth
density, the 3- and 4-slab levels agreed and the returned answer was wrong by **0.855** in
probability.  ``strict_convergence`` requires two *consecutive* agreements for that reason.

It is fair to ask why the gap is not converted into an error estimate by Richardson
extrapolation -- for refinement ratio :math:`r` and order :math:`p`, the finer level's
error is :math:`\text{gap}/(r^p - 1)` -- which is what the sibling NuOscProbExact does.
The answer is that the required :math:`p` is not available.  Fitting the observed order on
the Earth chord against a 4096-slab reference gives :math:`p` = 3.84, 5.62 and 4.06 at 1, 2
and 10 GeV for ``magnus_exp_order=4`` (nominal 4), but **1.15, 2.66 and 1.59** for
``magnus_exp_order=2`` (nominal 2) -- scattered by more than a factor of two, with one
sequence not even monotone.  On solar configurations under ``strategy='magnus'`` the error
sequence is frankly non-monotone at the slab counts the ladder visits, so no power law
holds at all.  Assuming :math:`p` equals the requested Magnus order would divide the gap by
too large a denominator wherever the true order is lower, reporting an error *smaller* than
the truth -- the dangerous direction, and the same shape as the false-certification bug
already on record in ``adiabatic.hybrid_propagator``.


**The oracle discipline.** ``solve_ivp``/DOP853 at ``rtol=1e-12, atol=1e-14`` is the only
accuracy oracle, and its convergence is verified per configuration by tightening to
``rtol=1e-13`` and confirming the movement is far below the error being quoted. Where the
profile makes ``scipy.linalg.expm`` exact -- a constant or declared-piecewise-constant ``H``
-- ``expm`` is used instead, because it is not an approximation at all. **No Magnus path is
ever scored against another**: that is the mistake this whole page's robustness work exists
to avoid, and a cross-check between two paths is reported as agreement, never as accuracy.

Measured distributions, against those oracles:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Population
     - Median error
     - p90
     - Silent misses
   * - 145 random smooth profiles, d ∈ {2,3,4,5}, N ∈ {1,3,12,30,80}
     - 6.08e-08
     - 7.14e-04
     - 6 (4.1 %)
   * - the same, restricted to N ≥ 30 (cumulative scan)
     - ~8e-09
     - --
     - **0**
   * - 150 random piecewise-constant profiles, edges declared
     - 1.34e-12
     - 1.40e-11
     - **0**
   * - the same, edges **not** declared
     - 7.76e-04
     - 2.96e-03
     - 2
   * - 164 Earth/solar configurations (not adversarial)
     - --
     - --
     - **0**

A *silent miss* is an answer outside the requested tolerance with no warning of any kind.
It is the only failure mode that matters; an inaccurate answer that says so is the warnings'
job. Every remaining silent miss sits below the seam -- single points and short scans
on random smooth profiles, overshooting a requested 1e-3 by a factor of one to three.

**Unitarity** is exact by construction (every engine composes unitary factors), and measured
on the package's own probability output it degrades only from ~3e-12 to 1.6e-11 across four
decades of N, at d = 2…5.

``tests/test_fuzz_statistics.py`` runs a CI-sized version of the fuzzing above and asserts on
the **distribution** -- silent-miss rate, median, worst case -- rather than on individual
cases, because a per-case assertion on random input is brittle and an aggregate one still
catches a regression that moves the distribution.


.. _safeguard-limits:

Robustness, and what each safeguard cannot do
-----------------------------------------------

Each safeguard below is stated with its limit, because the limits are what a user needs and
what a reviewer will not otherwise find.

**The probe-scale resolution test** (``magnus.adiabatic._profile_is_resolved``). Decides
whether ``H`` is continuous at the scale this package samples it on, by comparing how much of
the variation inside a probe interval falls in one half. *What it cannot do:* a jump smaller
than **1.33×** the local smooth variation is genuinely indistinguishable from steep smooth
behavior at that sampling density; see :data:`magnus.adiabatic.RESOLUTION_RATIO` for the
derivation of that factor from the threshold.

**γ-aware certification** (:data:`magnus.adiabatic.GAMMA_TO_ERROR`). When no non-adiabatic
window opens, successive refinements differ only in the transport grid, so they converge to
the same adiabatic limit and agree with each other whether or not that limit is right.
Certifying an empty window list therefore additionally requires γ itself to be small enough
for the requested tolerance. *What it cannot do:* the constant converts γ into an error
estimate good to about a factor of two, so certification near the bound is a closer call than
it looks.

**The patch budget** (``max_n_slabs = 32768`` in ``_local_evolution_operator``). A patch is
meant to be a short, local repair; one needing more slabs than a plain Magnus integration of
the whole trajectory means the non-adiabatic region is not narrow and the hybrid strategy has
no reason to exist for that request. Declining is the honest answer, and the general path is
70× faster there.

**Cross-method agreement** (:func:`magnus.oscprob.cross_check_strategies`). Runs whichever
engines apply and reports the pairwise spread. On the pre-fix package it reports the
disagreement on **seven of the eight** constructions where a method was silently wrong, each
at least four times the requested tolerance. *What it cannot do:* see the one below.

**The sampling report** (:func:`magnus.adiabatic.oscillation_sampling`). Answers a question no
engine asks itself: how coarsely does this request sample the oscillation it is computing?  A
solar trajectory is a few thousand oscillations long and a supernova ray tens of thousands, so a
scan of any ordinary size returns correct values that must not be read as a curve.  Surfaced as
``strategy_info['sampling']`` and **never warned about** -- the Nyquist criterion would fire on
44 of 45 realistic scan sizes, and a warning at that rate is noise.  Computed only when
``strategy_info`` was supplied, so the default path pays nothing.  See
:doc:`averaged_probability` for what to do when it says ``aliased``.

**The sub-probe feature scan** (:func:`magnus.adiabatic.find_hidden_features`). Looks at the
*profile* rather than at the answers, which is what lets it reach the one class no cross-check
can: within each interval of the refinement-ceiling grid, it compares the total variation a
denser grid sees inside that interval with the change its endpoints show, and reports the
largest excess as a fraction of the total. **Concentration, not size** -- an aliased sinusoid
hides variation in every interval, a narrow bump hides all of it in one. Measured at **0 false
positives over 67 smooth and resolvable profiles**, detecting 68-90 % of features in the
unresolvable band, for 0.37 ms once per call. *What it cannot do:* detection falls to ~0.73 for
features far below the dense sampling, and it **reports rather than cures** -- it names the
position and the exact ``t_breakpoints`` to pass, verified end to end -- warn, pass the
printed edges back, re-run: 3.0e-02 to 1.0e-04 on the width-3e-5 calibration case.

**The second irreducible limit: broadband roughness.** The sub-probe scan is a
*concentration* statistic, and that is exactly what makes it blind to structure spread over
every scale rather than piled into one place. Measured on Kolmogorov density fluctuations built
the way the supernova literature builds them -- a :math:`k^{-5/3}` spectrum with a 40-50 dB
dynamic range, so there is power below every grid this package lays down -- the finest reachable
grid sees only a third of the profile's total variation, and **neither structural test notices**:
:func:`magnus.adiabatic.find_hidden_features` returns a concentration of 0.002 against its 0.30
threshold, and ``_profile_is_resolved`` declares the profile resolved on 144 of 144
configurations. A power law spreads its sub-grid variation evenly over all 6400 reference
intervals, so each carries about :math:`1/6400` of it however large the total is; no threshold
reaches that, in the same way that no threshold reaches a feature which was never sampled.

What saves the answer is unrelated machinery: the errors such a profile produces (up to
1.4e-02 instantaneous at 45 MeV) are caught by the **convergence** checks, which watch the
refinement ladder rather than the profile. So the outcome is correct -- the caller is warned --
but by accident of mechanism rather than because anything recognized the profile. If you are
propagating through a turbulent or noisy medium, treat the structural diagnostics as silent by
construction and rely on the tolerance machinery, or supply ``t_breakpoints`` yourself. A cheap
statistic that *would* see this is described in ``docs/dev/FINDINGS_ROBUSTNESS_PROGRAMME.md``
§13.14; it is not shipped because it would need its own false-positive measurement first, and
the errors it would flag are already reported.

**The one irreducible limit: a feature narrower than the probe spacing.** A Gaussian
resonance of width :math:`10^{-5}` of the trajectory is not sampled by the probe grid
(spacing :math:`5\times10^{-3}`), nor by its refinement ceiling
(:math:`1.6\times10^{-4}`), nor by the cumulative scan's grid. Every engine reports a smooth
profile, small γ, and a resolved Hamiltonian -- correctly, given what any of them can see --
and all of them are wrong together by **2.9e-02 against a requested 1e-3**. Because they are
wrong *together*, the cross-check sees nothing either: it detects disagreement, so it finds a
wrong engine exactly when some other engine got it right.

The cure is caller-supplied ``t_breakpoints`` at the feature, and it is verified: with
edges placed by hand at the feature's own width the same case goes to 8.8e-04 at a single
point and 8.9e-04 over a 60-point scan, and with the set the warning itself prints (it
localizes the feature by re-sampling the flagged interval) to 1.0e-04. This is a property
of any fixed grid, not of any particular test, and no detector that pretends otherwise would
be honest. What *has* changed is that the condition is now usually **detected and reported**
rather than silent -- see the feature scan above.

**The scan is sized to the request.** It runs once per call whatever the point count, so its
share of the work falls as the request grows: 8 sub-steps (0.37 ms) for a single point, 32
(2.85 ms) for a scan of sixteen or more, holding it under about 7 % of the call at every size.
A single point keeps the cheapest scan by design -- the extra reach that finer sampling buys is
at widths of :math:`3\times10^{-6}` of the trajectory and below, narrower than anything
physically plausible in a density profile.

**A cross-check cannot close the rest, and this was measured rather than assumed.** Having
``strategy='auto'`` verify its own window-free results against the general Magnus ladder below
the seam (then at N = 25) was built, measured and removed: on 200 random smooth profiles the ladder agreed
with all 25 window-free results, and when :data:`magnus.adiabatic.GAMMA_TO_ERROR` was
deliberately mis-calibrated by 2x the check still fired zero times while three answers went
genuinely wrong. **What is left in that band is not engines disagreeing -- it is engines being
wrong together**, which a cross-check cannot see by construction. See
``docs/dev/FINDINGS_ROBUSTNESS_PROGRAMME.md`` §11.2.


.. _warning-catalogue:

Warnings: what each one means and what to do about it
-------------------------------------------------------

A warning here is an **instruction, not a disclaimer**. Each one below is held to four
things, in this order: what was detected, what it means for the answer (including *by how
much*, where the code knows), what to change, and when it is genuinely safe to ignore.

.. list-table::
   :header-rows: 1
   :widths: 24 30 22 24

   * - Warning
     - Condition
     - Is the answer affected?
     - What to change
   * - :class:`magnus.magnus.ScalarHamiltonianWarning`
     - ``H_func`` accepts only one position at a time.
     - No -- output is bit-identical.
     - Make ``H_func`` array-capable (``VCC[..., None, None]*e00``). Measured 4.6× faster.
   * - :class:`magnus.matter.DensityUnitWarning` (over-declared)
     - A density declared in g cm⁻³ is denser than a neutron star.
     - Yes -- catastrophically. The potential is inflated by ~18 orders; the tell is
       :math:`P_{ee} = 1`.
     - The density is already in natural units: leave
       ``density_matter_is_in_g_per_cm3`` at False.
   * - :class:`magnus.matter.DensityUnitWarning` (under-declared)
     - A density left in natural units is far too small to be one -- anything physical is
       4.3e18 or more, since that is what one g cm⁻³ becomes.
     - Yes, and this is the dangerous direction. The potential comes out ~19 orders too
       small, i.e. zero, so the call returns **exactly the vacuum probability** -- which
       looks like an ordinary answer rather than a missing one, and there is no tell in
       the numbers at all.
     - Pass ``density_matter_is_in_g_per_cm3=True``, or convert yourself (multiply by
       ``gd.UNIT_G_PER_CM3``).
   * - :class:`magnus.oscprob.UnmarkedDiscontinuityWarning`
     - The Hamiltonian is discontinuous at the grid scale and no ``t_breakpoints`` were
       given.
     - Yes, and refinement cannot help -- a straddling slab only gets narrower.
     - ``t_breakpoints`` at the jumps. Measured: median 7.8e-04 → 1.3e-12.
   * - :class:`magnus.oscprob.PhaseAveragingWarning`
     - ``average=True`` where the oscillation has not averaged.
     - The matrix is valid; the *question* does not apply there.
     - Use ``average=False``; the s.e.m. is reported.
   * - :class:`magnus.magnus.MagnusHighOrderCostWarning`
     - ``magnus_exp_order`` above 6 on ``'trapezoid'``/``'simpson'``.
     - No -- it is a cost trade, not an error.
     - Usually narrower slabs at order 4 or 6 instead.
   * - :class:`magnus.oscprob.ToleranceNotAchievedWarning`
     - A refinement cap was reached with the last two levels still disagreeing.
     - Unverified. The message reports **how far** from converged it stopped, as a multiple
       of the tolerance.
     - Raise the named cap; or loosen ``rtol``/``atol``; or add ``t_breakpoints``.
   * - :class:`magnus.oscprob.HybridCertificationWarning`
     - ``strategy='hybrid'`` was forced and a point did not self-certify.
     - **Unverified, which is not the same as wrong.** The result is still exactly unitary.
     - ``strategy='auto'`` (falls back automatically); or ``t_breakpoints`` at known
       structure; or a looser tolerance.
   * - :class:`magnus.oscprob.HiddenFeatureWarning`
     - The profile has structure too narrow for **any** grid here to sample.
     - Possibly wrong, and no strategy or tolerance helps -- every engine misses it together.
     - ``t_breakpoints`` at the position named in the message. A partial cure.
   * - :class:`magnus.magnus.MagnusConvergenceWarning`
     - :math:`\lVert\Omega\rVert_2 \geq \pi` on some slab.
     - **Unknown.** This reports a slab width, not an error.
     - Narrower slabs (smaller ``rtol``/``atol``, larger ``n_slabs``); ``t_breakpoints`` at
       any jump. Raising the order does not help.

**Measured false-positive rates** (``docs/dev/adversarial_batteries/warn_fp.py``, 168
configurations across the profile families this package serves, d = 2-5, scored against
``solve_ivp`` or -- for piecewise-constant profiles, where it is exact -- ``expm``):

.. list-table::
   :header-rows: 1
   :widths: 40 12 12 12 24

   * - Warning
     - Fired
     - TP
     - FP
     - FP rate
   * - :class:`magnus.magnus.MagnusConvergenceWarning`
     - 70
     - 17
     - 53
     - **76 %**
   * - :class:`magnus.oscprob.UnmarkedDiscontinuityWarning`
     - 56
     - 23
     - 33
     - 59 %
   * - :class:`magnus.oscprob.ToleranceNotAchievedWarning`
     - 37
     - 16
     - 21
     - 57 %

Silent misses across that whole population: **2 of 168 (1.2 %)**.

Read the 59 % carefully. ``UnmarkedDiscontinuityWarning`` reports a *condition about the input*,
not a prediction about the error, and on all 33 the condition was real -- there was an undeclared
discontinuity -- and the answer survived anyway. Declaring the edges would still have improved it
by orders of magnitude. A warning whose claim is true and whose advice is worth taking is not
made a false alarm by the answer surviving.

``MagnusConvergenceWarning``'s 76 % has a known cause and a **measured non-fix**. Of 66
single-point calls, some refinement level exceeded :math:`\pi` in 46, but **the level whose
answer was returned did so in only 7** -- so 85 % of its firings describe an intermediate grid
nobody receives. Keying it to the returned level therefore looks obviously right, and was
implemented. Re-measured over the same 168 configurations it made the warning *worse*: firings
fell 70 to 53, but **true positives fell 17 to 4** while false positives fell only 53 to 49.
"The ladder started far from convergence" predicts a bad answer better than "the final grid is
coarse" does. Nothing became silent either way (2 of 168 in both), because the cases it stopped
flagging are covered by :class:`magnus.oscprob.ToleranceNotAchievedWarning`. Reverted; the
mechanism and the numbers are kept in ``magnus._deferred_slab_norm``.

Two of these deserve their honesty spelled out rather than buried:

``MagnusConvergenceWarning`` **reports slab width, not accuracy.** The convergence bound it
checks is sufficient, not necessary, so exceeding it does not imply a wrong answer. It fires
on results accurate to 1.6e-06 and on results seven times outside a requested 1e-3, and
nothing available to it distinguishes the two. Until recently its message ended by telling
the reader that "if a target tolerance was requested … this warning can be ignored". That
clause was **false in exactly the cases where the warning matters**: on a sawtooth density
with ``rtol=atol=1e-3`` explicitly requested, under both ``strategy='auto'`` and
``strategy='magnus'``, the refinement ran and the answer was still 7.484e-03. It is gone.

``HybridCertificationWarning`` **means unverified, not wrong.** Every piece of the hybrid
propagator is unitary by construction, so the returned probabilities are a valid probability
matrix regardless; what is missing is the evidence that they are accurate to the tolerance
requested.
