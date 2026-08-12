Implementation Details
========================

.. contents::
   :local:
   :depth: 2

.. _the-engines:

The engines
-------------

Six independent engines can answer a request. None of them is a special case of another,
and each declines requests it cannot serve honestly.

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
   * - **Two-flavour interaction picture**
       (``_osc_prob_ip_exp_dispatch``)
     - A genuine exponential profile, built by
       :func:`magnus.matter.exp_density_profile`; exactly two flavours.
     - Single points and multi-energy scans at one baseline.
     - Non-exponential profiles, :math:`d > 2`, LIV, breakpoints, and whenever its own
       iteration fails to converge (typically near an MSW resonance).
   * - **Constant Hamiltonian**
       (``_osc_prob_scan_constant_h``)
     - ``V_CC`` does not depend on position, so neither does ``H``.
     - Vacuum and constant density, at any flavour count, for a single point or a scan,
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
   * - **Adiabatic + Magnus hybrid**
       (:func:`magnus.adiabatic.hybrid_propagator`)
     - ``H`` is smooth at the scale of a 200-point probe grid.
     - Any dimension, any smooth position-dependent profile, with a requested tolerance.
     - Breakpoints or slab edges supplied, a constant potential, no requested tolerance,
       a profile that fails the resolution test, or failure to self-certify.

A sixth reference, ``scipy.linalg.expm``, is not an engine but is used as an oracle by
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
* ``'interaction-picture'`` -- the two-flavour fast path. Same Magnus core, but the fast
  vacuum phase is factored out analytically first, so what it must resolve is a different
  function.
* ``'adiabatic'`` -- the hybrid strategy. A genuinely different method; its blind spots are
  the resonance detector's, not the quadrature's.
* ``'exact'`` -- ``expm``, independent of all of them.

Two engines in the same family can be wrong in the same way at the same time. Their
disagreement is informative; their agreement is not.


.. _dispatch-order:

Dispatch
----------

Every scenario wrapper (:func:`magnus.oscprob.osc_prob_matter_std_potential`,
:func:`magnus.oscprob.osc_prob_matter_nsi`, :func:`magnus.oscprob.osc_prob_liv`) tries the
engines in a fixed order, falling through on ``NotImplemented``:

.. mermaid::

   flowchart TD
       A["wrapper call"] --> AV{"average=True<br/>and H is position-independent?"}
       AV -- yes --> AVG["closed-form phase average<br/>(magnus.avgprob)"]
       AV -- no --> HY{"hybrid applies?<br/>smooth profile, tolerance requested,<br/>strategy != 'magnus',<br/>fewer than 25 scan points"}
       HY -- yes, and it certifies --> HYB["adiabatic + Magnus patch"]
       HY -- no --> IP{"exponential profile, d = 2?"}
       IP -- yes, and it converges --> IPX["interaction picture"]
       IP -- no --> SEP{"many energies, one baseline?"}
       SEP -- yes --> SEPE["energy-batched scan"]
       SEP -- no --> CUM{"baseline scan at one energy,<br/>2 or more points?"}
       CUM -- yes --> CUMS["cumulative scan"]
       CUM -- no --> GEN["general Magnus ladder"]

Two thresholds decide the seams, and both are constants with docstrings of their own:

* :data:`magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` = 8. Under
  ``strategy='auto'`` the hybrid strategy stands aside for a baseline scan of at least this
  many points, because the cumulative scan answers all of them from one traversal.  This was
  25; the constant's own docstring records why it moved, and why a later attempt to lower it
  to 1 was reverted.
* :data:`magnus.oscprob.CUMULATIVE_AUTO_MIN_POINTS` = 2. Below this there is no prefix to
  reuse.

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
    info['engine']      # 'hybrid', 'ip_exp', 'separable', 'cumulative', 'magnus', 'average'
    info['certified']   # for the hybrid strategy
    info['declined']    # [(engine, why it stood aside)]

This is the answer to "why did my result move?" and "why did this call get slow?", both of
which were previously unanswerable from outside the package.


Speed
-------

Measured by an alternating harness (``docs/dev/adversarial_batteries/timing.py``) that
interleaves the trees round-robin and carries two workloads the change cannot touch as
controls. **Ratios survive a loaded machine; absolute times do not** -- the controls came
back at exactly 1.00× while a fuzzer pegged a core, which is what makes the ratios below
usable at all.

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Workload
     - Versus the previous release
     - Note
   * - solar baseline scan, N = 400
     - **0.01×** (79× faster)
     - the cumulative scan; 8.1 s → 0.10 s
   * - single point, solar (hybrid)
     - **0.79×**
     - 21 % faster
   * - single point, 3ν solar
     - **0.83×**
     -
   * - single point, multi-resonance
     - **0.76×**
     -
   * - solar scan, N = 8 (hybrid)
     - **0.78×**
     -
   * - CONTROL: vacuum scan, N = 300
     - 1.00×
     - untouched by any change
   * - CONTROL: constant-density scan, N = 300
     - 1.00×
     - untouched by any change

Across 164 Earth and solar configurations spanning d = 2…5, standard/NSI/LIV, ν and ν̄, the
**median call is 2 ms** and the slowest is 0.90 s.

One lesson from that table is worth keeping, because it was nearly missed: an earlier
version of the γ sweep carried a comment claiming it *"reuses the eigendecomposition already
needed"* while in fact rebuilding it -- 600 extra Hamiltonian evaluations and a second
``eigh``, costing 1.4× at the entry point. **A claim in a comment about what code reuses is
not evidence that it reuses it.**


The palindrome, and what it is worth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A chord through a spherically symmetric Earth meets every radius twice, so its density
profile reads the same from either end.  :func:`magnus.magnus.magnus_expansion_multislab`
evaluates :math:`A` on the first half of such a slab chain and derives the rest by
reversal.  The saving is halved evaluations of the caller's Hamiltonian **and nothing
else** -- the matrix exponential is untouched, and so is the commutator algebra -- so it
is worth exactly what that Hamiltonian costs.  Measured through
:func:`magnus.oscprob.osc_prob_earth`, ``costhz = -0.9``, 2 GeV, against a vectorised
``H_func`` whose cost scales per position:

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Workload
     - Speed-up
     - Note
   * - single point, plain PREM
     - 0.905×
     - a density lookup is too cheap to be worth halving
   * - single point, expensive ``H_func``
     - **1.41×-1.67×**
     -
   * - 12- and 40-energy scan, expensive ``H_func``
     - **1.56×-1.64×**
     - falls to the general ladder, so the mirror applies
   * - energy scan, standard PREM
     - 1.00×
     - answered by the separable engine; see below

The ceiling is 1.67× rather than 2× because the refinement ladder and the unpaired middle
slab of an odd chain cut Hamiltonian evaluations from 159 positions to 93, not quite in
half.

**A standard PREM energy scan gains nothing, and that is correct rather than a gap.**  It
is answered by the separable engine, which already evaluates the profile once and shares
it across every energy -- the same saving, taken earlier and more completely.  Measured,
that engine spends a fraction :math:`f` = 0.001-0.026 of its time in the profile, which
caps any possible mirror gain at 1.001×-1.013×.

**Symmetry is declared, never detected.**  It cannot be detected where it would pay to
know it: the test needs the very evaluations the optimisation skips.  A test on the slab
*widths* is not a substitute -- a monotonic, solar-like profile on a uniform grid has
perfectly palindromic widths, and mirroring it is wrong by 3.3e-01.  So the declaration is
made by the Earth entry points, where a chord meeting every radius twice is geometry
rather than a claim, and it travels as the *interval* it holds over rather than as a flag:
a chord is symmetric over its full length and over no shorter prefix, so a request at a
shorter baseline fails the span check and takes the ordinary path with no extra
bookkeeping.

Set :data:`magnus.magnus.USE_PALINDROME` to ``False`` to evaluate every slab in full.  The
two routes agree to a few times 1e-15 rather than bitwise, because the mirrored slab's
nodes are reached as ``(L - b) + h*s`` on one route and ``a + h*s`` on the other -- two
floating-point expressions for the same real number.  On Earth single points that is worth
up to 8.6e-15 relative.


The matrix exponential, and which backend computes it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every slab ends in a matrix exponential, and ``np.linalg.eigh`` costs **about 1.27 µs per
3×3 whatever the stack size** -- measured 1.268 µs at N = 108 and 1.279 µs at N = 4096,
flat, because it loops over LAPACK internally instead of vectorising over the stack.
:data:`magnus.magnus.EXPM_BACKEND` selects between that and the compiled Cayley-Hamilton
kernel in :mod:`magnus.expmkernels`, which applies to :math:`K` the polynomial interpolating
:math:`\exp(-i\lambda)` on its spectrum -- no eigenvectors, and the eigenvalues in closed
form.

Interleaved round-robin, minima of many repetitions, with a control the change cannot
touch:

.. list-table:: The exponential alone, :math:`\exp(-iK)` for a stack of N matrices
   :header-rows: 1
   :widths: 10 12 20 20 20

   * - d
     - N
     - ``eigh``
     - ``numba``
     - Speed-up
   * - 3
     - 1
     - 14.2 µs
     - 7.4 µs
     - 1.9×
   * - 3
     - 108
     - 162.6 µs
     - 23.8 µs
     - **6.8×**
   * - 3
     - 4096
     - 6467 µs
     - 934 µs
     - **6.9×**
   * - 2
     - 108
     - 94.0 µs
     - 12.9 µs
     - **7.3×**
   * - 2
     - 1024
     - 716.9 µs
     - 54.5 µs
     - **13.2×**

.. list-table:: End to end, through ``osc_prob``
   :header-rows: 1
   :widths: 46 27 27

   * - Workload
     - Speed-up
     - Note
   * - 3ν PREM, 60-energy scan
     - **2.11×**
     - 9291 µs → 4409 µs (73.5 µs per energy)
   * - 3ν PREM chord, single point
     - 1.22×
     - dominated by the refinement ladder
   * - 3ν vacuum, single point
     - 1.11×
     - and see the constant-Hamiltonian engine below, which is the larger win here
   * - 3ν constant density, single point
     - 1.09×
     -
   * - CONTROL: 4ν vacuum
     - 1.00×
     - dimension 4 uses ``eigh`` on both settings

**A 6.8× exponential is a 2.1× call, and the gap is Amdahl's law rather than a
disappointment.** The exponential is roughly a third of a slab pass, so removing six
sevenths of a third is about what the table shows.  Anyone quoting the 6.8× as a package
speed-up is quoting the wrong number.

A caution about the PREM row, because the first version of this table got it wrong.
:func:`magnus.earth.distance_traveled_inside_earth` returns **kilometres**, while every
``osc_prob`` baseline is in natural units, and passing the raw value does not raise: it
returns a converged, unitary answer for a chord a few metres long, on which the refinement
ladder trivially agrees with itself at every tolerance.  Measured that way the PREM speed-up
reads 1.45× rather than 2.11×, because a metre-long chord needs almost no slabs and so hardly
exercises the exponential at all.

**At N = 1 the exponential is no longer the thing to optimise.** ``eigh`` on one 3×3 costs
3.5 µs, and reaching it through ``_expm_stack`` costs 14.2 µs -- the
difference is the anti-Hermiticity test and the temporaries around it, which do not shrink
with the stack.  That fixed cost, not the exponential, is what caps the single-point rows
above.

**Dimensions 4 and 5 keep ``eigh``, and always will.**  There is no practical closed form
for a 4×4 or 5×5 Hermitian eigenproblem, so 4ν and 5ν are correct and simply not
accelerated.  :func:`magnus.expmkernels.supports_dim` is the only place that decides this.

Neither backend is exactly unitary, and a previous version of ``_expm_stack``'s docstring
claimed the ``eigh`` one was.  It is not: :math:`U^\dagger U - I` measures 4e-16 for a
single 3×3 and 4e-15 for a stack of 4096, growing with stack size and never reaching zero.
Against a 40-digit reference the kernel is the same order or slightly better at every norm
from :math:`\lVert K \rVert` = 1 to 1e5 **on unclustered spectra**, and both degrade linearly
in that norm, which is the conditioning of the problem rather than a property of either route.
Probabilities sum to 1 to about 1e-15; they do not do so by construction.

That qualifier was missing from an earlier version of this page, and it mattered.  The closed
form was verified against random spectra at many norms, and separately at many eigenvalue
separations at norm ~1; where those two conditions hold *together* it reached 2.7e-07 against
``eigh``'s 3.0e-11, a factor of 7440, because :math:`\arccos` has infinite derivative at
:math:`u = \pm 1`.  Neither single-axis sweep visits that corner.  It is now closed by
:data:`magnus.expmkernels.SEV_TOL`, which hands such matrices to ``eigh``; the worst absolute
error over the whole separation-by-scale grid is 8.7e-14, and the fraction of matrices declined
on real work -- PREM chords, solar slab chains, constant density, NSI -- measures 0.00%.

Switching backend moves probabilities by at most 4.6e-15 across PREM chords, energy scans,
NSI resonances, constant density and vacuum -- except on a solar profile at
``strategy='magnus'``, which chains 33,575 slab exponentials and drifts 3.0e-12, within the
:math:`N\epsilon` = 7.4e-12 that an ordered product of that length allows.

numba is a required dependency, so ``'auto'`` reaches the compiled kernel on any
ordinary install.  It costs about 90 ms of ``import magnus``, and the first call to each
kernel pays a one-off ~0.7 s compile that is then cached to disk.

The ``'eigh'`` fallback is still there and still correct -- ``'auto'`` degrades to it if
the import fails for any reason, and nothing but speed changes, every result agreeing to
~1e-15.  numba was optional on exactly that argument, which held for the library and not
for its suite: a clean-room install of the published wheel produced twelve *failures*
rather than twelve skips, because ``tests/test_engines.py`` asks for
``expm_backend='numba'`` by name.  The trade is that numba lags new interpreters, so a
Python release it has no wheel for now makes the package uninstallable rather than merely
slower.


A constant Hamiltonian needs no ladder at all
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When :math:`V_\text{CC}` does not vary with position, neither does :math:`H`, and the Magnus
series **terminates at its first term**: :math:`\Omega_1 = -iH\Delta`, and every higher
:math:`\Omega_k` is a nested commutator of :math:`H` with itself, hence zero.  So
:math:`U = \exp(-iH\Delta)` is not an approximation to be refined but the exact answer, and an
entire energy scan is one stacked exponential over an ``(nE, d, d)`` array.

This case used to be turned away deliberately -- the separable dispatcher bailed out on a
non-callable potential, its docstring saying "a constant potential falls back to the generic
path" -- so the easiest Hamiltonian there is took the slowest route available: a 60-energy scan
made 18,000 ``osc_prob`` calls per 300 repetitions, each one rediscovering the same constancy.

.. list-table:: Against the per-point route it replaces (interleaved; control 1.00×)
   :header-rows: 1
   :widths: 16 22 22 20

   * - Flavours
     - Matter scan
     - Vacuum scan
     - Single point
   * - 2ν
     - **17.3×**
     - **24.7×**
     - 2.0×
   * - 3ν
     - **15.5×**
     - **18.9×**
     - 2.1×
   * - 4ν
     - 7.2×
     - 7.4×
     - 1.4×
   * - 5ν
     - 6.0×
     - 6.2×
     - 1.4×

4ν and 5ν gain less because they exponentiate through ``eigh``: the Cayley-Hamilton kernel
covers dimensions 2 and 3 only.  In absolute terms a 3ν constant-density scan costs 1.10 µs per
energy, against NuOscProbExact's 1.44 µs batched and 13.25 µs looped; a single point is 33.8 µs
against its 19.9 µs, and **what remains is wrapper parameter resolution rather than
arithmetic** -- the exponential itself is under a tenth of it.

Results are bit-identical to the per-point route on every flavour count and both neutrino signs.
``n_slabs``, ``n_tpts_per_slab``, ``t_breakpoints`` and ``rtol``/``atol`` are accepted and
ignored, because they can only ask for a refinement of something already exact.

**PREM and exponential profiles are untouched** -- their potential varies with position, so they
keep ``separable``, ``ip_exp`` or ``hybrid``.  A constant-H engine that captured one would
propagate a whole chord with a single exponential of a single Hamiltonian: wrong by O(1) and
still perfectly unitary, which is why ``tests/test_engines.py`` asserts the engine *identity*
for PREM and the Sun rather than only comparing numbers.

Two traps this engine paid for, both recorded because neither was visible in the answer.
``h_matt`` meant different things on different branches -- two of the three dispatch call sites
had folded :math:`V_\text{CC}` into it already, and the engine multiplied by
:math:`V_\text{CC}` again, giving :math:`V_\text{CC}^2 \sim` 1e-25 instead of 1e-13: the matter
term all but vanished and, because a square has no sign, the neutrino and antineutrino answers
came back *bit-identical*.  And a new engine absent from ``_CROSS_CHECK_FORCING``'s forbid lists
answers before the payload the independent ``expm`` oracle is built from is ever recorded, which
silently removed the only non-Magnus reference from the cross-check.


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
behaviour at that sampling density; see :data:`magnus.adiabatic.RESOLUTION_RATIO` for the
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
position and the ``t_breakpoints`` to pass, which is a partial fix (measured 3.9e-03 to
8.5e-05), not a complete one.

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
but by accident of mechanism rather than because anything recognised the profile. If you are
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

The cure is caller-supplied ``t_breakpoints`` at the feature, and it is verified: the same
case goes to 8.8e-04 at a single point and 8.9e-04 over a 60-point scan. This is a property
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
     - ``magnus_exp_order`` above 6.
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


.. _how-constants-were-set:

How the constants were set
----------------------------

This is the section no reader can reconstruct from the code. Every calibration constant is
listed with its provenance -- or with an explicit statement that it has none, which is the
honest entry for a number that has simply always been what it is.

Measured
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Constant
     - Value
     - Population it was measured on
   * - :data:`magnus.adiabatic.GAMMA_TO_ERROR`
     - 0.85
     - 149 configurations: resonance width over a decade, d = 2…5, 5–80 MeV, 0.5–2 density
       scale heights, pure adiabatic operator scored against ``solve_ivp``. **Restricted to
       the small-γ rows the rule actually governs** -- reading the unrestricted maximum was
       the over-correction that made an earlier value wrong.
   * - :data:`magnus.adiabatic.RESOLUTION_RATIO`
     - 0.70
     - 192 smooth configurations (ceiling 0.602) against 15 random piecewise-constant ones
       (1.000), plus a deliberately weak jump 4.7× smaller than the steepest smooth step
       (0.773).
   * - :data:`magnus.oscprob.BATCH_WORKING_ENTRIES`
     - 65 536
     - Fifteen workloads on three batched engines, d = 2…5, scans of 60 to 20 000 points,
       swept over 1 / 4.2 / 12.6 / 67 / 268 MB.  1 MB won eight of the eleven memory-bound
       rows and was never worse than the previous 67 MB: **1.19×-1.38×** on Earth energy
       scans, growing with both flavour count and scan length, 1.06×-1.16× on cumulative
       baseline scans, flat within 2 % on short scans.  The interaction-picture engine is
       flat at 1.00× -- it is compute-bound, so the constant does not reach it.  Every row
       was **bit-identical at every budget**, tiles being independent and only
       concatenated, so this is a pure performance knob.  Measured on one machine (13 MB
       L3, 6.5 MB L2), and note the optimum sits *below* the last-level cache, so sizing
       to a detected cache would land on a worse value than this fixed constant does.
   * - ``_local_evolution_operator`` ``max_n_slabs``
     - 32 768
     - Legitimate patches converge at 800–12 800 slabs; a patch covering 88 % of a solar
       trajectory needs 102 400 and should decline. 32 768 sits in the factor-of-eight gap.
   * - :data:`magnus.oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`
     - 8
     - Cost/accuracy crossover measured over scan sizes and re-measured over 42 workloads,
       which moved it from 25; the cumulative scan is cheaper on median at every size and
       three to six orders more accurate on the ones it serves.
   * - :data:`magnus.oscprob.CUMULATIVE_N_ACC_SAFETY`
     - 4
     - The longest baseline sets the grid; shorter ones in the same scan would have chosen a
       denser one for themselves.
   * - :data:`magnus.adiabatic.LOCAL_JUMP_RATIO`
     - 0.5
     - 79 flagged intervals over 1440 smooth configurations (ceiling **0.087**) against 348 over
       432 piecewise ones (floor **1.000**). Swept over *sub-intervals*, the axis the original
       ``RESOLUTION_RATIO`` measurement did not have.
   * - ``find_resonance_candidates`` ``fd_step_frac``
     - 1e-6
     - Scored against the **analytic** :math:`dH/dl`. The optimum moves with the profile's
       shortest length scale (1e-5 solar, 1e-6 sinusoid, 1e-7 for a narrow bump), but anywhere
       in 1e-8…1e-5 the relative error stays below 3e-09 -- six orders below anything that
       could move a probability here. **The band, not the value, is what to preserve.**
   * - ``hybrid_propagator`` ``threshold0``
     - 0.1
     - See :data:`magnus.adiabatic.THRESHOLD0_PROVENANCE`. Accuracy identical at every value in
       16 of 18 rows at a fixed baseline, and a lower start up to **6.5×** cheaper -- but a
       tolerance-derived rule built on that evidence made an **energy scan 20× worse**
       (2.5e-05 → 4.95e-04) and was reverted.
   * - :data:`magnus.adiabatic.HIDDEN_FEATURE_CONCENTRATION`
     - 0.3
     - 67 smooth and resolvable profiles (ceiling **0.060**) against features in the
       unresolvable band (0.91–1.00). **0 false positives at every threshold from 0.2 to 0.6**;
       0.3 maximises detection (68–90 %) at five times the measured ceiling.
   * - :data:`magnus.adiabatic.N_HIDDEN_FEATURE_SUBDIVISION`
     - 8
     - Chosen on cost, not on the statistic (which is flat in it): 0.37 ms against 2.85 ms at
       32, where the arrays stop fitting in cache.
   * - ``n_probe0``, ``n_points0``, ``patch_atol``, ``n_slabs0``, ``growth_factor_n_slabs``,
       ``min_n_tpts_per_slab``
     - 200, 201, 1e-7, 400, 1.5, 2
     - Swept across **18 workloads spanning single points, baseline scans and energy scans** ×
       3 profile families × d = 2, 3. The worst error is **4.49e-04 at essentially every value
       of every one of them**: these set where a doubling ladder starts, and the ladder reaches
       the same place regardless. ``patch_atol`` at 1e-9 is the one exception and is not really
       about this constant -- see :func:`magnus.adiabatic.hybrid_propagator`.
   * - ``min_threshold``
     - 1e-6
     - Identical at every value over 18 ordinary workloads, because the ladder stops long
       before the floor. **The regime it governs was then constructed rather than assumed**:
       the floor is reached only when :math:`\gamma_\max` is below it *and* the tolerance is
       tighter than ``GAMMA_TO_ERROR`` :math:`\times \gamma_\max`. There it does change
       behaviour (a window opens below :math:`\gamma_\max`) but not usefully --
       ``certified=False`` at every value, error three orders inside tolerance either way, and
       the window costs 2.4× the time.

**Why ``threshold0`` was measured, changed, and changed back.** The fixed-baseline sweep said a
tolerance-derived rule was safe and cheaper. It was built, and the package's bit-identity
workloads — which include an energy scan the sweep did not — said otherwise: one row 13711×
better, another 20× worse. **A population that does not contain the workload you are about to
change is not evidence about it**, which is the same mistake that made ``GAMMA_TO_ERROR`` wrong
twice, committed again while explicitly trying to avoid it. The measurement is kept; the default
is not changed.

Not measured
~~~~~~~~~~~~~~

The following carry no provenance beyond "it has always been that". They are listed rather
than quietly left out, because an unaudited constant that nobody has written down is
indistinguishable from an audited one.

``max_n_probe`` (6400), ``max_n_points`` (12864) and ``max_iters`` (12) in
:mod:`magnus.adiabatic`; ``max_num_loops`` (50) in :mod:`magnus.oscprob`.

All four are **cost ceilings rather than calibrations**: they bound work, and reaching one is
reported by :class:`magnus.oscprob.ToleranceNotAchievedWarning` rather than absorbed. So
"unmeasured" means something milder for them than for a threshold that silently decides an
outcome. Every constant that *does* silently decide an outcome now appears in the table above.


Reproducing any of this
-------------------------

Every measurement on this page comes from a script under
``docs/dev/adversarial_batteries/``, and their outputs are deliberately **not** committed, so
re-running is the only way to get them:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Script
     - What it measures
   * - ``crosscheck_acceptance.py``
     - Whether a cross-check between engines would have caught the known silent misses.
       Runs against either tree via ``PYTHONPATH``.
   * - ``invariants.py``
     - The oracle-free invariants, swept over a profile matrix.
   * - ``warn_fp.py``
     - Every warning's true- and false-positive rate.
   * - ``constants_audit.py``, ``constants_audit2.py``
     - Provenance for the calibration constants above; the second sweeps 18 workloads spanning
       points, baseline scans and energy scans.
   * - ``resolution_fp.py``
     - The resolution test's false-positive rate, swept over sub-intervals.
   * - ``weak_band.py``, ``crosscheck_benefit.py``
     - Where the hybrid path's self-certification is weak, and whether a default-path
       cross-check would earn its cost. It does not; see the robustness section.
   * - ``battery2.py`` … ``battery10_coverage.py``
     - The original adversarial batteries; see
       ``docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md``.

See also :doc:`adiabatic_strategy` for the hybrid strategy's derivation and validation, and
``docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md`` for the adversarial validation these
safeguards came out of.
