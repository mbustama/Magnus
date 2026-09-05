Performance
=============

.. contents::
   :local:
   :depth: 2

Where the time goes, what was tried and rejected, and the population every
tuned constant was measured on.

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
:func:`magnus.oscprob.osc_prob_earth`, ``costhz = -0.9``, 2 GeV, against a vectorized
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
know it: the test needs the very evaluations the optimization skips.  A test on the slab
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
flat, because it loops over LAPACK internally instead of vectorizing over the stack.
:data:`magnus.magnus.EXPM_BACKEND` selects between that and the compiled kernels in
:mod:`magnus.expmkernels` -- ``'numba'`` means the Cayley-Hamilton kernel at dimensions 2
and 3 and the Jacobi eigensolver at 4 and 5, a distinction that did not exist when this
paragraph was written.  The Cayley-Hamilton kernel which applies to :math:`K` the polynomial interpolating
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
     - dimension 4 used ``eigh`` on both settings *at the time of this measurement*,
       which is what made it a control; it now uses the Jacobi kernel

**A 6.8× exponential is a 2.1× call, and the gap is Amdahl's law rather than a
disappointment.** The exponential is roughly a third of a slab pass, so removing six
sevenths of a third is about what the table shows.  Anyone quoting the 6.8× as a package
speed-up is quoting the wrong number.

A caution about the PREM row, because the first version of this table got it wrong.
:func:`magnus.earth.distance_traveled_inside_earth` returns **kilometers**, while every
``osc_prob`` baseline is in natural units, and passing the raw value does not raise: it
returns a converged, unitary answer for a chord a few meters long, on which the refinement
ladder trivially agrees with itself at every tolerance.  Measured that way the PREM speed-up
reads 1.45× rather than 2.11×, because a meter-long chord needs almost no slabs and so hardly
exercises the exponential at all.

**At N = 1 the exponential is no longer the thing to optimize.** ``eigh`` on one 3×3 costs
3.5 µs, and reaching it through ``_expm_stack`` costs 14.2 µs -- the
difference is the anti-Hermiticity test and the temporaries around it, which do not shrink
with the stack.  That fixed cost, not the exponential, is what caps the single-point rows
above.

**Dimensions 4 and 5 no longer keep ``eigh``, and an earlier version of this paragraph
said they always would.**  The reasoning was that a 4×4 or 5×5 Hermitian eigenproblem has
no practical closed form -- true, and beside the point, because the closed form was never
what the speed-up came from: ``eigh``'s fixed per-matrix LAPACK overhead (~2.3 µs on a
4×4, two thirds of a d = 4 call) was, and a batched Jacobi eigensolver that warm-starts
each matrix from its predecessor's eigenvectors removes it with no closed form at all --
2.6× on the exponential stage at 4ν and 1.7× at 5ν, measured against ``eigh`` plus
reconstruction on a 13k-slab chain.  Unlike the d ≤ 3 kernels it is iterative, so the
backend swap is not bit-identical; it is held instead to ``eigh``'s accuracy class, within
6.4× of it at every norm, clustering and degeneracy measured, under the same 10× bar that
admits the closed forms.  :func:`magnus.expmkernels.supports_dim` is still the only place
that decides this.

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

Two paths are now compiled beyond the exponential itself.  The separable energy scan folds
its slab operators in a numba kernel that keeps the Python loop's association but
accumulates each matrix element as a compiled scalar sum, where BLAS orders the same
arithmetic its own way.  A numba-less install was bit-identical there and may now differ
at the 1e-14 level -- worst observed 1.28e-14 across 16 scan configurations, with every
refinement decision unchanged.  The commutators of the Magnus schemes run in a second
such kernel, which fuses the two batched matmuls of ``X @ Y - Y @ X`` -- nearly all
gufunc dispatch at these matrix sizes -- into one pass over the stack.  Measured on the
two benchmark profiles, that cuts the marginal cost per slab of the order-4 scheme by
2.0-2.2x at three flavors and 3.1-3.4x at two; orders 6 and 8, which build three and six
commutators per slab, gain 2.5-3.0x and 1.9-2.8x.  At four and five flavors the matmuls
carry enough arithmetic to amortize their dispatch, so the gain settles at 1.1x; on the
cumulative-quadrature methods, whose time goes to the integrals rather than the
commutators, it disappears into the noise.  Probabilities move by at most 6.7e-14 across
36 configurations, every refinement decision and warning unchanged; without numba the
kernel falls back to the expression it replaced, bit-identical.


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

   * - Flavors
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

4ν and 5ν gained less here because they exponentiated through ``eigh``: the
Cayley-Hamilton kernel covers dimensions 2 and 3 only, and it still does.  That is no
longer the whole story -- those dimensions now go to the Jacobi eigensolver instead of
``eigh``, worth a further 1.8-1.9× at 4ν and 1.5-1.6× at 5ν end to end.  The table above
predates it.  In absolute terms a 3ν constant-density scan costs 1.10 µs per
energy, against NuOscProbExact's 1.44 µs batched and 13.25 µs looped; a single point is 33.8 µs
against its 19.9 µs, and **what remains is wrapper parameter resolution rather than
arithmetic** -- the exponential itself is under a tenth of it.  Part of that resolution
cost has since been removed: the accepted pass-through keyword names were rebuilt by
``inspect.signature`` on every public call, twice per call, and are now cached.  The
figures in this paragraph predate that change and were not re-measured for it.

Results are bit-identical to the per-point route on every flavor count and both neutrino signs.
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
       scans, growing with both flavor count and scan length, 1.06×-1.16× on cumulative
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
       0.3 maximizes detection (68–90 %) at five times the measured ceiling.
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
       behavior (a window opens below :math:`\gamma_\max`) but not usefully --
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
