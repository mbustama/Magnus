Methodology
============

.. contents::
   :local:
   :depth: 2


This page documents the numerical machinery behind Magνs: the Magnus
expansion itself, the two families of integrators, the guarantees they
carry, the adaptive-refinement and performance engineering built around
them, and the evidence used to validate all of it.  See
:ref:`when-is-magnus-a-win` on the front page for the short version.

The Magnus expansion
-----------------------

Neutrino flavor evolution is governed by the Schrödinger-like equation

.. math::

   i \frac{d}{dl}\, |\psi(l)\rangle = H(l)\, |\psi(l)\rangle ,

with :math:`H(l)` the (possibly position-dependent) flavor Hamiltonian.
The evolution operator :math:`U(l_1, l_0)` satisfies the same equation with
:math:`U(l_0, l_0) = \mathbb{1}`.  When :math:`H` does not commute with
itself at different positions, :math:`U` is *not* simply
:math:`\exp\!\left[-i\int_{l_0}^{l_1} H(l)\, dl\right]`.

The Magnus expansion instead writes :math:`U(l_1, l_0) = \exp[\Omega(l_1)]`
exactly, where :math:`\Omega = \sum_k \Omega_k` is built order by order from
nested commutators of :math:`A(l) \equiv -i H(l)` :cite:p:`Blanes2009`:

.. math::

   \Omega_1(l) &= \int_{l_0}^{l} A(s)\, ds \\
   \Omega_n(l) &= \sum_{j=1}^{n-1} \frac{B_j}{j!} \int_{l_0}^{l} S_n^{(j)}(s)\, ds ,

with :math:`B_j` the Bernoulli numbers (:math:`B_1 = -1/2` convention) and
:math:`S_n^{(j)}` sums of :math:`j`-fold nested commutators of the
lower-order terms with :math:`A`.  Magνs implements this recursion through
:math:`n = 10` (odd Bernoulli numbers :math:`B_3 = B_5 = 0` vanish
identically, so only even-index commutator groups appear beyond
:math:`\Omega_3`); the coefficients and every term were verified
independently, term by term, against this recursion (see
:ref:`validation`).  Orders 1 to 6 are written out inline; beyond that the
terms are generated from the recursion, since their number roughly doubles
per order.  :doc:`expansion_terms` derives them symbolically at any order,
which is what the verification checks against.

**Truncating the series is exact for the group, not just approximate for
the answer.**  Whatever order the sum stops at, :math:`\Omega` remains
anti-Hermitian (since each :math:`\Omega_k` is a real combination of nested
commutators of anti-Hermitian matrices), so :math:`\exp(\Omega)` is
*exactly* unitary — probabilities computed from it are non-negative and sum
to one to machine precision, regardless of the truncation order or the
quadrature accuracy.  This is the central practical advantage over direct
ODE integration, whose iterates only approximately preserve unitarity (see
:ref:`accumulated-phase`).

The series converges absolutely whenever
:math:`\int_{l_0}^{l_1} \lVert A(l)\rVert_2\, dl < \pi` over the interval
in question.  Magνs partitions the trajectory into a chain of slabs and
evaluates the expansion independently in each one; a large accumulated
phase (a long baseline, a strong potential, or both) is handled by adding
more, narrower slabs rather than by raising the expansion order.  Since
:math:`\lVert\Omega_1\rVert_2` is a necessary (if not sufficient) proxy for
this criterion, and its value is obtained for free from the eigenvalues
already computed for the matrix exponential (see below), Magνs checks it
automatically and emits ``MagnusConvergenceWarning`` if a slab is
comfortably outside the guaranteed regime.

Two integration methods
--------------------------

Evaluating the nested integrals above requires sampling :math:`A(l)` inside
each slab.  Magνs offers two families, selected via
``integration_method``, which defaults to ``'gl'``:

**Gauss-Legendre commutator-free integrators (** ``'gl'`` **, the default).**
Following :cite:t:`Blanes2000`, orders 2, 4, and 6 can be reached from only
1, 2, or 3 evaluations of :math:`A` per slab, at the Gauss-Legendre nodes,
with no cumulative quadrature and no separate commutator bookkeeping:

.. math::

   \Omega^{(2)} &= h\, A_1 \\
   \Omega^{(4)} &= \frac{h}{2}(A_1 + A_2) + \frac{\sqrt{3}}{12} h^2\, [A_2, A_1] \\
   \Omega^{(6)} &= \ldots \quad \text{(three-node scheme; see the reference)}

with :math:`h` the slab width and :math:`A_i` the Hamiltonian sampled at
the corresponding node.  Because the quadrature order is matched exactly
to the truncation order, this method needs far fewer Hamiltonian
evaluations for the same accuracy -- it is simultaneously the fastest and
the most accurate choice whenever the Hamiltonian is smooth within a slab,
which is why it is the default.  Layer-aligned slabs (below) make that the
common case even across the Earth.

Because ``'gl'`` uses a fixed 1, 2, or 3 nodes per slab, ``n_tpts_per_slab``
plays no role for it: accuracy is controlled by the slab count alone, and the
adaptive refinement below grows only ``n_slabs``.  The physics-informed
starting slab count is likewise applied only for ``'gl'``, since for the
quadrature methods accuracy is governed jointly by ``n_slabs`` and
``n_tpts_per_slab``, and seeding only the slab count unbalances that ladder.

**Cumulative quadrature (** ``'trapezoid'`` **,** ``'simpson'`` **).**
Sample :math:`A` on a uniform grid of ``n_tpts_per_slab`` points and
integrate with cumulative trapezoid or Simpson's rule.  Slower for the same
accuracy on a smooth profile, but fully general, and so the safer choice if
:math:`A(l)` has a kink or a discontinuity *inside* a slab, where
Gauss-Legendre loses its order advantage.  The quadrature error
(:math:`O(h^2)` or :math:`O(h^4)` in the grid spacing :math:`h`) can dominate
the Magnus truncation error at high orders unless ``n_tpts_per_slab`` grows
accordingly.

Exact unitarity from the eigendecomposition
------------------------------------------------

Since :math:`\Omega` is anti-Hermitian, Magνs computes
:math:`\exp(\Omega)` from the eigendecomposition of the Hermitian matrix
:math:`K = i\Omega`:

.. math::

   \exp(\Omega) = V\, \mathrm{diag}\!\left(e^{-i\lambda}\right)\, V^\dagger ,
   \qquad K = V\, \mathrm{diag}(\lambda)\, V^\dagger .

This is both faster than a general (Padé-based) matrix exponential for
stacks of small matrices, and exactly unitary by construction — no residual
non-unitarity to track.  A general (non-anti-Hermitian) fallback based on
``scipy.linalg.expm`` remains available for exotic, non-physical uses of
the underlying :func:`magnus.magnus.magnus_expansion` engine.

Time-ordering
----------------

A neutrino traversing a chain of slabs accumulates the evolution operator
as a time-ordered product, with the *last* slab as the leftmost factor:

.. math::

   U_\mathrm{tot} = U_N \cdots U_2\, U_1 .

This matters physically whenever the Hamiltonians of different slabs do
not commute — e.g., an asymmetric density profile together with a nonzero
CP-violating phase — and is exercised directly in the test suite with an
exact two-constant-slab check (:math:`\exp(-iH_B L_2)\exp(-iH_A L_1)` from
matrix arithmetic alone, no quadrature).

.. _accumulated-phase:

Adaptive refinement and slab placement
-----------------------------------------

By default, ``osc_prob`` and its wrappers refine the number of slabs (and,
for the quadrature methods, the number of points per slab) until the
probability matrix stops changing within a requested tolerance
(``rtol``, ``atol``), doubling as the standard heuristic for an a
posteriori error estimate.  Three refinements make this efficient in
practice:

* **Physics-informed starting slab count.**  Rather than always starting
  from one slab, the refinement is seeded from an estimate of the
  accumulated (traceless) phase :math:`\lVert\Omega_1\rVert_2` over the
  whole trajectory, aiming for roughly :math:`2\pi` radians of phase per
  slab — enough for the Gauss-Legendre method to already be close to
  converged at the first attempt.
* **Warm starts across scan points.**  When computing many points (an
  energy scan, an oscillogram), each point's refinement is seeded from the
  previous point's converged slab count and point count, rather than
  reclimbing the same geometric ladder from scratch.
* **Slab edges aligned with density discontinuities.**  The PREM profile
  used for the Earth is piecewise-smooth, with density discontinuities at
  the boundaries between its ten shells :cite:p:`Dziewonski1981`.  A
  slab that straddles one of these boundaries locally degrades the
  quadrature to low order no matter how high ``magnus_exp_order`` is set.
  The Earth wrappers compute the exact chord positions where the
  trajectory crosses a PREM layer boundary (a closed-form quadratic in the
  zenith angle) and insert them as mandatory slab edges at every
  refinement level.
* **A caller-supplied floor.**  Passing ``n_slabs`` together with a
  tolerance sets a lower bound on the ladder: refinement starts at
  ``max(min_n_slabs, n_slabs)`` and only ever climbs from there (clipped at
  ``max_n_slabs``).  With the default ``n_slabs = 1`` the floor is inactive.

.. _refinement-blind-spot:

.. warning::

   The phase estimate that seeds the ladder is an *integral* of the
   Hamiltonian along the trajectory, and an integral is blind to structure
   that averages out.  A profile that oscillates rapidly about its mean --
   a castle wall, a periodically layered medium -- can accumulate very
   little net phase while still demanding many slabs to resolve, and will
   then be seeded with far too few.  The successive-iterate test is no
   protection here: refinements that all fail to see the profile can agree
   with each other while disagreeing with the truth, and a tighter ``rtol``
   only compares two answers that are both wrong.  Tightening the tolerance
   is the wrong lever; resolving the profile is the right one.

   If you know your profile's feature scale, say so, in either of two ways.
   Pass ``n_slabs`` (a floor, per the bullet above) so the ladder cannot
   start below it.  Better, where the features are discontinuities at known
   positions, pass those positions as ``t_breakpoints``: they become
   mandatory slab edges, which both resolves the profile and restores the
   quadrature's nominal order, and so costs less than the equivalent number
   of uniform slabs.  On a 50-wall castle-wall profile the two together
   reduce the worst-case error over a baseline scan from 0.855 to 1.9e-3,
   while running faster than the under-resolved version did.

The slab cap itself is method-aware.  ``max_n_slabs`` defaults to None,
meaning "use the cap appropriate to ``integration_method``": 20000 for
``'gl'`` and 2000 for the cumulative-quadrature methods (see
``magnus.oscprob.MAX_N_SLABS_DEFAULT``; an explicit value is always used as
given).  A single cap cannot serve both families, because their cost per
slab differs by more than an order of magnitude -- ``'gl'`` evaluates the
Hamiltonian 1 to 3 times per slab, the quadrature methods
``n_tpts_per_slab`` times.  With a shared cap of 2000, ``'gl'`` hit the
ceiling on problems it could resolve comfortably (eV-scale sterile
splittings over an Earth-crossing baseline need about 8,600 slabs) and
reported that it could not verify convergence, on answers that were in fact
far more accurate than the quadrature methods reached within the same cap.
Even at 20000 slabs, ``'gl'`` is the cheaper worst case: 40,000-60,000
Hamiltonian evaluations, against the ~200,000 that 2000 quadrature slabs at
100 points per slab already permit.

If a refinement cap (``max_n_slabs``, ``max_n_tpts_per_slab``,
``max_num_loops``) is reached before the tolerance is met, ``osc_prob``
returns its best available estimate but raises
``ToleranceNotAchievedWarning`` unconditionally (regardless of the
``verbose`` setting) — the returned probabilities remain exactly unitary,
so they can look entirely plausible while still being inaccurate.  This is
the practical manifestation of the convergence criterion above: it is the
expected behavior for extreme accumulated phases, such as low-energy solar
neutrinos traversing most of the Sun, where an adiabatic treatment is the
more natural tool — see :doc:`adiabatic_strategy` for the
``strategy='hybrid'``/``'auto'`` alternative that automates exactly this,
built directly on top of the machinery described on this page (its local
patches call the same :func:`magnus.magnus.magnus_expansion_multislab`
kernel).

Choosing the expansion order
-------------------------------

``magnus_exp_order`` defaults to 4, and for the tolerances most calculations
ask for that is the right choice.  The adaptive refinement already turns a
higher order into fewer slabs on its own, so the order and the requested
tolerance interact: raising the order pays only once the tolerance is tight
enough to make the extra work per slab worthwhile.

Measured wall time relative to order 4 on the same problem (greater than 1
means order 6 is faster):

.. list-table::
   :header-rows: 1
   :widths: 16 17 17 15 17 18

   * - Tolerance
     - Earth 3ν, 1 GeV
     - Earth 3ν, 10 GeV
     - Earth 5ν
     - Exponential density
     - 200-energy scan
   * - :math:`10^{-4}`
     - 0.89
     - 0.91
     - 0.89
     - 1.04
     - 0.70
   * - :math:`10^{-6}`
     - 0.98
     - 1.03
     - 1.07
     - 1.02
     - 1.05
   * - :math:`10^{-8}`
     - 1.25
     - 1.42
     - 1.08
     - 1.51
     - 1.93

So: leave the order alone for everyday work, and raise it to 6 if you are
asking for :math:`10^{-7}` or tighter, where it runs up to twice as fast.
Dropping to order 2 is almost never worthwhile -- at :math:`10^{-8}` on the
Earth cases it needs thousands of slabs where order 6 needs about a hundred,
and runs roughly twenty times slower.

Beyond order 6 the terms are generated rather than written out, the count
roughly doubles per order, and ``'gl'`` has no scheme at all (see
:doc:`expansion_terms`), so orders 7 to 10 require ``'trapezoid'`` or
``'simpson'`` and warn about their cost.  They are there for accuracy
studies rather than production runs.

.. note::
   How these numbers were obtained, since they are the basis for leaving the
   defaults alone.  Three measurements, all against a tight-tolerance
   reference computed at order 6 with the slab cap raised:

   #. **Cheapest configuration sweep.**  For each of seven cases -- Earth
      PREM 3ν at 0.5, 1 and 10 GeV; Earth PREM 5ν; an exponential density
      profile; the Sun at 100 MeV; and Earth 3ν with NSI -- and each of the
      targets :math:`10^{-4}`, :math:`10^{-6}`, :math:`10^{-8}`, the smallest
      slab count reaching that accuracy was found by explicit sweep at orders
      2, 4 and 6, with the adaptive loop switched off.  Counted in
      *Hamiltonian evaluations*, the optimal order rose monotonically with
      tolerance in every case.
   #. **Wall-time confirmation.**  Evaluation count turned out to be a poor
      proxy: the fixed per-slab overhead (array setup, the eigendecomposition
      for the matrix exponential, the slab product) outweighs the node count,
      so fewer slabs matters more than fewer evaluations.  Re-timing the same
      optima is what produced the table above, and it moved the crossover --
      order 2 wins on evaluations at :math:`10^{-4}` but loses on wall time.
   #. **Seed prototype, rejected.**  Because the starting slab count comes
      from a phase target that is order-independent (:math:`2\pi` radians per
      slab), an order-aware target was prototyped and A/B tested over 45
      configurations (five cases × three orders × three tolerances).  It gave
      no speed-up, and cost up to 20% on the energy scan: the final slab
      count is set by the refinement loop, not the seed, so starting coarser
      only adds an iteration.  The seed was left as it is.

Silent vectorization and the energy-batched scan engine
-------------------------------------------------------------

Two further layers of performance engineering do not change any physics
and require no change to user code *for correctness* -- though the first
of them rewards one:

* **Silent Hamiltonian vectorization.**  A user-supplied Hamiltonian or
  density-profile function is probed once: if it accepts an array of
  positions and returns a matching stack of matrices (verified against a
  scalar spot-check), that vectorized form is used for every subsequent
  evaluation; otherwise Magνs falls back transparently to evaluating it
  one point at a time.  Repeated evaluations of a density profile on
  identical position grids (common across an energy scan, where only the
  vacuum term of the Hamiltonian depends on energy) are additionally
  cached.

  **The fallback is correct but slow, and how slow is worth knowing.**
  The engine samples the Hamiltonian at every quadrature node of every
  slab -- a few hundred positions for a single probability, repeated at
  each level of the adaptive refinement -- so a scalar-only function
  turns that into a Python loop.  Measured on a three-flavor
  exponential-density profile, making the same ``H_func`` array-capable
  cut the time per :func:`~magnus.oscprob.osc_prob` call from 7.8 ms to
  1.7 ms, a factor of 4.6, with bit-identical output.  See
  :ref:`array-capable-hamiltonians` for how to write one.
* **Energy-batched scans.**  The standard, NSI, and LIV Hamiltonians all
  have the separable form :math:`H(E, l) = H_E(E) + V_\mathrm{CC}(l)\, M`,
  with :math:`H_E` collecting the energy-dependent (vacuum and LIV) terms
  and :math:`M` a fixed matrix.  When many energies share a single
  baseline, Magνs detects this and runs the *entire* scan as one batched
  pipeline: the potential is sampled once per refinement level and shared
  across all energies, and the quadrature, commutator algebra, matrix
  exponentials, and slab products all carry the energy axis as an
  additional batch dimension, with per-energy convergence masking so that
  energies that have already converged stop being recomputed.

.. _array-capable-hamiltonians:

Writing an array-capable Hamiltonian
--------------------------------------

If you pass your own ``H_func`` to :func:`~magnus.oscprob.osc_prob`, whether
it can be evaluated for many positions at once is the single largest factor
under your control.  The change is usually small: write the position
dependence with NumPy and let the matrix part broadcast.

.. code-block:: python

    # Slow: one position at a time
    def H_func(l):
        VCC = matter.VCC_func(l, num_density_e_func)
        return (1.0/energy)*h_vac + hamiltonians.hamiltonian_3nu_matter(VCC)

    # Fast: the same physics, all positions at once
    e00 = np.diag([1.0, 0.0, 0.0])
    def H_func(l):
        l = np.asarray(l, dtype=float)
        VCC = VCC_central*np.exp(-(l/gd.UNIT_KM)/l_scale)   # an array
        return (1.0/energy)*h_vac + VCC[..., None, None]*e00

The ``[..., None, None]`` is what does the work: it turns one potential per
position into a stack of matrices, so NumPy broadcasts where Python would
otherwise loop.  The function must still return a single ``(d, d)`` matrix
when handed a scalar -- the probe checks exactly that consistency before
trusting the vectorized form.

Note that this is a property of *your* function rather than of
:func:`~magnus.oscprob.osc_prob`, whose own inner loops are already
vectorized: the quadrature, the commutator algebra, the matrix exponentials
and the slab products all carry a batch dimension.

Two cases need no attention.  A Hamiltonian that **ignores** its argument --
constant density -- is detected separately and broadcast, so it is already on
a fast path.  And the ``osc_prob_{2,3,4,5}nu_*`` wrappers build their own
Hamiltonians, already array-capable, so this applies only when you supply one.

Since version 1.0.0 the fallback raises
:class:`~magnus.magnus.ScalarHamiltonianWarning` once per session, naming the
fix.  It was silent before, which is why the slow path is easy to sit on
without noticing -- the example notebooks shipped with it for years.

.. _validation:

Validation strategy
-----------------------

The `test suite <https://github.com/mbustama/Magnus/tree/main/tests>`_,
which runs in CI on every push (see the badge on :doc:`index`), validates
the methodology above directly:

* **The expansion terms** :math:`\Omega_1, \ldots, \Omega_6` are compared,
  term by term, to an independently coded implementation of the
  Bernoulli-number recursion, using a Hamiltonian with three independent,
  non-commuting generators — chosen specifically because a
  two-generator Hamiltonian causes one nested-commutator term of
  :math:`\Omega_4` to vanish identically, which would otherwise mask a
  coefficient error.
* **Convergence order** is checked against a high-accuracy
  ``scipy.integrate.solve_ivp`` (``DOP853``, ``rtol=1e-12``) solution of
  the same Schrödinger equation, confirming that each additional Magnus
  order improves the error, and that the Gauss-Legendre integrators
  achieve their nominal orders 2/4/6 (measured error reduction ratios of
  4.0/16.0/63.8 under slab halving, matching :math:`2^{\text{order}+1}`).
* **Physical probabilities** are cross-checked against closed-form
  expressions for 2ν and 3ν vacuum oscillations and 2ν constant-density
  matter oscillations (for both neutrinos and antineutrinos), and against
  ``solve_ivp`` for asymmetric, complex-valued profiles and for full
  PREM Earth crossings.
* **Time-ordering, unitarity, channel conventions, the silent
  vectorization path, and the energy-batched scan** each have dedicated
  regression tests, including a pure matrix-arithmetic check (no
  quadrature) that isolates the slab time-ordering from every other
  source of numerical error.

In practice, the default tolerance setting (``rtol = atol = 1e-3``, a
target for the difference between successive refinements rather than a
strict global error bound) delivers an actual accuracy of about
:math:`5\times10^{-4}` on Earth crossings, verified against
:math:`10^{-7}`-tolerance references.

See :doc:`references` for full citations of the works referred to above.


.. _conventions:

Conventions
------------

Everything below is a *choice*. None of it is forced by the physics, all of it
is forced by consistency, and a convention that is wrong **consistently** passes
every internal test — which is why they are written down here rather than left
in the code. Magνs has been bitten by exactly that: a reversed slab ordering, a
doubled antineutrino potential sign and a flipped two-flavour mass ordering were
all fixed on the same day, and each had been silently self-consistent.

Ordering of the probabilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every ``osc_prob_*`` function returns the probability matrix indexed
**initial flavour first**:

.. math::

   P[\nu_i][\nu_f] \;=\; P(\nu_i \to \nu_f) .

So ``P[1][0]`` is :math:`P(\nu_\mu \to \nu_e)`, not the reverse. Flavours are
in the standard order :math:`(e, \mu, \tau, s_1, s_2)`, so index 0 is always
:math:`\nu_e`.

Each **row** sums to one — a neutrino that started as :math:`\nu_i` ends as
something. Each column also sums to one, but that is a consequence of unitarity
rather than a separate statement. Passing ``nu_i`` and ``nu_f`` returns that one
entry instead of the matrix.

For a batched call the point index comes **first**: the shape is
``(n_points, d, d)``, so ``P[:, 1, 0]`` is :math:`P_{\mu e}` along a scan.

Sign of the matter potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The charged-current potential enters the electron-flavour diagonal entry,

.. math::

   H \;=\; H_\text{vac} \;+\; \mathrm{diag}(V_{CC},\, 0,\, \ldots) ,
   \qquad V_{CC} = +\sqrt{2}\, G_F n_e ,

and **for antineutrinos it changes sign**. That flip is applied once, inside
:func:`magnus.matter.vcc_func_from_rho_func`, so a caller passing
``nubar=True`` gets it automatically and code downstream must not apply it
again. It was applied twice once, which gave antineutrinos a positive potential
and answers that looked plausible.

Mass ordering
~~~~~~~~~~~~~

The ordering is carried by the **sign of** :math:`\Delta m^2_{31}`, not by a
flag: positive is normal, negative is inverted. ``OSC_PARAMS_DEFAULT`` is the
normal ordering, with :math:`\Delta m^2_{31} = +2.511 \times 10^{-3}`
eV\ :sup:`2`. It is NuFit 6.1, the same release
:func:`~magnus.globaldefs.load_nufit_params` returns by default, and is derived
from it rather than written out a second time.
``magnus.globaldefs.OSC_PARAMS_PREDEFINED`` also carries
``OSC_PARAMS_NU_FIT_6_1_SK_NO``, ``..._SK_IO`` and the 6.0 pair, if you want to
name the fit explicitly.

For two flavours the same rule applies to :math:`\Delta m^2`, which is what
makes the two-flavour case easy to get backwards: flipping its sign moves the
MSW resonance into the other channel, and the result is still a perfectly
ordinary-looking probability.

Mixing parameters
~~~~~~~~~~~~~~~~~

Angles are given as **sines** -- not as angles, and not as
:math:`\sin^2\theta`: ``s12`` is
:math:`\sin\theta_{12}`. Quoted fits usually give :math:`\sin^2\theta`, so
take the square root — ``gd.S12_NO_BF_NUFIT_6_0`` is ``np.sqrt(0.308)``.
Phases are in **radians**; the default :math:`\delta_{CP}` is 3.7001 rad, i.e.
212 degrees.

Two flavours take ``sth`` and ``Dm2`` rather than ``s12`` and ``D21``. Passing
the three-flavour names to a two-flavour call is not an error — the keys are
simply not recognised — so check the names if a two-flavour result looks
untouched by the parameters you set.

Units
~~~~~

Natural units throughout: energies in eV, baselines and positions in
eV\ :sup:`-1`, so that :math:`HL` is dimensionless.
:mod:`magnus.globaldefs` supplies the conversions — multiply by ``UNIT_KM``,
``UNIT_MEV``, ``UNIT_GEV``, ``UNIT_G_PER_CM3`` — and :ref:`units-table` lists
them.
