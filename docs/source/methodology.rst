Methodology
============

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
:math:`n = 6` (odd Bernoulli numbers :math:`B_3 = B_5 = 0` vanish
identically, so only even-index commutator groups appear beyond
:math:`\Omega_3`); the coefficients and every term were verified
independently, term by term, against this recursion (see
:ref:`validation`).

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
``integration_method``:

**Cumulative quadrature (** ``'trapezoid'`` **,** ``'simpson'`` **).**
Sample :math:`A` on a uniform grid of ``n_tpts_per_slab`` points and
integrate with cumulative trapezoid or Simpson's rule.  General-purpose,
but the quadrature error (:math:`O(h^2)` or :math:`O(h^4)` in the grid
spacing :math:`h`) can dominate the Magnus truncation error at high orders
unless ``n_tpts_per_slab`` grows accordingly.

**Gauss-Legendre commutator-free integrators (** ``'gl'`` **).**
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
evaluations for the same accuracy and is the recommended default whenever
the Hamiltonian is smooth within a slab (which layer-aligned slabs, below,
make the common case even across the Earth).

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

Silent vectorization and the energy-batched scan engine
-------------------------------------------------------------

Two further layers of performance engineering do not change any physics
and require no change to user code:

* **Silent Hamiltonian vectorization.**  A user-supplied Hamiltonian or
  density-profile function is probed once: if it accepts an array of
  positions and returns a matching stack of matrices (verified against a
  scalar spot-check), that vectorized form is used for every subsequent
  evaluation; otherwise Magνs falls back transparently to evaluating it
  one point at a time.  Repeated evaluations of a density profile on
  identical position grids (common across an energy scan, where only the
  vacuum term of the Hamiltonian depends on energy) are additionally
  cached.
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
