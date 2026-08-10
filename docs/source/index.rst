.. Magnus documentation master file

Magνs: Neutrino Oscillations via the Magnus Expansion
========================================================

.. image:: https://github.com/mbustama/Magnus/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/mbustama/Magnus/actions/workflows/tests.yml
   :alt: CI Tests

.. image:: https://github.com/mbustama/Magnus/actions/workflows/lint.yml/badge.svg
   :target: https://github.com/mbustama/Magnus/actions/workflows/lint.yml
   :alt: Code Quality

.. image:: https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg
   :target: https://mbustama.github.io/Magnus/
   :alt: Documentation

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: License: GPL v3

.. image:: https://img.shields.io/badge/python-3.10+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.10+

.. image:: https://codecov.io/gh/mbustama/Magnus/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/mbustama/Magnus
   :alt: codecov

.. image:: https://img.shields.io/pypi/v/magnuspy.svg
   :target: https://pypi.org/project/magnuspy/
   :alt: PyPI

.. image:: https://pepy.tech/badge/magnuspy
   :target: https://pepy.tech/project/magnuspy
   :alt: Downloads

.. image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Code style: ruff

.. hint::
   **How do I say that?** Just like the name **Magnus** — the Greek letter
   **ν** (nu), the neutrino's symbol, simply stands in for the "nu"
   syllable.  (And since most of this package was written while the author
   was based in Denmark, you are equally welcome to say it `the Danish way
   <https://translate.google.com/?sl=da&tl=en&text=Magnus&op=translate>`_.)

.. important::
   **Important Links:**

   * :doc:`What it can compute, with code <recipes>`
   * `GitHub Repository <https://github.com/mbustama/Magnus>`_
   * `Example Notebooks <https://github.com/mbustama/Magnus/tree/main/notebooks>`_ (see also :doc:`tutorials` for a guided tour)
   * :doc:`How to cite <citing>`
   * :doc:`changelog`

**Magνs** computes neutrino oscillation probabilities between an arbitrary
number of flavors, for any given Hamiltonian, time-dependent or
-independent.  Internally, it propagates the neutrino evolution operator
using the **Magnus expansion**: rather than integrating the Schrödinger
equation step by step, it exponentiates truncated time-ordered integrals of
the Hamiltonian over a chain of position slabs.  Any truncation of the
Magnus series lives in the Lie algebra, so the resulting evolution operator
is **exactly unitary by construction** — probabilities are non-negative and
sum to one at machine precision, at any accuracy setting.

.. _when-is-magnus-a-win:

.. _what-accuracy-means:

What "accurate" means here
---------------------------

Magνs is a numerical integrator, so unlike a closed-form method it has an error
that depends on how finely it discretises.  Two properties are exact regardless,
and the rest is measured rather than asserted.

**Exact at any order, by construction.** Every truncation of the Magnus series is
anti-Hermitian, so its exponential is exactly unitary --- not unitary to within a
tolerance.  Truncating early costs accuracy, never norm.

.. list-table::
   :header-rows: 1
   :widths: 62 38

   * - Property
     - Measured agreement
   * - Unitarity, :math:`U^\dagger U - \mathbb{1}`
     - 1e-16 to 1e-13
   * - Magnus terms vs an independently coded Bernoulli recursion, orders 1--6
     - machine precision
   * - Gauss--Legendre convergence rate under slab halving, orders 2/4/6
     - error ratios 4 / 16 / 64
   * - 2ν and 3ν vacuum vs the closed-form expression
     - machine precision
   * - 2ν constant-density matter vs the closed form, ν and ν̄
     - machine precision
   * - Earth crossing (PREM) at the default ``rtol = atol = 1e-3``, against a
       1e-7-tolerance reference
     - ~5e-4
   * - Asymmetric profiles with complex Hamiltonians vs ``solve_ivp``/DOP853 at
       ``rtol=1e-12``
     - 1e-4 to 1e-7
   * - Energy-batched scan vs the per-point path
     - exactly 0.0
   * - ``n_jobs > 1`` vs serial
     - exactly 0.0

The last two rows are the ones worth reading twice: they are *bit-identity*
assertions, not tolerances, so an optimisation that changed an answer would fail
them rather than pass quietly.

**And the honest caveat.** ``rtol``/``atol`` are a stopping criterion --- the
ladder halts when two successive refinement levels agree --- not a bound on the
error of what is returned.  Usually that is conservative.  It is not always:
:ref:`what-rtol-atol-control` gives the measured detail, including a case where
two levels agreed coincidentally and the answer was wrong by 0.855.  Magνs warns
loudly in that regime, and :doc:`implementation_details` reports the measured
false-alarm rate of each warning.

When is Magνs a win?
------------------------

Compared to solving the propagation ODE directly (e.g., with an adaptive
Runge–Kutta solver), Magνs wins when one or more of these apply:

#. **The matter profile varies slowly compared to the oscillation length.**
   A Magnus slab is *exact* for a constant Hamiltonian no matter how many
   oscillation cycles it spans, so the slab size is set by how fast the
   *profile* changes, not by how fast the phase winds.  An ODE solver must
   resolve every oscillation.  For a 1 GeV neutrino crossing the Earth
   (PREM profile), Magνs needs ~10 slabs plus the ~16 layer crossings,
   versus thousands of right-hand-side evaluations for ``solve_ivp`` —
   measured: **~2 ms vs ~360-700 ms per probability** at comparable
   accuracy.

#. **You scan over energy and/or direction** (spectra, oscillograms,
   sensitivity studies).  The Magnus kernel is built from fixed,
   data-independent matrix operations, so slabs — and, for the
   standard/NSI/LIV Hamiltonians, the *entire energy axis* — evaluate as
   batched NumPy/BLAS calls.  Adaptive ODE integration is inherently
   sequential and cannot share steps across energies.  Measured: a
   200-energy Earth-crossing scan takes **76 ms** (0.4 ms per energy); a
   100×100 oscillogram takes **~2 s**.

#. **Unitarity matters more than raw local error** — long baselines, small
   probabilities, CP/T asymmetries.  Runge–Kutta iterates drift off the
   unitary manifold (probability leaks of ~1e-6 at typical tolerances,
   growing with baseline); the Magnus route has no leakage to leak, ever
   (probability rows sum to 1 to ~1e-14).

#. **You want arbitrary physics with no per-model work**: any number of
   flavors, any Hermitian Hamiltonian — sterile neutrinos, non-standard
   interactions, Lorentz-invariance violation, or your own matrix function
   of energy and position.

When is it *not* the best tool?  For a single probability at a single
energy, any method is fast enough.  For **extreme accumulated phases** —
e.g., ~10 MeV neutrinos crossing most of the Sun (~1e4 rad of
matter-dominated phase) — the plain Magnus slab-refinement method can need
a very large slab count, and warns (``ToleranceNotAchievedWarning``)
instead of failing silently if it hits its caps first.  This regime is now
handled automatically by ``strategy='auto'`` (the default for
``osc_prob_matter_std_potential``, ``osc_prob_matter_nsi``, ``osc_prob_liv``,
and every wrapper built on them, including every ``osc_prob_*_sun*``
function): an adiabatic-transport-plus-Magnus-patch strategy that stays
exactly unitary and is 50-25,000x faster than direct integration across the
validation grid — see :doc:`adiabatic_strategy` for the full derivation.
And a tight-tolerance ODE solver remains the best *reference* for
validation regardless — Magνs's own test suite uses
``scipy.integrate.solve_ivp`` at ``rtol=1e-12`` as ground truth.
See :doc:`methodology` for the full numerical story, including how these
numbers were measured.

.. _when-is-magnus-not-the-right-tool:

When is Magνs not the right tool?
-------------------------------------

Magνs solves the **unitary** Schrödinger equation for a Hermitian
Hamiltonian: any truncation of the Magnus series lives in the Lie algebra, so
the package is architecturally committed to norm-preserving, reversible
evolution.  That rules out several classes of problems that show up in
neutrino phenomenology:

#. **Quantum decoherence.**  Wave-packet separation, quantum-gravity-induced
   decoherence, or any model where coherence between mass eigenstates is
   damped over the baseline requires evolving a density matrix under a
   non-unitary master equation (e.g., Lindblad/GKSL), not a state vector
   under a Hamiltonian.  Magνs has no dissipative term and cannot represent
   one.

#. **Open-system coupling to a bath.**  Any scenario where the neutrino
   exchanges energy or phase information with an environment --
   collisional decoherence, thermal baths, stochastic scattering beyond the
   mean-field matter potential -- needs a reduced density matrix with
   dissipators, which is again outside a Hermitian-Hamiltonian, pure-state
   framework.

#. **Neutrino decay.**  Invisible or visible decay into lighter states
   removes probability from the system, so the evolution is no longer
   norm-preserving.  A Hermitian effective Hamiltonian cannot encode a decay
   width -- that requires an anti-Hermitian term, which breaks the
   unitarity the whole method relies on.

#. **Self-consistent collective oscillations.**  Dense-environment (e.g.,
   supernova) neutrino self-interactions, where the effective Hamiltonian
   depends on the (unknown, evolving) neutrino/antineutrino flavor content
   itself, are a nonlinear, self-consistent problem.  Magνs assumes the
   Hamiltonian is a *known* function of energy and position supplied by the
   caller, not a functional of the solution.

If your problem needs any of the above, look instead at packages built
around density-matrix/Lindblad evolution (for decoherence or decay) or
dedicated collective-oscillation codes (for self-interaction problems).

.. _what-magnus-is-not:

What it is not
---------------

Separately from the physics above, and worth saying plainly so that nobody
evaluates Magνs for a job it was never meant to do:

* **Not a solver for constant Hamiltonians in a hurry.**  It will do them, but
  a closed form beats an integrator every time; see
  :ref:`use-nuoscprobexact-instead`.
* **Not a flux, cross-section or detector code.**  It computes oscillation
  probabilities and stops there.
* **Not a fitting framework.**  There is no likelihood machinery; the
  probabilities are meant to be handed to whatever does that.
* **Not an event generator, and not an unfolding tool.**

.. _use-nuoscprobexact-instead:

When to use NuOscProbExact instead
-----------------------------------

Magνs integrates a Hamiltonian *across* each slab, which is what lets it follow
a density that changes as the neutrino moves.  That machinery is wasted --- and
slower than the alternative --- when the Hamiltonian does not change at all.

Reach for `NuOscProbExact <https://github.com/mbustama/NuOscProbExact>`_ instead
when **the Hamiltonian is constant, or piecewise constant**.  It expands the
Hamiltonian and the evolution operator in the SU(2), SU(3) and SU(4) bases,
which gives a closed form rather than a numerical integration: exact up to
floating-point round-off, and with no slab count to choose.

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Situation
     - Use this
     - Because
   * - Constant density
     - **NuOscProbExact**
     - One closed form, no integration
   * - Piecewise constant, tens of layers --- the Earth through PREM
     - **NuOscProbExact**
     - Each layer solved exactly, operators multiplied
   * - Smoothly varying, slow against the oscillation
     - Either
     - Slabbing converges quickly
   * - Smoothly varying, fast against the oscillation --- the Sun, adiabatic MSW
     - **Magνs**
     - Slabbing needs :math:`\sim 10^4` steps per resonance crossing
   * - A front resolved across many slabs --- a shock from a simulation snapshot
     - **Magνs**
     - Smooth on the slab scale, so fourth order beats second
   * - A front thin against the oscillation length --- a real hydrodynamic shock
     - **NuOscProbExact**
     - To any sampling method that is a jump, which is a closed form's home ground
   * - A kink, a tabulated profile
     - **Magνs**
     - ``t_breakpoints`` puts a slab edge on the discontinuity
   * - More than four flavours
     - **Magνs**
     - The SU(N) expansions stop at SU(4); Magνs has no ceiling
   * - Genuinely open systems: decay, decoherence
     - Neither
     - Needs a Lindblad solver, not a unitary one

The two packages share conventions, units and parameter defaults deliberately,
so a calculation can be moved between them as a cross-check.  That is worth
doing: agreement between two methods with different failure modes is stronger
evidence than either one's internal convergence check.

.. _what-magnus-earns-its-place-on:

What Magνs earns its place on
------------------------------

Worth stating in one place, because the table above is a list of cases and not
a reason.  Where a closed form exists, an exact algebraic solution beats a
truncated series; that is arithmetic rather than a defect, and it is why the
first two rows go the way they do.  What is left is three axes, all of them
measured in :doc:`notebook 25 <tutorials>` with every code timed in one process
on one machine and refereed by a method that is neither code's.

**Reach --- accuracy past where a slab product stalls.**  Composing slabs is
second order in the slab width, so halving it buys a factor of four; the
Gauss--Legendre Magnus expansion is fourth order and buys sixteen.  More
importantly the slab product has a *floor*.  On a smooth exponential profile at
three flavours its error bottoms out at :math:`2.5\times10^{-11}` near 16 000
slabs and then **rises** --- past that point the round-off of composing that
many matrix products costs more than another halving of the width buys, so
32 768 slabs is worse than 16 384.  No setting reaches below that.  Magνs
continues to :math:`2.9\times10^{-13}`.

**Generality --- an arbitrary** :math:`H(t)`, **and five flavours.**  A custom
Hamiltonian, a BSM term nobody has diagonalised, a profile interpolated from a
simulation: none of these need per-model work, because nothing in the method
assumes a form for :math:`H`.  The SU(N) expansions stop at SU(4); at five
flavours there is no comparison to draw at all, which is the same point stated
at its limit.

**Pre-packaged observables --- the quantity an experiment measures.**  A 5-MeV
neutrino leaving the Sun accumulates some 12 800 radians of phase, so the
*instantaneous* survival probability at the surface is neither measurable nor
stable: neighbouring energies land anywhere between 0.15 and 0.9.  What a solar
experiment measures is the phase-averaged probability, and ``average=True``
returns it directly, transporting along the levels of the instantaneous
Hamiltonian instead of propagating.  On one BS2005-AGS,OP model file Magνs
returns 40 averaged energies in 0.66 s, against 131 s for 12 *instantaneous*
ones from nuSQuIDS --- and recovering the observable from those means averaging
many such evaluations on top.  Neither of the other codes offers an averaging
flag; this is a different algorithm for the question being asked, not the same
algorithm run faster.

.. _performance:

Performance
------------

A single 3ν Earth probability takes about 2 ms at the default tolerance, and the
median call across 164 Earth and solar configurations is 2 ms with the slowest at
0.90 s.  Scans are what the code mostly does, and three things make them much
faster without changing any answer.

**Pass arrays instead of looping.**  Every wrapper takes an array of energies, of
baselines, or both.  For a position-dependent Hamiltonian the matter profile is
then built once for the whole scan rather than once per point, which is what the
energy-batched engine exists to do.

**Write your ``H_func`` so it accepts an array of positions.**  The single largest
factor under a caller's control: measured at **4.6x** on a 3ν exponential-density
profile, with bit-identical output.  A scalar-only Hamiltonian raises
:class:`~magnus.magnus.ScalarHamiltonianWarning` once per session, naming the fix.
See :ref:`write-h-func-vectorised`.

**An Earth chord is a palindrome.**  A neutrino crossing a spherically symmetric
Earth meets every radius twice, so the Hamiltonian is evaluated on the first half
of the slab chain and the rest follows by reversal.  That halves the calls to
your ``H_func``, so it is worth what your Hamiltonian costs:

.. list-table::
   :header-rows: 1
   :widths: 55 22 23

   * - Workload
     - Speed-up
     - Note
   * - Single point, plain PREM
     - 0.91x
     - a density lookup is too cheap to halve
   * - Single point, expensive ``H_func``
     - **1.41x--1.67x**
     -
   * - 12- and 40-energy scan, expensive ``H_func``
     - **1.56x--1.64x**
     -
   * - Energy scan, standard PREM
     - 1.00x
     - the separable engine already shares the profile

:data:`magnus.magnus.USE_PALINDROME` switches it off.  Standard PREM scans are
unaffected because the batched engine already evaluates the profile once and
shares it across energies --- the same saving, taken earlier.

**And one cost that runs the other way.**  The adaptive ladder computes the
probability at several slab counts and stops when two agree, so a call at a tight
tolerance is doing real extra work rather than being slow.  ``rtol=atol=None``
runs once at the grid you specify.

:doc:`implementation_details` reports where the time goes, and what was tried and
rejected.

Salient Features
-----------------

* **Two ways to use it**: as an importable Python module (the full API --
  see :doc:`quickstart`) or as a ``magnus`` command-line calculator for a
  single probability with no Python required (see :doc:`cli`).
* **Any number of flavors, any Hamiltonian**: dedicated, validated wrappers
  for 2ν, 3ν, 4ν (3+1 sterile), and 5ν (3+2 sterile) systems (see
  :doc:`functions` for the full listing), plus a fully generic entry point
  (``osc_prob``) that accepts an arbitrary Hermitian Hamiltonian of any
  dimension.
* **Vacuum, matter, Earth, and Sun**: constant-density matter, exponentially
  falling density profiles, the Earth (`Preliminary Reference Earth Model
  <https://doi.org/10.1016/0031-9201(81)90046-7>`_, including chords between
  named detector sites), the Sun, or any density profile you supply.
* **Beyond the Standard Model**: non-standard neutrino interactions (NSI)
  and CPT-odd Lorentz-invariance violation (LIV), for every flavor count and
  environment above.
* **Magnus expansion to order 6**, with the term recursion verified
  term-by-term against the literature, and three integration methods.  The
  default, **Gauss-Legendre commutator-free integrators**, reaches orders
  2/4/6 from only 1/2/3 Hamiltonian evaluations per slab; cumulative
  trapezoid/Simpson quadrature remains available for Hamiltonians that are
  not smooth within a slab.
* **Exact unitarity**, adaptive refinement to a requested tolerance with
  physics-informed starting slab counts and warm starts across scans, slab
  edges aligned with density discontinuities, and an energy-batched scan
  engine for standard/NSI/LIV Hamiltonians.
* **Silent vectorization**: Hamiltonian and density-profile functions that
  accept position arrays are detected and used automatically, with a safe
  scalar fallback.

.. toctree::
   :maxdepth: 2
   :caption: Getting started:

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Using Magnus:

   recipes
   tutorials
   functions
   cli
   plotting

.. toctree::
   :maxdepth: 2
   :caption: How it works:

   methodology
   expansion_terms
   adiabatic_strategy
   averaged_probability
   architecture
   implementation_details

.. toctree::
   :maxdepth: 2
   :caption: Reference:

   api_reference
   citing
   references
   changelog

Author
-------

Magnus was written by Mauricio Bustamante (mbustamante@gmail.com).  Bug reports
and questions are best raised as `GitHub issues
<https://github.com/mbustama/Magnus/issues>`_, which leave a public record
others can find.

Citing
-------

If Magnus contributed to work you are publishing, please cite it, and say which
version you used -- results can depend on it.  :doc:`citing` has the BibTeX
entry and the two or three things worth stating in the text.

License
=========

Magνs is released under the `GNU General Public License v3.0 only
<https://www.gnu.org/licenses/gpl-3.0>`_ (``GPL-3.0-only``).  The full text
ships with the source, as ``LICENSE`` in the repository root, and inside the
installed distribution.

You are free to use, study, modify, and redistribute it, including for
commercial purposes, provided that derivative works are distributed under the
same license and with source available.  If you are unsure whether your
intended use is compatible, read the license itself rather than this summary.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
