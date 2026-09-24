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

**Flexible.**  The Hamiltonian is an argument, not an assumption.  Standard
oscillations, non-standard interactions, Lorentz-invariance violation, sterile
states, pseudo-Dirac pairs and a model of your own all go through the same call.
Two to five flavors ship ready-made; the generic entry points take any dimension
and any profile, given as a function of position.

**Fast.**  A scan over energy or arrival direction is one batched call rather than
a loop, worth one to two orders of magnitude per probability.  The median call
over 164 Earth and solar configurations is **2 ms**; a 200-energy Earth-crossing
scan takes 76 ms, and a 100x100 oscillogram about 2 s.

**Accurate.**  Probabilities are unitary by construction at every setting, not by
refinement, and on a smooth profile Magνs reaches **2.9e-13** where a composition
of constant slabs floors at 2.5e-11.  Where it cannot certify its own answer, it
says so.

What it can compute
--------------------

* Oscillations through a **varying profile**: the Earth's PREM layers, any of
  twelve tabulated standard solar models, a supernova shock front, or any density
  you supply.
* The **phase-averaged** probability a solar or astrophysical experiment actually
  measures, over its energy resolution, without resolving the oscillation.
* The **evolution operator** itself, alongside the probabilities, for observables
  built from amplitudes.
* The same probabilities **from a shell**, with no Python, through the ``magnus``
  command.

What it has been used for
--------------------------

Each of these is one call with a different Hamiltonian, profile or observable.

* **Beam experiments** — appearance probabilities along the DUNE, T2K, Hyper-K and
  ESS chords, from two named sites (`notebook 04
  <https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb>`_).
* **Atmospheric oscillograms** — probability over zenith angle and energy in a
  single batched call (`notebook 06
  <https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb>`_).
* **Solar neutrinos** — twelve standard solar models, taken by name, and the
  averaged probability an experiment sees (`notebook 13
  <https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb>`_).
* **Supernova shock fronts** — where a travelling discontinuity changes the
  conversion probability itself (`notebook 14
  <https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb>`_).
* **Astrophysical flavor composition**, including pseudo-Dirac pairs that stay
  coherent after everything else has averaged (`notebook 29
  <https://github.com/mbustama/Magnus/blob/main/notebooks/29_magnus_pseudo_dirac.ipynb>`_).
* **A Hamiltonian of your own** — a long-range :math:`L_e - L_\mu` interaction
  sourced by the Sun's electrons, a cavity in the Earth's crust, geoneutrinos, or
  a jet inside a collapsing star (`notebook 19
  <https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb>`_).

:doc:`recipes` gives the code for each in a few lines; :doc:`tutorials` is the
guided tour.

.. _what-accuracy-means:

What "accurate" means here
---------------------------

Magνs is a numerical integrator, so unlike a closed-form method it has an error
that depends on how finely it discretizes.  Two properties are exact regardless,
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
   * - Earth crossing (PREM) at the default ``rtol = atol = 1e-3``, against the
       same call at 1e-7
     - median 9e-7, worst 2e-3
   * - Asymmetric profiles with complex Hamiltonians vs ``solve_ivp``/DOP853 at
       ``rtol=1e-12``
     - 1e-4 to 1e-7
   * - Energy-batched scan vs the per-point path, grid and tolerances pinned
     - 1e-12
   * - ``n_jobs > 1`` vs serial
     - exactly 0.0

The last row is the one worth reading twice: it is a *bit-identity* assertion
rather than a tolerance, so an optimization that changed an answer would fail it
rather than pass quietly.  The batched scan is held to 1e-12, and only with the
grid pinned, which is what isolates the batching: left to refine on its own it
builds the matter profile once for the whole scan, moving the answer at the 1e-6
level.  Notebook 24 measures that comparison.

**And the honest caveat.** ``rtol``/``atol`` are a stopping criterion --- the
ladder halts when two successive refinement levels agree --- not a bound on the
error of what is returned.  Usually that is conservative.  It is not always:
:ref:`what-rtol-atol-control` gives the measured detail, including a case where
two levels agreed coincidentally and the answer was wrong by 0.855.  Magνs warns
loudly in that regime, and :doc:`diagnostics` reports the measured
false-alarm rate of each warning.

.. _when-is-magnus-a-win:

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

**And separately from the physics, worth saying plainly so that nobody
evaluates Magνs for a job it was never meant to do:**

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

Magνs answers a constant-density call exactly: the series terminates at its first
term, so the evolution operator is a single exponential carrying no discretization
at all, and a whole scan is one batched exponential.  What NuOscProbExact offers
there is not a better answer but a leaner route to it, being built for that case
alone rather than carrying a general solver's dispatch, validation and refinement
ladder.  The margin is narrow and runs both ways: on a 3ν constant-density scan
Magνs costs 1.10 µs per energy against NuOscProbExact's 1.44 µs batched, while at
a *single* point it costs 33.8 µs against 19.9 µs.  :doc:`comparison` carries the
decision table --- which of the two to reach for, case by case --- together with
the measured speed and accuracy behind it.


.. _what-magnus-earns-its-place-on:

What Magνs earns its place on
------------------------------

The lists above are cases, not a reason.  Where a closed form exists, an exact
algebraic solution beats a truncated series --- that is arithmetic, not a
defect in anybody's code.  What is left is three axes, every one of them
measured against the other codes in :doc:`comparison`:

**Reach.**  A slab product has a *floor*: on a smooth profile its error bottoms
out near :math:`2.5\times10^{-11}` and then **rises**, because the round-off of
composing that many matrix products outgrows what another halving of the width
buys.  No setting reaches below it.  Magνs continues to
:math:`2.9\times10^{-13}`.

**Generality.**  An arbitrary :math:`H(t)` --- a custom Hamiltonian, a BSM term
nobody has diagonalized, a profile interpolated from a simulation --- needs no
per-model work, because nothing in the method assumes a form for :math:`H`.
The SU(N) closed forms stop at SU(4); Magνs has no ceiling.

**Pre-packaged observables.**  ``average=True`` returns the phase-averaged
probability a solar experiment actually measures, without resolving some 13 000
radians of phase and averaging the result yourself.  Neither of the other codes offers it.  And every entry point
can hand back the converged evolution operator alongside the probabilities
(``return_evolution_operator=True``), for the observables that are built from
amplitudes rather than from probabilities.

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
See :ref:`write-h-func-vectorized`.

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

:doc:`performance` reports where the time goes, and what was tried and
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
  named detector sites), the Sun on an exponential fit or any of twelve standard
  solar models (:doc:`solar_models`), or any density profile you supply.
* **Beyond the Standard Model**: non-standard neutrino interactions (NSI)
  and CPT-odd Lorentz-invariance violation (LIV), for every flavor count and
  environment above.
* **Magnus expansion to order 10**, with the term recursion verified
  term-by-term against the literature, and three integration methods.  The
  default, **Gauss-Legendre collocation integrators**, reaches orders
  2/4/6/8 from only 1/2/3/4 Hamiltonian evaluations per slab.  Cumulative
  trapezoid and Simpson quadrature reach order 10 and serve Hamiltonians that
  are not smooth within a slab.
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
   comparison
   functions
   solar_models
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
   engines
   performance
   diagnostics

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
