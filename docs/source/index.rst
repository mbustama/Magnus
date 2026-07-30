.. Magnus documentation master file

Magνs: Neutrino Oscillations via the Magnus Expansion
========================================================

.. image:: https://github.com/mbustama/Magnus/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/mbustama/Magnus/actions
   :alt: CI Tests

.. important::
   **Important Links:**

   * `GitHub Repository <https://github.com/mbustama/Magnus>`_
   * `Example Notebooks <https://github.com/mbustama/Magnus/tree/main/notebooks>`_ (see also :doc:`tutorials` for a guided tour)

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
matter-dominated phase) — the required slab count can exceed the default
caps; Magνs then warns (``ToleranceNotAchievedWarning``) instead of
failing silently, and you should raise ``max_n_slabs`` (or use an adiabatic
approximation, the natural method in that regime).  And a tight-tolerance
ODE solver remains the best *reference* for validation — Magνs's own test
suite uses ``scipy.integrate.solve_ivp`` at ``rtol=1e-12`` as ground truth.
See :doc:`methodology` for the full numerical story, including how these
numbers were measured.

Salient Features
-----------------

* **Any number of flavors, any Hamiltonian**: dedicated, validated wrappers
  for 2ν, 3ν, 4ν (3+1 sterile), and 5ν (3+2 sterile) systems, plus a fully
  generic entry point (``osc_prob``) that accepts an arbitrary Hermitian
  Hamiltonian of any dimension.
* **Vacuum, matter, Earth, and Sun**: constant-density matter, exponentially
  falling density profiles, the Earth (`Preliminary Reference Earth Model
  <https://doi.org/10.1016/0031-9201(81)90046-7>`_, including chords between
  named detector sites), the Sun, or any density profile you supply.
* **Beyond the Standard Model**: non-standard neutrino interactions (NSI)
  and CPT-odd Lorentz-invariance violation (LIV), for every flavor count and
  environment above.
* **Magnus expansion to order 6**, with the term recursion verified
  term-by-term against the literature, and three integration methods —
  cumulative trapezoid/Simpson quadrature, and **Gauss-Legendre
  commutator-free integrators** that reach orders 2/4/6 from only 1/2/3
  Hamiltonian evaluations per slab.
* **Exact unitarity**, adaptive refinement to a requested tolerance with
  physics-informed starting slab counts and warm starts across scans, slab
  edges aligned with density discontinuities, and an energy-batched scan
  engine for standard/NSI/LIV Hamiltonians.
* **Silent vectorization**: Hamiltonian and density-profile functions that
  accept position arrays are detected and used automatically, with a safe
  scalar fallback.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   installation
   quickstart
   architecture
   methodology
   tutorials
   references

.. toctree::
   :maxdepth: 2
   :caption: API Reference:

   api/magnus/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
