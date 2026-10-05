.. Magnus documentation master file

Magνs: neutrino oscillations via the Magnus expansion
=====================================================

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

**Magνs** computes neutrino oscillation probabilities for two to five flavors, or
for any Hermitian Hamiltonian, in vacuum, in matter, through the Earth and through
the Sun.  Its evolution operator is unitary to round-off at any accuracy setting,
so every probability is non-negative and every row sums to one.

.. tip::
   **New here?**  Install it with ``pip install magnuspy``, then compute your first
   probability:

   .. code-block:: python

      import magnus.oscprob as oscprob
      import magnus.globaldefs as gd

      P = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM)
      print(P[gd.NUMU][gd.NUE])        # P(nu_mu -> nu_e) = 0.0313

   Energies and distances are in natural units, hence ``gd.UNIT_GEV`` and
   ``gd.UNIT_KM``.  The :doc:`quickstart` continues from here: matter, the Earth,
   the Sun, antineutrinos, new physics and a Hamiltonian of your own.

.. important::
   **Where to start**

   * :doc:`Quick start <quickstart>` and :doc:`installation`
   * :doc:`What it can compute, with code <recipes>`
   * `GitHub repository <https://github.com/mbustama/Magnus>`_
   * `Example notebooks <https://github.com/mbustama/Magnus/tree/main/notebooks>`_
     (:doc:`tutorials` gives a guided tour)
   * :doc:`How to cite <citing>`

**Flexible.**  The Hamiltonian is an argument, not an assumption.  Standard
oscillations, non-standard interactions, Lorentz-invariance violation, sterile
states, pseudo-Dirac pairs and a model of your own all go through the same call.
Two to five flavors ship ready-made; the generic entry points take any dimension
and any profile, given as a function of position.

**Fast.**  An energy scan is one batched call rather than a loop; batching makes
each probability one to two orders of magnitude cheaper.  The cost of a probability
follows how fast the density varies, not how many times the neutrino oscillates.
:ref:`performance` has the timings.

**Accurate.**  Magνs propagates the evolution operator with the **Magnus expansion**: it
exponentiates truncated integrals of the Hamiltonian over a chain of slabs, short consecutive
stretches of the path.  At a tight tolerance, it agrees with an independent integration to a few
parts in :math:`10^{12}` at two to five flavors.  Where it cannot certify its own answer, it
says so.

.. hint::
   **How do I say that?** It is pronounced like the name **Magnus**: the Greek
   letter **ν** (nu), the neutrino's symbol, stands in for the "nu" syllable.

What it can compute
-------------------

* Oscillations through a **varying profile**: the layers of the Preliminary Reference Earth
  Model (PREM), any of twelve tabulated standard solar models, a supernova shock front, or any
  density you supply.
* The **phase-averaged** probability a solar or astrophysical experiment
  measures, over its energy resolution, without resolving the oscillation.
* The **evolution operator** itself, alongside the probabilities, for observables
  built from amplitudes.
* The same probabilities **from a shell**, with no Python, through the ``magnus``
  command.

Examples
--------

Each of these is one call with a different Hamiltonian, profile or observable.

* **Beam experiments** — appearance probabilities along the DUNE, T2K, Hyper-K and
  ESS chords through the Earth, from two named sites (`notebook 04
  <https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb>`_).
* **Atmospheric oscillograms** — probability over zenith angle and energy, one
  batched energy scan per zenith angle (`notebook 06
  <https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb>`_).
* **Solar neutrinos** — standard solar models taken by name, and the averaged
  probability an experiment sees (`notebook 13
  <https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb>`_).
* **Supernova shock fronts** — where a traveling discontinuity changes the
  conversion probability itself (`notebook 14
  <https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb>`_).
* **Astrophysical flavor composition** — pseudo-Dirac pairs that stay coherent
  after everything else has averaged (`notebook 29
  <https://github.com/mbustama/Magnus/blob/main/notebooks/29_magnus_pseudo_dirac.ipynb>`_).
* **A Hamiltonian of your own** — a long-range :math:`L_e - L_\mu` interaction
  sourced by the Sun's electrons, a cavity in the Earth's crust, geoneutrinos, or
  a jet inside a collapsing star (`notebook 19
  <https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb>`_).

:doc:`recipes` gives the code for each in a few lines; :doc:`tutorials` is the
guided tour.

.. _what-accuracy-means:

What "accurate" means here
--------------------------

Magνs is a numerical integrator, so its error depends on how finely it discretizes.  One
property holds regardless: every truncation of the Magnus series is anti-Hermitian, so the
evolution operator is unitary at any order, any tolerance and any slab count.  The rest is
measured against checks that are independent of one another:

.. list-table::
   :header-rows: 1
   :widths: 58 42

   * - Checked against
     - Result
   * - The derivation: every expansion term, :math:`\Omega_1` to :math:`\Omega_{10}`, against
       terms generated independently from the recursion
     - agree to a relative 1e-11 at every order
   * - Closed forms: 2ν and 3ν vacuum, 2ν constant-density matter, ν and ν̄
     - agree to 1e-12
   * - An independent solver: DOP853 at ``rtol=1e-12``, ``atol=1e-14``; halving the slab
       width at orders 2, 4, 6
     - error divided by about 4, 16, 64
   * - Itself: repeated calls, and a baseline scan given in any order
     - identical, bit for bit
   * - Itself: a parallel run against a serial one
     - agree to the requested tolerance
   * - Itself: the energy-batched scan against the per-point path, grid pinned
     - agree to 1e-12 (asserted), 1e-14 (measured)
   * - A population: 40 random smooth profiles at the default tolerance
     - median error about 1e-8, one silent miss
   * - A population: 120 random piecewise-constant profiles, edges left undeclared
     - 19 answers outside the tolerance, every one warned

A *silent miss* is an answer outside the requested tolerance with no warning; it is the
failure that matters.  Larger batteries, run by hand, found about 4% silent misses on 145
random smooth profiles (each within three times the tolerance), none on 150 piecewise-constant
profiles with declared edges, and none on 164 Earth, solar, vacuum and constant-density
configurations.

The tolerances ``rtol`` and ``atol`` are a stopping rule, not a guarantee: the refinement ladder
(:ref:`glossary`) stops once two successive answers agree.  Two slab counts can be wrong by the
same amount and still agree, which is how a silent miss happens.  At the default
``rtol = atol = 1e-3``, a probability through the Earth is usually far more accurate than that:
on eight chords at six energies, the median difference from the same call at 1e-7 is about
1e-6, and the largest about 1e-3.  :ref:`what-rtol-atol-control` gives the details, and
:doc:`diagnostics` what each warning means.

.. _when-is-magnus-a-win:

When is Magνs the right tool?
-----------------------------

Magνs takes any density profile, any number of flavors and any Hermitian Hamiltonian.  Among
the public codes compared in :doc:`comparison`, it is the cheapest to reach high accuracy on the
core-crossing Earth chord measured there.  Where the accumulated phase is extreme and the profile
varies slowly, as in the Sun, it transports the state along the instantaneous eigenstates and
keeps the expansion for the narrow windows where that fails (:doc:`adiabatic_strategy`).  Magνs
decides the handoff itself.
:ref:`when-to-use-magnus` sets out, with measurements, when another code is cheaper.

.. _when-is-magnus-not-the-right-tool:

When is Magνs not the right tool?
---------------------------------

Some limits belong to the method, and no implementation would remove them:

* **Open systems.**  Decoherence, coupling to a bath and decay to invisible states remove
  probability or damp the coherence between mass eigenstates.  They need a density matrix
  under a non-unitary evolution equation, or an anti-Hermitian term, neither of which a
  unitary method admits.  nuSQuIDS is the tool for those.
* **Collective oscillations.**  In a dense neutrino gas the Hamiltonian depends on the
  flavor content of the neutrinos, so the problem is nonlinear.  In Magνs the Hamiltonian is
  fixed before the propagation.  Magνs could be the propagator inside a self-consistent
  iteration, but it does not ship one.
* **A feature narrower than every grid.**  Every engine (the algorithm that answers a call;
  :ref:`glossary`) samples the Hamiltonian on a grid of positions, so a feature narrower than
  the finest grid is missed by every engine.  Magνs scans the profile for such features
  and warns, naming the breakpoints to declare, but the scan does not catch every one.

Others belong to the implementation:

* **Constant density.**  The exponential comes wrapped in a general solver's dispatch and
  validation, so a code built for that case alone is cheaper.
* **Fluctuations at every scale.**  The structural diagnostics cannot see density
  fluctuations spread over every scale; for those, raise the slab-count floor ``n_slabs``.
* **Engine changes.**  Two nearly identical requests can be answered by different engines,
  so the error can change in a step between them (:doc:`engines`).

.. _what-magnus-is-not:

Magνs is also not a flux, cross-section or detector code, not a fitting framework and not an
event generator: it computes oscillation probabilities only.

.. _performance:

Performance
-----------

A single three-flavor probability through the Earth takes about 2 ms at the default
tolerance.  Over 164 Earth, solar, vacuum and constant-density configurations, the median
call takes 2 ms and the slowest under a second.  These times are per call, on one laptop,
and exclude the first call of a session.  That call also loads the compiled kernels, which
takes 0.1–0.3 s, or about 2 s on a machine's first run, when they compile.  The
configurations are ``docs/dev/adversarial_batteries/battery10_coverage.py`` and the timing
harness is ``timing.py`` beside it (in a source checkout).  Four things affect the cost; the
third is automatic:

* **Pass arrays.**  Every wrapper accepts arrays of energies, of baselines or both, and then
  shares work across the points, which makes it one to two orders of magnitude faster than a
  loop, at every number of flavors from two to five.
* **Write your** ``H_func`` **to accept an array of positions.**  It is then called once per
  refinement stage rather than once per quadrature node, which is several times faster, with
  identical output (:ref:`write-h-func-vectorized`).
* **Earth chords are symmetric.**  The Hamiltonian is evaluated on the first half of the chord
  and mirrored, which makes a call with an expensive Hamiltonian about 1.5× faster.
* **Ask for worker processes where no batched engine applies.**  Ten workers (``n_jobs=10``)
  make a per-point scan 2–3× faster.  Where a batched engine applies, ``n_jobs > 1`` sends
  the scan to the per-point path instead, about 10× slower on 5000 energies, so one
  process is the right default.

The refinement ladder works against these savings: it computes every slab count below the one
that converges.  On an Earth chord, that makes a call about 3–4× slower than one
given the right slab count in advance.  :doc:`performance` has the rest.

Salient features
----------------

* **Any number of flavors, any Hamiltonian**: validated wrappers for 2ν, 3ν,
  4ν (3+1) and 5ν (3+2) systems (:doc:`functions`), plus a
  generic entry point, ``osc_prob``, that takes a Hermitian Hamiltonian of any
  dimension.
* **Vacuum, matter, Earth and Sun**: constant-density matter, exponentially
  falling density profiles, the Earth (`Preliminary Reference Earth Model
  <https://doi.org/10.1016/0031-9201(81)90046-7>`_, including chords between
  named detector sites), the Sun on an exponential fit or a tabulated standard
  solar model (:doc:`solar_models`), or any density profile you supply.
* **Magnus expansion to order 10**, with every term checked against terms
  generated independently from the recursion, and three integration methods.
  The default **Gauss–Legendre collocation integrators** reach orders 2, 4, 6
  and 8 from only 1, 2, 3 and 4 Hamiltonian evaluations per slab.  Cumulative
  trapezoid and Simpson quadrature reach order 10.  A density jump or kink is
  declared with ``t_breakpoints`` and becomes a slab edge, where every method
  keeps its order.
* **Adaptive refinement** to a requested tolerance, starting from a slab count
  estimated from the accumulated phase and, across a scan, from the previous
  point's count; slab edges on density discontinuities; and an energy-batched
  scan engine for standard, NSI and LIV Hamiltonians.
* **Vectorized Hamiltonians**: a Hamiltonian or density function that accepts an
  array of positions is detected and evaluated on whole arrays; one that takes a
  single position still works, more slowly, and issues
  ``ScalarHamiltonianWarning``.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting started:

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Using Magνs:

   recipes
   examples
   tutorials
   functions
   conventions
   solar_models
   cli
   plotting

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: How it works:

   methodology
   expansion_terms
   adiabatic_strategy
   averaged_probability
   architecture
   engines
   performance
   diagnostics
   comparison

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference:

   api_reference
   citing
   references
   changelog

Author
------

Magνs was written by Mauricio Bustamante (mbustamante@gmail.com).  Bug reports
and questions are best raised as `GitHub issues
<https://github.com/mbustama/Magnus/issues>`_, which leave a public record
others can find.

Citing
------

If Magνs contributed to work you are publishing, please cite it and say which
version you used, since results can depend on it.  :doc:`citing` has the BibTeX
entry and what to state in the text: the version, the tolerance and the strategy.

License
-------

Magνs is released under the `GNU General Public License v3.0 only
<https://www.gnu.org/licenses/gpl-3.0>`_ (``GPL-3.0-only``).  The full text
ships with the source, as ``LICENSE`` in the repository root, and inside the
installed distribution.

You are free to use, study, modify and redistribute it, including for
commercial purposes, provided that derivative works are distributed under the
same license and with source available.  If you are unsure whether your
intended use is compatible, read the license itself rather than this summary.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
