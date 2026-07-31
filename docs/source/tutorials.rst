Tutorial Notebooks
===================

The ``notebooks/`` directory in the repository contains eleven runnable
Jupyter notebooks, numbered ``01``-``11`` in the order we'd suggest reading
them.

.. note::
   These pages are static documentation and do not execute the notebooks
   inline.  Follow a link below to view (or download and run) the notebook
   on GitHub.

.. list-table::
   :widths: 8 30 42 20
   :header-rows: 1

   * - #
     - Notebook
     - What it covers
     - Read this if...
   * - 01
     - `magnus_introduction.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/01_magnus_introduction.ipynb>`_
     - Basic usage of Magνs: computing a probability, single channels, arrays of energies and baselines.
     - **Start here.**  You will likely also want notebooks 02 and 03 before doing something useful.
   * - 02
     - `magnus_2nu_vacuum_matter.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/02_magnus_2nu_vacuum_matter.ipynb>`_
     - Two-neutrino probabilities in vacuum; constant-density matter; varying density (exponential, Gaussian); castle-wall and noisy density potentials; the Earth (PREM); the Sun (MSW resonance). Validated against the standard closed-form expressions.
     - You are working with a 2ν system, or want the fullest tour of the supported matter-density profiles.
   * - 03
     - `magnus_3nu_vacuum_matter.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/03_magnus_3nu_vacuum_matter.ipynb>`_
     - The same sequence of scenarios as notebook 02 (vacuum, constant/varying/castle-wall/noisy density, Earth, Sun), for the full three-flavor system.
     - You are working with the standard 3ν system -- the most common use case.
   * - 04
     - `magnus_long_baseline.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb>`_
     - Computing oscillation probabilities between two named points on the surface of the Earth (source and detector), as needed for long-baseline experiments (DUNE, Super-K, Hyper-K, T2K, ESS).
     - You need a probability between two real-world locations rather than a bare (direction, baseline) pair.
   * - 05
     - `magnus_biprobability.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/05_magnus_biprobability.ipynb>`_
     - Biprobability plots: the 3ν oscillation probability :math:`P(\nu_\alpha \to \nu_\beta)` vs. the antineutrino probability :math:`P(\bar\nu_\alpha \to \bar\nu_\beta)`, for different values of the CP-violating phase.
     - You are studying CP violation and want the standard biprobability visualization.
   * - 06
     - `magnus_oscillograms.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb>`_
     - Oscillograms: probability vs. neutrino direction (zenith angle) vs. energy, for neutrinos propagating inside the Earth. This is the workload the energy-batched scan engine (see :doc:`methodology`) targets directly.
     - You need a full (energy, direction) probability map, e.g. for atmospheric-neutrino analyses.
   * - 07
     - `magnus_bsm_sterile_nu.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/07_magnus_bsm_sterile_nu.ipynb>`_
     - Oscillation probabilities in systems of more than three neutrinos: 3+1 and 3+2 sterile-neutrino models.
     - Your model includes one or two additional (sterile) mass states.
   * - 08
     - `magnus_bsm_nsi.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/08_magnus_bsm_nsi.ipynb>`_
     - 2ν and 3ν probabilities under non-standard neutral-current interactions with matter, parametrized by the conventional :math:`\epsilon` couplings.
     - You need non-standard interactions (NSI) in matter.
   * - 09
     - `magnus_bsm_liv.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/09_magnus_bsm_liv.ipynb>`_
     - 2ν and 3ν probabilities with an additional effective, energy-dependent Hamiltonian representing Lorentz-invariance violation (LIV).
     - You need a CPT-odd LIV term in the Hamiltonian.
   * - 10
     - `magnus_matrix_exponential.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/10_magnus_matrix_exponential.ipynb>`_
     - Using Magνs as a general-purpose tool: computing the Magnus expansion and matrix exponential of an arbitrary time-dependent matrix :math:`A(t)`, i.e. :math:`\exp\!\left(\int_{t_i}^{t_f} A(t)\, dt\right)`, decoupled from any neutrino-physics interpretation.  Closes with the expansion's own coefficients: deriving them from the Bernoulli recursion at any order with :mod:`magnus.expansionterms` (see :doc:`expansion_terms`).
     - You want to call :func:`magnus.magnus.magnus_expansion` directly, outside of the oscillation-probability wrappers.
   * - 11
     - `magnus_adiabatic_hybrid_strategy.ipynb <https://github.com/mbustama/Magnus/blob/main/notebooks/11_magnus_adiabatic_hybrid_strategy.ipynb>`_
     - Live comparison of ``strategy='auto'``/``'hybrid'``/``'magnus'`` for 2-5 flavors, standard oscillations and an engineered BSM (NSI) resonance, each cross-checked against ``solve_ivp`` in both runtime and accuracy. Reproduces the validation in :doc:`adiabatic_strategy`.
     - You want to see the extreme-accumulated-phase problem (and its fix) reproduced live, or need a template for benchmarking your own Hamiltonian.
