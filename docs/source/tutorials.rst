Tutorial notebooks
==================

Twenty-nine worked notebooks live in `notebooks/
<https://github.com/mbustama/Magnus/tree/main/notebooks>`_, numbered in reading
order. Each carries its figures inline, so they can be read on GitHub without
being run, and each ends with a footer pointing at the previous notebook, the
next one, and the API reference.

They are the long form of :doc:`recipes`.  A recipe is a few lines and its
output; a notebook is the same calculation with the reasoning around it: why
the convention is what it is, what happens at the edges, and what the numbers
were checked against.

To run them rather than read them, clone the repository (the PyPI package does not
include the notebooks) and install from it, with a Jupyter front end::

   git clone https://github.com/mbustama/Magnus.git
   cd Magnus
   pip install -e ".[notebooks]" jupyterlab
   jupyter lab notebooks/

.. note::

   The notebooks are not built into this documentation, because executing all
   twenty-nine would take most of an hour.  The links below go to GitHub, which
   shows them with their stored outputs and figures.

   The notebooks are generated: ``notebooks/make_notebooks.py`` builds, executes
   and stores them, and CI executes them again on every change to the notebooks
   or the package.  To change a notebook, edit the generator, not the ``.ipynb``.


Start here
----------

The conventions everything else assumes, and the two systems every treatment of
oscillations opens with.

`01. Introduction <https://github.com/mbustama/Magnus/blob/main/notebooks/01_magnus_introduction.ipynb>`_
   The shortest path to a probability: single channels, arrays of energies and
   baselines, and what the returned matrix is indexed by.

`02. Two-neutrino probabilities <https://github.com/mbustama/Magnus/blob/main/notebooks/02_magnus_2nu_vacuum_matter.ipynb>`_
   Vacuum, constant density, exponential and Gaussian profiles, castle-wall and
   noisy potentials, the Earth and the Sun — each validated against the
   closed-form expression where one exists. The fullest tour of the supported
   matter profiles.

`03. Three-neutrino probabilities <https://github.com/mbustama/Magnus/blob/main/notebooks/03_magnus_3nu_vacuum_matter.ipynb>`_
   The same seven settings with three flavors and a CP-violating phase.
   Nothing about the method changes; the Hamiltonian is one dimension larger.


Geometry, and what experiments measure
--------------------------------------

Once the trajectory is a real one, the geometry starts to matter as much as the
Hamiltonian.

`04. Long baselines <https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb>`_
   Probabilities between two points on the Earth's surface — the geometry of
   DUNE, T2K, Hyper-K and ESS. Give the coordinates and the chord follows.

`05. Biprobability plots <https://github.com/mbustama/Magnus/blob/main/notebooks/05_magnus_biprobability.ipynb>`_
   Neutrino against antineutrino as the CP phase runs. The area enclosed is the
   CP violation an experiment is trying to measure.

`06. Oscillograms <https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb>`_
   Probability across zenith angle and energy at once. The workload that most
   rewards passing arrays rather than looping.


New physics
-----------

Each of these is a different Hermitian matrix in the same slot, so the
machinery is unchanged and only the Hamiltonian differs.

`07. Sterile neutrinos <https://github.com/mbustama/Magnus/blob/main/notebooks/07_magnus_bsm_sterile_nu.ipynb>`_
   Four- and five-flavor systems, where the extra states do not couple to the
   weak interaction.

`08. Non-standard interactions <https://github.com/mbustama/Magnus/blob/main/notebooks/08_magnus_bsm_nsi.ipynb>`_
   A new matter potential with off-diagonal couplings the Standard Model does
   not have.

`09. Lorentz-invariance violation <https://github.com/mbustama/Magnus/blob/main/notebooks/09_magnus_bsm_liv.ipynb>`_
   An energy dependence the vacuum term does not have.

`29. Pseudo-Dirac neutrinos <https://github.com/mbustama/Magnus/blob/main/notebooks/29_magnus_pseudo_dirac.ipynb>`_
   A sterile partner for each mass state, split by a :math:`\delta m^2` small
   enough that the pair stays coherent after everything else has averaged.  It
   is the physical case behind the coherent-block form of
   :doc:`averaged_probability`.


How the method works
--------------------

Three notebooks for readers who want to know why the answers are what they
are, rather than how to ask for them.

`10. Phase-averaged probabilities <https://github.com/mbustama/Magnus/blob/main/notebooks/10_magnus_averaged_probability.ipynb>`_
   What remains when the oscillation phase cannot be resolved, and why an
   error in the phase disappears under averaging while an error in the
   amplitude does not.

`11. The matrix exponential <https://github.com/mbustama/Magnus/blob/main/notebooks/11_magnus_matrix_exponential.ipynb>`_
   How :math:`\exp(\Omega)` is built, and why the method matters: the truncated
   series is anti-Hermitian, and its exponential stays unitary to round-off only
   if it is computed in a way that preserves that property.

`12. The strategy parameter <https://github.com/mbustama/Magnus/blob/main/notebooks/12_magnus_adiabatic_hybrid_strategy.ipynb>`_
   ``'auto'`` against ``'magnus'``, timed and scored against ``solve_ivp`` from
   two to five flavors.  ``'magnus'`` is more accurate in two of the seven
   cases.  At two flavors, though, it is slower than ``solve_ivp`` itself, and
   on the NSI cases it stops short of the tolerance.


Where the limits are
--------------------

Three notebooks that show where Magνs is inaccurate, and how to tell an error
in the phase from an error in the amplitude.

`13. Tabulated solar models <https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb>`_
   The twelve standard solar models of :doc:`solar_models`, taken by name. On
   BS2005-AGS,OP, it separates the instantaneous probability from the one an
   experiment measures.  Averaging a scan over a window does not converge to
   the second; ``average=True`` computes it in closed form, matching the
   adiabatic MSW expression to 3e-16.  On that observable, the twelve models
   agree to 2.4e-3, and the exponential fit is off by up to 0.1.

`14. A supernova shock front <https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb>`_
   Here, averaging barely reduces the error, because a shock changes the
   adiabaticity of the level crossing and so moves the conversion probability
   itself.  On this profile, Magνs warns that its answer is inaccurate.

`23. When averaging rescues you <https://github.com/mbustama/Magnus/blob/main/notebooks/23_magnus_when_averaging_helps.ipynb>`_
   The mechanism behind those two, isolated on a vacuum probability: averaging
   suppresses a phase error a hundredfold, but an amplitude error only by a
   fixed factor of about seven, however many cycles are averaged.


Conventions worth getting right
-------------------------------

Two notebooks about the places where a wrong answer looks exactly like a right
one.

`15. Antineutrinos, done properly <https://github.com/mbustama/Magnus/blob/main/notebooks/15_magnus_antineutrinos.ipynb>`_
   Conjugating the PMNS matrix and flipping the matter potential are two
   separate things, and doing one without the other returns a plausible wrong
   answer: the correct probability is 0.014, against 0.057 and 0.023 with only
   one of the two changes.

`18. Unusual density profiles <https://github.com/mbustama/Magnus/blob/main/notebooks/18_magnus_unusual_density_profiles.ipynb>`_
   Five profiles with the same mean density differing by up to 0.98 in
   probability — and the single rearrangement that changes nothing, exactly,
   whenever :math:`\delta_{\rm CP}` is 0 or :math:`\pi`.


Physics questions
-----------------

`16. Exact versus the approximations <https://github.com/mbustama/Magnus/blob/main/notebooks/16_magnus_exact_vs_approximations.ipynb>`_
   The textbook closed forms are exact, and Magνs reproduces them to
   :math:`10^{-14}`.  The error comes from replacing a varying density by its
   mean: up to 0.58 in probability on a core-crossing chord.

`17. Mass ordering and the octant <https://github.com/mbustama/Magnus/blob/main/notebooks/17_magnus_ordering_and_octant.ipynb>`_
   The ordering is carried entirely by the sign of ``D31``. Through the core it
   separates the two by 0.44; the octant, by about 0.015.


Using and diagnosing the machinery
----------------------------------

`19. Bring your own Hamiltonian <https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb>`_
   The interface is one callable returning a Hermitian matrix.  Covers writing
   it for arrays of positions, and what the Earth entry point declares for you.

`20. Numerical edge cases <https://github.com/mbustama/Magnus/blob/main/notebooks/20_magnus_numerical_edge_cases.ipynb>`_
   Exact degeneracies, zero baselines and empty requests all return numbers
   rather than ``NaN``.  The notebook also explains each warning class, and
   which ones to act on.

`21. What rtol and atol promise <https://github.com/mbustama/Magnus/blob/main/notebooks/21_magnus_what_tolerance_means.ipynb>`_
   A stopping criterion, not an error bound. Measured against an independent
   ``solve_ivp`` oracle: a request for :math:`10^{-2}` came back wrong by
   :math:`4.8\times10^{-2}` and reported success.

`22. Which engine answered, and why <https://github.com/mbustama/Magnus/blob/main/notebooks/22_magnus_which_engine_answered.ipynb>`_
   Seven engines (eight registered) in five families, and
   ``cross_check_strategies``, which flags a problem without a reference solution:
   two independent methods that disagree point to an inaccurate answer, although
   two that agree are not proven right.

`24. Performance <https://github.com/mbustama/Magnus/blob/main/notebooks/24_magnus_performance.ipynb>`_
   Which optimizations pay, measured as the notebook runs, and when each of them
   gains nothing.


Against other codes, and against the parameters
-----------------------------------------------

Where an independent method is the judge rather than Magνs itself.

`25. Against other codes <https://github.com/mbustama/Magnus/blob/main/notebooks/25_magnus_against_other_codes.ipynb>`_
   NuOscProbExact and nuSQuIDS on the same problems, each comparison checked
   against a third method: where a closed form is faster, where its accuracy
   stops improving, and a supernova shock where the *width of the front*, not
   the physics in it, decides which method suits the case.

`26. Fourteen years of NuFIT <https://github.com/mbustama/Magnus/blob/main/notebooks/26_magnus_nufit_evolution.ipynb>`_
   How the parameter likelihood, not just the best fit, moves the probability —
   and why the spread is comparable to effects other notebooks treat as
   signals.


Watching it happen
------------------

`27. Animated scenes <https://github.com/mbustama/Magnus/blob/main/notebooks/27_magnus_animations.ipynb>`_
   Nine parameter sweeps, drawn as filmstrips so that the notebook can be read
   without running it.  Four are the scenes of `NuOscProbExact's notebook 19
   <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/19_animations.ipynb>`_,
   computed here for comparison.  The other five need what a closed-form slab
   code does not have: a refinement ladder that checks convergence, a moving
   front, an averaged observable, and a Hamiltonian that varies along the path.

   Setting ``RENDER = True`` writes the scenes as GIFs, and
   ``tools/make_demo_video.py`` joins and compresses them.  The notebook describes
   the procedure, its cost and its pitfalls.


The paper
---------

`28. The paper's figures <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`_
   Every figure in the Magνs paper (``resources/paper/``), produced in one run.
   Magνs's numbers are computed as the notebook runs, so the figures always
   reflect the current code.  The other codes'
   numbers are read from the stored ``external_*.json`` datasets, so none of those
   codes has to be installed.  The notebook also shows that, on an Earth chord,
   the probability inherits the relative error of the matter potential almost one
   for one.  A comparison between codes there is therefore limited by their Earth
   models before it is limited by their solvers.  :doc:`examples` presents the
   paper's usage section, with its snippets and figures.
