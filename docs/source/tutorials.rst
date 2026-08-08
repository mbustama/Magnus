Tutorial notebooks
==================

Fourteen worked notebooks live in `notebooks/
<https://github.com/mbustama/Magnus/tree/main/notebooks>`_, numbered in reading
order. Each carries its figures inline, so they can be read on GitHub without
being run, and each ends with a footer pointing at the previous notebook, the
next one, and the API reference.

They are the long form of :doc:`recipes`. A recipe is a few lines and its
output; a notebook is the same calculation with the reasoning around it — why
the convention is what it is, what happens at the edges, and what the numbers
were checked against. Both call the same functions, so there is no third
version to drift out of step.

To run them rather than read them::

   pip install "magnus[notebooks]"
   jupyter lab notebooks/

.. note::

   The notebooks are not built into this documentation — executing fourteen of
   them on every docs build would take the better part of an hour, and they are
   more useful where their outputs are already stored. The links below go to
   GitHub, which renders them with their figures.

   They are generated. ``notebooks/make_notebooks.py`` builds all fourteen,
   executes them and stores their outputs, and CI runs the same execution on
   every change to the notebooks or to the package: a notebook is documentation
   that claims to work, and running it is what makes the claim checkable. Edit
   the generator, not the ``.ipynb``.


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
   The same seven settings with three flavours and a CP-violating phase.
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
   Four- and five-flavour systems, where the extra states do not couple to the
   weak interaction.

`08. Non-standard interactions <https://github.com/mbustama/Magnus/blob/main/notebooks/08_magnus_bsm_nsi.ipynb>`_
   A new matter potential with off-diagonal couplings the Standard Model does
   not have.

`09. Lorentz-invariance violation <https://github.com/mbustama/Magnus/blob/main/notebooks/09_magnus_bsm_liv.ipynb>`_
   An energy dependence the vacuum term does not have.


What the method actually does
-----------------------------

The two notebooks for readers who want to know why the answers are what they
are, rather than how to ask for them.

`10. Phase-averaged probabilities <https://github.com/mbustama/Magnus/blob/main/notebooks/10_magnus_averaged_probability.ipynb>`_
   What survives when the oscillation phase is unresolvable — and why an error
   that is a phase disappears under averaging while one that is an envelope
   does not.

`11. The matrix exponential <https://github.com/mbustama/Magnus/blob/main/notebooks/11_magnus_matrix_exponential.ipynb>`_
   How :math:`\exp(\Omega)` is built, and why the route matters: the truncated
   series is anti-Hermitian, so its exponential is exactly unitary only if the
   exponential itself preserves that.

`12. The strategy parameter <https://github.com/mbustama/Magnus/blob/main/notebooks/12_magnus_adiabatic_hybrid_strategy.ipynb>`_
   ``'auto'`` against ``'magnus'``, timed and scored against ``solve_ivp``. The
   headline is not the speed: for three or more flavours the old default can hit
   its refinement caps and return a plausible, exactly unitary, **wrong** answer.


Where the limits are
--------------------

The two notebooks that show what Magνs gets wrong, and how to tell the two
kinds of wrong apart.

`13. A tabulated solar model <https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb>`_
   A real BS05 profile rather than an exponential. A case that looks wrong by
   1.4e-03 and is not: the error is a phase, and averaging removes 53x of it.

`14. A supernova shock front <https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb>`_
   The contrast. Here averaging does essentially nothing, because a shock
   changes the adiabaticity of the level crossing and so moves the conversion
   probability itself. Wrong, and loud about it.
