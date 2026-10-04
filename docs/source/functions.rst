Available probability functions
===============================

.. contents::
   :local:
   :depth: 2


This page lists the ``osc_prob_*`` functions Magνs ships, grouped by
environment and scenario, with the name for each flavor count.  Use it when you
know what you want ("3-flavor, matter, with NSI") but not the function's name;
:doc:`api_reference` documents each one in full.

Two groups of functions are not listed here.  The layers below the wrappers (the
scenario functions, ``osc_prob_energy_baseline`` and ``osc_prob``) are described in
:doc:`architecture`.  The closed-form functions ``osc_prob_2nu_vacuum_std``,
``osc_prob_2nu_matter_std`` and ``osc_prob_3nu_vacuum_std``, used to validate the
numerical ones, are in :mod:`magnus.oscprob`.  :doc:`cli` is the command-line
calculator that wraps the same functions.

Every function below returns a full :math:`d \times d` probability matrix
(:math:`P[i][j] = P(\nu_i \to \nu_j)`), or a single channel if ``nu_i``
and ``nu_f`` are both given; every one also accepts ``nubar=True`` to
compute the antineutrino probability. Standard oscillation parameters
left unset default to the NuFIT 6.1 best fit (normal ordering)
:cite:p:`Esteban:2024eli`;
sterile-sector parameters (4th/5th flavor) default to zero mixing.

Vacuum
--------

No matter potential -- the flavor Hamiltonian is just the vacuum mass/mixing
term, evaluated once and scaled by :math:`1/E`.

.. list-table::
   :header-rows: 1
   :widths: 15 45 40

   * - Flavors
     - Standard
     - LIV
   * - 2
     - :py:func:`~magnus.oscprob.osc_prob_2nu_vacuum`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_vacuum_liv`
   * - 3
     - :py:func:`~magnus.oscprob.osc_prob_3nu_vacuum`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_vacuum_liv`
   * - 4 (3+1 sterile)
     - :py:func:`~magnus.oscprob.osc_prob_4nu_vacuum`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_vacuum_liv`
   * - 5 (3+2 sterile)
     - :py:func:`~magnus.oscprob.osc_prob_5nu_vacuum`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_vacuum_liv`

There is no "vacuum + NSI" family: NSI couplings scale the matter
potential, and vacuum has none to scale (the CLI rejects this combination
explicitly; see :doc:`cli`).

Pseudo-Dirac neutrinos in vacuum, with a sterile partner on any of the three
mass states, have their own function,
:py:func:`~magnus.oscprob.osc_prob_pseudo_dirac_vacuum`.  It takes the pairing
and the splittings, and returns :math:`(3 + n_\text{pairs})`-dimensional
probabilities.  Its ``average=True`` forms each pair phase from its splitting,
and stays exact where the generic route loses a small splitting; see
:ref:`the note on pseudo-Dirac pairs <avg-pseudo-dirac>`.

Matter, constant density
---------------------------

A user-supplied matter density, uniform along the trajectory.  ``rho`` is in
natural units, eV\ :sup:`4`: pass ``2.8*gd.UNIT_G_PER_CM3``, or pass ``2.8``
together with ``density_matter_is_in_g_per_cm3=True``, but not both.
``electron_fraction``, the electrons per atomic mass unit :math:`Y_e`
(:ref:`quickstart-conventions`), is 0.5 unless given.
The command line's ``--rho`` is in g cm\ :sup:`-3`.

.. list-table::
   :header-rows: 1
   :widths: 15 30 30 25

   * - Flavors
     - Standard
     - NSI
     - LIV
   * - 2
     - :py:func:`~magnus.oscprob.osc_prob_2nu_matter_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_matter_nsi_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_matter_liv_constant_density`
   * - 3
     - :py:func:`~magnus.oscprob.osc_prob_3nu_matter_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_matter_nsi_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_matter_liv_constant_density`
   * - 4
     - :py:func:`~magnus.oscprob.osc_prob_4nu_matter_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_matter_nsi_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_matter_liv_constant_density`
   * - 5
     - :py:func:`~magnus.oscprob.osc_prob_5nu_matter_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_matter_nsi_constant_density`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_matter_liv_constant_density`

Matter, exponential density
-------------------------------

A user-supplied matter density profile
:math:`\rho(l) = \rho_{\rm central}\, e^{-l/l_{\rm scale}}`, with ``rho_central``
in the same units as ``rho`` above and ``l_scale`` in eV\ :sup:`-1`.

.. list-table::
   :header-rows: 1
   :widths: 15 30 30 25

   * - Flavors
     - Standard
     - NSI
     - LIV
   * - 2
     - :py:func:`~magnus.oscprob.osc_prob_2nu_matter_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_matter_nsi_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_matter_liv_exp_density`
   * - 3
     - :py:func:`~magnus.oscprob.osc_prob_3nu_matter_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_matter_nsi_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_matter_liv_exp_density`
   * - 4
     - :py:func:`~magnus.oscprob.osc_prob_4nu_matter_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_matter_nsi_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_matter_liv_exp_density`
   * - 5
     - :py:func:`~magnus.oscprob.osc_prob_5nu_matter_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_matter_nsi_exp_density`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_matter_liv_exp_density`

Earth
-------

The Preliminary Reference Earth Model (PREM) density profile, along a
chord specified either by the cosine of the zenith angle (plus a
baseline) or by two named locations (``loc_ini``/``loc_fin``; see
:data:`magnus.earth.loc_coords_dms` for the predefined sites).
:math:`\cos\theta_z = -1` is straight up through the Earth's center and 0 is
horizontal; :func:`magnus.earth.distance_traveled_inside_earth` gives the chord
length for a zenith angle.  :math:`Y_e`, the electrons per atomic mass unit
(:ref:`quickstart-conventions`), is set per layer (0.4656 in the core, 0.4957 in the
mantle), and ``electron_fraction_core`` and its siblings override it.

Both named locations lie on the surface.  Either end of the trajectory can
instead be put underground with ``source_depth`` and ``detector_depth``,
which every function below accepts.  The zenith angle is measured at the
detector, so a buried detector also sees downward-going neutrinos
(:math:`\cos\theta_z > 0`), which cross its overburden; a detector on the
surface does not.  Naming ``detector_depth`` fixes where the trajectory ends, so
the baseline is computed rather than given.  A third keyword,
``density_matter_ocean``, replaces the density of PREM's outermost shell, a
global-average ocean that does not describe a detector under rock or ice.

.. list-table::
   :header-rows: 1
   :widths: 15 30 30 25

   * - Flavors
     - Standard
     - NSI
     - LIV
   * - 2
     - :py:func:`~magnus.oscprob.osc_prob_2nu_earth`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_earth_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_earth_liv`
   * - 3
     - :py:func:`~magnus.oscprob.osc_prob_3nu_earth`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_earth_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_earth_liv`
   * - 4
     - :py:func:`~magnus.oscprob.osc_prob_4nu_earth`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_earth_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_earth_liv`
   * - 5
     - :py:func:`~magnus.oscprob.osc_prob_5nu_earth`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_earth_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_earth_liv`

Sun
-----

The built-in exponentially-falling solar electron-density profile (see
:func:`magnus.oscprob.osc_prob_sun`), from an initial radial
position ``L0`` (0 is the center) to a final radial position ``L``, both in
eV\ :sup:`-1`; ``gd.SUN_RADIUS*gd.UNIT_KM`` is the surface.  A solar-neutrino
measurement is the phase-averaged probability: pass ``average=True``.
Every one takes ``density_profile`` to use one of twelve tabulated standard
solar models instead (see :doc:`solar_models`).

.. list-table::
   :header-rows: 1
   :widths: 15 30 30 25

   * - Flavors
     - Standard
     - NSI
     - LIV
   * - 2
     - :py:func:`~magnus.oscprob.osc_prob_2nu_sun`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_sun_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_2nu_sun_liv`
   * - 3
     - :py:func:`~magnus.oscprob.osc_prob_3nu_sun`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_sun_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_3nu_sun_liv`
   * - 4
     - :py:func:`~magnus.oscprob.osc_prob_4nu_sun`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_sun_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_4nu_sun_liv`
   * - 5
     - :py:func:`~magnus.oscprob.osc_prob_5nu_sun`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_sun_nsi`
     - :py:func:`~magnus.oscprob.osc_prob_5nu_sun_liv`

Generic entry points
------------------------

For anything the tables above don't cover -- any other number of flavors,
or a Hamiltonian that doesn't fit the vacuum/matter/NSI/LIV mold -- three
functions accept an arbitrary user-supplied Hamiltonian directly:

* :py:func:`~magnus.oscprob.osc_prob` -- the general Magnus ladder, the base layer:
  any Hamiltonian, any dimension, any environment you build yourself.
* :py:func:`~magnus.oscprob.osc_prob_earth` -- like ``osc_prob``,
  but handles the Earth-crossing geometry and PREM potential for you.
* :py:func:`~magnus.oscprob.osc_prob_sun` -- like ``osc_prob``,
  but handles the solar density profile for you.

The ``osc_prob_{N}nu_*`` functions above call these three internally; see
:doc:`architecture`.


Returning the evolution operator
-----------------------------------

Every function above returns probabilities, which is what most observables
need. Some need the amplitudes instead: the content of each mass eigenstate in
the state that leaves a star, the flavor composition at a detector so far away
that the phases have averaged, or any quantity built from a product of
operators. For those, pass ``return_evolution_operator=True`` to any of the
functions, and the call returns the pair ``(P, U)`` instead of ``P`` alone:

.. code-block:: python

    import numpy as np
    import magnus.oscprob as oscprob
    import magnus.globaldefs as gd
    import magnus.hamiltonians as hamiltonians

    osc = gd.load_nufit_params('NuFIT 6.1')
    P, U = oscprob.osc_prob_3nu_matter_exp_density(
        1.0*gd.UNIT_GEV, 5000.0*gd.UNIT_KM, 0.0, 10.0, 1000.0*gd.UNIT_KM,
        density_matter_is_in_g_per_cm3=True, return_evolution_operator=True,
        **osc)

``P`` is exactly what the call returns without the keyword, so ``nu_i``, ``nu_f``
and the batching over arrays keep their meaning. ``U`` is the evolution operator
over the same interval, in the flavor basis, complex and unitary, with
``U[final, initial]`` the amplitude from the initial to the final flavor, so that
``P == (abs(U)**2).T``; for arrays of points it has shape ``(n, d, d)``, and
``P == np.swapaxes(abs(U)**2, -1, -2)``.

The operator comes from the general Magnus ladder, the only engine that forms it.
With the keyword set, the ladder compares the operator itself between refinement
levels, at the same ``rtol`` and ``atol``, so its phases converge as well as its
moduli.  The specialized engines of :doc:`engines` are skipped for the call, and a
baseline scan is computed point by point.  Two combinations raise an error,
because they form no operator: ``average=True`` and ``strategy='hybrid'``.

The phase-averaged content at a distant detector, from ``U`` and ``osc`` of the call
above, is then two lines:

.. code-block:: python

    R = hamiltonians.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP'])
    content = abs(R.conj().T @ U)**2      # mass-state content, per initial flavor
    P_far = (abs(R)**2 @ content).T       # phases averaged on the way

with ``R`` the mixing matrix in vacuum.  The transpose puts ``P_far`` in the same
order as ``P``: ``P_far[i][f]`` is the probability from flavor ``i`` to flavor ``f``.

The phase average is also available on the direct route.  ``average=True`` on
``osc_prob_energy_baseline``, ``osc_prob_earth`` and ``osc_prob_sun`` returns what
it returns on a wrapper, by the same three routes: closed form, adiabatic
transport, or an energy-window average across declared discontinuities.  The same
keywords apply: ``average_spread`` sets the spread, ``average_n_samples`` the
number of energies in a window average, and ``average_initial_state`` the starting
state, the flavor state by default.  A Hamiltonian given as a matrix, or as a
function of position alone, does not depend on energy, so no spread acts on it;
pairs of levels with a large phase are averaged and the others stay coherent (see
:doc:`averaged_probability`).  ``osc_prob`` computes a single point, so it refuses
``average``, as it refuses ``cumulative``.
