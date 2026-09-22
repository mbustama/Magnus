Available Oscillation-Probability Functions
==============================================

.. contents::
   :local:
   :depth: 2


This page lists every user-facing ``osc_prob_*`` function Magνs ships,
grouped by environment and scenario, with the exact function name for each
flavor count. It complements the API reference (generated from the
docstrings directly, see the *API Reference* section in the sidebar) by
showing the *shape* of the whole family at a glance -- useful when you know
roughly what you want ("3-flavor, matter, with NSI") but not the exact name.

The internal middle layer these wrappers dispatch through
(``osc_prob_vacuum``, ``osc_prob_matter_std_potential``,
``osc_prob_matter_nsi``, ``osc_prob_liv``, and
``osc_prob_energy_baseline``) is deliberately not listed here; see
:doc:`architecture` for what it does, when you would call it directly, and how
the wrapper/middle/primordial layering fits together.  :doc:`cli` is the
command-line calculator that wraps the same functions.

Every function below returns a full :math:`d \times d` probability matrix
(:math:`P[i][j] = P(\nu_i \to \nu_j)`), or a single channel if ``nu_i``
and ``nu_f`` are both given; every one also accepts ``nubar=True`` to
compute the antineutrino probability. Standard oscillation parameters
left as ``None`` default to the NuFIT 6.1 best fit (normal ordering)
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

Matter, constant density
---------------------------

A user-supplied matter density, uniform along the trajectory (``rho``,
in :math:`\text{g cm}^{-3}` by default).

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
:math:`\rho(l) = \rho_{\rm central}\, e^{-l/l_{\rm scale}}`.

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

Both named locations lie on the surface.  Either end of the trajectory can
instead be put underground with ``source_depth`` and ``detector_depth``,
which every function below accepts.  The zenith angle is measured at the
detector, so a buried detector also sees downward-going neutrinos
(:math:`\cos\theta_z > 0`) through its overburden, which a detector on the
surface has no path for at all.  Naming ``detector_depth`` fixes where the
trajectory ends, so the baseline is computed rather than given.  A third
keyword, ``density_matter_ocean``, replaces the density of PREM's outermost
shell: that shell is a global-average ocean, and a detector under rock or
ice is not under one.

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
position ``L0`` (default: the center) to a final radial position ``L``.

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

* :py:func:`~magnus.oscprob.osc_prob` -- the primordial function:
  any Hamiltonian, any dimension, any environment you build yourself.
* :py:func:`~magnus.oscprob.osc_prob_earth` -- like ``osc_prob``,
  but handles the Earth-crossing geometry and PREM potential for you.
* :py:func:`~magnus.oscprob.osc_prob_sun` -- like ``osc_prob``,
  but handles the solar density profile for you.

See :doc:`architecture` for how these three relate to the ``osc_prob_{N}nu_*``
functions above (they are, in fact, what those functions call internally).


Returning the evolution operator
-----------------------------------

Every function above returns probabilities, which is what most observables
need. Some need the amplitudes instead: the content of each mass eigenstate in
the state that leaves a star, the flavor composition at a detector so far away
that the phases have averaged, or any quantity built from a product of
operators. For those, pass ``return_evolution_operator=True`` to any of the
functions, and the call returns the pair ``(P, U)`` instead of ``P`` alone:

.. code-block:: python

    P, U = oscprob.osc_prob_3nu_matter_exp_density(
        energy, L, 0.0, rho_central, l_scale,
        density_matter_is_in_g_per_cm3=True,
        return_evolution_operator=True)

``P`` is exactly what the call returns without the keyword, so ``nu_i``, ``nu_f``
and the batching over arrays keep their meaning. ``U`` is the evolution operator
over the same interval, in the flavor basis, complex and unitary, with
``U[final, initial]`` the amplitude from the initial to the final flavor, so that
``P == abs(U)**2.T``; for arrays of points it has shape ``(n, d, d)``.

The operator comes from the general Magnus ladder, the one engine that forms it.
With the keyword set, the ladder compares the operator itself between refinement
levels, at the same ``rtol`` and ``atol``, so what comes back is converged in its
phases and not only in its moduli; the specialized engines of :doc:`engines`
stand aside for the call, and a baseline scan that would otherwise take the
cumulative traversal takes the per-point path instead. Two combinations are
refused with an error, because no operator exists to return: ``average=True``,
and ``strategy='hybrid'``.

The phase-averaged content at a distant detector, from the operator, is then two
lines:

.. code-block:: python

    content = abs(R.conj().T @ U)**2      # mass-state content, per initial flavor
    P_far = abs(R)**2 @ content           # phases averaged on the way

with ``R`` the mixing matrix in vacuum (``magnus.hamiltonians.pmns_mixing_matrix``).

The phase-averaged limit is likewise available on the direct route:
``average=True`` on ``osc_prob_energy_baseline``, ``osc_prob_earth`` and
``osc_prob_sun`` returns what the same keyword returns on a wrapper, by the same
three routes (closed form, adiabatic transport, or an energy-window average
across declared discontinuities); see :doc:`averaged_probability`.
``osc_prob`` computes one point and refuses the keyword by name, as it refuses
``cumulative``.
