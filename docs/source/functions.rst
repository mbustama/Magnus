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
``osc_prob_energy_baseline``) is deliberately not listed
here; see :doc:`architecture` for what it does and when you would call it
directly.

See :doc:`architecture` for how these functions are organized internally
(the wrapper/middle/primordial layering), and :doc:`cli` for the
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
