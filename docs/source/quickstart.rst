Quick start guide
=================

This page takes you from installation to the probabilities most users need: vacuum,
constant-density matter, the Earth and the Sun, for neutrinos and antineutrinos.
Every block runs as written, in order.  Magνs also works from the shell, one
probability at a time, with no Python (:doc:`cli`).

Your first probability
----------------------

Install with ``pip install magnuspy`` (the import name is ``magnus``; see
:doc:`installation`), then:

.. code-block:: python

   import numpy as np
   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

   # 3 flavors, 1 GeV, 1300 km of vacuum
   P = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM)

   print(P[gd.NUMU][gd.NUE])        # P(nu_mu -> nu_e) = 0.0313

``P`` is the 3×3 matrix of probabilities, **initial flavor first**:
``P[i][f]`` is :math:`P(\nu_i \to \nu_f)`, and ``gd.NUE``, ``gd.NUMU`` and
``gd.NUTAU`` are 0, 1 and 2.  Every row sums to one.  The oscillation parameters
are the NuFIT 6.1 best fit, normal ordering, unless you pass others.

.. important::

   **Energies and distances are in natural units, not GeV and km.**  Multiply by
   ``gd.UNIT_GEV`` and ``gd.UNIT_KM`` on the way in, and divide by them on the way
   out.  ``osc_prob_3nu_vacuum(1.0, 1300.0)`` is accepted and means 1 eV and
   1300 eV\ :sup:`-1`.

.. _units-table:

Units
-----

Magνs works in natural units throughout: energies in eV, distances in
eV\ :sup:`-1`, so that the product :math:`HL` is dimensionless.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Quantity
     - Unit
     - Conversion constant (multiply by it)
   * - Neutrino energy
     - eV
     - ``UNIT_MEV`` = 1e6, ``UNIT_GEV`` = 1e9
   * - Baseline, position
     - eV\ :sup:`-1`
     - ``UNIT_KM`` = 5.068e9, ``UNIT_CM`` = 5.068e4; ``gd.SUN_RADIUS`` and
       ``gd.EARTH_RADIUS`` are in km
   * - Hamiltonian, matter potential
     - eV
     - none
   * - Mass density
     - eV\ :sup:`4`
     - ``UNIT_G_PER_CM3`` = 4.310e18
   * - Number density
     - eV\ :sup:`3`
     - ``UNIT_PER_CM3`` = 7.684e-15
   * - Mass-squared differences
     - eV\ :sup:`2`
     - none
   * - Mixing angles
     - :math:`\sin\theta`
     - none (``angles=`` accepts other forms; see below)
   * - CP phases
     - radian
     - none

.. _quickstart-conventions:

Conventions
-----------

They are the standard ones, stated here so that you can check them against other codes.

* **Mixing matrix**: the PDG parametrization,
  :math:`\mathbb{R} = \mathbb{R}_{23}(\theta_{23})\,\mathbb{R}_{13}(\theta_{13},\delta_{\rm CP})\,\mathbb{R}_{12}(\theta_{12})`,
  with :math:`\mathbb{R}_{e3} = \sin\theta_{13}\,e^{-i\delta_{\rm CP}}`.
* **Mass splittings**: ``D21`` :math:`= m_2^2 - m_1^2` and ``D31``
  :math:`= m_3^2 - m_1^2`.  The ordering is the sign of ``D31``: positive is
  normal, negative is inverted.
* **Antineutrinos**: ``nubar=True`` conjugates the mixing matrix and flips the sign
  of the matter potential, and of a Lorentz-violating term of even ``n_liv``
  (:doc:`conventions`).
* **Matter**: ``electron_fraction`` is :math:`Y_e`, the number of electrons per
  atomic mass unit of the material, so a mass density :math:`\rho` gives the electron
  number density :math:`n_e = \rho N_A Y_e`.  It is 0.5 unless given.  The Earth
  functions use one value per layer (0.4656 in the core, 0.4957 in the mantle);
  :doc:`functions` lists them.

:ref:`conventions` gives the details, and :ref:`coming-from-other-codes` sets them
beside those of GLoBES, Prob3++ and nuSQuIDS.

.. _glossary:

Terms used throughout
---------------------

* **Slab**: one step of the position grid; on each, Magνs exponentiates the Magnus
  expansion of :math:`H`, and the evolution operator is the product over slabs.
* **Ladder**: the refinement that raises the slab count until the answers at two
  successive counts agree to ``rtol`` and ``atol``.
* **Engine**: the algorithm that answers a call (the ladder, an energy-batched scan,
  a closed form, ...).  Magνs picks one from the shape of the request; :doc:`engines`
  lists them.
* **Strategy**: the ``strategy`` argument, which chooses between the Magnus-expansion
  engines only (``'magnus'``), the hybrid (``'hybrid'``) and letting Magνs decide
  (``'auto'``, the default).
* **Hybrid**: transport along the instantaneous eigenstates where the profile is
  adiabatic, with the Magnus expansion only across the windows where it is not
  (:doc:`adiabatic_strategy`).
* **Phase-averaged**: the probability averaged over oscillation phases too fast for a
  detector to resolve (``average=True``; :doc:`averaged_probability`).
* **Chord**: the straight path through the Earth between a source and a detector.
* **Oscillogram**: a map of probability over energy and zenith angle.

.. _nufit-parameters:

Choosing the oscillation parameters
-----------------------------------

Pass any parameter by name to override the default: ``s12``, ``s23``, ``s13``
(sines), ``dCP`` (radians), ``D21``, ``D31`` (eV\ :sup:`2`).
:func:`~magnus.globaldefs.load_nufit_params` returns exactly these six, for every
NuFIT release from 1.0 to 6.1:

.. code-block:: python

   energy = 1.0*gd.UNIT_GEV
   L = 1300.0*gd.UNIT_KM

   osc = gd.load_nufit_params('NuFIT 6.1', 'NO')        # the default set
   P = oscprob.osc_prob_3nu_vacuum(energy, L, **osc)

   # The inverted-ordering best fit (its theta23 and dCP differ too)
   osc_io_fit = gd.load_nufit_params('NuFIT 6.1', 'IO')

   # The same parameters with only the ordering flipped
   osc_io = {**osc, 'D31': -abs(osc['D31'])}

NuFIT quotes :math:`\Delta m^2_{32}` for the inverted ordering;
``load_nufit_params`` converts it, so ``D31`` is always :math:`m_3^2 - m_1^2`.
Older releases and their categories (with or without Super-Kamiokande data) are
listed in ``gd.NUFIT_GLOBAL_FITS``:

.. code-block:: python

   older = gd.load_nufit_params('NuFIT 5.2', 'NO', category='without_SK')

Vacuum: scans, one channel, antineutrinos
-----------------------------------------

Pass an array where one energy or one baseline would go.  ``nu_i`` and ``nu_f``
select one channel, and ``nubar=True`` gives antineutrinos:

.. code-block:: python

   energies = np.linspace(0.5, 5.0, 200)*gd.UNIT_GEV

   P_scan = oscprob.osc_prob_3nu_vacuum(energies, L)          # shape (200, 3, 3)
   P_mue = oscprob.osc_prob_3nu_vacuum(energies, L, nu_i=gd.NUMU, nu_f=gd.NUE)
   P_mue_bar = oscprob.osc_prob_3nu_vacuum(energies, L, nu_i=gd.NUMU, nu_f=gd.NUE,
                                           nubar=True)          # shape (200,)

   import matplotlib.pyplot as plt
   plt.plot(energies/gd.UNIT_GEV, P_mue, label=r'$\nu_\mu \to \nu_e$')
   plt.plot(energies/gd.UNIT_GEV, P_mue_bar, label=r'$\bar\nu_\mu \to \bar\nu_e$')
   plt.xlabel('Energy [GeV]')
   plt.ylabel('Probability')
   plt.legend()
   plt.show()

:doc:`plotting` has ready-made versions of this and other figures.  The same calls
work at two flavors, where there is no global fit to default to, so the mixing
and the splitting are required:

.. code-block:: python

   P2 = oscprob.osc_prob_2nu_vacuum(energy, L, sth=np.sqrt(0.5), Dm2=2.5e-3)

At four and five flavors (``osc_prob_4nu_vacuum``, ``osc_prob_5nu_vacuum``), the
sterile mixing defaults to zero.

Matter of constant density
--------------------------

A density can be given in natural units, or in g/cm³ with a flag.  Use one or the
other, never both:

.. code-block:: python

   rho = 2.848*gd.UNIT_G_PER_CM3                  # the average density along DUNE

   P_mue = oscprob.osc_prob_3nu_matter_constant_density(
       energies, L, rho, nu_i=gd.NUMU, nu_f=gd.NUE)

   # the same, with the density in g/cm^3
   P_mue = oscprob.osc_prob_3nu_matter_constant_density(
       energies, L, 2.848, density_matter_is_in_g_per_cm3=True,
       nu_i=gd.NUMU, nu_f=gd.NUE)

   # T2K: 295 km, 0.6 GeV, 2.6 g/cm^3
   P_t2k = oscprob.osc_prob_3nu_matter_constant_density(
       0.6*gd.UNIT_GEV, 295.0*gd.UNIT_KM, 2.6*gd.UNIT_G_PER_CM3,
       nu_i=gd.NUMU, nu_f=gd.NUE)                              # 0.0525 (0.0318 with nubar=True)

A density that falls exponentially, as in a supernova envelope, takes a central
density and a scale length:

.. code-block:: python

   P = oscprob.osc_prob_3nu_matter_exp_density(
       energy, L, L0=0.0, rho_central=1e3*gd.UNIT_G_PER_CM3,
       l_scale=100.0*gd.UNIT_KM)

The Earth
---------

A path through the Earth is fixed by two named sites, or by the zenith angle at
the detector and the length of the chord.  A value of :math:`\cos\theta_z = -1` is straight up
through the Earth's center, 0 is horizontal, and a neutrino with
:math:`\cos\theta_z > 0` crosses no Earth to reach a detector at the surface.

.. code-block:: python

   import magnus.earth as earth

   # Fermilab to the Sanford lab (DUNE): the chord is computed from the sites
   P = oscprob.osc_prob_3nu_earth(2.5*gd.UNIT_GEV, loc_ini='fermilab',
                                  loc_fin='homestake')     # P[1][0] = 0.0724

   # An upgoing atmospheric neutrino, cos(theta_z) = -0.8: a 10194 km chord
   costhz = -0.8
   L_chord = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM

   E_atm = np.logspace(0.0, 1.3, 200)*gd.UNIT_GEV
   P_atm = oscprob.osc_prob_3nu_earth(E_atm, costhz=costhz, L=L_chord,
                                      nu_i=gd.NUMU, nu_f=gd.NUE)

``magnus.earth.loc_coords_dms`` lists the named sites.  An energy scan is one
batched call; an oscillogram is one such call per zenith angle:

.. code-block:: python

   cos_grid = np.linspace(-1.0, -0.1, 50)
   E_grid = np.logspace(0.0, 1.5, 50)*gd.UNIT_GEV

   P_mumu = np.array([
       oscprob.osc_prob_3nu_earth(
           E_grid, costhz=c, L=earth.distance_traveled_inside_earth(c)*gd.UNIT_KM,
           nu_i=gd.NUMU, nu_f=gd.NUMU)
       for c in cos_grid])                                   # shape (50, 50)

Earth and Sun calls can print a ``MagnusConvergenceWarning``.  It reports that a
slab of the integration grid is wide, not that the result is wrong;
:doc:`diagnostics` explains every warning and what to do about it.

The Sun
-------

For solar neutrinos you almost always want the **phase-averaged** probability,
which is what a detector with a finite energy resolution measures.  Pass
``average=True``; without it you get the value at one exact energy and distance,
which oscillates rapidly:

.. code-block:: python

   R_sun = gd.SUN_RADIUS*gd.UNIT_KM        # gd.SUN_RADIUS is in km

   # nu_e survival at 8 MeV, produced at the center (L0 = 0)
   P_ee = oscprob.osc_prob_3nu_sun(8.0*gd.UNIT_MEV, R_sun, 0.0,
                                   nu_i=gd.NUE, nu_f=gd.NUE, average=True)   # 0.300

   # with a standard solar model instead of the built-in exponential fit
   P_ee = oscprob.osc_prob_3nu_sun(8.0*gd.UNIT_MEV, R_sun, 0.0,
                                   nu_i=gd.NUE, nu_f=gd.NUE, average=True,
                                   density_profile='B16-GS98')             # 0.332

Without ``average=True``, the exponential fit gives 0.175 at the same
energy.  :doc:`solar_models` lists the twelve standard solar models, and
:doc:`averaged_probability` explains the average.

New physics: NSI and LIV
------------------------

Non-standard interactions are given as couplings relative to the standard matter
potential.  Unset couplings are zero; the diagonal ones (``eps_ee``, ``eps_mm``,
``eps_tt``) are real and the off-diagonal ones may be complex:

.. code-block:: python

   P_nsi = oscprob.osc_prob_3nu_matter_nsi_constant_density(
       2.5*gd.UNIT_GEV, L, rho, eps_ee=0.1, eps_em=0.05j,
       nu_i=gd.NUMU, nu_f=gd.NUE)                              # 0.0848

   # Lorentz-invariance violation, in vacuum: b1, b2, b3 and Lambda in eV
   P_liv = oscprob.osc_prob_3nu_vacuum_liv(
       2.5*gd.UNIT_GEV, L, b1=1e-9, b2=1e-9, b3=2e-9, Lambda=1e12, n_liv=1,
       nu_i=gd.NUMU, nu_f=gd.NUE)                              # 0.0600

Every matter environment above has NSI and LIV versions, and vacuum has LIV ones
(``osc_prob_3nu_earth_nsi``, ``osc_prob_3nu_sun_liv``, and so on), at two to five
flavors; :doc:`functions` lists them all.

Your own Hamiltonian
--------------------

:func:`~magnus.oscprob.osc_prob_earth` and :func:`~magnus.oscprob.osc_prob_sun`
take the trajectory and the density profile from the package and the physics from
you: a function ``H(E, l, VCC)`` of the energy ``E``, the position ``l`` and ``VCC``,
the standard matter potential at ``l`` (already with the antineutrino sign).  Write it so that
it accepts an array of positions and returns a stack of matrices, which is several
times faster:

.. code-block:: python

   import magnus.hamiltonians as hamiltonians

   h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**osc)
   e_ee = np.diag([1.0, 0.0, 0.0])

   def H(E, l, VCC):
       # VCC[..., None, None] broadcasts over an array of positions
       return h_vac/E + np.asarray(VCC)[..., None, None]*e_ee

   P = oscprob.osc_prob_earth(H, 2.5*gd.UNIT_GEV, loc_ini='fermilab',
                              loc_fin='homestake')

Any Hermitian matrix function of position, of any dimension, goes through
:func:`~magnus.oscprob.osc_prob`, which every function above calls:

.. code-block:: python

   def H_of_l(l):
       vcc = 1.0e-13*np.exp(-np.asarray(l)/(500.0*gd.UNIT_KM))   # [eV]
       return h_vac/energy + vcc[..., None, None]*e_ee

   P = oscprob.osc_prob(H_of_l, t_ini=0.0, t_fin=L, rtol=1e-4, atol=1e-4)

.. note::

   Mixing angles are sines by default.  To pass them as published, use ``angles``:
   ``'sin2'`` for :math:`\sin^2\theta`, ``'rad'`` or ``'deg'`` for the angle (under
   ``'deg'`` the CP phase is in degrees too).  Every function that takes a mixing
   angle accepts it:

   .. code-block:: python

      h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
          s12=33.76, s23=43.28, s13=8.62, dCP=212.0, D21=7.537e-5, D31=2.511e-3,
          angles='deg')

   These are the NuFIT 6.1 values as published, rounded to two decimals, so this
   Hamiltonian agrees with the loader's to about 1e-4.  ``load_nufit_params``
   returns sines unless it is given the same ``angles``, so pass the same value
   to both.

Where next
----------

* :doc:`recipes`: a short, runnable snippet for each common task.
* :doc:`functions`: every entry point, grouped by environment.
* :doc:`tutorials`: the notebooks, with figures.
* :doc:`diagnostics`: what each warning means.
