Quick Start Guide
==================

Magνs works two ways: as an **importable Python module** (the full API --
this page) and as a **command-line calculator** (one probability, no Python
required -- see :doc:`cli`). Use the module for anything programmatic
(scans, plots, fitting); use the CLI for a quick one-off number or a shell
script.

.. note::
   All positions, baselines, and energies in Magνs are in **natural
   units** (inverse eV and eV, respectively).  The :mod:`magnus.globaldefs`
   module provides conversion constants (``UNIT_KM``, ``UNIT_MEV``,
   ``UNIT_GEV``, ``UNIT_G_PER_CM3``, ...): multiply a physical quantity by
   the matching constant to convert it, e.g. ``100.0*gd.UNIT_KM`` for a
   100 km baseline.

Install Magνs with ``pip install --pre magnus`` (see :doc:`installation`),
then:

.. code-block:: python

   import numpy as np
   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

Oscillation parameters that are not passed explicitly default to the
`NuFit 6.0 <http://www.nu-fit.org>`_ best fit (normal ordering); pass
``s12``, ``D31``, ``dCP``, etc., or ``nubar=True``, to change them.

1. Vacuum oscillations
------------------------

.. code-block:: python

   energy = 1.0*gd.UNIT_GEV       # [eV]
   L = 1000.0*gd.UNIT_KM          # [eV^-1]

   # Full 3x3 probability matrix, P[i][j] = P(nu_i -> nu_j)
   P = oscprob.osc_prob_3nu_vacuum(energy, L)

   # A single channel, and an array of energies
   energies = np.logspace(-1, 1, 50)*gd.UNIT_GEV
   P_emu = oscprob.osc_prob_3nu_vacuum(energies, L, nu_i=gd.NUE, nu_f=gd.NUMU)

The same pattern applies to :func:`~magnus.oscprob.osc_prob_2nu_vacuum`,
:func:`~magnus.oscprob.osc_prob_4nu_vacuum` (3+1 sterile), and :func:`~magnus.oscprob.osc_prob_5nu_vacuum`
(3+2 sterile).

2. Matter with constant or exponential density
------------------------------------------------

.. code-block:: python

   rho = 5.0*gd.UNIT_G_PER_CM3   # constant matter density [eV^4]
   P = oscprob.osc_prob_3nu_matter_constant_density(energy, L, rho)

   # Exponentially falling density profile, e.g. inside a supernova
   P = oscprob.osc_prob_3nu_matter_exp_density(
       energy, L, L0=0.0, rho_central=1e3*gd.UNIT_G_PER_CM3,
       l_scale=100.0*gd.UNIT_KM)

3. The Earth (PREM) and the Sun
----------------------------------

.. code-block:: python

   # By direction (cosine of the zenith angle) and baseline
   P = oscprob.osc_prob_3nu_earth(energy, costhz=-0.8,
                                   L=2.0*6371.0*0.8*gd.UNIT_KM)

   # By source and detector location (the chord through the Earth is
   # computed automatically; see magnus.earth.loc_coords_dms for the
   # predefined named locations)
   P = oscprob.osc_prob_3nu_earth(energy, loc_ini='fermilab',
                                   loc_fin='homestake')

   # A full energy scan is batched internally when the baseline is shared
   energies = np.logspace(-0.3, 1.3, 200)*gd.UNIT_GEV
   P_scan = oscprob.osc_prob_3nu_earth(
       energies, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
       nu_i=gd.NUE, nu_f=gd.NUMU)

   # The Sun, with its built-in exponential electron-density profile
   P = oscprob.osc_prob_2nu_sun(10.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM,
                                 L0=0.0, sth=np.sqrt(0.308), Dm2=7.5e-5)

None of these pass ``integration_method``, so they use the default,
``'gl'`` -- the Gauss-Legendre commutator-free integrators (see
:doc:`methodology`), which are both the fastest and the most accurate choice
whenever the Hamiltonian is smooth within each slab, the common case.  Pass
``integration_method='trapezoid'`` (or ``'simpson'``) for a Hamiltonian with
a kink or a discontinuity *inside* a slab.

4. Beyond the Standard Model: NSI and LIV
---------------------------------------------

.. code-block:: python

   # Non-standard neutrino interactions, in the Earth
   P = oscprob.osc_prob_3nu_earth_nsi(
       energy, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
       eps_ee=0.1, eps_em=0.05j, eps_et=0.0, eps_mm=0.0, eps_mt=0.02,
       eps_tt=0.0)

   # Lorentz-invariance violation, in vacuum
   P = oscprob.osc_prob_3nu_vacuum_liv(
       energy, L, b1=gd.B1, b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1)

5. Your own Hamiltonian through the Earth or the Sun
--------------------------------------------------------

:func:`~magnus.oscprob.osc_prob_earth` and :func:`~magnus.oscprob.osc_prob_sun` handle the trajectory
geometry and the built-in density profile for you, while leaving the
physics completely open: supply ``H_func(energy, l, VCC)`` (``VCC`` is the
charged-current potential at position ``l``, with the antineutrino sign
already applied) or the two-argument ``H_func(energy, l)`` to ignore the
built-in potential entirely.

.. code-block:: python

   import magnus.hamiltonians as hamiltonians

   h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
       s12=np.sqrt(0.308), s23=np.sqrt(0.470), s13=np.sqrt(2.215e-2),
       dCP=212./180.*np.pi, D21=7.49e-5, D31=2.513e-3)

   def H(energy, l, VCC):
       vcc = np.asarray(VCC)
       return (1.0/energy)*h_vac + vcc[..., None, None]*np.diag([1.0, 0.0, 0.0])

   P = oscprob.osc_prob_earth(H, energy, loc_ini='fermilab', loc_fin='homestake')

``H`` may accept an array of positions ``l`` and return a stack of
Hamiltonians (position axis leading) for extra speed; this is detected
automatically, with a safe per-point fallback if it is not supported.

6. Fully generic: any Hamiltonian, any environment
--------------------------------------------------------

:func:`~magnus.oscprob.osc_prob` is the primordial function that every wrapper above calls
internally.  It accepts any square, Hermitian-valued function of position
(or a constant matrix), for any number of flavors:

.. code-block:: python

   def H_func(l):
       # Your own position-dependent Hamiltonian, returning a (d, d) array
       ...

   P = oscprob.osc_prob(H_func, t_ini=0.0, t_fin=L,
                         magnus_exp_order=4,
                         rtol=1e-4, atol=1e-4)

Find a full worked example of using :func:`~magnus.oscprob.osc_prob` directly for a
time-dependent matrix exponential (not necessarily a physical Hamiltonian)
in notebook 10; see :doc:`tutorials`.
