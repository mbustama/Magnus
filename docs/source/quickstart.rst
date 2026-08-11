Quick Start Guide
==================

Magνs works two ways: as an **importable Python module** (the full API --
this page) and as a **command-line calculator** (one probability, no Python
required -- see :doc:`cli`). Use the module for anything programmatic
(scans, plots, fitting); use the CLI for a quick one-off number or a shell
script.

.. _units-table:

Units
------

Magνs works in **natural units** throughout: energies in eV, baselines and
positions in eV\ :sup:`-1`, so that the product :math:`HL` is dimensionless.
:mod:`magnus.globaldefs` supplies the conversions — multiply a physical
quantity by the matching constant, e.g. ``100.0*gd.UNIT_KM`` for a 100 km
baseline.

.. list-table::
   :header-rows: 1
   :widths: 38 26 36

   * - Quantity
     - Units
     - Constant
   * - Neutrino energy
     - eV
     - ``UNIT_MEV``, ``UNIT_GEV``
   * - Baseline, position
     - eV\ :sup:`-1`
     - ``UNIT_KM``, ``UNIT_CM``
   * - Hamiltonian
     - eV
     - ---
   * - Mass-squared differences
     - eV\ :sup:`2`
     - ---
   * - Matter potential
     - eV
     - ---
   * - Mass density
     - eV\ :sup:`4`
     - ``UNIT_G_PER_CM3``
   * - Number density
     - eV\ :sup:`3`
     - ``UNIT_PER_CM3``
   * - Mixing angles
     - given as :math:`\sin\theta`
     - ---
   * - CP phases
     - radian
     - ---

The last two are the ones to check first when a result looks untouched by the
parameters you set; see :ref:`conventions`.

Install Magνs with ``pip install magnuspy`` -- the distribution is
``magnuspy`` on PyPI, the import package is ``magnus`` (see
:doc:`installation`) -- then:

.. code-block:: python

   import numpy as np
   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

Oscillation parameters that are not passed explicitly default to the
`NuFit 6.0 <http://www.nu-fit.org>`_ best fit (normal ordering); pass
``s12``, ``D31``, ``dCP``, etc., or ``nubar=True``, to change them.

.. _nufit-parameters:

Choosing a global fit
---------------------

To use a different release, or the inverted ordering, ask
:func:`magnus.globaldefs.load_nufit_params` for it.  It returns **exactly the
six parameters** every ``osc_prob_3nu_*`` function takes -- ``s12``, ``s23``,
``s13``, ``dCP``, ``D21``, ``D31`` -- so the result can be passed straight
through:

.. code-block:: python

   energy = 1.0*gd.UNIT_GEV       # [eV]
   L = 1000.0*gd.UNIT_KM          # [eV^-1]

   osc = gd.load_nufit_params('NuFIT 6.1', 'NO')

   P = oscprob.osc_prob_3nu_vacuum(energy, L, **osc)

Every NuFit release from v1.0 to v6.1 is available, along with the
release-specific secondary category where one exists (``'with_SK'`` /
``'without_SK'`` from v4.0, ``'LEM'`` / ``'LID'`` for v2.1); omitting
``category`` takes the release's preferred one.  ``gd.NUFIT_GLOBAL_FITS.keys()``
lists what is available.

.. code-block:: python

   inverted = gd.load_nufit_params('NuFIT 6.1', 'IO')
   older = gd.load_nufit_params('NuFIT 5.2', 'NO', category='without_SK')

The mass ordering is carried by the **sign of** ``D31``, so the inverted set
differs from the normal one in that sign -- and, at the best fit, in the
:math:`\theta_{23}` octant and :math:`\delta_{\rm CP}` as well.  If you want to
vary the ordering alone, flip the sign of ``D31`` yourself rather than swapping
parameter sets; notebook 17 shows why.

Notebook 26 goes further and samples the :math:`\Delta\chi^2` profiles behind
these fits, to show how much of a predicted probability is really the
parameters.

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

   # Your own position-dependent Hamiltonian, returning a (d, d) array.
   # This one is the standard three-flavor vacuum term plus a matter
   # potential that falls off exponentially with position.
   h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
       osc['s12'], osc['s23'], osc['s13'], osc['dCP'], osc['D21'], osc['D31']))

   def H_func(l):
       vcc = 1.0e-13*np.exp(-l/(500.0*gd.UNIT_KM))        # [eV]
       return h_vac/energy + vcc*np.diag([1.0, 0.0, 0.0])

   P = oscprob.osc_prob(H_func, t_ini=0.0, t_fin=L,
                         magnus_exp_order=4,
                         rtol=1e-4, atol=1e-4)

Find a full worked example of using :func:`~magnus.oscprob.osc_prob` directly for a
time-dependent matrix exponential (not necessarily a physical Hamiltonian)
in notebook 11; see :doc:`tutorials`.
