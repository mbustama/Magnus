.. _ex-sec-wrapper-vs-direct:

Using wrappers vs. scenario calls vs. direct calls
--------------------------------------------------

.. contents:: On this page
   :local:
   :depth: 1

.. _ex-lst-constant:

**Vacuum and constant density, three ways.** The three-flavor probability, computed three ways. The wrappers take the oscillation parameters as keywords, here the NuFIT 6.1 values, and the density in g cm\ :math:`^{-3}`. The scenario function takes the density as a function of position. The hand-built form is the vacuum term divided by the energy plus the potential :math:`V_{\rm CC}` times the projector :math:`\mathbb{P}` (:ref:`conv-hamiltonian`); ``VCC_EARTH_CRUST`` is that potential at 3 g cm\ :math:`^{-3}` and :math:`Y_e = 1/2`. All three return the same probabilities; the times are measured under the protocol of :doc:`/performance`.

.. code-block:: python

   import numpy as np
   import magnus.globaldefs as gd
   import magnus.hamiltonians as hamiltonians
   import magnus.matter as matter
   import magnus.oscprob as oscprob

   osc = gd.load_nufit_params('NuFIT 6.1')
   E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
   Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV

   # --- 1. With the named wrappers
   P_vac = oscprob.osc_prob_3nu_vacuum(E, L, **osc)
   P_mat = oscprob.osc_prob_3nu_matter_constant_density(
       E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc)
   # Pme = 0.00377306 in vacuum, 0.01358909 in
   # matter, in 0.090 ms

   P_scan = oscprob.osc_prob_3nu_matter_constant_density(
       Es, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc)
   # P_scan.shape = (200, 3, 3), in 0.123 ms

   # --- 2. With a scenario function
   rho_func = lambda l: 3.0*gd.UNIT_G_PER_CM3
   P_mat = oscprob.osc_prob_matter_std_potential(
       3, rho_func, E, L, osc_params=osc)
   # The same number, to 4e-16, in 31.9 ms

   # --- 3. With a Hamiltonian built by hand
   H_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
       s12=osc['s12'], s23=osc['s23'], s13=osc['s13'], dCP=osc['dCP'],
       D21=osc['D21'], D31=osc['D31'])
   proj = matter.matter_potential_projector(3)

   P_vac = oscprob.osc_prob(H_vac/E, 0.0, L)
   P_mat = oscprob.osc_prob(
       H_vac/E + gd.VCC_EARTH_CRUST*proj, 0.0, L)
   # The same two numbers, in 0.013 ms

   P_scan = oscprob.osc_prob_energy_baseline(
       lambda e: H_vac/e + gd.VCC_EARTH_CRUST*proj, Es, L,
       H_func_is_function_only_of_energy=True)
   # The same 200 points, in 3.6 ms


Every wrapper function (e.g., ``osc_prob_3nu_vacuum``, ``osc_prob_3nu_matter_constant_density``) builds a Hamiltonian internally and hands it, through a scenario function, to the engines of :doc:`/engines`, the last of which is ``osc_prob``, the primitive probability function of Magνs. The wrappers are a convenience; a user can do the same by hand, or enter at the scenario layer between the two. Wrappers also contain the most safeguards and validation particular to the cases they are built for; scenario functions contain fewer; and the direct call to ``osc_prob``, the least. :doc:`/architecture` shows the wrappers, scenario functions, and their relation to one another and ``osc_prob``.

Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>` computes the same three-flavor probability all three ways, at one energy and over a scan of two hundred. The three agree: to round-off in matter, bit-for-bit in vacuum.

.. _ex-sec-wrappers:

Through a wrapper
~~~~~~~~~~~~~~~~~

.. _ex-tab-wrappers:

.. table:: The probability wrapper functions

   .. list-table::
      :header-rows: 1
      :widths: 22 22 56

      * - Environment
        - Physics
        - Wrapper, for N = 2, 3, 4, 5
      * - Vacuum
        - Standard
        - ``osc_prob_Nnu_vacuum``
      * -  
        - Lorentz violation
        - ``osc_prob_Nnu_vacuum_liv``
      * - Constant density
        - Standard
        - ``osc_prob_Nnu_matter_constant_density``
      * -  
        - Non-standard int.
        - ``osc_prob_Nnu_matter_nsi_constant_density``
      * -  
        - Lorentz violation
        - ``osc_prob_Nnu_matter_liv_constant_density``
      * - Exponential density
        - Standard
        - ``osc_prob_Nnu_matter_exp_density``
      * -  
        - Non-standard int.
        - ``osc_prob_Nnu_matter_nsi_exp_density``
      * -  
        - Lorentz violation
        - ``osc_prob_Nnu_matter_liv_exp_density``
      * - Earth (PREM)
        - Standard
        - ``osc_prob_Nnu_earth``
      * -  
        - Non-standard int.
        - ``osc_prob_Nnu_earth_nsi``
      * -  
        - Lorentz violation
        - ``osc_prob_Nnu_earth_liv``
      * - Sun
        - Standard
        - ``osc_prob_Nnu_sun``
      * -  
        - Non-standard int.
        - ``osc_prob_Nnu_sun_nsi``
      * -  
        - Lorentz violation
        - ``osc_prob_Nnu_sun_liv``

The 56 wrapper functions of Magνs, by environment and by the physics in the Hamiltonian. Each of the fourteen combinations exists at two, three, four, and five flavors. There is no wrapper for non-standard interactions in vacuum, since those are interactions with matter. Every function here builds a Hamiltonian and hands it to the engines of :doc:`/engines`. Three others, not shown here, do not: ``osc_prob_2nu_vacuum_std``, ``osc_prob_3nu_vacuum_std``, and ``osc_prob_2nu_matter_std`` evaluate the standard analytic expressions for the probability instead. One more wrapper, ``osc_prob_pseudo_dirac_vacuum``, covers a pseudo-Dirac spectrum in vacuum (:ref:`ex-sec-astro`).

Table :ref:`ex-tab-wrappers` lists the 56 wrappers of Magνs. A wrapper is the shortest call whenever the scenario is one Magνs ships with. It takes the oscillation parameters as keywords, fills in any left unset from a named set, converts the density from g cm\ :math:`^{-3}`, and assembles the Hamiltonian. That last step includes the vacuum term: ``osc_prob_3nu_matter_constant_density`` takes the same six oscillation parameters as ``osc_prob_3nu_vacuum`` and adds the potential to what it builds from them, which is why ``H_vac`` appears only in the third block of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`.

Every wrapper takes its arguments in the same four groups: the energy, the geometry of the environment, the oscillation parameters, and whatever the scenario adds. For instance, for non-standard interactions (NSI) at three flavors in constant density,

.. code-block:: python

   P = oscprob.\
    osc_prob_3nu_matter_nsi_constant_density(
       E, L,            # 1. Energy, baseline
       rho=3.0,         # 2. Density
       density_matter_is_in_g_per_cm3=True,
       **osc,           # 3. Mixing
       eps_em=0.05)     # 4. New physics
   # Pme = 0.017126.  Dropping eps_em gives
   # 0.013589, the standard wrapper's answer


The fourth group (“new physics”) exists only in the wrappers whose name carries ``_nsi`` or ``_liv``; Table :ref:`ex-tab-parameters` lists its parameters. Every one of them defaults to zero, so one of those wrappers called without them returns the standard result, bit-for-bit identical to what the standard wrapper of the same environment returns. That makes it the control for any new-physics scan.

Passing ``rho`` as a scalar makes the Hamiltonian independent of position, so the request reaches the constant-Hamiltonian engine of :doc:`/engines` and the whole scan becomes one stack of matrix exponentials. Two hundred energies cost 0.123 ms vs. 0.090 ms for one, so the scan is nearly free.

Through a scenario function
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _ex-tab-scenarios:

.. table:: The probability scenario functions

   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+
   | Scenario function                 | Covers                                 | Arguments after ``num_flavors``                                           |
   +===================================+========================================+===========================================================================+
   | ``osc_prob_vacuum``               | Vacuum                                 | ``energy``, ``L``, ``osc_params``                                         |
   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+
   | ``osc_prob_matter_std_potential`` | Standard oscillations in matter        | ``rho_func``, ``energy``, ``L``, ``osc_params``                           |
   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+
   | ``osc_prob_matter_nsi``           | Non-standard interactions              | ``rho_func``, ``energy``, ``L``, ``osc_params``, ``nsi_params``           |
   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+
   | ``osc_prob_liv``                  | Lorentz violation, in vacuum or matter | ``energy``, ``L``, ``osc_params``, ``liv_params``; ``rho_func`` optional  |
   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+
   | ``osc_prob_earth``                | A chord through PREM                   | ``H_func``, ``energy``, and ``costhz`` with ``L``, or two named locations |
   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+
   | ``osc_prob_sun``                  | A ray out of the Sun                   | ``H_func``, ``energy``, ``L``                                             |
   +-----------------------------------+----------------------------------------+---------------------------------------------------------------------------+

The scenario functions of Magνs, one per physics scenario. Each takes the flavor count as its first argument, so one function serves every flavor count where Table :ref:`ex-tab-wrappers` needs four wrappers. Magνs ships vacuum Hamiltonians for two through five flavors; past that you supply one through ``h_vac_energy_indep``, which the same functions then use. Each of the matter functions takes the density as ``rho_func``, a function of position, which is what makes this the layer for a profile that varies. The optional ``rho_func`` of ``osc_prob_liv`` is why Lorentz violation is the one non-standard scenario that also exists in vacuum. The refinement keywords of :doc:`/architecture` reach all of them through ``**kwargs``. The last two, ``osc_prob_earth`` and ``osc_prob_sun``, are not scenario functions: they take a Hamiltonian you supply and apply a geometry to it, so they pair either with a scenario function or with a Hamiltonian built by hand.

Table :ref:`ex-tab-scenarios` lists the four scenario functions shipped with Magνs (plus two other functions at the same level that handle geometry in the Earth and the Sun). The scenario functions sit one layer below the wrappers, one per physics scenario and each applicable at any flavor count. They take the density as a *function of position*, which is what makes them the entry point for a profile that varies; :doc:`/comparison` calls Magνs this way.

Every scenario function takes its arguments in the same groups, with the flavor count first and the physics handed over as dictionaries rather than as keywords. For instance, for NSI at three flavors in constant density,

.. code-block:: python

   rho_func = lambda l: 3.0*gd.UNIT_G_PER_CM3
   nsi = dict(eps_ee=0.0, eps_em=0.05,
              eps_et=0.0, eps_mm=0.0,
              eps_mt=0.0, eps_tt=0.0)

   P = oscprob.osc_prob_matter_nsi(
       3,               # 1. Flavor count
       rho_func,        # 2. The medium
       E, L,            # 3. Energy, baseline
       osc_params=osc,  # 4. Mixing
       nsi_params=nsi)  # 5. New physics
   # Pme = 0.017126, the same number the
   # wrapper above returns


Two features differ from the wrapper. First, the flavor count is an argument, not part of the function’s name, which is why one scenario function covers what Table :ref:`ex-tab-wrappers` needs four wrappers for. Second, the dictionaries have no defaults: every coupling has to be present, so the six entries of ``nsi_params`` are all written out even though five of them are zero. A wrapper fills those in.

In constant density, scenario functions are the wrong tool. A callable ``rho_func`` that returns the same number everywhere is not detected as constant by the scenario function, so the constant-Hamiltonian engine is never reached: a single probability falls through to the general machinery in ``osc_prob``, and a scan to the energy-batched engine. The answer is the same as the wrapper’s to round-off, at a far higher cost.

Through a direct call
~~~~~~~~~~~~~~~~~~~~~

Handing ``osc_prob`` a matrix is the only route on which you assemble the Hamiltonian yourself; the wrappers and the scenario functions build it from the parameters they are given. At the same time, it is the cheapest route per probability, about 0.01 ms, because the parameter resolution, input validation, and assembly have already been paid for by you.

The wrapper and scenario function above both reach the probability from named parameters, building the Hamiltonian internally on the way. At the base layer, you build it yourself. For the same case of NSI at three flavors in constant-density matter, this means computing and summing three matrices, i.e.,

.. code-block:: python

   H_vac = hamiltonians.hamiltonian_3nu_vacuum(
       E,
       s12=osc['s12'], s23=osc['s23'],
       s13=osc['s13'], dCP=osc['dCP'],
       D21=osc['D21'], D31=osc['D31'])

   # proj for three flavors is diag(1, 0, 0)
   proj = matter.matter_potential_projector(3)

   V = gd.VCC_EARTH_CRUST # About 1.145e-13 eV
   H_nsi = hamiltonians.hamiltonian_3nu_nsi(
       V, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0)

   # 1. mixing, 2. matter, 3. new physics
   H = H_vac + V*proj + H_nsi

   P = oscprob.osc_prob(H, 0.0, L)   # 4. path
   # Pme = 0.017126 again


``hamiltonian_3nu_nsi`` returns the non-standard term alone, :math:`V_{\rm CC}` times the matrix of couplings, Eq. :eq:`ex-equ-h-nsi-3nu`, without the standard charged-current term. The standard term is a separate summand, which means that leaving ``V*proj`` out of the sum does not raise an error: it returns a converged, unitary probability of 0.005962 for a medium with no ordinary matter effect in it. This is unphysical, since NSI coexist with standard interactions; Preventing this is up to you.

Assembling the Hamiltonian rarely means writing a matrix from scratch. Magνs ships forty-three Hamiltonian builders: vacuum, matter, non-standard interactions, and Lorentz violation at two to five flavors, in position-dependent and position-independent forms, and three for a pseudo-Dirac spectrum. Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>` uses one of them for the vacuum term and takes the projector from ``matter``. :ref:`ex-sec-hamiltonians` lists them all, so a direct call is usually a shipped Hamiltonian with something you added to it, e.g., a new non-standard contribution.

The direct route does not batch: ``osc_prob`` takes one Hamiltonian at one energy and one baseline, and returns one matrix of probabilities. To scan energies with a Hamiltonian built by hand, pass it to ``osc_prob_energy_baseline``, the third layer of :doc:`/architecture`, which calls ``osc_prob`` once per energy. For the 200 energies of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, that takes 3.6 ms. The wrapper computes the same 200 energies in 0.123 ms, about thirty times faster, because it hands them to the constant-Hamiltonian engine as one batch.

Calling ``osc_prob`` is the most general route: replacing the scalar ``V`` by a function of position gives a varying profile, as in :ref:`ex-sec-hamiltonians`; replacing the whole Hamiltonian matrix is illustrated in :ref:`ex-sec-building-new-hamiltonian` and :ref:`ex-sec-lri-sun`.

Which to use
~~~~~~~~~~~~

In brief, a user should use:

- **A wrapper** (Table :ref:`ex-tab-wrappers`), for any scenario Magνs already implements (vacuum, constant density, exponential density, the Earth, and the Sun, each at two to five flavors, with standard oscillations, non-standard interactions, or Lorentz violation). It is the shortest call and the one with the most safeguards: it fills in the parameters left unset, converts the density, and builds the Hamiltonian in the form the fastest applicable engine of :doc:`/engines` expects, e.g., a constant density as a number.

- **A scenario function** (Table :ref:`ex-tab-scenarios`), for a density profile that no wrapper provides, given as any function of position, or for more than five flavors, with a vacuum Hamiltonian supplied through ``h_vac_energy_indep``. Given the same inputs, it reaches the same engines as a wrapper.

- **A direct call** to ``osc_prob``, for a Hamiltonian that no shipped builder produces: a term from new physics beyond non-standard interactions and Lorentz violation, or a matrix that comes from somewhere else entirely.

.. _ex-sec-evolution-operator:

Returning the evolution operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every wrapper and scenario function above returns probabilities; for most purposes, that is the output wanted. But some observables are built from the amplitudes instead. Three examples are the mass-eigenstate content of a neutrino leaving a dense source, which is the factor the source contributes to :ref:`averaged limit on a varying profile <avg-varying>`; the flavor composition at a detector after an incoherent mixture of mass eigenstates has crossed the Earth; and the evolution through two media in sequence. None of them can be recovered from a probability matrix, since the phases are gone at that stage. For those cases, every probability function of Magνs, on any of the three routes above and at any flavor count, returns the evolution operator alongside the probabilities when asked via ``return_evolution_operator``:

.. code-block:: python

   P, U = oscprob.osc_prob_3nu_vacuum(
       E, L, return_evolution_operator=True,
       **osc)
   # P is the matrix the call returns without
   # the keyword; U is complex and unitary.
   # Both of these hold:
   np.allclose(P, abs(U.T)**2)
   np.allclose(U.conj().T @ U, np.eye(3))


The operator is indexed ``U[final, initial]``, so that its moduli squared, transposed, are the probability matrix, indexed the other way. Everything else keeps its meaning: ``nu_i`` and ``nu_f`` still select one channel of ``P`` while ``U`` stays the full operator; arrays of energies or baselines return one operator per point, of shape ``(n, d, d)``. ``return_evolution_operator`` is declared by ``osc_prob`` and forwarded by every wrapper, as ``t_breakpoints`` is.

Inside ``osc_prob``, the request changes two behaviors. First, the refinement ladder compares the operator itself between levels, at the same ``rtol`` and ``atol``, so the operator returned is converged in its phases, not only in its moduli. Second, the specialized engines of :doc:`/engines` stand aside for the call, since only the general Magnus ladder forms the evolution operator; a scan over baselines takes the per-point path instead of the cumulative traversal. Two settings are refused if asked together with ``return_evolution_operator``, with an error returned naming why: ``average=True``, since its routes decohere before any operator is formed; ``strategy='hybrid'``, since that engine answers with probabilities only.

The evolution operators returned can be composed. The evolution through two media in sequence is the product of their operators, the later one on the left:

.. code-block:: python

   # Two slabs of constant density in sequence,
   # the second denser than the first
   H1 = H_vac/E + gd.VCC_EARTH_CRUST*proj
   H2 = H_vac/E + 1.5*gd.VCC_EARTH_CRUST*proj
   _, U1 = oscprob.osc_prob(H1, 0.0, L,
       return_evolution_operator=True)
   _, U2 = oscprob.osc_prob(H2, 0.0, L,
       return_evolution_operator=True)
   U = U2 @ U1
   P = abs(U.T)**2       # Over both slabs, 2L


The same probabilities come from one call over a two-slab profile with the edge declared through ``t_breakpoints``, as in :doc:`/recipes`; the product above is what that call assembles internally. :ref:`ex-sec-jet` uses the ``return_evolution_operator`` flag and operator composition to first compute the mass-eigenstate content of the neutrino flux that leaves a star, and then its flavor content at a detector on Earth.

.. _ex-sec-average-keyword:

Asking for the averaged probability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every wrapper, scenario function, and ``osc_prob_energy_baseline`` return the phase average of :doc:`/averaged_probability` when called with ``average=True``. Two more keywords configure it: ``average_spread`` sets the relative energy spread, :math:`\sigma`, which is 10% by default, and ``average_initial_state`` sets the state the neutrino starts in, explained below. ``osc_prob`` itself, which computes a single point, does not accept ``average``. Nor can ``average=True`` be combined with ``return_evolution_operator`` (:ref:`ex-sec-evolution-operator`), since the averaged routes never form an evolution operator. On the wrappers and scenario functions, ``strategy_info['engine']`` reports ``'average'`` when the averaged route answered (:ref:`ex-sec-which-engine`).

The phase average is not the mean of the probability over an energy range: it damps each interference term by the spread of its phase, but keeps the eigenvectors at the central energy (:doc:`/averaged_probability`). Wherever the limit (:ref:`constant <avg-limit>` or :ref:`varying <avg-varying>` Hamiltonian) agrees with the phase average within :math:`10^{-4}`, or within the tighter of ``rtol`` and ``atol`` if that is smaller, both absolutely and relative to the probability, Magνs returns the limit. How the phase average is computed follows :doc:`/averaged_probability`.

*A constant Hamiltonian.—*\ For a Hamiltonian that does not depend on position, the phase average is the :ref:`averaged limit <avg-limit>`. It needs only the eigenvectors of :math:`\mathbb{H}`: one diagonalization per energy, no propagation. Over an astrophysical baseline, every interference term is damped, and the result is the :ref:`averaged limit <avg-limit>`. In the context of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`,

.. code-block:: python

   FAR = 1.0e8*gd.UNIT_KM   # Astrophysical L
   P = oscprob.osc_prob_3nu_vacuum(
       E, FAR, average=True, **osc)
   # P[0, 0] = 0.5481, P[1, 0] = 0.2155: the
   # familiar sum_i |U_ei|^2 |U_mi|^2, from the
   # eigenvectors alone


At the 1 000 km of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, the atmospheric phase is only about 6 rad, and a 10% spread damps its interference only partly. The result then depends on :math:`\sigma`, and Magνs warns about it:

.. code-block:: python

   P = oscprob.osc_prob_3nu_vacuum(
       E, L, average=True, **osc)   # 1000 km
   # PhaseAveragingWarning: the phase-averaged
   # probability depends on the energy spread
   # at 1 of 1 (energy, L) point(s) ...
   # P[0, 0] = 0.9844, P[1, 1] = 0.9088;
   # without the keyword: 0.9924, 0.9955


A pair of eigenstates whose phase stays small over the baseline, as in the coherent blocks of the :ref:`block form <avg-coherence>`, keeps its interference almost undamped, because a small phase also changes little across the spread.

*The initial state.—*\ For a Hamiltonian that varies along the path, the result depends on the state the neutrino starts in (:ref:`avg-initial-state`). With ``average_initial_state='flavor'``, the default, the neutrino starts as a flavor state, as one produced in the medium does, e.g., in a beam, in the atmosphere, in the Earth, or in the Sun. With ``'decohered'``, it starts as an incoherent mixture of the eigenstates at production, as one arriving from a distant source does. The two differ only by the interference present at production, so they agree wherever that interference is damped away, as it is for a neutrino made in the solar core:

.. code-block:: python

   for start in ('flavor', 'decohered'):
       P = oscprob.osc_prob_3nu_sun(
           10.0*gd.UNIT_MEV,
           gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
           nu_i=gd.NUE, nu_f=gd.NUE,
           average=True,
           average_initial_state=start)
       print(start, P)    # 0.297082 both times


Where it is not damped away, the two differ; on the chord of :ref:`ex-sec-solar-tomography`, they differ by 0.25. For a constant Hamiltonian, ``'decohered'`` returns the :ref:`averaged limit <avg-limit>`.

*Route 1: a smoothly varying Hamiltonian.—*\ This is the case of the Sun. Magνs carries the initial state along the instantaneous eigenstates of the Hamiltonian, applies the Magnus patch of :doc:`/adiabatic_strategy` across every crossing that is not adiabatic, and reads the result in the flavor basis at detection, as in :doc:`/methodology`. On a path without such a crossing, the result is the :ref:`averaged limit on a varying profile <avg-varying>` for a flavor start, and its decohered form for a decohered start. As in other calls, ``rtol`` and ``atol`` set the tolerance: the patches, and propagation along the adiabatic stretches between them, converge to the tighter of the two, :math:`10^{-3}` by default; if it is tighter than :math:`10^{-4}`, it replaces :math:`10^{-4}` in the choice between the limit and the phase average. Over chords through the core with impact parameters below :math:`0.15\,R_\odot`, from 30 GeV to 3 TeV, the probability at the default tolerance differs from its value at a tolerance of :math:`10^{-5}` by at most :math:`10^{-4}`. A discontinuity that is not declared, but could change the probability, is reported with ``UnmarkedDiscontinuityWarning``, not treated as smooth.

*Route 2: a profile with declared discontinuities.—*\ This is the case of the Earth, whose PREM layer boundaries are breakpoints, or of any profile with discontinuities declared through ``t_breakpoints`` or ``t_slab_edges``. Across a jump in the density, the eigenstates change abruptly, so there are none to carry the state along. Magνs instead propagates at energies across a window, averages the resulting probabilities, and warns that the result is the mean of the probability over that window, not the phase average. This route starts from a flavor state only; ``'decohered'`` raises an error.

*A Hamiltonian built by hand.—*\ Because ``osc_prob_energy_baseline`` accepts ``average``, a Hamiltonian built by hand is averaged in the same way as a shipped one:

.. code-block:: python

   H = H_vac/E + gd.VCC_EARTH_CRUST*proj
   P = oscprob.osc_prob_energy_baseline(
       H, E, FAR, average=True)
   # P[0, 0] = 0.8934: the decohered limit,
   # now from the eigenvectors in matter


*Where ``average`` is used.—*\ Several later sections use ``average=True``. For instance, in the Sun (:ref:`ex-sec-sun` and Listing :ref:`Averaged solar probabilities <ex-lst-sun>`), a ray from the center to the surface holds about :math:`1.4 \times 10^{5}` oscillation cycles at 5 MeV (:ref:`ex-sec-averaging-is-not-estimating`), far too many to resolve. For the flavor composition of TeV–PeV astrophysical fluxes (:ref:`ex-sec-astro`), the :ref:`averaged limit <avg-limit>` is the whole observable. For a long-range force in the Sun (:ref:`ex-sec-lri-sun`), the Hamiltonian is built by hand. :ref:`ex-sec-jet` does not use ``average=True``. There, the neutrino crosses the star coherently and decoheres only on its way to Earth, so the average is taken after it leaves the star, from the evolution operator across the star.
