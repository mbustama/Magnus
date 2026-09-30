.. _ex-sec-constant-density:

Constant density
----------------

In matter of constant density, the Hamiltonian does not depend on position. Every commutator in the :ref:`recursion <magnus-recursion>` vanishes, the Magnus series stops at its first term, and the evolution operator is a single matrix exponential. Like for vacuum, Magνs computes this without exercising any of the method of :doc:`/methodology`.

.. _ex-lst-const-density:

**Constant density at three and two flavors.** The density is swept from vacuum to about the value at the center of the Earth. The last call returns 200 energies at one density, which the constant-Hamiltonian engine of :doc:`/engines` answers as a single stack of matrix exponentials. See :ref:`ex-sec-constant-density` for details.

.. code-block:: python

   import numpy as np
   import magnus.globaldefs as gd
   import magnus.oscprob as oscprob

   # 'NuFIT 6.1' is the default; re-loaded for clarity only
   osc = gd.load_nufit_params('NuFIT 6.1')
   E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM

   # Three flavors, one energy and one baseline
   for rho in (0.0, 3.0, 8.0, 13.0):
       P = oscprob.osc_prob_3nu_matter_constant_density(
           E, L, rho=rho,
           density_matter_is_in_g_per_cm3=True, **osc)
       print('%5.1f g/cm3  Pme = %.6f'
             % (rho, P[gd.NUMU][gd.NUE]))
   #   0.0 g/cm3  Pme = 0.003773
   #   3.0 g/cm3  Pme = 0.013475
   #   8.0 g/cm3  Pme = 0.049843
   #  13.0 g/cm3  Pme = 0.111807

   # Two flavors: one angle and one splitting
   P = oscprob.osc_prob_2nu_matter_constant_density(
       E, L, rho=3.0, sth=osc['s13'], Dm2=osc['D31'],
       density_matter_is_in_g_per_cm3=True)
   # P[0][1] = 0.005651

   # Two hundred energies in one call
   Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
   P = oscprob.osc_prob_3nu_matter_constant_density(
       Es, L, rho=3.0,
       density_matter_is_in_g_per_cm3=True, **osc)
   # P.shape = (200, 3, 3)


Listing :ref:`Constant density at three and two flavors <ex-lst-const-density>` shows the constant-density wrapper at three and two flavors, over a range of densities, and also the whole calculation as one batched call. The number passed as ``rho`` is read in one of three ways, set by two keywords, both ``False`` by default:

- ``density_matter_is_in_g_per_cm3=True``, as in Listing :ref:`Constant density at three and two flavors <ex-lst-const-density>`, means ``rho`` is mass density in g cm\ :math:`^{-3}`.

- ``density_is_of_number_of_electrons=True`` means ``rho`` is the electron number density itself, in eV\ :math:`^3`.

- With neither, ``rho`` is a mass density in natural units, eV\ :math:`^4`.

A mass density and an electron number density differ by the composition of the medium, which enters through a third keyword, ``electron_fraction``, the number of electrons per nucleon, :math:`Y_e = n_e/(n_p + n_n)`, with :math:`n_e`, :math:`n_p`, and :math:`n_n` the number densities of electrons, protons, and neutrons. It defaults to 1/2 in the constant-density wrappers. (Inside the Earth, the wrappers instead take :math:`Y_e` layer by layer from the PREM, as a function of the distance from the center; see :ref:`ex-sec-prem` and Figure :ref:`The Earth’s density profile <ex-fig-prem>`.) Lowering :math:`Y_e` to 0.2 at 3 g cm\ :math:`^{-3}` takes :math:`P_{\nu_\mu \to \nu_e}` from 0.013589 to 0.006726, so its effect is small, but not insignificant. It is read only when ``rho`` is a mass density; a call that supplies the electron number density directly ignores it, since you have already done the conversion it governs.

Four and five flavors take the same form, with the active-sterile parameters added as further keywords. The call ``osc_prob_4nu_matter_constant_density`` with ``s14`` and ``D41`` returns a :math:`4 \times 4` matrix where the three-flavor call returns a :math:`3 \times 3` one; nothing else about the call changes.

One further keyword, ``ratio_number_neutrons_to_protons``, is the ratio :math:`r = n_n/n_p` that enters the projector :math:`\mathbb{P}` (:ref:`conv-hamiltonian`); since :math:`n_p = n_e`, it is also the number of neutrons per electron, :math:`Y_n`. It defaults to 1, the isoscalar value, which is the medium the default :math:`Y_e = 1/2` already describes. In neutral matter, :math:`Y_e = 1/(1 + r)`. Setting one of the two keywords does not change the other, so both should be set together. The ratio sets the sterile entry of the matter projector :math:`\mathbb{P}` (:ref:`conv-hamiltonian`): with :math:`\sin^2\theta_{14} = \sin^2\theta_{24} = 0.1` and :math:`\Delta m^2_{41} = 1` eV\ :math:`^2`, at the energy, baseline, and density of Listing :ref:`Constant density at three and two flavors <ex-lst-const-density>`, raising it to 1.5 moves :math:`P_{\nu_\mu \to \nu_\mu}` from 0.7764 to 0.7587, by 2.3%. It does not enter the conversion of a mass density into an electron number density, :math:`n_e = \rho N_A Y_e`, which the electron fraction alone sets (see :ref:`conv-hamiltonian`).
