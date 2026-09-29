.. _ex-sec-minimal:

A minimal example in vacuum, step-by-step
-----------------------------------------

.. _ex-lst-minimal:

**A minimal probability calculation.** The simplest complete calculation: the three-flavor vacuum probabilities for a 1-GeV neutrino over 1 300 km. The oscillation parameters—here and in the rest of :ref:`ex-sec-examples`— take their defaults, the NuFIT 6.1 best fit with Super-Kamiokande atmospheric data in normal ordering  :cite:p:`Esteban:2024eli`. See :ref:`ex-sec-minimal` for details.

.. code-block:: python

   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

   # Natural units: energies in eV, lengths in eV^-1
   E = 1.0*gd.UNIT_GEV
   L = 1300.0*gd.UNIT_KM

   P = oscprob.osc_prob_3nu_vacuum(E, L)

   print('Pme = %.5f, Pmm = %.5f, Pmt = %.5f'
    % (P[gd.NUMU][gd.NUE],
    P[gd.NUMU][gd.NUMU],
    P[gd.NUMU][gd.NUTAU]))

   # Pme = 0.03127, Pmm = 0.39231, Pmt = 0.57643


Listing :ref:`A minimal probability calculation <ex-lst-minimal>` is a complete calculation, built in four steps.

#. **Import the two modules a typical calculation needs.** ``oscprob`` holds the probability functions; ``globaldefs`` holds the physical constants, the unit conversions, and the flavor indices.

   .. code-block:: python

      import magnus.oscprob as oscprob
      import magnus.globaldefs as gd


#. **Give the energy and the baseline in natural units.** Every energy crossing the interface is in eV and every length in eV\ :math:`^{-1}` (:doc:`/conventions`). The factors ``UNIT_GEV`` and ``UNIT_KM`` convert from the units a user thinks in.

   .. code-block:: python

      E = 1.0*gd.UNIT_GEV
      L = 1300.0*gd.UNIT_KM


#. **Call one wrapper.** Its name fixes the flavor count and the environment. The oscillation parameters default to the NuFIT 6.1 best fit in normal ordering; any of them may be overridden by name, e.g., ``dCP=0.0``. Nothing in the call selects an engine or a numerical setting: Magνs chooses the engine and, where the profile varies, the slab count. [They can still be fixed via arguments (:doc:`/engines` and :doc:`/methodology`).]

   .. code-block:: python

      P = oscprob.osc_prob_3nu_vacuum(E, L)


#. **Read the matrix.** The call returns every flavor pair at once rather than one number, indexed by the constants ``NUE``, ``NUMU``, and ``NUTAU``.

   .. code-block:: python

      P[gd.NUMU][gd.NUE]


The returned matrix is indexed initial flavor first, :math:`P[\alpha][\beta] = P_{\nu_\alpha \to \nu_\beta}`. Each row and each column sums to one, because the evolution is unitary. The corresponding evolution operator element, which the call returns when asked (:ref:`ex-sec-evolution-operator`), is :math:`\mathbb{U}_{\beta\alpha}`, so :math:`\lvert \mathbb{U}\rvert^2` comes out indexed the other way:

.. code-block:: python

   # P is indexed initial flavor first
   P[gd.NUMU][gd.NUE]           # 0.03127

   # The amplitude is indexed the other way
   P, U = oscprob.osc_prob_3nu_vacuum(
       E, L, return_evolution_operator=True)
   abs(U[gd.NUE][gd.NUMU])**2   # 0.03127
   abs(U[gd.NUMU][gd.NUE])**2   # 0.00854


In vacuum, the Hamiltonian does not depend on position, so the Magnus series stops at its first term: the evolution operator is a single exponential, exact to round-off. Therefore, Listing :ref:`A minimal probability calculation <ex-lst-minimal>` exercises none of the method of :doc:`/methodology`. Everything else on this page is the same four steps with a different wrapper. Once the profile varies, the refinement of :doc:`/methodology` engages, at the default ``rtol`` :math:`=` ``atol`` :math:`= 10^{-3}` unless another tolerance is requested.
