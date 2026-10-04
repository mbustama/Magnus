.. _ex-sec-batched-calls:

Batched calls
-------------

Every wrapper and every scenario function accepts an array for the energy, for the baseline, or for both. The other arguments, such as the mixing parameters and the density, take a single value per call. An array of energies at one baseline returns one probability matrix per energy, stacked along an axis. The snippets on this page continue Listing :ref:`A minimal probability calculation <ex-lst-minimal>`, with its imports, ``E``, and ``L``:

.. code-block:: python

   import numpy as np

   Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
   P = oscprob.osc_prob_3nu_vacuum(Es, L)
   P.shape                    # (200, 3, 3)


An array of baselines at one energy works the same way, here at the 1 GeV of ``E``:

.. code-block:: python

   Ls = np.linspace(0, 1300, 200)*gd.UNIT_KM
   P = oscprob.osc_prob_3nu_vacuum(E, Ls)
   P.shape                    # (200, 3, 3)
   # P[-1] equals the single call at 1300 km


Arrays for both the energy and the baseline are paired element by element, not crossed into a grid, so they must have the same length; if they do not, Magνs raises a ``ValueError``. The first two pairs below have the same :math:`L/E`, so, in vacuum, they have the same probabilities, e.g., 0.048968 for :math:`\nu_\mu \to \nu_e`:

.. code-block:: python

   Es = np.array([1.0, 2.0, 5.0])*gd.UNIT_GEV
   Ls = np.array([500, 1000, 1300])*gd.UNIT_KM
   P = oscprob.osc_prob_3nu_vacuum(Es, Ls)
   P.shape          # (3, 3, 3): 3 pairs
   P[0, 1, 0], P[1, 1, 0]   # 0.048968 both


To compute a grid of energies and baselines, such as an oscillogram, make one call per baseline, each with the full array of energies.

A batched call is much faster than the same points computed one at a time, because Magνs shares the work across them (:doc:`/engines`). Where the Hamiltonian does not vary along the path, as in vacuum or at constant density, the whole scan is one stack of matrix exponentials. Where it varies, the energy-batched engine samples the density profile once for all the energies of a scan, and the cumulative engine computes all the baselines of a scan in a single pass. For scans of 200 points, at two to five flavors, batching is 77–139 times faster than one call per point in vacuum, 108–130 times faster at constant density, and 29–175 times faster along an Earth chord. The gain is smallest where the energies need the most different refinements, as for an eV-scale sterile neutrino along an Earth chord: :math:`1.7 \times`, for 40 energies at four flavors.

In vacuum and at constant density, batched and point-by-point results agree to within :math:`10^{-13}`, and exactly at two and three flavors. Along a varying profile, each route refines its answer to the requested tolerance on its own, so the two agree to within about that tolerance, but not exactly.
