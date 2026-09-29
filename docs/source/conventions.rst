.. _conventions:

Conventions
============

Everything below is a *choice*. None of it is forced by the physics, all of it
is forced by consistency, and a convention that is wrong **consistently** passes
every internal test — which is why they are written down here rather than left
in the code. Magνs has been bitten by exactly that: a reversed slab ordering, a
doubled antineutrino potential sign and a flipped two-flavor mass ordering were
all fixed on the same day, and each had been silently self-consistent.


At a glance
------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - What
     - Convention
   * - Probability matrix
     - ``P[i][f]`` is :math:`P(\nu_i \to \nu_f)`, initial flavor first; flavors in the
       order :math:`(e, \mu, \tau, s_1, s_2)`; a batched call puts the point index
       first, ``(n_points, d, d)``
   * - Mixing matrix
     - The PDG parametrization,
       :math:`U = R_{23}(\theta_{23})\,U_{13}(\theta_{13},\delta_{\rm CP})\,R_{12}(\theta_{12})`,
       with :math:`U_{e3} = \sin\theta_{13}\,e^{-i\delta_{\rm CP}}`
   * - Mixing angles
     - Given as sines by default; ``angles='sin2'``, ``'rad'`` or ``'deg'`` changes that
   * - CP phases
     - Radians (degrees under ``angles='deg'``)
   * - Mass splittings
     - ``D21`` :math:`= m_2^2 - m_1^2`, ``D31`` :math:`= m_3^2 - m_1^2`, in eV\ :sup:`2`;
       the ordering is the sign of ``D31``
   * - Inverted ordering from NuFIT
     - NuFIT quotes :math:`\Delta m^2_{32}`; ``load_nufit_params(..., 'IO')`` returns
       ``D31`` :math:`= \Delta m^2_{32} + \Delta m^2_{21}`
   * - Antineutrinos
     - ``nubar=True`` conjugates the mixing matrix and flips the sign of the matter
       potential
   * - Matter potential
     - :math:`V_{CC} = +\sqrt{2} G_F n_e` on the :math:`\nu_e` entry; :math:`Y_e = 0.5`
       by default, and one value per layer in the Earth (0.4656 core, 0.4957 mantle)
   * - Units
     - Natural units: energies in eV, distances in eV\ :sup:`-1`, densities in
       eV\ :sup:`4`; see :ref:`units-table`

.. _conv-hamiltonian:

The Hamiltonian and the probability
-----------------------------------

The state of :math:`n` flavors evolves under an :math:`n \times n` Hermitian
Hamiltonian :math:`\mathbb{H}`.  Its solution defines the evolution operator
:math:`\mathbb{U}(l_1, l_0)`, and the probability that a neutrino born as
:math:`\nu_\alpha` is detected as :math:`\nu_\beta` is

.. math::
   :label: conv-probability

   P_{\nu_\alpha \to \nu_\beta}(l_1, l_0)
   = \left\lvert \left[\mathbb{U}(l_1,l_0)\right]_{\beta\alpha} \right\rvert^2 .

In vacuum,

.. math::
   :label: conv-h-vacuum

   \mathbb{H}^{\rm vac}(E) = \frac{1}{2E}\, \mathbb{R}\, \mathbb{M}^2\, \mathbb{R}^\dagger ,
   \qquad \mathbb{M}^2 = {\rm diag}(0, \Delta m^2_{21}, \Delta m^2_{31}, \ldots) ,

with the mixing matrix :math:`\mathbb{R}` a product of rotations.  Each
:math:`\mathbb{R}_{ij}` rotates by :math:`\theta_{ij}` in the :math:`(i,j)` plane, with
:math:`\sin\theta_{ij}\,e^{-i\delta_{ij}}` in entry :math:`(i,j)`:

.. math::
   :label: conv-mixing-matrix

   \mathbb{R} =
   \begin{cases}
    \mathbb{R}_{12} , & n = 2 , \\
    \mathbb{R}_{23}\, \mathbb{R}_{13}\, \mathbb{R}_{12} , & n = 3 , \\
    \mathbb{R}_{34}\, \mathbb{R}_{24}\, \mathbb{R}_{14}\, \mathbb{R}_{23}\, \mathbb{R}_{13}\, \mathbb{R}_{12} , & n = 4 , \\
    \mathbb{R}_{35}\, \mathbb{R}_{25}\, \mathbb{R}_{15}\, \mathbb{R}_{34}\, \mathbb{R}_{24}\, \mathbb{R}_{14}\, \mathbb{R}_{23}\, \mathbb{R}_{13}\, \mathbb{R}_{12} , & n = 5 .
   \end{cases}

At three flavors this is the PMNS matrix, and the phase of :math:`\mathbb{R}_{13}` is
:math:`\delta_{\rm CP}`.  Of the sterile rotations only :math:`\mathbb{R}_{14}`,
:math:`\mathbb{R}_{24}`, :math:`\mathbb{R}_{15}` and :math:`\mathbb{R}_{35}` carry a phase.

In matter, the charged-current potential enters through a diagonal projector
:math:`\mathbb{P}`:

.. math::
   :label: conv-h-matter

   \mathbb{H}(E, l) = \mathbb{H}^{\rm vac}(E) + V_{\rm CC}(l)\, \mathbb{P}(l) ,
   \qquad V_{\rm CC} = \sqrt{2}\, G_F\, n_e(l) .

At three flavors :math:`\mathbb{P} = {\rm diag}(1, 0, 0)`.  With sterile states,
:math:`\mathbb{P} = {\rm diag}(1, 0, 0, r/2, \ldots)`, with :math:`r = n_n/n_p`
(``ratio_number_neutrons_to_protons``): the neutral-current potential, common to the
active flavors, is subtracted from the whole diagonal, and what it leaves on a sterile
state is :math:`(r/2)\,V_{\rm CC}`.  :func:`magnus.matter.matter_potential_projector`
builds :math:`\mathbb{P}` for every flavor count.

Ordering of the probabilities
-----------------------------

Every ``osc_prob_*`` function returns the probability matrix indexed
**initial flavor first**:

.. math::

   P[\nu_i][\nu_f] \;=\; P(\nu_i \to \nu_f) .

So ``P[1][0]`` is :math:`P(\nu_\mu \to \nu_e)`, not the reverse. Flavors are
in the standard order :math:`(e, \mu, \tau, s_1, s_2)`, so index 0 is always
:math:`\nu_e`.

Each **row** sums to one — a neutrino that started as :math:`\nu_i` ends as
something. Each column also sums to one, but that is a consequence of unitarity
rather than a separate statement. Passing ``nu_i`` and ``nu_f`` returns that one
entry instead of the matrix.

For a batched call the point index comes **first**: the shape is
``(n_points, d, d)``, so ``P[:, 1, 0]`` is :math:`P_{\mu e}` along a scan.

Sign of the matter potential
----------------------------

The charged-current potential enters the electron-flavor diagonal entry,

.. math::

   H \;=\; H_\text{vac} \;+\; \mathrm{diag}(V_{CC},\, 0,\, \ldots) ,
   \qquad V_{CC} = +\sqrt{2}\, G_F n_e ,

and **for antineutrinos it changes sign**. That flip is applied once, inside
:func:`magnus.matter.vcc_func_from_rho_func`, so a caller passing
``nubar=True`` gets it automatically and code downstream must not apply it
again. It was applied twice once, which gave antineutrinos a positive potential
and answers that looked plausible.

Mass ordering
-------------

The ordering is carried by the **sign of** :math:`\Delta m^2_{31}`, not by a
flag: positive is normal, negative is inverted. ``OSC_PARAMS_DEFAULT`` is the
normal ordering, with :math:`\Delta m^2_{31} = +2.511 \times 10^{-3}`
eV\ :sup:`2`. It is NuFIT 6.1 with Super-Kamiokande atmospheric data, the same
release :func:`~magnus.globaldefs.load_nufit_params` returns by default, and is
derived from it rather than written out a second time.
``magnus.globaldefs.OSC_PARAMS_PREDEFINED`` carries every NuFIT release from 1.0
on, in both orderings and, from 4.0 on, with and without that atmospheric data,
if you want to name the fit explicitly.

For two flavors the same rule applies to :math:`\Delta m^2`, which is what
makes the two-flavor case easy to get backwards: flipping its sign moves the
MSW resonance into the other channel, and the result is still a perfectly
ordinary-looking probability.

Mixing parameters
-----------------

Angles are given as **sines** by default (``angles='sin'``) -- not as angles,
and not as :math:`\sin^2\theta`; ``angles=`` also accepts ``'sin2'``,
``'rad'`` and ``'deg'``.  By default ``s12`` is
:math:`\sin\theta_{12}`. Quoted fits usually give :math:`\sin^2\theta`, so
take the square root — ``gd.S12_NO_BF_NUFIT_6_0`` is ``np.sqrt(0.308)``.
Phases are in **radians**; the default :math:`\delta_{CP}` is 3.7001 rad, i.e.
212 degrees.

Two flavors take ``sth`` and ``Dm2`` rather than ``s12`` and ``D21``. This is
one of the few convention errors here that cannot pass quietly: unrecognized
keywords are refused by name at the call site rather than forwarded down, so a
two-flavor call written with the three-flavor names raises instead of returning
a probability computed from the defaults.

Units
-----

Natural units throughout: energies in eV, baselines and positions in
eV\ :sup:`-1`, so that :math:`HL` is dimensionless.
:mod:`magnus.globaldefs` supplies the conversions — multiply by ``UNIT_KM``,
``UNIT_MEV``, ``UNIT_GEV``, ``UNIT_G_PER_CM3`` — and :ref:`units-table` lists
them.
