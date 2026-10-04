.. _conventions:

Conventions
============

These are the conventions Magνs follows: the parametrization, signs, units and
indexing that another code may choose differently.  A convention applied wrongly
but consistently passes every internal test, so each one is stated here.


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
       :math:`\mathbb{R} = \mathbb{R}_{23}(\theta_{23})\,\mathbb{R}_{13}(\theta_{13},\delta_{\rm CP})\,\mathbb{R}_{12}(\theta_{12})`,
       with :math:`\mathbb{R}_{e3} = \sin\theta_{13}\,e^{-i\delta_{\rm CP}}`
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
       potential.  A Lorentz-violating term is conjugated too, and changes sign when
       ``n_liv`` is even: an operator of dimension ``n_liv + 3`` is CPT-odd when that is
       odd, as in the Standard-Model Extension (Kostelecky & Mewes, Phys. Rev. D 85,
       096005 (2012))
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

   \mathbb{H} \;=\; \mathbb{H}^{\rm vac} \;+\; \mathrm{diag}(V_{CC},\, 0,\, \ldots) ,
   \qquad V_{CC} = +\sqrt{2}\, G_F n_e ,

and **for antineutrinos it changes sign**. That flip is applied once, inside
:func:`magnus.matter.vcc_func_from_rho_func`, so a user passing
``nubar=True`` gets it automatically, and code downstream must not apply it
again.  Applied twice, it would give antineutrinos a positive potential, and the
resulting probabilities would still look plausible.

Mass ordering
-------------

The ordering is carried by the **sign of** :math:`\Delta m^2_{31}`, not by a
flag: positive is normal, negative is inverted.
``OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']``, the default parameter set, is the
normal ordering, with :math:`\Delta m^2_{31} = +2.511 \times 10^{-3}`
eV\ :sup:`2`.  It is NuFIT 6.1 with Super-Kamiokande atmospheric data, the
release :func:`~magnus.globaldefs.load_nufit_params` returns by default.  To name
a fit explicitly, ``magnus.globaldefs.OSC_PARAMS_PREDEFINED`` holds every NuFIT
release from 1.0 on, in both orderings and, from 4.0 on, with and without that
atmospheric data.

For two flavors the same rule applies to :math:`\Delta m^2`, which is what
makes the two-flavor case easy to get backwards: flipping its sign moves the
MSW resonance into the other channel, and the result still looks like an
ordinary probability.

Mixing parameters
-----------------

Angles are given as **sines** by default (``angles='sin'``): ``s12`` is
:math:`\sin\theta_{12}`, not the angle and not :math:`\sin^2\theta_{12}`.
``angles='sin2'``, ``'rad'`` and ``'deg'`` select the other forms.  Under
``'rad'`` and ``'deg'``, an angle beyond :math:`\pm 90^\circ` is refused, because
the cosine is taken as :math:`+\sqrt{1 - s^2}`.  Fits usually quote
:math:`\sin^2\theta`; pass those with ``angles='sin2'``, or take the square root,
e.g., ``s12 = np.sqrt(0.308)``.
Phases are in **radians**; the default :math:`\delta_{CP}` is 3.7001 rad, i.e.
212 degrees.

Two flavors take ``sth`` and ``Dm2`` rather than ``s12`` and ``D21``.  A mistake
here cannot pass unnoticed: unrecognized keywords are refused by name, so a
two-flavor call written with the three-flavor names raises an error instead of
returning a probability computed from the defaults.

Units
-----

Natural units throughout: energies in eV, baselines and positions in
eV\ :sup:`-1`, so that :math:`\mathbb{H} L` is dimensionless.
:mod:`magnus.globaldefs` supplies the conversions — multiply by ``UNIT_KM``,
``UNIT_MEV``, ``UNIT_GEV``, ``UNIT_G_PER_CM3`` — and :ref:`units-table` lists
them.

.. _coming-from-other-codes:

Coming from GLoBES, Prob3++ or nuSQuIDS
----------------------------------------

The physics is the same; the bookkeeping differs.  Each entry for the other codes is
what the drivers behind :doc:`comparison` pass them (``resources/benchmarks/external_drivers``,
in a source checkout).

.. list-table::
   :header-rows: 1
   :widths: 16 21 21 21 21

   * -
     - Magνs
     - GLoBES
     - Prob3++
     - nuSQuIDS
   * - Energy
     - eV: multiply GeV by ``gd.UNIT_GEV``
     - GeV
     - GeV
     - its own units: ``E*units.GeV``, with ``units = nuSQuIDS.Const()``
   * - Baseline
     - eV\ :sup:`-1`: multiply km by ``gd.UNIT_KM``
     - km
     - km, or a zenith cosine through the Earth
     - its own units: ``L*units.km``
   * - Density
     - natural units, or g cm\ :sup:`-3` with ``density_matter_is_in_g_per_cm3=True``
     - g cm\ :sup:`-3`
     - g cm\ :sup:`-3`
     - g cm\ :sup:`-3`
   * - :math:`Y_e`
     - ``electron_fraction``, 0.5 by default; per layer through the Earth
     - fixed at 0.5 inside the code
     - from the profile file on the Earth path
     - an argument of the body, e.g., ``ConstantDensity(rho, 0.5)``
   * - Flavor indices
     - 0, 1, 2 = e, μ, τ; ``P[initial][final]``
     - 1, 2, 3; ``(initial, final, …)``
     - 1, 2, 3; ``GetProb(initial, final)``
     - 0, 1, 2; initial state set, then ``EvalFlavor(final)``
   * - Antineutrinos
     - ``nubar=True``
     - ``cp_sign = -1``
     - a negative neutrino type
     - ``NeutrinoType.antineutrino``
   * - Angles
     - sines by default; ``angles='sin2'``, ``'rad'`` or ``'deg'``
     - radians, in the order θ12, θ13, θ23
     - :math:`\sin^2\theta` (``kSquared=true``) or :math:`\sin^2 2\theta`, in the order
       θ12, θ13, θ23
     - radians, ``Set_MixingAngle(i, j, θ)``
   * - Splittings and ordering
     - ``D21``, ``D31``; the sign of ``D31`` is the ordering
     - Δm²21, Δm²31
     - Δm²21 and **Δm²32**, not Δm²31
     - ``Set_SquareMassDifference(1, Δm²21)``, ``(2, Δm²31)``
   * - :math:`\delta_{\rm CP}`
     - radians (degrees under ``angles='deg'``)
     - radians
     - radians
     - radians, ``Set_CPPhase(0, 2, δ)``

Two differences matter at the level of the comparison itself:

* **The matter potential for a given density.**  Every code converts g cm\ :sup:`-3` to an
  electron density as :math:`\rho N_A Y_e`, but with its own rounding of the constant
  :math:`\sqrt{2} G_F N_A`: Prob3++'s is 0.04% below it, and the comparison
  drivers rescale for it.  The stored benchmark data were measured with an older
  Magνs, whose potential was 0.8% lower; their rescaling factors (0.992 for GLoBES)
  include that difference (see the changelog for 1.2.0).
* **Prob3++ takes Δm²32.**  Passing Δm²31 in its place changes probabilities by up to 0.26.
