.. _ex-sec-hamiltonians:

The Hamiltonians
----------------

.. contents:: On this page
   :local:
   :depth: 1

A direct call to ``osc_prob`` takes the Hamiltonian as a matrix, or as a function of position that returns one (:ref:`ex-sec-wrapper-vs-direct`). Magνs ships builders for its standard and non-standard Hamiltonian terms (:ref:`ex-sec-shipped-hamiltonians`), and any other term can be written by hand (:ref:`ex-sec-building-new-hamiltonian`).

.. _ex-tab-hamiltonians:

.. table:: The Hamiltonian builders

   .. list-table::
      :header-rows: 1
      :widths: 40 60

      * - Returns
        - ``hamiltonian_`` + (N = 2, 3, 4, 5 where it appears)
      * - Vacuum term
        - ``Nnu_vacuum``
      * - Same, without the :math:`1/E`
        - ``Nnu_vacuum_energy_independent``
      * - Vacuum term at a position
        - ``Nnu_vacuum_td``
      * - Both of the above
        - ``Nnu_vacuum_energy_independent_td``
      * - :math:`V_{\rm CC}` times the projector
        - ``Nnu_matter``
      * - Same, from a potential function
        - ``Nnu_matter_td``
      * - Non-standard term alone
        - ``Nnu_nsi``
      * - Same, from a potential function
        - ``Nnu_nsi_td``
      * - Lorentz-violating term alone
        - ``Nnu_liv``
      * - Same, without the energy factor
        - ``Nnu_liv_energy_independent``
      * - Pseudo-Dirac vacuum term
        - ``pseudo_dirac_vacuum``
      * - Same, without the :math:`1/E`
        - ``pseudo_dirac_vacuum_energy_independent``
      * - Its matter term
        - ``pseudo_dirac_matter``

The :math:`43` Hamiltonian builders of ``magnus.hamiltonians``, each named with the prefix ``hamiltonian_``. The ten variants above the table break exist at two, three, four, and five flavors. The three below the break build a pseudo-Dirac spectrum from the active one, and take the pairs of states to split, so their names carry no number of flavors.

.. _ex-sec-shipped-hamiltonians:

Shipped Hamiltonians
~~~~~~~~~~~~~~~~~~~~

Table :ref:`ex-tab-hamiltonians` lists the 43 Hamiltonian builders shipped in ``magnus.hamiltonians``: vacuum, matter, non-standard interactions, and Lorentz violation at two to five flavors, plus three for a pseudo-Dirac spectrum.

A Hamiltonian is a sum of terms, as in the :ref:`matter Hamiltonian <conv-hamiltonian>`. Only the vacuum builders return a complete Hamiltonian; every other builder in Table :ref:`ex-tab-hamiltonians` returns its own term alone. For instance, at three flavors, the matter builder returns :math:`V_{\rm CC}` times the projector, which is zero except in the electron entry, ``[0][0]``. Added to a vacuum term, it gives the :ref:`matter Hamiltonian <conv-hamiltonian>`. Passed to ``osc_prob`` on its own, it returns the probabilities as the identity matrix, since a diagonal Hamiltonian does not mix flavors.

Every builder returns a :math:`d \times d` ``NumPy`` array, with :math:`d` the number of flavors, so terms add with ``+`` and scale with ``*``. The builders take their parameters under the same names as the wrappers, so a three-flavor set loaded with ``load_nufit_params`` can be passed to them with ``**``. Every physical parameter must be given, since a builder has no default set to fall back on. The one exception is the neutron-to-proton ratio of the four- and five-flavor matter builders, which is 1 by default.

The snippets on this page reuse ``E``, ``L``, ``Es``, and ``osc`` from Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, and import the builders as ``ham``.

*Vacuum.—*\ At three flavors, the vacuum Hamiltonian is

.. code-block:: python

   import magnus.hamiltonians as ham

   H = ham.hamiltonian_3nu_vacuum(E, **osc)
   H.shape                        # (3, 3)


The call returns the :ref:`vacuum Hamiltonian <conv-hamiltonian>` at three flavors, i.e.,

.. math::
   :label: ex-equ-h-vacuum-3nu

   \mathbb{H}_{3\nu}^{\rm vac}(E)
   =
   \frac{1}{2E}\,
   \mathbb{R}\,
   {\rm diag}\!\left(0,\, \Delta m^2_{21},\, \Delta m^2_{31}\right)
   \mathbb{R}^\dagger \;,

with :math:`\mathbb{R}` the PMNS matrix (:ref:`conv-hamiltonian`), built from :math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}`, and :math:`\delta_{\rm CP}`. At two flavors, the builder takes ``sth`` and ``Dm2``, and returns the :ref:`vacuum Hamiltonian <conv-hamiltonian>` with :math:`\mathbb{R} = \mathbb{R}_{12}` and :math:`\mathbb{M}^2 = {\rm diag}(0, \Delta m^2)`.

At four and five flavors, the builder also takes the mixing angles, phases, and mass-squared splitting of each sterile state, under the same names as in the wrappers:

.. code-block:: python

   st4 = dict(s14=0.1, s24=0.1, s34=0.1,
              d14=0.0, d24=0.0, D41=1.0)
   H = ham.hamiltonian_4nu_vacuum(
       E, **osc, **st4)
   H.shape                        # (4, 4)

   st5 = dict(st4, s15=0.1, s25=0.1, s35=0.1,
              d15=0.0, d35=0.0, D51=1.7)
   H = ham.hamiltonian_5nu_vacuum(
       E, **osc, **st5)
   H.shape                        # (5, 5)


The two calls return the :ref:`vacuum Hamiltonian <conv-hamiltonian>` at four and five flavors, i.e.,

.. math::
   :label: ex-equ-h-vacuum-45nu

   \begin{aligned}
   \mathbb{H}^{\rm vac}(E)
   &=
   \frac{1}{2E}\,
   \mathbb{R}\,
   \mathbb{M}^2\,
   \mathbb{R}^\dagger \;,
   \\
   \mathbb{M}^2_{3+1}
   &=
   {\rm diag}\!\left(0,\, \Delta m^2_{21},\, \Delta m^2_{31},\, \Delta m^2_{41}\right) ,
   \\
   \mathbb{M}^2_{3+2}
   &=
   {\rm diag}\!\left(0,\, \Delta m^2_{21},\, \Delta m^2_{31},\, \Delta m^2_{41},\, \Delta m^2_{51}\right) ,
   \end{aligned}

with :math:`\mathbb{R}` the :ref:`mixing matrix <conv-hamiltonian>` at :math:`n = 4` and 5, respectively.

A builder marked ``_energy_independent`` returns the same matrix without its :math:`1/E` factor. At one energy, the two builders are interchangeable. Across a scan in energy, the energy-independent one is faster, because the mixing rotation does not depend on the energy and is then computed only once. Over the 200 energies ``Es`` of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, rebuilding the Hamiltonian at each energy takes about 4 ms, and building it once and dividing it by each energy takes 0.2 ms:

.. code-block:: python

   # Interchangeable at one energy
   H1 = ham.hamiltonian_3nu_vacuum(E, **osc)
   H0 = ham.\
    hamiltonian_3nu_vacuum_energy_independent(
       **osc)
   H1 - H0/E                     # All zeros

   # Over 200 energies: 4.1 ms
   Hs = [ham.hamiltonian_3nu_vacuum(e, **osc)
         for e in Es]
   # ... and 0.16 ms
   Hs = [H0/e for e in Es]


*Matter.—*\ At three flavors, the matter term is

.. code-block:: python

   V = gd.VCC_EARTH_CRUST         # 1.145e-13 eV
   H = ham.hamiltonian_3nu_matter(V)
   H.shape                        # (3, 3)


The call returns the matter term of the :ref:`matter Hamiltonian <conv-hamiltonian>` at three flavors, i.e.,

.. math::
   :label: ex-equ-h-matter-3nu

   V_{\rm CC}\, \mathbb{P}
   =
   V_{\rm CC}\,
   {\rm diag}\!\left(1,\, 0,\, 0\right) \;,

the charged-current potential acting on the electron flavor alone. (``VCC_EARTH_CRUST`` is :math:`V_{\rm CC}` in matter of density 3 g cm\ :math:`^{-3}` with :math:`Y_e = 0.5`.) At two flavors, the call is the same, and the projector is :math:`{\rm diag}(1, 0)`.

At four and five flavors, the projector also has an entry for each sterile flavor, :math:`r/2`, where :math:`r` is the neutron-to-proton ratio of the medium (:doc:`/conventions`). The builder takes :math:`r` as a second parameter, 1 by default. In neutral matter, :math:`r = (1 - Y_e)/Y_e` (:ref:`ex-sec-constant-density`), so :math:`r` follows from the electron fraction; for instance, from that of the core of the Earth, which Magνs stores in ``magnus.earth``:

.. code-block:: python

   import magnus.earth as earth

   Ye = earth.Y_E_CORE_PREM       # 0.4656
   r = (1.0 - Ye)/Ye              # 1.15

   H = ham.hamiltonian_4nu_matter(V, r)
   H.shape                        # (4, 4)

   H = ham.hamiltonian_5nu_matter(V, r)
   H.shape                        # (5, 5)


The builder takes the potential and the ratio independently. The two calls return

.. math::
   :label: ex-equ-h-matter-45nu

   \begin{aligned}
   V_{\rm CC}\, \mathbb{P}_{3+1}
   &=
   V_{\rm CC}\,
   {\rm diag}\!\left(1,\, 0,\, 0,\, r/2\right) \;,
   \\
   V_{\rm CC}\, \mathbb{P}_{3+2}
   &=
   V_{\rm CC}\,
   {\rm diag}\!\left(1,\, 0,\, 0,\, r/2,\, r/2\right) \;.
   \end{aligned}

Four of the ten variants in Table :ref:`ex-tab-hamiltonians` end in ``_td`` (for *time-dependent*): the two vacuum forms, ``matter``, and ``nsi``. Each takes the position along the path, :math:`l`, as its first argument, just as the function of position passed to ``osc_prob`` does. The vacuum ones ignore it and return the same matrix at every position. ``matter_td`` and ``nsi_td`` take the potential as a function of position, ``f``, rather than as a number, and evaluate it at :math:`l`; e.g., ``matter_td(l, f)`` returns the same matrix as ``matter(f(l))``.

*Non-standard interactions.—*\ At three flavors, the NSI term is built as below. This block and the other NSI and LIV blocks show the signatures only: each argument stands for a number.

.. code-block:: python

   H = ham.hamiltonian_3nu_nsi(
       V, eps_ee, eps_em, eps_et,
       eps_mm, eps_mt, eps_tt)
   H.shape                        # (3, 3)


The call returns the charged-current potential times the matrix of NSI couplings, i.e.,

.. math::
   :label: ex-equ-h-nsi-3nu

   V_{\rm CC}
   \begin{pmatrix}
   \varepsilon_{ee} & \varepsilon_{e\mu} & \varepsilon_{e\tau} \\
   \varepsilon_{e\mu}^\ast & \varepsilon_{\mu\mu} & \varepsilon_{\mu\tau} \\
   \varepsilon_{e\tau}^\ast & \varepsilon_{\mu\tau}^\ast & \varepsilon_{\tau\tau}
   \end{pmatrix} \;.

The six couplings passed are the diagonal entries and those above it. The builder fills the entries below the diagonal with their complex conjugates, so the matrix is Hermitian even when the off-diagonal couplings are complex. At two flavors, the builder takes two couplings, ``eps_aa`` and ``eps_ab``, and returns

.. math::
   :label: ex-equ-h-nsi-2nu

   V_{\rm CC}
   \begin{pmatrix}
   \varepsilon_{aa} & \varepsilon_{ab} \\
   \varepsilon_{ab}^\ast & 0
   \end{pmatrix} \;,

where :math:`a` and :math:`b` are the two flavors, :math:`e` and :math:`\mu` by default (:ref:`ex-sec-two-flavors`). The lower diagonal entry is zero because only the difference between the two diagonal entries affects probabilities. Adding the same amount to both shifts every energy level of the Hamiltonian equally, which changes only an overall phase of the evolution, so :math:`\varepsilon_{aa}` stands for the difference :math:`\varepsilon_{aa} - \varepsilon_{bb}`.

At four and five flavors, the matrix of couplings also has rows and columns for the sterile flavors, and the builder takes the entries on and above its diagonal, ten and fifteen of them:

.. code-block:: python

   H = ham.hamiltonian_4nu_nsi(
       V, eps_ee, eps_em, eps_et, eps_es,
       eps_mm, eps_mt, eps_ms, eps_tt,
       eps_ts, eps_ss)
   H.shape                        # (4, 4)

   H = ham.hamiltonian_5nu_nsi(
       V, eps_ee, eps_em, eps_et, eps_es1,
       eps_es2, eps_mm, eps_mt, eps_ms1,
       eps_ms2, eps_tt, eps_ts1, eps_ts2,
       eps_s1s1, eps_s1s2, eps_s2s2)
   H.shape                        # (5, 5)


The two calls return the matrix whose entries are

.. math::
   :label: ex-equ-h-nsi-45nu

   \left[\mathbb{H}^{\rm NSI}\right]_{\alpha\beta}
   =
   V_{\rm CC}\, \varepsilon_{\alpha\beta} \;,
   \qquad
   \varepsilon_{\beta\alpha} = \varepsilon_{\alpha\beta}^\ast \;,

with :math:`\alpha, \beta \in \{e, \mu, \tau, s\}` at four flavors and :math:`\alpha, \beta \in \{e, \mu, \tau, s_1, s_2\}` at five; Eq. :eq:`ex-equ-h-nsi-3nu` is the same expression at three flavors.

*Lorentz-invariance violation.—*\ At three flavors, the Lorentz-invariance-violating (LIV) term is

.. code-block:: python

   H = ham.hamiltonian_3nu_liv(
       E, sxi12, sxi23, sxi13, dxiCP,
       b1, b2, b3, Lambda, n_liv)
   H.shape                        # (3, 3)


The call returns

.. math::
   :label: ex-equ-h-liv-3nu

   \mathbb{H}^{\rm LIV}(E)
   =
   \left(\frac{E}{\Lambda}\right)^{n_{\rm LIV}}
   \mathbb{V}_\xi\,
   {\rm diag}\!\left(b_1,\, b_2,\, b_3\right)
   \mathbb{V}_\xi^\dagger \;.

Here, :math:`b_1`, :math:`b_2`, and :math:`b_3` are the eigenvalues of the LIV operator, in eV. The unitary matrix :math:`\mathbb{V}_\xi` rotates them into the flavor basis; it is built from three angles, :math:`\xi_{12}`, :math:`\xi_{23}`, and :math:`\xi_{13}`, and one phase, :math:`\delta_\xi`, in the same way as the PMNS matrix (:ref:`conv-hamiltonian`). The prefactor sets how the term grows with energy: :math:`\Lambda` is the energy scale of the operator, in eV, and :math:`n_{\rm LIV}` is its dimension minus three, so that :math:`n_{\rm LIV} = 0` gives a term that does not depend on the energy.

The angles :math:`\xi_{ij}` and the phase :math:`\delta_\xi` are parameters of the LIV operator alone, independent of the PMNS angles and phase. They fix the basis in which the operator is diagonal, which need be neither the flavor basis nor the mass basis. For instance, with every :math:`\xi_{ij}` and :math:`\delta_\xi` set to zero, :math:`\mathbb{V}_\xi` is the identity and the operator is diagonal in the flavor basis; with them set equal to the PMNS values, it is diagonal in the mass basis.

At four and five flavors, the operator has one eigenvalue per flavor, and :math:`\mathbb{V}_\xi` is built as the :ref:`mixing matrix <conv-hamiltonian>` :math:`\mathbb{R}` at :math:`n = 4` and :math:`n = 5`, with angles :math:`\xi_{ij}` in place of :math:`\theta_{ij}` and phases :math:`\delta_{\xi, ij}` in place of :math:`\delta_{ij}`:

.. code-block:: python

   H = ham.hamiltonian_4nu_liv(
       E, sxi12, sxi23, sxi13, dxiCP,
       sxi14, dxi14, sxi24, dxi24, sxi34,
       b1, b2, b3, b4, Lambda, n_liv)
   H.shape                        # (4, 4)

   H = ham.hamiltonian_5nu_liv(
       E, sxi12, sxi23, sxi13, dxiCP,
       sxi14, dxi14, sxi15, dxi15, sxi24,
       dxi24, sxi25, sxi34, sxi35, dxi35,
       b1, b2, b3, b4, b5, Lambda, n_liv)
   H.shape                        # (5, 5)


The two calls return

.. math::
   :label: ex-equ-h-liv-45nu

   \begin{aligned}
   \mathbb{H}_{3+1}^{\rm LIV}(E)
   &=
   \left(\frac{E}{\Lambda}\right)^{n_{\rm LIV}}
   \mathbb{V}_\xi\,
   {\rm diag}\!\left(b_1,\, b_2,\, b_3,\, b_4\right)
   \mathbb{V}_\xi^\dagger \;,
   \\
   \mathbb{H}_{3+2}^{\rm LIV}(E)
   &=
   \left(\frac{E}{\Lambda}\right)^{n_{\rm LIV}}
   \mathbb{V}_\xi\,
   {\rm diag}\!\left(b_1,\, b_2,\, b_3,\, b_4,\, b_5\right)
   \mathbb{V}_\xi^\dagger \;.
   \end{aligned}

At two flavors, :math:`\mathbb{V}_\xi` is a single rotation, by the angle :math:`\xi`, and the operator has two eigenvalues.

The builders marked ``_energy_independent`` return these terms without the factor :math:`E^{n_{\rm LIV}}`, which you then multiply back in at each energy, as for the vacuum term.

*Pseudo-Dirac.—*\ A pseudo-Dirac neutrino is a mass state split into a nearly degenerate pair, one member active and the other sterile. The pseudo-Dirac builders take the mixing matrix and the mass-squared values of the active sector, and the mass states to split:

.. code-block:: python

   R = ham.pmns_mixing_matrix(
       osc['s12'], osc['s23'], osc['s13'],
       osc['dCP'])
   H = ham.hamiltonian_pseudo_dirac_vacuum(
       E, R, [0.0, osc['D21'], osc['D31']],
       {1: 1.0e-18})
   H.shape                        # (4, 4)


The last argument lists the mass states to split, each with the mass-squared splitting of its pair, :math:`\delta m^2_j`, in eV\ :math:`^2`. The mass states are counted from 0, so ``{1: 1.0e-18}`` splits :math:`\nu_2` by :math:`\delta m^2_2 = 10^{-18}` eV\ :math:`^2`. States not listed stay unsplit; with an empty mapping, the builder returns the three-flavor vacuum Hamiltonian, to round-off.

The call returns the :ref:`vacuum Hamiltonian <conv-hamiltonian>` for the enlarged set of states, i.e.,

.. math::
   :label: ex-equ-h-pseudo-dirac

   \mathbb{H}^{\rm vac}_{\rm PD}(E)
   =
   \frac{1}{2E}\,
   \mathbb{R}_{\rm PD}\,
   \mathbb{M}^2_{\rm PD}\,
   \mathbb{R}_{\rm PD}^\dagger \;.

Each split state :math:`\nu_j` is replaced by two mass states, :math:`(\nu_j \pm s_j)/\sqrt{2}`, where :math:`s_j` is its sterile partner, with masses squared :math:`m_j^2` and :math:`m_j^2 + \delta m^2_j`. For the example above, which splits :math:`\nu_2`, the rows of :math:`\mathbb{R}_{\rm PD}` are the flavors :math:`e`, :math:`\mu`, :math:`\tau`, and :math:`s_2`, and its columns are :math:`\nu_1`, the two members of the pair, and :math:`\nu_3`:

.. math::
   :label: ex-equ-mixing-pseudo-dirac

   \mathbb{R}_{\rm PD}
   =
   \begin{pmatrix}
   \mathbb{R}_{e1} & \mathbb{R}_{e2}/\sqrt{2} & \mathbb{R}_{e2}/\sqrt{2} & \mathbb{R}_{e3} \\
   \mathbb{R}_{\mu 1} & \mathbb{R}_{\mu 2}/\sqrt{2} & \mathbb{R}_{\mu 2}/\sqrt{2} & \mathbb{R}_{\mu 3} \\
   \mathbb{R}_{\tau 1} & \mathbb{R}_{\tau 2}/\sqrt{2} & \mathbb{R}_{\tau 2}/\sqrt{2} & \mathbb{R}_{\tau 3} \\
   0 & 1/\sqrt{2} & -1/\sqrt{2} & 0
   \end{pmatrix} ,

with :math:`\mathbb{R}` the PMNS matrix, and

.. math::
   :label: ex-equ-mass-pseudo-dirac

   \mathbb{M}^2_{\rm PD}
   =
   {\rm diag}\!\left(0,\, \Delta m^2_{21},\, \Delta m^2_{21} + \delta m^2_2,\, \Delta m^2_{31}\right) .

In general, there is one state per active state, plus one per split state. Both matrices are also available on their own, from ``pseudo_dirac_mixing_matrix`` and ``pseudo_dirac_mass_squared``.

The matter builder, ``hamiltonian_pseudo_dirac_matter``, takes the potential, the number of active states, the same mapping of split states, and the neutron-to-proton ratio, :math:`r`, 1 by default. It returns :math:`V_{\rm CC}` times the projector in the same enlarged flavor basis. For the example above, which splits :math:`\nu_2`,

.. math::
   :label: ex-equ-matter-pseudo-dirac

   V_{\rm CC}\, \mathbb{P}_{\rm PD}
   =
   V_{\rm CC}\,
   {\rm diag}\!\left(1,\, 0,\, 0,\, r/2\right) \;,

with rows :math:`e`, :math:`\mu`, :math:`\tau`, and :math:`s_2`. Since each partner :math:`s_j` is sterile, its entry is :math:`r/2`, as for the sterile states of Eq. :eq:`ex-equ-h-matter-45nu`. Split states other than :math:`\nu_2` add one such entry each.

.. _ex-sec-building-new-hamiltonian:

Building a new Hamiltonian from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two ways to build a Hamiltonian that Magνs does not ship. The first is to add shipped terms, as :ref:`ex-sec-wrapper-vs-direct` did with a vacuum, a matter, and a non-standard term. With a ``_td`` builder, a term can also vary along the path. For instance, a matter potential that grows linearly along the baseline is

.. code-block:: python

   H_vac = ham.hamiltonian_3nu_vacuum(E, **osc)
   VCC_func = lambda l: V*(1.0 + l/L)

   def H_func(l):
       H_mat = ham.hamiltonian_3nu_matter_td(
           l, VCC_func)
       return H_vac + H_mat

   P = oscprob.osc_prob(H_func, 0.0, L,
                        rtol=1e-8, atol=1e-8)
   # Pme = 0.021309


The second way is to write the matrix directly. ``osc_prob`` requires only a function that returns a Hermitian matrix, of any size, so an entirely new term reaches the solver in the same way as a shipped one. For instance, an off-diagonal term between :math:`\nu_e` and :math:`\nu_\tau` that grows along the path is

.. code-block:: python

   def H_own(l):
       H = H_vac.copy()
       H[0, 2] += 1.0e-13*(l/L)
       H[2, 0] += 1.0e-13*(l/L)
       return H

   P = oscprob.osc_prob(H_own, 0.0, L,
                        rtol=1e-8, atol=1e-8)
   # Pme = 0.004865


The Hermiticity of the matrix is yours to guarantee, since it is what keeps the evolution unitary (:doc:`/comparison`). Magνs neither checks nor enforces it, and a matrix that is not Hermitian returns probabilities that do not add up to one. Every entry of the matrix is in eV, and every position in eV\ :math:`^{-1}` (:doc:`/conventions`). :ref:`ex-sec-lri-sun` shows a complete example.

A custom Hamiltonian may have more than five flavors. The wrappers of Table :ref:`ex-tab-wrappers`, the builders of Table :ref:`ex-tab-hamiltonians`, and the compiled exponential kernels of :doc:`/methodology` stop at five. Above five, you supply the matrix, either to ``osc_prob`` or, as a vacuum term, to a scenario function through ``h_vac_energy_indep`` (:ref:`ex-sec-wrapper-vs-direct`), and the exponential of each slab goes through a general Hermitian eigendecomposition, ``numpy.linalg.eigh``. Each slab then costs more: on one exponential density profile, at the default tolerance, a probability takes about twice as long at six flavors as at five, and three and a half times as long at eight. The probabilities remain unitary to within :math:`10^{-13}`.
