.. _ex-sec-params:

Choosing the oscillation parameters
-----------------------------------

.. contents:: On this page
   :local:
   :depth: 1

Table :ref:`ex-tab-parameters` lists every parameter that the shipped wrappers take, at two to five flavors, together with its default. The subsections below show how they are set.

.. _ex-tab-parameters:

.. table:: The parameters of the wrapper functions

   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   |                                             | **Two flavors**        | **Three flavors**                   | **Four flavors (3+1)**                                       | **Five flavors (3+2)**                                              |
   +=============================================+========================+=====================================+==============================================================+=====================================================================+
   | *Standard oscillations, in every wrapper*                                                                                                                                                                                                       |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Mixing angles                               | ``sth``                | ``s12``, ``s23``, ``s13``           | as 3\ :math:`\nu`, :math:`+` ``s14``, ``s24``, ``s34``       | as 3\ :math:`\nu`, :math:`+` ``s14``, ``s15``, ``s24``,             |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   |                                             |                        |                                     |                                                              | ``s25``, ``s34``, ``s35``                                           |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | CP phases                                   | —                      | ``dCP``                             | as 3\ :math:`\nu`, :math:`+` ``d14``, ``d24``                | as 3\ :math:`\nu`, :math:`+` ``d14``, ``d15``, ``d24``, ``d35``     |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Mass splittings [eV\ :math:`^2`]            | ``Dm2``                | ``D21``, ``D31``                    | as 3\ :math:`\nu`, :math:`+` ``D41``                         | as 3\ :math:`\nu`, :math:`+` ``D41``, ``D51``                       |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Default                                     | required               | named set                           | 3\ :math:`\nu`: named set; sterile: 0                        | 3\ :math:`\nu`: named set; sterile: 0                               |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | *Non-standard interactions, in the ``_nsi`` wrappers; default 0*                                                                                                                                                                                |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Couplings :math:`\varepsilon_{\alpha\beta}` | ``eps_aa``, ``eps_ab`` | ``eps_ee``, ``eps_em``, ``eps_et``, | as 3\ :math:`\nu`, :math:`+` ``eps_es``, ``eps_ms``,         | as 3\ :math:`\nu`, :math:`+` ``eps_es1``, ``eps_es2``, ``eps_ms1``, |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   |                                             |                        | ``eps_mm``, ``eps_mt``, ``eps_tt``  | ``eps_ts``, ``eps_ss``                                       | ``eps_ms2``, ``eps_ts1``, ``eps_ts2``, ``eps_s1s1``,                |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   |                                             |                        |                                     |                                                              | ``eps_s1s2``, ``eps_s2s2``                                          |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | *Lorentz-invariance violation, in the ``_liv`` wrappers; default 0, except :math:`\Lambda = 1`*                                                                                                                                                 |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Mixing angles :math:`\xi_{ij}`              | ``sxi``                | ``sxi12``, ``sxi23``, ``sxi13``     | as 3\ :math:`\nu`, :math:`+` ``sxi14``, ``sxi24``, ``sxi34`` | as 3\ :math:`\nu`, :math:`+` ``sxi14``, ``sxi15``, ``sxi24``,       |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   |                                             |                        |                                     |                                                              | ``sxi25``, ``sxi34``, ``sxi35``                                     |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | CP phases                                   | —                      | ``dxiCP``                           | ``dxiCP``, ``dxi14``, ``dxi24``                              | ``dxiCP``, ``dxi14``, ``dxi15``,                                    |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   |                                             |                        |                                     |                                                              | ``dxi24``, ``dxi35``                                                |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Eigenvalues [eV]                            | ``b1``, ``b2``         | ``b1``–``b3``                       | ``b1``–``b4``                                                | ``b1``–``b5``                                                       |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+
   | Scale, power                                | ``Lambda``, ``n_liv``  | ``Lambda``, ``n_liv``               | ``Lambda``, ``n_liv``                                        | ``Lambda``, ``n_liv``                                               |
   +---------------------------------------------+------------------------+-------------------------------------+--------------------------------------------------------------+---------------------------------------------------------------------+

Every parameter the wrappers of Table :ref:`ex-tab-wrappers` take, by flavor count, beyond the energy, the baseline, and the description of the medium. The angles are read in the convention set by ``angles`` (:ref:`ex-sec-angle-conventions`); the phases in radians, or in degrees under ``angles='deg'``. The standard three-flavor parameters left unset are read from the named set ``OSC_PARAMS_DEFAULT``, or from another set named through ``default_osc_params_set_name``; the active-sterile ones default to zero, so that a four- or five-flavor call without them returns the three-flavor probabilities. At two flavors, ``sth`` and ``Dm2`` are required. The non-standard parameters default to zero, so a ``_nsi`` or ``_liv`` wrapper called without them returns the standard result; with :math:`n_{\rm LIV} = 0`, a term switched on through the eigenvalues is energy-independent.

Three flavors
~~~~~~~~~~~~~

Every ``osc_prob_3nu_*`` wrapper leaves its six standard parameters unset by default. An unset parameter is read from a named set, ``OSC_PARAMS_DEFAULT``, which is the NuFIT 6.1 best fit with Super-Kamiokande atmospheric data in normal ordering  :cite:p:`Esteban:2024eli`. Naming a different set changes all six at once.

.. code-block:: python

   P = oscprob.osc_prob_3nu_vacuum(E, L,
    default_osc_params_set_name=
     'OSC_PARAMS_NU_FIT_5_2_SK_IO')
   # Pme = 0.05166, against 0.03127 by default


Fifty-two NuFIT global-fit sets of best-fit three-flavor mixing parameters are predefined, together with ``OSC_PARAMS_DEFAULT``, as the keys of ``gd.OSC_PARAMS_PREDEFINED``; a set is passed by its name, as a string. The names follow the releases. From NuFIT 4.0 onward, a release splits its fits by whether Super-Kamiokande atmospheric data is included; both are included, e.g., ``OSC_PARAMS_NU_FIT_5_2_SK_IO`` and ``OSC_PARAMS_NU_FIT_5_2_NOSK_IO``. Earlier releases carry no such split, so their names drop that infix, as in ``OSC_PARAMS_NU_FIT_3_0_NO``. Every set can also be loaded by release, ordering, and category through ``load_nufit_params``, e.g.,

.. code-block:: python

   osc = gd.load_nufit_params('NuFIT 5.2',
    ordering='IO', category='without_SK')
   P = oscprob.osc_prob_3nu_vacuum(E, L, **osc)


Earlier releases divide their fits in other ways, by reactor-flux treatment before 2.0 and by MINOS event selection in 2.1. The named sets hold one category of each release; passing ``category`` to the loader reaches the others, e.g.,

.. code-block:: python

   # Additional set categories
   osc = gd.load_nufit_params('NuFIT 2.1',
    category='LID')   # MINOS event selection
   P = oscprob.osc_prob_3nu_vacuum(E, L,
    **osc)
   # Pme = 0.05579; 'LEM' gives 0.05558

   osc = gd.load_nufit_params('NuFIT 1.3',
    category='huber_fluxes_no_rsbl')
   P = oscprob.osc_prob_3nu_vacuum(E, L,
    **osc)
   # Pme = 0.04998; 'free_fluxes_rsbl'
   # gives 0.04552


:doc:`/conventions` describes the predefined sets. The full list of sets, and that of the releases and categories the loader reads, can be printed from Python, which keeps them current with the installed version:

.. code-block:: python

   # Every predefined set, with its description
   sets = gd.OSC_PARAMS_PREDEFINED
   for name in sorted(sets):
       desc = sets[name]['description']
       print(name, '|', desc)
   # 53 lines, e.g.,
   # OSC_PARAMS_NU_FIT_1_0_IO | NuFIT 1.0, IO

   # The releases the loader reads, and the
   # categories of each
   fits = gd.NUFIT_GLOBAL_FITS
   for v in fits:
       print(v, list(fits[v]['categories']))
   # 18 lines, e.g., NuFIT 2.1 ['LEM', 'LID']


Individual parameter values are passed by name. Anything left unset still comes from the default set, so one parameter can be moved without restating the other five, e.g.,

.. code-block:: python

   P = oscprob.osc_prob_3nu_vacuum(E, L,
    dCP=0.0)
   # Pme = 0.05098; the other five stay at
   # their defaults


The two mechanisms combine. A named set chooses the base, and a parameter passed alongside it replaces that one entry:

.. code-block:: python

   # A named set, with one value moved
   P = oscprob.osc_prob_3nu_vacuum(E, L,
    default_osc_params_set_name=
     'OSC_PARAMS_NU_FIT_6_1_SK_IO', dCP=0.0)
   # The other five stay at their NuFIT 6.1 IO
   # values.  Pme = 0.01800, against
   # 0.05241 with the set's own dCP


.. _ex-sec-angle-conventions:

Specifying angle conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mixing angles can be stated four ways, which the keyword ``angles`` selects: ``'sin'`` by default, ``'sin2'`` for the sines squared that global fits report, ``'rad'`` for the angles themselves, or ``'deg'``. With ``'deg'``, the CP phases are read as degrees as well. The keyword has to be given the same value in two places, on ``load_nufit_params`` and on the probability call, since a set loaded under one convention and read under another usually returns a plausible, wrong answer:

.. code-block:: python

   # The convention travels with the numbers
   for conv in ('sin', 'sin2', 'rad', 'deg'):
       osc = gd.load_nufit_params('NuFIT 6.1',
        angles=conv)
       P = oscprob.osc_prob_3nu_vacuum(E, L,
        angles=conv, **osc)
   # s12 reads 0.5557, 0.3088, 0.5892, 33.759
   # and Pme = 0.03127 in all four cases


.. _ex-sec-two-flavors:

Two flavors
~~~~~~~~~~~

Two flavors take a different route. A two-flavor system has one mixing angle and one splitting, which no global fit reports as such, so ``osc_prob_2nu_*`` has no parameter set to fall back on. It does not accept ``default_osc_params_set_name``. Both ``sth`` and ``Dm2`` are required arguments; omitting either raises an error. The ``angles`` keyword applies as it does at three flavors, so ``sth`` is read as a sine unless another convention is named:

.. code-block:: python

   # Two flavors: both values are required
   P = oscprob.osc_prob_2nu_vacuum(E, L,
    sth=0.5557, Dm2=7.537e-5)
   # P[0][1] = 0.013089

   # The same point, angles in degrees
   P = oscprob.osc_prob_2nu_vacuum(E, L,
    sth=33.759, Dm2=7.537e-5, angles='deg')


.. _ex-sec-four-flavors:

Four flavors
~~~~~~~~~~~~

A fourth, sterile flavor adds three mixing angles, two CP phases, and one mass splitting to the six three-flavor parameters, in the parametrization of Ref. :cite:p:`Kopp:2011qd`. The six three-flavor parameters behave as at three flavors: left unset, they are read from the named set. The six new ones default to zero, so a four-flavor call without them returns the three-flavor probabilities in its upper :math:`3 \times 3` block. All twelve can be given at once:

.. code-block:: python

   # Four flavors (3+1): twelve parameters
   P = oscprob.osc_prob_4nu_vacuum(E, L,
    s12=0.3088, s23=0.4700, s13=0.02248,
    dCP=3.700, D21=7.537e-5, D31=2.511e-3,
    s14=0.10, s24=0.10, s34=0.0,
    d14=0.0, d24=0.0, D41=1.0,
    angles='sin2')
   # Pme = 0.02801, Pmm = 0.39958


Here, the angles are given as their sines squared, as global fits report them; the first six values are the NuFIT 6.1 best fit of the default set, rounded. The call returns a :math:`4 \times 4` matrix, indexed by ``NUE``, ``NUMU``, ``NUTAU``, and ``NUS`` for the sterile flavor; at five flavors, the two sterile flavors are ``NUS1`` and ``NUS2``.

.. _ex-sec-five-flavors:

Five flavors
~~~~~~~~~~~~

A fifth flavor adds three further mixing angles, two further CP phases, and a second splitting, for 18 parameters in all:

.. code-block:: python

   # Five flavors (3+2): 18 parameters
   P = oscprob.osc_prob_5nu_vacuum(E, L,
    s12=0.3088, s23=0.4700, s13=0.02248,
    dCP=3.700, D21=7.537e-5, D31=2.511e-3,
    s14=0.10, s24=0.10, s34=0.0,
    s15=0.06, s25=0.06, s35=0.0,
    d14=0.0, d15=0.0, d24=0.0, d35=0.0,
    D41=1.0, D51=1.7, angles='sin2')
   # Pme = 0.03984, Pmm = 0.40684


As at four flavors, every active-sterile parameter defaults to zero. The mixing matrix has nine angles and five CP phases, ``dCP``, ``d14``, ``d15``, ``d24``, and ``d35``, where a general :math:`5 \times 5` unitary matrix has ten angles and six phases. The missing angle and phase belong to a rotation between the two sterile flavors. Since the two sterile states are indistinguishable, that rotation changes no probability among the active flavors, so it is left out.
