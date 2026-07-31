Code Architecture
===================

This page documents how Magνs's source code (under ``src/magnus/``) is
organized: the module layout, the three-layer call structure of the
oscillation-probability API, the contract between the layers, and a
worked walkthrough of how to add a new physics scenario yourself. See
:doc:`methodology` for the numerical machinery (the Magnus expansion
itself); this page is about the *code*, not the *math*.

Module layout
---------------

Magνs is split into eight modules under ``src/magnus/``, each with a
single, non-overlapping responsibility, all explicitly listed in
``magnus/__init__.py``'s ``submodules``/``__all__``:

.. mermaid::

   flowchart TD
       subgraph pure["Pure, self-contained modules (no internal dependencies)"]
           magnuscore["magnus.magnus<br/>Magnus expansion, GL integrators,<br/>matrix exponential"]
           ham["magnus.hamiltonians<br/>hamiltonians{2,3,4,5}nu.py<br/>mixing matrices, vacuum/matter/NSI/LIV H"]
       end
       gd["magnus.globaldefs<br/>physical constants, unit conversions,<br/>NuFit parameter sets"]
       earth["magnus.earth<br/>PREM density profile,<br/>chord/zenith-angle geometry"]
       matter["magnus.matter<br/>density profiles, electron number density,<br/>V_CC construction"]
       adiabatic["magnus.adiabatic<br/>adiabatic transport + Magnus-patch<br/>hybrid_propagator"]
       osc["magnus.oscprob<br/>the public API + oscprobstd.py<br/>(closed-form validation)"]

       gd --> earth
       gd --> matter
       magnuscore --> osc
       magnuscore --> adiabatic
       ham --> osc
       earth --> osc
       matter --> osc
       gd --> osc
       adiabatic --> osc

Two consequences of this layout are worth calling out because they are
easy to break by accident when adding code:

* ``magnus.magnus`` (the Magnus-expansion core) and ``magnus.hamiltonians``
  (the physics formulas) know nothing about each other, about
  ``oscprob``, or about ``earth``/``matter``. The core is pure numerical
  linear algebra on an arbitrary matrix function :math:`A(t)`; the
  Hamiltonians are pure algebra on mixing angles and potentials. Neither
  needs to change when the other does, and both are independently unit
  tested (``tests/test_magnus_expansion.py`` /
  ``tests/test_hamiltonians.py``) without touching ``oscprob`` at all.
  ``magnus.adiabatic`` follows the same rule: it depends only on
  ``magnus.magnus`` (for the local Magnus patch), never on ``oscprob``, so
  it is independently unit tested (``tests/test_adiabatic.py``) and usable
  directly on any Hamiltonian function, entirely outside the
  oscillation-probability API. See :doc:`adiabatic_strategy` for its
  numerical method.
* ``magnus.oscprob`` is the only module that imports everything else. It
  is where physics scenarios (vacuum/matter/NSI/LIV), environments
  (constant density/exponential density/Earth/Sun), and the Magnus core
  are wired together. This is deliberate: it keeps the wiring in one
  place instead of scattering it across the physics and numerical
  modules.

``earth.py``, ``globaldefs.py``, ``magnus.py``, ``adiabatic.py``,
``matter.py``, ``oscprob.py``, and ``oscprobstd.py`` are flat sibling
files directly under ``src/magnus/`` -- there is no subpackage directory
wrapping any of them. Only ``magnus.hamiltonians`` is a genuine
subpackage, since it holds four distinct, flavor-count-specific modules
(``hamiltonians2nu.py`` through ``hamiltonians5nu.py``); its
``__init__.py`` explicitly imports and re-exports each one's public
names (no ``from .module import *``). ``magnus/__init__.py`` itself
explicitly imports all eight top-level modules (again, no wildcard
imports) so that ``import magnus`` alone makes ``magnus.earth``,
``magnus.oscprob``, etc. immediately accessible. ``magnus.oscprob``
additionally imports and re-exports ``oscprobstd.py``'s five names (the
closed-form validation counterpart to the wrapper API), so both
``magnus.oscprob.osc_prob_3nu_vacuum_std`` and
``magnus.oscprobstd.osc_prob_3nu_vacuum_std`` work.

The three-layer structure of ``magnus.oscprob``
----------------------------------------------------

``magnus.oscprob`` is the largest module (~10,000 lines) because it exposes a
dedicated, explicitly-named function for every combination of
(flavor count) :math:`\times` (environment) :math:`\times` (BSM
scenario) — roughly 60 combinations. To keep that size from turning into
60 independent copies of the same logic (which is exactly what caused
several of the bugs this package's test suite now guards against — see
:ref:`layer-contract` below), every one of those 60 functions is a thin
call into a much smaller set of shared functions. There are three layers:

.. mermaid::

   flowchart TD
       subgraph L1["Layer 3 -- Wrappers (~60 functions)"]
           w1["osc_prob_3nu_earth(energy, costhz, s12=None, s23=None, ..., nubar=False, ...)"]
           w2["osc_prob_2nu_matter_nsi_constant_density(energy, L, rho, eps_aa, eps_ab, ...)"]
           w3["osc_prob_5nu_vacuum_liv(energy, L, s12=None, ..., b1=None, ..., Lambda=None, ...)"]
       end
       subgraph L2["Layer 2 -- Scenario (4 functions, generic in num_flavors)"]
           m1["osc_prob_vacuum(num_flavors, energy, L, osc_params, ...)"]
           m2["osc_prob_matter_std_potential(num_flavors, rho_func, energy, L, osc_params, ...)"]
           m3["osc_prob_matter_nsi(num_flavors, rho_func, energy, L, osc_params, nsi_params, ...)"]
           m4["osc_prob_liv(num_flavors, energy, L, osc_params, liv_params, ...)"]
       end
       subgraph L3["Layer 1 -- Primordial"]
           p1["osc_prob_energy_baseline(H_func, energy, L, ...)<br/>loops/parallelizes over (energy, L) points,<br/>warm-starts refinement across points"]
           p2["osc_prob(H_func, t_ini, t_fin, ...)<br/>owns adaptive slab refinement;<br/>calls the Magnus core once per point"]
       end
       core["magnus.magnus.compute_evolution_operator_multiple_slabs<br/>-> magnus_expansion_multislab"]

       w1 --> m2
       w2 --> m3
       w3 --> m4
       m1 --> p1
       m2 --> p1
       m3 --> p1
       m4 --> p1
       p1 --> p2 --> core
       direct["Your own H_func(l)"] -. bypasses layers 2 and 3 entirely .-> p2

**Layer 1 -- primordial.** ``osc_prob`` is the only function that calls
into the Magnus core. It owns the adaptive-refinement loop (grow
``n_slabs``/``n_tpts_per_slab``/``magnus_exp_order`` until ``rtol``/
``atol`` is met or a cap is hit), input validation, logging, and the
`~50`-line docstring documenting all of the refinement/logging keyword
arguments (see it directly in
:func:`~magnus.oscprob.osc_prob`). It is also a first-class
public entry point: pass it *any* callable ``H_func(l)`` (or
``H_func(enu, l)``, or a constant matrix) and it works with no wrapper at
all -- this is the escape hatch for Hamiltonians the package does not
already know about. ``osc_prob_energy_baseline`` sits just above it:
given arrays of ``energy`` and ``L``, it builds the right
energy-dependent closure over ``H_func``, decides whether to parallelize
over points (``joblib.Parallel``) or hand a single call straight to
``osc_prob``, and carries the *warm start* logic that seeds each point's
refinement from the previous point's converged (``n_slabs``,
``n_tpts_per_slab``).

**Layer 2 -- scenario.** ``osc_prob_vacuum``, ``osc_prob_matter_std_potential``,
``osc_prob_matter_nsi``, and ``osc_prob_liv`` are each generic in
``num_flavors`` (2, 3, 4, or 5): they unpack the relevant parameter dict
(via ``unpack_oscillation_params_from_dict``, ``unpack_nsi_params_from_dict``,
``unpack_liv_params_from_dict``), dispatch to the matching function in
``magnus.hamiltonians`` (e.g. ``hamiltonian_3nu_matter`` for
``num_flavors=3``), build the position-dependent matter potential where
relevant (via ``magnus.matter.vcc_func_from_rho_func``), and call
``osc_prob_energy_baseline`` with the resulting ``H_func``. This is where
"what physics scenario is this" is decided; it is *not* where "how many
flavors" or "which environment" is decided -- those come from the caller
(layer 3) and from ``rho_func``, respectively.

**Layer 3 -- wrappers.** Every ``osc_prob_{2,3,4,5}nu_{scenario}`` function
(e.g. ``osc_prob_3nu_matter_constant_density``,
``osc_prob_4nu_earth_nsi``, ``osc_prob_2nu_vacuum_liv``) exists purely so
that users get explicit, autocomplete-and-docs-friendly parameter names
(``s12``, ``eps_em``, ``rho``, ...) instead of having to build
``osc_params``/``nsi_params``/``liv_params`` dictionaries by hand. Its
entire job is: validate/repackage its named parameters into the right
dict(s), and forward everything else to the matching layer-2 function.
``osc_prob_earth``/``osc_prob_sun`` are a deliberate, bounded exception:
because they need to build a PREM-based (or solar-density-based)
``VCC_func`` and choose between the 2/3/4/5nu Hamiltonians themselves,
they sit one level below the per-flavor ``osc_prob_{N}nu_earth``/
``osc_prob_{N}nu_sun`` wrappers and forward a curated subset of
parameters *positionally* into ``_osc_prob_with_potential``, rather than
by name through ``**kwargs`` like every other wrapper.

.. _layer-contract:

The layer contract: what a wrapper must **not** do
------------------------------------------------------

Every wrapper function ends in ``**kwargs`` and must **not** redeclare
any of the refinement/logging keyword arguments that layers 1-2 own:

.. code-block:: text

   magnus_exp_order, n_jobs, integration_method, rtol, atol,
   growth_factor_n_slabs, growth_factor_n_tpts_per_slab, max_num_loops,
   min_n_slabs, max_n_slabs, min_n_tpts_per_slab, max_n_tpts_per_slab,
   iterate_over_magnus_exp_order, min_magnus_exp_order, max_magnus_exp_order,
   new_recursion_limit

This is not a style preference; it is a correctness requirement, and the
history of this package shows what happens when it is violated. Before
the refactor that introduced this contract (internally referred to as
"G1"), every wrapper declared its own copy of these ~15 parameters with
its own defaults. That duplication is exactly what let several bugs hide
for a long time: a wrapper with a silently different default tolerance
than its siblings, a wrapper missing ``nubar`` entirely, a wrapper with
an inconsistent validation bound. Fixing a default meant remembering to
fix it in ~60 places; inevitably, some were missed.

Two permanent tests in ``tests/test_oscprob.py`` enforce this contract in
CI, and will fail if it is ever violated again:

* ``test_no_wrapper_redeclares_standard_refinement_kwargs`` — inspects
  every ``osc_prob_{2,3,4,5}nu_*`` function's signature via
  :func:`inspect.signature` and fails if any of the 15 names above appear
  in it.
* ``test_nubar_present_across_all_flavor_counts_in_matter_families`` —
  fails if a matter/NSI/LIV wrapper family exposes ``nubar`` for some
  flavor counts but not others.

If you are adding a wrapper and find yourself typing
``rtol: Optional[float] = 1.e-3`` in its signature, that is a signal you
are working at the wrong layer: forward it through ``**kwargs`` instead.

Data flow: how the Hamiltonian and potential are built
-----------------------------------------------------------

The matter potential and the Hamiltonian are built once per call (not
once per slab), then passed down as a single callable:

.. mermaid::

   sequenceDiagram
       participant User
       participant W as osc_prob_3nu_matter_nsi_exp_density
       participant M as osc_prob_matter_nsi
       participant Matter as matter.vcc_func_from_rho_func
       participant Ham as hamiltonians.hamiltonian_3nu_nsi_td
       participant B as osc_prob_energy_baseline
       participant P as osc_prob
       participant K as magnus core

       User->>W: energy, L, rho_central, l_scale, eps_ee, ...
       W->>M: num_flavors=3, rho_func=exp profile, osc_params, nsi_params
       M->>Matter: build VCC_func(l) from the density profile
       M->>Ham: wrap VCC_func into H(l) = H_vac + H_matter(l) + H_NSI(l)
       M->>B: H_func=H(l), energy array, L array
       loop for each (energy, L) point (parallel if n_jobs != 1)
           B->>P: osc_prob(H_func, 0, L_i, ...)
           P->>K: adaptively refine slabs until rtol/atol met
           K-->>P: unitary evolution operator U
           P-->>B: probability matrix P = |U|^2
       end
       B-->>User: array of probability matrices

Every intermediate object here is a plain Python callable
(``VCC_func: l -> float``, ``H_func: l -> ndarray``); nothing is
precomputed on a grid before reaching ``osc_prob``, which is what lets
:func:`magnus.magnus.probe_eval_mode` decide, once, whether ``H_func``
can be evaluated on a vectorized array of positions (silent
vectorization -- see :doc:`methodology`) or must be called one position
at a time.

How to add your own wrapper
------------------------------

Suppose you want to add support for a new environment, e.g. a
user-supplied radial density profile for 3-flavor NSI oscillations,
``osc_prob_3nu_matter_nsi_custom_density``. The existing
``osc_prob_3nu_matter_nsi_exp_density`` (in ``magnus.oscprob``) is the
closest sibling to copy from. The recipe:

#. **Pick the right layer-2 function.** You are adding an environment
   (a new ``rho_func``), not a new physics scenario, so you call the
   existing ``osc_prob_matter_nsi`` — you do **not** need to touch
   ``magnus.hamiltonians`` or the Magnus core at all.

#. **Name only the parameters specific to your scenario.** Your
   function's signature should have: the standard positional physics
   inputs (``energy``, ``L``), whatever parametrizes *your* density
   profile (e.g. a ``density_func: Callable`` the user supplies
   directly), the standard oscillation parameters for 3 flavors
   (``s12, s23, s13, dCP, D21, D31``, all ``Optional[float] = None``),
   the standard NSI parameters (``eps_ee, eps_em, ...``), the standard
   trailing parameters every wrapper has
   (``ratio_number_neutrons_to_protons``, ``electron_fraction``,
   ``nubar``, ``nu_i``, ``nu_f``, ``validate_input``, ``save_log``,
   ``filename_log``, ``file_log``, ``close_file_log_upon_exit``,
   ``verbose``), and end with ``**kwargs``.

#. **Do not name any of the 15 refinement/logging kwargs listed in**
   :ref:`layer-contract` **above.** They flow through ``**kwargs``
   automatically. This is what the two permanent guard tests check.

#. **Write the body as a single call down**, packaging your named
   parameters into the ``osc_params``/``nsi_params`` dicts that
   ``osc_prob_matter_nsi`` expects:

   .. code-block:: python

       def osc_prob_3nu_matter_nsi_custom_density(
           energy, L, density_func,
           s12=None, s23=None, s13=None, dCP=None, D21=None, D31=None,
           eps_ee=0.0, eps_em=0.0j, eps_et=0.0j,
           eps_mm=0.0, eps_mt=0.0j, eps_tt=0.0,
           ratio_number_neutrons_to_protons=1.0, electron_fraction=0.5,
           nubar=False, nu_i=None, nu_f=None,
           validate_input=True, save_log=False, filename_log='./out.log',
           file_log=None, close_file_log_upon_exit=True, verbose=0,
           **kwargs
       ):
           r"""Compute the 3nu NSI oscillation probability for a
           user-supplied radial matter density profile.

           .. versionadded:: 0.10.0
           """
           return osc_prob_matter_nsi(
               num_flavors=3,
               rho_func=density_func,
               energy=energy, L=L,
               osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP,
                           'D21': D21, 'D31': D31},
               nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et,
                           'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_tt': eps_tt},
               ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
               electron_fraction=electron_fraction,
               nubar=nubar, nu_i=nu_i, nu_f=nu_f,
               validate_input=validate_input, save_log=save_log,
               filename_log=filename_log, file_log=file_log,
               close_file_log_upon_exit=close_file_log_upon_exit,
               verbose=verbose,
               **kwargs
           )

#. **Add it to the family-consistency tests.** ``test_oscprob.py``
   parametrizes several checks over "every osc_prob wrapper family" by
   name pattern; add your new function's family prefix alongside its
   3 siblings (2nu/4nu/5nu, if you are adding all four) so the same
   unitarity/``nubar``-sensitivity/API-shape checks cover it
   automatically instead of needing bespoke tests.

#. **Run the two guard tests** described in :ref:`layer-contract` before
   opening a pull request:

   .. code-block:: bash

      pytest tests/test_oscprob.py -k "redeclares_standard or nubar_present" -v

If your new function needs genuinely new *physics* (not just a new
environment) -- e.g. a Hamiltonian term that does not fit
vacuum/matter/NSI/LIV -- then the right layer to extend is layer 2: add
a new ``osc_prob_<scenario>`` function generic in ``num_flavors``,
following the pattern of ``osc_prob_liv``, and a matching
``hamiltonian_<n>nu_<scenario>`` in ``magnus.hamiltonians`` for each
flavor count you support. Only then add layer-3 wrappers on top of it.

Where things live: a quick lookup
-------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - If you need to change...
     - ...look in
   * - A default tolerance, the refinement/adaptive-slab-growth logic,
       or anything every scenario shares
     - ``osc_prob`` / ``osc_prob_energy_baseline``
       (``src/magnus/oscprob.py``)
   * - How a physics scenario's Hamiltonian is assembled from mixing
       angles/NSI epsilons/LIV coefficients
     - ``osc_prob_vacuum`` / ``osc_prob_matter_std_potential`` /
       ``osc_prob_matter_nsi`` / ``osc_prob_liv``
   * - The actual matrix form of a Hamiltonian (mixing matrix, vacuum
       term, matter term, NSI term, LIV term)
     - ``magnus.hamiltonians.hamiltonians{2,3,4,5}nu``
   * - A named parameter exposed to end users for one (flavor count,
       environment, scenario) combination
     - the matching ``osc_prob_{N}nu_{scenario}`` wrapper
   * - The Magnus term recursion, the Gauss-Legendre integrators, or the
       matrix exponential itself
     - ``magnus.magnus`` (:doc:`methodology`)
   * - The PREM density profile or Earth chord/zenith geometry
     - ``magnus.earth``
   * - A generic density profile, electron number density, or the
       :math:`V_{CC}` potential construction
     - ``magnus.matter``
   * - A physical constant, unit conversion, or a predefined oscillation
       parameter set (e.g. NuFit 6.0)
     - ``magnus.globaldefs``
