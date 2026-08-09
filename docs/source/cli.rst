Command-Line Calculator
=========================

.. contents::
   :local:
   :depth: 2


In addition to the Python API (:doc:`quickstart`), Magνs installs a
``magnus`` command (equivalently, ``python -m magnus``) for computing a
single oscillation probability directly from the shell, with no Python
required. It wraps the same ``osc_prob_{2,3,4,5}nu_*`` functions used by
the Python API: ``magnus prob`` picks the right one from
``--flavors``/``--environment``/``--scenario`` and calls it with the flags
you give (see :doc:`architecture` for how those functions themselves are
organized).

.. note::
   The CLI computes **one probability at a time** (a single energy and
   baseline). For scans, plots, or fitting, use the Python API -- it is
   the same underlying code, just called in a loop or with array
   arguments (which Magνs evaluates in a single batched, vectorized pass;
   see :doc:`methodology`).

Installation
--------------

.. code-block:: bash

   git clone https://github.com/mbustama/Magnus.git
   cd Magnus
   pip install -e .

This installs the ``magnus`` console script. If you would rather not
install the package, ``python -m magnus`` works identically from the
repository root once ``src/`` is on ``PYTHONPATH``.

``magnus --version`` (or ``-V``) prints the installed version and exits;
it reports the same number as ``magnus.__version__``, which is read
from the ``version`` field of ``pyproject.toml``.

Usage pattern
---------------

.. code-block:: text

   magnus prob --flavors {2,3,4,5} --environment {vacuum,matter,earth,sun}
               --scenario {std,nsi,liv} [environment- and scenario-specific flags]
               --energy ENERGY [--energy-unit UNIT] [--baseline BASELINE] ...

``--environment`` selects the propagation medium; ``--scenario`` selects
the physics on top of it. Not every combination exists: ``--scenario nsi``
requires a matter potential to modify, so it is not available with
``--environment vacuum`` (the CLI rejects this combination with a clear
error rather than silently ignoring the epsilon flags). The full dispatch
table:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 40

   * - ``--environment``
     - ``--scenario``
     - ``--density-profile``
     - Function called (per ``--flavors``)
   * - vacuum
     - std
     - --
     - ``osc_prob_{N}nu_vacuum``
   * - vacuum
     - liv
     - --
     - ``osc_prob_{N}nu_vacuum_liv``
   * - matter
     - std
     - constant / exp
     - ``osc_prob_{N}nu_matter_{constant,exp}_density``
   * - matter
     - nsi
     - constant / exp
     - ``osc_prob_{N}nu_matter_nsi_{constant,exp}_density``
   * - matter
     - liv
     - constant / exp
     - ``osc_prob_{N}nu_matter_liv_{constant,exp}_density``
   * - earth
     - std / nsi / liv
     - --
     - ``osc_prob_{N}nu_earth[_nsi|_liv]``
   * - sun
     - std / nsi / liv
     - --
     - ``osc_prob_{N}nu_sun[_nsi|_liv]``

Examples
----------

Three-flavor vacuum oscillation, full probability matrix (output captured
from this version):

.. code-block:: text

   $ magnus prob --flavors 3 --environment vacuum \
       --energy 1 --energy-unit GeV --baseline 1300 --baseline-unit km
   Magνs 1.0.0rc1 -- osc_prob_3nu_vacuum
   E = 1 GeV, L = 1300 km

               nu_e   nu_mu  nu_tau
   nu_e      0.9297  0.0085  0.0618
   nu_mu     0.0311  0.3885  0.5804
   nu_tau    0.0393  0.6029  0.3578

The same calculation, one channel only:

.. code-block:: text

   $ magnus prob --flavors 3 --environment vacuum --energy 1 --energy-unit GeV \
       --baseline 1300 --baseline-unit km --nu-i e --nu-f mu
   Magνs 1.0.0rc1 -- osc_prob_3nu_vacuum
   E = 1 GeV, L = 1300 km

   P = 0.0085

Earth crossing from the cosine of the zenith angle (equivalently, from two
named locations -- see ``--loc-ini``/``--loc-fin`` below):

.. code-block:: text

   $ magnus prob --flavors 3 --environment earth --energy 1 --energy-unit GeV \
       --costhz -0.8 --baseline 10193.6 --baseline-unit km
   Magνs 1.0.0rc1 -- osc_prob_3nu_earth
   E = 1 GeV, L = 10193.6 km

               nu_e   nu_mu  nu_tau
   nu_e      0.9128  0.0863  0.0009
   nu_mu     0.0629  0.6681  0.2690
   nu_tau    0.0243  0.2456  0.7301

Constant-density matter with non-standard interactions:

.. code-block:: text

   $ magnus prob --flavors 3 --environment matter --scenario nsi --rho 2.7 \
       --eps-ee 0.06 --eps-em -0.06 \
       --energy 1 --energy-unit GeV --baseline 1000 --baseline-unit km
   Magνs 1.0.0rc1 -- osc_prob_3nu_matter_nsi_constant_density
   E = 1 GeV, L = 1000 km

               nu_e   nu_mu  nu_tau
   nu_e      0.9898  0.0093  0.0009
   nu_mu     0.0093  0.9906  0.0001
   nu_tau    0.0009  0.0001  0.9990

Vacuum with a (deliberately large, for illustration) Lorentz-invariance-violating
term -- compare to the plain-vacuum result above at the same energy and baseline:

.. code-block:: text

   $ magnus prob --flavors 3 --environment vacuum --scenario liv \
       --sxi12 0.3 --b1 6e-13 --b2 1.2e-12 --liv-lambda 1e9 --n-liv 1 \
       --energy 1 --energy-unit GeV --baseline 1300 --baseline-unit km
   Magνs 1.0.0rc1 -- osc_prob_3nu_vacuum_liv
   E = 1 GeV, L = 1300 km

               nu_e   nu_mu  nu_tau
   nu_e      0.4971  0.0506  0.4523
   nu_mu     0.1341  0.7020  0.1639
   nu_tau    0.3688  0.2474  0.3838

A 3+2 sterile scenario (5 flavors), machine-readable output:

.. code-block:: text

   $ magnus prob --flavors 5 --environment earth --scenario liv \
       --costhz -0.8 --baseline 10193.6 --sxi12 0.2 --b1 1e-13 --liv-lambda 1e9 \
       --energy 1 --energy-unit GeV --json
   {
     "function": "osc_prob_5nu_earth_liv",
     "flavors": 5,
     "environment": "earth",
     "scenario": "liv",
     "nubar": false,
     "energy_eV": 1000000000.0,
     "baseline_eV-1": ...,
     "probability": [[...], [...], [...], [...], [...]]
   }

With ``s14 = s15 = s24 = s25 = s34 = s35 = 0`` (their defaults), the two
sterile states stay perfectly decoupled from the three active flavors and
from each other, as expected -- this exact check is one of the CLI's
regression tests (``tests/test_cli.py``).

Choosing a propagation strategy
------------------------------------

For a position-dependent Hamiltonian -- ``--environment sun``,
``--environment earth``, or ``--environment matter --density-profile exp`` --
``--strategy`` selects how the evolution operator is propagated, exactly as the
``strategy`` keyword does in the Python API (see :doc:`adiabatic_strategy` for
the full description of the three values).  It defaults to ``auto`` and is
ignored for vacuum and constant-density environments, whose Hamiltonians do not
depend on position at all.

This matters most for low-energy solar neutrinos, where the accumulated phase is
extreme.  ``magnus`` is not merely slower there -- it can hit its refinement caps
and return a confidently wrong number:

.. code-block:: bash

   magnus prob --flavors 3 --environment sun --energy 10 --energy-unit MeV \
       --baseline 626000 --nu-i e --nu-f e --strategy magnus

.. code-block:: text

   P = 0.6560

.. code-block:: bash

   magnus prob --flavors 3 --environment sun --energy 10 --energy-unit MeV \
       --baseline 626000 --nu-i e --nu-f e --strategy auto

.. code-block:: text

   P = 0.2905

The second value is the correct one.  Since ``auto`` is the default, you only
need this flag to *opt out* of the hybrid strategy (``--strategy magnus``, to
reproduce the older behavior) or to force it and be warned when it cannot
certify its own result (``--strategy hybrid``).

Errors are explicit rather than silent
------------------------------------------

Missing a required flag, or an invalid combination, produces a clear
message and a non-zero exit code instead of a wrong answer or a raw
traceback:

.. code-block:: text

   $ magnus prob --flavors 3 --environment vacuum --scenario nsi --energy 1 --baseline 1300
   magnus prob: --scenario nsi is not available with --environment vacuum (NSI couplings
   scale the matter potential, which vacuum has none of); use --environment matter/earth/sun
   instead.

   $ magnus prob --flavors 2 --environment vacuum --energy 1 --baseline 1300
   usage: magnus [-h] [-V] {prob} ...
   magnus: error: --sth and --dm2 are both required for --flavors 2.

Full flag reference
----------------------

The complete, current ``--help`` output (every flag is grouped by what it
configures):

.. code-block:: text

   usage: magnus prob [-h] [--flavors {2,3,4,5}] [--environment {vacuum,matter,earth,sun}]
                      [--scenario {std,nsi,liv}] [--density-profile {constant,exp}]
                      [--nubar] --energy ENERGY [--energy-unit {eV,keV,MeV,GeV,TeV,PeV}]
                      [--baseline BASELINE] [--l0 L0] [--baseline-unit {eV-1,km,cm}]
                      [--rho RHO] [--rho-central RHO_CENTRAL] [--l-scale L_SCALE]
                      [--density-unit {g/cm3,natural}] [--ratio-n-to-p RATIO_N_TO_P]
                      [--electron-fraction ELECTRON_FRACTION] [--costhz COSTHZ]
                      [--loc-ini LOC_INI] [--loc-fin LOC_FIN] [--sth STH] [--dm2 DM2]
                      [--s12 S12] [--s23 S23] [--s13 S13] [--dcp DCP] [--dm21 D21]
                      [--dm31 D31]
                      [--osc-params-set {OSC_PARAMS_DEFAULT,OSC_PARAMS_NU_FIT_6_0_SK_NO,OSC_PARAMS_NU_FIT_6_0_SK_IO}]
                      [--s14 S14] [--d14 D14] [--s24 S24] [--d24 D24] [--s34 S34]
                      [--dm41 D41] [--s15 S15] [--d15 D15] [--s25 S25] [--s35 S35]
                      [--d35 D35] [--dm51 D51] [--eps-aa EPS_AA] [--eps-ab EPS_AB]
                      [--eps-ee EPS_EE] [--eps-em EPS_EM] [--eps-et EPS_ET]
                      [--eps-mm EPS_MM] [--eps-mt EPS_MT] [--eps-tt EPS_TT]
                      [--eps-es EPS_ES] [--eps-ms EPS_MS] [--eps-ts EPS_TS]
                      [--eps-ss EPS_SS] [--eps-es1 EPS_ES1] [--eps-es2 EPS_ES2]
                      [--eps-ms1 EPS_MS1] [--eps-ms2 EPS_MS2] [--eps-ts1 EPS_TS1]
                      [--eps-ts2 EPS_TS2] [--eps-s1s1 EPS_S1S1] [--eps-s1s2 EPS_S1S2]
                      [--eps-s2s2 EPS_S2S2] [--sxi SXI] [--sxi12 SXI12] [--sxi23 SXI23]
                      [--sxi13 SXI13] [--dxicp DXICP] [--dxi13 DXI13] [--sxi14 SXI14]
                      [--dxi14 DXI14] [--sxi24 SXI24] [--dxi24 DXI24] [--sxi34 SXI34]
                      [--sxi15 SXI15] [--dxi15 DXI15] [--sxi25 SXI25] [--sxi35 SXI35]
                      [--dxi35 DXI35] [--b1 B1] [--b2 B2] [--b3 B3] [--b4 B4] [--b5 B5]
                      [--liv-lambda LAMBDA] [--n-liv N_LIV] [--nu-i NU_I] [--nu-f NU_F]
                      [--magnus-exp-order MAGNUS_EXP_ORDER]
                      [--integration-method {gl,trapezoid,simpson}] [--rtol RTOL]
                      [--atol ATOL] [--n-jobs N_JOBS] [--strategy {auto,hybrid,magnus}]
                      [--verbose {0,1,2}] [--json] [--precision PRECISION]

   options:
     -h, --help            show this help message and exit

   Environment:
     --flavors {2,3,4,5}   Number of neutrino flavors (default: 3).
     --environment {vacuum,matter,earth,sun}
                           Propagation environment (default: vacuum).
     --scenario {std,nsi,liv}
                           Physics scenario on top of the environment: 'std' (Standard
                           Model), 'nsi' (non-standard interactions), or 'liv' (Lorentz-
                           invariance violation). 'nsi' is not available with --environment
                           vacuum. Default: std.
     --density-profile {constant,exp}
                           Matter density profile, only used with --environment matter:
                           'constant' (requires --rho) or 'exp' (requires --rho-central and
                           --l-scale). Default: constant.
     --nubar               Compute the probability for antineutrinos instead of neutrinos.

   Energy and baseline:
     --energy ENERGY       Neutrino energy.
     --energy-unit {eV,keV,MeV,GeV,TeV,PeV}
                           Unit of --energy (default: GeV).
     --baseline BASELINE   Baseline / final position. Required for vacuum, matter, and sun,
                           and for earth when using --costhz. Only computed automatically
                           for earth when both --loc-ini and --loc-fin are given instead.
     --l0 L0               Initial position (used by --environment sun and --density-
                           profile exp). Default: 0.0.
     --baseline-unit {eV-1,km,cm}
                           Unit of --baseline, --l0, and --l-scale (default: km).

   Matter (--environment matter):
     --rho RHO             Matter density (constant profile).
     --rho-central RHO_CENTRAL
                           Matter density at the center of the profile, l=0 (exponential
                           profile).
     --l-scale L_SCALE     Length scale of the exponential density decrease (exponential
                           profile).
     --density-unit {g/cm3,natural}
                           Unit of --rho/--rho-central: g/cm3 (converted internally) or
                           natural units (eV^4). Default: g/cm3.
     --ratio-n-to-p RATIO_N_TO_P
                           Ratio of the number of neutrons to protons in matter. Default:
                           1.0.
     --electron-fraction ELECTRON_FRACTION
                           Electron fraction of matter. Default: 0.5.

   Earth (--environment earth):
     --costhz COSTHZ       Cosine of the neutrino zenith angle.
     --loc-ini LOC_INI     Initial location name (e.g. fermilab); see
                           magnus.earth.loc_coords_dms. Must be given together with --loc-
                           fin, as an alternative to --costhz.
     --loc-fin LOC_FIN     Final location name; see --loc-ini.

   Standard oscillation parameters (2-flavor):
     --sth STH             Sin(theta) (required for --flavors 2).
     --dm2 DM2             Mass-squared difference Delta m^2 (required for --flavors 2).

   Standard oscillation parameters (3+ flavors):
     --s12 S12             Sin(theta_12). Default: NuFit 6.0.
     --s23 S23             Sin(theta_23). Default: NuFit 6.0.
     --s13 S13             Sin(theta_13). Default: NuFit 6.0.
     --dcp DCP             delta_CP [radian]. Default: NuFit 6.0.
     --dm21 D21            Mass-squared difference Delta m^2_21. Default: NuFit 6.0.
     --dm31 D31            Mass-squared difference Delta m^2_31. Default: NuFit 6.0.
     --osc-params-set {OSC_PARAMS_DEFAULT,OSC_PARAMS_NU_FIT_6_0_SK_NO,OSC_PARAMS_NU_FIT_6_0_SK_IO}
                           Predefined set used to fill in any of s12/s23/s13/dCP/D21/D31
                           left unspecified: normal ordering (..._NO, the default) or
                           inverted ordering (..._IO).

   Additional sterile mixing (4+ flavors):
     --s14 S14             Sin(theta_14). Default: 0.0.
     --d14 D14             delta_14 [radian]. Default: 0.0.
     --s24 S24             Sin(theta_24). Default: 0.0.
     --d24 D24             delta_24 [radian]. Default: 0.0.
     --s34 S34             Sin(theta_34). Default: 0.0.
     --dm41 D41            Mass-squared difference Delta m^2_41. Default: 0.0.

   Additional sterile mixing (5 flavors):
     --s15 S15             Sin(theta_15). Default: 0.0.
     --d15 D15             delta_15 [radian]. Default: 0.0.
     --s25 S25             Sin(theta_25). Default: 0.0.
     --s35 S35             Sin(theta_35). Default: 0.0.
     --d35 D35             delta_35 [radian]. Default: 0.0.
     --dm51 D51            Mass-squared difference Delta m^2_51. Default: 0.0.

   NSI parameters (--scenario nsi):
     --eps-aa EPS_AA       2-flavor diagonal NSI coupling.
     --eps-ab EPS_AB       2-flavor off-diagonal NSI coupling.
     --eps-ee EPS_EE       Diagonal NSI coupling of nu_e.
     --eps-em EPS_EM       Off-diagonal (e-mu) NSI coupling.
     --eps-et EPS_ET       Off-diagonal (e-tau) NSI coupling.
     --eps-mm EPS_MM       Diagonal NSI coupling of nu_mu.
     --eps-mt EPS_MT       Off-diagonal (mu-tau) NSI coupling.
     --eps-tt EPS_TT       Diagonal NSI coupling of nu_tau.
     --eps-es EPS_ES       (4nu) Off-diagonal (e-s) NSI coupling.
     --eps-ms EPS_MS       (4nu) Off-diagonal (mu-s) NSI coupling.
     --eps-ts EPS_TS       (4nu) Off-diagonal (tau-s) NSI coupling.
     --eps-ss EPS_SS       (4nu) Diagonal NSI coupling of nu_s.
     --eps-es1 EPS_ES1     (5nu) Off-diagonal (e-s1) NSI coupling.
     --eps-es2 EPS_ES2     (5nu) Off-diagonal (e-s2) NSI coupling.
     --eps-ms1 EPS_MS1     (5nu) Off-diagonal (mu-s1) NSI coupling.
     --eps-ms2 EPS_MS2     (5nu) Off-diagonal (mu-s2) NSI coupling.
     --eps-ts1 EPS_TS1     (5nu) Off-diagonal (tau-s1) NSI coupling.
     --eps-ts2 EPS_TS2     (5nu) Off-diagonal (tau-s2) NSI coupling.
     --eps-s1s1 EPS_S1S1   (5nu) Diagonal NSI coupling of nu_s1.
     --eps-s1s2 EPS_S1S2   (5nu) Off-diagonal (s1-s2) NSI coupling.
     --eps-s2s2 EPS_S2S2   (5nu) Diagonal NSI coupling of nu_s2.

   LIV parameters (--scenario liv):
     --sxi SXI             2-flavor LIV mixing angle sine.
     --sxi12 SXI12         LIV mixing angle sine xi_12.
     --sxi23 SXI23         LIV mixing angle sine xi_23.
     --sxi13 SXI13         LIV mixing angle sine xi_13.
     --dxicp DXICP         (3nu) LIV CP-violation phase [radian].
     --dxi13 DXI13         (4/5nu) LIV CP-violation phase [radian] (replaces --dxicp).
     --sxi14 SXI14         (4/5nu) LIV mixing angle sine xi_14.
     --dxi14 DXI14         (4/5nu) LIV CP-violation phase [radian].
     --sxi24 SXI24         (4/5nu) LIV mixing angle sine xi_24.
     --dxi24 DXI24         (4/5nu) LIV CP-violation phase [radian].
     --sxi34 SXI34         (4/5nu) LIV mixing angle sine xi_34.
     --sxi15 SXI15         (5nu) LIV mixing angle sine xi_15.
     --dxi15 DXI15         (5nu) LIV CP-violation phase [radian].
     --sxi25 SXI25         (5nu) LIV mixing angle sine xi_25.
     --sxi35 SXI35         (5nu) LIV mixing angle sine xi_35.
     --dxi35 DXI35         (5nu) LIV CP-violation phase [radian].
     --b1 B1               LIV eigenvalue b1.
     --b2 B2               LIV eigenvalue b2.
     --b3 B3               LIV eigenvalue b3.
     --b4 B4               LIV eigenvalue b4.
     --b5 B5               LIV eigenvalue b5.
     --liv-lambda LAMBDA   LIV energy scale Lambda. Default: 1.0.
     --n-liv N_LIV         Power of the energy dependence of the LIV operator. Default: 0.

   Channel selection:
     --nu-i NU_I           Initial flavor (index or name: e, mu, tau, s, s1, s2). If given
                           with --nu-f, prints a single probability instead of the full
                           matrix.
     --nu-f NU_F           Final flavor; see --nu-i.

   Advanced numerics:
     --magnus-exp-order MAGNUS_EXP_ORDER
                           Highest order of the Magnus expansion (1-6). Default: 4.
     --integration-method {gl,trapezoid,simpson}
                           Quadrature method. 'gl' (Gauss-Legendre collocation) needs only
                           1-3 Hamiltonian evaluations per slab and matches its quadrature
                           order to the expansion order, so it is both the fastest and the
                           most accurate for a smooth Hamiltonian. 'trapezoid'/'simpson'
                           sample a uniform grid of --n-tpts-per-slab points instead, and
                           are the safer choice if the Hamiltonian is not smooth within a
                           slab. Default: gl.
     --rtol RTOL           Relative tolerance on the agreement between successive
                           refinement levels -- a stopping rule, not a guaranteed accuracy.
                           Default: 1e-3.
     --atol ATOL           Absolute tolerance on the same agreement; see --rtol. Default:
                           1e-3.
     --n-jobs N_JOBS       Number of parallel joblib workers. Default: 1.
     --strategy {auto,hybrid,magnus}
                           How to propagate a position-dependent Hamiltonian: 'magnus' uses
                           only the Magnus-expansion machinery; 'hybrid' also tries
                           adiabatic transport with a Magnus patch at each non-adiabatic
                           window, warning if it cannot certify the result; 'auto' tries
                           hybrid and falls back to magnus silently. Ignored for vacuum and
                           constant-density environments. Default: auto.
     --verbose {0,1,2}     Verbosity level. Default: 0.

   Output:
     --json                Print the result as JSON instead of a table.
     --precision PRECISION
                           Decimal digits shown in table output. Default: 4.

Implementation notes
-----------------------

``magnus prob`` does not reimplement any physics: it builds a keyword-argument
dictionary from the flags you passed and calls straight into the matching
``osc_prob_{N}nu_*`` function (see :func:`magnus.cli.main` and
``_wrapper_name``). Physics keyword arguments that a given
function does not explicitly accept (for example, ``default_osc_params_set_name``
is not defined on every LIV wrapper) are filtered out via
:func:`inspect.signature` before the call, rather than being silently
forwarded through ``**kwargs`` to a layer that does not expect them -- see
``_call``.
