# -*- coding: utf-8 -*-
"""oscprob.py

Contains routines to compute the neutrino oscillation probability.

Internally, the probability is computed using Magnus expansion, but the
user does not call the routines in the :py:mod:`magnus.magnus` module
directly. Instead, the user calls the :func:`osc_prob`, which calls the
Magnus expansion routines internally. The function :func:`osc_prob` is
generic, flexible, and computationally efficient. 

- :func:`osc_prob`: Primordial function to compute the oscillation
  probability, for any given Hamiltonian, either time-dependent or 
  -independent (or, equivalently, position-dependent or -independent). 
  Supports arbitrary number of neutrino flavors.

The module contains additional functions that are wrappers of 
:func:`osc_prob` to compute commonly studied cases.

Neutrino oscillations in **vacuum**:

- :func:`osc_prob_2nu_vacuum`: Two-neutrino oscillation probabilities.

- :func:`osc_prob_3nu_vacuum`: Three-flavor oscillation probabilities.

- :func:`osc_prob_4nu_vacuum`: One 
  additional flavor (i.e., 3+1 sterile neutrino model).

- :func:`osc_prob_5nu_vacuum`: Two 
  additional flavors (i.e., 3+2 sterile neutrino model).

Neutrino oscillations in **constant-density matter**:

- :func:`osc_prob_2nu_matter_constant_density`: Two-neutrino oscillation
  probabilities.

- :func:`osc_prob_3nu_matter_constant_density`: Three-neutrino
  oscillation probabilities.

- :func:`osc_prob_4nu_matter_constant_density`: One additional flavor 
  (i.e., 3+1 sterile model).

- :func:`osc_prob_5nu_matter_constant_density`: Two additional flavors 
  (i.e., 3+2 sterile model).

Neutrino oscillations in **exponentially falling matter density 
profile** (e.g., in a supernova or the Sun):

- :func:`osc_prob_2nu_matter_exp_density`: Two-neutrino oscillation
  probabilities.

- :func:`osc_prob_3nu_matter_exp_density`: Three-neutrino oscillation
  probabilities.

- :func:`osc_prob_4nu_matter_exp_density`: One additional flavor. Matter
  potential affects only :math:`\\nu_e`.

- :func:`osc_prob_5nu_matter_exp_density`: Two additional flavors.
  Matter potential affects only :math:`\\nu_e`.

Neutrino oscillations between any two locations on the surface of the
**Earth**, useful for long-baseline neutrino experiments:

- :func:`osc_prob_2nu_earth`: Two-neutrino oscillation probabilities.  

- :func:`osc_prob_3nu_earth`: Three-neutrino oscillation probabilities. 

- :func:`osc_prob_4nu_earth`: One additional flavor. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_5nu_earth`: Two additional flavors. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_earth`: Oscillation probabilities for arbitrary number
  of flavors and arbitrary Hamiltonian.  Does not assume standard 
  oscillations.

.. note::
   These routines use the `Preliminary Reference Earth Model 
   <https://www.cfa.harvard.edu/~lzeng/papers/PREM.pdfL>`_ for the 
   matter density profile inside Earth.  To use a different density 
   profile (including also profiles for bodies other than the Earth), 
   use the primordial function :func:`osc_prob` instead.

Neutrino oscillations in the **Sun**:

- :func:`osc_prob_2nu_sun`: Two-neutrino oscillation probabilities.  

- :func:`osc_prob_3nu_sun`: Three-neutrino oscillation probabilities. 

- :func:`osc_prob_4nu_sun`: One additional flavor. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_5nu_sun`: Two additional flavors. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_sun`: Oscillation probabilities for arbitrary number
  of flavors and arbitrary Hamiltonian.  Does not assume standard 
  oscillations.

.. note::
   These routines use a simple exponentially falling function of radial
   distance for the matter density inside the Sun: :math:`N_e(r) = 
   N_e(0) \\exp(-r/r_0)`, with 
   :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
   :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
   `Fundamentals of Neutrino Physics and Astrophysics 
   <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
   Wook Kim.

   To use a different density profile, use the primordial function 
   :func:`osc_prob` instead.

Functions designed for specific **beyond-the-Standard-Model** proposals:

- Non-standard neutrino interactions (NSI):

  - :func:`osc_prob_2nu_matter_nsi_constant_density`
  - :func:`osc_prob_3nu_matter_nsi_constant_density`
  - :func:`osc_prob_4nu_matter_nsi_constant_density`
  - :func:`osc_prob_5nu_matter_nsi_constant_density`
  - :func:`osc_prob_2nu_matter_nsi_exp_density`
  - :func:`osc_prob_3nu_matter_nsi_exp_density`
  - :func:`osc_prob_4nu_matter_nsi_exp_density`
  - :func:`osc_prob_5nu_matter_nsi_exp_density`
  - :func:`osc_prob_2nu_earth_nsi`
  - :func:`osc_prob_3nu_earth_nsi`
  - :func:`osc_prob_4nu_earth_nsi`
  - :func:`osc_prob_5nu_earth_nsi`
  - :func:`osc_prob_2nu_sun_nsi`
  - :func:`osc_prob_3nu_sun_nsi`
  - :func:`osc_prob_4nu_sun_nsi`
  - :func:`osc_prob_5nu_sun_nsi`

- Lorentz-invariance violation:

  - :func:`osc_prob_2nu_matter_liv_constant_density`
  - :func:`osc_prob_3nu_matter_liv_constant_density`
  - :func:`osc_prob_4nu_matter_liv_constant_density`
  - :func:`osc_prob_5nu_matter_liv_constant_density`
  - :func:`osc_prob_2nu_matter_liv_exp_density`
  - :func:`osc_prob_3nu_matter_liv_exp_density`
  - :func:`osc_prob_4nu_matter_liv_exp_density`
  - :func:`osc_prob_5nu_matter_liv_exp_density`
  - :func:`osc_prob_2nu_earth_liv`
  - :func:`osc_prob_3nu_earth_liv`
  - :func:`osc_prob_4nu_earth_liv`
  - :func:`osc_prob_5nu_earth_liv`
  - :func:`osc_prob_2nu_sun_liv`
  - :func:`osc_prob_3nu_sun_liv`
  - :func:`osc_prob_4nu_sun_liv`
  - :func:`osc_prob_5nu_sun_liv`

Examples
--------

.. seealso::
   Find many more examples, including advanced applications and plots,
   in the `Jupyter notebooks 
   <https://github.com/mbustama/Magnus/tree/main/notebooks>`_ that are 
   distributed with :math:`{\\rm Mag}{\\nu}s`.

The code blocks below run when these docs are built, so the output shown
is always current.

.. jupyter-execute::

    import numpy as np

    import magnus.oscprob as oscprob
    import magnus.globaldefs as gd

    np.set_printoptions(precision=3)

    # Warnings are normally prefixed with an ANSI-colored "Warning:", which is
    # meant for a terminal and renders as escape-code noise in HTML.  These docs
    # therefore switch to plain text; in a terminal, leave this alone.
    gd.set_color_output(False)

Calling :func:`osc_prob_3nu_vacuum` returns a :math:`3 \\times 3` NumPy array
of probabilities whose entry ``[i][j]`` is the probability of a neutrino
produced with flavor ``i`` being detected with flavor ``j``

For a single neutrino energy and baseline:

.. jupyter-execute::

    baseline = 10.0 * gd.UNIT_KM  # 10 km in natural units [eV^-1]
    energy = 1.0 * gd.UNIT_MEV    # [eV]

    oscprob.osc_prob_3nu_vacuum(energy, baseline)

The probabilities returned by :func:`osc_prob_3nu_vacuum` (and also
:func:`osc_prob_2nu_vacuum`, 
:func:`osc_prob_2nu_matter_constant_density`, and
:func:`osc_prob_3nu_matter_constant_density`) are returned with machine
(or NumPy) precision, since first-order Magnus expansion is enough to 
compute them.

Pick one channel only, e.g., :math:`\\nu_e \\to \\nu_\\mu`, by passing
an initial flavor, ``nu_i``, and a final flavor ``nu_f``:

.. jupyter-execute::

    oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU)

The flavor indices ``NUE``, ``NUMU``, and ``NUTAU`` are defined in the
:py:mod:`magnus.globaldefs` module. For anti-neutrinos, i.e.,
:math:`\\bar{\\nu}_e \\to \\bar{\\nu}_\\mu`:

.. jupyter-execute::

    oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU,
                                nubar=True)

Calling :func:`osc_prob_3nu_vacuum` without specifying the values of the
oscillation parameters will compute probabilities using the default 
values in :math:`{\\rm Mag}{\\nu}s` (see 
``gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']``.)

We can specify values of the oscillation parameters. Unspecified values
are set to their defaults (pass nonzero ``verbose`` to see this and 
other warnings):

.. jupyter-execute::

    oscprob.osc_prob_3nu_vacuum(energy, baseline, s12=0.0, verbose=1)

Fixed energy, multiple baselines:

.. jupyter-execute::

    baselines = gd.UNIT_KM * np.array([1.0, 10.0, 100.0])

    oscprob.osc_prob_3nu_vacuum(energy, baselines, nu_i=gd.NUE, nu_f=gd.NUMU)

Fixed baseline, multiple energies:

.. jupyter-execute::

    energies = gd.UNIT_MEV * np.array([1.0, 10.0, 100.0])

    oscprob.osc_prob_3nu_vacuum(energies, baseline, nu_i=gd.NUE, nu_f=gd.NUMU)

Multiple energies and baselines:

.. jupyter-execute::

    oscprob.osc_prob_3nu_vacuum(energies, baselines, nu_i=gd.NUE, nu_f=gd.NUMU)

To compute the oscillation probabilities in constant-density matter, we
need to specify the matter density, ``rho``, i.e.,

.. jupyter-execute::

    rho = 10.0 * gd.UNIT_G_PER_CM3

    oscprob.osc_prob_3nu_matter_constant_density(energy, baseline, rho,
                                                 nu_i=gd.NUE, nu_f=gd.NUMU)

To compute oscillation probabilities for a time-dependent Hamiltonian,
we need to call :func:`osc_prob` directly which, while still 
straightforward, requires us to pass a Hamiltonian function explicitly.

For instance, for density matter profile that is exponentially falling
with distance:

.. hint::
   There is a good chance that the scenario you are interested in 
   calculating was already developed in the :math:`{\\rm Mag}{\\nu}s` 
   `Jupyter 
   notebooks 
   <https://github.com/mbustama/Magnus/tree/main/notebooks>`_.  
   
   Worked-out examples include: oscillations in various matter density 
   profiles, in the Earth, and in the Sun, oscillograms, biprobability 
   plots, and new-physics models like additional neutrino flavors (3+1 
   and 3+2 sterile neutrino models), non-standard neutrino interactions,
   and Lorentz-invariance violation.

See :doc:`/architecture` for how the ``osc_prob_*`` functions listed above
are layered internally (primordial/middle/wrapper) and how to add a new one.
"""

__author__ = 'Mauricio Bustamante'


import numpy as np
import os
import sys
import warnings
from contextlib import contextmanager
from functools import reduce
from joblib import Parallel, delayed
from typing import Optional, Callable, Union, Tuple, Dict
from io import TextIOWrapper
from inspect import signature
# import numba as nb

# import numpy.typing

import magnus.magnus as magnus
import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.adiabatic as adiabatic
import magnus.avgprob as avgprob
from magnus import version
from magnus import authors


has_magnus_header_been_printed = False


MAX_N_SLABS_DEFAULT = {'gl': 20000, 'trapezoid': 2000, 'simpson': 2000}
r"""dict: Module-level constant

Default cap on the number of slabs, per ``integration_method``, used when
``max_n_slabs`` is left as None.

The cap exists to bound cost, and cost per slab differs by more than an order of
magnitude between the two families of integrators.  ``'gl'`` evaluates the Hamiltonian
1, 2, or 3 times per slab (set by the expansion order), while ``'trapezoid'`` and
``'simpson'`` evaluate it ``n_tpts_per_slab`` times (100 by default, up to 500).  A
single cap tuned for one family therefore starves the other: at 2000 slabs -- the
quadrature cap, unchanged here -- ``'gl'`` was hitting the ceiling on cases it could
resolve comfortably, and reporting that it could not verify convergence, on answers
that were in fact far more accurate than the quadrature methods managed within their
own cap.

The ``'gl'`` value is set from both directions.  The hardest case in the validation
suite (5 flavors, eV-scale sterile splittings, Earth-crossing baseline) converges at
about 8,600 slabs, so 20000 leaves better than a factor of two of headroom.  At the
same time, 20000 slabs at 2-3 nodes is roughly 40,000-60,000 Hamiltonian evaluations,
still well under the ~200,000 that 2000 quadrature slabs at 100 points per slab
already permit -- so the more generous cap is also the cheaper worst case.

Passing ``max_n_slabs`` explicitly always wins; this is only the fallback.

.. versionadded:: 1.0.0
"""


IP_EXP_N_SLABS_CAP = 2_000_000
r"""int: Module-level constant

Slab ceiling for the closed-form interaction-picture integrator
(``_osc_prob_ip_exp_core``), deliberately decoupled from the caller's
``max_n_slabs``, which is calibrated for the far more expensive quadrature slabs of
the general method.  Each slab here costs one 2x2 eigendecomposition, so the ceiling
is affordable in *time*.

It is not, on its own, affordable in *memory*: the integrator works on arrays of
shape ``(n_energies, n_slabs, d, d)``, so the working set scales with the number of
energies as well as the slab count, and at this ceiling it reached ~1.3 GB **per
energy** -- enough for a batched solar call to exhaust the machine.  The working set
is now tiled to :data:`BATCH_WORKING_ENTRIES` independently of both, so this ceiling
bounds only time.  See ``docs/dev/BUG_IP_EXP_MEMORY.md``.

Named rather than inlined so that the method's give-up behaviour at the ceiling can
be exercised by a test at a small cap: reaching two million slabs to check what
happens at the boundary would cost minutes, so with the value buried in the function
body those branches could not be tested at all.

.. versionadded:: 1.0.0
"""


BATCH_WORKING_ENTRIES = 4_194_304
r"""int: Module-level constant

Ceiling on the number of complex entries in any one temporary array of the batched
scan engines -- about 64 MB at 16 bytes each.  Both batched engines work on arrays
indexed by (energy, slab, ...), whose size is the product of quantities the caller
controls independently, so neither a slab cap nor an energy count bounds them on its
own.  Tiling against a fixed entry budget does, and it makes peak memory a property
of the library rather than of the call.

The value is a compromise: large enough that the tiles stay well inside the regime
where batched BLAS calls amortize their overhead, small enough to be invisible next
to the result arrays on any machine that can hold them.

.. versionadded:: 1.0.0
"""


CUMULATIVE_AUTO_MIN_POINTS = 2
r"""int: Module-level constant

Fewest baselines at which ``cumulative='auto'`` engages the cumulative scan in
:func:`osc_prob_energy_baseline`.

A single baseline has no prefix to reuse, and would pay for the adaptive probe that sizes the
inherited grid without getting anything back -- which matters because every single-point call
through the wrapper layer is served by ``osc_prob_energy_baseline``.  From two baselines
upward there is something to share.

The threshold is deliberately not set at the point where the cumulative scan becomes *faster*,
which is higher (measured against ``solve_ivp``, on a 5 MeV solar scan to one solar radius:
0.75x at N = 2, 0.87x at N = 10, 2.65x at N = 25, 84x at N = 1000).  Below that crossover the
cumulative scan is at most ~1.3x slower in wall time while being one to three orders of
magnitude more accurate -- and the per-point path it replaces returns answers *outside* the
requested 1e-3 there (9.7e-3 at N = 10, 5.6e-3 at N = 25, 2.6e-3 at N = 100).  Trading a few
milliseconds for that is the right way round.

.. versionadded:: 1.0.0
"""


# A weak-band self-cross-check was built here and REMOVED, on measurement.  The idea: below
# the N = 25 seam, where every remaining silent miss lives, have strategy='auto' verify a
# window-free hybrid result against the general Magnus ladder instead of taking its word --
# certification with no window rests on gamma alone, and GAMMA_TO_ERROR is good only to ~2x.
#
# It does not work, and the reason is worth keeping:
#
#   * with the shipped constant, over 200 random smooth profiles, 25 results were window-free
#     and certified and the ladder agreed with ALL 25.  Zero disagreements, so zero benefit
#     against a measured 9% of calls paying for a second engine.
#   * with GAMMA_TO_ERROR deliberately made optimistic by 2x, 41 results were certified,
#     3 of them genuinely outside tolerance -- and the check still fired zero times.  The
#     first version's trigger was computed FROM GAMMA_TO_ERROR, so mis-calibrating the
#     constant shrank the trigger in step with it: a check insured against a constant, keyed
#     on that same constant.  Self-referential, which is precisely the failure shape this
#     whole programme exists to find, reproduced by the person writing the fix.
#   * removing the circularity does not rescue it.  Verifying EVERY window-free result --
#     11% of random calls, and 100% of ordinary solar single points, which are window-free --
#     still fires zero times and still misses the same 3.
#
# The conclusion is structural: what is left in the weak band is not disagreement between
# engines, it is the engines being wrong TOGETHER.  A cross-check detects the former by
# construction and can never detect the latter -- the same reason cross_check_strategies
# cannot see a feature narrower than every grid.  The instrument that does reach that class is
# adiabatic.find_hidden_features, which looks at the profile rather than at the answers.
#
# Reproduce: docs/dev/adversarial_batteries/crosscheck_benefit.py and weak_band.py.

HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25
r"""int: Module-level constant

Fewest baselines at which ``_osc_prob_hybrid_dispatch`` stands aside, under
``strategy='auto'``, so that a single-energy baseline scan reaches the cumulative scan instead.

Deliberately **larger** than :data:`CUMULATIVE_AUTO_MIN_POINTS`, because the two thresholds
guard different trades.  On the primordial entry point the cumulative scan's alternative is the
general per-point path, which is not merely slower but returns answers outside the requested
tolerance on solar profiles (9.7e-3 at N = 10); taking the cumulative scan from N = 2 is right
there.  Reached through the wrapper layer the alternative is the *hybrid* strategy, which is
both accurate (~1e-5, well inside a requested 1e-3) and cheap per point (~20 ms), so the
cumulative scan's near-constant cost -- dominated by its strict probe -- only pays off once
there are enough points to amortise it.

Measured through ``osc_prob_2nu_sun`` on a solar profile, hybrid against the cumulative scan:

===== =============== ===============
N     5 MeV           10 MeV
===== =============== ===============
2     37 ms / 283 ms  39 ms / 297 ms
10    184 ms / 282 ms 197 ms / 454 ms
25    464 ms / 263 ms 489 ms / 432 ms
400   7459 ms / 292 ms 7904 ms / 368 ms
===== =============== ===============

The crossover sits near N = 14 at 5 MeV and N = 22 at 10 MeV; 25 clears both.  Erring high only
forgoes a speed-up, whereas erring low makes small scans several times slower for accuracy that
was already two orders inside what was asked for.

Being the larger of the two thresholds also keeps the fall-through safe: whenever the hybrid
dispatcher declines on this count, ``cumulative='auto'`` is guaranteed to engage, so a scan can
never decline both paths and land on the general per-point method.

**Accuracy steps at this threshold, and it is a large step.**  Because the two sides are
different methods rather than two settings of one method, adding a single baseline to a
24-point scan can change every answer in it.  Measured against ``solve_ivp``, ``err(N=24)``
against ``err(N=26)``:

========================= ========== ========== ==========
profile                   N = 24     N = 26     step
========================= ========== ========== ==========
solar exponential         3.30e-05   2.13e-08   1 546x
noisy                     6.27e-04   1.04e-08   60 418x
multi-resonance           1.58e-03   2.86e-09   552 945x
========================= ========== ========== ==========

The step is always *toward* the truth -- above the threshold the scan is more accurate, not
less -- so this is a discontinuity to know about rather than a defect.  But a caller sweeping N
and watching the answer move by five orders of magnitude at N = 25 is seeing the routing change,
not a numerical instability.  Pass ``cumulative=True`` to take the cumulative scan below the
threshold as well, or ``cumulative=False`` to stay off it entirely.

.. versionadded:: 1.0.0
"""


CUMULATIVE_N_ACC_SAFETY = 4
r"""int: Module-level constant

Multiple of the inherited slab count used for the accuracy grid of a cumulative baseline
scan (``osc_prob_energy_baseline(..., cumulative=True)``).

The grid is sized from one adaptive :func:`osc_prob` call at the longest baseline, which
reports the slab count *that* baseline needed.  Applied unmultiplied, the same uniform density
is thinner than what a per-point path would have chosen for the shorter baselines in the scan,
and the result -- while inside the requested tolerance -- comes out less accurate than the path
it replaces.  Measured on a 1000-point solar scan against ``solve_ivp``, where the per-point
path takes 12.0 s for an error of 5.6e-5:

===========  ==========  =========  ==========
safety       ``n_acc``   time       error
===========  ==========  =========  ==========
1            14 883      0.049 s    2.35e-04
2            29 766      0.097 s    5.10e-06
**4**        **59 532**  0.173 s    **3.34e-07**
8            119 064     0.346 s    1.80e-08
===========  ==========  =========  ==========

**Why four rather than two.**  Two was chosen when the cumulative scan's only alternative was
the general per-point path, against which it was already 124x faster and 11x more accurate.
Since the dispatch chain routes wrapper baseline scans here (see
``HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS``), the alternative is the *hybrid* strategy instead,
which is considerably more accurate than the per-point path -- so the bar moved.  At two, a
48-configuration sweep found three where the cumulative scan was less accurate than the hybrid
answer it replaced, all at high energy over a short baseline:

=========================  ==========  ==========  ==========
configuration              hybrid      safety 2    safety 4
=========================  ==========  ==========  ==========
60 MeV, N = 150, 0.4 R_sun  1.57e-05   5.03e-05    8.56e-07
100 MeV, N = 150, 0.4 R_sun 2.51e-05   3.77e-05    6.11e-07
100 MeV, N = 40, 0.4 R_sun  9.13e-06   1.10e-05    5.58e-07
=========================  ==========  ==========  ==========

Four removes all three and beats the hybrid answer on each, while improving the unaffected
configurations by roughly twenty times as well (5 MeV, N = 150: 8.4e-07 -> 3.7e-08).  It costs
about 1.4x in wall time -- 28 ms -> 40 ms, 260 ms -> 366 ms on the cases above -- against a
path it is still tens of times faster than.  Eight is better again but 2.4x, and buys accuracy
no longer needed to clear the bar.

Note that the error is **not** concentrated where the shape of this constant suggests: on the
60 MeV case it sits at the *longest* baselines (5.03e-05 there against 2.18e-06 over the
shortest third), and the grid density at the short end already matches what a probe there would
ask for to within 1%.  What the multiplier buys is total resolution, not better placement.

**The probe does not always converge**, and then this multiplier is doing more work than its
name suggests.  Over a full solar radius at 5 and 10 MeV the strict probe reaches
``max_n_slabs`` (20 000) without two successive levels agreeing, so the count it reports is the
cap rather than a converged requirement, and ``n_acc`` is 80 000 by way of a ceiling.  The
resulting scans are accurate (~5e-08 measured against ``solve_ivp``), but the safety margin is
what makes that so.  A caller who lowers ``max_n_slabs`` lowers the scan's resolution with it,
in proportion and without a separate warning.

.. versionadded:: 1.0.0
"""


OUTPUT_GUARD_MIN_BYTES = 64*1024*1024
r"""int: Module-level constant

Requested-output size below which :func:`osc_prob_energy_baseline` does not bother
checking whether the result will fit in memory.  The check itself costs one integer
multiply below this threshold and a single read of the operating system's free-memory
figure above it, so the floor exists to keep even that off the path of ordinary calls.

.. versionadded:: 1.0.0
"""


OUTPUT_GUARD_SAFETY = 2.0
r"""float: Module-level constant

Multiple of the requested result size that must fit in available memory before
:func:`osc_prob_energy_baseline` will attempt a scan.  Greater than one because the
adaptive engines hold at least the current and the previous probability matrices at
once, plus the caller's own input arrays.

Deliberately not larger: the guard exists to turn an out-of-memory kill into a
diagnosable error, not to second-guess a caller who knows their machine.  It refuses
only when the answer alone would claim more than half of what is free.

.. versionadded:: 1.0.0
"""


IP_EXP_LOOP_CAP = 30
r"""int: Module-level constant

Maximum number of refinement loops in ``_osc_prob_ip_exp_core``.

Note this ceiling is not what stops the loop in practice: the slab count doubles
each pass, so it reaches :data:`IP_EXP_N_SLABS_CAP` after about twenty passes and
returns there, well before a thirtieth pass could occur.  It is a backstop against a
future change to the growth factor or the slab ceiling, not a live limit.

.. versionadded:: 1.0.0
"""


def _n_required_params(func):
    r"""How many arguments ``func`` obliges its caller to supply positionally.

    Every ``H_func``/``rho_func`` arity check in this module used ``len(signature(f).parameters)``,
    which counts keyword parameters that already have defaults.  That breaks the ordinary Python
    idiom for binding a loop variable into a closure --

    .. code-block:: python

        def H(energy, l, VCC, _hvac=hvac, _proj=proj):   # 5 parameters, 3 required
            ...

    -- which the package's own documentation recommends the *factory* form of, precisely because
    this form used to fail.  With ``validate_input=True`` it raised "must be a function of either
    three arguments (energy, l, VCC) or two arguments (energy, l); the provided H_func takes 5";
    with ``validate_input=False`` it silently took the two-argument branch and died with a
    ``TypeError`` from inside the engine.  Neither is the user's fault.

    Counting required parameters instead makes both forms work and changes nothing for a function
    written without defaults.  A ``*args`` function is not counted this way -- it declares no
    required parameters at all, and the old total is the better guess there -- so those keep the
    previous behaviour.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    func : Callable
        The function to inspect.

    Returns
    -------
    int
        Number of positional parameters without defaults, or the total parameter count when
        ``func`` takes ``*args``.
    """
    params = list(signature(func).parameters.values())
    if any(q.kind is q.VAR_POSITIONAL for q in params):
        return len(params)
    return sum(1 for q in params
               if q.kind in (q.POSITIONAL_ONLY, q.POSITIONAL_OR_KEYWORD)
               and q.default is q.empty)


def _resolve_max_n_slabs(max_n_slabs, integration_method):
    """Fills in the per-method default cap when ``max_n_slabs`` is None.

    An explicitly passed value always wins.  Unknown method names fall back to the
    quadrature cap, so an invalid ``integration_method`` still fails in the validator that
    is meant to report it, rather than here with a KeyError.
    """
    if max_n_slabs is not None:
        return max_n_slabs
    return MAX_N_SLABS_DEFAULT.get(integration_method, MAX_N_SLABS_DEFAULT['trapezoid'])


def _tile_for_working_set(n_energies, n_inner, cell_entries, live_arrays=1,
                          max_entries=None):
    """Split an ``(n_energies, n_inner, ...)`` batch into tiles under a fixed entry budget.

    The batched engines build temporaries indexed by (energy, slab); their size is a
    product of two quantities the caller sets independently, so bounding either one alone
    does not bound the array.  This returns the tile to iterate in instead.

    Parameters
    ----------
    n_energies, n_inner : int
        Extent of the two batched axes.
    cell_entries : int
        Complex entries per (energy, inner) cell -- ``d*d`` for a stack of matrices,
        ``d*d*n_tpts_per_slab`` when each cell also carries quadrature samples.
    live_arrays : int, optional
        How many temporaries of this shape exist at once at the peak of the caller's
        loop.  The budget is divided by it, so that :data:`BATCH_WORKING_ENTRIES` bounds
        the engine's whole working set rather than one array of it -- the distinction is
        an eightfold one for the interaction-picture integrator, which holds the argument,
        the slab integral, ``Omega``, the slab operators and the matrix-exponential's own
        workspace simultaneously.  Default 1, which reproduces a budget stated per array.
    max_entries : int, optional
        Budget; defaults to :data:`BATCH_WORKING_ENTRIES`, read at call time rather than
        bound as a default argument so that the constant can be varied -- a test that
        monkeypatches a module attribute consumed as a default would pass trivially,
        comparing two identical runs.

    Returns
    -------
    (int, int)
        ``(energy_chunk, inner_block)``, both at least 1, whose product times
        ``cell_entries`` times ``live_arrays`` stays within the budget whenever that is
        possible at all.  A single cell larger than the whole budget cannot be split
        further and is returned as ``(1, 1)``: there is no tiling that helps, and refusing
        to proceed would be worse than a large allocation the caller can at least see.

    .. versionadded:: 1.0.0
    """
    if max_entries is None:
        max_entries = BATCH_WORKING_ENTRIES
    per_cell = max(1, int(cell_entries))*max(1, int(live_arrays))
    budget_cells = max(1, int(max_entries)//per_cell)
    n_energies = max(1, int(n_energies))
    n_inner = max(1, int(n_inner))
    inner_block = min(n_inner, max(1, budget_cells//n_energies))
    energy_chunk = min(n_energies, max(1, budget_cells//inner_block))
    return energy_chunk, inner_block


def _available_memory_bytes():
    """Best-effort free physical memory, or None if it cannot be had cheaply.

    ``MemAvailable`` is preferred where it exists because it accounts for reclaimable
    page cache, which the raw free-page count does not: on a machine with a warm cache
    the latter understates what a large allocation can actually get, and a guard built on
    it would refuse work that would have succeeded.

    Returns None rather than guessing on platforms that expose neither.  A guard that
    cannot measure must not block.

    .. versionadded:: 1.0.0
    """
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1])*1024
    except (OSError, ValueError, IndexError):
        pass
    try:
        return os.sysconf('SC_AVPHYS_PAGES')*os.sysconf('SC_PAGE_SIZE')
    except (AttributeError, ValueError, OSError):
        return None


def _check_output_fits(n_points, dim, source_func_name):
    """Refuse a scan whose *result* cannot fit in memory, before allocating anything.

    Tiling bounds the engines' working set, but nothing can shrink the answer itself: a
    scan of N points over d flavors returns ``N*d*d`` floats, and if that does not fit,
    no strategy helps.  Left unchecked, the failure arrives as an out-of-memory kill from
    somewhere deep in an engine -- or, on an overcommitting kernel, as the machine going
    down rather than the process.  Checking up front turns that into a message naming the
    number that is too large.

    Costs one multiply for ordinary calls: the free-memory figure is only consulted once
    the request passes :data:`OUTPUT_GUARD_MIN_BYTES`.

    Raises
    ------
    MemoryError
        If the result would claim more than 1/:data:`OUTPUT_GUARD_SAFETY` of available
        memory.  Never raised when free memory cannot be determined.

    .. versionadded:: 1.0.0
    """
    needed = int(n_points)*int(dim)*int(dim)*8          # float64 probability matrices
    if needed < OUTPUT_GUARD_MIN_BYTES:
        return
    available = _available_memory_bytes()
    if available is None:
        return
    if needed*OUTPUT_GUARD_SAFETY > available:
        raise MemoryError(
            gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": the requested "
            "scan would return " + f"{n_points:,}" + " probability matrices of size " +
            f"{dim}x{dim}" + " (" + f"{needed/2**30:.2f}" + " GiB), against " +
            f"{available/2**30:.2f}" + " GiB of available memory. The result alone would "
            "not fit, so no choice of method or tolerance can help. Split the scan into "
            "batches and concatenate the results.")


class ToleranceNotAchievedWarning(UserWarning):
    r"""Warns that a refinement ladder ran out of room before it converged.

    **What was detected.**  A cap was reached -- ``max_num_loops``, ``max_n_slabs``, or
    ``max_n_tpts_per_slab`` -- while the last two refinement levels still disagreed by more
    than the requested ``rtol``/``atol``.

    **What it means for the answer.**  The result is still exactly unitary, so it looks
    plausible; its accuracy is unverified.  The message says **how far from converged the ladder
    stopped**, as a multiple of the tolerance asked for, because that disagreement is computed
    anyway by the comparison that decides convergence and it is the difference between raising a
    cap as a guess and as a decision.  A ladder that stopped a few times outside tolerance is a
    different situation from one that stopped thirty times outside it.

    **What to change.**  Raise the cap the message names; or loosen ``rtol``/``atol`` to what
    the levels actually agreed to; or, if the profile has a density jump or a kink, pass
    ``t_breakpoints`` there -- no number of slabs fixes a slab that straddles one.

    **When it is safe to ignore.**  When the reported shortfall is smaller than the accuracy the
    result is used at.  Not otherwise: this reports a genuine failure of the convergence test,
    unlike :class:`magnus.magnus.MagnusConvergenceWarning`, which reports a slab width.

    Also raised by a **cumulative baseline scan** whose accuracy grid was sized from
    ``max_n_slabs`` rather than from a converged probe.  That instance is worth its own message
    because the consequence is larger than one point: the whole scan inherits the capped grid.

    Two subclasses narrow the diagnosis: :class:`HybridCertificationWarning` and
    :class:`UnmarkedDiscontinuityWarning`.  Code filtering on this class catches both.

    **Measured rates** (``docs/dev/adversarial_batteries/warn_fp.py``, 168 configurations):
    fired 37 times, **16 true positives and 21 false positives -- a 57 % false-positive rate**.
    A false positive here means the ladder genuinely ran out of room *and* the answer was
    nonetheless inside tolerance, which is the expected shape: a cap is reached before
    convergence has been *verified*, not before it has been *achieved*.

    .. versionadded:: 1.0.0
    """


class HybridCertificationWarning(ToleranceNotAchievedWarning):
    r"""Warns that ``strategy='hybrid'`` was requested (forcing
    :func:`magnus.adiabatic.hybrid_propagator`) but the self-certifying
    refinement of at least one requested (energy, L) point did not
    converge within its internal iteration/slab caps.

    The returned probabilities remain exactly unitary (every piece of the
    hybrid propagator -- the adiabatic transport and the local Magnus
    patches -- is unitary by construction) but their accuracy relative to
    the requested ``rtol``/``atol`` is not certified. This subclasses
    :class:`ToleranceNotAchievedWarning` so existing code that filters on
    the parent class also catches this warning; it is issued regardless
    of the verbosity setting. With the default ``strategy='auto'``, this
    situation instead falls back silently to the general slab-refinement
    method (which raises :class:`ToleranceNotAchievedWarning` itself if
    *it* also fails to converge), so this warning fires only when
    ``strategy='hybrid'`` was explicitly requested. See
    :doc:`/adiabatic_strategy`.

    **Uncertified means unverified, not wrong.**  The message says so, and names three things
    to change rather than leaving the reader with a disclaimer:

    * ``strategy='auto'`` -- the automatic fallback, which hands exactly these points to the
      general Magnus path and certifies there instead;
    * ``t_breakpoints`` at any known structure -- a density jump, a kink, or a feature narrower
      than 1/200 of the trajectory, which is the one cure for what the probe grid cannot
      resolve;
    * a looser ``rtol``/``atol``, when the accuracy needed is less than the accuracy requested.

    Which one applies is visible without guessing: pass ``strategy_info`` (see
    :func:`osc_prob_matter_std_potential`) and read ``'declined'``, or pass ``info`` to
    :func:`magnus.adiabatic.hybrid_propagator` directly for ``'resolved'`` and ``'gamma_max'``.

    .. versionadded:: 1.0.0
    """


class UnmarkedDiscontinuityWarning(ToleranceNotAchievedWarning):
    r"""Warns that a cumulative baseline scan was asked to integrate a Hamiltonian that is
    discontinuous at the scale of the grid it built, without being told where the
    discontinuities are.

    The cumulative scan lays a uniform accuracy grid over the trajectory (plus the requested
    baselines, plus any ``t_breakpoints``).  A slab that straddles a density jump degrades the
    quadrature to low order no matter how high ``magnus_exp_order`` is, and refining the grid
    does not fix it -- the straddling slab merely gets narrower.  The cure is to put an edge
    *on* the discontinuity, which is what ``t_breakpoints`` is for.

    Fuzzing 150 random piecewise-constant profiles, declaring the edges gave a median error of
    1.34e-12 and nothing outside tolerance; leaving them undeclared gave a median of 7.76e-04
    with 59 of 150 outside it.  Of those 59, all but **two** already warned for other reasons
    (usually :class:`ToleranceNotAchievedWarning` from the probe).  This warning exists for
    those two: measured at 1.36e-03 and 2.10e-03 against a requested 1e-3, silently, and
    2.33e-11 and 4.35e-14 once the edges were declared.

    Detection is a measurement, not a guess: the profile is sampled at two grid densities and
    the largest adjacent change in ``H`` is compared (see
    ``magnus.adiabatic._profile_is_resolved``).  A :math:`C^1` profile halves that change when
    the spacing halves; a jump does not.  On the profile families this package ships --- solar
    exponential, multi-resonance, noisy, sinusoidal --- the test reports "resolved" every time,
    and it flagged 12 of 12 random piecewise profiles.

    Subclasses :class:`ToleranceNotAchievedWarning` so that code already filtering on the parent
    also catches this.  Not raised when ``t_breakpoints`` was supplied: the caller has then said
    where the edges are, and the grid honours them.

    **Measured rates** (``docs/dev/adversarial_batteries/warn_fp.py``, 168 configurations
    including 48 random piecewise-constant profiles with the edges deliberately left
    undeclared): fired 56 times, **23 true positives and 33 false positives -- 59 %**.  Read
    that number carefully: this reports a *condition about the input*, not a prediction about
    the error, and on every one of those 33 the condition was real -- there was an undeclared
    discontinuity -- and the answer happened to come out inside tolerance anyway.  Declaring the
    edges would still have improved it (median 7.8e-04 to 1.3e-12 in ``FINDINGS`` §9.2).  A
    warning whose claim is true and whose advice is worth taking is not made a false alarm by
    the answer surviving.

    **Also raised on the hybrid path.**  The same detector already ran inside
    :func:`magnus.adiabatic.hybrid_propagator`, where failing it makes the strategy decline --
    silently, so the caller heard about slab widths from whichever engine answered instead,
    which is true and points at the wrong knob.  It now says what it found there too.  The
    detector, its two-stage protocol and its measured false-positive rate are unchanged; only
    the number of places that report it has grown.  On an unmarked density step the adiabatic
    answer was wrong by **0.54** in probability while reporting itself certified, and that is
    the case this instance exists for.

    .. versionadded:: 1.0.0
    """


def _shortfall_phrase(last_gap, rtol, atol) -> str:
    r"""How far from converged a refinement ladder stopped, as a fixed phrase.

    :class:`ToleranceNotAchievedWarning` used to say only *that* the ladder ran out of room.
    The disagreement between the last two levels is computed anyway, by the very comparison that
    decides convergence, so the warning can say by how much -- which is the difference between
    "raise the cap" as a guess and as a decision.

    Bucketed rather than numeric so the message stays one of four fixed strings and Python's
    default filter still shows each at most once per session.  The buckets are ratios to the
    requested tolerance, because that is the quantity the caller chose and can change.

    .. versionadded:: 1.0.0
    """
    tol = (atol or 0.0) + (rtol or 0.0)
    if (last_gap is None) or (tol <= 0.0):
        return "with no two levels to compare"
    ratio = last_gap/tol
    if ratio <= 3.0:
        return "the last two refinement levels still differing by a few times the tolerance"
    if ratio <= 30.0:
        return ("the last two refinement levels still differing by roughly ten times the "
                "tolerance")
    return ("the last two refinement levels still differing by more than thirty times the "
            "tolerance")


class HiddenFeatureWarning(ToleranceNotAchievedWarning):
    r"""Warns that the profile has structure too narrow for **any** grid this package lays down.

    **What was detected.**  A feature whose variation is concentrated between samples of even the
    finest grid the adaptive machinery reaches -- see
    :func:`magnus.adiabatic.find_hidden_features` and
    :data:`magnus.adiabatic.HIDDEN_FEATURE_CONCENTRATION`.  The message names the position.

    **What it means for the answer.**  Possibly wrong, and *no choice of strategy or tolerance
    helps*.  This is the one exposure the adversarial validation
    (``docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md`` §8.3) could not close: the hybrid
    strategy's probe grid, the general ladder's slabs and the cumulative scan's accuracy grid all
    miss the same feature, so they agree with each other and are wrong together -- which is also
    why :func:`cross_check_strategies` cannot see it either.  Measured on a Gaussian of width
    :math:`3\times10^{-5}` of the trajectory: **wrong by 2.9e-02 against a requested 1e-3, with
    no warning at all** before this existed.

    **What to change.**  Pass ``t_breakpoints`` bracketing the position in the message, padded by
    a few reference intervals.  Measured on the two constructions above, going from nothing to a
    padded breakpoint set: 3.9e-03 → 8.5e-05 and 1.3e-03 → 4.4e-04, and the answer stops being
    silent.  The cure is real but **partial** -- putting edges on a feature helps the quadrature,
    it does not conjure resolution that the sampling never had.  For a feature you know the width
    of, a denser grid there is better still.

    **When it is safe to ignore.**  When the narrow structure is an artefact of how the profile
    function was written rather than physics -- an interpolation kink, a rounding step in a
    tabulated density -- and you know the physical profile is smooth there.

    Subclasses :class:`ToleranceNotAchievedWarning`, so code already filtering on the parent
    catches it.  Not raised when ``t_breakpoints`` was supplied: the caller has then already said
    where the structure is.  The scan depends on the profile and the interval but not on energy,
    so it runs **once per call**, not once per (energy, L) point.

    .. versionadded:: 1.0.0
    """


class PhaseAveragingWarning(UserWarning):
    r"""Warns that ``average=True`` was requested at an (energy, L) point
    where the oscillation has not, in fact, averaged.

    The phase-averaged probability is the exact limit reached when every
    pair of eigenvalues has accumulated many cycles of relative phase (see
    :mod:`magnus.avgprob`).  A pair whose relative phase is neither much
    larger than :math:`2\pi` nor much smaller than one radian is in
    neither limit, and no averaged expression describes it -- the
    oscillation probability itself is the meaningful quantity there.

    This is not a statement about numerical accuracy: the returned matrix
    is still a valid, doubly stochastic probability matrix.  It is a
    statement that the *question* does not apply at that baseline, which
    is why it warns rather than refining anything.

    .. versionadded:: 1.0.0
    """


#-----------------------------------------------------------------------
# Helper functions
#-----------------------------------------------------------------------

def print_banner(file: TextIOWrapper=None):
    r"""Prints the Magnus ASCII banner, version, and author string.

    Prints an ASCII-art banner followed by the package version (``magnus.__version__``, resolved
    from ``pyproject.toml``) and author (``magnus.authors.__authors__``).  Both live in internal
    metadata modules that are excluded from the API reference, so they are shown as literals
    rather than as cross-references.  Uses ANSI color codes when printing to
    stdout (``file is None``); plain text otherwise (e.g., when writing to a log file).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    file : TextIOWrapper, optional
        File object to print to, in addition to (or instead of, depending on the caller) stdout.
        If None (default), print to stdout with color.

    Returns
    -------
    None
    """
    if file is None:
        print(gd.cstyle.CBLUEBG + ".----------------------------------------." + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|   __  __                               |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  |  \/  | __ _  __ _ _ __  _   _ ___   |" + gd.cstyle.CEND, 
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  | |  | | (_| | (_| | | | | |_| \__ \  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|                |___/                   |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "'----------------------------------------'" + gd.cstyle.CEND,
            file=file)
    else: 
        print(".----------------------------------------.", file=file)
        print("|   __  __                               |", file=file)
        print(r"|  |  \/  | __ _  __ _ _ __  _   _ ___   |", file=file)
        print(r"|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |", file=file)
        print(r"|  | |  | | (_| | (_| | | | | |_| \__ \  |", file=file)
        print(r"|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |", file=file)
        print("|                |___/                   |", file=file)
        print("'----------------------------------------'", file=file)
    print("Version: "+ version.__version__ + " | Author(s): " + authors.__authors__ + "\n",
        file=file)


def print_run_parameters(
    H_func: Union[Callable, np.ndarray], 
    t_ini: Union[int, float], 
    t_fin: Union[int, float], 
    n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='gl', 
    rtol: Optional[Union[int, float]]=None, 
    atol: Optional[Union[int, float]]=None, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=None, 
    min_n_tpts_per_slab: Optional[int]=2, 
    max_n_tpts_per_slab: Optional[int]=500,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log', 
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0, 
    file_log: Optional[TextIOWrapper]=None
):
    r"""Prints the banner (once per session) and the parameters passed to :func:`osc_prob`.

    Diagnostic/logging helper called from :func:`osc_prob` when ``verbose >= 1`` or ``save_log``
    is True.  Prints (to stdout, and additionally to ``file_log`` if ``save_log`` is True) the
    values of every refinement/logging parameter for the current call, to help reproduce or debug
    a specific run.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable or np.ndarray
        The Hamiltonian passed to :func:`osc_prob`.
    t_ini, t_fin : int or float
        Integration limits passed to :func:`osc_prob`.
    n_slabs : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    n_tpts_per_slab : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    t_slab_edges : list or np.ndarray, optional
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    magnus_exp_order : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    n_jobs : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    integration_method : str
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    rtol : int or float, optional
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    atol : int or float, optional
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    growth_factor_n_slabs : int or float
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    growth_factor_n_tpts_per_slab : int or float
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    max_num_loops : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    min_n_slabs : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    max_n_slabs : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    min_n_tpts_per_slab : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    max_n_tpts_per_slab : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    validate_input : bool
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    save_log : bool
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    filename_log : str
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    new_recursion_limit : int, optional
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    verbose : int
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.
    file_log : TextIOWrapper, optional
        Forwarded verbatim from the calling :func:`osc_prob`; see its docstring.

    Returns
    -------
    None
    """
    global has_magnus_header_been_printed

    for f in [None, file_log] if save_log else [None]:
        if not has_magnus_header_been_printed:
            print_banner(f)
            has_magnus_header_been_printed = True
        print("Parameters passed to function magnus.osc_prob in this run:", file=f)
        if callable(H_func):
            print("   H_func = " + H_func.__name__, file=f)
        else:
            print("   H_func = constant (time-independent)", file=f)
        print("   t_ini = " + str(t_ini), file=f)
        print("   t_fin = " + str(t_fin), file=f)
        print("   n_slabs = " + str(n_slabs), file=f)
        print("   n_tpts_per_slab = " + str(n_tpts_per_slab), file=f)
        if t_slab_edges is None:
            print("   t_slab_edges = None", file=f)
        else:
            print("   t_slab_edges = ", file=f)
            for i, t_slab in enumerate(t_slab_edges):
                print("      i" + ": " + str(t_slab), file=f)
        print("   magnus_exp_order = " + str(magnus_exp_order), file=f)
        print("   n_jobs = " + str(n_jobs), file=f)
        print("   integration_method = " + integration_method, file=f)
        print("   rtol = " + str(rtol), file=f)
        print("   atol = " + str(atol), file=f)
        print("   growth_factor_n_slabs = " + str(growth_factor_n_slabs), file=f)
        print("   growth_factor_n_tpts_per_slab = " + str(growth_factor_n_tpts_per_slab), file=f)
        print("   max_num_loops = " + str(max_num_loops), file=f)
        print("   min_n_slabs = " + str(min_n_slabs), file=f)
        print("   max_n_slabs = " + str(max_n_slabs), file=f)
        print("   min_n_tpts_per_slab = " + str(min_n_tpts_per_slab), file=f)
        print("   max_n_tpts_per_slab = " + str(max_n_tpts_per_slab), file=f)
        print("   validate_input = " + str(validate_input), file=f)
        print("   save_log = " + str(save_log), file=f)
        print("   filename_log = " + filename_log, file=f)
        print("   new_recursion_limit = " + str(new_recursion_limit), file=f)
        print("   verbose = " + str(verbose), file=f)

    return


def validate_input_battery(
    source_func_name: str, 
    energy: Optional[Union[int, float, list, np.ndarray]]=None, 
    L: Optional[Union[int, float, list, np.ndarray]]=None, 
    L0: Optional[Union[int, float]]=None,
    num_flavors: Optional[int]=None,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    osc_params: Optional[Union[list, np.ndarray]]=None,
    rho_func: Optional[Union[Callable, int, float]]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    validate_energy_and_L: Optional[bool]=True,
    validate_flavor_indices: Optional[bool]=True,
    validate_osc_params: Optional[bool]=True,
    validate_initial_position: Optional[bool]= False,
    validate_density: Optional[bool]=False
) -> None:
    r"""Validates the inputs common to the ``osc_prob_*`` family of functions.

    Runs a battery of type/shape/value checks (selected by the ``validate_*`` flags below) and
    raises :class:`ValueError` with a descriptive message identifying the offending argument and
    the calling function (via ``source_func_name``) if any check fails, rather than letting an
    invalid input propagate into a cryptic NumPy/linear-algebra error deep inside the Magnus
    core.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, used to build more informative error messages.
    energy : int, float, list, or np.ndarray, optional
        Neutrino energy/energies to validate (checked if ``validate_energy_and_L`` is True).
    L : int, float, list, or np.ndarray, optional
        Baseline(s) to validate (checked if ``validate_energy_and_L`` is True).
    L0 : int or float, optional
        Initial position to validate (checked if ``validate_initial_position`` is True).
    num_flavors : int, optional
        Number of neutrino flavors, used to validate ``nu_i``/``nu_f``/``osc_params``.
    nu_i : int, optional
        Initial flavor index to validate (checked if ``validate_flavor_indices`` is True).
    nu_f : int, optional
        Final flavor index to validate (checked if ``validate_flavor_indices`` is True).
    osc_params : list or np.ndarray, optional
        Unpacked oscillation parameters to validate (checked if ``validate_osc_params`` is True).
    rho_func : Callable, int, or float, optional
        Matter density (function or constant) to validate (checked if ``validate_density`` is
        True).
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter, validated alongside the density.
        Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction, validated alongside the density. Default: 0.5.
    validate_energy_and_L : bool, optional
        If True, validate ``energy`` and ``L``. Default: True.
    validate_flavor_indices : bool, optional
        If True, validate ``nu_i`` and ``nu_f`` against ``num_flavors``. Default: True.
    validate_osc_params : bool, optional
        If True, validate ``osc_params``. Default: True.
    validate_initial_position : bool, optional
        If True, validate ``L0``. Default: False.
    validate_density : bool, optional
        If True, validate ``rho_func``, ``ratio_number_neutrons_to_protons``, and
        ``electron_fraction``. Default: False.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If any requested check fails.  The message names the offending argument and the calling
        function.
    """
    if validate_energy_and_L:

        if ( (not isinstance(energy, int)) and (not isinstance(energy, float)) and \
            (not isinstance(energy, list)) and (not isinstance(energy, np.ndarray)) ):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": energy must be an int, a float, a 1D list, or a 1D NumPy array.")

        if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) and \
            (np.array(energy).ndim != 1) ):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": if energy is a list or NumPy array, it must be 1D.")

        if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) ):
            # (np.issubdtype is used instead of the np.float_/np.int_ aliases, which were
            # removed in NumPy 2.0)
            if not (np.issubdtype(np.asarray(energy).dtype, np.floating) or \
                np.issubdtype(np.asarray(energy).dtype, np.integer)):
                raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                    ": since energy is a list or NumPy array, all of its elements must be int" + \
                    " or float.")

        if ( (not isinstance(L, int)) and (not isinstance(L, float)) and \
            (not isinstance(L, list)) and (not isinstance(L, np.ndarray)) ):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": L must be an int, a float, a 1D list, or a 1D NumPy array.")

        if ( (isinstance(L, list) or isinstance(L, np.ndarray)) and \
            (np.array(L).ndim != 1) ):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": if L is a list or NumPy array, it must be 1D.")

        if ( (isinstance(L, list) or isinstance(L, np.ndarray)) ):
            if not (np.issubdtype(np.asarray(L).dtype, np.floating) or \
                np.issubdtype(np.asarray(L).dtype, np.integer)):
                raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                    ": since L is a list or NumPy array, all of its elements must be int or float.")

        if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) and \
            (isinstance(L, list) or isinstance(L, np.ndarray)) and \
            (len(energy) != len(L)) ):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": since the input energy and L are both lists or NumPy arrays, they must have " + \
                "the same length.")

        if (((nu_i is not None) and (nu_f is None)) or ((nu_i is None) and (nu_f is not None))):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": if either nu_i or nu_f is not None, the other flavor must also be not None.")

    if validate_flavor_indices:

        if ((nu_i is not None) and (nu_f is not None)):
            if (num_flavors <= gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
                if ((num_flavors == 2) or (num_flavors == 3)):
                    flavors = set([gd.NUE, gd.NUMU, gd.NUTAU])
                elif (num_flavors == 4):
                    flavors = set([gd.NUE, gd.NUMU, gd.NUTAU, gd.NUS])
                elif (num_flavors == 5):
                    flavors = set([gd.NUE, gd.NUMU, gd.NUTAU, gd.NUS1, gd.NUS2])
                if ((nu_i not in flavors) or (nu_f not in flavors)):
                    if ((num_flavors == 2) or (num_flavors == 3)):
                        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                            ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                            str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), or gd.NUTAU (" + \
                            str(gd.NUTAU) + ") only.")
                    elif (num_flavors == 4):
                        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                            ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                            str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), gd.NUTAU (" + \
                            str(gd.NUTAU) + "), or gd.NUS (" + str(gd.NUS) + ") only.")
                    elif (num_flavors == 5):
                        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                            ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                            str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), gd.NUTAU (" + \
                            str(gd.NUTAU) + "), gd.NUS1 (" + str(gd.NUS1) + "), or gd.NUS2 (" + \
                            str(gd.NUS2) + ") only.")
            else:
                print(gd.WARNING_MSG_IN_COLOR + " " + source_func_name + \
                    ": nu_i and nu_f are not None, but, since num_flavors = " + str(num_flavors) + \
                    " > globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
                    str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + ", input validation cannot " + \
                    "check if nu_e and nu_f are valid indices.")

    if validate_osc_params:

        ttest = [(isinstance(x, int) or isinstance(x, float) or (x is None)) 
            for x in osc_params]
        if (not np.all(ttest)):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ":"+\
                " the oscillation parameters must be int or float.")

    if validate_initial_position:

        if not ((isinstance(L0, int) or (isinstance(L0, float)))):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                " the initial neutrino position (L0) must be an int or float.")

    if validate_density:

        if (ratio_number_neutrons_to_protons < 0.0):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                " the ratio of neutrinos to protons (ratio_number_neutrons_to_protons) must" + \
                " be non-negative.")

        if (electron_fraction < 0.0):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                " the ratio of electrons to protons + neutrons (electron_fraction) must be " + \
                "non-negative.")

        if ((callable(rho_func)) and (_n_required_params(rho_func) > 1)):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                " the provided rho_func is a function of more than one parameter.")

        rho_test = rho_func(L0) if callable(rho_func) else rho_func

        if (rho_test < 0.0):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                " rho_func must be non-negative.")

        if not (isinstance(rho_test, int) or isinstance(rho_test, float)):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                " rho_func must be a float (or int) or must return a float (or int).")


def validate_input_osc_prob_earth(
    source_func_name: str,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    costhz: Optional[Union[int, float]]=None,
    L: Optional[Union[float, list, np.ndarray]]=None,
    verbose: Optional[int]=0,
    ) -> Tuple[float, np.ndarray]:
    r"""Resolves (costhz, L) for :func:`osc_prob_earth`, from either two locations or costhz+L.

    Implements the two mutually exclusive ways of specifying an Earth-crossing trajectory: either
    give both ``loc_ini`` and ``loc_fin`` (the chord's zenith angle and length are computed from
    their coordinates), or give ``costhz`` and ``L`` directly. Aborts with a descriptive error if
    exactly one location is given, or if neither locations nor (costhz, L) are given.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, used to build more informative error messages.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location, as (latitude, longitude) coordinates or a predefined location name (see
        :data:`magnus.earth.loc_coords_dms`). Must be given together with ``loc_fin``.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given together with ``loc_ini``.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used only if ``loc_ini``/``loc_fin`` are not
        given.
    L : float, list, or np.ndarray, optional
        Baseline(s). Used only if ``loc_ini``/``loc_fin`` are not given.
    verbose : int, optional
        Verbosity level: if > 0, print a note when the chord between the two given locations is
        used as the baseline. Default: 0.

    Returns
    -------
    (float, np.ndarray)
        The resolved ``(costhz, L)`` pair.
    """
    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given.
    # If only a single location is given, throw an exception.  If neither of the two locations are
    # given, use the given value of costhz and of baseline given (could be an array of baselines).
    
    if ( ((loc_ini is None) and (loc_fin is not None)) or \
        ((loc_ini is not None) and (loc_fin is None)) ):

        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": only one of the two " + \
            "locations on Earth (loc_ini or loc_fin) has been given. If one location is " + \
            "given (i.e., is not None), the other one must also be given.  Alternatively, " + \
            "both locations can be set to None, and the given value of costhz will be used " +\
            "(if it is not None).")

    elif ((loc_ini is not None) and (loc_fin is not None)):

        # Check that the location is a two-entry tuple, list, or array

        # Unpacking a sequence of the wrong length raises ValueError, and unpacking
        # something that is not iterable at all raises TypeError; neither is a KeyError,
        # which is what this used to catch, so the message below could never be reached.
        # A three-entry tuple reported "too many values to unpack (expected 2)", and an
        # int escaped as a TypeError -- breaking the convention that bad input to this
        # package raises ValueError with a message naming the parameter at fault.
        try:
            lat_ini, lon_ini = loc_ini
        except (TypeError, ValueError):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": if the initial " + \
                    "location (loc_ini) is given as coordinates, it must be a two-entry tuple," + \
                    " list, or NumPy array.")

        try:
            lat_fin, lon_fin = loc_fin
        except (TypeError, ValueError):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": if the final " + \
                    "location (loc_fin) is given as coordinates, it must be a two-entry tuple," + \
                    " list, or NumPy array.")

        # We use the function earth.costhz_between_points_on_surface to compute the cosine of the
        # zenith angle of the chord that joins two locations on the surface of the Earth, measured 
        # at one position (any of the two locations will give the same result).
        costhz = earth.costhz_between_points_on_surface(lat_ini, lon_ini, lat_fin, lon_fin)

        # Length of the chord is the baseline
        L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM # [eV^{-1}]

        if verbose > 0:
            print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": using as " + \
                "baseline the chord that joins the given initial and final locations on the " + \
                "surface of the Earth.")

        return costhz, L

    else: # (loc_ini is None) and (loc_fin is None)

        if costhz is None:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": no" + \
                " initial and final locations on the surface of the Earth given, and no " + \
                "value of costhz given.  This function requires either the two locations " + \
                "or, alternatively, the value of costhz.")

        if L is None:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": since two locations on the surface of the Earth have not been given, " + \
                "the value of costhz will be used to define the chord length, but the" + \
                " baseline, L, cannot be None.")

        return costhz, L


def valid_flavor_indices_2nu(nu_i: int, nu_f: int) -> Tuple[int, int]:
    r"""Remaps 3-flavor-style flavor indices onto valid 2-flavor indices (0 or 1).

    Two-flavor wrappers (e.g. :func:`osc_prob_2nu_matter_constant_density`) accept ``nu_i``/
    ``nu_f`` values from the same ``NUE``/``NUMU``/``NUTAU`` constants used by the 3/4/5-flavor
    wrappers, for interface consistency, even though a two-flavor system only has indices 0 and 1.
    This remaps the flavor not included in the two-flavor system (whichever of NUE/NUMU/NUTAU is
    not being used) onto the other valid index, so that, e.g., requesting the nu_e-nu_tau channel
    of a system parametrized by :math:`\theta_{13}` (which is really a nu_e-nu_x system) resolves correctly.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    nu_i : int
        Initial flavor index, as one of ``globaldefs.NUE``, ``NUMU``, ``NUTAU``, or already 0/1.
    nu_f : int
        Final flavor index, same convention as ``nu_i``.

    Returns
    -------
    (int, int)
        The remapped ``(nu_i, nu_f)``, each 0 or 1.
    """
    if ((nu_i == gd.NUE) and (nu_f == gd.NUTAU)):
        nu_f = 1
    elif ((nu_i == gd.NUTAU) and (nu_f == gd.NUE)):
        nu_i = 1
    elif ((nu_i == gd.NUMU) and (nu_f == gd.NUTAU)):
        nu_i, nu_f = 0, 1
    elif ((nu_i == gd.NUTAU) and (nu_f == gd.NUMU)):
        nu_i, nu_f = 1,0

    return nu_i, nu_f


def values_to_unspecified_osc_params(
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    verbose: Optional[int]=0
) -> Tuple[float, float, float, float, float, float]:
    r"""Return values of unspecified standard oscillation parameters

    If any of the oscillation parameters has not been given a value, assign to it the value from
    the specified parameter set with name default_osc_params_set_name.  When input validation is
    on (validate_input == True), the routine checks whether the parameter set name is among the
    predefined ones (see validation above).  Only the values of the parameters passed as None are
    assigned from the predefined set; other parameters are not modified.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`; if None, taken from the predefined set.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`; if None, taken from the predefined set.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`; if None, taken from the predefined set.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]; if None, taken from the predefined set.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`; if None, taken from the predefined set.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`; if None, taken from the predefined set.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set to draw missing values from (see
        ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    verbose : int, optional
        Verbosity level. Default: 0.

    Returns
    -------
    (float, float, float, float, float, float)
        ``(s12, s23, s13, dCP, D21, D31)``, with every previously-None entry filled in from the
        predefined set.
    """

    if default_osc_params_set_name not in list(gd.OSC_PARAMS_PREDEFINED.keys()):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.values_to_unspecified_osc_params"+ \
            ": the requested oscillation parameter set (default_osc_params_set_name = " + \
            default_osc_params_set_name + ") is not among the predefined sets in Magnus. " + \
            "Available sets are " + str(list(gd.OSC_PARAMS_PREDEFINED.keys())) + ".")

    global has_magnus_header_been_printed

    if ((s12 is None) or (s23 is None) or (s13 is None) or (s23 is None) or (dCP is None) or \
        (D21 is None) or (D31 is None)):

        default_osc_params = gd.OSC_PARAMS_PREDEFINED[default_osc_params_set_name]

        if verbose > 0:
            if verbose >= 2:
                if (not has_magnus_header_been_printed):
                    print_banner()
                    has_magnus_header_been_printed = True
            print(gd.WARNING_MSG_IN_COLOR + " Setting unspecified standard oscillation " + \
                "parameters to default values from the predefined set " + \
                default_osc_params['name'] + " (" + default_osc_params['description'] + "):\n" + \
                ("s12 = " + str(default_osc_params['s12']) + "\n" if (s12 is None) else '') + \
                ("s23 = " + str(default_osc_params['s23']) + "\n" if (s23 is None) else '') + \
                ("s13 = " + str(default_osc_params['s13']) + "\n" if (s13 is None) else '') + \
                ("dCP = " + str(default_osc_params['dCP']) + " rad\n" if (dCP is None) else '') + \
                ("D21 = " + str(default_osc_params['D21']) + " eV^2\n" if (D21 is None) else '') + \
                ("D31 = " + str(default_osc_params['D31']) + " eV^2\n" if (D31 is None) else ''))

        s12 = s12 if (s12 is not None) else default_osc_params['s12']
        s23 = s23 if (s23 is not None) else default_osc_params['s23']
        s13 = s13 if (s13 is not None) else default_osc_params['s13']
        dCP = dCP if (dCP is not None) else default_osc_params['dCP']
        D21 = D21 if (D21 is not None) else default_osc_params['D21']
        D31 = D31 if (D31 is not None) else default_osc_params['D31'] 

    return s12, s23, s13, dCP, D21, D31


def unpack_oscillation_params_from_dict(
    source_func_name: str,
    num_flavors: int,
    osc_params: Dict,
    h_vac_energy_indep: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Unpack oscillation parameters from the osc_params dict

    Extracts the standard oscillation parameters for ``num_flavors`` flavors from ``osc_params``
    (as built by each ``osc_prob_{N}nu_*`` wrapper), in the fixed order expected by the matching
    ``hamiltonians.hamiltonian_{N}nu_vacuum_energy_independent`` function. Aborts with a
    descriptive error if a required key is missing.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, used to build more informative error messages.
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_vac_energy_indep`` is given).
    osc_params : dict
        Dictionary of oscillation parameters. For ``num_flavors == 2``, must contain 'sth', 'Dm2'.
        For 3, 4, 5, must contain 's12', 's23', 's13', 'dCP', 'D21', 'D31', plus, for 4:
        's14', 'd14', 's24', 'd24', 's34', 'D41'; and for 5, additionally 's15', 'd15', 's25',
        's34', 's35', 'd35', 'D51'.
    h_vac_energy_indep : list or np.ndarray
        Precomputed energy-independent vacuum Hamiltonian, required (and used, instead of
        ``osc_params``) when ``num_flavors`` exceeds
        ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.

    Returns
    -------
    np.ndarray
        The unpacked oscillation parameters, in the order expected by the matching
        ``hamiltonian_{N}nu_vacuum_energy_independent`` function.
    """

    if (num_flavors == 2):
        try:
            sth = osc_params['sth']
            Dm2 = osc_params['Dm2']
            return np.array([sth, Dm2])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since "+ \
                    "num_flavors == 2, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 'sth' and 'Dm2'.")
    elif (num_flavors == 3):
        try:
            s12 = osc_params['s12']
            s23 = osc_params['s23']
            s13 = osc_params['s13']
            dCP = osc_params['dCP']
            D21 = osc_params['D21']
            D31 = osc_params['D31']
            return np.array([s12, s23, s13, dCP, D21, D31])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 3, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 's12', 's23', 's13', 'dCP', 'D21', and " + \
                    "'D31', even if they are None.")
    elif (num_flavors == 4):
        try:
            s12 = osc_params['s12']
            s23 = osc_params['s23']
            s13 = osc_params['s13']
            dCP = osc_params['dCP']
            s14 = osc_params['s14']
            d14 = osc_params['d14']
            s24 = osc_params['s24']
            d24 = osc_params['d24']
            s34 = osc_params['s34']
            D21 = osc_params['D21']
            D31 = osc_params['D31']
            D41 = osc_params['D41']
            return np.array([s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 4, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 's12', 's23', 's13', 'dCP', 'D21', and " + \
                    "'D31' (even if they are None); and 's14', 'd14', 's24', 'd24', 's34', and " + \
                    "'D41'.")
    elif (num_flavors == 5):
        try:
            s12 = osc_params['s12']
            s23 = osc_params['s23']
            s13 = osc_params['s13']
            dCP = osc_params['dCP']
            s14 = osc_params['s14']
            d14 = osc_params['d14']
            s15 = osc_params['s15']
            d15 = osc_params['d15']
            s24 = osc_params['s24']
            d24 = osc_params['d24']
            s25 = osc_params['s25']
            s34 = osc_params['s34']
            s35 = osc_params['s35']
            d35 = osc_params['d35']
            D21 = osc_params['D21']
            D31 = osc_params['D31']
            D41 = osc_params['D41']
            D51 = osc_params['D51']
            return np.array([s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, \
                D21, D31, D41, D51])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 5, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 's12', 's23', 's13', 'dCP', 'D21', and " + \
                    "'D31' (even if they are None); and 's14', 'd14', 's15', 'd15', 's24', " + \
                    "'d24', 's25', 's34', 's35', 'd35', 'D41', and 'D51'.")
    elif (num_flavors > gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
        print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": the number of " + \
            "flavors passed (num_flavors = " + str(num_flavors) + \
            ") exceeds the maximum number for which Magnus has predefined vacuum Hamiltonians " + \
            "(globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
            str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + "). Will use the Hamiltonian provided " + \
            "in h_vac_energy_indep.")
        if (h_vac_energy_indep is None):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": provided " + \
                "h_vac_energy_indep is None.")
    elif (num_flavors < 1):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": num_flavors must be " + \
            ">= 2.")


def unpack_nsi_params_from_dict(
    source_func_name: str,
    num_flavors: int,
    nsi_params: Dict,
    h_nsi: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Unpack NSI parameters from the nsi_params dict

    Extracts the NSI epsilon parameters for ``num_flavors`` flavors from ``nsi_params`` (as built
    by each ``osc_prob_{N}nu_*_nsi_*`` wrapper), in the fixed order expected by the matching
    ``hamiltonians.hamiltonian_{N}nu_nsi`` function. Aborts with a descriptive error if a required
    key is missing.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, used to build more informative error messages.
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_nsi`` is given).
    nsi_params : dict
        Dictionary of NSI parameters. For ``num_flavors == 2``, must contain 'eps_aa', 'eps_ab'.
        For 3, must contain 'eps_ee', 'eps_em', 'eps_et', 'eps_mm', 'eps_mt', 'eps_tt'. For 4,
        additionally 'eps_es', 'eps_ms', 'eps_ts', 'eps_ss'. For 5, instead of the sterile-flavor
        keys above, 'eps_es1', 'eps_es2', 'eps_ms1', 'eps_ms2', 'eps_ts1', 'eps_ts2', 'eps_s1s1',
        'eps_s1s2', 'eps_s2s2'.
    h_nsi : list or np.ndarray
        Precomputed NSI Hamiltonian, required when ``num_flavors`` exceeds
        ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS`` (in which case there is nothing to
        unpack from ``nsi_params`` and this function returns None).

    Returns
    -------
    np.ndarray or None
        The unpacked NSI parameters, in the order expected by the matching
        ``hamiltonian_{N}nu_nsi`` function; or None if ``num_flavors`` exceeds
        ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS`` (the caller uses ``h_nsi`` directly).
    """

    if (num_flavors == 2):
        try:
            eps_aa = nsi_params['eps_aa']
            eps_ab = nsi_params['eps_ab']
            return np.array([eps_aa, eps_ab])
        except KeyError :
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since "+ \
                    "num_flavors == 2, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_aa' and 'eps_ab'.")
    elif (num_flavors == 3):
        try:
            eps_ee = nsi_params['eps_ee']
            eps_em = nsi_params['eps_em']
            eps_et = nsi_params['eps_et']
            eps_mm = nsi_params['eps_mm']
            eps_mt = nsi_params['eps_mt']
            eps_tt = nsi_params['eps_tt']
            return np.array([eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 3, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_ee', 'eps_em', 'eps_et', 'eps_mm'," + \
                    " 'eps_mt', and 'eps_tt'.")
    elif (num_flavors == 4):
        try:
            eps_ee = nsi_params['eps_ee']
            eps_em = nsi_params['eps_em']
            eps_et = nsi_params['eps_et']
            eps_es = nsi_params['eps_es']
            eps_mm = nsi_params['eps_mm']
            eps_mt = nsi_params['eps_mt']
            eps_ms = nsi_params['eps_ms']
            eps_tt = nsi_params['eps_tt']
            eps_ts = nsi_params['eps_ts']
            eps_ss = nsi_params['eps_ss']
            return np.array([eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts,
                eps_ss])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 4, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_ee', 'eps_em', 'eps_et', 'eps_es'," + \
                    " 'eps_mm', 'eps_mt', 'eps_ms', 'eps_tt', 'eps_ts', and 'eps_ss'.")
    elif (num_flavors == 5):
        try:
            eps_ee = nsi_params['eps_ee']
            eps_em = nsi_params['eps_em']
            eps_et = nsi_params['eps_et']
            eps_es1 = nsi_params['eps_es1']
            eps_es2 = nsi_params['eps_es2']
            eps_mm = nsi_params['eps_mm']
            eps_mt = nsi_params['eps_mt']
            eps_ms1 = nsi_params['eps_ms1']
            eps_ms2 = nsi_params['eps_ms2']
            eps_tt = nsi_params['eps_tt']
            eps_ts1 = nsi_params['eps_ts1']
            eps_ts2 = nsi_params['eps_ts2']
            eps_s1s1 = nsi_params['eps_s1s1']
            eps_s1s2 = nsi_params['eps_s1s2']
            eps_s2s2 = nsi_params['eps_s2s2']
            return np.array([eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1,
                eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 5, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_ee', 'eps_em', 'eps_et', " + \
                    "'eps_es1', 'eps_es2', 'eps_mm', 'eps_mt', 'eps_ms1', 'eps_ms2', 'eps_tt', " + \
                    "'eps_ts1', 'eps_ts2', 'eps_s1s1', 'eps_s1s2', and 'eps_s2s2'.")
    elif (num_flavors > gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
        print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": the number of " + \
            "flavors passed (num_flavors = " + str(num_flavors) + \
            ") exceeds the maximum number for which Magnus has predefined vacuum Hamiltonians " + \
            "(globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
            str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + "). Will use the Hamiltonian provided " + \
            "in h_nsi.")
        if (h_nsi is None):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": provided " + \
                "h_nsi is None.")
        # num_flavors exceeds the predefined range: the caller builds its Hamiltonian directly from
        # h_nsi instead of from a flat parameter list, so there is nothing to unpack here.
        return None
    elif (num_flavors < 1):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": num_flavors must be " + \
            ">= 2.")


def unpack_liv_params_from_dict(
    source_func_name: str,
    num_flavors: int,
    liv_params: Dict,
    h_liv: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Unpack LIV parameters from the liv_params dict

    Extracts the LIV parameters for ``num_flavors`` flavors from ``liv_params`` (as built by each
    ``osc_prob_{N}nu_*_liv`` wrapper), in the fixed order expected by the matching
    ``hamiltonians.hamiltonian_{N}nu_liv_energy_independent`` function. Validates that ``Lambda``
    is positive and aborts with a descriptive error if a required key is missing.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, used to build more informative error messages.
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_liv`` is given).
    liv_params : dict
        Dictionary of LIV parameters. Always must contain 'Lambda' (LIV energy scale, must be
        positive) and 'n_liv' (power of the energy dependence). For ``num_flavors == 2``, must
        also contain 'sxi', 'b1', 'b2'. For 3, 'sxi12', 'sxi23', 'sxi13', 'dxiCP', 'b1', 'b2',
        'b3'. For 4, additionally 'dxi13' (replacing 'dxiCP'), 'sxi14', 'dxi14', 'sxi24', 'dxi24',
        'sxi34', 'b4'. For 5, additionally 'sxi15', 'dxi15', 'sxi25', 'sxi35', 'dxi35', 'b5'.
    h_liv : list or np.ndarray
        Precomputed LIV Hamiltonian, required when ``num_flavors`` exceeds
        ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS`` (in which case there is nothing to
        unpack from ``liv_params`` and this function returns None).

    Returns
    -------
    np.ndarray or None
        The unpacked LIV parameters, in the order expected by the matching
        ``hamiltonian_{N}nu_liv_energy_independent`` function; or None if ``num_flavors`` exceeds
        ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS`` (the caller uses ``h_liv`` directly).
    """

    if (num_flavors == 2):
        try:
            Lambda = liv_params['Lambda']
            if (Lambda <= 0.0):
                raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                    ": Lambda must be positive.")
            sxi = liv_params['sxi']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            n_liv = liv_params['n_liv']
            return np.array([sxi, b1, b2, Lambda, n_liv])
        except KeyError :
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since "+ \
                    "num_flavors == 2, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi', 'b1', 'b2', 'Lambda', 'n_liv'.")
    elif (num_flavors == 3):
        try:
            Lambda = liv_params['Lambda']
            if (Lambda <= 0.0):
                raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                    ": Lambda must be positive.")
            sxi12 = liv_params['sxi12']
            sxi23 = liv_params['sxi23']
            sxi13 = liv_params['sxi13']
            dxiCP = liv_params['dxiCP']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            b3 = liv_params['b3']
            n_liv = liv_params['n_liv']
            return np.array([sxi12, sxi23, sxi13, dxiCP, b1, b2, b3, Lambda, n_liv])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 3, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi12', 'sxi23', 'sxi13', 'dxiCP'," + \
                    " 'b1', 'b2', 'b3', 'Lambda', and 'n_liv'.")
    elif (num_flavors == 4):
        try:
            Lambda = liv_params['Lambda']
            if (Lambda <= 0.0):
                raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                    ": Lambda must be positive.")
            sxi12 = liv_params['sxi12']
            sxi23 = liv_params['sxi23']
            sxi13 = liv_params['sxi13']
            dxi13 = liv_params['dxi13']
            sxi14 = liv_params['sxi14']
            dxi14 = liv_params['dxi14']
            sxi24 = liv_params['sxi24']
            dxi24 = liv_params['dxi24']
            sxi34 = liv_params['sxi34']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            b3 = liv_params['b3']
            b4 = liv_params['b4']
            n_liv = liv_params['n_liv']
            return np.array([sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2,
                b3, b4, Lambda, n_liv])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 4, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi12', 'sxi23', 'sxi13', 'dxi13'," + \
                    " 'sxi14', 'dxi14', 'sxi24', 'dxi24', 'sxi34', 'b1', 'b2', 'b3', 'b4'," + \
                    " 'Lambda', and 'n_liv'.")
    elif (num_flavors == 5):
        try:
            Lambda = liv_params['Lambda']
            if (Lambda <= 0.0):
                raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                    ": Lambda must be positive.")
            sxi12 = liv_params['sxi12']
            sxi23 = liv_params['sxi23']
            sxi13 = liv_params['sxi13']
            dxi13 = liv_params['dxi13']
            sxi14 = liv_params['sxi14']
            dxi14 = liv_params['dxi14']
            sxi15 = liv_params['sxi15']
            dxi15 = liv_params['dxi15']            
            sxi24 = liv_params['sxi24']
            dxi24 = liv_params['dxi24']
            sxi25 = liv_params['sxi25']
            sxi34 = liv_params['sxi34']
            sxi35 = liv_params['sxi35']
            dxi35 = liv_params['dxi35']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            b3 = liv_params['b3']
            b4 = liv_params['b4']
            b5 = liv_params['b5']
            n_liv = liv_params['n_liv']
            return np.array([sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, 
                sxi25, sxi34, sxi35, dxi35, b1, b2, b3, b4, b5, Lambda, n_liv])
        except KeyError:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 5, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi12', 'sxi23', 'sxi13', 'dxi13'," + \
                    " 'sxi14', 'dxi14', 'sxi15', 'dxi15', 'sxi24', 'dxi24', 'sxi25' 'sxi34', " + \
                    " 'sxi35', 'dxi35', 'b1', 'b2', 'b3', 'b4', 'b5', 'Lambda', and 'n_liv'.")
    elif (num_flavors > gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
        print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": the number of " + \
            "flavors passed (num_flavors = " + str(num_flavors) + \
            ") exceeds the maximum number for which Magnus has predefined vacuum Hamiltonians " + \
            "(globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
            str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + "). Will use the Hamiltonian provided " + \
            "in h_liv.")
        if (h_liv is None):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": provided " + \
                "h_liv is None.")
        # num_flavors exceeds the predefined range: the caller builds its Hamiltonian directly from
        # h_liv instead of from a flat parameter list, so there is nothing to unpack here.
        return None
    elif (num_flavors < 1):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": num_flavors must be " + \
            ">= 2.")


# def chunkify(lst, n):
#     """Yield successive n-sized chunks from lst."""
#     for i in range(0, len(lst), n):
#         yield lst[i:i + n]


#-----------------------------------------------------------------------
# Primordial functions
#-----------------------------------------------------------------------

class _PositionProfileCache:
    r"""Memoizes a position-profile function on repeated position grids.

    Across an energy scan, the matter term of the Hamiltonian is evaluated on
    the same position grids for every energy (only the 1/E vacuum part
    changes), and across the adaptive refinement loops the same grids recur
    between neighboring points.  This tiny cache stores the profile values of
    the most recent grids, keyed by the exact grid contents, so the
    (comparatively expensive) density-profile chain runs once per distinct
    grid instead of once per Hamiltonian evaluation.  Scalar evaluations are
    passed through uncached.
    """

    def __init__(self, func: Callable, maxsize: Optional[int]=8):
        self.func = func
        self._cache = {}
        self._keys = []
        self._maxsize = maxsize
        # Preserve any profile tag set by the caller (e.g., is_exp_density_profile/l_scale from
        # matter.exp_density_profile, propagated through matter.vcc_func_from_rho_func), so that
        # wrapping in this cache does not hide it from _osc_prob_ip_exp_dispatch.
        if getattr(func, 'is_exp_density_profile', False):
            self.is_exp_density_profile = True
            self.l_scale = func.l_scale

    def __call__(self, l: Union[int, float, np.ndarray]):
        if np.ndim(l) == 0:
            return self.func(l)
        l = np.asarray(l, dtype=float)
        key = (l.shape, l.tobytes())
        val = self._cache.get(key)
        if val is None:
            val = np.asarray(self.func(l))
            self._cache[key] = val
            self._keys.append(key)
            if len(self._keys) > self._maxsize:
                self._cache.pop(self._keys.pop(0), None)
        return val


def compute_evolution_operator(
    H_func: Callable,
    t_slab: Union[list, np.ndarray],
    n_tpts_per_slab: int,
    magnus_exp_order: int,
    **kwargs
) -> np.ndarray:
    r"""Computes the evolution operator inside a given time slab.  This functions is not designed to
    be called directly by the user, but rather internally by :func:`osc_prob`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of time or position that returns a square matrix (NumPy array).
    t_slab : list or np.ndarray
        Start and end times or positions of the slab, ``[t0, t1]``.
    n_tpts_per_slab : int
        Number of time points inside the slab at which to evaluate ``H_func`` to numerically
        compute the integrals required by the Magnus expansion.
    magnus_exp_order : int
        Highest order of the Magnus expansion used to compute the evolution operator (should not
        exceed ``globaldefs.MAGNUS_EXP_ORDER_MAX``).
    \**kwargs
        Additional arguments passed to :func:`magnus.magnus.magnus_expansion` (e.g.,
        ``integration_method``).

    Returns
    -------
    np.ndarray
        The evolution operator for the given time slab.
    """
    if t_slab[1] > t_slab[0]:
        return magnus.magnus_expansion(
            lambda t: -1j * H_func(t),
            t0=t_slab[0],
            t1=t_slab[1],
            # t_slabs=[t_slab],
            order=magnus_exp_order,
            n_tpts=n_tpts_per_slab,
            **kwargs,
        )
    else:  # t_slab[1] == t_slab[0]
        n = H_func(t_slab[0]).shape[0]
        return np.eye(n)


def compute_evolution_operator_multiple_slabs(
    H_func: Callable,
    t_slabs: Union[list, np.ndarray],
    n_tpts_per_slab: int,
    magnus_exp_order: int,
    **kwargs
) -> np.ndarray:
    r"""Computes the evolution operators of a chain of time slabs.  This function is not designed
    to be called directly by the user, but rather internally by :func:`osc_prob`.

    All slabs are computed at once by :func:`magnus.magnus.magnus_expansion_multislab`, which
    batches the Hamiltonian evaluation, the quadrature, the commutator algebra, and the matrix
    exponentials over the slab axis.  Slabs of zero width yield identity operators.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of time or position that returns a square matrix (NumPy array).
        If it also accepts an array of times (returning a stack of matrices), the vectorized form
        is detected and used automatically for speed.
    t_slabs : list or np.ndarray
        Pairs specifying the start and end times or positions of each slab, i.e.,
        ``[[t0, t1], [t1, t2], ...]``.
    n_tpts_per_slab : int
        Number of time points inside each slab at which to evaluate ``H_func`` to numerically
        compute the integrals required by the Magnus expansion (ignored by the 'gl' integration
        method).
    magnus_exp_order : int
        Highest order of the Magnus expansion used to compute the evolution operator (should not
        exceed ``globaldefs.MAGNUS_EXP_ORDER_MAX``).
    \**kwargs
        Additional arguments passed to :func:`magnus.magnus.magnus_expansion_multislab` (e.g.,
        ``integration_method``).

    Returns
    -------
    np.ndarray
        Evolution operators, shape (n_slabs, dim, dim), ordered like ``t_slabs`` (earliest slab
        first). Note that the time-ordered product over the chain is
        ``U_total = U[-1] @ ... @ U[1] @ U[0]``, i.e., the last slab is the leftmost factor.
    """
    def hh(t):
        return -1j * H_func(t)

    return magnus.magnus_expansion_multislab(hh, t_slabs, n_tpts_per_slab=n_tpts_per_slab,
        order=magnus_exp_order, **kwargs)


def osc_prob(
    H_func: Union[Callable, np.ndarray], 
    t_ini: Union[int, float], 
    t_fin: Union[int, float], 
    n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='gl', 
    rtol: Optional[Union[int, float]]=1.e-3, 
    atol: Optional[Union[int, float]]=1.e-3, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=None, 
    min_n_tpts_per_slab: Optional[int]=2, 
    max_n_tpts_per_slab: Optional[int]=500, 
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0, 
    A_eval_mode: Optional[str]=None,
    convergence_info: Optional[Dict]=None,
    t_breakpoints: Optional[Union[list, np.ndarray]]=None,
    strict_convergence: Optional[bool]=False,
    **kwargs
) -> np.ndarray:
    r"""Computes and returns the neutrino oscillation probability.

    Computes the oscillation probability of neutrinos starting at time
    (or position) ``t_ini`` and ending at time (or position) ``t_fin``.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable or np.ndarray
        The Hamiltonian, which is a function of time or position that 
        returns a square matrix (a NumPy array). The Hamiltonian can 
        have complex-valued entries.
    t_ini : int or float
        Initial time or position of the neutrino.
    t_fin : int or float
        Final time or position of the neutrino.
    n_slabs : int, optional
        Number of slabs, or subintervals, into which the interval 
        [``t_ini``, ``t_fin``] is partitioned in order to compute the 
        neutrino evolution operators. A higher value of ``n_slabs`` 
        yields a more accurate probability.

        If no target tolerance is requested (i.e., if ``rtol`` and 
        ``atol`` are both ``None``), then the given value of ``n_slabs`` 
        is the final number of slabs used in the computation.

        If a target tolerance is requested (i.e., if either ``rtol`` or
        ``atol`` is not ``None``), then the given value of ``n_slabs``
        acts as a *floor*: the number of slabs is increased
        progressively, starting from ``max(min_n_slabs, n_slabs)``,
        until the tolerance is achieved or until we hit ``max_n_slabs``,
        whichever happens first.  The refinement never runs coarser than
        what was asked for, so a caller who knows the feature scale of
        their profile can state it here and have it respected.  With the
        default, ``n_slabs = 1``, the floor is inactive and refinement
        starts at ``min_n_slabs`` as before.
    n_tpts_per_slab : int, optional
        Number of time-points inside the slab at which to evaluate 
        H_func in order to numerically compute the integrals over time 
        required by the Magnus expansion. A higher value of 
        ``n_tpts_per_slab`` yields a more accurate probability.
    t_slab_edges : list or np.ndarray, optional
        Optional list of pairs [[t0, t1], [t1, t2], ...] with the edges
        of each time slab.  If given, it overrides ``n_slabs`` and the
        uniform partitioning of [``t_ini``, ``t_fin``]; the user must
        ensure that the slabs chain without gaps.  If a tolerance is
        requested, only ``n_tpts_per_slab`` is grown (the user-provided
        edges are kept fixed).
    magnus_exp_order : int, optional
        Order at which the Magnus expansion is truncated (1 to
        ``globaldefs.MAGNUS_EXP_ORDER_MAX``).
    n_jobs : int, optional
        Number of parallel joblib workers used to compute the per-slab
        evolution operators.  With the default, ``n_jobs = 1``, all
        slabs are computed in a single vectorized (batched) call, which
        is usually fastest; use ``n_jobs > 1`` only for very expensive
        Hamiltonian functions.
    integration_method : str, optional
        'gl' for Gauss-Legendre collocation, which needs only 1, 2, or 3
        Hamiltonian evaluations per slab for orders <= 2, <= 4, <= 6, and
        ignores ``n_tpts_per_slab``; or 'trapezoid'/'simpson' for cumulative
        quadrature over ``n_tpts_per_slab`` points per slab. Default: 'gl'.
    rtol : int or float, optional
        Target relative tolerance of the probability matrix between
        successive refinement loops.  Set both ``rtol`` and ``atol`` to
        ``None`` to run once with the given fixed parameters.  If only
        one of the two is ``None``, it is treated as 0.
    atol : int or float, optional
        Target absolute tolerance; see ``rtol``.
    growth_factor_n_slabs : int or float, optional
        Factor by which ``n_slabs`` is multiplied on each refinement
        loop (used only when a tolerance is requested).
    growth_factor_n_tpts_per_slab : int or float, optional
        Factor by which ``n_tpts_per_slab`` is multiplied on each
        refinement loop (used only when a tolerance is requested).
    max_num_loops : int, optional
        Maximum number of refinement loops.
    min_n_slabs : int, optional
        Number of slabs used in the first refinement loop.
    max_n_slabs : int, optional
        Maximum allowed number of slabs.  If None (default), a cap appropriate to
        ``integration_method`` is used: 20000 for 'gl', 2000 for the cumulative-quadrature
        methods (see :data:`MAX_N_SLABS_DEFAULT`).  'gl' costs 1-3 Hamiltonian evaluations
        per slab against the quadrature methods' ``n_tpts_per_slab``, so the same cost
        budget buys it far more slabs.  An explicit value is always used as given.
    min_n_tpts_per_slab : int, optional
        Number of time points per slab in the first refinement loop.
    max_n_tpts_per_slab : int, optional
        Maximum allowed number of time points per slab.
    validate_input : bool, optional
        If True, validate the input parameters (set to False for a
        small speed-up once a call is known to be well-formed).
    save_log : bool, optional
        If True, also write all messages to the log file.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no
        ``file_log`` object is given).
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning.
    new_recursion_limit : int, optional
        If not None, raise Python's recursion limit to this value.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the
        refinement loops).
    A_eval_mode : str, optional
        How the Hamiltonian can be evaluated: 'vector' (accepts an
        array of positions), 'constant', or 'scalar'.  Determined
        automatically when None; pass it explicitly (e.g., from
        :func:`magnus.magnus.probe_eval_mode`) to skip the probe.
    convergence_info : Dict, optional
        If a dict is passed, it is filled in place with the refinement
        parameters of the returned probability ('n_slabs',
        'n_tpts_per_slab'), which callers can use to warm-start
        neighboring computations.
    t_breakpoints : list or np.ndarray, optional
        Optional positions at which the Hamiltonian is known to be
        non-smooth (e.g., density discontinuities such as the PREM
        layer boundaries).  They are inserted as mandatory slab edges
        into the automatically generated slab grid at every refinement
        level, so that the quadrature never integrates across them.
        Ignored when ``t_slab_edges`` is given explicitly.
    strict_convergence : bool, optional
        Require the refinement ladder to agree **twice in a row** before
        declaring convergence, instead of once.  Default: False.

        *What the ladder normally does.*  With a tolerance requested,
        ``osc_prob`` computes the probability on a grid of ``n_slabs``
        slabs, then again on a finer one (``n_slabs`` grows by
        ``growth_factor_n_slabs`` each time), and returns as soon as two
        successive grids agree within ``rtol``/``atol``.  The assumption
        is that agreement between successive refinements means the answer
        has stopped changing because it has converged.

        *When that assumption fails.*  It is only safe while the sequence
        is settling down.  If the grid is still too coarse to resolve the
        Hamiltonian, successive refinements do not approach the answer
        smoothly -- they jump around it -- and two neighboring jumps can
        land close together by coincidence.  ``np.allclose`` cannot tell
        that apart from convergence, so the ladder stops early and returns
        a plausible, exactly unitary, *wrong* answer with no warning.
        Measured example (2 flavors, solar exponential profile, 10 MeV
        over one solar radius, default ``rtol=atol=1e-3``): the errors at
        successive levels run 5.9e-02, 3.8e-03, 1.6e-02, 1.7e-02, 8.1e-03,
        4.5e-03, 3.5e-06.  Levels 3 and 4 agree to 1.1e-03 -- inside the
        requested tolerance -- while both are wrong by ~1.6e-02, and the
        next level moves by 2.5e-02.

        *What this flag changes.*  Convergence is declared only after two
        consecutive agreements, so a lone coincidence is vetoed by the
        level that follows it.  On the example above the ladder continues
        to ``n_slabs = 20000`` and an error of 6.8e-08.

        *What it costs.*  One extra refinement level, each costing about
        ``growth_factor_n_slabs`` times the last: measured median **1.53x**
        (worst 1.75x) on calls whose first agreement was already genuine.
        On calls this actually rescues it costs 3.5-8.6x, because there the
        second agreement is several levels away -- that cost is paid only
        where the answer would otherwise have been wrong.

        *When you do not need it.*  If the quantity you care about is an
        average over many oscillations -- the usual case for solar
        neutrinos, where the survival probability oscillates thousands of
        times along the trajectory -- most of the error this guards
        against is in the *phase* and cancels in the average.  On the 10
        MeV example the pointwise error of 2.5e-02 becomes 1.9e-04 once
        averaged over 25 oscillations.  Prefer ``average=True`` on the
        wrapper functions (see :func:`osc_prob_matter_std_potential`),
        which computes the phase-averaged probability directly and far
        more cheaply.  Use ``strict_convergence`` when the oscillating
        probability itself is the answer you want -- a probability-versus-
        baseline or versus-energy curve, an oscillogram, or a fixed
        baseline and energy.

        *What it does not fix.*  A refinement ladder of any strictness is
        powerless against an **incomplete** ``t_breakpoints`` list.  If the
        Hamiltonian is discontinuous somewhere that is not marked as a slab
        edge, every level integrates across that discontinuity, successive
        levels can agree to machine precision, and the shared answer is
        simply wrong.  Measured on a 50-wall piecewise-constant profile
        whose first boundary was left unmarked: the error sat at 1.6e-02,
        bit-identical from ``n_slabs = 4`` through 32, and adding the one
        missing edge moved it to 3.6e-12 at every slab count.  When a
        profile is discontinuous, marking *every* discontinuity -- including
        where it switches on and off, which may lie inside the trajectory --
        is worth more than any amount of refinement.
    \**kwargs
        Additional arguments passed through to the Magnus-expansion
        routines

    Returns
    -------
    np.ndarray
        NumPy array containing the probability matrix of the same 
        dimensions as the Hamiltonian, ``H_func``.
    """

    # Checked before anything forwards **kwargs onwards, and regardless of validate_input:
    # these two keys are rejected several hops away otherwise, by a function the caller never
    # named (see _reject_parameter_set_metadata).
    _reject_parameter_set_metadata(kwargs, 'osc_prob')

    # Validate input; set validate_input to False for speed-up.
    # None means 'use the cap appropriate to this integration method'
    # (see MAX_N_SLABS_DEFAULT); an explicit value always wins.
    max_n_slabs = _resolve_max_n_slabs(max_n_slabs, integration_method)
    if validate_input:

        if (t_fin < t_ini): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: t_fin must be >=" + \
                " t_ini.")

        if (magnus_exp_order < 1): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: magnus_exp_order " + \
                "must be >= 1.")

        if ((rtol is not None) and (rtol <= 0.0)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: rtol must be None " + \
                "or > 0.0.")

        if ((atol is not None) and (atol <= 0.0)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: atol must be None " + \
                "or > 0.0.")

        if ((rtol is not None) and (atol is not None) and (growth_factor_n_slabs < 1.0)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: " + \
                "growth_factor_n_slabs must be >= 1.0.")

        if ((rtol is not None) and (atol is not None) and (growth_factor_n_tpts_per_slab < 1.0)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: " + \
                "growth_factor_n_tpts_per_slab must be >= 1.0.") 

        if ( ((rtol is not None) and (atol is not None)) and \
            ((growth_factor_n_slabs == 1.0) and (growth_factor_n_tpts_per_slab == 1.0)) ): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: since a target " + \
                "tolerance has been requested, either growth_factor_n_slabs, " + \
                "growth_factor_n_tpts_per_slab, or both must be > 1.")

        if ((rtol is not None) and (atol is not None) and (max_num_loops <= 1)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: max_num_loops must" + \
                " be > 1.")

        if ((rtol is not None) and (atol is not None) and (max_n_slabs <= 1)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: max_n_slabs must " + \
                "be > 1.")

        if ((rtol is not None) and (atol is not None) and (max_n_slabs <= 2)): 
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: max_n_tpts_per_slab" +\
                " must be > 2.")

        if ((callable(H_func)) and (_n_required_params(H_func) > 1)):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: the provided H_func" +\
                " is a function of more than one parameter")

        H_test = H_func(t_ini) if callable(H_func) else H_func

        if not isinstance(H_test, np.ndarray):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: H_func must be a " + \
                "NumPy (if the Hamiltonian is time-independent) or must return a NumPy array.")

        if H_test.shape[0] != H_test.shape[1]:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob: H_func must be a " + \
                "square matrix (if the Hamiltonian is time-independent) or must return a " + \
                "square matrix.")

    # If only one of rtol and atol was given (i.e., the other one is None), set the missing one to
    # 0.0, so that the requested tolerance is driven by the one that was given.  (Internally, the
    # code treats "a tolerance was requested" as both rtol and atol being not None.)
    if (rtol is None) != (atol is None):
        rtol = 0.0 if rtol is None else rtol
        atol = 0.0 if atol is None else atol

    # The Gauss-Legendre integration method ('gl') uses a fixed, small number of Hamiltonian
    # evaluations per slab (1, 2, or 3, depending on magnus_exp_order), so n_tpts_per_slab plays no
    # role: the accuracy is controlled by the number of slabs only.  Neutralize the growth of
    # n_tpts_per_slab so that the adaptive loop below grows only n_slabs.
    if integration_method == 'gl':
        growth_factor_n_tpts_per_slab = 1.0
        min_n_tpts_per_slab = max_n_tpts_per_slab = n_tpts_per_slab = 2

    # If there is no file object given (i.e., if file_log is None), open a log file if requested
    if file_log is None:
        file_log = open(filename_log, 'w') if save_log else None

    # Print a list of all the parameters passed to the osc_prob function and their values
    if (verbose > 1):
        print_run_parameters(H_func, t_ini, t_fin, n_slabs, n_tpts_per_slab, t_slab_edges,
            magnus_exp_order, n_jobs, integration_method, rtol, atol, growth_factor_n_slabs,
            growth_factor_n_tpts_per_slab, max_num_loops, min_n_slabs, max_n_slabs, 
            min_n_tpts_per_slab, max_n_tpts_per_slab, validate_input, save_log, filename_log, 
            new_recursion_limit, verbose, file_log)

    # Note: new_recursion_limit is accepted for backward compatibility but no longer used; the
    # probability calculation is fully iterative (nothing recurses), so there is no need to raise
    # Python's recursion limit.

    # osc_prob runs at a fixed magnus_exp_order; the requested tolerance is reached by
    # refining the number of slabs (and, for the quadrature methods, the points per
    # slab), never by raising the order.  See the note on choosing an order in
    # docs/source/methodology.rst.

    loop_count = 1 # Loop counter
    # Probability matrix of the current and the previous refinement loop; the two are compared
    # to decide convergence.  Both are None until the first loop has produced a matrix, which is
    # why the early-exit checks inside the loop are guarded on loop_count > 1.
    P = None
    P_old = None
    last_gap = None
    # Consecutive refinement levels that have agreed within (rtol, atol) so far.  The ladder
    # normally returns on the first agreement; strict_convergence requires two in a row, so that
    # a coincidental agreement between two levels of a sequence that is still jumping around is
    # vetoed by the level after it.  See the strict_convergence entry in the docstring above.
    n_agreements = 0
    agreements_required = 2 if strict_convergence else 1
    # Copy this to remember whether the function was originally called with predefine slab edges,
    # or whether we can increase the number of edges (n_slabs) progressively to reach tolerance
    t_slab_edges_original = t_slab_edges 
    # Flags to signal whether a loop has been run with n_slabs == max_n_slabs or 
    # n_tpts_per_slab = max_n_tpts_per_slab
    ran_with_max_n_slabs, ran_with_max_n_tpts_per_slab = False, False 
    # Flags to signal whether we have already printed the warning that we have reached 
    # n_slabs == max_n_slabs or n_tpts_per_slab = max_n_tpts_per_slab, so as not to print it again
    warned_reached_max_n_slabs, warned_reached_max_n_tpts_per_slab = False, False

    # If a tolerance is requested, start the iterations at the floor on the slab count.  That floor
    # used to be min_n_slabs alone, with the caller's n_slabs discarded outright; it is now the
    # larger of the two.  Discarding it was how a profile with 50 density walls, called with
    # n_slabs=150, came to be integrated on 4 slabs and declared converged.  The seed that replaced
    # the caller's number, magnus.suggest_n_slabs, measures the *integral* of the Hamiltonian along
    # the path, and an integral is blind to structure that averages out: that profile accumulates
    # only ~9 radians over the whole trajectory, so it was seeded with 2 slabs.  Below the count
    # that resolves the walls the ladder does not converge, it thrashes -- 0.43, 0.13, 0.13, 0.64,
    # 0.12 at 2, 3, 4, 5, 6 slabs -- and np.allclose fired on the accidental 3-vs-4 agreement,
    # returning an answer wrong by 0.855 in probability.  A stricter rtol is the wrong lever: it
    # tightens a comparison between two answers that both failed to see the profile.  Resolving the
    # profile is the right one, and the caller is who knows its feature scale (t_breakpoints is the
    # sharper tool still, where the features sit at known positions).  With the default n_slabs=1
    # the floor is inactive and nothing changes.  Clipped at max_n_slabs so that a floor above the
    # cap cannot make the ladder step *down* on its first growth; the usual "reached max_n_slabs"
    # warning then fires, as it should.
    if ((rtol is not None) and (atol is not None)):
        min_n_slabs = int(min(max(min_n_slabs, n_slabs), max_n_slabs))
        n_slabs = min_n_slabs
        n_tpts_per_slab = min_n_tpts_per_slab

    # The provided Hamiltonian, H_func, can be either a single-parameter function (of the neutrino
    # position) or, if time-independent, a constant NumPy array (e.g., for oscillations in vacuum
    # or in matter with constant density).  In the latter case, we use this constant Hamiltonian to
    # build a dummy one-parameter function of position that we will need later to call the function
    # compute_evolution_operator.  In this case, first-order Magnus expansion is enough, and so we
    # can overwrite the parameters provided to n_slabs = 1, n_tpts_per_slab = 2, rtol = None, 
    # atol = None for speed-up.
    if not callable(H_func):
        H = np.copy(H_func)
        def H_func(l: float) -> np.ndarray:
            return H
        magnus_exp_order = 1
        n_slabs = 1
        n_tpts_per_slab = 2
        rtol = None
        atol = None
        n_jobs = 1 # No need to parallelize for this simple computation in a single slab
        # A single slab is exact for a constant Hamiltonian, so drop any user-provided slab edges
        t_slab_edges = None
        t_slab_edges_original = None
        if verbose > 0:
            for f in [None, file_log] if save_log else [None]:
                warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                print("\n" + warn_msg + " The provided Hamiltonian is time-independent. " + \
                    "Overwriting the run parameters to magnus_exp_order = 1, n_slabs = 1, " + \
                    "n_tpts_per_slab = 2, rtol = None, atol = None, and n_jobs = 1 for speed-up.",
                    file=f)

    # If the user provided the slab edges explicitly, the number of slabs is set by them
    if t_slab_edges_original is not None:
        n_slabs = len(t_slab_edges_original)

    # Determine once how the Hamiltonian can be evaluated (vectorized over an array of positions,
    # constant, or scalar-only), so that the Magnus kernel does not have to re-probe it on every
    # refinement iteration below.
    if A_eval_mode is None:
        A_eval_mode = magnus.probe_eval_mode(lambda t: -1j*H_func(t), t_ini, t_fin)

    # Physics-informed starting number of slabs (Gauss-Legendre method only): rather than always
    # starting the refinement from min_n_slabs and climbing the geometric ladder, start from an
    # estimate based on the accumulated phase (see magnus.suggest_n_slabs).  min_n_slabs still
    # acts as a lower bound, so warm starts provided by the caller take precedence when they are
    # larger.  For the quadrature methods ('trapezoid', 'simpson') the accuracy is governed
    # jointly by n_slabs and n_tpts_per_slab, and seeding only the slab count unbalances that
    # ladder, so the seed is not applied there.
    if ((rtol is not None) and (atol is not None) and (t_slab_edges_original is None) and \
        (integration_method == 'gl')):
        n_slabs = int(np.clip(max(min_n_slabs,
            magnus.suggest_n_slabs(lambda t: -1j*H_func(t), t_ini, t_fin,
                A_eval_mode=A_eval_mode)), 1, max_n_slabs))

    while True:

        # These checks only apply when osc_prob is run with a requested tolerance (rtol, atol) that
        # should be achieved.
        if ((rtol is not None) and (atol is not None)):
            # Reached maximum allowed number of loops: exit loop, return the probability matrix.
            # Guarded on loop_count > 1 because these are refinement limits: there is nothing to
            # return until at least one loop has produced a probability matrix.  Without the
            # guard, max_num_loops < 1 with validate_input=False (the validator rejects it
            # otherwise) reached this return before P existed and raised UnboundLocalError.
            if (loop_count > 1) and (loop_count > max_num_loops):
                warnings.warn("osc_prob: requested tolerance not achieved "
                    "(max_num_loops reached), " + _shortfall_phrase(last_gap, rtol, atol) +
                    "; the returned probabilities may be inaccurate. Raise max_num_loops, or "
                    "loosen rtol/atol to what the last two levels actually agreed to. Shown "
                    "once per session.", ToleranceNotAchievedWarning, stacklevel=2)
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg + " Number of loops (loop_count = " + \
                            str(loop_count-1) + ") reached maximum allowed (max_num_loops = " + \
                            str(max_num_loops) + "). Requested tolerance not achieved. Try " + \
                            "increasing max_num_loops.\n",
                            file=f)
                if save_log and close_file_log_upon_exit: file_log.close()
                return P
            # Reached maximum allowed number of slabs: continue execution
            if (n_slabs == max_n_slabs):
                if ((verbose > 0) and not warned_reached_max_n_slabs):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg +  " Number of slabs (n_slabs) reached maximum " + \
                            "allowed (max_n_slabs = " + str(max_n_slabs) + ").", file=f)
                        warned_reached_max_n_slabs = True
            # Reached maximum allowed number of time-points per slab: continue execution
            if (n_tpts_per_slab == max_n_tpts_per_slab):
                if ((verbose > 0) and not warned_reached_max_n_tpts_per_slab):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg + " Number of time-points per slab " + \
                            "(n_tpts_per_slab) reached maximum allowed (max_n_tpts_per_slab = " + \
                            str(max_n_tpts_per_slab) + ").", file=f)
                        warned_reached_max_n_tpts_per_slab = True
            # Reached maximum allowed number of slabs and maximum allowed number of time-points per
            # slab: exit loop, return the probability matrix
            if (loop_count > 1) and ran_with_max_n_slabs and ran_with_max_n_tpts_per_slab:
                # 'gl' pins n_tpts_per_slab (it uses a fixed 1-3 nodes per slab), so only
                # max_n_slabs is a meaningful knob for it; naming max_n_tpts_per_slab in the
                # message would send the reader after a setting that cannot help them.
                knobs = ("max_n_slabs" if integration_method == 'gl'
                         else "max_n_slabs and max_n_tpts_per_slab")
                warnings.warn("osc_prob: requested tolerance not achieved (" + knobs +
                    " reached), " + _shortfall_phrase(last_gap, rtol, atol) + ", so "
                    "convergence could not be verified by successive refinement and the "
                    "returned probabilities may be inaccurate. Raise " + knobs + ". This can "
                    "happen for very large accumulated phases, e.g., low-energy neutrinos over "
                    "very long baselines, or eV-scale sterile splittings over Earth-crossing "
                    "baselines. If the profile has a density jump or a kink, pass "
                    "t_breakpoints there as well -- a slab straddling one is never fixed by "
                    "more slabs. Shown once per session.",
                    ToleranceNotAchievedWarning, stacklevel=2)
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg + " Number of slabs (n_slabs) and time-points " + \
                            "per slab (n_tpts_per_slab) reached maximum allowed (max_n_slabs = " + \
                            str(max_n_slabs) + ", max_n_tpts_per_slab = " + \
                            str(max_n_tpts_per_slab) + ").", file=f)
                        print("   " + warn_msg + " Returning probability, but requested " + \
                            "tolerance (rtol = " + str(rtol) + ", atol = " + str(atol) + \
                            ") not achieved. Try increasing max_n_slabs or max_n_tpts_per_slab.\n",
                            file=f)
                if save_log and close_file_log_upon_exit: file_log.close()
                return P

        # The array (or list) t_slab_edges contains user-provided pairs of start and end times, 
        # [ti, tf]_k, that define the initial and final times of each of the k-th time slab.  It is 
        # up to the user to ensure that the chain of time slabs covers the full range [t_ini, t_fin] 
        # without leaving gaps.  I.e., the user should ensure that ti_{k+1} = tf_k.  
        if (t_slab_edges_original is None):
            # If t_slab_edges == None, then divide the interval [t_ini, t_fin] evenly into a number
            # n_slabs of time slabs.  Any t_breakpoints inside the interval (e.g., density
            # discontinuities) are inserted as additional mandatory slab edges: high-order
            # quadrature converges at its nominal order only if the Hamiltonian is smooth inside
            # each slab.
            grid = np.linspace(t_ini, t_fin, n_slabs+1)
            if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
                bp = np.atleast_1d(np.asarray(t_breakpoints, dtype=float))
                bp = bp[(bp > t_ini) & (bp < t_fin)]
                grid = np.unique(np.concatenate([grid, bp]))
            t_slab_edges = np.column_stack([grid[:-1], grid[1:]])

        # Within each slab, t_slab, we use n_tpts_per_slab time-evaluations to compute the integrals
        # of the Magnus expansion, from t_slab[0] to t_slab[1].  U_chain contains the chain of time-
        # ordered evolution operators, each computed in one time slab.  All slabs are computed in a
        # single batched call.  (Note: n_jobs is accepted for backward compatibility, but the
        # per-slab parallelization it used to trigger here has been retired: the batched kernel is
        # faster than distributing the small per-slab tasks over joblib workers.  Parallelism over
        # (energy, L) points is available in osc_prob_energy_baseline instead.)
        U_chain = compute_evolution_operator_multiple_slabs(H_func, t_slab_edges,
            n_tpts_per_slab, magnus_exp_order, integration_method=integration_method,
            A_eval_mode=A_eval_mode, **kwargs)

        # Now compute the time-ordered product of all evolution operators across all slabs.  The
        # neutrino crosses the slabs in the order in which they appear in U_chain (earliest first),
        # so the total operator is U_tot = U_chain[-1] @ ... @ U_chain[1] @ U_chain[0]: the operator
        # of the *last* slab is the leftmost factor.  (functools.reduce is used instead of
        # np.linalg.multi_dot because all factors are square matrices of the same size, for which
        # multi_dot wastes time computing an optimal parenthesization that does not exist.)
        Utot = reduce(np.matmul, U_chain[::-1]) if len(U_chain) > 1 else U_chain[0]

        # Using Utot, compute all the survival and transition probabilities in a probability matrix
        # P = (|Utot|^2).T and return that matrix, so that P[nu_i][nu_f] = |Utot[nu_f][nu_i]|^2.
        P = np.transpose(Utot.real**2 + Utot.imag**2)

        # Record the refinement parameters of this (latest) computation, so that callers (e.g.,
        # osc_prob_energy_baseline) can warm-start neighboring points
        if convergence_info is not None:
            convergence_info['n_slabs'] = n_slabs
            convergence_info['n_tpts_per_slab'] = n_tpts_per_slab

        # If no target relative tolerance (rtol) or absolute tolerance (atol) of the probability is
        # requested, then return the result obtained already.  If, instead, a target tolerance is
        # requested, then increase the number of points per slab approximately by the factor
        # growth_factor_n_tpts_per_slab, and repeat the probability calculation until the desired
        # tolerance is achieved.
        if ((rtol is None) and (atol is None)): # No target tolerance requested: return right away
            if save_log and close_file_log_upon_exit: file_log.close()
            return P
        else: # Target tolerance requested: iterate until tolerance is achieved
            if (verbose > 1):
                for f in [None, file_log] if save_log else [None]:
                    if (loop_count == 1):
                        print("\nRunning loops until requested rtol and atol are achieved:", file=f)
                    print("   Loop #" + str(loop_count) + ":", file=f)
                    print("      magnus_exp_order = " + str(magnus_exp_order), file=f)                    
                    print("      n_slabs = " + str(n_slabs), file=f)
                    print("      n_tpts_per_slab = " + str(n_tpts_per_slab), file=f)
            if P_old is not None:
                # Compare the new and old probability matrices element-wise.  A run of agreements
                # is tracked rather than a single one: a disagreement resets it, so with
                # strict_convergence the two agreements must be genuinely consecutive.
                # Kept so the tolerance-not-achieved warnings below can say how far from
                # converged the refinement stopped, rather than only that it stopped.  The
                # comparison is being made anyway; this is the number it is made on.
                last_gap = float(np.max(np.abs(P - P_old)))
                if np.allclose(P, P_old, rtol=rtol, atol=atol):
                    n_agreements += 1
                else:
                    n_agreements = 0
                if n_agreements >= agreements_required:
                    if (verbose > 0):
                        for f in [None, file_log] if save_log else [None]:
                            tol_msg = gd.TOL_MSG_IN_COLOR if f is None else gd.TOL_MSG_NO_COLOR
                            print("   " + tol_msg + " (for fixed magnus_exp_order "+ \
                                "= " + str(magnus_exp_order) + "): rtol = " + str(rtol) + \
                                ", atol = " + str(atol) + ".\n", file=f)
                    if save_log and close_file_log_upon_exit: file_log.close()
                    return P
            P_old = np.ndarray.copy(P)
            # Increase the number of slabs approximately by growth_factor_n_slabs.  Do it only
            # if the slab edges have *not* been explicitly provided by the user in t_slab_edges.
            if t_slab_edges_original is None:
                ran_with_max_n_slabs = False if n_slabs < max_n_slabs else True
                n_slabs_old = n_slabs
                n_slabs = min(round(growth_factor_n_slabs*n_slabs), max_n_slabs)
                # Occasionally, the new number of slabs could be equal to the old number (i.e., if
                # growth_factor_n_slabs is too small or if n_slabs = 1).  If this happens, increase
                # the new number of slabs by 1.
                if ((growth_factor_n_slabs > 1.0) and (n_slabs < max_n_slabs) and \
                    (n_slabs == n_slabs_old)): n_slabs += 1
            # Increase the number of points per slab approximately by growth_factor_n_tpts_per_slab
            ran_with_max_n_tpts_per_slab = False if n_tpts_per_slab < max_n_tpts_per_slab else True
            n_tpts_per_slab_old = n_tpts_per_slab
            n_tpts_per_slab = min(int(growth_factor_n_tpts_per_slab*n_tpts_per_slab), 
                max_n_tpts_per_slab)
            if ((growth_factor_n_tpts_per_slab > 1.0) and \
                (n_tpts_per_slab < max_n_tpts_per_slab) and \
                (n_tpts_per_slab == n_tpts_per_slab_old)): n_tpts_per_slab += 1
            loop_count += 1


def _avg_prob_dispatch(
    htot: Callable,
    htot_is_function_only_of_energy: bool,
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    average: bool,
    source_func_name: str,
    smooth_profile: Optional[bool] = True,
    engine_kwargs: Optional[dict] = None
):
    r"""Phase-averaged probabilities, for the position-independent Hamiltonians.

    Returns ``NotImplemented`` when ``average`` is falsy, so a caller can place this ahead of
    its ordinary dispatch chain and fall through untouched in the default case.

    Averaging is exact here and costs one eigendecomposition per energy: with the Hamiltonian
    independent of position, the evolution is a fixed set of phases whose averages are known in
    closed form (see :mod:`magnus.avgprob`).  Which pairs of eigenvalues have actually averaged
    is decided from the baseline rather than assumed, so a request made where the oscillation
    has not decohered is warned about instead of being answered with an expression that does not
    describe it.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    htot : Callable
        Total Hamiltonian as a function of energy alone [eV].
    htot_is_function_only_of_energy : bool
        Whether ``htot`` is independent of position.  Averaging a position-dependent
        Hamiltonian needs the adiabatic treatment, which this does not yet implement.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies [eV].
    L : int, float, list, or np.ndarray
        Baseline(s) [:math:`\text{eV}^{-1}`], used to decide which eigenvalue pairs have
        decohered.
    nu_i, nu_f : int or None
        Initial and final flavor; if both are given, the single probability is returned.
    average : bool
        Whether the caller asked for the averaged probability.
    source_func_name : str
        Name of the calling function, for error messages.

    Returns
    -------
    np.ndarray, float, or NotImplemented
        The averaged probabilities, shaped as the caller's ordinary return value, or
        ``NotImplemented`` if ``average`` is falsy.
    """
    if not average:
        return NotImplemented

    sample_numerically = (not htot_is_function_only_of_energy) and (not smooth_profile)
    if sample_numerically and (engine_kwargs is None):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": average=True "
            "needs a Hamiltonian that is either constant along the trajectory or smooth enough "
            "to have an instantaneous eigenbasis, and this caller cannot fall back to sampling.")

    energy_arr, L_arr, return_float, ok = _normalize_energy_L(energy, L)
    if not ok:
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": the energy "
            "and L arrays must have the same length, or one of them must be a single value.")
    if len(energy_arr) == 1 and len(L_arr) > 1:
        energy_arr = np.repeat(energy_arr, len(L_arr))
    if len(L_arr) == 1 and len(energy_arr) > 1:
        L_arr = np.repeat(L_arr, len(energy_arr))

    n_pts = len(energy_arr)
    undecided_points = 0
    uncertified_points = 0

    if htot_is_function_only_of_energy:
        # Constant along the trajectory: the averaged limit is closed-form, one
        # eigendecomposition per energy.
        H = np.stack([np.asarray(htot(float(enu)), dtype=complex) for enu in energy_arr])
        eigenvalues, eigenvectors = np.linalg.eigh(H)
        d = H.shape[-1]
        P_out = np.empty((n_pts, d, d))
        for i in range(n_pts):
            blocks, undecided = avgprob.coherence_report(eigenvalues[i], float(L_arr[i]))
            if undecided: undecided_points += 1
            P_out[i] = avgprob.averaged_probabilities_from_eigenbasis(eigenvectors[i],
                blocks=blocks)
    elif sample_numerically:
        # No closed form: the profile steps through discontinuities (PREM layer boundaries),
        # so there is no instantaneous eigenbasis to decohere in.  The probability is instead
        # propagated for real across an energy window and averaged over it.
        #
        # This is a different quantity from the two branches above.  They return the exact
        # L/E -> infinity limit, which needs no window; this returns the average over one
        # particular window, and the answer depends on its width.  That is why it is warned
        # about below, and why the width is a named constant rather than a literal.
        d = np.asarray(htot(float(energy_arr[0]), float(L0)), dtype=complex).shape[-1]
        P_out = np.empty((n_pts, d, d))
        eng = dict(engine_kwargs)
        extra = eng.pop('kwargs', None) or {}
        worst_sem = 0.0
        for i in range(n_pts):
            L_i = float(L_arr[i])

            def prob_of_energy(enu, L_i=L_i):
                return osc_prob_energy_baseline(htot, enu, L_i, L0, None, None,
                    htot_is_function_only_of_energy, **eng, **extra)

            P_out[i], sem = avgprob.averaged_probabilities_numerically(prob_of_energy,
                float(energy_arr[i]))
            worst_sem = max(worst_sem, sem)

        warnings.warn(gd.WARNING_MSG_NO_COLOR + " oscprob." + source_func_name + ": average=True "
            "on a profile with discontinuities has no closed form, so the probability was "
            "propagated across an energy window of +/-" +
            str(100.0*avgprob.AVG_DEFAULT_ENERGY_SPREAD) + "% and averaged over " +
            str(avgprob.AVG_DEFAULT_N_SAMPLES) + " samples.  This is the average over that "
            "window, not the L/E -> infinity limit, and it depends on the window: the largest "
            "standard error of the mean here is " + format(worst_sem, '.2e') + ".  Pass an "
            "explicit width via magnus.avgprob.averaged_probabilities_numerically if the "
            "measurement has a known resolution.  Shown once per session.",
            PhaseAveragingWarning, stacklevel=3)

    else:
        # Position-dependent and smooth: decohere in the eigenbasis at production, transport
        # along the levels of the instantaneous Hamiltonian (with the exact crossing
        # probabilities wherever the evolution stops being adiabatic), and read out in the
        # eigenbasis at detection.  See magnus.avgprob.averaged_probabilities_adiabatic.
        d = np.asarray(htot(float(energy_arr[0]), float(L0)), dtype=complex).shape[-1]
        P_out = np.empty((n_pts, d, d))
        for i in range(n_pts):
            enu = float(energy_arr[i])

            def H_of_l(l, enu=enu):
                return htot(enu, l)

            P_out[i], report = avgprob.averaged_probabilities_adiabatic(H_of_l, float(L0),
                float(L_arr[i]))
            if report['undecided'] or report['undecided_between_crossings']:
                undecided_points += 1
            if not report['patches_converged']:
                uncertified_points += 1

    if uncertified_points > 0:
        warnings.warn(gd.WARNING_MSG_NO_COLOR + " oscprob." + source_func_name + ": the local "
            "Magnus patch across a non-adiabatic crossing did not converge at " +
            str(uncertified_points) + " of " + str(n_pts) + " (energy, L) point(s), so the "
            "level-crossing probabilities there are not trustworthy.  Shown once per session.",
            HybridCertificationWarning, stacklevel=3)

    if undecided_points > 0:
        warnings.warn(gd.WARNING_MSG_NO_COLOR + " oscprob." + source_func_name + ": the averaged "
            "probability was requested at " + str(undecided_points) + " of " +
            str(len(energy_arr)) + " (energy, L) point(s) where at least one pair of eigenvalues "
            "has neither decohered nor stayed coherent, so no averaged expression describes it.  "
            "The oscillation probability itself (average=False) is the meaningful quantity there. "
            "Shown once per session.", PhaseAveragingWarning, stacklevel=3)

    _note_engine('average')
    if (nu_i is not None) and (nu_f is not None):
        P_out = P_out[:, nu_i, nu_f]

    return P_out.__getitem__(0 if return_float else slice(None))


PARAMETER_SET_METADATA_KEYS = ('name', 'description')
r"""tuple: Module-level constant

Keys carried by the entries of :data:`magnus.globaldefs.OSC_PARAMS_PREDEFINED` that
label the parameter set rather than parameterize the physics.

They are the reason ``**OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']`` cannot be
splatted into a probability function: the two strings travel down the shared
``**kwargs`` chain and are rejected at the far end by
``magnus_expansion_multislab``, whose complaint names neither the caller nor the
parameter set it came from.  :func:`magnus.globaldefs.load_nufit_params` returns
the same numbers without them.

.. versionadded:: 1.0.0
"""


def _reject_parameter_set_metadata(kwargs: dict, source_func_name: str) -> None:
    r"""Rejects the label keys of a predefined parameter set, with the remedy.

    Splatting a whole ``OSC_PARAMS_PREDEFINED`` entry is the natural thing to
    write and it does not work: ``name`` and ``description`` are not oscillation
    parameters, so they flow through every ``**kwargs`` hop until the Magnus core
    raises ``TypeError: magnus_expansion_multislab() got an unexpected keyword
    argument 'name'`` -- a message that points at the one function in the chain
    that has nothing to do with the mistake.

    Caught here instead, where the caller and the fix can both be named.

    .. versionadded:: 1.0.0
    """
    found = [key for key in PARAMETER_SET_METADATA_KEYS if key in kwargs]
    if not found:
        return

    raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": received " +
        str(found) + ", which label a predefined parameter set rather than parameterize the "
        "physics.  This happens when an entry of globaldefs.OSC_PARAMS_PREDEFINED is passed "
        "whole, as **OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'].  Use "
        "globaldefs.load_nufit_params(...), which returns the same numbers without the labels, "
        "or drop the labelling keys.")


def _normalize_energy_L(
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, bool, bool]:
    r"""Normalize energy and L to same-length 1D arrays.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).

    Returns
    -------
    (np.ndarray, np.ndarray, bool, bool)
        ``(energy, L, return_float, ok)``: ``energy`` and ``L`` broadcast to the same length;
        ``return_float`` records whether both inputs were scalars (so that the caller returns a
        scalar-like result); ``ok`` records whether the input lengths were compatible (equal, or
        one of them of length 1).
    """
    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L
    return_float = isinstance(energy, float) and isinstance(L, float)
    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)
    L = np.array([L]) if isinstance(L, float) else np.array(L)
    ok = ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or
          (len(energy) > 1 and len(L) == 1))
    if ok:
        energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
        L = np.full(len(energy), L[0]) if (len(L) == 1) else L
    return energy, L, return_float, ok


#-----------------------------------------------------------------------
# Which engine answered: the shared instrument behind strategy_info and
# cross_check_strategies
#-----------------------------------------------------------------------

ENGINE_FAMILIES = {
    'hybrid': 'adiabatic',
    'ip_exp': 'interaction-picture',
    'magnus': 'magnus-ladder',
    'cumulative': 'magnus-ladder',
    'separable': 'magnus-ladder',
    'average': 'phase-average',
    'expm': 'exact',
}
r"""dict: Module-level constant

Which engines share machinery, and therefore which pairwise comparisons in
:func:`cross_check_strategies` carry information.  Two engines in the **same** family can be
wrong in the same way at the same time, so their agreement is not evidence; two in different
families fail for different reasons.

* ``'adiabatic'`` -- :func:`magnus.adiabatic.hybrid_propagator`: transport in the instantaneous
  eigenbasis, with Magnus patches only inside non-adiabatic windows.  Its blind spots are the
  detector's (a feature narrower than the probe grid, a profile the resolution test rejects).
* ``'magnus-ladder'`` -- the general per-point path, the cumulative baseline scan, and the
  energy-batched separable scan.  All three walk slabs with
  :func:`magnus.magnus.magnus_expansion_multislab`, and the cumulative scan additionally *sizes*
  its grid from an ordinary adaptive :func:`osc_prob` probe, so it inherits that path's stopping
  rule as well.  Grouping them is deliberate: the accuracy step at
  :data:`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` shows they are not interchangeable, but a
  quadrature that cannot see a feature will not see it on any of the three grids.
* ``'interaction-picture'`` -- the 2-flavour exponential-profile fast path.  It uses the same
  Magnus core, but factors the fast vacuum phase out analytically first, so what it must resolve
  is a different function; it is kept separate for that reason and not because the core differs.
* ``'exact'`` -- ``scipy.linalg.expm``, used only where it is the exact answer rather than an
  approximation (see :func:`cross_check_strategies`).
* ``'phase-average'`` -- :mod:`magnus.avgprob`'s closed form, which answers a different question
  and is never compared against the others.

.. versionadded:: 1.0.0
"""


_ENGINE_TRACE = None
r"""list or None: set to a list by ``_engine_probe`` while a diagnostic is watching; ``None``
(and therefore free) on every ordinary call."""


_ENGINES_DISABLED = frozenset()
r"""frozenset: engine labels that the dispatchers must decline, set by ``_engine_probe``.  Used
only by :func:`cross_check_strategies`, to reach an engine that a faster one would otherwise
answer for -- there is no user-facing way to ask for the general ladder specifically when the
interaction-picture path applies, and a cross-check that silently compared the same engine with
itself would be exactly the failure it exists to detect."""


def _scan_for_hidden_features(profile, l0, L, t_breakpoints=None) -> Optional[Dict]:
    r"""Run the sub-probe feature scan once for a whole call, and warn if it finds something.

    Placed at the entry points rather than inside any one engine because the blind spot is not
    one engine's: the hybrid probe grid, the general ladder's slabs and the cumulative scan's
    accuracy grid all miss the same feature.  See :class:`HiddenFeatureWarning`.

    Skipped when the caller supplied ``t_breakpoints`` (they have already said where the
    structure is) and when the profile is not a function of position (nothing to hide in).

    .. versionadded:: 1.0.0

    Returns
    -------
    dict or None
        The scan result from :func:`magnus.adiabatic.find_hidden_features`, or None if the scan
        did not apply.
    """
    if not isinstance(profile, Callable):
        return None
    if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
        return None
    L_arr = np.atleast_1d(np.asarray(L, dtype=float))
    l1 = float(np.max(L_arr))
    if l1 == float(l0):
        return None

    scan = adiabatic.find_hidden_features(profile, float(l0), l1)
    if not scan['hidden']:
        return scan

    # The position IS the actionable content here, so it goes in the message even though that
    # costs the static-message dedup other warnings keep.  It does not cost much: the scan
    # depends on the profile and the interval only, so every call over the same profile
    # produces the same string and Python's default filter still collapses them.
    pad = 3.0*(scan['l_hi'] - scan['l_lo'])
    warnings.warn(
        "osc_prob: the density profile has structure too narrow for any grid this package "
        "samples on, near l = " + format(scan['l_centre'], '.6e') + ". Every engine misses it "
        "-- the adiabatic probe grid, the general Magnus slabs and the cumulative scan alike -- "
        "so they agree with each other and may be wrong together, and no choice of strategy or "
        "tolerance helps. Measured on such a profile: wrong by 2.9e-02 against a requested "
        "1e-3. Pass t_breakpoints=["
        + format(scan['l_lo'] - pad, '.6e') + ", " + format(scan['l_centre'], '.6e') + ", "
        + format(scan['l_hi'] + pad, '.6e') + "] to put slab edges on it; that is a partial "
        "cure (measured 3.9e-03 -> 8.5e-05), not a complete one. If the narrow structure is an "
        "artefact of how the profile function was written rather than physics, this is safe to "
        "ignore. Shown once per profile per session.",
        HiddenFeatureWarning, stacklevel=3)
    return scan


def _note_engine(label: str, answered: bool = True, **detail) -> None:
    r"""Record that ``label`` answered (or declined), if anything is watching.

    One dict per dispatch decision, in the order the decisions were taken, so the trace reads as
    the route the request took: ``hybrid declined (uncertified) -> magnus answered``.

    .. versionadded:: 1.0.0
    """
    if _ENGINE_TRACE is not None:
        _ENGINE_TRACE.append(dict(engine=label, answered=answered, **detail))


@contextmanager
def _engine_probe(disabled=(), info=None, extra=None):
    r"""Watch which engine answers, and optionally forbid some of them.

    Restores both globals on the way out, including on an exception, so a raising call (which
    ``cumulative=True`` does by design when it cannot serve a request) cannot leave a dispatcher
    disabled for the rest of the session.  ``info``, if given, is the caller's ``strategy_info``
    dict and is filled on the way out -- also on an exception, since "which engine was I in when
    this raised" is exactly what a caller debugging one wants.

    .. versionadded:: 1.0.0
    """
    global _ENGINE_TRACE, _ENGINES_DISABLED
    prev_trace, prev_disabled = _ENGINE_TRACE, _ENGINES_DISABLED
    # Nested probes SHARE one trace, and a nested one can only add to the disabled set.  Both
    # matter because nesting is the normal case, not an edge case: cross_check_strategies
    # watches from outside the wrapper, and the wrapper opens its own probe for strategy_info.
    # A fresh list at the inner level collected every note and left the outer one empty, so the
    # cross-check reported that no engine had answered at all -- and re-assigning rather than
    # unioning the disabled set would have let the inner probe re-enable an engine the outer one
    # had switched off, which is how the cross-check reaches an engine a faster one shadows.
    trace = prev_trace if prev_trace is not None else []
    start = len(trace)
    _ENGINE_TRACE = trace
    _ENGINES_DISABLED = prev_disabled | frozenset(disabled)
    try:
        yield trace
    finally:
        _ENGINE_TRACE, _ENGINES_DISABLED = prev_trace, prev_disabled
        if info is not None:
            info.update(_summarize_engine_trace(trace[start:]))
            if extra is not None:
                info.update(extra)


def _summarize_engine_trace(trace) -> Dict:
    r"""Turn a raw trace into the ``strategy_info`` an ordinary caller wants.

    .. versionadded:: 1.0.0
    """
    answered = [e for e in trace if e['answered']]
    used = answered[-1] if answered else None
    return {
        'engine': used['engine'] if used else None,
        'family': ENGINE_FAMILIES.get(used['engine']) if used else None,
        'certified': used.get('certified') if used else None,
        'declined': [(e['engine'], e.get('reason', 'does not apply'))
                     for e in trace if not e['answered']],
        'trace': [{k: v for k, v in e.items() if not k.startswith('_')} for e in trace],
    }


def _osc_prob_scan_separable(
    H_E: np.ndarray,
    VCC_func: Callable,
    h_matt: np.ndarray,
    L0: float,
    L_val: float,
    t_breakpoints: Optional[np.ndarray],
    magnus_exp_order: int,
    integration_method: str,
    rtol: Optional[float],
    atol: Optional[float],
    growth_factor_n_slabs: float,
    growth_factor_n_tpts_per_slab: float,
    max_num_loops: int,
    min_n_slabs: int,
    max_n_slabs: int,
    min_n_tpts_per_slab: int,
    max_n_tpts_per_slab: int,
    n_slabs: int,
    n_tpts_per_slab: int
) -> np.ndarray:
    r"""Energy-batched probability scan for separable Hamiltonians.

    Computes the probabilities of many neutrino energies that share the same
    baseline [``L0``, ``L_val``] in one batched pipeline, for Hamiltonians of
    the separable form

        H(E, l) = H_E(E) + VCC(l) * h_matt ,

    where ``H_E`` (shape (nE, d, d)) collects all the position-independent,
    energy-dependent terms (vacuum, LIV, ...), ``VCC_func`` is the scalar
    matter potential along the trajectory, and ``h_matt`` (shape (d, d)) is
    the constant matter matrix it multiplies.  The position samples of the
    potential are computed once per refinement level and shared by all
    energies, and the Magnus kernel (quadrature, commutators, exponentials,
    slab products) runs with the energy axis batched in front of the slab
    axis.

    The adaptive refinement mirrors :func:`osc_prob`: the slab count (and,
    for the quadrature methods, the points per slab) grows geometrically
    until the probabilities of each energy agree between successive levels
    within (rtol, atol); converged energies drop out of the batch.  Energies
    are processed in chunks to bound the memory of the sample array.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_E : np.ndarray
        Position-independent, energy-dependent part of the Hamiltonian for each energy, shape
        (nE, d, d) (vacuum, LIV, ...).
    VCC_func : Callable
        Scalar matter potential along the trajectory, as a function of position (accepts an array).
    h_matt : np.ndarray
        Constant matrix multiplying ``VCC_func(l)``, shape (d, d).
    L0 : float
        Initial position.
    L_val : float
        Final position (baseline).
    t_breakpoints : np.ndarray, optional
        Mandatory slab edges (e.g., PREM layer boundaries) inserted into the grid at every
        refinement level.
    magnus_exp_order : int
        Highest order of the Magnus expansion.
    integration_method : str
        'gl', 'trapezoid', or 'simpson'.
    rtol, atol : float, optional
        Target relative/absolute tolerance between successive refinement levels. If both None,
        run once with the given fixed ``n_slabs``/``n_tpts_per_slab``.
    growth_factor_n_slabs, growth_factor_n_tpts_per_slab : float
        Factors by which ``n_slabs``/``n_tpts_per_slab`` are multiplied on each refinement loop.
    max_num_loops : int
        Maximum number of refinement loops.
    min_n_slabs, max_n_slabs : int
        Bounds on the number of slabs.  ``max_n_slabs=None`` selects the per-method cap;
        see :data:`MAX_N_SLABS_DEFAULT`.
    min_n_tpts_per_slab, max_n_tpts_per_slab : int
        Bounds on the number of time points per slab.
    n_slabs, n_tpts_per_slab : int
        Starting number of slabs/time points per slab.  Under a tolerance, ``n_slabs`` is a
        floor on the refinement ladder rather than a discarded argument; see :func:`osc_prob`.

    Returns
    -------
    np.ndarray
        Stacked probability matrices, shape (nE, d, d).
    """
    # None means 'use the cap appropriate to this integration method'
    # (see MAX_N_SLABS_DEFAULT); an explicit value always wins.
    max_n_slabs = _resolve_max_n_slabs(max_n_slabs, integration_method)
    nE, dim = H_E.shape[0], H_E.shape[-1]
    tol_requested = ((rtol is not None) and (atol is not None))

    if integration_method == 'gl':
        # The accuracy of the GL method is controlled by n_slabs only
        growth_factor_n_tpts_per_slab = 1.0
        min_n_tpts_per_slab = max_n_tpts_per_slab = n_tpts_per_slab = 2
        s_nodes = magnus.gl_nodes(magnus_exp_order)

    if tol_requested:
        # The caller's n_slabs is a floor on the refinement ladder, not something to discard; see
        # the corresponding note in osc_prob.
        min_n_slabs = int(min(max(min_n_slabs, n_slabs), max_n_slabs))
        n_tpts_per_slab = min_n_tpts_per_slab
        # Physics-informed starting number of slabs (see magnus.suggest_n_slabs):
        # integral of the traceless Hamiltonian over the trajectory, maximized
        # over the energies of the scan
        if integration_method == 'gl':
            ts = np.linspace(L0, L_val, 17)
            V17 = np.asarray(VCC_func(ts))
            I_V = (np.sum(V17) - 0.5*(V17[0] + V17[-1]))*(L_val - L0)/16.0
            M = (L_val - L0)*H_E + I_V*h_matt
            M = M - (np.trace(M, axis1=-2, axis2=-1)/dim)[:, None, None]*np.eye(dim)
            try:
                phase = np.max(np.linalg.svd(M, compute_uv=False))
            except np.linalg.LinAlgError:
                phase = 0.0
            n_slabs = int(np.clip(max(min_n_slabs,
                np.ceil(phase/(2.0*np.pi))), 1, max_n_slabs))
        else:
            n_slabs = min_n_slabs

    P_prev = np.full((nE, dim, dim), np.nan)
    P_out = np.empty((nE, dim, dim))
    active = np.arange(nE)
    mA = -1j*h_matt.astype(complex)
    HE_c = -1j*H_E.astype(complex)

    loop_count = 1
    while True:
        # Slab grid shared by all energies (PREM-layer breakpoints included)
        grid = np.linspace(L0, L_val, n_slabs + 1)
        if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
            bp = np.atleast_1d(np.asarray(t_breakpoints, dtype=float))
            bp = bp[(bp > L0) & (bp < L_val)]
            grid = np.unique(np.concatenate([grid, bp]))
        edges = np.column_stack([grid[:-1], grid[1:]])
        widths = edges[:, 1] - edges[:, 0]

        if integration_method == 'gl':
            s = s_nodes
        else:
            s = np.linspace(0.0, 1.0, n_tpts_per_slab)
        tgrid = edges[:, :1] + widths[:, None]*s              # (n_slabs, m)
        V = np.asarray(VCC_func(tgrid.ravel())).reshape(tgrid.shape)
        Vmat = V[:, :, None, None]*mA                         # (n_slabs, m, d, d)

        # Batched kernel over the active energies, chunked so that each
        # sample array At holds at most ~4M complex entries (~64 MB)
        chunk, _ = _tile_for_working_set(len(active), 1, tgrid.size*dim*dim)
        P_new = np.empty((len(active), dim, dim))
        for i0 in range(0, len(active), chunk):
            sel = active[i0:i0+chunk]
            At = HE_c[sel][:, None, None, :, :] + Vmat[None, :, :, :, :]
            U = magnus.evolution_operators_from_samples(At, widths,
                magnus_exp_order, integration_method, validate_input=False)
            Utot = U[:, -1]
            for k in range(U.shape[1] - 2, -1, -1):
                Utot = Utot @ U[:, k]
            P_new[i0:i0+chunk] = np.swapaxes(
                Utot.real**2 + Utot.imag**2, -1, -2)

        if not tol_requested:
            P_out[active] = P_new
            return P_out

        prev = P_prev[active]
        have_prev = ~np.isnan(prev[:, 0, 0])
        conv = have_prev & np.all(np.abs(P_new - prev) <= atol + rtol*np.abs(prev),
                                  axis=(-1, -2))
        P_out[active[conv]] = P_new[conv]
        P_prev[active] = P_new
        active = active[~conv]
        if active.size == 0:
            return P_out

        at_caps = ((n_slabs >= max_n_slabs) and
                   (n_tpts_per_slab >= max_n_tpts_per_slab))
        if (loop_count >= max_num_loops) or at_caps:
            warnings.warn("osc_prob (energy-batched scan): requested tolerance "
                "not achieved for some energies (refinement caps reached); the "
                "returned probabilities may be inaccurate. Try increasing "
                "max_n_slabs, max_n_tpts_per_slab, or max_num_loops. Shown "
                "once per session.", ToleranceNotAchievedWarning, stacklevel=2)
            P_out[active] = P_new[~conv]
            return P_out

        n_slabs_old = n_slabs
        n_slabs = min(round(growth_factor_n_slabs*n_slabs), max_n_slabs)
        if ((growth_factor_n_slabs > 1.0) and (n_slabs < max_n_slabs) and
                (n_slabs == n_slabs_old)):
            n_slabs += 1
        n_tpts_old = n_tpts_per_slab
        n_tpts_per_slab = min(int(growth_factor_n_tpts_per_slab*n_tpts_per_slab),
                              max_n_tpts_per_slab)
        if ((growth_factor_n_tpts_per_slab > 1.0) and
                (n_tpts_per_slab < max_n_tpts_per_slab) and
                (n_tpts_per_slab == n_tpts_old)):
            n_tpts_per_slab += 1
        loop_count += 1


def _osc_prob_scan_separable_dispatch(
    h_vac_energy_indep: np.ndarray,
    VCC_func: Union[Callable, float],
    h_matt: np.ndarray,
    h_liv_energy_indep: Optional[np.ndarray],
    n_liv: Optional[Union[int, float]],
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    scan_kwargs: Dict
):
    r"""Decide whether the energy-batched scan engine applies; run it if so.

    Returns NotImplemented when the request does not fit the engine (single
    point, per-point baselines, user-provided slab edges, parallel or logged
    runs, iteration over the expansion order, or unknown extra arguments), in
    which case the caller falls back to the generic per-point path.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    h_vac_energy_indep : np.ndarray
        Energy-independent part of the vacuum Hamiltonian.
    VCC_func : Callable or float
        Matter potential, as a function of position (required for the batched engine to apply;
        a constant potential falls back to the generic path).
    h_matt : np.ndarray
        Constant matrix multiplying ``VCC_func(l)``.
    h_liv_energy_indep : np.ndarray, optional
        Energy-independent part of the LIV Hamiltonian, if any.
    n_liv : int or float, optional
        Power of the energy dependence of the LIV operator, if ``h_liv_energy_indep`` is given.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s); the batched engine applies only when all requested baselines are equal.
    L0 : int or float
        Initial position.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    scan_kwargs : dict
        The refinement/logging keyword arguments of :func:`osc_prob_energy_baseline` (rtol, atol,
        magnus_exp_order, integration_method, growth factors, loop/slab/time-point bounds,
        t_slab_edges, n_jobs, save_log, file_log) plus a nested
        'kwargs' dict of any remaining, unrecognized keyword arguments.

    Returns
    -------
    np.ndarray or NotImplemented
        The oscillation probability (or single channel), computed via the batched engine; or the
        ``NotImplemented`` singleton if the request does not fit it.
    """
    if ('separable' in _ENGINES_DISABLED) or (scan_kwargs.get('cumulative') is True):
        return NotImplemented
    kwargs = dict(scan_kwargs.get('kwargs', {}))
    t_breakpoints = kwargs.pop('t_breakpoints', None)
    n_slabs = kwargs.pop('n_slabs', 1)
    n_tpts_per_slab = kwargs.pop('n_tpts_per_slab', 100)
    if len(kwargs) > 0:
        return NotImplemented
    if not isinstance(VCC_func, Callable):
        return NotImplemented
    if scan_kwargs['t_slab_edges'] is not None:
        return NotImplemented
    if (scan_kwargs['n_jobs'] != 1) or scan_kwargs['save_log'] or \
            (scan_kwargs['file_log'] is not None):
        return NotImplemented

    energy_arr, L_arr, return_float, ok = _normalize_energy_L(energy, L)
    if (not ok) or (len(energy_arr) < 2) or (not np.all(L_arr == L_arr[0])):
        return NotImplemented

    rtol, atol = scan_kwargs['rtol'], scan_kwargs['atol']
    if (rtol is None) != (atol is None):
        rtol = 0.0 if rtol is None else rtol
        atol = 0.0 if atol is None else atol

    # All the position-independent, energy-dependent terms of the Hamiltonian
    H_E = (1.0/energy_arr)[:, None, None]*np.asarray(h_vac_energy_indep)
    if h_liv_energy_indep is not None:
        H_E = H_E + (energy_arr**n_liv)[:, None, None]*np.asarray(h_liv_energy_indep)

    P = _osc_prob_scan_separable(H_E, VCC_func, np.asarray(h_matt), float(L0),
        float(L_arr[0]), t_breakpoints, scan_kwargs['magnus_exp_order'],
        scan_kwargs['integration_method'], rtol, atol,
        scan_kwargs['growth_factor_n_slabs'],
        scan_kwargs['growth_factor_n_tpts_per_slab'],
        scan_kwargs['max_num_loops'], scan_kwargs['min_n_slabs'],
        scan_kwargs['max_n_slabs'], scan_kwargs['min_n_tpts_per_slab'],
        scan_kwargs['max_n_tpts_per_slab'], n_slabs, n_tpts_per_slab)

    _note_engine('separable')
    if (nu_i is not None) and (nu_f is not None):
        P = P[:, nu_i, nu_f]
    return P.__getitem__(0 if return_float else slice(None))


def _osc_prob_ip_exp_core(
    H_E: np.ndarray,
    l_scale: float,
    VCC_func: Callable,
    h_matt: np.ndarray,
    L0: float,
    L_val: float,
    rtol: Optional[float],
    atol: Optional[float],
    growth_factor_n_slabs: float,
    max_num_loops: int,
    min_n_slabs: int,
    max_n_slabs: int,
    n_slabs: int
) -> Tuple[np.ndarray, bool]:
    r"""Interaction-picture Magnus integrator for an exponential matter profile.

    Computes the evolution operator for Hamiltonians of the separable form
    ``H(E, l) = H_E(E) + VCC_func(l) * h_matt``, with ``VCC_func(l) = VCC_func(0) * exp(-l/l_scale)``
    a genuine exponential profile, WITHOUT resolving the (possibly huge, at low energy) fast phase of
    ``H_E`` slab by slab. ``H_E`` is diagonalized once (it does not depend on position); in that
    eigenbasis, the free ("vacuum") evolution ``exp(-i H_E s)`` is an exactly known diagonal phase for
    any ``s``, so it is factored out analytically instead of being resolved by narrow slabs, leaving
    only the matter-potential envelope to be integrated. Within each slab ``[l0, l0+h]``, the
    first-order Magnus term of the resulting interaction-picture generator,

    .. math::

       \Omega_1(h) = -i \int_0^h e^{i H_E s}\, V_\text{CC}(l_0+s)\,
       h_\text{matt}\, e^{-i H_E s}\, ds ,

    has a closed form because ``VCC_func(l0+s) = VCC_func(l0) exp(-s/l_scale)`` is itself an
    exponential: in the eigenbasis of ``H_E`` (eigenvalues :math:`\lambda_j`, so
    :math:`\Delta_{jk} = \lambda_j - \lambda_k`),

    .. math::

       \left(\tilde{\Omega}_1\right)_{jk}(h) = -i \left(\tilde{h}_\text{matt}\right)_{jk}
       V_\text{CC}(l_0)\,
       \frac{e^{\left(i \Delta_{jk} - 1/l_\text{scale}\right) h} - 1}
            {i \Delta_{jk} - 1/l_\text{scale}} ,

    valid uniformly for :math:`j = k` too (the :math:`j = k` denominator,
    :math:`-1/l_\text{scale}`, is never zero for finite
    ``l_scale``). This is exact in the envelope (no local-constant or local-linear approximation of
    ``VCC_func`` is made); the only approximation is truncating the interaction-picture Magnus series
    at first order, which is accurate away from an MSW resonance (where the matter term becomes
    comparable to the *vacuum* splitting :math:`\Delta_{jk}`, rather than merely small) and improves,
    rather than worsens, at lower neutrino energy (:math:`\Delta_{jk}` grows as :math:`1/E`). Both
    factors (``exp(-i H_E h)`` and :math:`\exp(\Omega_1)`) are exactly unitary by construction (the
    former is a diagonal phase, the latter is exponentiated via ``_expm_stack`` from an
    anti-Hermitian generator), so the returned probabilities remain exactly unitary regardless of how
    good the approximation is. Accuracy is controlled the usual way: growing ``n_slabs`` shrinks the
    per-slab truncation error (which vanishes faster than the slab width itself), so successive
    refinements converge to the exact solution; the loop mirrors :func:`osc_prob`'s own
    successive-refinement comparison, batched over the leading energy axis of ``H_E``.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_E : np.ndarray
        Position-independent, energy-dependent part of the Hamiltonian, shape (nE, d, d).
    l_scale : float
        Length scale of the exponential density decrease.
    VCC_func : Callable
        Matter potential along the trajectory (accepts an array of positions), satisfying
        ``VCC_func(l) = VCC_func(0)*exp(-l/l_scale)``.
    h_matt : np.ndarray
        Constant matrix multiplying ``VCC_func(l)``, shape (d, d).
    L0 : float
        Initial position.
    L_val : float
        Final position (baseline).
    rtol, atol : float, optional
        Target relative/absolute tolerance between successive refinement levels. If both None, run
        once with the given fixed ``n_slabs``.
    growth_factor_n_slabs : float
        Factor by which ``n_slabs`` is multiplied on each refinement loop.
    max_num_loops : int
        Maximum number of refinement loops.
    min_n_slabs, max_n_slabs : int
        Bounds on the number of slabs.  ``max_n_slabs=None`` selects the 'gl' entry of
        :data:`MAX_N_SLABS_DEFAULT`, this integrator being closed-form per slab and so
        comparably cheap.
    n_slabs : int
        Starting number of slabs (or fixed count, if no tolerance is requested).

    Returns
    -------
    (np.ndarray, bool)
        Stacked probability matrices, shape (nE, d, d), and whether the requested tolerance was
        achieved (always True if no tolerance was requested). If False, the result should be
        discarded in favor of the general (slower, but unconditionally convergent) slab-refinement
        method; this happens when the matter term is not a small perturbation on the vacuum
        splitting anywhere along the trajectory (e.g., at an MSW resonance).
    """
    # This integrator has no `integration_method`: within each slab it is closed-form, with
    # no quadrature at all, so its cost per slab is comparable to 'gl' rather than to the
    # cumulative-quadrature methods.  It therefore takes the 'gl' cap when none is given.
    if max_n_slabs is None:
        max_n_slabs = MAX_N_SLABS_DEFAULT['gl']
    tol_requested = (rtol is not None) and (atol is not None)

    # Diagonalize the position-independent part of H once; this eigenbasis is shared by every slab,
    # so the fast phase of H_E is handled exactly, analytically, regardless of slab width.
    Lam, W = np.linalg.eigh(H_E)                                    # (nE, d), (nE, d, d)
    Wd = np.conj(np.swapaxes(W, -1, -2))
    Mt = Wd @ h_matt.astype(complex)[None, :, :] @ W                # h_matt in the H_E eigenbasis

    Delta = Lam[:, :, None] - Lam[:, None, :]                       # (nE, d, d)
    denom = 1j*Delta - (1.0/l_scale)                                # never zero (l_scale finite)

    if tol_requested:
        # The caller's n_slabs is a floor on the refinement ladder, not something to discard; see
        # the corresponding note in osc_prob.  Clipped at this method's own ceiling, not at
        # max_n_slabs: the slab budget here is decoupled from the caller's cap (see the note
        # below), so max_n_slabs is not the bound the growth step will respect.
        n_slabs = int(min(max(min_n_slabs, n_slabs), IP_EXP_N_SLABS_CAP))

    # The part of h_matt that is diagonal in the H_E eigenbasis commutes with H_E and does not
    # oscillate: it accumulates an ordinary, unsuppressed phase (proportional to the matter
    # potential integrated over the *whole* slab) that first-order Omega_1 does not shrink just
    # because Delta_jk is large. So max||Omega_1|| does not fall smoothly slab by slab the way it
    # would for a pure Magnus quadrature error -- for a slab still wide compared to the scale at
    # which this diagonal phase becomes O(1), successive refinements can land on essentially
    # uncorrelated (not just slowly converging) probabilities. The two safeguards below keep this
    # method from ever reporting a false convergence in that pre-asymptotic regime: (a) the
    # successive-refinement comparison is trusted only once max||Omega_1|| itself has dropped
    # below a conservative threshold (comfortably inside the regime where the neglected Omega_2 ~
    # O(||Omega_1||^2) term is genuinely small), and (b) even then, agreement is required twice in
    # a row. The slab budget is also decoupled from the caller's max_n_slabs/growth_factor_n_slabs
    # (calibrated for the much more expensive quadrature-based slabs of the general method): each
    # slab here costs one small eigendecomposition, so pushing to hundreds of thousands of slabs
    # when needed still completes in a couple of seconds.
    # The neglected Omega_2 ~ O(||Omega_1||^2) term sets the probability error at ~C*omega^2 for
    # some O(1) constant C (empirically ~0.005 for the sole off-diagonal pair of a genuine 2-level
    # H_E, which is all _osc_prob_ip_exp_dispatch admits here); tie the trust threshold to the
    # requested tolerance (with a safety factor) instead of a fixed value, so tighter requests
    # correctly demand more slabs rather than risking a plausible-looking but insufficiently
    # accurate "convergence". The slab budget is a fixed ceiling, independent of the caller's
    # max_n_slabs/growth_factor_n_slabs (calibrated for the much more expensive quadrature-based
    # slabs of the general method): each slab here costs one 2x2 eigendecomposition, so even the
    # full ceiling completes in a couple of seconds.
    omega_trust_threshold = min(0.1, np.sqrt((atol + rtol)/2.0)) if tol_requested else 0.1
    n_slabs_cap = IP_EXP_N_SLABS_CAP
    loop_cap = IP_EXP_LOOP_CAP
    growth = 2.0
    nE, dim = H_E.shape[0], H_E.shape[-1]

    # Can the trust gate above ever open, at any slab count this method is allowed?
    #
    # Certification requires max|Omega_t| < omega_trust_threshold, and that maximum is bounded
    # below by the largest *diagonal* entry, which has a closed form here.  On the diagonal
    # Delta_jj = 0, so denom_jj = -1/l_scale and the slab integral collapses to
    # l_scale*(1 - exp(-w/l_scale)) for slab width w, giving
    #
    #     max|Omega_jj| = max_ej|Mt[e,j,j]| * max_s|V(l_s)| * l_scale*(1 - exp(-w/l_scale)),
    #
    # decreasing in the slab count and independent of it otherwise.  If that alone still
    # exceeds the threshold at the ceiling, then max|Omega_t| does too at every reachable slab
    # count, no comparison is ever trusted, and the ladder is guaranteed to climb to the cap
    # and refuse.  Detecting it costs two evaluations of VCC_func and no allocation.
    #
    # This is a bound, not an estimate: it can only ever say "certification is impossible",
    # never "certification will succeed", so it cannot abandon a case that would have worked.
    # The pass below still runs once, at the starting slab count, because an uncertified result
    # is still required to be a genuine unitary probability matrix -- see the tests on this
    # function's give-up exits.
    certifiable = True
    if tol_requested:
        v_max = max(abs(float(np.real(np.asarray(VCC_func(L0))))),
                    abs(float(np.real(np.asarray(VCC_func(L_val))))))
        mt_diag_max = float(np.max(np.abs(np.diagonal(Mt, axis1=-2, axis2=-1))))
        scale = mt_diag_max*v_max
        if scale > 0.0:
            c = omega_trust_threshold/scale
            if c < l_scale:                       # otherwise the bound is below threshold at any w
                w_max = -l_scale*np.log1p(-c/l_scale)
                certifiable = np.ceil((L_val - L0)/w_max) <= n_slabs_cap

    P_prev = None
    n_slabs_prev = None
    consecutive_agreements = 0
    loop_count = 1
    while True:
        # The full edge grid is O(n_slabs) floats and independent of the energy count, so it
        # is affordable at any slab cap; it is built whole rather than per tile because
        # np.linspace's endpoint handling is not reproduced by arithmetic on a sub-range,
        # and the tiling below must not perturb a single slab edge.
        grid = np.linspace(L0, L_val, n_slabs + 1)

        # Everything below is tiled over (energy, slab).  The arrays this replaces --
        # (nE, n_slabs, d, d) complex, several live at once -- are the whole of the memory
        # bug: at the slab ceiling they reached ~1.3 GB per energy, so a batched solar call
        # could exhaust the machine.  See docs/dev/BUG_IP_EXP_MEMORY.md.
        #
        # The tiling is exact, not approximate.  Within a tile the arithmetic is elementwise,
        # so slicing changes no value; and the product is folded slab-by-slab in the same
        # descending order as before, with the accumulator on the left, so the parenthesis
        # nesting -- the only thing that could move a floating-point result -- is unchanged.
        # Blocks are therefore walked from the *last* slab backwards. A test pins the output
        # of a tiled run against an untiled one at exact equality.
        # live_arrays: arg, the exp() temporary, I, Omega_t, _expm_stack's eigenvectors and
        # its workspace, U_slab, and the accumulator's operand -- eight of this shape at the
        # peak, which is why the budget has to be divided rather than applied per array.
        e_chunk, blk = _tile_for_working_set(nE, n_slabs, dim*dim, live_arrays=8)
        Utot = np.empty((nE, dim, dim), dtype=complex)
        max_omega = 0.0

        for e0 in range(0, nE, e_chunk):
            esel = slice(e0, min(e0 + e_chunk, nE))
            acc = None
            for b1 in range(n_slabs, 0, -blk):                      # descending slab blocks
                b0 = max(0, b1 - blk)
                edges0 = grid[b0:b1]                                # (nb,) slab starts
                widths = grid[b0 + 1:b1 + 1] - edges0               # (nb,)
                V0 = np.asarray(VCC_func(edges0), dtype=complex)    # VCC at each slab's start

                arg = denom[esel, None, :, :]*widths[None, :, None, None]
                I = (np.exp(arg) - 1.0)/denom[esel, None, :, :]
                Omega_t = -1j*Mt[esel, None, :, :]*V0[None, :, None, None]*I
                max_omega = max(max_omega, float(np.max(np.abs(Omega_t))))

                U_free_diag = np.exp(-1j*Lam[esel, None, :]*widths[None, :, None])
                U_slab = U_free_diag[..., :, None]*magnus._expm_stack(
                    Omega_t, warn_wide=False)

                for k in range(U_slab.shape[1] - 1, -1, -1):
                    acc = U_slab[:, k] if acc is None else acc @ U_slab[:, k]
                del arg, I, Omega_t, U_free_diag, U_slab
            Utot[esel] = acc

        Utot = W @ Utot @ Wd                                        # back to the flavor basis

        P_new = np.swapaxes(Utot.real**2 + Utot.imag**2, -1, -2)

        if not tol_requested:
            return P_new, True

        # Certification is provably out of reach (see the bound above): refuse now rather
        # than doubling the slab count twenty more times to arrive at the same refusal.  The
        # caller's dispatcher discards this result and falls back to the general method, as
        # it would have anyway -- identically, and roughly a thousand times sooner.
        if not certifiable:
            return P_new, False

        at_cap = (n_slabs >= n_slabs_cap)
        # Once n_slabs is pinned at the cap, growth is a no-op: a further "refinement" would just
        # repeat this identical computation, which trivially agrees with itself and would falsely
        # look converged. So growing past the cap gives at most one genuine comparison (against
        # the last, truly smaller, n_slabs); if that one comparison does not already satisfy both
        # safeguards, there is no more evidence to be had, and the fast method must give up.
        is_repeat = (n_slabs == n_slabs_prev)
        if is_repeat:  # pragma: no cover - unreachable; see below
            # Unreachable as the loop currently stands, and kept as a guard rather than
            # deleted.  Below the cap the slab count strictly increases (the growth factor
            # is 2, and the clause at the foot of the loop forces progress even if it were
            # not); at the cap every branch below returns within the same iteration.  So
            # the loop never survives a pass at the cap to make a repeated comparison, and
            # n_slabs never equals n_slabs_prev.  It becomes live again the moment the
            # growth factor, the cap, or the returns below change.
            return P_new, False

        if (P_prev is not None) and (max_omega < omega_trust_threshold):
            if np.all(np.abs(P_new - P_prev) <= atol + rtol*np.abs(P_prev)):
                consecutive_agreements += 1
                if at_cap or (consecutive_agreements >= 2):
                    return P_new, True
            else:
                consecutive_agreements = 0
                if at_cap:
                    return P_new, False
        else:
            consecutive_agreements = 0
            if at_cap:
                return P_new, False
        P_prev = P_new
        n_slabs_prev = n_slabs

        if loop_count >= loop_cap:
            return P_new, False

        n_slabs_old = n_slabs
        n_slabs = min(round(growth*n_slabs), n_slabs_cap)
        if (n_slabs == n_slabs_old) and (n_slabs < n_slabs_cap):  # pragma: no cover
            # The no-progress guard, and unreachable while the growth factor is 2:
            # round(2n) == n has no solution for n >= 1.  It exists so that a smaller
            # growth factor -- 1.1, say, which rounds to no change at small n -- cannot
            # turn this into an infinite loop, which is precisely when it stops being
            # dead code.
            n_slabs += 1
        loop_count += 1


def _osc_prob_ip_exp_dispatch(
    h_vac_energy_indep: np.ndarray,
    VCC_func: Union[Callable, float],
    h_matt: np.ndarray,
    h_liv_energy_indep: Optional[np.ndarray],
    n_liv: Optional[Union[int, float]],
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    scan_kwargs: Dict
):
    r"""Decide whether the fast interaction-picture integrator applies; run it if so.

    Returns ``NotImplemented`` when the request does not fit the fast method (``VCC_func`` is not a
    genuine exponential profile built via :func:`magnus.matter.exp_density_profile`, user-provided
    slab edges or breakpoints, iteration over the expansion order, or logging requested) or when it
    fails to converge within the requested tolerance (signaling that the matter term is not a small
    perturbation on the vacuum splitting somewhere along the trajectory, e.g., an MSW resonance); the
    caller falls back to the general per-point path in either case. Unlike
    ``_osc_prob_scan_separable_dispatch``, this applies equally to a single (energy, L) point (the
    common case for :func:`osc_prob_sun`-family calls) and to a multi-energy scan at a shared baseline.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    h_vac_energy_indep : np.ndarray
        Energy-independent part of the vacuum (and, if present, LIV) Hamiltonian.
    VCC_func : Callable or float
        Matter potential, as a function of position (required for the fast method to apply; a
        constant potential, or one not tagged as exponential, falls back to the generic path).
    h_matt : np.ndarray
        Constant matrix multiplying ``VCC_func(l)``.
    h_liv_energy_indep : np.ndarray, optional
        Energy-independent part of the LIV Hamiltonian, if any.
    n_liv : int or float, optional
        Power of the energy dependence of the LIV operator, if ``h_liv_energy_indep`` is given.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s); the fast method applies only when all requested baselines are equal.
    L0 : int or float
        Initial position.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    scan_kwargs : dict
        The refinement/logging keyword arguments of :func:`osc_prob_energy_baseline` (rtol, atol,
        growth_factor_n_slabs, max_num_loops, min_n_slabs, max_n_slabs, t_slab_edges,
        save_log, file_log) plus a nested 'kwargs' dict of any
        remaining, unrecognized keyword arguments.

    Returns
    -------
    np.ndarray or NotImplemented
        The oscillation probability (or single channel), computed via the fast method; or the
        ``NotImplemented`` singleton if the request does not fit it or it failed to converge.
    """
    if ('ip_exp' in _ENGINES_DISABLED) or (scan_kwargs.get('cumulative') is True):
        return NotImplemented
    kwargs = dict(scan_kwargs.get('kwargs', {}))
    t_breakpoints = kwargs.pop('t_breakpoints', None)
    n_slabs0 = kwargs.pop('n_slabs', 1)
    kwargs.pop('n_tpts_per_slab', None)
    if len(kwargs) > 0:
        return NotImplemented
    if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
        return NotImplemented
    if not isinstance(VCC_func, Callable):
        return NotImplemented
    if not getattr(VCC_func, 'is_exp_density_profile', False):
        return NotImplemented
    if np.asarray(h_vac_energy_indep).shape[-1] != 2:
        # The neglected second-order term Omega_2 involves a sum over every off-diagonal pair in
        # the H_E eigenbasis; empirically (see the validation grid in the module docstring... see
        # tests/test_oscprob.py), its coefficient grows by three orders of magnitude already from
        # 2 to 3 flavors (one pair vs three, with sizable diagonal mixing-angle contributions in
        # each), which pushes the slab count needed for a certified answer far past what stays
        # fast. The fast path is therefore restricted to genuinely 2-level Hamiltonians for now,
        # where it has been validated against solve_ivp across the realistic solar-neutrino
        # energy range; 3+ flavor (and LIV, which adds further energy-dependent terms to H_E) fall
        # back to the general method unconditionally.
        return NotImplemented
    if scan_kwargs['t_slab_edges'] is not None:
        return NotImplemented
    if scan_kwargs['save_log'] or (scan_kwargs['file_log'] is not None):
        return NotImplemented

    l_scale = VCC_func.l_scale
    energy_arr, L_arr, return_float, ok = _normalize_energy_L(energy, L)
    if (not ok) or (not np.all(L_arr == L_arr[0])):
        return NotImplemented

    rtol, atol = scan_kwargs['rtol'], scan_kwargs['atol']
    if (rtol is None) != (atol is None):
        rtol = 0.0 if rtol is None else rtol
        atol = 0.0 if atol is None else atol

    H_E = (1.0/energy_arr)[:, None, None]*np.asarray(h_vac_energy_indep, dtype=complex)
    if h_liv_energy_indep is not None:
        H_E = H_E + (energy_arr**n_liv)[:, None, None]*np.asarray(h_liv_energy_indep, dtype=complex)

    P, converged = _osc_prob_ip_exp_core(H_E, l_scale, VCC_func, np.asarray(h_matt), float(L0),
        float(L_arr[0]), rtol, atol, scan_kwargs['growth_factor_n_slabs'],
        scan_kwargs['max_num_loops'], scan_kwargs['min_n_slabs'], scan_kwargs['max_n_slabs'], n_slabs0)
    if not converged:
        _note_engine('ip_exp', answered=False, reason='did not converge')
        return NotImplemented

    _note_engine('ip_exp')
    if (nu_i is not None) and (nu_f is not None):
        P = P[:, nu_i, nu_f]
    return P.__getitem__(0 if return_float else slice(None))


def _cumulative_scan_would_serve(energy_arr, L_arr, L0, min_points):
    r"""Whether the cumulative scan applies to this set of requested points.

    Used in two places with different ``min_points``, because the two guard different trades:
    ``osc_prob_energy_baseline`` resolves ``cumulative='auto'`` with
    :data:`CUMULATIVE_AUTO_MIN_POINTS`, while ``_osc_prob_hybrid_dispatch`` decides whether to
    stand aside with the larger :data:`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`.  Sharing the
    rest of the predicate keeps them from drifting apart on the conditions that are genuinely
    the same.

    That the dispatcher's threshold is the larger one is what makes the fall-through safe: when
    the hybrid dispatcher declines on this count, ``'auto'`` is guaranteed to accept, so a scan
    can never be declined by both and land on the general per-point path -- slower and less
    accurate than either, and silently so.

    The remaining ``'auto'`` conditions (a position-dependent Hamiltonian, no ``t_slab_edges``)
    are already guaranteed where the dispatcher calls this: it has returned ``NotImplemented``
    for a non-callable ``VCC_func`` and for user-supplied slab edges or breakpoints above.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy_arr : np.ndarray
        Requested energies, one per point.
    L_arr : np.ndarray
        Requested baselines, one per point.
    L0 : int or float
        Initial position.
    min_points : int
        Fewest points at which the cumulative scan is worth taking, for this caller.

    Returns
    -------
    bool
        True when the cumulative scan applies to this set of points.
    """
    return bool(len(L_arr) >= min_points
                and np.all(energy_arr == energy_arr[0])
                and np.all(np.asarray(L_arr, dtype=float) >= L0))


def _resolve_cumulative_kwarg(kwargs, strategy):
    r"""Pops a caller-supplied ``cumulative`` out of ``kwargs`` and decides what to forward.

    The three scenario wrappers (:func:`osc_prob_matter_std_potential`,
    :func:`osc_prob_matter_nsi`, :func:`osc_prob_liv`) each set ``cumulative`` themselves when
    calling :func:`osc_prob_energy_baseline`, so a caller who also passed it in ``**kwargs``
    used to get ``TypeError: got multiple values for keyword argument 'cumulative'`` --
    which made ``cumulative=False`` unreachable from the entire wrapper layer, and that is the
    one mitigation ``docs/dev/DECISION_CUMULATIVE_DEFAULT.md`` names for a caller who needs
    bit-reproducibility against pre-1.0.0 results.

    **An explicit value from the caller always wins**, including over the ``strategy='magnus'``
    opt-out below: naming ``cumulative`` is a specific request about the scan engine, and
    ``cumulative=True`` is documented to *raise* rather than fall back when it cannot be served,
    so honouring it cannot silently do the wrong thing.

    Otherwise ``strategy='magnus'`` resolves to ``False`` and everything else to ``'auto'``.
    That strategy promises the behaviour Mag(nu)s had before the adiabatic strategy existed,
    *unconditionally*; the cumulative scan is Magnus machinery but postdates the promise and
    builds a different grid, so the escape hatch would quietly stop being one for exactly the
    case -- a single-energy baseline scan -- where someone reproducing older numbers reaches
    for it.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    kwargs : dict
        The wrapper's ``**kwargs``, modified in place: ``'cumulative'`` is removed if present.
    strategy : str
        'auto', 'hybrid' or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential`.

    Returns
    -------
    bool or str
        The value to forward to :func:`osc_prob_energy_baseline` as ``cumulative``.
    """
    if 'cumulative' in kwargs:
        return kwargs.pop('cumulative')
    return 'auto' if strategy != 'magnus' else False


def _osc_prob_hybrid_dispatch(
    h_vac_energy_indep: np.ndarray,
    VCC_func: Union[Callable, float],
    h_matt: np.ndarray,
    h_liv_energy_indep: Optional[np.ndarray],
    n_liv: Optional[Union[int, float]],
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    scan_kwargs: Dict,
    strategy: str
):
    r"""Decide whether the adiabatic-transport-plus-Magnus-patch ("hybrid") strategy applies; run
    it if so.

    Unlike ``_osc_prob_ip_exp_dispatch``, this is not restricted to a genuine exponential
    profile or to two flavors: :func:`magnus.adiabatic.hybrid_propagator` locates non-adiabatic
    windows via an exact Hellmann-Feynman diagnostic that makes no assumption about the profile's
    functional form or the Hamiltonian's dimension (see :doc:`/adiabatic_strategy`). It is still
    restricted to a smooth ``VCC_func`` (no user-supplied slab edges or breakpoints -- a
    piecewise-discontinuous profile such as PREM breaks the finite-difference diagnostics this
    method relies on) and, since the method is fundamentally adaptive, to a genuinely requested
    tolerance (``rtol``/``atol`` not both ``None``).

    Each requested (energy, L) point is handled by an independent call to
    :func:`magnus.adiabatic.hybrid_propagator`, since the position of a resonance (if any) is
    generally energy-dependent -- unlike ``_osc_prob_ip_exp_dispatch``, this applies equally
    to a scan with per-point baselines, not only a shared one. If any requested point fails to
    self-certify, the whole batch is treated as not fitting this method with ``strategy='auto'``
    (returns ``NotImplemented``, so the caller falls back to the general per-point path);
    with ``strategy='hybrid'``, the best-effort result is returned together with
    ``HybridCertificationWarning``.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    h_vac_energy_indep : np.ndarray
        Energy-independent part of the vacuum Hamiltonian.
    VCC_func : Callable or float
        Matter potential, as a function of position (required for this method to apply; a
        constant potential falls back to the generic path, since there is then no position
        dependence for a resonance to hide in).
    h_matt : np.ndarray
        Constant matrix multiplying ``VCC_func(l)``.
    h_liv_energy_indep : np.ndarray, optional
        Energy-independent part of the LIV Hamiltonian, if any.
    n_liv : int or float, optional
        Power of the energy dependence of the LIV operator, if ``h_liv_energy_indep`` is given.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    scan_kwargs : dict
        The refinement/logging keyword arguments of :func:`osc_prob_energy_baseline` (rtol, atol,
        magnus_exp_order, integration_method, t_slab_edges,
        save_log, file_log) plus a nested 'kwargs' dict of any remaining, unrecognized keyword
        arguments.
    strategy : str
        'auto', 'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential`. If 'magnus', this function always returns
        ``NotImplemented`` without doing any work.

    Returns
    -------
    np.ndarray or NotImplemented
        The oscillation probability (or single channel), computed via the hybrid strategy; or the
        ``NotImplemented`` singleton if the request does not fit it, ``strategy == 'magnus'``, or
        (only with ``strategy == 'auto'``) it failed to self-certify for at least one point.
    """
    if (strategy == 'magnus') or ('hybrid' in _ENGINES_DISABLED):
        return NotImplemented
    if scan_kwargs.get('cumulative') is True:
        # An explicit cumulative=True is a request for one engine in particular, documented to
        # raise rather than fall back if it cannot be served.  Answering it here would be a
        # silent substitution -- the exact thing that flag exists to rule out.
        return NotImplemented

    kwargs = dict(scan_kwargs.get('kwargs', {}))
    t_breakpoints = kwargs.pop('t_breakpoints', None)
    kwargs.pop('n_slabs', None)
    kwargs.pop('n_tpts_per_slab', None)
    if len(kwargs) > 0:
        return NotImplemented
    if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
        return NotImplemented
    if not isinstance(VCC_func, Callable):
        return NotImplemented
    if scan_kwargs['t_slab_edges'] is not None:
        return NotImplemented
    if scan_kwargs['save_log'] or (scan_kwargs['file_log'] is not None):
        return NotImplemented

    rtol, atol = scan_kwargs['rtol'], scan_kwargs['atol']
    if (rtol is None) and (atol is None):
        # The hybrid method is fundamentally adaptive (self-certifying); "run once with a fixed,
        # non-adaptive number of slabs" has no hybrid analogue, so fall back to the general path,
        # which does support it.
        return NotImplemented
    rtol = 0.0 if rtol is None else rtol
    atol = 0.0 if atol is None else atol

    energy_arr, L_arr, return_float, ok = _normalize_energy_L(energy, L)
    if not ok:
        return NotImplemented

    # Stand aside for a *large enough* baseline scan at a single energy: the cumulative scan
    # answers all of those baselines from one traversal, where this method calls
    # hybrid_propagator once per point at its ~20 ms floor.  Measured through osc_prob_2nu_sun
    # on a solar profile at N = 400: 7.5 s -> 0.29 s, with the error improving from ~1e-5 to
    # ~1e-6 as well.
    #
    # The threshold is not 2.  This method is accurate and cheap per point, so below
    # HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS the cumulative scan's near-constant cost -- its
    # strict probe -- is not yet amortised, and yielding would make a small scan several times
    # slower (7.6x at N = 2) to buy accuracy that was already two orders inside what the caller
    # asked for.  See that constant for the measurements.
    #
    # Only under strategy='auto', which promises the best available answer rather than this
    # method in particular; strategy='hybrid' is an explicit request and still gets hybrid.
    # Declining here is enough to reach the cumulative path: ip_exp needs every baseline equal
    # and the separable engine needs a single shared baseline, so both decline a scan too, and
    # the caller falls through to osc_prob_energy_baseline, where cumulative='auto' engages --
    # guaranteed, since its threshold is the smaller one.
    if (strategy == 'auto') and _cumulative_scan_would_serve(
            energy_arr, L_arr, L0, HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS):
        return NotImplemented

    magnus_exp_order = scan_kwargs['magnus_exp_order']
    integration_method = scan_kwargs['integration_method']
    h_vac_energy_indep = np.asarray(h_vac_energy_indep, dtype=complex)
    h_matt = np.asarray(h_matt, dtype=complex)
    if h_liv_energy_indep is not None:
        h_liv_energy_indep = np.asarray(h_liv_energy_indep, dtype=complex)
    d = h_vac_energy_indep.shape[-1]

    def H_at_energy(enu):
        def H_of_l(l, enu=enu):
            vcc = np.asarray(VCC_func(l))
            H = (1.0/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt
            if h_liv_energy_indep is not None:
                H = H + (enu**n_liv)*h_liv_energy_indep
            return H
        return H_of_l

    P_out = _hybrid_propagator_scan(H_at_energy, energy_arr, L_arr, L0, rtol, atol,
        magnus_exp_order, integration_method, strategy, d)
    if P_out is NotImplemented:
        return NotImplemented

    if (nu_i is not None) and (nu_f is not None):
        P_out = P_out[:, nu_i, nu_f]
    return P_out.__getitem__(0 if return_float else slice(None))


def _warn_hybrid_unresolved() -> None:
    r"""The hybrid strategy declined because ``H_func`` is not resolved at the probe scale.

    One function rather than two call sites with the same string: this fires both when
    ``strategy='auto'`` declines (and the general path answers) and when ``strategy='hybrid'``
    was forced, and the two must not drift apart.

    .. versionadded:: 1.0.0
    """
    warnings.warn(
        "osc_prob (hybrid strategy): the Hamiltonian is not resolved at the scale this method "
        "samples it on -- a density jump, or a feature narrower than the probe grid can see -- "
        "and no t_breakpoints were given. The adiabatic strategy is built on finite differences "
        "of H between probe points, which mean nothing across a jump, so it declined; the "
        "answer comes from the general Magnus path instead, which is correct there but slower, "
        "and a slab straddling the same feature still limits its accuracy. Pass t_breakpoints "
        "at the feature: it is the cure in both cases. Measured on an unmarked density step, "
        "the adiabatic answer was wrong by 0.54 in probability while reporting itself "
        "certified. Shown once per session.",
        UnmarkedDiscontinuityWarning, stacklevel=4)


def _hybrid_propagator_scan(
    H_at_energy: Callable,
    energy_arr: np.ndarray,
    L_arr: np.ndarray,
    L0: Union[int, float],
    rtol: float,
    atol: float,
    magnus_exp_order: int,
    integration_method: str,
    strategy: str,
    d: int
):
    r"""Shared per-(energy, L)-point hybrid-propagator loop, used by both
    ``_osc_prob_hybrid_dispatch`` (separable vacuum + matter potential Hamiltonians) and
    ``_osc_prob_hybrid_dispatch_generic`` (an arbitrary user-supplied Hamiltonian, as accepted
    by :func:`osc_prob_earth`/:func:`osc_prob_sun`).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_at_energy : Callable
        Given a neutrino energy, returns ``H_of_l(l)``, the Hamiltonian at that energy as a
        function of position only, suitable for :func:`magnus.adiabatic.hybrid_propagator`.
    energy_arr : np.ndarray
        Neutrino energies, one per requested point.
    L_arr : np.ndarray
        Baselines, one per requested point (same length as ``energy_arr``).
    L0 : int or float
        Initial position (shared by every point).
    rtol, atol : float
        Target relative/absolute tolerance, already coerced to non-``None`` floats by the caller.
    magnus_exp_order : int
        Magnus expansion order used for the local patch inside each non-adiabatic window.
    integration_method : str
        Integration method used for the local patch.
    strategy : str
        'auto' or 'hybrid' (never 'magnus'; the caller already handles that case). If 'auto', any
        point that fails to self-certify aborts the whole scan (returns ``NotImplemented``); if
        'hybrid', the best-effort result is kept and ``HybridCertificationWarning`` is raised once
        at the end if at least one point was uncertified.
    d : int
        Hamiltonian dimension (number of flavors).

    Returns
    -------
    np.ndarray or NotImplemented
        Stacked probability matrices, shape ``(len(energy_arr), d, d)``; or ``NotImplemented`` if
        ``strategy == 'auto'`` and at least one point failed to self-certify.
    """
    n_pts = len(energy_arr)
    P_out = np.empty((n_pts, d, d))
    any_uncertified = False
    unresolved = False

    for i in range(n_pts):
        H_of_l = H_at_energy(energy_arr[i])

        info = {}
        U, _, certified = adiabatic.hybrid_propagator(H_of_l, float(L0), float(L_arr[i]),
            rtol=rtol, atol=atol, magnus_exp_order=magnus_exp_order,
            integration_method=integration_method, info=info)
        unresolved = unresolved or (not info.get('resolved', True))


        if not certified:
            if strategy == 'auto':
                _note_engine('hybrid', answered=False, certified=False,
                    reason=('the profile is not resolved at the probe scale'
                            if unresolved
                            else 'did not self-certify at the requested tolerance'))
                if unresolved:
                    _warn_hybrid_unresolved()
                return NotImplemented
            any_uncertified = True

        P_out[i] = np.swapaxes(U.real**2 + U.imag**2, -1, -2)

    _note_engine('hybrid', certified=not any_uncertified)
    if unresolved:
        # A new instance of an existing kind.  adiabatic._profile_is_resolved already runs
        # inside hybrid_propagator, where failing it causes a decline rather than a message --
        # so on an undeclared density jump the caller heard about slab widths from whichever
        # engine answered instead, which is true but points at the wrong knob.  Measured before
        # shipping, over the axis the earlier measurement did not have -- one call per baseline
        # means one call per sub-interval, not one per profile: 0 false positives over 1440
        # smooth configurations, and every sub-interval containing a jump caught.  Getting
        # there needed a repair to the detector; see adiabatic.LOCAL_JUMP_RATIO.
        _warn_hybrid_unresolved()

    if any_uncertified:
        warnings.warn("osc_prob (hybrid strategy): requested tolerance not achieved for at "
            "least one (energy, L) point; the returned probabilities remain exactly unitary "
            "but their accuracy is not certified -- unverified, which is not the same as "
            "wrong. To get a certified answer: use strategy='auto', which falls back to the "
            "general Magnus path for exactly these points; or, if the profile has known "
            "structure (a density jump, a kink, a feature narrower than 1/200 of the "
            "trajectory), pass t_breakpoints there, which is the one cure for a feature the "
            "probe grid cannot resolve; or request a looser rtol/atol, if the accuracy you "
            "need is less than the accuracy you asked for. Shown once per session.",
            HybridCertificationWarning, stacklevel=3)

    return P_out


def _osc_prob_hybrid_dispatch_generic(
    htot: Callable,
    VCC_func: Union[Callable, float],
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    t_breakpoints: Optional[np.ndarray],
    rtol: Optional[Union[int, float]],
    atol: Optional[Union[int, float]],
    magnus_exp_order: int,
    integration_method: str,
    strategy: str,
    kwargs: Dict
):
    r"""Decide whether the hybrid strategy applies to an arbitrary, user-supplied Hamiltonian (as
    accepted by :func:`osc_prob_earth`/:func:`osc_prob_sun`); run it if so.

    Same method and gating philosophy as ``_osc_prob_hybrid_dispatch`` (see its docstring and
    :doc:`/adiabatic_strategy`), adapted to ``htot(energy, l)`` -- the Hamiltonian already unified
    into a single two-argument function by ``_osc_prob_with_potential``, regardless of
    whether the user's own ``H_func`` takes ``(energy, l, VCC)`` or ``(energy, l)`` -- instead of
    a separable ``h_vac_energy_indep``/``VCC_func``/``h_matt`` decomposition, since no such
    decomposition is available (or needed: the resonance detector and adiabatic propagator make
    no assumption about the Hamiltonian's internal structure) for a fully generic ``H_func``.

    In practice, this means :func:`osc_prob_earth` almost always falls back to the ``'magnus'``
    strategies regardless of what ``strategy`` is requested, since ``t_breakpoints`` (the PREM
    layer-boundary crossings) is essentially always non-empty for a real Earth-crossing
    trajectory; :func:`osc_prob_sun` has no such restriction, since its density profile has no
    breakpoints.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    htot : Callable
        The Hamiltonian, as a function of ``(energy, l)`` -- already unified by
        ``_osc_prob_with_potential`` from the user's own ``H_func(energy, l, VCC)`` or
        ``H_func(energy, l)``.
    VCC_func : Callable or float
        The environment's matter potential, as a function of position (required for this method
        to apply; a constant potential falls back to the generic path).
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    t_breakpoints : np.ndarray, optional
        Mandatory slab edges (e.g., PREM layer boundaries); a non-empty array disables the hybrid
        dispatch (see above).
    rtol, atol : int or float, optional
        Target relative/absolute tolerance requested by the caller.
    magnus_exp_order : int
        Magnus expansion order used for the local patch inside each non-adiabatic window.
    integration_method : str
        Integration method used for the local patch.
    strategy : str
        'auto', 'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential`. If 'magnus', this function always returns
        ``NotImplemented`` without doing any work.
    kwargs : dict
        Any additional, unrecognized keyword arguments the caller received (e.g., forwarded from
        :func:`osc_prob_earth`/:func:`osc_prob_sun`'s own ``**kwargs``); any entry other than
        ``n_slabs``/``n_tpts_per_slab`` disables the hybrid dispatch, signaling that the caller
        wants low-level control of the general slab-refinement method specifically.

    Returns
    -------
    np.ndarray or NotImplemented
        The oscillation probability (or single channel), computed via the hybrid strategy; or the
        ``NotImplemented`` singleton if the request does not fit it, ``strategy == 'magnus'``, or
        (only with ``strategy == 'auto'``) it failed to self-certify for at least one point.
    """
    if (strategy == 'magnus') or ('hybrid' in _ENGINES_DISABLED):
        return NotImplemented

    if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
        return NotImplemented
    if not isinstance(VCC_func, Callable):
        return NotImplemented
    if len(set(kwargs) - {'n_slabs', 'n_tpts_per_slab'}) > 0:
        return NotImplemented
    if (rtol is None) and (atol is None):
        return NotImplemented
    rtol = 0.0 if rtol is None else rtol
    atol = 0.0 if atol is None else atol

    energy_arr, L_arr, return_float, ok = _normalize_energy_L(energy, L)
    if not ok:
        return NotImplemented

    def H_at_energy(enu):
        def H_of_l(l, enu=enu):
            return htot(enu, l)
        return H_of_l

    d = np.asarray(htot(energy_arr[0], L0)).shape[-1]

    P_out = _hybrid_propagator_scan(H_at_energy, energy_arr, L_arr, L0, rtol, atol,
        magnus_exp_order, integration_method, strategy, d)
    if P_out is NotImplemented:
        return NotImplemented

    if (nu_i is not None) and (nu_f is not None):
        P_out = P_out[:, nu_i, nu_f]
    return P_out.__getitem__(0 if return_float else slice(None))


def _cumulative_scan_grid(L_out, L0, n_acc, t_breakpoints):
    """Slab edges for a cumulative baseline scan, and where each output sits in them.

    The grid is the union of three things, and each is there for a different reason:

    * the **requested baselines**, so that every answer lands exactly on a slab edge and is
      read off the running product rather than interpolated;
    * a **uniform accuracy grid** of ``n_acc`` slabs, because the requested baselines are
      typically logarithmically spaced and so are dense where accuracy is cheap and sparse
      where it is expensive -- the opposite of what the integration needs;
    * the **breakpoints**, for the usual reason: a slab straddling a density discontinuity
      degrades the quadrature no matter how high ``magnus_exp_order`` is.

    Returns
    -------
    (np.ndarray, np.ndarray)
        Sorted, unique edges spanning ``[L0, max(L_out)]``, and the index into them of each
        entry of ``L_out``.

    .. versionadded:: 1.0.0
    """
    L_out = np.asarray(L_out, dtype=float)
    parts = [np.linspace(L0, float(L_out[-1]), int(n_acc) + 1), L_out, np.array([L0])]
    if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
        bp = np.atleast_1d(np.asarray(t_breakpoints, dtype=float))
        parts.append(bp[(bp > L0) & (bp < L_out[-1])])
    edges = np.unique(np.concatenate(parts))
    return edges, np.searchsorted(edges, L_out)


def _osc_prob_cumulative_scan(H_func, L_out, L0, n_acc, magnus_exp_order,
                              n_tpts_per_slab, integration_method, t_breakpoints,
                              A_eval_mode, **kwargs):
    r"""Every baseline in ``L_out`` from a single traversal of the profile.

    The evolution operator is a time-ordered product, so each requested answer is a *prefix*
    of the next: :math:`U(0 \to L_2) = U(L_1 \to L_2)\,U(0 \to L_1)`.  Computing N baselines
    independently therefore re-walks the profile N times over.  This walks it once and
    records the running product wherever an answer was asked for -- the ``reduce`` in
    :func:`osc_prob` with its intermediates kept rather than discarded.

    Two properties are requirements rather than optimisations, and both are about memory:

    * the traversal is **chunked**, so the slab operators are never all live at once;
    * each snapshot is converted to a probability **immediately**, so the recorded term
      collapses into the result the caller asked for instead of sitting beside it as N
      complex unitaries.

    Together they make peak memory ``O(block) + O(result)``, which above the block size is
    less than the per-point path uses for the same scan.

    Parameters
    ----------
    H_func : Callable
        Hamiltonian as a function of position alone; the energy is already bound.
    L_out : np.ndarray
        Requested baselines, strictly ascending and all greater than ``L0``.
    n_acc : int
        Slabs the accuracy grid would use over the whole path on its own; see
        :func:`_cumulative_scan_grid`.

    Returns
    -------
    np.ndarray
        Probability matrices, shape ``(len(L_out), d, d)``.

    .. versionadded:: 1.0.0
    """
    edges, out_idx = _cumulative_scan_grid(L_out, L0, n_acc, t_breakpoints)
    n_slabs = len(edges) - 1

    # A position-independent Hamiltonian arrives here as a bare array, as it does in
    # osc_prob; wrap it in the array-capable function the slab kernel expects.  A single slab
    # would be exact for it, but the caller asked for a scan and the extra slabs cost little
    # -- and keeping one code path avoids a second place where the time ordering could differ.
    if not callable(H_func):
        H_const = np.asarray(H_func)

        def H_func(l, _H=H_const):
            return np.broadcast_to(_H, np.shape(l) + _H.shape) if np.ndim(l) else _H

        A_eval_mode = 'vector'

    dim = np.asarray(H_func(L0)).shape[-1]
    P = np.empty((len(L_out), dim, dim))
    running = np.eye(dim, dtype=complex)

    # A baseline equal to L0 is the identity: no slab precedes it.
    for j in np.flatnonzero(out_idx == 0):
        P[j] = np.transpose(running.real**2 + running.imag**2)

    _, block = _tile_for_working_set(1, n_slabs, dim*dim, live_arrays=8)
    for start in range(0, n_slabs, block):
        stop = min(start + block, n_slabs)
        U = compute_evolution_operator_multiple_slabs(
            H_func, np.column_stack([edges[start:stop], edges[start + 1:stop + 1]]),
            n_tpts_per_slab, magnus_exp_order, integration_method=integration_method,
            A_eval_mode=A_eval_mode, **kwargs)
        # Outputs landing inside this block, in edge order, so the running product is
        # snapshotted at the right moment without a second pass.
        here = np.flatnonzero((out_idx > start) & (out_idx <= stop))
        order = here[np.argsort(out_idx[here], kind='stable')]
        pos = 0
        for k in range(stop - start):
            running = U[k] @ running
            while (pos < len(order)) and (out_idx[order[pos]] == start + k + 1):
                j = order[pos]
                P[j] = np.transpose(running.real**2 + running.imag**2)
                pos += 1
        del U
    return P


def osc_prob_energy_baseline(
    H_func: Union[Callable, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Optional[Union[int, float]]=0.0,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    H_func_is_function_only_of_energy: Optional[bool]=False,
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='gl', 
    rtol: Optional[Union[int, float]]=1.e-3, 
    atol: Optional[Union[int, float]]=1.e-3, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=None, 
    min_n_tpts_per_slab: Optional[int]=2, 
    max_n_tpts_per_slab: Optional[int]=500, 
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    cumulative: Optional[Union[bool, str]]='auto',
    **kwargs
) -> Union[int, float, np.ndarray]:
    r"""Compute and return oscillation probabilities for given arrays of
    neutrino energy and baseline, and an arbitrary Hamiltonian.

    Sits directly above :func:`osc_prob` in the primordial layer (see
    :doc:`/architecture`): given arrays of ``energy`` and ``L``, builds the
    right energy-dependent closure over ``H_func``, decides whether to parallelize over
    (energy, L) points or hand a single call straight to :func:`osc_prob`, and carries the warm
    start logic that seeds each point's refinement from the previous point's converged
    (``n_slabs``, ``n_tpts_per_slab``). Called directly by :func:`osc_prob_vacuum`,
    :func:`osc_prob_matter_std_potential`, :func:`osc_prob_matter_nsi`, and :func:`osc_prob_liv`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable or np.ndarray
        The Hamiltonian: a function of energy only, of position only, of both (in that
        parameter order), or a constant matrix.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s). Must have the same length as ``energy``, or either may be a single value
        broadcast against the other.
    L0 : int or float, optional
        Initial position. Default: 0.0.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    H_func_is_function_only_of_energy : bool, optional
        If True and ``H_func`` accepts a single argument, treat it as energy-only (returning a
        constant matrix per energy) rather than position-only. Default: False.
    t_slab_edges : list or np.ndarray, optional
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    magnus_exp_order : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    n_jobs : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    integration_method : str
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    rtol : int or float, optional
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    atol : int or float, optional
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    growth_factor_n_slabs : int or float
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.

        **Provenance.**  Swept over 1.2, 1.5, 2.0, 3.0 across 18 workloads spanning single
        points, baseline scans and energy scans
        (``docs/dev/adversarial_batteries/constants_audit2.py``): worst error 4.49e-04 at every
        value.  It sets how coarsely the ladder is sampled, not where it stops, so it trades
        wasted refinement against overshoot without moving the accepted answer.
    growth_factor_n_tpts_per_slab : int or float
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    max_num_loops : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    min_n_slabs : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    max_n_slabs : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    min_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.

        **Provenance.**  Swept over 2, 4, 8 on the same 18 workloads: worst error 4.49e-04 at
        every value.  With the default ``integration_method='gl'`` this is expected --
        Gauss-Legendre pins the node count per slab and ignores it -- so the sweep confirms the
        documented behaviour rather than calibrating anything.
    max_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    validate_input : bool
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    save_log : bool
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    filename_log : str
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    file_log : TextIOWrapper, optional
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    close_file_log_upon_exit : bool
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    new_recursion_limit : int, optional
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    verbose : int
        Forwarded to :func:`osc_prob` for each (energy, L) point; see its docstring.
    cumulative : bool or str, optional
        Compute a whole baseline scan from **one** traversal of the profile instead of one
        traversal per baseline.  The evolution operator is a time-ordered product, so
        :math:`U(0 \to L_2) = U(L_1 \to L_2)\,U(0 \to L_1)`: each requested answer is a
        prefix of the next, and recording the running product yields all of them at once.
        One of:

        * ``'auto'`` (default) -- use the cumulative scan whenever the request fits it, and
          the ordinary per-point path otherwise.  The request fits when the Hamiltonian varies
          with position, all the energies are equal, no ``t_slab_edges`` were given, every
          baseline is at or beyond ``L0``, and there are at least
          ``CUMULATIVE_AUTO_MIN_POINTS`` of them.  A position-independent Hamiltonian
          (vacuum, constant density) is excluded because :func:`osc_prob` integrates it
          exactly on a single slab, leaving no traversal to share.
        * ``True`` -- require it, and **raise** if the request does not fit.  Use this when
          the cumulative scan is what you want and silently getting the per-point path
          instead would be a problem.
        * ``False`` -- never use it.

        Applies to a **baseline scan at a single energy**.  The nesting it exploits belongs to
        the baseline axis alone -- :math:`P(E_1)` shares nothing with :math:`P(E_2)`, since
        each energy needs its own propagation through the whole profile -- so there is no
        energy-axis counterpart to this.

        Not compatible with ``t_slab_edges``, which it would have to override: the scan
        builds a grid that is the union of the requested baselines, an accuracy grid, and
        any ``t_breakpoints``.  The accuracy grid is sized by one ordinary adaptive
        :func:`osc_prob` call at the longest baseline, so the usual tolerance machinery and
        its warnings apply unchanged.

        The default became ``'auto'`` in 1.0.0, having been ``False``, because the cumulative
        scan measured **more accurate at every scan size tested**, not merely faster.  Against
        ``solve_ivp`` on a 5 MeV solar scan to one solar radius, the per-point path it replaces
        returns answers outside the requested 1e-3 at several sizes -- 9.7e-3 at N = 10, 5.6e-3
        at N = 25, 2.6e-3 at N = 100 -- where the cumulative scan stays near 5e-6 throughout.
        Speed follows from N ~ 25 upward (2.65x there, 84x at N = 1000); below it the
        cumulative scan can be ~1.3x slower in wall time, which is a few milliseconds.

        Because the two paths build different grids, results move -- within the requested
        tolerance, and generally toward the truth.  Pass ``cumulative=False`` to reproduce
        pre-1.0.0 numbers exactly.
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob`.

    Returns
    -------
    int, float, or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each (energy, L) point; a single value/matrix if both ``energy`` and ``L`` were floats.
    """

    if (isinstance(H_func, Callable) and (_n_required_params(H_func) > 2)):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_energy_baseline:"+\
            " H_func can be energy- and position-dependent, only energy-dependent, or only" + \
            " position-dependent. H_func cannot depend on more than two parameters. To vary" + \
            " the third parameter, call osc_prob_energy_baseline within a loop where it is" + \
            " varied.")

    # Turn int into float
    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L

    # Flag return_float remembers if energy and L were both floats.  If True,
    # osc_prob_energy_baseline returns a float, too.
    return_float = isinstance(energy, float) and isinstance(L, float)

    # If there is a single value of energy, make an array out of it.  Same for L.  This will allow
    # us to zip them later.
    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)
    L = np.array([L]) if isinstance(L, float) else np.array(L)

    # Either energy and L are both lists (or NumPy arrays) of the same length; or one is a float and
    # the other is a list (or NumPy array).  Any other possibility will generate an exception.  This
    # exception may be raised earlier in routines that call osc_prob_energy_baseline if they are
    # called wih validate_input == True, but we check below in case it osc_prob_energy_baseline was
    # set to False.
    if not ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or \
        (len(energy) > 1 and len(L) == 1)):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_energy_baseline: energy and L must be both " + \
            "int or float; or, if lists (or NumPy arrays), they must have the same length;" + \
            " or, if one is a float or single-entry list, the other must be a list with " + \
            "multiple entries.")

    # If energy is a single value, then transform it into an array containing the value energy
    # repeated a number of times equal to the length of the L, and vice versa, in order to zip them.
    energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
    L = np.full(len(energy), L[0]) if (len(L) == 1) else L

    n_points = len(energy)

    # When there are multiple (energy, L) points and n_jobs != 1, parallelize over the points, and
    # run each individual osc_prob call serially.  The per-point tasks are large enough for
    # process-based parallelism to pay off, unlike the much smaller per-slab tasks inside osc_prob.
    parallelize_over_points = (n_jobs != 1) and (n_points > 1)

    # Keyword arguments common to all the calls to osc_prob below.  Additional keyword arguments
    # received in **kwargs are passed through to osc_prob as well (e.g., n_slabs,
    # n_tpts_per_slab).
    osc_prob_kwargs = dict(
        t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order,
        n_jobs=1 if parallelize_over_points else n_jobs,
        integration_method=integration_method, rtol=rtol, atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit, verbose=verbose, **kwargs)

    # Build, for a given neutrino energy, the Hamiltonian to be passed to osc_prob: either a
    # one-parameter function of position or, if position-independent, a constant matrix (which
    # osc_prob detects and handles with internal speed-ups).
    if not isinstance(H_func, Callable):
        # H_func is position- and energy-independent
        def H_at_energy(enu: float) -> np.ndarray:
            return H_func
    elif (_n_required_params(H_func) == 2):
        # H_func is a function of two parameters; it is assumed that the first parameter is the
        # energy and the second one is the position
        def H_at_energy(enu: float) -> Callable:
            return lambda l: H_func(enu, l)
    elif H_func_is_function_only_of_energy:
        # H_func is a function only of energy: at fixed energy, it is a constant matrix
        def H_at_energy(enu: float) -> np.ndarray:
            return H_func(enu)
    else:
        # H_func is a function only of position
        def H_at_energy(enu: float) -> Callable:
            return H_func

    # Probe once how the Hamiltonian can be evaluated (vectorized over an array of positions,
    # constant, or scalar-only): the verdict is structural and holds for every (energy, L) point,
    # so probing here avoids re-probing inside every osc_prob call.
    H_first = H_at_energy(energy[0])
    if isinstance(H_first, Callable):
        osc_prob_kwargs['A_eval_mode'] = magnus.probe_eval_mode(
            lambda t: -1j*H_first(t), L0, np.max(L))

    # Refuse a request whose *result* cannot fit, before the scan allocates anything.  Every
    # batched engine either runs from here or falls back to here, and their working sets are
    # now tiled to a fixed budget, which leaves the result array as the only quantity still
    # free to grow without bound.  Placed here because this is the first point at which the
    # flavor count is known without evaluating the Hamiltonian specially for it: H_first has
    # been built already, and probe_eval_mode has just called it if it is a function of
    # position.  Costs one multiply for ordinary requests; see _check_output_fits.
    _check_output_fits(
        n_points,
        np.asarray(H_first(L0) if isinstance(H_first, Callable) else H_first).shape[-1],
        'osc_prob_energy_baseline')

    # Cumulative baseline scan: one traversal of the profile for every requested baseline,
    # instead of one traversal per baseline.  Applicable only when the Hamiltonian is the same
    # for every requested point -- i.e. a single energy, scanned over baselines -- because the
    # nesting it exploits is a property of the baseline axis alone: P(0->L1) is a prefix of
    # P(0->L2), while P(E1) shares nothing with P(E2).
    #
    # cumulative='auto' (the default) resolves to True exactly when the request fits, and to
    # False otherwise, so that the explicit cumulative=True can keep *raising* on a request it
    # cannot serve -- a caller who asked for it by name should hear that it did not happen.
    #
    # 'auto' adds two requirements beyond what cumulative=True checks:
    #
    #   - at least two baselines.  A single point has no prefix to reuse and would pay the
    #     inherited-grid probe for nothing, which matters because every single-point call
    #     through the wrapper layer is served from here.
    #   - a position-*dependent* Hamiltonian.  When H does not vary along the trajectory
    #     (vacuum, constant density), osc_prob integrates it exactly on one slab, so there is
    #     no traversal to share and the cumulative scan is strictly worse: it sizes a grid from
    #     an adaptive probe and then walks it.  H_first is the Hamiltonian at the first energy,
    #     already built above; it is a plain matrix exactly when the profile is constant.
    if isinstance(cumulative, str):
        if cumulative != 'auto':
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_energy_baseline: "
                "cumulative must be True, False, or 'auto'; got " + repr(cumulative) + ".")
        cumulative = bool(
            isinstance(H_first, Callable)
            and (t_slab_edges is None)
            and _cumulative_scan_would_serve(np.asarray(energy), np.asarray(L), L0,
                                             CUMULATIVE_AUTO_MIN_POINTS))

    if cumulative:
        if t_slab_edges is not None:
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_energy_baseline: "
                "cumulative=True builds its own slab grid and cannot also honour "
                "t_slab_edges.")
        if not np.all(energy == energy[0]):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_energy_baseline: "
                "cumulative=True scans baselines at one energy; the given energies differ. "
                "Baselines nest and energies do not, so there is nothing to reuse across "
                "energies.")
        if np.any(np.asarray(L, dtype=float) < L0):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_energy_baseline: "
                "cumulative=True requires every baseline to be at or beyond L0.")

        H_fixed = H_at_energy(energy[0])
        order = np.argsort(np.asarray(L, dtype=float), kind='stable')
        L_sorted = np.asarray(L, dtype=float)[order]

        # The accuracy grid is inherited, not invented.  One ordinary adaptive osc_prob call
        # at the longest baseline reports the slab count at which *it* converged -- which is
        # precisely "slabs needed for a uniform grid over the whole path", the definition of
        # n_acc -- and brings with it every safeguard that path already has, including the
        # n_slabs floor and the tolerance-not-achieved warning.  Getting this wrong is the
        # one way a cumulative scan goes silently wrong: on a solar profile an n_acc of 2000
        # is off by 1.6e-2 where 14883 is right, and nothing in the traversal itself notices.
        n_acc_from_ceiling = False
        if (rtol is None) and (atol is None):
            n_acc = kwargs.get('n_slabs', 1)
        else:
            probe_info = {}
            probe_kwargs = dict(osc_prob_kwargs)
            probe_kwargs['convergence_info'] = probe_info
            probe_kwargs.pop('A_eval_mode', None)
            # The probe is always strict, whatever the caller asked for their own points.  It is
            # the one call whose convergence decides the grid for the *whole* scan, so its
            # failure mode is not one bad point but N of them -- and the ladder's ordinary stop
            # rule can end on a coincidental agreement between two levels that are both wrong
            # (see the strict_convergence entry in osc_prob's docstring).  Measured on the solar
            # profile at 10 MeV, where that is exactly what the ordinary ladder does: the scan
            # came out at 5.2e-3 against a requested 1e-3 with a loose probe, and 1.0e-6 with a
            # strict one.  It costs one extra refinement level on a single call, amortised over
            # every baseline in the scan.
            probe_kwargs['strict_convergence'] = True
            # The probe's *probabilities* are discarded -- only its slab count is kept -- so a
            # MagnusConvergenceWarning about its intermediate refinement levels describes a
            # result nobody receives, and would be actively misleading: the grid this call
            # sizes produces no such warning when it is actually traversed.  Suppressed here
            # rather than globally, and only this one category: anything reporting that the
            # count itself is unreliable (ToleranceNotAchievedWarning) still reaches the caller,
            # because that does bear on the answer.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', magnus.MagnusConvergenceWarning)
                osc_prob(H_fixed, L0, float(L_sorted[-1]),
                         A_eval_mode=osc_prob_kwargs.get('A_eval_mode'), **probe_kwargs)
            # Scaled up because that count is what the *longest* baseline needed, and the
            # same uniform density is thinner than the shorter baselines in the scan would
            # have chosen for themselves; see CUMULATIVE_N_ACC_SAFETY for the measurement.
            n_acc = probe_info['n_slabs']*CUMULATIVE_N_ACC_SAFETY
            # Whether that count is a converged requirement or merely where the probe ran out
            # of room.  The probe's own ToleranceNotAchievedWarning says the *probe* did not
            # converge; it says nothing about the consequence that matters here, which is that
            # the grid for every baseline in the scan is now sized by a cap.  See the warning
            # a few lines below.
            n_acc_from_ceiling = bool(
                probe_info['n_slabs'] >= _resolve_max_n_slabs(max_n_slabs, integration_method))
            if n_acc_from_ceiling:
                # The probe raises ToleranceNotAchievedWarning of its own here, and that is a
                # statement about the probe -- one call whose probabilities are discarded.  The
                # consequence a caller needs is different and larger: the accuracy grid for
                # EVERY baseline in the scan is now sized by a cap rather than by a converged
                # requirement, so one capped call becomes N inaccurate answers.  Measured on a
                # solar profile, an n_acc of 2000 is off by 1.6e-2 where 14883 is right, and
                # nothing in the traversal itself notices.
                warnings.warn(
                    "osc_prob_energy_baseline (cumulative scan): the accuracy grid was sized "
                    "from max_n_slabs rather than from a converged requirement -- the probe "
                    "that sizes it ran out of slabs before it converged. Every baseline in "
                    "this scan inherits that grid, so the whole scan is affected, not one "
                    "point. Raise max_n_slabs (currently "
                    + str(_resolve_max_n_slabs(max_n_slabs, integration_method))
                    + "), or shorten the longest baseline, which is what sets the grid. Shown "
                      "once per session.",
                    ToleranceNotAchievedWarning, stacklevel=2)

        # A jump the caller did not declare is the one way this grid goes wrong that adding
        # slabs cannot fix: a slab straddling the discontinuity degrades the quadrature
        # regardless of magnus_exp_order, and refining only narrows the straddling slab.  Say
        # so rather than returning a quietly-wrong scan -- measured at 1.4e-03 and 2.1e-03,
        # silently, on two of 150 random piecewise profiles.  Skipped when the caller supplied
        # breakpoints, because then the grid already has edges on the discontinuities.
        # Only for a position-*dependent* Hamiltonian.  cumulative='auto' excludes a constant
        # one, but an explicit cumulative=True accepts it, and it then arrives here as a bare
        # matrix rather than a callable -- there is nothing to sample and nothing to be
        # discontinuous.  (Calling it anyway is a TypeError, which three existing
        # constant-Hamiltonian tests caught immediately.)
        if isinstance(H_fixed, Callable) and (
                (kwargs.get('t_breakpoints') is None)
                or (len(np.atleast_1d(kwargs.get('t_breakpoints'))) == 0)):
            if not (adiabatic._profile_is_resolved(H_fixed, float(L0), float(L_sorted[-1]), 200)
                    or adiabatic._profile_is_resolved(H_fixed, float(L0),
                                                      float(L_sorted[-1]), 6400)):
                warnings.warn(
                    "osc_prob_energy_baseline (cumulative scan): the Hamiltonian is "
                    "discontinuous at the scale of the grid this scan builds, and no "
                    "t_breakpoints were given. A slab straddling a density jump degrades the "
                    "quadrature no matter how many slabs are used; pass t_breakpoints at the "
                    "discontinuities. Shown once per session.",
                    UnmarkedDiscontinuityWarning, stacklevel=2)

        P_sorted = _osc_prob_cumulative_scan(
            H_fixed, L_sorted, L0, n_acc, magnus_exp_order,
            kwargs.get('n_tpts_per_slab', 100), integration_method,
            kwargs.get('t_breakpoints'), osc_prob_kwargs.get('A_eval_mode'),
            # strict_convergence is dropped rather than forwarded: the traversal walks a fixed
            # grid and runs no refinement ladder, so the flag has nothing to act on here, and
            # the engine below rejects unknown keywords.  The one adaptive step in a cumulative
            # scan is the probe above, which is unconditionally strict already.  Without this a
            # user who passes the flag to a baseline scan gets a TypeError out of
            # magnus_expansion_multislab.
            **{k: v for k, v in kwargs.items()
               if k not in ('n_slabs', 'n_tpts_per_slab', 't_breakpoints',
                            'strict_convergence')})

        _note_engine('cumulative', n_acc=int(n_acc), n_acc_from_ceiling=n_acc_from_ceiling)
        P_all = np.empty_like(P_sorted)
        P_all[order] = P_sorted
        if (nu_i is not None) and (nu_f is not None):
            P_all = P_all[:, nu_i, nu_f]
        return P_all

    # Warm starts: osc_prob reports the refinement parameters at which each point converged
    # (conv_info), and the next point starts its refinement from there (divided by one growth
    # factor, so that the comparison between successive refinements is still performed).
    # Neighboring points typically converge at (nearly) the same parameters, so this skips most
    # of the refinement ladder.
    conv_info = osc_prob_kwargs.get('convergence_info')
    if conv_info is None:
        conv_info = {}
    osc_prob_kwargs['convergence_info'] = conv_info
    warm_start = (t_slab_edges is None) and \
        ((rtol is not None) or (atol is not None))

    def apply_warm_start():
        # Seed the next point TWO growth steps below the last converged values.  One step below
        # reproduces exactly the pair of refinements at which the previous point was accepted;
        # starting one step lower than that lets the refinement scale decay geometrically across
        # points when the previous point was harder than the next ones (e.g., the lowest energy
        # of a scan), at the price of at most one extra refinement when it was not.
        if warm_start and conv_info:
            g1 = max(growth_factor_n_slabs, 1.0)**2
            g2 = max(growth_factor_n_tpts_per_slab, 1.0)**2
            osc_prob_kwargs['min_n_slabs'] = max(min_n_slabs,
                int(np.ceil(conv_info['n_slabs']/g1)))
            osc_prob_kwargs['min_n_tpts_per_slab'] = max(min_n_tpts_per_slab,
                int(np.ceil(conv_info['n_tpts_per_slab']/g2)))

    def compute_single_point(enu: float, baseline: float) -> Union[float, np.ndarray]:
        P = osc_prob(H_at_energy(enu), L0, baseline, **osc_prob_kwargs)
        # Select one oscillation channel if requested; otherwise return the full matrix
        if ((nu_i is not None) and (nu_f is not None)):
            return P[nu_i][nu_f]
        return P

    if parallelize_over_points:
        # Compute the first point serially to learn the refinement parameters, then distribute
        # the remaining points over the workers, warm-started from the first point.  (The shared
        # conv_info dict cannot be updated across processes, so it is dropped from the parallel
        # calls.)
        probs = [compute_single_point(energy[0], L[0])]
        apply_warm_start()
        osc_prob_kwargs.pop('convergence_info', None)
        probs += Parallel(n_jobs=n_jobs)(delayed(compute_single_point)(enu, baseline)
            for enu, baseline in zip(energy[1:], L[1:]))
    else:
        probs = []
        for enu, baseline in zip(energy, L):
            apply_warm_start()
            probs.append(compute_single_point(enu, baseline))

    # The private '_hamiltonian' payload is what lets cross_check_strategies build the 'expm'
    # reference without rebuilding any wrapper's physics: this is the one place every entry
    # point's Hamiltonian arrives already assembled.  Stripped from strategy_info.
    _note_engine('magnus', n_points=n_points, _hamiltonian=dict(
        H_at_energy=H_at_energy, L0=L0, energy=energy, L=L, nu_i=nu_i, nu_f=nu_f,
        t_breakpoints=kwargs.get('t_breakpoints')))
    # The call to __getitem__ below is a way to return a single float (or single probability
    # matrix) if both energy and L were given as floats.
    return np.array(probs).__getitem__(0 if return_float else slice(None))


#-----------------------------------------------------------------------
# Cross-method agreement
#-----------------------------------------------------------------------

_CROSS_CHECK_FORCING = {
    # label:      (strategy, cumulative, engines to forbid so that this one is reached)
    'hybrid':     ('hybrid', False, ('ip_exp', 'separable')),
    'ip_exp':     ('magnus', False, ('hybrid', 'separable')),
    'separable':  ('magnus', False, ('hybrid', 'ip_exp')),
    'cumulative': ('magnus', True, ('hybrid', 'ip_exp', 'separable')),
    'magnus':     ('magnus', False, ('hybrid', 'ip_exp', 'separable')),
}


def _expm_reference(payload: Dict) -> Tuple[Optional[np.ndarray], str]:
    r"""Oscillation probabilities from ``scipy.linalg.expm``, **only where it is exact**.

    :math:`U = \exp(-iH\,\Delta l)` solves :math:`dU/dl = -iH U` exactly when :math:`H` does not
    depend on :math:`l`; across a piecewise-constant profile whose edges are known, the
    time-ordered product of one exponential per piece is likewise exact.  Anywhere else it is a
    first-order approximation, which is worse than every engine it would be checking, so this
    declines instead of supplying a bad reference dressed as a good one.

    Constancy is *measured*, not assumed: ``H`` is sampled inside each piece and required to be
    constant to a relative :math:`10^{-12}`.  Sampling cannot prove constancy, but a profile that
    varies below that on 33 samples per piece and not between them is not one this test is
    protecting anybody from.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    payload : dict
        The ``_hamiltonian`` payload recorded by :func:`osc_prob_energy_baseline`:
        ``H_at_energy``, ``L0``, ``energy``, ``L``, ``nu_i``, ``nu_f``, ``t_breakpoints``.

    Returns
    -------
    (np.ndarray or None, str)
        Probabilities of shape ``(n_points, d, d)`` (or ``(n_points,)`` if a single channel was
        requested), and an empty reason; or ``None`` and the reason it declined.
    """
    from scipy.linalg import expm

    H_at_energy = payload['H_at_energy']
    L0 = float(payload['L0'])
    energy = np.asarray(payload['energy'], dtype=float)
    L = np.asarray(payload['L'], dtype=float)
    bp = payload['t_breakpoints']
    bp = (np.atleast_1d(np.asarray(bp, dtype=float)) if bp is not None
          else np.array([], dtype=float))

    P_out = []
    for enu, baseline in zip(energy, L):
        H_of_l = H_at_energy(float(enu))
        edges = np.unique(np.concatenate(
            [[L0, float(baseline)], bp[(bp > min(L0, baseline)) & (bp < max(L0, baseline))]]))
        U = None
        for a, b in zip(edges[:-1], edges[1:]):
            if isinstance(H_of_l, Callable):
                # Sample the open interval: the endpoints are exactly where a piecewise profile
                # is ambiguous, and a jump sitting on an edge is not a reason to decline.
                xs = np.linspace(a, b, 35)[1:-1]
                Hs = adiabatic._H_on_grid(H_of_l, xs)
                scale = np.max(np.abs(Hs))
                if scale > 0.0 and np.max(np.abs(Hs - Hs[0])) > 1.0e-12*scale:
                    return None, ('H varies with position on [%g, %g]; expm is exact only for '
                                  'a piecewise-constant H with declared edges' % (a, b))
                H_piece = Hs[len(Hs)//2]
            else:
                H_piece = np.asarray(H_of_l, dtype=complex)
            U_piece = expm(-1j*np.asarray(H_piece, dtype=complex)*(b - a))
            U = U_piece if U is None else U_piece @ U
        P_out.append(np.transpose(U.real**2 + U.imag**2))

    P_out = np.array(P_out)
    if (payload['nu_i'] is not None) and (payload['nu_f'] is not None):
        P_out = P_out[:, payload['nu_i'], payload['nu_f']]
    return P_out, ''


def cross_check_strategies(entry_point: Callable, *args, engines=None, **kwargs) -> Dict:
    r"""Answer the same request with every engine that applies, and report how far apart they are.

    **Why this exists.**  Every silently-wrong result found in
    ``docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md`` came from a method certifying itself by
    comparing itself with itself.  :func:`magnus.adiabatic.hybrid_propagator` refines its own
    knobs and checks the two answers agree; :func:`osc_prob`'s slab ladder does the same.  When
    the method has a blind spot, both sides of the comparison share it and the agreement carries
    no information -- that is not a bug in either comparison, it is a limit of the *shape* of the
    check.  This package contains genuinely different engines (see :data:`ENGINE_FAMILIES`), and
    running two of them needs **no oracle at all** while detecting exactly the class that
    self-certification cannot.

    **This is a diagnostic, not a safety net.**  It is never on by default: it multiplies the cost
    of a call by the number of engines that apply.  A large spread is *reported*, never raised --
    what it means depends on the request, and deciding that is the caller's job.

    **Acceptance.**  Measured by running the same comparison against the *pre-fix* package (a
    worktree at ``978663a``), on every construction of ``FINDINGS`` §3 that was silently wrong
    and reported ``certified=True`` there -- because a diagnostic validated only against code
    with no known defects has not been validated.  Reproduce with
    ``docs/dev/adversarial_batteries/crosscheck_acceptance.py``:

    ============================================ ============== ==============================
    construction                                 silent error   max cross-family spread
    ============================================ ============== ==============================
    step function, unmarked edge (§3.1)          5.395e-01      **5.399e-01**
    ten crossings (§3.2, worst found anywhere)   3.907e-02      **3.913e-02**
    sinusoid at span/7 (§3.2)                    1.672e-02      **1.687e-02**
    kink, :math:`C^0` but not :math:`C^1`        1.448e-02      **1.448e-02**
    singularity approached, not reached          8.625e-03      **8.613e-03**
    sub-threshold bump, w = 1e-2 span (§3.2)     7.701e-03      **7.768e-03**
    sub-threshold bump, w = 3e-2 span (§3.2)     4.388e-03      **4.594e-03**
    narrow bump, w = 3e-5 span (§3.3)            2.907e-02      3.5e-14 -- **not detected**
    ============================================ ============== ==============================

    Seven of eight, each at least four times the requested 1e-3.  The last row is the honest
    limit, and is why this table is here rather than a claim of coverage: a feature narrower than
    the probe spacing is invisible to *every* engine that samples the profile on a grid, so they
    agree -- correctly, given what they can see -- and are wrong together.  No cross-check
    between grid-based methods can find that; the cure is ``t_breakpoints`` at the feature (see
    :doc:`/adiabatic_strategy`).  Stated as a rule: **this sees a wrong engine exactly when some
    other engine got it right.**

    .. versionadded:: 1.0.0

    Parameters
    ----------
    entry_point : Callable
        The function to cross-check: :func:`osc_prob_matter_std_potential`,
        :func:`osc_prob_matter_nsi`, :func:`osc_prob_liv`, :func:`osc_prob_sun`,
        :func:`osc_prob_earth`, one of their fixed-flavour wrappers, or
        :func:`osc_prob_energy_baseline`.  ``strategy`` and ``cumulative`` are supplied by this
        function and dropped from ``kwargs`` if the caller passed them, since forcing them is
        how each engine is reached.
    \*args
        Positional arguments for ``entry_point``, exactly as in an ordinary call.
    engines : sequence of str, optional
        Restrict the check to these engines; see :data:`ENGINE_FAMILIES` for the labels.  Default:
        every engine that applies.
    \**kwargs
        Keyword arguments for ``entry_point``, exactly as in an ordinary call.

    Returns
    -------
    dict
        ``'answers'``
            ``{label: probabilities}`` for each engine that answered.
        ``'ran'``
            The labels that answered, in the order tried.
        ``'declined'``
            ``{label: reason}`` for each engine that did not.  Most engines decline on most
            requests, which is expected and not a finding.
        ``'spread'``
            ``{(label_a, label_b): max |P_a - P_b|}`` over every pair that ran.
        ``'max_spread'``, ``'max_spread_pair'``
            The largest spread and where it was.
        ``'max_spread_independent'``, ``'max_spread_independent_pair'``
            The same, restricted to pairs from **different** families (see
            :data:`ENGINE_FAMILIES`).  This is the number to read: two engines from the same
            family can be wrong in the same way, so their disagreement is informative but their
            *agreement* is not.
        ``'families'``
            ``{label: family}`` for the engines that ran.
        ``'warnings'``
            ``{label: [warning class names]}`` raised while that engine answered.
        ``'certified'``
            ``{label: bool}``, currently only for ``'hybrid'``.

    Examples
    --------
    >>> import numpy as np
    >>> import magnus.globaldefs as gd
    >>> import magnus.matter as matter
    >>> import magnus.oscprob as oscprob
    >>> ne = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    >>> params = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    >>> out = oscprob.cross_check_strategies(
    ...     oscprob.osc_prob_matter_std_potential, 2, ne, 10.0e6,
    ...     0.5*gd.SUN_RADIUS*gd.UNIT_KM,
    ...     {'s12': params['s12'], 'Dm2': params['D21']}, L0=0.0,
    ...     density_is_of_number_of_electrons=True)
    >>> sorted(out['ran'])                                    # doctest: +SKIP
    ['cumulative', 'hybrid', 'ip_exp', 'magnus']
    >>> out['max_spread_independent'] < 1.0e-3                # doctest: +SKIP
    True

    See Also
    --------
    ENGINE_FAMILIES : which engines share machinery, and therefore which pairs carry information.
    """
    wanted = tuple(_CROSS_CHECK_FORCING) + ('expm',) if engines is None else tuple(engines)
    unknown = set(wanted) - set(_CROSS_CHECK_FORCING) - {'expm'}
    if unknown:
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.cross_check_strategies: unknown "
            "engine label(s) " + repr(sorted(unknown)) + "; known labels are "
            + repr(sorted(set(_CROSS_CHECK_FORCING) | {'expm'})) + ".")

    call_kwargs = {k: v for k, v in kwargs.items() if k not in ('strategy', 'cumulative')}
    takes_strategy = 'strategy' in signature(entry_point).parameters

    answers, declined, warns, certified = {}, {}, {}, {}
    hamiltonian = None

    for label in wanted:
        if label == 'expm':
            continue
        strategy, cumulative, forbid = _CROSS_CHECK_FORCING[label]
        forced = dict(call_kwargs, cumulative=cumulative)
        if takes_strategy:
            forced['strategy'] = strategy
        elif strategy == 'hybrid':
            declined[label] = "entry point has no 'strategy' parameter"
            continue
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                with _engine_probe(forbid) as trace:
                    P = entry_point(*args, **forced)
            except Exception as exc:               # noqa: BLE001 -- a decline, not a failure
                declined[label] = type(exc).__name__ + ': ' + str(exc).strip().split('\n')[0]
                continue
        note = next((e for e in trace if e['engine'] == label and e['answered']), None)
        if note is None:
            other = next((e['engine'] for e in trace if e['answered']), 'nothing')
            refused = next((e.get('reason') for e in trace
                            if e['engine'] == label and not e['answered']), None)
            declined[label] = refused or ('does not apply to this request (answered by '
                                          + other + ')')
            continue
        answers[label] = np.asarray(P)
        warns[label] = sorted({w.category.__name__ for w in caught})
        if 'certified' in note:
            certified[label] = note['certified']
        if hamiltonian is None:
            hamiltonian = next((e['_hamiltonian'] for e in trace if '_hamiltonian' in e), None)

    if 'expm' in wanted:
        if hamiltonian is None:
            declined['expm'] = ('the Hamiltonian was not observed; run the general Magnus '
                                'engine as well (it is where the Hamiltonian is assembled)')
        else:
            P_expm, reason = _expm_reference(hamiltonian)
            if P_expm is None:
                declined['expm'] = reason
            else:
                answers['expm'] = P_expm
                warns['expm'] = []

    ran = [lab for lab in wanted if lab in answers]
    spread, best, best_pair, best_ind, best_ind_pair = {}, 0.0, None, 0.0, None
    for i, a in enumerate(ran):
        for b in ran[i + 1:]:
            Pa, Pb = np.ravel(answers[a]), np.ravel(answers[b])
            if Pa.shape != Pb.shape:
                # Different shapes mean the engines answered different questions; that is a
                # defect in this diagnostic's forcing, not a physics finding, so say so rather
                # than broadcasting them into a number.
                spread[(a, b)] = np.nan
                continue
            s = float(np.max(np.abs(Pa - Pb)))
            spread[(a, b)] = s
            if s > best:
                best, best_pair = s, (a, b)
            if ENGINE_FAMILIES[a] != ENGINE_FAMILIES[b] and s > best_ind:
                best_ind, best_ind_pair = s, (a, b)

    return {
        'answers': answers,
        'ran': tuple(ran),
        'declined': declined,
        'spread': spread,
        'max_spread': best,
        'max_spread_pair': best_pair,
        'max_spread_independent': best_ind,
        'max_spread_independent_pair': best_ind_pair,
        'families': {lab: ENGINE_FAMILIES[lab] for lab in ran},
        'warnings': warns,
        'certified': certified,
    }


#-----------------------------------------------------------------------
# General functions for vacuum, standard matter, NSI, LIV
#-----------------------------------------------------------------------

def osc_prob_vacuum(
    num_flavors: int,
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    average: Optional[bool]=False,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='gl', 
    rtol: Optional[Union[int, float]]=1.e-3, 
    atol: Optional[Union[int, float]]=1.e-3, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=None, 
    min_n_tpts_per_slab: Optional[int]=2, 
    max_n_tpts_per_slab: Optional[int]=500, 
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for
    oscillations in vacuum

    Middle (scenario) layer for the vacuum case, generic in ``num_flavors``: unpacks
    ``osc_params``, builds the energy-independent vacuum Hamiltonian via
    ``hamiltonians.hamiltonian_{num_flavors}nu_vacuum_energy_independent``, and calls
    :func:`osc_prob_energy_baseline`. Called by :func:`osc_prob_2nu_vacuum`,
    :func:`osc_prob_3nu_vacuum`, :func:`osc_prob_4nu_vacuum`, and :func:`osc_prob_5nu_vacuum`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_vac_energy_indep`` is given).
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    osc_params : dict
        Oscillation parameters; see :func:`unpack_oscillation_params_from_dict` for the required
        keys for each ``num_flavors``.
    h_vac_energy_indep : list or np.ndarray, optional
        Precomputed energy-independent vacuum Hamiltonian, used instead of ``osc_params`` when
        ``num_flavors`` exceeds ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any parameter left as
        None in ``osc_params``. Default: 'OSC_PARAMS_DEFAULT'.
    t_slab_edges : list or np.ndarray, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    magnus_exp_order : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    n_jobs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    integration_method : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    rtol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    atol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_slabs : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_tpts_per_slab : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_num_loops : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    validate_input : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    save_log : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    filename_log : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    file_log : TextIOWrapper, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    close_file_log_upon_exit : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    verbose : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`.

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each (energy, L) point.
    """

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name, 
        num_flavors, osc_params, h_vac_energy_indep)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list

    if validate_input:
        validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=0.0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            validate_energy_and_L=True, validate_flavor_indices=True, validate_osc_params=True, 
            validate_initial_position=False, validate_density=False)

    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar) 

    def htot(enu: Union[int, float]) -> np.ndarray:
        return (1/enu)*h_vac_energy_indep

    htot_is_function_only_of_energy = True

    # Checked here, rather than only in osc_prob: the averaged path returns before anything
    # forwards **kwargs onwards, so a check further down would see these keys on the ordinary
    # path and silently ignore them on the averaged one.
    _reject_parameter_set_metadata(kwargs, 'osc_prob_vacuum')

    # Phase-averaged limit, requested with average=True: exact and closed-form whenever the
    # Hamiltonian does not depend on position, so it is tried before any of the propagation
    # machinery below, all of which would resolve phases that the average discards (see
    # _avg_prob_dispatch and :mod:`magnus.avgprob`).
    P_avg = _avg_prob_dispatch(htot, htot_is_function_only_of_energy, energy, L, 0.0, nu_i, nu_f,
        average, 'osc_prob_vacuum')
    if P_avg is not NotImplemented:
        return P_avg

    # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).  (The
    # Hamiltonian is constant in position, so osc_prob computes each point exactly with a single
    # slab; the tolerance and refinement parameters play no role and are not forwarded.)
    return osc_prob_energy_baseline(htot, energy, L, 0.0, nu_i, nu_f,
        htot_is_function_only_of_energy, n_jobs=n_jobs, validate_input=validate_input,
        verbose=verbose, **kwargs)


def osc_prob_matter_std_potential(
    num_flavors: int,
    rho_func: Union[Callable, int, float],
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    L0: Optional[Union[int, float]]=0.0,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    average: Optional[bool]=False,
    strategy: Optional[str]='auto',
    strategy_info: Optional[Dict]=None,
    t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='gl',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5,
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5,
    max_num_loops: Optional[int]=50,
    min_n_slabs: Optional[int]=1,
    max_n_slabs: Optional[int]=None,
    min_n_tpts_per_slab: Optional[int]=2,
    max_n_tpts_per_slab: Optional[int]=500,
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    new_recursion_limit: Optional[int]=5000,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for
    standard oscillations in matter, i.e., the matter potential is only
    due to the coherent forward scattering of nu_e on electrons.

    Middle (scenario) layer for the standard-matter case, generic in ``num_flavors``: unpacks
    ``osc_params``, builds the vacuum + matter Hamiltonian (via
    ``hamiltonians.hamiltonian_{num_flavors}nu_vacuum_energy_independent`` and
    ``hamiltonian_{num_flavors}nu_matter_td``, with the potential from
    :func:`magnus.matter.vcc_func_from_rho_func`), and calls
    :func:`osc_prob_energy_baseline`. Called by every
    ``osc_prob_{2,3,4,5}nu_matter_{constant,exp}_density`` and
    ``osc_prob_{2,3,4,5}nu_earth``/``osc_prob_{2,3,4,5}nu_sun`` wrapper.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_vac_energy_indep`` is given).
    rho_func : Callable, int, or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is
        True), either as a function of position or as a constant.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    osc_params : dict
        Oscillation parameters; see :func:`unpack_oscillation_params_from_dict`.
    L0 : int or float, optional
        Initial position. Default: 0.0.
    h_vac_energy_indep : list or np.ndarray, optional
        Precomputed energy-independent vacuum Hamiltonian, used instead of ``osc_params`` when
        ``num_flavors`` exceeds ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos (flips the sign of the matter
        potential). Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``rho_func`` returns the density in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, ``rho_func`` directly returns the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any parameter left as
        None in ``osc_params``. Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid',
        or 'magnus'.

        * ``'magnus'`` uses only the traditional Magnus-expansion machinery (the closed-form
          two-flavor interaction-picture integrator when it applies, the energy-batched scan
          engine, or the general adaptive slab-refinement method) -- this reproduces the exact
          behavior of Magνs as it was before the adiabatic strategy was added,
          unconditionally.  It therefore also opts out of the cumulative baseline scan, which
          postdates that behavior: pass ``strategy='magnus'`` to reproduce older numbers
          exactly, on a baseline scan as well as at a single point.
        * ``'hybrid'`` additionally tries :func:`magnus.adiabatic.hybrid_propagator` (adiabatic
          transport, with a Magnus patch at any non-adiabatic window; see
          :doc:`/adiabatic_strategy`) for any requested (energy, L) point where ``rho_func`` is
          position-dependent and no ``t_slab_edges``/breakpoints are given, and a target
          tolerance (``rtol``/``atol``) is requested. If it fails to self-certify for at least
          one point, the best-effort result is still returned, together with
          :class:`HybridCertificationWarning`.
        * ``'auto'`` tries the hybrid strategy first, under the same conditions, but falls back
          silently to the ``'magnus'`` strategies above (no warning about the hybrid attempt
          itself) for any point where it does not apply or fails to self-certify.  It also
          stands aside for a **baseline scan at a single energy** of at least
          ``HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`` points, which the cumulative scan
          (see ``cumulative`` in :func:`osc_prob_energy_baseline`) answers from one traversal
          instead of one hybrid call per point -- measured on solar profiles as tens of times
          faster at equal or better accuracy.

        The hybrid strategy is the natural tool exactly where the plain Magnus refinement needs
        very many slabs (an extreme accumulated phase, e.g., low-energy solar neutrinos crossing
        an MSW resonance), and applies to any number of flavors and to genuinely complex
        Hamiltonians; see :doc:`/adiabatic_strategy` for the full derivation, validation, and
        performance comparison. Default: 'auto'.

        .. warning::

           **If your profile has a feature much narrower than the trajectory, say where it is.**
           The hybrid strategy locates resonances by sampling ``n_probe`` points (200, refined to
           at most 6400), and the general Magnus path seeds its grid from an integral along the
           path; neither can see a feature that falls between samples, and no refinement of
           either finds it, because refinement never puts a point inside it. Measured on a
           Gaussian resonance of width :math:`10^{-5}` of the trajectory, the returned
           probability was wrong by **2.9e-02 against a requested 1e-3** -- on the hybrid path,
           the general path, and the cumulative scan alike.

           It is **no longer silent**: :func:`magnus.adiabatic.find_hidden_features` scans the
           profile itself, once per call, and raises :class:`HiddenFeatureWarning` naming the
           position and the ``t_breakpoints`` to pass.  It reaches this class precisely because
           it looks at the profile rather than at the answers, which is what no comparison
           between engines can do when they are all wrong together.  Detection is 68-90 % over
           the unresolvable band with 0 false positives on 67 smooth profiles, so it is a
           report -- not a guarantee, and not a cure.

           Passing ``t_breakpoints`` at the feature fixes it, and is tested: the same case goes
           to 8.8e-04 at a single point and 8.9e-04 over a 60-point scan. It is the right tool
           twice over, since an edge placed *on* a sharp feature also stops a slab straddling it
           from degrading the quadrature. This is the one exposure the adversarial validation
           (``docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md``) could not close in the library
           itself: what a fixed grid never samples, it cannot report.

        .. versionadded:: 1.0.0
    strategy_info : dict, optional
        If given, filled in place with which engine actually answered, following the same
        out-parameter convention as ``convergence_info`` in :func:`osc_prob`.  Under
        ``strategy='auto'`` the fallbacks are silent by design -- that is right for ordinary
        calls, and wrong for anyone asking why a result moved or why a call got slow -- so this
        is how to see them without turning the fallbacks into warnings.  Keys:

        * ``'engine'`` -- ``'hybrid'``, ``'ip_exp'``, ``'separable'``, ``'cumulative'``,
          ``'magnus'`` or ``'average'``.
        * ``'family'`` -- the engine's family; see :data:`ENGINE_FAMILIES`.
        * ``'certified'`` -- for ``'hybrid'``, whether
          :func:`magnus.adiabatic.hybrid_propagator` self-certified.  ``None`` for engines
          that do not certify.  Under ``'auto'`` an uncertified hybrid result is never
          returned, so this is ``True`` whenever the engine is ``'hybrid'``; under
          ``'hybrid'`` it can be ``False``, and then it means the accuracy is **unverified**,
          not that the answer is wrong.
        * ``'declined'`` -- ``[(engine, reason)]`` for the engines that stood aside first.
          Most requests decline most engines, which is ordinary and not a finding.
        * ``'trace'`` -- every dispatch decision in order, with per-engine detail (for the
          cumulative scan, ``'n_acc'`` and whether it came from a ceiling).

        Costs nothing when omitted.  Default: None.

        .. versionadded:: 1.0.0
    t_slab_edges : list or np.ndarray, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    magnus_exp_order : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    n_jobs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    integration_method : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    rtol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    atol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_slabs : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_tpts_per_slab : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_num_loops : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    validate_input : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    save_log : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    filename_log : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    file_log : TextIOWrapper, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    close_file_log_upon_exit : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    verbose : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    new_recursion_limit : int, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`.

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each (energy, L) point.
    """

    if validate_input and (strategy not in ('auto', 'hybrid', 'magnus')):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_std_potential:" + \
            " strategy must be 'auto', 'hybrid', or 'magnus'.")

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, osc_params, h_vac_energy_indep)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list

    if validate_input:
        validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=L0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list,
            rho_func=rho_func, ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction, validate_energy_and_L=True,
            validate_flavor_indices=True, validate_osc_params=True, validate_initial_position=True,
            validate_density=True)

    # If any of the standard oscillation parameters has not been given a value, assign to it the
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21,
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)

    # Build the coherent forward potential function, VCC_func, from the density function, rho_func.
    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3].
    VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
        electron_fraction, nubar, density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons) # [eV]

    # Projector onto the nu_e--nu_e entry, multiplied below by the potential VCC.  Note that
    # VCC_func already carries the antineutrino sign flip (applied inside
    # matter.vcc_func_from_rho_func), so no extra sign is applied here.  [Previously, the sign was
    # applied twice, which gave the antineutrino matter potential the wrong (positive) sign.]
    h_matt_proj = np.zeros((num_flavors, num_flavors))
    h_matt_proj[0][0] = 1.0

    # Cache repeated evaluations of the potential on identical position grids (see
    # _PositionProfileCache)
    if isinstance(VCC_func, Callable):
        VCC_func = _PositionProfileCache(VCC_func)

    # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
    if isinstance(VCC_func, Callable):
        # VCC_func is a function of position, so the Hamiltonian is, too.  If l is an array, the
        # result is a stack of Hamiltonians with the position axis leading; this lets the Magnus
        # routines evaluate the Hamiltonian at all time points in a single vectorized call.
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            vcc = np.asarray(VCC_func(l))
            return (1/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt_proj
        htot_is_function_only_of_energy = False
    else:
        # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is passed to
        # osc_prob below, osc_prob will detect that VCC_func is constant and set parameters
        # internally for speed-up.
        h_matt = VCC_func*h_matt_proj
        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = True


    # Resolved -- and, crucially, POPPED out of kwargs -- before scan_kwargs is built, because
    # every dispatcher below declines outright on an unrecognized entry in kwargs.  Left in, an
    # explicitly-passed cumulative silently disabled the hybrid, interaction-picture and
    # separable engines: passing cumulative='auto', which is the documented default, changed
    # which engine answered (hybrid -> general ladder) and moved a 10 MeV solar single point by
    # 9.3e-06.  See _resolve_cumulative_kwarg.
    cumulative_resolved = _resolve_cumulative_kwarg(kwargs, strategy)

    scan_kwargs = dict(t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
        integration_method=integration_method, rtol=rtol, atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        save_log=save_log, file_log=file_log, cumulative=cumulative_resolved,
        kwargs=kwargs)

    # Checked here, rather than only in osc_prob: the averaged path returns before anything
    # forwards **kwargs onwards, so a check further down would see these keys on the ordinary
    # path and silently ignore them on the averaged one.
    _reject_parameter_set_metadata(kwargs, 'osc_prob_matter_std_potential')

    # Phase-averaged limit, requested with average=True: exact and closed-form whenever the
    # Hamiltonian does not depend on position, so it is tried before any of the propagation
    # machinery below, all of which would resolve phases that the average discards (see
    # _avg_prob_dispatch and :mod:`magnus.avgprob`).
    # `or` on an array would ask for its truth value, so the emptiness test is explicit.
    _breakpoints = kwargs.get('t_breakpoints')
    _profile_is_smooth = ((_breakpoints is None or len(np.atleast_1d(_breakpoints)) == 0)
                          and (t_slab_edges is None))
    # One scan of the profile for the whole call, before any engine sees it: a feature
    # narrower than every grid here is the one exposure none of them can detect for itself,
    # because they all miss it together.  Depends on the profile and the interval, never on
    # energy, so it is not repeated per point.  See _scan_for_hidden_features.
    _hidden = _scan_for_hidden_features(VCC_func, L0, L, kwargs.get('t_breakpoints'))

    # Everything below is dispatch: which engine gets the request.  Watched as a unit so
    # that strategy_info can report which one answered and, for the hybrid strategy,
    # whether it certified -- see _engine_probe.  Costs one list allocation per call when
    # nobody is watching.
    with _engine_probe(info=strategy_info, extra={'hidden_feature': _hidden}):
        P_avg = _avg_prob_dispatch(htot, htot_is_function_only_of_energy, energy, L, L0, nu_i, nu_f,
            average, 'osc_prob_matter_std_potential', smooth_profile=_profile_is_smooth, engine_kwargs=scan_kwargs)
        if P_avg is not NotImplemented:
            return P_avg

        # Hybrid strategy (adiabatic transport + Magnus patch at any non-adiabatic window; see
        # _osc_prob_hybrid_dispatch and :doc:`/adiabatic_strategy`): broader than the interaction-
        # picture fast path below (any number of flavors, any smooth position-dependent profile), and
        # tried first unless strategy == 'magnus'.  Falls back transparently (returns NotImplemented)
        # if it does not apply or (with strategy == 'auto' only) fails to self-certify.
        #
        # Tried *before* the fast path because that is what strategy='auto' has always been
        # documented to mean ("tries the hybrid strategy first ... but falls back silently to the
        # 'magnus' strategies", of which the interaction-picture integrator is one) -- and because
        # measurement says the documented order is also the better one.  On solar configurations,
        # across 50 (energy, baseline) points spanning the standard, NSI and LIV families and
        # 0.5-100 MeV, scored against solve_ivp/DOP853: the hybrid strategy certified 50/50 with a
        # worst error of 1.8e-04 against a requested 1e-3 and no warnings, while the fast path
        # certified 22/50 and took a mean of 13.2 s to decline the other 28.  Where both answer,
        # hybrid is 28-594x faster (median 397x).  The fast path is more accurate only at 40-100 MeV
        # and only at the default tolerance -- at rtol/atol <= 1e-5 it declines outright at every
        # energy measured.  See docs/dev/DECISION_DISPATCH_ORDER.md.
        P_hybrid = _osc_prob_hybrid_dispatch(h_vac_energy_indep, VCC_func, h_matt_proj, None, None,
            energy, L, L0, nu_i, nu_f, scan_kwargs, strategy)
        if P_hybrid is not NotImplemented:
            return P_hybrid

        # Fast path for a genuine exponential density profile (e.g., the Sun): factor out the
        # (possibly huge, at low energy) fast vacuum phase analytically in the interaction picture,
        # instead of resolving it slab by slab (see _osc_prob_ip_exp_dispatch).  Applies to a single
        # (energy, L) point as well as to a scan, and falls back transparently (returns
        # NotImplemented) if the profile is not exponential or if it fails to converge (e.g., near an
        # MSW resonance), in which case the general methods below are used instead.  Reached only
        # where the hybrid strategy declined, or with strategy == 'magnus'.
        P_ip = _osc_prob_ip_exp_dispatch(h_vac_energy_indep, VCC_func, h_matt_proj, None, None,
            energy, L, L0, nu_i, nu_f, scan_kwargs)
        if P_ip is not NotImplemented:
            return P_ip

        # Energy-batched fast path: when many energies share a single baseline and the Hamiltonian
        # is position-dependent, compute the whole scan in one batched pipeline, with the potential
        # samples shared across energies (see _osc_prob_scan_separable).  If the request does not fit
        # the engine, fall back to the generic per-point path below.
        P_scan = _osc_prob_scan_separable_dispatch(h_vac_energy_indep, VCC_func, h_matt_proj, None, None,
            energy, L, L0, nu_i, nu_f, scan_kwargs)
        if P_scan is not NotImplemented:
            return P_scan

        # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
        return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f,
            htot_is_function_only_of_energy, t_slab_edges=t_slab_edges,
            magnus_exp_order=magnus_exp_order, n_jobs=n_jobs, integration_method=integration_method,
            rtol=rtol, atol=atol, growth_factor_n_slabs=growth_factor_n_slabs,
            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
            min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
            new_recursion_limit=new_recursion_limit, verbose=verbose,
            # An explicit cumulative= from the caller wins; otherwise strategy='magnus' opts out
            # of the cumulative scan and everything else takes 'auto'.  Resolved near the top of
            # this function, which is also where it is removed from kwargs.
            cumulative=cumulative_resolved, **kwargs)


def osc_prob_matter_nsi(
    num_flavors: int,
    rho_func: Union[Callable, int, float],
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    nsi_params: Dict,
    L0: Optional[Union[int, float]]=0.0,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    h_nsi: Union[list, np.ndarray]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    average: Optional[bool]=False,
    strategy: Optional[str]='auto',
    strategy_info: Optional[Dict]=None,
    t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='gl',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5,
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5,
    max_num_loops: Optional[int]=50,
    min_n_slabs: Optional[int]=1,
    max_n_slabs: Optional[int]=None,
    min_n_tpts_per_slab: Optional[int]=2,
    max_n_tpts_per_slab: Optional[int]=500,
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    new_recursion_limit: Optional[int]=5000,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for
    oscillations in matter with non-standard interactions (NSI), i.e., the matter potential
    includes both the standard coherent-forward-scattering term and the NSI epsilon couplings.

    Middle (scenario) layer for the NSI case, generic in ``num_flavors``: unpacks ``osc_params``
    and ``nsi_params``, builds the vacuum + matter + NSI Hamiltonian (via
    ``hamiltonians.hamiltonian_{num_flavors}nu_vacuum_energy_independent``,
    ``hamiltonian_{num_flavors}nu_matter_td``, and ``hamiltonian_{num_flavors}nu_nsi_td``, with
    the potential from :func:`magnus.matter.vcc_func_from_rho_func`), and calls
    :func:`osc_prob_energy_baseline`. Called by every
    ``osc_prob_{2,3,4,5}nu_matter_nsi_{constant,exp}_density`` and
    ``osc_prob_{2,3,4,5}nu_earth_nsi``/``osc_prob_{2,3,4,5}nu_sun_nsi`` wrapper.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_vac_energy_indep``/``h_nsi``
        are given).
    rho_func : Callable, int, or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is
        True), either as a function of position or as a constant.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    osc_params : dict
        Oscillation parameters; see :func:`unpack_oscillation_params_from_dict`.
    nsi_params : dict
        NSI epsilon parameters; see :func:`unpack_nsi_params_from_dict`.
    L0 : int or float, optional
        Initial position. Default: 0.0.
    h_vac_energy_indep : list or np.ndarray, optional
        Precomputed energy-independent vacuum Hamiltonian, used instead of ``osc_params`` when
        ``num_flavors`` exceeds ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.
    h_nsi : list or np.ndarray, optional
        Precomputed NSI Hamiltonian, used instead of ``nsi_params`` when ``num_flavors`` exceeds
        ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos (flips the sign of the matter
        potential). Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``rho_func`` returns the density in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, ``rho_func`` directly returns the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any parameter left as
        None in ``osc_params``. Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid',
        or 'magnus'; see the ``strategy`` parameter of :func:`osc_prob_matter_std_potential` for
        the full description and :doc:`/adiabatic_strategy` for the derivation and validation of
        the ``'hybrid'``/``'auto'`` strategies (adiabatic transport with a Magnus patch at any
        non-adiabatic window, applicable to any number of flavors). Default: 'auto'.

        .. versionadded:: 1.0.0
    strategy_info : dict, optional
        If given, filled in place with which engine actually answered, following the same
        out-parameter convention as ``convergence_info`` in :func:`osc_prob`.  Under
        ``strategy='auto'`` the fallbacks are silent by design -- that is right for ordinary
        calls, and wrong for anyone asking why a result moved or why a call got slow -- so this
        is how to see them without turning the fallbacks into warnings.  Keys:

        * ``'engine'`` -- ``'hybrid'``, ``'ip_exp'``, ``'separable'``, ``'cumulative'``,
          ``'magnus'`` or ``'average'``.
        * ``'family'`` -- the engine's family; see :data:`ENGINE_FAMILIES`.
        * ``'certified'`` -- for ``'hybrid'``, whether
          :func:`magnus.adiabatic.hybrid_propagator` self-certified.  ``None`` for engines
          that do not certify.  Under ``'auto'`` an uncertified hybrid result is never
          returned, so this is ``True`` whenever the engine is ``'hybrid'``; under
          ``'hybrid'`` it can be ``False``, and then it means the accuracy is **unverified**,
          not that the answer is wrong.
        * ``'declined'`` -- ``[(engine, reason)]`` for the engines that stood aside first.
          Most requests decline most engines, which is ordinary and not a finding.
        * ``'trace'`` -- every dispatch decision in order, with per-engine detail (for the
          cumulative scan, ``'n_acc'`` and whether it came from a ceiling).

        Costs nothing when omitted.  Default: None.

        .. versionadded:: 1.0.0
    t_slab_edges : list or np.ndarray, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    magnus_exp_order : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    n_jobs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    integration_method : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    rtol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    atol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_slabs : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_tpts_per_slab : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_num_loops : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    validate_input : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    save_log : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    filename_log : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    file_log : TextIOWrapper, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    close_file_log_upon_exit : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    verbose : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    new_recursion_limit : int, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`.

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each (energy, L) point.
    """

    if validate_input and (strategy not in ('auto', 'hybrid', 'magnus')):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_matter_nsi: strategy" + \
            " must be 'auto', 'hybrid', or 'magnus'.")

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, osc_params, h_vac_energy_indep)
    nsi_params_list = unpack_nsi_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, nsi_params, h_nsi)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
        eps_aa, eps_ab = nsi_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
        eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = nsi_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
        eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss = \
            nsi_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list
        eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, \
            eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2 = nsi_params_list

    if validate_input:
        validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=L0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            rho_func=rho_func, ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction, validate_energy_and_L=True, 
            validate_flavor_indices=True, validate_osc_params=True, validate_initial_position=True,
            validate_density=True)

    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)

    # Compute the standard + NSI matter Hamiltonian *without* the multiplicative prefactor of VCC.
    # To do this we call the functions hamiltonians_Xnu_nsi(VCC, ...) with VCC = 1.0.  We add the
    # standard matter contribution to the NSI matter contribution by adding 1.0 to the eps_ee entry.
    # The overall antineutrino sign flip is carried by VCC_func (see
    # matter.vcc_func_from_rho_func); for antineutrinos, the NSI couplings are additionally
    # conjugated (H_matt -> -H_matt^* relative to neutrinos).
    if num_flavors == 2:
        h_matt = np.diag([1.0, 0.0]) + \
            hamiltonians.hamiltonian_2nu_nsi(1.0, eps_aa, eps_ab) # VCC = 1.0
    elif num_flavors == 3:
        h_matt = np.diag([1.0, 0.0, 0.0]) + \
            hamiltonians.hamiltonian_3nu_nsi(1.0, eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt)
    elif num_flavors == 4:
        h_matt = np.diag([1.0, 0.0, 0.0, 0.0]) + \
            hamiltonians.hamiltonian_4nu_nsi(1.0, eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt,
                eps_ms, eps_tt, eps_ts, eps_ss)
    elif num_flavors == 5:
        h_matt = np.diag([1.0, 0.0, 0.0, 0.0, 0.0]) + \
            hamiltonians.hamiltonian_5nu_nsi(1.0, eps_ee, eps_em, eps_et, eps_es1, eps_es2,
                eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2,
                eps_s2s2)

    if nubar:
        h_matt = np.conj(h_matt)

    # Build the coherent forward potential function, VCC_func, from the density function, rho_func.
    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3].
    VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
        electron_fraction, nubar, density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons) # [eV] 
    
    # Cache repeated evaluations of the potential on identical position grids (see
    # _PositionProfileCache)
    if isinstance(VCC_func, Callable):
        VCC_func = _PositionProfileCache(VCC_func)

    # Matter Hamiltonian function: (standard + NSI) matter matrix scaled by VCC
    if isinstance(VCC_func, Callable):
        # VCC_func is a function of position, so the Hamiltonian is, too.  If l is an array, the
        # result is a stack of Hamiltonians with the position axis leading; this lets the Magnus
        # routines evaluate the Hamiltonian at all time points in a single vectorized call.
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            vcc = np.asarray(VCC_func(l))
            return (1/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt
        htot_is_function_only_of_energy = False
    else:
        # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is passed to
        # osc_prob below, osc_prob will detect that VCC_func is constant and set parameters
        # internally for speed-up.
        h_matt = VCC_func*h_matt
        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = True


    # Resolved -- and, crucially, POPPED out of kwargs -- before scan_kwargs is built, because
    # every dispatcher below declines outright on an unrecognized entry in kwargs.  Left in, an
    # explicitly-passed cumulative silently disabled the hybrid, interaction-picture and
    # separable engines: passing cumulative='auto', which is the documented default, changed
    # which engine answered (hybrid -> general ladder) and moved a 10 MeV solar single point by
    # 9.3e-06.  See _resolve_cumulative_kwarg.
    cumulative_resolved = _resolve_cumulative_kwarg(kwargs, strategy)

    scan_kwargs = dict(t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
        integration_method=integration_method, rtol=rtol, atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        save_log=save_log, file_log=file_log, cumulative=cumulative_resolved,
        kwargs=kwargs)

    # Checked here, rather than only in osc_prob: the averaged path returns before anything
    # forwards **kwargs onwards, so a check further down would see these keys on the ordinary
    # path and silently ignore them on the averaged one.
    _reject_parameter_set_metadata(kwargs, 'osc_prob_matter_nsi')

    # Phase-averaged limit, requested with average=True: exact and closed-form whenever the
    # Hamiltonian does not depend on position, so it is tried before any of the propagation
    # machinery below, all of which would resolve phases that the average discards (see
    # _avg_prob_dispatch and :mod:`magnus.avgprob`).
    # `or` on an array would ask for its truth value, so the emptiness test is explicit.
    _breakpoints = kwargs.get('t_breakpoints')
    _profile_is_smooth = ((_breakpoints is None or len(np.atleast_1d(_breakpoints)) == 0)
                          and (t_slab_edges is None))
    # One scan of the profile for the whole call, before any engine sees it: a feature
    # narrower than every grid here is the one exposure none of them can detect for itself,
    # because they all miss it together.  Depends on the profile and the interval, never on
    # energy, so it is not repeated per point.  See _scan_for_hidden_features.
    _hidden = _scan_for_hidden_features(VCC_func, L0, L, kwargs.get('t_breakpoints'))

    # Everything below is dispatch: which engine gets the request.  Watched as a unit so
    # that strategy_info can report which one answered and, for the hybrid strategy,
    # whether it certified -- see _engine_probe.  Costs one list allocation per call when
    # nobody is watching.
    with _engine_probe(info=strategy_info, extra={'hidden_feature': _hidden}):
        P_avg = _avg_prob_dispatch(htot, htot_is_function_only_of_energy, energy, L, L0, nu_i, nu_f,
            average, 'osc_prob_matter_nsi', smooth_profile=_profile_is_smooth, engine_kwargs=scan_kwargs)
        if P_avg is not NotImplemented:
            return P_avg

        # Hybrid strategy, tried first: see _osc_prob_hybrid_dispatch and the matching comment in
        # osc_prob_matter_std_potential for why it precedes the interaction-picture fast path.
        P_hybrid = _osc_prob_hybrid_dispatch(h_vac_energy_indep, VCC_func, h_matt, None, None,
            energy, L, L0, nu_i, nu_f, scan_kwargs, strategy)
        if P_hybrid is not NotImplemented:
            return P_hybrid

        # Fast path for a genuine exponential density profile (e.g., the Sun), reached only where the
        # hybrid strategy declined: see _osc_prob_ip_exp_dispatch and the matching comment in
        # osc_prob_matter_std_potential.
        P_ip = _osc_prob_ip_exp_dispatch(h_vac_energy_indep, VCC_func, h_matt, None, None,
            energy, L, L0, nu_i, nu_f, scan_kwargs)
        if P_ip is not NotImplemented:
            return P_ip

        # Energy-batched fast path: when many energies share a single baseline and the Hamiltonian
        # is position-dependent, compute the whole scan in one batched pipeline, with the potential
        # samples shared across energies (see _osc_prob_scan_separable).  If the request does not fit
        # the engine, fall back to the generic per-point path below.
        P_scan = _osc_prob_scan_separable_dispatch(h_vac_energy_indep, VCC_func, h_matt, None, None,
            energy, L, L0, nu_i, nu_f, scan_kwargs)
        if P_scan is not NotImplemented:
            return P_scan

        # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
        return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f,
            htot_is_function_only_of_energy, t_slab_edges=t_slab_edges,
            magnus_exp_order=magnus_exp_order, n_jobs=n_jobs, integration_method=integration_method,
            rtol=rtol, atol=atol, growth_factor_n_slabs=growth_factor_n_slabs,
            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
            min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
            new_recursion_limit=new_recursion_limit, verbose=verbose,
            # An explicit cumulative= from the caller wins; otherwise strategy='magnus' opts out
            # of the cumulative scan and everything else takes 'auto'.  Resolved near the top of
            # this function, which is also where it is removed from kwargs.
            cumulative=cumulative_resolved, **kwargs)


def osc_prob_liv(
    num_flavors: int,
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    liv_params: Dict,
    rho_func: Optional[Union[Callable, int, float]]=0.0,
    L0: Optional[Union[int, float]]=0.0,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    h_liv_energy_indep: Union[list, np.ndarray]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    average: Optional[bool]=False,
    strategy: Optional[str]='auto',
    strategy_info: Optional[Dict]=None,
    t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='gl',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5,
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5,
    max_num_loops: Optional[int]=50,
    min_n_slabs: Optional[int]=1,
    max_n_slabs: Optional[int]=None,
    min_n_tpts_per_slab: Optional[int]=2,
    max_n_tpts_per_slab: Optional[int]=500,
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    new_recursion_limit: Optional[int]=5000,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for
    oscillations under (one form of) Lorentz-invariance violation, in
    vacuum or in matter.

    Middle (scenario) layer for the LIV case, generic in ``num_flavors``: unpacks ``osc_params``
    and ``liv_params``, builds the vacuum (+ matter, if ``rho_func`` is nonzero) + LIV
    Hamiltonian (via ``hamiltonians.hamiltonian_{num_flavors}nu_vacuum_energy_independent``,
    optionally ``hamiltonian_{num_flavors}nu_matter_td``, and
    ``hamiltonian_{num_flavors}nu_liv_energy_independent``), and calls
    :func:`osc_prob_energy_baseline`. Called by every ``osc_prob_{2,3,4,5}nu_vacuum_liv``,
    ``osc_prob_{2,3,4,5}nu_matter_liv_{constant,exp}_density``, and
    ``osc_prob_{2,3,4,5}nu_earth_liv``/``osc_prob_{2,3,4,5}nu_sun_liv`` wrapper.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    num_flavors : int
        Number of neutrino flavors (2, 3, 4, or 5; or higher, if ``h_vac_energy_indep``/
        ``h_liv_energy_indep`` are given).
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    osc_params : dict
        Oscillation parameters; see :func:`unpack_oscillation_params_from_dict`.
    liv_params : dict
        LIV parameters; see :func:`unpack_liv_params_from_dict`.
    rho_func : Callable, int, or float, optional
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is
        True). If 0.0 (default), the probability is for vacuum + LIV only, with no matter term.
    L0 : int or float, optional
        Initial position. Default: 0.0.
    h_vac_energy_indep : list or np.ndarray, optional
        Precomputed energy-independent vacuum Hamiltonian, used instead of ``osc_params`` when
        ``num_flavors`` exceeds ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.
    h_liv_energy_indep : list or np.ndarray, optional
        Precomputed energy-independent LIV Hamiltonian, used instead of ``liv_params`` when
        ``num_flavors`` exceeds ``globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS``.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``rho_func`` returns the density in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, ``rho_func`` directly returns the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any parameter left as
        None in ``osc_params``. Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid',
        or 'magnus'; see the ``strategy`` parameter of :func:`osc_prob_matter_std_potential` for
        the full description and :doc:`/adiabatic_strategy` for the derivation and validation of
        the ``'hybrid'``/``'auto'`` strategies (adiabatic transport with a Magnus patch at any
        non-adiabatic window, applicable to any number of flavors). Only relevant when
        ``rho_func`` is nonzero (there is no position dependence, hence no resonance, in pure
        vacuum + LIV). Default: 'auto'.

        .. versionadded:: 1.0.0
    strategy_info : dict, optional
        If given, filled in place with which engine actually answered, following the same
        out-parameter convention as ``convergence_info`` in :func:`osc_prob`.  Under
        ``strategy='auto'`` the fallbacks are silent by design -- that is right for ordinary
        calls, and wrong for anyone asking why a result moved or why a call got slow -- so this
        is how to see them without turning the fallbacks into warnings.  Keys:

        * ``'engine'`` -- ``'hybrid'``, ``'ip_exp'``, ``'separable'``, ``'cumulative'``,
          ``'magnus'`` or ``'average'``.
        * ``'family'`` -- the engine's family; see :data:`ENGINE_FAMILIES`.
        * ``'certified'`` -- for ``'hybrid'``, whether
          :func:`magnus.adiabatic.hybrid_propagator` self-certified.  ``None`` for engines
          that do not certify.  Under ``'auto'`` an uncertified hybrid result is never
          returned, so this is ``True`` whenever the engine is ``'hybrid'``; under
          ``'hybrid'`` it can be ``False``, and then it means the accuracy is **unverified**,
          not that the answer is wrong.
        * ``'declined'`` -- ``[(engine, reason)]`` for the engines that stood aside first.
          Most requests decline most engines, which is ordinary and not a finding.
        * ``'trace'`` -- every dispatch decision in order, with per-engine detail (for the
          cumulative scan, ``'n_acc'`` and whether it came from a ceiling).

        Costs nothing when omitted.  Default: None.

        .. versionadded:: 1.0.0
    t_slab_edges : list or np.ndarray, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    magnus_exp_order : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    n_jobs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    integration_method : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    rtol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    atol : int or float, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_slabs : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    growth_factor_n_tpts_per_slab : int or float
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_num_loops : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_slabs : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    min_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    max_n_tpts_per_slab : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    validate_input : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    save_log : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    filename_log : str
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    file_log : TextIOWrapper, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    close_file_log_upon_exit : bool
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    verbose : int
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    new_recursion_limit : int, optional
        Forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`; see their docstrings.
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`.

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each (energy, L) point.
    """

    if validate_input and (strategy not in ('auto', 'hybrid', 'magnus')):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob.osc_prob_liv: strategy must be" + \
            " 'auto', 'hybrid', or 'magnus'.")

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, osc_params, h_vac_energy_indep)
    liv_params_list = unpack_liv_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, liv_params, h_liv_energy_indep)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
        sxi, b1, b2, Lambda, n_liv = liv_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
        sxi12, sxi23, sxi13, dxiCP, b1, b2, b3, Lambda, n_liv = liv_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
        sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, \
            n_liv = liv_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list
        sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, \
            dxi35, b1, b2, b3, b4, b5, Lambda, n_liv = liv_params_list

    if validate_input:
        validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=L0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            rho_func=rho_func, ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction, validate_energy_and_L=True, 
            validate_flavor_indices=True, validate_osc_params=True, validate_initial_position=True,
            validate_density=True)
    
    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)
    
    # Compute the energy-independent part of the LIV Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_liv_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_liv_energy_indep = hamiltonians.hamiltonian_2nu_liv_energy_independent(sxi, b1, b2, 
            Lambda, n_liv)
    elif num_flavors == 3:
        h_liv_energy_indep = hamiltonians.hamiltonian_3nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxiCP, b1, b2, b3, Lambda, n_liv, nubar=nubar)
    elif num_flavors == 4:
        h_liv_energy_indep = hamiltonians.hamiltonian_4nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, n_liv,
            nubar=nubar)
    elif num_flavors == 5:
        h_liv_energy_indep = hamiltonians.hamiltonian_5nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, dxi35, b1,
            b2, b3, b4, b5, Lambda, n_liv, nubar=nubar)
   
    if (rho_func != 0.0): # Matter density is nonzero, include the matter term in the Hamiltonian

        # Projector onto the nu_e--nu_e entry, multiplied below by the potential VCC.  Note that
        # VCC_func already carries the antineutrino sign flip (applied inside
        # matter.vcc_func_from_rho_func), so no extra sign is applied here.
        h_matt = np.zeros((num_flavors, num_flavors))
        h_matt[0][0] = 1.0

        # Build the coherent forward potential function, VCC_func, from the density function,
        # rho_func. If the provided rho_func is the matter density (e.g., g cm^{-3}), convert
        # rho_func to a function that returns the electron number density [eV^3].
        VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
            electron_fraction, nubar, density_matter_is_in_g_per_cm3,
            density_is_of_number_of_electrons) # [eV]

        # Cache repeated evaluations of the potential on identical position grids (see
        # _PositionProfileCache)
        if isinstance(VCC_func, Callable):
            VCC_func = _PositionProfileCache(VCC_func)

        # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
        if isinstance(VCC_func, Callable):
            # VCC_func is a function of position, so the Hamiltonian is, too.  If l is an array,
            # the result is a stack of Hamiltonians with the position axis leading; this lets the
            # Magnus routines evaluate the Hamiltonian at all time points in a single call.
            def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
                vcc = np.asarray(VCC_func(l))
                return (1/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt + \
                    pow(enu,n_liv)*h_liv_energy_indep
            htot_is_function_only_of_energy = False
        else:
            # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is
            # passed to osc_prob below, osc_prob will detect that VCC_func is constant and set
            # parameters  internally for speed-up.
            h_matt = VCC_func*h_matt
            def htot(enu: Union[int, float]) -> np.ndarray:
                return (1/enu)*h_vac_energy_indep + h_matt + pow(enu,n_liv)*h_liv_energy_indep
            htot_is_function_only_of_energy = True

    else: # Matter density is zero; the only terms in the Hamiltonian are vacuum and LIV

        # Bound on this branch too, so that everything below can refer to it unconditionally.
        # Unlike the other two scenario wrappers, this one builds VCC_func only when the
        # density is nonzero -- and a LIV *vacuum* call therefore reached the profile scan with
        # the name unbound.  0.0 is the right value as well as a safe one: there is no
        # potential, so there is no position dependence for anything to hide in.
        VCC_func = 0.0

        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep + pow(enu,n_liv)*h_liv_energy_indep
        htot_is_function_only_of_energy = True


    # Built unconditionally: every entry is a parameter of this function, and the averaging
    # dispatch below needs them on the zero-density path too, where the rest of this block
    # does not apply.
    # Resolved -- and, crucially, POPPED out of kwargs -- before scan_kwargs is built, because
    # every dispatcher below declines outright on an unrecognized entry in kwargs.  Left in, an
    # explicitly-passed cumulative silently disabled the hybrid, interaction-picture and
    # separable engines: passing cumulative='auto', which is the documented default, changed
    # which engine answered (hybrid -> general ladder) and moved a 10 MeV solar single point by
    # 9.3e-06.  See _resolve_cumulative_kwarg.
    cumulative_resolved = _resolve_cumulative_kwarg(kwargs, strategy)

    scan_kwargs = dict(t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs, integration_method=integration_method, rtol=rtol, atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        save_log=save_log, file_log=file_log, cumulative=cumulative_resolved,
        kwargs=kwargs)

    # Checked here, rather than only in osc_prob: the averaged path returns before anything
    # forwards **kwargs onwards, so a check further down would see these keys on the ordinary
    # path and silently ignore them on the averaged one.
    _reject_parameter_set_metadata(kwargs, 'osc_prob_liv')

    # Phase-averaged limit, requested with average=True: exact and closed-form whenever the
    # Hamiltonian does not depend on position, so it is tried before any of the propagation
    # machinery below, all of which would resolve phases that the average discards (see
    # _avg_prob_dispatch and :mod:`magnus.avgprob`).
    # `or` on an array would ask for its truth value, so the emptiness test is explicit.
    _breakpoints = kwargs.get('t_breakpoints')
    _profile_is_smooth = ((_breakpoints is None or len(np.atleast_1d(_breakpoints)) == 0)
                          and (t_slab_edges is None))
    # One scan of the profile for the whole call, before any engine sees it: a feature
    # narrower than every grid here is the one exposure none of them can detect for itself,
    # because they all miss it together.  Depends on the profile and the interval, never on
    # energy, so it is not repeated per point.  See _scan_for_hidden_features.
    _hidden = _scan_for_hidden_features(VCC_func, L0, L, kwargs.get('t_breakpoints'))

    # Everything below is dispatch: which engine gets the request.  Watched as a unit so
    # that strategy_info can report which one answered and, for the hybrid strategy,
    # whether it certified -- see _engine_probe.  Costs one list allocation per call when
    # nobody is watching.
    with _engine_probe(info=strategy_info, extra={'hidden_feature': _hidden}):
        P_avg = _avg_prob_dispatch(htot, htot_is_function_only_of_energy, energy, L, L0, nu_i, nu_f,
            average, 'osc_prob_liv', smooth_profile=_profile_is_smooth, engine_kwargs=scan_kwargs)
        if P_avg is not NotImplemented:
            return P_avg

        # Energy-batched fast path: when many energies share a single baseline and the Hamiltonian
        # is position-dependent, compute the whole scan in one batched pipeline, with the potential
        # samples shared across energies (see _osc_prob_scan_separable).  If the request does not fit
        # the engine, fall back to the generic per-point path below.
        P_scan = NotImplemented
        if (rho_func != 0.0):  # VCC_func and h_matt exist only when there is matter
            # Hybrid strategy, tried first: see _osc_prob_hybrid_dispatch and the matching comment in
            # osc_prob_matter_std_potential for why it precedes the interaction-picture fast path.
            P_scan = _osc_prob_hybrid_dispatch(h_vac_energy_indep, VCC_func, h_matt,
                h_liv_energy_indep, n_liv, energy, L, L0, nu_i, nu_f, scan_kwargs, strategy)
            if P_scan is NotImplemented:
                # Fast path for a genuine exponential density profile (e.g., the Sun), reached only
                # where the hybrid strategy declined: see _osc_prob_ip_exp_dispatch and the matching
                # comment in osc_prob_matter_std_potential.
                P_scan = _osc_prob_ip_exp_dispatch(h_vac_energy_indep, VCC_func, h_matt,
                    h_liv_energy_indep, n_liv, energy, L, L0, nu_i, nu_f, scan_kwargs)
            if P_scan is NotImplemented:
                P_scan = _osc_prob_scan_separable_dispatch(h_vac_energy_indep, VCC_func, h_matt,
                    h_liv_energy_indep, n_liv, energy, L, L0, nu_i, nu_f, scan_kwargs)
        if P_scan is not NotImplemented:
            return P_scan

        # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
        return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f,
            htot_is_function_only_of_energy, t_slab_edges=t_slab_edges, 
            magnus_exp_order=magnus_exp_order, n_jobs=n_jobs, integration_method=integration_method,
            rtol=rtol, atol=atol, growth_factor_n_slabs=growth_factor_n_slabs,
            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
            min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
            new_recursion_limit=new_recursion_limit, verbose=verbose,
            # An explicit cumulative= from the caller wins; otherwise strategy='magnus' opts out
            # of the cumulative scan and everything else takes 'auto'.  Resolved near the top of
            # this function, which is also where it is removed from kwargs.
            cumulative=cumulative_resolved, **kwargs)


#-----------------------------------------------------------------------
# In vacuum
#-----------------------------------------------------------------------

def osc_prob_2nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    vacuum.

    By default, returns :math:`2 \times 2` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pab,
    Pba],[Pba,Pab]])``.  The matrix is symmetric, i.e., ``Pba == Pab``.

    If a single energy and baseline is given, the function returns a 
    single matrix.  If multiple energies and baselines are given, 
    function returns an NumPy array of matrices.  See examples below.

    If the probability needs to be computed multiple times, it is 
    recommended to pass the array of energies and the array of baselines
    to the function in a single call instead of calling the function
    separately for each combination of energy and baseline. The reason
    is that the function has an overhead that gets diluted when 
    computing when the input energies and baselines are many.

    Unlike :func:`osc_prob_3nu_vacuum` (and also 
    :func:`osc_prob_4nu_vacuum` and :func:`osc_prob_5nu_vacuum`), the
    oscillation parameters `sth` and `Dm2` are not optional, but must be
    passed.  Depending on the values passed, :func:`osc_prob_2nu_vacuum`
    will return probabilities for different two-neutrino systems: 

    - :math:`\nu_e-\nu_\mu` if ``sth`` is :math:`\sin \theta_{12}` and 
      ``Dm2`` is :math:`\Delta m_{21}^2`
    - :math:`\nu_\mu-\nu_\tau` if ``sth`` is :math:`\sin \theta_{23}`
      and ``Dm2`` is :math:`\Delta m_{32}^2`
    - :math:`\nu_e-\nu_\tau` if ``sth`` is :math:`\sin \theta_{13}` 
      and ``Dm2`` is :math:`\Delta m_{31}^2`.

    If the initial and final flavors, ``nu_i`` and ``nu_f``, are 
    specified (by setting them to ``NUE``, ``NUMU``, or ``NUTAU``
    from the :py:mod:`magnus.globaldefs` module), the function returns 
    instead a one-dimensional array of the probabilities computed for
    each value of energy and baseline requested. See examples below.

    Because this is a two-neutrino system, the flavor indices can only 
    be 0 or 1.  To prevent using other values, we convert the indices
    like this:
    
    - If ``nu_i == NUE`` (i.e., 0) and ``nu_f == NUTAU`` (i.e., 2), we 
      set ``nu_f = 1``
    - If ``nu_i == NUTAU`` (i.e., 2) and ``nu_f == NUE`` (i.e., 0), we 
      set ``nu_i = 1``
    - If ``nu_i == NUMU`` (i.e., 1) and ``nu_f == NUTAU`` (i.e., 2), we 
      set ``nu_i = 0`` and ``nu_f = 1``
    - If ``nu_i == NUTAU`` (i.e., 2) and ``nu_f == NUMU`` (i.e., 1), we
      set ``nu_i = 1`` and ``nu_f = 0``

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy, single value or array.
    L : int, float, list, or np.ndarray
        Neutrino baseline, single value or array.
    sth : int or float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2`.
    nu_i : int, optional
        Initial neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    nu_f : int, optional
        Final neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    validate_input : bool, optional
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose : int, optional
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

    Returns
    -------
    Union[float, np.ndarray]
        Neutrino oscillation probability matrix or probability for a 
        single oscillation channel, for the values of `energy` and `L`.

    Examples
    --------
    Single energy and baseline (the code below runs when these docs are
    built, so the output shown is always current):

    .. jupyter-execute::

        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd

        sth = gd.S12_NO_BF_NUFIT_6_0  # sin(theta) [adim]
        Dm2 = gd.D21_NO_BF_NUFIT_6_0  # [eV^2]
        baseline = 10.0 * gd.UNIT_KM  # 10 km in natural units [eV^-1]
        energy = 1.0 * gd.UNIT_MEV    # [eV]

        oscprob.osc_prob_2nu_vacuum(energy, baseline, sth, Dm2)

    .. seealso::
        :func:`osc_prob_3nu_vacuum`
            Three-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_4nu_vacuum`
            Four-flavor (3+1) oscillation probabilities in vacuum. 
        :func:`osc_prob_5nu_vacuum`
            Four-flavor (3+2) oscillation probabilities in vacuum. 
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_vacuum(
        num_flavors=2,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    vacuum.

    By default, returns :math:`3 \times 3` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pee,
    Pem,Pet],[Pme,Pmm,Pmt],[Pte,Ptm,Ptt]])``.  The matrix is symmetric, 
    i.e., ``Pme == Pee``, ``Pte == Pet``, and ``Ptm == Pmt``.  

    If a single energy and baseline is given, the function returns a 
    single matrix.  If multiple energies and baselines are given, 
    function returns an NumPy array of matrices.  See examples below.

    If the probability needs to be computed multiple times, it is 
    recommended to pass the array of energies and the array of baselines
    to the function in a single call instead of calling the function
    separately for each combination of energy and baseline. The reason
    is that the function has an overhead that gets diluted when 
    computing when the input energies and baselines are many.

    If the initial and final flavors, ``nu_i`` and ``nu_f``, are 
    specified (by setting them to ``NUE``, ``NUMU``, or ``NUTAU``
    from the :py:mod:`magnus.globaldefs` module), the function returns 
    instead a one-dimensional array of the probabilities computed for
    each value of energy and baseline requested. See examples below.

    If the function is called without specifying values of the standard
    oscillation parameters (``s12``, ``s23``, ``s13``, ``dCP``, ``D21``,
    ``D31``), the unspecified parameters are assigned default values 
    taken from a predefined parameter set.  The name of the default 
    parameter set can be changed by passing 
    ``default_osc_params_set_name``.  

    The names of the predefined parameter sets included in 
    :math:`\text{Mag}\nu\text{s}` can be seen by printing

    .. jupyter-execute::

        import magnus.globaldefs as gd

        list(gd.OSC_PARAMS_PREDEFINED.keys())

    And the default parameter values are from the set with name 
    ``'OSC_PARAMS_DEFAULT'``:

    .. jupyter-execute::

        gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy, single value or array.
    L : int, float, list, or np.ndarray
        Neutrino baseline, single value or array.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : int or float, optional
        CP-violation phase, :math:`\delta_\text{CP}`.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`.
    nubar : bool, optional
        False (default) for neutrinos; True for anti-neutrinos.
    nu_i : int, optional
        Initial neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    nu_f : int, optional
        Final neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    default_osc_params_set_name : str, optional
        Name of the predefined set of oscillation parameters to use when
        assigning default values to unspecified parameters.
    validate_input : bool, optional
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose : int, optional
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

    Returns
    -------
    Union[float, np.ndarray]
        Neutrino oscillation probability matrix or probability for a 
        single oscillation channel, for the values of `energy` and `L`.

    Examples
    --------
    If both ``energy`` and ``L`` are single values, this function returns
    the full :math:`3\times 3` probability matrix computed at those
    values, using the NuFit 6.0 (normal ordering) defaults for any
    oscillation parameter not passed explicitly:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd

        baseline = 10.0 * gd.UNIT_KM  # 10 km in natural units [eV^-1]
        energy = 1.0 * gd.UNIT_MEV    # [eV]

        oscprob.osc_prob_3nu_vacuum(energy, baseline)

    Pick one channel only, e.g., :math:`\nu_e \to \nu_\mu`, by passing an
    initial flavor, ``nu_i``, and a final flavor, ``nu_f`` (the flavor
    indices ``NUE``, ``NUMU``, ``NUTAU`` are defined in
    :py:mod:`magnus.globaldefs`); pass ``nubar=True`` for the
    antineutrino channel :math:`\bar\nu_e \to \bar\nu_\mu`:

    .. jupyter-execute::

        print(oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU))
        print(oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU,
                                           nubar=True))

    Any standard oscillation parameter can be overridden; the rest keep
    defaulting to the predefined set named by ``default_osc_params_set_name``
    (``'OSC_PARAMS_DEFAULT'`` unless changed -- see
    ``globaldefs.OSC_PARAMS_PREDEFINED`` for the values it uses):

    .. jupyter-execute::

        oscprob.osc_prob_3nu_vacuum(energy, baseline, s12=0.0)

    If a single energy and multiple baselines are passed, the result is
    an array of probabilities, one per baseline (and, conversely, one
    per energy for a single baseline and multiple energies; or a full
    grid for arrays of both -- paired index-by-index, not an outer
    product):

    .. jupyter-execute::

        baselines = gd.UNIT_KM * np.array([1.0, 10.0, 100.0])
        energies = gd.UNIT_MEV * np.array([1.0, 5.0, 20.0])

        print(oscprob.osc_prob_3nu_vacuum(energy, baselines, nu_i=gd.NUE, nu_f=gd.NUMU))
        print(oscprob.osc_prob_3nu_vacuum(energies, baseline, nu_i=gd.NUE, nu_f=gd.NUMU))
        print(oscprob.osc_prob_3nu_vacuum(energies, baselines, nu_i=gd.NUE, nu_f=gd.NUMU))

    .. seealso::
        :func:`osc_prob_2nu_vacuum`
            Two-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_4nu_vacuum`
            Four-flavor (3+1) oscillation probabilities in vacuum. 
        :func:`osc_prob_5nu_vacuum`
            Four-flavor (3+2) oscillation probabilities in vacuum. 
    """

    return osc_prob_vacuum(
        num_flavors=3,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in vacuum.

    By default, returns :math:`4 \times 4` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pee,
    Pem,Pet,Pes],[Pme,Pmm,Pmt,Pms],[Pte,Ptm,Ptt,Pts],
    [Pse,Psm,Pst,Pss]])``.  The matrix is symmetric, i.e., 
    ``Pme == Pee``, ``Pte == Pet``, ``Pse == Pes``, ``Ptm == Pmt``,
    ``Psm == Pms``, and ``Pst == Pts``.

    If a single energy and baseline is given, the function returns a 
    single matrix.  If multiple energies and baselines are given, 
    function returns an NumPy array of matrices.  See examples below.

    If the probability needs to be computed multiple times, it is 
    recommended to pass the array of energies and the array of baselines
    to the function in a single call instead of calling the function
    separately for each combination of energy and baseline. The reason
    is that the function has an overhead that gets diluted when 
    computing when the input energies and baselines are many.

    If the initial and final flavors, ``nu_i`` and ``nu_f``, are 
    specified (by setting them to ``NUE``, ``NUMU``, ``NUTAU``, or 
    ``NUS`` from the :py:mod:`magnus.globaldefs` module), the function
    returns instead a one-dimensional array of the probabilities 
    computed for each value of energy and baseline requested. See 
    examples below.

    If the function is called without specifying values of the standard
    oscillation parameters (``s12``, ``s23``, ``s13``, ``dCP``, ``D21``,
    ``D31``), the unspecified parameters are assigned default values 
    taken from a predefined parameter set.  The name of the default 
    parameter set can be changed by passing 
    ``default_osc_params_set_name``.  

    The names of the predefined parameter sets included in 
    :math:`\text{Mag}\nu\text{s}` can be seen by printing

    .. jupyter-execute::

        import magnus.globaldefs as gd

        list(gd.OSC_PARAMS_PREDEFINED.keys())

    And the default parameter values are from the set with name 
    ``'OSC_PARAMS_DEFAULT'``:

    .. jupyter-execute::

        gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy, single value or array.
    L : int, float, list, or np.ndarray
        Neutrino baseline, single value or array.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`.
    d14 : int or float, optional
        CP-violation phase, :math:`\delta_{14}`.
    d24 : int or float, optional
        CP-violation phase, :math:`\delta_{24}`.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : int or float, optional
        CP-violation phase, :math:`\delta_\text{CP}`.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`.
    nubar : bool, optional
        False (default) for neutrinos; True for anti-neutrinos.
    nu_i : int, optional
        Initial neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        or ``NUS`` from the :py:mod:`magnus.globaldefs` module.
    nu_f : int, optional
        Final neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        or ``NUS`` from the :py:mod:`magnus.globaldefs` module.
    default_osc_params_set_name : str, optional
        Name of the predefined set of standard oscillation parameters to
        use when assigning default values to unspecified parameters.
    validate_input : bool, optional
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose : int, optional
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

    Returns
    -------
    Union[float, np.ndarray]
        Neutrino oscillation probability matrix or probability for a 
        single oscillation channel, for the values of `energy` and `L`.

    Examples
    --------
    With the sterile-sector angles/phases given explicitly and the
    active-sector angles left at their NuFit 6.0 defaults:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd

        baseline = 10.0 * gd.UNIT_KM  # 10 km in natural units [eV^-1]
        energy = 1.0 * gd.UNIT_MEV    # [eV]
        s14, s24, s34 = 0.1, 0.2, 0.3
        d14, d24 = np.radians(10.0), np.radians(100.0)
        D41 = 0.1  # [eV^2]

        oscprob.osc_prob_4nu_vacuum(energy, baseline, s14=s14, s24=s24, s34=s34,
                                     d14=d14, d24=d24, D41=D41)

    Pick one channel only, e.g., :math:`\nu_e \to \nu_s`, and the
    antineutrino channel :math:`\bar\nu_e \to \bar\nu_\mu` (the flavor
    indices ``NUE``, ``NUMU``, ``NUTAU``, ``NUS`` are defined in
    :py:mod:`magnus.globaldefs`):

    .. jupyter-execute::

        common = dict(s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41)
        print(oscprob.osc_prob_4nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUS, **common))
        print(oscprob.osc_prob_4nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU,
                                           nubar=True, **common))

    Any active-sector parameter can still be overridden on top of the
    sterile-sector values:

    .. jupyter-execute::

        oscprob.osc_prob_4nu_vacuum(energy, baseline, s12=0.0, **common)

    Arrays of energies and/or baselines work the same way as for
    :func:`osc_prob_3nu_vacuum` (paired index-by-index for two arrays,
    not an outer product):

    .. jupyter-execute::

        baselines = gd.UNIT_KM * np.array([1.0, 10.0, 100.0])
        energies = gd.UNIT_MEV * np.array([1.0, 5.0, 20.0])

        print(oscprob.osc_prob_4nu_vacuum(energy, baselines, nu_i=gd.NUE, nu_f=gd.NUMU, **common))
        print(oscprob.osc_prob_4nu_vacuum(energies, baseline, nu_i=gd.NUE, nu_f=gd.NUMU, **common))
        print(oscprob.osc_prob_4nu_vacuum(energies, baselines, nu_i=gd.NUE, nu_f=gd.NUMU, **common))

    .. seealso::
        :func:`osc_prob_2nu_vacuum`
            Two-flavor oscillation probabilities in vacuum.
        :func:`osc_prob_3nu_vacuum`
            Three-flavor oscillation probabilities in vacuum.
        :func:`osc_prob_5nu_vacuum`
            Five-flavor (3+2) oscillation probabilities in vacuum.
    """

    return osc_prob_vacuum(
        num_flavors=4,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in vacuum.

    By default, returns :math:`5 \times 5` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pee,
    Pem,Pet,Pes1,Pes2],[Pme,Pmm,Pmt,Pms1,Pms2],[Pte,Ptm,Ptt,Pts1,Pts2],
    [Ps1e,Ps1m,Ps1t,Ps1s1,Ps1s2],[Ps2e,Ps2m,Ps2t,Ps2s1,Ps2s2]])``.  The 
    matrix is symmetric, i.e., ``Pme == Pee``, ``Pte == Pet``, 
    ``Ps1e == Pes1``, ``Ps2e == Pes2`` ``Ptm == Pmt``,
    ``Ps1m == Pms1``, ``Ps1t == Pts1``, ``Ps2t == Pts2``, and 
    ``Ps2s1 == Ps1s2``.

    If a single energy and baseline is given, the function returns a 
    single matrix.  If multiple energies and baselines are given, 
    function returns an NumPy array of matrices.  See examples below.

    If the probability needs to be computed multiple times, it is 
    recommended to pass the array of energies and the array of baselines
    to the function in a single call instead of calling the function
    separately for each combination of energy and baseline. The reason
    is that the function has an overhead that gets diluted when 
    computing when the input energies and baselines are many.

    If the initial and final flavors, ``nu_i`` and ``nu_f``, are 
    specified (by setting them to ``NUE``, ``NUMU``, ``NUTAU``, 
    ``NUS1``, or ``NUS2`` from the :py:mod:`magnus.globaldefs` module),
    the function returns instead a one-dimensional array of the 
    probabilities computed for each value of energy and baseline 
    requested. See examples below.

    If the function is called without specifying values of the standard
    oscillation parameters (``s12``, ``s23``, ``s13``, ``dCP``, ``D21``,
    ``D31``), the unspecified parameters are assigned default values 
    taken from a predefined parameter set.  The name of the default 
    parameter set can be changed by passing 
    ``default_osc_params_set_name``.  

    The names of the predefined parameter sets included in 
    :math:`\text{Mag}\nu\text{s}` can be seen by printing

    .. jupyter-execute::

        import magnus.globaldefs as gd

        list(gd.OSC_PARAMS_PREDEFINED.keys())

    And the default parameter values are from the set with name 
    ``'OSC_PARAMS_DEFAULT'``:

    .. jupyter-execute::

        gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy, single value or array.
    L : int, float, list, or np.ndarray
        Neutrino baseline, single value or array.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`.
    d14 : int or float, optional
        CP-violation phase, :math:`\delta_{14}`.
    d15 : int or float, optional
        CP-violation phase, :math:`\delta_{15}`.
    d24 : int or float, optional
        CP-violation phase, :math:`\delta_{24}`.
    d35 : int or float, optional
        CP-violation phase, :math:`\delta_{35}`.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : int or float, optional
        CP-violation phase, :math:`\delta_\text{CP}`.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`.
    nubar : bool, optional
        False (default) for neutrinos; True for anti-neutrinos.
    nu_i : int, optional
        Initial neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        ``NUS1``, or ``NUS2`` from the :py:mod:`magnus.globaldefs`
        module.
    nu_f : int, optional
        Final neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        ``NUS1``, or ``NUS2`` from the :py:mod:`magnus.globaldefs`
        module.
    default_osc_params_set_name : str, optional
        Name of the predefined set of standard oscillation parameters to
        use when assigning default values to unspecified parameters.
    validate_input : bool, optional
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose : int, optional
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

    Returns
    -------
    Union[float, np.ndarray]
        Neutrino oscillation probability matrix or probability for a 
        single oscillation channel, for the values of `energy` and `L`.

    Examples
    --------
    With the sterile-sector angles/phases given explicitly and the
    active-sector angles left at their NuFit 6.0 defaults:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd

        baseline = 10.0 * gd.UNIT_KM  # 10 km in natural units [eV^-1]
        energy = 1.0 * gd.UNIT_MEV    # [eV]
        s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 1.e-2, 1.e-2, 1.e-3, 1.e-3
        d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
        D41, D51 = 0.1, 0.001  # [eV^2]
        common = dict(s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35,
                      d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51)

        oscprob.osc_prob_5nu_vacuum(energy, baseline, **common).shape

    Pick one channel only -- :math:`\nu_e \to \nu_{s_1}`,
    :math:`\nu_e \to \nu_{s_2}`, and :math:`\nu_{s_1} \to \nu_{s_2}` (the
    flavor indices ``NUE``, ``NUMU``, ``NUTAU``, ``NUS1``, ``NUS2`` are
    defined in :py:mod:`magnus.globaldefs`); and the antineutrino channel
    :math:`\bar\nu_e \to \bar\nu_{s_1}`:

    .. jupyter-execute::

        print(oscprob.osc_prob_5nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUS1, **common))
        print(oscprob.osc_prob_5nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUS2, **common))
        print(oscprob.osc_prob_5nu_vacuum(energy, baseline, nu_i=gd.NUS1, nu_f=gd.NUS2, **common))
        print(oscprob.osc_prob_5nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUS1,
                                           nubar=True, **common))

    Any active-sector parameter can still be overridden on top of the
    sterile-sector values:

    .. jupyter-execute::

        oscprob.osc_prob_5nu_vacuum(energy, baseline, s12=0.0, **common)

    Arrays of energies and/or baselines work the same way as for
    :func:`osc_prob_3nu_vacuum` (paired index-by-index for two arrays,
    not an outer product):

    .. jupyter-execute::

        baselines = gd.UNIT_KM * np.array([1.0, 10.0, 100.0])
        energies = gd.UNIT_MEV * np.array([1.0, 5.0, 20.0])

        print(oscprob.osc_prob_5nu_vacuum(energy, baselines, nu_i=gd.NUE, nu_f=gd.NUMU, **common))
        print(oscprob.osc_prob_5nu_vacuum(energies, baseline, nu_i=gd.NUE, nu_f=gd.NUMU, **common))
        print(oscprob.osc_prob_5nu_vacuum(energies, baselines, nu_i=gd.NUE, nu_f=gd.NUMU, **common))

    .. seealso::
        :func:`osc_prob_2nu_vacuum`
            Two-flavor oscillation probabilities in vacuum.
        :func:`osc_prob_3nu_vacuum`
            Three-flavor oscillation probabilities in vacuum.
        :func:`osc_prob_4nu_vacuum`
            Four-flavor (3+1) oscillation probabilities in vacuum.
    """

    return osc_prob_vacuum(
        num_flavors=5,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


#-----------------------------------------------------------------------
# In matter, standard oscillations, constant density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    return osc_prob_matter_std_potential(
        num_flavors=3,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with a constant density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    return osc_prob_matter_std_potential(
        num_flavors=4,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with a constant density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, standard oscillations, exponentially falling density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    sth: Union[int, float],
    Dm2: Union[int, float],
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation 
    probability in matter with an exponentially falling density profile.

    .. versionadded:: 1.0.0

    .. note::
        Dispatches to a fast, closed-form interaction-picture Magnus
        integrator whenever the accumulated matter phase stays small enough
        to certify (see ``_osc_prob_ip_exp_dispatch``), giving warning-free
        results in a fraction of a second across the realistic solar-neutrino
        energy range for baselines up to a few e-folds of ``l_scale``. Longer
        baselines fall back transparently to the general slab-refinement
        method.

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_2nu_matter_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation 
    probability in matter with an exponentially falling density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_3nu_matter_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    return osc_prob_matter_std_potential(
        num_flavors=3,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with an exponentially falling density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_4nu_matter_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    return osc_prob_matter_std_potential(
        num_flavors=4,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with an exponentially falling density profile.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_5nu_matter_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, standard oscillations, in the Earth
#-----------------------------------------------------------------------

def osc_prob_2nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Standard two-neutrino oscillations through the Earth, specified by
    the cosine of the zenith angle:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        sth = gd.S12_NO_BF_NUFIT_6_0
        Dm2 = gd.D21_NO_BF_NUFIT_6_0
        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM  # chord length for this costhz
        energy = 1.0 * gd.UNIT_GEV

        # The small solar mass splitting Dm2 combined with this Earth baseline
        # means the adaptive refinement needs a few loops; this is the expected,
        # informational MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_2nu_earth(energy, sth, Dm2, costhz=costhz, L=baseline)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose) # L in eV^{-1}

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # [g cm^{-3}] as a function of radial distance, r, using the Preliminary Reference Earth Model 
    # (PREM). The function matter.num_density_e_func converts the matter density into electron 
    # number density [eV^3].

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L, # [eV^{-1}]
        t_breakpoints=t_breakpoints,
        osc_params={'sth': sth, 'Dm2': Dm2},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Standard three-neutrino oscillations through the Earth, using the
    NuFit 6.0 defaults for the oscillation parameters:

    .. jupyter-execute::

        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_3nu_earth(energy, costhz=costhz, L=baseline)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=3,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Four-neutrino (3+1 sterile) oscillations through the Earth, with a
    modest sterile mixing on top of the NuFit 6.0 active-sector defaults:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_4nu_earth(energy, costhz=costhz, L=baseline,
                                        s14=0.1, s24=0.05, s34=0.02,
                                        d14=np.radians(10.0), d24=np.radians(20.0), D41=0.1)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=4,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Five-neutrino (3+2 sterile) oscillations through the Earth, with
    modest sterile mixing on top of the NuFit 6.0 active-sector defaults:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_5nu_earth(energy, costhz=costhz, L=baseline,
                                        s14=0.1, s15=0.05, s24=0.05, s25=0.02, s34=0.02, s35=0.01,
                                        d14=np.radians(10.0), d15=np.radians(15.0),
                                        d24=np.radians(20.0), d35=np.radians(25.0),
                                        D41=0.1, D51=0.05)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_earth(
    H_func: Callable,
    energy: Union[int, float, list, np.ndarray],
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None,
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None,
    L: Optional[Union[float, list, np.ndarray]]=None,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='gl',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    validate_input: Optional[bool]=True,
    verbose: Optional[int]=0,
    strategy: Optional[str]='auto',
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the neutrino oscillation probability inside
    the Earth for a given arbitrary Hamiltonian.

    Does **not** assume standard oscillations nor a given number of
    neutrino flavors: the user supplies their own Hamiltonian function,
    ``H_func``, and this routine takes care of the geometry of the
    trajectory through the Earth and of the matter density along it.

    ``H_func`` must be a function of either three arguments,
    ``H_func(energy, l, VCC)``, or two arguments,
    ``H_func(energy, l)``, returning a square complex NumPy array (the
    Hamiltonian in the flavor basis, in eV).  In the three-argument
    form, ``VCC`` is the charged-current matter potential
    :math:`V_{\rm CC} = \sqrt{2} G_F N_e` [eV] at position ``l``
    along the chord, computed from the Preliminary Reference Earth
    Model; its sign is already flipped for antineutrinos
    (``nubar=True``).  The user is free to use it, scale it, or ignore
    it (e.g., to add non-standard matter potentials that affect flavors
    other than :math:`\nu_e`).  For extra speed, ``H_func`` may accept
    an array of positions ``l`` and return a stack of Hamiltonians with
    the position axis leading; this is detected automatically.

    The trajectory can be specified either by the cosine of the zenith
    angle (``costhz``) together with the baseline ``L`` [:math:`\text{eV}^{-1}`], or
    by an initial and a final location on the surface of the Earth
    (``loc_ini``, ``loc_fin``), given as (degree, minute, second)
    latitude/longitude tuples or as the names of predefined locations
    (see ``earth.loc_coords_dms``); in the latter case the neutrino
    travels the chord that joins the two locations.

    The slab edges used internally are aligned with the crossings of
    the PREM layer boundaries along the chord.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        The Hamiltonian, as ``H_func(energy, l, VCC)`` or ``H_func(energy, l)``; see above.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative
        to ``loc_ini``/``loc_fin``.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a
        predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``.
    L : float, list, or np.ndarray, optional
        Baseline(s) [:math:`\text{eV}^{-1}`]. Used together with ``costhz``, as an alternative to
        ``loc_ini``/``loc_fin``.
    nubar : bool, optional
        If True, compute the probability for antineutrinos (flips the sign of the PREM-based
        matter potential passed to ``H_func``). Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in Earth matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction of Earth matter. Default: 0.5.
    magnus_exp_order : int, optional
        Highest order of the Magnus expansion. Default: 4.
    n_jobs : int, optional
        Number of parallel joblib workers. Default: 1.
    integration_method : str, optional
        'gl', 'trapezoid', or 'simpson'. Default: 'gl'.
    rtol, atol : int or float, optional
        Target relative/absolute tolerance for the adaptive slab refinement. Default: 1e-3 each.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    verbose : int, optional
        Verbosity level. Default: 0.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid',
        or 'magnus'; see the ``strategy`` parameter of :func:`osc_prob_matter_std_potential` for
        the full description and :doc:`/adiabatic_strategy` for the derivation and validation. In
        practice, ``'hybrid'``/``'auto'`` rarely engage here: the PREM density profile has
        layer-boundary discontinuities (``t_breakpoints``), which this strategy does not support
        (see :doc:`/adiabatic_strategy`), so a real Earth-crossing trajectory almost always falls
        back to the ``'magnus'`` strategies regardless of what is requested. Default: 'auto'.

        .. versionadded:: 1.0.0
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`
        (e.g., the refinement-loop bounds).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each energy.

    Examples
    --------
    Standard three-neutrino oscillations, written by hand (the
    dedicated wrapper :func:`osc_prob_3nu_earth` does this internally):

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.hamiltonians as hamiltonians
        import magnus.globaldefs as gd

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        s12, s23, s13, dCP, D21, D31 = p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31']

        h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
            s12, s23, s13, dCP, D21, D31)

        # Written so that it accepts an array of positions: VCC[..., None, None]
        # turns one potential per position into a stack of matrices, which keeps
        # osc_prob on its vectorized path (see ScalarHamiltonianWarning).
        e00 = np.diag([1.0, 0.0, 0.0])
        def H(energy, l, VCC):
            return (1 / energy) * h_vac + np.asarray(VCC)[..., None, None] * e00

        oscprob.osc_prob_earth(H, energy=1.0 * gd.UNIT_GEV, loc_ini='fermilab',
                                loc_fin='homestake')
    """
    source_func_name = sys._getframe().f_code.co_name

    # If the location is given as a string, look it up among the predefined named locations
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # Resolve the trajectory: either the chord between two surface locations, or (costhz, L)
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]

    # Charged-current potential along the chord from the PREM electron density; the antineutrino
    # sign flip is applied inside matter.vcc_func_from_rho_func.  The profile evaluations are
    # cached on repeated position grids.
    VCC_func = matter.vcc_func_from_rho_func(
        rho_func=lambda l: matter.num_density_e_func(
            earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM),
            earth.density_matter_func_prem,
            ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction,
            density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        nubar=nubar,
        density_is_of_number_of_electrons=True) # [eV]
    VCC_func = _PositionProfileCache(VCC_func)

    return _osc_prob_with_potential(source_func_name, H_func, VCC_func, energy, L, 0.0, nu_i,
        nu_f, t_breakpoints, magnus_exp_order, n_jobs, integration_method, rtol, atol,
        validate_input, verbose, strategy=strategy, **kwargs)


def _osc_prob_with_potential(
    source_func_name: str,
    H_func: Callable,
    VCC_func: Callable,
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    t_breakpoints: Optional[np.ndarray],
    magnus_exp_order: int,
    n_jobs: int,
    integration_method: str,
    rtol: Optional[Union[int, float]],
    atol: Optional[Union[int, float]],
    validate_input: bool,
    verbose: int,
    strategy: Optional[str] = 'auto',
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Common machinery of :func:`osc_prob_earth` and
    :func:`osc_prob_sun`: wire a user-supplied Hamiltonian function --
    H_func(energy, l, VCC) or H_func(energy, l) -- to the environment
    potential ``VCC_func`` and hand it to
    :func:`osc_prob_energy_baseline`.

    .. versionadded:: 1.0.0

    .. note::
        With ``strategy='auto'`` (the default) or ``'hybrid'``, this also tries the
        adiabatic-transport-plus-Magnus-patch hybrid strategy (see
        ``_osc_prob_hybrid_dispatch_generic`` and :doc:`/adiabatic_strategy`) whenever
        ``t_breakpoints`` is empty and a target tolerance is requested, before falling back to
        the general slab-refinement method.

    Parameters
    ----------
    source_func_name : str
        Name of the calling function (``osc_prob_earth`` or ``osc_prob_sun``), used to build more
        informative error messages.
    H_func : Callable
        The Hamiltonian, as ``H_func(energy, l, VCC)`` or ``H_func(energy, l)``.
    VCC_func : Callable
        The environment's matter potential, as a function of position.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    t_breakpoints : np.ndarray, optional
        Mandatory slab edges (e.g., PREM layer boundaries).
    magnus_exp_order : int
        Highest order of the Magnus expansion.
    n_jobs : int
        Number of parallel joblib workers.
    integration_method : str
        'gl', 'trapezoid', or 'simpson'.
    rtol, atol : int or float, optional
        Target relative/absolute tolerance for the adaptive slab refinement.
    validate_input : bool
        If True, validate that ``H_func`` has the expected signature.
    verbose : int
        Verbosity level.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid',
        or 'magnus'; see the ``strategy`` parameter of :func:`osc_prob_matter_std_potential` for
        the full description and :doc:`/adiabatic_strategy` for the derivation and validation.
        Default: 'auto'.

        .. versionadded:: 1.0.0
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`.

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each (energy, L) point.
    """

    if validate_input:
        if not isinstance(H_func, Callable):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": H_func must be a function of (energy, l, VCC) or of (energy, l).")
        n_params_H = _n_required_params(H_func)
        if n_params_H not in (2, 3):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": H_func must be a function of either three arguments (energy, l, VCC) or" + \
                " two arguments (energy, l); the provided H_func takes " + \
                str(n_params_H) + " argument(s).")
        if strategy not in ('auto', 'hybrid', 'magnus'):
            raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + \
                ": strategy must be 'auto', 'hybrid', or 'magnus'.")

    n_params_H = _n_required_params(H_func)
    if n_params_H == 3:
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            return H_func(enu, l, VCC_func(l))
    else:
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            return H_func(enu, l)

    # Hybrid strategy (adiabatic transport + Magnus patch at any non-adiabatic window; see
    # _osc_prob_hybrid_dispatch_generic and :doc:`/adiabatic_strategy`). Falls back transparently
    # (returns NotImplemented) if it does not apply -- in particular, this is essentially always
    # the case for osc_prob_earth, since t_breakpoints (the PREM layer crossings) is virtually
    # never empty for a real trajectory -- or (with strategy == 'auto' only) fails to certify.
    #
    # cumulative is popped out of kwargs first, for the reason given in the three scenario
    # wrappers: every dispatcher declines on an unrecognized entry in kwargs, so leaving it
    # there made passing cumulative -- even the default 'auto' -- silently disable this
    # strategy instead of configuring the scan.  An explicit True still stands aside here,
    # since it names one engine and is documented to raise rather than be substituted for.
    cumulative = kwargs.pop('cumulative', 'auto')
    # The same one-per-call profile scan the scenario wrappers run; the blind spot is a property
    # of the grids, not of which entry point built the Hamiltonian.
    _scan_for_hidden_features(VCC_func, L0, L, t_breakpoints)
    P_hybrid = (NotImplemented if cumulative is True else
        _osc_prob_hybrid_dispatch_generic(htot, VCC_func, energy, L, L0, nu_i, nu_f,
            t_breakpoints, rtol, atol, magnus_exp_order, integration_method, strategy, kwargs))
    if P_hybrid is not NotImplemented:
        return P_hybrid

    return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f, False,
        t_breakpoints=t_breakpoints, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
        integration_method=integration_method, rtol=rtol, atol=atol,
        validate_input=validate_input, verbose=verbose, cumulative=cumulative, **kwargs)


#-----------------------------------------------------------------------
# In matter, standard oscillations, in the Sun
#-----------------------------------------------------------------------

def osc_prob_2nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    sth: Union[int, float],
    Dm2: Union[int, float],
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability
    for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\nu_e` on 
    electrons.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Standard two-neutrino oscillations from the center of the Sun to
    90% of its radius:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        sth = gd.S12_NO_BF_NUFIT_6_0
        Dm2 = gd.D21_NO_BF_NUFIT_6_0
        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_2nu_sun(energy, L, L0, sth, Dm2)
        P

    .. versionadded:: 1.0.0

    .. note::
        Dispatches to a fast, closed-form interaction-picture Magnus
        integrator whenever the accumulated matter phase stays small enough
        to certify (see ``_osc_prob_ip_exp_dispatch``), giving warning-free
        results in a fraction of a second across the realistic solar-neutrino
        energy range for baselines up to a few e-folds of ``l_scale``. Longer
        baselines fall back transparently to the general slab-refinement
        method.

    .. note::
        With the default ``strategy='auto'``, this also tries the more general
        adiabatic-transport-plus-Magnus-patch hybrid strategy (see
        :func:`magnus.adiabatic.hybrid_propagator` and :doc:`/adiabatic_strategy`) for baselines
        beyond the interaction-picture integrator's reach (e.g., low-energy neutrinos over most
        of the Sun's radius), before falling back to the general slab-refinement method.

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid', or
        'magnus'; see the ``strategy`` parameter of :func:`osc_prob_matter_std_potential` for the
        full description and :doc:`/adiabatic_strategy` for the derivation and validation.
        Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_2nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        sth=sth,
        Dm2=Dm2,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_3nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\nu_e` on 
    electrons.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Standard three-neutrino oscillations through the Sun, using the
    NuFit 6.0 defaults for the oscillation parameters:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_3nu_sun(energy, L, L0)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_3nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\nu_e` on 
    electrons.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Four-neutrino (3+1 sterile) oscillations through the Sun, with a
    modest sterile mixing on top of the NuFit 6.0 active-sector defaults:

    .. jupyter-execute::

        import warnings
        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_4nu_sun(energy, L, L0, s14=0.1, s24=0.05, s34=0.02,
                                      d14=np.radians(10.0), d24=np.radians(20.0), D41=0.1)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_4nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s24=s24,
        s34=s34,
        d14=d14,
        d24=d24,
        D41=D41,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\nu_e` on 
    electrons.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Five-neutrino (3+2 sterile) oscillations through the Sun, with
    modest sterile mixing on top of the NuFit 6.0 active-sector defaults:

    .. jupyter-execute::

        import warnings
        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_5nu_sun(energy, L, L0,
                                      s14=0.1, s15=0.05, s24=0.05, s25=0.02, s34=0.02, s35=0.01,
                                      d14=np.radians(10.0), d15=np.radians(15.0),
                                      d24=np.radians(20.0), d35=np.radians(25.0),
                                      D41=0.1, D51=0.05)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_5nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s15=s15,
        s24=s24,
        s25=s25,
        s34=s34,
        s35=s35,
        d14=d14,
        d15=d15,
        d24=d24,
        d35=d35,
        D41=D41,
        D51=D51,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_sun(
    H_func: Callable,
    energy: Union[int, float, list, np.ndarray],
    L: Union[float, list, np.ndarray],
    L0: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='gl',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    validate_input: Optional[bool]=True,
    verbose: Optional[int]=0,
    strategy: Optional[str]='auto',
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the neutrino oscillation probability inside
    the Sun for a given arbitrary Hamiltonian.

    Does **not** assume standard oscillations nor a given number of
    neutrino flavors: the user supplies their own Hamiltonian function,
    ``H_func``, and this routine provides the solar electron density
    along the (radial) trajectory.

    ``H_func`` must be a function of either three arguments,
    ``H_func(energy, l, VCC)``, or two arguments,
    ``H_func(energy, l)``, returning a square complex NumPy array (the
    Hamiltonian in the flavor basis, in eV).  In the three-argument
    form, ``VCC`` is the charged-current matter potential
    :math:`V_{\rm CC} = \sqrt{2} G_F N_e` [eV] at radial position
    ``l``; its sign is already flipped for antineutrinos
    (``nubar=True``).  For extra speed, ``H_func`` may accept an array
    of positions ``l`` and return a stack of Hamiltonians with the
    position axis leading; this is detected automatically.

    The neutrino travels radially outward from ``L0`` to ``L`` (both in
    :math:`\text{eV}^{-1}`, measured from the center of the Sun).

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`,
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung
    Wook Kim.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        The Hamiltonian, as ``H_func(energy, l, VCC)`` or ``H_func(energy, l)``; see above.
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Final radial position(s) [:math:`\text{eV}^{-1}`], measured from the center of the Sun.
    L0 : int or float, optional
        Initial radial position [:math:`\text{eV}^{-1}`]. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos (flips the sign of the solar matter
        potential passed to ``H_func``). Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned
        instead of the full probability matrix.
    nu_f : int, optional
        Final flavor index; see ``nu_i``.
    magnus_exp_order : int, optional
        Highest order of the Magnus expansion. Default: 4.
    n_jobs : int, optional
        Number of parallel joblib workers. Default: 1.
    integration_method : str, optional
        'gl', 'trapezoid', or 'simpson'. Default: 'gl'.
    rtol, atol : int or float, optional
        Target relative/absolute tolerance for the adaptive slab refinement. Default: 1e-3 each.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    verbose : int, optional
        Verbosity level. Default: 0.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default), 'hybrid',
        or 'magnus'; see the ``strategy`` parameter of :func:`osc_prob_matter_std_potential` for
        the full description and :doc:`/adiabatic_strategy` for the derivation and validation
        (adiabatic transport with a Magnus patch at any non-adiabatic window, applicable to any
        ``H_func`` regardless of its internal structure). Default: 'auto'.

        .. versionadded:: 1.0.0
    \**kwargs
        Additional arguments forwarded to :func:`osc_prob_energy_baseline`/:func:`osc_prob`
        (e.g., the refinement-loop bounds).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for
        each energy.

    Examples
    --------
    Standard two-neutrino oscillations, written by hand (the dedicated
    wrapper :func:`osc_prob_2nu_sun` does this internally):

    .. jupyter-execute::

        import warnings
        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.hamiltonians as hamiltonians
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        sth = gd.S12_NO_BF_NUFIT_6_0
        Dm2 = gd.D21_NO_BF_NUFIT_6_0
        h_vac = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)

        # Array-capable: VCC[..., None, None] broadcasts one potential per
        # position over a stack of matrices, keeping the vectorized path.
        e00 = np.diag([1.0, 0.0])
        def H(energy, l, VCC):
            return (1 / energy) * h_vac + np.asarray(VCC)[..., None, None] * e00

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_sun(H, energy=1.0 * gd.UNIT_GEV,
                                      L=0.9 * gd.SUN_RADIUS * gd.UNIT_KM)
        P
    """
    source_func_name = sys._getframe().f_code.co_name

    # Solar electron number density [eV^3] along the radial trajectory; the antineutrino sign
    # flip of the potential is applied inside matter.vcc_func_from_rho_func.  The profile
    # evaluations are cached on repeated position grids.
    VCC_func = matter.vcc_func_from_rho_func(
        rho_func=lambda l: matter.density_matter_func_exp(l, gd.NUM_DENSITY_E_SUN_CENTRAL,
            gd.L_SCALE_SUN), # [eV^3] (l in eV^{-1})
        L0=L0,
        nubar=nubar,
        density_is_of_number_of_electrons=True) # [eV]
    VCC_func = _PositionProfileCache(VCC_func)

    return _osc_prob_with_potential(source_func_name, H_func, VCC_func, energy, L, L0, nu_i,
        nu_f, None, magnus_exp_order, n_jobs, integration_method, rtol, atol,
        validate_input, verbose, strategy=strategy, **kwargs)


#-----------------------------------------------------------------------
# In matter, NSI, constant density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    eps_aa : int or float, optional
        Non-universal diagonal NSI coupling of nu_e (relative to nu_mu, whose diagonal coupling is fixed to 0 by convention); see ``hamiltonians.hamiltonian_2nu_nsi``. Default: 0.0.
    eps_ab : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_nsi(
        num_flavors=2,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nsi_params={'eps_aa': eps_aa, 'eps_ab': eps_ab},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    return osc_prob_matter_nsi(
        num_flavors=3,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_mm': eps_mm,
            'eps_mt': eps_mt, 'eps_tt': eps_tt},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s) NSI coupling. Default: 0.0.
    eps_ss : int or float, optional
        Diagonal NSI coupling of nu_s. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    return osc_prob_matter_nsi(
        num_flavors=4,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14,
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es': eps_es,
            'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms': eps_ms, 'eps_tt': eps_tt,
            'eps_ts': eps_ts, 'eps_ss': eps_ss},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es1 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s1) NSI coupling. Default: 0.0.
    eps_es2 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s2) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms1 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s1) NSI coupling. Default: 0.0.
    eps_ms2 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s2) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts1 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s1) NSI coupling. Default: 0.0.
    eps_ts2 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s2) NSI coupling. Default: 0.0.
    eps_s1s1 : int or float, optional
        Diagonal NSI coupling of nu_s1. Default: 0.0.
    eps_s1s2 : int or float, optional
        Flavor-off-diagonal (nu_s1-nu_s2) NSI coupling. Default: 0.0.
    eps_s2s2 : int or float, optional
        Diagonal NSI coupling of nu_s2. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    return osc_prob_matter_nsi(
        num_flavors=5,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14,
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35,
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es1': eps_es1,
            'eps_es2': eps_es2, 'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms1': eps_ms1,
            'eps_ms2': eps_ms2, 'eps_tt': eps_tt, 'eps_ts1': eps_ts1, 'eps_ts2': eps_ts2,
            'eps_s1s1': eps_s1s1, 'eps_s1s2': eps_s1s2, 'eps_s2s2': eps_s2s2},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, NSI, exponentially falling density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_nsi_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    sth: Union[int, float],
    Dm2: Union[int, float],
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with an exponentially falling density profile, including
    non-standard interactions (NSI).

    .. versionadded:: 1.0.0

    .. note::
        Dispatches to a fast, closed-form interaction-picture Magnus
        integrator whenever the accumulated matter phase stays small enough
        to certify (see ``_osc_prob_ip_exp_dispatch``), giving warning-free
        results in a fraction of a second across the realistic solar-neutrino
        energy range for baselines up to a few e-folds of ``l_scale`` (the NSI
        couplings are folded into the same fast path). Longer baselines fall
        back transparently to the general slab-refinement method.

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    eps_aa : int or float, optional
        Non-universal diagonal NSI coupling of nu_e (relative to nu_mu, whose diagonal coupling is fixed to 0 by convention); see ``hamiltonians.hamiltonian_2nu_nsi``. Default: 0.0.
    eps_ab : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_2nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
            "non-negative.")

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_nsi(
        num_flavors=2,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nsi_params={'eps_aa': eps_aa, 'eps_ab': eps_ab},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_nsi_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with an exponentially falling density profile, including
    non-standard interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_3nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
            "non-negative.")

    return osc_prob_matter_nsi(
        num_flavors=3,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_mm': eps_mm,
            'eps_mt': eps_mt, 'eps_tt': eps_tt},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_nsi_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0,
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with an exponentially falling density profile,
    including non-standard interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s) NSI coupling. Default: 0.0.
    eps_ss : int or float, optional
        Diagonal NSI coupling of nu_s. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_4nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
            "non-negative.")

    return osc_prob_matter_nsi(
        num_flavors=4,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es': eps_es, 
            'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms': eps_ms, 'eps_tt': eps_tt,
            'eps_ts': eps_ts, 'eps_ss': eps_ss},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_nsi_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with an exponentially falling density profile,
    including non-standard interactions (NSI).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es1 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s1) NSI coupling. Default: 0.0.
    eps_es2 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s2) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms1 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s1) NSI coupling. Default: 0.0.
    eps_ms2 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s2) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts1 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s1) NSI coupling. Default: 0.0.
    eps_ts2 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s2) NSI coupling. Default: 0.0.
    eps_s1s1 : int or float, optional
        Diagonal NSI coupling of nu_s1. Default: 0.0.
    eps_s1s2 : int or float, optional
        Flavor-off-diagonal (nu_s1-nu_s2) NSI coupling. Default: 0.0.
    eps_s2s2 : int or float, optional
        Diagonal NSI coupling of nu_s2. Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_5nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
            "non-negative.")

    return osc_prob_matter_nsi(
        num_flavors=5,
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14,
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35,
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es1': eps_es1,
            'eps_es2': eps_es2, 'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms1': eps_ms1,
            'eps_ms2': eps_ms2, 'eps_tt': eps_tt, 'eps_ts1': eps_ts1, 'eps_ts2': eps_ts2,
            'eps_s1s1': eps_s1s1, 'eps_s1s2': eps_s1s2, 'eps_s2s2': eps_s2s2},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


#-----------------------------------------------------------------------
# In matter, NSI, in the Earth
#-----------------------------------------------------------------------

def osc_prob_2nu_earth_nsi(
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Two-neutrino oscillations through the Earth with non-standard
    interactions:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        sth = gd.S12_NO_BF_NUFIT_6_0
        Dm2 = gd.D21_NO_BF_NUFIT_6_0
        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The small solar mass splitting Dm2 combined with this Earth baseline
        # means the adaptive refinement needs a few loops; this is the expected,
        # informational MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_2nu_earth_nsi(energy, sth, Dm2, eps_aa=0.05, eps_ab=0.02,
                                                costhz=costhz, L=baseline)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    eps_aa : int or float, optional
        Non-universal diagonal NSI coupling of nu_e (relative to nu_mu, whose diagonal coupling is fixed to 0 by convention); see ``hamiltonians.hamiltonian_2nu_nsi``. Default: 0.0.
    eps_ab : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nsi_params={'eps_aa': eps_aa, 'eps_ab': eps_ab},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_earth_nsi(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Three-neutrino oscillations through the Earth with non-standard
    interactions, using the NuFit 6.0 defaults for the standard
    oscillation parameters:

    .. jupyter-execute::

        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_3nu_earth_nsi(energy, costhz=costhz, L=baseline,
                                            eps_ee=0.05, eps_em=-0.03, eps_et=0.01)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=3,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_mm': eps_mm,
            'eps_mt': eps_mt, 'eps_tt': eps_tt},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_earth_nsi(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0,
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Four-neutrino (3+1 sterile) oscillations through the Earth with
    non-standard interactions, including the sterile-flavor couplings:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_4nu_earth_nsi(
                energy, costhz=costhz, L=baseline,
                s14=0.1, s24=0.05, s34=0.02, d14=np.radians(10.0), d24=np.radians(20.0), D41=0.1,
                eps_ee=0.05, eps_em=-0.03, eps_et=0.01, eps_es=0.02, eps_ms=0.01, eps_ts=0.01)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s) NSI coupling. Default: 0.0.
    eps_ss : int or float, optional
        Diagonal NSI coupling of nu_s. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=4,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es': eps_es, 
            'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms': eps_ms, 'eps_tt': eps_tt,
            'eps_ts': eps_ts, 'eps_ss': eps_ss},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_earth_nsi(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Five-neutrino (3+2 sterile) oscillations through the Earth with
    non-standard interactions, including the sterile-flavor couplings:

    .. jupyter-execute::

        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_5nu_earth_nsi(
                energy, costhz=costhz, L=baseline,
                s14=0.1, s15=0.05, s24=0.05, s25=0.02, s34=0.02, s35=0.01,
                d14=np.radians(10.0), d15=np.radians(15.0), d24=np.radians(20.0), d35=np.radians(25.0),
                D41=0.1, D51=0.05,
                eps_ee=0.05, eps_em=-0.03, eps_et=0.01, eps_es1=0.02, eps_es2=0.01,
                eps_ms1=0.01, eps_ms2=0.01, eps_ts1=0.01, eps_ts2=0.01)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es1 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s1) NSI coupling. Default: 0.0.
    eps_es2 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s2) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms1 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s1) NSI coupling. Default: 0.0.
    eps_ms2 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s2) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts1 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s1) NSI coupling. Default: 0.0.
    eps_ts2 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s2) NSI coupling. Default: 0.0.
    eps_s1s1 : int or float, optional
        Diagonal NSI coupling of nu_s1. Default: 0.0.
    eps_s1s2 : int or float, optional
        Flavor-off-diagonal (nu_s1-nu_s2) NSI coupling. Default: 0.0.
    eps_s2s2 : int or float, optional
        Diagonal NSI coupling of nu_s2. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=5,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es1': eps_es1,
            'eps_es2': eps_es2, 'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms1': eps_ms1, 
            'eps_ms2': eps_ms2, 'eps_tt': eps_tt, 'eps_ts1': eps_ts1, 'eps_ts2': eps_ts2, 
            'eps_s1s1': eps_s1s1, 'eps_s1s2': eps_s1s2, 'eps_s2s2': eps_s2s2},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, NSI, in the Sun
#-----------------------------------------------------------------------

def osc_prob_2nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    sth: Union[int, float],
    Dm2: Union[int, float],
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    for neutrinos inside the Sun, including non-standard interactions
    (NSI).

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Two-neutrino oscillations through the Sun with non-standard
    interactions:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        sth = gd.S12_NO_BF_NUFIT_6_0
        Dm2 = gd.D21_NO_BF_NUFIT_6_0
        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_2nu_sun_nsi(energy, L, L0, sth, Dm2, eps_aa=0.05, eps_ab=0.02)
        P

    .. versionadded:: 1.0.0

    .. note::
        Dispatches to a fast, closed-form interaction-picture Magnus
        integrator whenever the accumulated matter phase stays small enough
        to certify (see ``_osc_prob_ip_exp_dispatch``), giving warning-free
        results in a fraction of a second across the realistic solar-neutrino
        energy range for baselines up to a few e-folds of ``l_scale`` (the NSI
        couplings are folded into the same fast path). Longer baselines fall
        back transparently to the general slab-refinement method.

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    eps_aa : int or float, optional
        Non-universal diagonal NSI coupling of nu_e (relative to nu_mu, whose diagonal coupling is fixed to 0 by convention); see ``hamiltonians.hamiltonian_2nu_nsi``. Default: 0.0.
    eps_ab : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_2nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        sth=sth,
        Dm2=Dm2,
        eps_aa=eps_aa,
        eps_ab=eps_ab,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_3nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    for neutrinos inside the Sun.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Three-neutrino oscillations through the Sun with non-standard
    interactions, using the NuFit 6.0 defaults for the standard
    oscillation parameters:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_3nu_sun_nsi(energy, L, L0, eps_ee=0.05, eps_em=-0.03, eps_et=0.01)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_3nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        eps_ee=eps_ee,
        eps_em=eps_em,
        eps_et=eps_et,
        eps_mm=eps_mm,
        eps_mt=eps_mt,
        eps_tt=eps_tt,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability for neutrinos inside the Sun.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Four-neutrino (3+1 sterile) oscillations through the Sun with
    non-standard interactions, including the sterile-flavor couplings:

    .. jupyter-execute::

        import warnings
        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_4nu_sun_nsi(
                energy, L, L0,
                s14=0.1, s24=0.05, s34=0.02, d14=np.radians(10.0), d24=np.radians(20.0), D41=0.1,
                eps_ee=0.05, eps_em=-0.03, eps_et=0.01, eps_es=0.02, eps_ms=0.01, eps_ts=0.01)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s) NSI coupling. Default: 0.0.
    eps_ss : int or float, optional
        Diagonal NSI coupling of nu_s. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_4nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s24=s24,
        s34=s34,
        d14=d14,
        d24=d24,
        D41=D41,
        eps_ee=eps_ee,
        eps_em=eps_em,
        eps_et=eps_et,
        eps_es=eps_es,
        eps_mm=eps_mm,
        eps_mt=eps_mt,
        eps_ms=eps_ms,
        eps_tt=eps_tt,
        eps_ts=eps_ts,
        eps_ss=eps_ss,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability for neutrinos inside the Sun.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    Five-neutrino (3+2 sterile) oscillations through the Sun with
    non-standard interactions, including the sterile-flavor couplings:

    .. jupyter-execute::

        import warnings
        import numpy as np
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        energy = 15.0 * gd.UNIT_GEV  # high enough energy for the adaptive refinement
        L0 = 0.0                     # to converge cleanly under the default tolerance
        L = 0.9 * gd.SUN_RADIUS * gd.UNIT_KM

        # A trajectory through most of the Sun accumulates a large phase, so the
        # adaptive refinement needs a few loops to narrow the slabs; this is the
        # expected, informational MagnusConvergenceWarning discussed in the
        # package README, suppressed here to keep the example focused.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_5nu_sun_nsi(
                energy, L, L0,
                s14=0.1, s15=0.05, s24=0.05, s25=0.02, s34=0.02, s35=0.01,
                d14=np.radians(10.0), d15=np.radians(15.0), d24=np.radians(20.0), d35=np.radians(25.0),
                D41=0.1, D51=0.05,
                eps_ee=0.05, eps_em=-0.03, eps_et=0.01, eps_es1=0.02, eps_es2=0.01,
                eps_ms1=0.01, eps_ms2=0.01, eps_ts1=0.01, eps_ts2=0.01)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    eps_ee : int or float, optional
        Diagonal NSI coupling of nu_e. Default: 0.0.
    eps_em : int or float, optional
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling. Default: 0.0.
    eps_et : int or float, optional
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling. Default: 0.0.
    eps_es1 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s1) NSI coupling. Default: 0.0.
    eps_es2 : int or float, optional
        Flavor-off-diagonal (nu_e-nu_s2) NSI coupling. Default: 0.0.
    eps_mm : int or float, optional
        Diagonal NSI coupling of nu_mu. Default: 0.0.
    eps_mt : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling. Default: 0.0.
    eps_ms1 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s1) NSI coupling. Default: 0.0.
    eps_ms2 : int or float, optional
        Flavor-off-diagonal (nu_mu-nu_s2) NSI coupling. Default: 0.0.
    eps_tt : int or float, optional
        Diagonal NSI coupling of nu_tau. Default: 0.0.
    eps_ts1 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s1) NSI coupling. Default: 0.0.
    eps_ts2 : int or float, optional
        Flavor-off-diagonal (nu_tau-nu_s2) NSI coupling. Default: 0.0.
    eps_s1s1 : int or float, optional
        Diagonal NSI coupling of nu_s1. Default: 0.0.
    eps_s1s2 : int or float, optional
        Flavor-off-diagonal (nu_s1-nu_s2) NSI coupling. Default: 0.0.
    eps_s2s2 : int or float, optional
        Diagonal NSI coupling of nu_s2. Default: 0.0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_5nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s15=s15,
        s24=s24,
        s25=s25,
        s34=s34,
        s35=s35,
        d14=d14,
        d15=d15,
        d24=d24,
        d35=d35,
        D41=D41,
        D51=D51,
        eps_ee=eps_ee,
        eps_em=eps_em,
        eps_et=eps_et,
        eps_es1=eps_es1,
        eps_es2=eps_es2,
        eps_mm=eps_mm,
        eps_mt=eps_mt,
        eps_ms1=eps_ms1,
        eps_ms2=eps_ms2,
        eps_tt=eps_tt,
        eps_ts1=eps_ts1,
        eps_ts2=eps_ts2,
        eps_s1s1=eps_s1s1,
        eps_s1s2=eps_s1s2,
        eps_s2s2=eps_s2s2,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


#-----------------------------------------------------------------------
# In vacuum, LIV
#-----------------------------------------------------------------------

def osc_prob_2nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    sxi : int or float, optional
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_liv(
        num_flavors=2,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxiCP : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_liv(
        num_flavors=3,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_liv(
        num_flavors=4,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi15 : int or float, optional
        Sin(xi_15); see ``sxi12``. Default: 0.0.
    dxi15 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi25 : int or float, optional
        Sin(xi_25); see ``sxi12``. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    sxi35 : int or float, optional
        Sin(xi_35); see ``sxi12``. Default: 0.0.
    dxi35 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    b5 : int or float, optional
        Eigenvalue b5 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_liv(
        num_flavors=5,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, constant density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None,
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile, under (one form of)
    Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    sxi : int or float, optional
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_liv(
        num_flavors=2,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Optional[Union[int, float]], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxiCP : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_liv(
        num_flavors=3,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Optional[Union[int, float]],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_liv(
        num_flavors=4,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Optional[Union[int, float]],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    L : int, float, list, or np.ndarray
        Baseline(s).
    rho : int or float
        Matter density (or electron number density, if ``density_is_of_number_of_electrons`` is True).
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi15 : int or float, optional
        Sin(xi_15); see ``sxi12``. Default: 0.0.
    dxi15 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi25 : int or float, optional
        Sin(xi_25); see ``sxi12``. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    sxi35 : int or float, optional
        Sin(xi_35); see ``sxi12``. Default: 0.0.
    dxi35 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    b5 : int or float, optional
        Eigenvalue b5 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    
    return osc_prob_liv(
        num_flavors=5,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, exponentially falling density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_liv_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    sth: Union[int, float],
    Dm2: Union[int, float],
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    sxi : int or float, optional
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_2nu_matter_liv_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_liv(
        num_flavors=2,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_liv_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxiCP : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_3nu_matter_liv_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    return osc_prob_liv(
        num_flavors=3,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_liv_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_4nu_matter_liv_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    return osc_prob_liv(
        num_flavors=4,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_liv_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    rho_central : int or float
        Matter density (or electron number density) at the center of the exponential profile (l = 0).
    l_scale : int or float
        Length scale of the exponential density decrease.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi15 : int or float, optional
        Sin(xi_15); see ``sxi12``. Default: 0.0.
    dxi15 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi25 : int or float, optional
        Sin(xi_25); see ``sxi12``. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    sxi35 : int or float, optional
        Sin(xi_35); see ``sxi12``. Default: 0.0.
    dxi35 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    b5 : int or float, optional
        Eigenvalue b5 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    if (rho_central < 0.0 or l_scale <= 0.0):
        raise ValueError(gd.ERROR_MSG_NO_COLOR + \
            " oscprob.osc_prob_5nu_matter_liv_exp_density: rho_central must be non-negative" + \
            " and l_scale must be positive.")

    return osc_prob_liv(
        num_flavors=5,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=matter.exp_density_profile(rho_central, l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, in the Earth
#-----------------------------------------------------------------------

def osc_prob_2nu_earth_liv(
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Two-neutrino oscillations through the Earth under Lorentz-invariance
    violation:

    .. jupyter-execute::

        import warnings
        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        from magnus.magnus import MagnusConvergenceWarning

        sth = gd.S12_NO_BF_NUFIT_6_0
        Dm2 = gd.D21_NO_BF_NUFIT_6_0
        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The small solar mass splitting Dm2 combined with this Earth baseline
        # means the adaptive refinement needs a few loops; this is the expected,
        # informational MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_2nu_earth_liv(energy, sth, Dm2, sxi=0.1, b1=1.e-13, b2=2.e-13,
                                                Lambda=1.e9, n_liv=1, costhz=costhz, L=baseline)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    sxi : int or float, optional
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_earth_liv(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Three-neutrino oscillations through the Earth under Lorentz-invariance
    violation, using the NuFit 6.0 defaults for the standard oscillation
    parameters:

    .. jupyter-execute::

        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_3nu_earth_liv(
                energy, costhz=costhz, L=baseline,
                sxi12=0.1, sxi23=0.05, sxi13=0.02, dxiCP=0.3,
                b1=1.e-13, b2=2.e-13, b3=3.e-13, Lambda=1.e9, n_liv=1)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxiCP : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=3,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_earth_liv(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Four-neutrino (3+1 sterile) oscillations through the Earth under
    Lorentz-invariance violation:

    .. jupyter-execute::

        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_4nu_earth_liv(
                energy, costhz=costhz, L=baseline,
                sxi12=0.1, sxi23=0.05, sxi13=0.02, dxi13=0.2,
                sxi14=0.05, dxi14=0.4, sxi24=0.03, dxi24=0.5, sxi34=0.02,
                b1=1.e-13, b2=2.e-13, b3=3.e-13, b4=1.e-13, Lambda=1.e9, n_liv=1)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=4,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_earth_liv(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.

    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    .. jupyter-execute::

        import magnus.earth as earth

        list(earth.loc_coords_dms.keys())

    .. jupyter-execute::

        print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

    Examples
    --------
    Five-neutrino (3+2 sterile) oscillations through the Earth under
    Lorentz-invariance violation:

    .. jupyter-execute::

        import magnus.oscprob as oscprob
        import magnus.globaldefs as gd
        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        costhz = -0.8
        baseline = 2.0 * gd.EARTH_RADIUS * 0.8 * gd.UNIT_KM
        energy = 1.0 * gd.UNIT_GEV

        # The chosen baseline/energy combination needs a few adaptive-refinement
        # loops to converge; this is the expected, informational
        # MagnusConvergenceWarning discussed in the package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            P = oscprob.osc_prob_5nu_earth_liv(
                energy, costhz=costhz, L=baseline,
                sxi12=0.1, sxi23=0.05, sxi13=0.02, dxi13=0.2,
                sxi14=0.05, dxi14=0.4, sxi15=0.02, dxi15=0.6,
                sxi24=0.03, dxi24=0.5, sxi25=0.01, sxi34=0.02, sxi35=0.01, dxi35=0.7,
                b1=1.e-13, b2=2.e-13, b3=3.e-13, b4=1.e-13, b5=1.e-13, Lambda=1.e9, n_liv=1)
        P

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : int, float, list, or np.ndarray
        Neutrino energy/energies.
    costhz : int or float, optional
        Cosine of the zenith angle of the neutrino. Used together with ``L``, as an alternative to ``loc_ini``/``loc_fin``. Default: None.
    loc_ini : tuple, list, np.ndarray, or str, optional
        Initial location on the surface of the Earth, as (latitude, longitude) coordinates or a predefined location name (see ``earth.loc_coords_dms``). Must be given with ``loc_fin``. Default: None.
    loc_fin : tuple, list, np.ndarray, or str, optional
        Final location, same format as ``loc_ini``. Must be given with ``loc_ini``. Default: None.
    L : float, list, or np.ndarray, optional
        Baseline(s). Default: None.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi15 : int or float, optional
        Sin(xi_15); see ``sxi12``. Default: 0.0.
    dxi15 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi25 : int or float, optional
        Sin(xi_25); see ``sxi12``. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    sxi35 : int or float, optional
        Sin(xi_35); see ``sxi12``. Default: 0.0.
    dxi35 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    b5 : int or float, optional
        Eigenvalue b5 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    default_osc_params_set_name : str, optional
        Name of the predefined oscillation-parameter set used to fill in any oscillation parameter left as None (see ``globaldefs.OSC_PARAMS_PREDEFINED``). Default: 'OSC_PARAMS_DEFAULT'.
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=5,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, in the Sun
#-----------------------------------------------------------------------

def osc_prob_2nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    sth: Union[int, float],
    Dm2: Union[int, float],
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    sth : int or float
        Sine of the mixing angle :math:`\theta`, the single mixing angle of the two-flavor system.
    Dm2 : int or float
        Mass-squared difference :math:`\Delta m^2` of the two-flavor system.
    sxi : int or float, optional
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_2nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        sth=sth,
        Dm2=Dm2,
        sxi=sxi,
        b1=b1,
        b2=b2,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_3nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxiCP : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_3nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        sxi12=sxi12,
        sxi23=sxi23,
        sxi13=sxi13,
        dxiCP=dxiCP,
        b1=b1,
        b2=b2,
        b3=b3,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_4nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        s14=s14,
        s24=s24,
        s34=s34,
        d14=d14,
        d24=d24,
        D41=D41,
        sxi12=sxi12,
        sxi23=sxi23,
        sxi13=sxi13,
        dxi13=dxi13,
        sxi14=sxi14,
        dxi14=dxi14,
        sxi24=sxi24,
        dxi24=dxi24,
        sxi34=sxi34,
        b1=b1,
        b2=b2,
        b3=b3,
        b4=b4,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    strategy: Optional[str]='auto',
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float, list, or np.ndarray
        Neutrino energy/energies.
    L : float, list, or np.ndarray
        Baseline(s).
    L0 : int or float
        Initial position.
    s12 : int or float, optional
        Sine of the mixing angle :math:`\theta_{12}`. Default: None.
    s23 : int or float, optional
        Sine of the mixing angle :math:`\theta_{23}`. Default: None.
    s13 : int or float, optional
        Sine of the mixing angle :math:`\theta_{13}`. Default: None.
    dCP : int or float, optional
        :math:`\delta_\text{CP}` [radian]. Default: None.
    D21 : int or float, optional
        Mass-squared difference :math:`\Delta m_{21}^2`. Default: None.
    D31 : int or float, optional
        Mass-squared difference :math:`\Delta m_{31}^2`. Default: None.
    s14 : int or float, optional
        Sine of the mixing angle :math:`\theta_{14}`. Default: 0.0.
    s15 : int or float, optional
        Sine of the mixing angle :math:`\theta_{15}`. Default: 0.0.
    s24 : int or float, optional
        Sine of the mixing angle :math:`\theta_{24}`. Default: 0.0.
    s25 : int or float, optional
        Sine of the mixing angle :math:`\theta_{25}`. Default: 0.0.
    s34 : int or float, optional
        Sine of the mixing angle :math:`\theta_{34}`. Default: 0.0.
    s35 : int or float, optional
        Sine of the mixing angle :math:`\theta_{35}`. Default: 0.0.
    d14 : int or float, optional
        :math:`\delta_{14}` [radian]. Default: 0.0.
    d15 : int or float, optional
        :math:`\delta_{15}` [radian]. Default: 0.0.
    d24 : int or float, optional
        :math:`\delta_{24}` [radian]. Default: 0.0.
    d35 : int or float, optional
        :math:`\delta_{35}` [radian]. Default: 0.0.
    D41 : int or float, optional
        Mass-squared difference :math:`\Delta m_{41}^2`. Default: 0.0.
    D51 : int or float, optional
        Mass-squared difference :math:`\Delta m_{51}^2`. Default: 0.0.
    sxi12 : int or float, optional
        Sin(xi_12), one of the mixing angles between the space of the eigenvectors of the LIV operator and the flavor states. Default: 0.0.
    sxi23 : int or float, optional
        Sin(xi_23); see ``sxi12``. Default: 0.0.
    sxi13 : int or float, optional
        Sin(xi_13); see ``sxi12``. Default: 0.0.
    dxi13 : int or float, optional
        CP-violation phase of the LIV operator [radian] (replaces ``dxiCP`` for 4/5-flavor systems). Default: 0.0.
    sxi14 : int or float, optional
        Sin(xi_14); see ``sxi12``. Default: 0.0.
    dxi14 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi15 : int or float, optional
        Sin(xi_15); see ``sxi12``. Default: 0.0.
    dxi15 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi24 : int or float, optional
        Sin(xi_24); see ``sxi12``. Default: 0.0.
    dxi24 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    sxi25 : int or float, optional
        Sin(xi_25); see ``sxi12``. Default: 0.0.
    sxi34 : int or float, optional
        Sin(xi_34); see ``sxi12``. Default: 0.0.
    sxi35 : int or float, optional
        Sin(xi_35); see ``sxi12``. Default: 0.0.
    dxi35 : int or float, optional
        CP-violation phase of the LIV operator [radian]. Default: 0.0.
    b1 : int or float, optional
        Eigenvalue b1 of the LIV operator. Default: 0.0.
    b2 : int or float, optional
        Eigenvalue b2 of the LIV operator. Default: 0.0.
    b3 : int or float, optional
        Eigenvalue b3 of the LIV operator. Default: 0.0.
    b4 : int or float, optional
        Eigenvalue b4 of the LIV operator. Default: 0.0.
    b5 : int or float, optional
        Eigenvalue b5 of the LIV operator. Default: 0.0.
    Lambda : int or float, optional
        Energy scale of the LIV operator. Default: 1.0.
    n_liv : int, optional
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3). Default: 0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, compute the probability for antineutrinos. Default: False.
    nu_i : int, optional
        Initial flavor index. If given together with ``nu_f``, a single channel is returned instead of the full probability matrix. Default: None.
    nu_f : int, optional
        Final flavor index; see ``nu_i``. Default: None.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, the density is given in :math:`\text{g cm}^{-3}`. Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, the density parameter directly gives the electron number density [:math:`\text{eV}^{3}`]. Default: False.
    strategy : str, optional
        Numerical strategy used to compute the evolution operator: 'auto' (default),
        'hybrid', or 'magnus'; see the ``strategy`` parameter of
        :func:`osc_prob_matter_std_potential` for the full description and
        :doc:`/adiabatic_strategy` for the derivation and validation. Default: 'auto'.

        .. versionadded:: 1.0.0
    validate_input : bool, optional
        If True, validate the input parameters. Default: True.
    save_log : bool, optional
        If True, also write log messages to a file. Default: False.
    filename_log : str, optional
        Name of the log file (used if ``save_log`` is True and no ``file_log`` object is given). Default: './out.log'.
    file_log : TextIOWrapper, optional
        Optional file object to write log messages to. Default: None.
    close_file_log_upon_exit : bool, optional
        If True, close the log file before returning. Default: True.
    verbose : int, optional
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the refinement loops). Default: 0.
    \**kwargs
        Additional arguments forwarded to the underlying middle-layer function (e.g., the standard refinement/logging kwargs; see :func:`osc_prob`).

    Returns
    -------
    float or np.ndarray
        Oscillation probability matrix (or single channel, if ``nu_i``/``nu_f`` are given) for each (energy, L) point.
    """

    return osc_prob_5nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        s14=s14,
        s15=s15,
        s24=s24,
        s25=s25,
        s34=s34,
        s35=s35,
        d14=d14,
        d15=d15,
        d24=d24,
        d35=d35,
        D41=D41,
        D51=D51,
        sxi12=sxi12,
        sxi23=sxi23,
        sxi13=sxi13,
        dxi13=dxi13,
        sxi14=sxi14,
        dxi14=dxi14,
        sxi15=sxi15,
        dxi15=dxi15,
        sxi24=sxi24,
        dxi24=dxi24,
        sxi25=sxi25,
        sxi34=sxi34,
        sxi35=sxi35,
        dxi35=dxi35,
        b1=b1,
        b2=b2,
        b3=b3,
        b4=b4,
        b5=b5,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        strategy=strategy,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )



from .oscprobstd import (
    osc_prob_2nu_vacuum_std,
    osc_prob_2nu_matter_std,
    delta,
    J,
    osc_prob_3nu_vacuum_std,
)

__all__ = [
    # re-exported from oscprobstd.py
    'osc_prob_2nu_vacuum_std',
    'osc_prob_2nu_matter_std',
    'delta',
    'J',
    'osc_prob_3nu_vacuum_std',
    'MAX_N_SLABS_DEFAULT',
    'ToleranceNotAchievedWarning',
    'HybridCertificationWarning',
    'UnmarkedDiscontinuityWarning',
    'HiddenFeatureWarning',
    'ENGINE_FAMILIES',
    'cross_check_strategies',
    'print_banner',
    'print_run_parameters',
    'validate_input_battery',
    'validate_input_osc_prob_earth',
    'valid_flavor_indices_2nu',
    'values_to_unspecified_osc_params',
    'unpack_oscillation_params_from_dict',
    'unpack_nsi_params_from_dict',
    'unpack_liv_params_from_dict',
    'compute_evolution_operator',
    'compute_evolution_operator_multiple_slabs',
    'osc_prob',
    'osc_prob_energy_baseline',
    'osc_prob_vacuum',
    'osc_prob_matter_std_potential',
    'osc_prob_matter_nsi',
    'osc_prob_liv',
    'osc_prob_2nu_vacuum',
    'osc_prob_3nu_vacuum',
    'osc_prob_4nu_vacuum',
    'osc_prob_5nu_vacuum',
    'osc_prob_2nu_matter_constant_density',
    'osc_prob_3nu_matter_constant_density',
    'osc_prob_4nu_matter_constant_density',
    'osc_prob_5nu_matter_constant_density',
    'osc_prob_2nu_matter_exp_density',
    'osc_prob_3nu_matter_exp_density',
    'osc_prob_4nu_matter_exp_density',
    'osc_prob_5nu_matter_exp_density',
    'osc_prob_2nu_earth',
    'osc_prob_3nu_earth',
    'osc_prob_4nu_earth',
    'osc_prob_5nu_earth',
    'osc_prob_earth',
    'osc_prob_2nu_sun',
    'osc_prob_3nu_sun',
    'osc_prob_4nu_sun',
    'osc_prob_5nu_sun',
    'osc_prob_sun',
    'osc_prob_2nu_matter_nsi_constant_density',
    'osc_prob_3nu_matter_nsi_constant_density',
    'osc_prob_4nu_matter_nsi_constant_density',
    'osc_prob_5nu_matter_nsi_constant_density',
    'osc_prob_2nu_matter_nsi_exp_density',
    'osc_prob_3nu_matter_nsi_exp_density',
    'osc_prob_4nu_matter_nsi_exp_density',
    'osc_prob_5nu_matter_nsi_exp_density',
    'osc_prob_2nu_earth_nsi',
    'osc_prob_3nu_earth_nsi',
    'osc_prob_4nu_earth_nsi',
    'osc_prob_5nu_earth_nsi',
    'osc_prob_2nu_sun_nsi',
    'osc_prob_3nu_sun_nsi',
    'osc_prob_4nu_sun_nsi',
    'osc_prob_5nu_sun_nsi',
    'osc_prob_2nu_vacuum_liv',
    'osc_prob_3nu_vacuum_liv',
    'osc_prob_4nu_vacuum_liv',
    'osc_prob_5nu_vacuum_liv',
    'osc_prob_2nu_matter_liv_constant_density',
    'osc_prob_3nu_matter_liv_constant_density',
    'osc_prob_4nu_matter_liv_constant_density',
    'osc_prob_5nu_matter_liv_constant_density',
    'osc_prob_2nu_matter_liv_exp_density',
    'osc_prob_3nu_matter_liv_exp_density',
    'osc_prob_4nu_matter_liv_exp_density',
    'osc_prob_5nu_matter_liv_exp_density',
    'osc_prob_2nu_earth_liv',
    'osc_prob_3nu_earth_liv',
    'osc_prob_4nu_earth_liv',
    'osc_prob_5nu_earth_liv',
    'osc_prob_2nu_sun_liv',
    'osc_prob_3nu_sun_liv',
    'osc_prob_4nu_sun_liv',
    'osc_prob_5nu_sun_liv',
]
