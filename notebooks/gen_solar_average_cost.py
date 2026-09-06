# -*- coding: utf-8 -*-
r"""gen_solar_average_cost.py

What the *averaged* solar probability costs, per physics configuration, on the
real BS2005-AGS,OP model.  Writes ``external_solar_average_cost.json``.

    python notebooks/gen_solar_average_cost.py > notebooks/external_solar_average_cost.json

WHY COST AND NOT SPEED-AGAINST-ACCURACY.  On this profile the averaged
probability has no accuracy dial.  Every eigenvalue pair is fully decohered and
the evolution is adiabatic throughout -- the adiabaticity parameter stays a
factor of eleven or more under the detector's threshold at every flavour count
from 0.1 to 20 MeV -- so the crossing matrix is the identity and the answer
reduces to two eigendecompositions, exact to ~1e-15.  Sweeping ``n_points``,
``threshold``, ``n_probe``, ``magnus_exp_order`` or ``integration_method``
returns bit-identical output.  There is therefore nothing to put on an accuracy
axis, and what remains worth showing is the price of the answer across the
scenarios the package covers.  See ``docs/dev/SCOPE_FIG13.md``.

THE PROFILE IS THE TABULATED MODEL, NOT THE EXPONENTIAL FIT.  ``osc_prob_*_sun``
uses the fit; this script uses ``gen_profile_benchmarks.solar_profile()``, which
log-interpolates the BS2005-AGS,OP table.  The fit is high by a factor 2.4 inside
0.05 R_sun, so they are different problems rather than two versions of one -- and
they cost differently: 5-8 ms on the fit against 8-40 ms here.

Timings are best-of-N wall clock with the first call discarded (it carries the
numba compile), and an interleaved control workload is recorded so a later run
can say whether this machine still times the same way.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import json
import pathlib
import platform
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_profile_benchmarks as gpb                      # noqa: E402
import magnus.globaldefs as gd                            # noqa: E402
import magnus.oscprob as oscprob                          # noqa: E402

ENERGY = 5.0e6                                            # 5 MeV, the solar workhorse
REPEAT = 7

# Large enough to matter, small enough to stay a perturbation -- the same scale
# gen_shock_nsi.py uses, so the two figures' NSI mean the same thing.
_EPS_ACTIVE = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0,
                   eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
# The sterile sector carries no NSI of its own: these rows ask what the *same*
# active-sector NSI costs once there are more states to carry it, not what a
# richer NSI model costs.  Every added entry is therefore zero, and the names
# differ by flavour count because the wrappers unpack them positionally.
_EPS_STERILE = {4: ('eps_es', 'eps_ms', 'eps_ts', 'eps_ss'),
                5: ('eps_es1', 'eps_es2', 'eps_ms1', 'eps_ms2', 'eps_ts1',
                    'eps_ts2', 'eps_s1s1', 'eps_s1s2', 'eps_s2s2')}


def nsi_params_for(d):
    eps = dict(_EPS_ACTIVE)
    eps.update({k: 0.0 for k in _EPS_STERILE.get(d, ())})
    return eps
LIV = dict(sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0,
           b1=gd.B1, b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1)
STERILE_KEYS = ('s14', 's24', 's34', 'd14', 'd24', 'D41')


def timed(call):
    """Best of REPEAT, first pass discarded: it carries the one-off numba compile."""
    call()
    best = np.inf
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        call()
        best = min(best, time.perf_counter() - t0)
    return best


def control():
    a = np.random.default_rng(0).normal(size=(180, 180))
    return a @ a


def main():
    warnings.simplefilter('ignore')
    prof = gpb.solar_profile()
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
    L = prof['baseline']

    def ne_of(x):
        return prof['vcc'](x)/per_ne

    common = dict(L0=0.0, density_is_of_number_of_electrons=True,
                  average=True, validate_input=False)

    def std(d):
        return lambda: oscprob.osc_prob_matter_std_potential(
            d, ne_of, ENERGY, L, gpb.osc_params(d), **common)

    def nsi(d):
        return lambda: oscprob.osc_prob_matter_nsi(
            d, ne_of, ENERGY, L, gpb.osc_params(d), nsi_params_for(d), **common)

    def liv(d):
        # rho_func is the sixth positional argument; omitting it silently gives
        # vacuum LIV, which runs a hundred times faster and is a different case.
        return lambda: oscprob.osc_prob_liv(
            d, ENERGY, L, gpb.osc_params(d), LIV, ne_of, **common)

    cases = [('2nu', 2, std(2)), ('3nu', 3, std(3)),
             ('3+1', 4, std(4)), ('3+2', 5, std(5)),
             ('3nu + NSI', 3, nsi(3)), ('3nu + LIV', 3, liv(3)),
             ('3+1 + NSI', 4, nsi(4)), ('3+2 + NSI', 5, nsi(5))]

    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], timed(control))
    out = {'note': 'Cost of the averaged probability, BS2005-AGS,OP, per configuration',
           'profile': 'BS2005-AGS,OP, log-interpolated (NOT the exponential fit)',
           'energy_ev': ENERGY, 'baseline_inv_ev': L,
           'why_no_accuracy_axis': ('adiabatic throughout and fully decohered, so the '
                                    'crossing matrix is the identity and no knob moves '
                                    'the value; see docs/dev/SCOPE_FIG13.md'),
           'machine': platform.platform(), 'python': platform.python_version(),
           'numpy': np.__version__, 'magnus': getattr(__import__('magnus'), '__version__', '?'),
           'control_ratio': best['a']/best['b'], 'repeat': REPEAT, 'cases': []}

    for label, d, call in cases:
        try:
            P = np.asarray(call())
        except Exception as exc:                          # noqa: BLE001
            print('  %-11s SKIPPED: %s' % (label, exc), file=sys.stderr, flush=True)
            continue
        t = timed(call)
        row = dict(label=label, flavors=d, seconds=t, ms=1.0e3*t,
                   p_ee=float(np.asarray(P)[0, 0]),
                   unitarity=float(np.max(np.abs(np.asarray(P).sum(axis=1) - 1.0))))
        out['cases'].append(row)
        print('  %-11s %8.2f ms   P_ee = %.6f   unitarity %.1e'
              % (label, row['ms'], row['p_ee'], row['unitarity']),
              file=sys.stderr, flush=True)

    json.dump(out, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
