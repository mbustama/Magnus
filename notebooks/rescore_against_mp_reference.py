# -*- coding: utf-8 -*-
r"""rescore_against_mp_reference.py

Re-score every series in ``external_profile_benchmarks.json`` against the extended-
precision reference, replacing the double-precision DOP853 one.

WHY.  Measured against itself at a tenth of its tolerance, the DOP853 referee is
uncertain by 1.1e-13 (3nu), 2.3e-11 (4nu) and 5.2e-11 (5nu) -- 1.1x the floor every
series in Figure 11 settles onto.  The bottom of those curves was the ruler, not the
codes.  ``mp_reference_profile.json`` replaces it with a triple-Richardson mpmath ladder
converged to 1e-15 or better.

WHAT MOVES AND WHAT DOES NOT.  Only ``max_abs_error``.  Every ``us_per_probability``
stays exactly as measured -- the probabilities are recomputed in untimed calls, so no
timing is disturbed and the figure's single time axis remains one measurement session
per point.  The original errors are kept beside the new ones as ``max_abs_error_ode``
so the change is auditable rather than silent.

    python notebooks/rescore_against_mp_reference.py
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import json
import pathlib
import shutil
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import gen_profile_benchmarks as gpb                      # noqa: E402

from magnus import globaldefs as gd                       # noqa: E402

HERE = pathlib.Path(__file__).parent
BENCH = HERE/'external_profile_benchmarks.json'
MPREF = HERE/'mp_reference_profile.json'


def magnus_P(d, prof, rtol, order):
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
    return np.asarray(gpb.oscprob.osc_prob_matter_std_potential(
        d, lambda x: prof['vcc'](x)/per_ne, prof['energies'], prof['baseline'],
        gpb.osc_params(d), L0=0.0, nu_i=gd.NUMU, nu_f=gd.NUMU,
        density_is_of_number_of_electrons=True, rtol=rtol, atol=rtol*1.0e-2,
        magnus_exp_order=order, strategy='magnus'))


def main():
    bench = json.loads(BENCH.read_text())
    mpref = json.loads(MPREF.read_text())
    shutil.copy(BENCH, BENCH.with_suffix('.json.pre_mpref'))
    prof = gpb.exponential_profile()
    by_d = {c['flavours']: c for c in mpref['cases']}

    for case in bench['cases']:
        d = case['flavours']
        if d not in by_d:
            print('  %dnu: no mpmath reference, left alone' % d, file=sys.stderr)
            continue
        ref = np.array(by_d[d]['reference'])
        case['reference_ode'] = case['reference']
        case['reference'] = by_d[d]['reference']
        case['reference_self_convergence'] = by_d[d]['self_convergence']
        case['reference_ode_crosscheck'] = by_d[d]['ode_crosscheck']
        print('--- %d flavours ---' % d, file=sys.stderr, flush=True)
        for s in case['series']:
            if s['dial'] != 'rtol':
                continue                       # NuOscProbExact rescored below
            order = s.get('magnus_exp_order', 4)
            for p in s['points']:
                P = magnus_P(d, prof, p['rtol'], order)
                p['max_abs_error_ode'] = p['max_abs_error']
                p['max_abs_error'] = float(np.max(np.abs(P - ref)))
            print('  %-18s %s' % (s['name'], '  '.join(
                '%.2e' % p['max_abs_error'] for p in s['points'])),
                file=sys.stderr, flush=True)
        npe = next((s for s in case['series'] if s['dial'] == 'n_slabs'), None)
        if npe is not None:
            pts = gpb.npe_points(d, prof, [p['n_slabs'] for p in npe['points']], ref)
            for old, new in zip(npe['points'], pts):
                old['max_abs_error_ode'] = old['max_abs_error']
                old['max_abs_error'] = new['max_abs_error']
            print('  %-18s %s' % ('NuOscProbExact', '  '.join(
                '%.2e' % p['max_abs_error'] for p in npe['points'])),
                file=sys.stderr, flush=True)

    bench['referee'] = ('mp_reference_profile.json -- triple-Richardson mpmath; the '
                        'former DOP853 errors are kept as max_abs_error_ode')
    BENCH.write_text(json.dumps(bench, indent=1))
    print('wrote %s' % BENCH.name, file=sys.stderr)


if __name__ == '__main__':
    main()
