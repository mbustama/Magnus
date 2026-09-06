# -*- coding: utf-8 -*-
r"""append_order_series.py

Add Mag(nu)s order-6 and order-8 series to ``external_profile_benchmarks.json``,
which Figure 11 of the paper reads.

WHAT THIS DOES NOT DO.  It does not re-time anything already in the file.  The
Mag(nu)s order-4 points, the NuOscProbExact points and the DOP853 references stay
exactly as they were measured on 2026-08-10; the new points are scored against those
same stored references, so accuracy is compared on one ruler and only the clock is new.

WHY THAT IS LEGITIMATE.  The figure puts every series on one time axis and its caption
says the codes were timed in one process on one machine.  Appending timings from a later
session only keeps that claim honest if the machine still times the same way, which is a
measurable question rather than an assumption.  ``probe_commensurability.py`` asked it
first: the interleaved control moved 0.986 -> 1.003, one stored Mag(nu)s point re-timed
to within -0.5%, and its error reproduced exactly (2.929e-13).  This script records its
own control ratio alongside the appended points so the same check is auditable later.

    python notebooks/append_order_series.py

Run it on an idle machine.  It rewrites the JSON in place, after saving a backup beside
it, and prints progress to stderr.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import datetime
import json
import pathlib
import platform
import shutil
import sys
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import gen_profile_benchmarks as gpb                      # noqa: E402

HERE = pathlib.Path(__file__).parent
STORE = HERE/'external_profile_benchmarks.json'
NEW_ORDERS = (6, 8)


def magnus_points_at_order(d, prof, dials, ref, order):
    """(dial, time per probability, worst error) at a fixed magnus_exp_order.

    Deliberately not gen_profile_benchmarks.magnus_points(): that one recomputes the
    DOP853 referee, which the stored file already carries and which must not move.
    Everything else -- the call, the batching over energies, the timing protocol -- is
    the same, so the new points sit on the same axes as the old ones.
    """
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)

    def ne_of(x):
        return prof['vcc'](x)/per_ne

    E = prof['energies']
    osc = gpb.osc_params(d)
    out = []
    for rtol in dials:
        def call(r=rtol):
            return np.asarray(gpb.oscprob.osc_prob_matter_std_potential(
                d, ne_of, E, prof['baseline'], osc, L0=0.0,
                nu_i=gpb.gd.NUMU, nu_f=gpb.gd.NUMU,
                density_is_of_number_of_electrons=True, rtol=r, atol=r*1.0e-2,
                magnus_exp_order=order, strategy='magnus'))

        P = call()
        us = 1.0e6*gpb.timed(call)/len(E)
        err = float(np.max(np.abs(P - ref)))
        out.append(dict(label='%.0e' % rtol, rtol=rtol,
                        us_per_probability=us, max_abs_error=err))
        print('      rtol %.0e  %10.2f us/prob  err %.3e' % (rtol, us, err),
              file=sys.stderr, flush=True)
    return out


def main():
    store = json.loads(STORE.read_text())
    shutil.copy(STORE, STORE.with_suffix('.json.bak'))
    print('backup: %s' % STORE.with_suffix('.json.bak').name, file=sys.stderr)

    # The same interleaved control the original run used, so the appended points carry
    # their own contention figure rather than inheriting one measured a month earlier.
    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], gpb.timed(gpb.control, repeat=1))
    ratio = best['a']/best['b']
    print('control ratio for this run: %.3f  (stored run: %.3f)'
          % (ratio, store['control_ratio']), file=sys.stderr)

    prof = gpb.exponential_profile()
    with warnings.catch_warnings():
        # The convergence warning fires on this profile at every order, including the
        # stored order-4 points; it is a slab-width statement, not an error, and the
        # accuracy column is what says whether it mattered.
        warnings.simplefilter('once')
        for case in store['cases']:
            d = case['flavours']
            ref = np.array(case['reference'])
            have = {s['name'] for s in case['series']}
            print('  %d flavours' % d, file=sys.stderr, flush=True)
            for order in NEW_ORDERS:
                name = 'Magnus, order %d' % order
                if name in have:
                    print('    %s already present, skipping' % name, file=sys.stderr)
                    continue
                print('    %s' % name, file=sys.stderr, flush=True)
                pts = magnus_points_at_order(d, prof, gpb.MG_DIALS, ref, order)
                case['series'].append({'name': name, 'dial': 'rtol',
                                       'magnus_exp_order': order, 'points': pts})

    store['appended'] = {
        'what': ('Magnus order-6 and order-8 series, scored against the stored '
                 'references; nothing already in this file was re-timed'),
        'orders': list(NEW_ORDERS),
        'unlabelled_magnus_series_is_order': 4,
        'date': datetime.date.today().isoformat(),
        'machine': platform.platform(),
        'control_ratio': ratio,
        'commensurability_probe': ('control 0.986 -> 1.003; one stored order-4 point '
                                   're-timed to -0.5% with its error reproducing '
                                   'exactly at 2.929e-13'),
    }
    STORE.write_text(json.dumps(store, indent=1))
    print('wrote %s' % STORE.name, file=sys.stderr)


if __name__ == '__main__':
    main()
