# -*- coding: utf-8 -*-
r"""retime_magnus_series.py

**PARTLY UNVERIFIED -- read this before running.**  The Mag(nu)s half is correct
and was checked.  The NuOscProbExact half is NOT: it routes every series through
one helper, while the stored series were built through three different routes
(see ``npe_rtol_series.route`` in each file).  Re-timing through the wrong route
produced NuOscProbExact 1.2-2.7x slower on the smooth profile and 0.42-1.05x on
PREM, which is a different measurement rather than a re-measurement.  Fixing that
is point 1 of the blocked note in resources/paper/pending-edits.md.

Re-time Mag(nu)s's series in the two files Figure 11 reads, and nothing else.

WHY.  Between v1.0.6 and v1.0.12 five compiled kernels landed -- the separable
scan's slab product, the commutator, the Gauss-Legendre Omega at orders 4, 6 and
8, the anti-Hermiticity framing, a Jacobi eigensolver at four and five flavours,
and the interaction-picture and baseline-scan folds.  The Mag(nu)s points in
``external_profile_benchmarks.json`` and ``external_prem_chord_benchmarks.json``
predate all of them, so Figure 11 shows a code that no longer exists.

WHAT THIS DOES NOT DO.  It does not touch the DOP853 or mpmath references, the
self-convergence records, or the energy grids.  Every point -- Mag(nu)s and
NuOscProbExact alike -- is re-run at *the same dial setting* it was run at before
and scored against *the same stored reference*, so only the clock is new.  The
accuracy column moves only insofar as a code's own accuracy did.

WHY NuOscProbExact IS RE-TIMED TOO, having first been left alone.  The stored
points date from 2026-08-10.  Re-timing one code and not the other is only honest
if the machine still times the same way, and ``probe_commensurability.py`` reported
its control ratio at 0.986 -> 1.011, +2.5%, which looked like a pass.  It is not
the right test for these points.  That control is a 180x180 matmul and is
BLAS-bound; the loose-tolerance points here are dominated by Python-level overhead,
and the two drift independently.  Measured directly: the pre-kernel code at
a55b8a4, run today, is **12-20% slower than its own stored numbers** on this
configuration.  So the stored NuOscProbExact points cannot share a time axis with
anything measured now, and re-timing Mag(nu)s alone would have flattered the
closed form by that much -- comparable to the whole effect at the loose end of the
dial, where Mag(nu)s gains only about 1.2x.  Re-timing both costs about a minute.

COMMENSURABILITY.  Figure 11 puts every series on one time axis and its caption
says both codes were timed in one process on one machine.  Re-timing one code in
a later session keeps that honest only if the machine still times the same way,
which ``probe_commensurability.py`` asks directly.  Run it first.  Its control
ratio is the test: it exercises a workload none of these kernels touch.  Its
second check, a re-timed Mag(nu)s point, will be several times faster now and
that is the change being recorded, not drift -- do not read its verdict as a
veto without looking at which line moved.

Run on an idle machine.  Rewrites both files in place after backing each up::

    python notebooks/retime_magnus_series.py
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import json
import pathlib
import platform
import shutil
import sys
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import gen_profile_benchmarks as gpb                      # noqa: E402
import gen_prem_benchmarks as gpb_prem                    # noqa: E402
import append_order_series as aos                         # noqa: E402
import prem_chord_common as pcc                           # noqa: E402
import append_npe_rtol_series as anrs                     # noqa: E402
import append_npe_rtol_prem as anrp                       # noqa: E402

HERE = pathlib.Path(__file__).parent
PROFILE_STORE = HERE/'external_profile_benchmarks.json'
PREM_STORE = HERE/'external_prem_chord_benchmarks.json'
MAGNUS_ORDERS = {'Magnus': 4, 'Magnus, order 6': 6, 'Magnus, order 8': 8}


def control_ratio():
    """The stored files' own interleaved control, so this run carries its own figure."""
    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], gpb.timed(gpb.control, repeat=1))
    return best['a']/best['b']


def retime(store_path, point_fn, prof_fn, npe_fn, label):
    store = json.loads(store_path.read_text())
    backup = store_path.with_suffix('.json.pre_kernels')
    if backup.exists():
        print('%s: backup %s already exists, keeping it' % (label, backup.name),
              file=sys.stderr)
    else:
        shutil.copy(store_path, backup)
        print('%s: backup -> %s' % (label, backup.name), file=sys.stderr)

    prof = prof_fn()
    with warnings.catch_warnings():
        # The convergence warning fires on these profiles at every order, including for
        # the stored points; it is a slab-width statement and the accuracy column is what
        # says whether it mattered.
        warnings.simplefilter('once')
        for case in store['cases']:
            d = case['flavours']
            ref = np.array(case['reference'])
            print('  %s %d flavours' % (label, d), file=sys.stderr, flush=True)
            for series in case['series']:
                print('    %s' % series['name'], file=sys.stderr, flush=True)
                order = MAGNUS_ORDERS.get(series['name'])
                if order is not None:
                    dials = [p['rtol'] for p in series['points']]
                    series['points'] = point_fn(d, prof, dials, ref, order)
                    continue
                npe = npe_fn(d, prof, series, ref)
                if npe:
                    series['points'] = npe
    return store, backup


def main():
    ratio = control_ratio()
    print('control ratio for this run: %.3f' % ratio, file=sys.stderr)

    note = {'date': '2026-09-05',
            'why': ('Magnus re-timed after the compiled kernels of v1.0.6-v1.0.12; '
                    'NuOscProbExact, the references and the grids are untouched, and '
                    'every point was re-run at its own stored tolerance against its own '
                    'stored reference'),
            'control_ratio_this_run': ratio,
            'magnus_version': '1.0.12',
            'machine': platform.platform()}

    chord = pcc.chord()
    prem_energies = np.asarray(json.loads(
        (HERE/'prem_chord_reference.json').read_text())['energy_ev'], dtype=float)

    def prem_npe(d, prof, series, ref):
        dials = [p['n_slabs'] for p in series['points']] if series['dial'] == 'n_slabs' \
            else [p['rtol'] for p in series['points']]
        if series['dial'] == 'n_slabs':
            return gpb.npe_points(d, prof, dials, ref)
        return anrp.npe_rtol_points(d, prem_energies, dials, ref)

    def profile_npe(d, prof, series, ref):
        if series['dial'] == 'n_slabs':
            return gpb.npe_points(d, prof, [p['n_slabs'] for p in series['points']], ref)
        return anrs.npe_rtol_points(d, prof, [p['rtol'] for p in series['points']], ref)

    prem_store, _ = retime(
        PREM_STORE,
        lambda d, prof, dials, ref, order: gpb_prem.magnus_points(
            d, prof, dials, ref, order, chord['edges'][1:-1]),
        lambda: {'vcc': chord['vcc'], 'baseline': chord['baseline'],
                 'energies': prem_energies},
        prem_npe, 'PREM')
    prem_store['magnus_rerun_kernels'] = note
    PREM_STORE.write_text(json.dumps(prem_store, indent=1) + '\n')
    print('PREM: written', file=sys.stderr)

    prof_store, _ = retime(
        PROFILE_STORE,
        aos.magnus_points_at_order,
        gpb.exponential_profile,
        profile_npe, 'profile')
    prof_store['magnus_rerun_kernels'] = note
    PROFILE_STORE.write_text(json.dumps(prof_store, indent=1) + '\n')
    print('profile: written', file=sys.stderr)


if __name__ == '__main__':
    main()
