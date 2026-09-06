# -*- coding: utf-8 -*-
r"""append_npe_rtol_series.py

Add a tolerance-dialled NuOscProbExact series to the smooth-profile benchmarks, beside
the slab-count one already there.

WHY.  Figure 11 dials Mag(nu)s by a requested tolerance and NuOscProbExact by a slab
count, so the two curves are parameterised by different knobs and a reader has to take
on trust that they are comparable.  NuOscProbExact 1.14.0 refines to a tolerance and
accepts a batched Hamiltonian, so both codes can now answer the same request: "you asked
for 1e-8; here is what it cost and what you got."

WHAT IS KEPT.  The existing ``NuOscProbExact`` series, dialled by ``n_slabs``, stays in
the file untouched.  It is the measurement the published figure was drawn from, and the
new series is added beside it rather than over it.  Which one the figure plots is the
figure's business.

BATCHING MATTERS AND IS NOT DECORATION.  The slab-count series times twelve energies in
one call; the generator's own note says timing them one at a time "would flatter
Mag(nu)s roughly fivefold".  1.14.0 accepts ``(n_energies, n, d, d)`` from
``hamiltonian_of``, so the tolerance series is timed the same way and the two sit on one
axis honestly.  Refinement is then all-energies-at-once: the slab count is set by the
hardest energy in the stack, which is also how Mag(nu)s treats an energy array.

The Hamiltonian is built here, from Mag(nu)s's own V_CC, and handed to both codes -- so
this measures the methods and not the two packages' matter-potential conventions.

    python notebooks/append_npe_rtol_series.py

Run it on an idle machine; it writes a backup beside the file first.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import datetime
import json
import pathlib
import re
import shutil
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, '/home/mbustamante/Research/NuOscProb/NuOscProbExact/src')
import gen_profile_benchmarks as gpb                      # noqa: E402

from magnus import globaldefs as gd, matter               # noqa: E402

HERE = pathlib.Path(__file__).parent
STORE = HERE/'external_profile_benchmarks.json'
N_MAX = 131072          # the slab series needed 32768; exhaustion raises, so leave room
SERIES = 'NuOscProbExact, rtol'


def npe_rtol_points(d, prof, dials, ref):
    import slabs as npe
    fn = getattr(npe, 'probabilities_%dnu_profile' % d, None)
    if fn is None:
        return []                                   # no five-flavour route exists

    E = prof['energies']
    hv = gpb.h_vac(d)
    proj = gpb.matter.matter_potential_projector(d)

    def hamiltonian_of(x):
        """(n_energies, n_positions, d, d) -- the batch 1.14.0 accepts."""
        v = np.asarray(prof['vcc'](np.asarray(x, dtype=float)))
        H = np.broadcast_to((hv[None, None]/E[:, None, None, None]).astype(complex),
                            (len(E), len(v), d, d)).copy()
        H += v[None, :, None, None]*proj[None, None]
        return H

    col = gd.NUMU*d + gd.NUMU
    out = []
    for rtol in dials:
        n_used = {}

        def call(r=rtol, n_used=n_used):
            p, n = fn(hamiltonian_of, prof['baseline'], rtol=r, atol=r*1.0e-2,
                      n_max=N_MAX, return_n_slabs=True)
            n_used['n'] = n
            return np.asarray(p)

        try:
            P = call()
        except ValueError as exc:
            # Not a failure of the run: past some tolerance the round-off accumulated
            # over the slab product overtakes the discretisation error, and no slab count
            # reaches the request.  That floor is a property of the method and is worth
            # recording -- it is the same one the n_slabs series bottoms out on.
            msg = str(exc)
            best = re.search(r'lowest estimate reached was ([0-9.e-]+)', msg)
            print('    rtol %.0e  UNREACHABLE -- floor near %s'
                  % (rtol, best.group(1) if best else '?'), file=sys.stderr, flush=True)
            out.append(dict(label='%.0e' % rtol, rtol=rtol, unreachable=True,
                            best_error_estimate=float(best.group(1)) if best else None,
                            note=msg.split('.  ')[0]))
            continue
        us = 1.0e6*gpb.timed(call)/len(E)
        err = float(np.max(np.abs(P[..., col] - ref)))
        out.append(dict(label='%.0e' % rtol, rtol=rtol, n_slabs_used=int(n_used['n']),
                        us_per_probability=us, max_abs_error=err))
        print('    rtol %.0e  %9.2f us/prob  err %.3e  (%d slabs)'
              % (rtol, us, err, n_used['n']), file=sys.stderr, flush=True)
    return out


def main():
    store = json.loads(STORE.read_text())
    shutil.copy(STORE, STORE.with_suffix('.json.pre_npe_rtol'))
    prof = gpb.exponential_profile()

    for case in store['cases']:
        d = case['flavours']
        if any(s['name'] == SERIES for s in case['series']):
            continue
        print('--- %d flavours ---' % d, file=sys.stderr, flush=True)
        pts = npe_rtol_points(d, prof, gpb.MG_DIALS, np.array(case['reference']))
        if not pts:
            print('    no %dnu route' % d, file=sys.stderr, flush=True)
            continue
        case['series'].append({'name': SERIES, 'dial': 'rtol', 'points': pts})
        STORE.write_text(json.dumps(store, indent=1))

    store['npe_rtol_series'] = {
        'added': datetime.date.today().isoformat(),
        'why': ('so both codes answer the same request; the n_slabs series is kept '
                'beside it, unplotted'),
        'nuoscprobexact': '1.14.0, batched hamiltonian_of, refinement all-energies-at-once',
        'n_max': N_MAX,
    }
    STORE.write_text(json.dumps(store, indent=1))
    print('wrote %s' % STORE.name, file=sys.stderr)


if __name__ == '__main__':
    main()
