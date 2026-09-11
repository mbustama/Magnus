# -*- coding: utf-8 -*-
r"""Freezes the ``solve_ivp`` references for Figure 12's front-width scan.

``make_shock_reference.py`` freezes five widths for notebook 14.  This freezes twenty,
log-spaced across the same four decades, because five points cannot resolve the shape the
figure turned out to have: the closed form's cost peaks somewhere between $7$ and $70$~km
and five samples locate that only to within a decade.

**Why a separate store.**  The notebook's file keys its cases ``'%.0e' % width``, which is
unique for five widths a decade apart and collides for twenty a factor of $1.6$ apart --
$2.6 \cdot 10^{-6}$ and $4.3 \cdot 10^{-6}$ would both be ``'3e-06'`` and ``'4e-06'`` is
two different widths away.  Rather than re-key a file three other figures read, this writes
its own, at full precision.  The five already-frozen cases are copied in rather than
recomputed, so the scan costs seventeen solves and not twenty.

**Resumable.**  The store is written after every width, and a width already present is
skipped.  Interrupting this costs at most the width in flight, about three and a half
minutes.  Re-running it resumes.

    python notebooks/make_shock_scan_references.py

Each width costs two solves, about $105$~s at rtol $10^{-14}$ plus $85$~s at $10^{-13}$ for
the self-convergence bound, so a cold run is a little over an hour.
"""

import json
import pathlib
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'src'))

import gen_shock_cost as G                                   # noqa: E402

OUT = HERE/'shock_reference_scan.json'
OLD = HERE/'shock_reference.json'
N_FINGERPRINT = 40


def hexed(a):
    a = np.asarray(a, dtype=float)
    return [float(x).hex() for x in a.ravel()], list(a.shape)


def freeze(width, ns, L0, Ls, fingerprint_l):
    """One width: the reference, its own error, and the two fingerprints that guard it."""
    ne = ns['sn_shock_ne'](width)
    H = ns['make_H'](ne)

    def rhs(l, y, H=H):
        return (-1j*np.asarray(H(l)) @ y.reshape(3, 3)).ravel()

    def probabilities(rtol, atol):
        sol = solve_ivp(rhs, (float(L0), float(Ls[-1])), np.eye(3, dtype=complex).ravel(),
                        rtol=rtol, atol=atol, method='DOP853', t_eval=Ls)
        if not sol.success:
            raise SystemExit('solve_ivp failed at w=%.6e, rtol=%.0e: %s'
                             % (width, rtol, sol.message))
        U = np.array([sol.y[:, i].reshape(3, 3) for i in range(len(Ls))])
        return np.swapaxes(U.real**2 + U.imag**2, -1, -2)

    # The same settings, and for the same reasons, as make_shock_reference.py: SciPy
    # clamps rtol at 16*eps, so 1e-13 is the closest genuinely different solve and its
    # movement is a conservative bound on this one's error rather than an optimistic one.
    P = probabilities(1.0e-14, 1.0e-16)
    self_conv = float(np.max(np.abs(probabilities(1.0e-13, 1.0e-15) - P)))
    flat, shape = hexed(P)
    Hf = np.asarray(H(fingerprint_l), dtype=complex)
    return {'P': flat, 'shape': shape, 'self_convergence': self_conv,
            'fingerprint_ne': hexed(np.asarray(ne(fingerprint_l), dtype=float))[0],
            'fingerprint_h': hexed(np.concatenate([Hf.real.ravel(), Hf.imag.ravel()]))[0]}


def main():
    ns, L0, Ls = G.NS, G.L0, G.Ls
    fingerprint_l = np.linspace(float(L0), float(Ls[-1]), N_FINGERPRINT)

    if OUT.exists():
        store = json.loads(OUT.read_text())
    else:
        store = {'note': ('solve_ivp DOP853 rtol=1e-14 atol=1e-16, one case per front '
                          'width of Figure 12; keys are "%.6e" of the width, which twenty '
                          'log-spaced widths need and "%.0e" does not give; produced by '
                          'notebooks/make_shock_scan_references.py, do not edit by hand'),
                 'Ls': hexed(Ls)[0], 'fingerprint_l': hexed(fingerprint_l)[0], 'cases': {}}
        # The five already frozen for notebook 14 are the same physics at the same
        # baselines, so they are copied across rather than recomputed.  Their fingerprints
        # come with them and are checked on every read, so a stale one cannot ride in.
        if OLD.exists():
            old = json.loads(OLD.read_text())
            if old['Ls'] == store['Ls'] and old['fingerprint_l'] == store['fingerprint_l']:
                for key, case in old['cases'].items():
                    store['cases']['%.6e' % float(key)] = case
                print('  carried over %d already-frozen widths' % len(old['cases']),
                      flush=True)
            else:
                print('  notebook 14 store is on different baselines; not carried over',
                      flush=True)
        OUT.write_text(json.dumps(store, indent=1))

    todo = [w for w in G.WIDTHS if '%.6e' % w not in store['cases']]
    print('  %d of %d widths to freeze' % (len(todo), len(G.WIDTHS)), flush=True)
    for i, width in enumerate(todo):
        t0 = time.perf_counter()
        case = freeze(width, ns, L0, Ls, fingerprint_l)
        store['cases']['%.6e' % width] = case
        OUT.write_text(json.dumps(store, indent=1))          # checkpoint: every width
        print('  [%2d/%2d] w=%.4e (%7.2f km)  self-convergence %.2e  %5.0f s'
              % (i + 1, len(todo), width, width*G.RAY_KM, case['self_convergence'],
                 time.perf_counter() - t0), flush=True)
    print('wrote %s (%d cases)' % (OUT.name, len(store['cases'])), flush=True)


if __name__ == '__main__':
    main()
