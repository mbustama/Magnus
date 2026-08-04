# -*- coding: utf-8 -*-
"""False-positive and true-positive rates for every warning the package can raise.

Robustness programme, item 5.  The brief's bar for a warning is four-part -- what was
detected, what it means for the answer, what to change, and when it is genuinely safe to
ignore -- and the second part cannot be written without knowing how often the warning fires
on an answer that was fine.  ``DECISION_DISPATCH_ORDER.md`` §5 records
``MagnusConvergenceWarning`` firing on rows accurate to 1.6e-06; this measures that properly,
over the profile families the package actually serves.

Definitions used throughout, stated before running:

  * **true positive**  -- the warning fired AND the answer is outside the requested tolerance.
  * **false positive** -- the warning fired AND the answer is inside it.
  * **silent miss**    -- no warning fired AND the answer is outside it.  This is the failure
    that matters; a false positive is only noise.

``MagnusConvergenceWarning`` additionally gets a mechanism measurement: every call to
``magnus._warn_slab_norm`` is recorded, so that "some refinement level had a slab too wide"
can be separated from "the level whose answer was returned did".  If the first dominates, the
warning is reporting on results nobody receives, and the fix is mechanical rather than
editorial.

**Cost, and how it was brought down.**  The first version of this script was killed unfinished
at 100 of 180 configurations: it put ``solve_ivp`` on every case, including N = 40 at 3nu and
10 MeV, which is the exact trap the brief warns about (25 s for 2 of 80 baselines).  This one
splits the population by oracle -- ``expm`` composed across segments for piecewise-constant
profiles, which is **exact** and costs nothing, and ``solve_ivp`` only for the smooth families,
held to >= 30 MeV and N <= 8.  Same question, roughly a tenth of the time.

Run:  python warn_fp.py [n_random]
"""

import sys
import time
import warnings

import numpy as np

from scipy.linalg import expm

import harness as H
import magnus.magnus as magnuscore
import magnus.oscprob as oscprob
from battery2 import bump_profile, ne_res_for
from battery3 import noisy_ne

TOL = 1e-3
L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0


class NormSpy:
    """Every ``||Omega||_2`` the convergence check saw, in the order it saw them.

    The last entry belongs to the refinement level whose answer is returned, so
    ``last >= pi`` is "the returned grid was too coarse" while ``any >= pi`` is what the
    warning currently reports.
    """

    def __enter__(self):
        self.vals = []
        self._orig = magnuscore._warn_slab_norm

        def spy(nmax):
            self.vals.append(float(nmax))
            return self._orig(nmax)

        magnuscore._warn_slab_norm = spy
        return self

    def __exit__(self, *e):
        magnuscore._warn_slab_norm = self._orig
        return False


def piecewise_profile(edges, values):
    def ne(l):
        x = np.asarray(l, dtype=float)
        idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(values) - 1)
        return H.scalarize(values[idx])
    return ne


def exact_piecewise_P(H_func, edges, Ls, d):
    """EXACT probabilities: H is constant on each segment, so the exponentials compose."""
    stops = sorted(set(np.concatenate([edges, np.atleast_1d(Ls)]).tolist()))
    U, cursor, at = np.eye(d, dtype=complex), float(stops[0]), {}
    for nxt in stops:
        if nxt > cursor:
            Hm = np.asarray(H_func(0.5*(cursor + nxt)), dtype=complex)
            U = expm(-1j*Hm*(nxt - cursor)) @ U
            cursor = nxt
        at[nxt] = U.copy()
    return np.array([H.P_of(at[float(L)]) for L in np.atleast_1d(Ls)])


def build_population(n_random=8):
    """(label, ne, d, energy, N, extra kwargs).  Families the package actually serves."""
    rng = np.random.default_rng(20260804)
    p_by_d = {d: H.params_for(d) for d in (2, 3)}
    ner2 = ne_res_for(2, p_by_d[2], 10.0e6)

    def castle(n_walls=12):
        lo, hi = 0.02*H.NE0, 0.30*H.NE0
        l_ini, l_fin = 0.1*L1, 0.9*L1
        edges = np.linspace(l_ini, l_fin, n_walls + 1)

        def ne(l):
            x = np.asarray(l, dtype=float)
            u = (x - l_ini)/(l_fin - l_ini)
            idx = np.floor(np.clip(u, 0.0, 1.0 - 1e-15)*n_walls)
            return H.scalarize(np.where((x < l_ini) | (x > l_fin), lo,
                                        np.where(idx % 2 == 0, lo, hi)))
        return ne, edges

    cw, cw_edges = castle()
    families = [
        ('solar exponential', H.solar_ne(), {}),
        ('multi-resonance', H.modulated_ne(amp=0.9, n_cycles=6.0, span=SPAN), {}),
        ('noisy', noisy_ne(), {}),
        ('sinusoid, span/7', H.sine_ne(SPAN/7.0, base_ratio=3.0e-2), {}),
        ('resolvable bump w=1e-2', bump_profile(ner2, 0.45*SPAN, 1e-2*SPAN), {}),
        ('castle wall, edges declared', cw, {'t_breakpoints': cw_edges}),
        ('castle wall, edges NOT declared', cw, {}),
    ]
    for i in range(n_random):
        families.append(('random Fourier #%d' % i,
                         H.fourier_ne(rng, n_modes=int(rng.integers(2, 9)), span=SPAN,
                                      base_ratio=float(10.0**rng.uniform(-2.5, -1.0)),
                                      amp=float(rng.uniform(0.2, 0.9))), {}))

    cases = []
    for label, ne, kw in families:
        for d in (2, 3):
            # >= 30 MeV and N <= 8: solve_ivp is what made the first version of this
            # unfinishable, and it is cheapest at high energy and few baselines.
            for energy in (30.0e6, 100.0e6):
                for N in (1, 8):
                    cases.append((label, ne, d, p_by_d[d], energy, N, kw, None))

    # Piecewise arm: exact oracle, so d = 4 and 5 and larger scans are affordable here.
    rng2 = np.random.default_rng(4242)
    for i in range(2*n_random):
        d = int(rng2.choice([2, 3, 4, 5]))
        n_seg = int(rng2.integers(2, 9))
        cuts = np.sort(rng2.uniform(0.0, 1.0, n_seg - 1))*L1
        edges = np.concatenate([[0.0], cuts, [L1]])
        values = H.NE0*10.0**rng2.uniform(-2.5, -0.5, n_seg)
        ne = piecewise_profile(edges, values)
        params = H.params_for(d)
        for N in (1, 8, 30):
            cases.append(('piecewise #%d' % i, ne, d, params,
                          float(10.0**rng2.uniform(7.3, 8.0)), N, {}, edges))
    return cases


def score(label, ne, d, params, energy, N, kw, edges=None):
    Ls = np.linspace(0.05*L1, L1, N) if N > 1 else np.array([L1])
    H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
    try:
        Pref = (exact_piecewise_P(H_of_l, edges, Ls, d) if edges is not None
                else np.array([H.P_of(U) for U in H.exact_U_many(H_of_l, L0, Ls, d)]))
    except Exception as exc:                       # noqa: BLE001
        return dict(label=label, skipped=type(exc).__name__)

    with NormSpy() as spy, warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            P = np.asarray(oscprob.osc_prob_matter_std_potential(
                d, ne, energy, Ls if N > 1 else float(Ls[0]), params, L0=L0,
                density_is_of_number_of_electrons=True, **kw)).reshape(N, d, d)
        except Exception as exc:                   # noqa: BLE001
            return dict(label=label, skipped=type(exc).__name__)

    names = sorted({w.category.__name__ for w in caught})
    err = H.maxabs(P - Pref)
    return dict(label=label, d=d, energy=energy, N=N, err=err, warns=names,
                outside=bool(err > TOL),
                norm_any=bool(spy.vals and max(spy.vals) >= np.pi),
                norm_last=bool(spy.vals and spy.vals[-1] >= np.pi),
                n_norm_calls=len(spy.vals))


def main():
    n_random = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    cases = build_population(n_random)
    print('# Warning false-positive measurement.  %d configurations, tolerance %.0e' %
          (len(cases), TOL))
    print('# families x d in {2,3} x E in {10, 50} MeV x N in {1, 8, 40}; L1 = 1 solar '
          'scale height\n')
    rows, t0 = [], time.time()
    for i, c in enumerate(cases):
        r = score(*c)
        rows.append(r)
        if (i + 1) % 20 == 0:
            print('   ... %d/%d  (%.0f s)' % (i + 1, len(cases), time.time() - t0),
                  flush=True)
    scored = [r for r in rows if 'skipped' not in r]
    print('\nscored %d of %d (%d skipped: %s)\n'
          % (len(scored), len(rows), len(rows) - len(scored),
             ', '.join(sorted({r['skipped'] for r in rows if 'skipped' in r})) or '-'))

    classes = sorted({n for r in scored for n in r['warns']})
    print('%-32s %6s %6s %6s %8s' % ('warning', 'fired', 'TP', 'FP', 'FP rate'))
    for cls in classes:
        fired = [r for r in scored if cls in r['warns']]
        tp = [r for r in fired if r['outside']]
        fp = [r for r in fired if not r['outside']]
        print('%-32s %6d %6d %6d %7.0f%%'
              % (cls, len(fired), len(tp), len(fp), 100.0*len(fp)/max(len(fired), 1)))

    silent = [r for r in scored if r['outside'] and not r['warns']]
    print('\noutside tolerance, ANY warning : %d'
          % len([r for r in scored if r['outside'] and r['warns']]))
    print('outside tolerance, SILENT      : %d' % len(silent))
    for r in sorted(silent, key=lambda r: -r['err'])[:12]:
        print('    %-32s d=%d E=%4.0fMeV N=%-3d err=%.3e'
              % (r['label'], r['d'], r['energy']/1e6, r['N'], r['err']))

    print('\n--- MagnusConvergenceWarning mechanism (single points only) ---')
    sp = [r for r in scored if r['N'] == 1 and r['n_norm_calls'] > 0]
    any_ = [r for r in sp if r['norm_any']]
    last = [r for r in sp if r['norm_last']]
    only_intermediate = [r for r in sp if r['norm_any'] and not r['norm_last']]
    print('single-point cases with a norm check : %d' % len(sp))
    print('  some refinement level had ||Om||>=pi: %d   (what the warning reports today)'
          % len(any_))
    print('  the RETURNED level had ||Om||>=pi   : %d' % len(last))
    print('  fired only on an intermediate level : %d' % len(only_intermediate))
    for tag, pop in (('any level', any_), ('returned level', last)):
        fp = [r for r in pop if not r['outside']]
        print('  FP rate if keyed on %-15s: %3d/%3d = %.0f%%'
              % (tag, len(fp), len(pop), 100.0*len(fp)/max(len(pop), 1)))
    np.save('warn_fp_rows.npy', np.array(rows, dtype=object), allow_pickle=True)
    print('\nrows saved to warn_fp_rows.npy')


if __name__ == '__main__':
    main()
