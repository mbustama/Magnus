# -*- coding: utf-8 -*-
"""Battery 8 -- fuzzing PIECEWISE (discontinuous) profiles, and the quadrature methods.

The Battery 5.4 fuzzer drew smooth profiles only (random Fourier sums), so it could not reach
the class the worst defect lived in: an unmarked density discontinuity, which returned a
probability wrong by 0.54 while reporting certified=True.  This battery closes that gap.

ORACLE.  Profiles here are piecewise CONSTANT, which makes `scipy.linalg.expm` composed across
the segments the **exact** evolution operator -- not an approximation, and immune to the
"adaptive solver steps over the feature" failure that makes solve_ivp untrustworthy on a jump.
Every error below is measured against that.

PASS CRITERIA, stated before running:
  (G1) every answer meets the requested tolerance OR warns.  A silent miss is a failure.
  (G2) declaring the edges via t_breakpoints must be at least as accurate as not declaring
       them, and must never warn where the undeclared version did not.
  (G3) trapezoid/simpson: no silent miss (they are expected to be *inaccurate* on hard
       profiles -- see DECISION_CUMULATIVE_DEFAULT.md 4f -- but not silently so).
"""

import sys
import warnings

import numpy as np
from scipy.linalg import expm

import harness as H
import harness6 as H6
import magnus.oscprob as oscprob

TOL = 1e-3
WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning')


def piecewise_profile(edges, values):
    """Piecewise-constant electron density: `values[i]` on [edges[i], edges[i+1])."""
    edges = np.asarray(edges, dtype=float)
    values = np.asarray(values, dtype=float)

    def ne(l):
        x = np.asarray(l, dtype=float)
        idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(values) - 1)
        return H.scalarize(values[idx])
    return ne


def exact_piecewise_U(H_func, edges, Ls):
    """EXACT U at each baseline: H is constant on each segment, so expm composes exactly."""
    edges = np.asarray(edges, dtype=float)
    d = np.asarray(H_func(0.5*(edges[0] + edges[1]))).shape[-1]
    out = []
    for L in np.atleast_1d(np.asarray(Ls, dtype=float)):
        U = np.eye(d, dtype=complex)
        cuts = np.unique(np.concatenate([edges[(edges > 0.0) & (edges < L)], [0.0, L]]))
        for a, b in zip(cuts[:-1], cuts[1:]):
            Hm = np.asarray(H_func(0.5*(a + b)), dtype=complex)
            U = expm(-1j*Hm*(b - a)) @ U
        out.append(U)
    return np.array(out)


def fuzz(n_cases=150, seed=20260804):
    print('## 8.1  piecewise-constant profile fuzzer, %d cases, d from {2,3,4,5}' % n_cases)
    print('##      oracle: expm composed across segments (EXACT for these profiles)\n')
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_cases):
        d = int(rng.choice([2, 3, 4, 5]))
        p = H.params_for(d, sterile_scale=float(rng.uniform(0.3, 1.5)))
        E = float(10.0**rng.uniform(6.7, 8.3))
        span = float(rng.uniform(0.3, 1.5))*H.L_SCALE
        n_seg = int(rng.integers(2, 13))
        cuts = np.sort(rng.uniform(0.0, 1.0, n_seg - 1))*span
        edges = np.concatenate([[0.0], cuts, [span]])
        values = H.NE0*10.0**rng.uniform(-2.5, -0.3, n_seg)
        ne = piecewise_profile(edges, values)
        nubar = bool(rng.random() < 0.3)
        N = int(rng.choice([1, 3, 12, 30, 80]))
        Ls = np.linspace(0.05*span, span, N) if N > 1 else np.array([span])

        H_of_l = H6.H_family('std', d, E, H.vcc_of(ne, nubar=nubar), p, nubar=nubar)
        Pref = np.array([H.P_of(U) for U in exact_piecewise_U(H_of_l, edges, Ls)])

        res = {}
        for tag, kw in (('bare', {}), ('breakpoints', {'t_breakpoints': edges[1:-1]})):
            try:
                with H.Caught() as c:
                    P = np.asarray(oscprob.osc_prob_matter_std_potential(
                        d, ne, E, Ls if N > 1 else float(Ls[0]), p, L0=0.0, nubar=nubar,
                        density_is_of_number_of_electrons=True, **kw)).reshape(N, d, d)
                err = H.maxabs(P - Pref)
                res[tag] = (err, any(x in WARN_OK for x in c.names), c.names)
            except Exception as ex:                # noqa: BLE001
                print('   case %3d (%s) RAISED %s: %s'
                      % (i, tag, type(ex).__name__, str(ex)[:80]))
                res[tag] = (float('nan'), True, ['RAISED'])

        rows.append(dict(i=i, d=d, N=N, E=E, n_seg=n_seg, nubar=nubar,
                         bare=res['bare'], bp=res['breakpoints']))
        if (i + 1) % 25 == 0:
            print('   ... %d/%d' % (i + 1, n_cases), flush=True)

    ok = [r for r in rows if np.isfinite(r['bare'][0]) and np.isfinite(r['bp'][0])]
    for tag in ('bare', 'bp'):
        e = np.array([r[tag][0] for r in ok])
        silent = [r for r in ok if r[tag][0] > TOL and not r[tag][1]]
        name = 'edges NOT declared' if tag == 'bare' else 'edges declared (t_breakpoints)'
        print('\n  %-32s median %.2e  p90 %.2e  max %.2e' % (name, np.median(e),
                                                             np.percentile(e, 90), e.max()))
        print('  %-32s outside 1e-3: %d   SILENT: %d   <-- criterion (G1)'
              % ('', sum(x > TOL for x in e), len(silent)))
        for r in sorted(silent, key=lambda r: -r[tag][0])[:12]:
            print('     case %3d d=%d N=%-3d segs=%-3d E=%6.1f MeV err=%.3e'
                  % (r['i'], r['d'], r['N'], r['n_seg'], r['E']/1e6, r[tag][0]))

    worse = [r for r in ok if r['bp'][0] > 10.0*r['bare'][0] and r['bp'][0] > TOL]
    print('\n  (G2) declaring the edges is >10x WORSE in %d / %d cases' % (len(worse), len(ok)))
    for r in worse[:6]:
        print('     case %3d bare %.3e -> breakpoints %.3e' % (r['i'], r['bare'][0], r['bp'][0]))

    by_N = {}
    for r in ok:
        by_N.setdefault(r['N'], []).append(r)
    print('\n  by scan size (edges NOT declared -- the adversarial case):')
    for N in sorted(by_N):
        g = by_N[N]
        s = [r for r in g if r['bare'][0] > TOL and not r['bare'][1]]
        print('    N=%-3d %3d cases, %2d silent (%4.1f%%), median err %.2e'
              % (N, len(g), len(s), 100*len(s)/len(g),
                 float(np.median([r['bare'][0] for r in g]))))
    np.save('battery8_rows.npy', np.array(rows, dtype=object), allow_pickle=True)


def quadrature():
    """8.2 trapezoid / simpson, left where DECISION_CUMULATIVE_DEFAULT.md 4f left them."""
    print('\n## 8.2  quadrature methods (criterion G3: inaccurate is allowed, SILENT is not)')
    p = H.params_for(2)
    ne = H.solar_ne()
    # Deliberately modest: trapezoid/simpson exhaust max_n_slabs AND max_n_tpts_per_slab on a
    # full solar radius (DECISION_CUMULATIVE_DEFAULT.md 4f), so a wide grid here costs tens of
    # minutes and measures only that already-recorded fact.  What is under test is whether they
    # stay HONEST, which a short scan shows just as well.
    for L_frac, E in ((0.05, 50e6), (0.15, 50e6)):
        LM = L_frac*H.R_SUN
        Ls = np.linspace(0.05*LM, LM, 12)
        H_of_l = H6.H_family('std', 2, E, H.vcc_of(ne), p)
        Pref = np.array([H.P_of(U) for U in H.exact_U_many(H_of_l, 0.0, Ls, 2)])
        for method in ('gl', 'trapezoid', 'simpson'):
            for strat in ('auto', 'magnus'):
                with H.Caught() as c:
                    P = np.asarray(oscprob.osc_prob_matter_std_potential(
                        2, ne, E, Ls, p, L0=0.0, density_is_of_number_of_electrons=True,
                        integration_method=method, strategy=strat)).reshape(12, 2, 2)
                err = H.maxabs(P - Pref)
                silent = err > TOL and not any(x in WARN_OK for x in c.names)
                print('  L=%.2f R_sun E=%5.1f MeV %-10s %-7s err=%.3e  %-11s %s'
                      % (L_frac, E/1e6, method, strat, err,
                         'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                         ','.join(c.names) or '-'), flush=True)
    print()


if __name__ == '__main__':
    warnings.simplefilter('always')
    print('# Battery 8 -- piecewise profiles and quadrature.  tol=%.0e\n' % TOL)
    args = sys.argv[1:]
    if not args or '1' in args:
        fuzz(int(args[1]) if len(args) > 1 and args[1].isdigit() else 150)
    if not args or '2' in args:
        quadrature()
