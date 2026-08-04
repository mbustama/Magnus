# -*- coding: utf-8 -*-
"""Battery 6 -- flavor count as a first-class axis (REQUIRED by the handover).

Every measurement behind the eleven commits was 2nu or 3nu.  4nu/5nu were touched once,
structurally only (shape, finiteness, unitarity).  Neither was ever scored against
solve_ivp on the new path.  This battery supplies that evidence.

PASS CRITERIA, stated before running:
  (C1) every answer either meets the requested tolerance (1e-3) or raises a warning
       (ToleranceNotAchievedWarning / MagnusConvergenceWarning / HybridCertification).
       A silent miss is a failure regardless of magnitude.
  (C2) unitarity better than 1e-9 everywhere.
  (C3) flagged separately: accuracy degrading MONOTONICALLY with d, which would suggest a
       d-dependent term neglected somewhere.
Oracle: solve_ivp/DOP853, rtol=1e-12 atol=1e-14, via a SINGLE traversal per scan (t_eval).
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import harness6 as H6
import magnus.adiabatic as ad
import magnus.oscprob as oscprob

TOL = 1.0e-3
# 50 MeV over 0.5 R_sun keeps the D31 phase at ~4e4 rad: a genuine integration problem, and
# a DOP853 oracle that is affordable at d = 5.  A single solve serves a whole scan.
E_DEF = 50.0e6
L_MAX = 0.5*H.R_SUN
NE = H.solar_ne()
VCC = H.vcc_of(NE)

WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning')

ROWS = []


def call(family, d, energy, L, params, nsi=None, liv=None, nubar=False, **kw):
    if family == 'std':
        return oscprob.osc_prob_matter_std_potential(
            d, NE, energy, L, params, L0=0.0, nubar=nubar,
            density_is_of_number_of_electrons=True, **kw)
    if family == 'nsi':
        return oscprob.osc_prob_matter_nsi(
            d, NE, energy, L, params, nsi, L0=0.0, nubar=nubar,
            density_is_of_number_of_electrons=True, **kw)
    if family == 'liv':
        return oscprob.osc_prob_liv(
            d, energy, L, params, liv, rho_func=NE, L0=0.0, nubar=nubar,
            density_is_of_number_of_electrons=True, **kw)
    raise ValueError(family)


def score(label, family, d, N, energy=E_DEF, nubar=False, ne_func=None, l_max=None,
          nsi=None, liv=None, params=None, **kw):
    """One configuration: wrapper vs a single-traversal solve_ivp oracle."""
    params = params if params is not None else H.params_for(d)
    ne_func = ne_func if ne_func is not None else NE
    vcc = H.vcc_of(ne_func, nubar=nubar)
    l_max = l_max if l_max is not None else L_MAX
    Ls = np.linspace(0.05*l_max, l_max, N) if N > 1 else np.array([l_max])

    H_of_l = H6.H_family(family, d, energy, vcc, params, nsi=nsi, liv=liv, nubar=nubar)

    t0 = time.time()
    with H.Caught() as c:
        P = np.asarray(call(family, d, energy, Ls if N > 1 else float(Ls[0]), params,
                            nsi=nsi, liv=liv, nubar=nubar, **kw))
    t_pkg = time.time() - t0
    P = P.reshape((N, d, d)) if N > 1 else P.reshape((1, d, d))

    t0 = time.time()
    Us = H.exact_U_many(H_of_l, 0.0, Ls, d)
    t_ref = time.time() - t0
    Pref = np.array([H.P_of(U) for U in Us])

    err = H.maxabs(P - Pref)
    unit = max(H.unitarity(U) for U in Us)
    warned = any(n in WARN_OK for n in c.names)
    silent_miss = bool(err > TOL and not warned)

    ROWS.append(dict(label=label, family=family, d=d, N=N, err=err, warned=c.names,
                     silent_miss=silent_miss, t_pkg=t_pkg, t_ref=t_ref, unit=unit))
    print('%-38s fam=%-3s d=%d N=%-4d err=%9.3e  %-7s warns=%-28s %s'
          % (label, family, d, N, err, '%.2fs' % t_pkg,
             ','.join(c.names) or '-',
             'SILENT MISS' if silent_miss else 'ok'), flush=True)
    return err


# ----------------------------------------------------------------------
def sub1():
    """6.1 Accuracy against solve_ivp, 4nu and 5nu, on the new path -- either side of N=25."""
    print('## 6.1  standard matter, solar profile, N either side of '
          'HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25')
    for d in (2, 3, 4, 5):
        for N in (1, 8, 24, 25, 26, 60):
            score('6.1 std N=%d' % N, 'std', d, N)
    print()


def sub2():
    """6.2 4nu/5nu with NSI, and with LIV (n_liv = 0 and 1)."""
    print('## 6.2  NSI and LIV at every flavor count, either side of the seam')
    for d in (2, 3, 4, 5):
        for N in (8, 60):
            score('6.2 nsi N=%d' % N, 'nsi', d, N, nsi=H6.nsi_params_for(d))
            for nl in (0, 1):
                score('6.2 liv n_liv=%d N=%d' % (nl, N), 'liv', d, N,
                      liv=H6.liv_params_for(d, n_liv=nl))
    print()


def sub3():
    """6.3 Multi-resonance at 4nu and 5nu: per-pair window bookkeeping.

    Confirms no crossing is swallowed by a window belonging to a DIFFERENT pair -- the
    over-merging-across-pairs case the handover calls untested.
    """
    print('## 6.3  multi-resonance, per-pair windows vs a dense per-pair gamma scan')
    l1 = 1.0*H.L_SCALE
    fd = l1*1e-6
    for d in (2, 3, 4, 5):
        p = H.params_for(d, sterile_scale=1.0)
        ne = H.modulated_ne(amp=0.9, n_cycles=6.0, span=l1)
        vcc = H.vcc_of(ne)
        H_of_l = H6.H_family('std', d, E_DEF, vcc, p)
        win, cands = ad.find_nonadiabatic_windows(H_of_l, 0.0, l1, threshold=0.1, n_probe=200)

        ls = np.linspace(0.0, l1, 20000)
        total_over, covered_over = 0, 0
        per_pair = []
        for j in range(d):
            for k in range(j + 1, d):
                g = np.array([ad._point_adiabaticity(H_of_l, float(l), j, k, fd, (0.0, l1))
                              for l in ls])
                over = ls[g > 0.1]
                total_over += over.size
                if over.size:
                    inside = np.zeros(over.size, dtype=bool)
                    for (a, b) in win:
                        inside |= (over >= a) & (over <= b)
                    covered_over += int(inside.sum())
                    per_pair.append('(%d,%d):%d/%d' % (j, k, int(inside.sum()), over.size))
        U, w2, cert = ad.hybrid_propagator(H_of_l, 0.0, l1, rtol=TOL, atol=TOL)
        err = H.maxabs(H.P_of(U) - H.P_of(H.exact_U(H_of_l, 0.0, l1, d)))
        miss = total_over - covered_over
        print('  d=%d pairs=%2d windows=%2d cands=%3d  dense gamma>thr pts %6d, '
              'UNCOVERED %6d  hybrid err=%.3e cert=%s  %s'
              % (d, d*(d - 1)//2, len(win), len(cands), total_over, miss, err, cert,
                 ' '.join(per_pair)), flush=True)
        ROWS.append(dict(label='6.3 multires', family='std', d=d, N=1, err=err,
                         warned=[], silent_miss=bool(cert and err > TOL),
                         t_pkg=0.0, t_ref=0.0, unit=H.unitarity(U)))
    print()


def sub4():
    """6.4 Degenerate and near-degenerate levels: gap ~ 0 over a stretch."""
    print('## 6.4  degenerate / near-degenerate mass splittings (gap ~ 0)')
    for d in (4, 5):
        base = H.params_for(d)
        for name, tweak in [('D41 == D31 exactly', {'D41': base['D31']}),
                            ('D41 = D31*(1+1e-9)', {'D41': base['D31']*(1 + 1e-9)}),
                            ('D41 = D31*(1+1e-14)', {'D41': base['D31']*(1 + 1e-14)}),
                            ('D21 == 0', {'D21': 0.0})]:
            p = dict(base)
            p.update(tweak)
            try:
                score('6.4 %s' % name, 'std', d, 30, params=p)
            except Exception as ex:      # noqa: BLE001 -- an exception IS the finding
                print('  d=%d %-22s RAISED %s: %s' % (d, name, type(ex).__name__, ex),
                      flush=True)
                ROWS.append(dict(label='6.4 %s' % name, family='std', d=d, N=30,
                                 err=float('nan'), warned=['RAISED'], silent_miss=False,
                                 t_pkg=0.0, t_ref=0.0, unit=0.0))
        # and the raw detector on an exactly-degenerate pair, where _point_adiabaticity
        # returns inf
        p = dict(base)
        p['D41'] = base['D31']
        H_of_l = H6.H_family('std', d, E_DEF, VCC, p)
        w, c = ad.find_nonadiabatic_windows(H_of_l, 0.0, 0.3*H.L_SCALE, threshold=0.1)
        gam = [x.get('gamma') for x in c]
        print('  d=%d exact degeneracy: %d candidates, %d windows, gamma inf/nan count = %d'
              % (d, len(c), len(w), sum(1 for g in gam if g is not None
                                        and not np.isfinite(g))), flush=True)
    print()


def sub5():
    """6.5 Cost scaling of the gamma sweep: 1.31x at 2nu -- what at d=4,5 (10x the pairs)?"""
    print('## 6.5  cost of the gamma sweep vs the extrema-only detector, by dimension')
    l1 = 1.0*H.L_SCALE
    ne = H.modulated_ne(amp=0.9, n_cycles=6.0, span=l1)
    vcc = H.vcc_of(ne)
    for d in (2, 3, 4, 5):
        p = H.params_for(d)
        H_of_l = H6.H_family('std', d, E_DEF, vcc, p)
        ts = []
        for _ in range(3):
            t0 = time.time()
            ad.find_nonadiabatic_windows(H_of_l, 0.0, l1, threshold=0.1, n_probe=200)
            ts.append(time.time() - t0)
        tc = []
        for _ in range(3):
            t0 = time.time()
            ad.find_resonance_candidates(H_of_l, 0.0, l1, n_probe=200)
            tc.append(time.time() - t0)
        print('  d=%d pairs=%2d  find_windows %.3f s   find_candidates %.3f s   ratio %.2fx'
              % (d, d*(d - 1)//2, min(ts), min(tc), min(ts)/min(tc)), flush=True)
    print()


def sub6():
    """6.6 Antineutrinos at 4nu and 5nu -- two uncovered axes compounded."""
    print('## 6.6  antineutrinos at every flavor count')
    for d in (2, 3, 4, 5):
        for N in (8, 60):
            score('6.6 nubar std N=%d' % N, 'std', d, N, nubar=True)
        score('6.6 nubar nsi N=60', 'nsi', d, 60, nubar=True, nsi=H6.nsi_params_for(d))
    print()


def sub7():
    """6.7 ip_exp must still decline for d > 2 at every entry point."""
    print('## 6.7  ip_exp still declines for d > 2 (the reason 3/4/5nu were unaffected '
          'by PR #23)')
    seen = {}
    real = oscprob._osc_prob_ip_exp_dispatch

    def spy(*a, **k):
        r = real(*a, **k)
        seen.setdefault(key[0], []).append(r is not NotImplemented)
        return r

    oscprob._osc_prob_ip_exp_dispatch = spy
    try:
        for d in (2, 3, 4, 5):
            for fam in ('std', 'nsi', 'liv'):
                key = ['%s d=%d' % (fam, d)]
                nsi = H6.nsi_params_for(d) if fam == 'nsi' else None
                liv = H6.liv_params_for(d, n_liv=0) if fam == 'liv' else None
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    # strategy='magnus' so hybrid stands aside and ip_exp is actually reached
                    call(fam, d, E_DEF, 0.2*H.R_SUN, H.params_for(d), nsi=nsi, liv=liv,
                         strategy='magnus')
    finally:
        oscprob._osc_prob_ip_exp_dispatch = real
    for k in sorted(seen):
        answered = any(seen[k])
        print('  %-12s reached %d time(s), ANSWERED=%s %s'
              % (k, len(seen[k]), answered,
                 '' if (k.endswith('d=2') or not answered) else '<-- UNEXPECTED for d>2'))
    print()


def summary():
    print('\n=== BATTERY 6 SUMMARY ===')
    scored = [r for r in ROWS if np.isfinite(r['err'])]
    print('configurations scored      : %d' % len(scored))
    print('worst error                : %.3e' % max(r['err'] for r in scored))
    print('worst unitarity            : %.3e  (criterion: < 1e-9)'
          % max(r['unit'] for r in scored))
    sm = [r for r in scored if r['silent_miss']]
    print('SILENT MISSES (>1e-3, no warning): %d   <-- criterion (C1)' % len(sm))
    for r in sorted(sm, key=lambda r: -r['err']):
        print('   %-38s fam=%-3s d=%d N=%-4d err=%.3e' % (r['label'], r['family'], r['d'],
                                                          r['N'], r['err']))
    print('\n(C3) worst error by dimension, standard matter:')
    for d in (2, 3, 4, 5):
        e = [r['err'] for r in scored if r['d'] == d and r['family'] == 'std']
        if e:
            print('   d=%d  worst %.3e  median %.3e' % (d, max(e), float(np.median(e))))
    np.save('battery6_rows.npy', np.array(ROWS, dtype=object), allow_pickle=True)


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3', '4', '5', '6', '7']
    print('# Battery 6 -- flavor count.  E=%.0f MeV, L_max=%.3f R_sun, tol=%.0e\n'
          % (E_DEF/1e6, L_MAX/H.R_SUN, TOL))
    for w in which:
        {'1': sub1, '2': sub2, '3': sub3, '4': sub4, '5': sub5, '6': sub6,
         '7': sub7}[w]()
    summary()
