# -*- coding: utf-8 -*-
"""Battery 5 -- designed to break, not to confirm.

The goal is a reproducible wrong answer.

PASS CRITERIA, stated before running:
  (E1) every answer meets the requested tolerance OR warns.  A silent miss is a failure
       regardless of magnitude.  This is the criterion for the fuzzer.
  (E2) determinism: identical calls bit-identical; shuffled input order identical after
       unshuffling; n_jobs>1 == n_jobs=1.
  (E3) adversarial n_slabs / n_tpts_per_slab must not produce a silently wrong answer.
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import harness6 as H6
import magnus.adiabatic as ad
import magnus.oscprob as oscprob
from battery2 import bump_profile, ne_res_for

TOL = 1e-3
WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning')
E_DEF = 50.0e6
L1 = 1.0*H.L_SCALE


def gamma_max(H_of_l, l0, l1, d, n=20000):
    fd = (l1 - l0)*1e-6
    ls = np.linspace(l0, l1, n)
    g = 0.0
    for j in range(d):
        for k in range(j + 1, d):
            v = np.array([ad._point_adiabaticity(H_of_l, float(l), j, k, fd, (l0, l1))
                          for l in ls])
            g = max(g, float(np.nanmax(v[np.isfinite(v)])) if np.isfinite(v).any() else 0.0)
    return g


def sub1():
    """5.1 A profile that defeats the gamma sweep: gamma just UNDER threshold, so nothing
    opens, but the accumulated non-adiabaticity is large.

    Commit 8 raised the bar (gamma is now swept along the path, not only at gap extrema)
    without removing the mechanism: a SUB-threshold accumulation is still invisible."""
    print('## 5.1  gamma held just below the 0.1 threshold over a long stretch')
    print('%10s %11s %11s %8s %11s %s'
          % ('width/span', 'gamma_max', 'hybrid err', 'windows', 'magnus err', 'verdict'))
    p2 = H.params_for(2)
    ner = ne_res_for(2, p2, E_DEF)
    span = L1
    lc = 0.5*span
    worst = None
    for wf in [3e-1, 1e-1, 6e-2, 4e-2, 3e-2, 2e-2, 1.5e-2, 1e-2, 7e-3, 5e-3]:
        ne = bump_profile(ner, lc, wf*span)
        H_of_l = H6.H_family('std', 2, E_DEF, H.vcc_of(ne), p2)
        gm = gamma_max(H_of_l, 0.0, span, 2, n=8000)
        U, win, cert = ad.hybrid_propagator(H_of_l, 0.0, span, rtol=TOL, atol=TOL)
        Pref = H.P_of(H.exact_U(H_of_l, 0.0, span, 2))
        err = H.maxabs(H.P_of(U) - Pref)
        with H.Caught():
            Pm = np.asarray(oscprob.osc_prob_matter_std_potential(
                2, ne, E_DEF, span, p2, L0=0.0, density_is_of_number_of_electrons=True,
                strategy='magnus'))
        em = H.maxabs(Pm - Pref)
        v = ('SUB-THRESHOLD FAILURE' if (gm < 0.1 and err > TOL and cert)
             else ('ok' if err <= TOL else 'wrong (window opened)'))
        if gm < 0.1 and err > TOL and (worst is None or err > worst[1]):
            worst = (wf, err, gm)
        print('%10.0e %11.3e %11.3e %8d %11.3e %s'
              % (wf, gm, err, len(win), em, v), flush=True)
    if worst:
        print('  WORST sub-threshold failure: width=%.0e span, gamma_max=%.3e (<0.1), '
              'err=%.3e = %.0fx the requested %.0e'
              % (worst[0], worst[2], worst[1], worst[1]/TOL, TOL))
    print()


def sub2():
    """5.2 A profile where the cumulative probe stops early -- the "frozen grid" mode of
    NOTES_ADAPTIVE_REFINEMENT §4b, which survives strictness because the agreeing levels are
    BIT-IDENTICAL, not merely close.  A user-supplied t_breakpoints set that dominates the
    grid should reproduce it -- and it would misplace a whole scan, not one point."""
    print('## 5.2  frozen-grid mode: bit-identical refinement levels defeat strictness')
    # Piecewise-constant profile, exactly constant between breakpoints: Magnus is exact on
    # every interval, so every refinement level returns THE SAME number.
    n_walls = 30
    l_ini, l_fin = 0.1*L1, 0.9*L1
    edges = np.linspace(l_ini, l_fin, n_walls + 1)
    lo, hi = 0.02*H.NE0, 0.30*H.NE0

    def ne(l):
        x = np.asarray(l, dtype=float)
        u = (x - l_ini)/(l_fin - l_ini)
        idx = np.floor(np.clip(u, 0.0, 1.0 - 1e-15)*n_walls)
        return H.scalarize(np.where((x < l_ini) | (x > l_fin), lo,
                                    np.where(idx % 2 == 0, lo, hi)))

    p2 = H.params_for(2)
    H_of_l = H6.H_family('std', 2, E_DEF, H.vcc_of(ne), p2)
    Ls = np.linspace(0.05*L1, L1, 60)
    Pref = np.array([H.P_of(U) for U in H.exact_U_many(H_of_l, 0.0, Ls, 2)])

    for name, bps in (('breakpoints INTERIOR only (l_ini/l_fin missing)', edges[1:-1]),
                      ('breakpoints complete (incl. l_ini, l_fin)', edges),
                      ('no breakpoints at all', None)):
        kw = {} if bps is None else {'t_breakpoints': bps}
        with H.Caught() as c:
            P = np.asarray(oscprob.osc_prob_matter_std_potential(
                2, ne, E_DEF, Ls, p2, L0=0.0, density_is_of_number_of_electrons=True,
                **kw)).reshape(60, 2, 2)
        err = H.maxabs(P - Pref)
        silent = err > TOL and not any(x in WARN_OK for x in c.names)
        print('  %-46s err=%.3e  %-11s warns=%s'
              % (name, err, 'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                 ','.join(c.names) or '-'), flush=True)

    # And the probe itself: does the strict ladder freeze on this profile?
    probe = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        oscprob.osc_prob(H_of_l, 0.0, float(Ls[-1]), rtol=TOL, atol=TOL,
                         strict_convergence=True, convergence_info=probe,
                         t_breakpoints=edges[1:-1])
    print('  strict probe on the incomplete-breakpoint profile stopped at n_slabs=%s'
          % probe.get('n_slabs'))
    print()


def sub3():
    """5.3 Discontinuous and pathological profiles."""
    print('## 5.3  discontinuous / pathological profiles (no breakpoints supplied)')
    p2 = H.params_for(2)
    lo, hi = 0.02*H.NE0, 0.30*H.NE0
    mid = 0.5*L1

    def step(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(np.where(x < mid, lo, hi))

    def kink(l):                      # C0 but not C1: |x - mid|
        x = np.asarray(l, dtype=float)
        return H.scalarize(lo + (hi - lo)*np.abs(x - mid)/L1)

    def near_sing(l):                 # approached but not reached
        x = np.asarray(l, dtype=float)
        return H.scalarize(lo + (hi - lo)*0.001/(np.abs(x - 1.02*L1)/L1 + 1e-3))

    def sawtooth(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(lo + (hi - lo)*((x/(0.07*L1)) % 1.0))

    for name, ne in (('step function, edge unmarked', step),
                     ('kink, no jump (C0 not C1)', kink),
                     ('singularity approached, not reached', near_sing),
                     ('sawtooth (jump every 0.07 span)', sawtooth)):
        H_of_l = H6.H_family('std', 2, E_DEF, H.vcc_of(ne), p2)
        Pref = H.P_of(H.exact_U(H_of_l, 0.0, L1, 2))
        for strat in ('auto', 'magnus'):
            with H.Caught() as c:
                P = np.asarray(oscprob.osc_prob_matter_std_potential(
                    2, ne, E_DEF, L1, p2, L0=0.0,
                    density_is_of_number_of_electrons=True, strategy=strat))
            err = H.maxabs(P - Pref)
            silent = err > TOL and not any(x in WARN_OK for x in c.names)
            print('  %-38s %-7s err=%.3e  %-11s %s'
                  % (name, strat, err,
                     'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                     ','.join(c.names) or '-'), flush=True)
    print()


def sub4(n_cases=200):
    """5.4 Random-profile fuzzing -- the highest-yield sub-test in the battery.

    Random smooth profiles (random Fourier sums, controlled bandwidth), random flavour count
    drawn from {2,3,4,5} with random sterile mixings, run through every entry point and
    scored against solve_ivp.  Reports the empirical error distribution and EVERY case
    outside its requested tolerance THAT DOES NOT WARN."""
    print('## 5.4  random-profile fuzzing, %d cases, d drawn from {2,3,4,5}' % n_cases)
    rng = np.random.default_rng(20260803)
    rows = []
    t0 = time.time()
    for i in range(n_cases):
        d = int(rng.choice([2, 3, 4, 5]))
        sterile = float(rng.uniform(0.3, 1.5))
        p = H.params_for(d, sterile_scale=sterile)
        E = float(10.0**rng.uniform(6.7, 8.3))          # 5 - 200 MeV
        span = float(rng.uniform(0.3, 1.5))*H.L_SCALE
        n_modes = int(rng.integers(2, 10))
        ne = H.fourier_ne(rng, n_modes=n_modes, span=span,
                          base_ratio=float(10.0**rng.uniform(-2.5, -0.7)),
                          amp=float(rng.uniform(0.2, 0.95)))
        nubar = bool(rng.random() < 0.3)
        N = int(rng.choice([1, 3, 12, 30, 80]))
        Ls = np.linspace(0.05*span, span, N) if N > 1 else np.array([span])
        vcc = H.vcc_of(ne, nubar=nubar)
        H_of_l = H6.H_family('std', d, E, vcc, p, nubar=nubar)
        try:
            with H.Caught() as c:
                P = np.asarray(oscprob.osc_prob_matter_std_potential(
                    d, ne, E, Ls if N > 1 else float(Ls[0]), p, L0=0.0, nubar=nubar,
                    density_is_of_number_of_electrons=True)).reshape(N, d, d)
            Pref = np.array([H.P_of(U) for U in H.exact_U_many(H_of_l, 0.0, Ls, d)])
        except Exception as ex:            # noqa: BLE001
            print('   case %3d RAISED %s: %s' % (i, type(ex).__name__, str(ex)[:90]))
            continue
        err = H.maxabs(P - Pref)
        warned = any(x in WARN_OK for x in c.names)
        rows.append(dict(i=i, d=d, E=E, N=N, nubar=nubar, err=err, warned=warned,
                         names=c.names, span=span, n_modes=n_modes))
        if (i + 1) % 25 == 0:
            print('   ... %d/%d  (%.0f s)' % (i + 1, n_cases, time.time() - t0), flush=True)

    errs = np.array([r['err'] for r in rows])
    silent = [r for r in rows if r['err'] > TOL and not r['warned']]
    warned_bad = [r for r in rows if r['err'] > TOL and r['warned']]
    print('\n  cases scored          : %d' % len(rows))
    print('  error distribution    : median %.2e  p90 %.2e  p99 %.2e  max %.2e'
          % (np.median(errs), np.percentile(errs, 90), np.percentile(errs, 99), errs.max()))
    print('  outside 1e-3, WARNED  : %d' % len(warned_bad))
    print('  outside 1e-3, SILENT  : %d   <-- criterion (E1)' % len(silent))
    for r in sorted(silent, key=lambda r: -r['err'])[:25]:
        print('     case %3d d=%d N=%-3d E=%6.1f MeV nubar=%-5s err=%.3e'
              % (r['i'], r['d'], r['N'], r['E']/1e6, r['nubar'], r['err']))
    for d in (2, 3, 4, 5):
        e = [r['err'] for r in rows if r['d'] == d]
        s = [r for r in silent if r['d'] == d]
        if e:
            print('  d=%d: n=%-3d median %.2e  max %.2e  silent misses %d'
                  % (d, len(e), float(np.median(e)), max(e), len(s)))
    np.save('battery5_fuzz.npy', np.array(rows, dtype=object), allow_pickle=True)
    print()


def sub5():
    """5.5 Adversarial caller-supplied n_slabs / n_tpts_per_slab."""
    print('## 5.5  deliberately absurd n_slabs / n_tpts_per_slab')
    p2 = H.params_for(2)
    ne = H.solar_ne()
    L = 0.5*H.R_SUN
    H_of_l = H6.H_family('std', 2, E_DEF, H.vcc_of(ne), p2)
    Pref = H.P_of(H.exact_U(H_of_l, 0.0, L, 2))
    for ns in (1, 2, 10**6):
        for nt in (2, 500):
            try:
                with H.Caught() as c:
                    P = np.asarray(oscprob.osc_prob_matter_std_potential(
                        2, ne, E_DEF, L, p2, L0=0.0,
                        density_is_of_number_of_electrons=True,
                        n_slabs=ns, n_tpts_per_slab=nt))
                err = H.maxabs(P - Pref)
                silent = err > TOL and not any(x in WARN_OK for x in c.names)
                print('  n_slabs=%-8d n_tpts=%-4d err=%.3e  %-11s %s'
                      % (ns, nt, err,
                         'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                         ','.join(c.names) or '-'), flush=True)
            except Exception as ex:        # noqa: BLE001
                print('  n_slabs=%-8d n_tpts=%-4d RAISED %s: %s'
                      % (ns, nt, type(ex).__name__, str(ex)[:80]), flush=True)
    print()


def sub6():
    """5.6 Determinism: repeat, shuffled order, n_jobs."""
    print('## 5.6  determinism (criterion E2)')
    p2 = H.params_for(2)
    ne = H.solar_ne()
    Ls = np.linspace(0.05*H.R_SUN, 0.5*H.R_SUN, 60)

    def call(L, **kw):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.asarray(oscprob.osc_prob_matter_std_potential(
                2, ne, E_DEF, L, p2, L0=0.0, density_is_of_number_of_electrons=True, **kw))

    a, b = call(Ls), call(Ls)
    print('  same call twice, max|diff|         : %.3e  %s'
          % (H.maxabs(a - b), 'ok' if H.maxabs(a - b) == 0.0 else '<-- NOT BIT-IDENTICAL'))

    rng = np.random.default_rng(1)
    perm = rng.permutation(len(Ls))
    c = call(Ls[perm])
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    dif = H.maxabs(a - c[inv])
    print('  shuffled input order, max|diff|    : %.3e  %s'
          % (dif, 'ok' if dif == 0.0 else '<-- ORDER-DEPENDENT'))

    for nj in (2, 4):
        e = call(Ls, n_jobs=nj)
        dif = H.maxabs(a - e)
        print('  n_jobs=%d vs n_jobs=1, max|diff|    : %.3e  %s'
              % (nj, dif, 'ok' if dif == 0.0 else '<-- DIFFERS'))
    print()


if __name__ == '__main__':
    args = sys.argv[1:] or ['1', '2', '3', '4', '5', '6']
    print('# Battery 5 -- designed to break.  tol=%.0e\n' % TOL)
    subs = {'1': sub1, '2': sub2, '3': sub3, '5': sub5, '6': sub6}
    i = 0
    while i < len(args):
        w = args[i]
        if w == '4':
            # '4' takes an optional case count as the NEXT argument; consume it here so
            # the dispatch loop never sees it as a sub-test name.
            n = 200
            if i + 1 < len(args) and args[i + 1].isdigit() and args[i + 1] not in subs:
                n = int(args[i + 1])
                i += 1
            sub4(n)
        else:
            subs[w]()
        i += 1
