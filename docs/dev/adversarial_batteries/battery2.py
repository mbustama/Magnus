# -*- coding: utf-8 -*-
"""Battery 2 -- aimed at the detector's fixed probe grid.

PREMISE (from the handover, suspicion #1): find_resonance_candidates /
find_nonadiabatic_windows sample n_probe = 200 LINEAR points (doubling to 6400 under
refinement).  Structure narrower than (l1-l0)/200 can be stepped straight over.  The gamma
sweep added in 9c7945a inherits exactly this grid.

PASS CRITERION, stated before running:
  (A) HONESTY: hybrid_propagator returning certified=True must imply
      max|P_hybrid - P_solve_ivp| <= atol + rtol*1  (i.e. <= 1e-3 at the default).
      A certified answer outside tolerance is a FAILURE, and is the headline result.
  (B) DETECTION: informational -- how many windows are found, and whether refinement
      (n_probe doubling to 6400) rescues a case the initial grid misses.
  (A) is the pass/fail criterion.  (B) is reported for every configuration, not only
  failures.

Oracle: solve_ivp/DOP853 at rtol=1e-12, atol=1e-14, ONLY.  No Magnus path is ever scored
against another.
"""

import sys
import time

import numpy as np

import harness as H
import magnus.adiabatic as ad

RTOL = ATOL = 1.0e-3

# Domain and energy.  Chosen so the accumulated vacuum phase (~1.2e3 rad) is large enough
# to be a real integration problem but small enough that the DOP853 oracle is affordable.
L0, L1 = 0.0, 1.0*H.L_SCALE
E = 10.0e6  # eV

# Conversion from electron number density [eV^3] to V_CC [eV], read off the library itself
# rather than re-derived (the constant is exactly what vcc_func_from_rho_func applies).
_C_VCC = float(np.asarray(H.vcc_of(lambda l: np.full_like(np.asarray(l, dtype=float),
                                                          H.NE0))(0.0)))/H.NE0


def ne_res_for(d, params, energy):
    """Electron number density at which levels 0,1 are resonant, found by minimising the
    0-1 gap of the constant-density Hamiltonian over ne.  Works for any d."""
    h_vac = np.asarray(H.h_vac_for(d, params), dtype=complex)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0

    def gap(ne):
        lam = np.linalg.eigvalsh(h_vac/energy + ne*_C_VCC*proj)
        return lam[1] - lam[0]

    lo, hi = H.NE0*1e-6, H.NE0*10.0
    xs = np.geomspace(lo, hi, 4000)
    gs = np.array([gap(x) for x in xs])
    i = int(np.argmin(gs))
    a, b = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    for _ in range(200):
        m1, m2 = a + (b - a)/3.0, b - (b - a)/3.0
        if gap(m1) < gap(m2):
            b = m2
        else:
            a = m1
    return 0.5*(a + b)


def bump_profile(ne_res, l_center, width, base_frac=0.30, peak_frac=3.0):
    """Quiet background at base_frac*ne_res with ONE Gaussian bump peaking at
    peak_frac*ne_res: the trajectory crosses the resonance density exactly twice, and the
    crossing sharpness is set by `width` alone."""
    def ne(l):
        x = np.asarray(l, dtype=float)
        return ne_res*(base_frac + (peak_frac - base_frac)
                       * np.exp(-0.5*((x - l_center)/width)**2))
    return ne


def many_bumps_profile(ne_res, centers, width, base_frac=0.30, peak_frac=3.0):
    centers = np.asarray(centers, dtype=float)

    def ne(l):
        x = np.asarray(l, dtype=float)[..., None]
        b = np.exp(-0.5*((x - centers)/width)**2).sum(axis=-1)
        return ne_res*(base_frac + (peak_frac - base_frac)*np.minimum(b, 1.0))
    return ne


def run_case(label, ne_func, d, params, energy=E, l0=L0, l1=L1, nubar=False,
             oracle_check=False, quiet=False):
    """One configuration: hybrid_propagator against solve_ivp, with certification."""
    vcc = H.vcc_of(ne_func, nubar=nubar)
    H_of_l = H.H_factory(d, params, vcc, energy, nubar=nubar)

    t0 = time.time()
    U, windows, certified = ad.hybrid_propagator(H_of_l, l0, l1, rtol=RTOL, atol=ATOL)
    t_hyb = time.time() - t0

    t0 = time.time()
    U_ref = H.exact_U(H_of_l, l0, l1, d)
    t_ref = time.time() - t0

    err = H.maxabs(H.P_of(U) - H.P_of(U_ref))
    tol = ATOL + RTOL*1.0
    dishonest = bool(certified and err > tol)

    # How many windows the *initial* (unrefined) grid finds, vs the refined ceiling.
    w200, _ = ad.find_nonadiabatic_windows(H_of_l, l0, l1, threshold=0.1, n_probe=200)
    w6400, _ = ad.find_nonadiabatic_windows(H_of_l, l0, l1, threshold=0.1, n_probe=6400)

    oracle_moved = None
    if oracle_check:
        oracle_moved, _ = H.oracle_converged(H_of_l, l0, l1, d, max(err, 1e-12))

    row = dict(label=label, d=d, err=err, certified=certified, n_win=len(windows),
               n_win_200=len(w200), n_win_6400=len(w6400), unit=H.unitarity(U),
               dishonest=dishonest, t_hyb=t_hyb, t_ref=t_ref, oracle_moved=oracle_moved)
    if not quiet:
        print('%-42s d=%d  err=%9.3e  cert=%-5s  win(final/200/6400)=%2d/%2d/%2d  '
              'unit=%8.2e  %s%s'
              % (label, d, err, certified, len(windows), len(w200), len(w6400),
                 row['unit'], 'FAIL-DISHONEST' if dishonest else 'ok',
                 '' if oracle_moved is None else '  [oracle moved %.1e]' % oracle_moved),
              flush=True)
    return row


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    rows = []
    span = L1 - L0
    print('# Battery 2 -- detector probe grid.  domain=[%.3e, %.3e], probe spacing '
          '(n=200) = %.3e' % (L0, L1, span/199.0))
    print('# tolerance = %.1e; FAIL = certified=True with err > tolerance\n' % (ATOL + RTOL))

    p2 = H.params_for(2)
    ner2 = ne_res_for(2, p2, E)
    print('# 2nu resonance ne = %.4e eV^3  (= %.3f x solar central)\n' % (ner2, ner2/H.NE0))

    # ---- 2.1 single narrow resonance, width swept 1e-1 .. 1e-5 of the domain -----------
    if which in ('all', '1'):
        print('## 2.1  one Gaussian crossing at a random position, width swept')
        rng = np.random.default_rng(7)
        lc = L0 + (0.37 + 0.2*rng.random())*span   # "random position", pinned by seed
        for wf in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]:
            w = wf*span
            rows.append(run_case('2.1 w/(l1-l0)=%.0e @ l=%.3f' % (wf, lc/span),
                                 bump_profile(ner2, lc, w), 2, p2,
                                 oracle_check=(wf <= 1e-3)))
        print()

    # ---- 2.2 aliasing: sinusoid at/near the probe period -------------------------------
    if which in ('all', '2'):
        print('## 2.2  sinusoidal profile, period swept through the probe resonances')
        # period == span/199 is exactly the probe spacing; rational multiples alias too.
        for mult, name in [(1.0, 'span/199 (= probe spacing)'),
                           (2.0, '2 x probe spacing'),
                           (0.5, '1/2 probe spacing'),
                           (199/100.0, 'span/100'),
                           (1.01, '1.01 x probe spacing (beat)'),
                           (0.99, '0.99 x probe spacing (beat)'),
                           (199/7.0, 'span/7 (well resolved)')]:
            per = mult*span/199.0
            def ne(l, _p=per):
                return ner2*(1.0 + 0.9*np.sin(2.0*np.pi*np.asarray(l, float)/_p))
            rows.append(run_case('2.2 period = %-28s' % name, ne, 2, p2))
        print()

    # ---- 2.3 edge crossings ------------------------------------------------------------
    if which in ('all', '3'):
        print('## 2.3  crossings at/near the domain edges and on/off a probe point')
        dl = span/199.0
        w = 0.02*span
        for name, lc in [('within one probe spacing of l0', L0 + 0.5*dl),
                         ('within one probe spacing of l1', L1 - 0.5*dl),
                         ('exactly ON a probe point (i=100)', L0 + 100*dl),
                         ('exactly MIDWAY between probes', L0 + 100.5*dl),
                         ('exactly at l0', L0),
                         ('exactly at l1', L1)]:
            rows.append(run_case('2.3 %-34s' % name, bump_profile(ner2, lc, w), 2, p2))
        # and the same with a NARROW bump, where the position matters far more
        w = 1e-3*span
        for name, lc in [('narrow, ON a probe point', L0 + 100*dl),
                         ('narrow, MIDWAY between probes', L0 + 100.5*dl)]:
            rows.append(run_case('2.3 %-34s' % name, bump_profile(ner2, lc, w), 2, p2))
        print()

    # ---- 2.4 many crossings, more than the probe has points ----------------------------
    if which in ('all', '4'):
        print('## 2.4  many crossings (up to more than the probe has points)')
        for n in [10, 50, 100, 200, 400]:
            centers = L0 + span*(np.arange(n) + 0.5)/n
            w = 0.25*span/n
            rows.append(run_case('2.4 %3d crossings' % n,
                                 many_bumps_profile(ner2, centers, w), 2, p2))
        print()

    # ---- 2.5 clustered crossings -------------------------------------------------------
    if which in ('all', '5'):
        print('## 2.5  all crossings packed into 1%% of the range, rest quiet')
        for n in [5, 20, 50]:
            lo = L0 + 0.5*span
            centers = lo + 0.01*span*(np.arange(n) + 0.5)/n
            w = 0.25*0.01*span/n
            rows.append(run_case('2.5 %2d crossings in 1%% of range' % n,
                                 many_bumps_profile(ner2, centers, w), 2, p2))
        print()

    # ---- 2.7 the narrow-resonance sweep at d = 4 and 5 ---------------------------------
    if which in ('all', '7'):
        print('## 2.7  narrow-resonance sweep at d = 4 and 5 (more level pairs to mask)')
        for d in (4, 5):
            pp = H.params_for(d)
            nerd = ne_res_for(d, pp, E)
            print('#   d=%d resonance ne = %.4e (= %.3f x solar central)'
                  % (d, nerd, nerd/H.NE0))
            rng = np.random.default_rng(7)
            lc = L0 + (0.37 + 0.2*rng.random())*span
            for wf in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
                rows.append(run_case('2.7 d=%d w/(l1-l0)=%.0e' % (d, wf),
                                     bump_profile(nerd, lc, wf*span), d, pp))
        print()

    # ---- summary ------------------------------------------------------------------------
    print('\n=== BATTERY 2 SUMMARY ===')
    bad = [r for r in rows if r['dishonest']]
    print('configurations run          : %d' % len(rows))
    print('certified=True              : %d' % sum(r['certified'] for r in rows))
    print('CERTIFIED BUT OUT OF TOL    : %d   <-- pass criterion (A)' % len(bad))
    print('worst error overall         : %.3e' % max(r['err'] for r in rows))
    print('worst unitarity             : %.3e' % max(r['unit'] for r in rows))
    if bad:
        print('\nFAILURES (certified=True, err > %.1e):' % (ATOL + RTOL))
        for r in sorted(bad, key=lambda r: -r['err']):
            print('  %-44s d=%d err=%.3e  windows final/200/6400 = %d/%d/%d'
                  % (r['label'], r['d'], r['err'], r['n_win'], r['n_win_200'],
                     r['n_win_6400']))
    np.save('battery2_rows.npy', np.array(rows, dtype=object), allow_pickle=True)


if __name__ == '__main__':
    main()
