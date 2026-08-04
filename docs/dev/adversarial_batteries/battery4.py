# -*- coding: utf-8 -*-
"""Battery 4 -- extreme numerics.

PASS CRITERIA, stated before running:
  (F1) every answer meets its requested tolerance OR warns.
  (F2) memory is O(output) in N -- the chunking claim -- measured with tracemalloc.
  (F3) lowering max_n_slabs must not degrade the answer SILENTLY (handover suspicion #3).
  (F4) geometry edge cases return correct answers or raise; never silently wrong.
"""

import sys
import time
import tracemalloc
import warnings

import numpy as np

import harness as H
import harness6 as H6
import magnus.oscprob as oscprob

TOL = 1e-3
WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning')
NE = H.solar_ne()


def std(*a, **k):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        *a, L0=k.pop('L0', 0.0), density_is_of_number_of_electrons=True, **k))


def sub1():
    """4.1 Energy from 0.05 MeV to 100 GeV.  At the low end the ORACLE is checked first."""
    print('## 4.1  energy extremes (oracle convergence verified at the low end)')
    p = H.params_for(2)
    for E_MeV, Lfrac in ((0.05, 0.02), (0.1, 0.02), (0.5, 0.05), (10.0, 0.3),
                         (1e3, 0.5), (1e4, 0.5), (1e5, 0.5)):
        E = E_MeV*1e6
        L = Lfrac*H.R_SUN
        Hf = H6.H_family('std', 2, E, H.vcc_of(NE), p)
        Uref = H.exact_U(Hf, 0.0, L, 2)
        with H.Caught() as c:
            P = std(2, NE, E, L, p)
        err = H.maxabs(P - H.P_of(Uref))
        moved, ok = H.oracle_converged(Hf, 0.0, L, 2, max(err, 1e-14))
        silent = err > TOL and not any(x in WARN_OK for x in c.names)
        print('  E=%9.3f MeV L=%.2f R_sun  err=%.3e  oracle moved %.1e (%s)  %-11s %s'
              % (E_MeV, Lfrac, err, moved, 'trustworthy' if ok else 'ORACLE SUSPECT',
                 'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                 ','.join(c.names) or '-'), flush=True)
    print()


def sub2():
    """4.2 N = 1e4, 1e5, 1e6: wall time and PEAK MEMORY (expect O(output), flat in N)."""
    print('## 4.2  very large N: wall time and peak memory (criterion F2)')
    p = H.params_for(2)
    LM = 0.3*H.R_SUN
    E = 50e6
    base = None
    for N in (10**3, 10**4, 10**5, 10**6):
        Ls = np.linspace(0.05*LM, LM, N)
        tracemalloc.start()
        t0 = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                P = std(2, NE, E, Ls, p)
            dt = time.time() - t0
            cur, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            out_mb = P.nbytes/1e6
            if base is None:
                base = peak/1e6
            print('  N=%-8d %7.2f s  peak %8.1f MB  output %7.1f MB  '
                  'peak-minus-output %7.1f MB' % (N, dt, peak/1e6, out_mb,
                                                  peak/1e6 - out_mb), flush=True)
        except Exception as ex:            # noqa: BLE001
            tracemalloc.stop()
            print('  N=%-8d RAISED %s: %s' % (N, type(ex).__name__, str(ex)[:110]),
                  flush=True)
    print()


def sub3():
    """4.3 Tolerances 1e-1 .. 1e-12: met, or ToleranceNotAchievedWarning.

    The suppression added in commit e34cbeb (MagnusConvergenceWarning silenced on the
    cumulative probe) must not have swallowed the signal that matters."""
    print('## 4.3  tolerance sweep -- met, or warned?  (the commit-5 suppression check)')
    p = H.params_for(2)
    LM = 0.5*H.R_SUN
    E = 50e6
    Ls = np.linspace(0.05*LM, LM, 60)
    Hf = H6.H_family('std', 2, E, H.vcc_of(NE), p)
    Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, 0.0, Ls, 2)])
    for tol in (1e-1, 1e-3, 1e-6, 1e-9, 1e-12):
        with H.Caught() as c:
            P = std(2, NE, E, Ls, p, rtol=tol, atol=tol).reshape(60, 2, 2)
        err = H.maxabs(P - Pref)
        silent = err > tol and not any(x in WARN_OK for x in c.names)
        print('  tol=%.0e  err=%.3e  %-11s warns=%s'
              % (tol, err, 'SILENT MISS' if silent else ('warned' if err > tol else 'ok'),
                 ','.join(c.names) or '-'), flush=True)
    print()


def sub4():
    """4.4 max_n_slabs lowered -- handover suspicion #3: does the answer degrade SILENTLY?

    Over a full solar radius the strict probe reaches max_n_slabs without converging, so
    n_acc is ceiling-derived: lowering the cap lowers the whole scan's resolution in
    proportion.  Untested below the default until now."""
    print('## 4.4  max_n_slabs lowered below the default (criterion F3)')
    p = H.params_for(2)
    LM = 1.0*H.R_SUN
    E = 10e6
    Ls = np.linspace(0.05*LM, LM, 60)
    Hf = H6.H_family('std', 2, E, H.vcc_of(NE), p)
    print('  building the solve_ivp reference over a FULL solar radius at 10 MeV ...',
          flush=True)
    t0 = time.time()
    Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, 0.0, Ls, 2)])
    print('  ... %.0f s' % (time.time() - t0), flush=True)
    for cap in (500, 2000, 5000, 20000, None):
        probe = {}
        with H.Caught() as c:
            P = std(2, NE, E, Ls, p, max_n_slabs=cap).reshape(60, 2, 2)
        err = H.maxabs(P - Pref)
        # what n_acc did the probe report at this cap?
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            oscprob.osc_prob(Hf, 0.0, float(Ls[-1]), rtol=TOL, atol=TOL,
                             strict_convergence=True, convergence_info=probe,
                             max_n_slabs=cap)
        silent = err > TOL and not any(x in WARN_OK for x in c.names)
        print('  max_n_slabs=%-6s probe n_slabs=%-6s err=%.3e  %-11s warns=%s'
              % (cap, probe.get('n_slabs'), err,
                 'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                 ','.join(c.names) or '-'), flush=True)
    print()


def sub5():
    """4.5 Geometry: L0 != 0, huge spans, unsorted/duplicated/degenerate baselines."""
    print('## 4.5  geometry edge cases (criterion F4)')
    p = H.params_for(2)
    E = 50e6
    LM = 0.5*H.R_SUN

    def check(name, Ls, L0=0.0):
        Ls = np.asarray(Ls, dtype=float)
        Hf = H6.H_family('std', 2, E, H.vcc_of(NE), p)
        try:
            with H.Caught() as c:
                P = np.asarray(std(2, NE, E, Ls, p, L0=L0)).reshape(len(Ls), 2, 2)
        except Exception as ex:            # noqa: BLE001
            print('  %-44s RAISED %s: %s' % (name, type(ex).__name__, str(ex)[:70]))
            return
        # oracle: sort, solve once, map back (handles duplicates and unsorted input)
        order = np.argsort(Ls, kind='stable')
        Ls_s = Ls[order]
        uniq, inv = np.unique(Ls_s, return_inverse=True)
        keep = uniq[uniq > L0]
        Uu = H.exact_U_many(Hf, L0, keep, 2) if len(keep) else np.empty((0, 2, 2))
        Pu = {float(k): H.P_of(U) for k, U in zip(keep, Uu)}
        Pu[float(L0)] = np.eye(2)
        Pref_s = np.array([Pu[float(u)] for u in uniq])[inv]
        Pref = np.empty_like(Pref_s)
        Pref[order] = Pref_s
        err = H.maxabs(P - Pref)
        silent = err > TOL and not any(x in WARN_OK for x in c.names)
        print('  %-44s err=%.3e  %-11s %s'
              % (name, err, 'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
                 ','.join(c.names) or '-'), flush=True)

    check('L0 = 0, ordinary scan', np.linspace(0.05*LM, LM, 30))
    check('L0 mid-profile (0.2 R_sun)', np.linspace(0.25*H.R_SUN, 0.6*H.R_SUN, 30),
          L0=0.2*H.R_SUN)
    check('span of 8 orders of magnitude',
          np.geomspace(1e-8*LM, LM, 40))
    check('UNSORTED baselines', np.random.default_rng(2).permutation(
        np.linspace(0.05*LM, LM, 30)))
    check('DUPLICATED baselines',
          np.repeat(np.linspace(0.05*LM, LM, 15), 2))
    check('L == L0 exactly included',
          np.r_[0.0, np.linspace(0.05*LM, LM, 29)])
    check('all baselines identical', np.full(20, 0.3*LM))
    print()


if __name__ == '__main__':
    print('# Battery 4 -- extreme numerics.  tol=%.0e\n' % TOL)
    for w in (sys.argv[1:] or ['1', '2', '3', '4', '5']):
        {'1': sub1, '2': sub2, '3': sub3, '4': sub4, '5': sub5}[w]()
