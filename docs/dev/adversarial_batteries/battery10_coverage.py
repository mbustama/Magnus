# -*- coding: utf-8 -*-
"""Battery 10 -- do the Earth and solar cases run, for every flavour count, family and
density, in reasonable time and to reasonable accuracy?

This is a coverage question rather than an adversarial one: not "can I break it" but "does the
whole surface a user actually touches work".  Axes swept:

  * environment  : Sun (exponential) and Earth (PREM), plus constant density and vacuum
  * flavours     : d = 2, 3, 4, 5
  * family       : standard matter, NSI, LIV (n_liv = 0 and 1)
  * antineutrinos: both
  * shape        : single point and baseline scan (either side of the N = 25 seam)
  * entry point  : the scenario wrappers and the generic user-Hamiltonian ones

PASS CRITERIA, stated before running:
  (J1) ACCURACY  -- meets the requested 1e-3 against solve_ivp/DOP853, or warns.
  (J2) TIME      -- no single call over SLOW_S seconds; nothing hangs.
  (J3) COMPLETES -- no exceptions anywhere.
Every row is reported, not only failures, and every warning is named.
"""

import sys
import time
import traceback
import warnings

import numpy as np

import harness as H
import harness6 as H6
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.oscprob as op

TOL = 1e-3
SLOW_S = 5.0
WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning', 'UnmarkedDiscontinuityWarning')
ROWS = []


def record(env, d, fam, nubar, N, err, dt, names, exc=None):
    silent = (exc is None) and err > TOL and not any(x in WARN_OK for x in names)
    ROWS.append(dict(env=env, d=d, fam=fam, nubar=nubar, N=N, err=err, dt=dt,
                     names=names, silent=silent, exc=exc))
    if exc is not None:
        flag = 'RAISED'
    elif silent:
        flag = 'SILENT MISS'
    elif err > TOL:
        flag = 'warned'
    elif dt > SLOW_S:
        flag = 'SLOW'
    else:
        flag = 'ok'
    print('  %-8s d=%d %-11s %-6s N=%-3d err=%9.3e %7.2fs  %-11s %s'
          % (env, d, fam, 'nubar' if nubar else 'nu', N, err, dt, flag,
             (','.join(names) or '-') if exc is None else exc), flush=True)


def solar(d, fam, nubar, N, energy=50.0e6, l_max=None):
    l_max = l_max if l_max is not None else 0.4*H.R_SUN
    p = H.params_for(d)
    ne = H.solar_ne()
    Ls = np.linspace(0.05*l_max, l_max, N) if N > 1 else np.array([l_max])
    nsi = H6.nsi_params_for(d) if fam == 'nsi' else None
    liv = H6.liv_params_for(d, n_liv=int(fam[-1])) if fam.startswith('liv') else None
    famkey = 'liv' if fam.startswith('liv') else fam
    try:
        with H.Caught() as c:
            t0 = time.time()
            if famkey == 'std':
                P = op.osc_prob_matter_std_potential(
                    d, ne, energy, Ls if N > 1 else float(Ls[0]), p, L0=0.0, nubar=nubar,
                    density_is_of_number_of_electrons=True)
            elif famkey == 'nsi':
                P = op.osc_prob_matter_nsi(
                    d, ne, energy, Ls if N > 1 else float(Ls[0]), p, nsi, L0=0.0,
                    nubar=nubar, density_is_of_number_of_electrons=True)
            else:
                P = op.osc_prob_liv(
                    d, energy, Ls if N > 1 else float(Ls[0]), p, liv, rho_func=ne, L0=0.0,
                    nubar=nubar, density_is_of_number_of_electrons=True)
            dt = time.time() - t0
        P = np.asarray(P).reshape((N, d, d))
        Hf = H6.H_family(famkey, d, energy, H.vcc_of(ne, nubar=nubar), p,
                         nsi=nsi, liv=liv, nubar=nubar)
        Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, 0.0, Ls, d)])
        record('Sun', d, fam, nubar, N, H.maxabs(P - Pref), dt, c.names)
    except Exception as ex:                       # noqa: BLE001
        record('Sun', d, fam, nubar, N, float('nan'), 0.0, [],
               exc='%s: %s' % (type(ex).__name__, str(ex)[:70]))


def earth_case(d, costhz, E_GeV, nubar, N):
    E = E_GeV*gd.UNIT_GEV
    p = H.params_for(d)
    Lc = 2.0*gd.EARTH_RADIUS*abs(costhz)*gd.UNIT_KM
    Ls = np.linspace(0.05*Lc, Lc, N) if N > 1 else np.array([Lc])

    def rho(x, _c=costhz):
        return earth.density_matter_func_prem(x, _c)

    try:
        with H.Caught() as c:
            t0 = time.time()
            P = op.osc_prob_matter_std_potential(
                d, rho, E, Ls if N > 1 else float(Ls[0]), p, L0=0.0, nubar=nubar,
                density_matter_is_in_g_per_cm3=True,
                t_breakpoints=earth.prem_layer_edges_along_chord(costhz))
            dt = time.time() - t0
        P = np.asarray(P).reshape((N, d, d))
        vcc = matter.vcc_func_from_rho_func(rho, 0.0, 1.0, 0.5, nubar=nubar,
                                            density_matter_is_in_g_per_cm3=True,
                                            density_is_of_number_of_electrons=False)
        Hf = H6.H_family('std', d, E, vcc, p, nubar=nubar)
        Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, 0.0, Ls, d)])
        record('Earth', d, 'std c=%.1f' % costhz, nubar, N, H.maxabs(P - Pref), dt, c.names)
    except Exception as ex:                       # noqa: BLE001
        record('Earth', d, 'std c=%.1f' % costhz, nubar, N, float('nan'), 0.0, [],
               exc='%s: %s' % (type(ex).__name__, str(ex)[:70]))


def sub_sun():
    print('## 10.1  Sun: every flavour count x family x nu/nubar x scan size')
    for d in (2, 3, 4, 5):
        for fam in ('std', 'nsi', 'liv0', 'liv1'):
            for nubar in (False, True):
                for N in (1, 8, 40):
                    solar(d, fam, nubar, N)
    print()


def sub_earth():
    print('## 10.2  Earth: real PREM, several directions and energies')
    for d in (2, 3, 4, 5):
        for costhz in (-0.2, -0.6, -0.95):
            for E_GeV in (1.0, 10.0):
                for nubar in (False, True):
                    earth_case(d, costhz, E_GeV, nubar, 1)
    for d in (2, 3):
        for costhz in (-0.6, -0.95):
            earth_case(d, costhz, 5.0, False, 40)
    print()


def sub_other():
    print('## 10.3  vacuum and constant density, every flavour count')
    for d in (2, 3, 4, 5):
        for nubar in (False, True):
            for rho, nm in ((0.0, 'vacuum'), (0.05*H.NE0, 'const')):
                p = H.params_for(d)
                L = 0.4*H.R_SUN
                Ls = np.linspace(0.05*L, L, 40)
                try:
                    with H.Caught() as c:
                        t0 = time.time()
                        P = np.asarray(op.osc_prob_matter_std_potential(
                            d, rho, 50e6, Ls, p, L0=0.0, nubar=nubar,
                            density_is_of_number_of_electrons=True)).reshape(40, d, d)
                        dt = time.time() - t0
                    Hf = H6.H_family('std', d, 50e6,
                                     H.vcc_of(lambda x, _r=rho: np.full_like(
                                         np.asarray(x, dtype=float), _r), nubar=nubar),
                                     p, nubar=nubar)
                    Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, 0.0, Ls, d)])
                    record(nm, d, 'std', nubar, 40, H.maxabs(P - Pref), dt, c.names)
                except Exception as ex:           # noqa: BLE001
                    record(nm, d, 'std', nubar, 40, float('nan'), 0.0, [],
                           exc='%s: %s' % (type(ex).__name__, str(ex)[:70]))
    print()


def summary():
    print('=== BATTERY 10 SUMMARY ===')
    ok = [r for r in ROWS if r['exc'] is None]
    print('  configurations run       : %d' % len(ROWS))
    print('  raised an exception      : %d   <-- criterion (J3)'
          % sum(1 for r in ROWS if r['exc'] is not None))
    for r in ROWS:
        if r['exc'] is not None:
            print('     %s d=%d %s: %s' % (r['env'], r['d'], r['fam'], r['exc']))
    if ok:
        e = np.array([r['err'] for r in ok])
        t = np.array([r['dt'] for r in ok])
        print('  accuracy: median %.2e  p90 %.2e  max %.2e' % (np.median(e),
                                                               np.percentile(e, 90), e.max()))
        print('  time    : median %.3fs  p90 %.3fs  max %.3fs' % (np.median(t),
                                                                  np.percentile(t, 90), t.max()))
        print('  SILENT MISSES            : %d   <-- criterion (J1)'
              % sum(1 for r in ok if r['silent']))
        for r in ok:
            if r['silent']:
                print('     %s d=%d %s %s N=%d err=%.3e'
                      % (r['env'], r['d'], r['fam'], 'nubar' if r['nubar'] else 'nu',
                         r['N'], r['err']))
        slow = [r for r in ok if r['dt'] > SLOW_S]
        print('  calls over %.0fs           : %d   <-- criterion (J2)' % (SLOW_S, len(slow)))
        for r in sorted(slow, key=lambda r: -r['dt'])[:8]:
            print('     %s d=%d %s N=%d  %.2fs' % (r['env'], r['d'], r['fam'], r['N'], r['dt']))
    np.save('battery10_rows.npy', np.array(ROWS, dtype=object), allow_pickle=True)


if __name__ == '__main__':
    warnings.simplefilter('always')
    print('# Battery 10 -- Earth and solar coverage.  tol=%.0e, slow>%.0fs\n' % (TOL, SLOW_S))
    subs = {'1': sub_sun, '2': sub_earth, '3': sub_other}
    try:
        for w in (sys.argv[1:] or ['1', '2', '3']):
            subs[w]()
    except Exception:                             # noqa: BLE001
        traceback.print_exc()
    summary()
