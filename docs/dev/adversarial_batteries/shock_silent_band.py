# -*- coding: utf-8 -*-
"""Map the silent band: which shock widths are wrong, unflagged AND unwarned?

`FINDINGS_ROBUSTNESS_PROGRAMME.md` §13.19 left exactly one open question. Of everything the
physical population measured, only a supernova shock produces an error that survives phase
averaging -- and the sharp shocks that do it worst are **loud**: flagged 144/144 by the
resolution test, warned every time, cured to 8.75e-06 by ``t_breakpoints``.

But at w = 1e-3 (a 70 km front on a 7e4 km ray) the error is 9.773e-04 on the averaged
observable -- inside the requested 1e-3 by two per cent -- and it is **neither flagged nor
warned**. That matters more than the number suggests, because a real hydrodynamic shock is
mean-free-path thin (w ~ 1e-6, loud) and **nobody hands this package a real shock**: they hand it
a simulation snapshot, where the front is smeared across a few grid cells, i.e. tens of km. The
silent band is plausibly the most likely form a user's shock actually arrives in, and two per
cent of margin is luck rather than headroom.

This sweeps width x energy and reports, per configuration:

  * the error on the **averaged** observable -- the lesson of §13.17-13.19 is that the
    instantaneous error at one baseline is not what a user sees;
  * whether any warning fired;
  * what the resolution test says at ``n_probe`` = 200 and 6400, which is what decides
    ``UnmarkedDiscontinuityWarning``;
  * the error with ``t_breakpoints`` declared, which confirms the front is the cause rather
    than something else about the profile.

A configuration is **SILENT-WRONG** when the averaged error is outside tolerance and nothing
warned. Those are the only rows that matter.

Run:  python shock_silent_band.py
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.globaldefs as gd
import magnus.oscprob as op
import physical_profiles as pp

TOL = 1e-3
D = 3
WIDTHS = (3.0e-2, 1.0e-2, 5.0e-3, 3.0e-3, 2.0e-3, 1.0e-3, 5.0e-4, 3.0e-4, 1.0e-4)
ENERGIES = (15.0e6, 45.0e6)
N_AVG = 41


def averaged_error(ne, energy, l0, l1, params, l_osc, kw=None):
    """(averaged error, instantaneous error at the nominal baseline, warned)."""
    Ls = np.linspace(max(l1 - 6.0*l_osc, l0 + 0.5*(l1 - l0)), l1, N_AVG)
    Hf = H.H_factory(D, params, H.vcc_of(ne), energy)
    Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, l0, Ls, D)])
    got = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        for L in Ls:
            got.append(np.asarray(op.osc_prob_matter_std_potential(
                D, ne, energy, float(L), params, L0=l0,
                density_is_of_number_of_electrons=True, **(kw or {}))))
        warned = sorted({w.category.__name__ for w in caught})
    got = np.array(got)
    return (float(np.max(np.abs(got.mean(axis=0) - Pref.mean(axis=0)))),
            float(np.max(np.abs(got[-1] - Pref[-1]))), warned)


def main():
    params = H.params_for(D)
    p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    span_km = pp.SN_R1_KM - pp.SN_R0_KM
    print('# Silent-band map for the supernova shock.  d = %d, tolerance %.0e' % (D, TOL))
    print('# averaged over %d points spanning 6 oscillation lengths (D31)\n' % N_AVG)
    print('%9s %8s %7s %10s %10s %8s %8s %11s  %s'
          % ('w', 'w/km', 'E/MeV', 'avg err', 'inst err', 'res200', 'res6400', 'w/ breakpts',
             'verdict'))
    rows = []
    for w in WIDTHS:
        ne = pp.sn_shock_ne(w)
        bps = pp.sn_shock_breakpoints(w)
        for energy in ENERGIES:
            l_osc = 4.0*np.pi*energy/p['D31']
            t0 = time.time()
            try:
                avg, inst, warned = averaged_error(ne, energy, pp.SN_L0, pp.SN_L1, params, l_osc)
                avg_bp, _, _ = averaged_error(ne, energy, pp.SN_L0, pp.SN_L1, params, l_osc,
                                              kw={'t_breakpoints': bps})
            except Exception as exc:                       # noqa: BLE001
                print('%9.0e %8.2f %7.0f  %s' % (w, w*span_km, energy/1e6, type(exc).__name__))
                continue
            Hf = H.H_factory(D, params, H.vcc_of(ne), energy)
            res200 = ad._profile_is_resolved(Hf, pp.SN_L0, pp.SN_L1, 200)
            res6400 = ad._profile_is_resolved(Hf, pp.SN_L0, pp.SN_L1, 6400)
            outside, silent = avg > TOL, (avg > TOL and not warned)
            verdict = ('SILENT-WRONG' if silent else
                       ('outside, warned' if outside else
                        ('inside, warned' if warned else 'inside, quiet')))
            rows.append(dict(w=w, energy=energy, avg=avg, inst=inst, warned=bool(warned),
                             res200=bool(res200), res6400=bool(res6400), avg_bp=avg_bp,
                             silent=silent))
            print('%9.0e %8.2f %7.0f %10.3e %10.3e %8s %8s %11.3e  %s   (%.0fs)'
                  % (w, w*span_km, energy/1e6, avg, inst, res200, res6400, avg_bp, verdict,
                     time.time() - t0), flush=True)

    print('\n=== THE BAND ===')
    silent = [r for r in rows if r['silent']]
    if not silent:
        print('no configuration is outside tolerance on the averaged observable AND unwarned.')
        worst = max(rows, key=lambda r: r['avg'] if not r['warned'] else -1.0)
        print('closest approach among the quiet rows: w=%.0e E=%.0f MeV avg=%.3e (%.0f%% of tol)'
              % (worst['w'], worst['energy']/1e6, worst['avg'], 100.0*worst['avg']/TOL))
    else:
        print('SILENT-WRONG configurations: %d' % len(silent))
        for r in silent:
            print('   w=%.0e (%.1f km) E=%.0f MeV  avg err %.3e = %.1fx tolerance'
                  % (r['w'], r['w']*span_km, r['energy']/1e6, r['avg'], r['avg']/TOL))

    # The hypothesis in §13.19: the resolution test keys on the probe grid, so a front that the
    # 6400-point grid resolves is declared fine even when the transport grid straddles it.
    print('\n=== probe-grid hypothesis ===')
    print('rows where the 6400 probe grid says RESOLVED but the answer is outside tolerance:')
    bad = [r for r in rows if r['res6400'] and r['avg'] > TOL]
    for r in bad:
        print('   w=%.0e E=%.0f MeV  avg %.3e  warned=%s' % (r['w'], r['energy']/1e6, r['avg'],
                                                             r['warned']))
    print('   (%d such rows)' % len(bad))
    print('\nt_breakpoints cure, worst remaining error with them declared: %.3e'
          % max(r['avg_bp'] for r in rows))
    np.save('shock_band_rows.npy', np.array(rows, dtype=object), allow_pickle=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
