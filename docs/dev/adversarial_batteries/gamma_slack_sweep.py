# -*- coding: utf-8 -*-
"""Is GAMMA_SLACK = 2.0 justified, or did four points get lucky?

GAMMA_SLACK decides when hybrid_propagator may certify a result with no non-adiabatic window:
it does so when ``GAMMA_TO_ERROR * gamma_max <= GAMMA_SLACK * (atol + rtol)``.  With
GAMMA_TO_ERROR = 1.0 that is simply ``gamma_max <= 2 * tolerance``.

The constant was chosen from a gap between two numbers -- 1.42 (a case that must certify) and
2.16 (a case that must refine) -- which is a factor-of-1.5 margin resting on four measurements.
This measures the underlying quantity densely instead.

THE QUANTITY.  Write the pure adiabatic answer's error as

    |dP|  =  k * gamma_max .

Then certifying an empty window list is safe exactly when k * gamma_max <= tolerance, i.e. when
gamma_max <= tolerance / k.  So **GAMMA_SLACK is 1/k_max**, and the question "is 2.0 right?" is
the question "is k_max <= 0.5?".

METHOD.  For each configuration, compute the PURE adiabatic operator (adiabatic_propagator, no
windows, no patching -- this is exactly what gets returned when no window opens), score it
against solve_ivp/DOP853, and divide by gamma_max measured on the same probe grid.  Sweeping
resonance width, energy, domain length and flavour count moves gamma_max over several decades.

Reports the distribution of k, and the implied safe slack 1/max(k).
"""

import sys
import warnings

import numpy as np

import harness as H
import harness6 as H6
import magnus.adiabatic as ad
from battery2 import bump_profile, ne_res_for

warnings.simplefilter('ignore')


def k_for(d, energy, span, width_frac, base=0.30, peak=3.0, centre=0.5):
    """k = |dP_adiabatic| / gamma_max for one configuration, plus the pieces."""
    p = H.params_for(d)
    ne_res = ne_res_for(d, p, energy)
    ne = bump_profile(ne_res, centre*span, width_frac*span, base_frac=base, peak_frac=peak)
    H_of_l = H6.H_family('std', d, energy, H.vcc_of(ne), p)

    info = {}
    ad.find_nonadiabatic_windows(H_of_l, 0.0, span, threshold=0.1, n_probe=200, info=info)
    gmax = float(info['gamma_max'])
    if not np.isfinite(gmax) or gmax <= 0.0:
        return None

    # The pure adiabatic answer: what hybrid returns when no window opens.  n_points at the
    # ceiling so this measures the ADIABATIC APPROXIMATION's error, not quadrature error.
    U_ad = ad.adiabatic_propagator(H_of_l, 0.0, span, n_points=12864)
    P_ad = H.P_of(U_ad)
    P_ex = H.P_of(H.exact_U(H_of_l, 0.0, span, d))
    err = H.maxabs(P_ad - P_ex)
    return dict(d=d, E=energy, span=span, wf=width_frac, gamma_max=gmax, err=err,
                k=err/gmax)


def main():
    n_arg = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    rows = []
    print('# GAMMA_SLACK sweep: k = |dP_adiabatic| / gamma_max')
    print('# GAMMA_SLACK is 1/k_max, so 2.0 is right iff k_max <= 0.5\n')
    print('%2s %9s %8s %8s %11s %11s %7s' % ('d', 'E[MeV]', 'span/ls', 'width', 'gamma_max',
                                             '|dP|', 'k'))

    grid = []
    for d in (2, 3):
        for E_MeV in (5.0, 10.0, 30.0, 80.0):
            for span_ls in (0.5, 1.0, 2.0):
                for wf in (3e-1, 1e-1, 6e-2, 4e-2, 2e-2, 1e-2):
                    grid.append((d, E_MeV*1e6, span_ls*H.L_SCALE, wf))
    for d in (4, 5):
        for E_MeV in (10.0, 50.0):
            for span_ls in (0.5, 1.0):
                for wf in (1e-1, 4e-2, 2e-2):
                    grid.append((d, E_MeV*1e6, span_ls*H.L_SCALE, wf))
    if n_arg:
        grid = grid[:n_arg]

    for (d, E, span, wf) in grid:
        try:
            r = k_for(d, E, span, wf)
        except Exception as ex:                    # noqa: BLE001
            print('  d=%d E=%.0f wf=%.0e RAISED %s' % (d, E/1e6, wf, type(ex).__name__))
            continue
        if r is None:
            continue
        rows.append(r)
        print('%2d %9.1f %8.2f %8.0e %11.3e %11.3e %7.3f'
              % (r['d'], r['E']/1e6, r['span']/H.L_SCALE, r['wf'], r['gamma_max'],
                 r['err'], r['k']), flush=True)

    ks = np.array([r['k'] for r in rows])
    print('\n=== k distribution over %d configurations ===' % len(ks))
    for q in (50, 90, 95, 99, 100):
        print('  p%-3d  k = %.3f   -> implied safe slack 1/k = %.2f'
              % (q, np.percentile(ks, q), 1.0/max(np.percentile(ks, q), 1e-12)))
    print('  mean %.3f   min %.3f   max %.3f' % (ks.mean(), ks.min(), ks.max()))
    print('\n  GAMMA_SLACK currently 2.0, i.e. it assumes k <= 0.500')
    over = [r for r in rows if r['k'] > 0.5]
    print('  configurations with k > 0.5 (where 2.0 is optimistic): %d of %d'
          % (len(over), len(rows)))
    for r in sorted(over, key=lambda r: -r['k'])[:10]:
        print('     d=%d E=%6.1f MeV span=%.1f ls wf=%.0e  k=%.3f  (gamma_max %.2e, |dP| %.2e)'
              % (r['d'], r['E']/1e6, r['span']/H.L_SCALE, r['wf'], r['k'],
                 r['gamma_max'], r['err']))
    if len(ks):
        print('\n  SAFE SLACK from this population: %.2f  (= 1/k_max)' % (1.0/ks.max()))
    np.save('gamma_slack_rows.npy', np.array(rows, dtype=object), allow_pickle=True)


if __name__ == '__main__':
    main()
