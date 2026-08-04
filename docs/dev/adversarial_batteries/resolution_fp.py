# -*- coding: utf-8 -*-
"""False-positive rate of the probe-scale resolution test, over SUB-INTERVALS.

The original measurement behind ``RESOLUTION_RATIO`` swept profile families, dimensions and
energies, but always on the **full** trajectory ``[l0, l1]``.  An ordinary baseline scan calls
the test once per baseline, i.e. on thirty different sub-intervals of the same profile, and the
statistic it computes turns out to depend on where a smooth extremum happens to fall inside its
probe interval.  This sweeps that axis, which is the one the earlier measurement did not have.

Reports, per profile family:

  * how many (d, E, sub-interval) configurations are declared UNRESOLVED -- false positives for
    a smooth family, true positives for a piecewise one;
  * the distribution of the **local confirmation** statistic on the intervals that were flagged,
    which is what sets ``LOCAL_JUMP_RATIO``.

Run:  python resolution_fp.py
"""

import sys

import numpy as np

import harness as H
import magnus.adiabatic as ad
from battery2 import bump_profile, ne_res_for
from battery3 import noisy_ne

L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0


def local_statistic(H_func, l0, l1, n_probe):
    """The largest local-confirmation fraction over the intervals the cheap test flags.

    Returns None when nothing is flagged (the common case), so that the distribution reported
    below is over the population the constant actually governs -- the mistake the
    ``GAMMA_TO_ERROR`` audit was written up to avoid.
    """
    ls = np.linspace(l0, l1, n_probe)
    mids = 0.5*(ls[:-1] + ls[1:])
    Hc = ad._H_on_grid(H_func, ls)
    Hm = ad._H_on_grid(H_func, mids)
    first = np.max(np.abs(Hm - Hc[:-1]), axis=(1, 2))
    second = np.max(np.abs(Hc[1:] - Hm), axis=(1, 2))
    total = first + second
    live = total > 1.0e-12*np.max(np.abs(Hc))
    if not np.any(live):
        return None
    live &= total > 0.25*np.median(total[live])
    if not np.any(live):
        return None
    ratio = np.where(live, np.maximum(first, second)/np.where(total > 0.0, total, 1.0), 0.0)
    flagged = np.where(ratio > ad.RESOLUTION_RATIO)[0]
    if flagged.size == 0:
        return None
    best = 0.0
    for i in flagged[np.argsort(-total[flagged])][:ad.MAX_LOCAL_CONFIRMATIONS]:
        xs = np.linspace(ls[i], ls[i + 1], ad.N_LOCAL_CONFIRM)
        steps = np.max(np.abs(np.diff(ad._H_on_grid(H_func, xs), axis=0)), axis=(1, 2))
        best = max(best, float(steps.max()/total[i]))
    return best


def families():
    p2 = H.params_for(2)
    ner2 = ne_res_for(2, p2, 10.0e6)
    rng = np.random.default_rng(20260804)
    lo, hi = 0.02*H.NE0, 0.30*H.NE0

    def step(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(np.where(x < 0.5*L1, lo, hi))

    def two_steps(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(np.where(x < 0.33*L1, lo, np.where(x < 0.66*L1, hi, 0.5*(lo + hi))))

    def castle(l, n_walls=12):
        x = np.asarray(l, dtype=float)
        u = np.clip((x - 0.1*L1)/(0.8*L1), 0.0, 1.0 - 1e-15)
        return H.scalarize(np.where(np.floor(u*n_walls) % 2 == 0, lo, hi))

    smooth = [
        ('solar exponential', H.solar_ne()),
        ('multi-resonance', H.modulated_ne(amp=0.9, n_cycles=6.0, span=SPAN)),
        ('noisy', noisy_ne()),
        ('sinusoid span/7', H.sine_ne(SPAN/7.0, base_ratio=3.0e-2)),
        ('gaussian bump w=1e-2', bump_profile(ner2, 0.45*SPAN, 1e-2*SPAN)),
        ('gaussian bump w=1e-1', bump_profile(ner2, 0.45*SPAN, 1e-1*SPAN)),
    ]
    for i in range(4):
        smooth.append(('random Fourier #%d' % i,
                       H.fourier_ne(rng, n_modes=int(rng.integers(2, 9)), span=SPAN,
                                    base_ratio=float(10.0**rng.uniform(-2.5, -1.0)),
                                    amp=float(rng.uniform(0.2, 0.9)))))
    piecewise = [('single step', step), ('two steps', two_steps), ('castle wall', castle)]
    return smooth, piecewise


def sweep(fams, kind):
    print('\n%-24s %10s %10s   %s' % (kind, 'unresolved', 'of', 'local statistic on flagged'))
    stats, n_bad_total, n_total = [], 0, 0
    for label, ne in fams:
        bad, tot, loc = 0, 0, []
        for d in (2, 3, 4, 5):
            params = H.params_for(d)
            for energy in (10.0e6, 50.0e6, 200.0e6):
                H_func = H.H_factory(d, params, H.vcc_of(ne), energy)
                for frac in np.linspace(0.2, 1.0, 12):
                    l1 = L0 + frac*SPAN
                    tot += 1
                    ok = (ad._profile_is_resolved(H_func, L0, l1, 200)
                          or ad._profile_is_resolved(H_func, L0, l1, 6400))
                    if not ok:
                        bad += 1
                    s = local_statistic(H_func, L0, l1, 200)
                    if s is not None:
                        loc.append(s)
        stats += loc
        n_bad_total += bad
        n_total += tot
        print('%-24s %10d %10d   %s'
              % (label, bad, tot,
                 ('n=%d  median %.3f  max %.3f' % (len(loc), np.median(loc), max(loc)))
                 if loc else 'nothing flagged'))
    print('%-24s %10d %10d' % ('TOTAL', n_bad_total, n_total))
    return stats, n_bad_total, n_total


def main():
    smooth, piecewise = families()
    print('# Resolution-test false positives, swept over SUB-INTERVALS')
    print('# 12 sub-intervals x d = 2-5 x 3 energies per family')
    print('# RESOLUTION_RATIO = %.2f  LOCAL_JUMP_RATIO = %.2f  N_LOCAL_CONFIRM = %d'
          % (ad.RESOLUTION_RATIO, ad.LOCAL_JUMP_RATIO, ad.N_LOCAL_CONFIRM))
    s_stats, s_bad, s_tot = sweep(smooth, 'SMOOTH (want 0 unresolved)')
    p_stats, p_bad, p_tot = sweep(piecewise, 'PIECEWISE (want all unresolved)')
    print('\n=== SUMMARY ===')
    print('smooth    : %d / %d declared unresolved   <-- false positives' % (s_bad, s_tot))
    print('piecewise : %d / %d declared unresolved   <-- true positives' % (p_bad, p_tot))
    if s_stats:
        print('local statistic, smooth flagged intervals   : n=%d max %.3f'
              % (len(s_stats), max(s_stats)))
    if p_stats:
        print('local statistic, piecewise flagged intervals: n=%d min %.3f'
              % (len(p_stats), min(p_stats)))
    if s_stats and p_stats:
        print('separation: smooth ceiling %.3f  <  LOCAL_JUMP_RATIO %.2f  <  jump floor %.3f'
              % (max(s_stats), ad.LOCAL_JUMP_RATIO, min(p_stats)))


if __name__ == '__main__':
    sys.exit(main())
