# -*- coding: utf-8 -*-
"""Diagnose the Battery 2.1 failures, and run Battery 7.5's premise check.

Battery 7.5: "every stretch with gamma > threshold should lie inside a reported window."
This is commit 8's premise.  Scan gamma on a grid 500x denser than the detector's, and ask:

  (a) what is gamma_max, and is it above the 0.1 threshold at all?
  (b) does every dense-grid exceedance lie inside a window the detector reports?

Distinguishes the two failure mechanisms:
  * gamma_max > threshold but no window  -> a DETECTION miss (the probe grid stepped over it)
  * gamma_max < threshold, answer wrong  -> SUB-THRESHOLD ACCUMULATION (no threshold can help)
"""

import numpy as np

import harness as H
import magnus.adiabatic as ad
from battery2 import L0, L1, E, bump_profile, ne_res_for

span = L1 - L0
p2 = H.params_for(2)
ner2 = ne_res_for(2, p2, E)
rng = np.random.default_rng(7)
lc = L0 + (0.37 + 0.2*rng.random())*span
FD = span*1e-6


def dense_gamma(H_of_l, l0, l1, d, n, focus=None):
    """gamma_01 on a dense grid; `focus` adds extra density around a known feature."""
    ls = np.linspace(l0, l1, n)
    if focus is not None:
        c, w = focus
        ls = np.unique(np.concatenate([ls, np.linspace(max(l0, c - 30*w),
                                                       min(l1, c + 30*w), n)]))
    g = np.array([ad._point_adiabaticity(H_of_l, float(l), 0, 1, FD, (l0, l1)) for l in ls])
    return ls, g


print('%-8s %11s %11s %8s %9s %-26s %s'
      % ('w/span', 'hybrid err', 'gamma_max', 'n>thr', 'covered?', 'windows(final)',
         'MECHANISM'))

for wf in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]:
    w = wf*span
    H_of_l = H.H_factory(2, p2, H.vcc_of(bump_profile(ner2, lc, w)), E)
    U, win, cert = ad.hybrid_propagator(H_of_l, L0, L1, rtol=1e-3, atol=1e-3)
    err = H.maxabs(H.P_of(U) - H.P_of(H.exact_U(H_of_l, L0, L1, 2)))

    ls, g = dense_gamma(H_of_l, L0, L1, 2, 100_000, focus=(lc, w))
    gmax = float(np.nanmax(g))
    over = ls[g > 0.1]
    n_over = int(over.size)

    if n_over == 0:
        covered = 'n/a'
    else:
        inside = np.zeros(n_over, dtype=bool)
        for (a, b) in win:
            inside |= (over >= a) & (over <= b)
        covered = 'YES' if inside.all() else 'NO (%d/%d)' % ((~inside).sum(), n_over)

    if err <= 2e-3:
        mech = 'ok'
    elif gmax > 0.1 and covered != 'YES':
        mech = 'DETECTION MISS (probe stepped over it)'
    elif gmax <= 0.1:
        mech = 'SUB-THRESHOLD ACCUMULATION'
    else:
        mech = 'covered but still wrong'

    print('%-8.0e %11.3e %11.3e %8d %9s %-26s %s'
          % (wf, err, gmax, n_over, covered,
             '[' + ', '.join('%.4f-%.4f' % (a/span, b/span) for a, b in win) + ']', mech))
