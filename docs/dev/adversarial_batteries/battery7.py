# -*- coding: utf-8 -*-
"""Battery 7 -- cross-module and oracle diversity.

Not solve_ivp.  Each of these is an INDEPENDENT check of a different property:
  7.1 expm for constant H, d = 2..5      -- time-ordering and indexing, not accuracy
  7.2 the analytic two-flavour vacuum formula
  7.3 composition: U(0->L2) vs U(L1->L2) . U(0->L1) computed independently
  7.4 avgprob against a brute-force window average of the oscillating probability
  7.5 windows vs a dense independent gamma scan -- commit 8's premise (see diag_gamma.py
      and battery6 sub3 for the aimed versions; here on the ordinary solar profile)
"""

import warnings

import numpy as np
from scipy.linalg import expm

import harness as H
import harness6 as H6
import magnus.adiabatic as ad
import magnus.oscprob as oscprob

warnings.simplefilter('ignore')
E = 50e6
LM = 0.5*H.R_SUN
NE = H.solar_ne()

print('# Battery 7 -- oracle diversity\n')

# ---------------------------------------------------------------- 7.1
print('## 7.1  expm for a CONSTANT Hamiltonian, d = 2..5 (time-ordering + indexing)')
for d in (2, 3, 4, 5):
    p = H.params_for(d)
    ne_const = 0.05*H.NE0
    Hm = np.asarray(H6.H_family('std', d, E, H.vcc_of(lambda l: ne_const), p)(0.0))
    Ls = np.linspace(0.02*LM, LM, 60)
    P = np.asarray(oscprob.osc_prob_matter_std_potential(
        d, ne_const, E, Ls, p, L0=0.0,
        density_is_of_number_of_electrons=True)).reshape(60, d, d)
    Pref = np.array([H.P_of(expm(-1j*Hm*L)) for L in Ls])
    print('  d=%d  max|P - expm| over 60 baselines = %.3e' % (d, H.maxabs(P - Pref)))

# and through the cumulative scan explicitly (a position-dependent H is required for
# 'auto', so force it by calling the scan directly on a constant H)
for d in (2, 3):
    p = H.params_for(d)
    ne_const = 0.05*H.NE0
    Hm = np.asarray(H6.H_family('std', d, E, H.vcc_of(lambda l: ne_const), p)(0.0))
    Ls = np.linspace(0.02*LM, LM, 60)
    Pc = oscprob._osc_prob_cumulative_scan(Hm, Ls, 0.0, 400, 6, 100, 'gl', None, None)
    Pref = np.array([H.P_of(expm(-1j*Hm*L)) for L in Ls])
    print('  d=%d  cumulative scan directly vs expm            = %.3e'
          % (d, H.maxabs(Pc - Pref)))
print()

# ---------------------------------------------------------------- 7.2
print('## 7.2  analytic two-flavour VACUUM formula')
p2 = H.params_for(2)
sth, Dm2 = p2['sth'], p2['Dm2']
s2th2 = (2.0*sth*np.sqrt(1.0 - sth**2))**2
Ls = np.linspace(1e12, 1e15, 200)
P = np.asarray(oscprob.osc_prob_matter_std_potential(
    2, 0.0, E, Ls, p2, L0=0.0, density_is_of_number_of_electrons=True)).reshape(200, 2, 2)
P_analytic = 1.0 - s2th2*np.sin(Dm2*Ls/(4.0*E))**2
print('  max|P_ee - analytic| over 200 baselines = %.3e'
      % H.maxabs(P[:, 0, 0] - P_analytic))
print()

# ---------------------------------------------------------------- 7.3
print('## 7.3  composition law: U(0->L2) vs U(L1->L2) . U(0->L1)')
for d in (2, 3, 4):
    p = H.params_for(d)
    Hf = H6.H_family('std', d, E, H.vcc_of(NE), p)
    L1, L2 = 0.2*LM, 0.9*LM
    Pdirect = np.asarray(oscprob.osc_prob_matter_std_potential(
        d, NE, E, np.array([L2, 1.01*L2]), p, L0=0.0,
        density_is_of_number_of_electrons=True)).reshape(2, d, d)[0]
    # each leg computed independently, composed via the evolution operators
    Ua = H.exact_U(Hf, 0.0, L1, d)
    Ub = H.exact_U(Hf, L1, L2, d)
    Pcomp = H.P_of(Ub @ Ua)
    Pfull = H.P_of(H.exact_U(Hf, 0.0, L2, d))
    print('  d=%d  |U(L1->L2)U(0->L1) - U(0->L2)| = %.3e ;  package vs composed = %.3e'
          % (d, H.maxabs(Pcomp - Pfull), H.maxabs(Pdirect - Pcomp)))
print()

# ---------------------------------------------------------------- 7.4
print('## 7.4  average=True vs a brute-force window average of the oscillating probability')
p2 = H.params_for(2)
L_c = 0.9*H.R_SUN
Pavg = np.asarray(oscprob.osc_prob_matter_std_potential(
    2, NE, 10e6, L_c, p2, L0=0.0, density_is_of_number_of_electrons=True, average=True))
# brute force: average P over a window of ~25 oscillations near L_c
Dm2 = p2['Dm2']
osc_len = 4.0*np.pi*10e6/Dm2
Lw = np.linspace(L_c - 12.5*osc_len, L_c + 12.5*osc_len, 4001)
Hf = H6.H_family('std', 2, 10e6, H.vcc_of(NE), p2)
Uw = H.exact_U_many(Hf, 0.0, Lw, 2)
Pw = np.array([H.P_of(U)[0][0] for U in Uw])
print('  average=True P_ee   = %.6f' % Pavg[0][0])
print('  brute-force <P_ee>  = %.6f  (25 oscillations, 4001 samples, solve_ivp)' % Pw.mean())
print('  |difference|        = %.3e' % abs(Pavg[0][0] - Pw.mean()))
print()

# ---------------------------------------------------------------- 7.5
print('## 7.5  reported windows vs a dense independent gamma scan (commit 8 premise)')
for name, ne, l1 in (('solar exponential', NE, 1.0*H.L_SCALE),
                     ('multi-resonance', H.modulated_ne(amp=0.9, n_cycles=6.0,
                                                        span=1.0*H.L_SCALE),
                      1.0*H.L_SCALE)):
    for d in (2, 3):
        p = H.params_for(d)
        Hf = H6.H_family('std', d, E, H.vcc_of(ne), p)
        win, _ = ad.find_nonadiabatic_windows(Hf, 0.0, l1, threshold=0.1, n_probe=200)
        fd = l1*1e-6
        ls = np.linspace(0.0, l1, 30000)
        tot = unc = 0
        gmax = 0.0
        for j in range(d):
            for k in range(j + 1, d):
                g = np.array([ad._point_adiabaticity(Hf, float(x), j, k, fd, (0.0, l1))
                              for x in ls])
                gf = g[np.isfinite(g)]
                gmax = max(gmax, float(gf.max()) if gf.size else 0.0)
                over = ls[g > 0.1]
                tot += over.size
                if over.size:
                    ins = np.zeros(over.size, dtype=bool)
                    for (a, b) in win:
                        ins |= (over >= a) & (over <= b)
                    unc += int((~ins).sum())
        print('  %-18s d=%d  windows=%2d  gamma_max=%.3e  dense pts>thr=%6d  '
              'UNCOVERED=%6d %s'
              % (name, d, len(win), gmax, tot, unc, '' if unc == 0 else '<-- 7.5 FAILURE'))
print()
