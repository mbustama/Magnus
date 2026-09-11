# -*- coding: utf-8 -*-
r"""Does Figure 12's cost peak sit where the crossing stops being sudden?

RESULT, measured 2026-09-11 on the 23-width scan: yes, and to the width.  gamma^-1 runs
from 728 at a 0.07 km front to 0.076 at 701 km and passes through 1 at 62.1 km, which is
where the closed form's cost peaks.  The prediction uses the Hamiltonian alone -- no
solver, no fitting, nothing from the timings -- so the agreement is a test the explanation
passed and not a story fitted to the plot.

The explanation on offer is that a front much thinner than the local oscillation length is
crossed suddenly, a front much wider is crossed adiabatically, and the expensive case is
the crossover.  That is testable without running either code.  For each front width, take
the standard non-adiabaticity of the instantaneous eigenbasis,

    gamma^-1 = max over r in the front, over pairs i!=j, of
               |<i| dH/dr |j>| / (lambda_i - lambda_r_j)^2,

which is >> 1 for a sudden crossing and << 1 for an adiabatic one.  If the explanation is
right, gamma^-1 passes through 1 at about the width where the measured cost peaks.  If it
crosses somewhere else, the explanation is wrong and I should not publish it.
"""
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'src'))

import gen_shock_cost as G                                   # noqa: E402

NS = G.NS
KM, R_F = NS['KM'], NS['R_FORWARD_KM']

def non_adiabaticity(width):
    """Max |<i|dH/dr|j>| / (dlambda)^2 across the forward front, over eigenvalue pairs."""
    w_km = width*(NS['R1_KM'] - NS['R0_KM'])
    H = NS['make_H'](NS['sn_shock_ne'](width))
    r = np.linspace(R_F - 0.75*w_km, R_F + 0.75*w_km, 401)          # km
    dr = (r[1] - r[0])*KM                                            # eV^-1
    Hs = np.asarray(H(r*KM), dtype=complex)
    worst = 0.0
    for k in range(1, len(r) - 1):
        lam, V = np.linalg.eigh(Hs[k])
        dH = (Hs[k + 1] - Hs[k - 1])/(2.0*dr)
        M = V.conj().T @ dH @ V
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                gap = lam[i] - lam[j]
                worst = max(worst, abs(M[i, j])/gap**2)
    return worst

d = json.loads((HERE/'external_shock_cost.json').read_text())
rows = sorted((c for c in d['cases']
               if c['code'] == 'NuOscProbExact' and c['target'] == 1.0e-7),
              key=lambda c: c['width'])
peak = max(rows, key=lambda c: c['seconds'])
print('measured cost peak for the closed form: %.2f km\n' % peak['width_km'], flush=True)
print('  %9s %12s %12s' % ('width_km', 'gamma^-1', 'ms/probability'), flush=True)
prev = None
for c in rows:
    g = non_adiabaticity(c['width'])
    ms = 1.0e3*c['seconds']/d['n_baselines']
    mark = ''
    if prev is not None and (prev - 1.0)*(g - 1.0) < 0:
        mark = '   <-- gamma^-1 crosses 1 here'
    prev = g
    print('  %9.3f %12.3e %12.3f%s' % (c['width_km'], g, ms, mark), flush=True)
