"""average=True on a pseudo-Dirac pair against a 50-digit reference (issue #165).

Three flavors at the NuFIT defaults, one pair at a time, E = 1 TeV.  For each ratio
delta/max m^2 and each pair phase phi = delta L / 2E, two routes are compared with the
reference over the 3x3 active block:

  generic   osc_prob_energy_baseline(lambda e: H0/e, E, L, average=True), with H0 from
            hamiltonian_pseudo_dirac_vacuum_energy_independent: the phase average of an
            eigendecomposition of H0/E;
  wrapper   osc_prob_pseudo_dirac_vacuum(E, L, pairs, average=True): the phase average of
            the known eigensystem, avgprob.phase_averaged_probabilities_pseudo_dirac.

The reference is the phase average of the exact eigenbasis: W from pseudo_dirac_mixing_matrix,
every phase (M2_i - M2_j) L / 2E with the splitting added in 50-digit arithmetic, damped by the
default spread, exp(-sigma^2 phi^2 / 2).  Reported: the largest error over the three states and
the phases.

Run from the repository root:
    python docs/dev/measurements/issue165_pseudo_dirac_resolution/resolution.py
"""
import warnings

import mpmath as mp
import numpy as np

import magnus.avgprob as ap
import magnus.globaldefs as gd
import magnus.hamiltonians as hh
import magnus.oscprob as o
from magnus.hamiltonians.hamiltonians_pseudodirac import _pair_layout

mp.mp.dps = 50
p = gd.load_nufit_params()
U = hh.pmns_mixing_matrix(p['s12'], p['s23'], p['s13'], p['dCP'])
m2 = [0.0, p['D21'], p['D31']]
E = 1e3*gd.UNIT_GEV
sigma = ap.AVG_PHASE_SPREAD


def reference(pairs, L):
    W = hh.pseudo_dirac_mixing_matrix(U, pairs)
    columns, _ = _pair_layout(3, pairs)
    M2 = [mp.mpf(m2[j]) + (mp.mpf(pairs[j]) if sign == -1 else 0) for j, sign in columns]
    n = len(columns)
    Wm = [[mp.mpc(complex(W[a, i])) for i in range(n)] for a in range(n)]
    P = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            total = mp.mpc(0)
            for i in range(n):
                for j in range(n):
                    phi = (M2[i] - M2[j])*mp.mpf(L)/(2*mp.mpf(E))
                    damp = mp.e**(-sigma**2*phi**2/2)
                    if damp < mp.mpf('1e-30'):
                        continue
                    total += (mp.conj(Wm[a][i])*Wm[b][i]*Wm[a][j]*mp.conj(Wm[b][j])
                              *mp.e**(-1j*phi)*damp)
            P[a, b] = float(mp.re(total))
    return P


print('delta/max m2   phases [rad]          worst error: generic    wrapper')
for ratio in (1e-4, 1e-6, 1e-8, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 4e-16):
    delta = ratio*max(m2)
    worst_generic = worst_wrapper = 0.0
    for phase in (0.5, 1.0, 3.0, 10.0, 20.0, 30.0):
        for state in (0, 1, 2):
            pairs = {state: delta}
            L = phase*2*E/delta
            H0 = hh.hamiltonian_pseudo_dirac_vacuum_energy_independent(U, m2, pairs)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                generic = np.asarray(o.osc_prob_energy_baseline(
                    lambda e: H0/e, E, L, 0.0, None, None, True, average=True))
                wrapper = np.asarray(o.osc_prob_pseudo_dirac_vacuum(
                    E, L, pairs, average=True))
            ref = reference(pairs, L)
            worst_generic = max(worst_generic, float(np.abs(generic[:3, :3] - ref).max()))
            worst_wrapper = max(worst_wrapper, float(np.abs(wrapper[:3, :3] - ref).max()))
    print('%-14.0e 0.5, 1, 3, 10, 20, 30   %20.1e %10.1e' % (ratio, worst_generic,
                                                            worst_wrapper))
