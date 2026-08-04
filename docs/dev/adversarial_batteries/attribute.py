# -*- coding: utf-8 -*-
"""Attribution: are the Battery 2 / 5.3 failures PRE-EXISTING on main, or introduced by
commit 9c7945a (the gamma sweep)?

Run the same failing configurations under whichever magnus is on sys.path and dump the
probabilities plus the window count and certification.  Diff the two runs.
"""

import sys
import warnings

import numpy as np

warnings.simplefilter('ignore')

import magnus.adiabatic as ad             # noqa: E402
import magnus.globaldefs as gd            # noqa: E402
import magnus.matter as matter            # noqa: E402
import magnus.oscprob as oscprob          # noqa: E402

L_SCALE = gd.L_SCALE_SUN
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL
_p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
P2 = {'sth': _p['s12'], 'Dm2': _p['D21']}
E = 10.0e6
L0, L1 = 0.0, 1.0*L_SCALE
span = L1 - L0

_C = float(np.asarray(matter.vcc_func_from_rho_func(
    lambda l: NE0, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
    density_is_of_number_of_electrons=True)(0.0)))/NE0
h_vac = np.asarray(oscprob.hamiltonians.hamiltonian_2nu_vacuum_energy_independent(
    P2['sth'], P2['Dm2']), dtype=complex)
proj = np.diag([1.0, 0.0]).astype(complex)


def ne_res():
    def gap(ne):
        lam = np.linalg.eigvalsh(h_vac/E + ne*_C*proj)
        return lam[1] - lam[0]
    xs = np.geomspace(NE0*1e-6, NE0*10, 4000)
    i = int(np.argmin([gap(x) for x in xs]))
    a, b = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    for _ in range(200):
        m1, m2 = a + (b - a)/3, b - (b - a)/3
        if gap(m1) < gap(m2):
            b = m2
        else:
            a = m1
    return 0.5*(a + b)


NER = ne_res()
rng = np.random.default_rng(7)
lc = L0 + (0.37 + 0.2*rng.random())*span


def bump(w):
    def ne(l):
        x = np.asarray(l, dtype=float)
        y = NER*(0.30 + 2.70*np.exp(-0.5*((x - lc)/w)**2))
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a
    return ne


def sine(per):
    def ne(l):
        y = NER*(1.0 + 0.9*np.sin(2*np.pi*np.asarray(l, float)/per))
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a
    return ne


def step():
    mid = 0.5*L1
    lo, hi = 0.02*NE0, 0.30*NE0

    def ne(l):
        y = np.where(np.asarray(l, float) < mid, lo, hi)
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a
    return ne


CASES = {
    'B2.1 w=3e-2 span (sub-threshold)': bump(3e-2*span),
    'B2.1 w=1e-2 span (sub-threshold)': bump(1e-2*span),
    'B2.1 w=3e-5 span (detection miss)': bump(3e-5*span),
    'B2.1 w=1e-5 span (detection miss)': bump(1e-5*span),
    'B2.2 sinusoid period span/7': sine(span/7.0),
    'B5.3 step function, edge unmarked': step(),
}

out = {}
for name, ne in CASES.items():
    def H_of_l(l, _ne=ne):
        vcc = np.asarray(matter.vcc_func_from_rho_func(
            _ne, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
            density_is_of_number_of_electrons=True)(l))
        return (1.0/E)*h_vac + vcc[..., None, None]*proj

    U, win, cert = ad.hybrid_propagator(H_of_l, L0, L1, rtol=1e-3, atol=1e-3)
    Pw = np.asarray(oscprob.osc_prob_matter_std_potential(
        2, ne, E, L1, P2, L0=L0, density_is_of_number_of_electrons=True))
    out[name] = np.array([np.transpose(U.real**2 + U.imag**2)[0][0], float(len(win)),
                          float(cert), np.asarray(Pw)[0][0]])

np.savez(sys.argv[1], **out)
print('wrote', sys.argv[1])
