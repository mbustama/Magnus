# -*- coding: utf-8 -*-
"""Family-aware Hamiltonian builders: std / NSI / LIV, at any d, mirroring exactly what
osc_prob_matter_std_potential, osc_prob_matter_nsi and osc_prob_liv assemble internally.

Every closure returned takes EXACTLY ONE parameter (the trap in the handover)."""

import numpy as np

import harness as H
import magnus.hamiltonians as hams


def nsi_params_for(d, scale=1.0):
    if d == 2:
        return {'eps_aa': 0.05*scale, 'eps_ab': 0.02*scale}
    if d == 3:
        return {'eps_ee': 0.05*scale, 'eps_em': 0.02*scale, 'eps_et': 0.03*scale,
                'eps_mm': 0.01*scale, 'eps_mt': 0.02*scale, 'eps_tt': 0.0}
    if d == 4:
        return {'eps_ee': 0.05*scale, 'eps_em': 0.02*scale, 'eps_et': 0.03*scale,
                'eps_es': 0.04*scale, 'eps_mm': 0.01*scale, 'eps_mt': 0.02*scale,
                'eps_ms': 0.01*scale, 'eps_tt': 0.0, 'eps_ts': 0.015*scale,
                'eps_ss': 0.02*scale}
    if d == 5:
        return {'eps_ee': 0.05*scale, 'eps_em': 0.02*scale, 'eps_et': 0.03*scale,
                'eps_es1': 0.04*scale, 'eps_es2': 0.03*scale, 'eps_mm': 0.01*scale,
                'eps_mt': 0.02*scale, 'eps_ms1': 0.01*scale, 'eps_ms2': 0.012*scale,
                'eps_tt': 0.0, 'eps_ts1': 0.015*scale, 'eps_ts2': 0.011*scale,
                'eps_s1s1': 0.02*scale, 'eps_s1s2': 0.008*scale, 'eps_s2s2': 0.017*scale}
    raise ValueError(d)


def liv_params_for(d, n_liv=0, frac=0.10, energy=1.0e7):
    """LIV eigenvalues at `frac` of the vacuum splitting at `energy`, as the decision doc used."""
    p = H.params_for(3)
    scale = frac*p['D31']/(2.0*energy)
    Lam = 1.0e9
    common = {'Lambda': Lam, 'n_liv': n_liv, 'b1': 0.0, 'b2': scale, 'b3': 2.0*scale}
    if d == 2:
        return {'sxi': 0.3, 'b1': 0.0, 'b2': scale, 'Lambda': Lam, 'n_liv': n_liv}
    if d == 3:
        return {'sxi12': 0.3, 'sxi23': 0.2, 'sxi13': 0.1, 'dxiCP': 0.5, **common}
    if d == 4:
        return {'sxi12': 0.3, 'sxi23': 0.2, 'sxi13': 0.1, 'dxi13': 0.5,
                'sxi14': 0.15, 'dxi14': 0.3, 'sxi24': 0.12, 'dxi24': 0.7, 'sxi34': 0.1,
                'b4': 3.0*scale, **common}
    if d == 5:
        return {'sxi12': 0.3, 'sxi23': 0.2, 'sxi13': 0.1, 'dxi13': 0.5,
                'sxi14': 0.15, 'dxi14': 0.3, 'sxi15': 0.13, 'dxi15': 1.1,
                'sxi24': 0.12, 'dxi24': 0.7, 'sxi25': 0.1, 'sxi34': 0.1,
                'sxi35': 0.09, 'dxi35': 0.2,
                'b4': 3.0*scale, 'b5': 4.0*scale, **common}
    raise ValueError(d)


def _h_nsi_matt(d, nsi, nubar=False):
    """standard + NSI matter matrix at VCC = 1, exactly as osc_prob_matter_nsi builds it."""
    e = nsi
    if d == 2:
        m = np.diag([1.0, 0.0]) + hams.hamiltonian_2nu_nsi(1.0, e['eps_aa'], e['eps_ab'])
    elif d == 3:
        m = np.diag([1.0, 0.0, 0.0]) + hams.hamiltonian_3nu_nsi(
            1.0, e['eps_ee'], e['eps_em'], e['eps_et'], e['eps_mm'], e['eps_mt'], e['eps_tt'])
    elif d == 4:
        m = np.diag([1.0, 0.0, 0.0, 0.0]) + hams.hamiltonian_4nu_nsi(
            1.0, e['eps_ee'], e['eps_em'], e['eps_et'], e['eps_es'], e['eps_mm'], e['eps_mt'],
            e['eps_ms'], e['eps_tt'], e['eps_ts'], e['eps_ss'])
    elif d == 5:
        m = np.diag([1.0, 0.0, 0.0, 0.0, 0.0]) + hams.hamiltonian_5nu_nsi(
            1.0, e['eps_ee'], e['eps_em'], e['eps_et'], e['eps_es1'], e['eps_es2'], e['eps_mm'],
            e['eps_mt'], e['eps_ms1'], e['eps_ms2'], e['eps_tt'], e['eps_ts1'], e['eps_ts2'],
            e['eps_s1s1'], e['eps_s1s2'], e['eps_s2s2'])
    else:
        raise ValueError(d)
    return np.conj(m) if nubar else m


def _h_liv(d, liv, nubar=False):
    p = liv
    if d == 2:
        # The 2nu builder takes no `nubar` -- matching osc_prob_liv, which also omits it there.
        return hams.hamiltonian_2nu_liv_energy_independent(
            p['sxi'], p['b1'], p['b2'], p['Lambda'], p['n_liv'])
    if d == 3:
        return hams.hamiltonian_3nu_liv_energy_independent(
            p['sxi12'], p['sxi23'], p['sxi13'], p['dxiCP'], p['b1'], p['b2'], p['b3'],
            p['Lambda'], p['n_liv'], nubar=nubar)
    if d == 4:
        return hams.hamiltonian_4nu_liv_energy_independent(
            p['sxi12'], p['sxi23'], p['sxi13'], p['dxi13'], p['sxi14'], p['dxi14'], p['sxi24'],
            p['dxi24'], p['sxi34'], p['b1'], p['b2'], p['b3'], p['b4'], p['Lambda'], p['n_liv'],
            nubar=nubar)
    if d == 5:
        return hams.hamiltonian_5nu_liv_energy_independent(
            p['sxi12'], p['sxi23'], p['sxi13'], p['dxi13'], p['sxi14'], p['dxi14'], p['sxi15'],
            p['dxi15'], p['sxi24'], p['dxi24'], p['sxi25'], p['sxi34'], p['sxi35'], p['dxi35'],
            p['b1'], p['b2'], p['b3'], p['b4'], p['b5'], p['Lambda'], p['n_liv'], nubar=nubar)
    raise ValueError(d)


def H_family(family, d, energy, vcc_func, osc_params, nsi=None, liv=None, nubar=False):
    """ONE-argument H_of_l closure reproducing the wrapper's htot(E, l) for the family."""
    h_vac = np.asarray(H.h_vac_for(d, osc_params, nubar=nubar), dtype=complex)
    inv_E = 1.0/float(energy)

    if family == 'std':
        m = np.zeros((d, d), dtype=complex)
        m[0, 0] = 1.0
        h_extra = None
    elif family == 'nsi':
        m = np.asarray(_h_nsi_matt(d, nsi, nubar=nubar), dtype=complex)
        h_extra = None
    elif family == 'liv':
        m = np.zeros((d, d), dtype=complex)
        m[0, 0] = 1.0
        h_extra = np.asarray(_h_liv(d, liv, nubar=nubar), dtype=complex) \
            * float(energy)**liv['n_liv']
    else:
        raise ValueError(family)

    base = inv_E*h_vac if h_extra is None else inv_E*h_vac + h_extra

    def H_of_l(l):
        vcc = np.asarray(vcc_func(l))
        return base + vcc[..., None, None]*m

    return H_of_l
