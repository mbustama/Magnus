# -*- coding: utf-8 -*-
"""Tests of the adiabatic-transport-plus-Magnus-patch hybrid strategy (magnus.adiabatic)."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import magnus.adiabatic as ad
import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.matter as matter

RNG = np.random.default_rng(11)


def maxabs(x):
    return np.max(np.abs(np.asarray(x)))


def exact_U(H_func, l0, l1, dim):
    """Ground truth evolution operator via a tight-tolerance ODE solve, matching the sign/time-
    ordering convention used throughout magnus.oscprob (dU/dl = -i H(l) U(l))."""
    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(dim, dim)).ravel()
    sol = solve_ivp(rhs, (l0, l1), np.eye(dim, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853')
    return sol.y[:, -1].reshape(dim, dim)


# ----------------------------------------------------------------------
# adiabatic_propagator
# ----------------------------------------------------------------------

def test_adiabatic_propagator_is_unitary_and_matches_exact_away_from_resonance():
    """A slowly-varying, non-crossing 2-level real Hamiltonian: pure adiabatic transport (no
    patch needed) must be exactly unitary and match a tight-tolerance ODE solve.

    Kept deliberately small in magnitude relative to l_scale (unlike a "natural units" physical
    Hamiltonian, where couplings and l_scale differ by many orders of magnitude but the *product*
    of any constant diagonal term and the domain width still stays modest): a constant diagonal
    entry as large as O(1) times a domain of O(l_scale) would add a huge, physically irrelevant
    accumulated phase that only makes the solve_ivp cross-check itself slow, not anything more
    rigorous (this is exactly the trap hit and fixed during development of this module, with a
    synthetic multi-resonance Hamiltonian whose sector-separating baseline offset was needlessly
    large)."""
    l_scale = 10.0

    def H_func(l):
        v = 0.1*np.exp(-l/l_scale)  # stays below the constant 0.3: no crossing in range
        return np.array([[v, 0.02], [0.02, 0.3]])

    l0, l1 = 0.0, 2.0*l_scale
    U = ad.adiabatic_propagator(H_func, l0, l1, n_points=201)
    assert maxabs(U.conj().T @ U - np.eye(2)) < 1e-12

    U_exact = exact_U(H_func, l0, l1, 2)
    assert maxabs(U - U_exact) < 2e-2


def test_adiabatic_propagator_zero_width_is_identity():
    def H_func(l):
        return np.eye(3)
    U = ad.adiabatic_propagator(H_func, 5.0, 5.0)
    assert np.array_equal(U, np.eye(3, dtype=complex))


def test_adiabatic_propagator_handles_real_valued_hamiltonian():
    """A genuinely real (no CP violation) Hamiltonian must not crash on dtype casting inside the
    parallel-transport phase-fixing step (regression test for a bug caught during development:
    eigh on a real array returns real eigenvectors that cannot hold a complex phase)."""
    def H_func(l):
        v = np.exp(-l/50.0)
        return np.array([[v, 0.05, 0.0], [0.05, 0.5, 0.03], [0.0, 0.03, 0.2]])
    U = ad.adiabatic_propagator(H_func, 0.0, 100.0, n_points=101)
    assert maxabs(U.conj().T @ U - np.eye(3)) < 1e-12


# ----------------------------------------------------------------------
# find_resonance_candidates / find_nonadiabatic_windows
# ----------------------------------------------------------------------

def test_find_resonance_candidates_locates_known_crossing():
    """A 2-level H(l) = diag(V(l), v1) with V(l) monotonically decreasing crosses v1 exactly
    once, at a position solvable in closed form; the Hellmann-Feynman detector must find it to
    high precision without any assumption about the profile beyond smoothness."""
    l_scale = 1.0e5
    V0 = 4.0
    v1 = 1.5
    eps = 1e-3

    def H_func(l):
        v = V0*np.exp(-l/l_scale)
        return np.array([[v, eps], [eps, v1]])

    l0, l1 = 0.0, 4.0*l_scale
    candidates = ad.find_resonance_candidates(H_func, l0, l1)
    assert len(candidates) == 1
    l_star_expected = -l_scale*np.log(v1/V0)
    assert abs(candidates[0]['l'] - l_star_expected) < 1e-3*l_scale
    assert candidates[0]['j'] == 0 and candidates[0]['k'] == 1
    assert candidates[0]['gap'] == pytest.approx(2*eps, rel=1e-6)


def test_find_nonadiabatic_windows_empty_for_strong_coupling_present_for_weak():
    """The same crossing position (same gap-extremum location), with a coupling large enough to
    open a wide avoided crossing (adiabatic, no window) versus one small enough to make the
    crossing sharp and diabatic (one window straddling it).

    This is the standard Landau-Zener direction, and it is the opposite of what "weak/strong
    coupling" might naively suggest: the adiabaticity parameter is gamma ~ coupling/gap^2, and
    the gap itself scales with the coupling (gap = 2*eps at an exact 2-level crossing), so
    gamma ~ 1/eps -- a *smaller* coupling gives a *larger* (more non-adiabatic) gamma, matching
    the real LMA solar case (large mixing angle -> adiabatic conversion) versus an artificially
    tiny mixing angle (needed, in earlier development/validation of this module, to force a
    genuine non-adiabatic crossing at all)."""
    l_scale = 1.0e5
    V0, v1 = 0.3, 0.1

    def make_H(eps):
        def H_func(l):
            v = V0*np.exp(-l/l_scale)
            return np.array([[v, eps], [eps, v1]])
        return H_func

    l0, l1 = 0.0, 4.0*l_scale
    l_star = -l_scale*np.log(v1/V0)

    windows_strong, _ = ad.find_nonadiabatic_windows(make_H(1e-2), l0, l1, threshold=0.1)
    assert windows_strong == []

    windows_weak, cands = ad.find_nonadiabatic_windows(make_H(1e-4), l0, l1, threshold=0.1)
    assert len(windows_weak) == 1
    w = windows_weak[0]
    assert w[0] < l_star < w[1]
    assert cands[0]['gamma'] > 0.1


def test_find_nonadiabatic_windows_merges_close_resonances():
    """Two distinct crossings close enough together that their individually-grown windows
    overlap must come back as a single merged window from the real detect->grow->merge pipeline
    (not a hand-specified window list)."""
    l_scale = 1.0e5
    v1 = 1.5

    def H_func(l):
        v = 4.0*np.exp(-l/l_scale)
        # 3-level chain, state 1 shared between both crossings (0-1 and 1-2), with the second
        # crossing's position controlled by scaling state 2's coupling to V(l) by alpha != 1.
        alpha = 1.005
        eps = 1e-3
        return np.array([
            [v, eps, 0.0],
            [eps, v1, eps],
            [0.0, eps, alpha*v],
        ])

    l0, l1 = 0.0, 4.0*l_scale
    windows, cands = ad.find_nonadiabatic_windows(H_func, l0, l1, threshold=0.1)
    assert len(windows) == 1
    # Both bare crossings (0,1) and (1,2) must be genuine (non-adiabatic) candidates inside it.
    pairs_in_window = {(c['j'], c['k']) for c in cands
                       if c['gamma'] > 0.1 and windows[0][0] <= c['l'] <= windows[0][1]}
    assert (0, 1) in pairs_in_window
    assert (1, 2) in pairs_in_window


def test_find_nonadiabatic_windows_keeps_far_apart_resonances_separate():
    """Two independent, well-separated 2-level sectors (no shared level, no cross-coupling) with
    distinct crossing positions must come back as two separate windows, not merged and not
    conflated into one."""
    l_scale = 1.0e5

    def H_func(l):
        v = np.exp(-l/l_scale)
        eps_a, eps_b = 3e-4, 1e-4
        v_res_a, v_res_b = 0.535, 0.0535
        base_b = 1.3
        return np.array([
            [v, eps_a, 0.0, 0.0],
            [eps_a, v_res_a, 0.0, 0.0],
            [0.0, 0.0, base_b + v, eps_b],
            [0.0, 0.0, eps_b, base_b + v_res_b],
        ])

    l0, l1 = 0.0, 4.0*l_scale
    windows, _ = ad.find_nonadiabatic_windows(H_func, l0, l1, threshold=0.1)
    assert len(windows) == 2
    assert windows[0][1] < windows[1][0]  # disjoint and correctly ordered


# ----------------------------------------------------------------------
# hybrid_propagator: unitarity, certification, and accuracy vs solve_ivp
# ----------------------------------------------------------------------

def test_hybrid_propagator_pure_adiabatic_matches_exact():
    """No resonance in range: hybrid_propagator must reduce to (and certify) pure adiabatic
    transport, matching a tight-tolerance ODE solve."""
    osc = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
        osc['s12'], osc['s23'], osc['s13'], osc['dCP'], osc['D21'], osc['D31'])
    e00 = np.diag([1.0, 0.0, 0.0])
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)
    energy = 18.0*gd.UNIT_MEV
    l_scale = gd.L_SCALE_SUN

    def H_func(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))*e00

    l0, l1 = 0.0, 4.0*l_scale
    U, windows, certified = ad.hybrid_propagator(H_func, l0, l1)
    assert certified
    assert windows == []
    assert maxabs(U.conj().T @ U - np.eye(3)) < 1e-12

    U_exact = exact_U(H_func, l0, l1, 3)
    assert maxabs(U - U_exact) < 1e-3


def test_hybrid_propagator_patches_genuine_resonance():
    """An NSI coupling engineered to induce a genuine non-adiabatic crossing: hybrid_propagator
    must detect and patch the window, stay unitary, certify, and match solve_ivp."""
    osc = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
        osc['s12'], osc['s23'], osc['s13'], osc['dCP'], osc['D21'], osc['D31'])
    h_matt = np.diag([1.0, 0.0, 0.0]) + hams.hamiltonian_3nu_nsi(1.0, 0.0, 0.0j, 3.0, 0.0, 0.0j, 0.0)
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)
    energy = 18.0*gd.UNIT_MEV
    l_scale = gd.L_SCALE_SUN

    def H_func(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))*h_matt

    l0, l1 = 2.5*l_scale, 4.5*l_scale
    U, windows, certified = ad.hybrid_propagator(H_func, l0, l1)
    assert certified
    assert len(windows) == 1
    assert maxabs(U.conj().T @ U - np.eye(3)) < 1e-10

    U_exact = exact_U(H_func, l0, l1, 3)
    assert maxabs(U - U_exact) < 2e-2


@pytest.mark.parametrize("num_flavors", [4, 5])
def test_hybrid_propagator_handles_higher_dimensions(num_flavors):
    """The detector and propagator must work unmodified on 4- and 5-level Hamiltonians (3+1/3+2
    sterile neutrinos), not just the 2/3-flavor cases used elsewhere -- this is the core
    "arbitrary Hamiltonian, arbitrary dimension" requirement."""
    osc = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)
    energy = 18.0*gd.UNIT_MEV
    l_scale = gd.L_SCALE_SUN

    if num_flavors == 4:
        hvac = hams.hamiltonian_4nu_vacuum_energy_independent(
            osc['s12'], osc['s23'], osc['s13'], osc['dCP'], 0.15, 1.2, 0.10, 0.0, 0.05,
            osc['D21'], osc['D31'], 1.5*osc['D31'])
        e00 = np.diag([1.0, 0.0, 0.0, 0.0])
    else:
        hvac = hams.hamiltonian_5nu_vacuum_energy_independent(
            osc['s12'], osc['s23'], osc['s13'], osc['dCP'], 0.15, 1.2, 0.08, 0.5, 0.10, 0.0,
            0.05, 0.05, 0.03, 0.9, osc['D21'], osc['D31'], 1.5*osc['D31'], 2.5*osc['D31'])
        e00 = np.diag([1.0, 0.0, 0.0, 0.0, 0.0])

    def H_func(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))*e00

    l0, l1 = 0.0, 4.0*l_scale
    U, windows, certified = ad.hybrid_propagator(H_func, l0, l1)
    assert certified
    assert windows == []  # standard mixing stays adiabatic at this energy, as in the 3nu case
    assert maxabs(U.conj().T @ U - np.eye(num_flavors)) < 1e-10
