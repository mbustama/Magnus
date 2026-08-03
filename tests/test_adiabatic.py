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


# ----------------------------------------------------------------------
# hybrid_propagator: the uncertified path
# ----------------------------------------------------------------------

def _sun_like_3nu_H():
    """A standard 3nu solar Hamiltonian, adiabatic over the whole range (0 windows)."""
    osc = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
        osc['s12'], osc['s23'], osc['s13'], osc['dCP'], osc['D21'], osc['D31'])
    e00 = np.diag([1.0, 0.0, 0.0])
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)
    energy = 18.0*gd.UNIT_MEV

    def H_func(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))*e00

    return H_func, 0.0, 4.0*gd.L_SCALE_SUN


def test_hybrid_propagator_does_not_certify_when_every_knob_is_saturated():
    """Regression test: with all three refinement knobs pinned at their ceilings from the very
    first iteration, successive iterations recompute bit-identical inputs.  The agreement test
    would then compare a result with itself and pass trivially, so the loop must refuse to
    certify rather than report success on the strength of a comparison carrying no information.
    The returned operator is still exactly unitary -- that never depended on certification."""
    H_func, l0, l1 = _sun_like_3nu_H()

    U, _, certified = ad.hybrid_propagator(
        H_func, l0, l1,
        rtol=1e-300, atol=1e-300,          # unreachable, so it can never certify honestly
        threshold0=0.1, min_threshold=0.1,  # already at its floor
        n_probe0=200, max_n_probe=200,      # already at its ceiling
        n_points0=201, max_n_points=201)    # already at its ceiling

    assert certified is False
    assert maxabs(U.conj().T @ U - np.eye(3)) < 1e-12


def test_hybrid_propagator_does_not_certify_when_iterations_run_out():
    """An unreachable tolerance with room left to refine must exhaust ``max_iters`` and report
    ``certified=False``, again while staying exactly unitary."""
    H_func, l0, l1 = _sun_like_3nu_H()

    U, _, certified = ad.hybrid_propagator(
        H_func, l0, l1, rtol=1e-300, atol=1e-300, max_iters=2)

    assert certified is False
    assert maxabs(U.conj().T @ U - np.eye(3)) < 1e-12


# ----------------------------------------------------------------------
# The two remaining refusals: a local Magnus patch that does not converge.
#
# hybrid_propagator already refuses to certify when its knobs saturate and
# when it runs out of iterations (both tested above). It has two further
# refusals, for when an evaluation reports that the patch across a
# non-adiabatic window did not itself converge, and a coverage run found
# neither was ever taken. Reaching them honestly means driving
# _local_evolution_operator to its 500,000-slab ceiling, so instead the
# single evaluation is replaced: what is under test here is the
# certification contract -- one non-converged patch anywhere means the
# whole result is uncertified -- not the numerics of the patch, which the
# solve_ivp cross-checks above already cover.
# ----------------------------------------------------------------------

def _flat_2level_H():
    """Any smooth H will do; these tests never reach the real numerics."""
    def H_func(l):
        return np.array([[1.0e-12, 1.0e-13], [1.0e-13, -1.0e-12]], dtype=complex)
    return H_func, 0.0, 1.0e3


def test_hybrid_propagator_does_not_certify_when_the_first_patch_fails(monkeypatch):
    """A non-converged patch on the very first evaluation is fatal at once:
    there is no point tightening the knobs and paying for another pass when
    the answer already cannot be certified."""
    calls = []

    def fake_once(H_func, l0, l1, threshold, n_probe, n_points, fd_step_frac,
                  magnus_exp_order, integration_method):
        calls.append((threshold, n_probe, n_points))
        return np.eye(2, dtype=complex), [(1.0, 2.0)], False, 0.0

    monkeypatch.setattr(ad, '_hybrid_propagator_once', fake_once)
    H_func, l0, l1 = _flat_2level_H()

    U, windows, certified = ad.hybrid_propagator(H_func, l0, l1)

    assert certified is False
    assert len(calls) == 1, "it kept refining after a result it already knew it could not certify"
    assert windows == [(1.0, 2.0)], "the windows found must still be reported back to the caller"
    assert maxabs(U.conj().T @ U - np.eye(2)) < 1e-12


def test_hybrid_propagator_does_not_certify_when_a_later_patch_fails(monkeypatch):
    """The same refusal from inside the refinement loop. The first pass
    converged, so the loop tightened its knobs and tried again; that second
    pass reports a failed patch, and one failure is enough to withhold
    certification even though a comparison was available."""
    results = [True, False]
    calls = []

    def fake_once(H_func, l0, l1, threshold, n_probe, n_points, fd_step_frac,
                  magnus_exp_order, integration_method):
        calls.append((threshold, n_probe, n_points))
        return np.eye(2, dtype=complex), [], results[len(calls) - 1], 0.0

    monkeypatch.setattr(ad, '_hybrid_propagator_once', fake_once)
    H_func, l0, l1 = _flat_2level_H()

    U, _, certified = ad.hybrid_propagator(H_func, l0, l1)

    assert certified is False, \
        "identical operators from two passes must not certify when the second reported failure"
    assert len(calls) == 2
    assert calls[1] != calls[0], "the second pass must run at genuinely tightened knobs"
    assert maxabs(U.conj().T @ U - np.eye(2)) < 1e-12


def test_windows_are_found_away_from_gap_extrema():
    """The adiabaticity parameter is not largest where the gap is stationary.

    find_resonance_candidates locates gap *extrema* -- sign changes of the Hellmann-Feynman
    derivative difference. Evaluating gamma only there understates its maximum badly on a
    rapidly varying profile: measured on this one, 3.6e-04 at the extrema against 7.0e-02 along
    the path, a factor of 196.

    The consequence was not merely a loose estimate. No window ever opened, so successive
    refinements of hybrid_propagator differed only in the adiabatic-transport grid, converged to
    the same wrong adiabatic limit, agreed with each other to 8.3e-05, and certified a result
    that was wrong by 4.3e-02 against solve_ivp -- the false-convergence failure class recorded
    in docs/dev/NOTES_ADAPTIVE_REFINEMENT.md, on the adiabatic path.

    The profile is a solar exponential modulated by a strong, fast sine, with an NSI coupling; it
    is smooth, so nothing here is about discontinuities."""
    osc = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
        osc['s12'], osc['s23'], osc['s13'], osc['dCP'], osc['D21'], osc['D31'])
    h_matt = np.diag([1.0, 0.0, 0.0]) + hams.hamiltonian_3nu_nsi(1.0, 0.0, 0.0j, 3.0,
                                                                 0.0, 0.0j, 0.0)
    LS = gd.L_SCALE_SUN
    energy = 18.0*gd.UNIT_MEV

    def rho(l):
        l = np.asarray(l)
        return gd.NUM_DENSITY_E_SUN_CENTRAL*np.exp(-l/LS)*(
            1.0 + 0.9*np.sin(2.0*np.pi*l/(0.45*LS)))

    VCC_func = matter.vcc_func_from_rho_func(rho, 0.0, 1.0, 0.5, False, False, True)

    def H_func(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))[..., None, None]*h_matt

    l0, l1 = 0.5*LS, 1.54*LS

    # The gap extrema really are quiet here -- that is the point of the case, so assert it
    # rather than trusting it, and the test degrades loudly if the profile ever changes.
    fd_step = (l1 - l0)*1e-6
    cands = ad.find_resonance_candidates(H_func, float(l0), float(l1))
    gamma_at_extrema = max(
        ad._point_adiabaticity(H_func, c['l'], c['j'], c['k'], fd_step, (float(l0), float(l1)))
        for c in cands)
    assert gamma_at_extrema < 1e-3, \
        f"gap extrema are no longer quiet on this profile (gamma {gamma_at_extrema:.2e})"

    # Probed at a threshold between the two scales: below the ~7e-2 that gamma reaches along the
    # path, far above the ~3.6e-4 it reaches at the extrema. Only a sweep of the path can open a
    # window here; extrema-only evaluation cannot, at any threshold above 3.6e-4.
    #
    # Not the default 0.1: gamma never exceeds that anywhere on this profile, so nothing opens on
    # the first iteration and it is hybrid_propagator's own threshold refinement (0.1 -> 0.0333)
    # that reaches the regime where the sweep matters.
    windows, _ = ad.find_nonadiabatic_windows(H_func, float(l0), float(l1), threshold=0.03)
    assert windows, "no window opened at threshold 0.03, though gamma reaches ~7e-2 on the path"

    U, _, certified = ad.hybrid_propagator(H_func, float(l0), float(l1))
    assert certified
    U_exact = exact_U(H_func, float(l0), float(l1), 3)
    P = np.abs(U).T**2
    P_exact = np.abs(U_exact).T**2
    err = np.max(np.abs(P - P_exact))
    assert err < 1e-3, f"certified but wrong by {err:.2e}"


# ----------------------------------------------------------------------
# Two ways the hybrid strategy used to certify a wrong answer.  Both were found by the
# adversarial-validation batteries (docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md) and both are
# pre-existing rather than introduced by the gamma sweep; see that document for the full
# reproductions and for the third, undetectable, case these do not cover.
# ----------------------------------------------------------------------

def _solar_step_H(mid_frac=0.5, l1=None):
    """2nu Hamiltonian on a density profile with one unmarked step, and its exact answer.

    H is constant on each side of the step, so scipy.linalg.expm composed across the two
    pieces is the *exact* evolution operator, not an approximation -- a stronger oracle than
    solve_ivp, and one that cannot itself step over the discontinuity.
    """
    l1 = gd.L_SCALE_SUN if l1 is None else l1
    mid = mid_frac*l1
    lo, hi = 0.02*gd.NUM_DENSITY_E_SUN_CENTRAL, 0.30*gd.NUM_DENSITY_E_SUN_CENTRAL

    def ne(l):
        y = np.where(np.asarray(l, dtype=float) < mid, lo, hi)
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a

    vcc = matter.vcc_func_from_rho_func(ne, 0.0, 1.0, 0.5, nubar=False,
                                        density_matter_is_in_g_per_cm3=False,
                                        density_is_of_number_of_electrons=True)
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    energy = 50.0e6

    def H_func(l):
        v = np.asarray(vcc(l))
        return (1.0/energy)*h_vac + v[..., None, None]*proj

    return H_func, mid, l1


def _solar_bump_H(width_frac, energy=10.0e6, centre_frac=0.5):
    """2nu solar-scale Hamiltonian with one Gaussian resonance of the given fractional width."""
    l1 = gd.L_SCALE_SUN
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    ne0 = gd.NUM_DENSITY_E_SUN_CENTRAL
    c_vcc = float(np.asarray(matter.vcc_func_from_rho_func(
        lambda l: ne0, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)(0.0)))/ne0

    def gap(ne):
        lam = np.linalg.eigvalsh(h_vac/energy + ne*c_vcc*proj)
        return lam[1] - lam[0]

    xs = np.geomspace(ne0*1e-6, ne0*10.0, 2000)
    ne_res = float(xs[int(np.argmin([gap(x) for x in xs]))])
    width = width_frac*l1

    def ne(l):
        x = np.asarray(l, dtype=float)
        y = ne_res*(0.30 + 2.70*np.exp(-0.5*((x - centre_frac*l1)/width)**2))
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a

    vcc = matter.vcc_func_from_rho_func(ne, 0.0, 1.0, 0.5, nubar=False,
                                        density_matter_is_in_g_per_cm3=False,
                                        density_is_of_number_of_electrons=True)

    def H_func(l):
        v = np.asarray(vcc(l))
        return (1.0/energy)*h_vac + v[..., None, None]*proj

    return H_func, l1


def test_unmarked_discontinuity_is_detected_and_not_certified():
    """A density step whose edge the caller did not declare.

    The dispatcher's guard against a non-smooth profile used to be "did the caller pass
    t_breakpoints", which fails open -- it declined when told about the discontinuity and
    accepted when not.  Measured then: certified=True with the probability wrong by 0.54.
    """
    H_func, mid, l1 = _solar_step_H()

    assert not ad._profile_is_resolved(H_func, 0.0, l1, 200), \
        "the resolution test no longer sees an unmarked density step"
    assert not ad._profile_is_resolved(H_func, 0.0, l1, 6400), \
        "a genuine jump must stay unresolved at every probe density"

    _, _, certified = ad.hybrid_propagator(H_func, 0.0, l1)
    assert not certified, "hybrid certified an unmarked discontinuity"


def test_a_smooth_profile_is_still_reported_as_resolved():
    """The mirror of the test above: the resolution test must not fire on ordinary profiles,
    or every solar call would abandon the hybrid strategy."""
    l1 = gd.L_SCALE_SUN
    vcc = matter.vcc_func_from_rho_func(
        matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN),
        0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)

    def H_func(l):
        v = np.asarray(vcc(l))
        return (1.0/5.0e6)*h_vac + v[..., None, None]*proj

    assert ad._profile_is_resolved(H_func, 0.0, l1, 200)
    _, _, certified = ad.hybrid_propagator(H_func, 0.0, l1)
    assert certified, "an ordinary solar profile stopped certifying"

    # A constant Hamiltonian has nothing to resolve and must not be rejected either.
    const = np.diag([1.0e-12, 2.0e-12]).astype(complex)
    assert ad._profile_is_resolved(lambda l: const, 0.0, l1, 200)


def test_a_sharp_but_smooth_feature_is_not_mistaken_for_a_discontinuity():
    """Refinement must be allowed to rescue a feature that is merely sharp at the *starting*
    probe density.

    A Gaussian of width 1e-3 of the domain is unresolved at n_probe=200 and resolved at 6400,
    and this module answers it to ~1e-11 once it is.  Testing at the starting density alone
    abandoned it as though it were a step function, which is both a large accuracy regression
    and the wrong diagnosis -- a jump stays unresolved at *every* density, which is exactly
    what separates the two.
    """
    H_func, l1 = _solar_bump_H(1.0e-3, centre_frac=0.495)

    assert not ad._profile_is_resolved(H_func, 0.0, l1, 200)
    assert ad._profile_is_resolved(H_func, 0.0, l1, 6400)

    U, windows, certified = ad.hybrid_propagator(H_func, 0.0, l1)
    assert certified, "a feature the refinement resolves was abandoned as a discontinuity"
    assert windows
    P = np.abs(U).T**2
    P_exact = np.abs(exact_U(H_func, 0.0, l1, 2)).T**2
    assert np.max(np.abs(P - P_exact)) < 1e-6


def test_subthreshold_nonadiabaticity_is_not_certified_on_agreement_alone():
    """gamma just below threshold0 everywhere: no window opens, successive iterations differ
    only in the transport grid, so they agree with each other and used to certify.

    Measured then: 1.77e-02 against a requested 1e-3, certified, silent.  The fix requires
    gamma itself to be small enough for the tolerance before an empty window list may be
    certified, so the threshold keeps dropping until a window opens.
    """
    H_func, l1 = _solar_bump_H(0.04)

    # The premise: nothing opens at the default threshold, and gamma stays under it.
    info = {}
    windows0, _ = ad.find_nonadiabatic_windows(H_func, 0.0, l1, threshold=0.1, info=info)
    assert not windows0, "profile no longer sits below the default threshold"
    assert info['gamma_max'] < 0.1, f"gamma_max {info['gamma_max']:.2e} is above threshold0"
    # The property that makes this a hard case rather than an arbitrary number: gamma_max sits
    # ABOVE the requested tolerance, so the pure adiabatic answer cannot be good enough and
    # certifying on agreement alone would be certifying something wrong.
    assert info['gamma_max'] > 1e-3, \
        f"gamma_max {info['gamma_max']:.2e} no longer exceeds the tolerance; case has no teeth"

    U, windows, certified = ad.hybrid_propagator(H_func, 0.0, l1, rtol=1e-3, atol=1e-3)
    assert certified, "the case is recoverable by refinement and should still certify"
    assert windows, "certified with no window on a profile that needed one"

    P = np.abs(U).T**2
    P_exact = np.abs(exact_U(H_func, 0.0, l1, 2)).T**2
    err = np.max(np.abs(P - P_exact))
    assert err < 1e-3, f"certified but wrong by {err:.2e}"


def test_find_nonadiabatic_windows_reports_gamma_max_via_info():
    """The out-parameter hybrid_propagator's certification rests on, pinned directly."""
    H_func, _, l1 = _solar_step_H()
    info = {}
    windows, candidates = ad.find_nonadiabatic_windows(H_func, 0.0, l1, info=info)
    assert 'gamma_max' in info
    assert info['gamma_max'] >= 0.0
    # It must dominate every candidate's own gamma, since it is a max over the probe grid too.
    for c in candidates:
        assert info['gamma_max'] >= c['gamma'] - 1e-30
    # Omitting info must remain valid (backward compatibility of the public signature).
    ad.find_nonadiabatic_windows(H_func, 0.0, l1)
