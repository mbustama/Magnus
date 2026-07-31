# -*- coding: utf-8 -*-
"""Tests of the oscillation-probability engine (magnus.oscprob)."""

import warnings

import numpy as np
import pytest
import scipy as sp
from scipy.integrate import solve_ivp

import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.magnus as mg
import magnus.matter as matter
import magnus.oscprob as op
import magnus.oscprobstd as opstd

RNG = np.random.default_rng(7)

# Oscillation parameters used throughout (dCP far from 0 and pi so that the
# Hamiltonian is genuinely complex)
S12, S23, S13, DCP = 0.55, 0.68, 0.15, 3.7
D21, D31 = 7.5e-5, 2.5e-3
ENERGY = 1.0*gd.UNIT_GEV
BASELINE = 1000.0*gd.UNIT_KM


def random_hermitian(dim, rng=RNG):
    X = rng.standard_normal((dim, dim)) + 1j*rng.standard_normal((dim, dim))
    return 0.5*(X + X.conj().T)


def maxabs(x):
    return np.max(np.abs(np.asarray(x)))


# ----------------------------------------------------------------------
# Closed-form cross-checks
# ----------------------------------------------------------------------

def test_3nu_vacuum_matches_closed_form():
    U_pmns = hams.pmns_mixing_matrix(S12, S23, S13, DCP)
    for nubar in [False, True]:
        P = op.osc_prob_3nu_vacuum(ENERGY, BASELINE, s12=S12, s23=S23,
                                   s13=S13, dCP=DCP, D21=D21, D31=D31,
                                   nubar=nubar)
        P_std = opstd.osc_prob_3nu_vacuum_std(U_pmns, D21, D31, ENERGY,
                                              BASELINE, nubar=nubar)
        assert maxabs(P - P_std) < 1e-12


def test_2nu_vacuum_matches_closed_form():
    sth, Dm2 = 0.4, 2.5e-3
    P = op.osc_prob_2nu_vacuum(ENERGY, BASELINE, sth, Dm2)
    P_std = opstd.osc_prob_2nu_vacuum_std(sth, Dm2, ENERGY, BASELINE)
    assert maxabs(np.asarray(P) - np.asarray(P_std)) < 1e-12


def test_2nu_constant_matter_matches_closed_form_nu_and_nubar():
    """Regression test for two sign bugs: the doubled antineutrino sign of
    the matter potential, and the flipped 2nu mass-ordering convention
    (which put the MSW resonance in the wrong channel)."""
    sth, Dm2 = 0.4, 2.5e-3
    rho = 5.0  # [g cm^-3]
    energy, L = 2.0*gd.UNIT_GEV, 2000.0*gd.UNIT_KM
    ne = rho*gd.CONV_G_TO_EV/((gd.MASS_PROTON + gd.MASS_NEUTRON)/2.0)*0.5 \
        / gd.CONV_CM3_TO_INV_EV3   # [eV^3]
    VCC = np.sqrt(2.0)*gd.GF*ne    # [eV]
    for nubar, sign in [(False, +1.0), (True, -1.0)]:
        P = op.osc_prob_2nu_matter_constant_density(
            energy, L, rho*gd.UNIT_G_PER_CM3, sth, Dm2, nubar=nubar,
            validate_input=False)
        P_std = opstd.osc_prob_2nu_matter_std(sth, Dm2, sign*VCC, energy, L)
        assert maxabs(np.asarray(P) - np.asarray(P_std)) < 1e-12, \
            f"nubar={nubar}"


# ----------------------------------------------------------------------
# Slab ordering (time-ordered product)
# ----------------------------------------------------------------------

def test_slab_ordering_two_constant_slabs():
    """A neutrino crosses slab A then slab B: the exact propagator is
    U = expm(-i H_B L2) @ expm(-i H_A L1) (B leftmost). Regression test for
    the reversed product. Uses the 'gl' method, whose interior nodes never
    touch the piecewise-constant discontinuity, so the comparison is exact
    to machine precision."""
    H_A, H_B = random_hermitian(3), random_hermitian(3)
    L1 = L2 = 0.4

    def H_piecewise(l):
        if np.ndim(l) == 0:
            return H_A if l < L1 else H_B
        raise TypeError("scalar only")

    U_exact = sp.linalg.expm(-1j*H_B*L2) @ sp.linalg.expm(-1j*H_A*L1)
    U_reversed = sp.linalg.expm(-1j*H_A*L1) @ sp.linalg.expm(-1j*H_B*L2)
    P_exact = np.abs(U_exact.T)**2
    P_reversed = np.abs(U_reversed.T)**2

    P = op.osc_prob(H_piecewise, 0.0, L1 + L2,
                    t_slab_edges=[[0.0, L1], [L1, L1 + L2]],
                    magnus_exp_order=2, integration_method='gl',
                    rtol=None, atol=None, validate_input=False)

    assert maxabs(P - P_exact) < 1e-12
    # The reversed product must be clearly distinguishable (the two
    # Hamiltonians do not commute)
    assert maxabs(P_exact - P_reversed) > 1e-2


def test_asymmetric_profile_matches_ode_solution():
    """Smooth asymmetric profile with a complex Hamiltonian, against a
    high-accuracy ODE solution."""
    H0, H1 = random_hermitian(3), random_hermitian(3)
    Lc = 1.0

    def H_ramp(l):
        lr = np.asarray(l)/Lc
        return H0 + lr[..., None, None]*H1 if lr.ndim else H0 + float(lr)*H1

    def rhs(t, y):
        return (-1j*H_ramp(t) @ y.reshape(3, 3)).ravel()

    sol = solve_ivp(rhs, (0.0, Lc), np.eye(3, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853')
    P_exact = np.abs(sol.y[:, -1].reshape(3, 3).T)**2

    for method in ['trapezoid', 'gl']:
        P = op.osc_prob(H_ramp, 0.0, Lc, n_slabs=64, n_tpts_per_slab=33,
                        magnus_exp_order=4, integration_method=method,
                        rtol=None, atol=None, validate_input=False)
        assert maxabs(P - P_exact) < 1e-6, method


# ----------------------------------------------------------------------
# osc_prob interface behavior
# ----------------------------------------------------------------------

def test_user_t_slab_edges_are_all_used():
    """Regression test: with user-provided t_slab_edges, all slabs (not just
    the first) must enter the product."""
    H0, H1 = random_hermitian(3), random_hermitian(3)

    def H_ramp(l):
        lr = np.asarray(l)
        return H0 + lr[..., None, None]*H1 if lr.ndim else H0 + float(lr)*H1

    edges = [[0.0, 0.3], [0.3, 1.0]]
    P_edges = op.osc_prob(H_ramp, 0.0, 1.0, t_slab_edges=edges,
                          n_tpts_per_slab=101, magnus_exp_order=4,
                          rtol=None, atol=None, validate_input=False)
    P_auto = op.osc_prob(H_ramp, 0.0, 1.0, n_slabs=50, n_tpts_per_slab=101,
                         magnus_exp_order=4, rtol=None, atol=None,
                         validate_input=False)
    # Same interval, different slabbings: results must agree well, and in
    # particular P_edges must NOT equal the single-slab [0, 0.3] result
    assert maxabs(P_edges - P_auto) < 1e-3
    P_first_only = op.osc_prob(H_ramp, 0.0, 0.3, n_slabs=1,
                               n_tpts_per_slab=101, magnus_exp_order=4,
                               rtol=None, atol=None, validate_input=False)
    assert maxabs(P_edges - P_first_only) > 1e-2


def test_single_none_tolerance_is_accepted():
    """Regression test: rtol=None with atol set used to crash with a
    TypeError inside np.allclose."""
    H0, H1 = random_hermitian(3), random_hermitian(3)

    def H_ramp(l):
        lr = np.asarray(l)
        return H0 + lr[..., None, None]*H1 if lr.ndim else H0 + float(lr)*H1

    P = op.osc_prob(H_ramp, 0.0, 1.0, rtol=None, atol=1e-4,
                    validate_input=False)
    assert np.all(np.isfinite(P))
    P = op.osc_prob(H_ramp, 0.0, 1.0, rtol=1e-4, atol=None,
                    validate_input=False)
    assert np.all(np.isfinite(P))


def test_energy_baseline_shapes_and_channel_selection():
    sth, Dm2 = 0.3, 2.4e-3
    energies = np.array([0.5, 1.0, 2.0])*gd.UNIT_GEV
    # scalar energy, scalar L -> 2x2 matrix
    P = op.osc_prob_2nu_vacuum(ENERGY, BASELINE, sth, Dm2)
    assert np.shape(P) == (2, 2)
    # array energy, scalar L -> (3, 2, 2)
    P = op.osc_prob_2nu_vacuum(energies, BASELINE, sth, Dm2)
    assert np.shape(P) == (3, 2, 2)
    # channel selection -> (3,)
    P = op.osc_prob_2nu_vacuum(energies, BASELINE, sth, Dm2, nu_i=gd.NUE,
                               nu_f=gd.NUMU)
    assert np.shape(P) == (3,)
    # scalar + channel -> scalar float
    P = op.osc_prob_2nu_vacuum(ENERGY, BASELINE, sth, Dm2, nu_i=gd.NUE,
                               nu_f=gd.NUMU)
    assert np.ndim(P) == 0


def test_point_parallelism_matches_serial():
    # Serial and parallel scans may follow different (warm-started) adaptive
    # refinement paths, so the results agree at the level of the requested
    # tolerance (1e-3 by default), not bit for bit.
    energies = np.array([0.5, 1.0, 2.0])*gd.UNIT_GEV
    common = dict(costhz=-0.8, L=np.full(3, 2.0*6371.0*0.8)*gd.UNIT_KM,
                  nu_i=gd.NUE, nu_f=gd.NUMU, validate_input=False)
    P_serial = op.osc_prob_3nu_earth(energies, **common)
    P_parallel = op.osc_prob_3nu_earth(energies, n_jobs=2, **common)
    assert maxabs(np.asarray(P_serial) - np.asarray(P_parallel)) < 5e-3


# ----------------------------------------------------------------------
# Physics wrappers: unitarity and vectorized-Hamiltonian consistency
# ----------------------------------------------------------------------

def test_earth_probability_rows_sum_to_one():
    P = op.osc_prob_3nu_earth(1.0*gd.UNIT_GEV, costhz=-0.8,
                              L=2.0*6371.0*0.8*gd.UNIT_KM,
                              validate_input=False)
    assert np.allclose(np.sum(P, axis=1), 1.0, atol=1e-9)
    assert np.all((P >= 0.0) & (P <= 1.0))


def test_sun_probability_rows_sum_to_one():
    # 50 MeV over 0.3 R_sun: moderate accumulated phase, converges quickly
    # within the default refinement caps (see test_tolerance_cap_warns for
    # the extreme-phase behavior)
    P = op.osc_prob_2nu_sun(50.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM,
                            0.0, np.sqrt(0.308), 7.5e-5,
                            integration_method='gl', validate_input=False)
    assert np.allclose(np.sum(P, axis=1), 1.0, atol=1e-9)


def test_bsm_wrappers_run_and_are_unitary():
    P4 = op.osc_prob_4nu_vacuum(ENERGY, BASELINE, s14=0.1, d14=0.0, s24=0.1,
                                d24=0.0, s34=0.1, D41=1.0,
                                validate_input=False)
    assert np.allclose(np.sum(P4, axis=1), 1.0, atol=1e-9)
    Pnsi = op.osc_prob_3nu_matter_nsi_constant_density(
        ENERGY, BASELINE, 5.0*gd.UNIT_G_PER_CM3, eps_ee=0.1, eps_em=0.05j,
        eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0, validate_input=False)
    assert np.allclose(np.sum(np.asarray(Pnsi), axis=1), 1.0, atol=1e-9)
    Pliv = op.osc_prob_3nu_vacuum_liv(ENERGY, BASELINE, b1=gd.B1, b2=gd.B2,
                                      b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1,
                                      validate_input=False)
    assert np.allclose(np.sum(np.asarray(Pliv), axis=1), 1.0, atol=1e-9)


def _earth_hamiltonian_chain(costhz):
    """Independent construction of the standard 3nu PREM Hamiltonian,
    mirroring what osc_prob_3nu_earth builds internally."""
    import magnus.earth as earth
    import magnus.matter as matter
    params = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
        params['s12'], params['s23'], params['s13'], params['dCP'],
        params['D21'], params['D31'])
    e00 = np.diag([1.0, 0.0, 0.0])

    def H(enu, l):
        ne = matter.num_density_e_func(
            earth.earth_radial_distance_from_depth(costhz,
                np.asarray(l)/gd.UNIT_KM),
            earth.density_matter_func_prem, density_matter_is_in_g_per_cm3=True)
        vcc = np.sqrt(2.0)*gd.GF*np.asarray(ne)
        return (1.0/enu)*hvac + vcc[..., None, None]*e00

    return H


def test_earth_probability_matches_ode_solution():
    """End-to-end physics test: 3nu through PREM at 1 GeV against a
    high-accuracy direct integration of the Schrodinger equation."""
    import magnus.earth as earth
    costhz = -0.8
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    H = _earth_hamiltonian_chain(costhz)
    enu = 1.0*gd.UNIT_GEV

    def rhs(l, y):
        return (-1j*H(enu, l) @ y.reshape(3, 3)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(3, dtype=complex).ravel(),
                    rtol=1e-10, atol=1e-12, method='DOP853')
    P_ode = np.abs(sol.y[:, -1].reshape(3, 3).T)**2

    P = op.osc_prob_3nu_earth(enu, costhz=costhz, L=L, rtol=1e-5, atol=1e-5,
                              integration_method='gl', validate_input=False)
    assert maxabs(np.asarray(P) - P_ode) < 1e-3


def test_prem_breakpoints_improve_accuracy():
    """At a fixed slab count, aligning slab edges with the PREM layer
    boundaries must reduce the error (the density is discontinuous there,
    which otherwise spoils the high-order convergence of the quadrature)."""
    import magnus.earth as earth
    costhz = -0.9
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    H = _earth_hamiltonian_chain(costhz)
    enu = 2.0*gd.UNIT_GEV

    def rhs(l, y):
        return (-1j*H(enu, l) @ y.reshape(3, 3)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(3, dtype=complex).ravel(),
                    rtol=1e-10, atol=1e-12, method='DOP853')
    P_ode = np.abs(sol.y[:, -1].reshape(3, 3).T)**2

    bp = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM
    common = dict(n_slabs=24, magnus_exp_order=4, integration_method='gl',
                  rtol=None, atol=None, validate_input=False)
    P_plain = op.osc_prob(lambda l: H(enu, l), 0.0, L, **common)
    P_bp = op.osc_prob(lambda l: H(enu, l), 0.0, L, t_breakpoints=bp, **common)
    err_plain = maxabs(P_plain - P_ode)
    err_bp = maxabs(P_bp - P_ode)
    assert err_bp < err_plain
    assert err_bp < 1e-4


@pytest.mark.parametrize("family", ['std', 'nsi', 'liv'])
def test_energy_batched_scan_matches_per_point(family):
    """The energy-batched scan engine must reproduce the per-point path.

    With fixed refinement parameters (rtol=atol=None) both paths use
    identical grids, so they must agree to near machine precision; single
    energies (nE=1) always take the per-point path, which is how the
    reference is generated here."""
    energies = np.array([0.8, 1.5, 3.0])*gd.UNIT_GEV
    L = 2.0*6371.0*0.7*gd.UNIT_KM
    common = dict(costhz=-0.7, L=L, integration_method='gl',
                  magnus_exp_order=4, rtol=None, atol=None, n_slabs=32,
                  validate_input=False)
    if family == 'std':
        f = lambda e: op.osc_prob_3nu_earth(e, **common)
    elif family == 'nsi':
        f = lambda e: op.osc_prob_3nu_earth_nsi(
            e, eps_ee=0.1, eps_em=0.05j, eps_et=0.0, eps_mm=0.0, eps_mt=0.02,
            eps_tt=0.0, **common)
    else:
        f = lambda e: op.osc_prob_3nu_earth_liv(
            e, sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0, b1=gd.B1,
            b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1, **common)
    P_scan = np.asarray(f(energies))              # batched engine (nE > 1)
    P_pts = np.array([f(float(E)) for E in energies])  # per-point path
    assert P_scan.shape == (3, 3, 3)
    assert maxabs(P_scan - P_pts) < 1e-12
    assert np.allclose(np.sum(P_scan, axis=2), 1.0, atol=1e-9)


def test_energy_batched_scan_adaptive_matches_per_point():
    """With adaptive tolerances the two paths refine differently but must
    agree at the tolerance level, and channel selection must work."""
    energies = np.linspace(0.6, 4.0, 8)*gd.UNIT_GEV
    common = dict(costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM, nu_i=gd.NUE,
                  nu_f=gd.NUMU, integration_method='gl', validate_input=False)
    P_scan = np.asarray(op.osc_prob_3nu_earth(energies, **common))
    P_pts = np.array([op.osc_prob_3nu_earth(float(E), **common)
                      for E in energies])
    assert P_scan.shape == (8,)
    assert maxabs(P_scan - P_pts) < 5e-3


def test_generic_osc_prob_earth_matches_wrapper():
    """osc_prob_earth with a hand-written standard 3nu Hamiltonian must
    reproduce the dedicated osc_prob_3nu_earth wrapper."""
    params = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    for nubar in [False, True]:
        hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
            params['s12'], params['s23'], params['s13'], params['dCP'],
            params['D21'], params['D31'], nubar=nubar)
        e00 = np.diag([1.0, 0.0, 0.0])

        def H(energy, l, VCC):
            vcc = np.asarray(VCC)
            return (1.0/energy)*hvac + vcc[..., None, None]*e00

        common = dict(costhz=-0.6, L=2.0*6371.0*0.6*gd.UNIT_KM,
                      validate_input=False)
        P_gen = op.osc_prob_earth(H, 2.0*gd.UNIT_GEV, nubar=nubar, **common)
        P_wrap = op.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, nubar=nubar, **common)
        assert maxabs(np.asarray(P_gen) - np.asarray(P_wrap)) < 3e-3
        assert np.allclose(np.sum(P_gen, axis=1), 1.0, atol=1e-9)


def test_generic_osc_prob_earth_named_locations_and_two_arg_H():
    """Named locations resolve, and a two-argument H (no matter potential)
    reproduces the vacuum probability over the chord."""
    params = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(
        params['s12'], params['s23'], params['s13'], params['dCP'],
        params['D21'], params['D31'])

    def H_vacuum_only(energy, l):
        return (1.0/energy)*hvac

    import magnus.earth as earth
    lat1 = earth.loc_coords_dms['fermilab']['lat']
    lon1 = earth.loc_coords_dms['fermilab']['lon']
    lat2 = earth.loc_coords_dms['homestake']['lat']
    lon2 = earth.loc_coords_dms['homestake']['lon']
    Lchord = earth.chord_length_inside_earth(lat1, lon1, lat2, lon2)*gd.UNIT_KM

    P = op.osc_prob_earth(H_vacuum_only, 3.0*gd.UNIT_GEV, loc_ini='fermilab',
                          loc_fin='homestake', validate_input=False)
    P_vac = op.osc_prob_3nu_vacuum(3.0*gd.UNIT_GEV, Lchord)
    assert maxabs(np.asarray(P) - np.asarray(P_vac)) < 1e-6


def test_generic_osc_prob_sun_matches_wrapper():
    """osc_prob_sun with a hand-written standard 2nu Hamiltonian must
    reproduce the dedicated osc_prob_2nu_sun wrapper.

    Uses a moderate accumulated phase (50 MeV over 0.3 R_sun, ~800 rad)
    so that both paths converge within the default refinement caps; at,
    e.g., 10 MeV over 0.9 R_sun (~1.2e4 rad of phase), the default
    max_n_slabs is insufficient and osc_prob warns
    (ToleranceNotAchievedWarning) instead of failing silently."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    e00 = np.diag([1.0, 0.0])

    def H(energy, l, VCC):
        vcc = np.asarray(VCC)
        return (1.0/energy)*hvac + vcc[..., None, None]*e00

    energy, L = 50.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM
    common = dict(integration_method='gl', rtol=1e-4, atol=1e-4,
                  validate_input=False)
    P_gen = op.osc_prob_sun(H, energy, L, **common)
    P_wrap = op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2,
                                 magnus_exp_order=4, **common)
    assert maxabs(np.asarray(P_gen) - np.asarray(P_wrap)) < 1e-3
    assert np.allclose(np.sum(P_gen, axis=1), 1.0, atol=1e-9)


def test_generic_osc_prob_sun_hybrid_strategy_resolves_hard_case():
    """The fully generic osc_prob_sun (arbitrary user H_func, no separable
    vacuum/matter decomposition available to the dispatcher) must also pick up the hybrid
    strategy: at 10 MeV over 0.9 R_sun (the same hard case as test_tolerance_cap_warns /
    test_sun_2nu_default_strategy_avoids_tolerance_cap), strategy='magnus' should still hit
    ToleranceNotAchievedWarning with tight caps, while the default strategy='auto' resolves
    warning-free and matches solve_ivp -- confirming _osc_prob_hybrid_dispatch_generic (used by
    osc_prob_sun/osc_prob_earth) works independently of _osc_prob_hybrid_dispatch (used by the
    standard/NSI/LIV wrappers)."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    e00 = np.diag([1.0, 0.0])

    def H(energy, l, VCC):
        vcc = np.asarray(VCC)
        return (1.0/energy)*hvac + vcc[..., None, None]*e00

    energy, L = 10.0*gd.UNIT_MEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM

    with pytest.warns(op.ToleranceNotAchievedWarning):
        op.osc_prob_sun(H, energy, L, integration_method='gl', max_n_slabs=64,
                        max_num_loops=8, strategy='magnus', validate_input=False)

    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        P_auto = op.osc_prob_sun(H, energy, L, strategy='auto', validate_input=False)
    assert not any(issubclass(w.category, (op.ToleranceNotAchievedWarning,
                                           op.HybridCertificationWarning,
                                           mg.MagnusConvergenceWarning)) for w in wlist)
    assert np.allclose(np.sum(P_auto, axis=1), 1.0, atol=1e-9)

    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)

    def rhs(l, y):
        Hl = (1.0/energy)*hvac + VCC_func(l)*e00
        return (-1j*Hl @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-11, atol=1e-13, method='DOP853')
    P_exact = np.abs(sol.y[:, -1].reshape(2, 2)).T**2
    assert maxabs(np.asarray(P_auto) - P_exact) < 1e-2


def test_generic_osc_prob_earth_strategy_falls_back_to_magnus():
    """osc_prob_earth's PREM density profile always carries layer-boundary breakpoints, which the
    hybrid strategy does not support (see docs/source/adiabatic_strategy.rst): strategy='auto'
    and strategy='hybrid' must therefore give results identical to strategy='magnus' for a real
    Earth-crossing trajectory -- confirming that wiring strategy into osc_prob_earth introduced
    no behavior change, only a (here, inert) new code path."""
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(S12, S23, S13, DCP, D21, D31)
    e00 = np.diag([1.0, 0.0, 0.0])

    def H(energy, l, VCC):
        vcc = np.asarray(VCC)
        return (1.0/energy)*hvac + vcc[..., None, None]*e00

    common = dict(loc_ini='fermilab', loc_fin='homestake', validate_input=False)
    P_magnus = op.osc_prob_earth(H, 1.0*gd.UNIT_GEV, strategy='magnus', **common)
    P_auto = op.osc_prob_earth(H, 1.0*gd.UNIT_GEV, strategy='auto', **common)
    P_hybrid = op.osc_prob_earth(H, 1.0*gd.UNIT_GEV, strategy='hybrid', **common)
    assert maxabs(np.asarray(P_auto) - np.asarray(P_magnus)) == 0.0
    assert maxabs(np.asarray(P_hybrid) - np.asarray(P_magnus)) == 0.0


def test_tolerance_cap_warns():
    """Hitting the refinement caps must warn even at verbose=0 (the result
    can look plausible while being unconverged).

    strategy='magnus' is required here: with the default strategy='auto' (added in 0.11.0), this
    exact (energy, baseline) case is precisely the one the hybrid strategy exists to fix (see
    test_sun_2nu_default_strategy_avoids_tolerance_cap below) and no longer hits the general
    method's refinement caps at all, so forcing strategy='magnus' is what still exercises the
    warning path being tested here."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    with pytest.warns(op.ToleranceNotAchievedWarning):
        op.osc_prob_2nu_sun(10.0*gd.UNIT_MEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM,
                            0.0, sth, Dm2, integration_method='gl',
                            max_n_slabs=64, max_num_loops=8,
                            strategy='magnus', validate_input=False)


def test_sun_2nu_default_strategy_avoids_tolerance_cap():
    """The default strategy='auto' must resolve, warning-free, exactly the case that
    test_tolerance_cap_warns shows still hits the refinement caps under strategy='magnus': 10 MeV
    over 0.9 R_sun is deep enough into the accumulated-phase regime that even the 2-flavor
    interaction-picture fast path (_osc_prob_ip_exp_dispatch) does not certify, so this also
    confirms the hybrid strategy (_osc_prob_hybrid_dispatch) is the one resolving it, not the
    pre-existing fast path."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    energy = 10.0*gd.UNIT_MEV
    L = 0.9*gd.SUN_RADIUS*gd.UNIT_KM

    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        P = op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2, strategy='auto', validate_input=False)
    assert not any(issubclass(w.category, (op.ToleranceNotAchievedWarning,
                                           op.HybridCertificationWarning,
                                           mg.MagnusConvergenceWarning)) for w in wlist)
    assert np.allclose(np.sum(P, axis=1), 1.0, atol=1e-9)

    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    e00 = np.diag([1.0, 0.0])
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)

    def rhs(l, y):
        H = (1.0/energy)*hvac + VCC_func(l)*e00
        return (-1j*H @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-11, atol=1e-13, method='DOP853')
    P_exact = np.abs(sol.y[:, -1].reshape(2, 2)).T**2
    assert maxabs(np.asarray(P) - P_exact) < 1e-2


@pytest.mark.parametrize("energy_mev", [1.0, 5.0, 15.0])
def test_sun_2nu_fast_path_matches_solve_ivp(energy_mev):
    """The interaction-picture fast path for a genuine exponential density profile (Sun-like)
    must reproduce the exact (solve_ivp) probability at realistic, low solar-neutrino energies,
    without hitting the refinement caps or emitting the slab-width convergence warning -- this is
    the regime (large accumulated vacuum phase, far below the 1 GeV point that already saturates
    the general method's default max_n_slabs) the fast path exists to fix. A short baseline (a
    fraction of an e-fold of the density profile) keeps solve_ivp itself tractable at these
    energies while still exercising a genuinely varying matter potential."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    energy = energy_mev*gd.UNIT_MEV
    L = 0.3*gd.L_SCALE_SUN

    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        P = op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2, validate_input=False)
    assert not any(issubclass(w.category, (op.ToleranceNotAchievedWarning,
                                           mg.MagnusConvergenceWarning)) for w in wlist)
    assert np.allclose(np.sum(P, axis=1), 1.0, atol=1e-9)

    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    e00 = np.diag([1.0, 0.0])
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)

    def rhs(l, y):
        H = (1.0/energy)*hvac + VCC_func(l)*e00
        return (-1j*H @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853')
    P_exact = np.abs(sol.y[:, -1].reshape(2, 2)).T**2
    assert maxabs(np.asarray(P) - P_exact) < 1e-3


def test_vectorized_profile_matches_scalar_profile():
    """The silently vectorized Hamiltonian evaluation must give exactly the
    same probabilities as a scalar-only density profile."""
    def rho_vec(l):
        return 3.0 + 2.0*np.sin(np.asarray(l)*1e-12)**2

    def rho_scalar(l):
        if np.ndim(l) != 0:
            raise TypeError("scalar only")
        return float(rho_vec(l))

    common = dict(num_flavors=3, energy=1.0*gd.UNIT_GEV,
                  L=5000.0*gd.UNIT_KM,
                  osc_params={'s12': S12, 's23': S23, 's13': S13,
                              'dCP': DCP, 'D21': D21, 'D31': D31},
                  density_matter_is_in_g_per_cm3=True,
                  validate_input=False)
    P_vec = op.osc_prob_matter_std_potential(rho_func=rho_vec, **common)
    P_scl = op.osc_prob_matter_std_potential(rho_func=rho_scalar, **common)
    assert maxabs(np.asarray(P_vec) - np.asarray(P_scl)) < 1e-12


# ----------------------------------------------------------------------
# Regression tests for the wrapper-layer bugs (A4, B1-B4) found in the
# src/magnus audit and fixed alongside the G1 de-duplication refactor.
# ----------------------------------------------------------------------

RHO_C = 10.0*gd.UNIT_G_PER_CM3
L_SCALE = 100.0*gd.UNIT_KM
S5 = dict(s14=0.1, d14=0.0, s15=0.05, d15=0.0, s24=0.05, d24=0.0, s25=0.02,
         s34=0.02, s35=0.01, d35=0.0, D41=1.0, D51=2.0)
NSI5_ZERO = dict(eps_ee=0.0, eps_em=0.0, eps_et=0.0, eps_es1=0.0, eps_es2=0.0,
                 eps_mm=0.0, eps_mt=0.0, eps_ms1=0.0, eps_ms2=0.0, eps_tt=0.0,
                 eps_ts1=0.0, eps_ts2=0.0, eps_s1s1=0.0, eps_s1s2=0.0, eps_s2s2=0.0)
NSI5_NONZERO = dict(NSI5_ZERO, eps_ee=0.2, eps_em=0.1)


def test_5nu_matter_nsi_exp_density_applies_nsi():
    """Regression test for A4: this wrapper used to call
    osc_prob_matter_std_potential instead of osc_prob_matter_nsi, so its
    epsilon parameters had no effect (or crashed, in the current codebase's
    stricter kwargs handling)."""
    common = dict(energy=ENERGY, L=BASELINE, L0=0.0, rho_central=RHO_C,
                  l_scale=L_SCALE, **S5, validate_input=False)
    P0 = op.osc_prob_5nu_matter_nsi_exp_density(**common, **NSI5_ZERO)
    P1 = op.osc_prob_5nu_matter_nsi_exp_density(**common, **NSI5_NONZERO)
    assert maxabs(np.asarray(P1) - np.asarray(P0)) > 1e-6


def test_5nu_sun_nsi_applies_nsi():
    """osc_prob_5nu_sun_nsi delegates to the function fixed above."""
    common = dict(energy=ENERGY, L=0.5*gd.SUN_RADIUS*gd.UNIT_KM, L0=0.0,
                  **S5, validate_input=False)
    P0 = op.osc_prob_5nu_sun_nsi(**common, **NSI5_ZERO)
    P1 = op.osc_prob_5nu_sun_nsi(**common, **NSI5_NONZERO)
    assert maxabs(np.asarray(P1) - np.asarray(P0)) > 1e-6


@pytest.mark.parametrize("name", ['osc_prob_2nu_sun_liv', 'osc_prob_3nu_sun_liv',
                                  'osc_prob_4nu_sun_liv', 'osc_prob_5nu_sun_liv'])
def test_sun_liv_atol_independent_of_rtol(name):
    """Regression test for B1: these four functions used to pass
    atol=rtol to the inner call, silently discarding the requested atol.
    Not directly observable from outputs (both are now just **kwargs
    passthrough, with no per-family atol= line left to typo at all), so
    this test inspects the source directly -- the same check that would
    have caught the original bug."""
    import inspect
    fn = getattr(op, name)
    src = inspect.getsource(fn)
    assert 'atol=rtol' not in src


@pytest.mark.parametrize("name", ['osc_prob_3nu_matter_nsi_constant_density',
                                  'osc_prob_4nu_matter_nsi_constant_density',
                                  'osc_prob_5nu_matter_nsi_constant_density'])
def test_matter_nsi_constant_density_nubar_has_effect(name):
    """Regression test for B2: these three were missing nubar entirely."""
    fn = getattr(op, name)
    import inspect
    params = inspect.signature(fn).parameters
    assert 'nubar' in params
    kwargs = dict(energy=ENERGY, L=BASELINE, rho=RHO_C, eps_ee=0.1, eps_em=0.05j,
                 validate_input=False)
    if 's14' in params:
        kwargs.update(S5 if 's15' in params else
                      {k: v for k, v in S5.items() if k not in ('s15', 'd15', 's25', 's35', 'd35', 'D51')})
    P_nu = fn(**kwargs, nubar=False)
    P_nubar = fn(**kwargs, nubar=True)
    assert maxabs(np.asarray(P_nu) - np.asarray(P_nubar)) > 1e-6


def test_matter_liv_constant_density_2nu_nubar_has_effect():
    """Regression test for B3: osc_prob_2nu_matter_liv_constant_density was
    uniquely missing nubar among its family (the matter potential is
    nonzero here, so nubar has a real, physical effect unlike the pure-LIV
    or pure-vacuum 2nu cases)."""
    import inspect
    assert 'nubar' in inspect.signature(op.osc_prob_2nu_matter_liv_constant_density).parameters
    common = dict(energy=ENERGY, L=BASELINE, rho=RHO_C, sth=0.3, Dm2=2.5e-3,
                 sxi=0.2, b1=gd.B1, b2=gd.B2, Lambda=gd.LAMBDA, n_liv=1,
                 validate_input=False)
    P_nu = op.osc_prob_2nu_matter_liv_constant_density(**common, nubar=False)
    P_nubar = op.osc_prob_2nu_matter_liv_constant_density(**common, nubar=True)
    assert maxabs(np.asarray(P_nu) - np.asarray(P_nubar)) > 1e-6


@pytest.mark.parametrize("name", ['osc_prob_2nu_matter_nsi_exp_density',
                                  'osc_prob_3nu_matter_nsi_exp_density',
                                  'osc_prob_4nu_matter_nsi_exp_density',
                                  'osc_prob_5nu_matter_nsi_exp_density'])
def test_nsi_exp_density_accepts_zero_rho_central(name):
    """Regression test for B4: the NSI exp-density family rejected
    rho_central == 0.0 (`<= 0.0` check) while the plain and LIV exp-density
    families correctly allow it (`< 0.0`)."""
    fn = getattr(op, name)
    import inspect
    params = inspect.signature(fn).parameters
    kwargs = dict(energy=ENERGY, L=BASELINE, L0=0.0, rho_central=0.0,
                 l_scale=L_SCALE, validate_input=True)
    if 'sth' in params:
        kwargs.update(sth=0.3, Dm2=2.5e-3, eps_aa=0.1, eps_ab=0.0)
    else:
        # 3/4/5nu: NSI eps parameters default to 0.0 already
        pass
    # Should not raise / sys.exit
    P = fn(**kwargs)
    assert np.all(np.isfinite(np.asarray(P)))


# ----------------------------------------------------------------------
# Structural smoke test (G1): every physics wrapper family, called once
# per flavor count, must run without error and return a valid, unitary
# probability matrix. Exercises code paths the targeted tests above don't.
# ----------------------------------------------------------------------

@pytest.mark.parametrize("num_flavors", [2, 3, 4, 5])
@pytest.mark.parametrize("family", ['vacuum', 'matter_constant_density', 'matter_exp_density',
                                    'earth', 'sun'])
def test_every_standard_wrapper_runs_and_is_unitary(family, num_flavors):
    fn = getattr(op, f'osc_prob_{num_flavors}nu_{family}')
    kwargs = dict(validate_input=False)
    if num_flavors == 2:
        kwargs.update(sth=0.3, Dm2=2.5e-3)
    elif num_flavors == 4:
        kwargs.update({k: v for k, v in S5.items()
                      if k not in ('s15', 'd15', 's25', 's35', 'd35', 'D51')})
    elif num_flavors == 5:
        kwargs.update(S5)

    if family == 'vacuum':
        P = fn(ENERGY, BASELINE, **kwargs)
    elif family == 'matter_constant_density':
        P = fn(ENERGY, BASELINE, RHO_C, **kwargs)
    elif family == 'matter_exp_density':
        P = fn(ENERGY, BASELINE, 0.0, RHO_C, L_SCALE, **kwargs)
    elif family == 'earth':
        P = fn(ENERGY, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM, **kwargs)
    elif family == 'sun':
        P = fn(ENERGY, 0.5*gd.SUN_RADIUS*gd.UNIT_KM, 0.0, **kwargs)

    P = np.asarray(P)
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-6)
    assert np.all((P >= -1e-9) & (P <= 1.0 + 1e-9))


# ----------------------------------------------------------------------
# Permanent API-consistency guard (see the G1 write-up): every bug in
# category B was a sibling wrapper silently drifting from its family's
# convention (a stray default, a missing parameter, an inconsistent
# validation bound). This test inspects the *shape* of the whole wrapper
# family at once, the same way the audit that found those bugs did, so
# a future copy-pasted wrapper cannot silently reintroduce the pattern.
# ----------------------------------------------------------------------

def test_no_wrapper_redeclares_standard_refinement_kwargs():
    """None of the osc_prob_{2,3,4,5}nu_* physics wrappers should declare
    the standard refinement/logging kwargs in their own signature -- they
    must flow through **kwargs from the single middle-layer default.
    Guards against the 'fat wrapper' pattern (G1) recurring, which is
    where every B1/E6-E8 bug hid. (Restricted to the {N}nu_* wrapper
    functions themselves -- the middle-layer functions they delegate to,
    e.g. osc_prob_vacuum/osc_prob_matter_nsi/osc_prob_liv, and the
    primordial osc_prob/osc_prob_energy_baseline, are the canonical
    source of these defaults and are expected to declare them.)"""
    import inspect
    import re
    standard_params = {
        'magnus_exp_order', 'n_jobs', 'integration_method', 'rtol', 'atol',
        'growth_factor_n_slabs', 'growth_factor_n_tpts_per_slab', 'max_num_loops',
        'min_n_slabs', 'max_n_slabs', 'min_n_tpts_per_slab', 'max_n_tpts_per_slab',
        'iterate_over_magnus_exp_order', 'min_magnus_exp_order', 'max_magnus_exp_order',
        'new_recursion_limit',
    }
    wrapper_pattern = re.compile(r'^osc_prob_[2345]nu_')
    offenders = {}
    for name in dir(op):
        if not wrapper_pattern.match(name):
            continue
        fn = getattr(op, name)
        if not inspect.isfunction(fn):
            continue
        params = set(inspect.signature(fn).parameters)
        hit = params & standard_params
        if hit:
            offenders[name] = hit
    assert not offenders, offenders
    # Sanity check that the pattern actually matched a substantial number of
    # real wrapper functions (i.e., the test isn't vacuously passing)
    n_matched = sum(1 for name in dir(op) if wrapper_pattern.match(name))
    assert n_matched >= 50, n_matched


@pytest.mark.parametrize("family_prefix", ['osc_prob_{n}nu_matter_nsi_constant_density',
                                           'osc_prob_{n}nu_matter_nsi_exp_density',
                                           'osc_prob_{n}nu_matter_liv_constant_density',
                                           'osc_prob_{n}nu_matter_liv_exp_density',
                                           'osc_prob_{n}nu_earth_nsi', 'osc_prob_{n}nu_earth_liv',
                                           'osc_prob_{n}nu_sun_nsi', 'osc_prob_{n}nu_sun_liv'])
def test_nubar_present_across_all_flavor_counts_in_matter_families(family_prefix):
    """Every matter/NSI/LIV family (where the potential is generically
    nonzero, so nubar always has a real physical effect) must expose
    nubar for all four flavor counts, not just some (the exact shape of
    bugs B2/B3)."""
    import inspect
    missing = []
    for n in (2, 3, 4, 5):
        name = family_prefix.format(n=n)
        fn = getattr(op, name)
        if 'nubar' not in inspect.signature(fn).parameters:
            missing.append(name)
    assert not missing, missing
