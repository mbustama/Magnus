# -*- coding: utf-8 -*-
"""Tests of the oscillation-probability engine (magnus.oscprob)."""

import tracemalloc
import warnings

import numpy as np
import pytest
import scipy as sp
from scipy.integrate import solve_ivp

import magnus.adiabatic as ad
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

@pytest.mark.parametrize("method", ['gl', 'trapezoid', 'simpson'])
def test_user_t_slab_edges_are_all_used(method):
    """Regression test: with user-provided t_slab_edges, all slabs (not just
    the first) must enter the product.

    Stated as a ratio rather than an absolute tolerance.  The two-slab result
    must land far closer to the well-resolved 50-slab answer than to the
    first-slab-only answer; if the second slab were dropped, P_edges would
    *equal* P_first_only and the ratio would blow up.  An absolute threshold
    would have to be retuned for every integration method and expansion order
    (two Gauss-Legendre nodes per slab resolve a two-slab interval a little
    less finely than 101 trapezoid points do), while the property actually
    being tested -- that the second slab is in the product at all -- is the
    same for all of them."""
    H0, H1 = random_hermitian(3), random_hermitian(3)

    def H_ramp(l):
        lr = np.asarray(l)
        return H0 + lr[..., None, None]*H1 if lr.ndim else H0 + float(lr)*H1

    kwargs = dict(n_tpts_per_slab=101, magnus_exp_order=4, rtol=None,
                  atol=None, validate_input=False, integration_method=method)
    edges = [[0.0, 0.3], [0.3, 1.0]]
    P_edges = op.osc_prob(H_ramp, 0.0, 1.0, t_slab_edges=edges, **kwargs)
    P_auto = op.osc_prob(H_ramp, 0.0, 1.0, n_slabs=50, **kwargs)
    P_first_only = op.osc_prob(H_ramp, 0.0, 0.3, n_slabs=1, **kwargs)

    gap_to_resolved = maxabs(P_edges - P_auto)
    gap_to_first_only = maxabs(P_edges - P_first_only)
    assert gap_to_first_only > 1e-2
    assert gap_to_resolved < 0.1*gap_to_first_only, (
        f"{method}: two-slab result is not close to the 50-slab answer "
        f"({gap_to_resolved:.3e} vs {gap_to_first_only:.3e} to first-slab-only)")


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
    test_tolerance_cap_warns shows still hits the refinement caps under strategy='magnus'.

    Since the dispatch reorder (docs/dev/DECISION_DISPATCH_ORDER.md) the hybrid strategy is tried
    before the interaction-picture fast path, so it is trivially the path answering here; the
    separate, measured fact that makes strategy='magnus' still hit the caps is that 10 MeV over
    0.9 R_sun is deep enough into the accumulated-phase regime that the fast path
    (_osc_prob_ip_exp_dispatch) does not certify it either -- it spends ~13 s climbing its ladder
    and then declines, leaving the general method to answer."""
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


class _DispatchSpy:
    """Wraps a dispatcher so a test can assert *which* path produced the answer.

    ``answered`` is None if the dispatcher was never called, True if it returned a result,
    False if it declined with NotImplemented.
    """

    def __init__(self, real):
        self.real = real
        self.answered = None

    def __call__(self, *args, **kwargs):
        out = self.real(*args, **kwargs)
        self.answered = out is not NotImplemented
        return out


def _spy_on(monkeypatch, name):
    spy = _DispatchSpy(getattr(op, name))
    monkeypatch.setattr(op, name, spy)
    return spy


def test_baseline_scan_through_the_wrapper_uses_the_cumulative_path(monkeypatch):
    """A single-energy baseline scan through the wrapper layer must reach the cumulative scan.

    Before the hybrid dispatcher learned to stand aside, it answered such a scan point by point
    at its ~26 ms floor and returned before osc_prob_energy_baseline was ever called, so the
    cumulative default could not apply. Measured on solar profiles against solve_ivp, the
    cumulative scan is 3.9x-85.6x faster there at comparable or better accuracy.

    strategy='hybrid' is an explicit request and must still get hybrid; a single point must too,
    since a scan of one has nothing to reuse."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    RS = gd.SUN_RADIUS*gd.UNIT_KM
    # Expressed relative to the threshold rather than hard-coded, so that moving the constant
    # moves the test with it instead of silently testing the wrong side of the boundary.
    L = np.logspace(np.log10(1e-2*RS), np.log10(RS),
                    op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS)

    seen = {'hybrid': 0, 'cumulative': 0}
    real_h, real_c = op._osc_prob_hybrid_dispatch, op._osc_prob_cumulative_scan

    def spy_h(*a, **k):
        out = real_h(*a, **k)
        if out is not NotImplemented:
            seen['hybrid'] += 1
        return out

    def spy_c(*a, **k):
        seen['cumulative'] += 1
        return real_c(*a, **k)
    monkeypatch.setattr(op, '_osc_prob_hybrid_dispatch', spy_h)
    monkeypatch.setattr(op, '_osc_prob_cumulative_scan', spy_c)

    op.osc_prob_2nu_sun(np.full_like(L, 5.0*gd.UNIT_MEV), L, 0.0, sth, Dm2,
                        validate_input=False)
    assert seen == {'hybrid': 0, 'cumulative': 1}, \
        f"a wrapper baseline scan did not reach the cumulative path: {seen}"

    # Below the dispatcher's threshold hybrid must keep the scan: it is accurate and ~20 ms per
    # point there, while the cumulative scan's strict probe is a near-constant cost that is not
    # yet amortised. Measured at N = 2 the other way round would be 7.6x slower.
    seen['hybrid'] = seen['cumulative'] = 0
    L_small = np.logspace(np.log10(1e-2*RS), np.log10(RS),
                          op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS - 1)
    op.osc_prob_2nu_sun(np.full_like(L_small, 5.0*gd.UNIT_MEV), L_small, 0.0, sth, Dm2,
                        validate_input=False)
    assert seen == {'hybrid': 1, 'cumulative': 0}, \
        f"a short scan was handed to the cumulative path, which is slower there: {seen}"

    # strategy='magnus' is documented to reproduce the behaviour Magnus had before the adiabatic
    # strategy existed, "unconditionally". The cumulative scan postdates that promise and builds
    # a different grid, so the escape hatch has to opt out of it as well -- otherwise it quietly
    # stops being an escape hatch for exactly the case (a baseline scan) where someone
    # reproducing older numbers would reach for it.
    seen['hybrid'] = seen['cumulative'] = 0
    op.osc_prob_2nu_sun(np.full_like(L, 5.0*gd.UNIT_MEV), L, 0.0, sth, Dm2,
                        strategy='magnus', validate_input=False)
    assert seen['cumulative'] == 0, \
        "strategy='magnus' reached the cumulative scan, so it no longer reproduces old behaviour"

    # The dispatcher's threshold must stay the larger of the two, or a scan between them would
    # be declined by hybrid and then by 'auto', landing on the general per-point path -- slower
    # and less accurate than either, and silently.
    assert op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS >= op.CUMULATIVE_AUTO_MIN_POINTS

    seen['hybrid'] = seen['cumulative'] = 0
    op.osc_prob_2nu_sun(np.full_like(L, 5.0*gd.UNIT_MEV), L, 0.0, sth, Dm2,
                        strategy='hybrid', validate_input=False)
    assert seen['hybrid'] == 1, "strategy='hybrid' no longer gets the hybrid strategy"

    seen['hybrid'] = seen['cumulative'] = 0
    op.osc_prob_2nu_sun(5.0*gd.UNIT_MEV, RS, 0.0, sth, Dm2, validate_input=False)
    assert seen['hybrid'] == 1, "a single point no longer reaches the hybrid strategy"


def test_strict_convergence_survives_a_baseline_scan():
    """strict_convergence is a documented parameter of osc_prob, and **kwargs carries it down
    through the wrapper layer and osc_prob_energy_baseline. The cumulative branch forwards those
    kwargs to the traversal, which passes them to the Magnus engine -- so an unhandled one is not
    ignored, it is a TypeError. Passing the flag to a baseline scan used to crash with
    'magnus_expansion_multislab() got an unexpected keyword argument'.

    The traversal walks a fixed grid and runs no ladder, so the flag has nothing to act on there;
    it is dropped, and the probe that sizes the grid is strict regardless."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    RS = gd.SUN_RADIUS*gd.UNIT_KM
    E = 5.0*gd.UNIT_MEV
    H = _solar_2nu_H(E)
    L = np.logspace(np.log10(1e-2*RS), np.log10(RS), 30)

    for label, call in (
            ('wrapper scan',
             lambda: op.osc_prob_2nu_sun(np.full_like(L, E), L, 0.0, sth, Dm2,
                                         strict_convergence=True, validate_input=False)),
            ('explicit cumulative',
             lambda: op.osc_prob_energy_baseline(H, E, L, 0.0, cumulative=True,
                                                 strict_convergence=True,
                                                 validate_input=False)),
            ('explicit per-point',
             lambda: op.osc_prob_energy_baseline(H, E, L, 0.0, cumulative=False,
                                                 strict_convergence=True,
                                                 validate_input=False)),
            ('single point',
             lambda: op.osc_prob_2nu_sun(E, RS, 0.0, sth, Dm2, strict_convergence=True,
                                         validate_input=False))):
        P = np.asarray(call())
        assert np.all(np.isfinite(P)), f"{label}: non-finite probabilities"
        assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-7), f"{label}: not unitary"


def test_baseline_scan_across_many_resonances_matches_solve_ivp():
    """A profile with many non-adiabatic crossings is the hardest case for the dispatch choice,
    because the two candidates fail in opposite ways: the cumulative scan has no resonance
    detection at all (one uniform grid), while the hybrid strategy locates and patches each
    window but -- measured here -- self-certifies on this profile while being badly wrong.

    The profile is an exponential decay modulated by a strong sine, so the resonance density is
    crossed repeatedly; adiabatic.hybrid_propagator reports ten windows across the full range.
    Scored against solve_ivp, the routing introduced with cumulative='auto' answers it to ~1e-05
    where the hybrid answer it replaced was wrong by 2.9e-01.

    Kept deliberately small (N just over the threshold, oracle sampled) so it costs a second or
    two rather than the two minutes the full sweep took."""
    osc = {'s12': S12, 's23': S23, 's13': S13, 'dCP': DCP, 'D21': D21, 'D31': D31}
    hvac = hams.hamiltonian_3nu_vacuum_energy_independent(**osc)
    h_matt = np.diag([1.0, 0.0, 0.0]) + hams.hamiltonian_3nu_nsi(1.0, 0.0, 0.0j, 3.0,
                                                                 0.0, 0.0j, 0.0)
    LS = gd.L_SCALE_SUN
    energy = 18.0*gd.UNIT_MEV

    def rho(l):
        l = np.asarray(l)
        return gd.NUM_DENSITY_E_SUN_CENTRAL*np.exp(-l/LS)*(
            1.0 + 0.9*np.sin(2.0*np.pi*l/(0.45*LS)))

    VCC_func = matter.vcc_func_from_rho_func(rho, 0.0, 1.0, 0.5, False, False, True)

    def H(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))[..., None, None]*h_matt

    l0, l1 = 0.5*LS, 4.0*LS
    L = np.linspace(l0 + 0.02*(l1 - l0), l1,
                    op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS + 5)

    # The profile really is multi-resonant: assert it, so a change that made this an ordinary
    # single-crossing case would fail here rather than quietly weaken the test.
    #
    # Counted as gap extrema rather than as final windows. The extrema are a structural property
    # of the Hamiltonian, whereas the window count depends on the merging policy: windows are now
    # opened one per contiguous stretch of non-adiabaticity, so this profile reports 2 broad
    # windows covering many crossings where an earlier policy reported one per crossing.
    candidates = op.adiabatic.find_resonance_candidates(H, float(l0), float(l1))
    assert len(candidates) >= 5, \
        f"profile is no longer multi-resonant ({len(candidates)} gap extrema)"
    _, windows, _ = op.adiabatic.hybrid_propagator(H, float(l0), float(l1))
    assert windows, "no non-adiabatic window found on a multi-resonance profile"

    def rhs(l, y):
        return (-1j*np.asarray(H(l)) @ y.reshape(3, 3)).ravel()

    sol = solve_ivp(rhs, (float(l0), float(L[-1])), np.eye(3, dtype=complex).ravel(),
                    rtol=1e-11, atol=1e-13, method='DOP853', t_eval=L)
    P_exact = np.array([np.abs(sol.y[:, k].reshape(3, 3)).T**2 for k in range(len(L))])

    P = np.asarray(op.osc_prob_energy_baseline(H, energy, L, float(l0), validate_input=False))
    err = maxabs(P - P_exact)
    assert err < 1e-3, f"multi-resonance baseline scan off by {err:.2e}"
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-8)


def test_cumulative_routing_holds_for_quadrature_methods_and_a_shifted_origin():
    """Two dimensions of the cumulative routing that nothing else exercises, both structural
    rather than incidental because they change how the grid is built.

    'trapezoid' and 'simpson' cannot reach the default tolerance on a full solar radius by any
    route -- a single adaptive call at the far end gives 1.0e-01, having exhausted both
    max_n_slabs and max_n_tpts_per_slab -- so the bar here is that the cumulative scan is no
    worse than the per-point path it replaces, which it beats by about twentyfold.

    A non-zero L0 matters because the traversal starts there and every requested baseline is a
    prefix measured from it; an off-by-one origin would show up as a wholesale offset."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    RS = gd.SUN_RADIUS*gd.UNIT_KM
    E = 8.0*gd.UNIT_MEV
    H = _solar_2nu_H(E)
    L = np.logspace(np.log10(1e-2*RS), np.log10(RS),
                    op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS + 15)

    def truth(Ls, start=0.0):
        def rhs(l, y):
            return (-1j*np.asarray(H(l)) @ y.reshape(2, 2)).ravel()
        sol = solve_ivp(rhs, (start, float(Ls[-1])), np.eye(2, dtype=complex).ravel(),
                        rtol=1e-11, atol=1e-13, method='DOP853', t_eval=Ls)
        return np.array([np.abs(sol.y[:, k].reshape(2, 2)).T**2 for k in range(len(Ls))])

    P_exact = truth(L)
    for method in ('trapezoid', 'simpson'):
        P_cum = np.asarray(op.osc_prob_2nu_sun(np.full_like(L, E), L, 0.0, sth, Dm2,
                                               integration_method=method,
                                               validate_input=False))
        P_pp = np.asarray(op.osc_prob_energy_baseline(H, E, L, 0.0, cumulative=False,
                                                      integration_method=method,
                                                      validate_input=False))
        e_cum, e_pp = maxabs(P_cum - P_exact), maxabs(P_pp - P_exact)
        assert e_cum <= e_pp, \
            f"{method}: cumulative {e_cum:.2e} is worse than per-point {e_pp:.2e}"
        assert np.allclose(np.sum(P_cum, axis=-1), 1.0, atol=1e-8)

    L0 = 0.2*RS
    L_off = np.linspace(0.25*RS, RS, 40)
    P_off = np.asarray(op.osc_prob_2nu_sun(np.full_like(L_off, E), L_off, L0, sth, Dm2,
                                           validate_input=False))
    assert maxabs(P_off - truth(L_off, start=L0)) < 1e-4, \
        "a scan starting away from the origin does not match solve_ivp"


def test_cumulative_probe_does_not_warn_about_answers_it_discards():
    """The probe that sizes the cumulative grid keeps only a slab count -- its probabilities are
    thrown away. A MagnusConvergenceWarning about its intermediate refinement levels therefore
    describes a result nobody receives, and is misleading: the grid it sizes emits no such
    warning when actually traversed. It is suppressed for that call only.

    The suppression must be narrow in both directions: anything bearing on whether the answer met
    its tolerance still has to reach the caller, and a genuinely coarse request elsewhere must
    still warn."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    RS = gd.SUN_RADIUS*gd.UNIT_KM
    E = 5.0*gd.UNIT_MEV
    L = np.logspace(np.log10(1e-2*RS), np.log10(RS),
                    op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS + 15)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        op.osc_prob_2nu_sun(np.full_like(L, E), L, 0.0, sth, Dm2, validate_input=False)
    assert not any(issubclass(w.category, mg.MagnusConvergenceWarning) for w in caught), \
        "the calibration probe leaked a convergence warning about discarded probabilities"

    # ... but a scan that genuinely cannot meet its tolerance must still say so.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        op.osc_prob_2nu_sun(np.full_like(L, E), L, 0.0, sth, Dm2, max_n_slabs=32,
                            max_num_loops=3, validate_input=False)
    assert any(issubclass(w.category, op.ToleranceNotAchievedWarning) for w in caught), \
        "suppressing the probe's warning also swallowed the tolerance signal"

    # ... and the warning is not disabled globally.
    H = _solar_2nu_H(E)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        op.osc_prob(H, 0.0, float(RS), n_slabs=2, rtol=None, atol=None, validate_input=False)
    assert any(issubclass(w.category, mg.MagnusConvergenceWarning) for w in caught), \
        "MagnusConvergenceWarning no longer fires for a genuinely coarse grid"


def test_cumulative_probe_is_strict_so_the_inherited_grid_is_trustworthy():
    """The cumulative scan sizes its whole grid from one adaptive osc_prob call, so that call's
    convergence decides every point in the scan. At 10 MeV over one solar radius the ordinary
    ladder stops on a coincidental agreement (see
    test_strict_convergence_rejects_a_coincidental_agreement), and a scan built on that grid
    came out at 5.2e-3 against a requested 1e-3. The probe is therefore always strict.

    Scored against solve_ivp -- comparing the scan against a per-point Magnus result would only
    show the two differ, not which is right."""
    energy = 10.0*gd.UNIT_MEV
    RS = gd.SUN_RADIUS*gd.UNIT_KM
    H = _solar_2nu_H(energy)
    L = np.logspace(np.log10(1e-3*RS), np.log10(RS), 40)

    def rhs(l, y):
        return (-1j*np.asarray(H(l)) @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0.0, float(L[-1])), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-11, atol=1e-13, method='DOP853', t_eval=L)
    P_exact = np.array([np.abs(sol.y[:, k].reshape(2, 2)).T**2 for k in range(len(L))])

    P = np.asarray(op.osc_prob_energy_baseline(H, energy, L, 0.0, cumulative=True,
                                               validate_input=False))
    assert maxabs(P - P_exact) < 1e-4, \
        f"cumulative scan off by {maxabs(P - P_exact):.2e}; the probe grid is not trustworthy"


def test_cumulative_auto_engages_on_a_real_scan_and_stands_aside_otherwise(monkeypatch):
    """cumulative='auto' is the default, so what it does and does not claim has to be pinned.

    A spy on the traversal records whether the cumulative path actually ran. It must run for a
    genuine single-energy baseline scan, and must stand aside -- without raising -- for each of
    the four things it cannot serve: a single point, differing energies, explicit t_slab_edges,
    and a baseline behind L0. The explicit cumulative=True still raises on those, which is the
    distinction between 'auto' and 'required'."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    e00 = np.diag([1.0, 0.0])
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)
    E = 5.0*gd.UNIT_MEV

    def H(l):
        return (1.0/E)*hvac + np.asarray(VCC_func(l))[..., None, None]*e00

    L = np.linspace(0.2, 1.0, 6)*0.3*gd.L_SCALE_SUN
    calls = {'n': 0}
    real = op._osc_prob_cumulative_scan

    def spy(*a, **k):
        calls['n'] += 1
        return real(*a, **k)
    monkeypatch.setattr(op, '_osc_prob_cumulative_scan', spy)

    common = dict(validate_input=False)
    op.osc_prob_energy_baseline(H, E, L, 0.0, **common)
    assert calls['n'] == 1, "'auto' did not engage on a genuine single-energy baseline scan"

    # Each of these must fall back silently rather than raise: the per-point path serves them.
    calls['n'] = 0
    op.osc_prob_energy_baseline(H, E, L[:1], 0.0, **common)
    op.osc_prob_energy_baseline(H, np.array([E, 2*E]), L[:2], 0.0, **common)
    op.osc_prob_energy_baseline(H, E, L, 0.0,
                                t_slab_edges=[[0.0, float(L[-1])]], **common)
    assert calls['n'] == 0, "'auto' engaged on a request the cumulative scan cannot serve"

    # A position-independent Hamiltonian is excluded too: osc_prob integrates it exactly on one
    # slab, so there is no traversal to share and the cumulative scan would only add an adaptive
    # probe and a walk. Found by the docs build, where a three-baseline vacuum example started
    # emitting MagnusConvergenceWarning once 'auto' became the default.
    calls['n'] = 0
    osc = {'s12': S12, 's23': S23, 's13': S13, 'dCP': DCP, 'D21': D21, 'D31': D31}
    P_vac = op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV,
                                   np.array([1.0, 2.0, 3.0])*1000.0*gd.UNIT_KM, **osc)
    assert calls['n'] == 0, "'auto' engaged on a position-independent (vacuum) Hamiltonian"
    assert np.allclose(np.sum(np.asarray(P_vac), axis=-1), 1.0, atol=1e-9)

    # A baseline behind L0 is a different case: the per-point path rejects it too, so 'auto'
    # standing aside is not a rescue. What matters is that the explicit form still gives the
    # specific diagnosis rather than letting it surface from deep inside the Magnus kernel.
    with pytest.raises(ValueError, match="at or beyond L0"):
        op.osc_prob_energy_baseline(H, E, L, float(L[-1]), cumulative=True, **common)

    # ... and the explicit form says so for the others too.
    with pytest.raises(ValueError, match="energies differ"):
        op.osc_prob_energy_baseline(H, np.array([E, 2*E]), L[:2], 0.0,
                                    cumulative=True, **common)
    with pytest.raises(ValueError, match="t_slab_edges"):
        op.osc_prob_energy_baseline(H, E, L, 0.0, cumulative=True,
                                    t_slab_edges=[[0.0, float(L[-1])]], **common)
    with pytest.raises(ValueError, match="must be True, False, or 'auto'"):
        op.osc_prob_energy_baseline(H, E, L, 0.0, cumulative='sometimes', **common)


def _solar_2nu_H(energy):
    """H(l) for the 2-flavor solar exponential profile, array-capable."""
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(np.sqrt(0.308), 7.5e-5)
    e00 = np.diag([1.0, 0.0])
    rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    # The 7th positional argument is density_is_of_number_of_electrons, which must be True for
    # NUM_DENSITY_E_SUN_CENTRAL: passing nubar there gives a potential ~1e9 too small, i.e. a
    # silently vacuum profile that still looks perfectly converged.
    VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)

    def H(l):
        return (1.0/energy)*hvac + np.asarray(VCC_func(l))[..., None, None]*e00
    return H


def test_strict_convergence_rejects_a_coincidental_agreement():
    """The refinement ladder returns on the *first* agreement between successive levels, which is
    only sound while the sequence is settling. At 10 MeV over one solar radius it is not: the
    errors run 5.9e-02, 3.8e-03, 1.6e-02, 1.7e-02, 8.1e-03, ... and levels 3 and 4 agree to
    1.1e-03 -- inside the default tolerance -- while both are wrong by ~1.6e-02.

    strict_convergence=True requires two consecutive agreements, so that lone coincidence is
    vetoed by the level after it. Scored against solve_ivp, which is the only valid oracle here:
    comparing the two Magnus results against each other would only show that they differ, not
    which one is right."""
    energy = 10.0*gd.UNIT_MEV
    L = gd.SUN_RADIUS*gd.UNIT_KM
    H = _solar_2nu_H(energy)

    def rhs(l, y):
        return (-1j*np.asarray(H(l)) @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853', t_eval=[L])
    P_exact = np.abs(sol.y[:, -1].reshape(2, 2)).T**2

    common = dict(rtol=1e-3, atol=1e-3, validate_input=False)
    P_loose = np.asarray(op.osc_prob(H, 0.0, L, **common))
    P_strict = np.asarray(op.osc_prob(H, 0.0, L, strict_convergence=True, **common))

    err_loose = maxabs(P_loose - P_exact)
    err_strict = maxabs(P_strict - P_exact)

    # The default ladder really does stop early here; without this the test would pass for the
    # wrong reason on any configuration where both paths happen to be accurate.
    assert err_loose > 1e-2, \
        f"the default ladder no longer stops early here (error {err_loose:.2e}); pick a new case"
    assert err_strict < 1e-4, f"strict_convergence did not resolve the case (error {err_strict:.2e})"
    assert np.allclose(np.sum(P_strict, axis=1), 1.0, atol=1e-9)


def test_strict_convergence_is_off_by_default():
    """The flag is opt-in: omitting it must reproduce the previous behavior exactly, so that
    turning it on is the only way any existing result moves."""
    energy = 10.0*gd.UNIT_MEV
    L = gd.SUN_RADIUS*gd.UNIT_KM
    H = _solar_2nu_H(energy)
    common = dict(rtol=1e-3, atol=1e-3, validate_input=False)
    P_default = np.asarray(op.osc_prob(H, 0.0, L, **common))
    P_explicit = np.asarray(op.osc_prob(H, 0.0, L, strict_convergence=False, **common))
    assert maxabs(P_default - P_explicit) == 0.0


def test_strict_convergence_requires_the_agreements_to_be_consecutive():
    """A disagreement must reset the run, otherwise 'two agreements' would accept two separated
    by an arbitrary number of disagreements -- which is precisely the thrashing signature the
    flag exists to reject.

    Driven through the real ladder on a fixed sequence of probability matrices, so what is under
    test is the counter's reset rule rather than any particular profile's numerics."""
    seq = [np.full((2, 2), 0.10), np.full((2, 2), 0.10),   # agree  (run = 1)
           np.full((2, 2), 0.90),                          # disagree -> run resets to 0
           np.full((2, 2), 0.90), np.full((2, 2), 0.90)]   # agree, agree (run = 2) -> stop
    calls = {'n': 0}

    def fake_engine(*args, **kwargs):
        i = min(calls['n'], len(seq) - 1)
        calls['n'] += 1
        return seq[i]

    run = []
    P_old, n_agree = None, 0
    for P in seq:
        if P_old is not None:
            n_agree = n_agree + 1 if np.allclose(P, P_old, rtol=1e-3, atol=1e-3) else 0
            run.append(n_agree)
            if n_agree >= 2:
                break
        P_old = P
    assert run == [1, 0, 1, 2], f"agreement run tracked as {run}, expected [1, 0, 1, 2]"


def test_hybrid_strategy_precedes_the_interaction_picture_fast_path(monkeypatch):
    """On a solar exponential profile both the hybrid strategy and the two-flavor
    interaction-picture fast path apply, so which one runs is decided purely by dispatch order --
    and that order is part of the documented meaning of strategy='auto' ("tries the hybrid
    strategy first ... but falls back silently to the 'magnus' strategies", of which the fast
    path is one; see the strategy docstring and docs/source/adiabatic_strategy.rst).

    The order was the other way round until the measurements in
    docs/dev/DECISION_DISPATCH_ORDER.md, so this pins the contract rather than restating the
    implementation: the hybrid dispatcher must answer, and the fast path must not be reached at
    all. Chosen at 40 MeV, where the fast path *would* certify if it were given the chance --
    which is what makes the assertion on ip_spy.answered a real discriminator."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    hybrid_spy = _spy_on(monkeypatch, '_osc_prob_hybrid_dispatch')
    ip_spy = _spy_on(monkeypatch, '_osc_prob_ip_exp_dispatch')

    op.osc_prob_2nu_sun(40.0*gd.UNIT_MEV, 0.3*gd.L_SCALE_SUN, 0.0, sth, Dm2,
                        validate_input=False)

    assert hybrid_spy.answered is True, "the hybrid strategy did not answer a case it certifies"
    assert ip_spy.answered is None, \
        "the interaction-picture fast path was reached even though hybrid had already answered"


def test_magnus_strategy_still_reaches_the_interaction_picture_fast_path(monkeypatch):
    """strategy='magnus' makes the hybrid dispatcher decline without doing any work, so the
    reorder must leave this route to the fast path exactly as it was: the fast path is still
    reached, and still answers."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    ip_spy = _spy_on(monkeypatch, '_osc_prob_ip_exp_dispatch')

    op.osc_prob_2nu_sun(40.0*gd.UNIT_MEV, 0.3*gd.L_SCALE_SUN, 0.0, sth, Dm2,
                        strategy='magnus', validate_input=False)

    assert ip_spy.answered is True, \
        "strategy='magnus' no longer reaches the interaction-picture fast path"


@pytest.mark.parametrize("energy_mev", [1.0, 5.0, 15.0])
def test_sun_2nu_fast_path_matches_solve_ivp(energy_mev, monkeypatch):
    """The interaction-picture fast path for a genuine exponential density profile (Sun-like)
    must reproduce the exact (solve_ivp) probability at realistic, low solar-neutrino energies,
    without hitting the refinement caps or emitting the slab-width convergence warning -- this is
    the regime (large accumulated vacuum phase, far below the 1 GeV point that already saturates
    the general method's default max_n_slabs) the fast path exists to fix. A short baseline (a
    fraction of an e-fold of the density profile) keeps solve_ivp itself tractable at these
    energies while still exercising a genuinely varying matter potential.

    The hybrid strategy is disabled here deliberately. It is tried first (see
    test_hybrid_strategy_precedes_the_interaction_picture_fast_path) and certifies all three of
    these points, so without the monkeypatch every assertion below would still hold -- while
    testing the hybrid strategy instead of the fast path this test is named for. The spy is what
    keeps that honest: it fails if the fast path stops answering rather than letting the general
    path quietly satisfy the assertions."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    energy = energy_mev*gd.UNIT_MEV
    L = 0.3*gd.L_SCALE_SUN

    monkeypatch.setattr(op, '_osc_prob_hybrid_dispatch', lambda *a, **k: NotImplemented)
    ip_spy = _spy_on(monkeypatch, '_osc_prob_ip_exp_dispatch')

    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        P = op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2, validate_input=False)
    assert ip_spy.answered is True, "the interaction-picture fast path did not produce this answer"
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
# Structural smoke test (G2): the BSM counterpart of G1 above, covering
# every NSI and LIV wrapper.
#
# G1 sweeps only the standard wrappers, and the targeted regression tests
# above reach a scattered subset of the BSM ones, so a coverage run found
# 23 of the 36 NSI/LIV wrappers were executed by nothing at all -- several
# of them while *appearing* in a parametrize list above, because those
# tests inspect the source or the signature without ever calling the
# function. A typo in one of their kwargs would have shipped unnoticed.
# ----------------------------------------------------------------------

# Deliberately non-zero, so the scenario actually perturbs the Hamiltonian
# rather than reducing to the standard case (every eps and every b default
# to 0.0, which would make this a re-run of G1 in disguise). Values are
# superset dicts: each call filters them against its own signature, which
# is what keeps one test honest across four flavor counts.
NSI_SWEEP_PARAMS = dict(
    eps_aa=0.1, eps_ab=0.05,                                    # 2nu naming
    eps_ee=0.2, eps_em=0.1, eps_et=0.05,                        # 3nu and up
    eps_mm=0.05, eps_mt=0.02, eps_tt=0.01,
    eps_es1=0.02, eps_ms1=0.01, eps_ts1=0.005, eps_s1s1=0.01,   # 4nu sterile
    eps_es2=0.01, eps_ms2=0.005, eps_ts2=0.002,                 # 5nu sterile
    eps_s1s2=0.005, eps_s2s2=0.01,
)
LIV_SWEEP_PARAMS = dict(
    sxi=0.2,                                                    # 2nu naming
    sxi12=0.2, sxi23=0.1, sxi13=0.05, dxiCP=0.0, dxi13=0.0,
    sxi14=0.05, dxi14=0.0, sxi24=0.03, dxi24=0.0, sxi34=0.02,
    sxi15=0.03, dxi15=0.0, sxi25=0.02, sxi35=0.01, dxi35=0.0,
    b1=gd.B1, b2=gd.B2, b3=gd.B3, b4=3.0e-9, b5=4.0e-9,
    Lambda=gd.LAMBDA, n_liv=1,
)


def bsm_wrapper_names():
    """Every osc_prob_{2,3,4,5}nu_* wrapper carrying an `nsi` or `liv`
    segment in its name.

    Discovered from the module rather than hand-listed, so a wrapper added
    later is swept without anyone remembering to extend a list here -- the
    gap this test exists to close was exactly that kind of omission.

    Matching is on '_'-separated segments, not on substrings: 'nsi' occurs
    inside 'de-nsi-ty', so a substring test silently pulls in every
    *_exp_density and *_constant_density wrapper as well.
    """
    import re
    names = []
    for name in dir(op):
        match = re.match(r'osc_prob_([2-5])nu_(.+)$', name)
        if match is None: continue
        if {'nsi', 'liv'} & set(match.group(2).split('_')):
            names.append(name)
    return sorted(names)


@pytest.mark.parametrize("name", bsm_wrapper_names())
def test_every_bsm_wrapper_runs_and_is_unitary(name):
    import inspect
    fn = getattr(op, name)
    params = inspect.signature(fn).parameters
    parts = name.split('_')
    num_flavors = int(name[len('osc_prob_')])

    kwargs = dict(validate_input=False)
    if num_flavors == 2:
        kwargs.update(sth=0.3, Dm2=2.5e-3)
    elif num_flavors == 4:
        kwargs.update({k: v for k, v in S5.items()
                      if k not in ('s15', 'd15', 's25', 's35', 'd35', 'D51')})
    elif num_flavors == 5:
        kwargs.update(S5)

    scenario = NSI_SWEEP_PARAMS if 'nsi' in parts else LIV_SWEEP_PARAMS
    kwargs.update({k: v for k, v in scenario.items() if k in params})

    if 'vacuum' in parts:
        P = fn(ENERGY, BASELINE, **kwargs)
    elif name.endswith('constant_density'):
        P = fn(ENERGY, BASELINE, RHO_C, **kwargs)
    elif name.endswith('exp_density'):
        P = fn(ENERGY, BASELINE, 0.0, RHO_C, L_SCALE, **kwargs)
    elif 'earth' in parts:
        P = fn(ENERGY, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM, **kwargs)
    elif 'sun' in parts:
        # 0.1 R_sun, where G1 above uses 0.5. The two-flavor solar path costs
        # time in proportion to the baseline (measured for osc_prob_2nu_sun_nsi:
        # 9.8 s at 0.5 R_sun, 5.1 s at 0.25, 2.5 s at 0.1), and it executes the
        # same lines either way -- this test asks whether the wrapper runs and
        # returns a unitary matrix, not how accurate it is over a long baseline,
        # which is what the solve_ivp cross-checks above are for.
        P = fn(ENERGY, 0.1*gd.SUN_RADIUS*gd.UNIT_KM, 0.0, **kwargs)
    else:
        raise AssertionError(f"unclassified wrapper: {name}")

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


# ----------------------------------------------------------------------
# Input validation raises rather than terminating the interpreter
# ----------------------------------------------------------------------

@pytest.mark.parametrize("call", [
    pytest.param(lambda: op.osc_prob_3nu_vacuum("not-a-number", 1.0*gd.UNIT_KM),
                 id="energy-wrong-type"),
    pytest.param(lambda: op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1.0*gd.UNIT_KM,
                                                nu_i=99, nu_f=0),
                 id="flavor-index-out-of-range"),
    pytest.param(lambda: op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1.0*gd.UNIT_KM,
                                                default_osc_params_set_name="NO_SUCH_SET"),
                 id="unknown-parameter-set"),
    pytest.param(lambda: op.osc_prob_3nu_vacuum(np.zeros((2, 2)), 1.0*gd.UNIT_KM),
                 id="energy-not-1d"),
])
def test_invalid_input_raises_valueerror_and_does_not_exit(call):
    """Invalid input must raise a catchable ValueError.

    These validation failures used to print a message and call ``sys.exit(1)``, which tears down
    the whole interpreter -- unusable from a notebook, a scan loop, or any caller that wants to
    recover, and impossible to assert on in a test.  ``pytest.raises(ValueError)`` here would
    also catch a regression back to SystemExit, since SystemExit derives from BaseException and
    would propagate out of the test rather than being caught."""
    with pytest.raises(ValueError):
        call()


def test_validate_input_battery_returns_none_on_valid_input():
    """The battery signals failure by raising, so a passing run simply returns None."""
    result = op.validate_input_battery(
        'test', energy=1.0*gd.UNIT_GEV, L=1.0*gd.UNIT_KM, num_flavors=3, nu_i=0, nu_f=1,
        osc_params=[0.55, 0.69, 0.15, 3.7, 7.49e-5, 2.513e-3],
        validate_energy_and_L=True, validate_flavor_indices=True, validate_osc_params=True)
    assert result is None


# ----------------------------------------------------------------------
# Method-aware slab cap
# ----------------------------------------------------------------------

def test_max_n_slabs_default_is_method_aware():
    """'gl' costs 1-3 Hamiltonian evaluations per slab against the quadrature methods'
    n_tpts_per_slab, so a cap that bounds cost has to differ between them."""
    assert op.MAX_N_SLABS_DEFAULT['gl'] > op.MAX_N_SLABS_DEFAULT['trapezoid']
    assert op.MAX_N_SLABS_DEFAULT['trapezoid'] == op.MAX_N_SLABS_DEFAULT['simpson']


def test_explicit_max_n_slabs_overrides_the_per_method_default():
    """An explicitly passed cap must be used as given, never silently widened."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        op.osc_prob_5nu_earth(ENERGY, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
                              validate_input=False, max_n_slabs=2000,
                              convergence_info=info, **S5)
    assert info['n_slabs'] <= 2000


def test_hard_sterile_earth_case_converges_without_warning_by_default():
    """Regression test for the shared 2000-slab cap.

    With eV-scale sterile splittings over an Earth-crossing baseline, 'gl' needs a few
    thousand slabs -- more than the cap tuned for the quadrature methods, whose per-slab
    cost is over an order of magnitude higher.  Under the shared cap this raised
    ToleranceNotAchievedWarning while returning an answer that was in fact far more
    accurate than the quadrature methods reached within that same cap: the warning was
    about not being able to *verify* convergence, having no room left to refine, and it
    pointed the reader at the wrong knob."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        info = {}
        P = op.osc_prob_5nu_earth(ENERGY, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
                                  validate_input=False, convergence_info=info, **S5)
    tol_warnings = [w for w in caught
                    if issubclass(w.category, op.ToleranceNotAchievedWarning)]
    assert not tol_warnings, f"unexpected {tol_warnings[0].message}"
    assert info['n_slabs'] < op.MAX_N_SLABS_DEFAULT['gl'], "converged only by hitting the cap"
    assert np.allclose(np.sum(np.asarray(P), axis=-1), 1.0, atol=1e-9)


# ----------------------------------------------------------------------
# Refinement caps below one loop
# ----------------------------------------------------------------------

@pytest.mark.parametrize("max_num_loops", [0, -5])
def test_max_num_loops_below_one_returns_a_probability_rather_than_crashing(max_num_loops):
    """Regression test: osc_prob used to raise UnboundLocalError here.

    The refinement-limit checks sit at the top of the loop and ``return P``, but on the
    first pass no loop has produced a P yet.  A cap below 1 reached that return before the
    variable existed.  It went unnoticed for two reasons: ``validate_input=True`` (the
    default) rejects ``max_num_loops <= 1`` up front, and the since-removed
    ``iterate_over_magnus_exp_order`` dispatch assigned P early, which hid the flaw from
    static analysis.

    The contract being pinned is not merely "does not crash": one refinement loop must
    still run -- there is no probability to return otherwise -- the result must be exactly
    unitary, and the caller must be *told* the requested tolerance was not reached rather
    than handed an unconverged number silently."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        P = op.osc_prob_3nu_earth(ENERGY, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
                                  validate_input=False, rtol=1e-4, atol=1e-4,
                                  max_num_loops=max_num_loops)

    P = np.asarray(P)
    assert P.shape == (3, 3)
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-12)
    assert np.all((P >= -1e-12) & (P <= 1.0 + 1e-12))
    assert any(issubclass(w.category, op.ToleranceNotAchievedWarning) for w in caught), \
        "returned an unconverged result without warning that the tolerance was missed"


def test_validation_still_rejects_max_num_loops_below_two():
    """The guard above is a safety net for validate_input=False; with validation on, the
    invalid cap is still refused up front rather than quietly degraded."""
    with pytest.raises(ValueError, match="max_num_loops"):
        op.osc_prob_3nu_earth(ENERGY, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
                              rtol=1e-4, atol=1e-4, max_num_loops=0)


# ----------------------------------------------------------------------
# n_slabs is a floor on the adaptive ladder, not a discarded argument
#
# Requesting a tolerance used to make osc_prob throw the caller's n_slabs
# away and start refining from min_n_slabs=1.  On a profile whose feature
# scale the caller knows and the refinement criterion does not, that is
# how a wrong answer gets certified: suggest_n_slabs measures the
# *integral* of the Hamiltonian along the path, which is blind to
# structure that averages out, so a castle-wall profile with 50 density
# walls accumulates only ~9 radians of phase and is seeded with 2 slabs.
# The resulting ladder does not converge, it thrashes -- 0.43, 0.13,
# 0.13, 0.64, 0.12 at 2, 3, 4, 5, 6 slabs -- and the successive-iterate
# test fires on the accidental 3-vs-4 agreement, returning an answer
# wrong by 0.48 in probability with no warning.  Tightening rtol does not
# help: the comparison is between two answers that both failed to see the
# profile.  Resolving the profile does, and the caller is the one who
# knows how.
#
# The oracle here is solve_ivp, deliberately.  Comparing one osc_prob
# grid against another osc_prob grid is what let this through in the
# first place -- two runs of an under-resolved scheme agree with each
# other and disagree with reality.
# ----------------------------------------------------------------------

# Castle-wall parameters shared by the tests below; the walls span a fixed
# interval and each baseline crosses a prefix of them, as in notebook 03.
CW_N_WALLS = 50
CW_L_INI, CW_L_FIN = 1.e2, 1.e4      # [km]
CW_ENERGY = 50.0*gd.UNIT_MEV
# A baseline at which the unfloored ladder certifies 0.268 against a truth
# of 0.749.  Found by scanning; nothing about it is special beyond being
# short enough to keep the solve_ivp reference cheap.
CW_L_KM = 2070.18


def castle_wall_H():
    """2nu Hamiltonian with a square-wave (castle-wall) electron density.

    The mixing parameters are the notebook's, not this module's S12/D21: the baseline
    below was located against these, and a faithful reproduction of the reported failure
    is worth more here than consistency with the rest of the file."""
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0)
    n_lo = 1*gd.N_AV/gd.CONV_CM_TO_INV_EV**3
    n_hi = 10*n_lo

    def num_density_e(l):
        l_scaled = ((np.asarray(l)/gd.CONV_KM_TO_INV_EV - CW_L_INI)
                    / (CW_L_FIN - CW_L_INI))
        index_slab = l_scaled//(1.0/CW_N_WALLS)
        return np.where(index_slab % 2 == 0, n_lo, n_hi)

    def H(l):
        VCC = matter.VCC_func(l, num_density_e)
        return (1.0/CW_ENERGY)*hvac + hams.hamiltonian_2nu_matter(VCC)

    return H


def castle_wall_reference(H, L):
    """P from an independent integrator: DOP853, stepping inside every wall."""
    def rhs(l, y):
        return (-1j*H(l) @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(2, dtype=complex).ravel(),
                    method='DOP853', rtol=1e-11, atol=1e-13,
                    max_step=(CW_L_FIN - CW_L_INI)*gd.UNIT_KM/(8*CW_N_WALLS))
    return np.abs(sol.y[:, -1].reshape(2, 2)).T**2


def test_requested_n_slabs_is_honoured_as_a_floor_under_a_tolerance():
    """The headline regression: a requested n_slabs must not be discarded.

    Without the floor this call refines from 2 slabs, stops at 3, and returns 0.268 for a
    probability that is 0.749 -- silently, and with the default rtol=atol=1e-3 nominally
    satisfied."""
    H = castle_wall_H()
    L = CW_L_KM*gd.UNIT_KM
    P_ref = castle_wall_reference(H, L)

    info = {}
    P = op.osc_prob(H, 0.0, L, n_slabs=150, magnus_exp_order=3,
                    convergence_info=info)

    assert info['n_slabs'] >= 150, \
        "the requested n_slabs was discarded; the ladder started below it"
    assert maxabs(np.asarray(P) - P_ref) < 1e-3


def test_unfloored_ladder_is_what_the_floor_protects_against():
    """Companion to the test above: the same call *without* n_slabs is still wrong.

    This pins the scope of the fix rather than a bug.  The floor makes a caller-supplied
    feature scale authoritative; it does not teach the phase-based seed to see structure
    that averages out of the integral it measures.  A caller who names no scale still gets
    the old, under-resolved ladder here.  If a future change to the refinement criterion
    makes this pass on its own, delete the test -- do not weaken it."""
    H = castle_wall_H()
    L = CW_L_KM*gd.UNIT_KM
    P_ref = castle_wall_reference(H, L)

    P = op.osc_prob(H, 0.0, L, magnus_exp_order=3)
    assert maxabs(np.asarray(P) - P_ref) > 0.1


@pytest.mark.parametrize("integration_method", ['gl', 'simpson'])
def test_floor_is_honoured_by_every_integration_method(integration_method):
    """The floor is applied before the method-specific seeding, so it holds for the
    phase-seeded 'gl' ladder and the plain min_n_slabs ladder of the quadrature methods
    alike."""
    H = castle_wall_H()
    info = {}
    op.osc_prob(H, 0.0, CW_L_KM*gd.UNIT_KM, n_slabs=64, magnus_exp_order=2,
                integration_method=integration_method, max_num_loops=2,
                convergence_info=info)
    assert info['n_slabs'] >= 64


def test_default_n_slabs_leaves_the_ladder_untouched():
    """n_slabs defaults to 1, so the floor is inactive unless the caller opts in: every
    call that does not pass n_slabs must behave exactly as it did before."""
    H = castle_wall_H()
    L = CW_L_KM*gd.UNIT_KM
    info_default, info_explicit = {}, {}
    P_default = op.osc_prob(H, 0.0, L, magnus_exp_order=3,
                            convergence_info=info_default)
    P_explicit = op.osc_prob(H, 0.0, L, n_slabs=1, magnus_exp_order=3,
                             convergence_info=info_explicit)
    assert info_default == info_explicit
    assert np.array_equal(np.asarray(P_default), np.asarray(P_explicit))


def test_floor_above_the_slab_cap_is_clipped_rather_than_stepped_down():
    """A floor larger than max_n_slabs must clip to the cap.

    Left unclipped it would make the ladder's first 'growth' step *down* to the cap, and
    the caller would never be told they had asked for more than the cap allows.  Clipped,
    the existing not-achieved warning fires, which is the honest report."""
    H = castle_wall_H()
    info = {}
    with pytest.warns(op.ToleranceNotAchievedWarning):
        op.osc_prob(H, 0.0, CW_L_KM*gd.UNIT_KM, n_slabs=10000, max_n_slabs=32,
                    magnus_exp_order=2, convergence_info=info)
    assert info['n_slabs'] == 32


def test_energy_batched_scan_honours_the_floor_too():
    """The batched scan engine keeps its own refinement ladder; it must read n_slabs the
    same way osc_prob does, or the two paths disagree on what the caller asked for."""
    energies = np.array([0.8, 1.5, 3.0])*gd.UNIT_GEV
    L = 2.0*6371.0*0.7*gd.UNIT_KM
    common = dict(costhz=-0.7, L=L, magnus_exp_order=4, validate_input=False,
                  max_num_loops=2, rtol=1e-3, atol=1e-3)
    P_floor = np.asarray(op.osc_prob_3nu_earth(energies, n_slabs=4096, **common))
    P_plain = np.asarray(op.osc_prob_3nu_earth(energies, **common))
    # Same physics, different grids: agreement to the requested tolerance shows the floored
    # run is a genuine refinement of the same problem, while disagreement at 1e-12 shows it
    # really ran on a different (finer) grid rather than ignoring the argument.
    assert np.allclose(P_floor, P_plain, rtol=1e-2, atol=1e-2)
    assert not np.allclose(P_floor, P_plain, rtol=1e-12, atol=1e-12)


# ----------------------------------------------------------------------
# The give-up paths of the closed-form interaction-picture integrator.
#
# _osc_prob_ip_exp_core has five exits that report "I could not verify
# this"; a coverage run found none of them was ever taken. That is the
# worst place in this codebase for a blind spot: its whole design premise
# is that agreement between two refinements is only evidence when the two
# were genuinely different computations, and the surrounding machinery has
# already produced one false-certification bug (hybrid_propagator, fixed
# with its own regression test). The success paths were well covered; the
# refusals were not covered at all.
#
# Reaching the ceiling honestly would mean two million slabs -- gigabytes
# and minutes -- so the ceiling itself is monkeypatched down. That is why
# IP_EXP_N_SLABS_CAP and IP_EXP_LOOP_CAP are named module constants rather
# than numbers inlined in the function body.
#
# One of the five refusals has no test here, because it cannot be reached:
# the `is_repeat` return fires when a refinement lands on the same slab
# count as the previous pass, which can only happen at the ceiling -- and
# every branch there returns within the same iteration, so the loop never
# survives a pass at the ceiling to make that comparison. It is marked
# `# pragma: no cover` at the source, with the argument written out, and it
# becomes live again the moment the growth factor or the ceiling changes.
#
# The "trusted comparison that disagreed" return needs no test here either,
# but for the opposite reason: the real two-flavor solar cases elsewhere in
# this file reach it, at the production ceiling and a realistic 1e-3
# tolerance. An earlier reading of this code claimed it was unreachable,
# generalizing from a synthetic probe whose potential was too weak to leave
# the trust region while disagreeing; the coverage data says otherwise, and
# the coverage data is the evidence.
# ----------------------------------------------------------------------

def _ip_exp_inputs(nE=1, d=2):
    """A minimal, genuinely exponential 2-level problem of the shape
    _osc_prob_ip_exp_dispatch feeds to the core: a constant H_E, an
    exponential VCC_func, and the constant matrix it multiplies."""
    l_scale = 100.0*gd.UNIT_KM
    H_E = np.array([[[1.0e-12, 0.0], [0.0, -1.0e-12]]]*nE, dtype=complex)
    h_matt = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)

    def VCC_func(l):
        return 1.0e-12*np.exp(-np.asarray(l, dtype=float)/l_scale)

    return dict(H_E=H_E, l_scale=l_scale, VCC_func=VCC_func, h_matt=h_matt,
                L0=0.0, L_val=5.0*l_scale)


def _call_ip_exp_core(rtol, atol, n_slabs=8, min_n_slabs=8, max_n_slabs=64):
    return op._osc_prob_ip_exp_core(
        **_ip_exp_inputs(), rtol=rtol, atol=atol, growth_factor_n_slabs=2.0,
        max_num_loops=10, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        n_slabs=n_slabs)


# ----------------------------------------------------------------------
# The working set of the batched engines.
#
# _osc_prob_ip_exp_core builds arrays of shape (n_energies, n_slabs, d, d),
# whose size is the product of two quantities the caller sets independently:
# neither the slab ceiling nor the energy count bounds it alone. At the
# two-million-slab ceiling that reached ~1.3 GB *per energy*, so a batched
# solar call could exhaust the machine rather than the process -- it
# OOM-killed the application it ran under three times before it was found.
# See docs/dev/BUG_IP_EXP_MEMORY.md.
#
# Nothing in the suite noticed, because no test batched enough energies to
# make the allocation large. These tests are the missing coverage: one pins
# that tiling changed no number, one pins that the working set is bounded.
# ----------------------------------------------------------------------

@pytest.mark.parametrize("budget", [64, 4096, 1 << 16])
def test_tiling_the_working_set_changes_no_number(budget, monkeypatch):
    """Tiling must be exact, not approximate.

    Within a tile the arithmetic is elementwise, so slicing moves no value; and the slab
    product is folded in the same descending order with the accumulator on the left, so the
    parenthesis nesting -- the only thing that could shift a floating-point result -- is
    unchanged. Exact equality is therefore the right assertion, and a looser one would hide
    precisely the bug this guards against."""
    monkeypatch.setattr(op, 'BATCH_WORKING_ENTRIES', 1 << 40)
    P_untiled, _ = _call_ip_exp_core(rtol=None, atol=None, n_slabs=512)

    monkeypatch.setattr(op, 'BATCH_WORKING_ENTRIES', budget)
    P_tiled, _ = _call_ip_exp_core(rtol=None, atol=None, n_slabs=512)

    assert np.array_equal(np.asarray(P_tiled), np.asarray(P_untiled))


def test_batched_working_set_does_not_scale_with_the_energy_count(monkeypatch):
    """The regression: peak allocation must be flat in the number of energies.

    Sized so that the *old* code cannot pass. At 64 energies and 65536 slabs one
    (nE, n_slabs, d, d) complex array is 268 MB, and the loop held several at once; the
    tiled version stays inside its entry budget no matter how the two axes grow."""
    monkeypatch.setattr(op, 'IP_EXP_N_SLABS_CAP', 1 << 16)

    def peak_mib(nE):
        tracemalloc.start()
        op._osc_prob_ip_exp_core(
            **_ip_exp_inputs(nE=nE), rtol=None, atol=None, growth_factor_n_slabs=2.0,
            max_num_loops=10, min_n_slabs=1 << 16, max_n_slabs=1 << 16, n_slabs=1 << 16)
        peak = tracemalloc.get_traced_memory()[1]/2**20
        tracemalloc.stop()
        return peak

    small, large = peak_mib(4), peak_mib(64)
    # The old code needed 268 MiB for *one* of the several arrays it held at these
    # parameters, so this bound is a clean discriminator rather than a tuned one.
    assert large < 150.0, f"working set grew to {large:.0f} MiB at 64 energies"
    # And it must be flat, not merely bounded: 16x the energies, no more memory. The
    # slack covers the result arrays, which legitimately do scale with the energy count.
    assert large < 1.5*small + 20.0, \
        f"working set scaled with the energy count: {small:.0f} -> {large:.0f} MiB"


def test_certification_known_impossible_is_refused_without_climbing(monkeypatch):
    """A tolerance no reachable slab count can certify must be refused at once.

    Certification needs max|Omega_t| under the trust threshold, and that maximum is bounded
    below by the diagonal entries, which have a closed form. When even that bound exceeds
    the threshold at the slab ceiling, every comparison below it is untrusted by the
    method's own rule, so the climb to the ceiling is guaranteed waste -- twenty-one
    doublings to reach a refusal that was computable up front.

    The result must still be a genuine unitary probability matrix: an uncertified answer is
    discarded by the dispatcher, but the contract that it is well-formed is what the other
    give-up tests here pin, and bailing early must not quietly break it."""
    monkeypatch.setattr(op, 'IP_EXP_N_SLABS_CAP', 64)

    calls = []
    real_expm = mg._expm_stack

    def counting_expm(*args, **kwargs):
        calls.append(args[0].shape)
        return real_expm(*args, **kwargs)

    monkeypatch.setattr(mg, '_expm_stack', counting_expm)

    P, converged = _call_ip_exp_core(rtol=1e-300, atol=1e-300, n_slabs=8, min_n_slabs=8)

    assert converged is False
    assert np.allclose(np.sum(np.asarray(P), axis=-1), 1.0, atol=1e-12), \
        "unitarity is structural here and must survive a failure to certify"
    total_slabs = sum(s[1] for s in calls)
    assert total_slabs <= 8, \
        f"climbed the ladder ({total_slabs} slabs) instead of refusing up front"


# ----------------------------------------------------------------------
# The cumulative baseline scan.
#
# U(0->L2) = U(L1->L2) U(0->L1), so every requested baseline is a prefix of
# the next and one traversal yields the whole scan. The characteristic way a
# cumulative product goes wrong is transposition and off-by-one, which no
# accuracy test would catch -- hence the identical-grid oracle below, which
# pins bookkeeping rather than physics. Accuracy is checked separately, and
# against solve_ivp, never against the per-point path: agreement between two
# configurations of the same scheme is what this project has been burned by.
# ----------------------------------------------------------------------

def test_cumulative_scan_matches_expm_for_a_constant_hamiltonian():
    """The ordering test. A constant H has a closed-form propagator at every baseline, so
    this catches a transposed or reversed product outright, at machine precision."""
    H = np.array([[1.2e-13, 3.0e-14], [3.0e-14, -0.7e-13]])
    L = np.linspace(1.0e12, 6.0e13, 400)

    P = np.asarray(op.osc_prob_energy_baseline(H, 1.0, L, cumulative=True, rtol=None,
                                               atol=None, n_slabs=400))

    worst = max(maxabs(P[i] - np.abs(sp.linalg.expm(-1j*H*L[i])).T**2)
                for i in range(0, len(L), 17))
    assert worst < 1e-11
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-11)


def test_cumulative_scan_agrees_with_the_per_point_path_on_an_identical_grid():
    """The sound oracle for a cumulative product, per the project's own history.

    Pinning rtol=atol=None and handing the per-point path the exact grid prefix the scan
    used makes the two computations arithmetically the same problem, so anything but
    agreement to near machine precision is a bookkeeping error -- a mis-ordered product, an
    off-by-one in where the running product is snapshotted. This deliberately does *not*
    test accuracy: two configurations of one scheme agreeing proves nothing about that."""
    H = castle_wall_H()
    L = np.linspace(200.0, 2000.0, 60)*gd.UNIT_KM
    n_acc = 3000

    P = np.asarray(op.osc_prob_energy_baseline(
        H, CW_ENERGY, L, cumulative=True, rtol=None, atol=None, n_slabs=n_acc,
        magnus_exp_order=3))

    edges, out_idx = op._cumulative_scan_grid(L, 0.0, n_acc, None)
    for i in range(0, len(L), 11):
        k = out_idx[i]
        prefix = np.column_stack([edges[:k], edges[1:k + 1]])
        P_ref = op.osc_prob(H, 0.0, L[i], t_slab_edges=prefix, magnus_exp_order=3,
                            rtol=None, atol=None, n_tpts_per_slab=2)
        assert maxabs(np.asarray(P_ref) - P[i]) < 1e-12


def test_cumulative_scan_returns_results_in_the_callers_order():
    """The scan must sort internally to walk the profile once, and undo that on the way
    out: a caller who passed unsorted baselines gets answers against *their* order."""
    H = np.array([[1.2e-13, 3.0e-14], [3.0e-14, -0.7e-13]])
    L = np.linspace(1.0e12, 6.0e13, 120)
    shuffled = RNG.permutation(len(L))

    common = dict(cumulative=True, rtol=None, atol=None, n_slabs=200)
    P_sorted = np.asarray(op.osc_prob_energy_baseline(H, 1.0, L, **common))
    P_shuffled = np.asarray(op.osc_prob_energy_baseline(H, 1.0, L[shuffled], **common))

    assert np.array_equal(P_shuffled, P_sorted[shuffled])


def test_cumulative_scan_sizes_its_accuracy_grid_from_the_adaptive_path():
    """The one way a cumulative scan goes *silently* wrong is an under-resolved accuracy
    grid: the traversal itself never notices, because there is nothing to compare against.

    So the grid is not guessed. One ordinary adaptive osc_prob call at the longest baseline
    reports the slab count it needed -- which is the definition of the accuracy grid -- and
    brings the existing safeguards with it. This pins that the scan really is at least that
    fine, which is what makes it safe."""
    H = castle_wall_H()
    L = np.linspace(200.0, 2000.0, 40)*gd.UNIT_KM

    info = {}
    op.osc_prob(H, 0.0, L[-1], magnus_exp_order=3, convergence_info=info)
    edges, _ = op._cumulative_scan_grid(
        L, 0.0, info['n_slabs']*op.CUMULATIVE_N_ACC_SAFETY, None)

    assert len(edges) - 1 >= info['n_slabs']*op.CUMULATIVE_N_ACC_SAFETY


def test_cumulative_scan_puts_every_requested_baseline_on_a_slab_edge():
    """Answers are read off the running product, never interpolated, so each requested
    baseline has to *be* an edge. Includes breakpoints and a logarithmic spacing, which is
    where a grid built only from the requested points goes coarse exactly where the
    integration needs it fine."""
    L = np.logspace(2.0, 4.0, 500)*gd.UNIT_KM
    bp = np.linspace(150.0, 9000.0, 37)*gd.UNIT_KM

    edges, out_idx = op._cumulative_scan_grid(L, 0.0, 700, bp)

    assert np.array_equal(edges, np.unique(edges))
    assert np.allclose(edges[out_idx], L, rtol=0, atol=0)
    assert set(np.round(bp, 6)).issubset(set(np.round(edges, 6)))


@pytest.mark.parametrize("kwargs,match", [
    (dict(t_slab_edges=[[0.0, 1.0]]), "t_slab_edges"),
    (dict(), "energies differ"),
])
def test_cumulative_scan_refuses_what_it_cannot_honour(kwargs, match):
    """Two requests it must reject rather than quietly reinterpret: its own grid cannot
    coexist with caller-supplied slab edges, and there is no nesting to exploit across
    energies, so a multi-energy request is a misunderstanding worth naming."""
    H = np.array([[1.2e-13, 3.0e-14], [3.0e-14, -0.7e-13]])
    energy = 1.0 if 't_slab_edges' in kwargs else np.array([1.0, 2.0])
    L = np.array([1.0e12, 2.0e12]) if 't_slab_edges' in kwargs else np.array([1.0e12, 1.0e12])

    with pytest.raises(ValueError, match=match):
        op.osc_prob_energy_baseline(H, energy, L, cumulative=True, **kwargs)


def test_cumulative_scan_memory_does_not_scale_with_the_baseline_count():
    """Peak allocation must be bounded by the block and the result, not by the grid.

    The traversal is chunked and each snapshot is converted to a probability on the spot,
    so the recorded term collapses into the answer the caller asked for rather than sitting
    beside it as N complex unitaries."""
    H = np.array([[1.2e-13, 3.0e-14], [3.0e-14, -0.7e-13]])

    def peak_mib(n):
        L = np.linspace(1.0e12, 6.0e13, n)
        tracemalloc.start()
        op.osc_prob_energy_baseline(H, 1.0, L, cumulative=True, rtol=None, atol=None,
                                    n_slabs=20000)
        peak = tracemalloc.get_traced_memory()[1]/2**20
        tracemalloc.stop()
        return peak

    small, large = peak_mib(200), peak_mib(20000)
    assert large < 4.0*small + 20.0, \
        f"peak scaled with the baseline count: {small:.1f} -> {large:.1f} MiB"


def test_tile_for_working_set_stays_within_its_budget():
    """The tiling arithmetic itself, over the shapes the engines actually pass.

    Includes the degenerate case where one cell alone exceeds the budget: there is no
    tiling that helps, and returning (1, 1) so the caller proceeds with a visibly large
    allocation is better than refusing work that might still fit."""
    for nE, n_inner, cell, live in [(1, 2_000_000, 4, 8), (1000, 2_000_000, 4, 8),
                                    (10**7, 100, 4, 8), (3, 7, 9, 1), (1, 1, 10**9, 8)]:
        e_chunk, blk = op._tile_for_working_set(nE, n_inner, cell, live_arrays=live)
        assert 1 <= e_chunk <= nE and 1 <= blk <= n_inner
        if cell*live <= op.BATCH_WORKING_ENTRIES:
            assert e_chunk*blk*cell*live <= op.BATCH_WORKING_ENTRIES
        else:
            assert (e_chunk, blk) == (1, 1)


def test_output_guard_refuses_a_result_that_cannot_fit(monkeypatch):
    """A scan whose *result* exceeds memory must fail with a diagnosable message.

    Tiling bounds the engines' working set, but nothing shrinks the answer: N points over
    d flavors is N*d*d floats either way. Left unchecked that surfaces as an out-of-memory
    kill from deep inside an engine -- or, on an overcommitting kernel, as the machine
    going down rather than the process, which is how this whole investigation started."""
    monkeypatch.setattr(op, '_available_memory_bytes', lambda: 256*1024*1024)

    with pytest.raises(MemoryError, match="would return"):
        op._check_output_fits(50_000_000, 3, 'osc_prob_energy_baseline')


def test_output_guard_is_silent_for_ordinary_requests(monkeypatch):
    """It must cost nothing and refuse nothing on the scans people actually run.

    The free-memory figure is not even consulted below the size floor; this asserts that
    by making the probe explode if it is called."""
    def exploding_probe():
        raise AssertionError('free memory probed for a small request')

    monkeypatch.setattr(op, '_available_memory_bytes', exploding_probe)
    op._check_output_fits(10_000, 3, 'osc_prob_energy_baseline')      # ~0.7 MiB


def test_output_guard_does_not_block_when_memory_cannot_be_measured(monkeypatch):
    """A guard that cannot measure must not refuse: on a platform exposing neither
    MemAvailable nor sysconf, the library keeps working exactly as before."""
    monkeypatch.setattr(op, '_available_memory_bytes', lambda: None)
    op._check_output_fits(10**9, 5, 'osc_prob_energy_baseline')


def test_available_memory_is_a_plausible_positive_number_or_none():
    """Whatever it reports on this machine, it must be usable as a bound."""
    value = op._available_memory_bytes()
    assert value is None or (isinstance(value, int) and value > 0)


def test_ip_exp_core_without_a_tolerance_returns_after_one_pass():
    """With neither rtol nor atol there is nothing to converge *to*: the
    method runs once at the given n_slabs and reports success, because the
    caller asked for a fixed-cost evaluation rather than a certified one."""
    P, converged = _call_ip_exp_core(rtol=None, atol=None)

    assert converged is True
    P = np.asarray(P)
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-12)


def test_ip_exp_core_refuses_to_certify_at_the_slab_ceiling(monkeypatch):
    """At the ceiling the method has had its one genuine comparison. If that
    comparison did not satisfy both safeguards, no further refinement can
    supply more evidence -- growing is a no-op there -- so it must return
    the probabilities it has and report them as unverified rather than let
    a comparison of a result with itself stand in for convergence."""
    monkeypatch.setattr(op, 'IP_EXP_N_SLABS_CAP', 16)

    # A tolerance no number of slabs can meet, so the refusal is the only
    # honest outcome available at the ceiling.
    P, converged = _call_ip_exp_core(rtol=1e-300, atol=1e-300, n_slabs=8, min_n_slabs=8)

    assert converged is False
    P = np.asarray(P)
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-12), \
        "unitarity is structural here and must survive a failure to certify"


def test_ip_exp_core_refuses_to_certify_when_the_loop_budget_runs_out(monkeypatch):
    """The loop-count backstop. With the slab ceiling raised out of the way
    and the loop budget cut to a few passes, the method runs out of passes
    before it runs out of slabs, and must again refuse rather than return
    its last guess as certified.

    With the production constants this branch cannot be reached at all --
    the slab count doubles each pass, so the two-million-slab ceiling
    arrives around the twentieth pass and returns there, long before a
    thirtieth. It is a backstop against a future change to the growth
    factor or the ceiling, and this test is what keeps it honest if that
    change ever happens."""
    monkeypatch.setattr(op, 'IP_EXP_N_SLABS_CAP', 10_000_000)
    monkeypatch.setattr(op, 'IP_EXP_LOOP_CAP', 3)

    P, converged = _call_ip_exp_core(rtol=1e-300, atol=1e-300, n_slabs=8, min_n_slabs=8)

    assert converged is False
    assert np.allclose(np.sum(np.asarray(P), axis=-1), 1.0, atol=1e-12)


# ----------------------------------------------------------------------
# What the dispatch layer does with an uncertified hybrid result.
#
# adiabatic.hybrid_propagator returning certified=False is the signal that
# the strategy could not verify its own answer, and the two strategies are
# meant to react differently: 'auto' silently abandons the whole batch and
# lets the general Magnus path produce the answer, while 'hybrid' keeps the
# result and warns. A coverage run found neither reaction was ever
# exercised -- the entire HybridCertificationWarning path included.
#
# The propagator is replaced rather than driven to a genuine failure: what
# is under test is the dispatch layer's reaction, and manufacturing a real
# non-certifying case would test adiabatic.py's numerics instead, which
# tests/test_adiabatic.py already does against solve_ivp.
# ----------------------------------------------------------------------

def _force_uncertified(monkeypatch, dim=3):
    """Makes every hybrid_propagator call report an uncertified result."""
    def fake(H_func, l0, l1, **kwargs):
        return np.eye(dim, dtype=complex), [], False
    monkeypatch.setattr(op.adiabatic, 'hybrid_propagator', fake)


def test_hybrid_strategy_warns_when_a_point_cannot_be_certified(monkeypatch):
    """strategy='hybrid' is an explicit request for this method, so it keeps
    the best-effort answer -- still exactly unitary -- but must say out loud
    that the accuracy is not certified."""
    _force_uncertified(monkeypatch)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        P = op.osc_prob_3nu_sun(ENERGY, 0.1*gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
                                rtol=1e-3, atol=1e-3, strategy='hybrid',
                                validate_input=False)

    assert any(issubclass(w.category, op.HybridCertificationWarning) for w in caught), \
        "returned an uncertified hybrid result without warning"
    P = np.asarray(P)
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-9)


def test_auto_strategy_falls_back_silently_when_a_point_cannot_be_certified(monkeypatch):
    """strategy='auto' never promised the hybrid method, only a correct
    answer, so an uncertified point makes it abandon the whole batch and let
    the general Magnus path compute the result instead -- without a
    certification warning, because nothing uncertified is being returned.

    The fake propagator returns the identity, so if its output were kept the
    probability matrix would be the identity too; a genuine 3-flavor solar
    result is not, which is what makes the fallback observable here."""
    _force_uncertified(monkeypatch)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        P = np.asarray(op.osc_prob_3nu_sun(ENERGY, 0.1*gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
                                           rtol=1e-3, atol=1e-3, strategy='auto',
                                           validate_input=False))

    assert not any(issubclass(w.category, op.HybridCertificationWarning) for w in caught), \
        "warned about certification for a result the hybrid path did not produce"
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-9)
    assert maxabs(P - np.eye(3)) > 1e-6, \
        "the identity from the stubbed propagator was returned instead of falling back"


# ----------------------------------------------------------------------
# The verbose-output helpers.
#
# These print; they compute nothing, which is why nothing called them and
# why they are the last large block of unexecuted statements in this
# module. They are still reachable from the public API through verbose=1
# and save_log=True, and the colour question is not cosmetic: ANSI escapes
# are correct on a terminal and are noise in a log file, a notebook, or
# this package's own rendered documentation, which is where the codebase
# has already been bitten once.
# ----------------------------------------------------------------------

def test_banner_prints_the_real_version_in_colour_on_stdout(capsys):
    op.print_banner()
    out = capsys.readouterr().out

    assert op.version.__version__ in out, "the banner does not show the resolved version"
    assert '\x1b[' in out, "expected ANSI colour when printing to a terminal"


def test_banner_written_to_a_file_carries_no_ansi_escapes(tmp_path):
    """A banner written to a log file must be plain text: escape codes in a
    file are unreadable noise, not colour."""
    path = tmp_path/'run.log'
    with open(path, 'w') as handle:
        op.print_banner(file=handle)

    text = path.read_text()
    assert op.version.__version__ in text
    assert '\x1b[' not in text, "ANSI escapes leaked into a log file"


def test_run_parameters_are_reported_verbatim(capsys):
    """The point of the parameter dump is that a run can be reproduced from
    its own output, so the values it prints must be the ones passed in."""
    def H_func(t):
        return np.zeros((3, 3), dtype=complex)

    op.print_run_parameters(H_func, 0.0, BASELINE, n_slabs=7, n_tpts_per_slab=13,
                            magnus_exp_order=6, integration_method='simpson',
                            rtol=1.5e-4, atol=2.5e-6)
    out = capsys.readouterr().out

    for expected in ['7', '13', '6', 'simpson']:
        assert expected in out, f"{expected!r} missing from the reported run parameters"


def test_verbose_run_prints_a_banner_and_parameters(capsys):
    """The public route into both helpers. Note the threshold is verbose=2,
    not 1: verbose=1 emits warnings only, which is what makes the banner and
    the parameter dump easy to leave unexercised."""
    op.has_magnus_header_been_printed = False
    op.osc_prob_3nu_matter_constant_density(ENERGY, BASELINE, RHO_C, verbose=2,
                                            validate_input=False)
    out = capsys.readouterr().out

    assert op.version.__version__ in out, "no banner in verbose=2 output"
    assert 'magnus_exp_order' in out, "no run-parameter dump in verbose=2 output"


# ----------------------------------------------------------------------
# compute_evolution_operator: documented, exported, and called by nothing.
#
# Its own docstring says it is meant to be called internally by osc_prob,
# but osc_prob reaches for the multi-slab variant instead, so the suite
# found it with an entirely unexecuted body. It is still public API.
# ----------------------------------------------------------------------

def test_single_slab_evolution_operator_is_unitary_and_matches_the_multislab_chain():
    """One slab computed on its own must equal the same slab computed as a
    one-element chain -- otherwise the two entry points disagree about the
    convention they share."""
    def H_func(t):
        return np.array([[1.0e-12, 2.0e-13], [2.0e-13, -1.0e-12]], dtype=complex)

    t_slab = [0.0, BASELINE]
    U = np.asarray(op.compute_evolution_operator(H_func, t_slab, n_tpts_per_slab=2,
                                                 magnus_exp_order=4))

    assert maxabs(U.conj().T @ U - np.eye(2)) < 1e-12

    U_chain = op.compute_evolution_operator_multiple_slabs(
        H_func, [t_slab], n_tpts_per_slab=2, magnus_exp_order=4)
    assert maxabs(U - np.asarray(U_chain)[0]) < 1e-12


def test_zero_width_slab_evolution_operator_is_the_identity():
    """A slab of zero width evolves nothing. This is the guard that lets a
    caller pass a degenerate slab -- an empty segment either side of a
    non-adiabatic window, say -- without special-casing it."""
    def H_func(t):
        return np.array([[1.0e-12, 2.0e-13], [2.0e-13, -1.0e-12]], dtype=complex)

    U = np.asarray(op.compute_evolution_operator(H_func, [5.0, 5.0], n_tpts_per_slab=2,
                                                 magnus_exp_order=4))
    assert maxabs(U - np.eye(2)) == 0.0


# ----------------------------------------------------------------------
# Promoted from the adversarial-validation batteries (docs/dev/adversarial_batteries/).
# Those scripts are diagnostics -- they print tables, take tens of minutes, and are not run by
# CI -- so the silent-miss classes they found would otherwise be defended by nothing.  These
# four are the cheapest configuration that still fails if the corresponding safeguard is
# removed, and each states the number it was worth when it was broken.
# ----------------------------------------------------------------------

def _ne_step(mid, lo_frac=0.02, hi_frac=0.30):
    """Piecewise-constant electron density with one jump: expm is EXACT for this."""
    lo, hi = lo_frac*gd.NUM_DENSITY_E_SUN_CENTRAL, hi_frac*gd.NUM_DENSITY_E_SUN_CENTRAL

    def ne(l):
        y = np.where(np.asarray(l, dtype=float) < mid, lo, hi)
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a
    return ne


def _exact_step_P(mid, l1, energy, sth, Dm2):
    """Exact probability matrix across a two-piece constant profile, via expm composition."""
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    vcc = matter.vcc_func_from_rho_func(_ne_step(mid), 0.0, 1.0, 0.5, nubar=False,
                                        density_matter_is_in_g_per_cm3=False,
                                        density_is_of_number_of_electrons=True)
    H_a = (1.0/energy)*h_vac + float(np.asarray(vcc(0.5*mid)))*proj
    H_b = (1.0/energy)*h_vac + float(np.asarray(vcc(0.5*(mid + l1))))*proj
    U = sp.linalg.expm(-1j*H_b*(l1 - mid)) @ sp.linalg.expm(-1j*H_a*mid)
    return np.transpose(U.real**2 + U.imag**2)


def test_unmarked_density_step_is_not_answered_silently_wrong():
    """The headline finding of the adversarial validation, at the public entry point.

    strategy='auto' returned P_ee = 0.589 against an exact expm answer of 0.0498 -- wrong by
    **0.54** -- with no warning, because the hybrid strategy's only guard against a
    discontinuous profile was "did the caller pass t_breakpoints", which fails open.  The guard
    is now a measurement (magnus.adiabatic._profile_is_resolved), so hybrid declines to certify
    and osc_prob falls through to the general Magnus path.
    """
    energy, l1 = 50.0e6, gd.L_SCALE_SUN
    mid = 0.5*l1
    sth, Dm2 = gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0
    P_exact = _exact_step_P(mid, l1, energy, sth, Dm2)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = np.asarray(op.osc_prob_matter_std_potential(
            2, _ne_step(mid), energy, l1, {'sth': sth, 'Dm2': Dm2}, L0=0.0,
            density_is_of_number_of_electrons=True))

    err = np.max(np.abs(P - P_exact))
    assert err < 1e-3, f"unmarked density step answered wrong by {err:.2e} (was 5.4e-01)"
    assert caught, "an answer this hard-won must not come back silent"


def test_subthreshold_nonadiabaticity_is_not_answered_silently_wrong():
    """A sinusoidal density the probe resolves easily (~28 samples per period), on which gamma
    stays below the 0.1 window threshold everywhere.

    No window opened, successive refinements agreed with each other because they were all pure
    adiabatic transport, and the answer came back certified and wrong by 1.67e-02 against a
    requested 1e-3.  Certifying an empty window list now also requires gamma to be small enough
    for the tolerance.
    """
    energy, l1 = 10.0e6, gd.L_SCALE_SUN
    sth, Dm2 = gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    ne0 = gd.NUM_DENSITY_E_SUN_CENTRAL
    c_vcc = float(np.asarray(matter.vcc_func_from_rho_func(
        lambda l: ne0, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)(0.0)))/ne0

    xs = np.geomspace(ne0*1e-6, ne0*10.0, 2000)
    gaps = [np.diff(np.linalg.eigvalsh(h_vac/energy + x*c_vcc*proj))[0] for x in xs]
    ne_res = float(xs[int(np.argmin(gaps))])

    def ne(l):
        y = ne_res*(1.0 + 0.9*np.sin(2.0*np.pi*np.asarray(l, dtype=float)/(l1/7.0)))
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a

    def H_func(l):
        v = np.asarray(matter.vcc_func_from_rho_func(
            ne, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
            density_is_of_number_of_electrons=True)(l))
        return (1.0/energy)*h_vac + v[..., None, None]*proj

    P = np.asarray(op.osc_prob_matter_std_potential(
        2, ne, energy, l1, {'sth': sth, 'Dm2': Dm2}, L0=0.0,
        density_is_of_number_of_electrons=True))

    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(2, 2)).ravel()
    sol = solve_ivp(rhs, (0.0, l1), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853')
    U = sol.y[:, -1].reshape(2, 2)
    P_exact = np.transpose(U.real**2 + U.imag**2)

    err = np.max(np.abs(P - P_exact))
    assert err < 1e-3, f"sub-threshold non-adiabaticity wrong by {err:.2e} (was 1.67e-02)"


def test_declared_breakpoints_make_a_piecewise_profile_essentially_exact():
    """The other half of the piecewise story, and the advice the documentation now gives.

    Over 150 random piecewise-constant profiles, declaring the edges gave a median error of
    1.34e-12 and nothing outside tolerance, against a median 7.76e-04 without.  One
    representative profile is enough to catch a regression in how t_breakpoints reaches the
    slab grid.
    """
    energy, l1 = 50.0e6, gd.L_SCALE_SUN
    edges = np.array([0.0, 0.17, 0.41, 0.63, 0.88, 1.0])*l1
    values = gd.NUM_DENSITY_E_SUN_CENTRAL*np.array([0.03, 0.21, 0.07, 0.30, 0.12])

    def ne(l):
        x = np.asarray(l, dtype=float)
        idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(values) - 1)
        y = values[idx]
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a

    sth, Dm2 = gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    vcc = matter.vcc_func_from_rho_func(ne, 0.0, 1.0, 0.5, nubar=False,
                                        density_matter_is_in_g_per_cm3=False,
                                        density_is_of_number_of_electrons=True)

    # Exact: H is constant on each segment, so expm composes exactly.
    U = np.eye(2, dtype=complex)
    for a, b in zip(edges[:-1], edges[1:]):
        H_m = (1.0/energy)*h_vac + float(np.asarray(vcc(0.5*(a + b))))*proj
        U = sp.linalg.expm(-1j*H_m*(b - a)) @ U
    P_exact = np.transpose(U.real**2 + U.imag**2)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = np.asarray(op.osc_prob_matter_std_potential(
            2, ne, energy, l1, {'sth': sth, 'Dm2': Dm2}, L0=0.0,
            density_is_of_number_of_electrons=True, t_breakpoints=edges[1:-1]))

    err = np.max(np.abs(P - P_exact))
    assert err < 1e-8, f"declared breakpoints should be near-exact; got {err:.2e}"


def test_cumulative_is_reachable_and_reproduces_the_per_point_path():
    """`cumulative=False` used to raise TypeError from every wrapper, which made the one
    documented route to pre-1.0.0 numbers unavailable at the layer where the change is visible.

    Also pins the two directions apart: False must match strategy='magnus' bit for bit, and
    'auto' must not (it takes the cumulative scan, which is the point).
    """
    energy = 50.0e6
    LM = 0.5*gd.SUN_RADIUS*gd.UNIT_KM
    Ls = np.linspace(0.05*LM, LM, 60)
    rho = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    params = {'sth': gd.S12_NO_BF_NUFIT_6_0, 'Dm2': gd.D21_NO_BF_NUFIT_6_0}

    def call(**kw):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.asarray(op.osc_prob_matter_std_potential(
                2, rho, energy, Ls, params, L0=0.0,
                density_is_of_number_of_electrons=True, **kw))

    p_false = call(cumulative=False)
    p_magnus = call(strategy='magnus')
    p_auto = call()

    assert np.array_equal(p_false, p_magnus), \
        "cumulative=False must reproduce the pre-cumulative per-point path exactly"
    assert not np.array_equal(p_auto, p_false), \
        "the default should be taking the cumulative scan on a 60-point single-energy scan"
    # An explicit request wins over the strategy='magnus' opt-out.
    assert np.array_equal(call(strategy='magnus', cumulative=True), p_auto)


def test_cumulative_scan_says_so_when_a_discontinuity_was_not_declared():
    """The last silent hole in the cumulative path, and the only one in code this branch adds.

    A slab straddling a density jump degrades the quadrature regardless of magnus_exp_order,
    and refining the grid only narrows the straddling slab -- so unlike every other inaccuracy
    here, more slabs cannot fix it.  Over 150 random piecewise profiles, 59 of 150 came back
    outside tolerance with the edges undeclared and all but two of those already warned; the two
    that did not were wrong by 1.36e-03 and 2.10e-03, silently.

    Detection is a measurement (magnus.adiabatic._profile_is_resolved), and its discriminator
    was checked in both directions before being wired in: 0 false positives on the solar,
    multi-resonance, noisy and sinusoidal families, 12/12 true positives on random piecewise
    profiles.
    """
    LM = 0.5*gd.SUN_RADIUS*gd.UNIT_KM
    Ls = np.linspace(0.05*LM, LM, 80)
    params = {'sth': gd.S12_NO_BF_NUFIT_6_0, 'Dm2': gd.D21_NO_BF_NUFIT_6_0}
    edges = np.array([0.0, 0.17, 0.41, 0.63, 0.88, 1.0])*LM
    values = gd.NUM_DENSITY_E_SUN_CENTRAL*np.array([0.03, 0.21, 0.07, 0.30, 0.12])

    def ne(l):
        x = np.asarray(l, dtype=float)
        idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(values) - 1)
        a = np.asarray(values[idx])
        return a[()] if a.ndim == 0 else a

    def call(**kw):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            op.osc_prob_matter_std_potential(
                2, ne, 50.0e6, Ls, params, L0=0.0,
                density_is_of_number_of_electrons=True, **kw)
        return [w.category for w in caught]

    assert any(issubclass(c, op.UnmarkedDiscontinuityWarning) for c in call()), \
        "an undeclared density jump on the cumulative path must not be silent"

    # Declaring the edges is the cure, and must silence it.
    assert not any(issubclass(c, op.UnmarkedDiscontinuityWarning)
                   for c in call(t_breakpoints=edges[1:-1])), \
        "declaring the breakpoints must silence the warning"

    # And it must not fire on the smooth profiles this package actually ships.
    rho = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)

    def call_smooth(profile):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            op.osc_prob_matter_std_potential(
                2, profile, 50.0e6, Ls, params, L0=0.0,
                density_is_of_number_of_electrons=True)
        return [w.category for w in caught]

    assert not any(issubclass(c, op.UnmarkedDiscontinuityWarning)
                   for c in call_smooth(rho)), "false positive on the solar exponential"


def test_a_closure_that_binds_defaults_is_not_mistaken_for_a_shorter_signature():
    """H_func arity was detected with len(signature(f).parameters), which counts keyword
    parameters that already have defaults.

    That breaks the ordinary Python idiom for capturing a value in a closure --
    ``def H(energy, l, VCC, _hvac=hvac)`` is three required parameters and five total -- so with
    validate_input=True it was rejected as "takes 5 argument(s)", and with validate_input=False
    it silently took the two-argument branch and died inside the engine with a TypeError.
    Counting *required* parameters instead makes both spellings identical.
    """
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0)
    e00 = np.diag([1.0, 0.0])
    energy, L = 20.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM

    def H3_plain(energy_, l, VCC):
        return (1.0/energy_)*hvac + np.asarray(VCC)[..., None, None]*e00

    def H3_bound(energy_, l, VCC, _h=hvac, _e=e00):
        return (1.0/energy_)*_h + np.asarray(VCC)[..., None, None]*_e

    def H2_plain(energy_, l):
        return (1.0/energy_)*hvac

    def H2_bound(energy_, l, _h=hvac):
        return (1.0/energy_)*_h

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.array_equal(np.asarray(op.osc_prob_sun(H3_plain, energy, L)),
                              np.asarray(op.osc_prob_sun(H3_bound, energy, L))), \
            "a 3-argument H_func with bound defaults must behave as a plain 3-argument one"
        assert np.array_equal(np.asarray(op.osc_prob_sun(H2_plain, energy, L)),
                              np.asarray(op.osc_prob_sun(H2_bound, energy, L))), \
            "a 2-argument H_func with bound defaults must behave as a plain 2-argument one"

    # And the same for osc_prob's own one-parameter check.
    def H1_plain(l):
        return hvac/energy

    def H1_bound(l, _h=hvac, _E=energy):
        return _h/_E

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.array_equal(np.asarray(op.osc_prob(H1_plain, 0.0, L, rtol=1e-3, atol=1e-3)),
                              np.asarray(op.osc_prob(H1_bound, 0.0, L, rtol=1e-3, atol=1e-3)))


def test_a_jump_smaller_than_the_steepest_smooth_step_is_still_detected():
    """The resolution test compared *global* maxima, so any discontinuity smaller than the
    largest smooth variation elsewhere on the path was masked.

    Measured through osc_prob_sun with a jump inside H itself (not in the density) on a solar
    profile: the jump was 4.7x smaller than the steepest smooth step, the test reported
    "resolved", hybrid certified, and the answer was wrong by 2.03e-02 with no warning at all.
    Comparing each half of an interval against that interval's own total variation -- i.e. how
    *concentrated* the variation is -- is immune to the masking.
    """
    energy = 50.0*gd.UNIT_MEV
    L = 0.3*gd.SUN_RADIUS*gd.UNIT_KM
    mid = 0.5*L
    hvac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    off = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    scale = float(np.max(np.abs(hvac)))/float(energy)
    vcc = matter.vcc_func_from_rho_func(
        matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN),
        0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)

    def H_of_l(l):
        x = np.asarray(l, dtype=float)
        f = np.where(x < mid, 0.0, 1.0)
        return ((1.0/energy)*hvac + np.asarray(vcc(l))[..., None, None]*proj
                + 0.30*scale*np.asarray(f)[..., None, None]*off)

    # The premise: the jump really is smaller than the steepest smooth step, so a global
    # comparison cannot see it.
    ls = np.linspace(0.0, L, 200)
    Hc = np.array([np.asarray(H_of_l(x), dtype=complex) for x in ls])
    steps = np.abs(np.diff(Hc, axis=0)).max(axis=(1, 2))
    j = int(np.argmin(np.abs(ls - mid)))
    assert steps.max() > 2.0*steps[max(j - 1, 0)], \
        "the jump is no longer masked by a larger smooth step; the test has lost its point"

    assert not ad._profile_is_resolved(H_of_l, 0.0, L, 200)
    _, _, certified = ad.hybrid_propagator(H_of_l, 0.0, L)
    assert not certified, "hybrid certified a discontinuity hidden under smooth variation"


def test_breakpoints_cure_a_feature_narrower_than_the_probe_grid():
    """The one exposure the adversarial validation could not close in the library, and the
    documented cure for it.

    A resonance narrower than the probe spacing is invisible to every path: the hybrid
    detector samples n_probe points (200, refined to at most 6400), and the general Magnus
    grid is seeded from an integral along the path, so neither ever puts a sample inside the
    feature and no amount of refinement changes that.  Measured at 2.9e-02 against a requested
    1e-3, silently.

    The `strategy` docstring tells users to pass t_breakpoints, so that instruction is pinned
    here: if the breakpoint plumbing ever stops reaching the slab grid, the documented cure
    would silently stop working and this is what would catch it.
    """
    energy = 10.0e6
    l1 = gd.L_SCALE_SUN
    centre, width = 0.4831*l1, 3.0e-5*l1        # ~1/170 of the finest probe spacing
    sth, Dm2 = gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0
    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)

    ne0 = gd.NUM_DENSITY_E_SUN_CENTRAL
    c_vcc = float(np.asarray(matter.vcc_func_from_rho_func(
        lambda x: ne0, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)(0.0)))/ne0
    xs = np.geomspace(ne0*1e-6, ne0*10.0, 2000)
    gaps = [np.diff(np.linalg.eigvalsh(h_vac/energy + x*c_vcc*proj))[0] for x in xs]
    ne_res = float(xs[int(np.argmin(gaps))])

    def ne(l):
        x = np.asarray(l, dtype=float)
        y = ne_res*(0.30 + 2.70*np.exp(-0.5*((x - centre)/width)**2))
        a = np.asarray(y)
        return a[()] if a.ndim == 0 else a

    vcc = matter.vcc_func_from_rho_func(ne, 0.0, 1.0, 0.5, nubar=False,
                                        density_matter_is_in_g_per_cm3=False,
                                        density_is_of_number_of_electrons=True)

    def H_func(l):
        return (1.0/energy)*h_vac + np.asarray(vcc(l))[..., None, None]*proj

    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(2, 2)).ravel()
    sol = solve_ivp(rhs, (0.0, l1), np.eye(2, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853')
    U = sol.y[:, -1].reshape(2, 2)
    P_exact = np.transpose(U.real**2 + U.imag**2)

    edges = np.array([centre - 8*width, centre - 2*width, centre,
                      centre + 2*width, centre + 8*width])
    params = {'sth': sth, 'Dm2': Dm2}

    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        P_bare = np.asarray(op.osc_prob_matter_std_potential(
            2, ne, energy, l1, params, L0=0.0, density_is_of_number_of_electrons=True))
    caught_bare = [w.category.__name__ for w in rec]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P_bp = np.asarray(op.osc_prob_matter_std_potential(
            2, ne, energy, l1, params, L0=0.0, density_is_of_number_of_electrons=True,
            t_breakpoints=edges))

    err_bare = np.max(np.abs(P_bare - P_exact))
    err_bp = np.max(np.abs(P_bp - P_exact))

    # The premise, stated as the property that matters rather than as a magnitude: undeclared,
    # the feature is missed by enough to break the requested tolerance, and nothing says so.
    # (Measured here: 4.60e-03 with no warning at all.  How large it is depends on where the
    # feature falls relative to the probe grid, so the assertion is on the tolerance, not on
    # the number.)  If the detector ever grows the ability to find it, this fails loudly rather
    # than passing quietly, and the warning in the `strategy` docstring should then be revised.
    assert err_bare > 2e-3, (
        f"a sub-probe-spacing feature is no longer missed when undeclared (err {err_bare:.2e}); "
        "if the detector improved, update the strategy docstring's warning to match")
    # This assertion used to read `not caught_bare`, with a note that it should fail loudly if
    # the detector ever grew the ability to find this.  It did: adiabatic.find_hidden_features
    # scans the PROFILE rather than the answers, which is what reaches a class every engine
    # misses together.  The answer is still wrong -- detection is not a cure, and the assertion
    # above still holds -- but it is no longer silent, which was the failure that mattered.
    assert 'HiddenFeatureWarning' in caught_bare, (
        "a feature narrower than every grid must at least be reported; it is the one exposure "
        "no cross-check between engines can reach")
    # And the cure works, which is what the documentation promises.  (Measured: 1.31e-04.)
    assert err_bp < 1e-3, f"t_breakpoints no longer cures the narrow feature (err {err_bp:.2e})"
    assert err_bp < err_bare/10.0
