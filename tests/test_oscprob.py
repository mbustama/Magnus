# -*- coding: utf-8 -*-
"""Tests of the oscillation-probability engine (magnus.oscprob)."""

import numpy as np
import pytest
import scipy as sp
from scipy.integrate import solve_ivp

import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.oscprob.oscprob as op
import magnus.oscprob.oscprobstd as opstd

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


def test_tolerance_cap_warns():
    """Hitting the refinement caps must warn even at verbose=0 (the result
    can look plausible while being unconverged)."""
    sth, Dm2 = np.sqrt(0.308), 7.5e-5
    with pytest.warns(op.ToleranceNotAchievedWarning):
        op.osc_prob_2nu_sun(10.0*gd.UNIT_MEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM,
                            0.0, sth, Dm2, integration_method='gl',
                            max_n_slabs=64, max_num_loops=8,
                            validate_input=False)


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
