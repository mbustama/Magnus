# -*- coding: utf-8 -*-
"""Tests of the phase average (magnus.avgprob, issue #64).

The phase average keeps every interference term with its phase at the central energy and weights
it by exp(-sigma^2 phi'^2 / 2), phi' = d phi / d ln E.  Its ground truths are:

* at sigma = 0, the oscillation probability itself (a flavor start for a constant Hamiltonian, a
  start decohered in the local eigenbasis on a profile);
* where every phase runs fast with energy, the decohered limit the older functions return;
* in between, the definition computed by brute force: the Gaussian average over u = delta ln E of
  the evolution under H + u D_diag, D_diag being the part of dH/dlnE diagonal in the
  instantaneous eigenbasis.

The last is what makes the answer independent of where the windows are drawn, which is the
property the step on the solar disk (issue #62) lacked: it is tested directly, one window against
several.
"""

import numpy as np
import pytest
import scipy.linalg

import magnus.adiabatic as ad
import magnus.avgprob as ap
import magnus.globaldefs as gd
import magnus.hamiltonians as hams

S12, S23, S13, DCP = 0.55, 0.68, 0.15, 3.7
D21, D31 = 7.5e-5, 2.5e-3


def maxabs(x):
    return float(np.max(np.abs(x)))


def vacuum_H(energy):
    """Three-flavor vacuum Hamiltonian and its exact derivative in ln E (it goes as 1/E)."""
    H = np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(S12, S23, S13, DCP, D21, D31),
                   dtype=complex)/energy
    return H, -H


# ------------------------------------------------------------------ constant Hamiltonian

@pytest.mark.parametrize("L_km", [1.0, 300.0, 1.0e4, 3.0e5])
def test_constant_H_at_zero_spread_is_the_oscillation_probability(L_km):
    H, D = vacuum_H(1.0*gd.UNIT_GEV)
    L = L_km*gd.UNIT_KM
    P, sens = ap.phase_averaged_probabilities_constant_hamiltonian(H, D, L, spread=0.0)
    S = scipy.linalg.expm(-1j*H*L)
    assert maxabs(P - np.abs(S.T)**2) < 1e-12
    assert sens == 0.0


def test_constant_H_accepts_a_batch_of_energies_and_baselines():
    E = np.array([0.5, 1.0, 2.0])*gd.UNIT_GEV
    L = np.array([100.0, 1000.0, 5000.0])*gd.UNIT_KM
    HD = [vacuum_H(e) for e in E]
    P, sens = ap.phase_averaged_probabilities_constant_hamiltonian(
        np.array([h for h, _ in HD]), np.array([d for _, d in HD]), L)
    for k in range(3):
        Pk, sk = ap.phase_averaged_probabilities_constant_hamiltonian(HD[k][0], HD[k][1], L[k])
        assert maxabs(P[k] - Pk) < 1e-15 and abs(sens[k] - sk) < 1e-15


def test_constant_H_decoheres_where_every_phase_runs_fast():
    H, D = vacuum_H(1.0*gd.UNIT_GEV)
    P, sens = ap.phase_averaged_probabilities_constant_hamiltonian(H, D, 1.0e8*gd.UNIT_KM)
    decohered = ap.averaged_probabilities_from_eigenbasis(np.linalg.eigh(H)[1])
    assert maxabs(P - decohered) < 1e-12
    assert sens < 1e-12


@pytest.mark.parametrize("phase", [0.5, 2*np.pi, 10.0, 20.0, 40.0])
def test_two_flavor_vacuum_damps_the_oscillation_by_the_gaussian_of_its_phase(phase):
    """P_ee = 1 - sin^2(2 theta)/2 (1 - exp(-sigma^2 phi^2/2) cos phi), phi = Dm2 L/2E."""
    theta, dm2, E, sigma = 0.6, 7.5e-5, 1.0e6, 0.1
    c, s = np.cos(theta), np.sin(theta)
    U = np.array([[c, s], [-s, c]])
    H = U @ np.diag([0.0, dm2/(2*E)]) @ U.T
    L = phase/(dm2/(2*E))
    P, _ = ap.phase_averaged_probabilities_constant_hamiltonian(H, -H, L, spread=sigma)
    expected = 1 - 0.5*np.sin(2*theta)**2*(1 - np.exp(-0.5*(sigma*phase)**2)*np.cos(phase))
    assert abs(P[0, 0] - expected) < 1e-12
    assert maxabs(P.sum(axis=1) - 1) < 1e-12


def test_a_round_off_pseudo_dirac_pair_stays_coherent():
    """Issue #61: at 100 TeV over 100 Mpc a pair split by 1e-21 eV^2 has a true phase of 8e-5 rad,
    but its eigenvalue gap, and a finite-difference slope, are round-off.  The slope floor keeps
    it coherent: the result is the coherent-block one, and does not depend on the spread."""
    U = hams.pmns_mixing_matrix(S12, S23, S13, DCP)
    m2 = np.array([0.0, D21, D31])
    E, L, h = 100.0*gd.UNIT_TEV, 100.0*3.0857e19*gd.UNIT_KM, 1.0e-3

    def H_at(e):
        return np.asarray(hams.hamiltonian_pseudo_dirac_vacuum(e, U, m2, {1: 1.0e-21}), dtype=complex)

    H = H_at(E)
    D = (H_at(E*np.exp(h)) - H_at(E*np.exp(-h)))/(2*h)
    P, sens = ap.phase_averaged_probabilities_constant_hamiltonian(H, D, L, dH_dlnE_step=h)
    blocks = ap.coherence_blocks(np.linalg.eigvalsh(H), L)
    coherent = ap.averaged_probabilities_from_eigenbasis(np.linalg.eigh(H)[1], blocks)
    assert maxabs(P - coherent) < 1e-5
    assert sens < 1e-6


def test_constant_H_rejects_what_it_cannot_average():
    H, D = vacuum_H(1.0*gd.UNIT_GEV)
    with pytest.raises(ValueError, match='square'):
        ap.phase_averaged_probabilities_constant_hamiltonian(H[:2], D[:2], 1.0)
    with pytest.raises(ValueError, match='spread'):
        ap.phase_averaged_probabilities_constant_hamiltonian(H, D, 1.0, spread=-0.1)


# ------------------------------------------------------------------ smooth profiles

E0 = 1.0


def crossing_H(g=0.2, delta=1.0, v0=3.0, l_scale=1.0):
    """Two levels crossing where v(l) = v0 exp(-l/l_scale) passes delta, coupled by g.

    The constant part stands for vacuum oscillations and goes as 1/E, the falling diagonal for
    matter and does not depend on energy: so dH/dlnE = -(H - diag(v, 0)) exactly.  With the
    defaults the crossing sits at l = ln 3 and hops with probability ~0.8 (Landau-Zener), and
    past it the levels accumulate about 1 rad per unit length -- ~20 rad by l = 22, enough that
    a 10 % spread leaves part of the interference.
    """
    A = np.array([[0.0, g], [g, delta]], dtype=complex)

    def H(l, energy=E0):
        v = v0*np.exp(-np.asarray(l, dtype=float)/l_scale)
        return (E0/energy)*A + np.asarray(v)[..., None, None]*np.diag([1.0, 0.0])

    def D(l):
        return -(H(l) - np.asarray(v0*np.exp(-np.asarray(l, dtype=float)/l_scale))[..., None, None]
                 * np.diag([1.0, 0.0]))

    return H, D, 0.0, 22.0


def decohered_start_probability(H, U, l0):
    """P[alpha, beta] from a start decohered in the eigenbasis at l0, through the operator U."""
    V0 = np.linalg.eigh(np.asarray(H(l0), dtype=complex))[1]
    W0 = np.abs(V0)**2
    return W0 @ (np.abs(U @ V0)**2).T


def decohered_rho0(H, l0):
    """The decohered start, diag |V0[a, :]|^2, in the eigenbasis of H(l0) (issue #73)."""
    V0 = np.linalg.eigh(np.asarray(H(l0), dtype=complex))[1]
    W0 = V0.real**2 + V0.imag**2
    return np.array([np.diag(W0[a]) for a in range(W0.shape[0])], dtype=complex)


def brute_force(H, D, l0, l1, sigma, n_nodes=301, start='decohered'):
    """The definition, by quadrature over u of the exact evolution under H + u D_diag, from a
    decohered start or from the flavor state."""
    def Hu(u):
        def f(l):
            Hl = np.asarray(H(l), dtype=complex)
            Dl = np.asarray(D(l), dtype=complex)
            lam, V = np.linalg.eigh(Hl)
            sl = np.real(np.einsum('...ai,...ab,...bi->...i', V.conj(), Dl, V))
            return Hl + u*((V*sl[..., None, :]) @ np.swapaxes(V.conj(), -1, -2))
        return f
    x, w = np.polynomial.hermite_e.hermegauss(n_nodes)
    w = w/w.sum()
    P = 0.0
    for xk, wk in zip(x, w):
        if wk < 1e-16:
            continue
        U, ok = ad._local_evolution_operator(Hu(sigma*xk), l0, l1, 6, 'gl')
        assert ok
        P = P + wk*(decohered_start_probability(H, U, l0) if start == 'decohered'
                    else np.abs(U).T**2)
    return P


def test_a_profile_without_windows_returns_the_decohered_limit():
    """No window: adiabatic transport of a decohered start carries no interference, and the
    result is what averaged_probabilities_adiabatic returns."""
    H, D, l0, _ = crossing_H(g=2.0, l_scale=30.0)     # strongly coupled and slow: no hop
    P, report = ap.phase_averaged_probabilities_adiabatic(H, D, l0, 150.0,
                                                          rho0=decohered_rho0(H, l0))
    P_old, old = ap.averaged_probabilities_adiabatic(H, l0, 150.0)
    assert report['windows'] == [] and old['windows'] == []
    assert report['method'] == 'none'
    assert maxabs(P - P_old) < 1e-15


def test_at_zero_spread_it_is_the_probability_from_a_decohered_start():
    H, D, l0, l1 = crossing_H()
    P, report = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1, spread=0.0, threshold=0.01,
                                                          rho0=decohered_rho0(H, l0))
    assert report['windows']
    U, ok = ad._local_evolution_operator(H, l0, l1, 6, 'gl')
    assert ok
    assert maxabs(P - decohered_start_probability(H, U, l0)) < 1e-5


def test_it_is_the_definition_computed_by_brute_force():
    H, D, l0, l1 = crossing_H()
    P, report = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1, rho0=decohered_rho0(H, l0))
    assert report['sigma_sensitivity'] > 1e-2           # the spread matters here: a real test
    assert maxabs(P - brute_force(H, D, l0, l1, ap.AVG_PHASE_SPREAD)) < 1e-4


def test_from_the_flavor_state_it_is_the_definition_computed_by_brute_force():
    """Issue #73: the default start is the flavor state.  Its coherence is carried from l0 to the
    first window, through it, and on to detection; the brute force starts in flavor too."""
    H, D, l0, l1 = crossing_H()
    P, report = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1)
    assert report['windows']
    assert maxabs(P - brute_force(H, D, l0, l1, ap.AVG_PHASE_SPREAD, start='flavor')) < 1e-4
    # and it is a different answer from the decohered start's, on this path
    assert maxabs(P - brute_force(H, D, l0, l1, ap.AVG_PHASE_SPREAD)) > 1e-3


def test_it_does_not_depend_on_where_the_windows_are_drawn():
    """One window over the whole path, the windows the search finds, and those windows split
    by a stretch the search would have carried adiabatically, all give the same number."""
    H, D, l0, l1 = crossing_H()
    P_found, found = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1)
    P_one, _ = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1, windows=[(l0, l1)])
    (b, c), = found['windows']
    P_wide, _ = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1,
                                                          windows=[(b, c), (c + 2.0, c + 5.0)])
    assert maxabs(P_found - P_one) < 1e-4
    assert maxabs(P_found - P_wide) < 1e-4


def test_a_final_stretch_whose_phase_runs_fast_carries_no_interference():
    """Past about 100 rad at sigma = 10 % a stretch's phase is gone: lengthening the path from
    there moves the phase by tens of radians and the probability only by what the window search
    moves -- its probe grid spans the path, so the window edges shift, by the 1e-7 the patches
    are converged to."""
    H, D, l0, _ = crossing_H()
    P1, _ = ap.phase_averaged_probabilities_adiabatic(H, D, l0, 300.0)
    P2, _ = ap.phase_averaged_probabilities_adiabatic(H, D, l0, 340.0)
    assert maxabs(P1 - P2) < 1e-6


def test_rows_sum_to_one():
    H, D, l0, l1 = crossing_H()
    for spread in (0.0, 0.03, 0.1, 0.3):
        P, _ = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1, spread=spread)
        assert maxabs(P.sum(axis=1) - 1) < 1e-9
        assert np.all(P > -1e-12)


def test_the_transport_carries_the_geometric_phase_of_a_complex_hamiltonian():
    """With complex non-standard interactions the eigenvectors pick up a geometric phase along a
    stretch, and eigh gives each of them an arbitrary phase besides; the stretch's transport,
    expressed between arbitrary end bases, must be the package's adiabatic propagator."""
    hv = np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(S12, S23, S13, DCP, D21, D31),
                    dtype=complex)
    eps = np.array([[0.2, 0.3*np.exp(1.1j), 0.2*np.exp(-0.7j)],
                    [0.0, -0.1, 0.4*np.exp(0.5j)], [0.0, 0.0, 0.05]], dtype=complex)
    eps = np.triu(eps) + np.triu(eps, 1).conj().T
    E, scale = 5.0e6, gd.L_SCALE_SUN

    def H(l):
        v = 3.0e-12*np.exp(-np.asarray(l, dtype=float)/scale)
        return hv/E + np.asarray(v)[..., None, None]*(np.diag([1.0, 0.0, 0.0]) + eps)

    def D(l):
        return -hv/E + 0.0*np.asarray(H(l))

    a, z = 0.5*scale, 2.5*scale
    rng = np.random.default_rng(3)
    Va = np.linalg.eigh(H(a))[1]*np.exp(1j*rng.uniform(0.0, 2*np.pi, 3))
    Vz = np.linalg.eigh(H(z))[1]*np.exp(1j*rng.uniform(0.0, 2*np.pi, 3))
    st = ap._stretch(H, D, a, z, Va, Vz, 0.1)
    assert st['converged']
    T = Vz.conj().T @ ad.adiabatic_propagator(H, a, z, n_points=20001) @ Va
    assert maxabs(T - np.diag(np.diag(T))) < 1e-10
    assert maxabs(np.diag(T) - np.exp(-1j*st['phase'])) < 1e-6
    # and the geometric part is really there: not a multiple of pi on every level
    transport = ap._stretch_once(H, D, a, z, 801)[2]
    assert np.any(np.abs(np.sin(transport)) > 0.01)


def test_the_profile_route_rejects_a_negative_spread():
    H, D, l0, l1 = crossing_H()
    with pytest.raises(ValueError, match='spread'):
        ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1, spread=-1.0)


# ------------------------------------------------------------------ through the wrappers

import warnings  # noqa: E402

import magnus.oscprob as op  # noqa: E402
import magnus.solarmodels as solarmodels  # noqa: E402

OSC = dict(s12=S12, s23=S23, s13=S13, dCP=DCP, D21=D21, D31=D31)


def call(fn, *args, **kwargs):
    """The result and the categories of the warnings the call raised."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = np.asarray(fn(*args, **kwargs))
    return P, {x.category.__name__ for x in w}


def test_where_every_phase_has_decohered_the_old_value_comes_back_bit_for_bit():
    E = np.array([1.0, 10.0, 100.0])*gd.UNIT_TEV
    P, warned = call(op.osc_prob_3nu_vacuum, E, 1.0e20*gd.UNIT_KM, **OSC, average=True)
    H = np.array([vacuum_H(e)[0] for e in E])
    assert np.array_equal(P, ap.averaged_probabilities_constant_hamiltonian(H))
    assert 'PhaseAveragingWarning' not in warned


def test_a_solar_curve_without_windows_is_untouched():
    """Adiabatic transport of a decohered start carries no interference, so the solar MSW curve
    is what the decohered route returns, bit for bit, and it is not even recomputed."""
    model = 'BS05-AGS-OP'
    E = np.logspace(-1.0, np.log10(20.0), 12)*gd.UNIT_MEV
    R = solarmodels.table_edge(model)
    info = {}
    P, warned = call(op.osc_prob_3nu_sun, E, R, 0.0, **OSC, nu_i=0, nu_f=0, average=True,
                     density_profile=model, strategy_info=info)
    assert 'PhaseAveragingWarning' not in warned
    assert info['trace'][-1]['recomputed'] == 0


def test_a_short_baseline_returns_the_phase_average_and_says_it_depends_on_the_spread():
    E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
    P, warned = call(op.osc_prob_3nu_vacuum, E, L, **OSC, average=True)
    H, D = vacuum_H(E)
    expected, sens = ap.phase_averaged_probabilities_constant_hamiltonian(H, D, L)
    assert maxabs(P - expected) < 1e-6            # D here is exact, the wrapper's a difference
    assert sens > ap.PHASE_SPREAD_SENSITIVITY_THRESHOLD
    assert 'PhaseAveragingWarning' in warned
    # and it is not the old answer, which kept every pair here at zero phase: the identity
    assert maxabs(P - np.eye(3)) > 0.05


def test_average_spread_reaches_the_wrappers_through_their_keywords():
    E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
    H, D = vacuum_H(E)
    for spread in (0.0, 0.05, 0.3):
        P, _ = call(op.osc_prob_3nu_vacuum, E, L, **OSC, average=True, average_spread=spread)
        expected, _ = ap.phase_averaged_probabilities_constant_hamiltonian(H, D, L, spread=spread)
        assert maxabs(P - expected) < 1e-6
    P0, _ = call(op.osc_prob_3nu_vacuum, E, L, **OSC, average=True, average_spread=0.0)
    assert maxabs(P0 - np.abs(scipy.linalg.expm(-1j*H*L).T)**2) < 1e-6


@pytest.mark.parametrize("bad", [-0.1, float('nan'), 'wide', True])
def test_average_spread_must_be_a_non_negative_number(bad):
    with pytest.raises(ValueError, match='average_spread'):
        op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM, **OSC, average=True,
                               average_spread=bad)


# ------------------------------------------------------------ the energy-window route (#134)

_STEP_AT = 1500.0*gd.UNIT_KM


def _step_density(l):
    """Two layers with one declared jump: the energy-window route of average=True."""
    out = np.where(np.asarray(l, dtype=float) < _STEP_AT, 3.0, 5.0)
    return out[()] if np.ndim(out) == 0 else out


_STEP_CALL = dict(t_breakpoints=[_STEP_AT], nu_i=0, nu_f=0, density_matter_is_in_g_per_cm3=True)
_STEP_OSC = dict(sth=0.55, Dm2=2.5e-3)
_STEP_E, _STEP_L = 1.0*gd.UNIT_GEV, 3000.0*gd.UNIT_KM


def _step_window_average(**kw):
    return op.osc_prob_matter_std_potential(2, _step_density, _STEP_E, _STEP_L, _STEP_OSC,
                                            average=True, **_STEP_CALL, **kw)


def _step_by_hand(relative_spread, n_samples):
    def prob_of_energy(e):
        return op.osc_prob_matter_std_potential(2, _step_density, e, _STEP_L, _STEP_OSC,
                                                **_STEP_CALL)
    return ap.averaged_probabilities_numerically(prob_of_energy, _STEP_E,
                                                 relative_spread=relative_spread,
                                                 n_samples=n_samples)


def test_average_spread_and_n_samples_set_the_energy_window():
    """Issue #134: across declared discontinuities, average_spread was accepted and ignored, and
    the sample count could not be set.  Both now reach the window: the wrapper's average is
    the by-hand one to round-off.  Bit for bit on some platforms, one ulp apart on others
    (CI's Python 3.11 to 3.13 runners, since the #168 conversion changed the potential)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = _step_window_average(average_spread=0.05, average_n_samples=11)
        P_default = _step_window_average()
        expected, _ = _step_by_hand(0.05, 11)
        expected_default, _ = _step_by_hand(ap.AVG_DEFAULT_ENERGY_SPREAD,
                                            ap.AVG_DEFAULT_N_SAMPLES)
    assert abs(float(P) - float(expected)) < 1e-14
    assert abs(float(P_default) - float(expected_default)) < 1e-14
    assert abs(float(P) - float(P_default)) > 1e-3     # the keywords are not ignored


def test_the_window_warning_and_the_trace_quote_the_window_used():
    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        _step_window_average(average_spread=0.05, average_n_samples=11, strategy_info=info)
    text = ' '.join(str(w.message) for w in caught if w.category is op.PhaseAveragingWarning)
    assert '+/-5.0%' in text and 'over 11 samples' in text
    last = info['trace'][-1]
    assert last['engine'] == 'average'
    assert last['window_half_width'] == 0.05 and last['n_samples'] == 11
    assert last['largest_sem'] > 0.0


@pytest.mark.parametrize("bad", [1, 0, True, 2.0, 'many'])
def test_average_n_samples_must_be_an_integer_of_at_least_two(bad):
    with pytest.raises(ValueError, match='average_n_samples'):
        _step_window_average(average_n_samples=bad)


@pytest.mark.parametrize("bad", [0.0, 1.0, 1.5])
def test_the_window_half_width_must_lie_between_zero_and_one(bad):
    with pytest.raises(ValueError, match='average_spread'):
        _step_window_average(average_spread=bad)


def test_average_n_samples_is_refused_where_nothing_is_sampled():
    """The closed form and the adiabatic route sample no energies: the keyword would be ignored
    there, so it is refused rather than dropped."""
    with pytest.raises(ValueError, match='average_n_samples'):
        op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM, **OSC, average=True,
                               average_n_samples=11)
    with pytest.raises(ValueError, match='average_n_samples'):
        op.osc_prob_matter_std_potential(2, lambda l: 3.0 + 2.0*float(l)/_STEP_L, _STEP_E,
                                         _STEP_L, _STEP_OSC, average=True,
                                         average_n_samples=11, nu_i=0, nu_f=0,
                                         density_matter_is_in_g_per_cm3=True)


def test_a_profile_with_a_window_is_recomputed_and_reported():
    """The end of the solar disk at 10 TeV (issue #62): one window over the whole chord, whose
    readout the decohered route discarded; the phase average keeps it."""
    import magnus.matter  # noqa: F401 -- the chord below is a plain density function
    R = gd.SUN_RADIUS*gd.UNIT_KM
    ne = solarmodels.electron_density_profile('B16-GS98')
    b = 0.9*R
    hl = np.sqrt(R**2 - b**2)

    def ne_chord(l):
        return ne(np.sqrt((np.asarray(l, dtype=float) - hl)**2 + b**2))

    info = {}
    nufit = gd.load_nufit_params('NuFIT 6.1')          # the parameters of the reference below
    kw = dict(nu_i=0, nu_f=0, density_is_of_number_of_electrons=True, L0=0.0,
              average_initial_state='decohered')     # a neutrino crossing the Sun (issue #73)
    P, warned = call(op.osc_prob_matter_std_potential, 3, ne_chord, 1.0e4*gd.UNIT_GEV, 2*hl,
                     nufit, average=True, strategy_info=info, **kw)
    assert info['trace'][-1]['recomputed'] == 1
    # a brute-force average of the definition gives 0.54109; the decohered route gave 0.40432
    assert abs(float(P) - 0.54109) < 1e-4
    assert 'PhaseAveragingWarning' not in warned      # matter phases: the spread does not matter


def solar_chord(br):
    """The electron density along the chord at impact parameter br R_sun, B16-GS98, and its
    length."""
    R = gd.SUN_RADIUS*gd.UNIT_KM
    ne = solarmodels.electron_density_profile('B16-GS98')
    b = br*R
    hl = np.sqrt(R**2 - b**2)

    def ne_chord(l):
        return ne(np.sqrt((np.asarray(l, dtype=float) - hl)**2 + b**2))

    return ne_chord, 2*hl


# Chords through the Sun: neutrinos from a distant source, which arrive decohered (issue #73).
CHORD_KW = dict(nu_i=0, nu_f=0, density_is_of_number_of_electrons=True, L0=0.0,
                average_initial_state='decohered')


def test_rtol_and_atol_are_the_tolerance_of_the_phase_average(monkeypatch):
    """Issue #65: on a smooth profile they were accepted and dropped.  The window patches and
    the stretch phases now converge to the tighter of the two."""
    seen = []
    real = ap.phase_averaged_probabilities_adiabatic

    def spy(*args, **kwargs):
        seen.append((kwargs.get('patch_atol'), kwargs.get('phase_tol')))
        return real(*args, **kwargs)

    monkeypatch.setattr(ap, 'phase_averaged_probabilities_adiabatic', spy)
    ne_chord, L = solar_chord(0.9)
    nufit = gd.load_nufit_params('NuFIT 6.1')
    for kw, want in ((dict(), 1.0e-3), (dict(rtol=1.0e-4), 1.0e-4),
                     (dict(rtol=1.0e-2, atol=1.0e-6), 1.0e-6)):
        call(op.osc_prob_matter_std_potential, 3, ne_chord, 1.0e4*gd.UNIT_GEV, L, nufit,
             average=True, **CHORD_KW, **kw)
        assert seen[-1] == (want, want)


def test_the_decohered_limit_stands_only_within_the_tolerance_asked_for(monkeypatch):
    """A phase average 5e-5 from the limit: inside the gate of 1e-4 at the default tolerance,
    so the limit comes back bit for bit; outside it at atol = 1e-5, so the phase average does."""
    old = {}
    real_old = ap.averaged_probabilities_adiabatic

    def spy_old(*args, **kwargs):
        P, report = real_old(*args, **kwargs)
        old['P'] = P.copy()
        return P, report

    def near_new(*args, **kwargs):
        return old['P'] + 5.0e-5, dict(sigma_sensitivity=0.0, patches_converged=True,
                                       phases_converged=True)

    monkeypatch.setattr(ap, 'averaged_probabilities_adiabatic', spy_old)
    monkeypatch.setattr(ap, 'phase_averaged_probabilities_adiabatic', near_new)
    ne_chord, L = solar_chord(0.9)
    nufit = gd.load_nufit_params('NuFIT 6.1')
    args = (op.osc_prob_matter_std_potential, 3, ne_chord, 1.0e4*gd.UNIT_GEV, L, nufit)
    P, _ = call(*args, average=True, **CHORD_KW)
    assert float(P) == old['P'][0, 0]
    P, _ = call(*args, average=True, atol=1.0e-5, **CHORD_KW)
    assert float(P) == old['P'][0, 0] + 5.0e-5


def test_the_module_tolerances_are_read_at_each_call(monkeypatch):
    """PHASE_AVERAGE_PATCH_ATOL used to be bound as a default argument at import, so setting it
    did nothing."""
    seen = []
    real = ad._local_evolution_operator

    def spy(*args, **kwargs):
        seen.append(kwargs.get('patch_atol'))
        return real(*args, **kwargs)

    monkeypatch.setattr(ad, '_local_evolution_operator', spy)
    monkeypatch.setattr(ap, 'PHASE_AVERAGE_PATCH_ATOL', 3.0e-4)
    H, D, l0, l1 = crossing_H()
    _, report = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1)
    assert report['windows'] and seen and all(x == 3.0e-4 for x in seen)


def test_a_converged_call_does_not_warn_about_the_coarse_levels_of_its_windows():
    """Issue #66: the paper's chord (b = 0.3 R_sun, 100 GeV) printed MagnusConvergenceWarning
    four times, from the first levels of the window ladders, and returned a converged value."""
    import magnus.magnus as mm
    ne_chord, L = solar_chord(0.3)
    nufit = gd.load_nufit_params('NuFIT 6.1')
    args = (op.osc_prob_matter_std_potential, 3, ne_chord, 100.0*gd.UNIT_GEV, L, nufit)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with mm._deferred_slab_norm() as sink:
            call(*args, average=True, **CHORD_KW)
    # The four warnings of #66 came from the trace of Omega, a global phase: measured on its
    # traceless part (issue #155 §1), no level of these windows is wide at all.  Deferral is
    # still covered by test_only_the_phase_average_defers_the_slab_norm_check.
    assert max(sink) < np.pi
    P, warned = call(*args, average=True, **CHORD_KW)
    assert 'MagnusConvergenceWarning' not in warned
    # 0.304987 until 1.2.0: this chord is near the atmospheric resonance, where the 0.4% higher
    # solar electron density of issue #168 moves the average by 2e-2.
    assert abs(float(P) - 0.283618) < 1e-6


def test_a_window_ladder_still_warns_about_the_level_it_returns():
    """Deferred, not dropped: a ladder stopped at a coarse level warns, deferring or not."""
    import magnus.magnus as mm
    H, _, l0, l1 = crossing_H()
    for defer in (False, True):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            # Stopped at two slabs over 22 units of length: wide in the traceless part too.
            _, ok = ad._local_evolution_operator(H, l0, l1, 6, 'gl', n_slabs0=1, max_n_slabs=2,
                                                 defer_slab_norm=defer)
        assert not ok
        assert any(issubclass(x.category, mm.MagnusConvergenceWarning) for x in w)


def test_only_the_phase_average_defers_the_slab_norm_check():
    """A ladder that starts coarse and converges: the hybrid strategy's patch still warns about
    the coarse level, the phase average's does not."""
    import magnus.magnus as mm
    H, _, l0, l1 = crossing_H()
    for defer, expected in ((False, True), (True, False)):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            _, ok = ad._local_evolution_operator(H, l0, l1, 6, 'gl', n_slabs0=2,
                                                 defer_slab_norm=defer)
        assert ok
        assert any(issubclass(x.category, mm.MagnusConvergenceWarning) for x in w) == expected


def test_a_round_off_pseudo_dirac_pair_does_not_warn_through_the_dispatcher():
    """Issue #61 as reported: osc_prob_energy_baseline on a pseudo-Dirac pair at 100 TeV over
    100 Mpc warned that the pair was undecided below 1e-18 eV^2, on a round-off phase."""
    U = hams.pmns_mixing_matrix(S12, S23, S13, DCP)
    m2 = np.array([0.0, D21, D31])
    E, L = 100.0*gd.UNIT_TEV, 100.0*3.0857e19*gd.UNIT_KM
    for dm2 in (1.0e-19, 1.0e-21):
        P, warned = call(op.osc_prob_energy_baseline,
                         lambda e: hams.hamiltonian_pseudo_dirac_vacuum(e, U, m2, {1: dm2}), E, L,
                         0.0, H_func_is_function_only_of_energy=True, average=True)
        assert 'PhaseAveragingWarning' not in warned
        assert maxabs(P.sum(axis=-1) - 1) < 1e-12


def test_a_hamiltonian_without_energy_dependence_keeps_the_decohered_limit():
    """A fixed matrix has no slope for an energy spread to act on: average=True returns the
    L/E -> infinity limit for it, as before, and warns as before where that does not apply."""
    H, _ = vacuum_H(1.0*gd.UNIT_GEV)
    P, _ = call(op.osc_prob_energy_baseline, H, 1.0*gd.UNIT_GEV, 1.0e8*gd.UNIT_KM, 0.0,
                average=True)
    assert np.array_equal(P, ap.averaged_probabilities_constant_hamiltonian(H, 1.0e8*gd.UNIT_KM))


# ----------------------------------------------------------------------------------------------
# The initial state of average=True (issue #73)
# ----------------------------------------------------------------------------------------------

# The values moved with the electron density in 1.2.0, from 0.36263161 and 0.45279506 (issue #168).
@pytest.mark.parametrize('start, expected', [('flavor', 0.36263347), ('decohered', 0.45316241)])
def test_the_initial_state_means_the_same_on_every_route(start, expected):
    """The reproduction of issue #73: one medium through the constant route and, as an
    exponential profile of scale height 1e9 km, through the smooth route.  They used to give
    0.3626 (a flavor start) and 0.4528 (a decohered start); each start now gives one answer on
    both routes, the hand-computed closed forms of the issue."""
    nufit = gd.load_nufit_params('NuFIT 6.1')
    kw = dict(nu_i=gd.NUMU, nu_f=gd.NUMU, density_matter_is_in_g_per_cm3=True, average=True,
              average_initial_state=start, **nufit)
    E, L = np.array([1.0*gd.UNIT_GEV]), 300.0*gd.UNIT_KM
    P_const = float(np.ravel(op.osc_prob_3nu_matter_constant_density(E, L=L, rho=3.0, **kw))[0])
    P_flat = float(np.ravel(op.osc_prob_3nu_matter_exp_density(E, L=L, L0=0.0, rho_central=3.0,
                                                                l_scale=1e9*gd.UNIT_KM, **kw))[0])
    assert abs(P_const - expected) < 1e-7
    assert abs(P_flat - expected) < 1e-7


def test_a_decohered_start_on_a_constant_hamiltonian_is_the_limit_exactly():
    """Whatever the phases, even where none has averaged away (a short vacuum baseline), an
    incoherent mixture of eigenstates never interferes: the limit, sum_i |U_ai|^2 |U_bi|^2."""
    P = op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 10.0*gd.UNIT_KM, **OSC, average=True,
                               average_initial_state='decohered')
    U = np.asarray(hams.pmns_mixing_matrix(OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP']))
    W = np.abs(U)**2
    assert maxabs(np.asarray(P) - W @ W.T) < 1e-14


def test_for_neutrinos_made_in_the_solar_core_both_starts_agree_bit_for_bit():
    """Produced in the core, every interference term is damped from production on, so the
    default start in flavor changes nothing on the solar curves: the value is the decohered
    start's, bit for bit, and no transport is paid for."""
    R = gd.SUN_RADIUS*gd.UNIT_KM
    E = np.geomspace(0.1e6, 20e6, 40)
    kw = dict(nu_i=0, nu_f=0, density_profile='BS05-AGS-OP', average=True)
    P_flavor = op.osc_prob_3nu_sun(E, R, 0.0, **kw)
    P_decoh = op.osc_prob_3nu_sun(E, R, 0.0, average_initial_state='decohered', **kw)
    assert np.array_equal(np.asarray(P_flavor), np.asarray(P_decoh))


def test_the_flavor_start_on_a_path_with_no_crossing_is_the_definition():
    """No window: the closed form of the paper's Eq. (phase_average_nocross), against the brute
    force started in flavor.  The bound is the adiabatic transport's own error on this path,
    which no tolerance moves: 4.9e-6 from a flavor start, 8.2e-7 from a decohered one, since a
    coherent start feels the weak non-adiabaticity no window was opened for at first order."""
    H, D, l0, _ = crossing_H(g=2.0, l_scale=30.0)     # strongly coupled and slow: no hop
    l1 = 8.0
    P, report = ap.phase_averaged_probabilities_adiabatic(H, D, l0, l1)
    assert report['windows'] == [] and report['method'] == 'closed_form'
    assert maxabs(P - brute_force(H, D, l0, l1, ap.AVG_PHASE_SPREAD, start='flavor')) < 1e-5


def test_a_decohered_start_is_refused_on_declared_discontinuities():
    """That route averages over an energy window from the flavor state (paper Route 2); there is
    no eigenbasis to decohere in at a jump."""
    with pytest.raises(ValueError, match='average_initial_state'):
        op.osc_prob_3nu_earth(3.0*gd.UNIT_GEV, loc_ini='fermilab', loc_fin='homestake', nu_i=1,
                              nu_f=0, average=True, average_initial_state='decohered')


@pytest.mark.parametrize('bad', ['coherent', 'Flavor', 1, True])
def test_an_unknown_initial_state_raises(bad):
    with pytest.raises(ValueError, match='average_initial_state'):
        op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM, **OSC, average=True,
                               average_initial_state=bad)


def test_the_initial_state_is_refused_without_average():
    """It was ignored without average=True, as documented; by the author's decision on
    issue #160 it is refused by name, valid value or not."""
    for value in ('decohered', 'bogus'):
        with pytest.raises(ValueError, match='average_initial_state'):
            op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM, **OSC,
                                   average_initial_state=value)
