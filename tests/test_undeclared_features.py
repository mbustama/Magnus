# -*- coding: utf-8 -*-
"""Undeclared jumps and narrow features on the refinement ladder (issue #161).

Two failures, both at the default tolerance, both silent or nearly so on 1.1:

* A slab straddling a density jump nobody declared converges as 1/n, so the per-point ladder
  thrashes and stops wherever two levels happen to agree: the castle wall of the issue was off
  by 0.043 at one baseline, and 24 of 40 baselines were off by more than 1e-3.
  ``osc_prob_energy_baseline`` now notices it in the ladder's own samples, locates the jumps,
  declares them, and warns ``UnmarkedDiscontinuityWarning``.
* A spike narrower than the ladder's starting slabs is stepped over by every level, which then
  agree on the same wrong answer: the #154 spike was off by 5.4e-3 with no warning at all.
  ``strategy='auto'`` now starts the ladder at the spacing of the probe grid that flags it.

References are DOP853 at rtol=1e-12, split at the features.
"""

import warnings

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import magnus.adiabatic as ad
import magnus.globaldefs as gd
import magnus.hamiltonians as hm
import magnus.oscprob as op

KM = gd.UNIT_KM
GEV = gd.UNIT_GEV
E00 = np.diag([1.0, 0.0, 0.0])


def _vcc(rho):
    """sqrt(2) G_F rho N_A Y_e at Y_e = 1/2, the package's conversion since 1.2.0 (issue #168)."""
    return (np.sqrt(2.0)*gd.GF*rho*gd.UNIT_G_PER_CM3/gd.ATOMIC_MASS_UNIT*0.5)


def _dop853(H, L, splits):
    psi = np.eye(3, dtype=complex)
    pts = [0.0] + sorted(s for s in splits if 0.0 < s < L) + [L]
    for a, b in zip(pts[:-1], pts[1:]):
        sol = solve_ivp(lambda l, y: (-1j*H(l) @ y.reshape(3, 3)).ravel(), (a, b), psi.ravel(),
                        method='DOP853', rtol=1e-12, atol=1e-14)
        psi = sol.y[:, -1].reshape(3, 3)
    return (np.abs(psi)**2).T


# The castle wall of the issue: three layers, the Hamiltonian written by hand, no t_breakpoints.
HV_ISSUE = np.asarray(hm.hamiltonian_3nu_vacuum_energy_independent(
    0.55, 0.72, 0.15, 1.3, 7.5e-5, 2.5e-3))
E_CASTLE = 3.0e9
WALLS = [3000.0*KM, 9000.0*KM]


def _castle(e, l):
    return HV_ISSUE/e + _vcc(3.3 if (l < 3000.0*KM or l > 9000.0*KM) else 11.5)*E00


def _castle_ref(L):
    return _dop853(lambda l: _castle(E_CASTLE, l), L, WALLS)


def test_the_castle_wall_of_the_issue_is_found_declared_and_right():
    L = 6202.5641*KM
    with pytest.warns(op.UnmarkedDiscontinuityWarning, match='jumps at l = 1.52'):
        P = op.osc_prob_energy_baseline(_castle, E_CASTLE, L)
    assert np.max(np.abs(P - _castle_ref(L))) < 1e-10          # 0.043 before


def test_at_a_tight_tolerance_too():
    L = 6202.5641*KM
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = op.osc_prob_energy_baseline(_castle, E_CASTLE, L, rtol=1e-8, atol=1e-8)
    assert np.max(np.abs(P - _castle_ref(L))) < 1e-8           # 0.043 before


def test_a_loop_over_baselines_agrees_with_the_reference_at_every_one():
    """The tester's scan, one point at a time: 24 of 40 were off by more than 1e-3."""
    worst = 0.0
    for L in np.linspace(1000.0, 12000.0, 12)*KM:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            P = op.osc_prob_energy_baseline(_castle, E_CASTLE, L)
        worst = max(worst, float(np.max(np.abs(P - _castle_ref(L)))))
    assert worst < 1e-6


def test_a_scan_declares_the_jumps_for_every_point_and_recomputes_the_first():
    Es = np.array([1.0, 3.0, 10.0])*GEV
    L = 6202.5641*KM
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = op.osc_prob_energy_baseline(_castle, Es, L, cumulative=False)
    for i, e in enumerate(Es):
        ref = _dop853(lambda l, e=e: _castle(e, l), L, WALLS)
        assert np.max(np.abs(P[i] - ref)) < 1e-9


def test_declared_breakpoints_and_smooth_profiles_are_never_searched(monkeypatch):
    calls = []
    real = ad._located_jumps
    monkeypatch.setattr(ad, '_located_jumps', lambda *a, **k: calls.append(1) or real(*a, **k))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        op.osc_prob_energy_baseline(_castle, E_CASTLE, 6202.5641*KM, t_breakpoints=WALLS[:1])

        def smooth(e, l):
            return HV_ISSUE/e + _vcc(4.0 + 2.0*np.sin(l/(700.0*KM)))*E00
        op.osc_prob_energy_baseline(smooth, E_CASTLE, 6200.0*KM)
        op.osc_prob_energy_baseline(smooth, np.array([1.0, 3.0])*GEV, 6200.0*KM)
    assert calls == []


def test_the_cross_check_sees_only_the_breakpoints_the_caller_declared():
    """Notebook 22: the jumps found were handed to the expm reference as declared ones, which
    then ran on an undeclared step and answered in the wrong shape."""
    def step(l):
        return np.where(np.asarray(l) < 2000.0*KM, 3.0, 8.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = op.cross_check_strategies(op.osc_prob_matter_std_potential, 3, step, 3.0*GEV,
            5000.0*KM, {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k]
                        for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')},
            density_matter_is_in_g_per_cm3=True)
    assert 'expm' not in out['ran']
    for engine in out['ran']:
        assert np.shape(out['answers'][engine]) == (3, 3)


def test_located_jumps_are_where_the_hamiltonian_jumps_and_nowhere_else():
    def H(l):
        return _castle(E_CASTLE, l)
    jumps = ad._located_jumps(H, 0.0, 11000.0*KM, 200)
    assert np.allclose(jumps, WALLS, rtol=1e-12)

    def steep_but_smooth(l):
        return HV_ISSUE/E_CASTLE + _vcc(4.0 + 3.0*np.tanh((l - 3000.0*KM)/(3.0*KM)))*E00
    assert ad._located_jumps(steep_but_smooth, 0.0, 6000.0*KM, 200) == []


# The spike of #154: smooth, 5 km wide, stepped over by the ladder's starting slabs.
OSC = {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k]
       for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
HV_STD = np.asarray(hm.hamiltonian_3nu_vacuum_energy_independent(**OSC))


def _spike(width_km):
    def rho(x):
        return 3.0 + 50.0*np.exp(-((np.asarray(x, dtype=float) - 1000.0*KM)/(width_km*KM))**2)
    return rho


@pytest.mark.parametrize('energy_gev, width_km', [(1.0, 5.0), (0.5, 5.0), (3.0, 2.0)])
def test_a_narrow_spike_is_resolved_under_auto(energy_gev, width_km):
    """#154 §6: 5.4e-3 off at 1 GeV, width 5 km, with no warning."""
    rho = _spike(width_km)
    E = energy_gev*GEV
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = op.osc_prob_matter_std_potential(3, rho, E, 4000.0*KM, OSC,
                                             density_matter_is_in_g_per_cm3=True)
    ref = _dop853(lambda l: HV_STD/E + _vcc(float(rho(l)))*E00, 4000.0*KM,
                  [(1000.0 - 8*width_km)*KM, (1000.0 + 8*width_km)*KM])
    assert np.max(np.abs(P - ref)) < 1e-3


def test_the_probe_reports_sharp_features_separately_from_jumps():
    def H_of(rho):
        return lambda l: HV_STD/GEV + _vcc(np.asarray(rho(l)))[..., None, None]*E00
    smooth = H_of(lambda x: 3.0 + 2.0*np.sin(np.asarray(x)/(800.0*KM)))
    spike = H_of(_spike(5.0))
    assert ad._profile_resolution(smooth, 0.0, 4000.0*KM, 200) == (True, False)
    assert ad._profile_resolution(spike, 0.0, 4000.0*KM, 200)[1]
    for H in (smooth, spike):
        assert (ad._profile_resolution(H, 0.0, 4000.0*KM, 200)[0]
                == ad._profile_is_resolved(H, 0.0, 4000.0*KM, 200))


def test_a_smooth_profile_keeps_its_starting_grid_under_auto():
    """The floor is raised only where the probe flags a sharp feature: a smooth profile's
    request is handed to the ladder exactly as before."""
    def H_at_energy(e):
        return lambda l: HV_STD/e + _vcc(3.0 + 2.0*np.sin(np.asarray(l)/(800.0*KM)))[..., None, None]*E00
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        prefer = op._auto_prefers_ladder(H_at_energy, np.array([3.0*GEV]), np.array([4000.0*KM]),
                                         0.0, 1e-3, 1e-3, 20000)
    assert prefer is not None
    assert prefer.min_n_slabs < op.AUTO_SHARP_FEATURE_SLABS_PER_PROBE*(op.AUTO_PROBE_POINTS - 1)
