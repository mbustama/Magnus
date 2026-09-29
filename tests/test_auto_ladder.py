# -*- coding: utf-8 -*-
"""Tests of the ladder route of ``strategy='auto'`` (issue #70).

On a smooth profile, ``'auto'`` used to run the hybrid strategy at every tolerance.  The hybrid's
cost is its window search, which does not follow the tolerance, so at the default of 1e-3 it was
the slower route by one to two orders of magnitude: the four scans of the paper's Fig. 1 took
8 s where the ladder takes 40 ms.  ``'auto'`` now hands a request to the ladder when its
estimated phase is at most ``AUTO_LADDER_MAX_PHASE`` and its tolerance no tighter than
``AUTO_LADDER_MIN_TOLERANCE``, provided its starting slab count leaves room below the cap; the
ladder then runs at a tenth of the tolerance, above a slab floor, and without the
interaction-picture fast path.

Below ``AUTO_LADDER_MIN_TOLERANCE`` the route used to close.  It now stays open on
``integration_method='gl'`` for a phase within a limit that shrinks with the tolerance and the
order, capped at ``AUTO_LADDER_TIGHT_MAX_PHASE``, and the ladder runs at the tolerance itself
(issue #120).  That is what lets the paper's Listing 1 drop ``strategy='magnus'``.

The reference throughout is the hybrid strategy at a tolerance of 1e-12, which the paper checks
against independent codes.
"""

import warnings

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import magnus.adiabatic as ad
import magnus.globaldefs as gd
import magnus.magnus as mg
import magnus.oscprob as op

OSC = dict(s12=0.55, s23=0.68, s13=0.15, dCP=3.7, D21=7.5e-5, D31=2.5e-3)
# The exponential profile of the paper's Fig. 1: 3e3 g/cm^3 at the origin, 10 km scale height,
# 25 km.
PROFILE = dict(L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.e3, l_scale=10.0*gd.UNIT_KM,
               density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE)
ENERGIES = np.logspace(np.log10(0.002), np.log10(0.2), 12)*gd.UNIT_GEV


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x))))


def p3(energy, **kw):
    return np.asarray(op.osc_prob_3nu_matter_exp_density(energy, **PROFILE, **OSC, **kw))


def p2(energy, **kw):
    return np.asarray(op.osc_prob_2nu_matter_exp_density(energy, **PROFILE, sth=0.55, Dm2=7.5e-5,
                                                         **kw))


def declined(info):
    return dict(info.get('declined', []))


@pytest.fixture(scope='module')
def reference():
    return p3(ENERGIES, rtol=1e-12, atol=1e-14, strategy='hybrid')


def test_default_tolerance_takes_the_ladder_and_meets_it(reference):
    """The case the issue is about: within the tolerance asked for, with no convergence warning.
    The warning fired on every such request before the slab floor, about rungs the ladder went
    on to refine away."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('error', mg.MagnusConvergenceWarning)
        P = p3(ENERGIES, rtol=1e-3, atol=1e-3, strategy_info=info)
    assert declined(info).get('hybrid') == 'auto prefers the ladder'
    assert info['engine'] != 'hybrid'
    assert maxabs(P - reference) < 1e-3


def test_the_decision_is_reported():
    info = {}
    p3(ENERGIES[0], rtol=1e-3, atol=1e-3, strategy_info=info)
    note = [t for t in info['trace'] if t.get('reason') == 'auto prefers the ladder'][0]
    assert 0.0 < note['estimated_phase'] <= op.AUTO_LADDER_MAX_PHASE
    assert note['min_n_slabs'] >= 1
    assert note['tolerance_margin'] == op.AUTO_LADDER_TOLERANCE_MARGIN


@pytest.mark.parametrize('tol', [1e-9, 1e-12])
def test_a_tight_tolerance_takes_the_ladder_on_a_small_phase(tol):
    """Replaces the test that a tolerance below AUTO_LADDER_MIN_TOLERANCE keeps the hybrid.  The
    ladder had cost more there only because it ran at a tenth of the tolerance; at the tolerance
    itself it is the faster route on this profile (78 rad) at every tolerance measured for issue
    #120."""
    info = {}
    p3(ENERGIES[:2], rtol=tol, atol=tol, strategy_info=info)
    assert info['engine'] == 'separable'
    note = [t for t in info['trace'] if t.get('reason') == 'auto prefers the ladder'][0]
    assert note['tolerance_margin'] == 1.0
    assert note['estimated_phase'] <= note['phase_limit'] == op._auto_ladder_max_phase(tol, 4, True)


def _listing1_curves():
    """The four calls of the paper's Listing 1, four energies per curve instead of 140."""
    osc = gd.load_nufit_params('NuFIT 6.1')
    s14 = s24 = np.sqrt(0.10)
    s15 = s25 = np.sqrt(0.06)
    return (('2nu', op.osc_prob_2nu_matter_exp_density, 0.0005, 0.05,
             dict(sth=osc['s12'], Dm2=osc['D21'])),
            ('3nu', op.osc_prob_3nu_matter_exp_density, 0.002, 0.2, dict(osc)),
            ('4nu', op.osc_prob_4nu_matter_exp_density, 2.0, 20.0,
             dict(osc, s14=s14, s24=s24, D41=1.0)),
            ('5nu', op.osc_prob_5nu_matter_exp_density, 2.0, 20.0,
             dict(osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7)))


@pytest.mark.parametrize('order', [8, None])
def test_listing_1_without_a_strategy_meets_its_tolerance(monkeypatch, order):
    """The paper's Listing 1 without ``strategy='magnus'`` (issue #120), at rtol = 1e-12 and
    atol = 1e-14.  At order 8, as printed, every curve goes to the energy-batched ladder.  With
    the order left at its default, whatever the rule decides, the worst point (3nu, 2 MeV) is
    checked against DOP853 at 1e-13 on the Hamiltonian the wrapper builds: 2.3e-14 off at order
    8 and 3.1e-14 at order 4, against a tolerance of 7.1e-13 there."""
    osc_kw = dict(rtol=1e-12, atol=1e-14)
    if order is not None:
        osc_kw['magnus_exp_order'] = order
    for name, fn, lo, hi, ex in _listing1_curves():
        info = {}
        P = np.asarray(fn(np.logspace(np.log10(lo), np.log10(hi), 4)*gd.UNIT_GEV, **PROFILE,
                          **ex, **osc_kw, strategy_info=info)).ravel()
        if order == 8:
            assert info['engine'] == 'separable', name
            assert declined(info).get('hybrid') == 'auto prefers the ladder', name
        if name == '3nu':
            P_worst = P[0]
            osc3 = ex

    captured = {}
    real = ad.hybrid_propagator

    def hp(H_func, *args, **kwargs):
        captured['H'] = H_func
        return real(H_func, *args, **kwargs)

    monkeypatch.setattr(ad, 'hybrid_propagator', hp)
    op.osc_prob_3nu_matter_exp_density(np.array([0.002])*gd.UNIT_GEV, **PROFILE, **osc3,
                                       rtol=1e-3, atol=1e-3, strategy='hybrid')
    H_func = captured['H']
    sol = solve_ivp(lambda l, y: (-1j*np.asarray(H_func(l)) @ y.reshape(3, 3)).ravel(),
                    (0.0, 25.0*gd.UNIT_KM), np.eye(3, dtype=complex).ravel(), rtol=1e-13,
                    atol=1e-15, method='DOP853')
    P_ref = abs(sol.y[0, -1])**2
    assert abs(P_worst - P_ref) <= 1e-14 + 1e-12*P_ref


@pytest.mark.parametrize('strategy', ['hybrid', 'magnus'])
def test_an_explicit_strategy_is_not_rerouted_at_a_tight_tolerance(strategy):
    """The tightened route is 'auto''s alone as well."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', mg.MagnusConvergenceWarning)
        p3(ENERGIES[:2], rtol=1e-12, atol=1e-12, strategy=strategy, strategy_info=info)
    assert (info['engine'] == 'hybrid') == (strategy == 'hybrid')
    assert 'auto prefers the ladder' not in declined(info).values()


def test_a_tight_tolerance_runs_the_ladder_at_the_tolerance_itself(monkeypatch):
    """At a tenth of the tolerance the ladder was the slower route below
    AUTO_LADDER_MIN_TOLERANCE; at the tolerance itself its deep rungs already overestimate the
    error, by about 1.5**p - 1."""
    seen = {}
    real = op._osc_prob_scan_separable_dispatch

    def spy(*args):
        seen.update(args[-1])
        return real(*args)

    monkeypatch.setattr(op, '_osc_prob_scan_separable_dispatch', spy)
    p3(ENERGIES, rtol=1e-9, atol=2e-9)
    assert seen['rtol'] == 1e-9 and seen['atol'] == 2e-9
    assert op._PreferLadder(40, 1.0).request(1e-12, 1e-14, None, None, 'gl') == (1e-12, 1e-14, 40)


def test_the_tight_phase_limit():
    """AUTO_LADDER_MAX_PHASE at AUTO_LADDER_MIN_TOLERANCE, shrinking as tol**(1/p) below it and
    capped at AUTO_LADDER_TIGHT_MAX_PHASE: 1 000 rad for the paper's Listing 1 at order 8, 100 at
    the default order 4.  Only a scan at a loose tolerance goes without a limit (issue #84)."""
    f = op._auto_ladder_max_phase
    assert f(1e-3, 4, False) == op.AUTO_LADDER_MAX_PHASE
    assert f(1e-3, 4, True) == np.inf
    assert f(1e-14, 8, True) == pytest.approx(1000.0)
    assert f(1e-14, 4, False) == pytest.approx(100.0)
    assert f(1e-9, 8, False) == op.AUTO_LADDER_TIGHT_MAX_PHASE


def test_a_tight_energy_scan_obeys_the_phase_limit(monkeypatch):
    """At a loose tolerance an energy scan goes to the ladder whatever its phase (issue #84; see
    test_an_energy_scan_ignores_the_phase_threshold).  That was measured at loose tolerances
    only, so below AUTO_LADDER_MIN_TOLERANCE the tightened limit applies to a scan as well."""
    monkeypatch.setattr(op, 'AUTO_LADDER_TIGHT_MAX_PHASE', 1.0)
    info = {}
    p3(ENERGIES[:2], rtol=1e-9, atol=1e-9, strategy_info=info)
    assert info['engine'] == 'hybrid'


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
def test_a_tight_tolerance_off_gauss_legendre_is_not_rerouted(method):
    """Issue #120 measured the ladder on 'gl' only, so below AUTO_LADDER_MIN_TOLERANCE the other
    quadratures keep what 'auto' did before: the hybrid is tried first."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', mg.MagnusConvergenceWarning)
        p3(ENERGIES[:2], rtol=1e-9, atol=1e-9, integration_method=method, strategy_info=info)
    assert 'auto prefers the ladder' not in declined(info).values()


def test_the_full_sun_is_not_rerouted_at_a_tight_tolerance():
    """The slab-count condition keeps it off the ladder at any tolerance: the core density sets a
    starting count near the cap.  Two flavors at 10 MeV over 0.9 R_sun, the case of
    AUTO_LADDER_MAX_FLOOR_FRACTION.  At 1e-9 the hybrid is tried first, as before."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        op.osc_prob_2nu_sun(np.array([10.0e-3])*gd.UNIT_GEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
                            sth=0.55, Dm2=7.5e-5, nu_i=gd.NUE, nu_f=gd.NUE, rtol=1e-9, atol=1e-9,
                            strategy_info=info)
    assert 'hybrid' in declined(info) or info['engine'] == 'hybrid'
    assert 'auto prefers the ladder' not in declined(info).values()


@pytest.mark.parametrize('order', [0, 10])
def test_an_invalid_order_at_a_tight_tolerance_raises_as_before(order):
    """The tightened limit divides by the order, so it is taken only for the orders 'gl'
    supports; any other keeps the hybrid strategy's path and the error it raises there."""
    with pytest.raises(ValueError, match='order'):
        p3(ENERGIES[:2], rtol=1e-9, atol=1e-11, magnus_exp_order=order)


def test_a_phase_above_the_threshold_keeps_the_hybrid(monkeypatch):
    """One energy: an energy scan at a shared baseline is exempt (next test)."""
    info = {}
    monkeypatch.setattr(op, 'AUTO_LADDER_MAX_PHASE', 1.0)
    p3(ENERGIES[:1], rtol=1e-3, atol=1e-3, strategy_info=info)
    assert info['engine'] == 'hybrid'


def test_an_energy_scan_ignores_the_phase_threshold(monkeypatch):
    """The phase limit prices the ladder point by point; the energy-batched engine shares its
    slabs across the energies, so an energy scan at one baseline goes to it whatever its
    phase.  Only the slab-count condition, which keeps the full Sun on the hybrid, applies."""
    info = {}
    monkeypatch.setattr(op, 'AUTO_LADDER_MAX_PHASE', 1.0)
    p3(ENERGIES[:2], rtol=1e-3, atol=1e-3, strategy_info=info)
    assert info['engine'] == 'separable'
    assert declined(info).get('hybrid') == 'auto prefers the ladder'


@pytest.mark.parametrize('strategy', ['hybrid', 'magnus'])
def test_an_explicit_strategy_is_not_rerouted(strategy):
    """The preference is 'auto''s alone: 'hybrid' still runs the hybrid, and 'magnus' runs at
    the tolerance it was given rather than a tenth of it."""
    info = {}
    with warnings.catch_warnings():
        # 'magnus' keeps the loose seed of suggest_n_slabs, whose first rungs warn by design.
        warnings.simplefilter('ignore', mg.MagnusConvergenceWarning)
        P = p3(ENERGIES[:2], rtol=1e-3, atol=1e-3, strategy=strategy, strategy_info=info)
    assert (info['engine'] == 'hybrid') == (strategy == 'hybrid')
    assert 'auto prefers the ladder' not in declined(info).values()
    if strategy == 'magnus':
        Q = p3(ENERGIES[:2], rtol=1e-3, atol=1e-3)
        assert not np.array_equal(P, Q)


def test_the_ladder_runs_at_a_tenth_of_the_tolerance_above_the_floor(monkeypatch):
    seen = {}
    real = op._osc_prob_scan_separable_dispatch

    def spy(*args):
        seen.update(args[-1])
        return real(*args)

    monkeypatch.setattr(op, '_osc_prob_scan_separable_dispatch', spy)
    info = {}
    p3(ENERGIES, rtol=1e-3, atol=2e-3, strategy_info=info)
    floor = [t for t in info['trace'] if t.get('reason') == 'auto prefers the ladder'][0]
    assert seen['rtol'] == pytest.approx(1e-4) and seen['atol'] == pytest.approx(2e-4)
    assert seen['min_n_slabs'] == floor['min_n_slabs']


def test_a_larger_caller_floor_wins_and_the_cap_holds():
    prefer = op._PreferLadder(40)
    assert prefer.request(1e-3, None, 100, None, 'gl') == (1e-4, None, 100)
    assert prefer.request(1e-3, 1e-3, 1, 25, 'gl')[2] == 25
    assert prefer.request(1e-3, 1e-3, None, None, 'gl')[2] == 40


def test_two_flavors_skip_the_interaction_picture():
    """The two-flavor exponential profile is where the interaction-picture engine answers, and
    at a loose tolerance it was the slow route the ladder is taken to avoid."""
    info = {}
    P = p2(ENERGIES, rtol=1e-3, atol=1e-3, strategy_info=info)
    assert info['engine'] not in ('hybrid', 'ip_exp')
    ref = p2(ENERGIES, rtol=1e-12, atol=1e-14, strategy='hybrid')
    assert maxabs(P - ref) < 1e-3


def test_the_floor_keeps_every_slab_convergent():
    """The floor is the sufficient condition itself: the largest spectral radius of H met, times
    the span, over pi.  On the Fig. 1 profile no slab the ladder builds reaches pi."""
    norms = []
    with mg._deferred_slab_norm() as sink:
        p3(ENERGIES, rtol=1e-3, atol=1e-3)
        norms = list(sink)
    assert norms and max(norms) < np.pi


def test_a_floor_near_the_cap_keeps_the_hybrid():
    """The Sun is the case: its core density sets the floor for the whole path, and a ladder
    that starts at its cap cannot refine.  Here the cap is lowered instead of the Sun used."""
    info = {}
    p3(ENERGIES[:2], rtol=1e-3, atol=1e-3, max_n_slabs=60, strategy_info=info)
    assert info['engine'] == 'hybrid'


def test_an_undeclared_step_still_warns():
    """The hybrid's resolution test is what warns about a density jump nobody declared.  The
    ladder route skips the hybrid, so it runs the test itself."""
    l1 = gd.L_SCALE_SUN

    def ne(l):
        x = np.asarray(l, dtype=float)
        out = np.where(x < 0.5*l1, 0.02, 0.30)*gd.NUM_DENSITY_E_SUN_CENTRAL
        return out[()] if out.ndim == 0 else out

    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        op.osc_prob_matter_std_potential(2, ne, 50.0e6, l1, {'sth': 0.55, 'Dm2': 7.5e-5},
                                         L0=0.0, density_is_of_number_of_electrons=True,
                                         strategy_info=info)
    assert any(issubclass(w.category, op.UnmarkedDiscontinuityWarning) for w in caught)
    assert 'not resolved' in declined(info)['hybrid']
    assert info['engine'] != 'hybrid'


def test_the_phase_is_the_spread_of_the_eigenvalues():
    """Not the norm of the integrated Hamiltonian, which lets opposite signs cancel: an MSW
    region does that to the matter and vacuum terms, and on the Sun it read 2.2 to 2.9 times
    low.  Here it cancels completely, and the spread does not."""
    a, L = 2.0e-3, 1.0e4
    sz = np.diag([1.0, -1.0]).astype(complex)

    def H_at_energy(enu):
        return lambda l: a*np.cos(np.pi*l/L)*sz

    phase, n_floor = op._estimated_phase(H_at_energy, np.array([1.0]), np.array([L]), 0.0)
    assert phase == pytest.approx(4.0*a*L/np.pi, rel=1e-2)
    assert n_floor == int(np.ceil(L*a/np.pi))
