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

The reference throughout is the hybrid strategy at a tolerance of 1e-12, which the paper checks
against independent codes.
"""

import warnings

import numpy as np
import pytest

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
def test_a_tight_tolerance_keeps_the_hybrid(tol):
    """Below AUTO_LADDER_MIN_TOLERANCE the hybrid reaches the tolerance at no extra cost, and
    the ladder's cost grows with it."""
    info = {}
    p3(ENERGIES[:2], rtol=tol, atol=tol, strategy_info=info)
    assert info['engine'] == 'hybrid'
    assert 'hybrid' not in declined(info)


def test_a_phase_above_the_threshold_keeps_the_hybrid(monkeypatch):
    info = {}
    monkeypatch.setattr(op, 'AUTO_LADDER_MAX_PHASE', 1.0)
    p3(ENERGIES[:2], rtol=1e-3, atol=1e-3, strategy_info=info)
    assert info['engine'] == 'hybrid'


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
