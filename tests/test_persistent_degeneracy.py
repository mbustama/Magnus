# -*- coding: utf-8 -*-
"""average=True with sterile flavors nothing couples to (issue #148).

At the default sterile parameters of the 5-flavor Sun wrappers the two sterile states are
decoupled and degenerate along the whole path.  The adiabatic averaging route read the pair as a
crossing spanning the Sun, laid 338 193 energy nodes across it, and never returned.  The
decoupled flavors are now averaged out: they are eigenstates at every position and never
oscillate, and the rest is the 3-flavor problem.
"""

import signal
import warnings

import numpy as np
import pytest

import magnus.avgprob as ap
import magnus.globaldefs as gd
import magnus.oscprob as op

R_SUN = gd.SUN_RADIUS*gd.UNIT_KM
MEV = gd.UNIT_MEV


class _TookTooLong(Exception):
    pass


def _within(seconds, f):
    """Fails, rather than hangs, if the call does not return."""
    def alarm(signum, frame):
        raise _TookTooLong()
    old = signal.signal(signal.SIGALRM, alarm)
    signal.alarm(seconds)
    try:
        return f()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


@pytest.mark.parametrize('nubar', [False, True])
@pytest.mark.parametrize('wrapper, kw', [(op.osc_prob_5nu_sun, {}),
                                         (op.osc_prob_5nu_sun_nsi, {'eps_ee': 0.1}),
                                         (op.osc_prob_5nu_sun, {'density_profile': 'BP04'})])
def test_the_5nu_sun_at_its_defaults_returns_the_3nu_answer(wrapper, kw, nubar):
    E = np.array([1.0, 5.0, 20.0])*MEV
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P5 = np.asarray(_within(60, lambda: wrapper(E, R_SUN, 0.0, average=True, nubar=nubar,
                                                     **kw)))
        three = op.osc_prob_3nu_sun_nsi if wrapper is op.osc_prob_5nu_sun_nsi else op.osc_prob_3nu_sun
        P3 = np.asarray(three(E, R_SUN, 0.0, average=True, nubar=nubar, **kw))
    assert np.max(np.abs(P5[:, :3, :3] - P3)) < 1e-6
    assert np.array_equal(P5[:, 3, 3], np.ones(3)) and np.array_equal(P5[:, 4, 4], np.ones(3))
    assert not P5[:, :3, 3:].any() and not P5[:, 3:, :3].any()
    assert not P5[:, 3, 4].any() and not P5[:, 4, 3].any()


def test_a_single_channel_and_a_scalar_energy():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = _within(60, lambda: op.osc_prob_5nu_sun(5.0*MEV, R_SUN, 0.0, nu_i=0, nu_f=0,
                                                   average=True))
        P3 = op.osc_prob_3nu_sun(5.0*MEV, R_SUN, 0.0, nu_i=0, nu_f=0, average=True)
    assert np.ndim(P) == 0 and abs(float(P) - float(P3)) < 1e-6


def _H(diag, couple=None):
    def H(e, l):
        M = np.diag(np.asarray(diag, dtype=complex)*(1.0 + 0.1*np.sin(l)))
        if couple:
            i, j, c = couple
            M[i, j] = M[j, i] = c
        return M
    return H


def test_only_a_decoupled_flavor_degenerate_along_the_path_is_split_off():
    # Flavors 3 and 4 decoupled and degenerate with each other: every decoupled flavor, 2
    # included, is split off.
    H = _H([1.0, 2.0, 3.0, 5.0, 5.0], couple=(0, 1, 0.3))
    assert op._decoupled_degenerate_flavors(H, 1.0, 0.0, 10.0) == ([0, 1], [2, 3, 4])
    # Decoupled but not degenerate with anything: left as it is (the 4nu Sun at its defaults).
    H = _H([1.0, 2.0, 3.0, 5.0], couple=(0, 1, 0.3))
    assert op._decoupled_degenerate_flavors(H, 1.0, 0.0, 10.0) is None
    # Degenerate but coupled: left as it is.
    H = _H([1.0, 2.0, 5.0, 5.0], couple=(2, 3, 0.3))
    assert op._decoupled_degenerate_flavors(H, 1.0, 0.0, 10.0) is None


def test_the_4nu_sun_at_its_defaults_is_not_split():
    """Its sterile state crosses an active level at a point, which the windows handle (#59)."""
    seen = {}

    def spy(htot, only, *a, **k):
        seen['htot'], seen['E'] = htot, float(np.ravel(a[0])[0])
        raise KeyboardInterrupt
    real = op._avg_prob_dispatch
    op._avg_prob_dispatch = spy
    try:
        op.osc_prob_4nu_sun(5.0*MEV, R_SUN, 0.0, average=True)
    except KeyboardInterrupt:
        pass
    finally:
        op._avg_prob_dispatch = real
    assert op._decoupled_degenerate_flavors(seen['htot'], seen['E'], 0.0, R_SUN) is None


def test_the_grid_route_stops_at_its_node_cap(monkeypatch):
    """A backstop for a persistent degeneracy the split does not catch: past the cap the phase
    average raises, and the dispatch returns the decohered limit with PhaseAveragingWarning.
    The 70 km shock fronts of test_avgprob take the grid route, on 81 nodes."""
    import test_avgprob as ta
    P_full, _ = ta._shock_average(70.0)
    monkeypatch.setattr(ap, 'PHASE_AVERAGE_MAX_GRID_NODES', 80)
    P_capped, caught = ta._shock_average(70.0)
    assert op.PhaseAveragingWarning in caught
    assert np.isfinite(P_capped) and P_capped != P_full
