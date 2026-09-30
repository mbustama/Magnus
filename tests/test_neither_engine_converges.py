# -*- coding: utf-8 -*-
"""When neither engine of strategy='auto' reaches the tolerance (issue #167).

At 1 MeV across the Sun, rtol = atol = 1e-8 is beyond what the adiabatic hybrid certifies, so
'auto' declines it and the Magnus ladder answers instead; the ladder then needs about 2e6 slabs
and stops at max_n_slabs = 20000 on a single level, with nothing to compare it against.  On 1.1
that level was returned: 4.0e-3 off, 450 times worse than the default tolerance's 8.8e-6.  The
hybrid's own uncertified answer, computed and thrown away on the way, was 1.7e-6 off.

Now the per-point ladder returns the hybrid's answer when it has less to show for its own, and
says so.  The reference here is the hybrid certified at 1e-5 (7.5e-6 off the 2e6-slab ladder),
because the fine ladder takes 12 s; every assertion leaves room for that.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

MEV = gd.UNIT_MEV
R_SUN = gd.SUN_RADIUS*gd.UNIT_KM
OSC = {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k]
       for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}


def _sun(energy_mev, **kw):
    return np.asarray(op.osc_prob_3nu_sun(energy_mev*MEV, R_SUN, 0.0, **OSC, **kw))


def _certified_reference(energy_mev, **kw):
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = _sun(energy_mev, rtol=1e-5, atol=1e-5, strategy='hybrid', strategy_info=info, **kw)
    assert info['certified'] is True
    return P


def test_a_tight_tolerance_is_no_longer_worse_than_the_default():
    """The reproducer of the issue: 4.0e-3 off on 1.1."""
    info = {}
    with pytest.warns(op.ToleranceNotAchievedWarning, match='neither engine'):
        P = _sun(1.0, rtol=1e-8, atol=1e-8, strategy_info=info)
    assert np.max(np.abs(P - _certified_reference(1.0))) < 3e-5
    assert info['engine'] == 'hybrid'
    assert info['certified'] is False
    assert 0.0 < info['trace'][-1]['error_estimate'] < 1e-4


def test_the_default_tolerance_never_reaches_the_fallback(monkeypatch):
    """Read only after the ladder fails: a call that converges pays nothing and is unchanged."""
    calls = []
    real = op._hybrid_fallback
    monkeypatch.setattr(op, '_hybrid_fallback', lambda *a: calls.append(a) or real(*a))
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _sun(1.0, strategy_info=info)
        _sun(1.0, rtol=1e-8, atol=1e-8, strategy='magnus', min_n_slabs=4, max_n_slabs=8)
    assert info['engine'] == 'hybrid' and info['certified'] is True
    assert calls == [(1.0*MEV, R_SUN, 0.0)]      # the forced ladder asked, and found nothing


def test_a_forced_ladder_still_answers_with_the_ladder():
    """strategy='magnus' asked for the ladder: no hybrid ran, so there is nothing to return."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _sun(1.0, rtol=1e-8, atol=1e-8, strategy='magnus', min_n_slabs=4, max_n_slabs=8,
             strategy_info=info)
    assert info['engine'] == 'magnus'


def test_a_note_left_by_an_earlier_answer_is_not_read():
    """Nested probes share one trace: a decline behind an answer belongs to a finished call."""
    fallback = dict(energy=1.0, L=2.0, L0=0.0, P=np.eye(3), error_estimate=1e-6)
    with op._engine_probe() as trace:
        op._note_engine('hybrid', answered=False, _fallback=fallback)
        assert op._hybrid_fallback(1.0, 2.0, 0.0) is fallback
        assert op._hybrid_fallback(1.0, 3.0, 0.0) is None
        op._note_engine('magnus')
        assert op._hybrid_fallback(1.0, 2.0, 0.0) is None
    assert all(not k.startswith('_') for e in op._summarize_engine_trace(trace)['trace']
               for k in e)


def test_a_scan_that_starts_at_the_slab_cap_says_it_checked_nothing():
    """Issue #184: across the Sun at 1e-6 every energy's slab count starts at max_n_slabs, so
    the energy-batched ladder computes one level and returns it unverified.  The warning says
    that, and where to look instead, rather than only that a cap was reached."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        _sun(np.array([1.0, 2.0, 5.0, 10.0]), rtol=1e-6, atol=1e-6)
    messages = [str(w.message) for w in caught
                if issubclass(w.category, op.ToleranceNotAchievedWarning)]
    assert any('starts at max_n_slabs' in m and "strategy='hybrid'" in m for m in messages)
