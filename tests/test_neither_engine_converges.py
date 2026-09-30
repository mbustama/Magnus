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
    real = op._hybrid_rescue
    monkeypatch.setattr(op, '_hybrid_rescue', lambda *a: calls.append(a) or real(*a))
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
    kept = dict(P=np.eye(3), certified=False, error_estimate=1e-6)

    class Asked(Exception):
        pass

    def H_at_energy(energy):
        raise Asked(energy)
    rescue = dict(H_at_energy=H_at_energy, L0=0.0, rtol=1e-8, atol=1e-8, magnus_exp_order=4,
                  integration_method='gl', answers={(1.0, 2.0): dict(kept, full=True)})
    with op._engine_probe() as trace:
        op._note_engine('hybrid', answered=False, _rescue=rescue)
        assert op._hybrid_rescue(1.0, 2.0, 0.0)['P'] is kept['P']
        assert op._hybrid_rescue(1.0, 2.0, 5.0) is None          # another starting point
        with pytest.raises(Asked):                                # a point it never reached
            op._hybrid_rescue(1.0, 3.0, 0.0)
        op._note_engine('magnus')
        assert op._hybrid_rescue(1.0, 2.0, 0.0) is None
    assert all(not k.startswith('_') for e in op._summarize_engine_trace(trace)['trace']
               for k in e)


# Issue #184: a scan of several energies is answered by the energy-batched engine, which across
# the Sun seeds every energy at the 20000-slab cap and returned one unverifiable level for all
# of them: 2.0e-3 off at every rtol = atol from 1e-4 to 1e-8, over 0.5-20 MeV.
SCAN = np.array([1.0, 2.0, 5.0, 10.0])


def _hybrid_reference(energies_mev):
    """The hybrid, certified at 3e-4 at each energy: within 7.0e-5 of the 2e6-slab ladder."""
    out = []
    for e in energies_mev:
        info = {}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out.append(_sun(e, rtol=3e-4, atol=3e-4, strategy='hybrid', strategy_info=info))
        assert info['certified'] is True
    return np.array(out)


def test_a_tight_scan_is_answered_by_the_better_supported_engine_and_says_so():
    """1 MeV was 4.0e-3 off and 10 MeV 7.9e-4 off; 2 and 5 MeV, 3e-5 and 5e-5, keep the ladder."""
    info = {}
    with pytest.warns(op.ToleranceNotAchievedWarning, match='neither engine'):
        P = _sun(SCAN, rtol=1e-6, atol=1e-6, strategy_info=info)
    assert np.max(np.abs(P - _hybrid_reference(SCAN))) < 2e-4
    hybrid = [e for e in info['trace'] if e['engine'] == 'hybrid' and e['answered']]
    assert hybrid and hybrid[0]['certified'] is False and 1 <= hybrid[0]['n_points'] <= len(SCAN)


def test_where_the_capped_ladder_is_right_it_is_kept():
    """B16-GS98: the ladder's one capped level is 1e-6 to 8e-6 off, the hybrid's one-iteration
    answer 2e-3 to 8e-3, and the hybrid's own estimate says so.  Nothing changes there."""
    energies = np.array([1.5, 6.0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        auto = _sun(energies, rtol=1e-6, atol=1e-6, density_profile='B16-GS98')
        ladder = _sun(energies, rtol=1e-6, atol=1e-6, density_profile='B16-GS98',
                      strategy='magnus')
    assert np.array_equal(auto, ladder)


def test_energies_the_hybrid_certified_keep_their_certified_answer():
    """At 1e-4 the hybrid certifies 1 and 2 MeV and declines the scan at 10 MeV."""
    energies = np.array([1.0, 2.0, 10.0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = _sun(energies, rtol=1e-4, atol=1e-4)
        for k in (0, 1):
            info = {}
            alone = _sun(energies[k], rtol=1e-4, atol=1e-4, strategy='hybrid', strategy_info=info)
            assert info['certified'] is True
            assert np.array_equal(P[k], alone)
    assert np.max(np.abs(P - _hybrid_reference(energies))) < 2e-4


def test_a_scan_the_hybrid_certifies_never_reaches_the_rescue(monkeypatch):
    calls = []
    real = op._hybrid_rescue_note
    monkeypatch.setattr(op, '_hybrid_rescue_note', lambda *a: calls.append(a) or real(*a))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _sun(SCAN)
    assert calls == []


def test_the_hybrid_reports_the_gap_between_its_last_two_levels():
    import magnus.adiabatic as ad
    H = lambda l: np.diag([0.0, 1e-12, 3e-12]) + 1e-13*np.ones((3, 3))   # noqa: E731
    info = {}
    ad.hybrid_propagator(H, 0.0, 1e12, rtol=1e-6, atol=1e-6, max_iters=1, info=info)
    assert isinstance(info['last_gap'], float) and info['last_gap'] >= 0.0
