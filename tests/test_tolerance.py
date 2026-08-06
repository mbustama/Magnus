# -*- coding: utf-8 -*-
r"""Tests for what ``rtol``/``atol`` promise, and for the effective-refinement gate.

``rtol``/``atol`` are a stopping criterion -- the ladder halts when two successive levels
agree -- not an estimate of the error.  That is only sound if the two levels compared are
genuinely different grids, and with ``t_breakpoints`` re-inserted at every level they need
not be: on an Earth chord a nominal 2 -> 3 slab step is a 16 -> 17 edge step.  Two grids
differing by 6% agree for reasons unrelated to convergence, and the ladder was stopping on
that, returning answers outside the tolerance asked for without warning.
"""

import warnings

import numpy as np
import pytest

from magnus import oscprob, earth, globaldefs as gd


REF_SLABS = 6000


def _call(costhz, energy, tol=None, n_slabs=None):
    r"""One Earth probability, plus what the ladder reported about it."""
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    info = {}
    kw = dict(costhz=costhz, L=L, validate_input=False, strategy='magnus',
              convergence_info=info)
    if n_slabs is None:
        kw.update(rtol=tol, atol=tol)
    else:
        kw.update(rtol=None, atol=None, n_slabs=n_slabs)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = np.asarray(oscprob.osc_prob_3nu_earth(energy, **kw), dtype=float)
    return P, info, sorted({w.category.__name__ for w in caught})


# --------------------------------------------------------- the reported facts

def test_convergence_info_reports_the_two_levels_compared():
    r"""Not just the level returned: without the previous one a caller cannot see that a
    nominal 6 -> 9 slab step was really a 20 -> 23 edge step."""
    _, info, _ = _call(-0.8, 2.0*gd.UNIT_GEV, tol=1.0e-3)
    for key in ('n_slabs', 'n_slabs_previous', 'n_slab_edges',
                'n_slab_edges_previous', 'last_gap', 'n_agreements',
                'tolerance_achieved'):
        assert key in info, key
    assert info['n_slab_edges'] > info['n_slabs']          # breakpoints were inserted
    assert info['n_slab_edges_previous'] < info['n_slab_edges']
    assert info['tolerance_achieved'] is True
    assert info['last_gap'] is not None


def test_last_gap_is_none_when_only_one_level_was_computed():
    r"""``None`` means "nothing was compared", which is not the same as "the gap was small"."""
    _, info, _ = _call(-0.8, 2.0*gd.UNIT_GEV, n_slabs=40)
    assert info['last_gap'] is None
    assert info['n_slabs_previous'] is None
    assert info['tolerance_achieved'] is None              # none was requested


def test_tolerance_achieved_is_false_when_the_ladder_runs_out_of_room():
    r"""A tight tolerance is not enough to force this -- ``gl`` on an Earth chord reaches
    1e-12 comfortably -- so the room itself is capped."""
    L = earth.distance_traveled_inside_earth(-0.8)*gd.UNIT_KM
    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        oscprob.osc_prob_3nu_earth(
            2.0*gd.UNIT_GEV, costhz=-0.8, L=L, validate_input=False,
            strategy='magnus', rtol=1.0e-10, atol=1.0e-10, max_n_slabs=12,
            convergence_info=info)
    names = sorted({w.category.__name__ for w in caught})
    assert info['tolerance_achieved'] is False
    assert 'ToleranceNotAchievedWarning' in names


# ------------------------------------------------- the effective-refinement gate

def test_the_gate_constant_is_the_measured_one():
    assert oscprob.MIN_EFFECTIVE_REFINEMENT == 1.25


@pytest.mark.parametrize('costhz, energy_gev, tol', [
    (-0.45, 2.0, 1.0e-5),        # the configuration that was silently wrong
    (-0.8, 1.0, 1.0e-3),
    (-0.6, 3.0, 1.0e-3),
    (-0.95, 8.0, 1.0e-4),
])
def test_a_converged_answer_is_actually_within_the_requested_tolerance(costhz, energy_gev, tol):
    r"""The property ``rtol``/``atol`` is supposed to deliver, checked against a reference.

    Every configuration here missed its tolerance before the effective-refinement gate; the
    first missed it *silently*, claiming success with no warning, by 2.1x.
    """
    energy = energy_gev*gd.UNIT_GEV
    ref, _, _ = _call(costhz, energy, n_slabs=REF_SLABS)
    P, info, names = _call(costhz, energy, tol=tol)
    if not info['tolerance_achieved']:
        pytest.skip('ladder did not claim convergence; the warning is the contract there')
    allowed = tol + tol*float(np.max(np.abs(ref)))
    assert float(np.max(np.abs(P - ref))) <= allowed


def test_an_agreement_across_a_barely_refined_grid_does_not_count():
    r"""The mechanism itself: two levels whose *edge* counts barely differ must not be
    allowed to certify, however well they agree."""
    _, info, _ = _call(-0.45, 2.0*gd.UNIT_GEV, tol=1.0e-5)
    if info['tolerance_achieved']:
        ratio = info['n_slab_edges']/info['n_slab_edges_previous']
        assert ratio >= oscprob.MIN_EFFECTIVE_REFINEMENT


def test_the_gate_is_inert_without_breakpoints(monkeypatch):
    r"""No breakpoints means the edge count is the slab count, so the ratio is
    ``growth_factor_n_slabs`` and the gate can never fire.  Solar results must be
    bit-identical whatever the constant is set to."""
    from magnus import matter
    profile = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    params = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    p2 = {'sth': params['s12'], 'Dm2': params['D21']}

    def run():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.asarray(oscprob.osc_prob_matter_std_potential(
                2, profile, 5.0e6, 0.3*gd.SUN_RADIUS*gd.UNIT_KM, p2, L0=0.0,
                density_is_of_number_of_electrons=True, strategy='magnus'), dtype=float)

    before = run()
    monkeypatch.setattr(oscprob, 'MIN_EFFECTIVE_REFINEMENT', 1.0)
    assert np.array_equal(run(), before)
    monkeypatch.setattr(oscprob, 'MIN_EFFECTIVE_REFINEMENT', 1.25)
    assert np.array_equal(run(), before)
    assert oscprob.MIN_EFFECTIVE_REFINEMENT < 1.5
