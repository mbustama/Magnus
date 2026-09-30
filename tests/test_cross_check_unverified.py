# -*- coding: utf-8 -*-
"""The cross-check says which engines did not reach the tolerance (issue #166).

On the Sun at 1 MeV and the default tolerance, the forced Magnus ladder needs more than its
20000-slab cap and stops on one level, 1.9e-3 off, and warns ToleranceNotAchievedWarning; the
certified hybrid is 3.8e-5 off.  max_spread_independent named the pair (hybrid, magnus) and
left the reader to find out which of the two had already said it had failed.
"""

import warnings


import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.oscprob as op

OSC = {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k]
       for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
R_SUN = gd.SUN_RADIUS*gd.UNIT_KM


def _check(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return op.cross_check_strategies(*args, **kwargs)


def test_an_engine_that_ran_out_of_room_is_named_and_left_out_of_the_verified_spread():
    out = _check(op.osc_prob_3nu_sun, 1.0*gd.UNIT_MEV, R_SUN, 0.0,
                 density_profile='B16-GS98', **OSC)
    assert 'magnus' in out['unverified'] and 'hybrid' not in out['unverified']
    assert out['max_spread_independent'] > 1e-3                 # the spread itself is unchanged
    assert 'magnus' not in (out['max_spread_verified_pair'] or ())
    assert out['max_spread_verified'] < 1e-3
    assert set(out['unverified']) <= set(out['ran'])


def test_when_every_engine_converged_nothing_is_unverified_and_the_spreads_agree():
    ne = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    out = _check(op.osc_prob_matter_std_potential, 2, ne, 10.0e6, 0.5*R_SUN,
                 {'sth': OSC['s12'], 'Dm2': OSC['D21']}, L0=0.0,
                 density_is_of_number_of_electrons=True)
    assert out['unverified'] == ()
    assert out['max_spread_verified'] == out['max_spread_independent']
    assert out['max_spread_verified_pair'] == out['max_spread_independent_pair']
