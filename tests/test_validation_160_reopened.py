r"""The #160 checklist items #173 did not implement, one test per item of the reopened issue.

Each test reproduces the case the audit (``docs/dev/audit/audit160.py``) found failing on
``main`` at 0939598, and requires what the checklist asks: a Magnus error naming the argument,
a warning, or the documented behaviour.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

E = 1.0*gd.UNIT_GEV
L = 1000.0*gd.UNIT_KM
KM = gd.UNIT_KM
OSC = gd.load_nufit_params('NuFIT 6.1')


def _refused(match, fn, *args, **kwargs):
    with pytest.raises(ValueError) as err:
        fn(*args, **kwargs)
    msg = str(err.value)
    assert 'Error in magnus' in msg and match in msg, msg
    return msg


# §11 angles beyond 90 degrees -------------------------------------------------------------

@pytest.mark.parametrize('s12, angles', [(2.0, 'rad'), (-2.0, 'rad'), (120.0, 'deg'),
                                         (-100.0, 'deg')])
def test_angle_with_negative_cosine_is_refused(s12, angles):
    kw = dict(s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3)
    if angles == 'deg':
        kw.update(s23=40.0, s13=8.6, dCP=200.0)
    msg = _refused('s12', op.osc_prob_3nu_vacuum, E, L, s12=s12, angles=angles, **kw)
    assert 'negative cosine' in msg


@pytest.mark.parametrize('s12, angles', [(-1.0, 'rad'), (np.pi/2, 'rad'), (90.0, 'deg'),
                                         (-90.0, 'deg'), (33.8, 'deg')])
def test_angle_in_first_or_fourth_quadrant_is_accepted(s12, angles):
    kw = dict(s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3)
    if angles == 'deg':
        kw.update(s23=40.0, s13=8.6, dCP=200.0)
    P = op.osc_prob_3nu_vacuum(E, L, s12=s12, angles=angles, **kw)
    assert np.allclose(np.sum(P, axis=-1), 1.0)


# §11 angle-convention warning ---------------------------------------------------------------

def test_angle_warning_judges_by_the_active_angles():
    """Sines in the active slots are not hidden by 5-degree sterile angles."""
    with pytest.warns(gd.MixingAngleConventionWarning):
        op.osc_prob_4nu_vacuum(E, L, s12=0.55, s23=0.75, s13=0.15, dCP=1.0, s14=5.0, s24=5.0,
                               s34=5.0, d14=0.0, d24=0.0, D21=7.5e-5, D31=2.5e-3, D41=1.0,
                               angles='deg')


def test_angle_warning_says_how_to_silence_a_false_positive():
    with pytest.warns(gd.MixingAngleConventionWarning, match='filterwarnings'):
        op.osc_prob_2nu_vacuum(E, L, sth=0.5, Dm2=7.5e-5, angles='deg')


def test_angle_warning_is_quiet_for_real_degrees():
    with warnings.catch_warnings():
        warnings.simplefilter('error', gd.MixingAngleConventionWarning)
        op.osc_prob_4nu_vacuum(E, L, s12=33.7, s23=43.3, s13=8.6, dCP=212.0, s14=0.5, s24=0.5,
                               s34=0.5, d14=0.0, d24=0.0, D21=7.5e-5, D31=2.5e-3, D41=1.0,
                               angles='deg')
