# -*- coding: utf-8 -*-
"""``average=True`` on a constant Hamiltonian at small phase (issue #163).

The closed-form route computed the phase average and then kept the zero-phase limit wherever
the two differed by less than 1e-4 in absolute terms, whatever the tolerance: at 1 GeV over 1
to 3 km it returned about 1e-32 for P(numu -> nue) = 4.7e-7, in silence.  The gate is now
relative as well as absolute, and follows rtol and atol as on the profile route.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

OSC = dict(s12=0.55, s23=0.7, s13=0.15, dCP=1.0, D21=7.4e-5, D31=2.5e-3)
KW = dict(nu_i=gd.NUMU, nu_f=gd.NUE, **OSC)
E = 1.0*gd.UNIT_GEV
WRAPPERS = [lambda L, **k: op.osc_prob_3nu_vacuum(E, L, **KW, **k),
            lambda L, **k: op.osc_prob_3nu_matter_constant_density(
                E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **KW, **k)]


def P(f, L_km, **k):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return float(f(L_km*gd.UNIT_KM, **k))


@pytest.mark.parametrize('f', WRAPPERS, ids=['vacuum', 'constant density'])
def test_a_narrow_spread_returns_the_unaveraged_probability_at_small_phase(f):
    for L_km in (1.0, 3.0):                      # 1e-32 before issue #163
        assert P(f, L_km, average=True, average_spread=1e-3) == pytest.approx(P(f, L_km), rel=1e-3)


@pytest.mark.parametrize('f', WRAPPERS, ids=['vacuum', 'constant density'])
def test_the_default_spread_is_continuous_in_baseline(f):
    """At the default spread the average sits about 1% above the unaveraged value at these
    phases; the ratio has to stay smooth instead of dropping to zero below some baseline."""
    ratios = np.array([P(f, L, average=True)/P(f, L) for L in np.geomspace(1.0, 30.0, 25)])
    assert np.all(ratios > 0.99) and np.all(ratios < 1.03)
    assert np.max(np.abs(np.diff(ratios))) < 1e-3


def test_large_phase_averages_are_unchanged():
    """Where the limit was right, the gate still returns it bit for bit."""
    assert P(WRAPPERS[0], 1.0e4, average=True) == P(WRAPPERS[0], 1.0e4, average=True,
                                                     rtol=1e-3, atol=1e-3)
