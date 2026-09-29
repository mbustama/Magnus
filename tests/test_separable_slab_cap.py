# -*- coding: utf-8 -*-
r"""Tests of the energy-batched engine at its slab cap (issue #71).

Once ``'trapezoid'`` or ``'simpson'`` reach ``max_n_slabs``, each further level refines only the
points per slab.  Two such levels agreeing verifies the quadrature inside each slab, not the
slab count, and the engine accepted them silently: on the multi-resonance test profile at
rtol = atol = 1e-6 the answer was 2.5e-4 off while successive levels agreed to 7e-8.  Such an
acceptance now raises ``ToleranceNotAchievedWarning``, as the per-point ladder already did.
The warning changes no result and no work.

The reference is the same scan with ``'gl'`` at rtol = atol = 1e-10 and a cap far above what it
needs.

``'gl'`` fixes the points per slab, so at the cap it stops with the refinement-caps warning.  Its
last step onto the cap is clamped, though, and can refine the grid by a sliver.  An agreement
across such a step now counts only within the fraction of the tolerance that the step can vouch
for (issue #122); ``test_separable_breakpoints.py`` keeps a capped step that does certify.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

L = 250.0*gd.UNIT_KM
E = np.linspace(0.1, 0.4, 8)*gd.UNIT_GEV


def rho(l):
    return 3e3*np.exp(-np.asarray(l, dtype=float)/(100.0*gd.UNIT_KM))


def scan(tol, **kw):
    r"""The batched scan's probabilities and the tolerance warnings it raised."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = np.asarray(op.osc_prob_matter_std_potential(
            3, rho, E, L, gd.load_nufit_params('NuFIT 6.1'), L0=0.0,
            density_matter_is_in_g_per_cm3=True, rtol=tol, atol=tol, strategy='magnus', **kw))
    return P, [w for w in caught if issubclass(w.category, op.ToleranceNotAchievedWarning)]


@pytest.fixture(scope='module')
def reference():
    P, warned = scan(1e-10, max_n_slabs=200000)
    assert not warned
    return P


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
def test_an_acceptance_at_the_slab_cap_warns_once(reference, method):
    r"""At ``max_n_slabs=8`` the levels past the cap agree while the answer is thousands of
    tolerances off: the scan now says so, once per call."""
    tol = 1e-6
    P, warned = scan(tol, integration_method=method, max_n_slabs=8)
    assert np.max(np.abs(P - reference)) > 1000*tol
    assert len(warned) == 1
    assert 'slab cap' in str(warned[0].message)


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
def test_below_the_default_cap_nothing_warns(reference, method):
    r"""With the default cap the slab count keeps growing until the scan converges: no
    warning, and the answer is within tolerance."""
    tol = 1e-6
    P, warned = scan(tol, integration_method=method)
    assert not warned
    assert np.all(np.abs(P - reference) <= tol + tol*np.abs(reference))


def test_gl_at_its_cap_keeps_its_own_warning():
    r"""``'gl'`` pins the points per slab, so at the cap it stops with the existing
    refinement-caps warning rather than refining in place."""
    _, warned = scan(1e-6, max_n_slabs=8)
    assert len(warned) == 1
    assert 'refinement caps reached' in str(warned[0].message)


def test_gl_refuses_an_agreement_across_the_clamped_step_onto_the_cap():
    r"""On three energies of the paper's five-flavor Listing 1 scan at rtol = 5e-13, floors of 13
    and 19 slabs take the ladder to 19 926 and 18 346 slabs, then clamp it to the cap of 20 000.
    The two grids agreed to 0.4 and 0.8 of the tolerance, beyond the 0.015 and 0.41 that steps
    this short can vouch for, and the scan returned a probability 6.6 times outside the
    tolerance with no warning (issue #122); it now ends in the caps warning."""
    osc = gd.load_nufit_params('NuFIT 6.1')
    sterile = dict(osc, s14=np.sqrt(0.1), s24=np.sqrt(0.1), s15=np.sqrt(0.06),
                   s25=np.sqrt(0.06), D41=1.0, D51=1.7)
    energies = np.logspace(np.log10(2.0), np.log10(20.0), 26)[1:4]*gd.UNIT_GEV
    for floor in (13, 19):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            op.osc_prob_5nu_matter_exp_density(
                energies, L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.0e3,
                l_scale=10.0*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE,
                nu_f=gd.NUE, strategy='magnus', magnus_exp_order=4, rtol=5e-13, atol=5e-15,
                n_slabs=floor, **sterile)
        warned = [w for w in caught if issubclass(w.category, op.ToleranceNotAchievedWarning)]
        assert len(warned) == 1
        assert 'refinement caps reached' in str(warned[0].message)
