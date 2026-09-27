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
