r"""Energy scans accept every form of density function the per-point path accepts (issue #113).

The energy-batched engine samples the potential on arrays of positions.  A density written
as a constant, ``lambda l: 3.0*UNIT_G_PER_CM3``, returned one number for such an array and
crashed the scan with an IndexError; one written for a single position, with ``float(l)``,
crashed it with a TypeError.  Both work point by point, so a scan must return the same.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

OSC = gd.load_nufit_params('NuFIT 6.1')
ENERGIES = np.logspace(-1, 1, 5)*gd.UNIT_GEV
L = 1000.0*gd.UNIT_KM
RHO = 3.0*gd.UNIT_G_PER_CM3

FORMS = {
    'scalar_return': (lambda l: RHO, 'separable'),
    # Evaluated position by position by the scenario function since issue #144 §4, so the
    # Hamiltonian built on it is array-capable and the scan stays on the batched engine too.
    'scalar_only': (lambda l: RHO*(1.0 if float(l) < 0.5*L else 2.0/3.0), 'separable'),
}


@pytest.mark.parametrize('form', sorted(FORMS))
@pytest.mark.parametrize('tol', [{}, dict(rtol=None, atol=None)], ids=['default_tol', 'no_tol'])
def test_scan_matches_points(form, tol):
    rho_func, engine = FORMS[form]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        info = {}
        scan = np.asarray(op.osc_prob_matter_std_potential(
            3, rho_func, ENERGIES, L, osc_params=OSC, strategy_info=info, **tol))
        points = np.array([op.osc_prob_matter_std_potential(
            3, rho_func, e, L, osc_params=OSC, **tol) for e in ENERGIES])
    assert scan.shape == points.shape
    assert np.max(np.abs(scan - points)) < 1.0e-12
    # A constant written as a scalar stays on the batched engine, and so does one that takes a
    # single position at a time.
    assert info['engine'] == engine
