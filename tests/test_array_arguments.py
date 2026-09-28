r"""Arguments other than energy and L refuse arrays with a message that names them (issue #116).

Only ``energy`` and ``L`` take arrays.  An array anywhere else used to fail deep in the call
with a NumPy message ("The truth value of an array ... is ambiguous", "setting an array element
with a sequence", "unhashable type") that named neither the argument nor the rule.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

E = 1.0*gd.UNIT_GEV
L = 1300.0*gd.UNIT_KM
KM = gd.UNIT_KM
TWO = np.array([0.1, 0.2])

CASES = {
    's13': (lambda: op.osc_prob_3nu_vacuum(E, L, s13=TWO), 's13'),
    'sth': (lambda: op.osc_prob_2nu_vacuum(E, L, sth=TWO, Dm2=2.5e-3), 'sth'),
    's14': (lambda: op.osc_prob_4nu_vacuum(E, L, s14=TWO), 's14'),
    's15': (lambda: op.osc_prob_5nu_vacuum(E, L, s15=TWO), 's15'),
    'eps_em': (lambda: op.osc_prob_3nu_matter_nsi_constant_density(
        E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, eps_em=TWO), 'eps_em'),
    'b3': (lambda: op.osc_prob_3nu_vacuum_liv(E, L, b3=TWO*1.0e-22), 'b3'),
    'nu_i': (lambda: op.osc_prob_3nu_vacuum(E, L, nu_i=np.array([0, 1]), nu_f=0), 'nu_i'),
    'rho': (lambda: op.osc_prob_3nu_matter_constant_density(
        E, L, rho=np.array([2.0, 3.0]), density_matter_is_in_g_per_cm3=True), 'rho'),
    'rho_central': (lambda: op.osc_prob_3nu_matter_exp_density(
        E, L, 0.0, np.array([3.0, 4.0]), 500.0*KM, density_matter_is_in_g_per_cm3=True),
        'rho_central'),
    'costhz': (lambda: op.osc_prob_3nu_earth(
        E, costhz=np.array([-0.3, -0.5]), L=np.array([3822.6, 6371.0])*KM), 'costhz'),
    'electron_fraction': (lambda: op.osc_prob_3nu_earth(
        E, costhz=-0.5, L=6371.0*KM, electron_fraction=np.array([0.5, 0.49])),
        'electron_fraction'),
}


@pytest.mark.parametrize('case', sorted(CASES))
def test_array_argument_is_named(case):
    call, name = CASES[case]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises(ValueError, match=name + r'.* single number.*Only energy and L'):
            call()


def test_both_named_when_two_are_arrays():
    with pytest.raises(ValueError, match=r's13, dCP are not single numbers'):
        op.osc_prob_3nu_vacuum(E, L, s13=TWO, dCP=[0.0, 1.0])


def test_scalar_arguments_unaffected():
    # The checks sit on the failing path; a valid call returns what it always did.
    P = op.osc_prob_3nu_vacuum(E, L, s13=0.15)
    assert P.shape == (3, 3)
    assert np.allclose(P.sum(axis=1), 1.0)
