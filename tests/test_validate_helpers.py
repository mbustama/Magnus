"""The shared argument checks of magnus._validate (issue #160)."""

import numpy as np
import pytest

from magnus import _validate as v

W = 'oscprob.osc_prob_3nu_vacuum'


@pytest.mark.parametrize('x', [1, 1.5, np.float32(1.5), np.int64(3), np.array(2.0), np.float64(4)])
def test_real_accepts_every_real_scalar(x):
    assert v.check_real('s12', x, W) == float(x)


@pytest.mark.parametrize('x, exc', [(True, TypeError), (np.True_, TypeError), ('a', TypeError),
                                     (None, TypeError), (1 + 2j, TypeError),
                                     (np.nan, ValueError), (np.inf, ValueError),
                                     ([1.0], ValueError), (np.ones(2), ValueError)])
def test_real_refuses_and_names(x, exc):
    with pytest.raises(exc, match=r'osc_prob_3nu_vacuum: D31'):
        v.check_real('D31', x, W)


def test_real_bounds():
    assert v.check_real('Ye', 1.0, W, lo=0.0, hi=1.0) == 1.0
    with pytest.raises(ValueError, match='Ye must be in'):
        v.check_real('Ye', 1.5, W, lo=0.0, hi=1.0)
    with pytest.raises(ValueError, match='positive'):
        v.check_real('l_scale', 0.0, W, positive=True)
    assert v.check_real('rtol', None, W, positive=True, allow_none=True) is None


@pytest.mark.parametrize('x', [True, 2.5, '3', np.nan, 1 + 0j, None])
def test_int_refuses_non_integers(x):
    with pytest.raises(TypeError, match='n_slabs must be an integer'):
        v.check_int('n_slabs', x, W, lo=1)


def test_int_accepts_numpy_and_bounds():
    assert v.check_int('n_slabs', np.int64(4), W, lo=1) == 4
    with pytest.raises(ValueError, match='max_n_slabs must be a positive integer; got -1'):
        v.check_int('max_n_slabs', -1, W, lo=1)


@pytest.mark.parametrize('x', ['False', 'no', 2, 1.0, None, np.array([True, False])])
def test_bool_refuses_truthy_values(x):
    with pytest.raises(TypeError, match='average must be True or False'):
        v.check_bool('average', x, W)


def test_bool_accepts_numpy_bool():
    assert v.check_bool('nubar', np.True_, W) is True


def test_choice_compares_types():
    assert v.check_choice('cumulative', 'auto', W, ('auto', True, False)) == 'auto'
    with pytest.raises(ValueError, match='integration_method must be one of'):
        v.check_choice('integration_method', 'foo', W, ('gl', 'trapezoid', 'simpson'))
    with pytest.raises(ValueError):
        v.check_choice('x', 1, W, (True, False))


def test_real_array_names_the_entry():
    v.check_real_array('L', np.array([1.0, 2.0]), W, nonnegative=True)
    with pytest.raises(ValueError, match='L must be finite; entry 1 is nan'):
        v.check_real_array('L', [1.0, np.nan], W)
    with pytest.raises(ValueError, match='energy must be positive; entry 0'):
        v.check_real_array('energy', [-1.0, 2.0], W, positive=True)
    with pytest.raises(TypeError, match='masked'):
        v.check_real_array('L', np.ma.masked_array([1.0, 2.0], [0, 1]), W)
    with pytest.raises(TypeError, match='complex'):
        v.check_real_array('L', np.array([1j, 2]), W)
    with pytest.raises(ValueError, match='empty'):
        v.check_real_array('L', [], W)


@pytest.mark.parametrize('edges, text', [
    ([[0, 1], [1.5, 2]], 'a gap'),
    ([[0, 1.2], [1, 2]], 'an overlap'),
    ([[0, 1], [1, 1], [1, 2]], 'end <= start'),
    ([[0, 1], [1, 1.5]], 'span the whole path'),
    ([0, 2], 'pairs'),
    ([[0, np.nan]], 'finite'),
    ([], 'pairs'),
])
def test_slab_edges(edges, text):
    with pytest.raises(ValueError, match=text):
        v.check_slab_edges(edges, 0.0, 2.0, W)


def test_slab_edges_valid_partition():
    e = v.check_slab_edges([[0, 1], [1, 2]], 0.0, 2.0, W)
    assert e.shape == (2, 2)


def test_hamiltonian_sample():
    H = np.diag([0.0, 1.0, 2.0])
    assert v.check_hamiltonian_sample('H_func', H, W).dtype == np.complex128
    with pytest.raises(ValueError, match='not Hermitian'):
        v.check_hamiltonian_sample('H_func', np.triu(np.ones((3, 3))), W)
    with pytest.raises(ValueError, match='not finite'):
        v.check_hamiltonian_sample('H_func', np.full((2, 2), np.nan), W)
    with pytest.raises(ValueError, match='square'):
        v.check_hamiltonian_sample('H_func', np.ones((2, 3)), W)
    with pytest.raises(TypeError):
        v.check_hamiltonian_sample('H_func', 'abc', W)


def test_refinement_rules():
    v.check_refinement(W, dict(rtol=1e-3, atol=None, n_slabs=4, max_n_slabs=None, n_jobs=-1))
    for bad in (dict(rtol=0.0), dict(atol=-1.0), dict(rtol=np.nan), dict(rtol=True),
                dict(n_slabs=2.5), dict(max_num_loops=2.5), dict(growth_factor_n_slabs=1.0),
                dict(magnus_exp_order=True), dict(magnus_exp_order=11), dict(n_jobs=0),
                dict(n_jobs=2.5), dict(n_tpts_per_slab=np.inf), dict(cumulative=2),
                dict(integration_method='foo'), dict(strategy_info=[]),
                dict(n_slabs=10, max_n_slabs=5)):
        with pytest.raises((TypeError, ValueError), match=list(bad)[0]):
            v.check_refinement(W, bad)
