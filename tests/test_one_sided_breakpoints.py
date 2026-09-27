"""One-sided endpoint samples at declared breakpoints ('trapezoid'/'simpson').

Cumulative quadrature samples each slab at both of its ends.  At a density jump the profile
returns one side's value there, so without care the slab on the other side integrates a sample
from the wrong layer: an error proportional to 1/m, whatever the rule.  With the breakpoint
declared, that sample is taken just inside its own slab.  On a piecewise-constant profile whose
layers are slabs, every rule is then exact, which is what these tests check against: the matrix
exponential of each layer, or the same call on 'gl' (Gauss-Legendre samples no endpoints).
"""
import warnings

import numpy as np
import pytest
import scipy as sp

import magnus.globaldefs as gd
import magnus.magnus as mg
import magnus.oscprob as op

RNG = np.random.default_rng(71)


def _herm(d, scale):
    X = RNG.normal(size=(d, d)) + 1j*RNG.normal(size=(d, d))
    return scale*(X + X.conj().T)/2.0


H0, M = _herm(3, 0.4), _herm(3, 0.4)


def _A_step(jumps, values):
    """-i(H0 + v(t) M), v piecewise constant; at a jump it returns the right-hand value."""
    jumps, values = np.asarray(jumps, float), np.asarray(values, float)

    def A(t):
        v = values[np.searchsorted(jumps, np.asarray(t, float), side='right')]
        return -1j*(H0 + v[..., None, None]*M)
    return A


def _exact(edges, A):
    U = np.eye(3, dtype=complex)
    for a, b in edges:
        U = sp.linalg.expm(A(0.5*(a + b))*(b - a)) @ U
    return U


def _chain(U):
    return mg.ordered_product(U)


EDGES_1 = np.array([[0.0, 0.5], [0.5, 1.0]])
A_1 = _A_step([0.5], [1.0, 3.0])


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
@pytest.mark.parametrize('m', [3, 5, 9])
def test_core_is_exact_on_layers_once_the_jump_is_declared(method, m):
    kw = dict(n_tpts_per_slab=m, order=4, integration_method=method)
    exact = _exact(EDGES_1, A_1)
    without = _chain(mg.magnus_expansion_multislab(A_1, EDGES_1, **kw))
    with_bp = _chain(mg.magnus_expansion_multislab(A_1, EDGES_1, t_breakpoints=[0.5], **kw))
    assert np.max(np.abs(without - exact)) > 1e-4          # the 1/m error the fix removes
    assert np.max(np.abs(with_bp - exact)) < 1e-12


def test_core_without_breakpoints_and_on_gl_is_unchanged():
    kw = dict(n_tpts_per_slab=5, order=4, integration_method='simpson')
    a = mg.magnus_expansion_multislab(A_1, EDGES_1, **kw)
    b = mg.magnus_expansion_multislab(A_1, EDGES_1, t_breakpoints=None, **kw)
    assert np.array_equal(a, b)
    g = mg.magnus_expansion_multislab(A_1, EDGES_1, order=4)
    h = mg.magnus_expansion_multislab(A_1, EDGES_1, order=4, t_breakpoints=[0.5])
    assert np.array_equal(g, h)


def _counting(A):
    seen = [0]

    def f(t):
        seen[0] += np.size(t)
        return A(t)
    return f, seen


EDGES_SYM = np.array([[0.0, 0.3], [0.3, 0.5], [0.5, 0.7], [0.7, 1.0]])
A_SYM = _A_step([0.3, 0.7], [1.0, 3.0, 1.0])


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
def test_mirror_keeps_working_with_palindromic_breakpoints(method, monkeypatch):
    kw = dict(n_tpts_per_slab=5, order=4, integration_method=method, A_eval_mode='vector')
    f, seen = _counting(A_SYM)
    mirrored = mg.magnus_expansion_multislab(f, EDGES_SYM, symmetric_over=(0.0, 1.0),
                                             t_breakpoints=[0.3, 0.7], **kw)
    n_mirrored = seen[0]
    monkeypatch.setattr(mg, 'USE_PALINDROME', False)
    seen[0] = 0
    plain = mg.magnus_expansion_multislab(f, EDGES_SYM, symmetric_over=(0.0, 1.0),
                                          t_breakpoints=[0.3, 0.7], **kw)
    assert n_mirrored < seen[0]                            # the mirror still halves the samples
    assert np.max(np.abs(mirrored - plain)) < 1e-14
    assert np.max(np.abs(_chain(mirrored) - _exact(EDGES_SYM, A_SYM))) < 1e-12


def test_non_palindromic_breakpoints_take_the_unmirrored_path(monkeypatch):
    kw = dict(n_tpts_per_slab=5, order=4, integration_method='simpson', A_eval_mode='vector')
    one = mg.magnus_expansion_multislab(A_SYM, EDGES_SYM, symmetric_over=(0.0, 1.0),
                                        t_breakpoints=[0.3], **kw)
    monkeypatch.setattr(mg, 'USE_PALINDROME', False)
    plain = mg.magnus_expansion_multislab(A_SYM, EDGES_SYM, t_breakpoints=[0.3], **kw)
    assert np.array_equal(one, plain)


def _warned(bp):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        mg.magnus_expansion_multislab(A_1, EDGES_1, n_tpts_per_slab=5, order=4,
                                      integration_method='simpson', t_breakpoints=bp)
    return [c for c in caught if issubclass(c.category, op.UnmarkedDiscontinuityWarning)]


def test_a_breakpoint_inside_a_slab_warns():
    assert len(_warned([0.25])) == 1
    assert len(_warned([0.25, 0.5])) == 1


@pytest.mark.parametrize('bp', [[0.5], [np.nextafter(0.5, 1.0)], [np.nextafter(0.5, 0.0)],
                                [0.0, 1.0], [-1.0, 2.0], []])
def test_breakpoints_on_edges_or_outside_the_chain_are_quiet(bp):
    assert _warned(bp) == []


@pytest.mark.parametrize('bad', [[[0.5]], [np.nan], [0.5 + 0.1j]])
def test_breakpoints_are_validated(bad):
    with pytest.raises(ValueError):
        mg.magnus_expansion_multislab(A_1, EDGES_1, n_tpts_per_slab=5, order=4,
                                      integration_method='simpson', t_breakpoints=bad)


# --- through the wrappers: per-point, energy-batched and cumulative engines ----------------

L_TOT = 3000.0*gd.UNIT_KM
OSC = {'s12': 0.55, 's23': 0.68, 's13': 0.15, 'dCP': 3.7, 'D21': 7.5e-5, 'D31': 2.5e-3}


def _rho(l):
    r = np.where(np.asarray(l) < 0.5*L_TOT, 3.0, 8.0)
    return float(r) if r.ndim == 0 else r


def _call(energy, method, **kw):
    info = {}
    P = op.osc_prob_matter_std_potential(3, _rho, energy, kw.pop('L', L_TOT), OSC,
                                         density_matter_is_in_g_per_cm3=True, strategy='magnus',
                                         t_breakpoints=[0.5*L_TOT], integration_method=method,
                                         strategy_info=info, **kw)
    return np.asarray(P), info.get('engine')


FIXED = dict(rtol=None, atol=None, n_slabs=2, n_tpts_per_slab=5, magnus_exp_order=4)


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
@pytest.mark.parametrize('energy, engine', [(1.0*gd.UNIT_GEV, 'magnus'),
                                            (np.array([0.5, 1.0, 2.0, 4.0])*gd.UNIT_GEV,
                                             'separable')])
def test_per_point_and_batched_engines_are_exact_on_layers(method, energy, engine):
    P, eng = _call(energy, method, **FIXED)
    ref, _ = _call(energy, 'gl', **FIXED)
    assert eng == engine
    assert np.max(np.abs(P - ref)) < 1e-10


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
def test_cumulative_scan_is_exact_on_layers(method):
    L = np.array([0.25, 0.5, 0.75, 1.0])*L_TOT
    P, eng = _call(1.0*gd.UNIT_GEV, method, L=L, cumulative=True, n_tpts_per_slab=5)
    ref, _ = _call(1.0*gd.UNIT_GEV, 'gl', L=L, cumulative=True)
    assert eng == 'cumulative'
    assert np.max(np.abs(P - ref)) < 1e-10


@pytest.mark.parametrize('method', ['simpson', 'trapezoid'])
def test_batched_ladder_meets_its_tolerance_on_layers(method):
    """Exercises the phase seed on the quadrature rules and the nudge on every level."""
    E = np.array([0.5, 1.0, 2.0, 4.0])*gd.UNIT_GEV
    P, eng = _call(E, method, rtol=1e-6, atol=1e-6)
    ref, _ = _call(E, 'gl', **FIXED)
    assert eng == 'separable'
    assert np.max(np.abs(P - ref)) < 1e-5
