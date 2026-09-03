# -*- coding: utf-8 -*-
"""Tests of the symbolic Magnus-term generator (magnus.expansionterms).

The point of this module is to check the numerical core's hard-coded coefficients against
a derivation that shares no code with them, so most of these tests are cross-checks rather
than assertions about the generator alone.
"""

from fractions import Fraction
from math import comb

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid, solve_ivp

import magnus.expansionterms as et
import magnus.magnus as mg

DIM = 3
RNG = np.random.default_rng(17)


def maxabs(x):
    return np.max(np.abs(np.asarray(x)))


# ----------------------------------------------------------------------
# Bernoulli numbers
# ----------------------------------------------------------------------

def test_bernoulli_matches_classical_values():
    """B_1 = -1/2 convention, exact rationals, odd B_n above 1 vanishing."""
    assert et.bernoulli(0) == Fraction(1)
    assert et.bernoulli(1) == Fraction(-1, 2)
    assert et.bernoulli(2) == Fraction(1, 6)
    assert et.bernoulli(4) == Fraction(-1, 30)
    assert et.bernoulli(6) == Fraction(1, 42)
    assert et.bernoulli(8) == Fraction(-1, 30)
    assert all(et.bernoulli(n) == 0 for n in (3, 5, 7, 9, 11))


def test_bernoulli_rejects_negative_index():
    with pytest.raises(ValueError):
        et.bernoulli(-1)


def test_hard_coded_group_factors_match_the_derivation():
    """The whole reason this module exists: magnus.py's F1..F4 are typed-in constants, and
    nothing else checked them against the recursion they are supposed to come from."""
    assert mg.F1 == pytest.approx(float(et.bernoulli_factor(2)), rel=0, abs=0)
    assert mg.F2 == pytest.approx(float(et.bernoulli_factor(4)), rel=0, abs=0)
    assert mg.F3 == pytest.approx(float(et.bernoulli_factor(6)), rel=0, abs=0)
    assert mg.F4 == pytest.approx(float(et.bernoulli_factor(8)), rel=0, abs=0)
    # and the -1/2 that multiplies the first group
    assert et.bernoulli_factor(1) == Fraction(-1, 2)


# ----------------------------------------------------------------------
# Structure of the generated expansion
# ----------------------------------------------------------------------

@pytest.mark.parametrize("n", range(1, 11))
def test_term_count_matches_the_composition_formula(n):
    """The j-th group of Omega_n has one term per composition of n-1 into j parts, so
    C(n-2, j-1) terms, summed over the j with a nonzero Bernoulli number."""
    if n == 1:
        assert et.count_terms(1) == 1
        return
    expected = sum(comb(n - 2, j - 1) for j in range(1, n)
                   if et.bernoulli_factor(j) != 0)
    assert et.count_terms(n) == expected


def test_low_order_terms_are_the_published_ones():
    """Omega_2 and Omega_3 written out, as a plain-sight check on the machinery."""
    assert et.omega_terms(2) == ((Fraction(-1, 2), ('c', ('Om', 1), 'A')),)
    assert set(et.omega_terms(3)) == {
        (Fraction(-1, 2), ('c', ('Om', 2), 'A')),
        (Fraction(1, 12), ('c', ('Om', 1), ('c', ('Om', 1), 'A'))),
    }


def test_every_term_is_a_nested_chain_ending_in_A():
    """Each term is [Om_m1, [Om_m2, ... [Om_mj, A]]] with the indices summing to n-1.
    If that ever stops holding, the composition-driven path in magnus.py is wrong."""
    for n in range(2, 9):
        for _, word in et.omega_terms(n):
            indices = []
            while word != 'A':
                assert word[0] == 'c', f"unexpected node {word[0]!r}"
                left, word = word[1], word[2]
                assert left[0] == 'Om'
                indices.append(left[1])
            assert sum(indices) == n - 1, (n, indices)


def test_print_magnus_terms_runs_and_mentions_every_order(capsys):
    et.print_magnus_terms(5)
    out = capsys.readouterr().out
    for n in range(1, 6):
        assert f"Omega_{n}" in out
    assert "1/12" in out


# ----------------------------------------------------------------------
# Cross-check against the numerical core
# ----------------------------------------------------------------------

def _sampled_A(n_points=4001, width=0.6):
    """A smooth, non-commuting matrix function sampled on the normalized grid."""
    s = np.linspace(0.0, 1.0, n_points)
    H = [RNG.normal(size=(DIM, DIM)) for _ in range(3)]
    Bt = np.array([width*(H[0] + np.sin(3*x)*H[1] + x**2*H[2]) for x in s], dtype=complex)
    return s, Bt


def _terms_from_generator(s, Bt, order):
    """Evaluates the generated symbolic terms numerically, sharing no code with magnus.py
    beyond the commutator itself."""
    def integ(y):
        return cumulative_trapezoid(y, x=s, axis=0, initial=0)

    def comm(X, Y):
        return X @ Y - Y @ X

    def evaluate(word, om):
        if word == 'A':
            return Bt
        if word[0] == 'Om':
            return om[word[1]]
        return comm(evaluate(word[1], om), evaluate(word[2], om))

    om = {1: integ(Bt)}
    for n in range(2, order + 1):
        acc = np.zeros_like(Bt)
        for coeff, word in et.omega_terms(n):
            acc = acc + float(coeff)*evaluate(word, om)
        om[n] = integ(acc)
    return [om[k][-1] for k in range(1, order + 1)]


@pytest.mark.parametrize("order", range(1, 11))
def test_implementation_matches_the_generator(order):
    """Every Omega_k the numerical core produces must equal the independently generated
    one.  This covers the hand-written orders 1-6 and the composition-driven 7-10 in the
    same check, so neither can drift from the recursion unnoticed."""
    s, Bt = _sampled_A()
    from_code = mg._magnus_terms_quadrature(Bt, order, 'trapezoid')
    from_gen = _terms_from_generator(s, Bt, order)
    for k in range(order):
        scale = max(maxabs(from_code[k]), 1e-30)
        assert maxabs(from_code[k] - from_gen[k])/scale < 1e-11, f"Omega_{k+1}"


# ----------------------------------------------------------------------
# The orders above 6 actually converge faster
# ----------------------------------------------------------------------

def test_higher_orders_converge_faster():
    """A non-circular check on the new coefficients: comparing against an ODE ground truth,
    order 8 must reach a visibly smaller error than order 6, which must beat order 4.  A
    wrong coefficient at order 7 or 8 would break the ordering even though the symbolic
    cross-check above (which shares the recursion) would still pass."""
    H0 = RNG.normal(size=(DIM, DIM)) + 1j*RNG.normal(size=(DIM, DIM))
    H0 = H0 + H0.conj().T
    H1 = RNG.normal(size=(DIM, DIM)) + 1j*RNG.normal(size=(DIM, DIM))
    H1 = H1 + H1.conj().T

    def A(t):
        t = np.asarray(t)
        if t.ndim:
            return -1j*(H0 + np.sin(t)[..., None, None]*H1)
        return -1j*(H0 + np.sin(float(t))*H1)

    h = 0.35
    sol = solve_ivp(lambda t, y: (A(t) @ y.reshape(DIM, DIM)).ravel(), (0.0, h),
                    np.eye(DIM, dtype=complex).ravel(), rtol=3e-14, atol=1e-16,
                    method='DOP853')
    exact = sol.y[:, -1].reshape(DIM, DIM)

    errs = {}
    for order in (4, 6, 8):
        U = mg.magnus_expansion(A, 0.0, h, n_tpts=20001, order=order,
                                integration_method='trapezoid', validate_input=False)
        errs[order] = maxabs(U - exact)
    assert errs[6] < errs[4], errs
    assert errs[8] < errs[6], errs


# ----------------------------------------------------------------------
# Guards around the order ceiling
# ----------------------------------------------------------------------

def test_gl_refuses_orders_it_has_no_scheme_for():
    """'gl' tops out at order 8: the Gauss-Legendre collocation schemes are separately
    derived, not products of the Magnus recursion.  Silently computing order 8 for a
    higher request would be a quiet wrong answer."""
    def A(t):
        return -1j*np.eye(DIM, dtype=complex)

    for order in range(mg.MAGNUS_EXP_ORDER_MAX_GL + 1, mg.MAGNUS_EXP_ORDER_MAX + 1):
        with pytest.raises(ValueError, match="Gauss-Legendre|gl"):
            mg.magnus_expansion(A, 0.0, 1.0, order=order, integration_method='gl')


def test_gl_ceiling_holds_even_without_input_validation():
    """The guard in _validate() is skipped when validate_input=False, so _gl_nodes carries
    a backstop; without it that flag would reopen the silent-degradation path."""
    def A(t):
        return -1j*np.eye(DIM, dtype=complex)

    with pytest.raises(ValueError):
        mg.magnus_expansion(A, 0.0, 1.0, order=9, integration_method='gl',
                            validate_input=False)


def test_high_order_quadrature_warns_about_cost():
    """Orders above 6 cost 2.7x to 17x order 6 per slab, so the trade is worth flagging."""
    def A(t):
        return -1j*np.eye(DIM, dtype=complex)

    with pytest.warns(mg.MagnusHighOrderCostWarning):
        mg.magnus_expansion(A, 0.0, 1.0, n_tpts=51, order=7,
                            integration_method='trapezoid')


def test_order_six_and_below_does_not_warn():
    def A(t):
        return -1j*np.eye(DIM, dtype=complex)

    with warnings_as_errors():
        mg.magnus_expansion(A, 0.0, 1.0, n_tpts=51, order=6,
                            integration_method='trapezoid')


def warnings_as_errors():
    import warnings
    ctx = warnings.catch_warnings()
    ctx.__enter__()
    warnings.simplefilter("error", mg.MagnusHighOrderCostWarning)
    return _Closer(ctx)


class _Closer:
    def __init__(self, ctx):
        self.ctx = ctx

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.ctx.__exit__(*exc)
        return False


def test_magnus_exp_order_max_is_defined_once():
    """globaldefs re-exports magnus.py's constant rather than repeating it; the two used to
    be separate literals kept in step by hand."""
    import magnus.globaldefs as gd
    assert gd.MAGNUS_EXP_ORDER_MAX is mg.MAGNUS_EXP_ORDER_MAX
    assert mg.MAGNUS_EXP_ORDER_MAX == 10


def test_format_term_can_omit_the_coefficient():
    """``with_coeff=False`` prints the commutator word alone, which is what
    a caller assembling its own coefficient column wants. The default keeps
    the sign and magnitude in front."""
    term = et.omega_terms(2)[0]

    body = et.format_term(term, with_coeff=False)
    full = et.format_term(term)

    assert body, "the formatted word is empty"
    assert not body.startswith(('+', '-')), "the coefficient was not omitted"
    assert full.endswith(body), "the two forms disagree about the word itself"
    assert full[0] in '+-', "the default form dropped the sign"
