# -*- coding: utf-8 -*-
"""Tests of the Magnus-expansion core (magnus.magnus).

The reference implementation of the Magnus terms in this file is coded
independently from the production code, with explicit parentheses around
every commutator group, following the Bernoulli-number recursion of
Blanes, Casas, Oteo & Ros, Phys. Rep. 470, 151 (2009), Eq. (2.16).
"""

import warnings

import numpy as np
import pytest
import scipy as sp
from scipy.integrate import cumulative_trapezoid, solve_ivp

import magnus.magnus as mg
from magnus.magnus import (MagnusConvergenceWarning, commutator,
                           magnus_expansion, magnus_expansion_multislab)

RNG = np.random.default_rng(42)
DIM = 3


def random_hermitian(dim, rng=RNG):
    X = rng.standard_normal((dim, dim)) + 1j*rng.standard_normal((dim, dim))
    return 0.5*(X + X.conj().T)


# Three independent generators: for a two-generator Hamiltonian of the form
# A(t) = a + f(t) b, the term [Omega_2, [Omega_1, A]] vanishes identically,
# which would mask coefficient errors in Omega_4 (this is how the original
# operator-precedence bug escaped simple two-flavor checks).
H0 = random_hermitian(DIM)
H1 = random_hermitian(DIM)
H2 = random_hermitian(DIM)


def A_scalar(t):
    """Scalar-only matrix function (raises on array input)."""
    if np.ndim(t) != 0:
        raise TypeError("scalar input only")
    return -1j*(H0 + np.sin(2.0*t)*H1 + np.cos(3.0*t)*H2)


def A_vec(t):
    """Vectorized matrix function (accepts scalars and arrays)."""
    t = np.asarray(t)
    return -1j*(H0 + np.sin(2.0*t)[..., None, None]*H1
                + np.cos(3.0*t)[..., None, None]*H2)


def exact_U(A, t0, t1, dim=DIM):
    """High-accuracy solution of U' = A(t) U via adaptive ODE integration."""
    def rhs(t, y):
        return (A(t) @ y.reshape(dim, dim)).ravel()
    sol = solve_ivp(rhs, (t0, t1), np.eye(dim, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853')
    return sol.y[:, -1].reshape(dim, dim)


def maxabs(x):
    return np.max(np.abs(x))


def reference_magnus_terms(A, t0, t1, n_tpts, order):
    """Independent implementation of Omega_1..Omega_order (trapezoid)."""
    times = np.linspace(t0, t1, n_tpts)
    ic = lambda y: cumulative_trapezoid(y, x=times, axis=0, initial=0)
    At = np.array([A(t) for t in times])
    f1, f2 = 1.0/12.0, -1.0/720.0
    terms = []
    o1t = ic(At)
    terms.append(o1t[-1])
    if order >= 2:
        o2t = ic(-0.5*commutator(o1t, At))
        terms.append(o2t[-1])
    if order >= 3:
        o3t = ic(-0.5*commutator(o2t, At)
                 + f1*commutator(o1t, commutator(o1t, At)))
        terms.append(o3t[-1])
    if order >= 4:
        o4t = ic(-0.5*commutator(o3t, At)
                 + f1*(commutator(o1t, commutator(o2t, At))
                       + commutator(o2t, commutator(o1t, At))))
        terms.append(o4t[-1])
    if order >= 5:
        o5t = ic(-0.5*commutator(o4t, At)
                 + f1*(commutator(o1t, commutator(o3t, At))
                       + commutator(o2t, commutator(o2t, At))
                       + commutator(o3t, commutator(o1t, At)))
                 + f2*commutator(o1t, commutator(o1t, commutator(o1t,
                     commutator(o1t, At)))))
        terms.append(o5t[-1])
    if order >= 6:
        o6t = ic(-0.5*commutator(o5t, At)
                 + f1*(commutator(o1t, commutator(o4t, At))
                       + commutator(o2t, commutator(o3t, At))
                       + commutator(o3t, commutator(o2t, At))
                       + commutator(o4t, commutator(o1t, At)))
                 + f2*(commutator(o1t, commutator(o1t, commutator(o1t,
                           commutator(o2t, At))))
                       + commutator(o1t, commutator(o1t, commutator(o2t,
                           commutator(o1t, At))))
                       + commutator(o1t, commutator(o2t, commutator(o1t,
                           commutator(o1t, At))))
                       + commutator(o2t, commutator(o1t, commutator(o1t,
                           commutator(o1t, At))))))
        terms.append(o6t[-1])
    return terms


@pytest.mark.parametrize("order", [1, 2, 3, 4, 5, 6])
def test_terms_match_independent_reference(order):
    """Omega_k from the code must match the independently coded recursion.

    Same quadrature and same grid, so any difference is a term/coefficient
    bug (this is the regression test for the operator-precedence bug that
    dropped the 1/12 and -1/720 factors at orders >= 4, and for the
    order-6 array rebinding bug).

    ``integration_method`` is pinned to 'trapezoid' to match
    ``reference_magnus_terms``, which is a trapezoid implementation.  It also
    has to be a quadrature method at all: 'gl' (the package default) is
    commutator-free and never forms the individual Omega_k, returning only
    their sum, so there would be nothing term-by-term to compare."""
    _, terms_code = magnus_expansion(A_scalar, 0.0, 1.0, n_tpts=1001,
                                     order=order, return_magnus_terms=True,
                                     integration_method='trapezoid')
    terms_ref = reference_magnus_terms(A_scalar, 0.0, 1.0, 1001, order)
    for k in range(order):
        assert maxabs(terms_code[k] - terms_ref[k]) < 1e-10


def test_error_decreases_with_order():
    """Against an exact ODE solution, each order must improve on the last."""
    Uex = exact_U(A_vec, 0.0, 0.6)
    errs = []
    for order in range(1, 7):
        U = magnus_expansion(A_vec, 0.0, 0.6, n_tpts=1501, order=order,
                             integration_method='simpson')
        errs.append(maxabs(U - Uex))
    for k in range(5):
        assert errs[k+1] < 0.9*errs[k], f"order {k+2} did not improve: {errs}"
    assert errs[5] < 1e-3


def test_unitarity_for_hermitian_hamiltonian():
    U = magnus_expansion(A_vec, 0.0, 0.6, n_tpts=201, order=4)
    assert maxabs(U @ U.conj().T - np.eye(DIM)) < 1e-12


@pytest.mark.parametrize("method", ['trapezoid', 'simpson', 'gl'])
def test_vectorized_and_scalar_A_agree(method):
    """The silent vectorization of A must not change the result at all."""
    Uv = magnus_expansion(A_vec, 0.0, 0.6, n_tpts=201, order=4,
                          integration_method=method)
    Us = magnus_expansion(A_scalar, 0.0, 0.6, n_tpts=201, order=4,
                          integration_method=method)
    assert maxabs(Uv - Us) < 1e-14


@pytest.mark.parametrize("order,expected_ratio", [(2, 4.0), (4, 16.0),
                                                  (6, 64.0)])
def test_gl_convergence_rates(order, expected_ratio):
    """Gauss-Legendre integrators must converge at their nominal order."""
    Uex = exact_U(A_vec, 0.0, 2.0)
    errs = []
    for n_slabs in [8, 16, 32]:
        grid = np.linspace(0.0, 2.0, n_slabs + 1)
        edges = np.column_stack([grid[:-1], grid[1:]])
        Uc = magnus_expansion_multislab(A_vec, edges, order=order,
                                        integration_method='gl')
        Utot = Uc[0]
        for k in range(1, n_slabs):
            Utot = Uc[k] @ Utot
        errs.append(maxabs(Utot - Uex))
    for k in range(2):
        ratio = errs[k]/errs[k+1]
        assert 0.55*expected_ratio < ratio < 2.0*expected_ratio, \
            f"errors {errs}, ratio {ratio}"


@pytest.mark.parametrize("method", ['trapezoid', 'simpson', 'gl'])
def test_multislab_matches_per_slab_calls(method):
    grid = np.linspace(0.0, 1.0, 9)
    edges = np.column_stack([grid[:-1], grid[1:]])
    Um = magnus_expansion_multislab(A_vec, edges, n_tpts_per_slab=41,
                                    order=4, integration_method=method)
    Ul = np.array([magnus_expansion(A_scalar, e0, e1, n_tpts=41, order=4,
                                    integration_method=method)
                   for e0, e1 in edges])
    assert maxabs(Um - Ul) < 1e-13


def test_zero_width_slab_gives_identity():
    Um = magnus_expansion_multislab(A_vec, [[0.0, 0.5], [0.5, 0.5]], order=4)
    assert maxabs(Um[1] - np.eye(DIM)) < 1e-14


def test_constant_A_is_exact():
    """For constant A the result must equal expm(A T) to machine precision."""
    Aconst = -1j*H0

    def A(t):
        return Aconst

    T = 2.5
    U = magnus_expansion(A, 0.0, T, n_tpts=2, order=1)
    assert maxabs(U - sp.linalg.expm(Aconst*T)) < 1e-12


def test_non_antihermitian_A_falls_back_to_expm():
    def A_gen(t):
        return np.array([[1.0*t, 2.0 + 3j*t], [2.0 - 3j*t, 2.0]])
    U = magnus_expansion(A_gen, 0.0, 0.3, n_tpts=801, order=4)
    Uex = exact_U(A_gen, 0.0, 0.3, dim=2)
    assert maxabs(U - Uex) < 5e-3


@pytest.mark.parametrize("bad_kwargs", [dict(order=0), dict(order=7),
                                        dict(integration_method='nope')])
def test_invalid_input_raises_value_error(bad_kwargs):
    with pytest.raises(ValueError):
        magnus_expansion(A_vec, 0.0, 1.0, **bad_kwargs)


def test_convergence_warning_fires_for_wide_slab():
    with pytest.warns(MagnusConvergenceWarning):
        magnus_expansion(A_vec, 0.0, 10.0, n_tpts=301, order=2)


def test_no_convergence_warning_for_constant_A():
    """The series terminates at Omega_1 for constant A: no warning."""
    Aconst = -1j*10.0*H0

    def A(t):
        return Aconst

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", MagnusConvergenceWarning)
        magnus_expansion(A, 0.0, 5.0, n_tpts=2, order=1)


# ----------------------------------------------------------------------
# Vectorization of the Hamiltonian (ScalarHamiltonianWarning)
# ----------------------------------------------------------------------

def _h_pieces():
    """A constant part and a position-dependent potential, for the H below."""
    h0 = np.diag([0.0, 1.0e-12, 2.0e-12])
    e00 = np.diag([1.0, 0.0, 0.0])
    return h0, e00


def test_scalar_only_hamiltonian_warns_and_names_the_fix():
    """The scalar fallback is correct but slow, and silent without this."""
    h0, e00 = _h_pieces()

    def H_scalar(l):
        return h0 + float(np.exp(-l)) * e00      # float() rejects an array

    with pytest.warns(mg.ScalarHamiltonianWarning, match='accepts an array'):
        mg.magnus_expansion(lambda t: -1j * H_scalar(t), 0.0, 1.0,
                            n_tpts=8, order=2, integration_method='trapezoid')


def test_array_capable_hamiltonian_does_not_warn():
    h0, e00 = _h_pieces()

    def H_array(l):
        l = np.asarray(l, dtype=float)
        return h0 + np.exp(-l)[..., None, None] * e00

    with warnings.catch_warnings():
        warnings.simplefilter('error', mg.ScalarHamiltonianWarning)
        mg.magnus_expansion(lambda t: -1j * H_array(t), 0.0, 1.0,
                            n_tpts=8, order=2, integration_method='trapezoid')


def test_constant_hamiltonian_does_not_warn():
    """A Hamiltonian that ignores its argument is detected and broadcast, so
    it is already on a fast path and must not be told otherwise."""
    h0, _ = _h_pieces()

    with warnings.catch_warnings():
        warnings.simplefilter('error', mg.ScalarHamiltonianWarning)
        mg.magnus_expansion(lambda t: -1j * h0, 0.0, 1.0,
                            n_tpts=8, order=2, integration_method='trapezoid')


def test_scalar_and_array_hamiltonians_give_the_same_answer():
    """The warning is about speed only: the two paths must agree exactly."""
    h0, e00 = _h_pieces()

    def H_scalar(l):
        return h0 + float(np.exp(-l)) * e00

    def H_array(l):
        l = np.asarray(l, dtype=float)
        return h0 + np.exp(-l)[..., None, None] * e00

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', mg.ScalarHamiltonianWarning)
        a = mg.magnus_expansion(lambda t: -1j * H_scalar(t), 0.0, 1.0,
                                n_tpts=32, order=4, integration_method='trapezoid')
    b = mg.magnus_expansion(lambda t: -1j * H_array(t), 0.0, 1.0,
                            n_tpts=32, order=4, integration_method='trapezoid')
    assert np.allclose(a, b, rtol=0.0, atol=0.0)


def test_wrong_vector_mode_hint_falls_back_and_warns():
    """A caller-supplied A_eval_mode='vector' that turns out to be wrong must
    degrade safely rather than return a mis-shaped result."""
    h0, e00 = _h_pieces()

    def H_scalar(l):
        return h0 + float(np.exp(-l)) * e00

    with pytest.warns(mg.ScalarHamiltonianWarning):
        At, mode = mg._evaluate_A(lambda t: -1j * H_scalar(t),
                                  np.linspace(0.0, 1.0, 5), 'vector')
    assert mode == 'scalar'
    assert At.shape == (5, 3, 3)


# ----------------------------------------------------------------------
# ordered_product: the time-ordered chain, collapsed pairwise
# ----------------------------------------------------------------------

def test_ordered_product_matches_the_sequential_product():
    """The tree reduction must equal reduce(np.matmul, U[::-1]) exactly in
    value.  Associativity is what licenses it; commutativity is neither
    required nor assumed, so the operator ordering has to survive."""
    from functools import reduce
    rng = np.random.default_rng(20260809)
    for n in (1, 2, 3, 7, 8, 33, 108):
        U = (rng.normal(size=(n, 3, 3)) + 1j*rng.normal(size=(n, 3, 3)))/3.0
        expected = reduce(np.matmul, U[::-1]) if n > 1 else U[0]
        got = mg.ordered_product(U)
        scale = max(np.max(np.abs(expected)), 1.0e-300)
        assert np.max(np.abs(got - expected))/scale < 1.0e-12, n


def test_ordered_product_preserves_order_not_just_content():
    """A product of non-commuting factors pins the ordering down: if the tree
    combined them out of turn, the value would change."""
    A = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    B = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    # U[0] is earliest, so the product is B @ A, not A @ B (which differ here).
    got = mg.ordered_product(np.stack([A, B]))
    np.testing.assert_allclose(got, B @ A, atol=0.0)
    assert not np.allclose(B @ A, A @ B)


def test_ordered_product_keeps_a_unitary_chain_unitary():
    """The slab operators are unitary, and their product must stay so: this is
    the property the probabilities inherit."""
    from scipy.linalg import expm as sp_expm
    rng = np.random.default_rng(7)
    U = np.empty((256, 3, 3), dtype=complex)
    for i in range(U.shape[0]):
        A = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        U[i] = sp_expm(-1j*(A + A.conj().T)/2*0.2)
    P = mg.ordered_product(U)
    assert np.max(np.abs(P @ P.conj().T - np.eye(3))) < 1.0e-13


# ----------------------------------------------------------------------
# A Hamiltonian that is constant *within* each slab
# ----------------------------------------------------------------------

def test_constant_slab_samples_drop_the_vanishing_commutators():
    """Where a slab's nodes are identical the Magnus series terminates at its
    first term, every later one being a commutator of A with itself.  Taking
    Omega = h A there is exact, not an approximation, so it must agree with the
    full expression to the last bit -- the commutator it skips is zero."""
    rng = np.random.default_rng(3)
    n, d = 32, 3
    A = rng.normal(size=(n, d, d)) + 1j*rng.normal(size=(n, d, d))
    widths = np.full(n, 0.01)

    for order, n_nodes in ((4, 2), (6, 3)):
        An = np.repeat(A[:, None, :, :], n_nodes, axis=1)
        h = widths[:, None, None]
        got = mg._magnus_gl(An, widths, order)
        np.testing.assert_array_equal(got, h*A)


def test_a_varying_slab_still_gets_the_full_expansion():
    """The detection must be exact equality, not a tolerance: a nearly constant
    A has commutators that carry real information, and dropping them because two
    samples agreed to some epsilon would silently lower the order."""
    rng = np.random.default_rng(4)
    n, d = 16, 3
    An = rng.normal(size=(n, 2, d, d)) + 1j*rng.normal(size=(n, 2, d, d))
    widths = np.full(n, 0.01)
    h = widths[:, None, None]
    expected = (0.5*h*(An[:, 0] + An[:, 1])
                + (np.sqrt(3.0)/12.0)*h*h*mg.commutator(An[:, 1], An[:, 0]))
    np.testing.assert_allclose(mg._magnus_gl(An, widths, 4), expected,
                               rtol=0.0, atol=1.0e-15)

    # And a sample pair differing in the last bit is *not* treated as constant.
    An[3, 1, 0, 0] = np.nextafter(An[3, 0, 0, 0].real, 1.0) + 0j
    An[3, 0, 0, 0] = An[3, 0, 0, 0].real + 0j
    assert not mg._samples_identical(An[:, 0], An[:, 1])


def test_order_six_constant_shortcut_needs_all_three_nodes_equal():
    r"""``A1 == A2 != A3`` must take the full expansion, not the constant shortcut.

    ``_magnus_gl``'s order-6 branch fires only when *every* Gauss-Legendre node of
    the slab carries the same sample, and the conjunction is what makes that true:

    .. code-block:: python

        if _samples_identical(A1, A2) and _samples_identical(A2, A3):

    Nothing else in the suite distinguishes that ``and`` from an ``or``.  The two
    existing tests make all three nodes bit-identical, so both operands are true
    either way, and the varying-slab test is order 4 with two nodes.  Weakened to
    ``or``, this configuration returns ``h*A1`` -- a first-order answer for a slab
    that is not constant -- and measures 0.0122 against an ``h*A`` scale of 0.01,
    an O(100%) error.

    The pattern is physically reachable rather than contrived: a piecewise-constant
    profile with a jump between the second and third node is a
    ``t_breakpoints``-delimited region boundary, or a castle wall, which is the
    case the shortcut exists to serve.
    """
    A = np.array([[0.0, 1.0e-13], [1.0e-13, 3.0e-13]], dtype=complex)
    B = np.array([[0.0, 4.0e-13], [4.0e-13, 9.0e-13]], dtype=complex)
    h = 0.01

    An = np.stack([A, A, B])[None, ...]          # (1 slab, 3 nodes, 2, 2)
    Om = mg._magnus_gl(An, np.array([h]), 6)[0]
    shortcut = h*A

    # The shortcut must NOT have fired: the third node differs.
    assert np.max(np.abs(Om - shortcut)) > 1e-3*np.max(np.abs(shortcut)), (
        "order 6 took the constant shortcut on a slab whose third node differs; "
        "check that _magnus_gl still requires ALL nodes identical")

    # And where all three DO agree, it fires and is exact.
    Ac = np.stack([A, A, A])[None, ...]
    Om_c = mg._magnus_gl(Ac, np.array([h]), 6)[0]
    assert np.array_equal(Om_c, h*A)


# ----------------------------------------------------------------------
# The evaluation-mode cache and the interval it was learned on.
#
# `probe_eval_mode`'s verdict is a property of the function *on an interval*,
# not of the function alone: 'constant' means "sampling this A across [t0, t1]
# gave the same matrix every time".  A cache keyed on the function alone
# therefore answers a question it was never asked, and `_evaluate_A`'s
# 'constant' branch applies the hint by evaluating once and broadcasting, with
# no spot-check -- unlike the 'vector' and 'scalar' hints, which self-heal.
#
# The result is a silent wrong answer that depends on the order of a loop.
# ----------------------------------------------------------------------

def _layered_A():
    r"""A profile with the short-circuit any careful author would write.

    Two layers; if every requested position falls in one of them, return that
    layer's single matrix rather than a stack.  Entirely reasonable, and it makes
    `probe_eval_mode` answer 'constant' on a short interval and 'vector' on a
    long one -- for the same function object.
    """
    lo = np.array([[0.0, 1.0e-13], [1.0e-13, 2.0e-13]], dtype=complex)
    hi = np.array([[0.0, 8.0e-13], [8.0e-13, 5.0e-13]], dtype=complex)
    edge = 5.0e12

    def A(t):
        t = np.asarray(t, dtype=float)
        if t.ndim == 0:
            return -1j*(lo if t < edge else hi)
        if np.all(t < edge):
            return -1j*lo
        if np.all(t >= edge):
            return -1j*hi
        return -1j*np.where((t < edge)[..., None, None], lo, hi)

    return A


def test_probe_eval_mode_verdict_depends_on_the_interval():
    r"""The premise: the same function is 'constant' on one interval, not another."""
    A = _layered_A()
    assert mg.probe_eval_mode(A, 0.0, 4.0e12) == 'constant'
    assert mg.probe_eval_mode(A, 0.0, 2.0e13) != 'constant'


def test_cached_eval_mode_does_not_reuse_a_mode_across_intervals():
    r"""A mode learned on a short interval must not be served for a long one.

    Keyed on the function alone, the first call teaches the cache 'constant' and
    the second is answered from it -- so the whole trajectory is propagated from
    one Hamiltonian sample, and the answer is unitary, unwarned, and wrong.
    """
    A = _layered_A()
    short = mg.cached_eval_mode(A, 0.0, 4.0e12)
    long_cached = mg.cached_eval_mode(A, 0.0, 2.0e13)
    long_fresh = mg.probe_eval_mode(A, 0.0, 2.0e13)
    assert short == 'constant'
    assert long_cached == long_fresh, (
        "cached_eval_mode returned %r for [0, 2e13] after learning %r on "
        "[0, 4e12]; probe_eval_mode says %r. The cache key is missing the "
        "interval." % (long_cached, short, long_fresh))


def test_evaluate_A_is_not_poisoned_by_a_mode_from_a_shorter_interval():
    r"""The consequence, at the level a caller would see it.

    Order the two calls shortest-first and the long one inherits 'constant';
    the samples it then uses are one matrix broadcast over a profile that
    genuinely varies.
    """
    A = _layered_A()
    times = np.linspace(0.0, 2.0e13, 9)
    mg.cached_eval_mode(A, 0.0, 4.0e12)          # teach it the short interval
    mode = mg.cached_eval_mode(A, 0.0, 2.0e13)
    At, used = mg._evaluate_A(A, times, mode)
    ref, _ = mg._evaluate_A(A, times, None)      # probe honestly
    assert np.max(np.abs(At - ref)) < 1e-30, (
        "samples built from a cached mode differ from honestly probed ones by "
        "%.3e" % np.max(np.abs(At - ref)))


@pytest.mark.parametrize('kind', ['slots', 'dataclass', 'custom_eq', 'plain'])
def test_cached_eval_mode_accepts_hamiltonians_it_cannot_cache(kind):
    r"""A cache must not decide whether a call is allowed to succeed.

    ``WeakKeyDictionary`` raises ``TypeError`` for a key that cannot be weakly
    referenced (a ``__slots__`` class) or cannot be hashed (a ``@dataclass``,
    which sets ``__hash__ = None``, or anything defining ``__eq__`` without one).
    All three are ordinary ways to write a Hamiltonian, and all three worked
    before this cache existed, so the lookup is guarded and simply probes.

    The guard has to wrap the *lookup*, not only the store: an earlier copy of
    this cache in ``oscprob`` guarded the store alone and turned every such
    Hamiltonian into a hard ``TypeError`` from inside ``osc_prob``.
    """
    from dataclasses import dataclass

    H = np.diag([0.0, 2.5e-13, 7.5e-13]).astype(complex)

    class Slotted:
        __slots__ = ('H',)

        def __init__(self, H):
            self.H = H

        def __call__(self, l):
            return -1j*self.H

    @dataclass
    class Data:
        H: object

        def __call__(self, l):
            return -1j*self.H

    class CustomEq:
        def __init__(self, H):
            self.H = H

        def __eq__(self, other):
            return isinstance(other, CustomEq)

        __hash__ = None

        def __call__(self, l):
            return -1j*self.H

    A = {'slots': Slotted(H), 'dataclass': Data(H), 'custom_eq': CustomEq(H),
         'plain': (lambda l: -1j*H)}[kind]
    mode = mg.cached_eval_mode(A, 0.0, 1.0e13)
    assert mode == mg.probe_eval_mode(A, 0.0, 1.0e13)
    # and calling twice must not raise on the store path either
    assert mg.cached_eval_mode(A, 0.0, 1.0e13) == mode
