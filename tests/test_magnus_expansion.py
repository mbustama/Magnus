# -*- coding: utf-8 -*-
"""Tests of the Magnus-expansion core (magnus.magnus).

The reference implementation of the Magnus terms in this file is coded
independently from the production code, with explicit parentheses around
every commutator group, following the Bernoulli-number recursion of
Blanes, Casas, Oteo & Ros, Phys. Rep. 470, 151 (2009), Eq. (2.16).
"""

import numpy as np
import pytest
import scipy as sp
from scipy.integrate import cumulative_trapezoid, solve_ivp

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
