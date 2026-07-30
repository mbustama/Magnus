# -*- coding: utf-8 -*-
r"""magnus.py

Compute the time-evolution operator using the Magnus expansion.

This module contains the numerical core of Magnus: routines to compute
the matrix exponential of the Magnus expansion of a (possibly
time-dependent) matrix function A(t), i.e.,

    U(t1, t0) = exp[Omega_1 + Omega_2 + ... + Omega_k] ,

where the terms Omega_k are built from time-ordered integrals of nested
commutators of A(t).  For neutrino oscillations, A(t) = -i H(t), with
H(t) the Hamiltonian, but the routines below work for arbitrary
matrix-valued A(t).

The terms are generated with the standard recursion based on Bernoulli
numbers [1]_ (in the B_1 = -1/2 convention):

    Omega_1(t) = int_0^t A(s) ds
    Omega_n(t) = sum_{j=1}^{n-1} (B_j / j!) int_0^t S_n^(j)(s) ds

with S_n^(j) the sums of nested commutators of the lower-order terms
with A.  Orders 1--6 are implemented (B_3 = B_5 = 0, so those groups
vanish identically).

Two families of methods are available, selected via
``integration_method``:

* ``'trapezoid'`` / ``'simpson'``: sample A(t) on a uniform grid of
  ``n_tpts`` points and evaluate the nested integrals with cumulative
  quadrature.  Fully general, but the quadrature error (O(h^2) or
  O(h^4)) can dominate the Magnus truncation error at high orders.

* ``'gl'``: Gauss-Legendre commutator-free collocation [1]_ [2]_.  For a
  slab of width h it needs only 1, 2, or 3 evaluations of A to reach
  order 2, 4, or 6, respectively, with quadrature error matched to the
  truncation order.  ``n_tpts`` is ignored.  This is the recommended
  method when A(t) is smooth within each slab.

References
----------
.. [1] S. Blanes, F. Casas, J. A. Oteo & J. Ros, "The Magnus expansion
   and some of its applications", Phys. Rep. 470, 151 (2009).
.. [2] S. Blanes, F. Casas & J. Ros, "Improved high order integrators
   based on the Magnus expansion", BIT Numer. Math. 40, 434 (2000).

Routine listings
----------------

    * commutator - Returns [X, Y] = XY - YX
    * probe_eval_mode - Determines how a matrix function can be evaluated
    * suggest_n_slabs - Suggests a starting number of time slabs
    * magnus_expansion - Computes exp(Omega) for a single time slab
    * evolution_operators_from_samples - Evolution operators of a chain
           of slabs from precomputed samples of A
    * gl_nodes - Returns the Gauss-Legendre nodes used by the 'gl' method
    * magnus_expansion_multislab - Computes the evolution operators of
           all time slabs at once, from A directly
    * MagnusConvergenceWarning - Warning class for slabs too wide for
           guaranteed Magnus convergence
"""

__version__ = "2.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import warnings
from typing import Optional, Callable, Union, Tuple, List

import numpy as np
import scipy as sp


class MagnusConvergenceWarning(UserWarning):
    r"""Warns that a time slab may be too wide for the Magnus series.

    The Magnus series is guaranteed to converge when
    int_{t0}^{t1} ||A(t)||_2 dt < pi.  We use ||Omega||_2 >= pi as a
    cheap (necessary, not sufficient) proxy to flag slabs that are
    likely too wide; raising the expansion order will not help in that
    regime -- use more (narrower) slabs instead.  The norm comes for
    free from the eigenvalues already computed for the matrix
    exponential.

    .. versionadded:: 0.10.0
    """


# Bernoulli numbers B_k (negative-B_1 convention), kept for reference;
# only B_1, B_2, and B_4 enter at the orders implemented here (<= 6).
B = {
    0: 1.0, 1: -0.5, 2: 1.0/6.0, 3: 0.0, 4: -1.0/30.0, 5: 0.0, 6: 1.0/42.0,
}

# Multiplicative factors of the commutator groups in the recursion
F1 = 1.0 / 12.0    # B_2 / 2!
F2 = -1.0 / 720.0  # B_4 / 4!

# Backward-compatible aliases
f1 = F1
f2 = F2

# Highest order of the Magnus expansion implemented here
MAGNUS_EXP_ORDER_MAX = 6

# Valid values of integration_method
valid_integration_methods = ['trapezoid', 'simpson', 'gl']

# Gauss-Legendre nodes on [0, 1] used by the 'gl' method
_GL1_NODES = np.array([0.5])
_GL2_NODES = np.array([0.5 - np.sqrt(3.0)/6.0, 0.5 + np.sqrt(3.0)/6.0])
_GL3_NODES = np.array([0.5 - np.sqrt(15.0)/10.0, 0.5,
                       0.5 + np.sqrt(15.0)/10.0])

_HAS_CUMULATIVE_SIMPSON = hasattr(sp.integrate, 'cumulative_simpson')


def commutator(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    r"""Returns the commutator [X, Y] = X Y - Y X.

    Works on single matrices and on stacks of matrices (the matrix
    product broadcasts over all leading axes).

    .. versionadded:: 0.10.0

    Parameters
    ----------
    X : np.ndarray
        Left matrix (or stack of matrices).
    Y : np.ndarray
        Right matrix (or stack of matrices), broadcastable against X.

    Returns
    -------
    np.ndarray
        The commutator X @ Y - Y @ X.
    """
    return X @ Y - Y @ X


def _evaluate_A(A: Callable, times: np.ndarray,
                A_eval_mode: Optional[str] = None) -> Tuple[np.ndarray, str]:
    r"""Evaluate the matrix function A at all requested times.

    Silently tries a single vectorized call, A(times), which is much
    faster than evaluating point by point.  The vectorized result is
    accepted only if it has the expected shape and matches a scalar
    spot-check evaluation; otherwise (or if the vectorized call raises)
    the routine falls back to a per-point loop.  A constant A (one that
    ignores its argument) is detected and broadcast.

    Parameters
    ----------
    A : Callable
        Function of time returning a (d, d) matrix; may optionally
        accept an array of times and return a (..., d, d) stack.
    times : np.ndarray
        Times at which to evaluate A; any shape.
    A_eval_mode : str, optional
        If given ('vector', 'constant', or 'scalar', e.g., from a
        previous call or from :func:`probe_eval_mode`), skip the probe
        and evaluate directly in that mode.  This avoids re-probing A
        (two extra scalar evaluations) on every call.

    Returns
    -------
    (np.ndarray, str)
        Array of shape ``times.shape + (d, d)`` and complex dtype, and
        the evaluation mode that was used ('vector', 'constant', or
        'scalar').
    """
    times = np.asarray(times, dtype=float)
    flat = times.ravel()

    if A_eval_mode == 'vector':
        try:
            At = np.asarray(A(flat))
        except Exception:
            At = None
        if (At is None) or (At.ndim < 3) or (At.shape[0] != flat.shape[0]):
            # The mode hint was wrong for this A: fall back safely
            At = np.array([A(t) for t in flat])
            A_eval_mode = 'scalar'
        return (At.reshape(times.shape + At.shape[-2:])
                .astype(complex, copy=False), A_eval_mode)

    if A_eval_mode == 'constant':
        A0 = np.asarray(A(flat[0]))
        At = np.broadcast_to(A0, flat.shape + A0.shape)
        return (At.reshape(times.shape + A0.shape)
                .astype(complex, copy=False), A_eval_mode)

    if A_eval_mode == 'scalar':
        At = np.array([A(t) for t in flat])
        return (At.reshape(times.shape + At.shape[-2:])
                .astype(complex, copy=False), A_eval_mode)

    # No mode given: probe
    A0 = np.asarray(A(flat[0]))
    target_shape = flat.shape + A0.shape

    At = None
    mode = 'scalar'
    try:
        cand = np.asarray(A(flat))
    except Exception:
        cand = None

    if cand is not None:
        if cand.shape == target_shape:
            # Guard against silent mis-broadcasting: spot-check one point
            k = len(flat) // 2
            spot = np.asarray(A(flat[k]))
            if np.allclose(cand[k], spot, rtol=1.e-10, atol=0.0):
                At = cand
                mode = 'vector'
        elif cand.shape == A0.shape:
            # A returned a single matrix for an array argument: constant A
            k = len(flat) // 2
            spot = np.asarray(A(flat[k]))
            if np.allclose(cand, spot, rtol=1.e-10, atol=0.0):
                At = np.broadcast_to(cand, target_shape)
                mode = 'constant'

    if At is None:  # Fall back to the (slow but safe) per-point loop
        At = np.array([A(t) for t in flat])

    At = At.reshape(times.shape + A0.shape).astype(complex, copy=False)
    return At, mode


def probe_eval_mode(A: Callable, t0: float, t1: float,
                    n_probe: Optional[int] = 5) -> str:
    r"""Determine how the matrix function A can be evaluated.

    Returns 'vector' if A accepts an array of times (fast path),
    'constant' if A ignores its argument, and 'scalar' otherwise.  Use
    the result as the ``A_eval_mode`` argument of
    :func:`magnus_expansion` and :func:`magnus_expansion_multislab` to
    avoid re-probing A on every call.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    A : Callable
        Matrix function of time; see :func:`magnus_expansion`.
    t0, t1 : float
        Interval over which A is probed (t1 >= t0).
    n_probe : int, optional
        Number of sample times used for the probe. Default: 5.

    Returns
    -------
    str
        'vector', 'constant', or 'scalar'.
    """
    times = np.linspace(t0, t1, n_probe)
    _, mode = _evaluate_A(A, times, None)
    return mode


def suggest_n_slabs(
    A: Callable,
    t0: float,
    t1: float,
    A_eval_mode: Optional[str] = None,
    n_probe: Optional[int] = 17,
    phase_per_slab: Optional[float] = 2.0*np.pi
) -> int:
    r"""Suggest a starting number of time slabs for [t0, t1].

    Estimates the accumulated phase ||Omega_1||_2 over the whole
    interval from a coarse sample of A (with the trace removed, since a
    global phase does not affect the probabilities) and suggests enough
    slabs to keep roughly ``phase_per_slab`` (radians) of phase per
    slab.  Starting an adaptive refinement from this estimate skips
    most of the geometric ladder that would otherwise climb from a
    single slab.

    The default of 2 pi radians per slab is deliberately *looser* than
    the Magnus convergence guarantee (pi): empirically, for smooth
    profiles, order-4 methods reach ~1e-3 accuracy already at this slab
    width, and the adaptive refinement loop -- which remains the sole
    arbiter of accuracy -- grows the slab count from here when the
    requested tolerance demands it.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    A : Callable
        Matrix function of time; see :func:`magnus_expansion`.
    t0, t1 : float
        Interval over which the phase is estimated (t1 >= t0).
    A_eval_mode : str, optional
        Skip probing how A can be evaluated; see :func:`probe_eval_mode`.
    n_probe : int, optional
        Number of sample points used to estimate the accumulated phase. Default: 17.
    phase_per_slab : float, optional
        Target accumulated phase per slab, in radians. Default: 2*pi.

    Returns
    -------
    int
        Suggested starting number of slabs (at least 1).
    """
    if not (t1 > t0):
        return 1
    times = np.linspace(t0, t1, n_probe)
    At, _ = _evaluate_A(A, times, A_eval_mode)
    M = (float(t1) - float(t0))*_full_integral(At, 1.0/(n_probe - 1),
                                               'trapezoid')
    dim = M.shape[-1]
    M = M - (np.trace(M)/dim)*np.eye(dim)
    try:
        nrm = np.max(np.linalg.svd(M, compute_uv=False))
    except np.linalg.LinAlgError:
        return 1
    return int(max(1, np.ceil(nrm/phase_per_slab)))


def _cumulative_integral(y: np.ndarray, ds: float, method: str) -> np.ndarray:
    r"""Cumulative integral of y along axis -3 on a uniform grid.

    The result has the same shape as y, with the first entry equal to
    zero, matching scipy's ``cumulative_trapezoid(..., initial=0)``.
    Handles complex integrands.  For 'simpson', splits the integrand
    into real and imaginary parts because scipy's
    ``cumulative_simpson`` silently discards the imaginary part.

    Parameters
    ----------
    y : np.ndarray
        Integrand samples, shape (..., m, d, d), on a uniform grid along axis -3.
    ds : float
        Grid spacing.
    method : str
        'trapezoid' or 'simpson'.

    Returns
    -------
    np.ndarray
        Cumulative integral, same shape as y.
    """
    m = y.shape[-3]
    if method == 'simpson' and m >= 3 and _HAS_CUMULATIVE_SIMPSON:
        re = sp.integrate.cumulative_simpson(y.real, dx=ds, axis=-3, initial=0)
        im = sp.integrate.cumulative_simpson(y.imag, dx=ds, axis=-3, initial=0)
        return re + 1j*im
    # Trapezoid (also the fallback for simpson with very few points or
    # old scipy versions without cumulative_simpson)
    c = np.cumsum(0.5*(y[..., 1:, :, :] + y[..., :-1, :, :]), axis=-3)*ds
    return np.concatenate([np.zeros_like(y[..., :1, :, :]), c], axis=-3)


def _full_integral(y: np.ndarray, ds: float, method: str) -> np.ndarray:
    r"""Integral of y along axis -3 over the full grid (endpoint only).

    Cheaper than :func:`_cumulative_integral` when only the total
    integral is needed (i.e., for the highest requested Magnus order).

    Parameters
    ----------
    y : np.ndarray
        Integrand samples, shape (..., m, d, d), on a uniform grid along axis -3.
    ds : float
        Grid spacing.
    method : str
        'trapezoid' or 'simpson'.

    Returns
    -------
    np.ndarray
        Integral over the full grid, shape (..., d, d).
    """
    m = y.shape[-3]
    if method == 'simpson' and m >= 3:
        return sp.integrate.simpson(y, dx=ds, axis=-3)
    return (np.sum(y, axis=-3) - 0.5*(y[..., 0, :, :] + y[..., -1, :, :]))*ds


def _magnus_terms_quadrature(
    Bt: np.ndarray,
    order: int,
    integration_method: str
) -> np.ndarray:
    r"""Magnus terms Omega_1..Omega_order from samples of A.

    Parameters
    ----------
    Bt : np.ndarray
        Samples of the rescaled matrix function, shape (..., m, d, d):
        Bt = width * A(t(s)) on the uniform normalized grid s in [0, 1]
        with m points, so that all integrals run over [0, 1].  Any
        leading axes (e.g., a slab axis) broadcast through.
    order : int
        Highest Magnus order to compute (1 <= order <= 6).
    integration_method : str
        'trapezoid' or 'simpson'.

    Returns
    -------
    np.ndarray
        Stacked terms, shape (order, ..., d, d).

    Notes
    -----
    Implements the Bernoulli-number recursion (see module docstring).
    The commutators C_k = [Omega_k(s), A(s)] and the nested combinations
    that repeat across orders are computed once and reused.  For the
    highest requested order only the endpoint integral is computed.
    """
    m = Bt.shape[-3]
    ds = 1.0/(m - 1)

    def integ(y: np.ndarray, k: int) -> np.ndarray:
        # Cumulative integral if Omega_k(s) is needed by higher orders;
        # plain endpoint integral for the highest requested order.
        if k < order:
            return _cumulative_integral(y, ds, integration_method)
        return _full_integral(y, ds, integration_method)

    def last(ot: np.ndarray, k: int) -> np.ndarray:
        return ot[..., -1, :, :] if k < order else ot

    terms = []

    o1t = integ(Bt, 1)
    terms.append(last(o1t, 1))

    if order >= 2:
        C1 = commutator(o1t, Bt)                      # [Omega_1, A]
        o2t = integ(-0.5*C1, 2)
        terms.append(last(o2t, 2))

    if order >= 3:
        C2 = commutator(o2t, Bt)                      # [Omega_2, A]
        D11 = commutator(o1t, C1)                     # [Omega_1, [Omega_1, A]]
        o3t = integ(-0.5*C2 + F1*D11, 3)
        terms.append(last(o3t, 3))

    if order >= 4:
        C3 = commutator(o3t, Bt)                      # [Omega_3, A]
        D12 = commutator(o1t, C2)                     # [Omega_1, [Omega_2, A]]
        D21 = commutator(o2t, C1)                     # [Omega_2, [Omega_1, A]]
        o4t = integ(-0.5*C3 + F1*(D12 + D21), 4)
        terms.append(last(o4t, 4))

    if order >= 5:
        C4 = commutator(o4t, Bt)                      # [Omega_4, A]
        o5t = integ(
            -0.5*C4
            + F1*(commutator(o1t, C3) + commutator(o2t, C2)
                  + commutator(o3t, C1))
            + F2*commutator(o1t, commutator(o1t, D11)),
            5)
        terms.append(last(o5t, 5))

    if order >= 6:
        C5 = commutator(o5t, Bt)                      # [Omega_5, A]
        o6t = integ(
            -0.5*C5
            + F1*(commutator(o1t, C4) + commutator(o2t, C3)
                  + commutator(o3t, C2) + commutator(o4t, C1))
            + F2*(commutator(o1t, commutator(o1t, D12))
                  + commutator(o1t, commutator(o1t, D21))
                  + commutator(o1t, commutator(o2t, D11))
                  + commutator(o2t, commutator(o1t, D11))),
            6)
        terms.append(last(o6t, 6))

    return np.stack(terms, axis=0)


def _magnus_gl(
    An: np.ndarray,
    widths: Union[float, np.ndarray],
    order: int
) -> np.ndarray:
    r"""Magnus operator Omega from Gauss-Legendre collocation.

    Commutator-free Magnus integrators of order 2, 4, and 6 based on
    Gauss-Legendre nodes (Blanes, Casas & Ros 2000; Blanes et al. 2009,
    Sec. 5.4).  Exact quadrature order matched to the truncation order,
    using only 1, 2, or 3 evaluations of A per slab.

    Parameters
    ----------
    An : np.ndarray
        A evaluated at the GL nodes, shape (..., n_nodes, d, d).
    widths : float or np.ndarray
        Slab widths h, broadcastable against the leading axes of An.
    order : int
        Requested order; mapped to the smallest GL scheme with at least
        that order (1-2 -> GL1, 3-4 -> GL2, 5-6 -> GL3).

    Returns
    -------
    np.ndarray
        The total Magnus operator Omega, shape (..., d, d).
    """
    h = np.asarray(widths)[..., None, None]

    if order <= 2:
        # Midpoint rule: Omega = h A(t0 + h/2)
        return h*An[..., 0, :, :]

    if order <= 4:
        # Omega = (h/2)(A1 + A2) + (sqrt(3)/12) h^2 [A2, A1]
        A1 = An[..., 0, :, :]
        A2 = An[..., 1, :, :]
        return 0.5*h*(A1 + A2) + (np.sqrt(3.0)/12.0)*h*h*commutator(A2, A1)

    # Order 6 (Blanes, Casas & Ros 2000):
    A1 = An[..., 0, :, :]
    A2 = An[..., 1, :, :]
    A3 = An[..., 2, :, :]
    a1 = h*A2
    a2 = (np.sqrt(15.0)/3.0)*h*(A3 - A1)
    a3 = (10.0/3.0)*h*(A3 - 2.0*A2 + A1)
    C1 = commutator(a1, a2)
    C2 = (-1.0/60.0)*commutator(a1, 2.0*a3 + C1)
    return a1 + a3/12.0 + (1.0/240.0)*commutator(-20.0*a1 - a3 + C1, a2 + C2)


def _gl_nodes(order: int) -> np.ndarray:
    r"""Returns the Gauss-Legendre nodes on [0, 1] for the given Magnus order.

    Parameters
    ----------
    order : int
        Requested Magnus order; mapped to the smallest GL scheme with at least that order
        (1-2 -> 1 node, 3-4 -> 2 nodes, 5-6 -> 3 nodes).

    Returns
    -------
    np.ndarray
        GL nodes on [0, 1] (1, 2, or 3 of them).
    """
    if order <= 2:
        return _GL1_NODES
    if order <= 4:
        return _GL2_NODES
    return _GL3_NODES


def _warn_slab_norm(nmax: float):
    r"""Warn if the slab norm proxy nmax = max ||Omega||_2 is >= pi (see
    :class:`MagnusConvergenceWarning`).

    Parameters
    ----------
    nmax : float
        Largest ||Omega||_2 (or a proxy for it) encountered across the slab(s) just evaluated.

    Returns
    -------
    None
    """
    if nmax >= np.pi:
        # The message is intentionally static (no numbers) so that Python's
        # default warning filter shows it only once per session.
        warnings.warn(
            "at least one time slab is too wide for guaranteed convergence "
            "of the Magnus series (||Omega||_2 >= pi); raising the "
            "expansion order will not help there -- more (narrower) slabs "
            "are needed. If a target tolerance (rtol/atol) was requested, "
            "the adaptive refinement narrows the slabs automatically and "
            "this warning can be ignored. Shown once per session.",
            MagnusConvergenceWarning, stacklevel=4)


def _expm_stack(Om: np.ndarray, warn_wide: bool = False,
                A_is_const: bool = False) -> np.ndarray:
    r"""Matrix exponential of one matrix or a stack of matrices.

    If Om is anti-Hermitian (as is always the case for A = -i H with a
    Hermitian Hamiltonian H), the exponential is computed from the
    eigendecomposition of the Hermitian matrix K = i Om:
    exp(Om) = V diag(exp(-i lambda)) V^dagger.  This is faster than
    scipy's Pade-based expm for stacks of small matrices, and it yields
    an exactly unitary result (probabilities that sum to 1 by
    construction).  Otherwise it falls back to scipy.linalg.expm.

    If ``warn_wide`` is True, the eigenvalues (whose maximum modulus is
    ||Om||_2) are also used to warn about slabs too wide for the Magnus
    series to converge; for a constant A (``A_is_const``) the series
    terminates exactly and the check is skipped.

    Parameters
    ----------
    Om : np.ndarray
        Matrix (or stack of matrices), shape (..., d, d).
    warn_wide : bool, optional
        If True, check the slab norm and emit :class:`MagnusConvergenceWarning` if it is too
        large. Default: False.
    A_is_const : bool, optional
        If True, A is constant in time/position, so the Magnus series terminates exactly and the
        convergence check is skipped even if ``warn_wide`` is True. Default: False.

    Returns
    -------
    np.ndarray
        exp(Om), same shape as Om.
    """
    Om = np.asarray(Om)
    K = 1j*Om
    Kh = np.conj(np.swapaxes(K, -1, -2))
    scale = np.max(np.abs(K))
    if scale == 0.0:
        return np.broadcast_to(np.eye(Om.shape[-1], dtype=complex),
                               Om.shape).copy()
    if np.max(np.abs(K - Kh)) <= 1.e-12*scale:
        lam, V = np.linalg.eigh(0.5*(K + Kh))
        if warn_wide and not A_is_const:
            _warn_slab_norm(np.max(np.abs(lam)))  # ||Om||_2 = max |lambda|
        Vh = np.conj(np.swapaxes(V, -1, -2))
        return (V*np.exp(-1j*lam)[..., None, :]) @ Vh
    # General (non-anti-Hermitian) fallback
    if warn_wide and not A_is_const:
        try:
            _warn_slab_norm(np.max(np.linalg.svd(Om, compute_uv=False)))
        except np.linalg.LinAlgError:
            pass
    try:
        return np.asarray(sp.linalg.expm(Om))
    except Exception:
        # Very old scipy without stacked-input support
        shape = Om.shape
        flat = Om.reshape((-1,) + shape[-2:])
        out = np.array([sp.linalg.expm(w) for w in flat])
        return out.reshape(shape)


def _validate(order: int, integration_method: str):
    r"""Validates ``order`` and ``integration_method``.

    Parameters
    ----------
    order : int
        Requested Magnus order; must satisfy 1 <= order <= MAGNUS_EXP_ORDER_MAX.
    integration_method : str
        Must be one of ``valid_integration_methods``.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If ``order`` or ``integration_method`` is invalid.
    """
    if not (integration_method in valid_integration_methods):
        raise ValueError(
            "magnus.magnus_expansion: integration_method must be one of "
            + str(valid_integration_methods) + ", not '"
            + str(integration_method) + "'.")
    if not (1 <= order <= MAGNUS_EXP_ORDER_MAX):
        raise ValueError(
            "magnus.magnus_expansion: order must be between 1 and "
            + str(MAGNUS_EXP_ORDER_MAX) + ", not " + str(order) + ".")


def magnus_expansion(
    A: Callable,
    t0: float,
    t1: float,
    n_tpts: Optional[int] = 50,
    order: Optional[int] = 2,
    integration_method: Optional[str] = 'trapezoid',
    return_magnus_terms: Optional[bool] = False,
    validate_input: Optional[bool] = True,
    A_eval_mode: Optional[str] = None
) -> np.ndarray:
    r"""Compute exp(Omega_1 + ... + Omega_order) of A(t) from t0 to t1.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    A : Callable
        Matrix function of time; must return a (d, d) NumPy array for a
        scalar time.  If it also accepts an array of times (returning
        a (n, d, d) stack), the vectorized form is used automatically
        for speed; this is detected silently and verified against a
        scalar evaluation.
    t0, t1 : float
        Integration limits (t1 >= t0).
    n_tpts : int, optional
        Number of uniformly spaced time points used to evaluate the
        integrals ('trapezoid'/'simpson' methods only; >= 2).
    order : int, optional
        Highest Magnus order (1 to 6).
    integration_method : str, optional
        'trapezoid', 'simpson', or 'gl' (Gauss-Legendre collocation;
        ignores ``n_tpts`` and uses 1, 2, or 3 nodes for orders <= 2,
        <= 4, <= 6, respectively).
    return_magnus_terms : bool, optional
        If True, also return the individual Magnus terms.  For the
        'gl' method the terms are not separable, and a single-element
        list containing the total Omega is returned instead.
    validate_input : bool, optional
        If True, validate ``order`` and ``integration_method``
        (raises ValueError on invalid input).

    Returns
    -------
    np.ndarray, or (np.ndarray, np.ndarray)
        The evolution operator U = exp(sum_k Omega_k); if
        ``return_magnus_terms`` is True, also the stacked terms.
    """
    if validate_input:
        _validate(order, integration_method)

    if integration_method == 'gl':
        nodes = _gl_nodes(order)
        width = float(t1) - float(t0)
        tnodes = t0 + width*nodes
        An, used_mode = _evaluate_A(A, tnodes, A_eval_mode)
        Om = _magnus_gl(An, width, order)
        U = _expm_stack(Om, warn_wide=True, A_is_const=(used_mode == 'constant'))
        if not return_magnus_terms:
            return U
        return U, np.stack([Om], axis=0)

    times = np.linspace(t0, t1, n_tpts)
    At, used_mode = _evaluate_A(A, times, A_eval_mode)
    Bt = (float(t1) - float(t0))*At  # rescale to the unit interval
    magnus_terms = _magnus_terms_quadrature(Bt, order, integration_method)

    U = _expm_stack(np.sum(magnus_terms, axis=0), warn_wide=True,
                    A_is_const=(used_mode == 'constant'))
    if not return_magnus_terms:
        return U
    return U, magnus_terms


def evolution_operators_from_samples(
    At: np.ndarray,
    widths: Union[list, np.ndarray],
    order: Optional[int] = 2,
    integration_method: Optional[str] = 'trapezoid',
    A_is_const: Optional[bool] = False,
    validate_input: Optional[bool] = True
) -> np.ndarray:
    r"""Evolution operators of a chain of slabs from precomputed samples.

    Mid-level entry point for callers that build the samples of A
    themselves -- e.g., to batch extra axes (such as the neutrino
    energy) in front of the slab axis, which this routine broadcasts
    through all operations.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    At : np.ndarray
        Samples of A, shape (..., n_slabs, m, d, d).  For the
        quadrature methods ('trapezoid'/'simpson'), the m samples of
        each slab lie on the uniform grid spanning the slab (endpoints
        included).  For 'gl', they lie on the Gauss-Legendre nodes
        (m = 1, 2, or 3 for orders <= 2, <= 4, <= 6; see
        :func:`gl_nodes`).
    widths : list or np.ndarray
        Slab widths, shape (n_slabs,) (or broadcastable to the leading
        axes of ``At`` without the last three).
    order : int, optional
        Highest Magnus order (1 to 6).
    integration_method : str, optional
        'trapezoid', 'simpson', or 'gl'.
    A_is_const : bool, optional
        Set to True if A is constant in time to skip the (inapplicable)
        slab-width convergence warning.
    validate_input : bool, optional
        If True, validate order and integration_method.

    Returns
    -------
    np.ndarray
        Evolution operators, shape (..., n_slabs, d, d).
    """
    if validate_input:
        _validate(order, integration_method)
    w = np.asarray(widths, dtype=float)
    if integration_method == 'gl':
        Om = _magnus_gl(At, w, order)
        return _expm_stack(Om, warn_wide=True, A_is_const=A_is_const)
    Bt = w[..., None, None, None]*At        # rescale to the unit interval
    magnus_terms = _magnus_terms_quadrature(Bt, order, integration_method)
    return _expm_stack(np.sum(magnus_terms, axis=0), warn_wide=True,
                       A_is_const=A_is_const)


def gl_nodes(order: int) -> np.ndarray:
    r"""Returns the Gauss-Legendre nodes on [0, 1] used by the 'gl' method.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    order : int
        Requested Magnus order; mapped to the smallest GL scheme with at least that order
        (1-2 -> 1 node, 3-4 -> 2 nodes, 5-6 -> 3 nodes).

    Returns
    -------
    np.ndarray
        GL nodes on [0, 1] (1, 2, or 3 of them).
    """
    return _gl_nodes(order)


def magnus_expansion_multislab(
    A: Callable,
    t_slab_edges: Union[list, np.ndarray],
    n_tpts_per_slab: Optional[int] = 50,
    order: Optional[int] = 2,
    integration_method: Optional[str] = 'trapezoid',
    validate_input: Optional[bool] = True,
    A_eval_mode: Optional[str] = None
) -> np.ndarray:
    r"""Compute the evolution operators of all time slabs at once.

    Vectorized (batched) version of :func:`magnus_expansion` for a
    chain of time slabs: A is evaluated for all slabs in a single call
    (when it supports array input), and the quadrature, commutator
    algebra, and matrix exponentials are evaluated as batched NumPy
    operations with the slab axis leading.  This is much faster than
    calling :func:`magnus_expansion` slab by slab.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    A : Callable
        Matrix function of time (see :func:`magnus_expansion`).
    t_slab_edges : list or np.ndarray
        Slab edges, shape (n_slabs, 2): [[t0, t1], [t1, t2], ...].
        Slabs of zero width yield identity operators.
    n_tpts_per_slab : int, optional
        Number of time points per slab ('trapezoid'/'simpson' only).
    order : int, optional
        Highest Magnus order (1 to 6).
    integration_method : str, optional
        'trapezoid', 'simpson', or 'gl'.
    validate_input : bool, optional
        If True, validate input (raises ValueError on invalid input).

    Returns
    -------
    np.ndarray
        Stack of evolution operators, shape (n_slabs, d, d), ordered
        like ``t_slab_edges`` (i.e., earliest slab first).

    Notes
    -----
    The time-ordered product over the chain is
    ``U_total = U[n_slabs-1] @ ... @ U[1] @ U[0]``,
    i.e., the *last* slab is the leftmost factor.
    """
    edges = np.asarray(t_slab_edges, dtype=float)
    if edges.ndim == 1:
        edges = edges[None, :]
    widths = edges[:, 1] - edges[:, 0]

    if validate_input:
        _validate(order, integration_method)
        if np.any(widths < 0.0):
            raise ValueError(
                "magnus.magnus_expansion_multislab: all slabs must have "
                "t1 >= t0.")

    if integration_method == 'gl':
        s = _gl_nodes(order)                            # (k,) GL nodes
    else:
        s = np.linspace(0.0, 1.0, n_tpts_per_slab)      # normalized grid
    tgrid = edges[:, :1] + widths[:, None]*s            # (n_slabs, m)
    At, used_mode = _evaluate_A(A, tgrid, A_eval_mode)  # (n_slabs, m, d, d)
    return evolution_operators_from_samples(At, widths, order,
        integration_method, A_is_const=(used_mode == 'constant'),
        validate_input=False)


if __name__ == "__main__":

    def A(t):
        return np.array([[1.0*t, 2.0 + 3j*t], [2.0 - 3j*t, 2.0]])

    t0, t1 = 0.0, 1.0
    exp_Omega_1 = magnus_expansion(A, t0, t1, n_tpts=100, order=4,
                                   integration_method='trapezoid')
    print(exp_Omega_1)
    exp_Omega_2 = magnus_expansion(A, t0, t1, n_tpts=100, order=4,
                                   integration_method='simpson')
    print(exp_Omega_2)
    exp_Omega_3 = magnus_expansion(A, t0, t1, order=4,
                                   integration_method='gl')
    print(exp_Omega_3)
    print(exp_Omega_1 - exp_Omega_2)
    print(exp_Omega_1 - exp_Omega_3)


__all__ = [
    'MagnusConvergenceWarning',
    'B',
    'F1',
    'F2',
    'f1',
    'f2',
    'MAGNUS_EXP_ORDER_MAX',
    'valid_integration_methods',
    'commutator',
    'probe_eval_mode',
    'suggest_n_slabs',
    'magnus_expansion',
    'evolution_operators_from_samples',
    'gl_nodes',
    'magnus_expansion_multislab',
]
