# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""magnus.py

Compute the time-evolution operator using the Magnus expansion.

This module contains the numerical core of Magnus: routines to compute
the matrix exponential of the Magnus expansion of a (possibly
time-dependent) matrix function :math:`A(t)`, i.e.,

.. math::

   U(t_1, t_0) = \exp\!\left[\Omega_1 + \Omega_2 + \cdots + \Omega_k\right] ,

where the terms :math:`\Omega_k` are built from time-ordered integrals
of nested commutators of :math:`A(t)`.  For neutrino oscillations,
:math:`A(t) = -i H(t)`, with :math:`H(t)` the Hamiltonian, but the
routines below work for arbitrary matrix-valued :math:`A(t)`.

The terms are generated with the standard recursion based on Bernoulli
numbers [1]_ (in the :math:`B_1 = -1/2` convention):

.. math::

   \Omega_1(t) &= \int_0^t A(s)\, ds \\
   \Omega_n(t) &= \sum_{j=1}^{n-1} \frac{B_j}{j!} \int_0^t S_n^{(j)}(s)\, ds ,

with :math:`S_n^{(j)}` the sums of nested commutators of the lower-order
terms with :math:`A` (:math:`B_3 = B_5 = 0`, so those groups vanish
identically).  Orders 1--6 are written out inline; above that the terms
are generated from the same recursion, at any order.

Two families of methods are available, selected via
``integration_method``:

* ``'gl'`` (the default): Gauss-Legendre collocation
  [1]_ [2]_.  For a slab of width :math:`h` it needs only 1, 2, 3, or 4
  evaluations of :math:`A` to reach order 2, 4, 6, or 8, respectively, with
  quadrature error matched to the truncation order.  ``n_tpts`` is
  ignored.  Both faster and more accurate than the alternatives whenever
  :math:`A(t)` is smooth within each slab, which is the common case --
  and, for the Earth, is what aligning slab edges with the PREM layer
  boundaries is for.

* ``'trapezoid'`` / ``'simpson'``: sample :math:`A(t)` on a uniform grid
  of ``n_tpts`` points and evaluate the nested integrals with cumulative
  quadrature.  Fully general, and so the safer choice if :math:`A(t)`
  has a kink or a discontinuity *inside* a slab, where Gauss-Legendre
  loses its order advantage.  The quadrature error
  (:math:`\mathcal{O}(h^2)` or :math:`\mathcal{O}(h^4)`) can dominate
  the Magnus truncation error at high orders unless ``n_tpts`` grows
  accordingly.

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
    * cached_eval_mode - Context manager reusing one probe_eval_mode result
           for a callable that will be probed more than once
    * ordered_product - Time-ordered product of a stack of slab operators,
           earliest slab first
    * palindromic - Returns whether every array given reads the same both
           ways, the geometric precondition for the half-chord optimization
    * suggest_n_slabs - Suggests a starting number of time slabs
    * magnus_expansion - Computes :math:`\exp(\Omega)` for a single time slab
    * evolution_operators_from_samples - Evolution operators of a chain
           of slabs from precomputed samples of A
    * gl_nodes - Returns the Gauss-Legendre nodes used by the 'gl' method
    * magnus_expansion_multislab - Computes the evolution operators of
           all time slabs at once, from A directly
    * MagnusConvergenceWarning - Warning class for slabs too wide for
           guaranteed Magnus convergence
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import os
import warnings
import weakref
from contextlib import contextmanager
from typing import Optional, Callable, Union, Tuple

import numpy as np
import scipy as sp

from magnus import expmkernels


class MagnusHighOrderCostWarning(UserWarning):
    r"""Warns that a Magnus order above 6 costs substantially more per slab.

    Orders 1-6 are written out inline with their shared subexpressions named and reused.
    Above that the terms are generated from the recursion, and their number roughly doubles
    per order (9 terms at order 6; 17, 33, 65, 129 at orders 7-10), so the work per slab
    grows with it -- measured at roughly 2.7x order 6 at order 7, rising to about 17x at
    order 10, for the same grid.

    Higher order buys a genuinely faster convergence rate in the slab width, so this is a
    trade rather than a mistake.  But it is often the worse side of the trade: narrowing
    the slabs at order 4 or 6 usually reaches a given accuracy for less total work, and
    beyond the Magnus series' convergence radius no order helps at all (see
    :class:`MagnusConvergenceWarning`).

    .. versionadded:: 1.0.0
    """


class ScalarHamiltonianWarning(UserWarning):
    r"""Warns that ``H_func`` accepts only one position at a time.

    The engine evaluates the Hamiltonian at every quadrature node of every slab
    -- often a few hundred positions for a single probability, and the adaptive
    refinement repeats that at each level. ``_evaluate_A`` therefore tries a
    single vectorized call, ``A(times)``, and uses the result if it has the
    right shape and agrees with a scalar spot-check. If that fails it falls back
    to a Python loop, one call per position.

    That fallback is correct but typically several times slower, and it is
    *silent*: nothing about a scalar-only ``H_func`` looks wrong, so the slow
    path is easy to sit on indefinitely. Measured on a three-flavor
    exponential-density profile, making the same ``H_func`` array-capable cut
    the time per :func:`magnus.oscprob.osc_prob` call from 7.8 ms to 1.7 ms,
    a factor of 4.6, with bit-identical output.

    Making a Hamiltonian array-capable usually means no more than writing its
    position dependence with NumPy and letting the matrix part broadcast::

        # slow: one position at a time
        def H_func(l):
            VCC = matter.VCC_func(l, num_density_e_func)
            return (1.0/energy)*h_vac + hamiltonians.hamiltonian_3nu_matter(VCC)

        # fast: the same physics, evaluated for all positions at once
        e00 = np.diag([1.0, 0.0, 0.0])
        def H_func(l):
            l = np.asarray(l, dtype=float)
            VCC = vcc_of(l)                      # returns an array
            return (1.0/energy)*h_vac + VCC[..., None, None]*e00

    The trailing ``[..., None, None]`` is what lets one potential per position
    multiply a stack of matrices. A Hamiltonian that ignores its argument
    entirely is detected separately and costs nothing, so constant-density cases
    never trigger this.

    Pass ``A_eval_mode='scalar'`` to :func:`magnus_expansion` (or accept the
    warning) when a scalar-only Hamiltonian is genuinely unavoidable.

    .. versionadded:: 1.0.0
    """


class MagnusConvergenceWarning(UserWarning):
    r"""Warns that a time slab may be too wide for the Magnus series.

    **What was detected.**  The Magnus series is guaranteed to converge when
    :math:`\int_{t_0}^{t_1} \lVert A(t)\rVert_2\, dt < \pi`.  :math:`\lVert\Omega\rVert_2 \geq
    \pi` is used as a cheap proxy for that integral -- it comes free from the eigenvalues already
    computed for the matrix exponential -- so this fires when a *sufficient* condition for
    convergence was not met on at least one slab.  The message says how far past :math:`\pi`, in
    three buckets, which is the one quantity this check actually knows.

    **What it means for the answer: unknown, and that is the honest answer.**  This is a
    statement about the slab width, not about the error.  The condition is sufficient, not
    necessary, so exceeding it does not imply a wrong answer -- and it fires on results accurate
    to 1.6e-06 (``docs/dev/DECISION_DISPATCH_ORDER.md`` §5) as well as on results seven times
    outside a requested 1e-3.  Anything that claims to tell you which of those you have is
    claiming more than this check can support; :class:`magnus.oscprob.ToleranceNotAchievedWarning` is the one
    that reports a failed convergence *test*.

    **What to change.**  More, narrower slabs: request a smaller ``rtol``/``atol``, or raise
    ``n_slabs``.  Raising ``magnus_exp_order`` does **not** help in this regime -- beyond the
    series' radius no order converges.  If the profile has a density jump or a kink, pass
    ``t_breakpoints`` there as well: a slab straddling one is never fixed by more slabs, only
    narrowed.

    **When it is safe to ignore.**  When the answer has been checked another way -- a tighter
    tolerance giving the same result, or :func:`magnus.oscprob.cross_check_strategies` showing
    a different engine agreeing.  **Not** merely because a tolerance was requested.  That advice
    used to be in this message and it is false in exactly the cases where the warning matters:
    measured on a sawtooth density with ``rtol=atol=1e-3`` explicitly requested, under both
    ``strategy='auto'`` and ``strategy='magnus'``, the adaptive refinement ran and the answer was
    still **7.484e-03**, seven times outside the tolerance asked for, with this warning showing.

    **Measured rates** (``docs/dev/adversarial_batteries/warn_fp.py``, 168 configurations across
    the profile families this package serves, d = 2-5, scored against ``solve_ivp`` or, for
    piecewise profiles, against ``expm``): fired 70 times, of which **17 true positives and 53
    false positives -- a 76 % false-positive rate**, the highest of any warning here.  That is
    the price of reporting a *sufficient* condition, and it is why the text above refuses to
    translate the condition into a claim about the error.

    **Where that noise comes from, and what would fix it.**  Of 66 single-point calls, some
    refinement level exceeded :math:`\pi` in 46 -- but the level whose answer was actually
    returned did so in only **7**.  So **39 of 46 firings, 85 %, describe an intermediate grid
    that nobody receives**: the ladder started coarse, said so, then refined and never retracted
    it.  Keying the warning to the returned level alone would cut false alarms from 31 to 5 at a
    similar rate (67 % against 71 %).  That change is *mechanical* -- capture the norm per level
    and emit once the loop has decided -- and is **deliberately not made here**, because it
    touches the refinement loop and the warning plumbing several tests depend on.  It is written
    down with its numbers so it can be made deliberately rather than rediscovered.

    .. versionadded:: 1.0.0
    """


# Bernoulli numbers B_k (negative-B_1 convention), kept for reference;
# only B_1, B_2, and B_4 enter at the orders implemented here (<= 6).
B = {
    0: 1.0, 1: -0.5, 2: 1.0/6.0, 3: 0.0, 4: -1.0/30.0, 5: 0.0, 6: 1.0/42.0,
}

# Multiplicative factors B_j/j! of the commutator groups in the recursion.  Only j = 1 and
# the even j contribute, since B_j = 0 for every odd j >= 3.  These are the numbers
# magnus.expansionterms.bernoulli_factor() derives from the Bernoulli recursion in exact
# rational arithmetic; tests/test_expansionterms.py checks them against it.
F1 = 1.0 / 12.0          # B_2 / 2!
F2 = -1.0 / 720.0        # B_4 / 4!
F3 = 1.0 / 30240.0       # B_6 / 6!   (first needed at order 7)
F4 = -1.0 / 1209600.0    # B_8 / 8!   (first needed at order 9)

# Backward-compatible aliases
f1 = F1
f2 = F2

# Coefficient of each commutator group, keyed by the group index j.  Used by the
# composition-driven path for orders 7 and above; orders 1-6 spell these out inline.
_GROUP_FACTORS = {1: -0.5, 2: F1, 4: F2, 6: F3, 8: F4}

# Highest order of the Magnus expansion implemented here.  Re-exported by globaldefs so
# there is one definition rather than two that have to be kept in step by hand.
MAGNUS_EXP_ORDER_MAX = 10

# Highest order we implement on Gauss-Legendre nodes: one, two, three and four nodes
# give orders two, four, six and eight.  These are separately derived integrators, not
# products of the Magnus recursion, so they do not extend along with it -- each order is
# its own construction.  Orders two to six follow Blanes, Casas & Ros, BIT 40 (2000) 434;
# orders six and eight use the commutator-optimal forms of Blanes, Casas & Ros, BIT 42
# (2002) 262, which need three and six commutators, the fewest possible at each order.
MAGNUS_EXP_ORDER_MAX_GL = 8

# Valid values of integration_method
valid_integration_methods = ['gl', 'trapezoid', 'simpson']

# Gauss-Legendre nodes on [0, 1] used by the 'gl' method
_GL1_NODES = np.array([0.5])
_GL2_NODES = np.array([0.5 - np.sqrt(3.0)/6.0, 0.5 + np.sqrt(3.0)/6.0])
_GL3_NODES = np.array([0.5 - np.sqrt(15.0)/10.0, 0.5,
                       0.5 + np.sqrt(15.0)/10.0])
# Four-node Gauss-Legendre, for the order-8 scheme.  The offsets from the slab midpoint
# and the matching weights are kept as named constants because the order-8 expression
# needs the weights themselves, not only the node positions.
_GL4_V1 = 0.5*np.sqrt((3.0 + 2.0*np.sqrt(6.0/5.0))/7.0)
_GL4_V2 = 0.5*np.sqrt((3.0 - 2.0*np.sqrt(6.0/5.0))/7.0)
_GL4_W1 = 0.5 - np.sqrt(5.0/6.0)/6.0
_GL4_W2 = 0.5 + np.sqrt(5.0/6.0)/6.0
_GL4_NODES = np.array([0.5 - _GL4_V1, 0.5 - _GL4_V2,
                       0.5 + _GL4_V2, 0.5 + _GL4_V1])

_HAS_CUMULATIVE_SIMPSON = hasattr(sp.integrate, 'cumulative_simpson')


def commutator(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    r"""Returns the commutator [X, Y] = X Y - Y X.

    Works on single matrices and on stacks of matrices (the matrix
    product broadcasts over all leading axes).

    .. versionadded:: 1.0.0

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
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import magnus

        X = np.array([[0.0, 1.0], [0.0, 0.0]])
        Y = np.array([[0.0, 0.0], [1.0, 0.0]])

        print(magnus.commutator(X, Y))
        print('antisymmetric:',
              np.array_equal(magnus.commutator(X, Y), -magnus.commutator(Y, X)))
"""
    return X @ Y - Y @ X


def _commutator_batched_core(X, Y):  # pragma: no cover -- compiled below
    r"""Commutator [X, Y] of two equal-shaped matrix stacks, one matrix at a time.

    For ``X`` and ``Y`` of shape ``(nB, d, d)`` returns the ``(nB, d, d)`` stack
    of ``X[b] @ Y[b] - Y[b] @ X[b]``, both products fused into a single loop
    nest so each element is one accumulated scalar sum.  Written to be compiled
    by numba; the pure-Python form exists only as compilation input.
    """
    nB, d, _ = X.shape
    out = np.empty((nB, d, d), dtype=np.complex128)
    for b in range(nB):
        for i in range(d):
            for j in range(d):
                s = 0.0 + 0.0j
                for m in range(d):
                    s += X[b, i, m]*Y[b, m, j] - Y[b, i, m]*X[b, m, j]
                out[b, i, j] = s
    return out


if expmkernels.HAVE_NUMBA:
    _commutator_batched_kernel = expmkernels._jit(_commutator_batched_core)
else:                                                   # pragma: no cover
    _commutator_batched_kernel = None


def _commutator_batched(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    r"""The commutator [X, Y] for the Magnus hot paths, compiled when possible.

    Same mathematics as :func:`commutator`, which stays the public, pure-NumPy
    form.  On the stacks of small matrices the engine works with, the two
    batched matmuls of ``X @ Y - Y @ X`` cost mostly gufunc dispatch -- about
    185 ns per 3x3 matrix -- so the Magnus term recursion and the
    Gauss-Legendre schemes route their commutators here instead: with numba
    present and both operands complex128 stacks of one shape, a compiled
    kernel fuses the two products into one pass.  Anything else -- no numba,
    another dtype, operands that would need broadcasting -- falls through to
    the same expression :func:`commutator` computes, so a numba-less install
    is bit-identical to what these call sites always produced.

    With numba the two installs are no longer bit-identical to each other on
    these paths: the kernel accumulates each matrix element as one interleaved
    scalar sum where NumPy rounds the two products separately, so probabilities
    can move at the rounding level.  Worst observed shift 6.7e-14 across 36
    scan configurations (both benchmark profiles, d = 2-5, gl orders 4-8,
    simpson, ladders from rtol 1e-3 to 1e-11), every refinement decision and
    warning unchanged.  At three flavors and below the shift measures exactly
    zero at fixed slab counts: there the commutator enters :math:`\Omega`
    suppressed by the squared slab width, and the rounding difference falls
    below the ulp of :math:`\Omega`'s leading term.

    .. versionadded:: 1.0.7

    Parameters
    ----------
    X : np.ndarray
        Left stack of matrices, shape (..., d, d).
    Y : np.ndarray
        Right stack of matrices, same shape as X (the hot call sites never
        broadcast; a pair that would is handed to the NumPy expression).

    Returns
    -------
    np.ndarray
        The commutator X @ Y - Y @ X, shaped like X.
    """
    if ((_commutator_batched_kernel is not None) and (X.shape == Y.shape)
            and (X.dtype == np.complex128) and (Y.dtype == np.complex128)):
        d = X.shape[-1]
        return _commutator_batched_kernel(
            np.ascontiguousarray(X).reshape((-1, d, d)),
            np.ascontiguousarray(Y).reshape((-1, d, d))).reshape(X.shape)
    return X @ Y - Y @ X


def _warn_scalar_hamiltonian() -> None:
    r"""Warn that the Hamiltonian is being evaluated one position at a time.

    Raised where the vectorization probe fails, which is the only place the
    engine learns that ``H_func`` cannot take an array. See
    :class:`ScalarHamiltonianWarning` for why this matters and how to fix it.

    .. versionadded:: 1.0.0
    """
    warnings.warn(
        "magnus: the Hamiltonian could not be evaluated for several positions "
        "at once, so it is being called one position at a time. This is "
        "correct but slower -- measured 4.6x on a 3nu exponential-density "
        "profile -- because the engine samples the Hamiltonian at every "
        "quadrature node of every slab, and the adaptive refinement repeats "
        "that at each level. To take the fast path, write H_func so that it "
        "accepts an array of positions and returns a stack of matrices: turn "
        "the position dependence into NumPy operations and broadcast the "
        "matrix part, e.g. 'VCC[..., None, None]*e00' instead of "
        "'VCC*e00'. A Hamiltonian that ignores its argument is detected "
        "separately and never triggers this. Shown once per session.",
        ScalarHamiltonianWarning, stacklevel=3)


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
            _warn_scalar_hamiltonian()
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
        _warn_scalar_hamiltonian()
        At = np.array([A(t) for t in flat])

    At = At.reshape(times.shape + A0.shape).astype(complex, copy=False)
    return At, mode


# Whether a given H_func *accepts an array of positions* is a property of the
# function: a callable that does today will tomorrow.  Probing costs three calls
# into the user's Hamiltonian, cheap for a PREM lookup and emphatically not for
# an interpolated profile or a quadrature (notebook 19's long-range potential is
# about a third of its own call), so the answer is worth remembering.
#
# **But the interval is part of the question, not context for it.**  'constant'
# means "sampling A across [t0, t1] gave the same matrix every time", which a
# wider interval can falsify -- a two-layer profile that short-circuits when all
# requested positions fall in one layer probes as 'constant' on a short baseline
# and 'vector' on a long one, for the same function object.  Keyed on the
# function alone, a mode learned on the short one was served for the long one and
# `_evaluate_A` then broadcast a single sample over a profile that varies: a
# unitary, unwarned answer wrong by 5.8e-02, and wrong only for one order of the
# caller's loop.  The interval is in the key for that reason.
#
# Keyed weakly, so a closure rebuilt per call simply misses and re-probes, and
# nothing is kept alive that the caller has dropped.  Note that
# `WeakKeyDictionary` keys by the referent's *equality*, not by identity -- an
# earlier version of this comment claimed the opposite, which matters for a
# caller whose Hamiltonian defines `__eq__`: two objects comparing equal share an
# entry.  That is a documented property of the container, not a choice made here.
_EVAL_MODE_CACHE = weakref.WeakKeyDictionary()

_EVAL_MODE_CACHE_MAX = 256
r"""int: How many distinct intervals are remembered per Hamiltonian before the lot is dropped.

The weak keying bounds the *outer* dictionary but says nothing about the inner one: a
Hamiltonian defined at module scope never dies, so without a ceiling its span dict grows for
the lifetime of the process.  A direct ``osc_prob`` loop over distinct baselines is exactly
that shape -- the interval is what varies per point -- and 1000 baselines retained 1000
entries, about 184 KB.

Cleared wholesale rather than LRU-evicted, matching
``hamiltonians3nu._VACUUM_H_CACHE`` and ``matter._VCC_CONST_CACHE``.  The case this cache
exists for is the refinement ladder, which calls repeatedly at *one* interval and so holds a
single entry that no eviction can reach; the case that fills it is a scan, where each entry
is used once and evicting the wrong one costs nothing.  Neither population rewards a smarter
policy.

.. versionadded:: 1.0.0
"""


def cached_eval_mode(A: Callable, t0: float, t1: float, key=None) -> str:
    r"""``probe_eval_mode`` for a callable that will be probed more than once.

    Returns the same value :func:`probe_eval_mode` would for *this interval*, and
    remembers it against ``(key or A, t0, t1)`` so a repeated call on the same
    Hamiltonian over the same span does not evaluate it three more times.  Falls
    straight through for anything that cannot be weakly referenced or hashed.

    ``key`` exists because callers often have to wrap the object they want cached:
    ``probe_eval_mode`` needs :math:`A = -iH`, and a fresh ``lambda t: -1j*H(t)``
    per call would miss every time.  Passing ``key=H_func`` caches against the
    thing whose signature is actually being described.  Multiplying by a constant
    cannot change whether a function accepts an array, so the two share a verdict.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    A : Callable
        The matrix function to probe.
    t0, t1 : float
        The interval to probe over.  Part of the cache key: see the comment above
        ``_EVAL_MODE_CACHE`` for the wrong answer that omitting it produced.
    key : optional
        Object to cache against instead of ``A``.  Defaults to ``A``.

    Returns
    -------
    str
        'vector', 'scalar' or 'constant'; see :func:`probe_eval_mode`.
    """
    holder = A if key is None else key
    span = (float(t0), float(t1))
    try:
        by_span = _EVAL_MODE_CACHE.get(holder)
    except TypeError:
        # Not weak-referenceable (a __slots__ class) or not hashable (a dataclass,
        # which sets __hash__ = None, or anything defining __eq__).  All of those
        # are ordinary ways to write a Hamiltonian, so probe and move on rather
        # than letting the cache decide whether the call is allowed to succeed.
        return probe_eval_mode(A, t0, t1)
    if by_span is not None:
        hit = by_span.get(span)
        if hit is not None:
            return hit
    mode = probe_eval_mode(A, t0, t1)
    try:
        if by_span is None:
            _EVAL_MODE_CACHE[holder] = {span: mode}
        else:
            if len(by_span) >= _EVAL_MODE_CACHE_MAX:
                by_span.clear()
            by_span[span] = mode
    except TypeError:
        pass
    return mode


def probe_eval_mode(A: Callable, t0: float, t1: float,
                    n_probe: Optional[int] = 5) -> str:
    r"""Determine how the matrix function A can be evaluated.

    Returns 'vector' if A accepts an array of times (fast path),
    'constant' if A ignores its argument, and 'scalar' otherwise.  Use
    the result as the ``A_eval_mode`` argument of
    :func:`magnus_expansion` and :func:`magnus_expansion_multislab` to
    avoid re-probing A on every call.

    .. versionadded:: 1.0.0

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
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import magnus

        def vectorized(t):
            return -1j*np.eye(2)*np.asarray(t)[..., None, None]

        print('array-capable :', magnus.probe_eval_mode(vectorized, 0.0, 1.0))
        print('constant      :', magnus.probe_eval_mode(lambda t: -1j*np.eye(2),
                                                        0.0, 1.0))
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

    Estimates the accumulated phase :math:`\lVert\Omega_1\rVert_2` over
    the whole interval from a coarse sample of A (with the trace removed, since a
    global phase does not affect the probabilities) and suggests enough
    slabs to keep roughly ``phase_per_slab`` (radians) of phase per
    slab.  Starting an adaptive refinement from this estimate skips
    most of the geometric ladder that would otherwise climb from a
    single slab.

    The default of :math:`2\pi` radians per slab is deliberately
    *looser* than the Magnus convergence guarantee (:math:`\pi`): empirically, for smooth
    profiles, order-4 methods reach ~1e-3 accuracy already at this slab
    width, and the adaptive refinement loop -- which remains the sole
    arbiter of accuracy -- grows the slab count from here when the
    requested tolerance demands it.

    .. versionadded:: 1.0.0

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
        Target accumulated phase per slab, in radians. Default: :math:`2\pi`.

    Returns
    -------
    int
        Suggested starting number of slabs (at least 1).
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import magnus

        H = np.array([[0.0, 1.0], [1.0, 0.0]])
        print('slabs suggested:',
              magnus.suggest_n_slabs(lambda t: -1j*20.0*H, 0.0, 1.0))
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

    Cheaper than ``_cumulative_integral`` when only the total
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


def _compositions(total: int, parts: int):
    r"""Yields every ordered tuple of ``parts`` positive integers summing to ``total``.

    These index the terms of one commutator group: the :math:`j`-th group of
    :math:`\Omega_n` has one term per composition of :math:`n-1` into :math:`j` parts,
    so it holds :math:`\binom{n-2}{j-1}` terms.
    """
    if parts == 1:
        yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in _compositions(total - first, parts - 1):
            yield (first,) + rest


def _nested_chain(comp, om, Bt, cache):
    r"""Evaluates :math:`[\Omega_{m_1}, [\Omega_{m_2}, \ldots [\Omega_{m_j}, A] \ldots]]`.

    Memoized on the composition, which is what makes this affordable: distinct terms share
    long suffixes (every term of the :math:`j`-th group ending in the same tail reuses one
    stored array), so each distinct nested commutator is built once no matter how many
    terms contain it.  This is the same reuse the hand-written orders 1-6 get from naming
    ``C1``, ``D11`` and friends, done automatically.
    """
    hit = cache.get(comp)
    if hit is not None:
        return hit
    if len(comp) == 1:
        value = _commutator_batched(om[comp[0]], Bt)
    else:
        value = _commutator_batched(om[comp[0]], _nested_chain(comp[1:], om, Bt, cache))
    cache[comp] = value
    return value


def _omega_integrand(n: int, om: dict, Bt: np.ndarray, cache: dict) -> np.ndarray:
    r"""Integrand of :math:`\Omega_n`, summed over every commutator group.

    Uses the closed form of the Bernoulli recursion: each term is a right-nested chain of
    lower-order :math:`\Omega_m` around :math:`A`, with the indices running over the
    compositions of :math:`n-1`, and the whole :math:`j`-th group scaled by
    :math:`B_j/j!`.  Orders 1-6 are written out inline in
    :func:`_magnus_terms_quadrature` instead, both because that path is hot and because
    keeping the published low-order expressions literal makes them checkable by eye; the
    two agree exactly, which ``tests/test_expansionterms.py`` verifies.
    """
    total = None
    for j, factor in _GROUP_FACTORS.items():
        if j > n - 1:
            continue
        group = None
        for comp in _compositions(n - 1, j):
            chain = _nested_chain(comp, om, Bt, cache)
            group = chain if group is None else group + chain
        contribution = factor*group
        total = contribution if total is None else total + contribution
    return total


def _cgroup_headroom_bytes():
    """Headroom left by a cgroup memory limit, or None if there is no limit to find.

    **This is the figure that matters wherever the library actually runs at scale.**
    ``/proc/meminfo`` is not namespaced: inside a container or a batch-scheduler cgroup
    it reports the *host's* memory, so a guard built on it alone sees far more headroom
    than the process can use.  Measured inside a 3 GiB scope on an 8 GiB machine, the
    host figure read 8.27 GiB and a 2.2 GiB allocation was waved through and then killed
    by the cgroup -- which is the exact outcome :func:`_check_output_fits` exists to
    prevent.  Docker, Kubernetes, SLURM and HPC schedulers all impose limits this way.

    Both cgroup versions are read, and in both the effective limit is the **minimum over
    the whole ancestor chain**: a limit may be set on any ancestor rather than on the
    leaf, and the tightest one binds.  Headroom is ``limit - current`` rather than the
    limit itself, because the process is already using some of it.

    Returns None when unlimited, unreadable, or absent, so that a caller can fall back to
    the host figure.  A guard that cannot measure must not block.

    .. versionadded:: 1.0.0
    """
    # cgroup v2: one unified hierarchy on the "0::" line of /proc/self/cgroup.
    # cgroup v1: a "N:memory:/path" line, mounted under /sys/fs/cgroup/memory.
    v2_path, v1_path = None, None
    try:
        with open('/proc/self/cgroup') as f:
            for line in f:
                parts = line.strip().split(':', 2)
                if len(parts) != 3:
                    continue
                if parts[0] == '0' and not parts[1]:
                    v2_path = parts[2]
                elif 'memory' in parts[1].split(','):
                    v1_path = parts[2]
    except OSError:
        return None

    def read_int(path):
        """The file's contents as an int; None for absent, unreadable or 'max'."""
        try:
            with open(path) as f:
                text = f.read().strip()
        except (OSError, ValueError):
            return None
        if text == 'max':
            return None
        try:
            value = int(text)
        except ValueError:
            return None
        # cgroup v1 spells "unlimited" as a sentinel near 2**63 rather than as a word.
        return None if value >= 2**62 else value

    best = None
    for root, rel, limit_name, usage_name in (
            ('/sys/fs/cgroup', v2_path, 'memory.max', 'memory.current'),
            ('/sys/fs/cgroup/memory', v1_path,
             'memory.limit_in_bytes', 'memory.usage_in_bytes')):
        if rel is None:
            continue
        # Walk from the process's own cgroup up to the mount root; the tightest limit
        # anywhere on the chain is the one that will kill us.
        parts = [p for p in rel.split('/') if p]
        for depth in range(len(parts), -1, -1):
            base = os.path.join(root, *parts[:depth])
            limit = read_int(os.path.join(base, limit_name))
            if limit is None:
                continue
            used = read_int(os.path.join(base, usage_name)) or 0
            headroom = max(limit - used, 0)
            best = headroom if best is None else min(best, headroom)
    return best


def _available_memory_bytes():
    """Best-effort free memory for *this* process, or None if it cannot be had cheaply.

    ``MemAvailable`` is preferred over the raw free-page count because it accounts for
    reclaimable page cache: on a machine with a warm cache the latter understates what a
    large allocation can actually get, and a guard built on it would refuse work that
    would have succeeded.

    Whatever the host reports is then **capped by any cgroup limit** applying to this
    process -- see :func:`_cgroup_headroom_bytes` for why that is not optional.  The
    minimum of the two is what an allocation can actually claim.

    Returns None rather than guessing on platforms that expose none of these.  A guard
    that cannot measure must not block.

    .. versionadded:: 1.0.0
    """
    host = None
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    host = int(line.split()[1])*1024
                    break
    except (OSError, ValueError, IndexError):
        pass
    if host is None:
        try:
            host = os.sysconf('SC_AVPHYS_PAGES')*os.sysconf('SC_PAGE_SIZE')
        except (AttributeError, ValueError, OSError):
            host = None

    cgroup = _cgroup_headroom_bytes()
    if cgroup is None:
        return host
    return cgroup if host is None else min(host, cgroup)

# Cumulative number of commutator terms the recursion evaluates through each order, from
# the term counts 1, 1, 2, 3, 5, 9, 17, 33, 65, 129 at orders one to ten.  Doubled below,
# because every commutator forms X @ Y and Y @ X before subtracting them.
_CUMULATIVE_TERMS = (0, 1, 2, 4, 7, 12, 21, 38, 71, 136, 265)

WORKING_SET_SAFETY = 2.0
"""float: fraction of available memory the quadrature working set may claim.

Matches :data:`magnus.oscprob.OUTPUT_GUARD_SAFETY`; the two guards cover different
allocations and should refuse at the same point.
"""



def _probe_dim(A, t0):
    """Matrix dimension of A, for the guard, without committing to a full evaluation."""
    try:
        return int(np.asarray(A(t0)).shape[-1])
    except Exception:
        return 0


def _quadrature_working_set_bytes(n_slabs, n_tpts, dim, order):
    r"""Bytes the cumulative-quadrature recursion will hold at once.

    The recursion works on arrays of shape ``(n_slabs, n_tpts, dim, dim)`` and keeps one
    per commutator it has evaluated, so the working set is the cell count times twice the
    cumulative term count.  Measured peak resident set, in units of one such array: 35 at
    order six, 118 at order eight, 408 at order ten, each stable to better than a per cent
    across ``n_slabs`` and ``n_tpts``.  The estimate above gives 42, 142 and 530, so it
    runs about a quarter high -- deliberately, since a guard that under-estimates does not
    guard.

    .. versionadded:: 1.0.0
    """
    order = max(1, min(int(order), len(_CUMULATIVE_TERMS) - 1))
    cells = int(n_slabs)*int(n_tpts)*int(dim)*int(dim)
    return cells*16*2*_CUMULATIVE_TERMS[order]


def _working_set_chunk(n_lead, n_tpts, dim, order, integration_method):
    r"""How many slabs the quadrature may hold at once, so the intermediates fit.

    The guard tempers the run rather than refusing it.  Each slab's :math:`\Omega` depends
    only on its own samples, so evaluating the chain a chunk at a time is exact -- the same
    reasoning that lets :data:`magnus.oscprob.BATCH_WORKING_ENTRIES` tile the energy axis,
    applied to the axis that actually overflows here.

    Returns ``n_lead`` unchanged whenever the whole chain fits, which is the ordinary case
    and costs one multiply.  Returns a smaller chunk when it does not.

    Raises
    ------
    MemoryError
        Only when a *single* slab will not fit, where no chunking helps and the caller has
        to change the request.

    .. versionadded:: 1.0.0
    """
    needed = _quadrature_working_set_bytes(n_lead, n_tpts, dim, order)
    if needed < WORKING_SET_MIN_BYTES:
        return n_lead
    available = _available_memory_bytes()
    if available is None:
        return n_lead
    budget = available/WORKING_SET_SAFETY
    if needed <= budget:
        return n_lead
    per_slab = _quadrature_working_set_bytes(1, n_tpts, dim, order)
    chunk = int(budget//per_slab) if per_slab else n_lead
    if chunk < 1:
        raise MemoryError(
            "Error in magnus: magnus._working_set_chunk: order " + str(order) + " on '"
            + str(integration_method) + "' quadrature needs "
            + f"{per_slab/2**30:.2f}" + " GiB of intermediates for a *single* slab at "
            + f"{int(n_tpts):,}" + " points, against " + f"{available/2**30:.2f}"
            + " GiB available. Chunking the chain cannot help, because one slab is already "
            "too large. Lower magnus_exp_order, or cap max_n_tpts_per_slab: the cost is the "
            "product of the points per slab and the number of commutator terms at this "
            "order. Note that n_tpts_per_slab is refined upward unless min_n_tpts_per_slab "
            "and max_n_tpts_per_slab pin it.")
    return chunk


WORKING_SET_MIN_BYTES = 64*1024*1024
"""int: below this the working set is not worth a free-memory read.  Matches
:data:`magnus.oscprob.OUTPUT_GUARD_MIN_BYTES`."""


def _magnus_terms_quadrature(
    Bt: np.ndarray,
    order: int,
    integration_method: str
) -> np.ndarray:
    r"""Magnus terms :math:`\Omega_1 \ldots \Omega_\text{order}` from samples of :math:`A`.

    Parameters
    ----------
    Bt : np.ndarray
        Samples of the rescaled matrix function, shape (..., m, d, d):
        :math:`B(t) = \text{width} \times A(t(s))` on the uniform
        normalized grid :math:`s \in [0, 1]`
        with m points, so that all integrals run over :math:`[0, 1]`.  Any
        leading axes (e.g., a slab axis) broadcast through.
    order : int
        Index of the last term computed: returns Omega_1 ... Omega_order.
        This is the cumulative path's meaning of ``order``; see
        ``magnus_expansion`` for how it maps onto a delivered order.
    integration_method : str
        'trapezoid' or 'simpson'.

    Returns
    -------
    np.ndarray
        Stacked terms, shape (order, ..., d, d).

    Notes
    -----
    Implements the Bernoulli-number recursion (see module docstring).
    The commutators :math:`C_k = [\Omega_k(s), A(s)]` and the nested combinations
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
        C1 = _commutator_batched(o1t, Bt)             # [Omega_1, A]
        o2t = integ(-0.5*C1, 2)
        terms.append(last(o2t, 2))

    if order >= 3:
        C2 = _commutator_batched(o2t, Bt)             # [Omega_2, A]
        D11 = _commutator_batched(o1t, C1)            # [Omega_1, [Omega_1, A]]
        o3t = integ(-0.5*C2 + F1*D11, 3)
        terms.append(last(o3t, 3))

    if order >= 4:
        C3 = _commutator_batched(o3t, Bt)             # [Omega_3, A]
        D12 = _commutator_batched(o1t, C2)            # [Omega_1, [Omega_2, A]]
        D21 = _commutator_batched(o2t, C1)            # [Omega_2, [Omega_1, A]]
        o4t = integ(-0.5*C3 + F1*(D12 + D21), 4)
        terms.append(last(o4t, 4))

    if order >= 5:
        C4 = _commutator_batched(o4t, Bt)             # [Omega_4, A]
        o5t = integ(
            -0.5*C4
            + F1*(_commutator_batched(o1t, C3) + _commutator_batched(o2t, C2)
                  + _commutator_batched(o3t, C1))
            + F2*_commutator_batched(o1t, _commutator_batched(o1t, D11)),
            5)
        terms.append(last(o5t, 5))

    if order >= 6:
        C5 = _commutator_batched(o5t, Bt)             # [Omega_5, A]
        o6t = integ(
            -0.5*C5
            + F1*(_commutator_batched(o1t, C4) + _commutator_batched(o2t, C3)
                  + _commutator_batched(o3t, C2) + _commutator_batched(o4t, C1))
            + F2*(_commutator_batched(o1t, _commutator_batched(o1t, D12))
                  + _commutator_batched(o1t, _commutator_batched(o1t, D21))
                  + _commutator_batched(o1t, _commutator_batched(o2t, D11))
                  + _commutator_batched(o2t, _commutator_batched(o1t, D11))),
            6)
        terms.append(last(o6t, 6))

    if order >= 7:
        # Beyond order 6 the number of terms (17, 33, 65, 129 at orders 7-10) makes writing
        # them out unreadable, so they are generated from the same recursion instead.  The
        # cache is shared across orders: a nested chain built for Omega_7 is reused by
        # Omega_8 and beyond rather than rebuilt.
        om = {1: o1t, 2: o2t, 3: o3t, 4: o4t, 5: o5t, 6: o6t}
        chain_cache = {}
        for n in range(7, order + 1):
            ont = integ(_omega_integrand(n, om, Bt, chain_cache), n)
            om[n] = ont
            terms.append(last(ont, n))

    return np.stack(terms, axis=0)


def _samples_identical(X: np.ndarray, Y: np.ndarray) -> bool:
    r"""Whether two node samples are bit-identical for *every* slab.

    Used to detect a Hamiltonian that is constant within each slab, where the
    Magnus series terminates at its first term.  The test is exact equality
    rather than a tolerance, deliberately: a *nearly* constant A still has
    non-vanishing commutators that carry real information, and dropping them
    because two samples happened to agree to some epsilon would silently lower
    the order.  Exact equality is what a piecewise-constant profile actually
    produces -- the same lookup returning the same float -- so nothing is lost
    by refusing to guess.

    ``array_equal`` short-circuits on the first differing element, so on a
    smooth profile this costs one comparison and returns.
    """
    return X.shape == Y.shape and np.array_equal(X, Y)


def _magnus_gl(
    An: np.ndarray,
    widths: Union[float, np.ndarray],
    order: int
) -> np.ndarray:
    r"""Magnus operator :math:`\Omega` from Gauss-Legendre collocation.

    Gauss-Legendre collocation Magnus integrators of order 2, 4, 6 and 8 based on
    Gauss-Legendre nodes (Blanes, Casas & Ros 2000; Blanes, Casas & Ros 2002;
    Blanes et al. 2009, Sec. 5.4).  Exact quadrature order matched to the
    truncation order, using only 1, 2, 3, or 4 evaluations of A per slab.
    Orders 6 and 8 use the commutator-optimal forms of the 2002 paper, which
    need three and six commutators, the fewest possible at each order.

    Parameters
    ----------
    An : np.ndarray
        :math:`A` evaluated at the GL nodes, shape (..., n_nodes, d, d).
    widths : float or np.ndarray
        Slab widths :math:`h`, broadcastable against the leading axes of ``An``.
    order : int
        Requested order; mapped to the smallest GL scheme with at least
        that order (1-2 -> GL1, 3-4 -> GL2, 5-6 -> GL3).

    Returns
    -------
    np.ndarray
        The total Magnus operator :math:`\Omega`, shape (..., d, d).
    """
    h = np.asarray(widths)[..., None, None]

    if order <= 2:
        # Midpoint rule: Omega = h A(t0 + h/2)
        return h*An[..., 0, :, :]

    if order <= 4:
        # Omega = (h/2)(A1 + A2) + (sqrt(3)/12) h^2 [A2, A1]
        A1 = An[..., 0, :, :]
        A2 = An[..., 1, :, :]
        if _samples_identical(A1, A2):
            # A is constant across every slab's nodes, so [A2, A1] is identically zero and
            # Omega = h A.  Not an approximation: for constant A the Magnus series terminates
            # at the first term, every later one being a commutator of A with itself.  Skipping
            # the commutator also removes its round-off, so this is very slightly *more*
            # accurate as well as cheaper.  Fires on piecewise-constant profiles -- castle
            # walls, a t_breakpoints-delimited region of uniform density -- and not on a smooth
            # one like PREM, where the two nodes genuinely differ.
            return h*A1
        return 0.5*h*(A1 + A2) + (np.sqrt(3.0)/12.0)*h*h*_commutator_batched(A2, A1)

    if order <= 6:
        # Order 6 (Blanes, Casas & Ros 2002, Eqs. 3.5-3.7): three commutators, the fewest
        # with which sixth order can be reached.
        A1 = An[..., 0, :, :]
        A2 = An[..., 1, :, :]
        A3 = An[..., 2, :, :]
        if _samples_identical(A1, A2) and _samples_identical(A2, A3):
            # Same argument as at order 4, and worth more here: the order-6 expression
            # builds three nested commutators, all of which vanish for constant A.
            return h*A1
        a1 = h*A2
        a2 = (np.sqrt(15.0)/3.0)*h*(A3 - A1)
        a3 = (10.0/3.0)*h*(A3 - 2.0*A2 + A1)
        C1 = _commutator_batched(a1, a2)
        C2 = (-1.0/60.0)*_commutator_batched(a1, 2.0*a3 + C1)
        return a1 + a3/12.0 + (1.0/240.0)*_commutator_batched(-20.0*a1 - a3 + C1, a2 + C2)

    # Order 8 (Blanes, Casas & Ros 2002, Eqs. 3.8-3.10): six commutators, again the
    # fewest possible.  The four alpha are that paper's b_1..b_4, obtained from the
    # univariate integrals B^(i) of its Eq. (3.2) -- note those carry a 1/h^i prefactor,
    # one power of h more than the 1/h^(i+1) of the 2000 paper, so the B^(i) below are
    # h times a quadrature average rather than the average itself.
    A1 = An[..., 0, :, :]
    A2 = An[..., 1, :, :]
    A3 = An[..., 2, :, :]
    A4 = An[..., 3, :, :]
    if (_samples_identical(A1, A2) and _samples_identical(A2, A3)
            and _samples_identical(A3, A4)):
        # As at orders 4 and 6, and worth most here: six commutators all vanish.
        return h*A1
    S1 = A1 + A4
    S2 = A2 + A3
    R1 = A4 - A1
    R2 = A3 - A2
    hh = 0.5*h
    B0 = hh*(_GL4_W1*S1 + _GL4_W2*S2)
    B1 = hh*(_GL4_W1*_GL4_V1*R1 + _GL4_W2*_GL4_V2*R2)
    B2 = hh*(_GL4_W1*_GL4_V1**2*S1 + _GL4_W2*_GL4_V2**2*S2)
    B3 = hh*(_GL4_W1*_GL4_V1**3*R1 + _GL4_W2*_GL4_V2**3*R2)
    a1 = 0.75*(3.0*B0 - 20.0*B2)
    a2 = 15.0*(5.0*B1 - 28.0*B3)
    a3 = -15.0*(B0 - 12.0*B2)
    a4 = -140.0*(3.0*B1 - 20.0*B3)
    C1 = (-1.0/28.0)*_commutator_batched(a1 + a3/28.0, a2 + (3.0/28.0)*a4)
    C2 = (1.0/3.0)*_commutator_batched(a1, -a3/14.0 + C1)
    C3 = _commutator_batched(a1 + a3/28.0 + C1, a2 + (3.0/28.0)*a4 + C2)
    C4 = _commutator_batched(a2, C1)
    C5 = _commutator_batched(a1 + 1.25*C1, 2.0*a3 + C3 + 0.5*C4)
    C6 = _commutator_batched(a1 + a3/12.0 - (7.0/3.0)*C1 - C3/6.0,
                             -9.0*a2 - 2.25*a4 + 63.0*C2 + C5)
    return a1 + a3/12.0 - (7.0/120.0)*C3 + (1.0/360.0)*C6


def _gl_nodes(order: int) -> np.ndarray:
    r"""Returns the Gauss-Legendre nodes on [0, 1] for the given Magnus order.

    Parameters
    ----------
    order : int
        Requested Magnus order; mapped to the smallest GL scheme with at least that order
        (1-2 -> 1 node, 3-4 -> 2 nodes, 5-6 -> 3 nodes, 7-8 -> 4 nodes).

    Returns
    -------
    np.ndarray
        GL nodes on [0, 1] (1, 2, 3, or 4 of them).
    """
    if order > MAGNUS_EXP_ORDER_MAX_GL:
        # Backstop.  _validate() reports this with a fuller message, but it is skipped when
        # validate_input=False, and silently returning the 3-node (order-6) scheme for a
        # higher requested order would be exactly the kind of quiet wrong answer that is
        # worse than an exception.
        raise ValueError(
            "Error in magnus: magnus._gl_nodes: no Gauss-Legendre scheme of order " + str(order)
            + " exists (the highest is " + str(MAGNUS_EXP_ORDER_MAX_GL)
            + "); use integration_method='trapezoid' or 'simpson'.")
    if order <= 2:
        return _GL1_NODES
    if order <= 4:
        return _GL2_NODES
    if order <= 6:
        return _GL3_NODES
    return _GL4_NODES


_SLAB_NORM_SINK = None
r"""list or None: when a caller has opened ``_deferred_slab_norm``, every ``||Omega||_2`` the
convergence check computes is collected here instead of warned about immediately.  ``None`` (and
therefore free) otherwise."""


@contextmanager
def _deferred_slab_norm():
    r"""Collect slab norms instead of warning about them, for the duration of the block.

    :func:`magnus.oscprob.osc_prob` refines a slab ladder and returns **one** level's answer,
    so warning as each level is computed reports on grids nobody receives: measured over 66
    single-point calls, some level exceeded :math:`\pi` in 46 of them but the level actually
    returned did so in only **7**.  This exists so a caller can collect the norms and emit once,
    for the level it is about to return.

    **:func:`magnus.oscprob.osc_prob` deliberately does not use it**, and the measurement is why.
    Keying the warning to the returned level was implemented and then reverted: over 168
    configurations, firings fell 70 to 53 but **true positives fell 17 to 4** while false
    positives fell only 53 to 49.  "The ladder started far from convergence" predicts a bad
    answer better than "the final grid is coarse" does, so the suppression removed most of the
    signal to remove a twelfth of the noise.  Nothing became silent either way (2 of 168 in
    both), because the cases it stopped flagging are covered by
    :class:`magnus.oscprob.ToleranceNotAchievedWarning`.

    The honest way to use the discarded signal would be a *different* warning -- "this request
    needed many refinement levels" -- rather than a quieter version of this one.

    Private, and stays private: nothing in the package uses it, and shipping public API for a
    design that was measured and rejected would be worse than keeping the knowledge here.

    Nested blocks share the outermost sink, so an inner engine's slabs are attributed to the
    level being computed rather than starting a fresh collection.

    .. versionadded:: 1.0.0

    Yields
    ------
    list of float
        Every norm seen inside the block, in the order seen.
    """
    global _SLAB_NORM_SINK
    prev = _SLAB_NORM_SINK
    sink = prev if prev is not None else []
    _SLAB_NORM_SINK = sink
    try:
        yield sink
    finally:
        _SLAB_NORM_SINK = prev


def _warn_slab_norm(nmax: float):
    r"""Warn if the slab norm proxy ``nmax`` :math:`= \max \lVert\Omega\rVert_2` is
    :math:`\geq \pi` (see :class:`MagnusConvergenceWarning`).

    Parameters
    ----------
    nmax : float
        Largest :math:`\lVert\Omega\rVert_2` (or a proxy for it) encountered across the slab(s)
        just evaluated.

    Returns
    -------
    None
    """
    if _SLAB_NORM_SINK is not None:
        # A caller is running a refinement ladder and will decide, once it knows which level it
        # is returning, whether this is worth saying.  See deferred_slab_norm.
        _SLAB_NORM_SINK.append(float(nmax))
        return
    if nmax < np.pi:
        return
    # Bucketed rather than numeric, so that the message stays one of three fixed strings and
    # Python's default filter still shows each at most once per session -- while carrying the
    # one quantity this function actually knows.  How far past pi is not the error, but it
    # separates "one slab marginally over" from "the grid is nowhere near fine enough".
    if nmax < 2.0*np.pi:
        how_far = "marginally over"
    elif nmax < 10.0*np.pi:
        how_far = "over by up to a factor of ten"
    else:
        how_far = "over by more than a factor of ten"
    warnings.warn(
        "at least one time slab is too wide for guaranteed convergence of the Magnus "
        "series (||Omega||_2 >= pi, " + how_far + "). This is a statement about the slab "
        "width, not about the answer: it reports that a sufficient condition for "
        "convergence was not met somewhere, and the error may be anywhere from negligible "
        "to large. To act on it, use more (narrower) slabs -- request a smaller rtol/atol, "
        "or raise n_slabs; raising magnus_exp_order will not help in this regime. If the "
        "profile has a density jump or a kink, pass t_breakpoints there as well: a slab "
        "straddling one is not fixed by any number of slabs. Do NOT assume the adaptive "
        "refinement has already taken care of it -- measured on a sawtooth density with "
        "rtol=atol=1e-3 explicitly requested, the refinement ran and the answer was still "
        "7.5e-03, seven times outside the tolerance asked for. Shown once per session.",
        MagnusConvergenceWarning, stacklevel=4)


def ordered_product(U: np.ndarray) -> np.ndarray:
    r"""Time-ordered product of a stack of slab operators, earliest slab first.

    Returns :math:`U_{n-1} \cdots U_1 U_0` for ``U`` of shape ``(n, d, d)`` --
    the same quantity as ``functools.reduce(np.matmul, U[::-1])`` and, because
    matrix multiplication is associative, the same value.

    The difference is how it gets there.  ``reduce`` walks the stack one matrix
    at a time, which is :math:`n-1` separate Python-level calls into NumPy for
    matrices of size 3; the array was already materialized as a single
    ``(n, d, d)`` block, so nearly all of that time is call overhead rather than
    arithmetic.  Multiplying adjacent pairs instead collapses the stack in
    :math:`\lceil\log_2 n\rceil` **batched** matmuls, each of which does its
    whole level in one call.

    Measured on unitary 3x3 stacks: 92 -> 29 us at n = 108, 1764 -> 370 us at
    n = 2048, agreeing to 8e-16 with no systematic loss of unitarity.  The gain
    grows with the slab count, which is the direction the adaptive refinement
    moves in.

    Associativity is what makes this legitimate; commutativity is *not* required
    and is *not* assumed.  Adjacent pairs are combined in order, so the operator
    ordering is preserved exactly -- an odd element is carried forward untouched
    rather than being folded in out of turn.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    U : np.ndarray
        Stack of operators, shape ``(n, d, d)``, ordered earliest slab first.

    Returns
    -------
    np.ndarray
        The ordered product, shape ``(d, d)``.
    """
    M = np.asarray(U)
    if M.ndim == 2:
        return M
    if M.shape[0] == 1:
        return M[0]
    M = M[::-1]                       # leftmost factor first
    while M.shape[0] > 1:
        n = M.shape[0]
        half = n//2
        prod = M[:2*half:2] @ M[1:2*half:2]
        M = np.concatenate([prod, M[-1:]], axis=0) if n % 2 else prod
    return M[0]


def _ordered_product_batched_core(U):  # pragma: no cover -- compiled below
    r"""Left-fold slab product of a batched operator stack, one batch at a time.

    For ``U`` of shape ``(nB, n, d, d)`` returns the ``(nB, d, d)`` stack whose
    element ``b`` is ``U[b, n-1] @ ... @ U[b, 1] @ U[b, 0]`` -- the slab crossed
    first standing rightmost, exactly the association the Python loop it
    replaces used (``Utot = Utot @ U[:, k]``, k descending).  Written to be
    compiled by numba; the pure-Python form exists only as compilation input.
    """
    nB, n, d, _ = U.shape
    out = np.empty((nB, d, d), dtype=np.complex128)
    acc = np.empty((d, d), dtype=np.complex128)
    tmp = np.empty((d, d), dtype=np.complex128)
    for b in range(nB):
        for i in range(d):
            for j in range(d):
                acc[i, j] = U[b, n - 1, i, j]
        for k in range(n - 2, -1, -1):
            for i in range(d):
                for j in range(d):
                    s = 0.0 + 0.0j
                    for m in range(d):
                        s += acc[i, m]*U[b, k, m, j]
                    tmp[i, j] = s
            acc, tmp = tmp, acc
        for i in range(d):
            for j in range(d):
                out[b, i, j] = acc[i, j]
    return out


if expmkernels.HAVE_NUMBA:
    _ordered_product_batched_kernel = expmkernels._jit(_ordered_product_batched_core)
else:                                                   # pragma: no cover
    _ordered_product_batched_kernel = None


def _ordered_product_batched(U: np.ndarray) -> np.ndarray:
    r"""Time-ordered slab product for a stack with a leading batch axis.

    The batched sibling of :func:`ordered_product`: ``U`` has shape
    ``(nB, n, d, d)`` with the slab axis second, and the return is the
    ``(nB, d, d)`` product ``U[:, n-1] @ ... @ U[:, 0]`` -- earliest slab
    rightmost, because the operators act on the state to their right.

    With numba present the product runs in a compiled kernel that keeps the
    *same left-fold association* as the Python loop it replaces; without numba
    (or on a dtype the kernel was not built for) that loop itself runs, so a
    numba-less install is bit-identical to what it always computed.

    The two installs are no longer bit-identical to *each other* on this path,
    where before this kernel they were: the association is the same, but the
    compiled kernel accumulates each matrix element as a scalar sum where BLAS
    orders the same arithmetic its own way, so probabilities can move at the
    rounding level.  Worst observed shift 1.28e-14 across 16 scan
    configurations, with every refinement decision unchanged.

    .. versionadded:: 1.0.6

    Parameters
    ----------
    U : np.ndarray
        Stack of operators, shape ``(nB, n, d, d)``, slabs ordered earliest
        first along axis 1.

    Returns
    -------
    np.ndarray
        The ordered products, shape ``(nB, d, d)``.
    """
    if (_ordered_product_batched_kernel is not None) and (U.dtype == np.complex128):
        return _ordered_product_batched_kernel(U)
    Utot = U[:, -1]
    for k in range(U.shape[1] - 2, -1, -1):
        Utot = Utot @ U[:, k]
    return Utot


valid_expm_backends = ['auto', 'numba', 'eigh']
r"""list of str: The accepted values of ``EXPM_BACKEND`` and of every
``expm_backend`` parameter.
"""


EXPM_BACKEND = 'auto'
r"""str: Module-level switch selecting how :math:`\exp(\Omega)` is computed.

Which routine exponentiates each slab.  This is not a correctness switch: the two
backends agree to about 1e-15 wherever the kernel is used, which is the accuracy either one
has -- and where it would not, it is not used: the kernel reports the conditioning of its own
characteriztic cubic and ``eigh`` answers instead.  See :data:`magnus.expmkernels.SEV_TOL`.

* ``'auto'`` (the default): the compiled Cayley-Hamilton kernel of
  :mod:`magnus.expmkernels` for 2x2 and 3x3 matrices when numba is installed,
  and ``numpy.linalg.eigh`` for everything else.  Never fails: without numba, or
  at dimension 4 and above, it is silently ``'eigh'``.
* ``'numba'``: the same, except that a missing numba is an error rather than a
  fallback -- for a caller who means to be sure the fast path is the one
  running.  Dimensions 4 and 5 still use ``eigh`` even here, because there is no
  practical closed form for a 4x4 or 5x5 Hermitian eigenproblem; 4nu and 5nu
  stay correct and are simply not accelerated.
* ``'eigh'``: ``numpy.linalg.eigh`` always, ignoring numba.  The reference route,
  and what to set when comparing the two.

``eigh`` costs about 1.25 us per 3x3 whatever the stack size, because it loops
over LAPACK internally instead of vectorizing, which makes it roughly a quarter
of a 108-slab Magnus pass.  The kernel removes that.

Setting this is the way to reach the whole package, including every
:mod:`magnus.oscprob` wrapper; the ``expm_backend`` parameter on
:func:`magnus_expansion`, :func:`evolution_operators_from_samples` and
:func:`magnus_expansion_multislab` overrides it for one call.

That includes ``n_jobs != 1``, but only because it is carried across deliberately: a module
global does not survive a process boundary, and loky re-imports magnus in each worker with
this back at its default.  ``oscprob.osc_prob_energy_baseline`` reads the value in the parent
and re-applies it inside the worker.  Anything that adds a second parallel entry point has to
do the same, or that path silently runs ``'auto'`` whatever this says.

.. versionadded:: 1.0.0
"""


def _resolve_expm_backend(expm_backend: Optional[str]) -> str:
    r"""Validates a requested backend and falls back to the module default.

    Parameters
    ----------
    expm_backend : str or None
        One of ``valid_expm_backends``; None means use ``EXPM_BACKEND``.

    Returns
    -------
    str
        The backend to use, one of ``valid_expm_backends``.

    Raises
    ------
    ValueError
        If the name is not recognized, or if ``'numba'`` was asked for by name
        and numba is not installed.
    """
    backend = EXPM_BACKEND if expm_backend is None else expm_backend
    if backend not in valid_expm_backends:
        raise ValueError(
            "Error in magnus: magnus._expm_stack: expm_backend must be one of "
            + str(valid_expm_backends) + ", not '" + str(backend) + "'.")
    # Asked for by name, so a silent downgrade would be the wrong answer to give:
    # the caller wanted to know the compiled kernel was running.  'auto' is the
    # value that promises to work anywhere, and it is the default.
    if backend == 'numba' and not expmkernels.HAVE_NUMBA:
        raise ValueError(
            "Error in magnus: magnus._expm_stack: expm_backend='numba' was requested but numba is not "
            "installed. Install it (pip install 'magnuspy[fast]', or pip install numba), "
            "or use expm_backend='auto', which falls back to 'eigh' when numba is absent.")
    return backend


def _expm_stack(Om: np.ndarray, warn_wide: bool = False,
                A_is_const: bool = False,
                expm_backend: Optional[str] = None) -> np.ndarray:
    r"""Matrix exponential of one matrix or a stack of matrices.

    If ``Om`` is anti-Hermitian (as is always the case for
    :math:`A = -i H` with a Hermitian Hamiltonian :math:`H`), the
    exponential is computed from the spectrum of the Hermitian matrix
    :math:`K = i\Omega`, by one of two routes selected by
    ``EXPM_BACKEND``/``expm_backend``: the compiled Cayley-Hamilton kernel of
    :mod:`magnus.expmkernels`, or the eigendecomposition
    :math:`\exp(\Omega) = V\, \mathrm{diag}\!\left(e^{-i\lambda}\right)\, V^\dagger`.
    Either is faster than scipy's Pade-based expm for stacks of small matrices.
    Otherwise it falls back to scipy.linalg.expm.

    Neither route is *exactly* unitary, and an earlier version of this docstring
    claimed the ``eigh`` one was.  It is not: :math:`U^\dagger U - I` measures
    4e-16 for a single 3x3 and 4e-15 for a stack of 4096, growing with stack
    size and never reaching zero, because the reconstruction from eigenvectors
    rounds like any other floating-point product.  The Cayley-Hamilton kernel is
    the same order -- measured slightly better, not worse, at every norm tested
    from 1 to 1e5 against a 40-digit reference *on unclustered spectra*, and the
    ``SEV_TOL`` gate is what keeps that true for clustered ones -- ungated, the
    closed form reaches 2.7e-07 against ``eigh``'s 3.0e-11 where a clustered
    spectrum meets a large norm, which is a corner neither a norm sweep nor a
    degeneracy sweep alone visits.  Probabilities sum to 1 to about
    1e-15, which is worth relying on; they do not sum to 1 by construction,
    which is not.

    If ``warn_wide`` is True, the eigenvalues (whose maximum modulus is
    :math:`\lVert\Omega\rVert_2`) are also used to warn about slabs too
    wide for the Magnus series to converge; for a constant :math:`A`
    (``A_is_const``) the series terminates exactly and the check is
    skipped.  Both routes return the eigenvalues, so the check costs nothing
    either way.

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
    expm_backend : str, optional
        ``'auto'``, ``'numba'`` or ``'eigh'``; see ``EXPM_BACKEND``, which
        supplies the default when this is None.

    Returns
    -------
    np.ndarray
        :math:`\exp(\Omega)`, same shape as ``Om``.
    """
    backend = _resolve_expm_backend(expm_backend)
    Om = np.asarray(Om)
    K = 1j*Om
    Kh = np.conj(np.swapaxes(K, -1, -2))
    scale = np.max(np.abs(K))
    if scale == 0.0:
        return np.broadcast_to(np.eye(Om.shape[-1], dtype=complex),
                               Om.shape).copy()
    if np.max(np.abs(K - Kh)) <= 1.e-12*scale:
        if (backend != 'eigh' and expmkernels.HAVE_NUMBA
                and expmkernels.supports_dim(Om.shape[-1])):
            U, lam, sev = expmkernels.expm_herm_stack(K)
            # The kernel forecasts, from the conditioning of its own characteriztic cubic,
            # whether it has lost more digits than eigh would; where it has, eigh answers
            # instead.  It needs a clustered spectrum AND a large norm together -- measured
            # up to 7440x worse than eigh in that corner, and no worse at all outside it --
            # so this is a rare, exact repair rather than a routine second opinion.
            #
            # Comparing scalars in Python rather than reducing an array in numpy: a
            # np.any() here would cost more than the kernel saves on a single 3x3.
            if sev > expmkernels.SEV_TOL:
                lam, V = np.linalg.eigh(K)
                Vh = np.conj(np.swapaxes(V, -1, -2))
                U = (V*np.exp(-1j*lam)[..., None, :]) @ Vh
        else:
            # eigh reads a single triangle, so the explicit symmetrization it used to be handed
            # was doing nothing the routine does not already do -- and the branch condition has
            # just established that the two triangles agree to 1e-12 anyway.  The kernel reads
            # the same triangle eigh does (the lower; eigh's UPLO defaults to 'L'), so at the
            # edge of that 1e-12 tolerance the two routes still agree, rather than diverging by
            # it because each picked a different half of a not-quite-Hermitian matrix.
            lam, V = np.linalg.eigh(K)
            Vh = np.conj(np.swapaxes(V, -1, -2))
            U = (V*np.exp(-1j*lam)[..., None, :]) @ Vh
        if warn_wide and not A_is_const:
            _warn_slab_norm(np.max(np.abs(lam)))  # ||Om||_2 = max |lambda|
        return U
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
    if (order > 6) and (integration_method in ('trapezoid', 'simpson')):
        warnings.warn(
            "magnus: Magnus order " + str(order) + " costs substantially more per slab "
            "than order 6 (roughly 2.7x at order 7, rising to about 17x at order 10, for "
            "the same grid), because the number of commutator terms roughly doubles per "
            "order. It does converge faster in the slab width, so this may still be the "
            "right trade; but narrowing the slabs at order 4 or 6 often reaches a given "
            "accuracy for less total work. Shown once per session.",
            MagnusHighOrderCostWarning, stacklevel=3)

    if (integration_method == 'gl') and (order > MAGNUS_EXP_ORDER_MAX_GL):
        raise ValueError(
            "Error in magnus: magnus._validate: integration_method 'gl' supports orders up to "
            + str(MAGNUS_EXP_ORDER_MAX_GL) + ", not " + str(order) + ". The "
            "Gauss-Legendre collocation schemes are separately derived integrators, "
            "not products of the Magnus recursion, so they do not extend with it. Use "
            "integration_method='trapezoid' or 'simpson' for orders above "
            + str(MAGNUS_EXP_ORDER_MAX_GL) + ", or lower the order.")

    if integration_method not in valid_integration_methods:
        raise ValueError(
            "Error in magnus: magnus.magnus_expansion: integration_method must be one of "
            + str(valid_integration_methods) + ", not '"
            + str(integration_method) + "'.")
    if not (1 <= order <= MAGNUS_EXP_ORDER_MAX):
        raise ValueError(
            "Error in magnus: magnus.magnus_expansion: order must be between 1 and "
            + str(MAGNUS_EXP_ORDER_MAX) + ", not " + str(order) + ".")


def magnus_expansion(
    A: Callable,
    t0: float,
    t1: float,
    n_tpts: Optional[int] = 50,
    order: Optional[int] = 2,
    integration_method: Optional[str] = 'gl',
    return_magnus_terms: Optional[bool] = False,
    validate_input: Optional[bool] = True,
    A_eval_mode: Optional[str] = None,
    expm_backend: Optional[str] = None
) -> np.ndarray:
    r"""Compute :math:`\exp(\Omega_1 + \cdots + \Omega_\text{order})` of :math:`A(t)` from
    ``t0`` to ``t1``.

    .. versionadded:: 1.0.0

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
        Requested Magnus order.  Its meaning depends on
        ``integration_method``, and so does the order actually delivered.
        On ``'gl'`` (the default) it is the classical order of the method,
        reached by the smallest collocation scheme that attains it: 1-2 use
        one node, 3-4 two, 5-6 three, 7-8 four, and a request above 8 raises
        rather than quietly returning order 8.  On ``'simpson'`` and ``'trapz'``
        it is instead the index of the last term ``Omega_k`` retained, and
        the delivered order is ``2*(order//2) + 2`` because the truncation
        is symmetric about the slab midpoint.  Measured global rates:

        ==========  ==  ==  ==  ==  ==  ==  ==  ==
        ``order``    1   2   3   4   5   6   8  10
        ==========  ==  ==  ==  ==  ==  ==  ==  ==
        ``'gl'``     2   2   4   4   6   6   8   -
        cumulative   2   4   4   6   6   8  10  12
        ==========  ==  ==  ==  ==  ==  ==  ==  ==

        So ``order=6`` is a sixth-order method on ``'gl'`` and an
        eighth-order one on ``'simpson'``.  The two extra orders come from
        the three further ``Omega`` terms the cumulative path keeps under
        the same label, not from the quadrature rule: ``order=3`` on
        ``'simpson'`` retains the same ``Omega_1 + Omega_2 + Omega_3`` as
        ``order=6`` on ``'gl'`` and converges two orders more slowly.
    integration_method : str, optional
        'gl' (Gauss-Legendre collocation; ignores ``n_tpts`` and uses 1, 2,
        3, or 4 nodes for orders <= 2, <= 4, <= 6, <= 8, respectively), 'trapezoid',
        or 'simpson'. Default: 'gl'.
    return_magnus_terms : bool, optional
        If True, also return the individual Magnus terms.  For the
        'gl' method the terms are not separable, and a single-element
        list containing the total :math:`\Omega` is returned instead.
    validate_input : bool, optional
        If True, validate ``order`` and ``integration_method``
        (raises ValueError on invalid input).
    A_eval_mode : str, optional
        Skip probing how ``A`` can be evaluated by declaring it up front
        ('vector', 'scalar', or 'constant'); see :func:`probe_eval_mode`.
        If None (default), it is probed once and detected automatically.
    expm_backend : str, optional
        Which routine exponentiates the slab: ``'auto'``, ``'numba'`` or
        ``'eigh'``.  If None (default), the module-level ``EXPM_BACKEND``
        decides.

    Returns
    -------
    np.ndarray, or (np.ndarray, np.ndarray)
        The evolution operator :math:`U = \exp(\sum_k \Omega_k)`; if
        ``return_magnus_terms`` is True, also the stacked terms.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import magnus

        H = np.array([[0.0, 1.0], [1.0, 0.0]])
        U = magnus.magnus_expansion(lambda t: -1j*H, 0.0, np.pi/4, order=4)

        print(np.round(U, 6))
        print('unitary to %.1e' % np.max(np.abs(U.conj().T @ U - np.eye(2))))

    Unitary to rounding, and that is structural rather than lucky: the
    truncated series is anti-Hermitian at any order.
"""
    if validate_input:
        _validate(order, integration_method)

    if integration_method == 'gl':
        nodes = _gl_nodes(order)
        _working_set_chunk(1, len(nodes), _probe_dim(A, t0), order, integration_method)
        width = float(t1) - float(t0)
        tnodes = t0 + width*nodes
        An, used_mode = _evaluate_A(A, tnodes, A_eval_mode)
        Om = _magnus_gl(An, width, order)
        U = _expm_stack(Om, warn_wide=True, A_is_const=(used_mode == 'constant'),
                        expm_backend=expm_backend)
        if not return_magnus_terms:
            return U
        return U, np.stack([Om], axis=0)

    _working_set_chunk(1, n_tpts, _probe_dim(A, t0), order, integration_method)
    times = np.linspace(t0, t1, n_tpts)
    At, used_mode = _evaluate_A(A, times, A_eval_mode)
    Bt = (float(t1) - float(t0))*At  # rescale to the unit interval
    magnus_terms = _magnus_terms_quadrature(Bt, order, integration_method)

    U = _expm_stack(np.sum(magnus_terms, axis=0), warn_wide=True,
                    A_is_const=(used_mode == 'constant'),
                    expm_backend=expm_backend)
    if not return_magnus_terms:
        return U
    return U, magnus_terms


def evolution_operators_from_samples(
    At: np.ndarray,
    widths: Union[list, np.ndarray],
    order: Optional[int] = 2,
    integration_method: Optional[str] = 'gl',
    A_is_const: Optional[bool] = False,
    validate_input: Optional[bool] = True,
    expm_backend: Optional[str] = None
) -> np.ndarray:
    r"""Evolution operators of a chain of slabs from precomputed samples.

    Mid-level entry point for callers that build the samples of A
    themselves -- e.g., to batch extra axes (such as the neutrino
    energy) in front of the slab axis, which this routine broadcasts
    through all operations.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    At : np.ndarray
        Samples of A, shape (..., n_slabs, m, d, d).  For the
        quadrature methods ('trapezoid'/'simpson'), the m samples of
        each slab lie on the uniform grid spanning the slab (endpoints
        included).  For 'gl', they lie on the Gauss-Legendre nodes
        (m = 1, 2, 3, or 4 for orders <= 2, <= 4, <= 6, <= 8; see
        :func:`gl_nodes`).
    widths : list or np.ndarray
        Slab widths, shape (n_slabs,) (or broadcastable to the leading
        axes of ``At`` without the last three).
    order : int, optional
        Requested Magnus order.  Its meaning depends on
        ``integration_method``, and so does the order actually delivered.
        On ``'gl'`` (the default) it is the classical order of the method,
        reached by the smallest collocation scheme that attains it: 1-2 use
        one node, 3-4 two, 5-6 three, 7-8 four, and a request above 8 raises
        rather than quietly returning order 8.  On ``'simpson'`` and ``'trapz'``
        it is instead the index of the last term ``Omega_k`` retained, and
        the delivered order is ``2*(order//2) + 2`` because the truncation
        is symmetric about the slab midpoint.  Measured global rates:

        ==========  ==  ==  ==  ==  ==  ==  ==  ==
        ``order``    1   2   3   4   5   6   8  10
        ==========  ==  ==  ==  ==  ==  ==  ==  ==
        ``'gl'``     2   2   4   4   6   6   8   -
        cumulative   2   4   4   6   6   8  10  12
        ==========  ==  ==  ==  ==  ==  ==  ==  ==

        So ``order=6`` is a sixth-order method on ``'gl'`` and an
        eighth-order one on ``'simpson'``.  The two extra orders come from
        the three further ``Omega`` terms the cumulative path keeps under
        the same label, not from the quadrature rule: ``order=3`` on
        ``'simpson'`` retains the same ``Omega_1 + Omega_2 + Omega_3`` as
        ``order=6`` on ``'gl'`` and converges two orders more slowly.
    integration_method : str, optional
        'gl', 'trapezoid', or 'simpson'. Default: 'gl'.
    A_is_const : bool, optional
        Set to True if A is constant in time to skip the (inapplicable)
        slab-width convergence warning.
    validate_input : bool, optional
        If True, validate order and integration_method.
    expm_backend : str, optional
        Which routine exponentiates each slab: ``'auto'``, ``'numba'`` or
        ``'eigh'``.  If None (default), the module-level ``EXPM_BACKEND``
        decides.

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
        return _expm_stack(Om, warn_wide=True, A_is_const=A_is_const,
                           expm_backend=expm_backend)
    Bt = w[..., None, None, None]*At        # rescale to the unit interval
    lead = Bt.shape[:-3]
    n_lead = int(np.prod(lead)) if lead else 1
    chunk = _working_set_chunk(n_lead, Bt.shape[-3], Bt.shape[-1], order, integration_method)
    if chunk >= n_lead:
        magnus_terms = _magnus_terms_quadrature(Bt, order, integration_method)
        return _expm_stack(np.sum(magnus_terms, axis=0), warn_wide=True,
                           A_is_const=A_is_const, expm_backend=expm_backend)
    # Too large to hold at once.  Each slab's Omega depends only on its own samples, so
    # evaluating the chain a chunk at a time gives the same operators for less memory.
    flat = Bt.reshape((n_lead,) + Bt.shape[-3:])
    out = np.empty((n_lead,) + Bt.shape[-2:], dtype=complex)
    for a in range(0, n_lead, chunk):
        piece = _magnus_terms_quadrature(flat[a:a + chunk], order, integration_method)
        out[a:a + chunk] = _expm_stack(np.sum(piece, axis=0), warn_wide=True,
                                       A_is_const=A_is_const, expm_backend=expm_backend)
    return out.reshape(lead + Bt.shape[-2:])


def gl_nodes(order: int) -> np.ndarray:
    r"""Returns the Gauss-Legendre nodes on [0, 1] used by the 'gl' method.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    order : int
        Requested Magnus order; mapped to the smallest GL scheme with at least that order
        (1-2 -> 1 node, 3-4 -> 2 nodes, 5-6 -> 3 nodes, 7-8 -> 4 nodes).

    Returns
    -------
    np.ndarray
        GL nodes on [0, 1] (1, 2, 3, or 4 of them).
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import magnus

        for order in (2, 4, 6, 8):
            print('order %d -> %s' % (order, np.round(magnus.gl_nodes(order), 6)))

    One to four nodes: the scheme uses the fewest that reach the order.
"""
    return _gl_nodes(order)


USE_PALINDROME = True
r"""bool: Module-level switch.

Whether a slab chain whose profile reads the same from either end may be built by evaluating
:math:`A` on its first half only, the mirrored half following by reversal.  True by default:
every Earth chord qualifies, because a chord through a spherically symmetric Earth meets every
radius twice.

Set it to ``False`` to evaluate every slab in full.  **This is not a correctness switch**, but
neither is it a no-op: the two routes agree to a few times 1e-15 rather than bitwise, because the
mirrored slab's nodes are reached as ``(L - b) + h*s`` on one route and ``a + h*s`` on the other,
which are different floating-point expressions for the same real number.  Set it to False to ask
for the plain per-slab evaluation when a comparison needs one.

The saving is halved evaluations of the caller's Hamiltonian, so it is worth most where that
Hamiltonian is expensive: with ``f`` the share of slab time spent inside it, the speed-up is
about :math:`1/(1 - f/2)`.

.. versionadded:: 1.0.0
"""


def palindromic(*arrays: np.ndarray) -> bool:
    r"""Returns whether every array given reads the same both ways.

    .. versionadded:: 1.0.0

    The comparison is exact, deliberately.  The saving relies on the mirrored slab's inputs being
    *identical* to the reversal of its partner's, which follows from identical inputs and from
    nothing weaker; a tolerance here would silently return a different answer for a
    nearly-symmetric profile, which is the one thing an optimization must never do.  It is the
    producer's business to make a profile exactly symmetric rather than nearly so.

    This mirrors ``fastkernels.palindromic`` in NuOscProbExact, deliberately, down to treating an
    empty call and any array shorter than two entries as trivially palindromic.

    Parameters
    ----------
    arrays : np.ndarray
        Arrays to test, given as separate arguments and each reversed along its first axis.

    Returns
    -------
    bool
        Whether every array equals its own reverse exactly.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import magnus

        print(magnus.palindromic(np.array([1.0, 2.0, 1.0])))
        print(magnus.palindromic(np.array([1.0, 2.0, 3.0])))
    """
    for array in arrays:
        a = np.asarray(array)
        if a.shape[0] > 1 and not np.array_equal(a, a[::-1]):
            return False
    return True


_MIRROR_GRID_ULP = 4.0
r"""float: Bound, in units of ``n_slabs*eps``, on how far a declared-symmetric slab grid may
depart from its own mirror and still be mirrored.

The scaling is measured, not assumed.  ``edges[:, 1] - edges[:, 0]`` on a grid built as
``w0*arange(n+1)`` returns a relative asymmetry of **0.30*n*eps**, holding to two digits from 31
slabs to 1024; the same construction after a mirror-average of the widths returns exactly zero at
most counts, and the real PREM chord's ``t_breakpoints`` (15 slabs) returns 1.8e-16.  Four gives
roughly thirteen times headroom over the measured worst case while still admitting only
floating-point noise: even at 1024 slabs the bound is ~1e-13 relative, which no physically uneven
grid approaches.

See :func:`_mirror_applies` for why this one comparison is a tolerance when :func:`palindromic`
is exact."""


def _mirror_applies(edges: np.ndarray, widths: np.ndarray,
                    symmetric_over: Optional[tuple]) -> bool:
    r"""Whether the mirrored evaluation is valid for this slab chain.

    ``symmetric_over`` is the caller's declaration that :math:`A(t) = A(lo + hi - t)` on
    ``(lo, hi)``.  A declaration is not enough on its own: the engines above this layer call the
    Magnus routines on *sub*-ranges of a profile (the cumulative scan, the adiabatic and
    interaction-picture paths), and a sub-range of a symmetric profile is not itself symmetric.
    So the chain must be checked to span exactly the declared interval, and to be palindromic in
    its widths -- both exactly, never within a tolerance, for the reason given in
    :func:`palindromic`.

    A widths test alone would **not** be sufficient even with a symmetric profile, and is the
    trap this function exists to avoid: a monotonic (solar-like) profile on a uniform grid has
    palindromic widths, and mirroring it is wrong by 3.3e-01.  It is the conjunction of the
    declaration with the span check that carries the correctness.

    **Why the width test is a tolerance here, when :func:`palindromic` is exact.**  Magnus derives
    its widths as ``edges[:, 1] - edges[:, 0]``, and that subtraction does not preserve symmetry:
    a grid built from a single width ``w0`` as ``w0*arange(n+1)`` comes back with **six to nine
    distinct** width values, none of them bitwise palindromic.  An exact test on them therefore
    almost never passes, and the optimization would ship as a silent no-op.  NuOscProbExact can
    keep its test exact because it carries ``widths`` as an array its producer symmetrizes
    (``earth._earth_slabs_cached``); here the producer cannot reach past the subtraction.

    So the roles differ from that project's, deliberately: **the declaration is the correctness
    criterion**, and the width comparison is a consistency check on the declaration rather than
    the thing establishing it.  The bound is a few ulp -- tight enough that only floating-point
    noise passes, so a genuinely uneven grid still takes the ordinary path -- and the mirrored
    branch then *forces* the widths it uses to be exactly palindromic rather than trusting them,
    which is the same "make it exact rather than tolerate near-exactness" discipline applied at
    the only place that can apply it.
    """
    if not USE_PALINDROME or symmetric_over is None or edges.shape[0] < 2:
        return False
    lo, hi = symmetric_over
    if not (edges[0, 0] == lo and edges[-1, 1] == hi):
        return False
    scale = np.max(np.abs(widths))
    if not np.isfinite(scale) or scale == 0.0:
        return False
    bound = _MIRROR_GRID_ULP*widths.shape[0]*np.finfo(float).eps*scale
    return bool(np.max(np.abs(widths - widths[::-1])) <= bound)


def magnus_expansion_multislab(
    A: Callable,
    t_slab_edges: Union[list, np.ndarray],
    n_tpts_per_slab: Optional[int] = 50,
    order: Optional[int] = 2,
    integration_method: Optional[str] = 'gl',
    validate_input: Optional[bool] = True,
    A_eval_mode: Optional[str] = None,
    symmetric_over: Optional[tuple] = None,
    expm_backend: Optional[str] = None
) -> np.ndarray:
    r"""Compute the evolution operators of all time slabs at once.

    Vectorized (batched) version of :func:`magnus_expansion` for a
    chain of time slabs: A is evaluated for all slabs in a single call
    (when it supports array input), and the quadrature, commutator
    algebra, and matrix exponentials are evaluated as batched NumPy
    operations with the slab axis leading.  This is much faster than
    calling :func:`magnus_expansion` slab by slab.

    .. versionadded:: 1.0.0

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
        Requested Magnus order.  Its meaning depends on
        ``integration_method``, and so does the order actually delivered.
        On ``'gl'`` (the default) it is the classical order of the method,
        reached by the smallest collocation scheme that attains it: 1-2 use
        one node, 3-4 two, 5-6 three, 7-8 four, and a request above 8 raises
        rather than quietly returning order 8.  On ``'simpson'`` and ``'trapz'``
        it is instead the index of the last term ``Omega_k`` retained, and
        the delivered order is ``2*(order//2) + 2`` because the truncation
        is symmetric about the slab midpoint.  Measured global rates:

        ==========  ==  ==  ==  ==  ==  ==  ==  ==
        ``order``    1   2   3   4   5   6   8  10
        ==========  ==  ==  ==  ==  ==  ==  ==  ==
        ``'gl'``     2   2   4   4   6   6   8   -
        cumulative   2   4   4   6   6   8  10  12
        ==========  ==  ==  ==  ==  ==  ==  ==  ==

        So ``order=6`` is a sixth-order method on ``'gl'`` and an
        eighth-order one on ``'simpson'``.  The two extra orders come from
        the three further ``Omega`` terms the cumulative path keeps under
        the same label, not from the quadrature rule: ``order=3`` on
        ``'simpson'`` retains the same ``Omega_1 + Omega_2 + Omega_3`` as
        ``order=6`` on ``'gl'`` and converges two orders more slowly.
    integration_method : str, optional
        'gl', 'trapezoid', or 'simpson'. Default: 'gl'.
    validate_input : bool, optional
        If True, validate input (raises ValueError on invalid input).
    A_eval_mode : str, optional
        Skip probing how ``A`` can be evaluated by declaring it up front
        ('vector', 'scalar', or 'constant'); see :func:`probe_eval_mode`.
        If None (default), it is probed once and detected automatically.
    symmetric_over : tuple, optional
        Caller's declaration that ``A(t) == A(lo + hi - t)`` on ``(lo, hi)``.
        When given, and when the slab chain is found to span exactly that
        interval with exactly palindromic widths, ``A`` is evaluated on the
        first half of the slabs only and the rest follows by reversal --
        halving the calls to the caller's Hamiltonian.  Ignored when
        :data:`USE_PALINDROME` is False.

        This is a **declaration, not a test**: it is not checked, and cannot
        be cheaply, since testing it would require the evaluations it exists
        to avoid.  Declaring it of a profile that is not symmetric returns a
        silently wrong answer -- measured at 3.3e-01 on a monotonic profile.
        It is therefore not a user-facing knob: it is set by the Earth entry
        points, where the symmetry is a fact of chord geometry rather than a
        claim.  See ``docs/dev/PLAN_PALINDROMIC_PROFILES.md`` section 3d(ii).

    expm_backend : str, optional
        Which routine exponentiates each slab: ``'auto'``, ``'numba'`` or
        ``'eigh'``.  If None (default), the module-level ``EXPM_BACKEND``
        decides.

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
                "Error in magnus: magnus.magnus_expansion_multislab: all slabs must have "
                "t1 >= t0.")

    if integration_method == 'gl':
        s = _gl_nodes(order)                            # (k,) GL nodes
    else:
        s = np.linspace(0.0, 1.0, n_tpts_per_slab)      # normalized grid

    if _mirror_applies(edges, widths, symmetric_over):
        # Evaluate the first half only; the mirrored half is the same samples read backwards.
        # Both node sets are symmetric within their slab (Gauss-Legendre nodes are, and so is
        # linspace(0, 1, m)), so reversing the sample axis lands on the mirror slab's own nodes.
        #
        # n_half counts the middle slab in when the count is odd: that slab straddles the center,
        # is its own mirror, and is evaluated forward like any other.  Writing the mirrored block
        # as ``At[n_half:]`` rather than ``At[n_slabs-n_half:]`` is what keeps it from being
        # overwritten -- the two differ only for odd counts, where the latter aliases the middle
        # slab and, on an uninitialized array, returns whatever was in memory.
        n_slabs = edges.shape[0]
        n_half = (n_slabs + 1)//2
        tgrid = edges[:n_half, :1] + widths[:n_half, None]*s
        At_half, used_mode = _evaluate_A(A, tgrid, A_eval_mode)
        At = np.empty((n_slabs,) + At_half.shape[1:], dtype=At_half.dtype)
        At[:n_half] = At_half
        At[n_half:] = At_half[:n_slabs - n_half][::-1, ::-1]
        # Force the widths exactly palindromic rather than trusting them to be: the mirrored
        # slab is being given its partner's samples, so it must be given its partner's width
        # too, or the two halves are scaled by numbers differing in the last bits.  This is the
        # counterpart of NuOscProbExact's ``w = (w + w[::-1])/2``, applied here because the
        # subtraction that produces our widths is downstream of anything the producer controls.
        widths = np.concatenate([widths[:n_half],
                                 widths[:n_slabs - n_half][::-1]])
    else:
        tgrid = edges[:, :1] + widths[:, None]*s            # (n_slabs, m)
        At, used_mode = _evaluate_A(A, tgrid, A_eval_mode)  # (n_slabs, m, d, d)

    return evolution_operators_from_samples(At, widths, order,
        integration_method, A_is_const=(used_mode == 'constant'),
        validate_input=False, expm_backend=expm_backend)


__all__ = [
    'MagnusConvergenceWarning',
    'ordered_product',
    'MagnusHighOrderCostWarning',
    'ScalarHamiltonianWarning',
    'B',
    'F1',
    'F2',
    'f1',
    'f2',
    'MAGNUS_EXP_ORDER_MAX',
    'valid_integration_methods',
    'commutator',
    'probe_eval_mode',
    'cached_eval_mode',
    'suggest_n_slabs',
    'magnus_expansion',
    'evolution_operators_from_samples',
    'gl_nodes',
    'magnus_expansion_multislab',
    'USE_PALINDROME',
    'palindromic',
    'EXPM_BACKEND',
    'valid_expm_backends',
]
