# -*- coding: utf-8 -*-
r"""adiabatic.py

Contains the adiabatic-transport-plus-Magnus-patch ("hybrid") propagator
used as an alternative strategy to compute the neutrino evolution
operator when the Hamiltonian is position-dependent and develops an
extreme accumulated phase (e.g., low-energy solar neutrinos crossing an
MSW resonance), the regime in which the plain slab-refinement Magnus
engine in :mod:`magnus.oscprob` needs a very large number of slabs (and
may raise ``ToleranceNotAchievedWarning``).

Physical idea: away from an eigenvalue crossing (or narrowly-avoided
crossing) of the instantaneous Hamiltonian :math:`H(l)`, the adiabatic
theorem says the evolution operator is well approximated by transport in
the *instantaneous eigenbasis* of :math:`H(l)` -- a dynamical phase (the
integral of the eigenvalues) plus a geometric (Berry) phase, both cheap
to compute on a coarse grid regardless of how large the accumulated
phase is. Near a genuine MSW resonance, the adiabatic approximation
breaks down over a narrow window, which is patched with an exact,
short-baseline Magnus computation (:func:`magnus.magnus.magnus_expansion_multislab`,
the package's own, already-unitary integrator). The two pieces are
stitched together with the exact composition law of quantum evolution,
:math:`U(l_2, l_0) = U(l_2, l_1) U(l_1, l_0)`, so the result is exactly
unitary regardless of the approximation's accuracy.

Where a patch is needed is decided by an *exact* Hellmann-Feynman
diagnostic (no finite-differenced eigenvectors, which are gauge-
ambiguous), so this applies to any Hermitian Hamiltonian of any
dimension, with any number of simultaneous or sequential resonances --
see :doc:`/adiabatic_strategy` for the full derivation, validation, and
worked examples.

This module is self-contained: it depends only on :mod:`magnus.magnus`
(the Magnus-expansion core), not on :mod:`magnus.oscprob`, so it can be
used directly on any Hamiltonian function, independently of the rest of
the oscillation-probability API. :mod:`magnus.oscprob` calls
:func:`hybrid_propagator` internally when ``strategy='hybrid'`` or
``strategy='auto'`` (the default) is passed to
:func:`magnus.oscprob.osc_prob_matter_std_potential`,
:func:`magnus.oscprob.osc_prob_matter_nsi`, or
:func:`magnus.oscprob.osc_prob_liv` (and, transitively, to every
``osc_prob_*_sun``/``osc_prob_*_sun_nsi``/``osc_prob_*_sun_liv`` wrapper),
and also when it is passed to the fully generic user-Hamiltonian entry
points :func:`magnus.oscprob.osc_prob_sun` and
:func:`magnus.oscprob.osc_prob_earth` (via
``magnus.oscprob._osc_prob_with_potential``).  For ``osc_prob_earth``
the hybrid path is normally declined, since a real Earth trajectory
supplies PREM layer breakpoints; see :doc:`/adiabatic_strategy`.

Routine listings
----------------

    * adiabatic_propagator - Evolution operator via pure adiabatic
           (instantaneous-eigenbasis) transport, no resonance patching
    * find_resonance_candidates - Locates every exact eigenvalue-gap
           critical point of H(l) via the Hellmann-Feynman theorem
    * find_nonadiabatic_windows - Filters/grows/merges candidates into
           position windows that need a Magnus patch
    * hybrid_propagator - Adiabatic transport with Magnus patches at
           non-adiabatic windows, self-certified against successive
           refinement of every internal tolerance knob
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from functools import reduce
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import simpson

import magnus.magnus as magnuscore


def _eigs_along(H_func: Callable, ls: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    r"""Diagonalizes ``H_func`` on a grid, with eigenvectors phase-fixed by discrete parallel
    transport.

    Complex dtype is forced before every ``eigh`` call: a genuinely real-valued Hamiltonian (no
    CP violation) is a legitimate special case, not a reason to special-case the code path, and
    ``eigh`` on a real array returns real eigenvectors that later fail to hold a complex parallel-
    transport phase.

    At each step past the first, eigenvector ``k`` is multiplied by the complex phase that makes
    its overlap with the previous step's eigenvector ``k`` real and positive. This is the discrete
    analogue of parallel transport and implicitly captures the geometric (Berry) phase, exactly,
    with no separate formula: the dynamical phase (see :func:`adiabatic_propagator`) and this
    phase-fixing are jointly equivalent to solving the adiabatic evolution equation.

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square matrix.
    ls : np.ndarray
        Positions at which to diagonalize ``H_func``, shape ``(n,)``.

    Returns
    -------
    (np.ndarray, np.ndarray)
        Eigenvalues, shape ``(n, d)``, and phase-fixed eigenvectors, shape ``(n, d, d)`` (each
        ``W[i, :, k]`` is the ``k``-th eigenvector at ``ls[i]``).
    """
    n = len(ls)
    H0 = np.asarray(H_func(ls[0]), dtype=complex)
    d = H0.shape[-1]
    lam = np.empty((n, d))
    W = np.empty((n, d, d), dtype=complex)
    lam[0], W[0] = np.linalg.eigh(H0)
    for i in range(1, n):
        Hi = np.asarray(H_func(ls[i]), dtype=complex)
        li, Wi = np.linalg.eigh(Hi)
        for k in range(d):
            overlap = np.vdot(W[i - 1, :, k], Wi[:, k])
            phase = overlap / abs(overlap) if abs(overlap) > 1e-14 else 1.0
            Wi[:, k] *= np.conj(phase)
        lam[i], W[i] = li, Wi
    return lam, W


def adiabatic_propagator(H_func: Callable, l0: float, l1: float,
    n_points: Optional[int] = 201) -> np.ndarray:
    r"""Computes the evolution operator via pure adiabatic (instantaneous-eigenbasis) transport.

    Diagonalizes ``H_func`` on a grid of ``n_points`` positions between ``l0`` and ``l1``,
    integrates each eigenvalue's dynamical phase with Simpson's rule (trapezoidal quadrature
    leaves a spurious residual that can look like a physics limit but is pure quadrature error),
    and reassembles the evolution operator in the original (flavor) basis:

    .. math::

       U(l_1, l_0) \approx W(l_1)\, \mathrm{diag}\!\left(e^{-i\Phi_k}\right)\, W(l_0)^\dagger ,
       \qquad \Phi_k = \int_{l_0}^{l_1} \lambda_k(l)\, dl ,

    with :math:`W(l)` the (parallel-transported; see ``_eigs_along``) matrix of instantaneous
    eigenvectors of :math:`H(l)` and :math:`\lambda_k(l)` its eigenvalues. This is *exact* in the
    strict adiabatic limit (no eigenvalue crossing or narrowly-avoided crossing along the
    trajectory) and remains unitary by construction (a diagonal phase conjugated by unitary
    matrices) regardless of grid density -- the only thing ``n_points`` controls is how well the
    quadrature/parallel-transport approximate the continuum limit, not whether the result is
    unitary. See :func:`hybrid_propagator` for what to do when the trajectory does cross a
    resonance.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix. May be
        real- or complex-valued.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    n_points : int, optional
        Number of positions at which to diagonalize ``H_func`` between ``l0`` and ``l1``.
        Default: 201.

    Returns
    -------
    np.ndarray
        The evolution operator, exactly unitary.
    """
    if l1 == l0:
        d = np.asarray(H_func(l0)).shape[-1]
        return np.eye(d, dtype=complex)
    ls = np.linspace(l0, l1, n_points)
    lam, W = _eigs_along(H_func, ls)
    d = lam.shape[1]
    Phi = np.array([simpson(lam[:, k], x=ls) for k in range(d)])
    return W[-1] @ np.diag(np.exp(-1j * Phi)) @ W[0].conj().T


def _dH_dl(H_func: Callable, l: float, h: float,
    bounds: Optional[Tuple[float, float]] = None) -> np.ndarray:
    r"""Ordinary real central finite difference of ``H_func`` at ``l``, step ``h``.

    Deliberately **not** complex-step differentiation (``Im[H(l+ih)]/h``): that technique is
    valid only for functions that are real-valued at real input, and ``H_func`` here is routinely
    complex-valued at real ``l`` (e.g., any nonzero CP-violating phase). Applying it anyway
    divides the Hamiltonian's l-independent complex entries by the (tiny) step ``h``, which blows
    up to astronomically wrong derivatives -- a subtle, easy-to-miss trap, caught by comparing
    against this real, always-valid alternative. This real finite difference is robust because
    ``dH/dl`` is always smooth (tied to the smoothness of the underlying density/potential
    profile), independent of how sharp the resulting *eigenvalue* crossing is -- the sharpness
    lives entirely in the gap, in the denominator of ``_point_adiabaticity``, never in
    ``dH/dl`` itself.

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position.
    l : float
        Position at which to evaluate the derivative.
    h : float
        Absolute finite-difference step.
    bounds : (float, float), optional
        Interval ``(l0, l1)`` outside which ``H_func`` must never be evaluated. At a position
        within ``h`` of either end, the stencil is made one-sided so that it stays inside,
        rather than reaching past the boundary. Callers should pass this whenever ``H_func``
        is only defined on the physical domain -- a radial profile undefined for negative
        radius, or one that raises beyond a maximum radius, as
        :func:`magnus.earth.density_matter_func_prem` does. If None, an unclamped central
        difference is used.

    Returns
    -------
    np.ndarray
        Approximation to :math:`dH/dl` at ``l``.
    """
    lm, lp = l - h, l + h
    if bounds is not None:
        lo, hi = bounds
        if lm < lo:
            lm, lp = lo, min(lo + 2.0 * h, hi)
        elif lp > hi:
            lm, lp = max(hi - 2.0 * h, lo), hi
    span = lp - lm
    if span <= 0.0:
        return np.zeros_like(np.asarray(H_func(l), dtype=complex))
    Hp = np.asarray(H_func(lp), dtype=complex)
    Hm = np.asarray(H_func(lm), dtype=complex)
    return (Hp - Hm) / span


def find_resonance_candidates(H_func: Callable, l0: float, l1: float,
    n_probe: Optional[int] = 200, fd_step_frac: Optional[float] = 1e-6) -> List[Dict]:
    r"""Locates every exact eigenvalue-gap critical point of ``H_func`` between ``l0`` and ``l1``.

    For every pair of levels :math:`(j, k)`, scans for sign changes of

    .. math::

       f_{jk}(l) = \langle v_j(l)|\, dH/dl\, |v_j(l)\rangle - \langle v_k(l)|\, dH/dl\, |v_k(l)\rangle ,

    refined by bisection to machine precision in position. By the Hellmann-Feynman theorem,
    :math:`d\lambda_k/dl = \langle v_k|\, dH/dl\, |v_k\rangle` *exactly* (no eigenvector finite
    difference, which would be gauge-ambiguous and fragile), so a sign change of :math:`f_{jk}`
    is an exact critical point of the gap :math:`\lambda_j - \lambda_k` -- a genuine crossing or
    near-crossing candidate, for *any* Hermitian ``H_func`` of *any* dimension, with no assumption
    of a separable or otherwise special structure. Every pair is scanned, so any number of
    simultaneous or sequential resonances (between any pair of levels) are all found.

    A returned candidate is a structural fact about ``H_func`` (an extremum of that pair's gap);
    whether it is actually non-adiabatic (whether it needs a Magnus patch) is a separate question,
    answered by :func:`find_nonadiabatic_windows`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    n_probe : int, optional
        Number of positions on the initial scan grid used to bracket sign changes. Default: 200.
    fd_step_frac : float, optional
        Finite-difference step for ``_dH_dl``, as a fraction of ``l1 - l0``. Default: 1e-6.

    Returns
    -------
    list of dict
        One entry per candidate, with keys ``'l'`` (position), ``'j'``, ``'k'`` (the level
        indices, ``j < k``), and ``'gap'`` (:math:`\lambda_k - \lambda_j` at that position),
        sorted by position.
    """
    ls = np.linspace(l0, l1, n_probe)
    h = (l1 - l0) * fd_step_frac
    n = len(ls)
    Hs = np.array([np.asarray(H_func(l), dtype=complex) for l in ls])
    d = Hs.shape[-1]
    lam = np.empty((n, d))
    W = np.empty((n, d, d), dtype=complex)
    dH = np.empty((n, d, d), dtype=complex)
    bounds = (l0, l1)
    for i in range(n):
        lam[i], W[i] = np.linalg.eigh(Hs[i])
        dH[i] = _dH_dl(H_func, ls[i], h, bounds)

    def f_pair(l: float, j: int, k: int) -> float:
        H = np.asarray(H_func(l), dtype=complex)
        _, Wi = np.linalg.eigh(H)
        dHl = _dH_dl(H_func, l, h, bounds)
        vj, vk = Wi[:, j], Wi[:, k]
        return float(np.real(np.vdot(vj, dHl @ vj) - np.vdot(vk, dHl @ vk)))

    candidates = []
    for j in range(d):
        for k in range(j + 1, d):
            fjk = np.real(np.einsum('ni,nij,nj->n', np.conj(W[:, :, j]), dH, W[:, :, j])
                - np.einsum('ni,nij,nj->n', np.conj(W[:, :, k]), dH, W[:, :, k]))
            sgn = np.sign(fjk)
            changes = np.where(np.diff(sgn) != 0)[0]
            for idx in changes:
                a, b = ls[idx], ls[idx + 1]
                fa, fb = f_pair(a, j, k), f_pair(b, j, k)
                if fa == 0.0:
                    l_star = a
                elif fb == 0.0:
                    l_star = b
                else:
                    for _ in range(60):
                        m = 0.5 * (a + b)
                        fm = f_pair(m, j, k)
                        if np.sign(fm) == np.sign(fa):
                            a, fa = m, fm
                        else:
                            b, fb = m, fm
                    l_star = 0.5 * (a + b)
                H_star = np.asarray(H_func(l_star), dtype=complex)
                lam_star = np.linalg.eigvalsh(H_star)
                gap = float(lam_star[k] - lam_star[j])
                candidates.append({'l': l_star, 'j': j, 'k': k, 'gap': gap})
    candidates.sort(key=lambda c: c['l'])
    return candidates


def _point_adiabaticity(H_func: Callable, l: float, j: int, k: int, fd_step: float,
    bounds: Optional[Tuple[float, float]] = None) -> float:
    r"""Adiabaticity parameter :math:`\gamma_{jk}(l) = |\langle v_j|\, dH/dl\, |v_k\rangle| / (\lambda_k - \lambda_j)^2`
    (Landau-Zener form), computed exactly from the Hellmann-Feynman off-diagonal matrix element --
    no eigenvector finite difference. Large :math:`\gamma` signals a narrowly-avoided (or exact)
    crossing where the adiabatic approximation breaks down; ``fd_step`` is an *absolute* step
    (unlike ``fd_step_frac`` elsewhere), since callers evaluate this at positions found by
    bisection, arbitrarily close together. ``bounds``, if given, keeps the finite-difference
    stencil inside the physical domain (see ``_dH_dl``).
    """
    H = np.asarray(H_func(l), dtype=complex)
    lam, W = np.linalg.eigh(H)
    dH = _dH_dl(H_func, l, fd_step, bounds)
    vj, vk = W[:, j], W[:, k]
    coupling = np.abs(np.vdot(vj, dH @ vk))
    gap = abs(lam[k] - lam[j])
    return coupling / gap**2 if gap > 0 else np.inf


def _estimate_window_bounds(H_func: Callable, l_star: float, j: int, k: int, l0: float, l1: float,
    threshold: float, fd_step: float, safety_factor: Optional[float] = 2.0,
    max_doublings: Optional[int] = 60) -> Tuple[float, float]:
    r"""Grows a window outward from a candidate position until the adiabaticity parameter drops
    below ``threshold``, then pads it by ``safety_factor``.

    Growing by doubling (rather than tying the window width to the search-grid spacing used by
    :func:`find_resonance_candidates`) makes the window a property of the physical transition
    width alone: the same physical case gives the same window regardless of ``n_probe``.
    """
    def grow(sign: float) -> float:
        step = fd_step
        l_edge = l_star
        for _ in range(max_doublings):
            l_try = l_star + sign * step
            if sign > 0 and l_try >= l1:
                return l1
            if sign < 0 and l_try <= l0:
                return l0
            if _point_adiabaticity(H_func, l_try, j, k, fd_step, (l0, l1)) < threshold:
                l_edge = l_try
                break
            step *= 2.0
        else:
            return l1 if sign > 0 else l0
        width = abs(l_edge - l_star)
        l_pad = l_star + sign * safety_factor * width
        return min(l_pad, l1) if sign > 0 else max(l_pad, l0)
    return grow(-1.0), grow(+1.0)


def find_nonadiabatic_windows(H_func: Callable, l0: float, l1: float,
    threshold: Optional[float] = 0.1, n_probe: Optional[int] = 200,
    fd_step_frac: Optional[float] = 1e-6) -> Tuple[List[Tuple[float, float]], List[Dict]]:
    r"""Finds every position window along ``[l0, l1]`` where ``H_func`` needs a Magnus patch.

    Calls :func:`find_resonance_candidates`, evaluates the adiabaticity parameter
    :math:`\gamma_{jk}` (see ``_point_adiabaticity``) at each candidate, grows a window around
    every candidate with :math:`\gamma_{jk} > \text{threshold}` (see ``_estimate_window_bounds``),
    and merges any windows that overlap or touch -- so two (or more) resonances close enough
    together are correctly folded into a single patch, rather than either double-counted or
    (worse) silently dropped. This works for any number of simultaneous or sequential resonances,
    between any pair of levels.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    threshold : float, optional
        Adiabaticity parameter above which a candidate is treated as non-adiabatic. Default: 0.1.
    n_probe : int, optional
        Forwarded to :func:`find_resonance_candidates`. Default: 200.
    fd_step_frac : float, optional
        Forwarded to :func:`find_resonance_candidates`. Default: 1e-6.

    Returns
    -------
    (list of (float, float), list of dict)
        The merged, non-overlapping windows (each a ``(l_b, l_c)`` pair, sorted by position), and
        the candidate list from :func:`find_resonance_candidates`, each entry additionally
        carrying its evaluated ``'gamma'``.
    """
    candidates = find_resonance_candidates(H_func, l0, l1, n_probe=n_probe,
        fd_step_frac=fd_step_frac)
    fd_step = (l1 - l0) * fd_step_frac
    windows = []
    for c in candidates:
        gamma = _point_adiabaticity(H_func, c['l'], c['j'], c['k'], fd_step, (l0, l1))
        c['gamma'] = gamma
        if gamma > threshold:
            l_b, l_c = _estimate_window_bounds(H_func, c['l'], c['j'], c['k'], l0, l1, threshold,
                fd_step)
            windows.append([l_b, l_c])

    # Sweep the probe grid as well, not only the gap extrema above.  A gap extremum is where the
    # *gap* is stationary, which is not where gamma = |<v_j|dH/dl|v_k>| / gap^2 peaks: on a
    # rapidly varying profile the coupling can be large between the extrema, and evaluating only
    # at them understates the maximum badly -- measured at 196x on a 3nu NSI profile modulated by
    # a strong sine (3.6e-04 at the extrema against 7.0e-02 along the path).
    #
    # Without this the failure is silent rather than merely inaccurate.  No window ever opens, so
    # successive refinements differ only in the adiabatic-transport grid, converge to the same
    # wrong adiabatic limit, agree with each other, and hybrid_propagator certifies a result that
    # was off by 4.3e-02 against solve_ivp.  Lowering the threshold cannot rescue it, because the
    # threshold is only ever compared against values sampled where gamma happens to be small.
    #
    # Reuses the eigendecomposition already needed for the sweep, so the cost is one extra pass
    # over the probe grid rather than a new one per pair per point.
    ls_probe = np.linspace(l0, l1, n_probe)
    Hs = np.array([np.asarray(H_func(l), dtype=complex) for l in ls_probe])
    lam_p, W_p = np.linalg.eigh(Hs)
    dH_p = np.array([_dH_dl(H_func, l, fd_step, (l0, l1)) for l in ls_probe])
    d_p = Hs.shape[-1]
    for j in range(d_p):
        for k in range(j + 1, d_p):
            vj, vk = W_p[:, :, j], W_p[:, :, k]
            coupling = np.abs(np.einsum('ni,nij,nj->n', np.conj(vj), dH_p, vk))
            gap = np.abs(lam_p[:, k] - lam_p[:, j])
            gamma_p = np.where(gap > 0.0, coupling/np.where(gap > 0.0, gap, 1.0)**2, np.inf)
            over = np.where(gamma_p > threshold)[0]
            if over.size == 0:
                continue
            # One window per *contiguous run* of exceedance, grown from that run's peak --
            # not one per exceeding point.  Growing from every point and merging pads each
            # window outward independently, and enough of them bridge the quiet stretch
            # between two genuinely separate crossings: on the two-crossing fixture in
            # tests/test_avgprob.py that collapsed 2 windows into 1 spanning almost the whole
            # trajectory, which destroys the crossing structure the averaged-probability
            # report is built on.  Runs also make this cheap: two growth searches there
            # rather than forty-two.
            for run in np.split(over, np.where(np.diff(over) != 1)[0] + 1):
                peak = int(run[np.argmax(gamma_p[run])])
                l_b, l_c = _estimate_window_bounds(H_func, float(ls_probe[peak]), j, k, l0, l1,
                    threshold, fd_step)
                windows.append([l_b, l_c])

    if not windows:
        return [], candidates
    windows.sort()
    merged = [windows[0]]
    for w in windows[1:]:
        if w[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], w[1])
        else:
            merged.append(w)
    return [tuple(w) for w in merged], candidates


def _local_evolution_operator(H_func: Callable, l_b: float, l_c: float, magnus_exp_order: int,
    integration_method: str, n_slabs0: Optional[int] = 400, max_n_slabs: Optional[int] = 500_000,
    patch_atol: Optional[float] = 1e-7) -> Tuple[np.ndarray, bool]:
    r"""Computes the (exact, not adiabatic) evolution operator across a single non-adiabatic
    window, via the package's own Magnus kernel, doubling the slab count until convergence.

    Uses :func:`magnus.magnus.magnus_expansion_multislab` directly (not
    :func:`magnus.oscprob.compute_evolution_operator_multiple_slabs`, to keep this module free of
    any dependency on :mod:`magnus.oscprob`).

    Returns
    -------
    (np.ndarray, bool)
        The evolution operator across ``[l_b, l_c]``, and whether it converged (agreed with the
        previous, half-as-fine slab count) within ``max_n_slabs``.
    """
    def U_at(n: int) -> np.ndarray:
        edges_lin = np.linspace(l_b, l_c, n + 1)
        edges = np.column_stack([edges_lin[:-1], edges_lin[1:]])

        def hh(t):
            return -1j * np.asarray(H_func(t))

        U_chain = magnuscore.magnus_expansion_multislab(hh, edges, n_tpts_per_slab=2,
            order=magnus_exp_order, integration_method=integration_method)
        return reduce(np.matmul, U_chain[::-1]) if len(U_chain) > 1 else U_chain[0]

    n_slabs = n_slabs0
    U_prev = U_at(n_slabs)
    while n_slabs < max_n_slabs:
        n_slabs *= 2
        U_next = U_at(n_slabs)
        if np.max(np.abs(U_next - U_prev)) <= patch_atol:
            return U_next, True
        U_prev = U_next
    return U_prev, False


def _hybrid_propagator_once(H_func: Callable, l0: float, l1: float, threshold: float,
    n_probe: int, n_points: int, fd_step_frac: float, magnus_exp_order: int,
    integration_method: str) -> Tuple[np.ndarray, List[Tuple[float, float]], bool]:
    r"""One evaluation of the hybrid propagator at a fixed set of internal tolerance knobs (see
    :func:`hybrid_propagator` for the self-certifying refinement built on top of this)."""
    windows, _ = find_nonadiabatic_windows(H_func, l0, l1, threshold=threshold, n_probe=n_probe,
        fd_step_frac=fd_step_frac)
    if not windows:
        return adiabatic_propagator(H_func, l0, l1, n_points=n_points), windows, True
    d = np.asarray(H_func(l0), dtype=complex).shape[-1]
    U_total = np.eye(d, dtype=complex)
    cursor = l0
    all_patches_converged = True
    for (l_b, l_c) in windows:
        U_total = adiabatic_propagator(H_func, cursor, l_b, n_points=n_points) @ U_total
        U_patch, ok = _local_evolution_operator(H_func, l_b, l_c, magnus_exp_order,
            integration_method)
        all_patches_converged = all_patches_converged and ok
        U_total = U_patch @ U_total
        cursor = l_c
    U_total = adiabatic_propagator(H_func, cursor, l1, n_points=n_points) @ U_total
    return U_total, windows, all_patches_converged


def hybrid_propagator(H_func: Callable, l0: float, l1: float, rtol: Optional[float] = 1.e-3,
    atol: Optional[float] = 1.e-3, magnus_exp_order: Optional[int] = 6,
    integration_method: Optional[str] = 'gl', threshold0: Optional[float] = 0.1,
    min_threshold: Optional[float] = 1.e-6, n_probe0: Optional[int] = 200,
    max_n_probe: Optional[int] = 6400, n_points0: Optional[int] = 201,
    max_n_points: Optional[int] = 12864, fd_step_frac: Optional[float] = 1.e-6,
    max_iters: Optional[int] = 12) -> Tuple[np.ndarray, List[Tuple[float, float]], bool]:
    r"""Computes the evolution operator via adiabatic transport, with a Magnus patch at every
    non-adiabatic window, self-certified against successive refinement of every internal
    tolerance knob.

    This is the main entry point of this module (see :doc:`/adiabatic_strategy` for the full
    derivation and validation). Given any Hermitian ``H_func`` of any dimension:

    #. Locates every non-adiabatic window along ``[l0, l1]`` (see
       :func:`find_nonadiabatic_windows`).
    #. If there are none, returns the pure adiabatic-transport operator (see
       :func:`adiabatic_propagator`).
    #. Otherwise, composes adiabatic transport between windows with an exact local Magnus patch
       *inside* each window (see ``_local_evolution_operator``), using the exact composition
       law of quantum evolution, :math:`U(l_2, l_0) = U(l_2, l_1)\, U(l_1, l_0)`, so the result is
       exactly unitary regardless of any approximation's accuracy.
    #. Self-certifies the result: a single fixed adiabaticity ``threshold`` (deciding which
       candidates count as non-adiabatic) is not safe in general -- too loose, and a genuine
       resonance is patched too narrowly or missed; too tight, and windows are patched
       needlessly, at some (still usually modest) extra cost. Rather than trust one fixed value,
       the whole computation (window threshold, adiabatic-transport grid density, and the probe
       grid used to *locate* candidates) is repeated with the knobs tightened together
       (threshold divided by 3, ``n_points``/``n_probe`` doubled) until two successive results
       agree within ``rtol``/``atol``, mirroring the successive-refinement discipline
       :func:`magnus.oscprob.osc_prob` already uses for the number of slabs.

       Each knob stops at its own ceiling (``min_threshold``, ``max_n_probe``,
       ``max_n_points``), which they reach at different iterations, so the later iterations
       tighten fewer knobs than the earlier ones. Once *all* of them have saturated, a further
       iteration would recompute bit-identical inputs and the agreement test would pass
       trivially, comparing a result with itself; the loop therefore stops at that point and
       reports ``certified=False`` rather than certifying on the strength of a comparison that
       carries no information.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian, a function of position returning a square (Hermitian) matrix. May be
        real- or complex-valued, of any dimension.
    l0 : float
        Initial position.
    l1 : float
        Final position.
    rtol : float, optional
        Target relative tolerance between successive refinement levels. Default: 1e-3.
    atol : float, optional
        Target absolute tolerance between successive refinement levels. Default: 1e-3.
    magnus_exp_order : int, optional
        Magnus expansion order used for the local patch inside each non-adiabatic window.
        Default: 6.
    integration_method : str, optional
        Integration method used for the local patch ('gl', 'trapezoid', or 'simpson').
        Default: 'gl'.
    threshold0 : float, optional
        Starting adiabaticity threshold. Default: 0.1.
    min_threshold : float, optional
        Floor below which the threshold is not tightened further. Default: 1e-6.
    n_probe0 : int, optional
        Starting number of positions used to locate resonance candidates. Default: 200.
    max_n_probe : int, optional
        Ceiling on the probe grid density. Default: 6400.
    n_points0 : int, optional
        Starting number of positions used for adiabatic-transport quadrature. Default: 201.
    max_n_points : int, optional
        Ceiling on the adiabatic-transport grid density. Default: 12864.
    fd_step_frac : float, optional
        Finite-difference step for the Hellmann-Feynman diagnostics, as a fraction of
        ``l1 - l0``. Default: 1e-6.
    max_iters : int, optional
        Maximum number of refinement iterations. Default: 12.

    Returns
    -------
    (np.ndarray, list of (float, float), bool)
        The evolution operator (exactly unitary regardless of ``certified``), the non-adiabatic
        windows used in the last iteration, and whether the result is certified (``True``).
        ``certified`` is ``False`` if the refinement exhausted ``max_iters``, if every knob
        reached its ceiling before two successive results agreed, or if a local patch failed to
        converge within its own slab cap -- in all three cases the returned operator is the best
        available estimate, still exactly unitary, but its accuracy is not certified to the
        requested tolerance.
    """
    threshold, n_probe, n_points = threshold0, n_probe0, n_points0
    U_prev, windows_prev, ok_prev = _hybrid_propagator_once(H_func, l0, l1, threshold, n_probe,
        n_points, fd_step_frac, magnus_exp_order, integration_method)
    if not ok_prev:
        return U_prev, windows_prev, False

    for _ in range(max_iters):
        knobs_prev = (threshold, n_probe, n_points)
        threshold = max(threshold / 3.0, min_threshold)
        n_probe = min(n_probe * 2, max_n_probe)
        n_points = min(n_points * 2, max_n_points)
        if (threshold, n_probe, n_points) == knobs_prev:
            # Every knob has hit its ceiling. _hybrid_propagator_once is deterministic, so
            # rerunning it here would reproduce U_prev exactly and the agreement test below
            # would pass on a comparison of a result with itself -- which is no evidence of
            # convergence at all. Stop and report the result as uncertified instead.
            break
        U_next, windows_next, ok_next = _hybrid_propagator_once(H_func, l0, l1, threshold,
            n_probe, n_points, fd_step_frac, magnus_exp_order, integration_method)
        if not ok_next:
            return U_next, windows_next, False
        if np.max(np.abs(U_next - U_prev)) <= atol + rtol * np.max(np.abs(U_prev)):
            return U_next, windows_next, True
        U_prev, windows_prev, ok_prev = U_next, windows_next, ok_next

    return U_prev, windows_prev, False


__all__ = [
    'adiabatic_propagator',
    'find_resonance_candidates',
    'find_nonadiabatic_windows',
    'hybrid_propagator',
]
