# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""avgprob.py

Contains the *phase-averaged* (fully decohered) oscillation
probabilities, the exact :math:`L/E \to \infty` limit reached by
high-energy astrophysical neutrinos.

Physical idea: a neutrino produced at a cosmological distance arrives
with an oscillation phase :math:`\Delta m^2 L / 2E` of order
:math:`10^{15}` or more, and neither the source distance, nor the
production region, nor the detector's energy resolution is known to
anything close to that precision.  Every oscillatory term is therefore
averaged over many cycles and vanishes, leaving only the incoherent sum

.. math::

   P(\nu_\alpha \to \nu_\beta) = \sum_i |V_{\alpha i}|^2 |V_{\beta i}|^2 ,

where :math:`V` diagonalizes the Hamiltonian.  This is not an
approximation to be refined: it is the exact limit, and it costs one
matrix product rather than an integration.  For standard vacuum
oscillations the result does not depend on energy or baseline at all, so
a single matrix serves an entire flux calculation.

Coherence is decided physically, not numerically
------------------------------------------------

The formula above assumes every *relative* phase averages away.  That is
a statement about pairs of eigenvalues, not about the spectrum as a
whole: the pair :math:`(i,j)` decoheres only if
:math:`(\lambda_i - \lambda_j) L` sweeps through many cycles across the
averaging window.  Two eigenvalues that are close enough to keep their
relative phase fixed stay *coherent*, and their cross term survives.

This module therefore groups the spectrum into blocks of mutually
coherent eigenvalues and sums coherently inside each block,

.. math::

   P(\nu_\alpha \to \nu_\beta) = \sum_{b} \Big|
   \sum_{i \in b} V^*_{\alpha i} V_{\beta i} \Big|^2 ,

which reduces to the familiar expression when every block is a singleton.
The distinction is not academic here: a sterile state with a small
:math:`\Delta m^2_{41}`, or any degenerate spectrum, makes the naive sum
quietly wrong.

The same per-pair phase decides whether the averaged limit applies at
all.  A pair whose phase spread is neither much larger than
:math:`2\pi` (decohered) nor much smaller than one (coherent) sits in
between, where no closed form is valid; :func:`coherence_report` names
those pairs, and the callers in :mod:`magnus.oscprob` warn rather than
return a number the physics does not support.

This module is self-contained: it depends only on ``numpy``, not on
:mod:`magnus.oscprob`, so it can be applied to any Hermitian Hamiltonian
of any dimension independently of the rest of the API.

Routine listings
----------------

    coherence_blocks
    coherence_report
    averaged_probabilities_from_eigenbasis
    averaged_probabilities_constant_hamiltonian
    adiabatic_phase_differences
    level_crossing_matrix
    averaged_probabilities_adiabatic
    averaged_probabilities_numerically
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from typing import Callable, List, Optional, Sequence, Tuple, Union

import numpy as np

import magnus.adiabatic as adiabatic


DECOHERENCE_PHASE_THRESHOLD = 2.0*np.pi
r"""float: Module-level constant

Accumulated phase spread, in radians, above which a pair of eigenvalues is
treated as fully decohered.  One full cycle is the point at which the average of
:math:`\cos\Delta\phi` over the window has collapsed to a small fraction of its
coherent value, and every further cycle only reduces it.

.. versionadded:: 1.0.0
"""


COHERENCE_PHASE_THRESHOLD = 1.0e-2
r"""float: Module-level constant

Accumulated phase spread, in radians, below which a pair of eigenvalues is
treated as fully coherent, so that its cross term is kept in full.

The gap between this and :data:`DECOHERENCE_PHASE_THRESHOLD` is deliberate and
is not a tolerance to be tightened away: a pair falling between the two is in
neither limit, and no averaged expression describes it.  Such pairs are reported
by :func:`coherence_report` rather than silently assigned to one side.

.. versionadded:: 1.0.0
"""


def coherence_blocks(
    eigenvalues: Union[Sequence[float], np.ndarray],
    phase_scale: float,
    decoherence_threshold: Optional[float] = DECOHERENCE_PHASE_THRESHOLD
) -> List[List[int]]:
    r"""Groups eigenvalues into blocks that stay mutually coherent.

    Two eigenvalues belong to the same block when the phase they accumulate
    relative to each other, :math:`|\lambda_i - \lambda_j| \times`
    ``phase_scale``, stays below ``decoherence_threshold``, so that their cross
    term in the probability is not averaged away.

    Grouping is by transitive closure over that relation, which is the
    conservative choice: a chain of individually-close eigenvalues is kept in one
    block rather than split at an arbitrary point.  A spectrum whose spacings are
    all comparable to the threshold therefore collapses into a single block, and
    is exactly the case :func:`coherence_report` flags as having no valid
    averaged limit.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    eigenvalues : list or np.ndarray
        Eigenvalues of the Hamiltonian [eV].  Need not be sorted.
    phase_scale : float
        Baseline over which the phase accumulates [:math:`\text{eV}^{-1}`], so
        that ``(lambda_i - lambda_j)*phase_scale`` is a phase in radians.
    decoherence_threshold : float, optional
        Phase above which a pair is treated as decohered.  Default:
        :data:`DECOHERENCE_PHASE_THRESHOLD`.

    Returns
    -------
    list of list of int
        Indices of ``eigenvalues``, grouped into blocks and sorted within each
        block.  The blocks themselves are ordered by their smallest index, so the
        result is deterministic.

    Examples
    --------
    A spectrum whose splittings are all large is fully decohered, one index per
    block; two eigenvalues sharing a value stay together.

    .. jupyter-execute::

        import magnus.avgprob as ap

        ap.coherence_blocks([0.0, 1.0, 2.0], phase_scale=1.0e3)
    """
    lam = np.asarray(eigenvalues, dtype=float).ravel()
    n = lam.size
    if n == 0:
        return []

    # Union-find over "this pair is still coherent", so the blocks are the
    # connected components of that relation rather than an order-dependent
    # sweep.
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj: parent[max(ri, rj)] = min(ri, rj)

    phases = np.abs(lam[:, None] - lam[None, :])*abs(phase_scale)
    for i in range(n):
        for j in range(i + 1, n):
            if phases[i, j] < decoherence_threshold:
                union(i, j)

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    return [sorted(g) for _, g in sorted(groups.items())]


def coherence_report(
    eigenvalues: Union[Sequence[float], np.ndarray],
    phase_scale: float,
    decoherence_threshold: Optional[float] = DECOHERENCE_PHASE_THRESHOLD,
    coherence_threshold: Optional[float] = COHERENCE_PHASE_THRESHOLD
) -> Tuple[List[List[int]], List[Tuple[int, int, float]]]:
    r"""Reports the coherence structure of a spectrum, and which pairs sit in
    neither limit.

    Every pair of eigenvalues is in one of three regimes, set by the phase it
    accumulates relative to the others over ``phase_scale``:

    * far above ``decoherence_threshold``, the cross term has averaged away and
      the pair contributes incoherently;
    * far below ``coherence_threshold``, the relative phase has barely advanced
      and the pair is still fully coherent;
    * in between, neither statement holds, and *no* averaged expression is a
      valid description -- the honest answer there is the full oscillation
      probability, not an average.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    eigenvalues : list or np.ndarray
        Eigenvalues of the Hamiltonian [eV].
    phase_scale : float
        Baseline over which the phase accumulates [:math:`\text{eV}^{-1}`].
    decoherence_threshold : float, optional
        Phase above which a pair counts as decohered.  Default:
        :data:`DECOHERENCE_PHASE_THRESHOLD`.
    coherence_threshold : float, optional
        Phase below which a pair counts as fully coherent.  Default:
        :data:`COHERENCE_PHASE_THRESHOLD`.

    Returns
    -------
    (list of list of int, list of (int, int, float))
        The coherence blocks, and the list of ``(i, j, phase)`` triples for pairs
        that are in neither limit.  An empty second element means the averaged
        result is exact for this spectrum and baseline.
    """
    lam = np.asarray(eigenvalues, dtype=float).ravel()
    blocks = coherence_blocks(lam, phase_scale, decoherence_threshold)

    undecided = []
    for i in range(lam.size):
        for j in range(i + 1, lam.size):
            phase = abs(lam[i] - lam[j])*abs(phase_scale)
            if coherence_threshold <= phase <= decoherence_threshold:
                undecided.append((i, j, float(phase)))

    return blocks, undecided


def averaged_probabilities_from_eigenbasis(
    eigenvectors: Union[Sequence, np.ndarray],
    blocks: Optional[List[List[int]]] = None
) -> np.ndarray:
    r"""Phase-averaged oscillation probabilities from the eigenbasis of the
    Hamiltonian.

    Computes

    .. math::

       P_{\alpha\beta} = \sum_b \Big| \sum_{i \in b}
       V^*_{\alpha i} V_{\beta i} \Big|^2 ,

    the sum over coherence blocks ``b`` of the squared modulus of the coherent
    amplitude within each block.  With one index per block this is the familiar
    :math:`\sum_i |V_{\alpha i}|^2 |V_{\beta i}|^2`.

    The result is symmetric, so the averaged probability is the same in both
    directions, and identical for neutrinos and antineutrinos: conjugating
    :math:`V` leaves every term unchanged.  CP violation does not survive the
    average, even though the mixing angles and phases do enter through
    :math:`|V_{\alpha i}|`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    eigenvectors : list or np.ndarray
        Matrix whose *columns* are the eigenvectors of the Hamiltonian, shape
        ``(..., d, d)``.  A leading batch axis is allowed and is broadcast over,
        so an array of energies costs one contraction.
    blocks : list of list of int, optional
        Coherence blocks, as returned by :func:`coherence_blocks`.  If None
        (default), every eigenvalue is assumed to have decohered from every
        other, which is the astrophysical case.

    Returns
    -------
    np.ndarray
        Averaged probability matrix, shape ``(..., d, d)``, with the initial
        flavor as the row index, so each row sums to one.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.avgprob as ap
        import magnus.hamiltonians as hams

        U = hams.pmns_mixing_matrix(0.55, 0.68, 0.15, 3.7)
        P = ap.averaged_probabilities_from_eigenbasis(U)
        np.round(P, 4)
    """
    V = np.asarray(eigenvectors, dtype=complex)
    if V.ndim < 2 or V.shape[-1] != V.shape[-2]:
        raise ValueError("magnus.avgprob.averaged_probabilities_from_eigenbasis: eigenvectors "
            "must be square, of shape (..., d, d), not " + str(V.shape) + "Error in magnus: .")

    d = V.shape[-1]
    if blocks is None:
        blocks = [[i] for i in range(d)]

    seen = sorted(i for b in blocks for i in b)
    if seen != list(range(d)):
        raise ValueError("magnus.avgprob.averaged_probabilities_from_eigenbasis: the blocks must "
            "partition the " + str(d) + " eigenvalue indices exactly once each; got " + str(blocks)
            + "Error in magnus: .")

    P = np.zeros(V.shape[:-2] + (d, d), dtype=float)
    for block in blocks:
        # Amplitude summed coherently inside the block:
        #     A[alpha, beta] = sum_{i in block} conj(V[alpha, i]) V[beta, i]
        V_block = V[..., :, block]
        A = np.einsum('...ai,...bi->...ab', V_block.conj(), V_block)
        P += A.real**2 + A.imag**2

    return P


def averaged_probabilities_constant_hamiltonian(
    hamiltonian: Union[Sequence, np.ndarray],
    baseline: Optional[float] = None
) -> np.ndarray:
    r"""Phase-averaged oscillation probabilities for a constant Hamiltonian.

    Diagonalizes ``hamiltonian`` and applies
    :func:`averaged_probabilities_from_eigenbasis`.  This covers every
    position-independent case -- vacuum, matter of constant density, and their
    NSI and LIV variants -- exactly, at the cost of one eigendecomposition.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    hamiltonian : list or np.ndarray
        Hermitian Hamiltonian [eV], shape ``(..., d, d)``.  A leading batch axis
        (energies, say) is allowed.
    baseline : float, optional
        Baseline [:math:`\text{eV}^{-1}`], used only to decide which eigenvalues
        have decohered from each other.  If None (default), every pair is taken
        to be decohered, which is the astrophysical limit and makes the result
        independent of distance.

    Returns
    -------
    np.ndarray
        Averaged probability matrix, shape ``(..., d, d)``, rows summing to one.
    """
    H = np.asarray(hamiltonian, dtype=complex)
    if H.ndim < 2 or H.shape[-1] != H.shape[-2]:
        raise ValueError("magnus.avgprob.averaged_probabilities_constant_hamiltonian: the "
            "Hamiltonian must be square, of shape (..., d, d), not " + str(H.shape) + "Error in magnus: .")

    eigenvalues, eigenvectors = np.linalg.eigh(H)

    if baseline is None:
        return averaged_probabilities_from_eigenbasis(eigenvectors)

    if H.ndim > 2:
        raise ValueError("magnus.avgprob.averaged_probabilities_constant_hamiltonian: a baseline "
            "can only be given for a single Hamiltonian, not for a batch of shape "
            + str(H.shape) + "Error in magnus: , since the coherence structure may differ from one to the next.")

    blocks = coherence_blocks(eigenvalues, baseline)

    return averaged_probabilities_from_eigenbasis(eigenvectors, blocks=blocks)


AVG_DEFAULT_ENERGY_SPREAD = 0.1
r"""float: Module-level constant

Half-width of the energy window, as a fraction of the energy, used when the
averaged probability has to be obtained by sampling rather than in closed form.

Ten per cent is the order of a real detector's energy resolution, and it is the
*smearing* that does the averaging: the physical statement is that the
oscillation phase varies by many cycles across whatever window the measurement
integrates over.  It is a default, not a property of the physics, so it is named
here rather than buried, every use of it is warned about, and callers with an
actual resolution should pass theirs.

.. versionadded:: 1.0.0
"""


AVG_DEFAULT_N_SAMPLES = 41
r"""int: Module-level constant

Number of samples across the window used by
:func:`averaged_probabilities_numerically`.

The sampled phases are effectively independent when the accumulated phase is
large, so the error of the mean falls only as :math:`1/\sqrt{N}` -- 41 samples
give a few per cent.  Raising it buys accuracy slowly and costs a full
propagation each; the closed-form paths in this module exist precisely to avoid
this trade.

.. versionadded:: 1.0.0
"""


def averaged_probabilities_numerically(
    prob_of_energy: Callable,
    energy: float,
    relative_spread: Optional[float] = AVG_DEFAULT_ENERGY_SPREAD,
    n_samples: Optional[int] = AVG_DEFAULT_N_SAMPLES
) -> Tuple[np.ndarray, float]:
    r"""Averages a probability by sampling it across an energy window.

    The fallback for cases with no closed form -- a profile with discontinuities, say, where
    there is no instantaneous eigenbasis to decohere in.  Unlike the closed forms in this
    module, **this is not the** :math:`L/E \to \infty` **limit**: it is the average over a
    particular window, and the answer depends on that window.  Its width is therefore an
    argument, and callers that leave it at the default should say so to their own callers.

    Samples are uniform in :math:`1/E`, in which the oscillation phase is linear, so they are
    spread evenly in phase rather than bunched.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    prob_of_energy : Callable
        
    energy : int, float or np.ndarray
        Neutrino energy [eV].
    relative_spread : float, optional
        Width of the sampling window, as a fraction of ``energy``.
    n_samples : int, optional
        Number of samples across the window.

Returns the probability matrix at a given energy; called once per sample.
    energy : float
        Central energy [eV].
    relative_spread : float, optional
        Half-width of the window as a fraction of ``energy``.  Default:
        :data:`AVG_DEFAULT_ENERGY_SPREAD`.
    n_samples : int, optional
        Number of samples.  Default: :data:`AVG_DEFAULT_N_SAMPLES`.

    Returns
    -------
    (np.ndarray, float)
        The mean probability matrix, and the largest standard error of the mean across its
        entries -- the honest uncertainty of the result, which a closed form would not have.
    """
    if not (0.0 < relative_spread < 1.0):
        raise ValueError("magnus.avgprob.averaged_probabilities_numerically: relative_spread "
            "must be between 0 and 1, not " + str(relative_spread) + "Error in magnus: .")
    if int(n_samples) < 2:
        raise ValueError("magnus.avgprob.averaged_probabilities_numerically: n_samples must be "
            "at least 2, not " + str(n_samples) + "Error in magnus: .")

    e_low = float(energy)*(1.0 - relative_spread)
    e_high = float(energy)*(1.0 + relative_spread)
    energies = 1.0/np.linspace(1.0/e_low, 1.0/e_high, int(n_samples))

    samples = np.array([np.asarray(prob_of_energy(float(e)), dtype=float) for e in energies])
    mean = samples.mean(axis=0)
    sem = float(np.max(samples.std(axis=0)/np.sqrt(len(energies))))

    return mean, sem


def adiabatic_phase_differences(
    H_func: Callable,
    l0: float,
    l1: float,
    n_points: Optional[int] = 201
) -> np.ndarray:
    r"""Relative phases accumulated between instantaneous eigenvalues.

    In the adiabatic regime a neutrino stays on one level and accumulates the dynamical phase
    :math:`\int \lambda_i(l)\, dl`, so the phase that decides whether levels :math:`i` and
    :math:`j` still interfere is
    :math:`\Delta\phi_{ij} = \int_{l_0}^{l_1} [\lambda_i(l) - \lambda_j(l)]\, dl`.
    That integral, not the eigenvalue gap at any single point, is what the coherence tests in
    this module need for a position-dependent Hamiltonian.

    Integrated with Simpson's rule: the trapezoid leaves a residual here that is easily mistaken
    for a physical effect (the same error, in the same integral, once looked like a floor on the
    accuracy of adiabatic transport in :mod:`magnus.adiabatic`).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian as a function of position, ``H_func(l)`` [eV].
    l0, l1 : float
        Start and end of the trajectory [:math:`\text{eV}^{-1}`].
    n_points : int, optional
        Number of sampling points; forced to be odd for Simpson's rule.  Default: 201.

    Returns
    -------
    np.ndarray
        Matrix of accumulated phase differences, shape ``(d, d)``, antisymmetric.
    """
    n_points = int(n_points)
    if n_points < 3: n_points = 3
    if n_points % 2 == 0: n_points += 1

    grid = np.linspace(float(l0), float(l1), n_points)
    lam = np.array([np.linalg.eigvalsh(np.asarray(H_func(l), dtype=complex)) for l in grid])

    # Simpson weights, times the uniform spacing
    h = (grid[-1] - grid[0])/(n_points - 1)
    weights = np.ones(n_points)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    integral = (h/3.0)*(weights @ lam)                          # (d,), int lambda_i dl

    return integral[:, None] - integral[None, :]


def level_crossing_matrix(
    H_func: Callable,
    l0: float,
    l1: float,
    threshold: Optional[float] = 0.1,
    n_probe: Optional[int] = 200,
    fd_step_frac: Optional[float] = 1.0e-6,
    magnus_exp_order: Optional[int] = 6,
    integration_method: Optional[str] = 'gl'
) -> Tuple[np.ndarray, List[Tuple[float, float]], bool]:
    r"""Probability of ending on level :math:`j` having started on level :math:`i`.

    Adiabatic evolution keeps a neutrino on the level it was produced on, so this matrix is the
    identity wherever the adiabatic approximation holds.  It departs from the identity only
    across a non-adiabatic window -- a resonance sharp enough for levels to exchange character
    faster than the state can follow -- and it is exactly there that the averaged probability
    needs it.

    The window is located with the Hellmann-Feynman diagnostic in :mod:`magnus.adiabatic`, and
    the transfer across it is computed with that module's own convergence-checked Magnus patch
    rather than with a Landau-Zener formula, so it inherits an exact treatment of the crossing
    instead of an asymptotic approximation to it.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian as a function of position, ``H_func(l)`` [eV].
    l0, l1 : float
        Start and end of the trajectory [:math:`\text{eV}^{-1}`].
    threshold : float, optional
        Adiabaticity threshold passed to
        :func:`magnus.adiabatic.find_nonadiabatic_windows`.  Default: 0.1.
    n_probe : int, optional
        Density of the search grid for the same.  Default: 200.
    fd_step_frac : float, optional
        Finite-difference step, as a fraction of the domain, for the same.  Default: 1e-6.
    magnus_exp_order : int, optional
        Magnus order for the local patch.  Default: 6.
    integration_method : str, optional
        Integration method for the local patch.  Default: 'gl'.

    Returns
    -------
    (np.ndarray, list of (float, float), bool)
        The level-to-level probability matrix, with the starting level as the row index; the
        non-adiabatic windows found; and whether every local patch converged.  A False in the
        last position means the crossing probabilities are not trustworthy, not that they are
        merely imprecise.
    """
    d = np.asarray(H_func(l0), dtype=complex).shape[-1]

    windows, _ = adiabatic.find_nonadiabatic_windows(H_func, float(l0), float(l1),
        threshold=threshold, n_probe=n_probe, fd_step_frac=fd_step_frac)

    crossing = np.eye(d)
    converged = True
    for (l_b, l_c) in windows:
        U_patch, ok = adiabatic._local_evolution_operator(H_func, l_b, l_c, magnus_exp_order,
            integration_method)
        converged = converged and ok

        V_b = np.linalg.eigh(np.asarray(H_func(l_b), dtype=complex))[1]
        V_c = np.linalg.eigh(np.asarray(H_func(l_c), dtype=complex))[1]

        # M[j, i] is the amplitude to arrive on level j having entered on level i, so the
        # probability matrix indexed by the starting level is the transpose of |M|^2.
        M = V_c.conj().T @ U_patch @ V_b
        crossing = crossing @ (M.real**2 + M.imag**2).T

    return crossing, windows, converged


def averaged_probabilities_adiabatic(
    H_func: Callable,
    l0: float,
    l1: float,
    n_points: Optional[int] = 201,
    threshold: Optional[float] = 0.1,
    n_probe: Optional[int] = 200,
    fd_step_frac: Optional[float] = 1.0e-6,
    magnus_exp_order: Optional[int] = 6,
    integration_method: Optional[str] = 'gl'
) -> Tuple[np.ndarray, dict]:
    r"""Phase-averaged probabilities for a position-dependent Hamiltonian.

    A neutrino produced at :math:`l_0` decoheres in the eigenbasis *there*, is carried along the
    levels of the instantaneous Hamiltonian, and is detected in the eigenbasis at :math:`l_1`:

    .. math::

       P_{\alpha\beta} = \sum_{ij} |V_{\alpha i}(l_0)|^2\, P^\text{cross}_{ij}\,
       |V_{\beta j}(l_1)|^2 ,

    with :math:`P^\text{cross}` from :func:`level_crossing_matrix` -- the identity wherever the
    evolution is adiabatic.  This is the standard MSW-plus-decoherence result, generalized to any
    number of levels and any number of crossings.

    Two things have to hold for the expression to mean anything, and both are checked rather
    than assumed.  The levels must have decohered from each other by the time of detection, and
    if there is more than one crossing they must also have decohered *between* crossings, since
    otherwise composing the crossings as probabilities -- rather than as amplitudes -- discards
    interference that is still there.  Both are reported.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    H_func : Callable
        Hamiltonian as a function of position, ``H_func(l)`` [eV].
    l0, l1 : float
        Production and detection positions [:math:`\text{eV}^{-1}`].
    n_points : int, optional
        Sampling density for the accumulated-phase integrals.  Default: 201.
    threshold, n_probe, fd_step_frac : float, int, float, optional
        Passed to :func:`level_crossing_matrix`.
    magnus_exp_order : int, optional
        Magnus order for the local patches.  Default: 6.
    integration_method : str, optional
        Integration method for the local patches.  Default: 'gl'.

    Returns
    -------
    (np.ndarray, dict)
        The averaged probability matrix, rows summing to one, and a report with keys
        ``'windows'`` (the non-adiabatic windows), ``'patches_converged'`` (bool),
        ``'undecided'`` (pairs that are in neither the coherent nor the decohered limit over the
        whole trajectory) and ``'undecided_between_crossings'`` (the same, over each adiabatic
        stretch separating two crossings).
    """
    H0 = np.asarray(H_func(l0), dtype=complex)
    H1 = np.asarray(H_func(l1), dtype=complex)

    V0 = np.linalg.eigh(H0)[1]
    V1 = np.linalg.eigh(H1)[1]

    crossing, windows, converged = level_crossing_matrix(H_func, l0, l1, threshold=threshold,
        n_probe=n_probe, fd_step_frac=fd_step_frac, magnus_exp_order=magnus_exp_order,
        integration_method=integration_method)

    W0 = V0.real**2 + V0.imag**2
    W1 = V1.real**2 + V1.imag**2
    P = W0 @ crossing @ W1.T

    # Has everything decohered by detection?
    dphi = adiabatic_phase_differences(H_func, l0, l1, n_points=n_points)
    undecided = []
    for i in range(dphi.shape[0]):
        for j in range(i + 1, dphi.shape[0]):
            phase = abs(dphi[i, j])
            if COHERENCE_PHASE_THRESHOLD <= phase <= DECOHERENCE_PHASE_THRESHOLD:
                undecided.append((i, j, float(phase)))

    # And between successive crossings, which is what composing crossings as probabilities
    # rather than as amplitudes assumes.
    undecided_between = []
    for (l_end_prev, l_start_next) in zip([w[1] for w in windows[:-1]],
                                          [w[0] for w in windows[1:]]):
        gap = adiabatic_phase_differences(H_func, l_end_prev, l_start_next, n_points=n_points)
        for i in range(gap.shape[0]):
            for j in range(i + 1, gap.shape[0]):
                phase = abs(gap[i, j])
                if phase <= DECOHERENCE_PHASE_THRESHOLD:
                    undecided_between.append((float(l_end_prev), float(l_start_next), i, j,
                                              float(phase)))

    report = {
        'windows': windows,
        'patches_converged': bool(converged),
        'undecided': undecided,
        'undecided_between_crossings': undecided_between,
    }

    return P, report


__all__ = [
    'DECOHERENCE_PHASE_THRESHOLD',
    'COHERENCE_PHASE_THRESHOLD',
    'coherence_blocks',
    'coherence_report',
    'averaged_probabilities_from_eigenbasis',
    'averaged_probabilities_constant_hamiltonian',
    'AVG_DEFAULT_ENERGY_SPREAD',
    'AVG_DEFAULT_N_SAMPLES',
    'adiabatic_phase_differences',
    'level_crossing_matrix',
    'averaged_probabilities_adiabatic',
    'averaged_probabilities_numerically',
]
