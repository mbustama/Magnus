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

This module stands apart from :mod:`magnus.oscprob`, so it can be applied to any
Hermitian Hamiltonian of any dimension independently of the rest of the API.  It
depends on ``numpy`` and on :mod:`magnus.adiabatic`; everything except
:func:`level_crossing_matrix` and :func:`averaged_probabilities_adiabatic` needs
``numpy`` alone.

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
    phase_averaged_probabilities_constant_hamiltonian
    phase_averaged_probabilities_adiabatic
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
neither limit, and no averaged expression describes it.  Such a pair is still
placed in a block -- coherent below the decoherence threshold, decohered at or
above it -- so the accompanying number is a definite choice; what
:func:`coherence_report` adds is that the choice is not made silently.

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
        that are in neither limit.  An empty second element means no pair sits
        between the two thresholds; the averaged result is then exact up to the
        residual the thresholds themselves allow, not exactly exact.
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
        raise ValueError("Error in magnus: magnus.avgprob.averaged_probabilities_from_eigenbasis: eigenvectors "
            "must be square, of shape (..., d, d), not " + str(V.shape) + ".")

    d = V.shape[-1]
    if blocks is None:
        blocks = [[i] for i in range(d)]

    seen = sorted(i for b in blocks for i in b)
    if seen != list(range(d)):
        raise ValueError("Error in magnus: magnus.avgprob.averaged_probabilities_from_eigenbasis: the blocks must "
            "partition the " + str(d) + " eigenvalue indices exactly once each; got " + str(blocks)
            + ".")

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
        independent of distance.  Only for a single Hamiltonian: giving a baseline
        for a batch raises, since the coherence structure may differ from one entry
        to the next.

    Returns
    -------
    np.ndarray
        Averaged probability matrix, shape ``(..., d, d)``, rows summing to one.
    """
    H = np.asarray(hamiltonian, dtype=complex)
    if H.ndim < 2 or H.shape[-1] != H.shape[-2]:
        raise ValueError("Error in magnus: magnus.avgprob.averaged_probabilities_constant_hamiltonian: the "
            "Hamiltonian must be square, of shape (..., d, d), not " + str(H.shape) + ".")

    eigenvalues, eigenvectors = np.linalg.eigh(H)

    if baseline is None:
        return averaged_probabilities_from_eigenbasis(eigenvectors)

    if H.ndim > 2:
        raise ValueError("Error in magnus: magnus.avgprob.averaged_probabilities_constant_hamiltonian: a baseline "
            "can only be given for a single Hamiltonian, not for a batch of shape "
            + str(H.shape) + ", since the coherence structure may differ from one to the next.")

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
here rather than buried, every use of it through the :mod:`magnus.oscprob` entry
points is warned about, and callers with an actual resolution should pass theirs.
Calling this module directly warns nobody: the width and the standard error come
back in the result instead.

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
        raise ValueError("Error in magnus: magnus.avgprob.averaged_probabilities_numerically: relative_spread "
            "must be between 0 and 1, not " + str(relative_spread) + ".")
    if int(n_samples) < 2:
        raise ValueError("Error in magnus: magnus.avgprob.averaged_probabilities_numerically: n_samples must be "
            "at least 2, not " + str(n_samples) + ".")

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
        Number of sampling points.  Truncated to an integer, raised to 3 if smaller,
        then raised to the next odd number for Simpson's rule -- all three silently.
        Default: 201.

    Returns
    -------
    np.ndarray
        Matrix of accumulated phase differences, shape ``(d, d)``, antisymmetric.
    """
    n_points = int(n_points)
    if n_points < 3: n_points = 3
    if n_points % 2 == 0: n_points += 1

    grid = np.linspace(float(l0), float(l1), n_points)
    # One vectorized Hamiltonian call and one batched eigendecomposition: the same eigenvalues,
    # bit for bit, without a Python call per grid point (issue #64).
    lam = np.linalg.eigvalsh(adiabatic._H_on_grid(H_func, grid))

    # Simpson weights, times the uniform spacing
    h = (grid[-1] - grid[0])/(n_points - 1)
    weights = np.ones(n_points)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    integral = (h/3.0)*(weights @ lam)                          # (d,), int lambda_i dl

    return integral[:, None] - integral[None, :]


SUDDEN_TRANSFER_THRESHOLD = 1.0e-3
r"""float: Module-level constant

How much probability a feature must be able to move between levels before
:func:`averaged_probabilities_adiabatic` stops trusting its own 200-point search to have seen it.

That search looks for non-adiabatic windows once, on a fixed probe grid.  A density front
narrower than the probe spacing falls between two probes and is never examined: no window
opens, :math:`P^\text{cross}` is the identity, and the answer is the fully adiabatic one,
returned without a warning (issue #60).  Such a front can be *seen* cheaply -- one half of a
probe interval carries nearly all of that interval's change -- but seeing it is not enough,
because a solar-model table interpolated in log-density shows the same shape at every one of
its grid points in the core, where nothing happens.  What separates the two is whether the
feature could move a neutrino at all.  An instantaneous change from :math:`H(l_a)` to
:math:`H(l_b)` moves at most :math:`\max_{i\ne j}|\langle v_i(l_a)|v_j(l_b)\rangle|^2`
between levels, and a monotone passage between the two positions moves less; when even that
bound is below this threshold, the feature cannot change the averaged probability by more than
the default tolerance, and today's answer stands.

Measured (``docs/dev/adversarial_batteries/sudden_transfer_sweep.py``) as the largest bound over
the intervals the probe grid finds concentrated:

=======================================================  =================  ==============================
population                                               largest bound      averaged answer today
=======================================================  =================  ==============================
BS05 solar model, cubic and linear, 1-30 MeV (core)      4.8e-07            escalating moves it <= 1.7e-16
issue #60's shock ray, the 16 of 24 fronts it got wrong  0.20 - 0.74        off by 0.08 - 0.56
the same ray, broad fronts it got right                  9.3e-03 - 0.70     right
supernova turbulence, 5-30 MeV                           0.043 - 0.93       off by up to 0.25
Earth crust with undeclared layer edges, 5-30 MeV        7.3e-06 - 4.7e-03  not scored; 1 of 9 escalates
=======================================================  =================  ==============================

1e-3, the default tolerance, sits three orders above the solar model and two below the
smallest front the engine gets wrong.  Across the roughly 5,000 averaged calls the notebooks
make, it escalates 16 pixels of paper Figure 5f and nothing else; three of those change, each
to within 0.001 of a decohered reference.

.. versionadded:: 1.1.1
"""


def _sudden_transfer(H_func: Callable, l_a: float, l_b: float) -> float:
    r"""The most probability an instantaneous change from ``H(l_a)`` to ``H(l_b)`` moves between levels.

    See :data:`SUDDEN_TRANSFER_THRESHOLD`.

    .. versionadded:: 1.1.1
    """
    V_a = np.linalg.eigh(np.asarray(H_func(l_a), dtype=complex))[1]
    V_b = np.linalg.eigh(np.asarray(H_func(l_b), dtype=complex))[1]
    M = np.abs(V_a.conj().T @ V_b)**2
    np.fill_diagonal(M, 0.0)
    return float(np.max(M))


def _unseen_features(H_func: Callable, l0: float, l1: float, n_probe: int) -> List[Tuple[float, float]]:
    r"""Probe intervals too sharp for a grid of ``n_probe`` points that could still move probability.

    See :data:`SUDDEN_TRANSFER_THRESHOLD`.

    .. versionadded:: 1.1.1
    """
    ls, flagged, _ = adiabatic._concentrated_intervals(H_func, float(l0), float(l1), n_probe)
    return [(float(ls[i]), float(ls[i + 1])) for i in flagged
            if _sudden_transfer(H_func, ls[i], ls[i + 1]) > SUDDEN_TRANSFER_THRESHOLD]


def _crossing_from_windows(
    H_func: Callable,
    d: int,
    windows: List[Tuple[float, float]],
    magnus_exp_order: int,
    integration_method: str
) -> Tuple[np.ndarray, bool]:
    r"""Level-to-level probabilities across ``windows``, each patched exactly; see :func:`level_crossing_matrix`.

    .. versionadded:: 1.1.1
    """
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
    return crossing, converged


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
    crossing, converged = _crossing_from_windows(H_func, d, windows, magnus_exp_order,
        integration_method)
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

    The windows come from a single search on ``n_probe`` points, which cannot see a front
    narrower than their spacing.  So the profile is first checked for features that sharp and
    able to move probability (see :data:`SUDDEN_TRANSFER_THRESHOLD`); where there is one, the
    windows are taken from :func:`magnus.adiabatic.hybrid_propagator` instead, which refines its
    search until it certifies, and which reports a profile it cannot resolve at all.  Everywhere
    else the result is what it was, bit for bit.

    Two things have to hold for the expression to mean anything, and both are checked rather
    than assumed.  The levels must have decohered from each other by the time of detection, and
    if there is more than one crossing they must also have decohered *between* crossings, since
    otherwise composing the crossings as probabilities -- rather than as amplitudes -- discards
    interference that is still there.  Both are reported.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.1
       Checks for features narrower than the probe spacing that could move probability, and
       takes the windows from :func:`magnus.adiabatic.hybrid_propagator` where it finds one
       (issue #60).  The report gains ``'escalated'``, ``'resolved'`` and ``'certified'``.

    Parameters
    ----------
    H_func : Callable
        Hamiltonian as a function of position, ``H_func(l)`` [eV].
    l0, l1 : float
        Production and detection positions [:math:`\text{eV}^{-1}`].
    n_points : int, optional
        Sampling density for the accumulated-phase integrals.  Default: 201.
    threshold, n_probe, fd_step_frac : float, int, float, optional
        Passed to :func:`level_crossing_matrix`. Defaults: 0.1, 200 and 1e-6.
    magnus_exp_order : int, optional
        Magnus order for the local patches.  Default: 6.
    integration_method : str, optional
        Integration method for the local patches.  Default: 'gl'.

    Returns
    -------
    (np.ndarray, dict)
        The averaged probability matrix, rows summing to one, and a report with keys
        ``'windows'`` (the non-adiabatic windows), ``'patches_converged'`` (bool),
        ``'undecided'`` (pairs that are in neither the coherent nor the decohered limit over
        the whole trajectory, as ``(i, j, phase)`` triples) and
        ``'undecided_between_crossings'`` (every pair that has *not* decohered over an
        adiabatic stretch separating two crossings, coherent pairs included, since composing
        crossings as probabilities fails for those too; entries are
        ``(l_start, l_end, i, j, phase)``).  The returned matrix is always the fully
        decohered form, so these entries qualify a number that was computed regardless --
        unlike the constant-Hamiltonian route, which keeps coherent pairs coherent.

        Three more keys say which search the windows came from.  ``'escalated'`` is True when
        the profile has a feature the ``n_probe`` grid cannot see and that could move
        probability.  Then ``'resolved'`` is whether :func:`magnus.adiabatic.hybrid_propagator`
        could resolve it -- False means a discontinuity, the matrix is the unescalated one, and
        declaring the feature through ``t_breakpoints`` is the cure -- and ``'certified'`` is
        whether that refinement certified.  Both are None when nothing escalated.
    """
    H0 = np.asarray(H_func(l0), dtype=complex)
    H1 = np.asarray(H_func(l1), dtype=complex)

    V0 = np.linalg.eigh(H0)[1]
    V1 = np.linalg.eigh(H1)[1]

    # One search on n_probe points is all the windows usually need, and it is exactly what
    # this function did before 1.1.1.  It cannot see a front narrower than the probe spacing:
    # on a supernova shock ray every such front was missed, P^cross came out the identity, and
    # the fully adiabatic answer was returned wrong by up to 0.56, silently (issue #60).  Where
    # the profile has a feature that sharp and able to move probability, take the windows from
    # the refinement the instantaneous route already certifies with -- or learn that no
    # refinement resolves it, which the caller turns into a warning.
    escalated = bool(_unseen_features(H_func, l0, l1, n_probe))
    resolved = certified = None
    windows = None
    if escalated:
        h_info = {}
        _, h_windows, certified = adiabatic.hybrid_propagator(H_func, float(l0), float(l1),
            info=h_info)
        resolved, certified = bool(h_info.get('resolved', True)), bool(certified)
        if resolved:
            windows = [tuple(w) for w in h_windows]
            crossing, converged = _crossing_from_windows(H_func, H0.shape[-1], windows,
                magnus_exp_order, integration_method)
    if windows is None:
        crossing, windows, converged = level_crossing_matrix(H_func, l0, l1,
            threshold=threshold, n_probe=n_probe, fd_step_frac=fd_step_frac,
            magnus_exp_order=magnus_exp_order, integration_method=integration_method)

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
        'escalated': escalated,
        'resolved': resolved,
        'certified': certified,
    }

    return P, report



# ------------------------------------------------------------------------------------------------
# The phase average (issue #64)
# ------------------------------------------------------------------------------------------------

AVG_PHASE_SPREAD = 0.1
r"""float: Module-level constant

Default relative energy spread :math:`\sigma` of the phase average returned by
:func:`phase_averaged_probabilities_constant_hamiltonian` and
:func:`phase_averaged_probabilities_adiabatic`, and by ``average=True`` in
:mod:`magnus.oscprob`.

The phase average keeps every interference term with its phase at the central energy and
multiplies it by :math:`e^{-\sigma^2\phi'^2/2}`, where :math:`\phi' = d\phi/d\ln E` is how fast
that phase runs with energy.  A term whose phase runs through many cycles across a spread
:math:`\sigma` is dropped, as the :math:`L/E \to \infty` limit drops it; a term whose phase barely
moves is kept with its real value; the ones in between are damped smoothly.  Ten per cent is a
typical resolution of neutrino detectors and telescopes.  Mixing, crossing amplitudes and the
eigenbases at the two ends of the path stay at the central energy: this averages phases, not
probabilities, so a result without interference is returned unchanged.

.. versionadded:: 1.1.1
"""


PHASE_AVERAGE_WINDOW_THRESHOLD = 0.01
r"""float: Module-level constant

Adiabaticity threshold at which :func:`phase_averaged_probabilities_adiabatic` looks for
non-adiabatic windows, lower than the 0.1 of :func:`averaged_probabilities_adiabatic`.

Outside a window the evolution is carried as adiabatic, so whatever small transfer between levels
happens there is lost; the decohered limit hides that loss, and the phase average does not,
because it keeps the interference such a transfer carries.  Measured against a brute-force average
of the same definition on five solar chords from 10 GeV to 10 TeV: at 0.1 the error reaches
2.6e-03, at 0.03 and at 0.01 every chord is within 4.2e-05.  On the two-level crossing of
``tests/test_phase_average.py`` 0.03 leaves 3.4e-04 and 0.01 leaves 5.1e-07, hence 0.01; the cost
falls only on calls whose phases survive the spread, since the others never reach this search.

.. versionadded:: 1.1.1
"""


PHASE_SPREAD_SENSITIVITY_THRESHOLD = 1.0e-3
r"""float: Module-level constant

Largest :math:`|\sigma\, \partial P / \partial\sigma|` a phase-averaged probability may have
before :mod:`magnus.oscprob` warns that it depends on the spread.  That derivative is the change
per e-fold of :math:`\sigma`; the threshold is the default tolerance of the package.

.. versionadded:: 1.1.1
"""


_SLOPE_FLOOR = 10.0
_PRUNE_Z = 9.0
_HERMITE_MAX = 31
_MAX_TERMS = 200_000
_PHASE_TOL = 1.0e-5


def _pair_slopes(slope_diff: np.ndarray, phase_diff: np.ndarray, scale: float,
                 dH_dlnE_step: Optional[float]) -> np.ndarray:
    r"""Pair slopes, with those below their round-off floor replaced by minus the pair phase.

    A slope is a difference of two Hellmann-Feynman derivatives, each carrying the round-off of
    :math:`dH/d\ln E`: :math:`\epsilon |H| L`, times :math:`1/h` when the derivative is a finite
    difference of step :math:`h`.  Below that floor the computed slope is noise -- a pseudo-Dirac
    pair at 100 Mpc reads hundreds of radians and would be averaged away although it is coherent
    (the case of issue #61) -- so the pair is treated as vacuum-like, where the slope is exactly
    minus the phase.
    """
    floor = _SLOPE_FLOOR*np.finfo(float).eps*scale/(dH_dlnE_step if dH_dlnE_step else 1.0)
    return np.where(np.abs(slope_diff) < floor, -phase_diff, slope_diff)


def phase_averaged_probabilities_constant_hamiltonian(
    hamiltonian: Union[Sequence, np.ndarray],
    dH_dlnE: Union[Sequence, np.ndarray],
    baseline: Union[float, np.ndarray],
    spread: Optional[float] = AVG_PHASE_SPREAD,
    dH_dlnE_step: Optional[float] = None
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Phase-averaged probabilities for a constant Hamiltonian, from a flavor state at the start.

    .. math::

       P_{\alpha\beta} = \sum_{ij} V^*_{\alpha i} V_{\beta i} V_{\alpha j} V^*_{\beta j}\,
       e^{-i\phi_{ij}}\, e^{-\sigma^2 \phi_{ij}'^2/2} ,
       \qquad \phi_{ij} = (\lambda_i - \lambda_j) L ,

    with :math:`\phi'_{ij} = d\phi_{ij}/d\ln E` from the Hellmann-Feynman derivatives
    :math:`d\lambda_i/d\ln E = \langle v_i|\, dH/d\ln E\, |v_i\rangle`.  At :math:`\sigma = 0`
    this is the oscillation probability itself; for :math:`\sigma|\phi'| \gg 1` on every pair it
    is the decohered sum :math:`\sum_i |V_{\alpha i}|^2 |V_{\beta i}|^2`.  In vacuum
    :math:`\phi' = -\phi`, so a phase of :math:`2\pi` keeps 82 per cent of its interference at
    :math:`\sigma = 10\%`, a phase of 30 rad about one per cent, and a phase of 40 rad
    :math:`3\times10^{-4}`.  See :data:`AVG_PHASE_SPREAD`.

    .. versionadded:: 1.1.1

    Parameters
    ----------
    hamiltonian : list or np.ndarray
        Hermitian Hamiltonian [eV], shape ``(..., d, d)``; a leading batch axis is allowed.
    dH_dlnE : list or np.ndarray
        Its derivative with respect to :math:`\ln E` [eV], same shape.
    baseline : float or np.ndarray
        Length of the path [:math:`\text{eV}^{-1}`], broadcast against the batch axes.
    spread : float, optional
        Relative energy spread :math:`\sigma`.  Default: :data:`AVG_PHASE_SPREAD`.
    dH_dlnE_step : float, optional
        The step in :math:`\ln E` of the finite difference ``dH_dlnE`` came from, if it came
        from one: it sets the round-off floor below which a pair's slope is taken from its
        phase instead (see issue #61).  None (default) means the derivative is exact.

    Returns
    -------
    (np.ndarray, np.ndarray)
        The probability matrix, shape ``(..., d, d)``, the initial flavor as the row index; and
        :math:`\max |\sigma\, \partial P/\partial\sigma|` over its entries, shape ``(...)``.
    """
    H = np.asarray(hamiltonian, dtype=complex)
    D = np.asarray(dH_dlnE, dtype=complex)
    if H.ndim < 2 or H.shape[-1] != H.shape[-2] or D.shape != H.shape:
        raise ValueError("Error in magnus: magnus.avgprob.phase_averaged_probabilities_constant_hamiltonian: "
            "the Hamiltonian must be square, of shape (..., d, d), and its derivative the same shape; "
            "got " + str(H.shape) + " and " + str(D.shape) + ".")
    if spread is None or spread < 0.0:
        raise ValueError("Error in magnus: magnus.avgprob.phase_averaged_probabilities_constant_hamiltonian: "
            "the spread must be a non-negative number, not " + repr(spread) + ".")
    lam, V = np.linalg.eigh(H)
    slope = np.real(np.einsum('...ai,...ab,...bi->...i', V.conj(), D, V))
    L = np.asarray(baseline, dtype=float)[..., None, None]
    phi = (lam[..., :, None] - lam[..., None, :])*L
    scale = np.max(np.abs(lam), axis=-1)[..., None, None]*np.abs(L)
    dphi = _pair_slopes((slope[..., :, None] - slope[..., None, :])*L, phi, scale, dH_dlnE_step)
    x2 = (spread*dphi)**2
    w = np.exp(-0.5*x2)
    rot = np.exp(-1j*phi)
    # P_ab = Re sum_ij X_abi conj(X_abj) K_ij, X_abi = conj(V_ai) V_bi, as one batched product
    # over j for the weights and their sigma-derivative together.
    X = V.conj()[..., :, None, :]*V[..., None, :, :]
    K = np.stack([rot*w, rot*(-x2*w)], axis=-3)                       # (..., 2, i, j)
    Xc = X.conj()
    Kt = np.swapaxes(K, -1, -2)                                         # (..., 2, j, i)
    Y = Xc[..., None, :, :, :] @ Kt[..., :, None, :, :]                # (..., 2, a, b, i)
    PS = np.real(np.sum(X[..., None, :, :, :]*Y, axis=-1))             # (..., 2, a, b)
    P, S = PS[..., 0, :, :], PS[..., 1, :, :]
    return P, np.max(np.abs(S), axis=(-2, -1))


def _stretch_once(H_func: Callable, D_func: Callable, a: float, z: float, n: int):
    """Dynamical phase, its slope, and the transported eigenvectors at both ends, on n points."""
    xs = np.linspace(a, z, n)
    xs[0], xs[-1] = a, z
    Hs = adiabatic._H_on_grid(H_func, xs)
    Ds = adiabatic._H_on_grid(D_func, xs)
    lam, V = np.linalg.eigh(Hs)
    sl = np.real(np.einsum('nai,nab,nbi->ni', V.conj(), Ds, V))
    h = (z - a)/(n - 1)
    w = np.ones(n)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    ov = np.einsum('nai,nai->ni', V[:-1].conj(), V[1:])
    transport = -np.sum(np.angle(ov), axis=0)
    return ((h/3.0)*(w @ lam), (h/3.0)*(w @ sl), transport, V[0], V[-1],
            float(np.max(np.abs(lam)))*abs(z - a), float(np.min(np.abs(ov))))


def _stretch(H_func: Callable, D_func: Callable, a: float, z: float, V_start: np.ndarray,
             V_end: np.ndarray, spread: float, n0: int = 801, n_max: int = 102_401) -> dict:
    r"""Adiabatic transport from ``a`` to ``z`` as per-level phases and slopes in ``ln E``.

    The phase of level :math:`i` carries the dynamical phase :math:`\int\lambda_i`, the
    parallel-transport phase along the grid, and the phase differences between the grid's
    eigenvectors at the two ends and the bases ``V_start``, ``V_end`` the neighbouring windows
    and the readout use -- so the transport composes with them whatever phase ``eigh`` gave each
    eigenvector.  Simpson's rule on a grid doubled until the pair phases that can still matter
    (weight above :math:`10^{-12}`) move by less than ``_PHASE_TOL`` = 1e-5 rad, which moves a
    probability by at most as much.  A tabulated profile, interpolated with kinks at its rows,
    converges slowly: measured on a solar chord at 100 GeV, phases of 3e3 rad move by 1e-3 rad
    between 801 and 1601 points and by 1e-6 between 25 601 and 51 201.
    """
    d = V_start.shape[0]
    if z <= a:
        # Zero length: the only transport is the change of basis, which must be diagonal.
        ph = np.angle(np.einsum('ai,ai->i', V_end.conj(), V_start))
        return dict(phase=-ph, slope=np.zeros(d), scale=0.0, converged=True, min_overlap=1.0, n=1)
    n, prev, converged = n0, None, False
    while True:
        Phi, dPhi, transport, V0, V1, scale, min_ov = _stretch_once(H_func, D_func, a, z, n)
        match = (np.angle(np.einsum('ai,ai->i', V0.conj(), V_start))
                 + np.angle(np.einsum('ai,ai->i', V_end.conj(), V1)))
        phase = Phi - transport - match
        if prev is not None:
            dp = (phase - prev[0])
            ds = (dPhi - prev[1])
            wgt = np.exp(-0.5*(spread*(dPhi[:, None] - dPhi[None, :]))**2) > 1.0e-12
            change = max(float(np.max(np.abs(dp[:, None] - dp[None, :])[wgt], initial=0.0)),
                         spread*float(np.max(np.abs(ds[:, None] - ds[None, :]), initial=0.0)))
            if change < _PHASE_TOL:
                converged = True
                break
        if 2*n - 1 > n_max:
            break
        prev = (phase, dPhi)
        n = 2*n - 1
    return dict(phase=phase, slope=dPhi, scale=scale, converged=converged, min_overlap=min_ov, n=n)


def _window_amplitudes(H_func: Callable, D_func: Callable, l_b: float, l_c: float,
                       u_nodes: np.ndarray, V_b: np.ndarray, V_c: np.ndarray, magnus_exp_order: int,
                       integration_method: str, n_slabs0: int = 400, max_n_slabs: int = 32_768,
                       patch_atol: float = 1.0e-7) -> Tuple[np.ndarray, bool]:
    r"""The amplitude matrix :math:`V(l_c)^\dagger U_u V(l_b)` across a window, at every node.

    :math:`U_u` evolves with :math:`H + u\,D_\text{diag}`, where :math:`D_\text{diag}` is the part
    of :math:`dH/d\ln E` diagonal in the instantaneous eigenbasis: an energy offset :math:`u`
    moves every eigenvalue by :math:`u\, d\lambda_i/d\ln E` and leaves the eigenvectors alone,
    which is the definition of the phase average carried inside the window.  The node
    :math:`u = 0` is the patch of :func:`level_crossing_matrix` itself; the others share one slab
    count, converged at the largest :math:`|u|`, and one evaluation of the Hamiltonian, its
    derivative and its eigenbasis per quadrature position.
    """
    d = V_b.shape[0]
    cache = {}

    def HD(t):
        key = (np.shape(t), np.asarray(t, dtype=float).tobytes())
        if key not in cache:
            ta = np.atleast_1d(np.asarray(t, dtype=float))
            H = adiabatic._H_on_grid(H_func, ta)
            D = adiabatic._H_on_grid(D_func, ta)
            _, V = np.linalg.eigh(H)
            sl = np.real(np.einsum('nai,nab,nbi->ni', V.conj(), D, V))
            Dd = (V*sl[:, None, :]) @ V.conj().transpose(0, 2, 1)
            cache[key] = (H[0], Dd[0]) if np.ndim(t) == 0 else (H, Dd)
        return cache[key]

    def U_at(u, n):
        def A(t):
            H, Dd = HD(t)
            return -1j*(H + u*Dd)
        e = np.linspace(l_b, l_c, n + 1)
        chain = adiabatic.magnuscore.magnus_expansion_multislab(A, np.column_stack([e[:-1], e[1:]]),
            n_tpts_per_slab=2, order=magnus_exp_order, integration_method=integration_method,
            A_eval_mode='vector')
        return adiabatic.magnuscore.ordered_product(chain)

    u_max = float(np.max(np.abs(u_nodes))) if len(u_nodes) else 0.0
    n, converged = n_slabs0, True
    if u_max > 0.0:
        converged = False
        prev = U_at(u_max, n)
        while n < max_n_slabs:
            n *= 2
            nxt = U_at(u_max, n)
            if np.max(np.abs(nxt - prev)) <= patch_atol:
                converged = True
                break
            prev = nxt
        cache.clear()
    M = np.empty((len(u_nodes), d, d), dtype=complex)
    for k, u in enumerate(u_nodes):
        if u == 0.0:
            U, ok = adiabatic._local_evolution_operator(H_func, l_b, l_c, magnus_exp_order,
                integration_method)
            converged = converged and ok
        else:
            U = U_at(float(u), n)
        M[k] = V_c.conj().T @ U @ V_b
    return M, converged


def _hermite_order(x: float, tol: float = 1.0e-9) -> int:
    r"""Odd Gauss-Hermite order whose truncation of :math:`e^{-i\kappa u}`, :math:`\sigma|\kappa| \le x`, is below ``tol``."""
    if x < 1.0e-12:
        return 1
    K = 3
    while K <= _HERMITE_MAX:
        # |c_m| <~ x^m/m!, times max_t t^m e^{-t^2/2} = (m/e)^(m/2)
        bound = np.exp(K*np.log(x) + 0.5*K*(np.log(K) - 1.0) - np.sum(np.log(np.arange(1, K + 1))))
        if bound < tol:
            return K
        K += 2
    return K


def phase_averaged_probabilities_adiabatic(
    H_func: Callable,
    dH_dlnE_func: Callable,
    l0: float,
    l1: float,
    spread: Optional[float] = AVG_PHASE_SPREAD,
    threshold: Optional[float] = PHASE_AVERAGE_WINDOW_THRESHOLD,
    windows: Optional[List[Tuple[float, float]]] = None,
    n_probe: Optional[int] = 200,
    fd_step_frac: Optional[float] = 1.0e-6,
    magnus_exp_order: Optional[int] = 6,
    integration_method: Optional[str] = 'gl',
    dH_dlnE_step: Optional[float] = None
) -> Tuple[np.ndarray, dict]:
    r"""Phase-averaged probabilities on a smooth position-dependent Hamiltonian.

    The neutrino starts decohered in the eigenbasis at :math:`l_0`, as in
    :func:`averaged_probabilities_adiabatic`, and is read out in the flavor basis at :math:`l_1`.
    In between, every interference term is kept with its phase and weighted by the spread of
    that phase across a relative energy spread :math:`\sigma` (see :data:`AVG_PHASE_SPREAD`):
    formally, the Gaussian average over :math:`u = \delta\ln E` of the evolution under
    :math:`H + u\,D_\text{diag}`, with :math:`D_\text{diag}` the part of :math:`dH/d\ln E`
    diagonal in the instantaneous eigenbasis -- an energy offset moves the eigenvalues and
    leaves the eigenvectors.

    It is computed without sampling energies across the adiabatic stretches.  Each non-adiabatic
    window is an amplitude matrix, evaluated at a few Gauss-Hermite nodes in :math:`u` (a uniform
    grid when its internal phase runs too fast for that); each stretch between windows is a
    diagonal phase with an exact slope in :math:`u`.  The density matrix is carried as terms
    labelled by accumulated slope, whose Gaussian average is analytic, and a term is dropped only
    once no later stretch can bring its slope back within reach.  So the answer does not depend on
    where the windows are drawn: one window over a stretch or two windows with the stretch between
    them give the same number.

    Where there is no window the evolution is adiabatic, a decohered start carries no
    interference, and the result is the decohered expression of
    :func:`averaged_probabilities_adiabatic`.

    .. versionadded:: 1.1.1

    Parameters
    ----------
    H_func : Callable
        Hamiltonian as a function of position, ``H_func(l)`` [eV]; arrays of positions are used
        where it accepts them.
    dH_dlnE_func : Callable
        Its derivative with respect to :math:`\ln E`, as a function of position [eV].
    l0, l1 : float
        Production and detection positions [:math:`\text{eV}^{-1}`].
    spread : float, optional
        Relative energy spread :math:`\sigma`.  Default: :data:`AVG_PHASE_SPREAD`.
    threshold : float, optional
        Adiabaticity threshold of the window search.  Default:
        :data:`PHASE_AVERAGE_WINDOW_THRESHOLD`.
    windows : list of (float, float), optional
        Windows to use instead of searching.  Default: None.
    n_probe, fd_step_frac : int, float, optional
        Passed to the window search.  Defaults: 200 and 1e-6.
    magnus_exp_order : int, optional
        Magnus order of the window patches.  Default: 6.
    integration_method : str, optional
        Integration method of the window patches.  Default: 'gl'.
    dH_dlnE_step : float, optional
        See :func:`phase_averaged_probabilities_constant_hamiltonian`.  Default: None.

    Returns
    -------
    (np.ndarray, dict)
        The probability matrix, initial flavor as the row index; and a report with keys
        ``'windows'``, ``'escalated'``, ``'resolved'``, ``'certified'`` (as in
        :func:`averaged_probabilities_adiabatic`), ``'patches_converged'``,
        ``'phases_converged'``, ``'n_nodes'`` and ``'method'`` (``'hermite'``, ``'grid'``, or
        ``'none'`` without windows), ``'n_terms'``, and ``'sigma_sensitivity'``, the largest
        :math:`|\sigma\, \partial P/\partial\sigma|`.

    Raises
    ------
    RuntimeError
        If the number of terms would exceed an internal bound (``_MAX_TERMS``): many windows
        with many flavors whose phases never decohere.
    """
    if spread is None or spread < 0.0:
        raise ValueError("Error in magnus: magnus.avgprob.phase_averaged_probabilities_adiabatic: the "
            "spread must be a non-negative number, not " + repr(spread) + ".")
    l0, l1 = float(l0), float(l1)
    _, V0 = np.linalg.eigh(np.asarray(H_func(l0), dtype=complex))
    _, V1 = np.linalg.eigh(np.asarray(H_func(l1), dtype=complex))
    d = V0.shape[0]
    W0 = V0.real**2 + V0.imag**2
    W1 = V1.real**2 + V1.imag**2
    report = dict(escalated=False, resolved=None, certified=None)

    if windows is None:
        escalated = bool(_unseen_features(H_func, l0, l1, n_probe))
        report['escalated'] = escalated
        if escalated:
            h_info = {}
            _, h_windows, certified = adiabatic.hybrid_propagator(H_func, l0, l1, info=h_info)
            report['resolved'] = bool(h_info.get('resolved', True))
            report['certified'] = bool(certified)
            if report['resolved']:
                windows = [tuple(w) for w in h_windows]
        if windows is None:
            windows, _ = adiabatic.find_nonadiabatic_windows(H_func, l0, l1, threshold=threshold,
                n_probe=n_probe, fd_step_frac=fd_step_frac)
    windows = [(float(b), float(c)) for b, c in windows]
    report['windows'] = windows

    if not windows:
        report.update(patches_converged=True, phases_converged=True, n_nodes=0, method='none',
                      n_terms=1, sigma_sensitivity=0.0)
        return W0 @ W1.T, report

    # Eigenbases at every window edge: the windows' amplitudes and the stretches' transport
    # are both expressed in these, so they compose whatever phase eigh gave each vector.
    V_b = [np.linalg.eigh(np.asarray(H_func(b), dtype=complex))[1] for b, _ in windows]
    V_c = [np.linalg.eigh(np.asarray(H_func(c), dtype=complex))[1] for _, c in windows]
    ends = [b for b, _ in windows][1:] + [l1]
    V_ends = V_b[1:] + [V1]
    stretches = [_stretch(H_func, dH_dlnE_func, c, z, V_c[i], V_ends[i], spread)
                 for i, ((_, c), z) in enumerate(zip(windows, ends))]
    report['phases_converged'] = all(s['converged'] for s in stretches)

    # How fast the phases inside each window run with u: this sets the nodes.
    drift = []
    for (b, c) in windows:
        _, dPhi, _, _, _, _, _ = _stretch_once(H_func, dH_dlnE_func, b, c, 801)
        drift.append(float(np.max(dPhi) - np.min(dPhi)))
    D_W = float(sum(drift))
    K = _hermite_order(spread*D_W)
    if K <= _HERMITE_MAX:
        x, wq = np.polynomial.hermite_e.hermegauss(K)
        wq = wq/wq.sum()
        u = spread*x
        method = 'hermite'
    else:
        # Frequencies up to twice the windows' drift plus what survives pruning, alias-free
        du = 2*np.pi/(2*D_W + _PRUNE_Z/spread + 12.0/spread)
        m = int(np.ceil(6.0*spread/du))
        u = du*np.arange(-m, m + 1)
        wq = np.exp(-0.5*(u/spread)**2)
        wq = wq/wq.sum()
        method = 'grid'
    report.update(n_nodes=len(u), method=method)

    Ms, conv = [], True
    for (b, c), Vb, Vc in zip(windows, V_b, V_c):
        M, ok = _window_amplitudes(H_func, dH_dlnE_func, b, c, u, Vb, Vc, magnus_exp_order,
                                   integration_method)
        Ms.append(M)
        conv = conv and ok
    report['patches_converged'] = conv

    # Propagate rho(u) = sum_t R_t(u) exp(-i s_t u), R_t at every node, in the level bases.
    span = [float(np.max(s['slope']) - np.min(s['slope'])) for s in stretches]
    reach_after = [D_W + sum(span[j] for j in range(i + 1, len(windows)))
                   for i in range(len(windows))]
    R0 = np.zeros((len(u), d, d, d), dtype=complex)
    for a in range(d):
        R0[:, a] = np.diag(W0[a])[None]
    terms = {0: (0.0, R0)}
    for i in range(len(windows)):
        M = Ms[i]
        terms = {key: (s, np.einsum('kij,kajm,klm->kail', M, R, M.conj()))
                 for key, (s, R) in terms.items()}
        st = stretches[i]
        ph = np.exp(-1j*st['phase'])
        rot = ph[:, None]*ph.conj()[None, :]
        dsl = _pair_slopes(st['slope'][:, None] - st['slope'][None, :],
                           st['phase'][:, None] - st['phase'][None, :], st['scale'], dH_dlnE_step)
        new = {}
        for key, (s, R) in terms.items():
            Rr = R*rot[None, None]
            for p in range(d):
                for q in range(d):
                    snew = s + (0.0 if p == q else float(dsl[p, q]))
                    if spread*(abs(snew) - reach_after[i]) > _PRUNE_Z:
                        continue
                    kk = int(round(snew*spread*1.0e9)) if spread > 0.0 else 0
                    if kk not in new:
                        if len(new) >= _MAX_TERMS:
                            raise RuntimeError("Error in magnus: magnus.avgprob."
                                "phase_averaged_probabilities_adiabatic: more than "
                                + str(_MAX_TERMS) + " interference terms survive across "
                                + str(len(windows)) + " windows at " + str(d) + " flavors.")
                        new[kk] = (snew, np.zeros_like(R))
                    new[kk][1][:, :, p, q] += Rr[:, :, p, q]
        terms = new
    report['n_terms'] = len(terms)

    P = np.zeros((d, d), dtype=complex)
    dP = np.zeros((d, d), dtype=complex)
    if method == 'hermite':
        He = np.polynomial.hermite_e.hermevander(x, K - 1)
        fact = np.cumprod(np.r_[1.0, np.arange(1, K)])
        mm = np.arange(K)
    for s, R in terms.values():
        F = np.einsum('bi,kaij,bj->kab', V1, R, V1.conj())
        if method == 'hermite':
            # R(u) = sum_m c_m He_m(u/sigma); E[He_m(x) e^{-itx}] = (-it)^m e^{-t^2/2}, and
            # sigma d/dsigma brings in (x^2 - 1) = He_2, with He_2 He_m = He_{m+2} + 2m He_m
            # + m(m-1) He_{m-2}.
            c = np.einsum('k,km,kab->mab', wq, He, F)/fact[:, None, None]
            t = spread*s
            z = -1j*t
            g = np.exp(-0.5*t*t)*z**mm
            zm2 = np.where(mm >= 2, z**np.maximum(mm - 2, 0), 0.0)
            g2 = np.exp(-0.5*t*t)*(z**(mm + 2) + 2*mm*z**mm + mm*(mm - 1)*zm2)
            P += np.einsum('m,mab->ab', g, c)
            dP += np.einsum('m,mab->ab', g2, c)
        else:
            gk = np.exp(-1j*s*u)[:, None, None]*F
            P += np.einsum('k,kab->ab', wq, gk)
            dP += np.einsum('k,kab->ab', wq*((u/spread)**2 - 1.0), gk)
    report['sigma_sensitivity'] = float(np.max(np.abs(dP.real)))
    return P.real, report

__all__ = [
    'DECOHERENCE_PHASE_THRESHOLD',
    'COHERENCE_PHASE_THRESHOLD',
    'coherence_blocks',
    'coherence_report',
    'averaged_probabilities_from_eigenbasis',
    'averaged_probabilities_constant_hamiltonian',
    'SUDDEN_TRANSFER_THRESHOLD',
    'AVG_DEFAULT_ENERGY_SPREAD',
    'AVG_DEFAULT_N_SAMPLES',
    'adiabatic_phase_differences',
    'level_crossing_matrix',
    'averaged_probabilities_adiabatic',
    'averaged_probabilities_numerically',
    'AVG_PHASE_SPREAD',
    'PHASE_AVERAGE_WINDOW_THRESHOLD',
    'PHASE_SPREAD_SENSITIVITY_THRESHOLD',
    'phase_averaged_probabilities_constant_hamiltonian',
    'phase_averaged_probabilities_adiabatic',
]
