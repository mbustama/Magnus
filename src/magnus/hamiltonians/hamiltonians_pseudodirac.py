# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""hamiltonians_pseudodirac.py

Compute pseudo-Dirac neutrino Hamiltonians for an arbitrary pairing pattern.

A Dirac neutrino may in fact be two Majorana states separated by a tiny
mass-squared splitting.  Each mass eigenstate that carries such a partner
becomes two eigenstates, :math:`(|\nu_j\rangle \pm |s_j\rangle)/\sqrt{2}`, with
masses :math:`m_j^2` and :math:`m_j^2 + \delta m_j^2`; each eigenstate without a
partner stays single.  The pairing is chosen **per mass state**, so a
three-flavor spectrum with partners on two of its three states is a
five-dimensional problem, and is supported.

The physical content is a separation of scales.  The pair splittings are far
smaller than the standard ones, so the standard phases average away over an
astrophysical baseline while each pair stays mutually coherent -- which is
exactly the regime the coherent-block averaging form is for, and the regime in
which the naive sum over eigenstates is wrong.  See
:func:`magnus.avgprob.coherence_blocks` and
:func:`magnus.avgprob.coherence_report`.

**No partial averaging is provided here, deliberately.**  This package's
contract is that a fully coherent pair is handled by the block form, a fully
decohered pair by the ordinary sum, and anything in between is refused rather
than approximated.  These routines build the Hamiltonian; the existing engines
propagate it.

Routine listings
----------------

    * PseudoDiracSplittingWarning - Raised when a splitting is not small
    * pseudo_dirac_mixing_matrix - Returns the extended mixing matrix
    * pseudo_dirac_mass_squared - Returns the extended mass-squared values
    * hamiltonian_pseudo_dirac_vacuum_energy_independent - Returns H_vac (no
           1/E)
    * hamiltonian_pseudo_dirac_vacuum - Returns H_vac
    * hamiltonian_pseudo_dirac_matter - Returns H_matter
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import warnings
from typing import Mapping, Optional, Sequence, Tuple, Union

import numpy as np

import magnus.matter as matter


__all__ = [
    'PseudoDiracSplittingWarning',
    'pseudo_dirac_mixing_matrix',
    'pseudo_dirac_mass_squared',
    'hamiltonian_pseudo_dirac_vacuum_energy_independent',
    'hamiltonian_pseudo_dirac_vacuum',
    'hamiltonian_pseudo_dirac_matter',
]


class PseudoDiracSplittingWarning(UserWarning):
    r"""A pseudo-Dirac splitting is not small against the standard ones.

    The physics of the pseudo-Dirac case rests on a separation of scales: the
    pair splittings must be far below the standard mass-squared differences, so
    that the standard phases decohere while the pairs stay coherent.  A
    splitting comparable with :math:`\Delta m^2_{21}` describes a different
    system -- effectively a sterile state with an ordinary splitting -- for
    which the four- and five-flavor routines are the appropriate tools.

    The calculation still proceeds; nothing is clamped.

    .. versionadded:: 1.0.5
    """


def _pair_layout(n_active: int, pairs: Mapping[int, float]) -> Tuple[list, list]:
    r"""Validates ``pairs`` and returns the eigenstate and sterile-row layout.

    Returns ``(columns, sterile_of)``.  ``columns`` lists, for each output
    eigenstate in order, the pair ``(j, sign)`` it comes from: ``sign`` is
    ``0`` for an unpaired state, ``+1`` for the symmetric combination and
    ``-1`` for the antisymmetric one.  ``sterile_of`` lists the paired mass
    states in the order their sterile rows are appended.
    """
    if pairs is None:
        pairs = {}
    if not isinstance(pairs, Mapping):
        raise TypeError(
            "Error in magnus: hamiltonians_pseudodirac: `pairs` must be a "
            "mapping from mass-state index to splitting, e.g. {0: 1.0e-18}; "
            "got %s" % type(pairs).__name__)

    for index, splitting in pairs.items():
        if not isinstance(index, (int, np.integer)) or isinstance(index, bool):
            raise TypeError(
                "Error in magnus: hamiltonians_pseudodirac: pairing keys must "
                "be integer mass-state indices; got %r" % (index,))
        if not 0 <= int(index) < n_active:
            raise ValueError(
                "Error in magnus: hamiltonians_pseudodirac: pairing index %d "
                "is outside the range of mass states [0, %d)"
                % (int(index), n_active))
        value = float(splitting)
        if not np.isfinite(value):
            raise ValueError(
                "Error in magnus: hamiltonians_pseudodirac: the splitting for "
                "mass state %d is not finite" % int(index))
        if value <= 0.0:
            raise ValueError(
                "Error in magnus: hamiltonians_pseudodirac: the splitting for "
                "mass state %d must be positive; got %g.  A pair with zero "
                "splitting is not a pseudo-Dirac pair -- omit the state from "
                "`pairs` to leave it unpaired." % (int(index), value))

    columns, sterile_of = [], []
    for j in range(n_active):
        if j in pairs:
            columns.append((j, +1))
            columns.append((j, -1))
            sterile_of.append(j)
        else:
            columns.append((j, 0))
    return columns, sterile_of


def pseudo_dirac_mixing_matrix(
    mixing_matrix: Union[Sequence, np.ndarray],
    pairs: Optional[Mapping[int, float]] = None
) -> np.ndarray:
    r"""Returns the extended mixing matrix for a pseudo-Dirac spectrum.

    Each paired mass eigenstate :math:`j` is replaced by the two combinations
    :math:`(|\nu_j\rangle \pm |s_j\rangle)/\sqrt{2}`, which contribute two
    columns; each unpaired state contributes one.  Rows are ordered as the
    ``n_active`` active flavors first, then one sterile partner per paired mass
    state, in ascending order of that state's index.

    The result is unitary whenever the input is, for any pairing pattern.

    .. versionadded:: 1.0.5

    Parameters
    ----------
    mixing_matrix : list or np.ndarray
        The ``n_active`` x ``n_active`` mixing matrix of the active sector,
        e.g. the output of :func:`magnus.hamiltonians.hamiltonians3nu.pmns_mixing_matrix`.
    pairs : dict, optional
        Mapping from mass-state index to its pseudo-Dirac splitting
        :math:`\delta m^2_j`, in eV^2.  States absent from the mapping are
        unpaired.  An empty mapping (the default) returns the input unchanged,
        so the Dirac case is recovered exactly.

    Returns
    -------
    np.ndarray
        The ``n`` x ``n`` complex mixing matrix, with
        ``n = n_active + len(pairs)``.

    Raises
    ------
    ValueError
        If ``mixing_matrix`` is not square, or a pairing index is out of range,
        or a splitting is not positive and finite.

    See Also
    --------
    pseudo_dirac_mass_squared : the matching mass-squared values.

    Examples
    --------
    Three active flavors with partners on the first and third mass states is a
    five-dimensional problem, and the extended matrix is still unitary.

    .. jupyter-execute::

        import numpy as np
        import magnus.hamiltonians as hamiltonians

        U = hamiltonians.pmns_mixing_matrix(0.5558, 0.6856, 0.1499, 3.7001)
        W = hamiltonians.pseudo_dirac_mixing_matrix(U, {0: 1.0e-18, 2: 4.0e-18})

        print(W.shape, bool(np.allclose(W @ W.conj().T, np.eye(len(W)))))
    """
    U = np.asarray(mixing_matrix, dtype=complex)
    if U.ndim != 2 or U.shape[0] != U.shape[1]:
        raise ValueError(
            "Error in magnus: hamiltonians_pseudodirac: `mixing_matrix` must "
            "be square; got shape %s" % (U.shape,))
    n_active = U.shape[0]
    columns, sterile_of = _pair_layout(n_active, pairs)

    n = n_active + len(sterile_of)
    W = np.zeros((n, n), dtype=complex)
    root_half = 1.0/np.sqrt(2.0)
    for column, (j, sign) in enumerate(columns):
        if sign == 0:
            W[:n_active, column] = U[:, j]
        else:
            W[:n_active, column] = U[:, j]*root_half
            W[n_active + sterile_of.index(j), column] = sign*root_half
    return W


def pseudo_dirac_mass_squared(
    mass_squared: Union[Sequence[float], np.ndarray],
    pairs: Optional[Mapping[int, float]] = None
) -> np.ndarray:
    r"""Returns the extended mass-squared values for a pseudo-Dirac spectrum.

    A paired state :math:`j` contributes :math:`m_j^2` and
    :math:`m_j^2 + \delta m_j^2`, in that order; an unpaired state contributes
    :math:`m_j^2` alone.  The ordering matches
    :func:`pseudo_dirac_mixing_matrix` column for column.

    .. versionadded:: 1.0.5

    Parameters
    ----------
    mass_squared : list or np.ndarray
        The ``n_active`` mass-squared values of the active sector, in eV^2.
        Only differences matter, so these are usually
        ``[0, Dm21, Dm31]``.
    pairs : dict, optional
        Mapping from mass-state index to its splitting, as in
        :func:`pseudo_dirac_mixing_matrix`.

    Returns
    -------
    np.ndarray
        The ``n`` mass-squared values, ``n = n_active + len(pairs)``.

    Warns
    -----
    PseudoDiracSplittingWarning
        If a splitting is not small against the standard splittings.  The
        pseudo-Dirac regime rests on that separation of scales; a comparable
        splitting is a sterile state with an ordinary mass, better described by
        the four- or five-flavor routines.

    Examples
    --------
    .. jupyter-execute::

        import magnus.hamiltonians as hamiltonians

        hamiltonians.pseudo_dirac_mass_squared([0.0, 7.5e-5, 2.511e-3],
                                               {0: 1.0e-18, 2: 4.0e-18})
    """
    m2 = np.asarray(mass_squared, dtype=float).ravel()
    n_active = m2.size
    columns, _ = _pair_layout(n_active, pairs)
    pairs = {} if pairs is None else pairs

    if pairs:
        spread = float(np.max(m2) - np.min(m2))
        if spread > 0.0:
            largest = max(float(v) for v in pairs.values())
            if largest > 1.0e-3*spread:
                warnings.warn(
                    "Pseudo-Dirac splitting %.3g eV^2 is not small against the "
                    "standard mass-squared splittings (spread %.3g eV^2).  The "
                    "pseudo-Dirac regime assumes the pair phase is still "
                    "developing while the standard phases have averaged away; "
                    "at this size the two scales overlap and the spectrum is "
                    "better described as a sterile state with an ordinary "
                    "splitting, for which the four- and five-flavor routines "
                    "are appropriate.  Proceeding unchanged."
                    % (largest, spread), PseudoDiracSplittingWarning,
                    stacklevel=2)

    out = np.empty(len(columns), dtype=float)
    for column, (j, sign) in enumerate(columns):
        out[column] = m2[j] + (float(pairs[j]) if sign == -1 else 0.0)
    return out


def hamiltonian_pseudo_dirac_vacuum_energy_independent(
    mixing_matrix: Union[Sequence, np.ndarray],
    mass_squared: Union[Sequence[float], np.ndarray],
    pairs: Optional[Mapping[int, float]] = None,
    nubar: Optional[bool] = False
) -> np.ndarray:
    r"""Returns the pseudo-Dirac vacuum Hamiltonian, without the 1/E factor.

    :math:`H_{\rm vac} E = \tfrac{1}{2} W M^2 W^\dagger`, with :math:`W` the
    extended mixing matrix and :math:`M^2` the extended mass-squared values.
    Because the energy factors out, this is what the energy-batched engine
    reuses across a scan; see :func:`hamiltonian_pseudo_dirac_vacuum`.

    .. versionadded:: 1.0.5

    Parameters
    ----------
    mixing_matrix : list or np.ndarray
        The active-sector mixing matrix.
    mass_squared : list or np.ndarray
        The active-sector mass-squared values, in eV^2.
    pairs : dict, optional
        The pairing specification; see :func:`pseudo_dirac_mixing_matrix`.  An
        empty mapping reproduces the ordinary Dirac Hamiltonian exactly.
    nubar : bool, optional
        If True, compute the antineutrino Hamiltonian, which conjugates the
        mixing matrix.  Default: False.

    Returns
    -------
    np.ndarray
        The ``n`` x ``n`` complex Hamiltonian, multiplied by the energy.

    Examples
    --------
    With no pairs this is the ordinary three-flavor vacuum Hamiltonian, to
    machine precision.

    .. jupyter-execute::

        import numpy as np
        import magnus.hamiltonians as hamiltonians

        U = hamiltonians.pmns_mixing_matrix(0.5558, 0.6856, 0.1499, 3.7001)
        H = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(
                U, [0.0, 7.5e-5, 2.511e-3], {})
        H3 = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
                0.5558, 0.6856, 0.1499, 3.7001, 7.5e-5, 2.511e-3)

        print(float(np.max(np.abs(H - H3))))
    """
    W = pseudo_dirac_mixing_matrix(mixing_matrix, pairs)
    M2 = pseudo_dirac_mass_squared(mass_squared, pairs)
    if nubar:
        W = np.conj(W)
    return 0.5 * W @ np.diag(M2).astype(complex) @ np.conj(W.T)


def hamiltonian_pseudo_dirac_vacuum(
    energy: float,
    mixing_matrix: Union[Sequence, np.ndarray],
    mass_squared: Union[Sequence[float], np.ndarray],
    pairs: Optional[Mapping[int, float]] = None,
    nubar: Optional[bool] = False
) -> np.ndarray:
    r"""Returns the pseudo-Dirac Hamiltonian for oscillations in vacuum.

    :math:`H_{\rm vac} = W M^2 W^\dagger / (2E)`.

    .. versionadded:: 1.0.5

    Parameters
    ----------
    energy : float
        Neutrino energy, in eV.
    mixing_matrix : list or np.ndarray
        The active-sector mixing matrix.
    mass_squared : list or np.ndarray
        The active-sector mass-squared values, in eV^2.
    pairs : dict, optional
        The pairing specification; see :func:`pseudo_dirac_mixing_matrix`.
    nubar : bool, optional
        If True, compute the antineutrino Hamiltonian.  Default: False.

    Returns
    -------
    np.ndarray
        The ``n`` x ``n`` complex Hamiltonian, in eV.

    Examples
    --------
    .. jupyter-execute::

        import magnus.globaldefs as gd
        import magnus.hamiltonians as hamiltonians

        U = hamiltonians.pmns_mixing_matrix(0.5558, 0.6856, 0.1499, 3.7001)
        H = hamiltonians.hamiltonian_pseudo_dirac_vacuum(
                1.0*gd.UNIT_GEV, U, [0.0, 7.5e-5, 2.511e-3], {1: 1.0e-18})

        print(H.shape)
    """
    return hamiltonian_pseudo_dirac_vacuum_energy_independent(
        mixing_matrix, mass_squared, pairs, nubar=nubar)/energy


def hamiltonian_pseudo_dirac_matter(
    VCC: Union[float, np.ndarray],
    n_active: int,
    pairs: Optional[Mapping[int, float]] = None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]] = 1.0
) -> np.ndarray:
    r"""Returns the matter part of the pseudo-Dirac Hamiltonian.

    :math:`H_{\rm matter} = V_{\rm CC}\,P`, with :math:`P` the flavor structure
    returned by :func:`magnus.matter.matter_potential_projector`.  The sterile
    partners feel neither charged nor neutral current, so once the actives'
    common neutral-current term is removed they are left carrying
    :math:`-V_{\rm NC} = (r/2)\,V_{\rm CC}`, exactly as the sterile states of a
    3+N spectrum do.  The projector is taken from ``matter`` rather than
    written out again here: writing that structure a second time by hand is
    what once gave the sterile states of the NSI route no matter at all.

    Add this to :func:`hamiltonian_pseudo_dirac_vacuum` to obtain the full
    Hamiltonian.  The split is deliberate and preserves separability -- the
    vacuum part carries all of the energy dependence and the matter part all of
    the position dependence -- so the energy-batched engine applies unchanged.

    .. versionadded:: 1.0.5

    Parameters
    ----------
    VCC : float or np.ndarray
        The charged-current matter potential, in eV, already carrying its sign
        for neutrinos or antineutrinos as
        :func:`magnus.matter.vcc_func_from_rho_func` returns it; do not negate
        it again.  An array of potentials, one per position, returns a stack of
        Hamiltonians with the position axis leading.
    n_active : int
        The number of active flavors.  Must be 3: the shared projector places
        the charged-current entry on the first flavor and treats states beyond
        the third as sterile, so any other value would silently mislabel a
        sterile partner as an active flavor.
    pairs : dict, optional
        The pairing specification; only its length is used here, since the
        matter term does not depend on the splittings.
    ratio_number_neutrons_to_protons : int or float, optional
        :math:`r = n_n/n_p` of the medium.  Default: 1.0 (isoscalar).  The same
        limitation the Earth wrappers carry applies here: this is one scalar for
        the whole trajectory, while a layered profile has one per layer.

    Returns
    -------
    np.ndarray
        The ``n`` x ``n`` matter Hamiltonian, or a stack of them.

    Raises
    ------
    ValueError
        If ``n_active`` is not 3.

    Examples
    --------
    Three active flavors and two partners: the charged-current entry on
    :math:`\nu_e`, nothing on the other actives, and :math:`r/2` on each
    sterile.

    .. jupyter-execute::

        import numpy as np
        import magnus.hamiltonians as hamiltonians

        H = hamiltonians.hamiltonian_pseudo_dirac_matter(
                1.0, 3, {0: 1.0e-18, 2: 4.0e-18},
                ratio_number_neutrons_to_protons=1.2)

        print(np.diag(H).real)
    """
    if int(n_active) != 3:
        raise ValueError(
            "Error in magnus: hamiltonians_pseudodirac: the matter term is "
            "available for three active flavors only; got n_active=%d.  "
            "magnus.matter.matter_potential_projector places the "
            "charged-current entry on the first flavor and treats every state "
            "beyond the third as sterile, so for any other number of active "
            "flavors it would label a sterile partner as active and return a "
            "silently wrong matter potential.  Build the projector explicitly "
            "if you need another active-flavor count." % int(n_active))
    _, sterile_of = _pair_layout(int(n_active), pairs)
    n = int(n_active) + len(sterile_of)
    projector = matter.matter_potential_projector(
        n, ratio_number_neutrons_to_protons)
    VCC = np.asarray(VCC, dtype=float)
    return VCC[..., None, None] * projector
