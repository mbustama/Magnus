# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""hamiltonians4nu.py

Compute four-neutrino (3+1 sterile) Hamiltonians for selected scenarios.

This module contains the routines to compute the four-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * mixing_matrix_4x4 - Returns 4x4 PMNS-like mixing matrix (3+1)
    * hamiltonian_4nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_4nu_vacuum_energy_independent_td - Returns H_vac (no
           1/E), as a function of position
    * hamiltonian_4nu_vacuum - Returns H_vac
    * hamiltonian_4nu_vacuum_td - Returns H_vac, as a function of position
    * hamiltonian_4nu_matter - Returns H_matter
    * hamiltonian_4nu_matter_td - Returns H_matter, as a function of position
    * hamiltonian_4nu_nsi - Returns H_NSI
    * hamiltonian_4nu_nsi_td - Returns H_NSI, as a function of position
    * hamiltonian_4nu_liv - Returns H_LIV
    * hamiltonian_4nu_liv_energy_independent - Returns H_LIV (no energy
           dependence)
"""


__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# from numpy import *
import numpy as np

from magnus.hamiltonians import _angles

import magnus.matter as matter
from typing import Optional, Callable, Union


def mixing_matrix_4x4(s12: float, s23: float, s13:float, d13: float, s14: float, d14: float,
    s24: float, d24: float, s34: float,
    compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the 4x4 (3+1 sterile) mixing matrix.

    Computes and returns the 4x4 complex mixing matrix for a 3+1 sterile-neutrino scenario,
    parametrized by the three standard mixing angles (:math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}`) and CP phase
    (:math:`\delta_{13}`), plus three additional mixing angles (:math:`\theta_{14}`, :math:`\theta_{24}`, :math:`\theta_{34}`) and two
    additional CP phases (:math:`\delta_{14}`, :math:`\delta_{24}`) coupling the sterile state.  Follows the
    parametrization :math:`U = R_{34} \tilde R_{24} \tilde R_{14} R_{23} \tilde R_{13} R_{12}` of
    Kopp, Machado, Maltoni & Schwetz, arXiv:1103.4570 (see also arXiv:1105.3911).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    s12 : float
        Mixing angle :math:`\theta_{12}`, in the convention set by ``angles`` (default: its sine).
    s23 : float
        Mixing angle :math:`\theta_{23}`, in the convention set by ``angles`` (default: its sine).
    s13 : float
        Mixing angle :math:`\theta_{13}`, in the convention set by ``angles`` (default: its sine).
    d13 : float
        :math:`\delta_{13}` [radian, or degree if ``angles='deg'``].
    s14 : float
        Mixing angle :math:`\theta_{14}`, in the convention set by ``angles`` (default: its sine).
    d14 : float
        :math:`\delta_{14}` [radian, or degree if ``angles='deg'``].
    s24 : float
        Mixing angle :math:`\theta_{24}`, in the convention set by ``angles`` (default: its sine).
    d24 : float
        :math:`\delta_{24}` [radian, or degree if ``angles='deg'``].
    s34 : float
        Mixing angle :math:`\theta_{34}`, in the convention set by ``angles`` (default: its sine).
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed closed-form expressions for each entry;
        otherwise, build the matrix by multiplying the five rotation matrices live. Both paths
        must (and do, see ``tests/test_hamiltonians.py``) agree to machine precision.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        4x4 mixing matrix.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus.hamiltonians import hamiltonians4nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        U = np.asarray(hamiltonians4nu.mixing_matrix_4x4(
            p['s12'], p['s23'], p['s13'], p['dCP'], 0.15, 0.0, 0.10, 0.0, 0.05))

        print('shape', U.shape)
        print('unitary to %.1e' % np.max(np.abs(U.conj().T @ U - np.eye(4))))
"""
    # arXiv:1105.3911

    _r, _p = _angles.resolve(
        'hamiltonians.mixing_matrix_4x4', angles,
        {'s12': s12, 's23': s23, 's13': s13, 's14': s14, 's24': s24, 's34': s34},
        {'d13': d13, 'd14': d14, 'd24': d24})
    s12, s23, s13 = _r['s12'], _r['s23'], _r['s13']
    s14, s24, s34 = _r['s14'], _r['s24'], _r['s34']
    d13, d14, d24 = _p['d13'], _p['d14'], _p['d24']

    c12 = np.sqrt(1.0-s12*s12)
    c23 = np.sqrt(1.0-s23*s23)
    c13 = np.sqrt(1.0-s13*s13)
    c14 = np.sqrt(1.0-s14*s14)
    c24 = np.sqrt(1.0-s24*s24)
    c34 = np.sqrt(1.0-s34*s34)
    cd13 = np.cos(d13)
    sd13 = np.sin(d13)
    exp_d13_p = complex(cd13, sd13)
    exp_d13_m = np.conj(exp_d13_p)
    cd14 = np.cos(d14)
    sd14 = np.sin(d14)
    exp_d14_p = complex(cd14, sd14)
    exp_d14_m = np.conj(exp_d14_p)
    cd24 = np.cos(d24)
    sd24 = np.sin(d24)
    exp_d24_p = complex(cd24, sd24)
    exp_d24_m = np.conj(exp_d24_p)

    if not compute_matrix_multiplication:

        U00 = c12*c13*c14
        U01 = c13*c14*s12
        U02 = c14*s13*exp_d13_m
        U03 = s14*exp_d14_m

        f1 = -c24*s13*s23*exp_d13_p-c13*s14*s24*exp_d14_p*exp_d24_m
        U10 = -c23*c24*s12 + c12*f1
        U11 = c12*c23*c24 + s12*f1
        U12 = c13*c24*s23 - s13*s14*s24*exp_d13_m*exp_d14_p*exp_d24_m
        U13 = c14*s24*exp_d24_m

        f2 = -c34*s23 - c23*s24*s34*exp_d24_p
        f3 = -c13*c24*s14*s34*exp_d14_p - s13*exp_d13_p*(c23*c34-s23*s24*s34*exp_d24_p)
        U20 = -s12*f2 + c12*f3
        U21 = c12*f2 + s12*f3
        U22 = -c24*s13*s14*s34*exp_d13_m*exp_d14_p + c13*(c23*c34-s23*s24*s34*exp_d24_p)
        U23 = c14*c24*s34

        f4 = -c23*c34*s24*exp_d24_p + s23*s34
        f5 = -c13*c24*c34*s14*exp_d14_p - s13*exp_d13_p*(-c34*s23*s24*exp_d24_p - c23*s34)
        U30 = -s12*f4 + c12*f5
        U31 = c12*f4 + s12*f5
        U32 = -c24*c34*s13*s14*exp_d13_m*exp_d14_p + c13*(-c34*s23*s24*exp_d24_p - c23*s34)
        U33 = c14*c24*c34

        return np.array([[U00,U01,U02,U03],[U10,U11,U12,U13],[U20,U21,U22,U23],[U30,U31,U32,U33]])

    else:

        # U = R34.~R24.~R14.R23.~R13.R12
        R12 = np.array([[c12, s12, 0, 0], [-s12, c12, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        R13 = np.array([[c13, 0, s13*exp_d13_m, 0], [0, 1, 0, 0], [-s13*exp_d13_p, 0, c13, 0],
            [0, 0, 0, 1]])
        R23 = np.array([[1, 0, 0, 0], [0, c23, s23, 0], [0, -s23, c23, 0], [0, 0, 0, 1]])
        R14 = np.array([[c14, 0, 0, s14*exp_d14_m], [0, 1, 0, 0], [0, 0, 1, 0],
            [-s14*exp_d14_p, 0, 0, c14]])
        R24 = np.array([[1, 0, 0, 0], [0, c24, 0, s24*exp_d24_m], [0, 0, 1, 0],
            [0, -s24*exp_d24_p, 0, c24]])
        R34 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, c34, s34], [0, 0, -s34, c34]])

        return np.linalg.multi_dot([R34, R24, R14, R23, R13, R12])


def hamiltonian_4nu_vacuum_energy_independent(s12: float, s23: float, s13:float, d13: float,
    s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float, D41: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the four-neutrino (3+1) Hamiltonian for vacuum oscillations.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations in vacuum,
    parametrized by the six 3+1 mixing angles and two CP phases of :func:`mixing_matrix_4x4`, and
    three mass-squared differences (:math:`\Delta m_{21}^2`, :math:`\Delta m_{31}^2`, :math:`\Delta m_{41}^2`).  The Hamiltonian is
    H = (1/2)*R.M2.R^dagger, with R the 4x4 mixing matrix and M2 the mass matrix.  The
    multiplicative factor 1/E is not applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    s12 : float
        Mixing angle :math:`\theta_{12}`, in the convention set by ``angles`` (default: its sine).
    s23 : float
        Mixing angle :math:`\theta_{23}`, in the convention set by ``angles`` (default: its sine).
    s13 : float
        Mixing angle :math:`\theta_{13}`, in the convention set by ``angles`` (default: its sine).
    d13 : float
        :math:`\delta_{13}` [radian, or degree if ``angles='deg'``].
    s14 : float
        Mixing angle :math:`\theta_{14}`, in the convention set by ``angles`` (default: its sine).
    d14 : float
        :math:`\delta_{14}` [radian, or degree if ``angles='deg'``].
    s24 : float
        Mixing angle :math:`\theta_{24}`, in the convention set by ``angles`` (default: its sine).
    d24 : float
        :math:`\delta_{24}` [radian, or degree if ``angles='deg'``].
    s34 : float
        Mixing angle :math:`\theta_{34}`, in the convention set by ``angles`` (default: its sine).
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    D41 : float
        Mass-squared difference :math:`\Delta m_{41}^2`.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos (conjugates the mixing matrix,
        equivalent to negating every CP phase). Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`mixing_matrix_4x4`. If False (default), use the pre-computed
        expressions; otherwise, multiply R.M2.R^dagger live.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    # Converted here rather than forwarded to mixing_matrix_4x4, so the message names THIS
    # function: the mistake this guard exists for is a phase landing in a sine slot when
    # thirteen arguments are passed positionally, and naming the mixing matrix instead would
    # point the reader at the wrong signature.
    _r, _p = _angles.resolve(
        'hamiltonians.hamiltonian_4nu_vacuum_energy_independent', angles,
        {'s12': s12, 's23': s23, 's13': s13, 's14': s14, 's24': s24, 's34': s34},
        {'d13': d13, 'd14': d14, 'd24': d24})
    s12, s23, s13 = _r['s12'], _r['s23'], _r['s13']
    s14, s24, s34 = _r['s14'], _r['s24'], _r['s34']
    d13, d14, d24 = _p['d13'], _p['d14'], _p['d24']
    # 4x4 mixing matrix
    R = mixing_matrix_4x4(s12, s23, s13, d13, s14, d14, s24, d24, s34,
        compute_matrix_multiplication=compute_matrix_multiplication) if not nubar else \
            np.conj(mixing_matrix_4x4(s12, s23, s13, d13, s14, d14, s24, d24, s34,
                compute_matrix_multiplication=compute_matrix_multiplication))
    # Mass matrix
    M2 = np.diag([0.0, D21, D31, D41])

    return 0.5 * np.linalg.multi_dot([R, M2, np.conj(R.T)])
    # return 0.5 * R @ M2 @ np.conj(R.T)


def hamiltonian_4nu_vacuum_energy_independent_td(l: float, s12: float, s23: float, s13:float,
    d13: float, s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float,
    D41: float, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Same as :func:`hamiltonian_4nu_vacuum_energy_independent`, included for interface parity with
    the other, genuinely position-dependent Hamiltonians.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    s12, s23, s13, d13, s14, d14, s24, d24, s34 : float
        3+1 mixing angles (sines) and CP phases; see :func:`mixing_matrix_4x4`.
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    D41 : float
        Mass-squared difference :math:`\Delta m_{41}^2`.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos. Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`mixing_matrix_4x4`.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, d13, s14, d14, s24, d24, s34,
        D21, D31, D41, nubar=nubar, compute_matrix_multiplication=compute_matrix_multiplication,
        angles=angles)


def hamiltonian_4nu_vacuum(energy: float, s12: float, s23: float, s13:float, d13: float,
    s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float, D41: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for vacuum oscillations.

    Same as :func:`hamiltonian_4nu_vacuum_energy_independent`, but with the 1/E factor applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float
        Neutrino energy.
    s12, s23, s13, d13, s14, d14, s24, d24, s34 : float
        3+1 mixing angles (sines) and CP phases; see :func:`mixing_matrix_4x4`.
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    D41 : float
        Mass-squared difference :math:`\Delta m_{41}^2`.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos. Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`mixing_matrix_4x4`.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    return (1/energy)*hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, d13, s14, d14, s24,
        d24, s34, D21, D31, D41, nubar=nubar,
        compute_matrix_multiplication=compute_matrix_multiplication, angles=angles)


def hamiltonian_4nu_vacuum_td(l: float, energy: float, s12: float, s23: float, s13:float, d13: float,
    s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float, D41: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Same as :func:`hamiltonian_4nu_vacuum`, included for interface parity with the other,
    genuinely position-dependent Hamiltonians.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    energy : float
        Neutrino energy.
    s12, s23, s13, d13, s14, d14, s24, d24, s34 : float
        3+1 mixing angles (sines) and CP phases; see :func:`mixing_matrix_4x4`.
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    D41 : float
        Mass-squared difference :math:`\Delta m_{41}^2`.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos. Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`mixing_matrix_4x4`.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_vacuum(energy, s12, s23, s13, d13, s14, d14, s24, d24, s34, D21, D31,
        D41, nubar=nubar, compute_matrix_multiplication=compute_matrix_multiplication,
        angles=angles)


def hamiltonian_4nu_matter(VCC: float,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]] = 1.0
) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 4x4 real four-neutrino Hamiltonian for
    oscillations in matter with constant density.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.
    ratio_number_neutrons_to_protons : int or float, optional
        :math:`r = n_n/n_p` of the medium, which sets the sterile states' entry
        via :math:`-V_{\rm NC} = (r/2) V_{\rm CC}`.  Must match the value given
        to :func:`magnus.matter.vcc_func_from_rho_func`.  Default: 1.0
        (isoscalar matter).

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus.hamiltonians import hamiltonians4nu

        print(np.asarray(hamiltonians4nu.hamiltonian_4nu_matter(1.0e-13)))

    The sterile state feels neither the charged- nor the neutral-current
    potential -- it is a Standard Model gauge singlet, so there is no W and no
    Z for it to scatter off.  That is *why* its entry here is not zero.
    Oscillations depend only on differences between the diagonal entries, so a
    multiple of the identity is unobservable and the actives' common
    :math:`V_{\rm NC}` is subtracted.  What that leaves on the sterile state is
    :math:`-V_{\rm NC} = (r/2) V_{\rm CC}`, with :math:`r = n_n/n_p`, and it is
    that residue -- not any interaction of the sterile state -- that makes a 3+1
    scenario more than a relabeling.  Setting it to zero instead is worth 0.29
    in probability on a PREM chord.  The derivation is in
    :func:`magnus.matter.matter_potential_projector`.
"""
    # Built by broadcasting rather than np.diag so that VCC may be an array of
    # positions: VCC[..., None, None] turns one potential per position into a
    # stack of matrices, which is what lets a caller's H_func take the engine's
    # vectorized path (see magnus.magnus.ScalarHamiltonianWarning). A scalar VCC
    # still returns a plain (4, 4) matrix.
    VCC = np.asarray(VCC, dtype=float)
    # NOT e_ee: a sterile state feels neither current, so once the actives' common V_NC
    # is removed it carries -V_NC = (r/2) V_CC.  The docstring above has always said the
    # sterile state feels neither potential; the code used to implement only half of that,
    # which costs 0.29 in probability on a 3+1 PREM chord.  See
    # matter.matter_potential_projector, which is the one definition of this structure.
    proj = matter.matter_potential_projector(4, ratio_number_neutrons_to_protons)
    return VCC[..., None, None] * proj


def hamiltonian_4nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for matter oscillations, as a function of distance.

    Computes and returns the 4x4 real four-neutrino Hamiltonian for oscillations in matter with a
    given density as a function of position.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    VCC_func : Callable
        Potential due to charged-current interactions of nu_e with electrons, as a function of
        position, l.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_matter(VCC_func(l))


def hamiltonian_4nu_nsi(
    VCC: float,
    eps_ee: float,
    eps_em: complex,
    eps_et: complex,
    eps_es: complex,
    eps_mm: float,
    eps_mt: complex,
    eps_ms: complex,
    eps_tt: float,
    eps_ts: complex,
    eps_ss: float
) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for oscillations w/ NSI.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations with
    non-standard interactions (NSI) in matter with constant density.  The additional 's' subscript
    denotes the sterile flavor.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    VCC : float
        Potential due to charged-current interactions of nu_e with electrons.
    eps_ee : float
        Diagonal NSI coupling of nu_e.
    eps_em : complex
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling.
    eps_et : complex
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling.
    eps_es : complex
        Flavor-off-diagonal (nu_e-nu_s) NSI coupling.
    eps_mm : float
        Diagonal NSI coupling of nu_mu.
    eps_mt : complex
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling.
    eps_ms : complex
        Flavor-off-diagonal (nu_mu-nu_s) NSI coupling.
    eps_tt : float
        Diagonal NSI coupling of nu_tau.
    eps_ts : complex
        Flavor-off-diagonal (nu_tau-nu_s) NSI coupling.
    eps_ss : float
        Diagonal NSI coupling of nu_s.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    return VCC * np.array([
        [eps_ee, eps_em, eps_et, eps_es],
        [np.conj(eps_em), eps_mm, eps_mt, eps_ms],
        [np.conj(eps_et), np.conj(eps_mt), eps_tt, eps_ts],
        [np.conj(eps_es), np.conj(eps_ms), np.conj(eps_ts), eps_ss]
        ], dtype=np.complex128)


def hamiltonian_4nu_nsi_td(l: float, VCC_func: Callable, eps_ee: float, eps_em: complex,
    eps_et: complex, eps_es: complex, eps_mm: float, eps_mt: complex, eps_ms: complex,
    eps_tt: float, eps_ts: complex, eps_ss: float) -> np.ndarray:
    r"""Returns the four-neutrino NSI Hamiltonian as a function of position.

    Same as :func:`hamiltonian_4nu_nsi`, but evaluates the position-dependent matter potential
    ``VCC_func(l)`` first.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    VCC_func : Callable
        Potential due to charged-current interactions of nu_e with electrons, as a function of
        position, l.
    eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss :
        NSI coupling parameters; see :func:`hamiltonian_4nu_nsi`.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_nsi(VCC_func(l), eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms,
        eps_tt, eps_ts, eps_ss)


def hamiltonian_4nu_liv(energy: float, sxi12: float, sxi23: float, sxi13: float, dxi13: float,
    sxi14: float, dxi14: float, sxi24: float, dxi24: float, sxi34: float, b1: float, b2: float,
    b3: float, b4: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations in a CPT-odd
    Lorentz invariance-violating background.  Same as
    :func:`hamiltonian_4nu_liv_energy_independent`, but with the
    :math:`E^{n_{\rm liv}}` energy dependence of the LIV operator applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float
        Neutrino energy.
    sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 : float
        Sines of the mixing angles between the space of the eigenvectors of the LIV operator B4
        and the flavor states, parametrized as in :func:`mixing_matrix_4x4`.
    dxi13, dxi14, dxi24 : float
        CP-violation phases of the LIV operator B4 [radian].
    b1 : float
        Eigenvalue b1 of the LIV operator B4.
    b2 : float
        Eigenvalue b2 of the LIV operator B4.
    b3 : float
        Eigenvalue b3 of the LIV operator B4.
    b4 : float
        Eigenvalue b4 of the LIV operator B4.
    Lambda : float
        Energy scale of the LIV operator B4.
    n_liv : int
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3).
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos (conjugates the LIV mixing matrix).
        Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`mixing_matrix_4x4`.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """

    return pow(energy, n_liv) * hamiltonian_4nu_liv_energy_independent(sxi12, sxi23, sxi13, dxi13,
        sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, n_liv, nubar=nubar,
        compute_matrix_multiplication=compute_matrix_multiplication, angles=angles)


def hamiltonian_4nu_liv_energy_independent(sxi12: float, sxi23: float, sxi13: float, dxi13: float,
    sxi14: float, dxi14: float, sxi24: float, dxi24: float, sxi34: float, b1: float, b2: float,
    b3: float, b4: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False,
    angles: Optional[str]='sin') -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations in a CPT-odd
    Lorentz invariance-violating background, without the energy-dependent prefactor.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 : float
        Sines of the mixing angles between the space of the eigenvectors of the LIV operator B4
        and the flavor states, parametrized as in :func:`mixing_matrix_4x4`.
    dxi13, dxi14, dxi24 : float
        CP-violation phases of the LIV operator B4 [radian].
    b1 : float
        Eigenvalue b1 of the LIV operator B4.
    b2 : float
        Eigenvalue b2 of the LIV operator B4.
    b3 : float
        Eigenvalue b3 of the LIV operator B4.
    b4 : float
        Eigenvalue b4 of the LIV operator B4.
    Lambda : float
        Energy scale of the LIV operator B4.
    n_liv : int
        Power of the energy dependence of the LIV operator; enters through the
        :math:`\Lambda^{-n_{\rm liv}}` normalization of the eigenvalues.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos (conjugates the LIV mixing matrix).
        Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`mixing_matrix_4x4`.
    angles : str, optional
        How the mixing angles are stated: ``'sin'`` (default) their sines, ``'sin2'``
        their sines *squared* -- which is what global fits report -- ``'rad'`` the angles
        themselves in radians, or ``'deg'`` in degrees.  Any other value raises.  Under
        ``'deg'`` the CP phases are read as degrees too; under the other three
        they stay in radians, a sine being no way to state a phase.

    Returns
    -------
    np.ndarray
        Hamiltonian 4x4 matrix.
    """
    # The LIV angles went through no guard at all before this.
    _r, _p = _angles.resolve(
        'hamiltonians.hamiltonian_4nu_liv_energy_independent', angles,
        {'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'sxi14': sxi14,
         'sxi24': sxi24, 'sxi34': sxi34},
        {'dxi13': dxi13, 'dxi14': dxi14, 'dxi24': dxi24})
    sxi12, sxi23, sxi13 = _r['sxi12'], _r['sxi23'], _r['sxi13']
    sxi14, sxi24, sxi34 = _r['sxi14'], _r['sxi24'], _r['sxi34']
    dxi13, dxi14, dxi24 = _p['dxi13'], _p['dxi14'], _p['dxi24']

    # 4x4 mixing matrix
    R = mixing_matrix_4x4(sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34,
        compute_matrix_multiplication=compute_matrix_multiplication) if not nubar else \
            np.conj(mixing_matrix_4x4(sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34,
                compute_matrix_multiplication=compute_matrix_multiplication))

    return pow(1.0/Lambda, n_liv) * R @ np.diag([b1, b2, b3, b4]) @ np.conj(R.T)
