# -*- coding: utf-8 -*-
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


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# from numpy import *
import numpy as np
from typing import Optional, Callable, Union


def mixing_matrix_4x4(s12: float, s23: float, s13:float, d13: float, s14: float, d14: float,
    s24: float, d24: float, s34: float,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the 4x4 (3+1 sterile) mixing matrix.

    Computes and returns the 4x4 complex mixing matrix for a 3+1 sterile-neutrino scenario,
    parametrized by the three standard mixing angles (:math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}`) and CP phase
    (:math:`\delta_{13}`), plus three additional mixing angles (:math:`\theta_{14}`, :math:`\theta_{24}`, :math:`\theta_{34}`) and two
    additional CP phases (:math:`\delta_{14}`, :math:`\delta_{24}`) coupling the sterile state.  Follows the
    parametrization :math:`U = R_{34} \tilde R_{24} \tilde R_{14} R_{23} \tilde R_{13} R_{12}` of
    Kopp, Machado, Maltoni & Schwetz, arXiv:1103.4570 (see also arXiv:1105.3911).

    .. versionadded:: 0.10.0

    Parameters
    ----------
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    d13 : float
        :math:`\delta_{13}` [radian].
    s14 : float
        Sine of the mixing angle :math:`\theta_{14}`.
    d14 : float
        :math:`\delta_{14}` [radian].
    s24 : float
        Sine of the mixing angle :math:`\theta_{24}`.
    d24 : float
        :math:`\delta_{24}` [radian].
    s34 : float
        Sine of the mixing angle :math:`\theta_{34}`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed closed-form expressions for each entry;
        otherwise, build the matrix by multiplying the five rotation matrices live. Both paths
        must (and do, see ``tests/test_hamiltonians.py``) agree to machine precision.

    Returns
    -------
    np.ndarray
        4x4 mixing matrix.
    """
    # arXiv:1105.3911

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
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino (3+1) Hamiltonian for vacuum oscillations.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations in vacuum,
    parametrized by the six 3+1 mixing angles and two CP phases of :func:`mixing_matrix_4x4`, and
    three mass-squared differences (:math:`\Delta m_{21}^2`, :math:`\Delta m_{31}^2`, :math:`\Delta m_{41}^2`).  The Hamiltonian is
    H = (1/2)*R.M2.R^dagger, with R the 4x4 mixing matrix and M2 the mass matrix.  The
    multiplicative factor 1/E is not applied.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    d13 : float
        :math:`\delta_{13}` [radian].
    s14 : float
        Sine of the mixing angle :math:`\theta_{14}`.
    d14 : float
        :math:`\delta_{14}` [radian].
    s24 : float
        Sine of the mixing angle :math:`\theta_{24}`.
    d24 : float
        :math:`\delta_{24}` [radian].
    s34 : float
        Sine of the mixing angle :math:`\theta_{34}`.
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

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """
    # 4x4 mixing matrix
    R = mixing_matrix_4x4(s12, s23, s13, d13, s14, d14, s24, d24, s34,
        compute_matrix_multiplication=compute_matrix_multiplication) if nubar == False else \
            np.conj(mixing_matrix_4x4(s12, s23, s13, d13, s14, d14, s24, d24, s34,
                compute_matrix_multiplication=compute_matrix_multiplication))
    # Mass matrix
    M2 = np.diag([0.0, D21, D31, D41])

    return 0.5 * np.linalg.multi_dot([R, M2, np.conj(R.T)])
    # return 0.5 * R @ M2 @ np.conj(R.T)


def hamiltonian_4nu_vacuum_energy_independent_td(l: float, s12: float, s23: float, s13:float,
    d13: float, s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float,
    D41: float, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Same as :func:`hamiltonian_4nu_vacuum_energy_independent`, included for interface parity with
    the other, genuinely position-dependent Hamiltonians.

    .. versionadded:: 0.10.0

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

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, d13, s14, d14, s24, d24, s34,
        D21, D31, D41, nubar=nubar, compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_4nu_vacuum(energy: float, s12: float, s23: float, s13:float, d13: float,
    s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float, D41: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for vacuum oscillations.

    Same as :func:`hamiltonian_4nu_vacuum_energy_independent`, but with the 1/E factor applied.

    .. versionadded:: 0.10.0

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

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """
    return (1/energy)*hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, d13, s14, d14, s24,
        d24, s34, D21, D31, D41, nubar=nubar,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_4nu_vacuum_td(l: float, energy: float, s12: float, s23: float, s13:float, d13: float,
    s14: float, d14: float, s24: float, d24: float, s34: float, D21: float, D31: float, D41: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Same as :func:`hamiltonian_4nu_vacuum`, included for interface parity with the other,
    genuinely position-dependent Hamiltonians.

    .. versionadded:: 0.10.0

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

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_vacuum(energy, s12, s23, s13, d13, s14, d14, s24, d24, s34, D21, D31,
        D41, nubar=nubar, compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_4nu_matter(VCC: float) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 4x4 real four-neutrino Hamiltonian for
    oscillations in matter with constant density.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """
    return np.diag([VCC, 0.0, 0.0, 0.0])


def hamiltonian_4nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for matter oscillations, as a function of distance.

    Computes and returns the 4x4 real four-neutrino Hamiltonian for oscillations in matter with a
    given density as a function of position.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    VCC_func : Callable
        Potential due to charged-current interactions of nu_e with electrons, as a function of
        position, l.

    Returns
    -------
    list
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

    .. versionadded:: 0.10.0

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
    list
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

    .. versionadded:: 0.10.0

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
    list
        Hamiltonian 4x4 matrix.
    """
    return hamiltonian_4nu_nsi(VCC_func(l), eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms,
        eps_tt, eps_ts, eps_ss)


def hamiltonian_4nu_liv(energy: float, sxi12: float, sxi23: float, sxi13: float, dxi13: float,
    sxi14: float, dxi14: float, sxi24: float, dxi24: float, sxi34: float, b1: float, b2: float,
    b3: float, b4: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations in a CPT-odd
    Lorentz invariance-violating background.  Same as
    :func:`hamiltonian_4nu_liv_energy_independent`, but with the
    :math:`E^{n_{\rm liv}}` energy dependence of the LIV operator applied.

    .. versionadded:: 0.10.0

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

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """

    return pow(energy, n_liv) * hamiltonian_4nu_liv_energy_independent(sxi12, sxi23, sxi13, dxi13,
        sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, n_liv, nubar=nubar,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_4nu_liv_energy_independent(sxi12: float, sxi23: float, sxi13: float, dxi13: float,
    sxi14: float, dxi14: float, sxi24: float, dxi24: float, sxi34: float, b1: float, b2: float,
    b3: float, b4: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 4x4 complex four-neutrino Hamiltonian for oscillations in a CPT-odd
    Lorentz invariance-violating background, without the energy-dependent prefactor.

    .. versionadded:: 0.10.0

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

    Returns
    -------
    list
        Hamiltonian 4x4 matrix.
    """
    # 4x4 mixing matrix
    R = mixing_matrix_4x4(sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34,
        compute_matrix_multiplication=compute_matrix_multiplication) if nubar == False else \
            np.conj(mixing_matrix_4x4(sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34,
                compute_matrix_multiplication=compute_matrix_multiplication))

    return pow(1.0/Lambda, n_liv) * R @ np.diag([b1, b2, b3, b4]) @ np.conj(R.T)
