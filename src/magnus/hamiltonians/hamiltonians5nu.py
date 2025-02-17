# -*- coding: utf-8 -*-
r"""Compute three-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the three-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * mixing_matrix_2nu - Returns 2x2 rotation matrix
    * hamiltonian_2nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * delta - Kronecker delta
    * J - Product of four elements of PMNS matrix
    * probabilities_3nu_vacuum_std - Vacuum probability, std. formula
    * hamiltonian_2nu_matter - Returns H_matter
    * hamiltonian_2nu_nsi - Returns H_NSI
    * hamiltonian_2nu_liv - Returns H_LIV

Created: 2019/04/17 17:14
Last modified: 2019/04/30 01:03
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# from numpy import *
import numpy as np
from typing import Optional, Callable, Union


def mixing_matrix_5x5(s12: float, s23: float, s13:float, d13: float, s14: float, d14: float,
    s15: float, d15: float, s24: float, d24: float, s25: float, s34: float, s35: float, d35: float,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    # arXiv:1105.3911

    c12 = np.sqrt(1.0-s12*s12)
    c23 = np.sqrt(1.0-s23*s23)
    c13 = np.sqrt(1.0-s13*s13)
    c14 = np.sqrt(1.0-s14*s14)
    c15 = np.sqrt(1.0-s15*s15)
    c24 = np.sqrt(1.0-s24*s24)
    c25 = np.sqrt(1.0-s25*s25)
    c34 = np.sqrt(1.0-s34*s34)
    c35 = np.sqrt(1.0-s35*s35)
    cd13 = np.cos(d13)
    sd13 = np.sin(d13)
    exp_d13_p = complex(cd13, sd13)
    exp_d13_m = np.conj(exp_d13_p)
    cd14 = np.cos(d14)
    sd14 = np.sin(d14)
    exp_d14_p = complex(cd14, sd14)
    exp_d14_m = np.conj(exp_d14_p)
    cd15 = np.cos(d15)
    sd15 = np.sin(d15)
    exp_d15_p = complex(cd15, sd15)
    exp_d15_m = np.conj(exp_d15_p)
    cd24 = np.cos(d24)
    sd24 = np.sin(d24)
    exp_d24_p = complex(cd24, sd24)
    exp_d24_m = np.conj(exp_d24_p)
    cd35 = np.cos(d35)
    sd35 = np.sin(d35)
    exp_d35_p = complex(cd35, sd35)
    exp_d35_m = np.conj(exp_d35_p)

    if not compute_matrix_multiplication:

        U00 = c12*c13*c14*c15
        U01 = c13*c14*c15*s12
        U02 = c14*c15*s13*exp_d13_m
        U03 = s14*c15*exp_d14_m
        U04 = s15*exp_d15_m

        f1 = -c25*s14*s24*exp_d14_p*exp_d24_m-c14*s15*s25*exp_d15_p
        f2 = -c24*c25*s13*s23*exp_d13_p + c13*f1
        U10 = -c23*c24*c25*s12 + c12*f2
        U11 = c12*c23*c24*c25 + s12*f2
        U12 = c13*c24*c25*s23 + s13*exp_d13_m*f1
        U13 = c14*c25*s24*exp_d24_m - s14*s15*s25*exp_d14_m*exp_d15_p
        U14 = c15*s25
       
        f3 = -c34*c35*s23 + c23*(-c35*s24*s34*exp_d24_p-c24*s25*s35*exp_d35_m)
        f4 = -c14*c25*s15*s35*exp_d15_p*exp_d35_m \
                - s14*exp_d14_p*(c24*c35*s34-s24*s25*s35*exp_d24_m*exp_d35_m)
        f5 = c23*c34*c35 + s23*(-c35*s24*s34*exp_d24_p-c24*s25*s35*exp_d35_m)
        f6 = -s13*exp_d13_p*f5
        U20 = -s12*f3 + c12*(f6 + c13*f4)
        U21 = c12*f3 + s12*(f6 + c13*f4)
        U22 = c13*(f5 + s13*exp_d13_m*f4)
        U23 = -c25*s14*s15*s35*exp_d14_m*exp_d15_p*exp_d35_m \
                + c14*(c24*c35*s34-s24*s25*s35*exp_d24_m*exp_d35_m)
        U24 = c15*c25*s35*exp_d35_m

        f7 = -c23*c34*s24*exp_d24_p + s23*s34
        f8 = -c34*s23*s24*exp_d24_p - c23*s34
        f9 = -c13*c24*c34*s14*exp_d14_p - s13*exp_d13_p*f8
        U30 = -s12*f7 + c12*f9
        U31 = c12*f7 + s12*f9
        U32 = -c24*c34*s13*s14*exp_d13_m*exp_d14_p + c13*f8
        U33 = c14*c24*c34
        U34 = 0

        f10 = c34*s23*s35*exp_d35_p + c23*(-c24*c35*s25+s24*s34*s35*exp_d24_p*exp_d35_p)
        f11 = -c35*s24*s25*exp_d24_m
        f12 = -c14*c25*c35*s15*exp_d15_p - s14*exp_d14_p*(f11-c24*s34*s35*exp_d35_p)
        f13 = -c23*c34*s35*exp_d35_p + s23*(-c24*c35*s25+s24*s34*s35*exp_d24_p*exp_d35_p)
        f14 = c13*f12 - s13*exp_d13_p*f13
        U40 = -s12*f10 + c12*f14
        U41 = c12*f10 + s12*f14
        U42 = s13*exp_d13_m*f12 + c13*f13
        U43 = -c25*c35*s14*s15*exp_d14_m*exp_d15_p + c14*(f11-c24*s34*s35*exp_d35_p)
        U44 = c15*c25*c35

        return np.array([
            [U00,U01,U02,U03,U04],
            [U10,U11,U12,U13,U14],
            [U20,U21,U22,U23,U24],
            [U30,U31,U32,U33,U34],
            [U40,U41,U42,U43,U44]])

    else:

        # U = ~R35.R34.R25.~R24.R23.~R15.~R14.~R13.R12
        R12 = np.array([
            [c12, s12, 0, 0, 0], 
            [-s12, c12, 0, 0, 0], 
            [0, 0, 1, 0, 0], 
            [0, 0, 0, 1, 0], 
            [0, 0, 0, 0, 1]])
        R13 = np.array([
            [c13, 0, s13*exp_d13_m, 0, 0], 
            [0, 1, 0, 0, 0], 
            [-s13*exp_d13_p, 0, c13, 0, 0], 
            [0, 0, 0, 1, 0], 
            [0, 0, 0, 0, 1]])
        R14 = np.array([
            [c14, 0, 0, s14*exp_d14_m, 0], 
            [0, 1, 0, 0, 0], 
            [0, 0, 1, 0, 0], 
            [-s14*exp_d14_p, 0, 0, c14, 0],
            [0, 0, 0, 0, 1]])
        R15 = np.array([
            [c15, 0, 0, 0, s15*exp_d15_m], 
            [0, 1, 0, 0, 0], 
            [0, 0, 1, 0, 0], 
            [0, 0, 0, 1, 0],
            [-s15*exp_d15_p, 0, 0, 0, c15]])
        R23 = np.array([
            [1, 0, 0, 0, 0], 
            [0, c23, s23, 0, 0], 
            [0, -s23, c23, 0, 0], 
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]])
        R24 = np.array([
            [1, 0, 0, 0, 0], 
            [0, c24, 0, s24*exp_d24_m, 0], 
            [0, 0, 1, 0, 0], 
            [0, -s24*exp_d24_p, 0, c24, 0],
            [0, 0, 0, 0, 1]])
        R25 = np.array([
            [1, 0, 0, 0, 0], 
            [0, c25, 0, 0, s25], 
            [0, 0, 1, 0, 0], 
            [0, 0, 0, 1, 0],
            [0, -s25, 0, 0, c25]])
        R34 = np.array([
            [1, 0, 0, 0, 0], 
            [0, 1, 0, 0, 0], 
            [0, 0, c34, s34, 0], 
            [0, 0, -s34, c34, 0],
            [0, 0, 0, 0, 1]])
        R35 = np.array([
            [1, 0, 0, 0, 0], 
            [0, 1, 0, 0, 0], 
            [0, 0, c35, 0, s35*exp_d35_p], 
            [0, 0, 0, 1, 0],
            [0, 0, -s35*exp_d35_p, 0, c35]])

        return np.linalg.multi_dot([R35, R34, R25, R24, R23, R15, R14, R13, R12])


def hamiltonian_5nu_vacuum_energy_independent(s12: float, s23: float, s13:float, d13: float, 
    s14: float, d14: float, s15: float, d15: float, s24: float, d24: float, s25: float, s34: float, 
    s35: float, d35: float, D21: float, D31: float, D41: float, D51: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the four-neutrino (3+1) Hamiltonian for vacuum oscillations.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations in vacuum, 
    parametrized by three mixing angles (theta_12, theta_23, theta_13), one CP-violation phase 
    (delta_CP), and two mass-squared difference (Delta m^2_21, Delta m^2_31).  The Hamiltonian is
    H = (1/2)*R.M2.R^dagger, with R the 3x3 PMNS matrix and M2 the mass matrix.  The multiplicative
    factor 1/E is not applied.

    Parameters
    ----------
    s12 : float
        Sin(theta_12).
    s23 : float
        Sin(theta_23).
    s13 : float
        Sin(theta_13).
    D21 : float
        Mass-squared difference Delta m^2_21.
    D31 : float
        Mass-squared difference Delta m^2_31.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger 
        live.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """    
    # 5x5 mixing matrix
    R = mixing_matrix_5x5(s12, s23, s13, d13, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, 
        compute_matrix_multiplication=compute_matrix_multiplication) if nubar == False else \
            np.conj(mixing_matrix_5x5(s12, s23, s13, d13, s14, d14, s15, d15, s24, d24, s25, s34, 
                s35, d35, compute_matrix_multiplication=compute_matrix_multiplication))
    # Mass matrix
    M2 = np.diag([0.0, D21, D31, D41, D51])
    
    return 0.5 * np.linalg.multi_dot([R, M2, np.conj(R.T)])
    # return 0.5 * R @ M2 @ np.conj(R.T)


def hamiltonian_5nu_vacuum_energy_independent_td(l: float, s12: float, s23: float, s13:float, 
    d13: float, s14: float, d14: float, s15: float, d15: float, s24: float, d24: float, s25: float,
    s34: float, s35: float, d35: float, D21: float, D31: float, D41: float, D51: float,
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.
    """
    return hamiltonian_5nu_vacuum_energy_independent(s12, s23, s13, d13, s14, d14, s15, d15, s24, 
        d24, s25, s34, s35, d35, D21, D31, D41, D51, nubar=nubar, 
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_5nu_vacuum(energy: float, s12: float, s23: float, s13:float, d13: float, s14: float,
    d14: float, s15: float, d15: float, s24: float, d24: float, s25: float, s34: float, s35: float, 
    d35: float, D21: float, D31: float, D41: float, D51: float, nubar: Optional[bool]=False, 
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.
    """
    return (1/energy)*hamiltonian_5nu_vacuum_energy_independent(s12, s23, s13, d13, s14, d14, s15,
        d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51, nubar=nubar, 
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_5nu_vacuum_td(l: float, energy: float, s12: float, s23: float, s13:float, 
    d13: float, s14: float, d14: float, s15: float, d15: float, s24: float, d24: float, s25: float,
    s34: float, s35: float, d35: float, D21: float, D31: float, D41: float, D51: float, 
    nubar: Optional[bool]=False, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.
    """
    return hamiltonian_5nu_vacuum(energy, s12, s23, s13, d13, s14, d14, s15, d15, s24, d24, s25, 
        s34, s35, d35, D21, D31, D41, D51, nubar=nubar, 
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_5nu_matter(VCC: float) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 3x3 real three-neutrino Hamiltonian for
    oscillations in matter with constant density.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the three-neutrino Hamiltonian for
        oscillations in vacuum.  This is computed by the routine
        hamiltonian_3nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """
    return np.diag([VCC, 0.0, 0.0, 0.0, 0.0]) 


def hamiltonian_5nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    return hamiltonian_5nu_matter(VCC_func(l))


def hamiltonian_5nu_nsi(
    VCC: float,
    eps_ee: float, 
    eps_em: complex, 
    eps_et: complex, 
    eps_es1: complex, 
    eps_es2: complex, 
    eps_mm: float, 
    eps_mt: complex, 
    eps_ms1: complex, 
    eps_ms2: complex, 
    eps_tt: float,
    eps_ts1: complex, 
    eps_ts2: complex, 
    eps_s1s1: float, 
    eps_s1s2: complex, 
    eps_s2s2: float
) -> np.ndarray:
    r"""Returns the five-neutrino Hamiltonian for oscillations w/ NSI.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations with 
    non-standard interactions (NSI) in matter with constant density.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for oscillations in vacuum.  This is
        computed by the routine hamiltonian_2nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    VCC : float
        Potential due to charged-current interactions of nu_e with electrons.
    eps : list
        Vector of NSI strength parameters: eps = eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """
    return VCC * np.array([
        [eps_ee, eps_em, eps_et, eps_es1, eps_es2], 
        [np.conj(eps_em), eps_mm, eps_mt, eps_ms1, eps_ms2],
        [np.conj(eps_et), np.conj(eps_mt), eps_tt, eps_ts1, eps_ts2],
        [np.conj(eps_es1), np.conj(eps_ms1), np.conj(eps_ts1), eps_s1s1, eps_s1s2],
        [np.conj(eps_es2), np.conj(eps_ms2), np.conj(eps_ts2), np.conj(eps_s1s2), eps_s2s2]
        ], dtype=np.complex128)


def hamiltonian_5nu_liv(energy: float, sxi12: float, sxi23: float, sxi13:float, dxi13: float, 
    sxi14: float, dxi14: float, sxi15: float, dxi15: float, sxi24: float, dxi24: float, 
    sxi25: float, sxi34: float, sxi35: float, dxi35: float, b1: float, b2: float, b3: float, 
    b4: float, b5: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the five-neutrino Hamiltonian for oscillations w/ LIV.
    """

    return pow(energy, n_liv) * hamiltonian_5nu_liv_energy_independent(sxi12, sxi23, sxi13, dxi13,
        sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, dxi35, b1, b2, b3, b4, b5,
        Lambda, n_liv, nubar=nubar, compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_5nu_liv_energy_independent(sxi12: float, sxi23: float, sxi13:float, dxi13: float, 
    sxi14: float, dxi14: float, sxi15: float, dxi15: float, sxi24: float, dxi24: float, 
    sxi25: float, sxi34: float, sxi35: float, dxi35: float, b1: float, b2: float, b3: float, 
    b4: float, b5: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the five-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 5x5 complex four-neutrino Hamiltonian for oscillations in a CPT-odd 
    Lorentz invariance-violating background.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for oscillations in vacuum.  This is
        computed by the routine hamiltonian_2nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    sxi12 : float
        Sin(xi_12), with xi_12 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    sxi23 : float
        Sin(xi_23), with xi_23 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    sxi13 : float
        Sin(xi_12), with xi_13 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    dciCP : float
        CP-violation angle of the LIV operator B3 [radian].
    b1 : float
        Eigenvalue b1 of the LIV operator B3.
    b2 : float
        Eigenvalue b2 of the LIV operator B3.
    b3 : float
        Eigenvalue b3 of the LIV operator B3.
    Lambda : float
        Energy scale of the LIV operator B2.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """
    # 5x5 mixing matrix

    if not nubar:
        R = mixing_matrix_5x5(sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24,
            sxi25, sxi34, sxi35, dxi35, compute_matrix_multiplication=compute_matrix_multiplication)
    else:
        R = np.conj(mixing_matrix_5x5(sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24,
            dxi24, sxi25, sxi34, sxi35, dxi35, 
            compute_matrix_multiplication=compute_matrix_multiplication))
    
    return pow(1.0/Lambda, n_liv) * R @ np.diag([b1, b2, b3, b4, b5]) @ np.conj(R.T)

