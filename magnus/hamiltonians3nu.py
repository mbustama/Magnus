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

# import cmath
# import cmath as cmath
# import copy as cp

# import oscprob3nu
# from globaldefs import *


def pmns_mixing_matrix(s12: float, s23: float, s13:float, dCP: float) -> np.ndarray:
    r"""Returns the 3x3 PMNS mixing matrix.

    Computes and returns the 3x3 complex PMNS mixing matrix parametrized by three rotation angles, 
    theta_12, theta_23, theta_13, and one CP-violation phase, delta_CP.

    Parameters
    ----------
    s12 : float
        Sin(theta_12).
    s23 : float
        Sin(theta_23).
    s13 : float
        Sin(theta_13).
    dCP : float
        delta_CP [radian].

    Returns
    -------
    list
        3x3 PMNS mixing matrix.
    """
    c12 = np.sqrt(1.0-s12*s12)
    c23 = np.sqrt(1.0-s23*s23)
    c13 = np.sqrt(1.0-s13*s13)
    cdCP = np.cos(dCP)
    sdCP = np.sqrt(1.0-cdCP*cdCP)
    exp_dCP_p = complex(cdCP, sdCP)
    exp_dCP_m = np.conj(exp_dCP_p)

    U00 = c12*c13
    U01 = s12*c13
    U02 = s13*exp_dCP_m
    U10 = -s12*c23 - c12*s23*s13*exp_dCP_p
    U11 = c12*c23 - s12*s23*s13*exp_dCP_p
    U12 = s23*c13
    U20 = s12*s23 - c12*c23*s13*exp_dCP_p
    U21 = -c12*s23 - s12*c23*s13*exp_dCP_p
    U22 = c23*c13

    return np.array([[U00,U01,U02],[U10,U11,U12],[U20,U21,U22]])


def hamiltonian_3nu_vacuum_energy_independent(s12: float, s23: float, s13: float, dCP: float, 
    D21: float, D31: float, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.

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

    # f = 0.5

    if not compute_matrix_multiplication:

        c12 = np.sqrt(1.0-s12*s12)
        c23 = np.sqrt(1.0-s23*s23)
        c13 = np.sqrt(1.0-s13*s13)
        cdCP = np.cos(dCP)
        sdCP = np.sqrt(1.0-cdCP*cdCP)
        exp_dCP_p = complex(cdCP, sdCP)
        exp_dCP_m = np.conj(exp_dCP_p)

        # All Hij have units of [eV^2]
        H00 = c13*c13*D21*s12*s12 + D31*s13*s13
        H01 = c12*c13*c23*D21*s12 + c13*(D31-D21*s12*s12)*s13*s23*exp_dCP_m
        H02 = c13*c23*(D31-D21*s12*s12)*s13*exp_dCP_m - c12*c13*D21*s12*s23
        H10 = c12*c13*c23*D21*s12 + c13*(D31-D21*s12*s12)*s13*s23*exp_dCP_p
        H11 = c12*c12*c23*c23*D21 + (c13*c13*D31 + D21*s12*s12*s13*s13)*s23*s23 - \
                2.0*c12*c23*D21*s12*s13*s23*cdCP
        H12 = c13*c13*c23*D31*s23 + (c23*s12*s13*exp_dCP_m + c12*s23) * \
                (-c12*c23*D21 + D21*s12*s13*s23*exp_dCP_p)
        H20 = c13*c23*(D31-D21*s12*s12)*s13*exp_dCP_p - c12*c13*D21*s12*s23
        H21 = c13*c13*c23*D31*s23 - D21*(c23*s12*s13*exp_dCP_p + c12*s23) * \
                (c12*c23 - s12*s13*s23*exp_dCP_m)
        H22 = c23*c23*(c13*c13*D31 + D21*s12*s12*s13*s13) + c12*c12*D21*s23*s23 + \
                2.0*c12*c23*D21*s12*s13*s23*cdCP

        return 0.5*np.array([[H00,H01,H02], [H10,H11,H12], [H20,H21,H22]])

    else:

        # PMNS matrix
        R = pmns_mixing_matrix(s12, s23, s13, dCP)
        # Mass matrix
        M2 = np.diag([0.0, D21, D31])
        # Hamiltonian
        return 0.5 * R @ M2 @ R.T 


def hamiltonian_3nu_vacuum_energy_independent_td(l: float, s12: float, s23: float, s13: float, 
    dCP: float, D21: float, D31: float, 
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.
    """
    return hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_vacuum(energy: float, s12: float, s23: float, s13: float, dCP: float, 
    D21: float, D31: float, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.
    """
    return (1/energy)*hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_vacuum(l: float, energy: float, s12: float, s23: float, s13: float, dCP: float, 
    D21: float, D31: float, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.
    """
    return hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, 
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_matter(VCC: float) -> np.ndarray:
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
    return np.diag([VCC, 0.0, 0.0]) 


def hamiltonian_3nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    return hamiltonian_3nu_matter(VCC_func(l))


def hamiltonian_3nu_nsi(VCC: float, eps: Union[list, np.ndarray]) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ NSI.

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
    eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = eps
    return VCC * np.array([
        [1.0+eps_ee, eps_em, eps_et], 
        [np.conj(eps_em), eps_mm, eps_mt],
        [np.conj(eps_et), np.conj(eps_mt), eps_tt],
        ], dtype=np.complex128)


def hamiltonian_3nu_nsi_td(l: float, VCC_func: Callable, 
    eps: Union[list, np.ndarray]) -> np.ndarray:

    return hamiltonian_3nu_nsi(VCC_func(l), eps)


def hamiltonian_3nu_liv(sxi12: float, sxi23: float, sxi13: float, dxiCP: float, b1: float, 
    b2: float, b3: float, Lambda: float) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations in a CPT-odd 
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
    R = pmns_mixing_matrix(sxi12, sxi23, sxi13, dxiCP)
    return (energy/Lambda) * R @ np.diag([b1, b2, b3]) @ np.conj(R.T)
