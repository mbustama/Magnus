# -*- coding: utf-8 -*-
r"""Compute two-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the two-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * mixing_matrix_2nu - Returns 2x2 rotation matrix
    * hamiltonian_2nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_2nu_matter - Returns H_matter
    * hamiltonian_2nu_nsi - Returns H_NSI
    * hamiltonian_2nu_liv - Returns H_LIV

Created: 2019/04/21 15:00
Last modified: 2019/04/23 21:04
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Callable, Union
# from globaldefs import *


def mixing_matrix_2nu(sth: float) -> np.ndarray:
    r"""Returns the 2x2 rotation matrix.

    Computes and returns a 2x2 real rotation matrix parametrized by a single rotation angle theta.

    Parameters
    ----------
    sth : float
        Sin(theta).

    Returns
    -------
    list
        Rotation matrix [[cth, sth], [-sth, cth]], with cth = cos(theta) and sth = sin(theta).
    """
    cth = np.sqrt(1.0-sth*sth)

    return np.array([[cth,sth],[-sth,cth]])


def hamiltonian_2nu_vacuum_energy_independent(sth: float, Dm2: float, 
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in vacuum,
    parametrized by a single mixing angle theta and a single mass-squared difference Dm2.  The
    Hamiltonian is H = (1/2)*R.M2.R^dagger, with R the 2x2 rotation matrix and M2 the mass matrix.
    The multiplicative factor 1/E is not applied.

    Parameters
    ----------
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    cth = np.sqrt(1.0-sth*sth)
    c2th = cth*cth-sth*sth
    s2th = 2.0*cth*sth

    if not compute_matrix_multiplication:

        H = (Dm2/4.0)*np.array([[c2th,-s2th], [-s2th,-c2th]])

    else:

        # PMNS matrix
        R = mixing_matrix_2nu(sth)
        # Mass matrix
        M2 = np.array([[1.0, 0.0], [0.0, -1.0]])
        # Hamiltonian
        H = (Dm2/4.0)*np.matmul(R, np.matmul(M2, np.transpose(R)))

    return H


def hamiltonian_2nu_vacuum_energy_independent_td(l: float, sth: float, Dm2: float, 
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in vacuum, as a 
    function of distance, parametrized by a single mixing angle theta and a single mass-squared 
    difference Dm2.  The Hamiltonian is H = (1/2)*R.M2.R^dagger, with R the 2x2 rotation matrix and
    M2 the mass matrix.  The multiplicative factor 1/E is not applied.  The vacuum Hamiltonian does
    not depend on distance in reality, but we include the dependence here as a way to validate the
    routine to compute probabilities for time-dependent Hamiltonians.

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger 
        live.

    Returns 
    -------
    list
        Hamiltonian 2x2 matrix.
    """

    H = hamiltonian_2nu_vacuum_energy_independent(sth, Dm2, 
        compute_matrix_multiplication=compute_matrix_multiplication)

    return H


def hamiltonian_2nu_vacuum(energy: float, sth: float, Dm2: float, 
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:

    h_vac = hamiltonian_2nu_vacuum_energy_independent(sth, Dm2, 
        compute_matrix_multiplication=compute_matrix_multiplication)

    return (1/energy)*h_vac


def hamiltonian_2nu_vacuum_td(l: float, energy: float, sth: float, Dm2: float, 
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in vacuum, as a
    function of distance, parametrized by a single mixing angle theta and a single mass-squared
    difference Dm2.  The Hamiltonian is H = (1/2)*R.M2.R^dagger, with R the 2x2 rotation matrix and 
    M2 the mass matrix.  The multiplicative factor 1/E is not applied.  The vacuum Hamiltonian does 
    not depend on distance in reality, but we include the dependence here as a way to validate the 
    routine to compute probabilities for time-dependent Hamiltonians.

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns 
    -------
    list
        Hamiltonian 2x2 matrix.
    """

    H = hamiltonian_2nu_vacuum(energy, sth, Dm2, 
        compute_matrix_multiplication=compute_matrix_multiplication)

    return H


def hamiltonian_2nu_matter(VCC: float) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in matter with
    constant density.

    Parameters
    ----------
    energy : float
        Neutrino energy.
    VCC : float
        Potential due to charged-current interactions of nu_e with electrons.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    h_matter = np.zeros((2,2))

    # Add the matter potential to the ee term to find the matter Hamiltonian
    h_matter[0][0] = VCC

    return h_matter


def hamiltonian_2nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for matter oscillations, as a function of distance.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in matter with a given density as a function of
    position.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for oscillations in vacuum.  This is
        computed by the routine hamiltonian_2nu_vacuum_energy_independent.
    VCC_func : float
        Potential due to charged-current interactions of nu_e with electrons. This is a function 
        only of the position, l.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    h_matter = hamiltonian_2nu_matter(VCC_func(l))

    return h_matter


def hamiltonian_2nu_nsi(VCC: float, eps: Union[list, np.ndarray]) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with NSI.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations with non-standard 
    interactions (NSI) in matter with constant density.

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
        Vector of NSI strength parameters: eps = eps_ee, eps_em, eps_mm.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    h_nsi = np.zeros((2,2))

    eps_ee, eps_em, eps_mm = eps

    h_nsi[0][0] = VCC*(1.0+eps_ee)
    h_nsi[0][1] = VCC*eps_em
    h_nsi[1][0] = VCC*np.conj(eps_em)
    h_nsi[1][1] = VCC*eps_mm

    return h_nsi


def hamiltonian_2nu_nsi_td(l: float, VCC_func: Callable, 
    eps: Union[list, np.ndarray]) -> np.ndarray:

    return hamiltonian_2nu_nsi(VCC_func(l), eps)


def hamiltonian_2nu_liv(energy: float, sxi: float, b1: float, b2: float, 
    Lambda: float) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with LIV.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in a CPT-odd Lorentz
    invariance-violating background.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for oscillations in vacuum.  This is
        computed by the routine hamiltonian_2nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    sxi : float
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of B2 and the 
        flavor states.
    b1 : float
        Eigenvalue b1 of the LIV operator B2.
    b2 : float
        Eigenvalue b2 of the LIV operator B2.
    Lambda : float
        Energy scale of the LIV operator B2.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    h_liv = np.zeros((2,2))

    cxi = np.sqrt(1.0-sxi-sxi)
    h_liv[0][0] = (b1*cxi*cxi + b2*sxi*sxi)
    h_liv[0][1] = ((-b1+b2)*cxi*sxi)
    h_liv[1][0] = ((-b1+b2)*cxi*sxi)
    h_liv[1][1] = (b2*cxi*cxi + b1*sxi*sxi)
    h_liv = (energy/Lambda)*h_liv

    return h_liv
