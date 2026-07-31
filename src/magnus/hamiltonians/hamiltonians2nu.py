# -*- coding: utf-8 -*-
r"""hamiltonians2nu.py

Compute two-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the two-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * mixing_matrix_2nu - Returns 2x2 rotation matrix
    * hamiltonian_2nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_2nu_vacuum_energy_independent_td - Returns H_vac (no
           1/E), as a function of position
    * hamiltonian_2nu_vacuum - Returns H_vac
    * hamiltonian_2nu_vacuum_td - Returns H_vac, as a function of position
    * hamiltonian_2nu_matter - Returns H_matter
    * hamiltonian_2nu_matter_td - Returns H_matter, as a function of position
    * hamiltonian_2nu_nsi - Returns H_NSI
    * hamiltonian_2nu_nsi_td - Returns H_NSI, as a function of position
    * hamiltonian_2nu_liv - Returns H_LIV
    * hamiltonian_2nu_liv_energy_independent - Returns H_LIV (no energy
           dependence)
"""


__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Callable
# from globaldefs import *


def mixing_matrix_2nu(sth: float) -> np.ndarray:
    r"""Returns the 2x2 rotation matrix.

    Computes and returns a 2x2 real rotation matrix parametrized by a single rotation angle theta.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    sth : float
        Sine of the mixing angle :math:`\theta`.

    Returns
    -------
    np.ndarray
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

    .. versionadded:: 1.0.0

    Parameters
    ----------
    sth : float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    cth = np.sqrt(1.0-sth*sth)
    c2th = cth*cth-sth*sth
    s2th = 2.0*cth*sth

    if not compute_matrix_multiplication:

        # Standard convention, Dm2 = m2^2 - m1^2 > 0 for m2 > m1 (matches the
        # three-neutrino convention M2 = diag(0, D21, D31) and the closed-form
        # matter probability in oscprobstd): H = (Dm2/4E) [[-c2th, s2th], [s2th, c2th]]
        return (Dm2/4.0)*np.array([[-c2th,s2th], [s2th,c2th]])

    else:

        # 2D mixing matrix
        R = mixing_matrix_2nu(sth)
        # Mass matrix
        M2 = np.diag([-1.0, 1.0]) # (m1^2 - m2^2, m2^2 - m1^2)/Dm2, traceless form
        # Hamiltonian
        return (Dm2/4.0) * R @ M2 @ R.T  # Use matrix multiplication operator
        # H = (Dm2/4.0)*np.matmul(R, np.matmul(M2, np.transpose(R)))

    # return H


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

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    sth : float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger
        live.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """

    return hamiltonian_2nu_vacuum_energy_independent(sth, Dm2,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_2nu_vacuum(energy: float, sth: float, Dm2: float,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations.

    Same as :func:`hamiltonian_2nu_vacuum_energy_independent`, but with the 1/E factor applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float
        Neutrino energy.
    sth : float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    return (1/energy)*hamiltonian_2nu_vacuum_energy_independent(sth, Dm2,
        compute_matrix_multiplication=compute_matrix_multiplication)


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

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    energy : float
        Neutrino energy.
    sth : float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """

    return hamiltonian_2nu_vacuum(energy, sth, Dm2,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_2nu_matter(VCC: float) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in matter with
    constant density.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    VCC : float
        Potential due to charged-current interactions of nu_e with electrons.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    # The matter Hamiltonian is [[VCC,0],[0,0]]
    return np.diag([VCC, 0.0])
    # h_matter = np.zeros((2,2))

    # # Add the matter potential to the ee term to find the matter Hamiltonian
    # h_matter[0][0] = VCC

    # return h_matter


def hamiltonian_2nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for matter oscillations, as a function of distance.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in matter with a given density as a function of
    position.

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
        Hamiltonian 2x2 matrix.
    """
    return hamiltonian_2nu_matter(VCC_func(l))


def hamiltonian_2nu_nsi(VCC: float, eps_aa: float, eps_ab: complex) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with NSI.

    Computes and returns the 2x2 complex two-neutrino Hamiltonian for oscillations with
    non-standard interactions (NSI) in matter with constant density.

    Two flavors admit only one physically meaningful diagonal NSI degree of freedom: an overall
    (flavor-universal) diagonal shift is proportional to the identity matrix, so it commutes with
    everything, contributes only an unobservable global phase to the evolution operator, and cannot
    affect any oscillation probability.  ``eps_aa`` is therefore defined here as the non-universal
    (flavor-off-diagonal-*difference*) coupling, following the convention eps_mumu = 0, i.e., it
    parametrizes the coupling of :math:`\nu_e` alone, relative to :math:`\nu_\mu`.  [Earlier
    versions of this function placed eps_aa on *both* diagonal entries, making it a pure multiple
    of the identity and therefore a no-op on every oscillation probability -- this was a bug, not a
    convention choice, confirmed by direct calculation.]

    .. versionadded:: 1.0.0

    Parameters
    ----------
    VCC : float
        Potential due to charged-current interactions of nu_e with electrons.
    eps_aa : float
        Non-universal diagonal NSI coupling of nu_e (relative to nu_mu, whose diagonal coupling is
        fixed to 0 by this convention).
    eps_ab : complex
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    return VCC * np.array([[eps_aa, eps_ab], [np.conj(eps_ab), 0.0]], dtype=np.complex128)


def hamiltonian_2nu_nsi_td(l: float, VCC_func: Callable, eps_aa: float,
    eps_ab: complex) -> np.ndarray:
    r"""Returns the two-neutrino NSI Hamiltonian as a function of position.

    Same as :func:`hamiltonian_2nu_nsi`, but evaluates the position-dependent matter potential
    ``VCC_func(l)`` first.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    VCC_func : Callable
        Potential due to charged-current interactions of nu_e with electrons, as a function of
        position, l.
    eps_aa : float
        Non-universal diagonal NSI coupling of nu_e; see :func:`hamiltonian_2nu_nsi`.
    eps_ab : complex
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    return hamiltonian_2nu_nsi(VCC_func(l), eps_aa, eps_ab)


def hamiltonian_2nu_liv(energy: float, sxi: float, b1: float, b2: float, Lambda: float, n_liv: int,
    nubar: Optional[bool]=False) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with LIV.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in a CPT-odd Lorentz
    invariance-violating background.  Same as
    :func:`hamiltonian_2nu_liv_energy_independent`, but with the
    :math:`E^{n_{\rm liv}}` energy dependence of the LIV operator applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
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
    n_liv : int
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3).
    nubar : bool, optional
        Accepted for interface parity with :func:`magnus.hamiltonians.hamiltonians3nu.hamiltonian_3nu_liv` and its
        4nu/5nu siblings, which conjugate their (complex) LIV mixing matrix for antineutrinos.  The
        2-flavor LIV rotation has no CP-violating phase (only the real angle ``sxi``), so there is
        nothing to conjugate and this parameter currently has no effect. Default: False.

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    return pow(energy, n_liv) * hamiltonian_2nu_liv_energy_independent(sxi, b1, b2, Lambda, n_liv)


def hamiltonian_2nu_liv_energy_independent(sxi: float, b1: float, b2: float,
    Lambda: float, n_liv: int) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with LIV.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for oscillations in a CPT-odd Lorentz
    invariance-violating background, without the energy-dependent prefactor.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    sxi : float
        Sin(xi), with xi the rotation angle between the space of the eigenvectors of B2 and the
        flavor states.
    b1 : float
        Eigenvalue b1 of the LIV operator B2.
    b2 : float
        Eigenvalue b2 of the LIV operator B2.
    Lambda : float
        Energy scale of the LIV operator B2.
    n_liv : int
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3).

    Returns
    -------
    np.ndarray
        Hamiltonian 2x2 matrix.
    """
    # H = R . diag(b1, b2) . R^T, with R = mixing_matrix_2nu(sxi) -- the same convention used by
    # every sibling Hamiltonian (2nu vacuum's slow path, and the 3/4/5nu LIV Hamiltonians).  The
    # off-diagonal sign below was previously flipped relative to this convention (a confirmed bug).
    cxi = np.sqrt(1.0 - sxi * sxi)
    delta_b = b2 - b1

    return pow(1.0 / Lambda, n_liv) * np.array([
        [b1 * cxi * cxi + b2 * sxi * sxi, delta_b * cxi * sxi],
        [delta_b * cxi * sxi, b2 * cxi * cxi + b1 * sxi * sxi]
    ], dtype=np.float64)
