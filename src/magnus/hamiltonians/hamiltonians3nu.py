# -*- coding: utf-8 -*-
r"""hamiltonians3nu.py

Compute three-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the three-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * pmns_mixing_matrix - Returns the 3x3 PMNS mixing matrix
    * mixing_matrix_3x3 - Alias of pmns_mixing_matrix
    * hamiltonian_3nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_3nu_vacuum_energy_independent_td - Returns H_vac (no
           1/E), as a function of position
    * hamiltonian_3nu_vacuum - Returns H_vac
    * hamiltonian_3nu_vacuum_td - Returns H_vac, as a function of position
    * hamiltonian_3nu_matter - Returns H_matter
    * hamiltonian_3nu_matter_td - Returns H_matter, as a function of position
    * hamiltonian_3nu_nsi - Returns H_NSI
    * hamiltonian_3nu_nsi_td - Returns H_NSI, as a function of position
    * hamiltonian_3nu_liv - Returns H_LIV
    * hamiltonian_3nu_liv_energy_independent - Returns H_LIV (no energy
           dependence)
"""


__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# from numpy import *
import numpy as np
from typing import Optional, Callable

# import cmath
# import cmath as cmath
# import copy as cp

# import oscprob3nu
# from globaldefs import *


def pmns_mixing_matrix(s12: float, s23: float, s13:float, dCP: float) -> np.ndarray:
    r"""Returns the 3x3 PMNS mixing matrix.

    Computes and returns the 3x3 complex PMNS mixing matrix parametrized by three rotation angles,
    :math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}`, and one CP-violation phase, :math:`\delta_\text{CP}`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : float
        :math:`\delta_\text{CP}` [radian].

    Returns
    -------
    np.ndarray
        3x3 PMNS mixing matrix.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus.hamiltonians import hamiltonians3nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        U = np.asarray(hamiltonians3nu.pmns_mixing_matrix(
            p['s12'], p['s23'], p['s13'], p['dCP']))

        print('|U_e2|^2 = %.4f   (sin^2 th12 = %.3f)' % (abs(U[0][1])**2,
                                                         p['s12']**2))
        print('unitary to %.1e' % np.max(np.abs(U.conj().T @ U - np.eye(3))))
"""
    c12 = np.sqrt(1.0-s12*s12)
    c23 = np.sqrt(1.0-s23*s23)
    c13 = np.sqrt(1.0-s13*s13)
    cdCP = np.cos(dCP)
    # sdCP = np.sqrt(1.0-cdCP*cdCP)
    sdCP = np.sin(dCP)
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


def mixing_matrix_3x3(s12: float, s23: float, s13:float, dCP: float) -> np.ndarray:
    r"""Returns the 3x3 PMNS mixing matrix.

    Alias of :func:`pmns_mixing_matrix`, kept for naming parity with
    :func:`magnus.hamiltonians.hamiltonians4nu.mixing_matrix_4x4` and
    :func:`magnus.hamiltonians.hamiltonians5nu.mixing_matrix_5x5`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : float
        :math:`\delta_\text{CP}` [radian].

    Returns
    -------
    np.ndarray
        3x3 PMNS mixing matrix.
    """
    return pmns_mixing_matrix(s12, s23, s13, dCP)


def hamiltonian_3nu_vacuum_energy_independent(s12: float, s23: float, s13: float, dCP: float,
    D21: float, D31: float, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations in vacuum,
    parametrized by three mixing angles (:math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}`), one CP-violation phase
    (:math:`\delta_\text{CP}`), and two mass-squared difference (:math:`\Delta m_{21}^2`, :math:`\Delta m_{31}^2`).  The Hamiltonian is
    H = (1/2)*R.M2.R^dagger, with R the 3x3 PMNS matrix and M2 the mass matrix.  The multiplicative
    factor 1/E is not applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : float
        :math:`\delta_\text{CP}` [radian].
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos (conjugates the PMNS matrix, equivalent
        to :math:`\delta_\text{CP}` -> -:math:`\delta_\text{CP}`). Default: False.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger
        live.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import magnus.globaldefs as gd
        from magnus.hamiltonians import hamiltonians3nu

        p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31']))

        print('Hermitian to %.1e' % np.max(np.abs(h - h.conj().T)))
        print('eigenvalues [eV^2]:', np.round(np.linalg.eigvalsh(h), 6))

    The eigenvalues are :math:`(0, \Delta m^2_{21}, \Delta m^2_{31})`: only
    mass-squared *differences* appear, which is why the first is zero.
"""

    # f = 0.5

    if not compute_matrix_multiplication:

        c12 = np.sqrt(1.0-s12*s12)
        c23 = np.sqrt(1.0-s23*s23)
        c13 = np.sqrt(1.0-s13*s13)
        cdCP = np.cos(dCP)
        # sdCP = np.sqrt(1.0-cdCP*cdCP)
        sdCP = np.sin(dCP)
        # exp_dCP_p = complex(cdCP, sdCP)
        exp_dCP_p = complex(cdCP, sdCP) if not nubar else complex(cdCP, -sdCP)
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
        # if not nubar:
        #     R = pmns_mixing_matrix(s12, s23, s13, dCP)
        # else:
        #     R = np.conj(pmns_mixing_matrix(s12, s23, s13, dCP))
        R = pmns_mixing_matrix(s12, s23, s13, dCP) if not nubar \
                else np.conj(pmns_mixing_matrix(s12, s23, s13, dCP))
        # Mass matrix
        M2 = np.diag([0.0, D21, D31])
        # Hamiltonian
        return 0.5 * R @ M2 @ np.conj(R.T)


def hamiltonian_3nu_vacuum_energy_independent_td(l: float, s12: float, s23: float, s13: float,
    dCP: float, D21: float, D31: float,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Same as :func:`hamiltonian_3nu_vacuum_energy_independent`, included for interface parity with
    the other, genuinely position-dependent Hamiltonians (see, e.g., :func:`hamiltonian_3nu_matter_td`).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : float
        :math:`\delta_\text{CP}` [radian].
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger
        live.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    """
    return hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_vacuum(energy: float, s12: float, s23: float, s13: float, dCP: float,
    D21: float, D31: float, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.

    Same as :func:`hamiltonian_3nu_vacuum_energy_independent`, but with the 1/E factor applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float
        Neutrino energy.
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : float
        :math:`\delta_\text{CP}` [radian].
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos. Default: False.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger
        live.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    """
    return (1/energy)*hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31,
        nubar=nubar, compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_vacuum_td(l: float, energy: float, s12: float, s23: float, s13: float, dCP: float,
    D21: float, D31: float, compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations, as a function of distance,
    even if it does not depend on it.

    Same as :func:`hamiltonian_3nu_vacuum`, included for interface parity with the other,
    genuinely position-dependent Hamiltonians.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    energy : float
        Neutrino energy.
    s12 : float
        Sine of the mixing angle :math:`\theta_{12}`.
    s23 : float
        Sine of the mixing angle :math:`\theta_{23}`.
    s13 : float
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP : float
        :math:`\delta_\text{CP}` [radian].
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise, multiply R.M2.R^dagger
        live.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    """
    return hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_matter(VCC: float) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 3x3 real three-neutrino Hamiltonian for
    oscillations in matter with constant density.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus.hamiltonians import hamiltonians3nu

        print(np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(1.0e-13)))

    Add it to the vacuum term divided by the energy to get the full
    Hamiltonian.  For antineutrinos the potential arrives already negated by
    :func:`magnus.matter.vcc_func_from_rho_func`; do not negate it again.
"""
    # Built by broadcasting rather than np.diag so that VCC may be an array of
    # positions: VCC[..., None, None] turns one potential per position into a
    # stack of matrices, which is what lets a caller's H_func take the engine's
    # vectorized path (see magnus.magnus.ScalarHamiltonianWarning). A scalar VCC
    # still returns a plain (3, 3) matrix.
    VCC = np.asarray(VCC, dtype=float)
    e00 = np.zeros((3, 3))
    e00[0, 0] = 1.0
    return VCC[..., None, None] * e00


def hamiltonian_3nu_matter_td(l: float, VCC_func: Callable) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for matter oscillations, as a function of distance.

    Computes and returns the 3x3 real three-neutrino Hamiltonian for oscillations in matter with a
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
        Hamiltonian 3x3 matrix.
    """
    return hamiltonian_3nu_matter(VCC_func(l))


def hamiltonian_3nu_nsi(
    VCC: float,
    eps_ee: float,
    eps_em: complex,
    eps_et: complex,
    eps_mm: float,
    eps_mt: complex,
    eps_tt: float
) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ NSI.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations with
    non-standard interactions (NSI) in matter with constant density.

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
    eps_mm : float
        Diagonal NSI coupling of nu_mu.
    eps_mt : complex
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling.
    eps_tt : float
        Diagonal NSI coupling of nu_tau.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus.hamiltonians import hamiltonians3nu

        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_nsi(
            1.0e-13, 0.1, 0.05, 0.0, 0.0, 0.0, 0.0))

        print(np.round(h/1e-13, 4), ' [1e-13 eV]')

    The couplings are dimensionless and multiply the same :math:`V_{CC}`, so
    ``eps_ee = 0.1`` is a ten-per-cent correction to the standard potential and
    ``eps_em`` is an off-diagonal one the Standard Model does not have.
"""
    return VCC * np.array([
        [eps_ee, eps_em, eps_et],
        [np.conj(eps_em), eps_mm, eps_mt],
        [np.conj(eps_et), np.conj(eps_mt), eps_tt],
        ], dtype=np.complex128)


def hamiltonian_3nu_nsi_td(l: float, VCC_func: Callable, eps_ee: float, eps_em: complex,
    eps_et: complex, eps_mm: float, eps_mt: complex, eps_tt: float) -> np.ndarray:
    r"""Returns the three-neutrino NSI Hamiltonian as a function of position.

    Same as :func:`hamiltonian_3nu_nsi`, but evaluates the position-dependent matter potential
    ``VCC_func(l)`` first.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    VCC_func : Callable
        Potential due to charged-current interactions of nu_e with electrons, as a function of
        position, l.
    eps_ee : float
        Diagonal NSI coupling of nu_e.
    eps_em : complex
        Flavor-off-diagonal (nu_e-nu_mu) NSI coupling.
    eps_et : complex
        Flavor-off-diagonal (nu_e-nu_tau) NSI coupling.
    eps_mm : float
        Diagonal NSI coupling of nu_mu.
    eps_mt : complex
        Flavor-off-diagonal (nu_mu-nu_tau) NSI coupling.
    eps_tt : float
        Diagonal NSI coupling of nu_tau.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    """
    return hamiltonian_3nu_nsi(VCC_func(l), eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt)


def hamiltonian_3nu_liv(energy: float, sxi12: float, sxi23: float, sxi13: float, dxiCP: float, b1: float,
    b2: float, b3: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations in a CPT-odd
    Lorentz invariance-violating background.  Same as
    :func:`hamiltonian_3nu_liv_energy_independent`, but with the
    :math:`E^{n_{\rm liv}}` energy dependence of the LIV operator applied.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energy : float
        Neutrino energy.
    sxi12 : float
        Sin(xi_12), with xi_12 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    sxi23 : float
        Sin(xi_23), with xi_23 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    sxi13 : float
        Sin(xi_13), with xi_13 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    dxiCP : float
        CP-violation angle of the LIV operator B3 [radian].
    b1 : float
        Eigenvalue b1 of the LIV operator B3.
    b2 : float
        Eigenvalue b2 of the LIV operator B3.
    b3 : float
        Eigenvalue b3 of the LIV operator B3.
    Lambda : float
        Energy scale of the LIV operator B3.
    n_liv : int
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3).
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos (conjugates the LIV mixing matrix,
        equivalent to dxiCP -> -dxiCP). Default: False.
    compute_matrix_multiplication : bool, optional
        Forwarded to :func:`hamiltonian_3nu_liv_energy_independent` (currently unused there; kept
        for interface parity with the vacuum Hamiltonian).

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    """
    return pow(energy, n_liv) * hamiltonian_3nu_liv_energy_independent(sxi12, sxi23, sxi13, dxiCP,
        b1, b2, b3, Lambda, n_liv, nubar=nubar,
        compute_matrix_multiplication=compute_matrix_multiplication)


def hamiltonian_3nu_liv_energy_independent(sxi12: float, sxi23: float, sxi13: float, dxiCP: float,
    b1: float, b2: float, b3: float, Lambda: float, n_liv: int, nubar: Optional[bool]=False,
    compute_matrix_multiplication: Optional[bool]=False) -> np.ndarray:
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for oscillations in a CPT-odd
    Lorentz invariance-violating background, without the energy-dependent prefactor.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    sxi12 : float
        Sin(xi_12), with xi_12 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    sxi23 : float
        Sin(xi_23), with xi_23 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    sxi13 : float
        Sin(xi_13), with xi_13 the one of the mixing angles between the space of the eigenvectors of
        B3 and the flavor states.
    dxiCP : float
        CP-violation angle of the LIV operator B3 [radian].
    b1 : float
        Eigenvalue b1 of the LIV operator B3.
    b2 : float
        Eigenvalue b2 of the LIV operator B3.
    b3 : float
        Eigenvalue b3 of the LIV operator B3.
    Lambda : float
        Energy scale of the LIV operator B3.
    n_liv : int
        Power of the energy dependence of the LIV operator (dimension of the operator minus 3);
        enters here through the :math:`\Lambda^{-n_{\rm liv}}` normalization of the eigenvalues.
    nubar : bool, optional
        If True, compute the Hamiltonian for antineutrinos (conjugates the LIV mixing matrix,
        equivalent to dxiCP -> -dxiCP). Default: False.
    compute_matrix_multiplication : bool, optional
        Currently unused; accepted for interface parity with the vacuum Hamiltonian.

    Returns
    -------
    np.ndarray
        Hamiltonian 3x3 matrix.
    """
    R = pmns_mixing_matrix(sxi12, sxi23, sxi13, dxiCP) if not nubar \
            else np.conj(pmns_mixing_matrix(sxi12, sxi23, sxi13, dxiCP))

    return pow(1.0/Lambda, n_liv) * R @ np.diag([b1, b2, b3]) @ np.conj(R.T)
