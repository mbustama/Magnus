# -*- coding: utf-8 -*-
r"""oscprobstd.py

Closed-form (non-Magnus) two- and three-neutrino oscillation
probabilities, computed from the standard analytical expressions rather
than from the Magnus expansion.  Used by the test suite to validate
oscprob.py's Magnus-based results against an independent method, not
intended as a general-purpose replacement for it (it does not cover
matter with non-constant density, NSI, LIV, or more than three flavors).

Routine listings
----------------

    * osc_prob_2nu_vacuum_std - Returns 2nu vacuum probabilities, closed form
    * osc_prob_2nu_matter_std - Returns 2nu constant-density matter
           probabilities, closed form
    * delta - Kronecker delta
    * J - Returns U*_ak * U_bk * U_aj * U*_bj, a building block of the
           3nu vacuum probability
    * osc_prob_3nu_vacuum_std - Returns 3nu vacuum probabilities, closed form
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Union


def osc_prob_2nu_vacuum_std(sth: float, Dm2: float, energy: float, L: float) -> np.ndarray:
    r"""Returns 2nu oscillation vacuum probabilities, std. computation.

    Returns the probabilities for two-neutrino oscillations in vacuum, computed using the standard
    analytical expression of the probabilities.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    sth : float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2`.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    np.ndarray
        List of probabilities [Pee, Pem, Pme, Pmm].
    """
    # arg = 1.27*Dm2*L/energy#/4.0
    cth = np.sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth

    Pem = s2th*s2th * pow(np.sin(Dm2*L/energy/4.0), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    prob = np.array([[Pee, Pem], [Pme, Pmm]])

    return prob


def osc_prob_2nu_matter_std(sth: float, Dm2: float, VCC: float, energy: float, 
    L: float) -> np.ndarray:
    r"""Returns 2nu oscillation matter probabilities, std. computation.

    Returns the probabilities for two-neutrino oscillations in matter, computed using the standard
    analytical expression of the probabilities.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    sth : float
        Sine of the mixing angle :math:`\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2`.
    VCC : float
        Potential due to charged-current interactions of :math:`\nu_e` with electrons.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    np.ndarray
        List of probabilities [Pee, Pem, Pme, Pmm].
    """
    # x = 2.0*VCC*(energy*1.e9)/Dm2
    x = 2.0*VCC*(energy)/Dm2
    cth = np.sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth
    s2thsq = s2th*s2th
    c2th = np.sqrt(1.0-s2thsq)

    Dm2m = Dm2*np.sqrt(s2thsq+pow(c2th-x, 2.0))
    s2thmsq = s2thsq / (s2thsq+pow(c2th-x, 2.0))

    # arg = 1.27*Dm2m*L/energy#/4.0
    Pem = s2thmsq * pow(np.sin(Dm2m*L/energy/4.0), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    prob = np.array([[Pee, Pem], [Pme, Pmm]])

    return prob


def delta(a: int, b: int) -> int:
    r"""Returns the Kronecker delta function.

    Returns the delta function delta(a, b) = 1 if a == b and 0 if a != b.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    a : int
        First index.
    b : int
        Second index.

    Returns
    -------
    int
        delta(a, b).
    """
    # if (a == b):
    #     return 1
    # else:
    #     return 0
    return 1 if (a == b) else 0


def J(U: Union[list, np.ndarray], alpha: int, beta: int, k: int, j: int) -> complex:
    r"""Returns U*_ak * U_bk * U_aj * U*_bj, with U the PMNS matrix.

    Returns the product U*_ak * U_bk * U_aj * U*_bj, where U is the PMNS mixing matrix.  This
    product appears in the standard expression for the three-neutrino oscillation probability in
    vacuum.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    U : list or np.ndarray
        3x3 PMNS complex mixing matrix.
    alpha : int
        Index of the initial flavor (0: e, 1: mu, 2: tau).
    beta : int
        Index of the final flavor (0: e, 1: mu, 2: tau).
    k : int
        First index of the sum over mass eigenstates (k = 0, 1, 2).
    j : int
        Second index of the sum over mass eigenstates (j = 0, 1, 2).

    Returns
    -------
    complex
        J(U, alpha, beta, k, j)
    """
    return np.conj(U[alpha][k])*U[beta][k]*U[alpha][j]*np.conj(U[beta][j])


def osc_prob_3nu_vacuum_std(U: Union[list, np.ndarray], D21: float, D31: float, energy: float, 
    L: float, nubar: Optional[bool]=False) -> np.ndarray:
    r"""Returns 3nu oscillation vacuum probabilities, std. computation.

    Returns the probabilities for three-neutrino oscillations in vacuum, computed using the standard
    analytical expression of the probabilities.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    U : list or np.ndarray
        3x3 PMNS complex mixing matrix.
    D21 : float
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31 : float
        Mass-squared difference :math:`\Delta m_{31}^2`.
    energy : float
        Neutrino energy.
    L : float
        Baseline.
    nubar : bool, optional
        If True, compute the probability for antineutrinos (flips the
        sign of the CP-violating term). Default: False.

    Returns
    -------
    np.ndarray
        List of probabilities [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt].
    """
    D32 = D31-D21
    # arg21 = 2.54*D21*L/energy
    # arg31 = 2.54*D31*L/energy
    # arg32 = 2.54*D32*L/energy
    arg21 = D21*L/energy/2.0
    arg31 = D31*L/energy/2.0
    arg32 = D32*L/energy/2.0
    s21 = np.sin(arg21)
    s31 = np.sin(arg31)
    s32 = np.sin(arg32)
    ss21 = pow(np.sin(arg21/2.0), 2.0)
    ss31 = pow(np.sin(arg31/2.0), 2.0)
    ss32 = pow(np.sin(arg32/2.0), 2.0)
    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    prob = [delta(alpha, beta) \
            - 4.0 * ( J(U, alpha, beta, 1, 0).real*ss21
                    + J(U, alpha, beta, 2, 0).real*ss31
                    + J(U, alpha, beta, 2, 1).real*ss32 ) \
            + (2.0 if (not nubar) else -2.0) * \
                    ( J(U, alpha, beta, 1, 0).imag*s21
                    + J(U, alpha, beta, 2, 0).imag*s31
                    + J(U, alpha, beta, 2, 1).imag*s32 ) \
            for alpha in [0,1,2] for beta in [0,1,2]]

    # [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
    return np.array(prob).reshape((3,3))


__all__ = [
    'osc_prob_2nu_vacuum_std',
    'osc_prob_2nu_matter_std',
    'delta',
    'J',
    'osc_prob_3nu_vacuum_std',
]
