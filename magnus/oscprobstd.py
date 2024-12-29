import numpy as np
from typing import Union


def osc_prob_2nu_vacuum_std(sth: float, Dm2: float, energy: float, L: float) -> np.ndarray:
    r"""Returns 2nu oscillation vacuum probabilities, std. computation.

    Returns the probabilities for two-neutrino oscillations in vacuum, computed using the standard
    analytical expression of the probabilities.

    Parameters
    ----------
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    list
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

    Parameters
    ----------
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    VCC : float
        Potential due to charged-current interactions of nu_e with electrons.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    list
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

    Parameters
    ----------
    U : list
        3x3 PMNS complex mixing matrix.
    alpha : int
        Index of the initial flavor (0: e, 1: mu, 2: tau).
    beta : int
        Index of the final flavor (0: e, 1: mu, 2: tau).
    k : int
        First index of the sum over mass eigenstates (k = 0, 1, 2).
    j : int
        First index of the sum over mass eigenstates (k = 0, 1, 2).

    Returns
    -------
    float
        J(U, alpha, beta, j, j)
    """
    return np.conj(U[alpha][k])*U[beta][k]*U[alpha][j]*np.conj(U[beta][j])


def osc_prob_3nu_vacuum_std(U: Union[list, np.ndarray], D21: float, D31: float, energy: float, 
    L: float, nubar: Optional[bool]=False) -> np.ndarray:
    r"""Returns 3nu oscillation vacuum probabilities, std. computation.

    Returns the probabilities for three-neutrino oscillations in vacuum, computed using the standard
    analytical expression of the probabilities.

    Parameters
    ----------
    U : list
        3x3 PMNS complex mixing matrix.
    D21 : float
        Mass-squared difference Delta m^2_21.
    D31 : float
        Mass-squared difference Delta m^2_31.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    list
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
    prob = np.array(prob).reshape((3,3))

    return prob
