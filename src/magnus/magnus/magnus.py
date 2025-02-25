import numpy as np
import scipy as sp
import sys
from typing import Optional, Callable

# from scipy.integrate import quad
# from scipy.linalg import expm
# from scipy.special import factorial
# from numba import njit, prange
# from joblib import Parallel, delayed
import numba as nb


# Bernoulli numbers
B = {
    0:  1,
    1:  -0.5, # Use the negative-number convention
    2:  1./6.,
    3:  0,
    4:  -1./30.,
    5:  0,
    6:  1./42.,
    7:  0,
    8:  -1./30.,
    9:  0,
    10: 5./66.,
    11: 0,
    12: -691./2730.,
    13: 0,
    14: 7./6.,
    15: 0,
    16: -3617./510.,
    17: 0,
    18: 43867./798.,
    19: 0,
    20: -174611./330.
}

# Multiplicative factors
f1 = 1.0/12.0
f2 = -1.0/720.0

# Validate input
valid_integration_methods = ['trapezoid', 'simpson']


# Function to compute the matrix exponential using Magnus expansion
# @jit(nopython=True)
def magnus_expansion(A: np.ndarray, t0: float, t1: float, n_tpts: Optional[int]=50, 
    order:Optional[int]=2, integration_method:Optional[str]='trapezoid',
    return_magnus_terms: Optional[bool]=False, 
    validate_input: Optional[bool]=True) -> np.ndarray:
    """
    Compute the matrix exponential of A(t) from t0 to t1 using the Magnus expansion.
    """

    # Validate input
    if validate_input:
        try:
            if not (integration_method in valid_integration_methods): 
                raise ValueError("Error in magnus: magnus.compute_magnus_terms: " + \
                    "integration_method must be one of " + str(valid_integration_methods) + ".")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

    # Precompute time points and weights
    if t0 > 0.0:
        times = np.logspace(np.log10(t0), np.log10(t1), n_tpts)
    else:
        times = np.linspace(t0, t1, n_tpts)
    # delta = (t1 - t0) / (n_tpts - 1)

    def commutator(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        return X @ Y - Y @ X

    # @nb.njit(parallel=True, fastmath=True, nogil=True)
    def integral_cumulative_simpson(matrices: np.ndarray, x: np.ndarray, **kwargs) -> np.ndarray:
        # We need to write our custom routine to compute cumulative matrix integrals because the
        # scipy routine `integrate.cumulative_simpson` does not handle complex numbers; it casts
        # them into real numbers.  This is not a problem of the `integrate.cumulative_trapezoid`, 
        # which we do use.
        return np.array([sp.integrate.simpson(matrices[:i+1], x=x[:i+1], axis=0) \
            for i in range(len(times))])

    if integration_method == 'trapezoid':
        integral_cumulative = sp.integrate.cumulative_trapezoid 
    elif integration_method == 'simpson':
        integral_cumulative = integral_cumulative_simpson

    magnus_terms = []

    # Precompute the A(t) terms
    At = np.array([A(t) for t in times])

    matrix_dim = At[0].shape[0]
    magnus_terms = np.empty((order, matrix_dim, matrix_dim), dtype=complex) 

    # Precompute the Omega_1(t) terms, integrating from t0 to t = t0, ..., t1
    o1t = integral_cumulative(At, x=times, axis=0, initial=0)
    magnus_terms[0] = o1t[-1] # Integral from t0 to t1

    # Precompute the Omega_2(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 2:
        o2t_integrand = -0.5 * np.stack([commutator(o1t[i], At[i]) for i in range(n_tpts)], axis=0)
        o2t = integral_cumulative(o2t_integrand, x=times, axis=0, initial=0)
        magnus_terms[1] = o2t[-1]

    # Precompute the Omega_3(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 3:
        t1 = -0.5 * np.stack([commutator(o2t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        # o3t_integrand = t1 + t2
        # o3t = integral_cumulative(o3t_integrand, x=times, axis=0, initial=0)
        o3t = integral_cumulative(t1+t2, x=times, axis=0, initial=0)
        magnus_terms[2] = o3t[-1]

    # Precompute the Omega_4(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 4:
        t1 = -0.5 * np.stack([commutator(o3t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o2t[i], At[i])) \
                + commutator(o2t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        # o4t_integrand = t1 + t2
        # o4t = integral_cumulative(o4t_integrand, x=times, axis=0, initial=0)
        o4t = integral_cumulative(t1+t2, x=times, axis=0, initial=0)
        magnus_terms[3] = o4t[-1]

    # Precompute the Omega_5(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 5:
        t1 = -0.5 * np.stack([commutator(o4t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o3t[i], At[i])) \
            + commutator(o2t[i], commutator(o2t[i], At[i])) \
            + commutator(o3t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        t3 = f2 * np.stack(
            [commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], 
            commutator(o1t[i], At[i])))) for i in range(n_tpts)], axis=0
        )
        # o5t_integrand = t1 + t2 + t3
        # o5t = integral_cumulative(o5t_integrand, x=times, axis=0, initial=0)
        o5t = integral_cumulative(t1+t2+t3, x=times, axis=0, initial=0)
        magnus_terms[4] = o5t[-1]

    # Precompute the Omega_6(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 6:
        t1 = -0.5 * np.stack([commutator(o5t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o4t[i], At[i]))  \
            + commutator(o2t[i], commutator(o3t[i], At[i])) \
            + commutator(o3t[i], commutator(o2t[i], At[i])) \
            + commutator(o4t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        t3 = f2 * np.stack(
            [commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], commutator(o2t[i], At[i])))) \
            + commutator(o1t[i], commutator(o1t[i], commutator(o2t[i], commutator(o1t[i], At[i]))))\
            + commutator(o1t[i], commutator(o2t[i], commutator(o1t[i], commutator(o1t[i], At[i]))))\
            + commutator(o2t[i], commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], At[i]))))\
            for i in range(n_tpts)], axis=0)
        # o6t_integrand = t1 + t2 + t3
        # o6t = integral_cumulative(o6t_integrand, x=times, axis=0, initial=0)
        o6t = integral_cumulative(t1+t2+t3, x=times, axis=0, initial=0)
        magnus_terms = o6t[-1]

    # Sum of the Magnus terms: Omega = np.sum(magnus_terms, axis=0)
    # Return the exponential matrix (expm) of Omega and the terms of the expansion if requested
    if not return_magnus_terms:
        return np.array(sp.linalg.expm(np.sum(magnus_terms, axis=0)))
    else:
        return np.array(sp.linalg.expm(np.sum(magnus_terms, axis=0))), magnus_terms


if __name__ == "__main__":
    def A(t):
        return np.array([[1.0*t, 2.0+3j*t],[2.0-3j*t, 2.0]])
    t0, t1 = 0.0, 1.0
    exp_Omega_1 = magnus_expansion(A, t0, t1, n_tpts=100, order=6, integration_method='trapezoidd')
    print(exp_Omega_1)
    exp_Omega_2 = magnus_expansion(A, t0, t1, n_tpts=100, order=6, integration_method='simpson')
    print(exp_Omega_2)
    print(exp_Omega_1 == exp_Omega_2)
    print(exp_Omega_1-exp_Omega_2)


