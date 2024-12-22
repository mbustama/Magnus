import numpy as np
import scipy as sp
# from scipy.integrate import quad
# from scipy.linalg import expm
# from scipy.special import factorial
from numba import njit, prange
# from joblib import Parallel, delayed


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


# Function to compute the Magnus expansion terms
def compute_magnus_terms(A, t0, t1, n_tpts=50, order=2):
    """
    Compute the Magnus expansion terms up to a given order over the range [t0, t1].
    """
    nA = A(t0).shape[0]
    zero_matrix = np.zeros((nA, nA), dtype=complex)

    # Precompute time points and weights
    # n_tpts = 10  # Predefined number of points for evaluation
    if t0 > 0.0:
        times = np.logspace(np.log10(t0), np.log10(t1), n_tpts)
    else:
        times = np.linspace(t0, t1, n_tpts)
    delta = (t1 - t0) / (n_tpts - 1)

    def commutator(X, Y):
        return X @ Y - Y @ X

    # def matrix_integral_simpson(M, ta, tb):
    #     """
    #     Numerically integrate M(t) from ta to tb element-wise using Simpson's rule.
        
    #     Parameters:
    #         M (function): A function that takes a scalar t and returns a square complex-valued matrix.
    #         ta (float): Start of the integration interval.
    #         tb (float): End of the integration interval.

    #     Returns:
    #         np.ndarray: The integrated matrix.
    #     """
    #     # Precompute time points and weights
    #     n_points = 100  # Predefined number of points for evaluation
    #     times = np.linspace(ta, tb, n_points)
    #     delta = (tb - ta) / (n_points - 1)

    #     # Evaluate the matrix function at all time points in one go
    #     matrices = np.stack([M(t) for t in times], axis=0)  # Shape: (n_points, nA, nA)
    #     # # Parallelize the evaluation of the matrix function
    #     # matrices = np.stack(Parallel(n_jobs=-1)(delayed(M)(t) for t in times), axis=0)  # Shape: (n_points, nA, nA)

    #     # Use Simpson's rule with precomputed weights for faster integration
    #     result = delta / 3 * (matrices[0] + 4 * np.sum(matrices[1:-1:2], axis=0) + 2 * np.sum(matrices[2:-2:2], axis=0) + matrices[-1])

    #     return result

    # def matrix_integral(matrices):
    #     """
    #     Numerically integrate M(t) from ta to tb element-wise using Simpson's rule.
        
    #     Parameters:
    #         M (function): A function that takes a scalar t and returns a square complex-valued matrix.
    #         ta (float): Start of the integration interval.
    #         tb (float): End of the integration interval.

    #     Returns:
    #         np.ndarray: The integrated matrix.
    #     """
    #     # Precompute time points and weights
    #     # n_points = 100  # Predefined number of points for evaluation
    #     # times = np.linspace(ta, tb, n_points)
    #     # delta = (tb - ta) / (n_points - 1)

    #     # Evaluate the matrix function at all time points in one go
    #     # matrices = np.stack([M(t) for t in times], axis=0)  # Shape: (n_points, nA, nA)
    #     # # Parallelize the evaluation of the matrix function
    #     # matrices = np.stack(Parallel(n_jobs=-1)(delayed(M)(t) for t in times), axis=0)  # Shape: (n_points, nA, nA)

    #     # Use Simpson's rule with precomputed weights for faster integration
    #     result = delta / 3 * (matrices[0] + 4 * np.sum(matrices[1:-1:2], axis=0) + 2 * np.sum(matrices[2:-2:2], axis=0) + matrices[-1])

    #     return result

    # @njit(parallel=True)
    def matrix_integral(times, matrices):
        # n = matrices.shape[2]
        n_times = len(times)
        if (n_times == 1): # Integral from t0 to t0
            return zero_matrix
        elif (n_times == 2): # Integral from t0 to t1
            dt = times[1]-times[0]
            return 0.5*dt*sum(matrices)
        else:
            return sp.integrate.simpson(matrices, x=times, axis=0) 

    magnus_terms = []

    # Precompute the A(t) terms
    At = np.array([A(t) for t in times])

    # Precompute the Omega_1(t) terms, integrating from t0 to t = t0, ..., t1
    o1t = sp.integrate.cumulative_trapezoid(At, x=times, axis=0, initial=0)
    # o1t = [matrix_integral(times[0:i+1], At[0:i+1]) for i in range(n_tpts)]
    magnus_terms.append(o1t[-1]) # Integral from t0 to t1

    # Precompute the Omega_2(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 2:
        o2t_integrand = -0.5 * np.stack([commutator(o1t[i], At[i]) for i in range(n_tpts)], axis=0)
        o2t = sp.integrate.cumulative_trapezoid(o2t_integrand, x=times, axis=0, initial=0)
        # o2t = [matrix_integral(times[0:i+1], o2t_integrand[0:i+1]) for i in range(n_tpts)]
        magnus_terms.append(o2t[-1])

    # Precompute the Omega_3(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 3:
        # Vectorized computation of t1
        t1 = -0.5 * np.stack([commutator(o2t[i], At[i]) for i in range(n_tpts)], axis=0)
        # Vectorized computation of t2
        t2 = (1.0 / 12.0) * np.stack(
            [commutator(o1t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        o3t_integrand = t1 + t2
        # Efficient integration using cumulative trapezoidal rule
        o3t = sp.integrate.cumulative_trapezoid(o3t_integrand, x=times, axis=0, initial=0)
        # The term in the Magnus expansion is the last one, which goes from t0 to t1
        magnus_terms.append(o3t[-1])
        # ------------------------------------------------------------------------------------------
        # Deprecated (slower):
        # t1 = -0.5 * np.array([commutator(o2t[i], At[i]) for i in range(n_tpts)])
        # t2 = (1.0/12.0) * np.array([commutator(o1t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)])
        # o3t_integrand = t1+t2
        # o3t = [matrix_integral(times[0:i+1], o3t_integrand[0:i+1]) for i in range(n_tpts)]
        # magnus_terms.append(o3t[-1])
        # ------------------------------------------------------------------------------------------
        # Deprecated (slower, not sure numba is handled properly):
        # t1, t2, o3t = [], [], []
        # for i in prange(n_tpts):
        #     t1.append(-0.5*commutator(o2t[i], At[i]))
        #     t2.append((1.0/12.0)*commutator(o1t[i], commutator(o1t[i], At[i])))
        # o3t_integrand = t1+t2
        # for i in prange(n_tpts):
        #     o3t.append(matrix_integral(times[0:i+1], o3t_integrand[0:i+1]))
        # magnus_terms.append(o3t[-1])
        # ------------------------------------------------------------------------------------------

    # Precompute the Omega_4(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 4:
        t1 = -0.5 * np.stack([commutator(o3t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = (1.0 / 12.0) * np.stack(
            [commutator(o1t[i], commutator(o2t[i], At[i])) \
                + commutator(o2t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        o4t_integrand = t1 + t2
        o4t = sp.integrate.cumulative_trapezoid(o4t_integrand, x=times, axis=0, initial=0)
        # o4t = [matrix_integral(times[0:i+1], o4t_integrand[0:i+1]) for i in range(n_tpts)]
        magnus_terms.append(o4t[-1])

    # Precompute the Omega_5(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 5:
        t1 = -0.5 * np.stack([commutator(o4t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = (1.0/12.0) * np.stack(
            [commutator(o1t[i], commutator(o3t[i], At[i])) \
            + commutator(o2t[i], commutator(o2t[i], At[i])) \
            + commutator(o3t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        t3 = (-1.0/720.0) * np.stack(
            [commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], 
            commutator(o1t[i], At[i])))) for i in range(n_tpts)], axis=0
        )
        o5t_integrand = t1 + t2 + t3
        o5t = sp.integrate.cumulative_trapezoid(o5t_integrand, x=times, axis=0, initial=0)
        # o5t = [matrix_integral(times[0:i+1], o5t_integrand[0:i+1]) for i in range(n_tpts)]
        magnus_terms.append(o5t[-1])

    # Precompute the Omega_6(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 6:
        t1 = -0.5 * np.stack([commutator(o5t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = (1.0/12.0) * np.stack(
            [commutator(o1t[i], commutator(o4t[i], At[i]))  \
            + commutator(o2t[i], commutator(o3t[i], At[i])) \
            + commutator(o3t[i], commutator(o2t[i], At[i])) \
            + commutator(o4t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        t3 = (-1.0/720.0) * np.stack(
            [commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], commutator(o2t[i], At[i])))) \
            + commutator(o1t[i], commutator(o1t[i], commutator(o2t[i], commutator(o1t[i], At[i]))))\
            + commutator(o1t[i], commutator(o2t[i], commutator(o1t[i], commutator(o1t[i], At[i]))))\
            + commutator(o2t[i], commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], At[i]))))\
            for i in range(n_tpts)], axis=0)
        o6t_integrand = t1 + t2 + t3
        o6t = sp.integrate.cumulative_trapezoid(o6t_integrand, x=times, axis=0, initial=0)
        # o6t = [matrix_integral(times[0:i+1], o6t_integrand[0:i+1]) for i in range(n_tpts)]
        magnus_terms.append(o6t[-1])

    return magnus_terms

# Function to compute the matrix exponential using Magnus expansion
def magnus_expansion(A, t0, t1, n_tpts=50, order=2):
    """
    Compute the matrix exponential of A(t) from t0 to t1 using the Magnus expansion.
    """
    magnus_terms = compute_magnus_terms(A, t0, t1, n_tpts=n_tpts, order=order)
    # print(magnus_terms)
    Omega = sum(magnus_terms)  # Sum the Magnus terms
    # return sp.linalg.cosm(Omega) + 1j*sp.linalg.sinm(Omega)
    return sp.linalg.expm(Omega)


# Tests
# def A(t):
#     return np.array([[1.0*t, 2.0+3j*t],[2.0-3j*t, 2.0]])
# t0, t1 = 0.0, 1.e4
# exp_Omega = magnus_expansion(A, t0, t1, n_tpts=100, order=2)
# print(exp_Omega)

# def A(times):
#     return [np.array([[1.0*t, 2.0+3j*t],[2.0-3j*t, 2.0]]) for t in times]
# times = np.array([0.0, 1.0])
# print(A(times))

