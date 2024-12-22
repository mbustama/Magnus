import numpy as np
import scipy as sp
# from scipy.integrate import quad
# from scipy.linalg import expm
# from scipy.special import factorial
# from numba import njit, prange
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
# @njit(parallel=True)
def compute_magnus_terms(A, t0, t1, order, epsabs=1e-8, epsrel=1e-8, limit=50):
    """
    Compute the Magnus expansion terms up to a given order over the range [t0, t1].
    """
    nA = A(t0).shape[0]
    zero_matrix = np.zeros((nA, nA), dtype=complex)

    def commutator(X, Y):
        return X @ Y - Y @ X

    # Use vectorized integration for improved performance
    def matrix_integral_quad(M, ta, tb):
        """
        Numerically integrate A(t) from ta to tb element-wise.
        """
        result = zero_matrix #np.zeros((nA, nA), dtype=complex)
        
        # Vectorized computation using parallel loops
        for i in prange(nA):
            for j in prange(nA):
                result[i, j] = sp.integrate.quad(lambda t: M(t)[i, j], ta, tb, 
                    epsabs=epsabs, epsrel=epsrel, limit=limit, complex_func=True)[0]
        return np.array(result)

    # def matrix_integral_simpson(M, ta, tb):
    #     """
    #     Numerically integrate A(t) from ta to tb element-wise.
    #     """
    #     result = zero_matrix #np.zeros((nA, nA), dtype=complex)
        
    #     # Vectorized computation using parallel loops
    #     times = np.linspace(ta, tb, 200)
    #     y = np.array([M(t) for t in times])
    #     for i in prange(nA):
    #         for j in prange(nA):
    #             # Sample function
    #             result[i, j] = sp.integrate.simpson(y[:, i, j], x=times)
    #     return np.array(result)

    # def matrix_integral_simpson(M, ta, tb):
    #     """
    #     Numerically integrate A(t) from ta to tb element-wise.
    #     """
    #     npts = 10
    #     times = np.linspace(ta, tb, npts)
        
    #     # Evaluate the matrix function at all time points (shape: n_points x nA x nA)
    #     matrices = np.array([M(t) for t in times])

    #     # Integrate each matrix element over time
    #     result = sp.integrate.simpson(matrices, x=times, axis=0)

    #     return result

    def matrix_integral_simpson(M, ta, tb):
        """
        Numerically integrate M(t) from ta to tb element-wise using Simpson's rule.
        
        Parameters:
            M (function): A function that takes a scalar t and returns a square complex-valued matrix.
            ta (float): Start of the integration interval.
            tb (float): End of the integration interval.

        Returns:
            np.ndarray: The integrated matrix.
        """
        # Precompute time points and weights
        n_points = 100  # Predefined number of points for evaluation
        times = np.linspace(ta, tb, n_points)
        delta = (tb - ta) / (n_points - 1)

        # Evaluate the matrix function at all time points in one go
        matrices = np.stack([M(t) for t in times], axis=0)  # Shape: (n_points, nA, nA)
        # # Parallelize the evaluation of the matrix function
        # matrices = np.stack(Parallel(n_jobs=-1)(delayed(M)(t) for t in times), axis=0)  # Shape: (n_points, nA, nA)

        # Use Simpson's rule with precomputed weights for faster integration
        result = delta / 3 * (matrices[0] + 4 * np.sum(matrices[1:-1:2], axis=0) + 2 * np.sum(matrices[2:-2:2], axis=0) + matrices[-1])

        return result

    magnus_terms = []

    matrix_integral = matrix_integral_simpson

    def omega_2(ta, tb):
        def integrand_omega_2(t):
            return -0.5 * commutator(matrix_integral(A, ta, t), A(t))
        return matrix_integral(integrand_omega_2, ta, tb)

    def omega_3(ta, tb):
        def integrand_omega_3(t):
            At = A(t)
            o1t = matrix_integral(A, t0, t)
            o2t = omega_2(t0, t)
            t1 = -0.5 * commutator(o2t, At)
            t2 = (1.0 / 12.0) * commutator(o1t, commutator(o1t, At))
            return t1 + t2
        return matrix_integral(integrand_omega_3, ta, tb)

    def omega_4(ta, tb):
        def integrand_omega_4(t):
            At = A(t)
            o1t = matrix_integral(A, t0, t)
            o2t = omega_2(t0, t)
            o3t = omega_3(t0, t)
            t1 = -0.5*commutator(o3t, At)
            t2 = commutator(o1t, commutator(o2t, At)) + commutator(o2t, commutator(o1t, At))
            t2 = (1.0/12.0)*t2
            return t1+t2
        return matrix_integral(integrand_omega_4, ta, tb)

    def omega_5(ta, tb):
        def integrand_omega_5(t):
            At = A(t)
            o1t = matrix_integral(A, t0, t)
            o2t = omega_2(t0, t)
            o3t = omega_3(t0, t)
            o4t = omega_4(t0, t)
            t1 = -0.5*commutator(o4t, At)
            t2 = commutator(o1t, commutator(o3t, At))
            t2 += commutator(o2t, commutator(o2t, At))
            t2 += commutator(o3t, commutator(o1t, At))
            t2 = (1./12.)*t2
            t3 = (-1./720.)*commutator(o1t, commutator(o1t, commutator(o1t, commutator(o1t, At))))
            return t1+t2+t3
        return matrix_integral(integrand_omega_5, ta, tb)

    def omega_6(ta, tb):
        def integrand_omega_6(t):
            At = A(t)
            o1t = matrix_integral(A, t0, t)
            o2t = omega_2(t0, t)
            o3t = omega_3(t0, t)
            o4t = omega_4(t0, t)
            o5t = omega_5(t0, t)
            t1 = -0.5*commutator(o5t, At)
            t2 = commutator(o1t, commutator(o4t, At))
            t2 += commutator(o2t, commutator(o3t, At))
            t2 += commutator(o3t, commutator(o2t, At))
            t2 += commutator(o4t, commutator(o1t, At))
            t2 = (1./12.)*t2
            t3 = commutator(o1t, commutator(o1t, commutator(o1t, commutator(o2t, At))))
            t3 += commutator(o1t, commutator(o1t, commutator(o2t, commutator(o2t, At))))
            t3 += commutator(o1t, commutator(o2t, commutator(o1t, commutator(o1t, At))))
            t3 += commutator(o2t, commutator(o1t, commutator(o1t, commutator(o1t, At))))
            t3 = (-1./720)*t3
            return t1+t2+t3
        return matrix_integral(integrand_omega_6, ta, tb)


    # First-order term: integral of A(t) over [t0, t1]
    if order >= 1:
        magnus_terms.append(matrix_integral(A, t0, t1))

    if order >= 2:
        magnus_terms.append(omega_2(t0, t1))

    if order >= 3:
        magnus_terms.append(omega_3(t0, t1))

    if order >= 4:
        magnus_terms.append(omega_4(t0, t1))

    if order >= 5:
        magnus_terms.append(omega_5(t0, t1))

    if order >= 6:
        magnus_terms.append(omega_6(t0, t1))


    # def ad(k, t):
    #     if (k > 0):
    #         return commutator(matrix_integral(A, t0, t), ad(k-1, t))
    #     else:
    #         return A(t)

    # def S(n, j, t):
    #     # print("b")
    #     if (j == 1):
    #         return commutator(omega_n(n-1, t0, t), A(t))
    #     elif (j == n-1):
    #         return ad(n-1, t)
    #     else: # (2 <= j <= (n-1)):
    #         terms = []
    #         for m in prange(1, n-j+1):
    #             terms.append(commutator(omega_n(m, t0, t), S(n-m, j-1, t)))
    #         return sum(terms)

    # def omega_n(n, ta, tb):
    #     if (n == 1):
    #         return matrix_integral(A, ta, tb)
    #     else:
    #         def integrand_omega_n(t):
    #             terms = []
    #             for j in prange(1, n):
    #                 terms.append(B[j]/sp.special.factorial(j)*S(n, j, t) if B[j] != 0 else zero_matrix)
    #             return sum(terms)
    #         return matrix_integral(integrand_omega_n, ta, tb)

    # if order >= 7:
    #     for n in range(7, order+1):
    #         print("n = ", n)
    #         magnus_terms.append(omega_n(n, t0, t1))
    #         print(magnus_terms)

    return magnus_terms

# Function to compute the matrix exponential using Magnus expansion
def magnus_expansion(A, t0, t1, order, epsabs=1e-8, epsrel=1e-8, limit=100):
    """
    Compute the matrix exponential of A(t) from t0 to t1 using the Magnus expansion.
    """
    magnus_terms = compute_magnus_terms(A, t0, t1, order, epsabs, epsrel, limit)
    # print(magnus_terms)
    Omega = sum(magnus_terms)  # Sum the Magnus terms
    return sp.linalg.expm(Omega)
    # return sp.linalg.exp_multiply(Omega, np.eye(nA))


# Tests
# def A(t):
#     return np.array([[1.0*t, 2.0+3j*t],[2.0-3j*t, 2.0]])
# t0, t1, order = 0.0, 0.5, 2
# exp_Omega = magnus_expansion(A, t0, t1, order, epsabs=1e0, epsrel=1e0, limit=100)
# print(exp_Omega)

# def A(times):
#     return [np.array([[1.0*t, 2.0+3j*t],[2.0-3j*t, 2.0]]) for t in times]
# times = np.array([0.0, 1.0])
# print(A(times))

