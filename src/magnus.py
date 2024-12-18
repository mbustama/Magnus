import numpy as np
import scipy as sp
from scipy.integrate import quad
import mpmath as mp

# Function to compute the Magnus expansion terms
def compute_magnus_terms(A, t0, t1, order, epsabs=1.e-8, epsrel=1.e-8, limit=50):
    """
    Compute the Magnus expansion terms up to a given order over the range [t0, t1].

    Parameters:
        A (callable): A function returning the n*n matrix A(t) at time t.
        t0 (float): The initial time.
        t1 (float): The final time.
        order (int): The order of the Magnus expansion.

    Returns:
        list: A list of Magnus expansion terms up to the given order.
    """
    def commutator(X, Y):
        return X @ Y - Y @ X

    def integral_of_A(ta, tb):
        """
        Numerically integrate A(t) from ta to tb element-wise.
        """
        n = A(0).shape[0]
        result = np.zeros((n, n), dtype=complex)
        for i in range(n):
            for j in range(n):
                result[i, j] = quad(lambda s: A(s)[i, j], ta, tb, 
                    epsabs=epsabs, epsrel=epsrel, limit=limit, complex_func=True)[0]
        return result

    magnus_terms = []

    # First-order term is integral of A(t) over [t0, t1]
    if order >= 1:
        omega_1 = integral_of_A(t0, t1)
        magnus_terms.append(omega_1)

    # Compute higher-order terms recursively
    for k in range(2, order + 1):
        term = 0
        for j in range(1, k):
            term += (-1)**(j+1) / j * commutator(magnus_terms[j-1], integral_of_A(t0, t1))
        magnus_terms.append(term / k)

    return magnus_terms

# Function to compute the matrix exponential using Magnus expansion
def magnus_expansion(A, t0, t1, order, epsabs=1.e-8, epsrel=1.e-8, limit=100):
    """
    Compute the matrix exponential of A(t) from t0 to t1 using the Magnus expansion.

    Parameters:
        A (callable): A function returning the n*n matrix A(t) at time t.
        t0 (float): The initial time.
        t1 (float): The final time.
        order (int): The order of the Magnus expansion.

    Returns:
        np.ndarray: The matrix exponential exp(Omega(t0, t1)) where Omega is the Magnus series.
    """
    # Compute Magnus terms
    magnus_terms = compute_magnus_terms(A, t0, t1, order, epsabs=epsabs, epsrel=epsrel, limit=limit)

    # Sum the Magnus terms to form Omega
    Omega = sum(magnus_terms)
    # print(Omega, sp.linalg.expm(Omega))
    # print()

    # mp.dps = 50
    # mp.trap_complex = False
    # return np.array(mp.expm(mp.matrix(Omega), method='pade').tolist(), dtype=complex)

    # Compute the matrix exponential of Omega
    # return sp.linalg.cosm(Omega)+1j*sp.linalg.sinm(Omega)
    return sp.linalg.expm(Omega)

# # Example usage
# if __name__ == "__main__":
#     # Example matrix function A(t)
#     def A(t):
#         return np.array([[0, -t], [t, 0]], dtype=complex)

#     t0 = 0.0
#     t1 = 1.0
#     order = 3

#     result = magnus_expansion(A, t0, t1, order)
#     print("Matrix exponential from t0 to t1:")
#     print(result)
