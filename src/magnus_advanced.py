import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import expm
from joblib import Parallel, delayed

# Function to compute the Magnus expansion terms
def compute_magnus_terms(A, t0, t1, magnus_exp_order=2):
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

    def integral_of_A(t0, t1):
        """
        Numerically integrate A(t) from t0 to t1 as a full matrix.
        """
        def A_flat(t):
            return A(t).flatten()

        n = A(0).shape[0]
        result, _ = quad_vec(A_flat, t0, t1)
        return result.reshape((n, n))

    # Cache for intermediate results
    integrals = {}
    integrals[(t0, t1)] = integral_of_A(t0, t1)

    magnus_terms = []

    # First-order term
    if magnus_exp_order >= 1:
        magnus_terms.append(integrals[(t0, t1)])

    # Higher-order terms
    for k in range(2, magnus_exp_order + 1):
        term = 0
        for j in range(1, k):
            term += (-1)**(j+1) / j * commutator(magnus_terms[j-1], integrals[(t0, t1)])
        term = term / k
        magnus_terms.append(term)

    return magnus_terms

# Function to compute the matrix exponential using Magnus expansion
def magnus_expansion(A, t0, t1, magnus_exp_order):
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
    magnus_terms = compute_magnus_terms(A, t0, t1, magnus_exp_order)

    # Sum the Magnus terms to form Omega
    Omega = sum(magnus_terms)

    # Compute the matrix exponential of Omega
    return expm(Omega)

# # Example usage
# if __name__ == "__main__":
#     # Example matrix function A(t)
#     def A(t):
#         return np.array([[0, -t], [t, 0]], dtype=complex)

#     t0 = 0.0
#     t1 = 1.0
#     order = 3

#     result = magnus_exponential(A, t0, t1, order)
#     print("Matrix exponential from t0 to t1:")
#     print(result)
