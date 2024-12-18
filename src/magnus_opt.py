import numpy as np
from scipy.linalg import expm
from multiprocessing import Pool
from functools import partial

def commutator(X, Y):
    return X @ Y - Y @ X


def magnus_expansion(A_func, t_ini, t_fin, t_npts=100, magnus_exp_order=4, num_procs=None, **kwargs):
    """
    Compute the exponential of a square complex matrix A(t) using the Magnus expansion.

    Parameters:
        A_func (callable): A function that takes a scalar t and returns an (n x n) complex matrix A(t).
        order (int): The order of the Magnus expansion to use (default is 4).
        time (float): The time at which to evaluate the exponential (default is 1.0).
        steps (int): The number of steps to use in numerical integration (default is 100).
        num_procs (int): The number of processors to use for multiprocessing (default is None, uses all available processors).

    Returns:
        np.ndarray: The matrix exponential of A(t) over the interval [0, time].
    """
    dt = (t_fin-t_ini) / (t_npts-1)
    times = np.linspace(t_ini, t_fin, t_npts)

    n = A_func(t_ini).shape[0]

    # Recursive computation of nested commutators
    def compute_nested_commutators(A_func, time_points, depth):
        if depth == 1:
            return np.array([A_func(t) for t in time_points])
        prev_commutators = compute_nested_commutators(A_func, time_points, depth - 1)
        new_commutators = []
        for t1_idx, t1 in enumerate(time_points):
            commutator_sum = np.zeros((n, n), dtype=np.complex128)
            for t2_idx, t2 in enumerate(time_points[:t1_idx + 1]):
                commutator_sum += commutator(prev_commutators[t1_idx], prev_commutators[t2_idx]) * dt
            new_commutators.append(commutator_sum)
        return np.array(new_commutators)

    Omega = np.zeros((n, n), dtype=np.complex128)
    for k in range(1, magnus_exp_order + 1):
        coeff = dt ** k / np.math.factorial(k)
        nested_commutators = compute_nested_commutators(A_func, times[:-1], k)
        integral = np.sum(nested_commutators, axis=0)
        Omega += coeff * integral

    # def compute_nested_commutators(A_func, times, depth, t1_idx, dt):
    #     if depth == 1:
    #         return A_func(times[t1_idx])
    #     prev = compute_nested_commutators(A_func, times, depth - 1, t1_idx, dt)
    #     return sum(commutator(prev, A_func(t)) * dt for t in times[:t1_idx + 1])

    # Omega = np.zeros((n, n), dtype=np.complex128)
    # with Pool(processes=num_procs) as pool:
    #     for k in range(1, magnus_exp_order + 1):
    #         coeff = dt ** k / np.math.factorial(k)
    #         args = [(A_func, times, k, t_idx, dt) for t_idx in range(len(times) - 1)]
    #         nested_commutators = pool.map(partial(compute_nested_commutators, A_func, times), args)
    #         Omega += coeff * sum(nested_commutators)

    return expm(Omega)


# # Example usage
# if __name__ == "__main__":
#     def A_func(t):
#         return np.array([[1 + 1j * t, 2 * t], [3 * t, 4 - 1j * t]], dtype=np.complex128)

#     exp_A = magnus_expansion(A_func, order=4, time=2.0, steps=100, num_procs=4)
#     print("Exponential of A(t) over [0, 2]:")
#     print(exp_A)
