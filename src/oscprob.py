import numpy as np
import magnus as magnus


def osc_prob(H_func, t_ini, t_fin, t_npts=100, magnus_exp_order=4, **kwargs):

    # Validate input
    if (t_fin < t_ini):
    	raise ValueError("osc_prob: t_fin must be >= t_ini.")
    if (magnus_exp_order < 1):
    	raise ValueError("osc_prob: magnus_exp_order must be >= 1.")
    H_ini = H_func(t_ini)
    if not isinstance(H_func(t_ini), np.ndarray):
        raise ValueError("osc_prob: H_func must return a numpy array.")
    if H_func(t_ini).shape[0] != H_func(t_ini).shape[1]:
        raise ValueError("osc_prob: H_func must return a square matrix.")

    # Compute the evolution operator from t=t_ini to t=t_fin using Magnus expansion
    U = magnus.magnus_expansion(lambda t: -1j*H_func(t), 
        t0=t_ini, t1=t_fin, order=magnus_exp_order, **kwargs)
    # U = magnus.magnus_expansion(H_func, t_ini=t_ini, t_fin=t_fin, t_npts=t_npts, 
    # 	magnus_exp_order=magnus_exp_order, **kwargs)
    # U = 1.0/np.sqrt(2.0)*np.array([[1.0, 1.0j], [1.0j, 1.0]])

    # Matrix dimension
    n = H_ini.shape[0]

	# Create an identity matrix representing all standard basis vectors
    I = np.eye(n)

    # Using U, compute all the survival and transition probabilities; save them in the matrix P
    if t_fin > t_ini:
        P = (np.abs(U)**2).T
		# # Compute the intermediate matrix U @ I
        # U_I = U @ I  # Each column of U_I is U @ nu_a for a specific i
		# # Compute the outer product of each column of U_I with itself
        # P = np.abs(U_I.T @ U_I)**2  # This computes the probability matrix P
        # # P = np.abs(I.T @ U_I)**2  # This computes the probability matrix P
    else: # t_fin = t_ini
        P = I

    return P


# if __name__ == "__main__":
#     def H_func(t):
#         return np.array([[1 + 1j * t, 2 * t], [3 * t, 4 - 1j * t]], dtype=np.complex128)

#     t_ini, t_fin, t_npts, magnus_exp_order = 0.0, 2.0, 10, 4
#     prob = osc_prob(H_func, t_ini, t_fin, t_npts=t_npts, magnus_exp_order=magnus_exp_order)

#     print(prob)
