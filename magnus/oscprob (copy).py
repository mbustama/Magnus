import numpy as np
import magnus as magnus
# import magnus_dev as magnus


def osc_prob(H_func, t_ini, t_fin, n_slabs=1, n_tpts_per_slab=100, t_slab_edges=None,
    magnus_exp_order=4, **kwargs):

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

    nA = H_ini.shape[0] # Matrix dimension

    # The array (or list) t_slab_edges contains user-provided pairs of start and end times, 
    # [ti, tf]_k, that define the initial and final times of each of the k-th time slab.  It is up
    # to the user to ensure that the chain of time slabs covers the full range [t_ini, t_fin] 
    # without leaving gaps.  I.e., the user should ensure that ti_{k+1} = tf_k.  
    if (t_slab_edges is None):
        # If t_slab_edges == None, then divide the internval [t_ini, t_fin] evenly into a number
        # n_slabs of time slabs.  
        dt = (t_fin-t_ini)/n_slabs # Size of one time slab
        t_slab_edges = [[t_ini+dt*i, t_ini+dt*(i+1)] for i in range(n_slabs)]

    U_chain = [] # Chain of time-ordered evolution operators, each computed in one time slab 
    for t_slab in t_slab_edges:
        if t_slab[1] > t_slab[0]:
            # Within each slab, t_slab, we use n_tpts_per_slab time-evaluations to compute the 
            # integrals of the Magnus expansion, from t_slab[0] to t_slab[1].
            U = magnus.magnus_expansion(lambda t: -1j*H_func(t), 
                t0=t_slab[0], t1=t_slab[1], order=magnus_exp_order, n_tpts=n_tpts_per_slab, 
                **kwargs)
            U_chain.append(U)
        else: # t_slab[1] == t_slab[0]
            U_chain.append(np.eye((nA, nA)))
    # Now compute the time-ordered product of all evolution operators across all slabs
    if n_slabs > 1:
        Utot = np.linalg.multi_dot(U_chain)
    else:
        Utot = U_chain[0]

    # Using Utot, compute all the survival and transition probabilities; save them in the matrix P
    P = (np.abs(Utot)**2).T

    return P


if __name__ == "__main__":
    def H_func(t):
        return np.array([[1 + 1j * t, 2 * t], [3 * t, 4 - 1j * t]], dtype=np.complex128)

    t_ini, t_fin = 0.0, 1.0
    prob = osc_prob(H_func, t_ini, t_fin, n_slabs=2, n_tpts_per_slab=100, magnus_exp_order=2)

    print(prob)

