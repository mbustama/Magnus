import numpy as np
import sys
from joblib import Parallel, delayed
from typing import Optional, Callable, Union

# TO-DO: remove this once setup.py and pip are working
import os
sys.path.append(os.path.split(os.path.split(os.getcwd())[0])[0])

import magnus.magnus as magnus

def compute_evolution_operator(H_func: Callable, t_slab: Union[list, np.ndarray], 
    n_tpts_per_slab: int, magnus_exp_order: int, **kwargs) -> np.ndarray:
    """Compute the evolution operator for a given time slab."""
    if t_slab[1] > t_slab[0]:
        return magnus.magnus_expansion(
            lambda t: -1j * H_func(t),
            t0=t_slab[0],
            t1=t_slab[1],
            order=magnus_exp_order,
            n_tpts=n_tpts_per_slab,
            **kwargs,
        )
    else:  # t_slab[1] == t_slab[0]
        n = H_func(t_slab[0]).shape[0]
        return np.eye(n, n)


def osc_prob(H_func: Callable, t_ini: float, t_fin: float, n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4, n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[float]=None, atol: Optional[float]=None, 
    growth_factor_n_tpts_per_slab: Optional[float]=1.5, validate_input: Optional[bool]=True,
    verbose: Optional[bool]=False, **kwargs) -> np.ndarray:

    # Validate input; set validate_input to False for speed-up.
    if validate_input:

        try:
            if (t_fin < t_ini): 
                raise ValueError("Error in magnus: oscprob.osc_prob: t_fin must be >= t_ini.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if (magnus_exp_order < 1): 
                raise ValueError("Error in magnus: oscprob.osc_prob: magnus_exp_order must be >= 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (rtol <= 0.0)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: rtol must be None or > 0.0.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((atol is not None) and (atol <= 0.0)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: atol must be None or > 0.0.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (growth_factor_n_tpts_per_slab < 1.0)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: growth_factor_n_tpts_per_slab" + 
                    " must be >= 1.0.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        H_ini = H_func(t_ini)

        try:
            if not isinstance(H_func(t_ini), np.ndarray):
                raise ValueError("Error in magnus: oscprob.osc_prob: H_func must return a numpy array.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if H_func(t_ini).shape[0] != H_func(t_ini).shape[1]:
                raise ValueError("Error in magnus: oscprob.osc_prob: H_func must return a square matrix.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

    # The array (or list) t_slab_edges contains user-provided pairs of start and end times, 
    # [ti, tf]_k, that define the initial and final times of each of the k-th time slab.  It is up
    # to the user to ensure that the chain of time slabs covers the full range [t_ini, t_fin] 
    # without leaving gaps.  I.e., the user should ensure that ti_{k+1} = tf_k.  
    if (t_slab_edges is None):
        # If t_slab_edges == None, then divide the internval [t_ini, t_fin] evenly into a number
        # n_slabs of time slabs.  
        dt = (t_fin-t_ini)/n_slabs # Size of one time slab
        t_slab_edges = [[t_ini+dt*i, t_ini+dt*(i+1)] for i in range(n_slabs)]

    first_iteration = True # Flag to know if this is the first iteration of the calculation
    while True:
        # Within each slab, t_slab, we use n_tpts_per_slab time-evaluations to compute the integrals
        # of the Magnus expansion, from t_slab[0] to t_slab[1].  U_chain contains the chain of time-
        # ordered evolution operators, each computed in one time slab 
        if (n_jobs == 1): # No parallelization
            U_chain = [compute_evolution_operator(H_func, t_slab, n_tpts_per_slab, magnus_exp_order,
                integration_method=integration_method, **kwargs) for t_slab in t_slab_edges]
        else: # Run n_jobs jobs in parallel
            U_chain = Parallel(n_jobs=n_jobs)(  
                delayed(compute_evolution_operator)(
                    H_func, t_slab, n_tpts_per_slab, magnus_exp_order, 
                    integration_method=integration_method, **kwargs
                )
                for t_slab in t_slab_edges
            )

        # Now compute the time-ordered product of all evolution operators across all slabs
        Utot = np.linalg.multi_dot(U_chain) if n_slabs > 1 else U_chain[0]

        # Using Utot, compute all the survival and transition probabilities in a probability matrix
        # P = (np.abs(Utot)**2).T and return that matrix.
        P = (np.abs(Utot)**2).T

        # If no target relative tolerance (rtol) or absolute tolerance (atol) of the probability is
        # requested, then return the result obtained already.  If, instead, a target tolerance is
        # requested, then increase the number of points per slab approximately by the factor
        # growth_factor_n_tpts_per_slab, and repeat the probability calculation until the desired
        # tolerance is reached.
        if ((rtol is None) and (atol is None)): # No target tolerance requested: return right away
            return P
        else: # Target tolerance requested: iterate until tolerance is reached
            if first_iteration == False:
                if np.allclose(P, P_old, rtol=rtol, atol=atol): # Compare old and new P matrices
                    return P
            else:
                P_old = P
                first_iteration = False
            # Increase the number of points per slab approximately by growth_factor_n_tpts_per_slab
            n_tpts_per_slab = int(growth_factor_n_tpts_per_slab*n_tpts_per_slab) 


if __name__ == "__main__":
    def H_2nu_func(t):
        return np.array([[1+1j*t, 2*t], [2*t, 4-1j*t]], dtype=np.complex128)
    def H_3nu_func(t):
        return np.array([[1+1j*t, 2*t, 3j*t], [2*t, 4-1j*t, 5+2j*t], [-3j*t, 5-2j*t, 1]], 
            dtype=np.complex128)

    t_ini, t_fin = 0.0, 1.0
    # prob = osc_prob(H_2nu_func, t_ini, t_fin, n_slabs=100, n_tpts_per_slab=100, magnus_exp_order=6,
    #     integration_method='simpson', n_jobs=1)
    # print(prob)
    # prob = osc_prob(H_3nu_func, t_ini, t_fin, n_slabs=100, n_tpts_per_slab=100, magnus_exp_order=6,
    #     integration_method='simpson', n_jobs=1)
    # print(prob)
    prob = osc_prob(H_3nu_func, t_ini, t_fin, n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=2,
        integration_method='simpson', n_jobs=1, rtol=1.e-8, atol=1.e-8, 
        growth_factor_n_tpts_per_slab=1.5)
    print(prob)

