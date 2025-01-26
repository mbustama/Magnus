import numpy as np
import sys
from joblib import Parallel, delayed
from typing import Optional, Callable, Union
from io import TextIOWrapper

# TO-DO: remove this once setup.py and pip are working
import os
sys.path.append(os.path.split(os.path.split(os.getcwd())[0])[0])
sys.path.append(os.path.split(os.getcwd())[0])

import magnus.magnus as magnus
import version as version


def print_run_parameters(H_func: Callable, t_ini: float, t_fin: float, n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4, n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[float]=None, atol: Optional[float]=None, 
    growth_factor_n_slabs: Optional[float]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[float]=1.5, 
    max_num_loops: Optional[int]=50, max_n_slabs: Optional[float]=2000, 
    max_n_tpts_per_slab: Optional[int]=500,
    validate_input: Optional[bool]=True, save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log', verbose: Optional[int]=0, 
    file_log: Optional[TextIOWrapper]=None):

    for f in [None, file_log] if save_log else [None]:
        print(".----------------------------------------.", file=f)
        print("|   __  __                               |", file=f)
        print("|  |  \/  | __ _  __ _ _ __  _   _ ___   |", file=f)
        print("|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |", file=f)
        print("|  | |  | | (_| | (_| | | | | |_| \__ \  |", file=f)
        print("|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |", file=f)
        print("|                |___/                   |", file=f)
        print("'----------------------------------------'", file=f)
        print("Version: "+ version.__version__+"\n", file=f)
        print("Parameters passed to function magnus.osc_prob in this run:", file=f)
        print("   H_func = " + H_func.__name__, file=f)
        print("   t_ini = " + str(t_ini), file=f)
        print("   t_fin = " + str(t_fin), file=f)
        print("   n_slabs = " + str(n_slabs), file=f)
        print("   n_tpts_per_slab = " + str(n_slabs), file=f)
        if t_slab_edges is None:
            print("   n_tpts_per_slab = None", file=f)
        else:
            print("   n_tpts_per_slab = ", file=f)
            for i, t_slab in enumerate(t_slab_edges):
                print("      i" + ": " + str(t_slab), file=f)
        print("   magnus_exp_order = " + str(magnus_exp_order), file=f)
        print("   n_jobs = " + str(n_jobs), file=f)
        print("   integration_method = " + integration_method, file=f)
        print("   rtol = " + str(rtol), file=f)
        print("   atol = " + str(atol), file=f)
        print("   growth_factor_n_slabs = " + str(growth_factor_n_slabs), file=f)
        print("   growth_factor_n_tpts_per_slab = " + str(growth_factor_n_tpts_per_slab), file=f)
        print("   max_num_loops = " + str(max_num_loops), file=f)
        print("   max_n_slabs = " + str(max_n_slabs), file=f)
        print("   max_n_tpts_per_slab = " + str(max_n_tpts_per_slab), file=f)
        print("   validate_input = " + str(validate_input), file=f)
        print("   save_log = " + str(save_log), file=f)
        print("   filename_log = " + filename_log, file=f)
        print("   verbose = " + str(verbose), file=f)

    return


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
    growth_factor_n_slabs: Optional[float]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[float]=1.5, 
    max_num_loops: Optional[int]=50, max_n_slabs: Optional[float]=2000, 
    max_n_tpts_per_slab: Optional[int]=500, validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, filename_log: Optional[str]='./out.log',
    verbose: Optional[int]=0, **kwargs) -> np.ndarray:

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
            if ((rtol is not None) and (atol is not None) and (growth_factor_n_slabs < 1.0)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: growth_factor_n_slabs" + 
                    " must be >= 1.0.")
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

        try:
            if ((rtol is not None) and (atol is not None) and (max_num_loops <= 1)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: max_num_loops must be > 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (max_n_slabs <= 1)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: max_n_slabs must be > 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (max_n_slabs <= 2)): 
                raise ValueError("Error in magnus: oscprob.osc_prob: max_n_tpts_per_slab" + \
                    " must be > 2.")
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

    # Open a log file if requested
    file_log = open(filename_log, 'w') if save_log else None

    # Print a list of all the parameters passed to the osc_prob function and their values
    if (verbose > 1):
        print_run_parameters(H_func, t_ini, t_fin, n_slabs, n_tpts_per_slab, t_slab_edges,
            magnus_exp_order, n_jobs, integration_method, rtol, atol, growth_factor_n_slabs,
            growth_factor_n_tpts_per_slab, max_num_loops, max_n_slabs, max_n_tpts_per_slab,
            validate_input, save_log, filename_log, verbose, file_log)

    loop_count = 1 # Loop counter
    # Copy this to remember whether the function was originally called with predefine slab edges, 
    # or whether we can increase the number of edges (n_slabs) progressively to reach tolerance
    t_slab_edges_original = t_slab_edges 
    # Flags to signal whether a loop has been run with n_slabs == max_n_slabs or 
    # n_tpts_per_slab = max_n_tpts_per_slab
    ran_with_max_n_slabs, ran_with_max_n_tpts_per_slab = False, False 
    # Flags to signal whether we have already printed the warning that we have reached 
    # n_slabs == max_n_slabs or n_tpts_per_slab = max_n_tpts_per_slab, so as not to print it again
    warned_reached_max_n_slabs, warned_reached_max_n_tpts_per_slab = False, False

    while True:

        # These checks only apply when osc_prob is run with a requested tolerance (rtol, atol) that
        # should be achieved.
        if ((rtol is not None) and (atol is not None)):
            # Reached maximum allowed number of loops: exit loop, return the probability matrix
            if (loop_count > max_num_loops):
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        print("   Warning: Number of loops (loop_count = " + str(loop_count-1) + \
                            ") reached maximum allowed (max_num_loops = " + str(max_num_loops) + \
                            "). Requested tolerance not achieved. Try increasing max_num_loops.\n",
                            file=f)
                if save_log: file_log.close()
                return P
            # Reached maximum allowed number of slabs: continue execution
            if (n_slabs == max_n_slabs):
                if ((verbose > 0) and not warned_reached_max_n_slabs):
                    for f in [None, file_log] if save_log else [None]:
                        print("   Warning: Number of slabs (n_slabs)" + \
                            " reached maximum allowed (max_n_slabs = " + str(max_n_slabs) + ").",
                            file=f)
                        warned_reached_max_n_slabs = True
            # Reached maximum allowed number of time-points per slab: continue execution
            if (n_tpts_per_slab == max_n_tpts_per_slab):
                if ((verbose > 0) and not warned_reached_max_n_tpts_per_slab):
                    for f in [None, file_log] if save_log else [None]:
                        print("   Warning: Number of time-points per slab (n_tpts_per_slab)" + \
                            " reached maximum allowed (max_n_tpts_per_slab = " + \
                            str(max_n_tpts_per_slab) + ").", file=f)
                        warned_reached_max_n_tpts_per_slab = True
            # Reached maximum allowed number of slabs and maximum allowed number of time-points per
            # slab: exit loop, return the probability matrix
            if (ran_with_max_n_slabs and ran_with_max_n_tpts_per_slab):
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        print("   Warning: Number of slabs (n_slabs) and time-points per slab" + \
                            " (n_tpts_per_slab) reached maximum allowed (max_n_slabs = " + \
                            str(max_n_slabs) + ", max_n_tpts_per_slab = " + \
                            str(max_n_tpts_per_slab) + ").", file=f)
                        print("   Warning: Returning probability, but requested tolerance (rtol = " + \
                            str(rtol) + ", atol = " + str(atol) + ") not achieved." + \
                            " Try increasing max_n_slabs or max_n_tpts_per_slab.\n", file=f)
                if save_log: file_log.close()
                return P

        # The array (or list) t_slab_edges contains user-provided pairs of start and end times, 
        # [ti, tf]_k, that define the initial and final times of each of the k-th time slab.  It is 
        # up to the user to ensure that the chain of time slabs covers the full range [t_ini, t_fin] 
        # without leaving gaps.  I.e., the user should ensure that ti_{k+1} = tf_k.  
        if (t_slab_edges_original is None):
            # If t_slab_edges == None, then divide the internval [t_ini, t_fin] evenly into a number
            # n_slabs of time slabs.  
            dt = (t_fin-t_ini)/n_slabs # Size of one time slab
            t_slab_edges = [[t_ini+dt*i, t_ini+dt*(i+1)] for i in range(n_slabs)]

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
        # tolerance is achieved.
        if ((rtol is None) and (atol is None)): # No target tolerance requested: return right away
            if save_log: file_log.close()
            return P
        else: # Target tolerance requested: iterate until tolerance is achieved
            if (verbose > 1):
                for f in [None, file_log] if save_log else [None]:
                    if (loop_count == 1):
                        print("\nRunning loops until requested rtol and atol are achieved:", file=f)
                    print("   Loop #" + str(loop_count) + ":", file=f)
                    print("      n_slabs = " + str(n_slabs), file=f)
                    print("      n_tpts_per_slab = " + str(n_tpts_per_slab), file=f)
            if loop_count > 1:
                # Compare the new and old probability matrices element-wise
                if np.allclose(P, P_old, rtol=rtol, atol=atol):
                    if (verbose > 0):
                        for f in [None, file_log] if save_log else [None]:
                            print("   Requested tolerance achieved\n", file=f)
                    if save_log: file_log.close()
                    return P
                else:
                    P_old = np.ndarray.copy(P)
            else: # loop_count == 1
                P_old = np.ndarray.copy(P)
            # Increase the number of slabs approximately by growth_factor_n_slabs.  Do it only
            # if the slab edges have *not* been explicitly provided by the user in t_slab_edges.
            if t_slab_edges_original is None:
                ran_with_max_n_slabs = False if n_slabs < max_n_slabs else True
                n_slabs = min(int(growth_factor_n_slabs*n_slabs), max_n_slabs)
            # Increase the number of points per slab approximately by growth_factor_n_tpts_per_slab
            ran_with_max_n_tpts_per_slab = False if n_tpts_per_slab < max_n_tpts_per_slab else True
            n_tpts_per_slab = min(int(growth_factor_n_tpts_per_slab*n_tpts_per_slab), 
                max_n_tpts_per_slab)
            loop_count += 1


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
    prob = osc_prob(H_3nu_func, t_ini, t_fin, n_slabs=10, n_tpts_per_slab=20, magnus_exp_order=4,
        integration_method='simpson', n_jobs=10, rtol=1e-5, atol=1.e-5, 
        growth_factor_n_slabs=1.5, growth_factor_n_tpts_per_slab=1.5, 
        max_num_loops=50, max_n_slabs=200, max_n_tpts_per_slab=150, 
        save_log=True, filename_log='./out.log', verbose=2)
    print(prob)
