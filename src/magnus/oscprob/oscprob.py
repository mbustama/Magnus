"""Module oscprob

Contains routines to compute the neutrino oscillation probability.
"""

__version__ = '0.1'
__author__ = 'Mauricio Bustamante'

import numpy as np
import sys
import platform
from joblib import Parallel, delayed
from typing import Optional, Callable, Union
from io import TextIOWrapper
from inspect import signature

# TO-DO: remove this once setup.py and pip are working
import os
sys.path.append(os.path.split(os.path.split(os.getcwd())[0])[0])
sys.path.append(os.path.split(os.getcwd())[0])

import magnus.magnus as magnus
import magnus.globaldefs as gd
import magnus.hamiltonians.hamiltonians2nu as hamiltonians2nu
import magnus.hamiltonians.hamiltonians3nu as hamiltonians3nu
import magnus.matter as matter
import version as version
import authors as authors


def print_run_parameters(
    H_func: Union[Callable, np.ndarray], 
    t_ini: float, 
    t_fin: float, 
    n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[float]=None, 
    atol: Optional[float]=None, 
    growth_factor_n_slabs: Optional[float]=1.5, 
    growth_factor_n_tpts_per_slab: 
    Optional[float]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[float]=1, 
    max_n_slabs: Optional[float]=2000, 
    min_n_tpts_per_slab: Optional[int]=2, 
    max_n_tpts_per_slab: Optional[int]=500,
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX, 
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log', 
    verbose: Optional[int]=0, 
    file_log: Optional[TextIOWrapper]=None
):

    def print_banner(file: TextIOWrapper=None):
        if file is None:
            print(gd.cstyle.CBLUEBG + ".----------------------------------------." + gd.cstyle.CEND,
                file=f)
            print(gd.cstyle.CBLUEBG + "|   __  __                               |" + gd.cstyle.CEND,
                file=f)
            print(gd.cstyle.CBLUEBG + "|  |  \/  | __ _  __ _ _ __  _   _ ___   |" + gd.cstyle.CEND, 
                file=f)
            print(gd.cstyle.CBLUEBG + "|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |" + gd.cstyle.CEND,
                file=f)
            print(gd.cstyle.CBLUEBG + "|  | |  | | (_| | (_| | | | | |_| \__ \  |" + gd.cstyle.CEND,
                file=f)
            print(gd.cstyle.CBLUEBG + "|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |" + gd.cstyle.CEND,
                file=f)
            print(gd.cstyle.CBLUEBG + "|                |___/                   |" + gd.cstyle.CEND,
                file=f)
            print(gd.cstyle.CBLUEBG + "'----------------------------------------'" + gd.cstyle.CEND,
                file=f)
        else: 
            print(".----------------------------------------.", file=f)
            print("|   __  __                               |", file=f)
            print("|  |  \/  | __ _  __ _ _ __  _   _ ___   |", file=f)
            print("|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |", file=f)
            print("|  | |  | | (_| | (_| | | | | |_| \__ \  |", file=f)
            print("|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |", file=f)
            print("|                |___/                   |", file=f)
            print("'----------------------------------------'", file=f)

    for f in [None, file_log] if save_log else [None]:
        print_banner(f)
        print("Version: "+ version.__version__ + " | Author(s): " + authors.__authors__ + "\n", 
            file=f)
        print("Parameters passed to function magnus.osc_prob in this run:", file=f)
        if callable(H_func):
            print("   H_func = " + H_func.__name__, file=f)
        else:
            print("   H_func = constant (time-independent)", file=f)
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
        print("   min_n_slabs = " + str(min_n_slabs), file=f)
        print("   max_n_slabs = " + str(max_n_slabs), file=f)
        print("   min_n_tpts_per_slab = " + str(min_n_tpts_per_slab), file=f)
        print("   max_n_tpts_per_slab = " + str(max_n_tpts_per_slab), file=f)
        print("   iterate_over_magnus_exp_order = " + str(iterate_over_magnus_exp_order), file=f)
        print("   min_magnus_exp_order = " + str(min_magnus_exp_order), file=f)
        print("   max_magnus_exp_order = " + str(max_magnus_exp_order), file=f)
        print("   validate_input = " + str(validate_input), file=f)
        print("   save_log = " + str(save_log), file=f)
        print("   filename_log = " + filename_log, file=f)
        print("   verbose = " + str(verbose), file=f)

    return


def validate_input_battery(
    source_func_name: str, 
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None
) -> int:

    try:
        if ( (not isinstance(energy, int)) and (not isinstance(energy, float)) and \
            (not isinstance(energy, list)) and (not isinstance(energy, np.ndarray)) ):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                ": energy must be an int, a float, a 1D list, or a 1D NumPy array.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    try:
        if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) and \
            (np.array(energy).ndim != 1) ):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                ": if energy is a list or NumPy array, it must be 1D.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    try:
        if ( (not isinstance(L, int)) and (not isinstance(L, float)) and \
            (not isinstance(L, list)) and (not isinstance(L, np.ndarray)) ):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                ": L must be an int, a float, a 1D list, or a 1D NumPy array.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    try:
        if ( (isinstance(L, list) or isinstance(L, np.ndarray)) and \
            (np.array(L).ndim != 1) ):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                ": if L is a list or NumPy array, it must be 1D.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    try:
        if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) and \
            (isinstance(L, list) or isinstance(L, np.ndarray)) and \
            (len(energy) != len(L)) ):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                ": since the input energy and L are both lists or NumPy arrays, they must have " + \
                "the same length.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    try:
        if (((nu_i is not None) and (nu_f is None)) or ((nu_i is None) and (nu_f is not None))):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                ": if either nu_i or nu_f is not None, the other flavor must also be not None.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    try:
        if ((nu_i is not None) and (nu_f is not None)):
            flavors = set([gd.NUE, gd.NUMU, gd.NUTAU])
            if ( (not (nu_i in flavors)) or (not (nu_f in flavors))):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                    str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), or gd.NUTAU (" + \
                    str(gd.NUTAU) + ") only.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        return 1
        # sys.exit(1)

    return 0


def compute_evolution_operator(
    H_func: Callable, 
    t_slab: Union[list, np.ndarray], 
    n_tpts_per_slab: int, 
    magnus_exp_order: int, 
    **kwargs
) -> np.ndarray:
    r"""Computes the evolution operator inside a given time slab.  This functions is not designed to
    be called directly by the user, but rather internally by :func:`osc_prob`.

    :param H_func: Hamiltonian, which is a function of time or position that returns a square matrix
        in the form of NumPy array
    :param t_slab: List or Numpy Array specifying the start and end times or positions of the slab,
        i.e., [t0, t1]
    :param n_tpts_per_slab: Number of time-points inside the slab at which to evaluate H_func in 
        order to numerically compute the integrals over time required by the Magnus expansion
    :param magnus_exp_order: Maximum order of Matnus expansion used to compute the evolution
        operator (should not exceed :func:`magnus.globaldefs.MAGNUS_EXP_ORDER_MAX`)
    :param \**kwargs: Additional unspecified arguments

    :return: An NumPy array containing the evolution operator for the given time-slab.

    """
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


def osc_prob(
    H_func: Union[Callable, np.ndarray], 
    t_ini: float, 
    t_fin: float, 
    n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[float]=1.e-3, 
    atol: Optional[float]=1.e-3, 
    growth_factor_n_slabs: Optional[float]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[float]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[float]=1, 
    max_n_slabs: Optional[float]=2000, 
    min_n_tpts_per_slab: Optional[int]=2, 
    max_n_tpts_per_slab: Optional[int]=500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX, 
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0, 
    **kwargs
) -> np.ndarray:
    r"""Computes and returns the neutrino oscillation probability.

    Computes the oscillation probability of neutrinos starting at time 
    (or position) ``t_ini`` and ending at time (or position) ``t_fin``.

    Parameters
    ----------
    H_func
        The Hamiltonian, which is a function of time or position that 
        returns a square matrix (a NumPy array). The Hamiltonian can 
        have complex-valued entries.
    t_ini
        Initial time or position of the neutrino.
    t_fin
        Final time or position of the neutrino.
    n_slabs
        Number of slabs, or subintervals, into which the interval 
        [t_ini, t_fin] is partitioned in order to compute the neutrino 
        evolution operators. A higher value of ``n_slabs`` yields a 
        more accurate probability.

        If no target tolerance is requested (i.e., if ``rtol`` and 
        ``atol`` are both ``None``), then the given value of `n_slabs` 
        is the final number of slabs used in the computation.

        If a target tolerance is requested (i.e., if either ``rtol`` or 
        ``atol`` is not ``None``), then the given value of `n_slabs` is 
        ignored. Instead, the number of slabs is increased 
        progressively, starting from ``min_n_slabs``, until the 
        tolerance is achieved or until we hit ``max_n_slabs``, whichever
        happens first.
    n_tpts_per_slab
        Number of time-points inside the slab at which to evaluate 
        H_func in order to numerically compute the integrals over time 
        required by the Magnus expansion. A higher value of 
        ``n_tpts_per_slab`` yields a more accurate probability.
    t_slab_edges
        XXX
    magnus_exp_order
        XXX
    n_jobs
        XXX
    integration_method
        XXX
    rtol
        XXX
    atol
        XXX
    growth_factor_n_slabs
        XXX
    growth_factor_n_tpts_per_slab
        XXX
    max_num_loops
        XXX
    min_n_slabs
        XXX
    max_n_slabs
        XXX
    min_n_tpts_per_slab
        XXX
    max_n_tpts_per_slab
        XXX
    iterate_over_magnus_exp_order
        XXX
    min_magnus_exp_order
        XXX
    max_magnus_exp_order
        XXX
    validate_input
        XXX
    save_log
        XXX
    filename_log
        XXX
    file_log
        XXX
    close_file_log_upon_exit
        XXX
    verbose
        XXX
    \**kwargs
        Additional unspecified arguments

    Returns
    -------
    np.ndarray
        NumPy array containing the probability matrix.
    """

    # Validate input; set validate_input to False for speed-up.
    if validate_input:

        try:
            if (t_fin < t_ini): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: t_fin must be >=" + \
                    " t_ini.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if (magnus_exp_order < 1): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: magnus_exp_order " + \
                    "must be >= 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (rtol <= 0.0)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: rtol must be None " + \
                    "or > 0.0.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((atol is not None) and (atol <= 0.0)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: atol must be None " + \
                    "or > 0.0.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (growth_factor_n_slabs < 1.0)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: " + \
                    "growth_factor_n_slabs must be >= 1.0.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (growth_factor_n_tpts_per_slab < 1.0)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: " + \
                    "growth_factor_n_tpts_per_slab must be >= 1.0.") 
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ( ((rtol is not None) and (atol is not None)) and \
                ((growth_factor_n_slabs == 1.0) and (growth_factor_n_tpts_per_slab == 1.0)) ): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: since a target " + \
                    "tolerance has been requested, either growth_factor_n_slabs, " + \
                    "growth_factor_n_tpts_per_slab, or both must be > 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (max_num_loops <= 1)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: max_num_loops must" + \
                    " be > 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (max_n_slabs <= 1)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: max_n_slabs must " + \
                    "be > 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((rtol is not None) and (atol is not None) and (max_n_slabs <= 2)): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: max_n_tpts_per_slab" +\
                    " must be > 2.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if ((callable(H_func)) and (len(signature(H_func).parameters) > 1)):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: the provided H_func" +\
                    " is a function of more than one parameter")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        H_test = H_func(t_ini) if callable(H_func) else H_func

        try:
            if not isinstance(H_test, np.ndarray):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: H_func must be a " + \
                    "NumPy (if the Hamiltonian is time-independent) or must return a NumPy array.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if H_test.shape[0] != H_test.shape[1]:
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob: H_func must be a " + \
                    "square matrix (if the Hamiltonian is time-independent) or must return a " + \
                    "square matrix.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

    # If there is no file object given (i.e., if file_log is None), open a log file if requested
    if file_log is None:
        file_log = open(filename_log, 'w') if save_log else None

    # Print a list of all the parameters passed to the osc_prob function and their values
    if (verbose > 1):
        print_run_parameters(H_func, t_ini, t_fin, n_slabs, n_tpts_per_slab, t_slab_edges,
            magnus_exp_order, n_jobs, integration_method, rtol, atol, growth_factor_n_slabs,
            growth_factor_n_tpts_per_slab, max_num_loops, min_n_slabs, max_n_slabs, 
            min_n_tpts_per_slab, max_n_tpts_per_slab, iterate_over_magnus_exp_order,
            min_magnus_exp_order, max_magnus_exp_order, validate_input, save_log, filename_log, 
            verbose, file_log)

    # By default, osc_prob is run using a fixed order of the Magnus expansion (magnus_exp_order), 
    # and the tolerance is achieved (see below) only by changing the number of slabs (n_slabs), of
    # time-points per slab (n_tpts_per_slab), or both, but not by changing the expansion order, 
    # since doing that can be computationally taxing. However, if iterate_over_magnus_exp_order is
    # True, then magnus_exp_order will be progressively increased, from min_magnus_exp_order to
    # max_magnus_exp_order, until the requested tolerance (rtol, atol) is achieved.  This is done
    # by calling the function osc_prob_iterate_over_magnus_exp_order, which in turn calls osc_prob
    # with varying values of magnus_exp_order.
    if ((rtol is not None) and (atol is not None) and iterate_over_magnus_exp_order):
        if (max_magnus_exp_order == min_magnus_exp_order):
            magnus_exp_order = min_magnus_exp_order
            if verbose > 0:
                for f in [None, file_log] if save_log else [None]:
                    warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                    print("\n" + warn_msg + "The flag iterate_over_magnus_exp_order has been " + \
                        "set to True, but with min_magnus_exp_order = max_magnus_exp_order. " + \
                        "Bypassing iteration over magnus_exp_order and calling osc_prob with " + \
                        "fixed magnus_exp_order = " + str(magnus_exp_order) + ".", file=f)
        else: # max_magnus_exp_order == min_magnus_exp_order (further input validation in function)
            if verbose > 0:
                for f in [None, file_log] if save_log else [None]:
                    warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                    print("\n" + warn_msg + " The flag iterate_over_magnus_exp_order has been " + \
                        "set to True, so the calculation of the probability will increase the" + \
                        " value of magnus_exp_order progressively from min_magnus_exp_order = " + \
                        str(magnus_exp_order) + " to max_magnus_exp_order = " + \
                        str(max_magnus_exp_order) + " until the requested tolerance is achieved.", 
                        file=f)
            P = osc_prob_iterate_over_magnus_exp_order(H_func, t_ini, t_fin, n_slabs, 
                n_tpts_per_slab, t_slab_edges, magnus_exp_order, n_jobs, integration_method,
                rtol, atol, growth_factor_n_slabs, growth_factor_n_tpts_per_slab,
                max_num_loops, min_n_slabs, max_n_slabs, min_n_tpts_per_slab, 
                max_n_tpts_per_slab, min_magnus_exp_order, 
                max_magnus_exp_order, validate_input, save_log, filename_log, verbose, 
                file_log, **kwargs)
            if save_log and close_file_log_upon_exit: file_log.close()
            return P

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

    # If a tolerance is requested, start the iterations with a number of slabs equal to the given
    # value of min_n_slabs.
    if ((rtol is not None) and (atol is not None)):
        n_slabs = min_n_slabs
        n_tpts_per_slab = min_n_tpts_per_slab

    # The provided Hamiltonian, H_func, can be either a single-parameter function (of the neutrino
    # position) or, if time-independent, a constant NumPy array (e.g., for oscillations in vacuum
    # or in matter with constant density).  In the latter case, we use this constant Hamiltonian to
    # build a dummy one-parameter function of position that we will need later to call the function
    # compute_evolution_operator.  In this case, first-order Magnus expansion is enough, and so we
    # can overwrite the parameters provided to n_slabs = 1, n_tpts_per_slab = 2, rtol = None, 
    # atol = None for speed-up.
    if not callable(H_func): 
        H = np.copy(H_func)
        def H_func(l: float) -> np.ndarray:
            return H
        magnus_exp_order = 1
        n_slabs = 1
        n_tpts_per_slab = 2
        rtol = None
        atol = None
        n_jobs = 1 # No need to parallelize for this simple computation in a single slab
        if verbose > 0:
            for f in [None, file_log] if save_log else [None]:
                warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                print("\n" + warn_msg + " The provided Hamiltonian is time-independent. " + \
                    "Overwriting the run parameters to magnus_exp_order = 1, n_slabs = 1, " + \
                    "n_tpts_per_slab = 2, rtol = None, atol = None, and n_jobs = 1 for speed-up.",
                    file=f)

    while True:

        # These checks only apply when osc_prob is run with a requested tolerance (rtol, atol) that
        # should be achieved.
        if ((rtol is not None) and (atol is not None)):
            # Reached maximum allowed number of loops: exit loop, return the probability matrix
            if (loop_count > max_num_loops):
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg + " Number of loops (loop_count = " + \
                            str(loop_count-1) + ") reached maximum allowed (max_num_loops = " + \
                            str(max_num_loops) + "). Requested tolerance not achieved. Try " + \
                            "increasing max_num_loops.\n",
                            file=f)
                if save_log and close_file_log_upon_exit: file_log.close()
                return P
            # Reached maximum allowed number of slabs: continue execution
            if (n_slabs == max_n_slabs):
                if ((verbose > 0) and not warned_reached_max_n_slabs):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg +  " Number of slabs (n_slabs) reached maximum " + \
                            "allowed (max_n_slabs = " + str(max_n_slabs) + ").", file=f)
                        warned_reached_max_n_slabs = True
            # Reached maximum allowed number of time-points per slab: continue execution
            if (n_tpts_per_slab == max_n_tpts_per_slab):
                if ((verbose > 0) and not warned_reached_max_n_tpts_per_slab):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg + " Number of time-points per slab " + \
                            "(n_tpts_per_slab) reached maximum allowed (max_n_tpts_per_slab = " + \
                            str(max_n_tpts_per_slab) + ").", file=f)
                        warned_reached_max_n_tpts_per_slab = True
            # Reached maximum allowed number of slabs and maximum allowed number of time-points per
            # slab: exit loop, return the probability matrix
            if (ran_with_max_n_slabs and ran_with_max_n_tpts_per_slab):
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                        print("   " + warn_msg + " Number of slabs (n_slabs) and time-points " + \
                            "per slab (n_tpts_per_slab) reached maximum allowed (max_n_slabs = " + \
                            str(max_n_slabs) + ", max_n_tpts_per_slab = " + \
                            str(max_n_tpts_per_slab) + ").", file=f)
                        print("   " + warn_msg + " Returning probability, but requested " + \
                            "tolerance (rtol = " + str(rtol) + ", atol = " + str(atol) + \
                            ") not achieved. Try increasing max_n_slabs or max_n_tpts_per_slab.\n",
                            file=f)
                if save_log and close_file_log_upon_exit: file_log.close()
                return P

        # The array (or list) t_slab_edges contains user-provided pairs of start and end times, 
        # [ti, tf]_k, that define the initial and final times of each of the k-th time slab.  It is 
        # up to the user to ensure that the chain of time slabs covers the full range [t_ini, t_fin] 
        # without leaving gaps.  I.e., the user should ensure that ti_{k+1} = tf_k.  
        if (t_slab_edges_original is None):
            # If t_slab_edges == None, then divide the interval [t_ini, t_fin] evenly into a number
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
            if save_log and close_file_log_upon_exit: file_log.close()
            return P
        else: # Target tolerance requested: iterate until tolerance is achieved
            if (verbose > 1):
                for f in [None, file_log] if save_log else [None]:
                    if (loop_count == 1):
                        print("\nRunning loops until requested rtol and atol are achieved:", file=f)
                    print("   Loop #" + str(loop_count) + ":", file=f)
                    print("      magnus_exp_order = " + str(magnus_exp_order), file=f)                    
                    print("      n_slabs = " + str(n_slabs), file=f)
                    print("      n_tpts_per_slab = " + str(n_tpts_per_slab), file=f)
            if loop_count > 1:
                # Compare the new and old probability matrices element-wise
                if np.allclose(P, P_old, rtol=rtol, atol=atol):
                    if (verbose > 0):
                        for f in [None, file_log] if save_log else [None]:
                            tol_msg = gd.TOL_MSG_IN_COLOR if f is None else gd.TOL_MSG_NO_COLOR
                            print("   " + tol_msg + " (for fixed magnus_exp_order "+ \
                                "= " + str(magnus_exp_order) + ").\n", file=f)
                    if save_log and close_file_log_upon_exit: file_log.close()
                    return P
                else:
                    P_old = np.ndarray.copy(P)
            else: # loop_count == 1
                P_old = np.ndarray.copy(P)
            # Increase the number of slabs approximately by growth_factor_n_slabs.  Do it only
            # if the slab edges have *not* been explicitly provided by the user in t_slab_edges.
            if t_slab_edges_original is None:
                ran_with_max_n_slabs = False if n_slabs < max_n_slabs else True
                n_slabs_old = n_slabs
                n_slabs = min(round(growth_factor_n_slabs*n_slabs), max_n_slabs)
                # Occasionally, the new number of slabs could be equal to the old number (i.e., if
                # growth_factor_n_slabs is too small or if n_slabs = 1).  If this happens, increase
                # the new number of slabs by 1.
                if ((growth_factor_n_slabs > 1.0) and (n_slabs < max_n_slabs) and \
                    (n_slabs == n_slabs_old)): n_slabs += 1
            # Increase the number of points per slab approximately by growth_factor_n_tpts_per_slab
            ran_with_max_n_tpts_per_slab = False if n_tpts_per_slab < max_n_tpts_per_slab else True
            n_tpts_per_slab_old = n_tpts_per_slab
            n_tpts_per_slab = min(int(growth_factor_n_tpts_per_slab*n_tpts_per_slab), 
                max_n_tpts_per_slab)
            if ((growth_factor_n_tpts_per_slab > 1.0) and \
                (n_tpts_per_slab < max_n_tpts_per_slab) and \
                (n_tpts_per_slab == n_tpts_per_slab_old)): n_tpts_per_slab += 1
            loop_count += 1


def osc_prob_iterate_over_magnus_exp_order(
    H_func: Union[Callable, np.ndarray],
    t_ini: float, 
    t_fin: float,
    n_slabs: Optional[int]=1,
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[float]=1.e-3,
    atol: Optional[float]=1.e-3, 
    growth_factor_n_slabs: Optional[float]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[float]=1.5, 
    max_num_loops: Optional[int]=50,
    min_n_slabs: Optional[float]=1, 
    max_n_slabs: Optional[float]=2000, 
    min_n_tpts_per_slab: Optional[int]=2,
    max_n_tpts_per_slab: Optional[int]=500, 
    min_magnus_exp_order: Optional[int]=1, 
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX, 
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False,
    filename_log: Optional[str]='./out.log',
    verbose: Optional[int]=0,
    file_log: Optional[TextIOWrapper]=None,
    **kwargs
) -> np.ndarray:
    r"""Computes the neutrino oscillation probability until a requested
    tolerance is achieved, including progressively increasing the order
    of the Magnus expansion.
    """

    # Validate input; set validate_input to False for speed-up.
    if validate_input:

        try:
            if (max_magnus_exp_order > gd.MAGNUS_EXP_ORDER_MAX): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                    " oscprob.osc_prob_iterate_over_magnus_exp_order: max_magnus_exp_order must" + \
                    " be <= globaldefs.MAGNUS_EXP_ORDER_MAX = " + str(gd.MAGNUS_EXP_ORDER_MAX) + \
                    ".")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if (min_magnus_exp_order < 1): 
                raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                    " oscprob.osc_prob_iterate_over_magnus_exp_order: max_magnus_exp_order must" + \
                    " be >= 1.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

    # Do this to prevent printing the Magnus header multiple times
    verbose = 1 if verbose > 0 else verbose

    # Call osc_prob to compute the probabilities using increasing order of the Magnus expansion
    # (magnus_exp_order) until, ideally, the requested tolerance (rtol, atol) is reached.  See
    # additional comments in the osc_prob function.  (We call osc_prob below with 
    # iterate_over_magnus_exp_order=False to prevent infinite recursion.)
    iterate_over_magnus_exp_order = False
    # close_file_log_upon_exit = False
    for magnus_exp_order in range(min_magnus_exp_order, max_magnus_exp_order+1):
        if (verbose > 0):
            for f in [None, file_log] if save_log else [None]:
                print("\nComputing probabilities using magnus_exp_order = " + \
                    str(magnus_exp_order) + ": \n", file=f)
        P = osc_prob(H_func, t_ini, t_fin, n_slabs, n_tpts_per_slab, t_slab_edges, magnus_exp_order,
            n_jobs, integration_method, rtol, atol, growth_factor_n_slabs, 
            growth_factor_n_tpts_per_slab, max_num_loops, min_n_slabs, max_n_slabs, 
            min_n_tpts_per_slab, max_n_tpts_per_slab, iterate_over_magnus_exp_order,
            min_magnus_exp_order, max_magnus_exp_order, validate_input, save_log, filename_log,
            file_log=file_log, close_file_log_upon_exit=False, verbose=verbose)
        if (magnus_exp_order == min_magnus_exp_order):
            P_old = np.ndarray.copy(P)
        else: # magnus_exp_order > min_magnus_exp_order
            if np.allclose(P, P_old, rtol=rtol, atol=atol):
                if (verbose > 0):
                    for f in [None, file_log] if save_log else [None]:
                        tol_msg = gd.TOL_MSG_IN_COLOR if f is None else gd.TOL_MSG_NO_COLOR
                        print(tol_msg + " using magnus_exp_order = " + str(magnus_exp_order) + \
                            "\n", file=f)
                if save_log: file_log.close()
                return P
    
    # If the for loop finishes, then it means that the requested tolerance could achieved using the
    # maximum test magnus_exp_order allowed for the run.  Return the probability matrix, but show
    # a warning (if verbose).
    if (verbose > 0):
        for f in [None, file_log] if save_log else [None]:
            warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
            print(warn_msg + " returning probability, but requested tolerance not achieved using" +\
                " even the maximum allowed order of the Magnus expansion for this run " + \
                "(max_magnus_exp_order = " + str(max_magnus_exp_order) + ").  Try increasing " + \
                "max_n_slabs, max_n_tpts_per_slab, or max_num_loops.\n", file=f)
    return P

    
def osc_prob_2nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    sth: float, 
    Dm2: float, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    vacuum.

    Parameters
    ----------
    energy
        XXX
    L
        XXX
    sth
        XXX
    Dm2
        XXX
    nu_i
        XXX
    nu_f
        XXX
    validate_input
        XXX

    Returns
    -------
    Union[float, np.narray]
        XXX

    Examples
    --------
    >>> import magnus.oscprob as oscprob
    >>> import magnus.globaldefs as gd
    >>> sth = gd.S12_NO_BF_NUFIT_6_0 # sin(theta) [adim]
    >>> Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]

    Single energy and baseline:

    >>> baseline = 10.*gd.UNIT_KM # 10 km natural units [eV^{-1}]
    >>> energy = 1.*gd.UNIT_MEV # [eV]
    >>> oscprob.osc_prob_2nu_vacuum(energy, baseline, sth, Dm2)
    array([[0.43678029, 0.56321971],
       [0.56321971, 0.43678029]])
    """

    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L

    if validate_input:
        # The function name is sys._getframe().f_code.co_name
        if validate_input_battery(sys._getframe().f_code.co_name,energy, L, nu_i, nu_f) == 1:
            sys.exit(1)

    # Compute the energy-independent part of the Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor below when calling osc_prob.
    h_vac_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)

    # Flag return_float remembers if energy and L were both floats.  If True, return a float, too.
    return_float = isinstance(energy, float) and isinstance(L, float)

    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)  
    L = np.array([L]) if isinstance(L, float) else np.array(L) 

    # Either energy and L are both lists (or NumPy arrays) of the same length; or one is a float and
    # the other is a list (or NumPy array).  Any other possibility will generate an exception.  This
    # exception is raised earlier if validate_input == True, but we check below in case it has been
    # set to False.
    try:
        if not ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or \
            (len(energy) > 1 and len(L) == 1)):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_2nu_vacuum: energy and " + \
                "L must be both int or float; or, if lists (or NumPy arrays), they must have " + \
                "the same length; or, if one is a float or single-entry list, the other must " + \
                "be a list with multiple entries.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If energy is a single value, then transform it into an array containing the value energy 
    # repeated a number of times equal to the length of the L, and vice versa, in order to zip them.
    energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
    L = np.full(len(energy), L[0]) if (len(L) == 1) else L

    # The call to __getitem__ below is a way to return a float if both energy and L were floats
    if ((nu_i is not None) and (nu_f is not None)):
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep, 0.0, xy[1])[nu_i][nu_f]
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
    else:
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep, 0.0, xy[1])
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))


def osc_prob_2nu_matter_constant_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray], 
    sth: float, 
    Dm2: float, 
    rho: float, 
    ratio_number_neutrons_to_protons: Optional[float]=1.0, 
    electron_fraction: Optional[float]=0.5, 
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True
) -> Union[float, np.ndarray]:

    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L

    if validate_input:
        # The function name is sys._getframe().f_code.co_name
        if validate_input_battery(sys._getframe().f_code.co_name, energy, L, nu_i, nu_f) == 1:
            sys.exit(1)

    # Flag return_float remembers if energy and L were both floats.  If True, return a float, too.
    return_float = isinstance(energy, float) and isinstance(L, float)

    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)  
    L = np.array([L]) if isinstance(L, float) else np.array(L) 

    # Either energy and L are both lists (or NumPy arrays) of the same length; or one is a float and
    # the other is a list (or NumPy array).  Any other possibility will generate an exception.  This
    # exception is raised earlier if validate_input == True, but we check below in case it has been
    # set to False.
    try:
        if not ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or \
            (len(energy) > 1 and len(L) == 1)):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_2nu_matter_constant_density: energy and L must be both " + \
                "int or float; or, if lists (or NumPy arrays), they must have the same length;" + \
                " or, if one is a float or single-entry list, the other must be a list with " + \
                "multiple entries.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    s = 1.0 if not nubar else -1.0

    # Electron number density [eV^3]
    num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, 
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction) 

    # Coherent forward potential, VCC [eV]
    VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) 

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor below when calling osc_prob.
    h_vac_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 

    # Compute the matter Hamiltonian only once, to save time.
    h_matt = s*hamiltonians2nu.hamiltonian_2nu_matter(VCC)

    # If energy is a single value, then transform it into an array containing the value energy 
    # repeated a number of times equal to the length of the L, and vice versa, in order to zip them.
    energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
    L = np.full(len(energy), L[0]) if (len(L) == 1) else L

    # The call to __getitem__ below is a way to return a float if both energy and L were floats
    if ((nu_i is not None) and (nu_f is not None)):
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep+h_matt, 0.0, xy[1])[nu_i][nu_f]
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
    else:
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep+h_matt, 0.0, xy[1])
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))


def osc_prob_3nu_vacuum(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray], 
    s12: Optional[float]=None, 
    s23: Optional[float]=None, 
    s13: Optional[float]=None, 
    dCP: Optional[float]=None, 
    D21: Optional[float]=None, 
    D31: Optional[float]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0, 
    **kwargs
) -> Union[float, np.ndarray]:

    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L

    if validate_input:
        # The function name is sys._getframe().f_code.co_name
        if validate_input_battery(sys._getframe().f_code.co_name, energy, L, nu_i, nu_f) == 1:
            sys.exit(1)

    # Flag return_float remembers if energy and L were both floats.  If True, return a float, too.
    return_float = isinstance(energy, float) and isinstance(L, float)

    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)  
    L = np.array([L]) if isinstance(L, float) else np.array(L) 

    # Either energy and L are both lists (or NumPy arrays) of the same length; or one is a float and
    # the other is a list (or NumPy array).  Any other possibility will generate an exception.  This
    # exception is raised earlier if validate_input == True, but we check below in case it has been
    # set to False.
    try:
        if not ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or \
            (len(energy) > 1 and len(L) == 1)):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_3nu_vacuum: energy and " + \
                "L must be both int or float; or, if lists (or NumPy arrays), they must have " + \
                "the same length; or, if one is a float or single-entry list, the other must " + \
                "be a list with multiple entries.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If any of the oscillation parameters has not been given a value, assign to it the value from
    # the specified parameter set with name default_osc_params_set_name.  When input validation is
    # on (validate_input == True), the routine checks whether the parameter set name is among the 
    # predefined ones (see validation above).  Only the values of the parameters passed as None are
    # assigned from the predefined set; other parameters are not modified.
    if ((s12 is None) or (s23 is None) or (s13 is None) or (s23 is None) or (dCP is None) or \
        (D21 is None) or (D31 is None)):

        default_osc_params = gd.OSC_PARAMS_PREDEFINED[default_osc_params_set_name]

        if verbose > 0:
            print(gd.WARNING_MSG_IN_COLOR + " Setting unspecified oscillation parameters to " + \
                "default values from the predefined set " + default_osc_params['name'] + " (" + \
                default_osc_params['description'] + "):\n" + \
                ("s12 = " + str(default_osc_params['s12']) + "\n" if (s12 is None) else '') + \
                ("s23 = " + str(default_osc_params['s23']) + "\n" if (s23 is None) else '') + \
                ("s13 = " + str(default_osc_params['s13']) + "\n" if (s13 is None) else '') + \
                ("dCP = " + str(default_osc_params['dCP']) + " rad\n" if (dCP is None) else '') + \
                ("D21 = " + str(default_osc_params['D21']) + " eV^2\n" if (D21 is None) else '') + \
                ("D31 = " + str(default_osc_params['D31']) + " eV^2\n" if (D31 is None) else ''))

        s12 = s12 if (s12 is not None) else default_osc_params['s12']
        s23 = s23 if (s23 is not None) else default_osc_params['s23']
        s13 = s13 if (s13 is not None) else default_osc_params['s13']
        dCP = dCP if (dCP is not None) else default_osc_params['dCP']
        D21 = D21 if (D21 is not None) else default_osc_params['D21']
        D31 = D31 if (D31 is not None) else default_osc_params['D31']            

    # Compute the energy-independent part of the Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor below when calling osc_prob.
    h_vac_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13,
        dCP, D21, D31, nubar=nubar) 

    # If energy is a single value, then transform it into an array containing the value energy 
    # repeated a number of times equal to the length of the L, and vice versa, in order to zip them.
    energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
    L = np.full(len(energy), L[0]) if (len(L) == 1) else L

    # The call to __getitem__ below is a way to return a float if both energy and L were floats
    if ((nu_i is not None) and (nu_f is not None)):
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep, 0.0, xy[1])[nu_i][nu_f]
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
    else:
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep, 0.0, xy[1])
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))


def osc_prob_3nu_matter_constant_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray], 
    rho: float, 
    s12: Optional[float]=None, 
    s23: Optional[float]=None, 
    s13: Optional[float]=None, 
    dCP: Optional[float]=None, 
    D21: Optional[float]=None, 
    D31: Optional[float]=None, 
    ratio_number_neutrons_to_protons: Optional[float]=1.0, 
    electron_fraction: Optional[float]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:

    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L

    if validate_input:
        # The function name is sys._getframe().f_code.co_name
        if validate_input_battery(sys._getframe().f_code.co_name, energy, L, nu_i, nu_f) == 1:
            sys.exit(1)

    # Flag return_float remembers if energy and L were both floats.  If True, return a float, too.
    return_float = isinstance(energy, float) and isinstance(L, float)

    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)  
    L = np.array([L]) if isinstance(L, float) else np.array(L) 

    # Either energy and L are both lists (or NumPy arrays) of the same length; or one is a float and
    # the other is a list (or NumPy array).  Any other possibility will generate an exception.  This
    # exception is raised earlier if validate_input == True, but we check below in case it has been
    # set to False.
    try:
        if not ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or \
            (len(energy) > 1 and len(L) == 1)):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_3nu_matter_constant_density: energy and L must be both " + \
                "int or float; or, if lists (or NumPy arrays), they must have the same length;" + \
                " or, if one is a float or single-entry list, the other must be a list with " + \
                "multiple entries.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If any of the oscillation parameters has not been given a value, assign to it the value from
    # the specified parameter set with name default_osc_params_set_name.  When input validation is
    # on (validate_input == True), the routine checks whether the parameter set name is among the 
    # predefined ones (see validation above).  Only the values of the parameters passed as None are
    # assigned from the predefined set; other parameters are not modified. 
    if ((s12 is None) or (s23 is None) or (s13 is None) or (s23 is None) or (dCP is None) or \
        (D21 is None) or (D31 is None)):

        default_osc_params = gd.OSC_PARAMS_PREDEFINED[default_osc_params_set_name]

        if verbose > 0:
            print(gd.WARNING_MSG_IN_COLOR + " Setting unspecified oscillation parameters to " + \
                "default values from the predefined set " + default_osc_params['name'] + " (" + \
                default_osc_params['description'] + "):\n" + \
                ("s12 = " + str(default_osc_params['s12']) + "\n" if (s12 is None) else '') + \
                ("s23 = " + str(default_osc_params['s23']) + "\n" if (s23 is None) else '') + \
                ("s13 = " + str(default_osc_params['s13']) + "\n" if (s13 is None) else '') + \
                ("dCP = " + str(default_osc_params['dCP']) + " rad\n" if (dCP is None) else '') + \
                ("D21 = " + str(default_osc_params['D21']) + " eV^2\n" if (D21 is None) else '') + \
                ("D31 = " + str(default_osc_params['D31']) + " eV^2\n" if (D31 is None) else ''))

        s12 = s12 if (s12 is not None) else default_osc_params['s12']
        s23 = s23 if (s23 is not None) else default_osc_params['s23']
        s13 = s13 if (s13 is not None) else default_osc_params['s13']
        dCP = dCP if (dCP is not None) else default_osc_params['dCP']
        D21 = D21 if (D21 is not None) else default_osc_params['D21']
        D31 = D31 if (D31 is not None) else default_osc_params['D31']            

    s = 1.0 if not nubar else -1.0

    # Electron number density [eV^3]
    num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, 
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction, 
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3) 

    # Coherent forward potential, VCC [eV]
    VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) 

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor below when calling osc_prob.
    h_vac_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, 
        dCP, D21, D31, nubar=nubar) 

    # Compute the matter Hamiltonian only once, to save time.
    h_matt = s*hamiltonians3nu.hamiltonian_3nu_matter(VCC)

    # If energy is a single value, then transform it into an array containing the value energy 
    # repeated a number of times equal to the length of the L, and vice versa, in order to zip them.
    energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
    L = np.full(len(energy), L[0]) if (len(L) == 1) else L

    # The call to __getitem__ below is a way to return a float if both energy and L were floats
    if ((nu_i is not None) and (nu_f is not None)):
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep+h_matt, 0.0, xy[1])[nu_i][nu_f]
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
    else:
        return np.array([osc_prob((1/xy[0])*h_vac_energy_indep+h_matt, 0.0, xy[1])
            for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))


if __name__ == "__main__":
    def H_2nu_func(t):
        return np.array([[1+1j*t, 2*t], [2*t, 4-1j*t]], dtype=np.complex128)
    def H_3nu_func(t):
        return np.array([[1+1j*t, 2*t, 3j*t], [2*t, 4-1j*t, 5+2j*t], [-3j*t, 5-2j*t, 1]], 
            dtype=np.complex128)

    t_ini, t_fin = 0.0, 1.0

    # prob = osc_prob(H_2nu_func, t_ini, t_fin, n_slabs=100, n_tpts_per_slab=100, magnus_exp_order=6,
    #     integration_method='simpson', n_jobs=1, save_log=True, verbose=2)
    # print(prob)
    # prob = osc_prob(H_3nu_func, t_ini, t_fin, n_slabs=100, n_tpts_per_slab=100, magnus_exp_order=6,
    #     integration_method='simpson', n_jobs=1)
    # print(prob)

    # prob = osc_prob(H_3nu_func, t_ini, t_fin, n_slabs=10, n_tpts_per_slab=20, magnus_exp_order=4,
    #     integration_method='simpson', n_jobs=10, rtol=1e-5, atol=1.e-5, 
    #     growth_factor_n_slabs=1.5, growth_factor_n_tpts_per_slab=1.5, 
    #     max_num_loops=50, max_n_slabs=200, max_n_tpts_per_slab=150, 
    #     save_log=True, filename_log='./out.log', verbose=2)
    # print(prob)

    # Test iteration over magnus_exp_order
    prob = osc_prob(H_3nu_func, t_ini, t_fin, 
        integration_method='simpson', n_jobs=10, rtol=1e-3, atol=1.e-3, 
        max_n_slabs=10, max_n_tpts_per_slab=10, 
        iterate_over_magnus_exp_order=True, min_magnus_exp_order=1, 
        max_magnus_exp_order=gd.MAGNUS_EXP_ORDER_MAX,
        save_log=True, filename_log='./out.log', verbose=2)
    print(prob)

    # Test use of default values of oscillation parameters: vacuum
    # print(osc_prob_3nu_vacuum(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, 
    #     validate_input=True, verbose=0), end='\n\n')
    # print(osc_prob_3nu_vacuum(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUMU,
    #     validate_input=True, verbose=0), end='\n\n')
    # print(osc_prob_3nu_vacuum(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUMU,
    #     validate_input=True, verbose=1), end='\n\n')
    # print(osc_prob_3nu_vacuum(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUMU,
    #     s13=0.0, s23=0.0, dCP=0.0, D31=0.0, validate_input=True, verbose=1), end='\n\n')
    # print(osc_prob_3nu_vacuum(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUMU,
    #     default_osc_params_set_name='xxx', validate_input=True, verbose=1), end='\n\n')

    # Test use of default values of oscillation parameters: constant-density matter
    # print(osc_prob_3nu_matter_constant_density(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, 
    #     10.0*gd.UNIT_G_PER_CM3, validate_input=True, verbose=0), end='\n\n')
    # print(osc_prob_3nu_matter_constant_density(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, 
    #     10.0*gd.UNIT_G_PER_CM3, nu_i=gd.NUE, nu_f=gd.NUMU, validate_input=True, verbose=0), 
    # end='\n\n')
    # print(osc_prob_3nu_matter_constant_density(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, 
    #     10.0*gd.UNIT_G_PER_CM3, nu_i=gd.NUE, nu_f=gd.NUMU, validate_input=True, verbose=1), 
    # end='\n\n')
    # print(osc_prob_3nu_matter_constant_density(1.*gd.UNIT_MEV, 100*gd.UNIT_KM,
    #     10.0*gd.UNIT_G_PER_CM3, nu_i=gd.NUE, nu_f=gd.NUMU,
    #     s13=0.0, s23=0.0, dCP=0.0, D31=0.0, validate_input=True, verbose=1), end='\n\n')
    # print(osc_prob_3nu_matter_constant_density(1.*gd.UNIT_MEV, 100*gd.UNIT_KM, 
    #     10.0*gd.UNIT_G_PER_CM3, nu_i=gd.NUE, nu_f=gd.NUMU,
    #     default_osc_params_set_name='xxx', validate_input=True, verbose=1), end='\n\n')

    # prob = osc_prob_2nu_vacuum('x', 1.0*gd.UNIT_KM, 0.5, 1.e-4)
    # print(prob)
