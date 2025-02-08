"""Contains routines to compute the neutrino oscillation probability.

Internally, the probability is computed using Magnus expansion, but the
user does not call the routines in the :py:mod:`magnus.magnus` module
directly. Instead, the user calls the :func:`osc_prob`, which calls the
Magnus expansion routines internally. The function :func:`osc_prob` is
generic, flexible, and computationally efficient. 

- :func:`osc_prob`: Primordial function to compute the oscillation
  probability, for any given Hamiltonian, either time-dependent or 
  -independent (or, equivalently, position-dependent or -independent). 
  Supports arbitrary number of neutrino flavors.

The module contains additional functions that are wrappers of 
:func:`osc_prob` to compute commonly studied cases.

Neutrino oscillations in **vacuum**:

- :func:`osc_prob_2nu_vacuum`: Two-neutrino oscillation probabilities.

- :func:`osc_prob_3nu_vacuum`: Three-flavor oscillation probabilities.

- :func:`osc_prob_4nu_vacuum`: One 
  additional flavor (i.e., 3+1 sterile neutrino model).

- :func:`osc_prob_5nu_vacuum`: Two 
  additional flavors (i.e., 3+2 sterile neutrino model).

Neutrino oscillations in **constant-density matter**:

- :func:`osc_prob_2nu_matter_constant_density`: Two-neutrino oscillation
  probabilities.

- :func:`osc_prob_3nu_matter_constant_density`: Three-neutrino
  oscillation probabilities.

- :func:`osc_prob_4nu_matter_constant_density`: One additional flavor 
  (i.e., 3+1 sterile model).

- :func:`osc_prob_5nu_matter_constant_density`: Two additional flavors 
  (i.e., 3+2 sterile model).

Neutrino oscillations in **exponentially falling matter density 
profile** (e.g., in a supernova or the Sun):

- :func:`osc_prob_2nu_matter_exp_density`: Two-neutrino oscillation
  probabilities.

- :func:`osc_prob_3nu_matter_exp_density`: Three-neutrino oscillation
  probabilities.

- :func:`osc_prob_4nu_matter_exp_density`: One additional flavor. Matter
  potential affects only :math:`\\nu_e`.

- :func:`osc_prob_5nu_matter_exp_density`: Two additional flavors. 
  Matter potential affects only :math:`\\nu_e`.

- :func:`osc_prob_matter_exp_density`: Oscillation probabilities for
  arbitrary number of flavors and arbitrary Hamiltonian.  Does not 
  assume standard oscillations.

Neutrino oscillations between any two locations on the surface of the
**Earth**, useful for long-baseline neutrino experiments:

- :func:`osc_prob_2nu_earth`: Two-neutrino oscillation probabilities.  

- :func:`osc_prob_3nu_earth`: Three-neutrino oscillation probabilities. 

- :func:`osc_prob_4nu_earth`: One additional flavor. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_5nu_earth`: Two additional flavors. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_earth`: Oscillation probabilities for arbitrary number
  of flavors and arbitrary Hamiltonian.  Does not assume standard 
  oscillations.

.. note::
   These routines use the `Preliminary Reference Earth Model 
   <https://www.cfa.harvard.edu/~lzeng/papers/PREM.pdfL>`_ for the 
   matter density profile inside Earth.  To use a different density 
   profile (including also profiles for bodies other than the Earth), 
   use the primordial function :func:`osc_prob` instead.

Neutrino oscillations in the **Sun**:

- :func:`osc_prob_2nu_sun`: Two-neutrino oscillation probabilities.  

- :func:`osc_prob_3nu_sun`: Three-neutrino oscillation probabilities. 

- :func:`osc_prob_4nu_sun`: One additional flavor. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_5nu_sun`: Two additional flavors. Matter potential 
  affects only :math:`\\nu_e`.

- :func:`osc_prob_sun`: Oscillation probabilities for arbitrary number
  of flavors and arbitrary Hamiltonian.  Does not assume standard 
  oscillations.

.. note::
   These routines use a simple exponentially falling function of radial
   distance for the matter density inside the Sun: :math:`N_e(r) = 
   N_e(0) \\exp(-r/r_0)`, with 
   :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
   :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
   `Fundamentals of Neutrino Physics and Astrophysics 
   <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
   Wook Kim.

   To use a different density profile, use the primordial function 
   :func:`osc_prob` instead.

Functions designed for specific **beyond-the-Standard-Model** proposals:

- Non-standard neutrino interactions (NSI):

  - :func:`osc_prob_2nu_matter_constant_density_nsi`
  - :func:`osc_prob_3nu_matter_constant_density_nsi`
  - :func:`osc_prob_2nu_matter_exp_density_nsi`
  - :func:`osc_prob_3nu_matter_exp_density_nsi`
  - :func:`osc_prob_2nu_matter_earth_nsi`
  - :func:`osc_prob_3nu_matter_earth_nsi`
  - :func:`osc_prob_2nu_matter_sun_nsi`
  - :func:`osc_prob_3nu_matter_sun_nsi`

- Lorentz-invariance violation:

  - :func:`osc_prob_2nu_vacuum_liv`
  - :func:`osc_prob_3nu_vacuum_liv`
  - :func:`osc_prob_2nu_matter_constant_density_liv`
  - :func:`osc_prob_3nu_matter_constant_density_liv`
  - :func:`osc_prob_2nu_matter_exp_density_liv`
  - :func:`osc_prob_3nu_matter_exp_density_liv`
  - :func:`osc_prob_2nu_earth_liv`
  - :func:`osc_prob_3nu_earth_liv`
  - :func:`osc_prob_2nu_sun_liv`
  - :func:`osc_prob_3nu_sun_liv`

Examples
--------

.. seealso::
   Find many more examples, including advanced applications and plots,
   in the `Jupyter notebooks 
   <https://github.com/mbustama/Magnus/tree/main/notebooks>`_ that are 
   distributed with :math:`{\\rm Mag}{\\nu}s`.

>>> import magnus.oscprob as oscprob
>>> import magnus.globaldefs as gd

Calling :func:`osc_prob_3nu_vacuum` returns a :math:`3 \\times 3` NumPy array
with entries XXX

For a single neutrino energy and baseline:

>>> baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
>>> energy = 1.*gd.UNIT_MEV # [eV]
>>> np.set_printoptions(precision=3)
>>> oscprob.osc_prob_3nu_vacuum(energy, baseline)
[[0.445 0.299 0.257]
 [0.251 0.639 0.11 ]
 [0.304 0.062 0.634]]

The probabilities returned by :func:`osc_prob_3nu_vacuum` (and also
:func:`osc_prob_2nu_vacuum`, 
:func:`osc_prob_2nu_matter_constant_density`, and
:func:`osc_prob_3nu_matter_constant_density`) are returned with machine
(or NumPy) precision, since first-order Magnus expansion is enough to 
compute them.

Pick one channel only, e.g., :math:`\\nu_e \\to \\nu_\\mu`, by passing
an initial flavor, ``nu_i``, and a final flavor ``nu_f``:

>>> oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, \
nu_f=gd.NUMU)

The flavor indices ``NUE``, ``NUMU``, and ``NUMU`` are defined in the 
:py:mod:`magnus.globaldefs` module. For anti-neutrinos, i.e., 
:math:`\\bar{\\nu}_e \\to \\bar{\\nu}_\\mu`:

>>> oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, \
nu_f=gd.NUMU, nubar=True)

Calling :func:`osc_prob_3nu_vacuum` without specifying the values of the
oscillation parameters will compute probabilities using the default 
values in :math:`{\\rm Mag}{\\nu}s` (see 
``gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']``.)

We can specify values of the oscillation parameters. Unspecified values
are set to their defaults (pass nonzero ``verbose`` to see this and 
other warnings):

>>> oscprob.osc_prob_3nu_vacuum(energy, baseline, s12=0.0, verbose=1)
Warning: Setting unspecified oscillation parameters to default values \
from the predefined set\n OSC_PARAMS_NU_FIT_6_0_NO (NuFit 6.0, NO, with \
SK atmospheric data):
s23 = 0.6855654600401044
s13 = 0.14882876066137216
dCP = 3.7000980142279785 rad
D21 = 7.49e-05 eV^2
D31 = 0.002513 eV^2
<BLANKLINE>
[[0.985 0.007 0.008]
 [0.007 0.736 0.257]
 [0.008 0.257 0.735]]

Fixed energy, multiple baselines:

>>> baselines = gd.UNIT_KM*np.array([1.0, 10.0 100.0])
>>> oscprob.osc_prob_3nu_vacuum(energy, baselines, nu_i=gd.NUE, \
nu_f=gd.NUMU)

Fixed baseline, multiple energies:

>>> energies = gd.UNIT_MEV*np.array([1.0, 10.0, 100.0])
>>> oscprob.osc_prob_3nu_vacuum(energies, baseline, nu_i=gd.NUE, \
nu_f=gd.NUMU)

Multiple energies and baselines:

>>> oscprob.osc_prob_3nu_vacuum(energies, baselines, nu_i=gd.NUE, \
nu_f=gd.NUMU)

To compute the oscillation probabilities in constant-density matter, we
need to specify the matter density, ``rho``, i.e.,

>>> rho = 10.0*gd.UNIT_G_PER_CM3
>>> osc_prob_3nu_matter_constant_density(energy, baseline, rho, \
nu_i=gd.NUE, nu_f=gd.NUMU)

To compute oscillation probabilities for a time-dependent Hamiltonian,
we need to call :func:`osc_prob` directly which, while still 
straightforward, requires us to pass a Hamiltonian function explicitly.

For instance, for density matter profile that is exponentially falling
with distance:

.. hint::
   There is a good chance that the scenario you are interested in 
   calculating was already developed in the :math:`{\\rm Mag}{\\nu}s` 
   `Jupyter 
   notebooks 
   <https://github.com/mbustama/Magnus/tree/main/notebooks>`_.  
   
   Worked-out examples include: oscillations in various matter density 
   profiles, in the Earth, and in the Sun, oscillograms, biprobability 
   plots, and new-physics models like additional neutrino flavors (3+1 
   and 3+2 sterile neutrino models), non-standard neutrino interactions,
   and Lorentz-invariance violation.

"""

__version__ = '0.1'
__author__ = 'Mauricio Bustamante'

import numpy as np
import sys
import platform
from joblib import Parallel, delayed
from typing import Optional, Callable, Union, Tuple, List
from io import TextIOWrapper
from inspect import signature

# import numpy.typing

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
    t_ini: Union[int, float], 
    t_fin: Union[int, float], 
    n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=None, 
    atol: Optional[Union[int, float]]=None, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=2000, 
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


def values_to_unspecified_osc_params(
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    verbose: Optional[int]=0
) -> Tuple[float, float, float, float, float, float]:
    r"""Return values of unspecified standard oscillation parameters

    If any of the oscillation parameters has not been given a value, assign to it the value from
    the specified parameter set with name default_osc_params_set_name.  When input validation is
    on (validate_input == True), the routine checks whether the parameter set name is among the 
    predefined ones (see validation above).  Only the values of the parameters passed as None are
    assigned from the predefined set; other parameters are not modified.
    """

    try:
        if not (default_osc_params_set_name in list(gd.OSC_PARAMS_PREDEFINED.keys()) ):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprobvalues_to_unspecified_osc_params:"+ \
                ": the requested oscillation parameter set (default_osc_params_set_name = " + \
                default_osc_params_set_name + ") is not among the predefined sets in Magnus. " + \
                "Available sets are " + str(list(gd.OSC_PARAMS_PREDEFINED.keys())) + ".")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

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

    return s12, s23, s13, dCP, D21, D31


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
    :param magnus_exp_order: Maximum order of Magnus expansion used to compute the evolution
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
    t_ini: Union[int, float], 
    t_fin: Union[int, float], 
    n_slabs: Optional[int]=1, 
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None, 
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-3, 
    atol: Optional[Union[int, float]]=1.e-3, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=2000, 
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
        [``t_ini``, ``t_fin``] is partitioned in order to compute the 
        neutrino evolution operators. A higher value of ``n_slabs`` 
        yields a more accurate probability.

        If no target tolerance is requested (i.e., if ``rtol`` and 
        ``atol`` are both ``None``), then the given value of ``n_slabs`` 
        is the final number of slabs used in the computation.

        If a target tolerance is requested (i.e., if either ``rtol`` or 
        ``atol`` is not ``None``), then the given value of ``n_slabs`` 
        is ignored. Instead, the number of slabs is increased 
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
        NumPy array containing the probability matrix of the same 
        dimensions as the Hamiltonian, ``H_func``.
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
    t_ini: Union[int, float], 
    t_fin: Union[int, float],
    n_slabs: Optional[int]=1,
    n_tpts_per_slab: Optional[int]=100, 
    t_slab_edges: Optional[Union[list, np.ndarray]]=None,
    magnus_exp_order: Optional[int]=4, 
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3, 
    growth_factor_n_slabs: Optional[Union[int, float]]=1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=1.5, 
    max_num_loops: Optional[int]=50,
    min_n_slabs: Optional[int]=1, 
    max_n_slabs: Optional[int]=2000, 
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
    sth: Union[int, float], 
    Dm2: Union[int, float], 
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

    .. seealso::
        :func:`osc_prob_3nu_vacuum`
            Three-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_4nu_vacuum`
            Four-flavor (3+1) oscillation probabilities in vacuum. 
        :func:`osc_prob_5nu_vacuum`
            Four-flavor (3+2) oscillation probabilities in vacuum. 
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


def osc_prob_3nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    vacuum.

    By default returns :math:`3 \times 3` probability matrices for all 
    the  oscillation channels. Each matrix has shape ``np.ndarray([Pee,
    Pem,Pet],[Pme,Pmm,Pmt],[Pte,Ptm,Ptt])``.  The matrix is symmetric, 
    i.e., ``Pme == Pee``, ``Pte == Pet``, and ``Ptm == Pmt``.  

    If a single energy and baseline is given, the function returns a 
    single matrix.  If multiple energies and baselines are given, 
    function returns an NumPy array of matrices.  See examples below.

    If the probability needs to be computed multiple times, it is 
    recommended to pass the array of energies and the array of baselines
    to the function in a single call instead of calling the function
    separately for each combination of energy and baseline. The reason
    is that the function has an overhead that gets diluted when 
    computing when the input energies and baselines are many.

    If the initial and final flavors, ``nu_i`` and ``nu_f``, are 
    specified (by setting them to ``NUE``, ``NUMU``, or ``NUTAU``
    from the :py:mod:`magnus.globaldefs` module), the function returns 
    instead a one-dimensional array of the probabilities computed for
    each value of energy and baseline requested. See examples below.

    If the function is called without specifying values of the standard
    oscillation parameters (``s12``, ``s23``, ``s13``, ``dCP``, ``D21``,
    ``D31``), the unspecified parameters are assigned default values 
    taken from a predefined parameter set.  The name of the default 
    parameter set can be changed by passing 
    ``default_osc_params_set_name``.  

    The names of the predefined parameter sets included in 
    :math:`\text{Mag}\nu\text{s}` can be seen by printing

    >>> import magnus.globaldefs as gd
    >>> list(gd.OSC_PARAMS_PREDEFINED.keys())

    And the default parameter values are from the set with name 
    ``'OSC_PARAMS_DEFAULT'``:

    >>> gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']

    If ``validate_input`` is set to True, the function validated the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    Parameters
    ----------
    energy
        Neutrino energy, single value or array.
    L
        Neutrino baseline, single value or array.
    s12
        Sine of the mixing angle :math:`\theta_{12}`.
    s23
        Sine of the mixing angle :math:`\theta_{23}`.
    s13
        Sine of the mixing angle :math:`\theta_{13}`.
    dCP
        CP-violation phase, :math:`\delta_\text{CP}`.
    D21
        Mass-squared difference :math:`\Delta m_{21}^2`.
    D31
        Mass-squared difference :math:`\Delta m_{31}^2`.
    nubar
        False (default) for neutrinos; True for anti-neutrinos.
    nu_i
        Initial neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    nu_f
        Final neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    default_osc_params_set_name
        Name of the predefined set of oscillation parameters to use when
        assigning default values to unspecified parameters.
    validate_input
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose
        0 not to print warnings and errors; 1 to print them.

    Returns
    -------
    Union[float, np.narray]
        Neutrino oscillation probability matrix or probability for a 
        single oscillation channel, for the values of `energy` and `L`.

    Examples
    --------
    >>> import magnus.oscprob as oscprob
    >>> import magnus.globaldefs as gd

    If both ``energy`` and ``L`` are single values (``float`` or 
    ``int``), this function returns the probability computed at these
    values.  

    >>> baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    >>> energy = 1.*gd.UNIT_MEV # [eV]
    >>> oscprob.osc_prob_3nu_vacuum(energy, baseline)
 
    Pick one channel only, e.g., :math:`\nu_e \to \nu_\mu`, by 
    passing an initial flavor, ``nu_i``, and a final flavor ``nu_f``:

    >>> oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU)

    The flavor indices ``NUE``, ``NUMU``, and ``NUMU`` are defined in 
    the :py:mod:`magnus.globaldefs` module. For anti-neutrinos, i.e., 
    :math:`\bar{\nu}_e \to \bar{\nu}_\mu`:

    >>> oscprob.osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUMU, nubar=True)

    We can specify values of the oscillation parameters. Unspecified 
    values are set to their defaults (pass nonzero ``verbose`` to see 
    this and other warnings):

    >>> oscprob.osc_prob_3nu_vacuum(energy, baseline, s12=0.0, verbose=1)
    Warning: Setting unspecified oscillation parameters to default 
    values from the predefined set\n OSC_PARAMS_NU_FIT_6_0_NO (NuFit \
    6.0, NO, with SK atmospheric data):
    s23 = 0.6855654600401044
    s13 = 0.14882876066137216
    dCP = 3.7000980142279785 rad
    D21 = 7.49e-05 eV^2
    D31 = 0.002513 eV^2
    <BLANKLINE>
    [[0.985 0.007 0.008]
     [0.007 0.736 0.257]
     [0.008 0.257 0.735]]

    If a single energy value and multiple baselines are passed, this
    function returns an array containing the probabilities computed for
    this fixed energy and each value of the baseline:
    
    >>> baselines = gd.UNIT_KM*np.array([1.0, 10.0 100.0])
    >>> oscprob.osc_prob_3nu_vacuum(energy, baselines, nu_i=gd.NUE, nu_f=gd.NUMU)

    Conversely, if a single baseline and multiple energies are passed,
    this function returns an array containing the probabilities computed
    for this fixed baseline and each value of the energy:

    >>> energies = gd.UNIT_MEV*np.array([1.0, 10.0, 100.0])
    >>> oscprob.osc_prob_3nu_vacuum(energies, baseline, nu_i=gd.NUE, nu_f=gd.NUMU)

    And, for multiple energies and baselines:

    >>> oscprob.osc_prob_3nu_vacuum(energies, baselines, nu_i=gd.NUE, nu_f=gd.NUMU)

    .. seealso::
        :func:`osc_prob_2nu_vacuum`
            Two-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_4nu_vacuum`
            Four-flavor (3+1) oscillation probabilities in vacuum. 
        :func:`osc_prob_5nu_vacuum`
            Four-flavor (3+2) oscillation probabilities in vacuum. 
    """
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
    # the specified parameter set with name default_osc_params_set_name.  Only the values of the 
    # parameters passed as None are assigned from the predefined set; others are not modified.
    s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, D31, 
        default_osc_params_set_name, verbose)

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


def osc_prob_4nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in vacuum.
    """

    pass 

    return


def osc_prob_5nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True,
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in vacuum.
    """

    pass 
    
    return


def osc_prob_2nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    rho: Union[int, float], 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile.
    """

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


def osc_prob_3nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile.
    """

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
    # the specified parameter set with name default_osc_params_set_name.  Only the values of the 
    # parameters passed as None are assigned from the predefined set; others are not modified.
    s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, D31, 
        default_osc_params_set_name, verbose)

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


def osc_prob_4nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with a constant density profile.
    """
    pass

    return


def osc_prob_5nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with a constant density profile.
    """
    pass

    return


def osc_prob_2nu_matter_exp_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation 
    probability in matter with an exponentially falling density profile.
    """
    pass

    return


def osc_prob_3nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation 
    probability in matter with an exponentially falling density profile.
    """
    pass

    return


def osc_prob_4nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with an exponentially falling density profile.
    """
    pass

    return


def osc_prob_5nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with an exponentially falling density profile.
    """
    pass

    return


def osc_prob_matter_exp_density(
    H_func: Union[Callable, np.ndarray],
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    r"""Compute and return the neutrino oscillation probability in 
    matter with an exponentially falling density profile for a given
    arbitrary Hamiltonian.
    """
    pass

    return


def osc_prob_2nu_earth(
    ra_dec_ini: Union[Tuple[float, float], list, np.ndarray, str], 
    ra_dec_fin: Union[Tuple[float, float], list, np.ndarray, str], 
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the two-neutrino oscillation probability 
    between two locations on the surface of the Earth.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass

    return 


def osc_prob_3nu_earth(
    ra_dec_ini: Union[Tuple[float, float], list, np.ndarray, str], 
    ra_dec_fin: Union[Tuple[float, float], list, np.ndarray, str], 
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the three-neutrino oscillation probability 
    between two locations on the surface of the Earth.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass

    return 


def osc_prob_4nu_earth(
    ra_dec_ini: Union[Tuple[float, float], list, np.ndarray, str], 
    ra_dec_fin: Union[Tuple[float, float], list, np.ndarray, str], 
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the four-neutrino (3+1) oscillation 
    probability between two locations on the surface of the Earth.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass


def osc_prob_5nu_earth(
    ra_dec_ini: Union[Tuple[float, float], list, np.ndarray, str], 
    ra_dec_fin: Union[Tuple[float, float], list, np.ndarray, str], 
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the five-neutrino (3+2) oscillation 
    probability between two locations on the surface of the Earth.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass


def osc_prob_earth(
    H_func: Union[Callable, np.ndarray],
    ra_dec_ini: Union[Tuple[float, float], list, np.ndarray, str], 
    ra_dec_fin: Union[Tuple[float, float], list, np.ndarray, str], 
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the neutrino oscillation probability between 
    two locations on the surface of the Earth for a given arbitrary 
    Hamiltonian.

    Does **not** assume standard oscillations nor a given number of 
    neutrino flavors: the user must supply their own Hamiltonian 
    function, ``H_func``.  The Hamiltonian can include matter potentials
    due not only to electrons and that affect not only :math:`\\nu_e`.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass

    return 


def osc_prob_2nu_sun(
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the two-neutrino oscillation probability 
    for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass

    return 


def osc_prob_3nu_sun(
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the three-neutrino oscillation probability 
    for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass

    return 


def osc_prob_4nu_sun(
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the four-neutrino (3+1) oscillation 
    probability for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass


def osc_prob_5nu_sun(
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the five-neutrino (3+2) oscillation 
    probability for neutrinos inside the Sun.

    Assumes that the matter potential is due only to the standard 
    charged-current coherent forward scattering of :math:`\\nu_e` on 
    electrons.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """

    pass


def osc_prob_sun(
    H_func: Union[Callable, np.ndarray],
    costhz: Union[int, float],
    energy: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    verbose: Optional[int]=0
) -> Union[float, np.ndarray]:
    """Compute and return the neutrino oscillation probability for 
    neutrinos inside the Sun for a given arbitrary Hamiltonian.

    Does **not** assume standard oscillations nor a given number of 
    neutrino flavors: the user must supply their own Hamiltonian 
    function, ``H_func``.  The Hamiltonian can include matter potentials
    due not only to electrons and that affect not only :math:`\\nu_e`.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    The location of each point on the surface can either be given as a 
    tuple (ra, dec) of right ascension and declination (i.e., 
    ``ra_dec_ini`` and ``ra_dec_fin``)---including the possibility of 
    using locations predefined in Magnus---or as the cosine of the
    zenith angle between the initial and final positions (i.e.,
    ``costhz``).

    Examples
    --------
    """


def osc_prob_matter_std_potential():

    pass

    return



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

    # # Test iteration over magnus_exp_order
    # prob = osc_prob(H_3nu_func, t_ini, t_fin, 
    #     integration_method='simpson', n_jobs=10, rtol=1e-3, atol=1.e-3, 
    #     max_n_slabs=10, max_n_tpts_per_slab=10, 
    #     iterate_over_magnus_exp_order=True, min_magnus_exp_order=1, 
    #     max_magnus_exp_order=gd.MAGNUS_EXP_ORDER_MAX,
    #     save_log=True, filename_log='./out.log', verbose=2)
    # print(prob)

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

    baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    energy = 1.*gd.UNIT_MEV # [eV]
    np.set_printoptions(precision=3)
    prob = osc_prob_3nu_vacuum(energy, baseline)
    print(prob)
    print(gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'])
    np.set_printoptions(precision=3)
    print(osc_prob_3nu_vacuum(energy, baseline, s12=0.0, verbose=1))

