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

  - :func:`osc_prob_2nu_matter_nsi_constant_density`
  - :func:`osc_prob_3nu_matter_nsi_constant_density`
  - :func:`osc_prob_4nu_matter_nsi_constant_density`
  - :func:`osc_prob_5nu_matter_nsi_constant_density`
  - :func:`osc_prob_2nu_matter_nsi_exp_density`
  - :func:`osc_prob_3nu_matter_nsi_exp_density`
  - :func:`osc_prob_4nu_matter_nsi_exp_density`
  - :func:`osc_prob_5nu_matter_nsi_exp_density`
  - :func:`osc_prob_2nu_earth_nsi`
  - :func:`osc_prob_3nu_earth_nsi`
  - :func:`osc_prob_4nu_earth_nsi`
  - :func:`osc_prob_5nu_earth_nsi`
  - :func:`osc_prob_2nu_sun_nsi`
  - :func:`osc_prob_3nu_sun_nsi`
  - :func:`osc_prob_4nu_sun_nsi`
  - :func:`osc_prob_5nu_sun_nsi`

- Lorentz-invariance violation:

  - :func:`osc_prob_2nu_matter_liv_constant_density`
  - :func:`osc_prob_3nu_matter_liv_constant_density`
  - :func:`osc_prob_4nu_matter_liv_constant_density`
  - :func:`osc_prob_5nu_matter_liv_constant_density`
  - :func:`osc_prob_2nu_matter_liv_exp_density`
  - :func:`osc_prob_3nu_matter_liv_exp_density`
  - :func:`osc_prob_4nu_matter_liv_exp_density`
  - :func:`osc_prob_5nu_matter_liv_exp_density`
  - :func:`osc_prob_2nu_earth_liv`
  - :func:`osc_prob_3nu_earth_liv`
  - :func:`osc_prob_4nu_earth_liv`
  - :func:`osc_prob_5nu_earth_liv`
  - :func:`osc_prob_2nu_sun_liv`
  - :func:`osc_prob_3nu_sun_liv`
  - :func:`osc_prob_4nu_sun_liv`
  - :func:`osc_prob_5nu_sun_liv`

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
>>> import numpy as np
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

__version__ = '0.10'
__author__ = 'Mauricio Bustamante'


import numpy as np
import sys
import platform
from joblib import Parallel, delayed
from typing import Optional, Callable, Union, Tuple, List, Dict
from io import TextIOWrapper
from inspect import signature
# from numba import jit

# import numpy.typing

# TO-DO: remove this once setup.py and pip are working
import os
sys.path.append(os.path.split(os.path.split(os.getcwd())[0])[0])
sys.path.append(os.path.split(os.getcwd())[0])

import magnus.magnus as magnus
import magnus.globaldefs as gd
import magnus.hamiltonians.hamiltonians2nu as hamiltonians2nu
import magnus.hamiltonians.hamiltonians3nu as hamiltonians3nu
import magnus.hamiltonians.hamiltonians4nu as hamiltonians4nu
import magnus.hamiltonians.hamiltonians5nu as hamiltonians5nu
import magnus.matter as matter
import magnus.earth as earth
import version as version
import authors as authors


has_magnus_header_been_printed = False


#-----------------------------------------------------------------------
# Helper functions
#-----------------------------------------------------------------------

def print_banner(file: TextIOWrapper=None):
    if file is None:
        print(gd.cstyle.CBLUEBG + ".----------------------------------------." + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|   __  __                               |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|  |  \/  | __ _  __ _ _ __  _   _ ___   |" + gd.cstyle.CEND, 
            file=file)
        print(gd.cstyle.CBLUEBG + "|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|  | |  | | (_| | (_| | | | | |_| \__ \  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|                |___/                   |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "'----------------------------------------'" + gd.cstyle.CEND,
            file=file)
    else: 
        print(".----------------------------------------.", file=file)
        print("|   __  __                               |", file=file)
        print("|  |  \/  | __ _  __ _ _ __  _   _ ___   |", file=file)
        print("|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |", file=file)
        print("|  | |  | | (_| | (_| | | | | |_| \__ \  |", file=file)
        print("|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |", file=file)
        print("|                |___/                   |", file=file)
        print("'----------------------------------------'", file=file)
    print("Version: "+ version.__version__ + " | Author(s): " + authors.__authors__ + "\n",
        file=file)


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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0, 
    file_log: Optional[TextIOWrapper]=None
):

    global has_magnus_header_been_printed

    for f in [None, file_log] if save_log else [None]:
        if not has_magnus_header_been_printed:
            print_banner(f)
            has_magnus_header_been_printed = True
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
            print("   t_slab_edges = None", file=f)
        else:
            print("   t_slab_edges = ", file=f)
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
        print("   new_recursion_limit = " + str(new_recursion_limit), file=f)
        print("   verbose = " + str(verbose), file=f)

    return


def validate_input_battery(
    source_func_name: str, 
    energy: Optional[Union[int, float, list, np.ndarray]]=None, 
    L: Optional[Union[int, float, list, np.ndarray]]=None, 
    L0: Optional[Union[int, float]]=None,
    num_flavors: Optional[int]=None,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    osc_params: Optional[Union[list, np.ndarray]]=None,
    rho_func: Optional[Union[Callable, int, float]]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    validate_energy_and_L: Optional[bool]=True,
    validate_flavor_indices: Optional[bool]=True,
    validate_osc_params: Optional[bool]=True,
    validate_initial_position: Optional[bool]= False,
    validate_density: Optional[bool]=False
) -> int:

    if validate_energy_and_L:

        try:
            if ( (not isinstance(energy, int)) and (not isinstance(energy, float)) and \
                (not isinstance(energy, list)) and (not isinstance(energy, np.ndarray)) ):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": energy must be an int, a float, a 1D list, or a 1D NumPy array.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) and \
                (np.array(energy).ndim != 1) ):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": if energy is a list or NumPy array, it must be 1D.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if ( (isinstance(energy, list) or isinstance(energy, np.ndarray)) ):
                if not (np.all(np.array(energy).dtype == np.float_) or \
                    np.all(np.array(energy).dtype == np.int_)):
                    raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                        ": since energy is a list or NumPy array, all of its elements must be int" + \
                        " or float.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if ( (not isinstance(L, int)) and (not isinstance(L, float)) and \
                (not isinstance(L, list)) and (not isinstance(L, np.ndarray)) ):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": L must be an int, a float, a 1D list, or a 1D NumPy array.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if ( (isinstance(L, list) or isinstance(L, np.ndarray)) and \
                (np.array(L).ndim != 1) ):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": if L is a list or NumPy array, it must be 1D.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if ( (isinstance(L, list) or isinstance(L, np.ndarray)) ):
                if not (np.all(np.array(L).dtype == np.float_) or np.all(np.array(L).dtype == np.int_)):
                    raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                        ": since L is a list or NumPy array, all of its elements must be int or float.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

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

        try:
            if (((nu_i is not None) and (nu_f is None)) or ((nu_i is None) and (nu_f is not None))):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": if either nu_i or nu_f is not None, the other flavor must also be not None.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

    if validate_flavor_indices:

        try:
            if ((nu_i is not None) and (nu_f is not None)):
                if (num_flavors <= gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
                    if ((num_flavors == 2) or (num_flavors == 3)):
                        flavors = set([gd.NUE, gd.NUMU, gd.NUTAU])
                    elif (num_flavors == 4):
                        flavors = set([gd.NUE, gd.NUMU, gd.NUTAU, gd.NUS])
                    elif (num_flavors == 5):
                        flavors = set([gd.NUE, gd.NUMU, gd.NUTAU, gd.NUS1, gd.NUS2])
                    if ((not (nu_i in flavors)) or (not (nu_f in flavors))):
                        if ((num_flavors == 2) or (num_flavors == 3)):
                            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                                ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                                str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), or gd.NUTAU (" + \
                                str(gd.NUTAU) + ") only.")
                        elif (num_flavors == 4):
                            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                                ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                                str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), gd.NUTAU (" + \
                                str(gd.NUTAU) + "), or gd.NUS (" + str(gd.NUS) + ") only.")
                        elif (num_flavors == 5):
                            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                                ": if nu_i and nu_f are not None, they must be either gd.NUE (" + \
                                str(gd.NUE) + "), gd.NUMU (" + str(gd.NUMU) + "), gd.NUTAU (" + \
                                str(gd.NUTAU) + "), gd.NUS1 (" + str(gd.NUS1) + "), or gd.NUS2 (" + \
                                str(gd.NUS2) + ") only.")
                else:
                    print(gd.WARNING_MSG_IN_COLOR + " " + source_func_name + \
                        ": nu_i and nu_f are not None, but, since num_flavors = " + str(num_flavors) + \
                        " > globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
                        str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + ", input validation cannot " + \
                        "check if nu_e and nu_f are valid indices.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

    if validate_osc_params:

        try:
            ttest = [(isinstance(x, int) or isinstance(x, float) or (x is None)) 
                for x in osc_params]
            if (not np.all(ttest)):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ":"+\
                    " the oscillation parameters must be int or float.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

    if validate_initial_position:

        try:
            if not ((isinstance(L0, int) or (isinstance(L0, float)))):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                    " the initial neutrino position (L0) must be an int or float.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

    if validate_density:

        try:
            if (ratio_number_neutrons_to_protons < 0.0):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                    " the ratio of neutrinos to protons (ratio_number_neutrons_to_protons) must" + \
                    " be non-negative.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if (electron_fraction < 0.0):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                    " the ratio of electrons to protons + neutrons (electron_fraction) must be " + \
                    "non-negative.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if ((callable(rho_func)) and (len(signature(rho_func).parameters) > 1)):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                    " the provided rho_func is a function of more than one parameter.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        rho_test = rho_func(L0) if callable(rho_func) else rho_func

        try:
            if (rho_test < 0.0):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                    " rho_func must be non-negative.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

        try:
            if not (isinstance(rho_test, int) or isinstance(rho_test, float)):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_matter_std_potential:"+\
                    " rho_func must be a float (or int) or must return a float (or int).")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            return 1

    return 0


def validate_input_osc_prob_earth(
    source_func_name: str,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    costhz: Optional[Union[int, float]]=None,
    L: Optional[Union[float, list, np.ndarray]]=None,
    verbose: Optional[int]=0,
    ) -> Tuple[float, np.ndarray]:

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given.
    # If only a single location is given, throw an exception.  If neither of the two locations are
    # given, use the given value of costhz and of baseline given (could be an array of baselines).
    
    if ( ((loc_ini is None) and (loc_fin is not None)) or \
        ((loc_ini is not None) and (loc_fin is None)) ):

        print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": only of the two " + \
            "locations on Earth (loc_ini or loc_fin) has been given. If one location is " + \
            "given (i.e., is not None), the other one must also be given.  Alternatively, " + \
            "both locations can be set to None, and the given value of costhz will be used " +\
            "(if it is not None).")
        print("Aborting execution...")
        sys.exit(1)

    elif ((loc_ini is not None) and (loc_fin is not None)):

        # Check that the location is a two-entry tuple, list, or array

        try:
            lat_ini, lon_ini = loc_ini
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": if the initial " + \
                    "location (loc_ini) is given as coordinates, it must be a two-entry tuple," + \
                    " list, or NumPy array.")
            print("Aborting execution...")
            sys.exit(1)
    
        try:
            lat_fin, lon_fin = loc_fin
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": if the final " + \
                    "location (loc_fin) is given as coordinates, it must be a two-entry tuple," + \
                    " list, or NumPy array.")
            print("Aborting execution...")
            sys.exit(1)

        # We use the function earth.costhz_between_points_on_surface to compute the cosine of the
        # zenith angle of the chord that joins two locations on the surface of the Earth, measured 
        # at one position (any of the two locations will give the same result).
        costhz = earth.costhz_between_points_on_surface(lat_ini, lon_ini, lat_fin, lon_fin)

        # Length of the chord is the baseline
        L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM # [eV^{-1}]

        if verbose > 0:
            print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": using as " + \
                "baseline the chord that joins the given initial and final locations on the " + \
                "surface of the Earth.")

        return costhz, L

    else: # (loc_ini is None) and (loc_fin is None)

        try:
            if costhz is None:
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": no" + \
                    " initial and final locations on the surface of the Earth give, and no " + \
                    "value of costhz given.  This function requires either the two locations " + \
                    "or, alternatively, the value of costhz.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        try:
            if L is None:
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": since two locations on the surface of the Earth have not been given, " + \
                    "the value of costhz wil be used to define the chord lengh, but the" + \
                    " baseline, L, cannot be None.")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

        return costhz, L


def valid_flavor_indices_2nu(nu_i: int, nu_f: int) -> Tuple[int, int]:

    if ((nu_i == gd.NUE) and (nu_f == gd.NUTAU)):
        nu_f = 1
    elif ((nu_i == gd.NUTAU) and (nu_f == gd.NUE)):
        nu_i = 1
    elif ((nu_i == gd.NUMU) and (nu_f == gd.NUTAU)):
        nu_i, nu_f = 0, 1
    elif ((nu_i == gd.NUTAU) and (nu_f == gd.NUMU)):
        nu_i, nu_f = 1,0

    return nu_i, nu_f


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

    global has_magnus_header_been_printed

    if ((s12 is None) or (s23 is None) or (s13 is None) or (s23 is None) or (dCP is None) or \
        (D21 is None) or (D31 is None)):

        default_osc_params = gd.OSC_PARAMS_PREDEFINED[default_osc_params_set_name]

        if verbose > 0:
            if verbose >= 2:
                if (not has_magnus_header_been_printed):
                    print_banner()
                    has_magnus_header_been_printed = True
            print(gd.WARNING_MSG_IN_COLOR + " Setting unspecified standard oscillation " + \
                "parameters to default values from the predefined set " + \
                default_osc_params['name'] + " (" + default_osc_params['description'] + "):\n" + \
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


def unpack_oscillation_params_from_dict(
    source_func_name: str,
    num_flavors: int,
    osc_params: Dict,
    h_vac_energy_indep: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Unpack oscillation parameters from the osc_params dict
    """

    if (num_flavors == 2):
        try:
            sth = osc_params['sth']
            Dm2 = osc_params['Dm2']
            return np.array([sth, Dm2])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since "+ \
                    "num_flavors == 2, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 'sth' and 'Dm2'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 3):
        try:
            s12 = osc_params['s12']
            s23 = osc_params['s23']
            s13 = osc_params['s13']
            dCP = osc_params['dCP']
            D21 = osc_params['D21']
            D31 = osc_params['D31']
            return np.array([s12, s23, s13, dCP, D21, D31])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 3, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 's12', 's23', 's13', 'dCP', 'D21', and " + \
                    "'D31', even if they are None.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 4):
        try:
            s12 = osc_params['s12']
            s23 = osc_params['s23']
            s13 = osc_params['s13']
            dCP = osc_params['dCP']
            s14 = osc_params['s14']
            d14 = osc_params['d14']
            s24 = osc_params['s24']
            d24 = osc_params['d24']
            s34 = osc_params['s34']
            D21 = osc_params['D21']
            D31 = osc_params['D31']
            D41 = osc_params['D41']
            return np.array([s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 4, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 's12', 's23', 's13', 'dCP', 'D21', and " + \
                    "'D31' (even if they are None); and 's14', 'd14', 's24', 'd24', 's34', and " + \
                    "'D41'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 5):
        try:
            s12 = osc_params['s12']
            s23 = osc_params['s23']
            s13 = osc_params['s13']
            dCP = osc_params['dCP']
            s14 = osc_params['s14']
            d14 = osc_params['d14']
            s15 = osc_params['s15']
            d15 = osc_params['d15']
            s24 = osc_params['s24']
            d24 = osc_params['d24']
            s25 = osc_params['s25']
            s34 = osc_params['s34']
            s35 = osc_params['s35']
            d35 = osc_params['d35']
            D21 = osc_params['D21']
            D31 = osc_params['D31']
            D41 = osc_params['D41']
            D51 = osc_params['D51']
            return np.array([s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, \
                D21, D31, D41, D51])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 5, the dictionary of oscillation parameters " + \
                    "(osc_params) must contain the keys 's12', 's23', 's13', 'dCP', 'D21', and " + \
                    "'D31' (even if they are None); and 's14', 'd14', 's15', 'd15', 's24', " + \
                    "'d24', 's25', 's34', 's35', 'd35', 'D41', and 'D51'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors > gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
        print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": the number of " + \
            "flavors passed (num_flavors = " + str(num_flavors) + \
            ") exceeds the maximum number for which Magnus has predefined vacuum Hamiltonians " + \
            "(globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
            str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + "). Will use the Hamiltonian provided " + \
            "in h_vac_energy_indep.")
        if (h_vac_energy_indep is None):
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": provided " + \
                "h_vac_energy_indep is None.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors < 1):
        print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": num_flavors must be " + \
            ">= 2.")
        print("Aborting execution...")
        sys.exit(1)


def unpack_nsi_params_from_dict(
    source_func_name: str,
    num_flavors: int,
    nsi_params: Dict,
    h_nsi: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Unpack NSI parameters from the nsi_params dict
    """

    if (num_flavors == 2):
        try:
            eps_aa = nsi_params['eps_aa']
            eps_ab = nsi_params['eps_ab']
            return np.array([eps_aa, eps_ab])
        except KeyError :
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since "+ \
                    "num_flavors == 2, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_aa' and 'eps_ab'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 3):
        try:
            eps_ee = nsi_params['eps_ee']
            eps_em = nsi_params['eps_em']
            eps_et = nsi_params['eps_et']
            eps_mm = nsi_params['eps_mm']
            eps_mt = nsi_params['eps_mt']
            eps_tt = nsi_params['eps_tt']
            return np.array([eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 3, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_ee', 'eps_em', 'eps_et', 'eps_mm'," + \
                    " 'eps_mt', and 'eps_tt'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 4):
        try:
            eps_ee = nsi_params['eps_ee']
            eps_em = nsi_params['eps_em']
            eps_et = nsi_params['eps_et']
            eps_es = nsi_params['eps_es']
            eps_mm = nsi_params['eps_mm']
            eps_mt = nsi_params['eps_mt']
            eps_ms = nsi_params['eps_ms']
            eps_tt = nsi_params['eps_tt']
            eps_ts = nsi_params['eps_ts']
            eps_ss = nsi_params['eps_ss']
            return np.array([eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts,
                eps_ss])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 4, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_ee', 'eps_em', 'eps_et', 'eps_es'," + \
                    " 'eps_mm', 'eps_mt', 'eps_ms', 'eps_tt', 'eps_ts', and 'eps_ss'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 5):
        try:
            eps_ee = nsi_params['eps_ee']
            eps_em = nsi_params['eps_em']
            eps_et = nsi_params['eps_et']
            eps_es1 = nsi_params['eps_es1']
            eps_es2 = nsi_params['eps_es2']
            eps_mm = nsi_params['eps_mm']
            eps_mt = nsi_params['eps_mt']
            eps_ms1 = nsi_params['eps_ms1']
            eps_ms2 = nsi_params['eps_ms2']
            eps_tt = nsi_params['eps_tt']
            eps_ts1 = nsi_params['eps_ts1']
            eps_ts2 = nsi_params['eps_ts2']
            eps_s1s1 = nsi_params['eps_s1s1']
            eps_s1s2 = nsi_params['eps_s1s2']
            eps_s2s2 = nsi_params['eps_s2s2']
            return np.array([eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1,
                eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 5, the dictionary of NSI parameters " + \
                    "(nsi_params) must contain the keys 'eps_ee', 'eps_em', 'eps_et', " + \
                    "'eps_es1', 'eps_es2', 'eps_mm', 'eps_mt', 'eps_ms1', 'eps_ms2', 'eps_tt', " + \
                    "'eps_ts1', 'eps_ts2', 'eps_s1s1', 'eps_s1s2', and 'eps_s2s2'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors > gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
        print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": the number of " + \
            "flavors passed (num_flavors = " + str(num_flavors) + \
            ") exceeds the maximum number for which Magnus has predefined vacuum Hamiltonians " + \
            "(globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
            str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + "). Will use the Hamiltonian provided " + \
            "in h_nsi.")
        if (h_nsi is None):
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": provided " + \
                "h_nsi is None.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors < 1):
        print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": num_flavors must be " + \
            ">= 2.")
        print("Aborting execution...")
        sys.exit(1)


def unpack_liv_params_from_dict(
    source_func_name: str,
    num_flavors: int,
    liv_params: Dict,
    h_liv: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Unpack LIV parameters from the liv_params dict
    """

    if (num_flavors == 2):
        try:
            Lambda = liv_params['Lambda']
            try:
                if (Lambda <= 0.0):
                    raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                        ": Lambda must be positive.")
            except ValueError as error:
                print(error)
                print("Aborting execution...")
                sys.exit(1)            
            sxi = liv_params['sxi']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            n_liv = liv_params['n_liv']
            return np.array([sxi, b1, b2, Lambda, n_liv])
        except KeyError :
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since "+ \
                    "num_flavors == 2, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi', 'b1', 'b2', 'Lambda', 'n_liv'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 3):
        try:
            Lambda = liv_params['Lambda']
            try:
                if (Lambda <= 0.0):
                    raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                        ": Lambda must be positive.")
            except ValueError as error:
                print(error)
                print("Aborting execution...")
                sys.exit(1)            
            sxi12 = liv_params['sxi12']
            sxi23 = liv_params['sxi23']
            sxi13 = liv_params['sxi13']
            dxiCP = liv_params['dxiCP']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            b3 = liv_params['b3']
            n_liv = liv_params['n_liv']
            return np.array([sxi12, sxi23, sxi13, dxiCP, b1, b2, b3, Lambda, n_liv])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 3, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi12', 'sxi23', 'sxi13', 'dxiCP'," + \
                    " 'b1', 'b2', 'b3', 'Lambda', and 'n_liv'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 4):
        try:
            Lambda = liv_params['Lambda']
            try:
                if (Lambda <= 0.0):
                    raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                        ": Lambda must be positive.")
            except ValueError as error:
                print(error)
                print("Aborting execution...")
                sys.exit(1)            
            sxi12 = liv_params['sxi12']
            sxi23 = liv_params['sxi23']
            sxi13 = liv_params['sxi13']
            dxi13 = liv_params['dxi13']
            sxi14 = liv_params['sxi14']
            dxi14 = liv_params['dxi14']
            sxi24 = liv_params['sxi24']
            dxi24 = liv_params['dxi24']
            sxi34 = liv_params['sxi34']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            b3 = liv_params['b3']
            b4 = liv_params['b4']
            n_liv = liv_params['n_liv']
            return np.array([sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2,
                b3, b4, Lambda, n_liv])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 4, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi12', 'sxi23', 'sxi13', 'dxi13'," + \
                    " 'sxi14', 'dxi14', 'sxi24', 'dxi24', 'sxi34', 'b1', 'b2', 'b3', 'b4'," + \
                    " 'Lambda', and 'n_liv'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors == 5):
        try:
            Lambda = liv_params['Lambda']
            try:
                if (Lambda <= 0.0):
                    raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                        ": Lambda must be positive.")
            except ValueError as error:
                print(error)
                print("Aborting execution...")
                sys.exit(1)            
            sxi12 = liv_params['sxi12']
            sxi23 = liv_params['sxi23']
            sxi13 = liv_params['sxi13']
            dxi13 = liv_params['dxi13']
            sxi14 = liv_params['sxi14']
            dxi14 = liv_params['dxi14']
            sxi15 = liv_params['sxi15']
            dxi15 = liv_params['dxi15']            
            sxi24 = liv_params['sxi24']
            dxi24 = liv_params['dxi24']
            sxi25 = liv_params['sxi25']
            sxi34 = liv_params['sxi34']
            sxi35 = liv_params['sxi35']
            dxi35 = liv_params['dxi35']
            b1 = liv_params['b1']
            b2 = liv_params['b2']
            b3 = liv_params['b3']
            b4 = liv_params['b4']
            b5 = liv_params['b5']
            n_liv = liv_params['n_liv']
            return np.array([sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, 
                sxi25, sxi34, sxi35, dxi35, b1, b2, b3, b4, b5, Lambda, n_liv])
        except KeyError:
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": since " + \
                    "num_flavors == 5, the dictionary of LIV parameters " + \
                    "(liv_params) must contain the keys 'sxi12', 'sxi23', 'sxi13', 'dxi13'," + \
                    " 'sxi14', 'dxi14', 'sxi15', 'dxi15', 'sxi24', 'dxi24', 'sxi25' 'sxi34', " + \
                    " 'sxi35', 'dxi35', 'b1', 'b2', 'b3', 'b4', 'b5', 'Lambda', and 'n_liv'.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors > gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS):
        print(gd.WARNING_MSG_IN_COLOR + " oscprob." + source_func_name + ": the number of " + \
            "flavors passed (num_flavors = " + str(num_flavors) + \
            ") exceeds the maximum number for which Magnus has predefined vacuum Hamiltonians " + \
            "(globaldefs.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = " + \
            str(gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS) + "). Will use the Hamiltonian provided " + \
            "in h_liv.")
        if (h_liv is None):
            print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": provided " + \
                "h_liv is None.")
            print("Aborting execution...")
            sys.exit(1)
    elif (num_flavors < 1):
        print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": num_flavors must be " + \
            ">= 2.")
        print("Aborting execution...")
        sys.exit(1)


# def chunkify(lst, n):
#     """Yield successive n-sized chunks from lst."""
#     for i in range(0, len(lst), n):
#         yield lst[i:i + n]


#-----------------------------------------------------------------------
# Primordial functions
#-----------------------------------------------------------------------

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
        return np.eye(n)


def compute_evolution_operator_multiple_slabs(
    H_func: Callable,
    t_slabs: Union[list, np.ndarray],
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
    t_slabs = np.asarray(t_slabs)  # Ensure t_slabs is a NumPy array
    n_slabs = t_slabs.shape[0]

    # Pre-allocate the output array.  Assumes all evolution operators are the same size.
    # Get the shape from a sample Hamiltonian evaluation. This is generally faster
    # than calling H_func repeatedly.
    sample_t = t_slabs[0, 0] if n_slabs > 0 else 0 # Handle empty t_slabs case
    sample_H = H_func(sample_t)
    matrix_dim = sample_H.shape[0]
    U_chain = np.empty((n_slabs, matrix_dim, matrix_dim), dtype=complex) 

    for i, t_slab in enumerate(t_slabs):
        if t_slab[1] > t_slab[0]:
            U_chain[i] = magnus.magnus_expansion(
                lambda t: -1j * H_func(t),
                t0=t_slab[0],
                t1=t_slab[1],
                order=magnus_exp_order,
                n_tpts=n_tpts_per_slab,
                **kwargs,
            )
        else:  # t1 == t0
            U_chain[i] = np.eye(matrix_dim, dtype=complex)  # Use pre-allocated identity

    return U_chain


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
    new_recursion_limit: Optional[int]=5000,
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
            new_recursion_limit, verbose, file_log)

    # Set the maximum recursion limit higher if requested (i.e., if new_recursion_limit is not None)
    if new_recursion_limit is not None:
        old_recursion_limit = sys.getrecursionlimit() 
        if new_recursion_limit > old_recursion_limit:
            sys.setrecursionlimit(new_recursion_limit)  
            if verbose > 0:
                for f in [None, file_log] if save_log else [None]:
                    print("\n" + gd.WARNING_MSG_IN_COLOR + " oscprob.osc_prob: raising recursion"+ \
                            " limit to new_recursion_limit = " + str(new_recursion_limit) + \
                            " (was " + str(old_recursion_limit) + ").\n", file=f)

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
            # U_chain = [compute_evolution_operator(H_func, t_slab, n_tpts_per_slab, magnus_exp_order,
            #     integration_method=integration_method, **kwargs) for t_slab in t_slab_edges]
            U_chain = compute_evolution_operator_multiple_slabs(H_func, t_slab_edges, 
                n_tpts_per_slab, magnus_exp_order, integration_method=integration_method, **kwargs)
        else: # Run n_jobs jobs in parallel
            # batch_size = max(1, n_slabs // (n_jobs * 2))  # Adjust the factor as needed.
            U_chain = Parallel(n_jobs=n_jobs, batch_size='auto')(  
                delayed(compute_evolution_operator)(
                    H_func, t_slab, n_tpts_per_slab, magnus_exp_order, 
                    integration_method=integration_method, **kwargs
                )
                for t_slab in t_slab_edges
            )

            # if n_slabs > 1:
            #     chunk_size = 4  # Experiment with this value
            #     U_chain = []
            #     t_slab_chunks = list(chunkify(t_slab_edges, chunk_size))
            #     U_chain_chunk = Parallel(n_jobs=n_jobs, batch_size=1)(
            #         delayed(compute_evolution_operator_multiple_slabs)(
            #             H_func, t_slabs[0], n_tpts_per_slab, magnus_exp_order,
            #             integration_method=integration_method, **kwargs
            #         ) for t_slabs in t_slab_chunks
            #     )
            #     U_chain.extend(U_chain_chunk)  # Combine results from chunks
            #     U_chain = np.array(U_chain) # Convert back to a NumPy array if needed
            # else:
            #     U_chain = compute_evolution_operator_multiple_slabs(H_func, t_slab_edges, 
            #         n_tpts_per_slab, magnus_exp_order, integration_method=integration_method, 
            #         **kwargs)


            # chunk_size = 4  # Experiment with this value
            # U_chain = []
            # for chunk in chunkify(t_slab_edges, chunk_size):
            #     U_chain_chunk = Parallel(n_jobs=n_jobs)(
            #         delayed(compute_evolution_operator)(
            #             H_func, t_slab, n_tpts_per_slab, magnus_exp_order,
            #             integration_method=integration_method, **kwargs
            #         ) for t_slab in chunk
            #     )
            #     U_chain.extend(U_chain_chunk)  # Combine results from chunks
            # U_chain = np.array(U_chain) # Convert back to a NumPy array if needed

        # Now compute the time-ordered product of all evolution operators across all slabs
        Utot = np.linalg.multi_dot(U_chain) if n_slabs > 1 else U_chain[0]

        # Using Utot, compute all the survival and transition probabilities in a probability matrix
        # P = (np.abs(Utot)**2).T and return that matrix.
        P = np.transpose(np.power(np.abs(Utot), 2))

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
                                "= " + str(magnus_exp_order) + "): rtol = " + str(rtol) + \
                                ", atol = " + str(atol) + ".\n", file=f)
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
    # verbose = 1 if verbose > 0 else verbose

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


def osc_prob_energy_baseline(
    H_func: Union[Callable, np.ndarray],
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    L0: Optional[Union[int, float]]=0.0,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    H_func_is_function_only_of_energy: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[int, float, np.ndarray]:
    r"""Compute and return oscillation probabilities for given arrays of
    neutrino energy and baseline, and an arbitrary Hamiltonian.

    Serves as primordial directly for osc_prob_vacuum,
    osc_prob_matter_std_potential, etc.
    """

    try:
        if (isinstance(H_func, Callable) and (len(signature(H_func).parameters) > 2)):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_energy_baseline:"+\
                " H_func can be energy- and position-dependent, only energy-dependent, or only" + \
                " position-dependent. H_func cannot depend on more than two parameters. To vary" + \
                " the third parameter, call osc_prob_energy_baseline within a loop where it is" + \
                " varied.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # Turn into into float
    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L

    # Flag return_float remembers if energy and L were both floats.  If True, 
    # osc_prob_energy_baseline returns a float, too.
    return_float = isinstance(energy, float) and isinstance(L, float)

    # If there is a single value of energy, make an array out of it.  Same for L.  This will allow
    # us to zip them later.
    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)  
    L = np.array([L]) if isinstance(L, float) else np.array(L) 

    # Either energy and L are both lists (or NumPy arrays) of the same length; or one is a float and
    # the other is a list (or NumPy array).  Any other possibility will generate an exception.  This
    # exception may be raised earlier in routines that call osc_prob_energy_baseline if they are 
    # called wih validate_input == True, but we check below in case it osc_prob_energy_baseline was
    # set to False.
    try:
        if not ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or \
            (len(energy) > 1 and len(L) == 1)):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_energy_baseline: energy and L must be both " + \
                "int or float; or, if lists (or NumPy arrays), they must have the same length;" + \
                " or, if one is a float or single-entry list, the other must be a list with " + \
                "multiple entries.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If energy is a single value, then transform it into an array containing the value energy 
    # repeated a number of times equal to the length of the L, and vice versa, in order to zip them.
    energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
    L = np.full(len(energy), L[0]) if (len(L) == 1) else L

    # The function osc_prob checks whether the Hamiltonian, H_func, is a one-dimensional function or
    # not (if H_func is a function, osc_prob modifies internal run parameters for speed-up).  
    # However, in a physical setting, H_func might still be a one-dimensional function, not of 
    # position, but of energy. Below we deal with these cases before calling osc_prob.

    # The call to __getitem__ below is a way to return a float if both energy and L were floats.
    
    # In the zip: xy[0]: energy, xy[1]: baseline
    
    try:
        if not isinstance(H_func, Callable): # H_func is position- and energy-independent
            if ((nu_i is not None) and (nu_f is not None)): # Select one oscillation channel
                return np.array([osc_prob(H_func, L0, xy[1],
                    t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                    integration_method=integration_method, rtol=rtol, atol=atol,
                    growth_factor_n_slabs=growth_factor_n_slabs,
                    growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                    max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
                    min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
                    iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                    min_magnus_exp_order=min_magnus_exp_order, 
                    max_magnus_exp_order=max_magnus_exp_order,
                    validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                    file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                    new_recursion_limit=new_recursion_limit, verbose=verbose)[nu_i][nu_f]
                for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
            else: # Select the full probability matrix
                return np.array([osc_prob(H_func, L0, xy[1], 
                    t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                    integration_method=integration_method, rtol=rtol, atol=atol,
                    growth_factor_n_slabs=growth_factor_n_slabs,
                    growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                    max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
                    min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
                    iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                    min_magnus_exp_order=min_magnus_exp_order, 
                    max_magnus_exp_order=max_magnus_exp_order,
                    validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                    file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                    new_recursion_limit=new_recursion_limit, verbose=verbose)
                for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
        else: # H_func is a function of one or more parameters 
            if (len(signature(H_func).parameters) == 2): # H_func is a function of two parameters
                # It is assumed that the first parameter is energy and the second is position.
                if ((nu_i is not None) and (nu_f is not None)): # Select one oscillation channel
                    return np.array([osc_prob(lambda l: H_func(xy[0], l), L0, xy[1], 
                        t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                        integration_method=integration_method, rtol=rtol, atol=atol,
                        growth_factor_n_slabs=growth_factor_n_slabs,
                        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
                        min_n_tpts_per_slab=min_n_tpts_per_slab,
                        max_n_tpts_per_slab=max_n_tpts_per_slab,
                        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                        min_magnus_exp_order=min_magnus_exp_order, 
                        max_magnus_exp_order=max_magnus_exp_order,
                        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                        new_recursion_limit=new_recursion_limit, verbose=verbose)[nu_i][nu_f]
                    for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
                else: # Select the full probability matrix
                    return np.array([osc_prob(lambda l: H_func(xy[0], l), L0, xy[1], 
                        t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                        integration_method=integration_method, rtol=rtol, atol=atol,
                        growth_factor_n_slabs=growth_factor_n_slabs,
                        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
                        min_n_tpts_per_slab=min_n_tpts_per_slab,
                        max_n_tpts_per_slab=max_n_tpts_per_slab,
                        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                        min_magnus_exp_order=min_magnus_exp_order, 
                        max_magnus_exp_order=max_magnus_exp_order,
                        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                        new_recursion_limit=new_recursion_limit, verbose=verbose)
                    for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
            elif (len(signature(H_func).parameters) == 1): # H_func is a function of one parameter
                if H_func_is_function_only_of_energy: # H_func is a function only of energy
                    if ((nu_i is not None) and (nu_f is not None)): # Select one oscillation channel
                        return np.array([osc_prob(H_func(xy[0]), L0, xy[1], 
                            t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                            integration_method=integration_method, rtol=rtol, atol=atol,
                            growth_factor_n_slabs=growth_factor_n_slabs,
                            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, 
                            max_n_slabs=max_n_slabs,
                            min_n_tpts_per_slab=min_n_tpts_per_slab,
                            max_n_tpts_per_slab=max_n_tpts_per_slab,
                            iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                            min_magnus_exp_order=min_magnus_exp_order, 
                            max_magnus_exp_order=max_magnus_exp_order,
                            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                            new_recursion_limit=new_recursion_limit, verbose=verbose)[nu_i][nu_f]
                        for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
                    else: # Select the full probability matrix
                        return np.array([osc_prob(H_func(xy[0]), L0, xy[1], 
                            t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                            integration_method=integration_method, rtol=rtol, atol=atol,
                            growth_factor_n_slabs=growth_factor_n_slabs,
                            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, 
                            max_n_slabs=max_n_slabs,
                            min_n_tpts_per_slab=min_n_tpts_per_slab,
                            max_n_tpts_per_slab=max_n_tpts_per_slab,
                            iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                            min_magnus_exp_order=min_magnus_exp_order, 
                            max_magnus_exp_order=max_magnus_exp_order,
                            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                            new_recursion_limit=new_recursion_limit, verbose=verbose)
                        for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
                else: # H_func is a function only of position
                    if ((nu_i is not None) and (nu_f is not None)): # Select one oscillation channel
                        return np.array([osc_prob(lambda l: H_func(l), L0, xy[1], 
                            t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                            integration_method=integration_method, rtol=rtol, atol=atol,
                            growth_factor_n_slabs=growth_factor_n_slabs,
                            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, 
                            max_n_slabs=max_n_slabs,
                            min_n_tpts_per_slab=min_n_tpts_per_slab,
                            max_n_tpts_per_slab=max_n_tpts_per_slab,
                            iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                            min_magnus_exp_order=min_magnus_exp_order, 
                            max_magnus_exp_order=max_magnus_exp_order,
                            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                            new_recursion_limit=new_recursion_limit, verbose=verbose)[nu_i][nu_f]
                        for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
                    else: # Select the full probability matrix
                        return np.array([osc_prob(lambda l: H_func(l), L0, xy[1], 
                            t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                            integration_method=integration_method, rtol=rtol, atol=atol,
                            growth_factor_n_slabs=growth_factor_n_slabs,
                            growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                            max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, 
                            max_n_slabs=max_n_slabs,
                            min_n_tpts_per_slab=min_n_tpts_per_slab,
                            max_n_tpts_per_slab=max_n_tpts_per_slab,
                            iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                            min_magnus_exp_order=min_magnus_exp_order, 
                            max_magnus_exp_order=max_magnus_exp_order,
                            validate_input=validate_input, save_log=save_log, filename_log=filename_log,
                            file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
                            new_recursion_limit=new_recursion_limit, verbose=verbose)
                        for xy in zip(energy, L)]).__getitem__(0 if return_float else slice(None))
    except RecursionError:
        print(gd.ERROR_MSG_IN_COLOR + " oscprob.osc_prob_energy_baseline: error improvement too" + \
            " slow and maximum recursion reached. Consider calling osc_prob_energy_baseline " + \
            "with a higher value of new_recursion_limit (current value of maximum recursion is " + \
            "sys.getrecursionlimit = " + str(sys.getrecursionlimit()) + ". [Failing that, " + \
            "consider running with a higher value of min_n_slabs (currently, min_n_slabs = " + \
            str(min_n_slabs) + ") or a lower requested tolerance (currently, rtol = " + \
            str(rtol) + ", atol = " + str(atol) + ").  Also, consider running using multiple " + \
            "cores by inreasing n_jobs (currently, n_jobs = " + str(n_jobs) + ").]")
        print("Aborting execution...")
        sys.exit(1)


#-----------------------------------------------------------------------
# General functions for vacuum, standard matter, NSI, LIV
#-----------------------------------------------------------------------

def osc_prob_vacuum(
    num_flavors: int,
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for 
    oscillations in vacuum
    """

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name, 
        num_flavors, osc_params, h_vac_energy_indep)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list

    if validate_input:
        if validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=0.0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            validate_energy_and_L=True, validate_flavor_indices=True, validate_osc_params=True, 
            validate_initial_position=False, validate_density=False) == 1:
            sys.exit(1)

    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians5nu.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar) 

    def htot(enu: Union[int, float]) -> np.ndarray:
        return (1/enu)*h_vac_energy_indep

    htot_is_function_only_of_energy = True

    # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
    return osc_prob_energy_baseline(htot, energy, L, 0.0, nu_i, nu_f, 
        htot_is_function_only_of_energy, new_recursion_limit=None)


def osc_prob_matter_std_potential(
    num_flavors: int,
    rho_func: Union[Callable, int, float],
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    L0: Optional[Union[int, float]]=0.0,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for 
    standard oscillations in matter, i.e., the matter potential is only
    due to the coherent forward scattering of nu_e on electrons.
    """

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, osc_params, h_vac_energy_indep)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list

    if validate_input:
        if validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=L0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            rho_func=rho_func, ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction, validate_energy_and_L=True, 
            validate_flavor_indices=True, validate_osc_params=True, validate_initial_position=True,
            validate_density=True) == 1:
            sys.exit(1)

    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians5nu.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)

    # Build the coherent forward potential function, VCC_func, from the density function, rho_func.
    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3].
    VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
        electron_fraction, nubar, density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons) # [eV] 
    
    s = -1.0 if nubar else 1.0

    # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
    if isinstance(VCC_func, Callable):
        # VCC_func is a function of position, so the Hamiltonian is, too
        def htot(enu: Union[int, float], l: Union[int, float]) -> np.ndarray:
            h_matt = np.zeros((num_flavors, num_flavors))
            h_matt[0][0] = s*VCC_func(l)
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = False
    else:
        # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is passed to
        # osc_prob below, osc_prob will detect that VCC_func is constant and set parameters 
        # internally for speed-up.
        h_matt = np.zeros((num_flavors, num_flavors))
        h_matt[0][0] = s*VCC_func
        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = True

    # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
    return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f,
        htot_is_function_only_of_energy, t_slab_edges=t_slab_edges, 
        magnus_exp_order=magnus_exp_order, n_jobs=n_jobs, integration_method=integration_method,
        rtol=rtol, atol=atol, growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order, max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit, verbose=verbose, **kwargs)


def osc_prob_matter_nsi(
    num_flavors: int,
    rho_func: Union[Callable, int, float],
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    nsi_params: Dict,
    L0: Optional[Union[int, float]]=0.0,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    h_nsi: Union[list, np.ndarray]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for 
    standard oscillations in matter, i.e., the matter potential is only
    due to the coherent forward scattering of nu_e on electrons.
    """

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, osc_params, h_vac_energy_indep)
    nsi_params_list = unpack_nsi_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, nsi_params, h_nsi)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
        eps_aa, eps_ab = nsi_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
        eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = nsi_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
        eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss = \
            nsi_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list
        eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, \
            eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2 = nsi_params_list

    if validate_input:
        if validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=L0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            rho_func=rho_func, ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction, validate_energy_and_L=True, 
            validate_flavor_indices=True, validate_osc_params=True, validate_initial_position=True,
            validate_density=True) == 1:
            sys.exit(1)

    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians5nu.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)

    s = -1.0 if nubar else 1.0

    # Compute the standard + NSI matter Hamiltonian *without* the multiplicative prefactor of VCC.
    # To do this we call the functions hamiltonians_Xnu_nsi(VCC, ...) with VCC = 1.0.  We add the
    # standard matter contribution to the NSI matter contribution by adding 1.0 to the eps_ee entry.
    if num_flavors == 2:
        h_matt = s*hamiltonians2nu.hamiltonian_2nu_nsi(1.0, 1.0+eps_aa, eps_ab) # VCC = 1.0
    elif num_flavors == 3:
        h_matt = s*hamiltonians3nu.hamiltonian_3nu_nsi(1.0, 1.0+eps_ee, eps_em, eps_et, eps_mm, 
            eps_mt, eps_tt)
    elif num_flavors == 4:
        h_matt = s*hamiltonians4nu.hamiltonian_4nu_nsi(1.0, 1.0+eps_ee, eps_em, eps_et, eps_es, 
            eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_tt)
    elif num_flavors == 5:
        h_matt = s*hamiltonians5nu.hamiltonian_5nu_nsi(1.0, 1.0+eps_ee, eps_em, eps_et, eps_es1, 
            eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2,
            eps_s2s2)

    # Build the coherent forward potential function, VCC_func, from the density function, rho_func.
    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3].
    VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
        electron_fraction, nubar, density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons) # [eV] 
    
    # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
    if isinstance(VCC_func, Callable):
        # VCC_func is a function of position, so the Hamiltonian is, too
        def htot(enu: Union[int, float], l: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+VCC_func(l)*h_matt
        htot_is_function_only_of_energy = False
    else:
        # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is passed to
        # osc_prob below, osc_prob will detect that VCC_func is constant and set parameters 
        # internally for speed-up.
        h_matt = VCC_func*h_matt
        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = True

    # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
    return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f,
        htot_is_function_only_of_energy, t_slab_edges=t_slab_edges, 
        magnus_exp_order=magnus_exp_order, n_jobs=n_jobs, integration_method=integration_method,
        rtol=rtol, atol=atol, growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order, max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit, verbose=verbose, **kwargs)


def osc_prob_liv(
    num_flavors: int,
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    osc_params: Dict,
    liv_params: Dict,
    rho_func: Optional[Union[Callable, int, float]]=0.0,
    L0: Optional[Union[int, float]]=0.0,
    h_vac_energy_indep: Union[list, np.ndarray]=None,
    h_liv_energy_indep: Union[list, np.ndarray]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Computes and returns neutrino oscillation probabilities for 
    oscillations under (one form of) Lorentz-invariance violation, in 
    vacuum or in matter.
    """

    # Unpack oscillation parameters from the osc_params dict, check if all values are available
    # The function name is sys._getframe().f_code.co_name
    osc_params_list = unpack_oscillation_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, osc_params, h_vac_energy_indep)
    liv_params_list = unpack_liv_params_from_dict(sys._getframe().f_code.co_name,
        num_flavors, liv_params, h_liv_energy_indep)
    if num_flavors == 2:
        sth, Dm2 = osc_params_list
        sxi, b1, b2, Lambda, n_liv = liv_params_list
    elif num_flavors == 3:
        s12, s23, s13, dCP, D21, D31 = osc_params_list
        sxi12, sxi23, sxi13, dxiCP, b1, b2, b3, Lambda, n_liv = liv_params_list
    elif num_flavors == 4:
        s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41 = osc_params_list
        sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, \
            n_liv = liv_params_list
    elif num_flavors == 5:
        s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51 = \
            osc_params_list
        sxi12, sxi23, sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, \
            dxi35, b1, b2, b3, b4, b5, Lambda, n_liv = liv_params_list

    if validate_input:
        if validate_input_battery(sys._getframe().f_code.co_name, energy=energy, L=L, L0=L0,
            num_flavors=num_flavors, nu_i=nu_i, nu_f=nu_f, osc_params=osc_params_list, 
            rho_func=rho_func, ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction, validate_energy_and_L=True, 
            validate_flavor_indices=True, validate_osc_params=True, validate_initial_position=True,
            validate_density=True) == 1:
            sys.exit(1)
    
    # If any of the standard oscillation parameters has not been given a value, assign to it the 
    # value from the specified parameter set with name default_osc_params_set_name.  Only the values
    # of the parameters passed as None are assigned from the predefined set; others are not 
    # modified.
    if num_flavors > 2:
        s12, s23, s13, dCP, D21, D31 = values_to_unspecified_osc_params(s12, s23, s13, dCP, D21, 
            D31, default_osc_params_set_name, verbose)

    # Compute the energy-independent part of the vacuum Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_vac_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_vac_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians5nu.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)
    
    # Compute the energy-independent part of the LIV Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_liv_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_liv_energy_indep = hamiltonians2nu.hamiltonian_2nu_liv_energy_independent(sxi, b1, b2, 
            Lambda, n_liv)
    elif num_flavors == 3:
        h_liv_energy_indep = hamiltonians3nu.hamiltonian_3nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxiCP, b1, b2, b3, Lambda, n_liv, nubar=nubar)
    elif num_flavors == 4:
        h_liv_energy_indep = hamiltonians4nu.hamiltonian_4nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, n_liv,
            nubar=nubar)
    elif num_flavors == 5:
        h_liv_energy_indep = hamiltonians5nu.hamiltonian_5nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, dxi35, b1,
            b2, b3, b4, b5, Lambda, n_liv, nubar=nubar)
   
    if (rho_func != 0.0): # Matter density is nonzero, include the matter term in the Hamiltonian

        # Compute the standard matter Hamiltonian *without* the multiplicative prefactor of VCC.
        h_matt = np.zeros((num_flavors, num_flavors))
        h_matt[0][0] = -1.0 if nubar else 1.0

        # Build the coherent forward potential function, VCC_func, from the density function, 
        # rho_func. If the provided rho_func is the matter density (e.g., g cm^{-3}), convert 
        # rho_func to a function that returns the electron number density [eV^3].
        VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
            electron_fraction, nubar, density_matter_is_in_g_per_cm3,
            density_is_of_number_of_electrons) # [eV] 
        
        # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
        if isinstance(VCC_func, Callable):
            # VCC_func is a function of position, so the Hamiltonian is, too
            def htot(enu: Union[int, float], l: Union[int, float]) -> np.ndarray:
                return (1/enu)*h_vac_energy_indep + VCC_func(l)*h_matt + \
                    pow(enu,n_liv)*h_liv_energy_indep
            htot_is_function_only_of_energy = False
        else:
            # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is 
            # passed to osc_prob below, osc_prob will detect that VCC_func is constant and set 
            # parameters  internally for speed-up.
            h_matt[0][0] *= VCC_func
            def htot(enu: Union[int, float]) -> np.ndarray:
                return (1/enu)*h_vac_energy_indep + h_matt + pow(enu,n_liv)*h_liv_energy_indep
            htot_is_function_only_of_energy = True

    else: # Matter density is zero; the only terms in the Hamiltonian are vacuum and LIV

        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep + pow(enu,n_liv)*h_liv_energy_indep
        htot_is_function_only_of_energy = True

    # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).
    return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f,
        htot_is_function_only_of_energy, t_slab_edges=t_slab_edges, 
        magnus_exp_order=magnus_exp_order, n_jobs=n_jobs, integration_method=integration_method,
        rtol=rtol, atol=atol, growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order, max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit, verbose=verbose, **kwargs)


#-----------------------------------------------------------------------
# In vacuum
#-----------------------------------------------------------------------

def osc_prob_2nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True,
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    vacuum.

    By default, returns :math:`2 \times 2` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pab,
    Pba],[Pba,Pab]])``.  The matrix is symmetric, i.e., ``Pba == Pab``.

    If a single energy and baseline is given, the function returns a 
    single matrix.  If multiple energies and baselines are given, 
    function returns an NumPy array of matrices.  See examples below.

    If the probability needs to be computed multiple times, it is 
    recommended to pass the array of energies and the array of baselines
    to the function in a single call instead of calling the function
    separately for each combination of energy and baseline. The reason
    is that the function has an overhead that gets diluted when 
    computing when the input energies and baselines are many.

    Unlike :func:`osc_prob_3nu_vacuum` (and also 
    :func:`osc_prob_4nu_vacuum` and :func:`osc_prob_5nu_vacuum`), the
    oscillation parameters `sth` and `Dm2` are not optional, but must be
    passed.  Depending on the values passed, :func:`osc_prob_2nu_vacuum`
    will return probabilities for different two-neutrino systems: 

    - :math:`\nu_e-\nu_\mu` if ``sth`` is :math:`\sin \theta_{12}` and 
      ``Dm2`` is :math:`\Delta m_{21}^2`
    - :math:`\nu_\mu-\nu_\tau` if ``sth`` is :math:`\sin \theta_{23}`
      and ``Dm2`` is :math:`\Delta m_{32}^2`
    - :math:`\nu_e-\nu_\tau` if ``sth`` is :math:`\sin \theta_{13}` 
      and ``Dm2`` is :math:`\Delta m_{31}^2`.

    If the initial and final flavors, ``nu_i`` and ``nu_f``, are 
    specified (by setting them to ``NUE``, ``NUMU``, or ``NUTAU``
    from the :py:mod:`magnus.globaldefs` module), the function returns 
    instead a one-dimensional array of the probabilities computed for
    each value of energy and baseline requested. See examples below.

    Because this is a two-neutrino system, the flavor indices can only 
    be 0 or 1.  To prevent using other values, we convert the indices
    like this:
    
    - If ``nu_i == NUE`` (i.e., 0) and ``nu_f == NUTAU`` (i.e., 2), we 
      set ``nu_f = 1``
    - If ``nu_i == NUTAU`` (i.e., 2) and ``nu_f == NUE`` (i.e., 0), we 
      set ``nu_i = 1``
    - If ``nu_i == NUMU`` (i.e., 1) and ``nu_f == NUTAU`` (i.e., 2), we 
      set ``nu_i = 0`` and ``nu_f = 1``
    - If ``nu_i == NUTAU`` (i.e., 2) and ``nu_f == NUMU`` (i.e., 1), we
      set ``nu_i = 1`` and ``nu_f = 0``

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    Parameters
    ----------
    energy
        Neutrino energy, single value or array.
    L
        Neutrino baseline, single value or array.
    sth
        Sine of the mixing angle :math:`\theta`.
    Dm2
        Mass-squared difference :math:`\Delta m^2`.
    nu_i
        Initial neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    nu_f
        Final neutrino flavor, either ``NUE``, ``NUMU``, or ``NUTAU``
        from the :py:mod:`magnus.globaldefs` module.
    validate_input
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

    Returns
    -------
    Union[float, np.narray]
        Neutrino oscillation probability matrix or probability for a 
        single oscillation channel, for the values of `energy` and `L`.

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
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_vacuum(
        num_flavors=2,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


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
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    vacuum.

    By default, returns :math:`3 \times 3` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pee,
    Pem,Pet],[Pme,Pmm,Pmt],[Pte,Ptm,Ptt]])``.  The matrix is symmetric, 
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

    If ``validate_input`` is set to True, the function validates the 
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
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

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

    return osc_prob_vacuum(
        num_flavors=3,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
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
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in vacuum.

    By default, returns :math:`4 \times 4` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pee,
    Pem,Pet,Pes],[Pme,Pmm,Pmt,Pms],[Pte,Ptm,Ptt,Pts],
    [Pse,Psm,Pst,Pss]])``.  The matrix is symmetric, i.e., 
    ``Pme == Pee``, ``Pte == Pet``, ``Pse == Pes``, ``Ptm == Pmt``,
    ``Psm == Pms``, and ``Pst == Pts``.

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
    specified (by setting them to ``NUE``, ``NUMU``, ``NUTAU``, or 
    ``NUS`` from the :py:mod:`magnus.globaldefs` module), the function
    returns instead a one-dimensional array of the probabilities 
    computed for each value of energy and baseline requested. See 
    examples below.

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

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    Parameters
    ----------
    energy
        Neutrino energy, single value or array.
    L
        Neutrino baseline, single value or array.
    s14
        Sine of the mixing angle :math:`\theta_{14}`.
    s24
        Sine of the mixing angle :math:`\theta_{24}`.
    s34
        Sine of the mixing angle :math:`\theta_{34}`.
    d14
        CP-violation phase, :math:`\delta_{14}`.
    d24
        CP-violation phase, :math:`\delta_{24}`.
    D41
        Mass-squared difference :math:`\Delta m_{41}^2`.
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
        Initial neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        or ``NUS`` from the :py:mod:`magnus.globaldefs` module.
    nu_f
        Final neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        or ``NUS`` from the :py:mod:`magnus.globaldefs` module.
    default_osc_params_set_name
        Name of the predefined set of standard oscillation parameters to
        use when assigning default values to unspecified parameters.
    validate_input
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

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
    >>> s14, s24, s34 = 0.1, 0.2, 0.3
    >>> d14, d24 = np.radians(10.0), np.radians(100.0)
    >>> D41 = 0.1  # [eV^2]
    >>> oscprob.osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41)
 
    Pick one channel only, e.g., :math:`\nu_e \to \nu_s`, by 
    passing an initial flavor, ``nu_i``, and a final flavor ``nu_f``:

    >>> oscprob.osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUS)

    The flavor indices ``NUE``, ``NUMU``, ``NUMU``, and ``NUS`` are 
    defined in  the :py:mod:`magnus.globaldefs` module. For 
    anti-neutrinos, i.e., :math:`\bar{\nu}_e \to \bar{\nu}_\mu`:

    >>> oscprob.osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUMU, nubar=True)

    We can specify values of the oscillation parameters. Unspecified 
    values are set to their defaults (pass nonzero ``verbose`` to see 
    this and other warnings):

    >>> oscprob.osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41, s12=0.0, verbose=1)
 
    If a single energy value and multiple baselines are passed, this
    function returns an array containing the probabilities computed for
    this fixed energy and each value of the baseline:
    
    >>> baselines = gd.UNIT_KM*np.array([1.0, 10.0 100.0])
    >>> oscprob.osc_prob_4nu_vacuum(energy, baselines, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUMU)

    Conversely, if a single baseline and multiple energies are passed,
    this function returns an array containing the probabilities computed
    for this fixed baseline and each value of the energy:

    >>> energies = gd.UNIT_MEV*np.array([1.0, 10.0, 100.0])
    >>> oscprob.osc_prob_4nu_vacuum(energies, baseline, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUMU)

    And, for multiple energies and baselines:

    >>> oscprob.osc_prob_4nu_vacuum(energies, baselines, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUMU)

    .. seealso::
        :func:`osc_prob_2nu_vacuum`
            Two-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_3nu_vacuum`
            Three-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_5nu_vacuum`
            Five-flavor (3+2) oscillation probabilities in vacuum. 
    """

    return osc_prob_vacuum(
        num_flavors=4,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_vacuum(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
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
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in vacuum.

    By default, returns :math:`5 \times 5` probability matrices for all 
    the oscillation channels. Each matrix has shape ``np.ndarray([[Pee,
    Pem,Pet,Pes1,Pes2],[Pme,Pmm,Pmt,Pms1,Pms2],[Pte,Ptm,Ptt,Pts1,Pts2],
    [Ps1e,Ps1m,Ps1t,Ps1s1,Ps1s2],[Ps2e,Ps2m,Ps2t,Ps2s1,Ps2s2]])``.  The 
    matrix is symmetric, i.e., ``Pme == Pee``, ``Pte == Pet``, 
    ``Ps1e == Pes1``, ``Ps2e == Pes2`` ``Ptm == Pmt``,
    ``Ps1m == Pms1``, ``Ps1t == Pts1``, ``Ps2t == Pts2``, and 
    ``Ps2s1 == Ps1s2``.

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
    specified (by setting them to ``NUE``, ``NUMU``, ``NUTAU``, 
    ``NUS1``, or ``NUS2`` from the :py:mod:`magnus.globaldefs` module),
    the function returns instead a one-dimensional array of the 
    probabilities computed for each value of energy and baseline 
    requested. See examples below.

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

    If ``validate_input`` is set to True, the function validates the 
    input arguments before calculating the probability, by calling the
    function :func:`validate_input_battery`.

    Parameters
    ----------
    energy
        Neutrino energy, single value or array.
    L
        Neutrino baseline, single value or array.
    s14
        Sine of the mixing angle :math:`\theta_{14}`.
    s15
        Sine of the mixing angle :math:`\theta_{15}`.
    s24
        Sine of the mixing angle :math:`\theta_{24}`.
    s25
        Sine of the mixing angle :math:`\theta_{25}`.
    s34
        Sine of the mixing angle :math:`\theta_{34}`.
    s35
        Sine of the mixing angle :math:`\theta_{35}`.
    d14
        CP-violation phase, :math:`\delta_{14}`.
    d15
        CP-violation phase, :math:`\delta_{15}`.
    d24
        CP-violation phase, :math:`\delta_{24}`.
    d35
        CP-violation phase, :math:`\delta_{35}`.
    D41
        Mass-squared difference :math:`\Delta m_{41}^2`.
    D51
        Mass-squared difference :math:`\Delta m_{51}^2`.
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
        Initial neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        ``NUS1``, or ``NUS2`` from the :py:mod:`magnus.globaldefs`
        module.
    nu_f
        Final neutrino flavor, either ``NUE``, ``NUMU``, ``NUTAU``,
        ``NUS1``, or ``NUS2`` from the :py:mod:`magnus.globaldefs`
        module.
    default_osc_params_set_name
        Name of the predefined set of standard oscillation parameters to
        use when assigning default values to unspecified parameters.
    validate_input
        True to validate input (default); False not to, which is faster
        but riskier.
    verbose
        0 not to print warnings and errors; 1 to print them; 2 to print
        progress.

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
    >>> s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 1.e-2, 1.e-2, 1.e-3, 1.e-3
    >>> d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    >>> D41, D51 = 0.1, 0.001  # [eV^2]
    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51)
 
    Pick one channel only, e.g., :math:`\nu_e \to \nu_{s_1}`, by 
    passing an initial flavor, ``nu_i``, and a final flavor ``nu_f``:

    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUS1)

    And, for :math:`\nu_e \to \nu_{s_2},

    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUS2)

    and :math:`\nu_{s_1} \to \nu_{s_2},

    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUS1, nu_f=gd.NUS2)

    The flavor indices ``NUE``, ``NUMU``, ``NUMU``, ``NUS1``, and 
    ``NUS2`` are defined in  the :py:mod:`magnus.globaldefs` module. For 
    anti-neutrinos, i.e., :math:`\bar{\nu}_e \to \bar{\nu}_{s_1}`:

    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUS1)
    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUS2)
    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUS1, nu_f=gd.NUS2)

    We can specify values of the oscillation parameters. Unspecified 
    values are set to their defaults (pass nonzero ``verbose`` to see 
    this and other warnings):

    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, s12=0.0, verbose=1)
 
    If a single energy value and multiple baselines are passed, this
    function returns an array containing the probabilities computed for
    this fixed energy and each value of the baseline:
    
    >>> baselines = gd.UNIT_KM*np.array([1.0, 10.0 100.0])
    >>> oscprob.osc_prob_5nu_vacuum(energy, baselines, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUMU)

    Conversely, if a single baseline and multiple energies are passed,
    this function returns an array containing the probabilities computed
    for this fixed baseline and each value of the energy:

    >>> energies = gd.UNIT_MEV*np.array([1.0, 10.0, 100.0])
    >>> oscprob.osc_prob_5nu_vacuum(energies, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUMU)

    And, for multiple energies and baselines:

    >>> oscprob.osc_prob_5nu_vacuum(energies, baselines, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUMU)

    .. seealso::
        :func:`osc_prob_2nu_vacuum`
            Two-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_3nu_vacuum`
            Three-flavor oscillation probabilities in vacuum. 
        :func:`osc_prob_4nu_vacuum`
            Four-flavor (3+1) oscillation probabilities in vacuum. 
    """

    return osc_prob_vacuum(
        num_flavors=5,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


#-----------------------------------------------------------------------
# In matter, standard oscillations, constant density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


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
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile.
    """
    return osc_prob_matter_std_potential(
        num_flavors=3,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
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
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with a constant density profile.
    """
    return osc_prob_matter_std_potential(
        num_flavors=4,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
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
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with a constant density profile.
    """
    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, standard oscillations, exponentially falling density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
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
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation 
    probability in matter with an exponentially falling density profile.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_2nu_matter_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


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
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation 
    probability in matter with an exponentially falling density profile.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_3nu_matter_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_std_potential(
        num_flavors=3,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  

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
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with an exponentially falling density profile.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_4nu_matter_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_std_potential(
        num_flavors=4,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  

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
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with an exponentially falling density profile.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_5nu_matter_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  

    return


#-----------------------------------------------------------------------
# In matter, standard oscillations, in the Earth
#-----------------------------------------------------------------------

def osc_prob_2nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

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
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)
    
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5),
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

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
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=3,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

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
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=4,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_earth(
    energy: Union[int, float, list, np.ndarray], 
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior.

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
    # If the location is given as a string, check if it is one of the predefined named locations in
    # Magnus.  The method sys._getframe().f_code.co_name returns the function name.  If the name
    # is one of the predefined ones, coordinates_of_named_location returns the coordinates as 
    # np.array([lat, lon]).  The latitude and longitude are each returned in day-minute-second 
    # format, (dd, mm, ss)

    source_func_name = sys._getframe().f_code.co_name
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # If the initial and final locations are given (i.e., if they are not None), then the neutrino 
    # travels the chord joining them through the Earth, overriding any given value of costhz given,
    # and using the chord length as the baseline. If only a single location is given, throw an 
    # exception.  If neither of the two locations are given, use the given value of costhz and of 
    # baseline given (could be an array of baselines).
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


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


#-----------------------------------------------------------------------
# In matter, standard oscillations, in the Sun
#-----------------------------------------------------------------------

def osc_prob_2nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    sth: Union[int, float],
    Dm2: Union[int, float],
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    magnus_exp_order: Optional[int]=3, #4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, #1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, #2000, 
    min_n_tpts_per_slab: Optional[int]=100, #10, 
    max_n_tpts_per_slab: Optional[int]=400, #500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
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

    Examples
    --------
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_2nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        sth=sth,
        Dm2=Dm2,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_3nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
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
    magnus_exp_order: Optional[int]=3, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, 
    min_n_tpts_per_slab: Optional[int]=100, 
    max_n_tpts_per_slab: Optional[int]=400, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
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

    Examples
    --------
    """

    return osc_prob_3nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
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
    magnus_exp_order: Optional[int]=3, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, 
    min_n_tpts_per_slab: Optional[int]=100, 
    max_n_tpts_per_slab: Optional[int]=400, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
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

    Examples
    --------
    """

    return osc_prob_4nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s24=s24,
        s34=s34,
        d14=d14,
        d24=d24,
        D41=D41,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_sun(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
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
    magnus_exp_order: Optional[int]=3, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, 
    min_n_tpts_per_slab: Optional[int]=100, 
    max_n_tpts_per_slab: Optional[int]=400, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
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

    Examples
    --------
    """

    return osc_prob_5nu_matter_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s15=s15,
        s24=s24,
        s25=s25,
        s34=s34,
        s35=s35,
        d14=d14,
        d15=d15,
        d24=d24,
        d35=d35,
        D41=D41,
        D51=D51,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


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


#-----------------------------------------------------------------------
# In matter, NSI, constant density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_nsi(
        num_flavors=2,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nsi_params={'eps_aa': eps_aa, 'eps_ab': eps_ab},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).
    """
    return osc_prob_matter_nsi(
        num_flavors=3,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_mm': eps_mm,
            'eps_mt': eps_mt, 'eps_tt': eps_tt},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).
    """
    return osc_prob_matter_nsi(
        num_flavors=4,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es': eps_es, 
            'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms': eps_ms, 'eps_tt': eps_tt,
            'eps_ts': eps_ts, 'eps_ss': eps_ss},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_nsi_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with a constant density profile, including non-standard
    interactions (NSI).
    """
    return osc_prob_matter_nsi(
        num_flavors=5,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es1': eps_es1,
            'eps_es2': eps_es2, 'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms1': eps_ms1, 
            'eps_ms2': eps_ms2, 'eps_tt': eps_tt, 'eps_ts1': eps_ts1, 'eps_ts2': eps_ts2, 
            'eps_s1s1': eps_s1s1, 'eps_s1s2': eps_s1s2, 'eps_s2s2': eps_s2s2},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        new_recursion_limit=None,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, NSI, exponentially falling density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_nsi_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    sth: Union[int, float],
    Dm2: Union[int, float],
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with an exponentially falling density profile, including
    non-standard interactions (NSI).
    """

    try:
        if (rho_central <= 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_2nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
                "non-negative.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_matter_nsi(
        num_flavors=2,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nsi_params={'eps_aa': eps_aa, 'eps_ab': eps_ab},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_nsi_exp_density(
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
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with an exponentially falling density profile, including
    non-standard interactions (NSI).
    """

    try:
        if (rho_central <= 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_3nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
                "non-negative.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_nsi(
        num_flavors=3,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_mm': eps_mm,
            'eps_mt': eps_mt, 'eps_tt': eps_tt},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  

    return


def osc_prob_4nu_matter_nsi_exp_density(
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
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0,
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with an exponentially falling density profile,
    including non-standard interactions (NSI).
    """

    try:
        if (rho_central <= 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_4nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
                "non-negative.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_nsi(
        num_flavors=4,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es': eps_es, 
            'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms': eps_ms, 'eps_tt': eps_tt,
            'eps_ts': eps_ts, 'eps_ss': eps_ss},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  

    return


def osc_prob_5nu_matter_nsi_exp_density(
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
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with an exponentially falling density profile,
    including non-standard interactions (NSI).
    """

    try:
        if (rho_central <= 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_5nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
                "non-negative.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_std_potential(
        num_flavors=5,
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es1': eps_es1,
            'eps_es2': eps_es2, 'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms1': eps_ms1, 
            'eps_ms2': eps_ms2, 'eps_tt': eps_tt, 'eps_ts1': eps_ts1, 'eps_ts2': eps_ts2, 
            'eps_s1s1': eps_s1s1, 'eps_s1s2': eps_s1s2, 'eps_s2s2': eps_s2s2},
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  

    return


#-----------------------------------------------------------------------
# In matter, NSI, in the Sun
#-----------------------------------------------------------------------

def osc_prob_2nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    sth: Union[int, float],
    Dm2: Union[int, float],
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    magnus_exp_order: Optional[int]=3, #4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, #1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, #2000, 
    min_n_tpts_per_slab: Optional[int]=100, #10, 
    max_n_tpts_per_slab: Optional[int]=400, #500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    """Compute and return the two-neutrino oscillation probability 
    for neutrinos inside the Sun, including non-standard interactions
    (NSI).
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    """

    return osc_prob_2nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        sth=sth,
        Dm2=Dm2,
        eps_aa=eps_aa,
        eps_ab=eps_ab,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_3nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    magnus_exp_order: Optional[int]=3, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, 
    min_n_tpts_per_slab: Optional[int]=100, 
    max_n_tpts_per_slab: Optional[int]=400, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    """Compute and return the three-neutrino oscillation probability 
    for neutrinos inside the Sun.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    """

    return osc_prob_3nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        eps_ee=eps_ee,
        eps_em=eps_em,
        eps_et=eps_et,
        eps_mm=eps_mm,
        eps_mt=eps_mt,
        eps_tt=eps_tt,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts: Optional[Union[int, float]]=0.0,
    eps_ss: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    magnus_exp_order: Optional[int]=3, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, 
    min_n_tpts_per_slab: Optional[int]=100, 
    max_n_tpts_per_slab: Optional[int]=400, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    """Compute and return the four-neutrino (3+1) oscillation 
    probability for neutrinos inside the Sun.

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    """

    return osc_prob_4nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s24=s24,
        s34=s34,
        d14=d14,
        d24=d24,
        D41=D41,
        eps_ee=eps_ee,
        eps_em=eps_em,
        eps_et=eps_et,
        eps_es=eps_es,
        eps_mm=eps_mm,
        eps_mt=eps_mt,
        eps_ms=eps_ms,
        eps_tt=eps_tt,
        eps_ts=eps_ts,
        eps_ss=eps_ss,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_sun_nsi(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    eps_ee: Optional[Union[int, float]]=0.0,
    eps_em: Optional[Union[int, float]]=0.0,
    eps_et: Optional[Union[int, float]]=0.0,
    eps_es1: Optional[Union[int, float]]=0.0,
    eps_es2: Optional[Union[int, float]]=0.0,
    eps_mm: Optional[Union[int, float]]=0.0,
    eps_mt: Optional[Union[int, float]]=0.0,
    eps_ms1: Optional[Union[int, float]]=0.0,
    eps_ms2: Optional[Union[int, float]]=0.0,
    eps_tt: Optional[Union[int, float]]=0.0,
    eps_ts1: Optional[Union[int, float]]=0.0,
    eps_ts2: Optional[Union[int, float]]=0.0,
    eps_s1s1: Optional[Union[int, float]]=0.0,
    eps_s1s2: Optional[Union[int, float]]=0.0,
    eps_s2s2: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    magnus_exp_order: Optional[int]=3, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, 
    min_n_tpts_per_slab: Optional[int]=100, 
    max_n_tpts_per_slab: Optional[int]=400, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    """Compute and return the five-neutrino (3+2) oscillation 
    probability for neutrinos inside the Sun.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \\exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\\text{Av}~\\text{cm}^{-3}` and 
    :math:`r_0 = R_\\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.

    Examples
    --------
    """

    return osc_prob_5nu_matter_nsi_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s14=s14,
        s15=s15,
        s24=s24,
        s25=s25,
        s34=s34,
        s35=s35,
        d14=d14,
        d15=d15,
        d24=d24,
        d35=d35,
        D41=D41,
        D51=D51,
        eps_ee=eps_ee,
        eps_em=eps_em,
        eps_et=eps_et,
        eps_es1=eps_es1,
        eps_es2=eps_es2,
        eps_mm=eps_mm,
        eps_mt=eps_mt,
        eps_ms1=eps_ms1,
        eps_ms2=eps_ms2,
        eps_tt=eps_tt,
        eps_ts1=eps_ts1,
        eps_ts2=eps_ts2,
        eps_s1s1=eps_s1s1,
        eps_s1s2=eps_s1s2,
        eps_s2s2=eps_s2s2,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


#-----------------------------------------------------------------------
# In vacuum, LIV
#-----------------------------------------------------------------------

def osc_prob_2nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_liv(
        num_flavors=2,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.
    """

    return osc_prob_liv(
        num_flavors=3,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.
    """

    return osc_prob_liv(
        num_flavors=4,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_vacuum_liv(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    vacuum under (one form of) Lorentz-invariance violation.
    """
    
    return osc_prob_liv(
        num_flavors=5,
        rho_func=0.0,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, constant density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Union[int, float], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.
    """
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_liv(
        num_flavors=2,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Optional[Union[int, float]], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.
    """

    return osc_prob_liv(
        num_flavors=3,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Optional[Union[int, float]],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.
    """

    return osc_prob_liv(
        num_flavors=4,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_liv_constant_density(
    energy: Union[int, float, list, np.ndarray], 
    L: Union[int, float, list, np.ndarray], 
    rho: Optional[Union[int, float]],
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with a constant density profile, under (one form of) 
    Lorentz-invariance violation.
    """
    
    return osc_prob_liv(
        num_flavors=5,
        rho_func=rho,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, exponentially falling density
#-----------------------------------------------------------------------

def osc_prob_2nu_matter_liv_exp_density(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    rho_central: Union[int, float], 
    l_scale: Union[int, float],
    sth: Union[int, float],
    Dm2: Union[int, float],
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_2nu_matter_liv_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    return osc_prob_liv(
        num_flavors=2,
        energy=energy,
        L=L,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_matter_liv_exp_density(
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
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_3nu_matter_liv_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_liv(
        num_flavors=3,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_matter_liv_exp_density(
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
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_4nu_matter_liv_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_liv(
        num_flavors=4,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_matter_liv_exp_density(
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
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
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
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability in
    matter with an exponentially falling density profile, under (one 
    form of) Lorentz-invariance violation.
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_5nu_matter_liv_exp_density: rho_central must be non-negative" + \
                " and l_scale must be positive.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_liv(
        num_flavors=5,
        energy=energy,
        L=L,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        rho_func=lambda r: rho_central*np.exp(-r/l_scale),
        L0=L0,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, in the Sun
#-----------------------------------------------------------------------

def osc_prob_2nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    sth: Union[int, float],
    Dm2: Union[int, float],
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    magnus_exp_order: Optional[int]=3, #4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, #1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, #2000, 
    min_n_tpts_per_slab: Optional[int]=100, #10, 
    max_n_tpts_per_slab: Optional[int]=400, #500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.
    """

    return osc_prob_2nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        sth=sth,
        Dm2=Dm2,
        sxi=sxi,
        b1=b1,
        b2=b2,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=rtol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_3nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxiCP: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    magnus_exp_order: Optional[int]=3, #4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, #1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, #2000, 
    min_n_tpts_per_slab: Optional[int]=100, #10, 
    max_n_tpts_per_slab: Optional[int]=400, #500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.
    """

    return osc_prob_3nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        sxi12=sxi12,
        sxi23=sxi23,
        sxi13=sxi13,
        dxiCP=dxiCP,
        b1=b1,
        b2=b2,
        b3=b3,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=rtol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_4nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    magnus_exp_order: Optional[int]=3, #4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, #1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, #2000, 
    min_n_tpts_per_slab: Optional[int]=100, #10, 
    max_n_tpts_per_slab: Optional[int]=400, #500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.
    """

    return osc_prob_4nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        s14=s14,
        s24=s24,
        s34=s34,
        d14=d14,
        d24=d24,
        D41=D41,
        sxi12=sxi12,
        sxi23=sxi23,
        sxi13=sxi13,
        dxi13=dxi13,
        sxi14=sxi14,
        dxi14=dxi14,
        sxi24=sxi24,
        dxi24=dxi24,
        sxi34=sxi34,
        b1=b1,
        b2=b2,
        b3=b3,
        b4=b4,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=rtol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_5nu_sun_liv(
    energy: Union[float, list, np.ndarray], 
    L: Union[float, list, np.ndarray],
    L0: Union[int, float], 
    s12: Optional[Union[int, float]]=None, 
    s23: Optional[Union[int, float]]=None, 
    s13: Optional[Union[int, float]]=None, 
    dCP: Optional[Union[int, float]]=None, 
    D21: Optional[Union[int, float]]=None, 
    D31: Optional[Union[int, float]]=None, 
    s14: Optional[Union[int, float]]=0.0,
    s15: Optional[Union[int, float]]=0.0,
    s24: Optional[Union[int, float]]=0.0,
    s25: Optional[Union[int, float]]=0.0,
    s34: Optional[Union[int, float]]=0.0,
    s35: Optional[Union[int, float]]=0.0,
    d14: Optional[Union[int, float]]=0.0,
    d15: Optional[Union[int, float]]=0.0,
    d24: Optional[Union[int, float]]=0.0,
    d35: Optional[Union[int, float]]=0.0,
    D41: Optional[Union[int, float]]=0.0, 
    D51: Optional[Union[int, float]]=0.0, 
    sxi12: Optional[Union[int, float]]=0.0,
    sxi23: Optional[Union[int, float]]=0.0,
    sxi13: Optional[Union[int, float]]=0.0,
    dxi13: Optional[Union[int, float]]=0.0,
    sxi14: Optional[Union[int, float]]=0.0,
    dxi14: Optional[Union[int, float]]=0.0,
    sxi15: Optional[Union[int, float]]=0.0,
    dxi15: Optional[Union[int, float]]=0.0,
    sxi24: Optional[Union[int, float]]=0.0,
    dxi24: Optional[Union[int, float]]=0.0,
    sxi25: Optional[Union[int, float]]=0.0,
    sxi34: Optional[Union[int, float]]=0.0,
    sxi35: Optional[Union[int, float]]=0.0,
    dxi35: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    b3: Optional[Union[int, float]]=0.0,
    b4: Optional[Union[int, float]]=0.0,
    b5: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0, 
    electron_fraction: Optional[Union[int, float]]=0.5, 
    nubar: Optional[bool]=False, 
    nu_i: Optional[int]=None, 
    nu_f: Optional[int]=None,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
    magnus_exp_order: Optional[int]=3, #4, 
    n_jobs: Optional[int]=1, 
    integration_method: Optional[str]='trapezoid', 
    rtol: Optional[Union[int, float]]=1.e-2, 
    atol: Optional[Union[int, float]]=1.e-2, 
    growth_factor_n_slabs: Optional[Union[int, float]]=2.0, #1.5, 
    growth_factor_n_tpts_per_slab: Optional[Union[int, float]]=2.0, 
    max_num_loops: Optional[int]=50, 
    min_n_slabs: Optional[int]=100, 
    max_n_slabs: Optional[int]=400, #2000, 
    min_n_tpts_per_slab: Optional[int]=100, #10, 
    max_n_tpts_per_slab: Optional[int]=400, #500, 
    iterate_over_magnus_exp_order: Optional[bool]=False,
    min_magnus_exp_order: Optional[int]=1,
    max_magnus_exp_order: Optional[int]=gd.MAGNUS_EXP_ORDER_MAX,
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    new_recursion_limit: Optional[int]=5000,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    for neutrinos inside the Sun, under (one form of) Lorentz-invariance
    violation.
    
    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`, 
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and 
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics 
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung 
    Wook Kim.
    """

    return osc_prob_5nu_matter_liv_exp_density(
        energy=energy,
        L=L,
        L0=L0,
        rho_central=gd.NUM_DENSITY_E_SUN_CENTRAL,
        l_scale=gd.L_SCALE_SUN,
        s12=s12,
        s23=s23,
        s13=s13,
        dCP=dCP,
        D21=D21,
        D31=D31,
        s14=s14,
        s15=s15,
        s24=s24,
        s25=s25,
        s34=s34,
        s35=s35,
        d14=d14,
        d15=d15,
        d24=d24,
        d35=d35,
        D41=D41,
        D51=D51,
        sxi12=sxi12,
        sxi23=sxi23,
        sxi13=sxi13,
        dxi13=dxi13,
        sxi14=sxi14,
        dxi14=dxi14,
        sxi15=sxi15,
        dxi15=dxi15,
        sxi24=sxi24,
        dxi24=dxi24,
        sxi25=sxi25,
        sxi34=sxi34,
        sxi35=sxi35,
        dxi35=dxi35,
        b1=b1,
        b2=b2,
        b3=b3,
        b4=b4,
        b5=b5,
        Lambda=Lambda,
        n_liv=n_liv,
        ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
        electron_fraction=electron_fraction,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons=density_is_of_number_of_electrons,
        magnus_exp_order=magnus_exp_order,
        n_jobs=n_jobs,
        integration_method=integration_method,
        rtol=rtol,
        atol=rtol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops,
        min_n_slabs=min_n_slabs,
        max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab,
        max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order,
        max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit,
        verbose=verbose,
        **kwargs
    )



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

    # np.set_printoptions(precision=3)
    # prob = osc_prob_3nu_vacuum(energy, baseline)
    # print(prob)

    # print(gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'])

    # np.set_printoptions(precision=3)
    # print(osc_prob_3nu_vacuum(energy, baseline, s12=0.0, verbose=1))

    # Two-neutrino oscillations in vacuum
    # np.set_printoptions(precision=3)
    # sth = 0.1
    # Dm2 = 0.1 # [eV^2]
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, verbose=1))
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU, verbose=1))
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUTAU, verbose=1))
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, nu_i=gd.NUMU, nu_f=gd.NUTAU, verbose=1))

    # Three-neutrino oscillations in vacuum
    # np.set_printoptions(precision=3)
    # print(osc_prob_3nu_vacuum(energy, [baseline, baseline], verbose=1))

    # Four-neutrino oscillations in vacuum
    # np.set_printoptions(precision=3)
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # print(osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41, verbose=1))

    # # Five-neutrino oscillations in vacuum
    # np.set_printoptions(precision=3)
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # print(osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35,
    #     D41, D51, verbose=1))

    # Two-neutrino oscillations in constant-density matter
    # np.set_printoptions(precision=3)
    # rho = 10.0*gd.UNIT_G_PER_CM3
    # sth = 0.1
    # Dm2 = 0.1 # [eV^2]
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho, sth, Dm2, verbose=1))
    # print(osc_prob_2nu_matter_constant_density(energy, gd.UNIT_KM*np.array([1.0, 10.0]), 
    #     rho, sth, Dm2, verbose=1))
    # print(osc_prob_2nu_matter_constant_density(gd.UNIT_MEV*np.array([1.0, 10.0]), baseline, 
    #     rho, sth, Dm2, verbose=1))
    # print(osc_prob_2nu_matter_constant_density(gd.UNIT_MEV*np.array([1.0, 5.0]), 
    #     gd.UNIT_KM*np.array([1.0, 10.0]), rho, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU, verbose=1))
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho, sth, Dm2, nu_i=gd.NUE, 
    #     nu_f=gd.NUMU, verbose=1))
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho, sth, Dm2, nu_i=gd.NUE, 
    #     nu_f=gd.NUTAU, verbose=1))
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho, sth, Dm2, nu_i=gd.NUMU, 
    #     nu_f=gd.NUTAU, verbose=1))

    # Three-neutrino oscillations in constant-density matter
    # np.set_printoptions(precision=3)
    # rho = 10.0*gd.UNIT_G_PER_CM3
    # print(osc_prob_3nu_matter_constant_density(energy, baseline, rho, verbose=1))

    # # Four-neutrino oscillations in constant-density matter
    # np.set_printoptions(precision=3)
    # rho = 10.0*gd.UNIT_G_PER_CM3
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # print(osc_prob_4nu_matter_constant_density(energy, baseline, rho, s14, s24, s34, d14, d24, D41,
    #     verbose=1))

    # # Five-neutrino oscillations in constant-density matter
    # np.set_printoptions(precision=3)
    # rho = 10.0*gd.UNIT_G_PER_CM3
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # print(osc_prob_5nu_matter_constant_density(energy, baseline, rho, s14, s15, s24, s25, s34, s35,
    #     d14, d15, d24, d35, D41, D51, verbose=2))
    # print(osc_prob_5nu_matter_constant_density(energy, baseline, rho, s14, s15, s24, s25, s34, s35,
    #     d14, d15, d24, d35, D41, D51, nubar=True, verbose=1))

    # Two-neutrino oscillations in exponentially falling matter density profile
    # np.set_printoptions(precision=3)
    # sth = 0.5
    # Dm2 = 1.e-3 # [eV^2]
    # rho_central = 10.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # print(osc_prob_2nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, sth, Dm2,
    #     verbose=2))
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho_central, sth, Dm2, nu_i=gd.NUE, 
    #     nu_f=gd.NUE, verbose=0))
    # print(osc_prob_2nu_matter_exp_density(energy, baseline, 0.0, rho_central, 1.0*gd.UNIT_KM,
    #     sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_2nu_matter_exp_density(energy, baseline, 10.0*gd.UNIT_KM, rho_central, 
    #     1.0*gd.UNIT_KM, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_2nu_matter_exp_density(energy, baseline, 0.0, rho_central, 100.0*gd.UNIT_KM,
    #     sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))

    # # Three-neutrino oscillations in exponentially falling matter density profile
    # np.set_printoptions(precision=3)
    # rho_central = 10.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # print(osc_prob_3nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, verbose=2))
    # print(osc_prob_3nu_vacuum(energy, baseline, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_3nu_matter_constant_density(energy, baseline, rho_central, nu_i=gd.NUE, 
    #     nu_f=gd.NUE, verbose=0))
    # print(osc_prob_3nu_matter_exp_density(energy, baseline, 0.0, rho_central, 
    #     l_scale=1.0*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_3nu_matter_exp_density(energy, baseline, 10.0*gd.UNIT_KM, rho_central, 
    #     l_scale=1.0*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_3nu_matter_exp_density(energy, baseline, 0.0, rho_central, 
    #     l_scale=100.0*gd.UNIT_KM, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))

    # # Four-neutrino oscillations in exponentially falling matter density profile
    # np.set_printoptions(precision=3)
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # rho_central = 10.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # # print(osc_prob_4nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, 
    # #     s14, s24, s34, d14, d24, D41, verbose=2))
    # # print(osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, 
    # #     nu_f=gd.NUE, verbose=0))
    # print(osc_prob_4nu_matter_constant_density(energy, baseline, rho_central, 
    #     s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_4nu_matter_exp_density(energy, baseline, 0.0, rho_central, 
    #     1.0*gd.UNIT_KM, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_4nu_matter_exp_density(energy, baseline, 10.0*gd.UNIT_KM, rho_central, 
    #     1.0*gd.UNIT_KM, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_4nu_matter_exp_density(energy, baseline, 0.0, rho_central, 
    #     100.0*gd.UNIT_KM, s14, s24, s34, d14, d24, D41, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))

    # # Five-neutrino oscillations in exponentially falling matter density profile
    # np.set_printoptions(precision=3)
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # rho_central = 10.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # # print(osc_prob_5nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, 
    # #     s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, verbose=2))
    # # print(osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, 
    # #     D41, D51, nu_i=gd.NUE, nu_f=gd.NUE, verbose=0))
    # print(osc_prob_5nu_matter_constant_density(energy, baseline, rho_central, 
    #     s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUE, 
    #     verbose=0))
    # print(osc_prob_5nu_matter_exp_density(energy, baseline, 0.0, rho_central, 
    #     1.0*gd.UNIT_KM, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, 
    #     nu_f=gd.NUE, verbose=0))
    # print(osc_prob_5nu_matter_exp_density(energy, baseline, 10.0*gd.UNIT_KM, rho_central, 
    #     1.0*gd.UNIT_KM, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, 
    #     nu_f=gd.NUE, verbose=0))
    # print(osc_prob_5nu_matter_exp_density(energy, baseline, 0.0, rho_central, 
    #     100.0*gd.UNIT_KM, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, 
    #     nu_f=gd.NUE, verbose=0))

    # # Two-neutrino oscillations in the Sun
    # np.set_printoptions(precision=3)
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.0*gd.UNIT_MEV # [eV]
    # print(osc_prob_2nu_sun(energy, baseline, L0, sth, Dm2, n_jobs=10, verbose=1))
    # energy = (10.0+1.e-4)*gd.UNIT_MEV # [eV]
    # print(osc_prob_2nu_sun(energy, baseline, L0, sth, Dm2, n_jobs=10, verbose=2))

    # Three-neutrino oscillations in the Sun
    # np.set_printoptions(precision=3)
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.0*gd.UNIT_MEV # [eV]
    # print(osc_prob_3nu_sun(energy, baseline, L0, n_jobs=10, verbose=1))

    # # Four-neutrino oscillations in exponentially falling matter density profile
    # np.set_printoptions(precision=3)
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.0*gd.UNIT_MEV # [eV]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # print(osc_prob_4nu_sun(energy, baseline, L0, s14, s24, s34, d14, d24, D41, 
    #     n_jobs=10, verbose=1))

    # # Five-neutrino oscillations in exponentially falling matter density profile
    # np.set_printoptions(precision=3)
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.0*gd.UNIT_MEV # [eV]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # print(osc_prob_5nu_sun(energy, baseline, L0, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, 
    #     D41, D51, n_jobs=10, verbose=1))

    # # Two-neutrino oscillations in constant-density matter, NSI
    # np.set_printoptions(precision=3)
    # rho = 100.0*gd.UNIT_G_PER_CM3
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # energy = 1.0*gd.UNIT_MEV
    # eps_aa = 0.1 #gd.EPS_EE
    # eps_ab = 0.2 #gd.EPS_EM
    # print(osc_prob_2nu_matter_nsi_constant_density(energy, baseline, rho, sth, Dm2, eps_aa, eps_ab,
    #     verbose=0))
    # print(osc_prob_2nu_matter_nsi_constant_density(energy, baseline, rho, sth, Dm2, 0.0, 0.0,
    #     verbose=0))
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho, sth, Dm2, verbose=0))
    # # print(osc_prob_2nu_matter_nsi_constant_density(energy, gd.UNIT_KM*np.array([1.0, 10.0]), 
    # #     rho, sth, Dm2, eps_aa, eps_ab, verbose=1))
    # # print(osc_prob_2nu_matter_nsi_constant_density(gd.UNIT_MEV*np.array([1.0, 10.0]), baseline, 
    # #     rho, sth, Dm2, eps_aa, eps_ab, verbose=1))
    # # print(osc_prob_2nu_matter_nsi_constant_density(gd.UNIT_MEV*np.array([1.0, 5.0]), 
    # #     gd.UNIT_KM*np.array([1.0, 10.0]), rho, sth, Dm2, eps_aa, eps_ab, nu_i=gd.NUE, nu_f=gd.NUMU,
    # #     verbose=1))
    # # print(osc_prob_2nu_matter_nsi_constant_density(energy, baseline, rho, sth, Dm2, eps_aa, eps_ab,
    # #     nu_i=gd.NUE, nu_f=gd.NUMU, verbose=1))
    # # print(osc_prob_2nu_matter_nsi_constant_density(energy, baseline, rho, sth, Dm2, eps_aa, eps_ab,
    # #     nu_i=gd.NUE, nu_f=gd.NUTAU, verbose=1))
    # # print(osc_prob_2nu_matter_nsi_constant_density(energy, baseline, rho, sth, Dm2, eps_aa, eps_ab,
    # #     nu_i=gd.NUMU, nu_f=gd.NUTAU, verbose=1))

    # # Three-neutrino oscillations in constant-density matter, NSI
    # np.set_printoptions(precision=3)
    # rho = 100.0*gd.UNIT_G_PER_CM3
    # energy = 1.0*gd.UNIT_MEV
    # eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = 0.1, 0.1, 1.1, 0.1, 0.1, 0.1
    # print(osc_prob_3nu_matter_nsi_constant_density(energy, baseline, rho, eps_ee, eps_em, eps_et, 
    #     eps_mm, eps_mt, eps_tt, verbose=0))
    # print(osc_prob_3nu_matter_constant_density(energy, baseline, rho, verbose=1))

    # # Four-neutrino oscillations in constant-density matter, NSI
    # np.set_printoptions(precision=3)
    # rho = 10.0*gd.UNIT_G_PER_CM3
    # energy = 1.0*gd.UNIT_GEV
    # baseline = 100.*gd.UNIT_KM
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss \
    #     = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # print(osc_prob_4nu_matter_nsi_constant_density(energy, baseline, rho, s14, s24, s34, d14, d24, 
    #     D41, eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss,
    #     verbose=0))
    # print(osc_prob_4nu_matter_constant_density(energy, baseline, rho, s14, s24, s34, d14, d24, D41,
    #     verbose=0))

    # # Five-neutrino oscillations in constant-density matter, NSI
    # np.set_printoptions(precision=3)
    # rho = 10.0*gd.UNIT_G_PER_CM3
    # energy = 1.0*gd.UNIT_GEV
    # baseline = 1000.*gd.UNIT_KM
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, \
    #     eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2 \
    #     = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # print(osc_prob_5nu_matter_nsi_constant_density(energy, baseline, rho, s14, s15, s24, s25, s34, 
    #     s35, d14, d15, d24, d35, D41, D51, eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt,
    #     eps_ms1, eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2, verbose=0))
    # print(osc_prob_5nu_matter_constant_density(energy, baseline, rho, s14, s15, s24, s25, s34, s35,
    #     d14, d15, d24, d35, D41, D51, verbose=0))

    # # Two-neutrino oscillations in exponentially falling matter density profile, NSI
    # np.set_printoptions(precision=3)
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # eps_aa = 0.1 #gd.EPS_EE
    # eps_ab = 0.2 #gd.EPS_EM
    # rho_central = 12.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # print(osc_prob_2nu_matter_nsi_exp_density(energy, baseline, 0.0, rho_central, l_scale, sth, Dm2,
    #     eps_aa, eps_ab, verbose=0))
    # print(osc_prob_2nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, sth, Dm2,
    #     verbose=0))
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, verbose=0))

    # Three-neutrino oscillations in exponentially falling matter density profile, NSI
    # np.set_printoptions(precision=3)
    # rho_central = 1000.0*gd.UNIT_G_PER_CM3
    # l_scale = 100.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = 0.1, 0.1, 1.1, 0.1, 0.1, 0.1
    # print(osc_prob_3nu_matter_nsi_exp_density(energy, baseline, 0.0, rho_central, l_scale, 
    #     eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt, verbose=0))
    # quit()
    # print(osc_prob_3nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, verbose=0))
    # print(osc_prob_3nu_vacuum(energy, baseline, verbose=0))
        # print(osc_prob_3nu_matter_nsi_exp_density(energy, baseline, 0.0, rho_central, l_scale, 
        # eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt, n_jobs=5, verbose=0))
    
    # import time
    # duration = 0
    # n_loops = 30
    # for i in range(n_loops):
    #     start = time.time()
    #     osc_prob_3nu_matter_nsi_exp_density(gd.UNIT_MEV*np.linspace(1,10000,100), 
    #         gd.UNIT_KM*np.linspace(0.1,1000,100), 0.0, rho_central, l_scale, 
    #         eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt, rtol=1.e-3, atol=1.e-3, n_jobs=1, 
    #         verbose=0)
    #     end = time.time()
    #     duration += end-start
    # print("Average time, n_jobs = 1: " + str(duration/n_loops) + " s")

    # duration = 0
    # n_loops = 30
    # for i in range(n_loops):
    #     start = time.time()
    #     osc_prob_3nu_matter_nsi_exp_density(gd.UNIT_MEV*np.linspace(1,10000,100), 
    #         gd.UNIT_KM*np.linspace(0.1,1000,100), 0.0, rho_central, l_scale, 
    #         eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt, rtol=1.e-3, atol=1.e-3, n_jobs=4, 
    #         verbose=0)
    #     end = time.time()
    #     duration += end-start
    # print("Average time, n_jobs = 10: " + str(duration/n_loops) + " s")


    # # Four-neutrino oscillations in exponentially falling matter density profile, NSI
    # np.set_printoptions(precision=3)
    # rho_central = 10.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 50.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss \
    #     = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # print(osc_prob_4nu_matter_nsi_exp_density(energy, baseline, 0.0, rho_central, l_scale, 
    #     s14, s24, s34, d14, d24, D41, eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, 
    #     eps_tt, eps_ts, eps_ss, validate_input=False, verbose=0))
    # print(osc_prob_4nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, s14, s24, 
    #     s34, d14, d24, D41,verbose=0))
    # print(osc_prob_4nu_vacuum(energy, baseline, s14, s24, s34, d14, d24, D41, verbose=0))

    # # Five-neutrino oscillations in exponentially falling matter density profile, NSI
    # np.set_printoptions(precision=3)
    # rho_central = 10.0*gd.UNIT_G_PER_CM3
    # l_scale = 200.0*gd.UNIT_KM
    # baseline = 1000.*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, \
    #     eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2 \
    #     = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # print(osc_prob_5nu_matter_nsi_exp_density(energy, baseline, 0.0, rho_central, l_scale, 
    #     s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, eps_ee, eps_em, eps_et, eps_es1,
    #     eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2,
    #     eps_s2s2, verbose=1))
    # print(osc_prob_5nu_matter_exp_density(energy, baseline, 0.0, rho_central, l_scale, s14, s15, 
    #     s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, verbose=1))
    # print(osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35,
    #     D41, D51, verbose=1))

    # # Two-neutrino oscillations in the Sun, NSI
    # np.set_printoptions(precision=3)
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # eps_aa = 0.1 #gd.EPS_EE
    # eps_ab = 0.2 #gd.EPS_EM
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.0*gd.UNIT_MEV # [eV]
    # print(osc_prob_2nu_sun_nsi(energy, baseline, L0, sth, Dm2, eps_aa, eps_ab,
    #     n_jobs=10, verbose=1))
    # # print(osc_prob_2nu_sun_nsi(energy, baseline, L0, sth, Dm2, eps_aa, eps_ab,
    # #     rtol=1.e-2, atol=1.e-2, magnus_exp_order=4, max_n_slabs=2000, n_jobs=10, verbose=1))

    # # Three-neutrino oscillations in the Sun, NSI
    # np.set_printoptions(precision=3)
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = 0.1, 0.1, 1.1, 0.1, 0.1, 0.1
    # print(osc_prob_3nu_sun_nsi(energy, baseline, L0, eps_ee, eps_em, eps_et, eps_mm, 
    #     eps_mt, eps_tt, n_jobs=10, verbose=1))

    # # Four-neutrino oscillations in the Sun, NSI
    # np.set_printoptions(precision=3)
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss \
    #     = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # print(osc_prob_4nu_sun_nsi(energy, baseline, L0, s14, s24, s34, d14, d24, D41, eps_ee, eps_em,
    #     eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss, n_jobs=10, verbose=1))

    # # Five-neutrino oscillations in the Sun, NSI
    # np.set_printoptions(precision=3)
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # baseline = 1.0*gd.SUN_RADIUS*gd.UNIT_KM # km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, \
    #     eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2 \
    #     = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # # print(osc_prob_5nu_sun_nsi(energy, baseline, L0, s14, s15, s24, s25, s34, s35, d14, d15, d24, 
    # #     d35, D41, D51, eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2,
    # #     eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2, n_jobs=10, verbose=1))
    # print(osc_prob_5nu_sun_nsi(energy, baseline, L0, s14=s14, s15=s15, D41=D41, D51=D51, 
    #     eps_ee=eps_ee, eps_em=eps_em, eps_es1=eps_es1, eps_es2=eps_es2, n_jobs=10, verbose=1))

    # # Two-neutrino oscillations, vacuum, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # sxi = 0.01
    # b1 = 1.e-2
    # b2 = 1.e-3
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 3
    # print(osc_prob_2nu_vacuum(energy, baseline, sth, Dm2, verbose=1))
    # print(osc_prob_2nu_vacuum_liv(energy, baseline, sth, Dm2, sxi, b1, b2, Lambda, n_liv, 
    #     verbose=1))

    # # Three-neutrino oscillations, vacuum, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # sxi12, sxi23, sxi13, dxiCP = 0.001, 0.002, 0.003, np.radians(10.0)
    # b1, b2, b3 = 2.e-10, 1.e-10, 3.e-10
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 0
    # print(osc_prob_3nu_vacuum(energy, baseline, verbose=0))
    # print(osc_prob_3nu_vacuum_liv(energy, baseline, sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, 
    #     dxiCP=dxiCP, b1=b1, b2=b2, b3=b3, Lambda=Lambda, n_liv=n_liv, verbose=0))

    # # Four-neutrino oscillations, vacuum, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 = 0.01, 0.05, 0.06, 0.1, 0.9, 0.02
    # dxi13, dxi14, dxi24 = np.radians([10,20,30])
    # b1, b2, b3, b4 = 1.e-2, 1.e-3, 5.e-3, 4.e-2
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 3
    # print(osc_prob_4nu_vacuum(energy, baseline, s14=s14, s24=s24, s34=s34, d14=d14, d24=d24,
    #     D41=D41, verbose=0))
    # print(osc_prob_4nu_vacuum_liv(energy, baseline, s14=s14, s24=s24, s34=s34, d14=d14, d24=d24,
    #     D41=D41, sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, sxi34=sxi34,
    #     dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, Lambda=Lambda,
    #     n_liv=n_liv, verbose=0))

    # # Five-neutrino oscillations, vacuum, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi15, sxi24, sxi25, sxi34, sxi35 = \
    #     0.01, 0.05, 0.06, 0.1, 0.9, 0.02, 0.04, 0.1, 0.05
    # dxi13, dxi14, dxi15, dxi24, dxi35 = np.radians([10, 20, 30, 40, 50])
    # b1, b2, b3, b4, b5 = 1.e-2, 1.e-3, 5.e-3, 4.e-2, 1.e-1
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 2
    # print(osc_prob_5nu_vacuum(energy, baseline, s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, 
    #     s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51, verbose=0))
    # print(osc_prob_5nu_vacuum_liv(energy, baseline, s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, 
    #     s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51, sxi12=sxi12, sxi23=sxi23,
    #     sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25, sxi34=sxi34, sxi35=sxi35,
    #     dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, b1=b1, b2=b2, b3=b3, b4=b4,
    #     b5=b5, Lambda=Lambda, n_liv=n_liv, verbose=0))

    # # Two-neutrino oscillations in constant-density matter, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 10.*gd.UNIT_GEV # [eV]
    # rho = 10.0*gd.UNIT_G_PER_CM3 # [g cm^{-3}]
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # sxi = 0.5
    # b1 = 1.e-2
    # b2 = 1.e-3
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 3
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho, sth, Dm2, verbose=0))
    # print(osc_prob_2nu_matter_liv_constant_density(energy, baseline, rho, sth, Dm2, sxi=sxi, 
    #     b1=b1, b2=b2, Lambda=Lambda, n_liv=n_liv, verbose=0))

    # # Three-neutrino oscillations in constant-density matter, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_GEV # [eV]
    # rho = 10.0*gd.UNIT_G_PER_CM3 # [g cm^{-3}]
    # sxi12, sxi23, sxi13, dxiCP = 0.01, 0.02, 0.03, np.radians(10.0)
    # b1, b2, b3 = 1.e-2, 1.e-3, 5.e-3
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 3
    # print(osc_prob_3nu_matter_constant_density(energy, baseline, rho, verbose=0))
    # print(osc_prob_3nu_matter_liv_constant_density(energy, baseline, rho, sxi12=sxi12, sxi23=sxi23,
    #     sxi13=sxi13, dxiCP=dxiCP, b1=b1, b2=b2, b3=b3, Lambda=Lambda, n_liv=n_liv, verbose=0))

    # # Four-neutrino oscillations in constant-density matter, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_GEV # [eV]
    # rho = 10.0*gd.UNIT_G_PER_CM3 # [g cm^{-3}]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 = 0.01, 0.05, 0.06, 0.1, 0.9, 0.02
    # dxi13, dxi14, dxi24 = np.radians([10,20,30])
    # b1, b2, b3, b4 = 1.e-2, 1.e-3, 5.e-3, 4.e-2
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 3
    # print(osc_prob_4nu_matter_constant_density(energy, baseline, rho, s14=s14, s24=s24, s34=s34,
    #     d14=d14, d24=d24, D41=D41, verbose=0))
    # print(osc_prob_4nu_matter_liv_constant_density(energy, baseline, rho, s14=s14, s24=s24, s34=s34,
    #     d14=d14, d24=d24, D41=D41, sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, 
    #     sxi34=sxi34, dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, 
    #     Lambda=Lambda, n_liv=n_liv, verbose=0))

    # # Five-neutrino oscillations in constant-density matter, LIV
    # np.set_printoptions(precision=3)
    # baseline = 10.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_GEV # [eV]
    # rho = 10.0*gd.UNIT_G_PER_CM3 # [g cm^{-3}]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi15, sxi24, sxi25, sxi34, sxi35 = \
    #     0.01, 0.05, 0.06, 0.1, 0.9, 0.02, 0.04, 0.1, 0.05
    # dxi13, dxi14, dxi15, dxi24, dxi35 = np.radians([10, 20, 30, 40, 50])
    # b1, b2, b3, b4, b5 = 1.e-2, 1.e-3, 5.e-3, 4.e-2, 1.e-1
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 2
    # print(osc_prob_5nu_matter_constant_density(energy, baseline, rho, s14=s14, s15=s15, s24=s24, 
    #     s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51, verbose=0))
    # print(osc_prob_5nu_matter_liv_constant_density(energy, baseline, rho, s14=s14, s15=s15, s24=s24,
    #     s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51, 
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25, 
    #     sxi34=sxi34, sxi35=sxi35, dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, 
    #     b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, Lambda=Lambda, n_liv=n_liv, verbose=0))

    # # Two-neutrino oscillations in exponentially falling matter density profile, LIV
    # np.set_printoptions(precision=3)
    # baseline = 20.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # L0 = 0.0
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # rho_central = 12.0*gd.UNIT_G_PER_CM3
    # l_scale = 4.0*gd.UNIT_KM
    # sxi = 0.5
    # b1 = 1.e-3#1.e-2
    # b2 = 1.e-4#1.e-3
    # Lambda = 10.*gd.UNIT_GEV
    # n_liv = 1
    # print(osc_prob_2nu_matter_constant_density(energy, baseline, rho_central, sth, Dm2, 
    #     verbose=0))
    # print(osc_prob_2nu_matter_exp_density(energy, baseline, L0, rho_central, l_scale, sth, Dm2, 
    #     verbose=0))
    # print(osc_prob_2nu_matter_liv_constant_density(energy, baseline, rho_central, sth, Dm2,
    #     sxi=sxi, b1=b1, b2=b2, Lambda=Lambda, n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))
    # print(osc_prob_2nu_matter_liv_exp_density(energy, baseline, L0, rho_central, l_scale, sth, Dm2,
    #     sxi=sxi, b1=b1, b2=b2, Lambda=Lambda, n_liv=n_liv, rtol=1.e-3, atol=1.e-3, n_jobs=1, 
    #     verbose=0))


    # # Three-neutrino oscillations in exponentially falling matter density profile, LIV
    # np.set_printoptions(precision=3)
    # baseline = 20.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # L0 = 0.0
    # rho_central = 12.0*gd.UNIT_G_PER_CM3
    # l_scale = 4.0*gd.UNIT_KM
    # sxi12, sxi23, sxi13, dxiCP = 0.001, 0.002, 0.003, np.radians(10.0)
    # b1, b2, b3 = 1.e-6, 2.e-6, 3.e-6 #1.e-2, 1.e-3, 5.e-3
    # Lambda = 100.*gd.UNIT_GEV
    # n_liv = 1
    # print(osc_prob_3nu_matter_constant_density(energy, baseline, rho_central, verbose=0))
    # print(osc_prob_3nu_matter_exp_density(energy, baseline, L0, rho_central, l_scale, verbose=0))
    # print(osc_prob_3nu_matter_liv_constant_density(energy, baseline, rho_central, 
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, dxiCP=dxiCP, b1=b1, b2=b2, b3=b3, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))
    # print(osc_prob_3nu_matter_liv_exp_density(energy, baseline, L0, rho_central, l_scale,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, dxiCP=dxiCP, b1=b1, b2=b2, b3=b3, Lambda=Lambda, 
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, n_jobs=1, verbose=0))


    # # Four-neutrino oscillations in exponentially falling matter density profile, LIV
    # np.set_printoptions(precision=3)
    # baseline = 20.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # L0 = 0.0
    # rho_central = 12.0*gd.UNIT_G_PER_CM3
    # l_scale = 4.0*gd.UNIT_KM
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 = 0.01, 0.05, 0.06, 0.01, 0.01, 0.02
    # dxi13, dxi14, dxi24 = np.radians([10,20,30])
    # b1, b2, b3, b4 = 2.e-8, 1.e-8, 5.e-8, 4.e-8 #1.e-2, 1.e-3, 5.e-3, 4.e-2
    # Lambda = 100.*gd.UNIT_GEV
    # n_liv = 1
    # print(osc_prob_4nu_matter_constant_density(energy, baseline, rho_central, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41, verbose=0))
    # print(osc_prob_4nu_matter_exp_density(energy, baseline, L0, rho_central, l_scale, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41, verbose=0))
    # print(osc_prob_4nu_matter_liv_constant_density(energy, baseline, rho_central, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, sxi34=sxi34, 
    #     dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))
    # print(osc_prob_4nu_matter_liv_exp_density(energy, baseline, L0, rho_central, l_scale,
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, sxi34=sxi34, 
    #     dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))


    # # Five-neutrino oscillations in exponentially falling matter density profile, LIV
    # np.set_printoptions(precision=3)
    # baseline = 20.*gd.UNIT_KM # 10 km in natural units [eV^{-1}]
    # energy = 1.*gd.UNIT_MEV # [eV]
    # L0 = 0.0
    # rho_central = 12.0*gd.UNIT_G_PER_CM3
    # l_scale = 4.0*gd.UNIT_KM
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi15, sxi24, sxi25, sxi34, sxi35 = \
    #     0.01, 0.05, 0.06, 0.1, 0.9, 0.02, 0.04, 0.1, 0.05
    # dxi13, dxi14, dxi15, dxi24, dxi35 = np.radians([10, 20, 30, 40, 50])
    # b1, b2, b3, b4, b5 = 1.e-2, 1.e-3, 5.e-3, 4.e-2, 1.e-1
    # Lambda = 100.*gd.UNIT_GEV
    # n_liv = 1
    # print(osc_prob_5nu_matter_constant_density(energy, baseline, rho_central, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51, verbose=0))
    # print(osc_prob_5nu_matter_exp_density(energy, baseline, L0, rho_central, l_scale, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51, verbose=0))
    # print(osc_prob_5nu_matter_liv_constant_density(energy, baseline, rho_central, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25,
    #     sxi34=sxi34, sxi35=sxi35, 
    #     dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, 
    #     b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))
    # print(osc_prob_5nu_matter_liv_exp_density(energy, baseline, L0, rho_central, l_scale,
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25,
    #     sxi34=sxi34, sxi35=sxi35, 
    #     dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, 
    #     b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))


    # # Two-neutrino oscillations in the Sun, LIV
    # np.set_printoptions(precision=3)
    # baseline = gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # sxi = 0.001
    # b1 = 2.e-4#1.e-2
    # b2 = 1.e-4#1.e-3
    # Lambda = 1000.*gd.UNIT_GEV
    # n_liv = 2
    # print(osc_prob_2nu_sun(energy, baseline, L0, sth, Dm2, verbose=0))
    # print(osc_prob_2nu_sun_liv(energy, baseline, L0, sth, Dm2, 
    #     sxi=sxi, b1=b1, b2=b2, Lambda=Lambda, n_liv=n_liv, verbose=0))


    # # Three-neutrino oscillations in the Sun, LIV
    # np.set_printoptions(precision=3)
    # baseline = gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # sxi12, sxi23, sxi13, dxiCP = 0.001, 0.002, 0.003, np.radians(10.0)
    # b1, b2, b3 = 1.e-6, 2.e-6, 3.e-6 #1.e-2, 1.e-3, 5.e-3
    # Lambda = 100.*gd.UNIT_GEV
    # n_liv = 1
    # print(osc_prob_3nu_sun(energy, baseline, L0, verbose=0))
    # print(osc_prob_3nu_sun_liv(energy, baseline, L0,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, dxiCP=dxiCP, b1=b1, b2=b2, b3=b3, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=0))


    # # Four-neutrino oscillations in the Sun, LIV
    # np.set_printoptions(precision=3)
    # baseline = gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 = 0.01, 0.05, 0.06, 0.01, 0.01, 0.02
    # dxi13, dxi14, dxi24 = np.radians([10,20,30])
    # b1, b2, b3, b4 = 2.e-8, 1.e-8, 5.e-8, 4.e-8 #1.e-2, 1.e-3, 5.e-3, 4.e-2
    # Lambda = 100.*gd.UNIT_GEV
    # n_liv = 2
    # print(osc_prob_4nu_sun(energy, baseline, L0, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41, verbose=2))
    # print(osc_prob_4nu_sun_liv(energy, baseline, L0, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, sxi34=sxi34, 
    #     dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-3, atol=1.e-3, verbose=2))


    # # Five-neutrino oscillations in the Sun, LIV
    # np.set_printoptions(precision=3)
    # baseline = gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # energy = 10.*gd.UNIT_MEV # [eV]
    # L0 = 0.1*gd.SUN_RADIUS*gd.UNIT_KM # [eV^{-1}]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # sxi12, sxi23, sxi13, sxi14, sxi15, sxi24, sxi25, sxi34, sxi35 = \
    #     0.01, 0.05, 0.06, 0.05, 0.01, 0.02, 0.04, 0.01, 0.05
    # dxi13, dxi14, dxi15, dxi24, dxi35 = np.radians([10, 20, 30, 40, 50])
    # b1, b2, b3, b4, b5 = 1.e-3, 2.e-3, 3.e-3, 4.e-3, 5.e-3
    # Lambda = 100.*gd.UNIT_GEV
    # n_liv = 1
    # print(osc_prob_5nu_sun(energy, baseline, L0, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51, verbose=0))
    # print(osc_prob_5nu_sun_liv(energy, baseline, L0, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51,
    #     sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25,
    #     sxi34=sxi34, sxi35=sxi35, 
    #     dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, 
    #     b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, Lambda=Lambda,
    #     n_liv=n_liv, rtol=1.e-2, atol=1.e-2, verbose=0))


    # Two-neutrino oscillations in Earth
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    ###
    # L = earth.distance_traveled_inside_earth(-0.05)
    # print(osc_prob_2nu_earth(energy, sth, Dm2, costhz=-0.05, L=L*gd.UNIT_KM, verbose=0))
    # print(osc_prob_2nu_earth(energy, sth, Dm2, costhz=-1, L=L*gd.UNIT_KM, verbose=0))
    ###
    # for costhz in np.linspace(-0.5, -1.0, 2):
    #     print(earth.distance_traveled_inside_earth(costhz))
    #     print(osc_prob_2nu_earth(energy, sth, Dm2, costhz=costhz,
    #         L=earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM,
    #         # L=6371*gd.UNIT_KM, 
    #         verbose=0))
    ###
    # loc_fin = 'fermilab'
    # for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    #     print(loc_ini)
    #     print(osc_prob_2nu_earth(energy, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
    #     print(osc_prob_3nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
    #     print()


    # # Three-neutrino oscillations in Earth
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # ###
    # # L = earth.distance_traveled_inside_earth(-0.05)
    # # print(osc_prob_3nu_earth(energy, costhz=-0.05, L=L*gd.UNIT_KM, verbose=0))
    # # print(osc_prob_3nu_earth(energy, costhz=-1, L=L*gd.UNIT_KM, verbose=0))
    # ###
    # # L = earth.distance_traveled_inside_earth(0)
    # # print(osc_prob_3nu_earth(energy, costhz=0, L=L*gd.UNIT_KM, verbose=0))
    # ###
    # # for costhz in np.linspace(-0.5, -1.0, 2):
    # #     print(earth.distance_traveled_inside_earth(costhz))
    # #     print(osc_prob_3nu_earth(energy, costhz=costhz,
    # #         L=earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM,
    # #         # L=6371*gd.UNIT_KM, 
    # #         verbose=0))
    # ###
    # # costhz = -0.8
    # # print(osc_prob_3nu_earth(energy, costhz=costhz, L=1000*gd.UNIT_KM, magnus_exp_order=4, 
    # #     verbose=0))
    # # print(osc_prob_3nu_earth(energy, costhz=costhz, L=1200*gd.UNIT_KM, magnus_exp_order=4, 
    # #     verbose=0))
    # # print(osc_prob_3nu_earth(energy, costhz=costhz, L=1000*gd.UNIT_KM, magnus_exp_order=5, 
    # #     verbose=0))
    # # print(osc_prob_3nu_earth(energy, costhz=costhz, L=1000*gd.UNIT_KM, magnus_exp_order=4, 
    # #     rtol=1.e-4, atol=1.e-4, verbose=0))
    # ###
    # loc_fin = 'fermilab'
    # for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    #     print(loc_ini)
    #     print(osc_prob_3nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))


    # # Four-neutrino oscillations in Earth
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s24, s34 = 0.1, 0.2, 0.3
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.1 # [eV^2]
    # ###
    # L = earth.distance_traveled_inside_earth(-0.05)
    # print(osc_prob_4nu_earth(energy, s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     costhz=-0.05, L=L*gd.UNIT_KM, verbose=0))
    # print(osc_prob_4nu_earth(energy, s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     costhz=-1, L=L*gd.UNIT_KM, verbose=0))
    # ###
    # loc_fin = 'fermilab'
    # for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    #     print(loc_ini)
    #     print(osc_prob_3nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
    #     print(osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))


    # Five-neutrino oscillations in Earth
    np.set_printoptions(precision=3)
    energy = 10.*gd.UNIT_MEV # [eV]
    s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    D41, D51 = 0.1, 0.01 # [eV^2]
    ###
    L = earth.distance_traveled_inside_earth(-0.05)
    print(osc_prob_5nu_earth(energy, 
        s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
        D41=D41, D51=D51, costhz=-0.05, L=L*gd.UNIT_KM, verbose=0))
    print(osc_prob_5nu_earth(energy, 
        s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
        D41=D41, D51=D51, costhz=-1, L=L*gd.UNIT_KM, verbose=0))
    ###
    loc_fin = 'fermilab'
    for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
        print(loc_ini)
        print(osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
             s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
            loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
        print(osc_prob_5nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
            s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24,
            d35=d35, D41=D41, D51=D51,
            loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
