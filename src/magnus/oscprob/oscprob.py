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
of probabilities whose entry ``[i][j]`` is the probability of a neutrino
produced with flavor ``i`` being detected with flavor ``j``

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
import warnings
from functools import reduce
from joblib import Parallel, delayed
from typing import Optional, Callable, Union, Tuple, List, Dict
from io import TextIOWrapper
from inspect import signature
# import numba as nb

# import numpy.typing

import magnus.magnus as magnus
import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
from magnus import version
from magnus import authors


has_magnus_header_been_printed = False


class ToleranceNotAchievedWarning(UserWarning):
    r"""Warns that the probability returned by :func:`osc_prob` did not
    reach the requested tolerance because a refinement cap was hit
    (max_num_loops, max_n_slabs, or max_n_tpts_per_slab).

    The result may look plausible (it is still exactly unitary) while
    being inaccurate, so this warning is issued regardless of the
    verbosity setting.  Raise the caps, loosen the tolerance, or use
    wider applicability methods for extreme-phase problems (e.g., many
    more slabs for low-energy solar neutrinos).
    """


#-----------------------------------------------------------------------
# Helper functions
#-----------------------------------------------------------------------

def print_banner(file: TextIOWrapper=None):
    if file is None:
        print(gd.cstyle.CBLUEBG + ".----------------------------------------." + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|   __  __                               |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  |  \/  | __ _  __ _ _ __  _   _ ___   |" + gd.cstyle.CEND, 
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  | |  | | (_| | (_| | | | | |_| \__ \  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + r"|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "|                |___/                   |" + gd.cstyle.CEND,
            file=file)
        print(gd.cstyle.CBLUEBG + "'----------------------------------------'" + gd.cstyle.CEND,
            file=file)
    else: 
        print(".----------------------------------------.", file=file)
        print("|   __  __                               |", file=file)
        print(r"|  |  \/  | __ _  __ _ _ __  _   _ ___   |", file=file)
        print(r"|  | |\/| |/ _` |/ _` | '_ \| | | / __|  |", file=file)
        print(r"|  | |  | | (_| | (_| | | | | |_| \__ \  |", file=file)
        print(r"|  |_|  |_|\__,_|\__, |_| |_|\__,_|___/  |", file=file)
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
        print("   n_tpts_per_slab = " + str(n_tpts_per_slab), file=f)
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
                # (np.issubdtype is used instead of the np.float_/np.int_ aliases, which were
                # removed in NumPy 2.0)
                if not (np.issubdtype(np.asarray(energy).dtype, np.floating) or \
                    np.issubdtype(np.asarray(energy).dtype, np.integer)):
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
                if not (np.issubdtype(np.asarray(L).dtype, np.floating) or \
                    np.issubdtype(np.asarray(L).dtype, np.integer)):
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
        # num_flavors exceeds the predefined range: the caller builds its Hamiltonian directly from
        # h_nsi instead of from a flat parameter list, so there is nothing to unpack here.
        return None
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
        # num_flavors exceeds the predefined range: the caller builds its Hamiltonian directly from
        # h_liv instead of from a flat parameter list, so there is nothing to unpack here.
        return None
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

class _PositionProfileCache:
    r"""Memoizes a position-profile function on repeated position grids.

    Across an energy scan, the matter term of the Hamiltonian is evaluated on
    the same position grids for every energy (only the 1/E vacuum part
    changes), and across the adaptive refinement loops the same grids recur
    between neighboring points.  This tiny cache stores the profile values of
    the most recent grids, keyed by the exact grid contents, so the
    (comparatively expensive) density-profile chain runs once per distinct
    grid instead of once per Hamiltonian evaluation.  Scalar evaluations are
    passed through uncached.
    """

    def __init__(self, func: Callable, maxsize: Optional[int]=8):
        self.func = func
        self._cache = {}
        self._keys = []
        self._maxsize = maxsize

    def __call__(self, l: Union[int, float, np.ndarray]):
        if np.ndim(l) == 0:
            return self.func(l)
        l = np.asarray(l, dtype=float)
        key = (l.shape, l.tobytes())
        val = self._cache.get(key)
        if val is None:
            val = np.asarray(self.func(l))
            self._cache[key] = val
            self._keys.append(key)
            if len(self._keys) > self._maxsize:
                self._cache.pop(self._keys.pop(0), None)
        return val


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
            # t_slabs=[t_slab],
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
    r"""Computes the evolution operators of a chain of time slabs.  This function is not designed
    to be called directly by the user, but rather internally by :func:`osc_prob`.

    All slabs are computed at once by :func:`magnus.magnus.magnus_expansion_multislab`, which
    batches the Hamiltonian evaluation, the quadrature, the commutator algebra, and the matrix
    exponentials over the slab axis.  Slabs of zero width yield identity operators.

    :param H_func: Hamiltonian, which is a function of time or position that returns a square matrix
        in the form of NumPy array.  If it also accepts an array of times (returning a stack of
        matrices), the vectorized form is detected and used automatically for speed.
    :param t_slabs: List or NumPy array of pairs specifying the start and end times or positions of
        each slab, i.e., [[t0, t1], [t1, t2], ...]
    :param n_tpts_per_slab: Number of time-points inside each slab at which to evaluate H_func in
        order to numerically compute the integrals over time required by the Magnus expansion
        (ignored by the 'gl' integration method)
    :param magnus_exp_order: Maximum order of Magnus expansion used to compute the evolution
        operator (should not exceed :func:`magnus.globaldefs.MAGNUS_EXP_ORDER_MAX`)
    :param \**kwargs: Additional arguments passed to
        :func:`magnus.magnus.magnus_expansion_multislab` (e.g., integration_method)

    :return: A NumPy array of shape (n_slabs, dim, dim) containing the evolution operators, ordered
        like ``t_slabs`` (earliest slab first).  Note that the time-ordered product over the chain
        is U_total = U[-1] @ ... @ U[1] @ U[0], i.e., the last slab is the leftmost factor.

    """
    def hh(t):
        return -1j * H_func(t)

    return magnus.magnus_expansion_multislab(hh, t_slabs, n_tpts_per_slab=n_tpts_per_slab,
        order=magnus_exp_order, **kwargs)


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
    A_eval_mode: Optional[str]=None,
    convergence_info: Optional[Dict]=None,
    t_breakpoints: Optional[Union[list, np.ndarray]]=None,
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
        Optional list of pairs [[t0, t1], [t1, t2], ...] with the edges
        of each time slab.  If given, it overrides ``n_slabs`` and the
        uniform partitioning of [``t_ini``, ``t_fin``]; the user must
        ensure that the slabs chain without gaps.  If a tolerance is
        requested, only ``n_tpts_per_slab`` is grown (the user-provided
        edges are kept fixed).
    magnus_exp_order
        Order at which the Magnus expansion is truncated (1 to
        ``globaldefs.MAGNUS_EXP_ORDER_MAX``).
    n_jobs
        Number of parallel joblib workers used to compute the per-slab
        evolution operators.  With the default, ``n_jobs = 1``, all
        slabs are computed in a single vectorized (batched) call, which
        is usually fastest; use ``n_jobs > 1`` only for very expensive
        Hamiltonian functions.
    integration_method
        'trapezoid' or 'simpson' for cumulative quadrature over
        ``n_tpts_per_slab`` points per slab, or 'gl' for Gauss-Legendre
        collocation, which needs only 1, 2, or 3 Hamiltonian
        evaluations per slab for orders <= 2, <= 4, <= 6, and ignores
        ``n_tpts_per_slab``.
    rtol
        Target relative tolerance of the probability matrix between
        successive refinement loops.  Set both ``rtol`` and ``atol`` to
        ``None`` to run once with the given fixed parameters.  If only
        one of the two is ``None``, it is treated as 0.
    atol
        Target absolute tolerance; see ``rtol``.
    growth_factor_n_slabs
        Factor by which ``n_slabs`` is multiplied on each refinement
        loop (used only when a tolerance is requested).
    growth_factor_n_tpts_per_slab
        Factor by which ``n_tpts_per_slab`` is multiplied on each
        refinement loop (used only when a tolerance is requested).
    max_num_loops
        Maximum number of refinement loops.
    min_n_slabs
        Number of slabs used in the first refinement loop.
    max_n_slabs
        Maximum allowed number of slabs.
    min_n_tpts_per_slab
        Number of time points per slab in the first refinement loop.
    max_n_tpts_per_slab
        Maximum allowed number of time points per slab.
    iterate_over_magnus_exp_order
        If True, additionally increase ``magnus_exp_order`` from
        ``min_magnus_exp_order`` to ``max_magnus_exp_order`` until the
        requested tolerance is achieved.
    min_magnus_exp_order
        Lowest expansion order tried when iterating over the order.
    max_magnus_exp_order
        Highest expansion order tried when iterating over the order.
    validate_input
        If True, validate the input parameters (set to False for a
        small speed-up once a call is known to be well-formed).
    save_log
        If True, also write all messages to the log file.
    filename_log
        Name of the log file (used if ``save_log`` is True and no
        ``file_log`` object is given).
    file_log
        Optional file object to write log messages to.
    close_file_log_upon_exit
        If True, close the log file before returning.
    new_recursion_limit
        If not None, raise Python's recursion limit to this value.
    verbose
        Verbosity level: 0 (silent), 1 (warnings), 2 (progress of the
        refinement loops).
    A_eval_mode
        How the Hamiltonian can be evaluated: 'vector' (accepts an
        array of positions), 'constant', or 'scalar'.  Determined
        automatically when None; pass it explicitly (e.g., from
        :func:`magnus.magnus.probe_eval_mode`) to skip the probe.
    convergence_info
        If a dict is passed, it is filled in place with the refinement
        parameters of the returned probability ('n_slabs',
        'n_tpts_per_slab'), which callers can use to warm-start
        neighboring computations.
    t_breakpoints
        Optional positions at which the Hamiltonian is known to be
        non-smooth (e.g., density discontinuities such as the PREM
        layer boundaries).  They are inserted as mandatory slab edges
        into the automatically generated slab grid at every refinement
        level, so that the quadrature never integrates across them.
        Ignored when ``t_slab_edges`` is given explicitly.
    \**kwargs
        Additional arguments passed through to the Magnus-expansion
        routines

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

    # If only one of rtol and atol was given (i.e., the other one is None), set the missing one to
    # 0.0, so that the requested tolerance is driven by the one that was given.  (Internally, the
    # code treats "a tolerance was requested" as both rtol and atol being not None.)
    if (rtol is None) != (atol is None):
        rtol = 0.0 if rtol is None else rtol
        atol = 0.0 if atol is None else atol

    # The Gauss-Legendre integration method ('gl') uses a fixed, small number of Hamiltonian
    # evaluations per slab (1, 2, or 3, depending on magnus_exp_order), so n_tpts_per_slab plays no
    # role: the accuracy is controlled by the number of slabs only.  Neutralize the growth of
    # n_tpts_per_slab so that the adaptive loop below grows only n_slabs.
    if integration_method == 'gl':
        growth_factor_n_tpts_per_slab = 1.0
        min_n_tpts_per_slab = max_n_tpts_per_slab = n_tpts_per_slab = 2

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

    # Note: new_recursion_limit is accepted for backward compatibility but no longer used; the
    # probability calculation is fully iterative (nothing recurses), so there is no need to raise
    # Python's recursion limit.

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
        # A single slab is exact for a constant Hamiltonian, so drop any user-provided slab edges
        t_slab_edges = None
        t_slab_edges_original = None
        if verbose > 0:
            for f in [None, file_log] if save_log else [None]:
                warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
                print("\n" + warn_msg + " The provided Hamiltonian is time-independent. " + \
                    "Overwriting the run parameters to magnus_exp_order = 1, n_slabs = 1, " + \
                    "n_tpts_per_slab = 2, rtol = None, atol = None, and n_jobs = 1 for speed-up.",
                    file=f)

    # If the user provided the slab edges explicitly, the number of slabs is set by them
    if t_slab_edges_original is not None:
        n_slabs = len(t_slab_edges_original)

    # Determine once how the Hamiltonian can be evaluated (vectorized over an array of positions,
    # constant, or scalar-only), so that the Magnus kernel does not have to re-probe it on every
    # refinement iteration below.
    if A_eval_mode is None:
        A_eval_mode = magnus.probe_eval_mode(lambda t: -1j*H_func(t), t_ini, t_fin)

    # Physics-informed starting number of slabs (Gauss-Legendre method only): rather than always
    # starting the refinement from min_n_slabs and climbing the geometric ladder, start from an
    # estimate based on the accumulated phase (see magnus.suggest_n_slabs).  min_n_slabs still
    # acts as a lower bound, so warm starts provided by the caller take precedence when they are
    # larger.  For the quadrature methods ('trapezoid', 'simpson') the accuracy is governed
    # jointly by n_slabs and n_tpts_per_slab, and seeding only the slab count unbalances that
    # ladder, so the seed is not applied there.
    if ((rtol is not None) and (atol is not None) and (t_slab_edges_original is None) and \
        (integration_method == 'gl')):
        n_slabs = int(np.clip(max(min_n_slabs,
            magnus.suggest_n_slabs(lambda t: -1j*H_func(t), t_ini, t_fin,
                A_eval_mode=A_eval_mode)), 1, max_n_slabs))

    while True:

        # These checks only apply when osc_prob is run with a requested tolerance (rtol, atol) that
        # should be achieved.
        if ((rtol is not None) and (atol is not None)):
            # Reached maximum allowed number of loops: exit loop, return the probability matrix
            if (loop_count > max_num_loops):
                warnings.warn("osc_prob: requested tolerance not achieved "
                    "(max_num_loops reached); the returned probabilities may be "
                    "inaccurate. Try increasing max_num_loops. Shown once per "
                    "session.", ToleranceNotAchievedWarning, stacklevel=2)
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
                warnings.warn("osc_prob: requested tolerance not achieved "
                    "(max_n_slabs and max_n_tpts_per_slab reached); the returned "
                    "probabilities may be inaccurate. Try increasing max_n_slabs "
                    "or max_n_tpts_per_slab. This can happen for very large "
                    "accumulated phases, e.g., low-energy neutrinos over very "
                    "long baselines. Shown once per session.",
                    ToleranceNotAchievedWarning, stacklevel=2)
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
            # n_slabs of time slabs.  Any t_breakpoints inside the interval (e.g., density
            # discontinuities) are inserted as additional mandatory slab edges: high-order
            # quadrature converges at its nominal order only if the Hamiltonian is smooth inside
            # each slab.
            grid = np.linspace(t_ini, t_fin, n_slabs+1)
            if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
                bp = np.atleast_1d(np.asarray(t_breakpoints, dtype=float))
                bp = bp[(bp > t_ini) & (bp < t_fin)]
                grid = np.unique(np.concatenate([grid, bp]))
            t_slab_edges = np.column_stack([grid[:-1], grid[1:]])

        # Within each slab, t_slab, we use n_tpts_per_slab time-evaluations to compute the integrals
        # of the Magnus expansion, from t_slab[0] to t_slab[1].  U_chain contains the chain of time-
        # ordered evolution operators, each computed in one time slab.  All slabs are computed in a
        # single batched call.  (Note: n_jobs is accepted for backward compatibility, but the
        # per-slab parallelization it used to trigger here has been retired: the batched kernel is
        # faster than distributing the small per-slab tasks over joblib workers.  Parallelism over
        # (energy, L) points is available in osc_prob_energy_baseline instead.)
        U_chain = compute_evolution_operator_multiple_slabs(H_func, t_slab_edges,
            n_tpts_per_slab, magnus_exp_order, integration_method=integration_method,
            A_eval_mode=A_eval_mode, **kwargs)

        # Now compute the time-ordered product of all evolution operators across all slabs.  The
        # neutrino crosses the slabs in the order in which they appear in U_chain (earliest first),
        # so the total operator is U_tot = U_chain[-1] @ ... @ U_chain[1] @ U_chain[0]: the operator
        # of the *last* slab is the leftmost factor.  (functools.reduce is used instead of
        # np.linalg.multi_dot because all factors are square matrices of the same size, for which
        # multi_dot wastes time computing an optimal parenthesization that does not exist.)
        Utot = reduce(np.matmul, U_chain[::-1]) if len(U_chain) > 1 else U_chain[0]

        # Using Utot, compute all the survival and transition probabilities in a probability matrix
        # P = (|Utot|^2).T and return that matrix, so that P[nu_i][nu_f] = |Utot[nu_f][nu_i]|^2.
        P = np.transpose(Utot.real**2 + Utot.imag**2)

        # Record the refinement parameters of this (latest) computation, so that callers (e.g.,
        # osc_prob_energy_baseline) can warm-start neighboring points
        if convergence_info is not None:
            convergence_info['n_slabs'] = n_slabs
            convergence_info['n_tpts_per_slab'] = n_tpts_per_slab

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
    
    # If the for loop finishes, then it means that the requested tolerance could not be achieved
    # using even the maximum magnus_exp_order allowed for the run.  Return the probability matrix,
    # but show a warning (if verbose).
    if (verbose > 0):
        for f in [None, file_log] if save_log else [None]:
            warn_msg = gd.WARNING_MSG_IN_COLOR if f is None else gd.WARNING_MSG_NO_COLOR
            print(warn_msg + " returning probability, but requested tolerance not achieved using" +\
                " even the maximum allowed order of the Magnus expansion for this run " + \
                "(max_magnus_exp_order = " + str(max_magnus_exp_order) + ").  Try increasing " + \
                "max_n_slabs, max_n_tpts_per_slab, or max_num_loops.\n", file=f)
    if save_log and (file_log is not None): file_log.close()
    return P


def _normalize_energy_L(
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, bool, bool]:
    r"""Normalize energy and L to same-length 1D arrays.

    Returns (energy, L, return_float, ok): return_float records whether both
    inputs were scalars (so that the caller returns a scalar-like result),
    and ok whether the input lengths were compatible.
    """
    energy = float(energy) if isinstance(energy, int) else energy
    L = float(L) if isinstance(L, int) else L
    return_float = isinstance(energy, float) and isinstance(L, float)
    energy = np.array([energy]) if isinstance(energy, float) else np.array(energy)
    L = np.array([L]) if isinstance(L, float) else np.array(L)
    ok = ((len(energy) == len(L)) or (len(energy) == 1 and len(L) > 1) or
          (len(energy) > 1 and len(L) == 1))
    if ok:
        energy = np.full(len(L), energy[0]) if (len(energy) == 1) else energy
        L = np.full(len(energy), L[0]) if (len(L) == 1) else L
    return energy, L, return_float, ok


def _osc_prob_scan_separable(
    H_E: np.ndarray,
    VCC_func: Callable,
    h_matt: np.ndarray,
    L0: float,
    L_val: float,
    t_breakpoints: Optional[np.ndarray],
    magnus_exp_order: int,
    integration_method: str,
    rtol: Optional[float],
    atol: Optional[float],
    growth_factor_n_slabs: float,
    growth_factor_n_tpts_per_slab: float,
    max_num_loops: int,
    min_n_slabs: int,
    max_n_slabs: int,
    min_n_tpts_per_slab: int,
    max_n_tpts_per_slab: int,
    n_slabs: int,
    n_tpts_per_slab: int
) -> np.ndarray:
    r"""Energy-batched probability scan for separable Hamiltonians.

    Computes the probabilities of many neutrino energies that share the same
    baseline [``L0``, ``L_val``] in one batched pipeline, for Hamiltonians of
    the separable form

        H(E, l) = H_E(E) + VCC(l) * h_matt ,

    where ``H_E`` (shape (nE, d, d)) collects all the position-independent,
    energy-dependent terms (vacuum, LIV, ...), ``VCC_func`` is the scalar
    matter potential along the trajectory, and ``h_matt`` (shape (d, d)) is
    the constant matter matrix it multiplies.  The position samples of the
    potential are computed once per refinement level and shared by all
    energies, and the Magnus kernel (quadrature, commutators, exponentials,
    slab products) runs with the energy axis batched in front of the slab
    axis.

    The adaptive refinement mirrors :func:`osc_prob`: the slab count (and,
    for the quadrature methods, the points per slab) grows geometrically
    until the probabilities of each energy agree between successive levels
    within (rtol, atol); converged energies drop out of the batch.  Energies
    are processed in chunks to bound the memory of the sample array.

    Returns the stacked probability matrices, shape (nE, d, d).
    """
    nE, dim = H_E.shape[0], H_E.shape[-1]
    tol_requested = ((rtol is not None) and (atol is not None))

    if integration_method == 'gl':
        # The accuracy of the GL method is controlled by n_slabs only
        growth_factor_n_tpts_per_slab = 1.0
        min_n_tpts_per_slab = max_n_tpts_per_slab = n_tpts_per_slab = 2
        s_nodes = magnus.gl_nodes(magnus_exp_order)

    if tol_requested:
        n_tpts_per_slab = min_n_tpts_per_slab
        # Physics-informed starting number of slabs (see magnus.suggest_n_slabs):
        # integral of the traceless Hamiltonian over the trajectory, maximized
        # over the energies of the scan
        if integration_method == 'gl':
            ts = np.linspace(L0, L_val, 17)
            V17 = np.asarray(VCC_func(ts))
            I_V = (np.sum(V17) - 0.5*(V17[0] + V17[-1]))*(L_val - L0)/16.0
            M = (L_val - L0)*H_E + I_V*h_matt
            M = M - (np.trace(M, axis1=-2, axis2=-1)/dim)[:, None, None]*np.eye(dim)
            try:
                phase = np.max(np.linalg.svd(M, compute_uv=False))
            except np.linalg.LinAlgError:
                phase = 0.0
            n_slabs = int(np.clip(max(min_n_slabs,
                np.ceil(phase/(2.0*np.pi))), 1, max_n_slabs))
        else:
            n_slabs = min_n_slabs

    P_prev = np.full((nE, dim, dim), np.nan)
    P_out = np.empty((nE, dim, dim))
    active = np.arange(nE)
    mA = -1j*h_matt.astype(complex)
    HE_c = -1j*H_E.astype(complex)

    loop_count = 1
    while True:
        # Slab grid shared by all energies (PREM-layer breakpoints included)
        grid = np.linspace(L0, L_val, n_slabs + 1)
        if (t_breakpoints is not None) and (len(np.atleast_1d(t_breakpoints)) > 0):
            bp = np.atleast_1d(np.asarray(t_breakpoints, dtype=float))
            bp = bp[(bp > L0) & (bp < L_val)]
            grid = np.unique(np.concatenate([grid, bp]))
        edges = np.column_stack([grid[:-1], grid[1:]])
        widths = edges[:, 1] - edges[:, 0]

        if integration_method == 'gl':
            s = s_nodes
        else:
            s = np.linspace(0.0, 1.0, n_tpts_per_slab)
        tgrid = edges[:, :1] + widths[:, None]*s              # (n_slabs, m)
        V = np.asarray(VCC_func(tgrid.ravel())).reshape(tgrid.shape)
        Vmat = V[:, :, None, None]*mA                         # (n_slabs, m, d, d)

        # Batched kernel over the active energies, chunked so that each
        # sample array At holds at most ~4M complex entries (~64 MB)
        chunk = max(1, int(4_194_304 // max(1, tgrid.size*dim*dim)))
        P_new = np.empty((len(active), dim, dim))
        for i0 in range(0, len(active), chunk):
            sel = active[i0:i0+chunk]
            At = HE_c[sel][:, None, None, :, :] + Vmat[None, :, :, :, :]
            U = magnus.evolution_operators_from_samples(At, widths,
                magnus_exp_order, integration_method, validate_input=False)
            Utot = U[:, -1]
            for k in range(U.shape[1] - 2, -1, -1):
                Utot = Utot @ U[:, k]
            P_new[i0:i0+chunk] = np.swapaxes(
                Utot.real**2 + Utot.imag**2, -1, -2)

        if not tol_requested:
            P_out[active] = P_new
            return P_out

        prev = P_prev[active]
        have_prev = ~np.isnan(prev[:, 0, 0])
        conv = have_prev & np.all(np.abs(P_new - prev) <= atol + rtol*np.abs(prev),
                                  axis=(-1, -2))
        P_out[active[conv]] = P_new[conv]
        P_prev[active] = P_new
        active = active[~conv]
        if active.size == 0:
            return P_out

        at_caps = ((n_slabs >= max_n_slabs) and
                   (n_tpts_per_slab >= max_n_tpts_per_slab))
        if (loop_count >= max_num_loops) or at_caps:
            warnings.warn("osc_prob (energy-batched scan): requested tolerance "
                "not achieved for some energies (refinement caps reached); the "
                "returned probabilities may be inaccurate. Try increasing "
                "max_n_slabs, max_n_tpts_per_slab, or max_num_loops. Shown "
                "once per session.", ToleranceNotAchievedWarning, stacklevel=2)
            P_out[active] = P_new[~conv]
            return P_out

        n_slabs_old = n_slabs
        n_slabs = min(round(growth_factor_n_slabs*n_slabs), max_n_slabs)
        if ((growth_factor_n_slabs > 1.0) and (n_slabs < max_n_slabs) and
                (n_slabs == n_slabs_old)):
            n_slabs += 1
        n_tpts_old = n_tpts_per_slab
        n_tpts_per_slab = min(int(growth_factor_n_tpts_per_slab*n_tpts_per_slab),
                              max_n_tpts_per_slab)
        if ((growth_factor_n_tpts_per_slab > 1.0) and
                (n_tpts_per_slab < max_n_tpts_per_slab) and
                (n_tpts_per_slab == n_tpts_old)):
            n_tpts_per_slab += 1
        loop_count += 1


def _osc_prob_scan_separable_dispatch(
    h_vac_energy_indep: np.ndarray,
    VCC_func: Union[Callable, float],
    h_matt: np.ndarray,
    h_liv_energy_indep: Optional[np.ndarray],
    n_liv: Optional[Union[int, float]],
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    scan_kwargs: Dict
):
    r"""Decide whether the energy-batched scan engine applies; run it if so.

    Returns NotImplemented when the request does not fit the engine (single
    point, per-point baselines, user-provided slab edges, parallel or logged
    runs, iteration over the expansion order, or unknown extra arguments), in
    which case the caller falls back to the generic per-point path.
    """
    kwargs = dict(scan_kwargs.get('kwargs', {}))
    t_breakpoints = kwargs.pop('t_breakpoints', None)
    n_slabs = kwargs.pop('n_slabs', 1)
    n_tpts_per_slab = kwargs.pop('n_tpts_per_slab', 100)
    if len(kwargs) > 0:
        return NotImplemented
    if not isinstance(VCC_func, Callable):
        return NotImplemented
    if scan_kwargs['t_slab_edges'] is not None:
        return NotImplemented
    if scan_kwargs['iterate_over_magnus_exp_order']:
        return NotImplemented
    if (scan_kwargs['n_jobs'] != 1) or scan_kwargs['save_log'] or \
            (scan_kwargs['file_log'] is not None):
        return NotImplemented

    energy_arr, L_arr, return_float, ok = _normalize_energy_L(energy, L)
    if (not ok) or (len(energy_arr) < 2) or (not np.all(L_arr == L_arr[0])):
        return NotImplemented

    rtol, atol = scan_kwargs['rtol'], scan_kwargs['atol']
    if (rtol is None) != (atol is None):
        rtol = 0.0 if rtol is None else rtol
        atol = 0.0 if atol is None else atol

    # All the position-independent, energy-dependent terms of the Hamiltonian
    H_E = (1.0/energy_arr)[:, None, None]*np.asarray(h_vac_energy_indep)
    if h_liv_energy_indep is not None:
        H_E = H_E + (energy_arr**n_liv)[:, None, None]*np.asarray(h_liv_energy_indep)

    P = _osc_prob_scan_separable(H_E, VCC_func, np.asarray(h_matt), float(L0),
        float(L_arr[0]), t_breakpoints, scan_kwargs['magnus_exp_order'],
        scan_kwargs['integration_method'], rtol, atol,
        scan_kwargs['growth_factor_n_slabs'],
        scan_kwargs['growth_factor_n_tpts_per_slab'],
        scan_kwargs['max_num_loops'], scan_kwargs['min_n_slabs'],
        scan_kwargs['max_n_slabs'], scan_kwargs['min_n_tpts_per_slab'],
        scan_kwargs['max_n_tpts_per_slab'], n_slabs, n_tpts_per_slab)

    if (nu_i is not None) and (nu_f is not None):
        P = P[:, nu_i, nu_f]
    return P.__getitem__(0 if return_float else slice(None))


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

    # Turn int into float
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

    n_points = len(energy)

    # When there are multiple (energy, L) points and n_jobs != 1, parallelize over the points, and
    # run each individual osc_prob call serially.  The per-point tasks are large enough for
    # process-based parallelism to pay off, unlike the much smaller per-slab tasks inside osc_prob.
    parallelize_over_points = (n_jobs != 1) and (n_points > 1)

    # Keyword arguments common to all the calls to osc_prob below.  Additional keyword arguments
    # received in **kwargs are passed through to osc_prob as well (e.g., n_slabs,
    # n_tpts_per_slab).
    osc_prob_kwargs = dict(
        t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order,
        n_jobs=1 if parallelize_over_points else n_jobs,
        integration_method=integration_method, rtol=rtol, atol=atol,
        growth_factor_n_slabs=growth_factor_n_slabs,
        growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
        max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
        min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
        iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
        min_magnus_exp_order=min_magnus_exp_order, max_magnus_exp_order=max_magnus_exp_order,
        validate_input=validate_input, save_log=save_log, filename_log=filename_log,
        file_log=file_log, close_file_log_upon_exit=close_file_log_upon_exit,
        new_recursion_limit=new_recursion_limit, verbose=verbose, **kwargs)

    # Build, for a given neutrino energy, the Hamiltonian to be passed to osc_prob: either a
    # one-parameter function of position or, if position-independent, a constant matrix (which
    # osc_prob detects and handles with internal speed-ups).
    if not isinstance(H_func, Callable):
        # H_func is position- and energy-independent
        def H_at_energy(enu: float) -> np.ndarray:
            return H_func
    elif (len(signature(H_func).parameters) == 2):
        # H_func is a function of two parameters; it is assumed that the first parameter is the
        # energy and the second one is the position
        def H_at_energy(enu: float) -> Callable:
            return lambda l: H_func(enu, l)
    elif H_func_is_function_only_of_energy:
        # H_func is a function only of energy: at fixed energy, it is a constant matrix
        def H_at_energy(enu: float) -> np.ndarray:
            return H_func(enu)
    else:
        # H_func is a function only of position
        def H_at_energy(enu: float) -> Callable:
            return H_func

    # Probe once how the Hamiltonian can be evaluated (vectorized over an array of positions,
    # constant, or scalar-only): the verdict is structural and holds for every (energy, L) point,
    # so probing here avoids re-probing inside every osc_prob call.
    H_first = H_at_energy(energy[0])
    if isinstance(H_first, Callable):
        osc_prob_kwargs['A_eval_mode'] = magnus.probe_eval_mode(
            lambda t: -1j*H_first(t), L0, np.max(L))

    # Warm starts: osc_prob reports the refinement parameters at which each point converged
    # (conv_info), and the next point starts its refinement from there (divided by one growth
    # factor, so that the comparison between successive refinements is still performed).
    # Neighboring points typically converge at (nearly) the same parameters, so this skips most
    # of the refinement ladder.
    conv_info = osc_prob_kwargs.get('convergence_info')
    if conv_info is None:
        conv_info = {}
    osc_prob_kwargs['convergence_info'] = conv_info
    warm_start = (t_slab_edges is None) and \
        ((rtol is not None) or (atol is not None))

    def apply_warm_start():
        # Seed the next point TWO growth steps below the last converged values.  One step below
        # reproduces exactly the pair of refinements at which the previous point was accepted;
        # starting one step lower than that lets the refinement scale decay geometrically across
        # points when the previous point was harder than the next ones (e.g., the lowest energy
        # of a scan), at the price of at most one extra refinement when it was not.
        if warm_start and conv_info:
            g1 = max(growth_factor_n_slabs, 1.0)**2
            g2 = max(growth_factor_n_tpts_per_slab, 1.0)**2
            osc_prob_kwargs['min_n_slabs'] = max(min_n_slabs,
                int(np.ceil(conv_info['n_slabs']/g1)))
            osc_prob_kwargs['min_n_tpts_per_slab'] = max(min_n_tpts_per_slab,
                int(np.ceil(conv_info['n_tpts_per_slab']/g2)))

    def compute_single_point(enu: float, baseline: float) -> Union[float, np.ndarray]:
        P = osc_prob(H_at_energy(enu), L0, baseline, **osc_prob_kwargs)
        # Select one oscillation channel if requested; otherwise return the full matrix
        if ((nu_i is not None) and (nu_f is not None)):
            return P[nu_i][nu_f]
        return P

    if parallelize_over_points:
        # Compute the first point serially to learn the refinement parameters, then distribute
        # the remaining points over the workers, warm-started from the first point.  (The shared
        # conv_info dict cannot be updated across processes, so it is dropped from the parallel
        # calls.)
        probs = [compute_single_point(energy[0], L[0])]
        apply_warm_start()
        osc_prob_kwargs.pop('convergence_info', None)
        probs += Parallel(n_jobs=n_jobs)(delayed(compute_single_point)(enu, baseline)
            for enu, baseline in zip(energy[1:], L[1:]))
    else:
        probs = []
        for enu, baseline in zip(energy, L):
            apply_warm_start()
            probs.append(compute_single_point(enu, baseline))

    # The call to __getitem__ below is a way to return a single float (or single probability
    # matrix) if both energy and L were given as floats.
    return np.array(probs).__getitem__(0 if return_float else slice(None))


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
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar) 

    def htot(enu: Union[int, float]) -> np.ndarray:
        return (1/enu)*h_vac_energy_indep

    htot_is_function_only_of_energy = True

    # Generate the probabilities for all pairs of energy and baseline in zip(energy, L).  (The
    # Hamiltonian is constant in position, so osc_prob computes each point exactly with a single
    # slab; the tolerance and refinement parameters play no role and are not forwarded.)
    return osc_prob_energy_baseline(htot, energy, L, 0.0, nu_i, nu_f,
        htot_is_function_only_of_energy, n_jobs=n_jobs, validate_input=validate_input,
        verbose=verbose, **kwargs)


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
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)

    # Build the coherent forward potential function, VCC_func, from the density function, rho_func.
    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3].
    VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
        electron_fraction, nubar, density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons) # [eV]

    # Projector onto the nu_e--nu_e entry, multiplied below by the potential VCC.  Note that
    # VCC_func already carries the antineutrino sign flip (applied inside
    # matter.vcc_func_from_rho_func), so no extra sign is applied here.  [Previously, the sign was
    # applied twice, which gave the antineutrino matter potential the wrong (positive) sign.]
    h_matt_proj = np.zeros((num_flavors, num_flavors))
    h_matt_proj[0][0] = 1.0

    # Cache repeated evaluations of the potential on identical position grids (see
    # _PositionProfileCache)
    if isinstance(VCC_func, Callable):
        VCC_func = _PositionProfileCache(VCC_func)

    # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
    if isinstance(VCC_func, Callable):
        # VCC_func is a function of position, so the Hamiltonian is, too.  If l is an array, the
        # result is a stack of Hamiltonians with the position axis leading; this lets the Magnus
        # routines evaluate the Hamiltonian at all time points in a single vectorized call.
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            vcc = np.asarray(VCC_func(l))
            return (1/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt_proj
        htot_is_function_only_of_energy = False
    else:
        # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is passed to
        # osc_prob below, osc_prob will detect that VCC_func is constant and set parameters
        # internally for speed-up.
        h_matt = VCC_func*h_matt_proj
        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = True

    # Energy-batched fast path: when many energies share a single baseline and the Hamiltonian
    # is position-dependent, compute the whole scan in one batched pipeline, with the potential
    # samples shared across energies (see _osc_prob_scan_separable).  If the request does not fit
    # the engine, fall back to the generic per-point path below.
    P_scan = _osc_prob_scan_separable_dispatch(h_vac_energy_indep, VCC_func, h_matt_proj, None, None,
        energy, L, L0, nu_i, nu_f,
        dict(t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
             integration_method=integration_method, rtol=rtol, atol=atol,
             growth_factor_n_slabs=growth_factor_n_slabs,
             growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
             max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
             min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
             iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
             save_log=save_log, file_log=file_log, kwargs=kwargs))
    if P_scan is not NotImplemented:
        return P_scan

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
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)

    # Compute the standard + NSI matter Hamiltonian *without* the multiplicative prefactor of VCC.
    # To do this we call the functions hamiltonians_Xnu_nsi(VCC, ...) with VCC = 1.0.  We add the
    # standard matter contribution to the NSI matter contribution by adding 1.0 to the eps_ee entry.
    # The overall antineutrino sign flip is carried by VCC_func (see
    # matter.vcc_func_from_rho_func); for antineutrinos, the NSI couplings are additionally
    # conjugated (H_matt -> -H_matt^* relative to neutrinos).
    if num_flavors == 2:
        h_matt = np.diag([1.0, 0.0]) + \
            hamiltonians.hamiltonian_2nu_nsi(1.0, eps_aa, eps_ab) # VCC = 1.0
    elif num_flavors == 3:
        h_matt = np.diag([1.0, 0.0, 0.0]) + \
            hamiltonians.hamiltonian_3nu_nsi(1.0, eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt)
    elif num_flavors == 4:
        h_matt = np.diag([1.0, 0.0, 0.0, 0.0]) + \
            hamiltonians.hamiltonian_4nu_nsi(1.0, eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt,
                eps_ms, eps_tt, eps_ts, eps_ss)
    elif num_flavors == 5:
        h_matt = np.diag([1.0, 0.0, 0.0, 0.0, 0.0]) + \
            hamiltonians.hamiltonian_5nu_nsi(1.0, eps_ee, eps_em, eps_et, eps_es1, eps_es2,
                eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2,
                eps_s2s2)

    if nubar:
        h_matt = np.conj(h_matt)

    # Build the coherent forward potential function, VCC_func, from the density function, rho_func.
    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3].
    VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
        electron_fraction, nubar, density_matter_is_in_g_per_cm3,
        density_is_of_number_of_electrons) # [eV] 
    
    # Cache repeated evaluations of the potential on identical position grids (see
    # _PositionProfileCache)
    if isinstance(VCC_func, Callable):
        VCC_func = _PositionProfileCache(VCC_func)

    # Matter Hamiltonian function: (standard + NSI) matter matrix scaled by VCC
    if isinstance(VCC_func, Callable):
        # VCC_func is a function of position, so the Hamiltonian is, too.  If l is an array, the
        # result is a stack of Hamiltonians with the position axis leading; this lets the Magnus
        # routines evaluate the Hamiltonian at all time points in a single vectorized call.
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            vcc = np.asarray(VCC_func(l))
            return (1/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt
        htot_is_function_only_of_energy = False
    else:
        # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is passed to
        # osc_prob below, osc_prob will detect that VCC_func is constant and set parameters
        # internally for speed-up.
        h_matt = VCC_func*h_matt
        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep+h_matt
        htot_is_function_only_of_energy = True

    # Energy-batched fast path: when many energies share a single baseline and the Hamiltonian
    # is position-dependent, compute the whole scan in one batched pipeline, with the potential
    # samples shared across energies (see _osc_prob_scan_separable).  If the request does not fit
    # the engine, fall back to the generic per-point path below.
    P_scan = _osc_prob_scan_separable_dispatch(h_vac_energy_indep, VCC_func, h_matt, None, None,
        energy, L, L0, nu_i, nu_f,
        dict(t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
             integration_method=integration_method, rtol=rtol, atol=atol,
             growth_factor_n_slabs=growth_factor_n_slabs,
             growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
             max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
             min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
             iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
             save_log=save_log, file_log=file_log, kwargs=kwargs))
    if P_scan is not NotImplemented:
        return P_scan

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
        h_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) 
    elif num_flavors == 3:
        h_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, D21, D31, nubar=nubar) 
    elif num_flavors == 4:
        h_vac_energy_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, 
            s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, nubar=nubar) 
    elif num_flavors == 5:
        h_vac_energy_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23,
            s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, D21, D31, D41, D51,
            nubar=nubar)
    
    # Compute the energy-independent part of the LIV Hamiltonian, i.e., everything but the 1/E 
    # prefactor, only once, to save time.  Multiply by the 1/E factor later when calling osc_prob.
    # If num_flavors > MAGNUS_MAX_PREDEFINED_NUM_FLAVORS, we use the h_liv_energy_indep that was
    # passed to the function.
    if num_flavors == 2:
        h_liv_energy_indep = hamiltonians.hamiltonian_2nu_liv_energy_independent(sxi, b1, b2, 
            Lambda, n_liv)
    elif num_flavors == 3:
        h_liv_energy_indep = hamiltonians.hamiltonian_3nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxiCP, b1, b2, b3, Lambda, n_liv, nubar=nubar)
    elif num_flavors == 4:
        h_liv_energy_indep = hamiltonians.hamiltonian_4nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxi13, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, n_liv,
            nubar=nubar)
    elif num_flavors == 5:
        h_liv_energy_indep = hamiltonians.hamiltonian_5nu_liv_energy_independent(sxi12, sxi23,
            sxi13, dxi13, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, dxi35, b1,
            b2, b3, b4, b5, Lambda, n_liv, nubar=nubar)
   
    if (rho_func != 0.0): # Matter density is nonzero, include the matter term in the Hamiltonian

        # Projector onto the nu_e--nu_e entry, multiplied below by the potential VCC.  Note that
        # VCC_func already carries the antineutrino sign flip (applied inside
        # matter.vcc_func_from_rho_func), so no extra sign is applied here.
        h_matt = np.zeros((num_flavors, num_flavors))
        h_matt[0][0] = 1.0

        # Build the coherent forward potential function, VCC_func, from the density function,
        # rho_func. If the provided rho_func is the matter density (e.g., g cm^{-3}), convert
        # rho_func to a function that returns the electron number density [eV^3].
        VCC_func = matter.vcc_func_from_rho_func(rho_func, L0, ratio_number_neutrons_to_protons,
            electron_fraction, nubar, density_matter_is_in_g_per_cm3,
            density_is_of_number_of_electrons) # [eV]

        # Cache repeated evaluations of the potential on identical position grids (see
        # _PositionProfileCache)
        if isinstance(VCC_func, Callable):
            VCC_func = _PositionProfileCache(VCC_func)

        # Matter Hamiltonian function: diagonal matrix with VCC in the top-left (ee) entry
        if isinstance(VCC_func, Callable):
            # VCC_func is a function of position, so the Hamiltonian is, too.  If l is an array,
            # the result is a stack of Hamiltonians with the position axis leading; this lets the
            # Magnus routines evaluate the Hamiltonian at all time points in a single call.
            def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
                vcc = np.asarray(VCC_func(l))
                return (1/enu)*h_vac_energy_indep + vcc[..., None, None]*h_matt + \
                    pow(enu,n_liv)*h_liv_energy_indep
            htot_is_function_only_of_energy = False
        else:
            # VCC_func is a constant in position, so the Hamiltonian is, too. When VCC_func is
            # passed to osc_prob below, osc_prob will detect that VCC_func is constant and set
            # parameters  internally for speed-up.
            h_matt = VCC_func*h_matt
            def htot(enu: Union[int, float]) -> np.ndarray:
                return (1/enu)*h_vac_energy_indep + h_matt + pow(enu,n_liv)*h_liv_energy_indep
            htot_is_function_only_of_energy = True

    else: # Matter density is zero; the only terms in the Hamiltonian are vacuum and LIV

        def htot(enu: Union[int, float]) -> np.ndarray:
            return (1/enu)*h_vac_energy_indep + pow(enu,n_liv)*h_liv_energy_indep
        htot_is_function_only_of_energy = True

    # Energy-batched fast path: when many energies share a single baseline and the Hamiltonian
    # is position-dependent, compute the whole scan in one batched pipeline, with the potential
    # samples shared across energies (see _osc_prob_scan_separable).  If the request does not fit
    # the engine, fall back to the generic per-point path below.
    P_scan = NotImplemented
    if (rho_func != 0.0):  # VCC_func and h_matt exist only when there is matter
        P_scan = _osc_prob_scan_separable_dispatch(h_vac_energy_indep, VCC_func, h_matt,
            h_liv_energy_indep, n_liv, energy, L, L0, nu_i, nu_f,
            dict(t_slab_edges=t_slab_edges, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                 integration_method=integration_method, rtol=rtol, atol=atol,
                 growth_factor_n_slabs=growth_factor_n_slabs,
                 growth_factor_n_tpts_per_slab=growth_factor_n_tpts_per_slab,
                 max_num_loops=max_num_loops, min_n_slabs=min_n_slabs, max_n_slabs=max_n_slabs,
                 min_n_tpts_per_slab=min_n_tpts_per_slab, max_n_tpts_per_slab=max_n_tpts_per_slab,
                 iterate_over_magnus_exp_order=iterate_over_magnus_exp_order,
                 save_log=save_log, file_log=file_log, kwargs=kwargs))
    if P_scan is not NotImplemented:
        return P_scan

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

    And, for :math:`\nu_e \to \nu_{s_2}`:

    >>> oscprob.osc_prob_5nu_vacuum(energy, baseline, s14, s15, s24, s25, s34, s35, d14, d15, d24, d35, D41, D51, nu_i=gd.NUE, nu_f=gd.NUS2)

    and :math:`\nu_{s_1} \to \nu_{s_2}`:

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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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
        verbose=verbose) # L in eV^{-1}

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]

    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # [g cm^{-3}] as a function of radial distance, r, using the Preliminary Reference Earth Model 
    # (PREM). The function matter.num_density_e_func converts the matter density into electron 
    # number density [eV^3].

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_std_potential(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L, # [eV^{-1}]
        t_breakpoints=t_breakpoints,
        osc_params={'sth': sth, 'Dm2': Dm2},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
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
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
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
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
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
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_earth(
    H_func: Callable,
    energy: Union[int, float, list, np.ndarray],
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None,
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None,
    L: Optional[Union[float, list, np.ndarray]]=None,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='trapezoid',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    validate_input: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the neutrino oscillation probability inside
    the Earth for a given arbitrary Hamiltonian.

    Does **not** assume standard oscillations nor a given number of
    neutrino flavors: the user supplies their own Hamiltonian function,
    ``H_func``, and this routine takes care of the geometry of the
    trajectory through the Earth and of the matter density along it.

    ``H_func`` must be a function of either three arguments,
    ``H_func(energy, l, VCC)``, or two arguments,
    ``H_func(energy, l)``, returning a square complex NumPy array (the
    Hamiltonian in the flavor basis, in eV).  In the three-argument
    form, ``VCC`` is the charged-current matter potential
    :math:`V_{\rm CC} = \sqrt{2} G_F N_e` [eV] at position ``l``
    along the chord, computed from the Preliminary Reference Earth
    Model; its sign is already flipped for antineutrinos
    (``nubar=True``).  The user is free to use it, scale it, or ignore
    it (e.g., to add non-standard matter potentials that affect flavors
    other than :math:`\nu_e`).  For extra speed, ``H_func`` may accept
    an array of positions ``l`` and return a stack of Hamiltonians with
    the position axis leading; this is detected automatically.

    The trajectory can be specified either by the cosine of the zenith
    angle (``costhz``) together with the baseline ``L`` [eV^{-1}], or
    by an initial and a final location on the surface of the Earth
    (``loc_ini``, ``loc_fin``), given as (degree, minute, second)
    latitude/longitude tuples or as the names of predefined locations
    (see ``earth.loc_coords_dms``); in the latter case the neutrino
    travels the chord that joins the two locations.

    The slab edges used internally are aligned with the crossings of
    the PREM layer boundaries along the chord.

    Examples
    --------
    Standard three-neutrino oscillations, written by hand (the
    dedicated wrapper :func:`osc_prob_3nu_earth` does this internally):

    >>> h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent( \
    ...     s12, s23, s13, dCP, D21, D31)
    >>> def H(energy, l, VCC):
    ...     return (1/energy)*h_vac + VCC*np.diag([1.0, 0.0, 0.0])
    >>> osc_prob_earth(H, energy=1.e9, loc_ini='fermilab', \
    ...     loc_fin='homestake')
    """
    source_func_name = sys._getframe().f_code.co_name

    # If the location is given as a string, look it up among the predefined named locations
    if isinstance(loc_ini, str):
        loc_ini = earth.coordinates_of_named_location(source_func_name, loc_name=loc_ini)
    if isinstance(loc_fin, str):
        loc_fin = earth.coordinates_of_named_location(source_func_name, loc_name=loc_fin)

    # Resolve the trajectory: either the chord between two surface locations, or (costhz, L)
    costhz, L = validate_input_osc_prob_earth(source_func_name, loc_ini, loc_fin, costhz, L,
        verbose=verbose)

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]

    # Charged-current potential along the chord from the PREM electron density; the antineutrino
    # sign flip is applied inside matter.vcc_func_from_rho_func.  The profile evaluations are
    # cached on repeated position grids.
    VCC_func = matter.vcc_func_from_rho_func(
        rho_func=lambda l: matter.num_density_e_func(
            earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM),
            earth.density_matter_func_prem,
            ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
            electron_fraction=electron_fraction,
            density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        nubar=nubar,
        density_is_of_number_of_electrons=True) # [eV]
    VCC_func = _PositionProfileCache(VCC_func)

    return _osc_prob_with_potential(source_func_name, H_func, VCC_func, energy, L, 0.0, nu_i,
        nu_f, t_breakpoints, magnus_exp_order, n_jobs, integration_method, rtol, atol,
        validate_input, verbose, **kwargs)


def _osc_prob_with_potential(
    source_func_name: str,
    H_func: Callable,
    VCC_func: Callable,
    energy: Union[int, float, list, np.ndarray],
    L: Union[int, float, list, np.ndarray],
    L0: Union[int, float],
    nu_i: Optional[int],
    nu_f: Optional[int],
    t_breakpoints: Optional[np.ndarray],
    magnus_exp_order: int,
    n_jobs: int,
    integration_method: str,
    rtol: Optional[Union[int, float]],
    atol: Optional[Union[int, float]],
    validate_input: bool,
    verbose: int,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Common machinery of :func:`osc_prob_earth` and
    :func:`osc_prob_sun`: wire a user-supplied Hamiltonian function --
    H_func(energy, l, VCC) or H_func(energy, l) -- to the environment
    potential ``VCC_func`` and hand it to
    :func:`osc_prob_energy_baseline`."""

    if validate_input:
        try:
            if not isinstance(H_func, Callable):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": H_func must be a function of (energy, l, VCC) or of (energy, l).")
            n_params_H = len(signature(H_func).parameters)
            if n_params_H not in (2, 3):
                raise ValueError(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + \
                    ": H_func must be a function of either three arguments (energy, l, VCC) or" + \
                    " two arguments (energy, l); the provided H_func takes " + \
                    str(n_params_H) + " argument(s).")
        except ValueError as error:
            print(error)
            print("Aborting execution...")
            sys.exit(1)

    n_params_H = len(signature(H_func).parameters)
    if n_params_H == 3:
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            return H_func(enu, l, VCC_func(l))
    else:
        def htot(enu: Union[int, float], l: Union[int, float, np.ndarray]) -> np.ndarray:
            return H_func(enu, l)

    return osc_prob_energy_baseline(htot, energy, L, L0, nu_i, nu_f, False,
        t_breakpoints=t_breakpoints, magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
        integration_method=integration_method, rtol=rtol, atol=atol,
        validate_input=validate_input, verbose=verbose, **kwargs)


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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


def osc_prob_sun(
    H_func: Callable,
    energy: Union[int, float, list, np.ndarray],
    L: Union[float, list, np.ndarray],
    L0: Optional[Union[int, float]]=0.0,
    nubar: Optional[bool]=False,
    nu_i: Optional[int]=None,
    nu_f: Optional[int]=None,
    magnus_exp_order: Optional[int]=4,
    n_jobs: Optional[int]=1,
    integration_method: Optional[str]='trapezoid',
    rtol: Optional[Union[int, float]]=1.e-3,
    atol: Optional[Union[int, float]]=1.e-3,
    validate_input: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the neutrino oscillation probability inside
    the Sun for a given arbitrary Hamiltonian.

    Does **not** assume standard oscillations nor a given number of
    neutrino flavors: the user supplies their own Hamiltonian function,
    ``H_func``, and this routine provides the solar electron density
    along the (radial) trajectory.

    ``H_func`` must be a function of either three arguments,
    ``H_func(energy, l, VCC)``, or two arguments,
    ``H_func(energy, l)``, returning a square complex NumPy array (the
    Hamiltonian in the flavor basis, in eV).  In the three-argument
    form, ``VCC`` is the charged-current matter potential
    :math:`V_{\rm CC} = \sqrt{2} G_F N_e` [eV] at radial position
    ``l``; its sign is already flipped for antineutrinos
    (``nubar=True``).  For extra speed, ``H_func`` may accept an array
    of positions ``l`` and return a stack of Hamiltonians with the
    position axis leading; this is detected automatically.

    The neutrino travels radially outward from ``L0`` to ``L`` (both in
    eV^{-1}, measured from the center of the Sun).

    For the electron density inside the Sun, it assumes an exponentially
    falling density profile: :math:`N_e(r) = N_e(0) \exp(-r/r_0)`,
    with :math:`N_e(0) = 245 N_\text{Av}~\text{cm}^{-3}` and
    :math:`r_0 = R_\odot/10.54`.  See Eq. (10.62) in
    `Fundamentals of Neutrino Physics and Astrophysics
    <https://academic.oup.com/book/3490>`_ by Carlo Giunti and Chung
    Wook Kim.

    Examples
    --------
    Standard two-neutrino oscillations, written by hand (the dedicated
    wrapper :func:`osc_prob_2nu_sun` does this internally):

    >>> h_vac = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    >>> def H(energy, l, VCC):
    ...     return (1/energy)*h_vac + VCC*np.diag([1.0, 0.0])
    >>> osc_prob_sun(H, energy=1.e7, L=0.9*gd.SUN_RADIUS*gd.UNIT_KM)
    """
    source_func_name = sys._getframe().f_code.co_name

    # Solar electron number density [eV^3] along the radial trajectory; the antineutrino sign
    # flip of the potential is applied inside matter.vcc_func_from_rho_func.  The profile
    # evaluations are cached on repeated position grids.
    VCC_func = matter.vcc_func_from_rho_func(
        rho_func=lambda l: matter.density_matter_func_exp(l, gd.NUM_DENSITY_E_SUN_CENTRAL,
            gd.L_SCALE_SUN), # [eV^3] (l in eV^{-1})
        L0=L0,
        nubar=nubar,
        density_is_of_number_of_electrons=True) # [eV]
    VCC_func = _PositionProfileCache(VCC_func)

    return _osc_prob_with_potential(source_func_name, H_func, VCC_func, energy, L, L0, nu_i,
        nu_f, None, magnus_exp_order, n_jobs, integration_method, rtol, atol,
        validate_input, verbose, **kwargs)


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
        nubar=nubar,
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
        nubar=nubar,
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
        nubar=nubar,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the two-neutrino oscillation probability in
    matter with an exponentially falling density profile, including
    non-standard interactions (NSI).
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability in
    matter with an exponentially falling density profile, including
    non-standard interactions (NSI).
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino (3+1) oscillation 
    probability in matter with an exponentially falling density profile,
    including non-standard interactions (NSI).
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino (3+2) oscillation 
    probability in matter with an exponentially falling density profile,
    including non-standard interactions (NSI).
    """

    try:
        if (rho_central < 0.0 or l_scale <= 0.0):
            raise ValueError(gd.ERROR_MSG_IN_COLOR + \
                " oscprob.osc_prob_5nu_matter_nsi_exp_density: rho_central and l_scale must be " + \
                "non-negative.")
    except ValueError as error:
        print(error)
        print("Aborting execution...")
        sys.exit(1)

    return osc_prob_matter_nsi(
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )


#-----------------------------------------------------------------------
# In matter, NSI, in the Earth
#-----------------------------------------------------------------------

def osc_prob_2nu_earth_nsi(
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    eps_aa: Optional[Union[int, float]]=0.0,
    eps_ab: Optional[Union[int, float]]=0.0,
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
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
    r"""Compute and return the two-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'sth': sth, 'Dm2': Dm2},
        nsi_params={'eps_aa': eps_aa, 'eps_ab': eps_ab},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_earth_nsi(
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=3,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_mm': eps_mm,
            'eps_mt': eps_mt, 'eps_tt': eps_tt},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_earth_nsi(
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=4,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es': eps_es, 
            'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms': eps_ms, 'eps_tt': eps_tt,
            'eps_ts': eps_ts, 'eps_ss': eps_ss},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_earth_nsi(
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, including
    non-standard interactions (NSI).
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_matter_nsi(
        num_flavors=5,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        nsi_params={'eps_ee': eps_ee, 'eps_em': eps_em, 'eps_et': eps_et, 'eps_es1': eps_es1,
            'eps_es2': eps_es2, 'eps_mm': eps_mm, 'eps_mt': eps_mt, 'eps_ms1': eps_ms1, 
            'eps_ms2': eps_ms2, 'eps_tt': eps_tt, 'eps_ts1': eps_ts1, 'eps_ts2': eps_ts2, 
            'eps_s1s1': eps_s1s1, 'eps_s1s2': eps_s1s2, 'eps_s2s2': eps_s2s2},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    r"""Compute and return the five-neutrino oscillation probability in
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
        nubar=nubar,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


#-----------------------------------------------------------------------
# In matter, LIV, in the Earth
#-----------------------------------------------------------------------

def osc_prob_2nu_earth_liv(
    energy: Union[int, float, list, np.ndarray], 
    sth: Union[int, float], 
    Dm2: Union[int, float], 
    sxi: Optional[Union[int, float]]=0.0,
    b1: Optional[Union[int, float]]=0.0,
    b2: Optional[Union[int, float]]=0.0,
    Lambda: Optional[Union[int, float]]=1.0,
    n_liv: Optional[int]=0,
    costhz: Optional[Union[int, float]]=None,
    loc_ini: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    loc_fin: Optional[Union[Tuple[float, float], list, np.ndarray, str]]=None, 
    L: Optional[Union[float, list, np.ndarray]]=None,
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
    r"""Compute and return the two-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # If any of the flavor indices is > 1, fix it (read the docstring above).
    nu_i, nu_f = valid_flavor_indices_2nu(nu_i, nu_f)

    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=2,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'sth': sth, 'Dm2': Dm2},
        liv_params={'sxi': sxi, 'b1': b1, 'b2': b2, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_3nu_earth_liv(
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
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the three-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=3,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 'D21': D21, 'D31': D31},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxiCP': dxiCP, 'b1': b1, 
            'b2': b2, 'b3': b3, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_4nu_earth_liv(
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
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the four-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=4,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's24': s24, 'd24': d24, 's34': s34, 'D21': D21, 'D31': D31, 'D41': D41},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi24': sxi24, 'dxi24': dxi24, 'sxi34': sxi34, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
        default_osc_params_set_name=default_osc_params_set_name,
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
        verbose=verbose,
        **kwargs
    )  


def osc_prob_5nu_earth_liv(
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
    default_osc_params_set_name: Optional[str]='OSC_PARAMS_DEFAULT',
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
    verbose: Optional[int]=0,
    **kwargs
) -> Union[float, np.ndarray]:
    r"""Compute and return the five-neutrino oscillation probability 
    inside the Earth, either between two locations on the surface of the
    Earth, or between the surface and a point in the interior, under
    (one form of) Lorentz-invariance violation.
    
    For the matter density inside the Earth, it uses the Preliminary 
    Reference Earth Model.

    If the initial location (``loc_ini``) and final location 
    (``loc_fin``) on the surface of the Earth are given (i.e., if they
    are not ``None``), then the neutrino travels the chord joining them 
    through the Earth, overriding any given value of costhz given, and 
    using the chord length as the baseline. 

    The initial and final location can be given as a three-entry tuple
    of coordinates in the (degree, minute, second) format.  Alternatively,
    any of the two locations can be given as a predefined named 
    location.  The predefined locations are in the earth.loc_coords_dms
    dictionary:

    >>> import magnus.earth as earth
    >>> list(earth.loc_coords_dms.keys())
    >>> print(earth.loc_coords_dms['fermilab'])

    See the example below.

    [If only a single location is given (i.e., if either ``loc_ini`` or
    ``loc_fin`` are ``None``), the function throws an exception.]

    If neither of the two locations is given, the function uses the 
    given value of ``costhz`` as direction and of ``L`` as baseline.
    (And ``L`` can be an array of baselines.)

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

    # Align the slab edges with the crossings of the PREM layer boundaries along the chord: the
    # matter density is discontinuous there, and the high-order quadrature of the Magnus kernel
    # converges at its nominal order only if the Hamiltonian is smooth inside each slab.
    t_breakpoints = earth.prem_layer_edges_along_chord(costhz)*gd.UNIT_KM # [eV^{-1}]
    
    # The function earth.density_matter_func_prem returns the internal matter density of the Earth
    # as a function of radial distance, r, using the Preliminary Reference Earth Model (PREM). The
    # function matter.num_density_e_func converts the matter density into electron number density.

    # The function earth.earth_radial_distance_from_depth returns the radial distance, measured from
    # the center of the Earth, given a neutrino direction (cosine of zenith angle, costhz) and the 
    # distance of the neutrino, or depth (l), measured from the surface of the Earth.

    return osc_prob_liv(
        num_flavors=5,
        rho_func=lambda l: matter.num_density_e_func(earth.earth_radial_distance_from_depth(costhz, 
            l/gd.UNIT_KM), earth.density_matter_func_prem, ratio_number_neutrons_to_protons=1.0,
            electron_fraction=0.5, density_matter_is_in_g_per_cm3=True), # [eV^3] (l in eV^{-1})
        energy=energy,
        L=L,
        t_breakpoints=t_breakpoints,
        osc_params={'s12': s12, 's23': s23, 's13': s13, 'dCP': dCP, 's14': s14, 'd14': d14, 
            's15': s15, 'd15': d15, 's24': s24, 'd24': d24, 's25': s25, 's34': s34, 's35': s35, 
            'd35': d35, 'D21': D21, 'D31': D31, 'D41': D41, 'D51': D51},
        liv_params={'sxi12': sxi12, 'sxi23': sxi23, 'sxi13': sxi13, 'dxi13': dxi13, 'sxi14': sxi14,
            'dxi14': dxi14, 'sxi15': sxi15, 'dxi15': dxi15, 'sxi24': sxi24, 'dxi24': dxi24, 
            'sxi25': sxi25, 'sxi34': sxi34, 'sxi35': sxi35, 'dxi35': dxi35, 'b1': b1, 'b2': b2, 
            'b3': b3, 'b4': b4, 'b5': b5, 'Lambda': Lambda, 'n_liv': n_liv},
        L0=0.0,
        nubar=nubar,
        nu_i=nu_i,
        nu_f=nu_f,
        density_is_of_number_of_electrons=True,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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
    validate_input: Optional[bool]=True, 
    save_log: Optional[bool]=False, 
    filename_log: Optional[str]='./out.log',
    file_log: Optional[TextIOWrapper]=None, 
    close_file_log_upon_exit: Optional[bool]=True,
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
        validate_input=validate_input,
        save_log=save_log,
        filename_log=filename_log,
        file_log=file_log,
        close_file_log_upon_exit=close_file_log_upon_exit,
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


    # # Five-neutrino oscillations in Earth
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.1, 0.01 # [eV^2]
    # ###
    # L = earth.distance_traveled_inside_earth(-0.05)
    # print(osc_prob_5nu_earth(energy, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51, costhz=-0.05, L=L*gd.UNIT_KM, verbose=0))
    # print(osc_prob_5nu_earth(energy, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
    #     D41=D41, D51=D51, costhz=-1, L=L*gd.UNIT_KM, verbose=0))
    # ###
    # loc_fin = 'fermilab'
    # for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    #     print(loc_ini)
    #     print(osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #          s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
    #     print(osc_prob_5nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24,
    #         d35=d35, D41=D41, D51=D51,
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))


    # # Two-neutrino oscillations in Earth, NSI
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # sth = gd.S12_NO_BF_NUFIT_6_0
    # Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    # eps_aa, eps_ab = 0.05, 0.02
    # ###
    # # costhz = -1.0
    # # L = earth.distance_traveled_inside_earth(costhz)
    # # print(L)
    # # print(osc_prob_2nu_vacuum(energy, L*gd.UNIT_KM, sth, Dm2, verbose=0))
    # # print(osc_prob_2nu_earth(energy, sth, Dm2, 
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, rtol=1.e-3, atol=1.e-3, n_jobs=1))
    # # print(osc_prob_2nu_earth_nsi(energy, sth, Dm2, eps_aa=0, eps_ab=0, 
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0))
    # # print(osc_prob_2nu_earth_nsi(energy, sth, Dm2, eps_aa=eps_aa, eps_ab=eps_ab, 
    # #     costhz=-1, L=L*gd.UNIT_KM, rtol=1.e-3, atol=1.e-3, verbose=0))
    # ###
    # loc_fin = 'fermilab'
    # for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    #     print(loc_ini)
    #     print(osc_prob_2nu_earth(energy, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
    #     print(osc_prob_2nu_earth_nsi(energy, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         eps_aa=eps_aa, eps_ab=eps_ab, loc_ini=loc_ini, loc_fin=loc_fin, verbose=0))
    #     print()


    # # Three-neutrino oscillations in Earth, NSI
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = 0.1, 0, 0, 0, 0, 0# 0.2, 1.1, 0.4, 0.6, 0.5
    # ###
    # costhz = -1.0
    # L = earth.distance_traveled_inside_earth(costhz)
    # print(L)
    # print("3nu, std: " + str(osc_prob_3nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=10, magnus_exp_order=3)))
    # print("3nu, nsi: " + str(osc_prob_3nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_mm=eps_mm, eps_mt=eps_mt, 
    #     eps_tt=eps_tt, costhz=costhz, L=L*gd.UNIT_KM, verbose=0, magnus_exp_order=3, n_jobs=10)))
    # print()
    # print("3nu, std: " + str(osc_prob_3nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=10, magnus_exp_order=4)))
    # print("3nu, nsi: " + str(osc_prob_3nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_mm=eps_mm, eps_mt=eps_mt, 
    #     eps_tt=eps_tt, costhz=costhz, L=L*gd.UNIT_KM, verbose=0, magnus_exp_order=4, n_jobs=10)))
    # ###
    # # loc_fin = 'fermilab'
    # # for loc_ini in ['South Pole']: #['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    # #     print(loc_ini)
    # #     print("3nu, std: " + str(osc_prob_3nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, n_jobs=10, max_magnus_exp_order=5,
    # #         iterate_over_magnus_exp_order=True)))
    # #     print("3nu, nsi: " + str(osc_prob_3nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #         eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_mm=eps_mm, eps_mt=eps_mt, 
    # #         eps_tt=eps_tt, loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, max_magnus_exp_order=5,
    # #         iterate_over_magnus_exp_order=True, n_jobs=10)))
    # #     print()

    # a = np.array([[1., 2.], [3., 4.]])
    # b = np.array([[5., 6.], [7., 8.]])
    # c = np.array([[9., 10.], [11., 12.]])
    # d = np.array([[13., 14.], [15., 16.]])
    # e = np.array([[17., 18.], [19., 20.]])
    # f = np.array([[21., 22.], [23., 24.]])
    # X = np.array([a,b,c])
    # Y = np.array([d,e,f])

    # def commutator(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    #     return X @ Y - Y @ X

    # Z = commutator(X, Y)
    # print(Z)

    # Z = np.stack([commutator(X[i], Y[i]) for i in range(3)], axis=0)
    # print(Z)

    # quit()

    # # Four-neutrino oscillations in Earth, NSI
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s24, s34 = 0.01, 0.02, 0.0
    # d14, d24 = np.radians(10.0), np.radians(100.0)
    # D41 = 0.02 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss \
    #     = 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.02, 0.03
    # print("4nu")
    # ###
    # costhz = -0.20#-1.0
    # L = earth.distance_traveled_inside_earth(costhz)
    # print(L)
    # print("4nu, std: " + str(osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=1, magnus_exp_order=3)))
    # # print("4nu, nsi: " + str(osc_prob_4nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    # #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_es=eps_es, eps_mm=eps_mm, eps_mt=eps_mt, 
    # #     eps_ms=eps_ms, eps_tt=eps_tt, eps_ts=eps_ts, eps_ss=eps_ss,
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=12, magnus_exp_order=3)))
    # # print()
    # # print("4nu, std: " + str(osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=10, magnus_exp_order=4)))
    # # print("4nu, nsi: " + str(osc_prob_4nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    # #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_es=eps_es, eps_mm=eps_mm, eps_mt=eps_mt, 
    # #     eps_ms=eps_ms, eps_tt=eps_tt, eps_ts=eps_ts, eps_ss=eps_ss,
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=10, magnus_exp_order=4)))
    # print()
    # ###
    # loc_fin = 'fermilab'
    # import time
    # for integration_method in ['trapezoid', 'simpson']:
    #     print(integration_method)
    #     for loc_ini in ['SNOLAB', 'Homestake', 'CERN']:#, "South Pole"]:
    #         print(loc_ini)
    #         start = time.time()
    #         for i in range(5):
    #             osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #                 s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #                 loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, n_jobs=1, max_magnus_exp_order=3,
    #                 integration_method=integration_method)
    #             # print("4nu, std: " + str(osc_prob_4nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #             #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #             #     loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, n_jobs=12, max_magnus_exp_order=3,
    #             #     integration_method=integration_method)))
    #         # print("4nu, nsi: " + str(osc_prob_4nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #         #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_es=eps_es, eps_mm=eps_mm,
    #         #     eps_mt=eps_mt, eps_ms=eps_ms, eps_tt=eps_tt, eps_ts=eps_ts, eps_ss=eps_ss,
    #         #     loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, max_magnus_exp_order=3, n_jobs=12,
    #         #     integration_method=integration_method)))
    #         print((time.time()-start)/5)
    #         print()


    # # Five-neutrino oscillations in Earth, NSI
    # np.set_printoptions(precision=3)
    # energy = 10.*gd.UNIT_MEV # [eV]
    # s14, s15, s24, s25, s34, s35 = 0.01, 0.0, 0.0, 0.02, 0.03, 0.0
    # d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    # D41, D51 = 0.01, 0.02 # [eV^2]
    # eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, \
    #     eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2 \
    #     = 0.0, 0.02, 0.0, 0.01, 0.0, 0.0, 0.0, 0.02, 0.02, 0.0, 0.01, 0.01, 0.02, 0.0, 0.01
    # print("5nu")
    # ###
    # costhz = -0.20
    # L = earth.distance_traveled_inside_earth(costhz)
    # print(L)
    # print("5nu, std: " + str(osc_prob_5nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, 
    #     d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51,
    #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=12, magnus_exp_order=3)))
    # print("5nu, nsi: " + str(osc_prob_5nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #     s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
    #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_es1=eps_es1, eps_es2=eps_es2, 
    #     eps_mm=eps_mm, eps_mt=eps_mt, eps_ms1=eps_ms1, eps_ms2=eps_ms2, eps_tt=eps_tt, 
    #     eps_ts1=eps_ts1, eps_ts2=eps_ts2, eps_s1s1=eps_s1s1, eps_s1s2=eps_s1s2, eps_s2s2=eps_s2s2,
    #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=12, magnus_exp_order=3)))
    # # print()
    # # print("5nu, std: " + str(osc_prob_5nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, 
    # #     d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51,
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=10, magnus_exp_order=4)))
    # # print("5nu, nsi: " + str(osc_prob_5nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    # #     s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, 
    # #     d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51,
    # #     eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_es1=eps_es1, eps_es2=eps_es2, 
    # #     eps_mm=eps_mm, eps_mt=eps_mt, eps_ms1=eps_ms1, eps_ms2=eps_ms2, eps_tt=eps_tt, 
    # #     eps_ts1=eps_ts1, eps_ts2=eps_ts2, eps_s1s1=eps_s1s1, eps_s1s2=eps_s1s2, eps_s2s2=eps_s2s2,
    # #     costhz=costhz, L=L*gd.UNIT_KM, verbose=0, n_jobs=10, magnus_exp_order=4)))
    # print()
    # ###
    # loc_fin = 'fermilab'
    # for loc_ini in ['SNOLAB', 'Homestake', 'CERN', "South Pole"]:
    #     print(loc_ini)
    #     print("5nu, std: " + str(osc_prob_5nu_earth(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, 
    #         d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51,
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, n_jobs=12, max_magnus_exp_order=3)))
    #     print("5nu, nsi: " + str(osc_prob_5nu_earth_nsi(energy, nu_i=gd.NUE, nu_f=gd.NUMU, 
    #         s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, 
    #         d14=d14, d15=d15, d24=d24, d35=d35, D41=D41, D51=D51,
    #         eps_ee=eps_ee, eps_em=eps_em, eps_et=eps_et, eps_es1=eps_es1, eps_es2=eps_es2, 
    #         eps_mm=eps_mm, eps_mt=eps_mt, eps_ms1=eps_ms1, eps_ms2=eps_ms2, eps_tt=eps_tt, 
    #         eps_ts1=eps_ts1, eps_ts2=eps_ts2, eps_s1s1=eps_s1s1, eps_s1s2=eps_s1s2, 
    #         eps_s2s2=eps_s2s2,
    #         loc_ini=loc_ini, loc_fin=loc_fin, verbose=0, max_magnus_exp_order=3, n_jobs=12)))
    #     print()

    # LIV in the Earth
    np.set_printoptions(precision=3)
    energy = 1.*gd.UNIT_GEV # [eV]
    costhz = -1.0
    L = earth.distance_traveled_inside_earth(costhz)
    # Two-neutrino
    sth = gd.S12_NO_BF_NUFIT_6_0
    Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
    sxi = 0.1
    b1 = 2.e-4
    b2 = 1.e-4
    Lambda = 1000.*gd.UNIT_GEV
    n_liv = 2
    print("Two-neutrino:")
    print(osc_prob_2nu_earth_liv(energy, sth, Dm2, 
        costhz=costhz, L=L,
        sxi=sxi, b1=b1, b2=b2, Lambda=Lambda, n_liv=n_liv, verbose=0))
    print(osc_prob_2nu_earth_liv(energy, sth, Dm2, 
        loc_ini='fermilab', loc_fin='homestake',
        sxi=sxi, b1=b1, b2=b2, Lambda=Lambda, n_liv=n_liv, verbose=0))
    print()
    # Three-neutrino
    sxi12, sxi23, sxi13, dxiCP = 0.1, 0.02, 0.03, np.radians(10.0)
    b1, b2, b3 = 1.e-4, 2.e-4, 3.e-4 
    Lambda = 1000.*gd.UNIT_GEV
    n_liv = 2
    print("Three-neutrino:")
    print(osc_prob_3nu_earth_liv(energy, costhz=costhz, L=L,
        sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, dxiCP=dxiCP, 
        b1=b1, b2=b2, b3=b3, Lambda=Lambda, n_liv=n_liv, verbose=0))
    print(osc_prob_3nu_earth_liv(energy, 
        loc_ini='fermilab', loc_fin='homestake',
        sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, dxiCP=dxiCP, 
        b1=b1, b2=b2, b3=b3, Lambda=Lambda, n_liv=n_liv, verbose=0))
    print()
    # Four-neutrino
    s14, s24, s34 = 0.1, 0.2, 0.3
    d14, d24 = np.radians(10.0), np.radians(100.0)
    D41 = 0.1 # [eV^2]
    sxi12, sxi23, sxi13, sxi14, sxi24, sxi34 = 0.1, 0.05, 0.06, 0.01, 0.01, 0.02
    dxi13, dxi14, dxi24 = np.radians([10,20,30])
    b1, b2, b3, b4 = 2.e-8, 1.e-8, 5.e-8, 4.e-8 
    Lambda = 1000.*gd.UNIT_GEV
    n_liv = 2
    print("Four-neutrino:")
    print(osc_prob_4nu_earth_liv(energy, costhz=costhz, L=L,
        s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
        sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, sxi34=sxi34, 
        dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, Lambda=Lambda,
        n_liv=n_liv, verbose=0))
    print(osc_prob_4nu_earth_liv(energy,
        loc_ini='fermilab', loc_fin='homestake',
        s14=s14, s24=s24, s34=s34, d14=d14, d24=d24, D41=D41,
        sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi24=sxi24, sxi34=sxi34, 
        dxi13=dxi13, dxi14=dxi14, dxi24=dxi24, b1=b1, b2=b2, b3=b3, b4=b4, Lambda=Lambda,
        n_liv=n_liv, verbose=0))
    print()
    # Five-neutrino
    s14, s15, s24, s25, s34, s35 = 0.1, 0.1, 0.2, 0.2, 0.3, 0.3
    d14, d15, d24, d35 = np.radians([10.0, 20.0, 30.0, 40.0])
    D41, D51 = 0.1, 0.01 # [eV^2]
    sxi12, sxi23, sxi13, sxi14, sxi15, sxi24, sxi25, sxi34, sxi35 = \
        0.1, 0.05, 0.6, 0.05, 0.01, 0.2, 0.4, 0.01, 0.05
    dxi13, dxi14, dxi15, dxi24, dxi35 = np.radians([10, 20, 30, 40, 50])
    b1, b2, b3, b4, b5 = 1.e-3, 2.e-3, 3.e-3, 4.e-3, 5.e-3
    Lambda = 1000.*gd.UNIT_GEV
    n_liv = 2
    print("Five-neutrino:")
    print(osc_prob_5nu_earth_liv(energy, costhz=costhz, L=L,
        s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
        D41=D41, D51=D51,
        sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25,
        sxi34=sxi34, sxi35=sxi35, 
        dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, 
        b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, Lambda=Lambda,
        n_liv=n_liv, verbose=0))
    print(osc_prob_5nu_earth_liv(energy, 
        loc_ini='fermilab', loc_fin='homestake',
        s14=s14, s15=s15, s24=s24, s25=s25, s34=s34, s35=s35, d14=d14, d15=d15, d24=d24, d35=d35, 
        D41=D41, D51=D51,
        sxi12=sxi12, sxi23=sxi23, sxi13=sxi13, sxi14=sxi14, sxi15=sxi15, sxi24=sxi24, sxi25=sxi25,
        sxi34=sxi34, sxi35=sxi35, 
        dxi13=dxi13, dxi14=dxi14, dxi15=dxi15, dxi24=dxi24, dxi35=dxi35, 
        b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, Lambda=Lambda,
        n_liv=n_liv, verbose=0))
