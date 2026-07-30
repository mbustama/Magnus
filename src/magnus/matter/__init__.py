# -*- coding: utf-8 -*-
r"""__init__.py

Contains helper functions to compute the oscillation probability in
matter.

This module contains routines to common matter density profiles (e.g.,
constant, exponentially decreasing), electron number density, and
coherent forward scattering potential.

Routine listings
----------------

    * density_matter_func_const - Returns the density for a constant
           matter density profile
    * density_matter_func_exp - Returns the density for an exponentially
           decreasing matter density profile
    * num_density_e_func - Converts a matter density to an electron
           number density
    * VCC_func - Returns the potential for coherent forward electron
           scattering
    * vcc_func_from_rho_func - Builds a VCC function (or constant) from
           a density profile, handling neutrino/antineutrino sign and
           unit conversion
"""

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Callable, Union, Tuple, List, Dict

import magnus.globaldefs as gd


def density_matter_func_const(l: float, 
    density_matter_const: Optional[float]=gd.DENSITY_MATTER_CRUST_G_PER_CM3) -> float:
    r"""Returns the matter density as a function of position, assuming a 
    constant density. Used for testing purposes.

    Returns the matter density as a function of position, assuming a
    constant density. Used for testing purposes.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated (in this
        case, the profile is uniform, so any value of l returns the same
        constant density).
    density_matter_const : float
        Matter density [g cm^{-3}]

    Returns
    -------
    float
        Matter density [g cm^{-3}]
    """

    return density_matter_const


def density_matter_func_exp(l: float, density_matter_central:float , l_scale: float) -> float:
    r"""Returns the matter density as a function of position, assuming  
    an exponentially decreasing density profile.

    Returns the matter density as a function of position, assuming  
    an exponentially decreasing density profile of the form
    rho(l) = density_matter_central*exp(-l/l_scale), for given values
    of density_matter_central and l_scale.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated.
    density_matter_central : float
        Matter density at the center of the profile (l = 0) [g cm^{-3}]
    l_scale : float
        Length scale of the exponential density decrease.

    Returns
    -------
    float
        Matter density [g cm^{-3}]
    """

    return density_matter_central*np.exp(-l/l_scale)


def num_density_e_func(l: float, density_matter_func: Callable, 
    ratio_number_neutrons_to_protons: Optional[float]=1.0,
    electron_fraction: Optional[float]=0.5,
    density_matter_is_in_g_per_cm3=False) -> float:
    r"""Converts matter density [g cm^{-3}] to electron number density
    [eV^3], for a given matter density profile and position.

    Converts the matter density [g cm^{-3}] to electron number density
    [eV^3], for a given matter density profile, density_matter_func,
    and position, l. Matter is assumed to be isoscalar, with the
    fraction of electrons given by electron_fraction.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated.
    density_matter_func : Callable
        Matter density as a function of l [g cm^{-3}] (or, if
        ``density_matter_is_in_g_per_cm3`` is False, already in natural units).
    ratio_number_neutrons_to_protons : float, optional
        Ratio of the number of neutrons to protons in matter, used to compute the average
        nucleon mass. Default: 1.0.
    electron_fraction : float, optional
        Electron fraction. Default: 0.5.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``density_matter_func`` returns the density in g cm^{-3} and it is converted to
        natural units internally; if False, it is assumed to already be in natural units.
        Default: False.

    Returns
    -------
    float
        Number density of electrons [eV^3]
    """
    avg_mass_nucleon = (gd.MASS_PROTON+gd.MASS_NEUTRON*ratio_number_neutrons_to_protons) \
                        / (1.0+ratio_number_neutrons_to_protons)

    # num_density_e = density_matter_func(l) * gd.CONV_G_TO_EV \
    #                     / avg_mass_nucleon * electron_fraction \
    #                     / gd.CONV_CM3_TO_INV_EV3 # [eV^3]

    # If the matter density is given in g cm^{-3} (density_matter_in_g_per_cm3 == True), convert it
    # natural units of eV^4.  Otherwise, it is assumed that the matter density is in natural units
    # already.
    
    return density_matter_func(l) / avg_mass_nucleon * electron_fraction * \
        (gd.UNIT_G_PER_CM3 if density_matter_is_in_g_per_cm3 else 1.0) # num_density_e [eV^3]


def VCC_func(l: float, num_density_e_func: Callable) -> float:
    r"""Computes and returns the coherent forward electron potential, 
    V_CC, at position l, for a given electron number density, 
    num_density_e_func.

    Computes and returns the coherent forward electron potential, V_CC,
    at position l, for a given electron number density profile,
    num_density_e_func.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated.
    num_density_e_func : Callable
        Electron number density as a function of l [eV^3].

    Returns
    -------
    float
        Coherent forward electron potntial, V_CC [eV]
    """

    return gd.SQRT_OF_2 * gd.GF * num_density_e_func(l) # VCC [eV]


def vcc_func_from_rho_func(
    rho_func: Union[Callable, int, float],
    L0: Optional[Union[int, float]]=0.0,
    ratio_number_neutrons_to_protons: Optional[Union[int, float]]=1.0,
    electron_fraction: Optional[Union[int, float]]=0.5,
    nubar: Optional[bool]=False,
    density_matter_is_in_g_per_cm3: Optional[bool]=False,
    density_is_of_number_of_electrons: Optional[bool]=False,
) -> Union[int, float, Callable]:
    r"""Builds a V_CC function (or constant) from a density profile.

    Builds the coherent forward-scattering potential, V_CC, from a matter- or electron-density
    profile, handling the neutrino/antineutrino sign of V_CC and the unit conversion in one place.
    This is the function that :func:`magnus.oscprob.osc_prob_matter_std_potential` and
    its NSI/LIV counterparts call to turn a user- or environment-supplied density (e.g., a
    constant, an exponential profile, or the Earth's PREM profile) into the ``VCC_func`` consumed
    by the ``hamiltonian_*nu_matter_td``/``hamiltonian_*nu_nsi_td`` functions.

    .. versionadded:: 0.10.0

    Parameters
    ----------
    rho_func : Callable or int or float
        Matter density (or, if ``density_is_of_number_of_electrons`` is True, electron number
        density directly), either as a function of position, l, or as a single constant value.
    L0 : int or float, optional
        Reference position at which to evaluate a constant ``rho_func`` (irrelevant when
        ``rho_func`` is a genuine constant, since the returned V_CC is then position-independent
        by construction). Default: 0.0.
    ratio_number_neutrons_to_protons : int or float, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.
    electron_fraction : int or float, optional
        Electron fraction. Default: 0.5.
    nubar : bool, optional
        If True, flip the sign of V_CC (electrons couple to nu_e and nu_e-bar with opposite-sign
        weak charge). Default: False.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``rho_func`` returns the matter density in g cm^{-3}; if False, it is assumed to
        already be in natural units. Ignored if ``density_is_of_number_of_electrons`` is True.
        Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, ``rho_func`` directly returns the electron number density [eV^3], skipping the
        matter-density-to-electron-density conversion. Default: False.

    Returns
    -------
    int, float, or Callable
        V_CC [eV], as a function of position if ``rho_func`` is a function, or as a constant
        (evaluated once, at ``L0``) if ``rho_func`` is a constant.
    """
    s = 1.0 if not nubar else -1.0

    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a 
    # function that returns the electron number density [eV^3]
    if not density_is_of_number_of_electrons: 
        # if isinstance(rho_func, Callable):
        #     density_matter_func = rho_func
        # else:
        #     density_matter_func = lambda r: rho_func # If rho_func is constant, pass a dummy function 
        # Number density of electrons [eV^3]
        num_density_e = lambda l: num_density_e_func(l, 
            # density_matter_func=density_matter_func,
            density_matter_func=rho_func if isinstance(rho_func, Callable) \
                else (lambda r: rho_func), # If rho_func is constant, pass a dummy function 
            ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons, 
            electron_fraction=electron_fraction, 
            density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3) 
        # Coherent forward potential, VCC [eV]
        if isinstance(rho_func, Callable):
            # Return VCC as a function, since the density is a function
            return lambda l: s*VCC_func(l, num_density_e_func=num_density_e)
        else:
            # Return VCC as a constant, since the density is a constant. Its value when evaluated at
            # L0 is the same at any other l.
            return s*VCC_func(l=L0, num_density_e_func=num_density_e)
    else: # rho_func is directly the electron number density [eV^3]
        if isinstance(rho_func, Callable):
            return lambda l: s*VCC_func(l, num_density_e_func=rho_func) 
        else:
            return s*VCC_func(l=L0, num_density_e_func=rho_func) 


if __name__ == "__main__":

    pass


name = 'matter'

__all__ = [s for s in dir() if not s.startswith('_')]
