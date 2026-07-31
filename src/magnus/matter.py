# -*- coding: utf-8 -*-
r"""matter.py

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
    * exp_density_profile - Builds an exponential density-profile
           callable tagged for the fast interaction-picture integrator
    * num_density_e_func - Converts a matter density to an electron
           number density
    * VCC_func - Returns the potential for coherent forward electron
           scattering
    * vcc_func_from_rho_func - Builds a VCC function (or constant) from
           a density profile, handling neutrino/antineutrino sign and
           unit conversion
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Callable, Union

import magnus.globaldefs as gd


def density_matter_func_const(l: float, 
    density_matter_const: Optional[float]=gd.DENSITY_MATTER_CRUST_G_PER_CM3) -> float:
    r"""Returns the matter density as a function of position, assuming a 
    constant density. Used for testing purposes.

    Returns the matter density as a function of position, assuming a
    constant density. Used for testing purposes.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated (in this
        case, the profile is uniform, so any value of l returns the same
        constant density).
    density_matter_const : float
        Matter density [:math:`\text{g cm}^{-3}`]

    Returns
    -------
    float
        Matter density [:math:`\text{g cm}^{-3}`]
    """

    return density_matter_const


def density_matter_func_exp(l: float, density_matter_central:float , l_scale: float) -> float:
    r"""Returns the matter density as a function of position, assuming  
    an exponentially decreasing density profile.

    Returns the matter density as a function of position, assuming  
    an exponentially decreasing density profile of the form

    .. math::

       \rho(l) = \rho_0\, e^{-l/l_\text{scale}} ,

    for given values of the central density :math:`\rho_0`
    (``density_matter_central``) and the length scale
    :math:`l_\text{scale}` (``l_scale``).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated.
    density_matter_central : float
        Matter density at the center of the profile (l = 0) [:math:`\text{g cm}^{-3}`]
    l_scale : float
        Length scale of the exponential density decrease.

    Returns
    -------
    float
        Matter density [:math:`\text{g cm}^{-3}`]
    """

    return density_matter_central*np.exp(-l/l_scale)


def exp_density_profile(density_matter_central: float, l_scale: float) -> Callable:
    r"""Builds an exponential density-profile callable tagged for the fast interaction-picture
    integrator.

    Same functional form as :func:`density_matter_func_exp` (curried over ``density_matter_central``
    and ``l_scale`` so it can be passed directly as ``rho_func``), but the returned callable also
    carries an ``l_scale`` attribute and an ``is_exp_density_profile`` marker set to ``True``.
    :func:`magnus.oscprob.osc_prob_matter_std_potential`, :func:`magnus.oscprob.osc_prob_matter_nsi`,
    and :func:`magnus.oscprob.osc_prob_liv` look for this marker (propagated through
    :func:`vcc_func_from_rho_func`) to detect a genuine exponential profile and automatically switch
    to the much faster interaction-picture Magnus integrator (see
    ``_osc_prob_ip_exp_dispatch``), with a transparent fallback to the general
    slab-refinement method whenever the fast method does not converge (e.g., near an MSW resonance).
    A plain lambda with the same functional form would work numerically but, lacking the marker,
    would silently skip the fast path -- always build exponential profiles through this function (or
    ``osc_prob_*_exp_density``/``osc_prob_*_sun*``, which already do) to get the speed-up.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    density_matter_central : float
        Matter density (or electron number density) at the center of the profile (l = 0).
    l_scale : float
        Length scale of the exponential density decrease.

    Returns
    -------
    Callable
        Function of position, l, tagged with ``is_exp_density_profile = True`` and
        ``l_scale = l_scale``.
    """
    def rho_func(l: Union[int, float, np.ndarray]) -> Union[float, np.ndarray]:
        return density_matter_func_exp(l, density_matter_central, l_scale)
    rho_func.is_exp_density_profile = True
    rho_func.l_scale = float(l_scale)
    return rho_func


def num_density_e_func(l: float, density_matter_func: Callable,
    ratio_number_neutrons_to_protons: Optional[float]=1.0,
    electron_fraction: Optional[float]=0.5,
    density_matter_is_in_g_per_cm3: Optional[bool]=False) -> float:
    r"""Converts matter density [:math:`\text{g cm}^{-3}`] to electron number density
    [:math:`\text{eV}^{3}`], for a given matter density profile and position.

    Converts the matter density [:math:`\text{g cm}^{-3}`] to electron number density
    [:math:`\text{eV}^{3}`], for a given matter density profile, density_matter_func,
    and position, l. Matter is assumed to be isoscalar, with the
    fraction of electrons given by electron_fraction.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated.
    density_matter_func : Callable
        Matter density as a function of l [:math:`\text{g cm}^{-3}`] (or, if
        ``density_matter_is_in_g_per_cm3`` is False, already in natural units).
    ratio_number_neutrons_to_protons : float, optional
        Ratio of the number of neutrons to protons in matter, used to compute the average
        nucleon mass. Default: 1.0.
    electron_fraction : float, optional
        Electron fraction. Default: 0.5.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``density_matter_func`` returns the density in :math:`\text{g cm}^{-3}` and it is converted to
        natural units internally; if False, it is assumed to already be in natural units.
        Default: False.

    Returns
    -------
    float
        Number density of electrons [:math:`\text{eV}^{3}`]
    """
    avg_mass_nucleon = (gd.MASS_PROTON+gd.MASS_NEUTRON*ratio_number_neutrons_to_protons) \
                        / (1.0+ratio_number_neutrons_to_protons)

    # num_density_e = density_matter_func(l) * gd.CONV_G_TO_EV \
    #                     / avg_mass_nucleon * electron_fraction \
    #                     / gd.CONV_CM3_TO_INV_EV3 # [eV^3]

    # If the matter density is given in g cm^{-3} (density_matter_is_in_g_per_cm3 == True), convert it
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

    .. versionadded:: 1.0.0

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated.
    num_density_e_func : Callable
        Electron number density as a function of l [:math:`\text{eV}^{3}`].

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

    .. versionadded:: 1.0.0

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
        If True, flip the sign of :math:`V_\text{CC}` (electrons couple to :math:`\nu_e` and
        :math:`\bar{\nu}_e` with opposite-sign weak charge). Default: False.
    density_matter_is_in_g_per_cm3 : bool, optional
        If True, ``rho_func`` returns the matter density in :math:`\text{g cm}^{-3}`; if False, it is assumed to
        already be in natural units. Ignored if ``density_is_of_number_of_electrons`` is True.
        Default: False.
    density_is_of_number_of_electrons : bool, optional
        If True, ``rho_func`` directly returns the electron number density [:math:`\text{eV}^{3}`], skipping the
        matter-density-to-electron-density conversion. Default: False.

    Returns
    -------
    int, float, or Callable
        V_CC [eV], as a function of position if ``rho_func`` is a function, or as a constant
        (evaluated once, at ``L0``) if ``rho_func`` is a constant.
    """
    s = 1.0 if not nubar else -1.0

    # If rho_func is a genuine exponential profile (tagged by exp_density_profile), propagate the
    # tag to the returned VCC_func: every conversion below (unit conversion, electron fraction,
    # sqrt(2)*G_F, the antineutrino sign s) is a plain scalar rescaling of rho_func(l), which leaves
    # the exponential functional form and l_scale unchanged.  Callers (osc_prob_matter_std_potential,
    # osc_prob_matter_nsi, osc_prob_liv) use this tag to detect the profile and switch to the fast
    # interaction-picture integrator.
    l_scale_tag = getattr(rho_func, 'l_scale', None) if isinstance(rho_func, Callable) else None

    def _tag(vcc: Callable) -> Callable:
        if l_scale_tag is not None:
            vcc.is_exp_density_profile = True
            vcc.l_scale = l_scale_tag
        return vcc

    # If the provided rho_func is the matter density (e.g., g cm^{-3}), convert rho_func to a
    # function that returns the electron number density [eV^3]
    if not density_is_of_number_of_electrons:
        # If rho_func is a constant rather than a callable, wrap it in a dummy function so that
        # num_density_e_func always receives a density it can evaluate at a position.
        if isinstance(rho_func, Callable):
            density_matter_func = rho_func
        else:
            def density_matter_func(r):
                return rho_func

        # Number density of electrons [eV^3]
        def num_density_e(l):
            return num_density_e_func(l,
                density_matter_func=density_matter_func,
                ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
                electron_fraction=electron_fraction,
                density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3)
        # Coherent forward potential, VCC [eV]
        if isinstance(rho_func, Callable):
            # Return VCC as a function, since the density is a function
            def vcc(l):
                return s*VCC_func(l, num_density_e_func=num_density_e)
            return _tag(vcc)
        else:
            # Return VCC as a constant, since the density is a constant. Its value when evaluated at
            # L0 is the same at any other l.
            return s*VCC_func(l=L0, num_density_e_func=num_density_e)
    else: # rho_func is directly the electron number density [eV^3]
        if isinstance(rho_func, Callable):
            def vcc(l):
                return s*VCC_func(l, num_density_e_func=rho_func)
            return _tag(vcc)
        else:
            return s*VCC_func(l=L0, num_density_e_func=rho_func)


__all__ = [
    'density_matter_func_const',
    'density_matter_func_exp',
    'exp_density_profile',
    'num_density_e_func',
    'VCC_func',
    'vcc_func_from_rho_func',
]
