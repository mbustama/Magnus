# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
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
    * matter_potential_projector - Returns the flavor structure of the
           matter term, diag(1, 0, ..., 0, r/2, ...), with the sterile
           states' neutral-current entry included
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import warnings

import numpy as np
from typing import Optional, Callable, Union

import magnus.globaldefs as gd


class DensityUnitWarning(UserWarning):
    r"""Warns that a matter density declared to be in
    :math:`\text{g cm}^{-3}` is too large to be one, and was most likely
    already converted to natural units.

    The conversion factor is :math:`4.3 \times 10^{18}`, so converting a
    second time inflates the matter potential far beyond anything physical.
    What makes it worth a warning is that the result does not look wrong:
    the matter term dominates every other scale, :math:`\nu_e` becomes an
    exact eigenstate of the Hamiltonian, and the calculation returns a
    perfectly self-consistent :math:`P_{ee} = 1`.  That reads as a broken
    formula rather than as a bad input.

    .. versionadded:: 1.0.0
    """


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
    
    Examples
    --------
    .. jupyter-execute::

        import magnus.globaldefs as gd
        from magnus import matter

        profile = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL,
                                             gd.L_SCALE_SUN)

        for frac in (0.0, 0.1, 0.5):
            print('l = %.1f R_sun -> n_e = %.3e eV^3'
                  % (frac, profile(frac*gd.SUN_RADIUS*gd.UNIT_KM)))
"""
    def rho_func(l: Union[int, float, np.ndarray]) -> Union[float, np.ndarray]:
        return density_matter_func_exp(l, density_matter_central, l_scale)
    rho_func.is_exp_density_profile = True
    rho_func.l_scale = float(l_scale)
    return rho_func


IMPLAUSIBLE_DENSITY_G_PER_CM3 = 1.0e16
r"""float: Module-level constant

Matter density [:math:`\text{g cm}^{-3}`] above which a value declared to be in
g cm^-3 is almost certainly already in natural units, and about to be converted a
second time.

The two scales do not overlap.  The densest matter anyone models is a neutron
star interior at some :math:`10^{15}\ \text{g cm}^{-3}`, while *any* density from
water upwards becomes :math:`4.3 \times 10^{18}` or more once multiplied by
``gd.UNIT_G_PER_CM3`` -- so a converted value re-declared as g cm^-3 lands at
least three orders of magnitude above anything physical.

Double conversion is silent and its consequences do not look like a unit
mistake: it inflates the matter potential by ~18 orders of magnitude, which makes
:math:`\nu_e` an exact eigenstate everywhere and returns a perfectly
self-consistent :math:`P_{ee} = 1`.  That reads as a broken formula, not as a
bad input, so it is worth catching where it happens.

.. versionadded:: 1.0.0
"""


IMPLAUSIBLE_DENSITY_NATURAL_UNITS = 1.0e10
r"""float: Module-level constant

Matter density in natural units below which a value declared *not* to be in
g cm^-3 is almost certainly a g cm^-3 number whose
``density_matter_is_in_g_per_cm3`` flag was left at its default of False.

This is the mirror of :data:`IMPLAUSIBLE_DENSITY_G_PER_CM3`, and it guards the
commoner mistake.  The flag defaults to False, so a density read straight off a
table -- 2.848 for the Earth's crust, 13 for its core, 150 for the Sun's center
-- is taken as already converted unless the caller says otherwise.  Anything
physical is :math:`4.3 \times 10^{18}` or more in natural units, so a table
value lands nine orders of magnitude below this threshold, which itself sits at
:math:`2 \times 10^{-9}\ \text{g cm}^{-3}` -- six orders more tenuous than air,
and far below any medium in which anyone computes oscillations.

Under-conversion is quieter than double conversion.  It does not inflate
anything; it makes the matter potential vanish, and the call returns *exactly*
the vacuum probability.  That is a perfectly ordinary-looking number, in the
right range, of the right shape, and nothing about it suggests matter was left
out -- which is why it is worth catching where it happens.

Deliberate vacuum is not caught: a density of exactly zero is left alone.

.. versionadded:: 1.0.0
"""


def _warn_if_density_is_probably_in_g_per_cm3(
    density: Union[int, float, np.ndarray],
    source_func_name: str
) -> None:
    r"""Warns when a density declared *not* to be in g cm^-3 is too small to be
    anything else.

    Called only when ``density_matter_is_in_g_per_cm3`` is False, so a nonzero
    value below :data:`IMPLAUSIBLE_DENSITY_NATURAL_UNITS` means the caller has
    very likely passed a g cm^-3 number and left the flag at its default.

    A warning rather than an error, for the same reason as its mirror: the
    threshold sits far below anything physical, but that is a statement about
    the media people currently model, not a law.

    .. versionadded:: 1.0.0
    """
    largest = float(np.max(np.abs(np.asarray(density, dtype=float))))
    if largest == 0.0 or largest >= IMPLAUSIBLE_DENSITY_NATURAL_UNITS:
        return

    warnings.warn(gd.WARNING_MSG_NO_COLOR + " matter." + source_func_name + ": a matter density "
        "of " + format(largest, '.3e') + " was declared to be in natural units "
        "(density_matter_is_in_g_per_cm3 is False, its default), but it is far too small to be "
        "one -- anything physical is 4.3e18 or more, since that is what one g cm^-3 becomes.  "
        "It was most likely read off a table in g cm^-3 and the flag left unset, in which case "
        "the matter potential is about to come out ~19 orders of magnitude too small, i.e. zero: "
        "the call will return exactly the vacuum probability, which looks like an ordinary "
        "answer rather than a missing one.  Either pass density_matter_is_in_g_per_cm3=True, or "
        "convert the density yourself (multiply by gd.UNIT_G_PER_CM3).  Shown once per "
        "session.", DensityUnitWarning, stacklevel=3)


def matter_potential_projector(
    num_flavors: int,
    ratio_number_neutrons_to_protons: Optional[Union[int, float, Callable]] = 1.0
) -> Union[np.ndarray, Callable]:
    r"""Returns the matrix that multiplies :math:`V_{\rm CC}` in the matter Hamiltonian.

    :math:`\mathrm{diag}(1, 0, \ldots, 0, r/2, \ldots, r/2)`: one on :math:`\nu_e`, zero on
    the other active flavors, and :math:`r/2` on every sterile state, with
    :math:`r = n_n/n_p`.  When :math:`r` is a function of position, so is the returned
    projector.

    **The sterile entries are the whole reason this function exists.**  In matter the active
    flavors all feel the same neutral-current potential :math:`V_{\rm NC}`, so it is
    proportional to the identity across them and drops out as an overall phase -- which is
    why two- and three-flavor codes can write the matter term as :math:`V_{\rm CC}` on
    :math:`\nu_e` alone.  **A sterile state feels neither current**, so once the common
    :math:`V_{\rm NC}` is removed it is left carrying :math:`-V_{\rm NC}`, and that is
    physical: it is the whole reason a 3+1 scenario is more than a relabeling.  With
    :math:`V_{\rm NC} = -G_F n_n/\sqrt{2}` and :math:`V_{\rm CC} = \sqrt{2} G_F n_e`, and
    :math:`n_e = n_p`,

    .. math::

       -V_{\rm NC} = \frac{r}{2} V_{\rm CC} , \qquad r = \frac{n_n}{n_p} .

    Omitting it is equivalent to giving the sterile state the actives' neutral-current
    potential.  Measured on a PREM chord at :math:`\cos\theta_z = -0.9` with
    :math:`\sin^2\theta_{14} = \sin^2\theta_{24} = 0.1` and
    :math:`\Delta m^2_{41} = 1\,{\rm eV}^2`, against a converged external reference, that
    omission costs **0.29 in probability** -- and it is flat in the requested tolerance, so
    no amount of refinement reveals it.  Restoring the term brings the same comparison to
    2.3e-04.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       ``ratio_number_neutrons_to_protons`` may be a callable of position, in which case
       the projector is returned as a callable of the same position.  This is how the
       Earth entry points now build it by default, so that the sterile entries follow the
       same layer-by-layer :math:`Y_e` the density uses.

    Parameters
    ----------
    num_flavors : int
        Number of neutrino flavors; the first three are active, the rest sterile.
    ratio_number_neutrons_to_protons : int, float, or Callable, optional
        :math:`r = n_n/n_p`.  1.0 (isoscalar matter) by default, matching the default of
        :func:`vcc_func_from_rho_func`, which should be given the same value.

        A callable is read as :math:`r(l)`, a function of the position the density
        profile takes -- pass the same composition to both, since :math:`r` and
        :math:`Y_e` are the same statement about the medium,
        :math:`r = (1 - Y_e)/Y_e`.  It must accept an array of positions and return an
        array of the same shape.  The Earth entry points build exactly this by default,
        with :math:`Y_e` resolved per PREM layer; a *scalar* passed there instead
        describes a single medium the layered density does not, which is worth up to
        ~0.4 in probability at 3+1 near the sterile matter resonance on a core-crossing
        chord -- see :class:`magnus.globaldefs.SterileMatterCompositionWarning`.

    Returns
    -------
    np.ndarray or Callable
        Real ``(num_flavors, num_flavors)`` matrix; or, for a callable :math:`r`, a
        function of position returning that matrix -- shape
        ``(num_flavors, num_flavors)`` at a scalar position, and a stack with the
        position axis leading at an array of positions.  At three flavors or fewer the
        constant matrix is returned even for a callable :math:`r`, since the sterile
        block it would act on is empty.
    """
    proj = np.zeros((num_flavors, num_flavors))
    proj[0][0] = 1.0
    if (num_flavors <= 3) or (not callable(ratio_number_neutrons_to_protons)):
        # Active flavors share V_NC and it cancels; sterile states do not, and keep -V_NC.
        for k in range(3, num_flavors):
            proj[k][k] = 0.5*float(ratio_number_neutrons_to_protons)
        return proj

    # A position-resolved ratio: the active entries stay constant, so only the sterile
    # mask is scaled per position.  Vectorized over a position array, with the position
    # axis leading, matching the convention of every Hamiltonian closure in oscprob.
    ratio_func = ratio_number_neutrons_to_protons
    sterile_mask = np.zeros((num_flavors, num_flavors))
    for k in range(3, num_flavors):
        sterile_mask[k][k] = 0.5

    def proj_of_l(l: Union[int, float, np.ndarray]) -> np.ndarray:
        r = np.asarray(ratio_func(l), dtype=float)
        return proj + r[..., None, None]*sterile_mask
    return proj_of_l


def _density_would_trip_a_unit_guard(
    density: Union[int, float, np.ndarray],
    density_matter_is_in_g_per_cm3: bool,
    density_is_of_number_of_electrons: bool
) -> bool:
    r"""Whether converting this density would emit a :class:`DensityUnitWarning`.

    The two guards above are the only reason :func:`vcc_func_from_rho_func` cannot simply
    return a cached constant: they live inside the conversion that a cache hit skips.  This
    reproduces their *conditions* -- and only their conditions -- so that a density which
    would be warned about is never cached in the first place, and every such call therefore
    takes the route that warns.

    That is the way round it has to be done.  Re-emitting the warning from the cache-hit site
    does not work: ``warnings.warn`` is called there with a ``stacklevel`` that attributes it
    to a different frame, and the frame is part of the interpreter's warning registry key, so
    the imitation registers separately and prints a *second* time under the default filter
    where an uncached call prints once.  Declining to cache costs nothing, because a density
    that trips a guard is a mistake being reported rather than a hot path.

    Kept adjacent to the two guards deliberately: it duplicates their thresholds, so it has to
    be read and changed with them.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    density : int, float or np.ndarray
        The density as the caller supplied it.
    density_matter_is_in_g_per_cm3 : bool
        Which of the two guards applies.
    density_is_of_number_of_electrons : bool
        When True the conversion is skipped entirely and neither guard runs.

    Returns
    -------
    bool
        Whether a warning would be emitted, and so whether caching must be declined.
    """
    if density_is_of_number_of_electrons:
        # rho_func is already an electron number density; no conversion, so no guard.
        return False
    try:
        largest = float(np.max(np.abs(np.asarray(density, dtype=float))))
    except (TypeError, ValueError):     # pragma: no cover -- float() already rejected it
        return True
    if density_matter_is_in_g_per_cm3:
        return largest > IMPLAUSIBLE_DENSITY_G_PER_CM3
    return 0.0 < largest < IMPLAUSIBLE_DENSITY_NATURAL_UNITS


_VCC_CONST_CACHE = {}
r"""dict: Memo for the constant-density branches of :func:`vcc_func_from_rho_func`.

Holds plain floats keyed on the seven scalars that determine them.  The callable branches are
deliberately not cached; see the comment at the lookup.
"""

_VCC_CONST_CACHE_MAX = 256
r"""int: How many entries :data:`_VCC_CONST_CACHE` holds before being cleared wholesale."""


def _remember_const_vcc(key, value):
    r"""Stores a constant V_CC against ``key`` and returns it.

    ``key`` is None when the caller's density was a function, in which case there is nothing to
    remember and the value passes straight through.

    Parameters
    ----------
    key : tuple or None
        Cache key from :func:`vcc_func_from_rho_func`, or None to skip caching.
    value : float
        The constant V_CC [eV].

    Returns
    -------
    float
        ``value``, unchanged.
    """
    if key is not None:
        if len(_VCC_CONST_CACHE) >= _VCC_CONST_CACHE_MAX:
            _VCC_CONST_CACHE.clear()
        _VCC_CONST_CACHE[key] = value
    return value


def _warn_if_density_was_probably_already_converted(
    density: Union[int, float, np.ndarray],
    source_func_name: str
) -> None:
    r"""Warns when a density declared as g cm^-3 is too large to be one.

    Called only when ``density_matter_is_in_g_per_cm3`` is True, so a value above
    :data:`IMPLAUSIBLE_DENSITY_G_PER_CM3` means the caller has very likely
    multiplied by ``gd.UNIT_G_PER_CM3`` already.

    A warning rather than an error: the threshold sits above anything physical,
    but "above anything physical" is a statement about the matter people
    currently model, not a law, and this should not stand in the way of someone
    deliberately exploring past it.

    .. versionadded:: 1.0.0
    """
    # The scalar case is the overwhelmingly common one and does not need numpy: asarray + abs +
    # max cost three dispatches to compare one number, on a guard that runs on every call.
    if type(density) is float or type(density) is int:
        if -IMPLAUSIBLE_DENSITY_G_PER_CM3 <= density <= IMPLAUSIBLE_DENSITY_G_PER_CM3:
            return
        largest = abs(float(density))
    else:
        largest = float(np.max(np.abs(np.asarray(density, dtype=float))))
        if largest <= IMPLAUSIBLE_DENSITY_G_PER_CM3:
            return

    warnings.warn(gd.WARNING_MSG_NO_COLOR + " matter." + source_func_name + ": a matter density "
        "of " + format(largest, '.3e') + " g cm^-3 was declared to be in g cm^-3, which is "
        "far denser than a neutron star (~1e15).  It was most likely already converted to "
        "natural units, in which case it is about to be converted a second time -- inflating "
        "the matter potential by ~18 orders of magnitude and returning a self-consistent but "
        "meaningless result (nu_e becomes an exact eigenstate, so P_ee = 1).  Either pass the "
        "density in g cm^-3, or leave density_matter_is_in_g_per_cm3 at False.  Shown once per "
        "session.", DensityUnitWarning, stacklevel=3)


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
    
    density = density_matter_func(l)

    if density_matter_is_in_g_per_cm3:
        _warn_if_density_was_probably_already_converted(density, 'num_density_e_func')
    else:
        _warn_if_density_is_probably_in_g_per_cm3(density, 'num_density_e_func')

    return density / avg_mass_nucleon * electron_fraction * \
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
    
    Examples
    --------
    .. jupyter-execute::

        import magnus.globaldefs as gd
        from magnus import matter

        profile = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL,
                                             gd.L_SCALE_SUN)
        print('V_CC at the center of the Sun: %.3e eV'
              % matter.VCC_func(0.0, profile))
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
    ratio_number_neutrons_to_protons : int, float, or Callable, optional
        Ratio of the number of neutrons to protons in matter. Default: 1.0.  A callable is
        read as :math:`r(l)` and evaluated at each position when converting a matter
        density -- it sets the average nucleon mass, a property of the local composition --
        so a position-resolved composition feeds the potential and the sterile projector
        consistently (see :func:`matter_potential_projector`).  A constant ``rho_func``
        with a callable ratio therefore still yields a position-*dependent* V_CC.  Ignored
        when ``density_is_of_number_of_electrons`` is True, where no conversion happens.

        .. versionchanged:: 1.1.0
           A callable is accepted; it used to have to be a scalar.
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
        V_CC [eV], as a function of position if ``rho_func`` is a function (or if a callable
        ``ratio_number_neutrons_to_protons`` makes the conversion position-dependent), or as
        a constant (evaluated once, at ``L0``) if both are constants.
    """
    s = 1.0 if not nubar else -1.0

    # A constant density makes this a pure function of seven scalars returning one float, and a
    # scan or a fit calls it once per point with the same six.  Only the constant branch is
    # cached: the callable branches return a *closure over rho_func*, which is neither reliably
    # hashable nor safe to hand out twice (callers attach is_exp_density_profile to it), so
    # caching those would trade a few microseconds for an aliasing bug.
    const_key = None
    if not callable(rho_func):
        try:
            const_key = (float(rho_func), float(L0),
                         float(ratio_number_neutrons_to_protons), float(electron_fraction),
                         bool(nubar), bool(density_matter_is_in_g_per_cm3),
                         bool(density_is_of_number_of_electrons))
        except (TypeError, ValueError):
            # An array-valued density is a legitimate input the arithmetic below handles, and
            # float() rejects it.  Not being cacheable is not a reason to refuse the call:
            # fall through uncached rather than raising where main returned an answer.
            const_key = None
        # A density that would be warned about is not cached at all, so the guard keeps firing
        # from the place it has always fired from.  See _density_would_trip_a_unit_guard for
        # why the warning cannot instead be re-emitted here, and note that the *dangerous*
        # case is the one where the flag was left unset: that returns exactly the vacuum
        # probability, so the warning is the only thing distinguishing it from an answer.
        if (const_key is not None) and _density_would_trip_a_unit_guard(
                rho_func, density_matter_is_in_g_per_cm3,
                density_is_of_number_of_electrons):
            const_key = None
        if const_key is not None:
            hit = _VCC_CONST_CACHE.get(const_key)
            if hit is not None:
                return hit

    # If rho_func is a genuine exponential profile (tagged by exp_density_profile), propagate the
    # tag to the returned VCC_func: every conversion below (unit conversion, electron fraction,
    # sqrt(2)*G_F, the antineutrino sign s) is a plain scalar rescaling of rho_func(l), which leaves
    # the exponential functional form and l_scale unchanged.  Callers (osc_prob_matter_std_potential,
    # osc_prob_matter_nsi, osc_prob_liv) use this tag to detect the profile and switch to the fast
    # interaction-picture integrator.
    l_scale_tag = getattr(rho_func, 'l_scale', None) if callable(rho_func) else None

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
        if callable(rho_func):
            density_matter_func = rho_func
        else:
            def density_matter_func(r):
                return rho_func

        # Number density of electrons [eV^3].  A callable ratio is resolved here, at each
        # position, rather than inside num_density_e_func: that function's arithmetic is a
        # plain elementwise formula, so handing it the local value (scalar or array, matching
        # l) keeps it untouched.
        ratio_is_of_l = callable(ratio_number_neutrons_to_protons)
        def num_density_e(l):
            return num_density_e_func(l,
                density_matter_func=density_matter_func,
                ratio_number_neutrons_to_protons=(
                    ratio_number_neutrons_to_protons(l) if ratio_is_of_l
                    else ratio_number_neutrons_to_protons),
                electron_fraction=electron_fraction,
                density_matter_is_in_g_per_cm3=density_matter_is_in_g_per_cm3)
        # Coherent forward potential, VCC [eV].  A callable ratio makes the average nucleon
        # mass -- and so the conversion -- position-dependent even over a constant matter
        # density, which is why it takes the function branch alongside a callable rho_func.
        if callable(rho_func) or ratio_is_of_l:
            # Return VCC as a function, since the density is a function
            def vcc(l):
                return s*VCC_func(l, num_density_e_func=num_density_e)
            return _tag(vcc)
        else:
            # Return VCC as a constant, since the density is a constant. Its value when evaluated at
            # L0 is the same at any other l.
            return _remember_const_vcc(
                const_key, s*VCC_func(l=L0, num_density_e_func=num_density_e))
    else: # rho_func is directly the electron number density [eV^3]
        if callable(rho_func):
            def vcc(l):
                return s*VCC_func(l, num_density_e_func=rho_func)
            return _tag(vcc)
        else:
            # Wrapped for the same reason as the constant density above: VCC_func evaluates
            # what it is given at a position, so handing it the bare number raised
            # TypeError: 'float' object is not callable.  The branch above wraps its
            # constant and this one did not, so a constant electron number density -- the
            # documented way to use density_is_of_number_of_electrons -- could not be used
            # at all.
            def num_density_e_const(l):
                return rho_func

            return _remember_const_vcc(
                const_key, s*VCC_func(l=L0, num_density_e_func=num_density_e_const))


__all__ = [
    'DensityUnitWarning',
    'IMPLAUSIBLE_DENSITY_G_PER_CM3',
    'IMPLAUSIBLE_DENSITY_NATURAL_UNITS',
    'density_matter_func_const',
    'density_matter_func_exp',
    'exp_density_profile',
    # Exported because it is the one definition of the matter term's structure, and
    # because every place that rebuilt that structure by hand instead got it wrong: four
    # inline copies in oscprob, the NSI route's literal diagonal, and notebook 12's
    # solve_ivp reference all gave the sterile states a zero where they carry
    # -V_NC = (r/2) V_CC.  A function that public docstrings point at, and that callers
    # writing their own H_func need, has no business being unexported -- while it was,
    # autoapi did not document it and every :func: reference to it failed to resolve.
    'matter_potential_projector',
    'num_density_e_func',
    'VCC_func',
    'vcc_func_from_rho_func',
]
