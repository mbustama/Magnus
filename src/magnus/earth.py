# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""earth.py

Contains helper functions related to the Earth: its internal matter
density and the geometry of neutrino trajectories through it.

Routine listings
----------------

    * density_matter_func_prem - Returns the density inside the Earth
           using the Preliminary Reference Earth Model (PREM)
    * prem_layer_edges_along_chord - Returns the positions at which a
           chord through the Earth crosses the PREM layer boundaries
    * distance_traveled_inside_earth - Returns the chord length for a
           given neutrino direction
    * earth_radial_distance_from_depth - Converts position along a
           chord to radial distance from the center of the Earth
    * dms_to_decimal - Converts (degree, minute, second) coordinates to
           decimal degrees
    * chord_length_inside_earth - Returns the chord length between two
           locations on the surface of the Earth
    * costhz_between_points_on_surface - Returns the zenith angle of
           the chord between two locations on the surface of the Earth
    * coordinates_of_named_location - Returns the coordinates of a
           predefined location (e.g., a neutrino detector site)
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Union

import magnus.globaldefs as gd

# Predefined locations in ISO 6709:
# North latitudes are positive, South latitudes are negative
# East longitudes are positive, West longitudes are negative
loc_coords_dms = {
    'baikal':      {'lat': (51, 45, 54),    'lon': (104, 24, 54)},
    'cern':        {'lat': (46, 14, 1.80),  'lon': (6, 3, 11.40)},
    'desy':        {'lat': (53, 34, 19.79), 'lon': (9, 52, 27.59)},
    'ess':         {'lat': (55, 44, 6),     'lon': (13, 15, 5.04)},
    'fermilab':    {'lat': (41, 49, 55),    'lon': (-88, 15, 26)},
    'gran_sasso':  {'lat': (42, 25, 15.8),  'lon': (13, 30, 58.43)},
    'homestake':   {'lat': (44, 21, 5.76),  'lon': (-103, 45, 4.68)},
    'kamioka':     {'lat': (36, 25, 50.05), 'lon': (137, 18, 41.15)}, # Mozumi mine
    'km3net_arca': {'lat': (36, 16, 0), 'lon': (16, 6, 0)}, 
    'km3net_orca': {'lat': (42, 48, 0), 'lon': (6, 2, 0)}, 
    'north_pole':  {'lat': (90, 0, 0),      'lon': (0, 0, 0)},
    'pyhaasalmi':  {'lat': (63, 39, 31),    'lon': (26, 2, 28)},
    'snolab':      {'lat': (46, 28, 18),    'lon': (-81, 11, 12)},
    'south_pole':  {'lat': (-90, 0, 0),     'lon': (0, 0, 0)},
    'tokai':       {'lat': (36, 27, 59),    'lon': (140, 36, 24)},
}


# PREM layers: inner radial boundary of each shell [km] (the last shell ends
# at the surface, gd.EARTH_RADIUS), and the coefficients (c0, c1, c2, c3) of
# the density polynomial rho(x) = c0 + c1*x + c2*x^2 + c3*x^3, with
# x = r/EARTH_RADIUS, inside each shell (Dziewonski & Anderson 1981).
PREM_BOUNDARIES = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0,
                            6346.6, 6356.0, 6368.0])
_PREM_COEFFS = np.array([
    [13.0885,  0.0,    -8.8381,  0.0],
    [12.5815, -1.2638, -3.6426, -5.5281],
    [ 7.9565, -6.4761,  5.5283, -3.0807],
    [ 5.3197, -1.4836,  0.0,     0.0],
    [11.2494, -8.0298,  0.0,     0.0],
    [ 7.1089, -3.8045,  0.0,     0.0],
    [ 2.6910,  0.6924,  0.0,     0.0],
    [ 2.900,   0.0,     0.0,     0.0],
    [ 2.600,   0.0,     0.0,     0.0],
    [ 1.020,   0.0,     0.0,     0.0],
])


def density_matter_func_prem(r: Union[float, np.ndarray],
    tol: Optional[float]=1.e-8) -> Union[float, np.ndarray]:
    r"""Returns the matter density inside the Earth according to the
    Preliminary Reference Earth Model (PREM) [1]_.

    Returns the matter density inside the Earth according to the PREM,
    for a given radial distance measured from the center of the Earth.
    Accepts a single radial distance or an array of radial distances;
    array input is evaluated in a single vectorized pass.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    r : float or np.ndarray
        Radial distance(s) measured from the center of the Earth [km].
    tol : float, optional
        Relative tolerance by which a radial distance may exceed
        ``globaldefs.EARTH_RADIUS`` before a ValueError is raised;
        radii within the tolerance are clamped onto the surface.
        Default: 1e-8.

    Returns
    -------
    float or np.ndarray
        Matter density [:math:`\text{g cm}^{-3}`].

    Raises
    ------
    ValueError
        If any radial distance exceeds globaldefs.EARTH_RADIUS by more
        than the relative tolerance tol.

    References
    ----------

    .. [1] Adam M. Dziewonski & Don L. Anderson, "Preliminary Reference
        Earth Model", Physics of the Earth and Planetary Interiors, 25,
        297 (1981).
    
    Examples
    --------
    .. jupyter-execute::

        from magnus import earth

        for r in (0.0, 3000.0, 5000.0, 6371.0):
            print('r = %6.0f km -> %6.2f g/cm^3'
                  % (r, earth.density_matter_func_prem(r)))
"""
    scalar_input = (np.ndim(r) == 0)
    r = np.asarray(r, dtype=float)

    x = r/gd.EARTH_RADIUS

    if np.any(x - 1.0 > tol):
        raise ValueError('Error in magnus: earth.density_matter_func_prem: value of r cannot exceed ' + \
            'globaldefs.EARTH_RADIUS = ' + str(gd.EARTH_RADIUS) + ' km by more than the ' + \
            'desired tolerance of tol = ' + str(tol))

    # Clamp radii within tolerance of the surface onto the surface
    r = np.minimum(r, gd.EARTH_RADIUS)
    x = np.minimum(x, 1.0)

    # Look up the PREM shell of each radius (side='left' reproduces the
    # right-closed bins of the piecewise definition, e.g., r <= 1221.5) and
    # evaluate the density polynomial via Horner's rule.  This is ~10x
    # faster than an np.select over the ten shells.
    c = _PREM_COEFFS[np.searchsorted(PREM_BOUNDARIES, r, side='left')]
    density = c[..., 0] + x*(c[..., 1] + x*(c[..., 2] + x*c[..., 3]))

    return float(density) if scalar_input else density


def distance_traveled_inside_earth(costhz: float) -> float:
    r"""Returns the distance traveled by a neutrino inside the Earth,
    traveling with a cosine of zenith angle costhz.
    
    Returns the length of the path traveled by a neutrino from the 
    surface ot the Earth, through it, until it reaches a detector. The
    direction of the neutrino is parametrized by the zenith angle of the
    neutrino. Assumes that the neutrino detector is on the surface of 
    the Earth, not underground. As a result, the distance is zero for
    all values of costhz > 0.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino.

    Returns
    -------
    float
        Path length inside the Earth [km].
    
    Examples
    --------
    .. jupyter-execute::

        from magnus import earth

        for costhz in (-0.2, -0.5, -1.0):
            print('costhz = %5.2f -> %8.1f km'
                  % (costhz, earth.distance_traveled_inside_earth(costhz)))
"""
    return 0.0 if costhz > 0.0 else -2.0 * gd.EARTH_RADIUS * costhz


def earth_radial_distance_from_depth(costhz: float, l: Union[float, np.ndarray],
    tol: Optional[float]=1.e-8) -> Union[float, np.ndarray]:
    r"""Returns the radial distance measured from the center of the
    Earth to a position inside the Earth, given by costhz and l.

    A neutrino with direction given by the cosine of the zenith angle,
    costhz, travels from l=0 to l=distance_traveled_inside_earth,
    computed below. The routine returns the radial distance to the
    neutrino when its distance from its point of entry into the Earth is
    l.  Accepts a single distance or an array of distances; array input
    is evaluated in a single vectorized pass.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino.
    l : float or np.ndarray
        Distance(s) of the neutrino from its point of entry into the
        Earth [km].
    tol : float, optional
        Absolute tolerance by which ``l`` may exceed the distance
        traveled inside the Earth before a ValueError is raised;
        distances within the tolerance are clamped onto the exit point.
        Default: 1e-8.

    Returns
    -------
    float or np.ndarray
        Radial distance to the neutrino [km].

    Raises
    ------
    ValueError
        If any l exceeds the distance traveled inside the Earth for
        this value of costhz by more than the tolerance tol.
    """
    scalar_input = (np.ndim(l) == 0)
    l = np.asarray(l, dtype=float)

    d = distance_traveled_inside_earth(costhz)

    if np.any(l - d > tol):
        raise ValueError('Error in magnus: earth_radial_distance_from_depth: value of ' + \
                'l cannot be larger than the distance traveled ' + \
                'inside Earth for this value of costhz')

    # Clamp values of l within tolerance of the exit point onto the exit point
    l = np.minimum(l, d)

    r2 = gd.EARTH_RADIUS*gd.EARTH_RADIUS
    r2 = r2 + (d-l)**2
    r2 = r2 + 2.0*gd.EARTH_RADIUS*(d-l)*costhz
    r = np.sqrt(np.abs(r2))

    return float(r) if scalar_input else r


def prem_layer_edges_along_chord(costhz: float) -> np.ndarray:
    r"""Returns the positions along a chord through the Earth at which
    the chord crosses the PREM layer boundaries.

    A neutrino entering the Earth with direction ``costhz`` travels along
    a chord from :math:`l = 0` to
    :math:`l =` :func:`distance_traveled_inside_earth` (``costhz``).
    The matter density along the chord is piecewise-smooth, with
    discontinuities (or kinks) where the chord crosses the boundaries
    between PREM shells.  This routine returns those crossing positions,
    which are useful as mandatory slab edges for the Magnus expansion:
    high-order quadrature converges at its nominal order only if the
    Hamiltonian is smooth within each slab.

    The crossing positions solve :math:`r(l) = r_b` for each boundary
    radius :math:`r_b`, which is a quadratic equation in :math:`l`: with
    :math:`u = d - l` and :math:`d = -2 R \cos\theta_z`, one has

    .. math::

       u^2 + 2 R \cos\theta_z\, u + \left(R^2 - r_b^2\right) = 0 .

    .. versionadded:: 1.0.0

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino (crossings exist
        only for costhz < 0).

    Returns
    -------
    np.ndarray
        Sorted crossing positions l [km], each strictly inside (0, d).
        Empty if the chord crosses no boundary.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        from magnus import earth

        edges = earth.prem_layer_edges_along_chord(-0.8)
        d = earth.distance_traveled_inside_earth(-0.8)

        print('%d crossings; the first three at %s km'
              % (len(edges), np.round(edges[:3], 1)))
        print('symmetric about the midpoint:',
              np.allclose(edges + edges[::-1], d))
"""
    if costhz >= 0.0:
        return np.array([])

    R = gd.EARTH_RADIUS
    d = -2.0*R*costhz                    # chord length [km]
    rmin2 = R*R*(1.0 - costhz*costhz)    # (squared) closest approach to the center

    crossings = []
    for rb in PREM_BOUNDARIES:
        disc = rb*rb - rmin2
        if disc <= 0.0:                  # chord never reaches this depth
            continue
        s = np.sqrt(disc)
        for u in (-R*costhz - s, -R*costhz + s):
            if 0.0 < u < d:
                crossings.append(d - u)

    return np.unique(np.array(sorted(crossings)))


def dms_to_decimal(degrees: float, minutes: float, seconds: float) -> float:
    r"""Converts (degree, minute, second) coordinates to decimal degrees.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    degrees : float
        Degree part of the coordinate.
    minutes : float
        Minute part of the coordinate.
    seconds : float
        Second part of the coordinate.

    Returns
    -------
    float
        Coordinate in decimal degrees.
    """
    return degrees + minutes / 60 + seconds / 3600


def chord_length_inside_earth(lat1_dms: tuple[float, float, float],
    lon1_dms: tuple[float, float, float], lat2_dms: tuple[float, float, float],
    lon2_dms: tuple[float, float, float]) -> float:
    r"""Returns the chord length between two locations on the surface of
    the Earth.

    Computes the straight-line (chord) distance between two locations on the surface of the
    Earth, assumed spherical, using the haversine formula for the central angle between the two
    locations and converting it to a chord length.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    lat1_dms : tuple of float
        Latitude of the first location, as (degrees, minutes, seconds).
    lon1_dms : tuple of float
        Longitude of the first location, as (degrees, minutes, seconds).
    lat2_dms : tuple of float
        Latitude of the second location, as (degrees, minutes, seconds).
    lon2_dms : tuple of float
        Longitude of the second location, as (degrees, minutes, seconds).

    Returns
    -------
    float
        Chord length between the two locations [km].
    
    Examples
    --------
    .. jupyter-execute::

        from magnus import earth

        fermilab = ((41.0, 49.0, 55.0), (-88.0, -15.0, -26.0))
        sanford = ((44.0, 21.0, 12.0), (-103.0, -45.0, -5.0))

        print('Fermilab to Sanford: %.1f km'
              % earth.chord_length_inside_earth(fermilab[0], fermilab[1],
                                                sanford[0], sanford[1]))
"""

    # Convert DMS to decimal degrees
    lat1 = dms_to_decimal(*lat1_dms)
    lon1 = dms_to_decimal(*lon1_dms)
    lat2 = dms_to_decimal(*lat2_dms)
    lon2 = dms_to_decimal(*lon2_dms)

    # Convert decimal degrees to radians
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    # Differences in coordinates
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad

    # Haversine formula to calculate the central angle
    a = np.sin(delta_lat / 2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(delta_lon / 2)**2
    central_angle = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Straight-line distance (chord length)
    distance = 2 * gd.EARTH_RADIUS * np.sin(central_angle / 2)

    return distance


def costhz_between_points_on_surface(lat1_dms: tuple[float, float, float],
    lon1_dms: tuple[float, float, float], lat2_dms: tuple[float, float, float],
    lon2_dms: tuple[float, float, float]) -> float:
    r"""Returns the zenith angle of the chord between two locations on
    the surface of the Earth.

    Computes the cosine of the zenith angle at which a neutrino would need to travel in a
    straight chord through the Earth to reach the second location from the first (e.g., a source
    and a detector both on the surface).  Assumes a spherical Earth and a detector on the
    surface, not underground, so the returned value is always non-positive (an upward- or
    horizontally-traveling neutrino, i.e. costhz > 0, would not cross the Earth's interior at all).

    .. versionadded:: 1.0.0

    Parameters
    ----------
    lat1_dms : tuple of float
        Latitude of the first location, as (degrees, minutes, seconds).
    lon1_dms : tuple of float
        Longitude of the first location, as (degrees, minutes, seconds).
    lat2_dms : tuple of float
        Latitude of the second location, as (degrees, minutes, seconds).
    lon2_dms : tuple of float
        Longitude of the second location, as (degrees, minutes, seconds).

    Returns
    -------
    float
        Cosine of the zenith angle of the chord connecting the two locations.
    """
    chord_length = chord_length_inside_earth(lat1_dms, lon1_dms, lat2_dms, lon2_dms) # [km]

    return -0.5 * chord_length / gd.EARTH_RADIUS


def coordinates_of_named_location(source_func_name: str, loc_name: str) -> np.ndarray:
    r"""Returns the coordinates of a predefined location (e.g., a
    neutrino detector site).

    Looks up ``loc_name`` (case-insensitively, spaces treated as underscores) in the
    ``loc_coords_dms`` dictionary of predefined locations (neutrino telescopes/detector sites and
    a few reference points) and returns its latitude and longitude.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, used only to build a more informative error message if
        ``loc_name`` is not found.
    loc_name : str
        Name of the predefined location (e.g., ``'kamioka'``, ``'south_pole'``). See
        ``earth.loc_coords_dms`` for the full list.

    Returns
    -------
    np.ndarray
        Array ``[lat, lon]``, with ``lat`` and ``lon`` each a (degree, minute, second) tuple.
    """
    # The latitude and longitude are each returned in day-minute-second format, (dd, mm, ss)

    try:
        lat = loc_coords_dms[loc_name.lower().replace(" ", "_")]['lat']
        lon = loc_coords_dms[loc_name.lower().replace(" ", "_")]['lon']
    except KeyError:
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " oscprob." + source_func_name + ": the given name of the" + \
                " location (" + loc_name + ") is not one of the predefined named locations" + \
                " in Magnus.  The available predefined named locations (in" + \
                " earth.loc_coords_dms)" + " are: " + str(list(loc_coords_dms.keys())) + ".")

    return np.array([lat, lon])


__all__ = [
    'loc_coords_dms',
    'PREM_BOUNDARIES',
    'density_matter_func_prem',
    'distance_traveled_inside_earth',
    'earth_radial_distance_from_depth',
    'prem_layer_edges_along_chord',
    'dms_to_decimal',
    'chord_length_inside_earth',
    'costhz_between_points_on_surface',
    'coordinates_of_named_location',
]
