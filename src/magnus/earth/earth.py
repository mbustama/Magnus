# -*- coding: utf-8 -*-
r"""Contains helper functions to compute the oscillation probability in
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
    * density_matter_prem - Returns the density inside the Earth using
           the Preliminary Reference Earth Model
    * num_density_e_func - Converts a matter density to an electron
           number density
    * VCC_func - Returns the potential for coherent forward electron
           scattering

Created: 2024/11/30 15:42
Last modified: 2024/11/30 21:23
"""

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np
from typing import Optional, Callable, Union

# TO-DO: remove this once setup.py and pip are working
import os, sys
sys.path.append(os.path.split(os.path.split(os.getcwd())[0])[0])
# sys.path.append('/home/mbustamante/Research/magnus/src/')
# print(os.path.split(os.path.split(os.getcwd())[0])[0])

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


def density_matter_func_prem(r: Union[float, np.ndarray],
    tol: Optional[float]=1.e-8) -> Union[float, np.ndarray]:
    r"""Returns the matter density inside the Earth according to the
    Preliminary Reference Earth Model (PREM).

    Returns the matter density inside the Earth according to the PREM,
    for a given radial distance measured from the center of the Earth.
    Accepts a single radial distance or an array of radial distances;
    array input is evaluated in a single vectorized pass.

    Parameters
    ----------
    r : float or np.ndarray
        Radial distance(s) measured from the center of the Earth [km].

    Returns
    -------
    float or np.ndarray
        Matter density [g cm^{-3}].

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
    """
    scalar_input = (np.ndim(r) == 0)
    r = np.asarray(r, dtype=float)

    x = r/gd.EARTH_RADIUS

    if np.any(x - 1.0 > tol):
        raise ValueError('earth.density_matter_func_prem: value of r cannot exceed ' + \
            'globaldefs.EARTH_RADIUS = ' + str(gd.EARTH_RADIUS) + ' km by more than the ' + \
            'desired tolerance of tol = ' + str(tol))

    # Clamp radii within tolerance of the surface onto the surface
    r = np.minimum(r, gd.EARTH_RADIUS)
    x = np.minimum(x, 1.0)

    density = np.select(
        [r <= 1221.5,
         r <= 3480.0,
         r <= 5701.0,
         r <= 5771.0,
         r <= 5971.0,
         r <= 6151.0,
         r <= 6346.6,
         r <= 6356.0,
         r <= 6368.0,
         r <= gd.EARTH_RADIUS],
        [13.0885-8.8381*x*x,
         12.5815-1.2638*x-3.6426*x*x-5.5281*x*x*x,
         7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x,
         5.3197-1.4836*x,
         11.2494-8.0298*x,
         7.1089-3.8045*x,
         2.6910+0.6924*x,
         np.full_like(x, 2.900),
         np.full_like(x, 2.600),
         np.full_like(x, 1.020)])

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

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino.

    Returns
    -------
    float
        Path length inside the Earth [km].
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

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino.
    l : float or np.ndarray
        Distance(s) of the neutrino from its point of entry into the
        Earth [km].

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
        raise ValueError('earth_radial_distance_from_depth: value of ' + \
                'l cannot be larger than the distance traveled ' + \
                'inside Earth for this value of costhz')

    # Clamp values of l within tolerance of the exit point onto the exit point
    l = np.minimum(l, d)

    r2 = gd.EARTH_RADIUS*gd.EARTH_RADIUS
    r2 = r2 + (d-l)**2
    r2 = r2 + 2.0*gd.EARTH_RADIUS*(d-l)*costhz
    r = np.sqrt(np.abs(r2))

    return float(r) if scalar_input else r


def dms_to_decimal(degrees: float, minutes: float, seconds: float) -> float:
    """
    Convert coordinates from degrees, minutes, and seconds to decimal degrees.
    
    Parameters:
        degrees: The degree part of the coordinate
        minutes: The minute part of the coordinate
        seconds: The second part of the coordinate
    
    Returns:
        Decimal degrees
    """
    return degrees + minutes / 60 + seconds / 3600


def chord_length_inside_earth(lat1_dms: tuple[float, float, float], 
    lon1_dms: tuple[float, float, float], lat2_dms: tuple[float, float, float], 
    lon2_dms: tuple[float, float, float]) -> float:
    """
    Calculate the straight-line distance between two points through the Earth.
    
    Parameters:
        lat1_dms, lon1_dms: Tuple of (degrees, minutes, seconds) for the first point
        lat2_dms, lon2_dms: Tuple of (degrees, minutes, seconds) for the second point
    
    Returns:
        Straight-line distance in kilometers
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

    # Assumes spherical Earth and detector on the surface, not underground.

    chord_length = chord_length_inside_earth(lat1_dms, lon1_dms, lat2_dms, lon2_dms) # [km]

    return -0.5 * chord_length / gd.EARTH_RADIUS


def coordinates_of_named_location(source_func_name: str, loc_name: str) -> np.ndarray:

    # The latitude and longitude are each returned in day-minute-second format, (dd, mm, ss)

    try:
        lat = loc_coords_dms[loc_name.lower().replace(" ", "_")]['lat']
        lon = loc_coords_dms[loc_name.lower().replace(" ", "_")]['lon']
    except KeyError:
        print(gd.ERROR_MSG_IN_COLOR + " oscprob." + source_func_name + ": the given name of the" + \
                " the location (" + loc_name + ") is not one of the predefined named locations" + \
                " in Magnus.  The available predefined named locations (in" + \
                " earth.loc_coords_dms)" + " are: " + str(list(loc_coords_dms.keys())) + ".")
        print("Aborting execution...")
        sys.exit(1)

    return np.array([lat, lon])


if __name__ == "__main__":

    lat1_dms = (52, 31, 12)  # Berlin latitude: 52°31'12"
    lon1_dms = (13, 24, 18)  # Berlin longitude: 13°24'18"
    lat2_dms = (48, 51, 24)  # Paris latitude: 48°51'24"
    lon2_dms = (2, 21, 7)    # Paris longitude: 2°21'7"

    distance = chord_length_inside_earth(lat1_dms, lon1_dms, lat2_dms, lon2_dms)
    costhz = costhz_between_points_on_surface(lat1_dms, lon1_dms, lat2_dms, lon2_dms)
    print(distance)
    print(costhz)

    # print(coord_cern_dms['lon'])