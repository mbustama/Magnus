# -*- coding: utf-8 -*-
"""Tests of the earth and matter helper modules."""

import numpy as np
import pytest

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.matter as matter


def prem_reference_scalar(r):
    """Original (pre-vectorization) piecewise PREM profile, kept here as an
    independent reference."""
    x = r/gd.EARTH_RADIUS
    if 0 <= r <= 1221.5:
        return 13.0885 - 8.8381*x*x
    elif r <= 3480.0:
        return 12.5815 - 1.2638*x - 3.6426*x*x - 5.5281*x*x*x
    elif r <= 5701.0:
        return 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x
    elif r <= 5771.0:
        return 5.3197 - 1.4836*x
    elif r <= 5971.0:
        return 11.2494 - 8.0298*x
    elif r <= 6151.0:
        return 7.1089 - 3.8045*x
    elif r <= 6346.6:
        return 2.6910 + 0.6924*x
    elif r <= 6356.0:
        return 2.900
    elif r <= 6368.0:
        return 2.600
    return 1.020


def test_prem_vectorized_matches_scalar_reference():
    # Dense grid plus every layer boundary from both sides
    boundaries = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0,
                           6346.6, 6356.0, 6368.0, gd.EARTH_RADIUS])
    r = np.concatenate([np.linspace(0.0, gd.EARTH_RADIUS, 20011),
                        boundaries, boundaries - 1e-9,
                        np.minimum(boundaries + 1e-9, gd.EARTH_RADIUS)])
    rho_vec = earth.density_matter_func_prem(r)
    rho_ref = np.array([prem_reference_scalar(x) for x in r])
    assert np.max(np.abs(rho_vec - rho_ref)) < 1e-12


def test_prem_scalar_in_scalar_out():
    rho = earth.density_matter_func_prem(3000.0)
    assert isinstance(rho, float)
    assert rho == pytest.approx(prem_reference_scalar(3000.0))


def test_prem_clamps_radius_within_tolerance():
    rho = earth.density_matter_func_prem(gd.EARTH_RADIUS*(1.0 + 1e-10))
    assert rho == pytest.approx(1.020)


def test_prem_raises_beyond_tolerance():
    with pytest.raises(ValueError):
        earth.density_matter_func_prem(gd.EARTH_RADIUS*1.1)


def test_radial_distance_geometry():
    R = gd.EARTH_RADIUS
    # Vertical crossing (costhz = -1): enters at r=R, passes through the
    # center at l=R, exits at r=R after l=2R
    assert earth.earth_radial_distance_from_depth(-1.0, 0.0) \
        == pytest.approx(R)
    assert earth.earth_radial_distance_from_depth(-1.0, R) \
        == pytest.approx(0.0, abs=1e-6)
    assert earth.earth_radial_distance_from_depth(-1.0, 2.0*R) \
        == pytest.approx(R)
    # Horizontal ray (costhz = 0) at its entry point sits on the surface
    # (regression test: this used to return 0, the center of the Earth)
    assert earth.earth_radial_distance_from_depth(0.0, 0.0) \
        == pytest.approx(R)


def test_radial_distance_vectorized_matches_scalar():
    costhz = -0.7
    d = earth.distance_traveled_inside_earth(costhz)
    ls = np.linspace(0.0, d, 1001)
    r_vec = earth.earth_radial_distance_from_depth(costhz, ls)
    r_scl = np.array([earth.earth_radial_distance_from_depth(costhz, l)
                      for l in ls])
    # Vectorized and scalar sqrt may differ by ~1 ulp; a micrometer is plenty
    assert np.max(np.abs(r_vec - r_scl)) < 1e-9


def test_radial_distance_raises_beyond_path():
    with pytest.raises(ValueError):
        earth.earth_radial_distance_from_depth(-0.5, 2.0*gd.EARTH_RADIUS)


def test_chord_length_haversine():
    # Independent haversine computation for Berlin -> Paris
    lat1, lon1 = (52, 31, 12), (13, 24, 18)
    lat2, lon2 = (48, 51, 24), (2, 21, 7)
    def dec(d, m, s):
        return d + m/60.0 + s/3600.0
    p1 = np.radians([dec(*lat1), dec(*lon1)])
    p2 = np.radians([dec(*lat2), dec(*lon2)])
    a = np.sin((p2[0]-p1[0])/2)**2 \
        + np.cos(p1[0])*np.cos(p2[0])*np.sin((p2[1]-p1[1])/2)**2
    central = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    expected = 2*gd.EARTH_RADIUS*np.sin(central/2)
    got = earth.chord_length_inside_earth(lat1, lon1, lat2, lon2)
    assert got == pytest.approx(expected, rel=1e-12)
    # ~878 km chord, sanity
    assert 850.0 < got < 900.0


def test_num_density_e_vectorizes():
    rho_func = lambda l: np.full_like(np.asarray(l, dtype=float), 3.0)
    ls = np.linspace(0.0, 10.0, 11)
    ne = matter.num_density_e_func(ls, rho_func,
                                   density_matter_is_in_g_per_cm3=True)
    ne0 = matter.num_density_e_func(0.0, lambda l: 3.0,
                                    density_matter_is_in_g_per_cm3=True)
    assert np.allclose(ne, ne0)


def test_vcc_sign_flips_for_antineutrinos():
    vcc_nu = matter.vcc_func_from_rho_func(3.0, density_matter_is_in_g_per_cm3=True)
    vcc_nubar = matter.vcc_func_from_rho_func(3.0, nubar=True,
                                              density_matter_is_in_g_per_cm3=True)
    assert vcc_nu > 0.0
    assert vcc_nubar == pytest.approx(-vcc_nu)
