# -*- coding: utf-8 -*-
"""Tests of the earth and matter helper modules."""

import warnings

import numpy as np
import pytest

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.oscprob as op


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


def test_prem_layer_edges_along_chord():
    # Every returned crossing must sit exactly on a PREM boundary radius,
    # and crossings must come in symmetric pairs about the chord midpoint
    costhz = -0.8
    d = earth.distance_traveled_inside_earth(costhz)
    ls = earth.prem_layer_edges_along_chord(costhz)
    assert len(ls) > 0
    assert np.all((ls > 0.0) & (ls < d))
    rs = earth.earth_radial_distance_from_depth(costhz, ls)
    dist_to_boundary = np.min(np.abs(rs[:, None] - earth.PREM_BOUNDARIES),
                              axis=1)
    assert np.max(dist_to_boundary) < 1e-6  # [km]
    assert np.allclose(np.sort(d - ls), np.sort(ls))  # symmetric pairs
    # A vertical chord (costhz = -1) crosses every boundary twice
    ls_vert = earth.prem_layer_edges_along_chord(-1.0)
    assert len(ls_vert) == 2*len(earth.PREM_BOUNDARIES)
    # A down-going direction has no crossings
    assert len(earth.prem_layer_edges_along_chord(0.5)) == 0
    # A shallow chord crosses only the outermost layers
    ls_shallow = earth.prem_layer_edges_along_chord(-0.05)
    rs_shallow = earth.earth_radial_distance_from_depth(-0.05, ls_shallow)
    assert np.min(rs_shallow) > 6151.0


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


# ----------------------------------------------------------------------
# Error and convenience paths a coverage run found unexercised.
# ----------------------------------------------------------------------

def test_named_location_lookup_is_case_and_space_insensitive():
    """The lookup normalizes case and spaces, so the same detector can be
    named the way a person would write it."""
    canonical = earth.coordinates_of_named_location('test', loc_name='gran_sasso')
    for variant in ['Gran Sasso', 'GRAN SASSO', 'gran sasso', 'Gran_Sasso']:
        assert np.array_equal(earth.coordinates_of_named_location('test', loc_name=variant),
                              canonical), variant


def test_unknown_named_location_raises_and_lists_the_known_ones():
    """An unrecognized site name must fail loudly, and the message must name
    the alternatives -- this is the one error a user hits by typing a
    detector name slightly wrong, so the list is the whole remedy."""
    with pytest.raises(ValueError) as excinfo:
        earth.coordinates_of_named_location('test', loc_name='not_a_detector')

    message = str(excinfo.value)
    assert 'not_a_detector' in message, "the rejected name is not echoed back"
    for known in earth.loc_coords_dms:
        assert known in message, f"the available location {known!r} is not offered"


def test_constant_density_profile_is_constant():
    rho = 4.5
    assert matter.density_matter_func_const(0.0, rho) == rho
    assert matter.density_matter_func_const(1.0e9, rho) == rho


def test_constant_density_profile_defaults_to_the_crust():
    assert matter.density_matter_func_const(0.0) == gd.DENSITY_MATTER_CRUST_G_PER_CM3


def test_constant_electron_number_density_gives_a_constant_potential():
    """When the profile is already an electron number density *and* a
    constant, the potential is a number rather than a callable: there is
    nothing for it to depend on. This is the one route through
    vcc_func_from_rho_func that neither converts units nor wraps a
    function, and it is what a caller passing a fixed n_e gets."""
    n_e = 1.0e-6   # [eV^3]
    vcc = matter.vcc_func_from_rho_func(n_e, density_is_of_number_of_electrons=True)

    assert not callable(vcc)
    assert vcc > 0.0

    # The antineutrino potential is the same number with the opposite sign.
    vcc_bar = matter.vcc_func_from_rho_func(n_e, nubar=True,
                                            density_is_of_number_of_electrons=True)
    assert vcc_bar == pytest.approx(-vcc)

    # And it is the constant that the callable form returns at every position.
    vcc_of_l = matter.vcc_func_from_rho_func(lambda l: n_e,
                                             density_is_of_number_of_electrons=True)
    assert callable(vcc_of_l)
    assert vcc_of_l(0.0) == pytest.approx(vcc)
    assert vcc_of_l(1.0e9) == pytest.approx(vcc)


# ----------------------------------------------------------------------
# The double-conversion guard
# ----------------------------------------------------------------------

def test_a_density_already_in_natural_units_is_flagged():
    """Passing an already-converted density while declaring it to be in
    g/cm^3 converts it twice, inflating the matter potential by some
    eighteen orders of magnitude.

    What makes this worth catching is that the result does not look like a
    unit error: the matter term swamps everything, nu_e becomes an exact
    eigenstate of the Hamiltonian, and the calculation returns a perfectly
    self-consistent P_ee = 1 -- which reads as a broken formula rather
    than as bad input."""
    rho_internal = 100.0*gd.UNIT_G_PER_CM3

    with pytest.warns(matter.DensityUnitWarning, match="natural units|neutron star"):
        vcc = matter.vcc_func_from_rho_func(matter.exp_density_profile(rho_internal, 1.0e5),
                                            density_matter_is_in_g_per_cm3=True)
        np.asarray(vcc(0.0))


@pytest.mark.parametrize("rho, in_g_per_cm3", [
    (100.0, True),                       # a real density, correctly declared
    (100.0*gd.UNIT_G_PER_CM3, False),    # already converted, correctly declared
    (13.0, True),                        # Earth's core
    (150.0, True),                       # the centre of the Sun
])
def test_plausible_densities_are_not_flagged(rho, in_g_per_cm3):
    """The threshold sits three orders of magnitude above the densest
    matter anyone models, so nothing physical trips it."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", matter.DensityUnitWarning)
        vcc = matter.vcc_func_from_rho_func(matter.exp_density_profile(rho, 1.0e5),
                                            density_matter_is_in_g_per_cm3=in_g_per_cm3)
        np.asarray(vcc(0.0))


def test_the_guard_only_applies_when_g_per_cm3_is_declared():
    """A large number in natural units is perfectly ordinary; it is only
    suspicious when the caller says it is a density in g/cm^3."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", matter.DensityUnitWarning)
        vcc = matter.vcc_func_from_rho_func(matter.exp_density_profile(1.0e30, 1.0e5))
        np.asarray(vcc(0.0))


def test_a_g_per_cm3_density_left_undeclared_is_flagged():
    """The mirror of the double-conversion guard, and the commoner mistake:
    density_matter_is_in_g_per_cm3 defaults to False, so a density read off a
    table is taken as already converted and the matter potential comes out ~19
    orders of magnitude too small.

    Under-conversion is quieter than double conversion -- it does not inflate
    anything, it makes the matter term vanish, and the call returns exactly the
    vacuum probability.  That is an ordinary-looking number of the right shape,
    which is why it needs saying out loud."""
    with pytest.warns(matter.DensityUnitWarning, match="too small|natural units"):
        vcc = matter.vcc_func_from_rho_func(matter.exp_density_profile(2.848, 1.0e5))
        np.asarray(vcc(0.0))


def test_an_undeclared_g_per_cm3_density_returns_exactly_the_vacuum_answer():
    """Why the guard above is worth having: the wrong answer is not merely
    close to the vacuum one, it *is* the vacuum one."""
    import magnus.oscprob as oscprob

    osc = gd.load_nufit_params('NuFIT 6.1', 'NO')
    energy, baseline = 2.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM

    vacuum = oscprob.osc_prob_3nu_vacuum(energy, baseline, **osc)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", matter.DensityUnitWarning)
        undeclared = oscprob.osc_prob_3nu_matter_constant_density(
            energy, baseline, 2.848, **osc)
    declared = oscprob.osc_prob_3nu_matter_constant_density(
        energy, baseline, 2.848, **osc, density_matter_is_in_g_per_cm3=True)

    np.testing.assert_allclose(np.asarray(undeclared), np.asarray(vacuum),
                               rtol=0.0, atol=1.0e-14)
    assert not np.allclose(np.asarray(declared), np.asarray(vacuum), atol=1.0e-3)


@pytest.mark.parametrize("rho", [
    2.848*gd.UNIT_G_PER_CM3,     # the Earth's crust, correctly converted
    150.0*gd.UNIT_G_PER_CM3,     # the centre of the Sun, correctly converted
    0.0,                         # deliberate vacuum
])
def test_natural_unit_densities_are_not_flagged(rho):
    """Anything physical is 4.3e18 or more once converted, nine orders above
    the threshold; and an explicit zero is a deliberate vacuum, not a mistake."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", matter.DensityUnitWarning)
        vcc = matter.vcc_func_from_rho_func(matter.exp_density_profile(rho, 1.0e5))
        np.asarray(vcc(0.0))


def test_the_electron_density_path_is_never_flagged():
    """Electron number densities are ~1e9-1e12 eV^3 -- straddling the
    threshold -- but they do not go through the matter-density conversion at
    all, so the guard must never see them."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", matter.DensityUnitWarning)
        vcc = matter.vcc_func_from_rho_func(
            matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN),
            density_is_of_number_of_electrons=True)
        np.asarray(vcc(0.0))


def test_the_undeclared_density_guard_survives_a_repeated_call():
    """The constant-V_CC memo skips the conversion that both unit guards live inside.

    Mirroring only the declared-in-g/cm^3 arm on a cache hit left the commoner and
    quieter mistake -- the undeclared one, which returns exactly the vacuum answer --
    warned about once and then silent, which is precisely the shape a scan or a fit
    has.  A density that trips either guard is now not memoised at all, so the guard
    keeps firing from where it always fired from.
    """
    matter._VCC_CONST_CACHE.clear()
    for _ in range(3):
        with pytest.warns(matter.DensityUnitWarning, match="too small|natural units"):
            matter.vcc_func_from_rho_func(2.848, L0=0.0)


def test_a_density_that_trips_a_guard_is_not_memoised():
    """The mechanism behind the test above, and the reason it is done this way round.

    The warning cannot simply be re-emitted from the cache-hit site: ``warnings.warn``
    is called there with a ``stacklevel`` that attributes it to a different frame, the
    frame is part of the interpreter's warning registry key, and the imitation
    therefore printed a *second* warning under the default filter where an uncached
    call printed one.  Declining to cache costs nothing, because a density that trips
    a guard is a mistake being reported rather than a hot path.
    """
    matter._VCC_CONST_CACHE.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", matter.DensityUnitWarning)
        matter.vcc_func_from_rho_func(2.848, L0=0.0)
    assert len(matter._VCC_CONST_CACHE) == 0

    # ...while an ordinary density is still memoised, which is the point of the cache.
    matter._VCC_CONST_CACHE.clear()
    first = matter.vcc_func_from_rho_func(2.848, L0=0.0,
                                          density_matter_is_in_g_per_cm3=True)
    assert len(matter._VCC_CONST_CACHE) == 1
    assert matter.vcc_func_from_rho_func(
        2.848, L0=0.0, density_matter_is_in_g_per_cm3=True) == first


# ----------------------------------------------------------------------
# Composition: Y_e by PREM layer
# ----------------------------------------------------------------------

def test_each_prem_layer_gets_its_own_electron_fraction():
    """PREM is a density model and carries no composition, so Y_e has to be
    supplied.  The library assumed 0.5 everywhere -- exactly isoscalar matter,
    which nothing in the Earth is."""
    cases = [(1000.0, earth.Y_E_CORE_PREM), (3480.0, earth.Y_E_CORE_PREM),
             (3480.1, earth.Y_E_MANTLE_PREM), (6346.6, earth.Y_E_MANTLE_PREM),
             (6346.7, earth.Y_E_CRUST_PREM), (6368.0, earth.Y_E_CRUST_PREM),
             (6368.1, earth.Y_E_OCEAN_PREM), (6371.0, earth.Y_E_OCEAN_PREM)]
    for r, expected in cases:
        got = float(earth.electron_fraction_func_prem(r))
        assert got == pytest.approx(expected), 'r = %.1f km gave Y_e = %r' % (r, got)


def test_the_neutron_ratio_is_derived_from_the_electron_fraction():
    """They are the same statement about composition: with charge neutrality,
    r = n_n/n_p = (1 - Y_e)/Y_e.  Carried as independent arguments, a caller
    could describe an iron core with isoscalar neutrons -- and silently, since
    r only shows up in the sterile sector."""
    assert earth.neutron_to_proton_ratio_from_electron_fraction(0.5) == pytest.approx(1.0)
    for ye in (earth.Y_E_CORE_PREM, earth.Y_E_MANTLE_PREM, earth.Y_E_CRUST_PREM, earth.Y_E_OCEAN_PREM):
        assert earth.neutron_to_proton_ratio_from_electron_fraction(ye) == \
            pytest.approx((1.0 - ye)/ye)


def test_the_uniform_override_reproduces_the_old_isoscalar_answer():
    """`electron_fraction=0.5` is how a result computed before the layered
    model is reproduced, so it has to keep meaning exactly that."""
    costhz = -0.9
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    energy = 1.0*gd.UNIT_GEV
    uniform = np.asarray(op.osc_prob_3nu_earth(energy, costhz=costhz, L=L,
                                               electron_fraction=0.5))
    layered = np.asarray(op.osc_prob_3nu_earth(energy, costhz=costhz, L=L))
    np.testing.assert_allclose(uniform.sum(axis=1), 1.0, atol=1.0e-12)
    assert not np.allclose(uniform, layered), \
        'the layered default should differ from uniform 0.5'


def test_a_core_crossing_chord_moves_most():
    """The core is iron (Y_e = 0.4656) and the mantle is rock (0.4957), so the
    size of the correction tracks how much core the chord crosses.  A shallow
    chord sees almost none of it."""
    energy = 1.0*gd.UNIT_GEV
    def channel(costhz, **kw):
        L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
        return np.asarray(op.osc_prob_3nu_earth(energy, costhz=costhz, L=L,
                                                **kw))[gd.NUMU][gd.NUE]
    deep = abs(channel(-1.0) - channel(-1.0, electron_fraction=0.5))
    shallow = abs(channel(-0.4) - channel(-0.4, electron_fraction=0.5))
    assert deep > 10.0*shallow, 'deep %r vs shallow %r' % (deep, shallow)


def test_uniform_and_per_layer_together_are_refused():
    """Any precedence rule would be a rule the caller has to know, and this
    package has already shipped two bugs of that shape."""
    costhz = -0.9
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    with pytest.raises(ValueError, match="one or the other"):
        op.osc_prob_3nu_earth(1.0*gd.UNIT_GEV, costhz=costhz, L=L,
                              electron_fraction=0.5, electron_fraction_core=0.46)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5])
def test_an_electron_fraction_outside_zero_to_one_is_refused(bad):
    """Y_e = <Z/A> is a fraction.  0.0 and 5.0 used to be accepted, returning
    answers 0.51 and 0.74 away from the default."""
    costhz = -0.9
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    with pytest.raises(ValueError, match="electron fraction"):
        op.osc_prob_3nu_earth(1.0*gd.UNIT_GEV, costhz=costhz, L=L,
                              electron_fraction_core=bad)


def test_the_sun_takes_no_electron_fraction():
    """The solar profile is an electron NUMBER density -- the standard
    exponential fit, in which Y_e is already folded in -- so the mass-density
    conversion never runs and Y_e does not enter.  `sun_liv` used to expose
    four such parameters and ignore all of them."""
    import inspect
    for name in ('osc_prob_2nu_sun', 'osc_prob_3nu_sun_nsi', 'osc_prob_3nu_sun_liv',
                 'osc_prob_5nu_sun_liv'):
        params = inspect.signature(getattr(op, name)).parameters
        for dead in ('electron_fraction', 'density_matter_is_in_g_per_cm3',
                     'density_is_of_number_of_electrons'):
            assert dead not in params, '%s still exposes %s' % (name, dead)


def test_the_neutron_to_proton_ratio_is_inert_on_the_sun_only_below_four_flavours():
    """It has two jobs, and only one of them is dead on a solar profile.

    Converting a mass density to an electron number density needs the average nucleon
    mass, which depends on r -- and that conversion never runs for the Sun, which is
    where the test above comes from.  But r also sets the sterile states' entry in the
    matter projector, r/2, and that has nothing to do with the density conversion.

    So removing it everywhere went one wrapper too far.  At two and three flavours the
    projector's sterile block is empty and r genuinely has nowhere to act.  At four and
    five it does: on ``osc_prob_5nu_sun_liv``, which previously pinned r = 1.0 in its
    own delegation with no way for a caller to reach it, moving to the Sun's own
    r ~ 0.29 changes the averaged survival probability by 4.5e-03 -- above the default
    tolerance.  The Sun is hydrogen-rich, so isoscalar r = 1.0 is not merely one choice
    among several: Y_e = (1 + X)/2 puts the true r between about 0.47 and 0.14, and 1.0
    is outside that range entirely.
    """
    import inspect
    for name in ('osc_prob_2nu_sun', 'osc_prob_3nu_sun', 'osc_prob_3nu_sun_nsi',
                 'osc_prob_3nu_sun_liv'):
        assert 'ratio_number_neutrons_to_protons' not in inspect.signature(
            getattr(op, name)).parameters, name
    for name in ('osc_prob_4nu_sun', 'osc_prob_5nu_sun', 'osc_prob_4nu_sun_nsi',
                 'osc_prob_5nu_sun_nsi', 'osc_prob_4nu_sun_liv', 'osc_prob_5nu_sun_liv'):
        assert 'ratio_number_neutrons_to_protons' in inspect.signature(
            getattr(op, name)).parameters, name
