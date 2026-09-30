r"""The #160 checklist items #173 did not implement, one test per item of the reopened issue.

Each test reproduces the case the audit (``docs/dev/audit/audit160.py``) found failing on
``main`` at 0939598, and requires what the checklist asks: a Magnus error naming the argument,
a warning, or the documented behaviour.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

E = 1.0*gd.UNIT_GEV
L = 1000.0*gd.UNIT_KM
KM = gd.UNIT_KM
OSC = gd.load_nufit_params('NuFIT 6.1')


def _refused(match, fn, *args, **kwargs):
    with pytest.raises(ValueError) as err:
        fn(*args, **kwargs)
    msg = str(err.value)
    assert 'Error in magnus' in msg and match in msg, msg
    return msg


# §11 angles beyond 90 degrees -------------------------------------------------------------

@pytest.mark.parametrize('s12, angles', [(2.0, 'rad'), (-2.0, 'rad'), (120.0, 'deg'),
                                         (-100.0, 'deg')])
def test_angle_with_negative_cosine_is_refused(s12, angles):
    kw = dict(s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3)
    if angles == 'deg':
        kw.update(s23=40.0, s13=8.6, dCP=200.0)
    msg = _refused('s12', op.osc_prob_3nu_vacuum, E, L, s12=s12, angles=angles, **kw)
    assert 'negative cosine' in msg


@pytest.mark.parametrize('s12, angles', [(-1.0, 'rad'), (np.pi/2, 'rad'), (90.0, 'deg'),
                                         (-90.0, 'deg'), (33.8, 'deg')])
def test_angle_in_first_or_fourth_quadrant_is_accepted(s12, angles):
    kw = dict(s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3)
    if angles == 'deg':
        kw.update(s23=40.0, s13=8.6, dCP=200.0)
    P = op.osc_prob_3nu_vacuum(E, L, s12=s12, angles=angles, **kw)
    assert np.allclose(np.sum(P, axis=-1), 1.0)


# §11 angle-convention warning ---------------------------------------------------------------

def test_angle_warning_judges_by_the_active_angles():
    """Sines in the active slots are not hidden by 5-degree sterile angles."""
    with pytest.warns(gd.MixingAngleConventionWarning):
        op.osc_prob_4nu_vacuum(E, L, s12=0.55, s23=0.75, s13=0.15, dCP=1.0, s14=5.0, s24=5.0,
                               s34=5.0, d14=0.0, d24=0.0, D21=7.5e-5, D31=2.5e-3, D41=1.0,
                               angles='deg')


def test_angle_warning_says_how_to_silence_a_false_positive():
    with pytest.warns(gd.MixingAngleConventionWarning, match='filterwarnings'):
        op.osc_prob_2nu_vacuum(E, L, sth=0.5, Dm2=7.5e-5, angles='deg')


def test_angle_warning_is_quiet_for_real_degrees():
    with warnings.catch_warnings():
        warnings.simplefilter('error', gd.MixingAngleConventionWarning)
        op.osc_prob_4nu_vacuum(E, L, s12=33.7, s23=43.3, s13=8.6, dCP=212.0, s14=0.5, s24=0.5,
                               s34=0.5, d14=0.0, d24=0.0, D21=7.5e-5, D31=2.5e-3, D41=1.0,
                               angles='deg')


# §11 per-node matter helpers --------------------------------------------------------------

@pytest.mark.parametrize('layer', ['core', 'mantle', 'crust', 'ocean'])
@pytest.mark.parametrize('value', [-0.1, 1.5, np.nan, np.inf, True])
def test_electron_fraction_func_prem_checks_its_fractions(layer, value):
    import magnus.earth as earth
    name = 'electron_fraction_' + layer
    with pytest.raises(ValueError, match=name):
        earth.electron_fraction_func_prem(1000.0, **{name: value})


@pytest.mark.parametrize('value', [0.0, -0.1, 1.5, np.nan, np.inf, True])
def test_neutron_to_proton_ratio_checks_the_electron_fraction(value):
    import magnus.earth as earth
    with pytest.raises(ValueError, match='electron_fraction'):
        earth.neutron_to_proton_ratio_from_electron_fraction(value)


def test_matter_helpers_still_take_arrays():
    import magnus.earth as earth
    import magnus.matter as matter
    ye = np.array([0.4656, 0.4957])
    assert np.allclose(earth.neutron_to_proton_ratio_from_electron_fraction(ye), (1 - ye)/ye)
    n = matter.num_density_e_func(np.array([0.0, 1.0]), lambda l: 3.0 + 0*np.asarray(l),
                                  ratio_number_neutrons_to_protons=(1 - ye)/ye,
                                  electron_fraction=ye, density_matter_is_in_g_per_cm3=True)
    assert np.all(np.isfinite(n))
    with pytest.raises(ValueError, match='electron_fraction'):
        earth.neutron_to_proton_ratio_from_electron_fraction(np.array([0.5, 1.5]))


@pytest.mark.parametrize('value', [1.5, -0.1, np.nan, np.inf, 1j, True])
def test_num_density_e_func_checks_the_electron_fraction(value):
    import magnus.matter as matter
    with pytest.raises(ValueError, match='electron_fraction'):
        matter.num_density_e_func(0.0, lambda l: 3.0, electron_fraction=value)


def test_num_density_e_func_checks_ratio_and_density_function():
    import magnus.matter as matter
    with pytest.raises(ValueError, match='ratio_number_neutrons_to_protons'):
        matter.num_density_e_func(0.0, lambda l: 3.0, ratio_number_neutrons_to_protons=-1.0)
    with pytest.raises(ValueError, match='density_matter_func'):
        matter.num_density_e_func(0.0, 3)


# §11 hot-path notes ---------------------------------------------------------------------------

def test_unchecked_per_node_functions_say_so():
    import magnus.earth as earth
    import magnus.hamiltonians as hams
    import magnus.matter as matter
    funcs = [matter.density_matter_func_const, matter.density_matter_func_exp,
             earth.density_matter_func_prem]
    funcs += [getattr(hams, 'hamiltonian_%dnu_%s' % (n, k)) for n in (2, 3, 4, 5)
              for k in ('matter', 'matter_td')]
    for f in funcs:
        assert 'hot path' in f.__doc__, f.__name__


# §10 commutator ------------------------------------------------------------------------------

def test_commutator_with_mismatched_shapes_names_its_arguments():
    import magnus.magnus as mm
    _refused('X and Y', mm.commutator, np.eye(2), np.eye(3))
    X = np.array([[0.0, 1.0], [0.0, 0.0]])
    assert np.array_equal(mm.commutator(X, X.T), X @ X.T - X.T @ X)


# §3 Earth helpers -------------------------------------------------------------------------------

@pytest.mark.parametrize('args, name', [((500, 0, 0), 'degrees'), ((-361, 0, 0), 'degrees'),
                                        ((10, -1, 0), 'minutes'), ((10, 5, -3), 'seconds')])
def test_dms_to_decimal_refuses_out_of_range_and_ambiguous_signs(args, name):
    import magnus.earth as earth
    _refused(name, earth.dms_to_decimal, *args)


@pytest.mark.parametrize('fn', ['chord_length_inside_earth', 'costhz_between_points_on_surface'])
def test_chord_helpers_name_the_longitude_beyond_360(fn):
    import magnus.earth as ea
    lat = (41, 50, 0)
    _refused('lon1_dms', getattr(ea, fn), lat, (400, 0, 0), (44, 21, 0), (-103, 45, 0))


@pytest.mark.parametrize('args, value', [((-46, 12, 0), -46.2), ((0, -30, 0), -0.5),
                                         ((359, 30, 0), 359.5), ((-88, -15, -36), -88.26)])
def test_dms_to_decimal_still_reads_signed_triples(args, value):
    import magnus.earth as earth
    assert abs(earth.dms_to_decimal(*args) - value) < 1e-12


@pytest.mark.parametrize('l', [np.nan, -10.0, True, np.array([10.0, np.nan])])
def test_earth_radial_distance_from_depth_checks_the_position(l):
    import magnus.earth as earth
    with pytest.raises(ValueError, match=r'\bl must be'):
        earth.earth_radial_distance_from_depth(-0.5, l)


def test_earth_radial_distance_from_depth_still_accepts_positions():
    import magnus.earth as earth
    assert np.all(np.isfinite(earth.earth_radial_distance_from_depth(-0.5, np.array([0.0, 100.0]))))


def test_earth_wrappers_document_the_partial_path():
    import re
    for name in dir(op):
        if re.match(r'osc_prob_\w*earth', name):
            assert 'partial path' in ' '.join(getattr(op, name).__doc__.split()), name


# §1 physics arguments ------------------------------------------------------------------------------

@pytest.mark.parametrize('arg, value', [('l_scale', np.nan), ('l_scale', 0.0),
                                        ('rho_central', np.nan), ('rho_central', -3.0)])
def test_exponential_profile_arguments_are_named(arg, value):
    kw = dict(rho_central=3.0, l_scale=300.0*KM)
    kw[arg] = value
    msg = _refused(arg, op.osc_prob_3nu_matter_exp_density, E, L, 0.0,
                   density_matter_is_in_g_per_cm3=True, **kw)
    assert 'osc_prob_3nu_matter_exp_density' in msg


RS = gd.SUN_RADIUS*KM


def _sun_kw(name):
    n = int(name[9])
    kw = dict(sth=0.55, Dm2=7.5e-5) if n == 2 else {}
    return kw


@pytest.mark.parametrize('name', [n for n in dir(op) if n.startswith('osc_prob_') and '_sun' in n
                                  and n[9].isdigit()])
def test_every_sun_wrapper_refuses_a_negative_start(name):
    f = getattr(op, name)
    _refused('L0', f, 1.0e7, RS, -1.0e5*KM, **_sun_kw(name))


def test_custom_sun_refuses_a_negative_start():
    H0v = hams_vacuum()
    _refused('L0', op.osc_prob_sun, lambda En, l, V: H0v/En + np.asarray(V)[..., None, None]*np.diag([1., 0, 0]),
             1.0e7, RS, -1.0e5*KM)


def hams_vacuum():
    import magnus.hamiltonians as hams
    return hams.hamiltonian_3nu_vacuum_energy_independent(**OSC)


@pytest.mark.parametrize('profile', ['exponential', 'B16-GS98'])
def test_sun_refuses_arguments_with_no_effect(profile):
    _refused('electron_fraction', op.osc_prob_3nu_sun, 1.0e7, RS, 0.0,
             density_profile=profile, electron_fraction=0.5)
    _refused('ratio_number_neutrons_to_protons', op.osc_prob_3nu_sun, 1.0e7, RS, 0.0,
             density_profile=profile, ratio_number_neutrons_to_protons=1.0)


def test_sun_start_at_zero_and_four_flavor_ratio_still_work():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        op.osc_prob_3nu_sun(1.0e7, RS, 0.0, nu_i=0, nu_f=0)
        op.osc_prob_3nu_sun(1.0e7, RS, 0, nu_i=0, nu_f=0)
        op.osc_prob_4nu_sun(1.0e7, RS, 0.0, s14=0.1, s24=0.05, s34=0.02, D41=1.0,
                            ratio_number_neutrons_to_protons=0.5, nu_i=0, nu_f=0)


def test_sun_docs_say_L0_is_a_single_radius():
    import os
    import re
    for name in dir(op):
        if re.match(r'osc_prob_(\dnu_)?sun', name):
            assert 'loop over' in ' '.join(getattr(op, name).__doc__.split()), name
    rst = os.path.join(os.path.dirname(__file__), '..', 'docs', 'source', 'solar_models.rst')
    if os.path.exists(rst):
        text = ' '.join(open(rst).read().split())
        assert 'may not, so for several production points call once per point' in text


@pytest.mark.parametrize('energy, L', [
    ([1.0e9], [100.*KM, 300.*KM, 800.*KM]),
    ([1.0e9, 2.0e9, 3.0e9], np.array([500.*KM])),
])
def test_wrappers_broadcast_a_single_entry_like_osc_prob_energy_baseline(energy, L):
    got = op.osc_prob_3nu_vacuum(energy, L, **OSC)
    n = max(len(energy), len(L))
    expected = op.osc_prob_3nu_vacuum(np.broadcast_to(energy, n).copy(),
                                      np.broadcast_to(L, n).copy(), **OSC)
    assert np.array_equal(got, expected)


def test_wrappers_still_refuse_unequal_lengths_above_one():
    _refused('single entry', op.osc_prob_3nu_vacuum, [1.0e9, 2.0e9], [1.*KM, 2.*KM, 3.*KM], **OSC)


@pytest.mark.parametrize('fn', ['osc_prob_3nu_earth', 'osc_prob_3nu_earth_nsi', 'osc_prob_3nu_earth_liv'])
@pytest.mark.parametrize('ratio', [0.5, lambda l: 1.0])
def test_three_flavor_earth_refuses_a_ratio_it_ignores(fn, ratio):
    _refused('ratio_number_neutrons_to_protons', getattr(op, fn), 1.0e9, costhz=-0.5,
             L=1000.*KM, ratio_number_neutrons_to_protons=ratio)


def test_two_flavor_earth_refuses_a_ratio_it_ignores():
    _refused('ratio_number_neutrons_to_protons', op.osc_prob_2nu_earth, 1.0e9, costhz=-0.5,
             L=1000.*KM, sth=0.5, Dm2=2.5e-3, ratio_number_neutrons_to_protons=1.0)


STERILE = dict(OSC, s14=0.1, s24=0.1, s34=0.1, D41=1.0)


def test_four_flavor_earth_checks_a_callable_ratio():
    _refused('ratio_number_neutrons_to_protons', op.osc_prob_4nu_earth, 1.0e9, costhz=-0.5,
             L=1000.*KM, ratio_number_neutrons_to_protons=lambda l: np.nan, **STERILE)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = op.osc_prob_4nu_earth(1.0e9, costhz=-0.5, L=1000.*KM,
                                  ratio_number_neutrons_to_protons=lambda l: 1.0 + 0.0*np.asarray(l),
                                  **STERILE)
    assert np.all(np.isfinite(a))


@pytest.mark.parametrize('scenario', ['osc_prob_vacuum', 'osc_prob_matter_std_potential'])
def test_scenario_functions_refuse_a_vacuum_hamiltonian_they_would_ignore(scenario):
    H4 = np.diag([0.0, 7.4e-5, 2.5e-3, 1.0]).astype(complex)
    params = dict(OSC, s14=0.1, d14=0.3, s24=0.15, d24=0.5, s34=0.2, D41=1.0)
    args = (4,) + ((3.0*gd.UNIT_G_PER_CM3,) if scenario != 'osc_prob_vacuum' else ()) + \
        (1.0e9, 1000.*KM, params)
    _refused('h_vac_energy_indep', getattr(op, scenario), *args, h_vac_energy_indep=H4)


def test_sample_count_without_average_is_refused():
    _refused('average_n_samples', op.osc_prob_3nu_vacuum, 1.0e9, 1000.*KM, average_n_samples=11,
             **OSC)
    _refused('average_n_samples', op.osc_prob_3nu_earth, 1.0e9, costhz=-0.5, L=5000.*KM,
             average_n_samples=11, **OSC)


def test_spread_and_initial_state_stay_ignored_without_average():
    plain = op.osc_prob_3nu_vacuum(1.0e9, 1000.*KM, **OSC)
    for kw in (dict(average_spread=0.1), dict(average_initial_state='decohered')):
        assert np.array_equal(op.osc_prob_3nu_vacuum(1.0e9, 1000.*KM, **kw, **OSC), plain)
