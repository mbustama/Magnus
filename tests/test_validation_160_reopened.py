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
