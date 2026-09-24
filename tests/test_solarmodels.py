# -*- coding: utf-8 -*-
"""Tests of the tabulated standard solar models (magnus.solarmodels) and of the Sun
wrappers' ``density_profile`` and ``stop_at_table_edge``.

Three things are pinned here.  The tables are what their headers say they are: every model
loads, carries its provenance, and gives the electron density the paper's notebooks compute
from the same columns.  The profile does what the documentation promises outside the table:
flat below the first row, the last interval's logarithmic slope past the last.  And the
wrappers route a named model to exactly that profile, with the model's own composition for
the sterile states unless a ratio is given, while the default -- the exponential fit -- is
the call it always was.

The phase-averaged probability is used wherever a probability is needed and the test is not
about the oscillation itself: it costs one propagation, where a coherent probability through
most of the Sun costs thousands of slabs.
"""

import inspect
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op
from magnus import solarmodels as sm

R_SUN = gd.SUN_RADIUS*gd.UNIT_KM
E = 10.0*gd.UNIT_MEV
OSC = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
OSC3 = {k: OSC[k] for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
M_N = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)

SUN_WRAPPERS = [f'osc_prob_{n}nu_sun{s}' for n in (2, 3, 4, 5) for s in ('', '_nsi', '_liv')]
STERILE_WRAPPERS = [w for w in SUN_WRAPPERS if w.startswith(('osc_prob_4nu', 'osc_prob_5nu'))]


def _required(name):
    """The arguments beyond (energy, L, L0) a wrapper cannot do without."""
    return dict(sth=0.55, Dm2=7.5e-5) if name.startswith('osc_prob_2nu') else {}


def _sterile(name):
    """A sterile state that mixes, so that the neutral-current term matters."""
    if name.startswith('osc_prob_4nu'):
        return dict(s14=0.3, s24=0.1, D41=1.0)
    if name.startswith('osc_prob_5nu'):
        return dict(s14=0.3, s24=0.1, D41=1.0, s15=0.2, D51=2.0)
    return {}


# ---------------------------------------------------------------------------------------------
# The tables
# ---------------------------------------------------------------------------------------------

def test_the_catalogue_is_the_twelve_shipped_models():
    assert sm.available_solar_models() == sm.SOLAR_MODELS
    assert len(sm.SOLAR_MODELS) == 12


@pytest.mark.parametrize('name', sm.SOLAR_MODELS)
def test_every_model_loads_and_is_physical(name):
    t = sm.load_solar_model(name)
    r, rho, X = t['r_over_r_sun'], t['rho_g_per_cm3'], t['x_hydrogen']
    assert len(r) == len(rho) == len(X) > 800
    assert np.all(np.diff(r) > 0.0)
    assert 0.0 <= r[0] < 0.01 and 0.94 < r[-1] <= 1.0
    # The centre of the Sun holds about 150 g/cm^3, and hydrogen is burnt down to ~35% there.
    assert 1.4e2 < rho[0] < 1.6e2
    assert np.all(rho > 0.0) and np.all((X > 0.3) & (X < 0.8))


@pytest.mark.parametrize('name', sm.SOLAR_MODELS)
def test_every_model_carries_its_provenance(name):
    info = sm.solar_model_info(name)
    assert info['name'] == name
    for key in ('reference', 'source', 'original', 'retrieved', 'terms'):
        assert info[key], key
    assert len(info['sha256']) == 64 and int(info['sha256'], 16) >= 0
    t = sm.load_solar_model(name)
    assert info['rows'] == len(t['r_over_r_sun'])
    assert (info['r_min'], info['r_max']) == (t['r_over_r_sun'][0], t['r_over_r_sun'][-1])
    # Computed as the profile computes its last row, so that 'past the edge' and 'past the
    # last row' are the same test.
    assert sm.table_edge(name) == info['r_max']*gd.SUN_RADIUS*gd.UNIT_KM


def test_names_match_in_any_case_and_unknown_ones_are_refused():
    assert sm.canonical_name('b16-gs98') == 'B16-GS98'
    assert sm.canonical_name(' BS05-ags-op ') == 'BS05-AGS-OP'
    with pytest.raises(ValueError, match='B23-MB22p'):
        sm.canonical_name('BS99')


# ---------------------------------------------------------------------------------------------
# The profile
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize('name', sm.SOLAR_MODELS)
def test_electron_density_at_the_rows_is_the_tabulated_one(name):
    """n_e = rho (1 + X)/(2 m_N), the formula notebooks 13 and 28 use."""
    t = sm.load_solar_model(name)
    expected = t['rho_g_per_cm3']*gd.UNIT_G_PER_CM3/M_N*0.5*(1.0 + t['x_hydrogen'])
    got = sm.electron_density_profile(name)(t['r_over_r_sun']*R_SUN)
    np.testing.assert_allclose(got, expected, rtol=1e-12)


@pytest.mark.parametrize('name', ['BP04', 'BS05-AGS-OP', 'B16-GS98'])
def test_between_rows_the_logarithm_is_interpolated_linearly(name):
    r = sm.load_solar_model(name)['r_over_r_sun']*R_SUN
    ne = sm.electron_density_profile(name)
    i = len(r)//2
    mid = ne(0.5*(r[i] + r[i + 1]))
    assert mid == pytest.approx(np.sqrt(ne(r[i])*ne(r[i + 1])), rel=1e-12)


@pytest.mark.parametrize('name', sm.SOLAR_MODELS)
def test_outside_the_table_flat_core_and_continued_slope(name):
    r = sm.load_solar_model(name)['r_over_r_sun']*R_SUN
    ne = sm.electron_density_profile(name)
    assert ne(0.0) == ne(r[0])
    slope = np.log(ne(r[-1])/ne(r[-2]))/(r[-1] - r[-2])
    assert slope < 0.0
    for d in (1e-3*R_SUN, 0.05*R_SUN):
        assert ne(r[-1] + d) == pytest.approx(ne(r[-1])*np.exp(slope*d), rel=1e-10)


def test_profile_accepts_scalars_and_arrays():
    ne = sm.electron_density_profile('B23-GS98')
    ratio = sm.neutron_to_proton_ratio_profile('B23-GS98')
    l = np.array([0.1, 0.5, 0.9])*R_SUN
    assert np.ndim(ne(l[0])) == 0 and ne(l).shape == (3,)
    assert np.ndim(ratio(l[0])) == 0 and ratio(l).shape == (3,)
    assert ne(l[1]) == ne(l)[1] and ratio(l[1]) == ratio(l)[1]


@pytest.mark.parametrize('name', ['BP2000', 'B16-AGSS09met', 'B23-C11'])
def test_neutron_to_proton_ratio_is_the_tabulated_composition(name):
    t = sm.load_solar_model(name)
    X = t['x_hydrogen']
    got = sm.neutron_to_proton_ratio_profile(name)(t['r_over_r_sun']*R_SUN)
    np.testing.assert_allclose(got, (1.0 - X)/(1.0 + X), rtol=1e-14)
    # Held at the last tabulated value past the edge.
    assert sm.neutron_to_proton_ratio_profile(name)(2.0*R_SUN) == pytest.approx(
        (1.0 - X[-1])/(1.0 + X[-1]), rel=1e-14)


# ---------------------------------------------------------------------------------------------
# The wrappers
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize('name', SUN_WRAPPERS + ['osc_prob_sun'])
def test_every_sun_entry_point_takes_the_two_keywords(name):
    params = inspect.signature(getattr(op, name)).parameters
    assert params['density_profile'].default == 'exp'
    assert params['stop_at_table_edge'].default is False


@pytest.mark.parametrize('name', SUN_WRAPPERS)
def test_naming_exp_is_the_default_call(name):
    f = getattr(op, name)
    kw = dict(energy=E, L=0.5*R_SUN, L0=0.0, average=True, **_required(name), **_sterile(name))
    assert np.array_equal(np.asarray(f(**kw)), np.asarray(f(density_profile='EXP', **kw)))


@pytest.mark.parametrize('name', SUN_WRAPPERS)
def test_every_wrapper_propagates_through_a_table(name):
    f = getattr(op, name)
    kw = dict(energy=E, L=0.9*R_SUN, L0=0.0, average=True, **_required(name), **_sterile(name))
    P_table = np.asarray(f(density_profile='B23-GS98', **kw))
    P_exp = np.asarray(f(**kw))
    assert np.all(P_table >= -1e-12)
    np.testing.assert_allclose(P_table.sum(axis=-1), 1.0, atol=1e-9)
    assert not np.allclose(P_table, P_exp, atol=1e-6)


@pytest.mark.parametrize('name', sm.SOLAR_MODELS)
def test_a_named_model_is_its_profile_through_the_general_routine(name):
    """Bit for bit: the wrapper builds nothing of its own around the table."""
    kw = dict(nu_i=gd.NUE, nu_f=gd.NUE, L0=0.0, average=True)
    P_wrapper = op.osc_prob_3nu_sun(energy=E, L=0.9*R_SUN, density_profile=name.lower(),
                                    **OSC3, **kw)
    P_general = op.osc_prob_matter_std_potential(3, sm.electron_density_profile(name), E,
                                                 0.9*R_SUN, OSC3,
                                                 density_is_of_number_of_electrons=True, **kw)
    assert float(np.asarray(P_wrapper)) == float(np.asarray(P_general))


@pytest.mark.parametrize('name', sm.SOLAR_MODELS)
def test_averaged_probability_raises_no_warning_on_any_table(name):
    """The tables are smooth enough that the averaged route neither escalates nor warns."""
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        for e in (1.0, 10.0):
            op.osc_prob_3nu_sun(energy=e*gd.UNIT_MEV, L=0.9*R_SUN, L0=0.0, density_profile=name,
                                average=True)


@pytest.mark.parametrize('name', STERILE_WRAPPERS)
def test_sterile_wrappers_take_the_model_composition_unless_told(name):
    f = getattr(op, name)
    kw = dict(energy=E, L=0.9*R_SUN, L0=0.0, average=True, **_sterile(name))
    P_default = np.asarray(f(density_profile='B16-GS98', **kw))
    P_table = np.asarray(f(density_profile='B16-GS98',
        ratio_number_neutrons_to_protons=sm.neutron_to_proton_ratio_profile('B16-GS98'), **kw))
    P_one = np.asarray(f(density_profile='B16-GS98', ratio_number_neutrons_to_protons=1.0, **kw))
    assert np.array_equal(P_default, P_table)
    assert np.max(np.abs(P_default - P_one)) > 1e-5
    # The exponential fit has no composition, and keeps the 1.0 it always had.
    assert np.array_equal(np.asarray(f(**kw)),
                          np.asarray(f(ratio_number_neutrons_to_protons=1.0, **kw)))


def test_an_unknown_model_is_refused_by_the_wrapper():
    with pytest.raises(ValueError, match='osc_prob_3nu_sun: density_profile: unknown solar'):
        op.osc_prob_3nu_sun(energy=E, L=0.5*R_SUN, L0=0.0, density_profile='BS99')


def test_osc_prob_sun_takes_a_model_too():
    osc = {'sth': OSC['s12'], 'Dm2': OSC['D21']}

    def H(energy, l, VCC):
        h = np.array([[0.0, 0.0], [0.0, osc['Dm2']/(2.0*energy)]])
        c, s = np.sqrt(1.0 - osc['sth']**2), osc['sth']
        U = np.array([[c, s], [-s, c]])
        return U @ h @ U.T + np.asarray(VCC)[..., None, None]*np.diag([1.0, 0.0])

    kw = dict(energy=E, L=0.9*R_SUN, average=True, nu_i=0, nu_f=0)
    P_general = float(np.asarray(op.osc_prob_sun(H, density_profile='B16-GS98', **kw)))
    P_wrapper = float(np.asarray(op.osc_prob_2nu_sun(L0=0.0, density_profile='B16-GS98', **osc,
                                                     **kw)))
    assert P_general == pytest.approx(P_wrapper, abs=1e-6)
    assert np.array_equal(np.asarray(op.osc_prob_sun(H, **kw)),
                          np.asarray(op.osc_prob_sun(H, density_profile='exp', **kw)))


# ---------------------------------------------------------------------------------------------
# stop_at_table_edge
# ---------------------------------------------------------------------------------------------

def test_stop_at_table_edge_blanks_only_what_lies_past_it():
    edge = sm.table_edge('BP04')
    L = np.array([0.5*R_SUN, 0.99*R_SUN, 0.9*R_SUN])
    kw = dict(energy=E, L0=0.0, nu_i=0, nu_f=0, average=True, density_profile='BP04')
    with pytest.warns(op.SolarModelRangeWarning, match=r'1 of 3 baseline\(s\).*0\.9468 R_sun'):
        P = np.asarray(op.osc_prob_3nu_sun(L=L, stop_at_table_edge=True, **kw))
    assert np.isnan(P[1]) and not np.any(np.isnan(P[[0, 2]]))
    # What is kept is what the same call computes with the refused point at the origin.
    P_same = np.asarray(op.osc_prob_3nu_sun(L=np.where(L > edge, 0.0, L), **kw))
    assert np.array_equal(P[[0, 2]], P_same[[0, 2]])


def test_stop_at_table_edge_is_silent_inside_the_table():
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        P = op.osc_prob_3nu_sun(energy=E, L=0.9*R_SUN, L0=0.0, average=True,
                                density_profile='BP04', stop_at_table_edge=True)
    assert not np.any(np.isnan(np.asarray(P)))


@pytest.mark.parametrize('kw', [dict(), dict(nu_i=0, nu_f=0), dict(average=True)],
                         ids=['matrix', 'channel', 'averaged'])
def test_stop_at_table_edge_blanks_every_shape(kw):
    with pytest.warns(op.SolarModelRangeWarning):
        P = op.osc_prob_4nu_sun(energy=E, L=0.99*R_SUN, L0=0.0, s14=0.1, D41=1.0,
                                density_profile='BS05-OP', stop_at_table_edge=True, **kw)
    assert np.all(np.isnan(np.asarray(P)))


def test_stop_at_table_edge_blanks_the_evolution_operator_too():
    with pytest.warns(op.SolarModelRangeWarning):
        P, U = op.osc_prob_3nu_sun(energy=E, L=0.99*R_SUN, L0=0.0, density_profile='BP2000',
                                   stop_at_table_edge=True, return_evolution_operator=True)
    assert np.all(np.isnan(P)) and np.all(np.isnan(U))


def test_stop_at_table_edge_reaches_osc_prob_sun():
    def H(energy, l, VCC):
        return np.diag([0.0, 7.5e-5/(2.0*energy)]) + np.asarray(VCC)[..., None, None]*np.diag([1.0, 0.0])

    with pytest.warns(op.SolarModelRangeWarning):
        P = op.osc_prob_sun(H, energy=E, L=0.99*R_SUN, density_profile='BP04',
                            stop_at_table_edge=True, average=True)
    assert np.all(np.isnan(np.asarray(P)))


def test_stop_at_table_edge_refuses_what_it_cannot_serve():
    with pytest.raises(ValueError, match='exponential fit has no last row'):
        op.osc_prob_3nu_sun(energy=E, L=0.5*R_SUN, L0=0.0, stop_at_table_edge=True)
    with pytest.raises(ValueError, match='path starts past the last tabulated radius'):
        op.osc_prob_3nu_sun(energy=E, L=0.99*R_SUN, L0=0.97*R_SUN, density_profile='BP04',
                            stop_at_table_edge=True)
