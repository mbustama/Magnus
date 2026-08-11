# -*- coding: utf-8 -*-
"""Tests of the input-validation guards (magnus.oscprob's validators).

These are the branches a user meets by getting an argument wrong, and a
coverage run found essentially none of them executed: 53 `raise` sites in
oscprob.py alone had never run. A guard that has never been taken is not
obviously a working guard -- writing these found one that could not fire at
all, because it caught the wrong exception type (see
test_earth_locations_of_the_wrong_shape_are_rejected_clearly).

What is asserted here is deliberately narrow: that the wrong input is
rejected, as ValueError rather than as whatever the interpreter happened to
raise, and that the message names the parameter at fault. Asserting the
full text would make these tests a transcription of the source, and they
would then fail on every reworded message rather than on a real change of
behaviour.
"""

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

ENERGY = 1.0*gd.UNIT_GEV
BASELINE = 1000.0*gd.UNIT_KM


# ----------------------------------------------------------------------
# unpack_*_params_from_dict: the dictionary routes into the Hamiltonians
# ----------------------------------------------------------------------

@pytest.mark.parametrize("num_flavors", [2, 3, 4, 5])
def test_incomplete_oscillation_parameter_dict_names_the_missing_keys(num_flavors):
    """Every flavor count needs its own set of keys, and the message has to
    say which -- the caller cannot guess that 4 flavors wants D41 as well."""
    with pytest.raises(ValueError, match="osc_params"):
        op.unpack_oscillation_params_from_dict('t', num_flavors, {'sth': 0.3}, None)


@pytest.mark.parametrize("num_flavors", [2, 3, 4, 5])
def test_incomplete_nsi_parameter_dict_names_the_missing_keys(num_flavors):
    with pytest.raises(ValueError, match="nsi_params"):
        op.unpack_nsi_params_from_dict('t', num_flavors, {'eps_ee': 0.1}, None)


@pytest.mark.parametrize("num_flavors", [2, 3, 4, 5])
def test_incomplete_liv_parameter_dict_names_the_missing_keys(num_flavors):
    with pytest.raises(ValueError, match="liv_params"):
        op.unpack_liv_params_from_dict('t', num_flavors, {'Lambda': 1.0}, None)


@pytest.mark.parametrize("num_flavors", [2, 3, 4, 5])
def test_non_positive_liv_scale_is_rejected(num_flavors):
    """Lambda is a scale that divides the LIV term, so zero or negative is
    not a physical choice; it is checked before the other keys are read, so
    the diagnostic is about Lambda rather than about a missing key."""
    liv_params = dict(Lambda=0.0, sxi=0.1, b1=1e-9, b2=1e-9, b3=1e-9, b4=1e-9, b5=1e-9,
                      n_liv=1, sxi12=0.1, sxi23=0.1, sxi13=0.1, dxiCP=0.0, dxi13=0.0,
                      sxi14=0.1, dxi14=0.0, sxi24=0.1, dxi24=0.0, sxi34=0.1,
                      sxi15=0.1, dxi15=0.0, sxi25=0.1, sxi35=0.1, dxi35=0.0)
    with pytest.raises(ValueError, match="Lambda"):
        op.unpack_liv_params_from_dict('t', num_flavors, liv_params, None)


@pytest.mark.parametrize("unpack, kind", [
    (op.unpack_oscillation_params_from_dict, 'oscillation'),
    (op.unpack_nsi_params_from_dict, 'nsi'),
    (op.unpack_liv_params_from_dict, 'liv'),
])
def test_more_flavors_than_predefined_requires_an_explicit_hamiltonian(unpack, kind):
    """Above five flavors there are no built-in Hamiltonians, so the caller
    has to supply the matrix itself; omitting it must be refused rather than
    silently producing a wrong-sized array."""
    with pytest.raises(ValueError):
        unpack('t', gd.MAGNUS_MAX_PREDEFINED_NUM_FLAVORS + 1, {}, None)


@pytest.mark.parametrize("unpack", [
    op.unpack_oscillation_params_from_dict,
    op.unpack_nsi_params_from_dict,
    op.unpack_liv_params_from_dict,
])
@pytest.mark.parametrize("num_flavors", [0, -3])
def test_fewer_than_one_flavor_is_rejected(unpack, num_flavors):
    with pytest.raises(ValueError):
        unpack('t', num_flavors, {}, None)


# ----------------------------------------------------------------------
# validate_input_osc_prob_earth: locations, zenith angle, baseline
# ----------------------------------------------------------------------

@pytest.mark.parametrize("loc_ini, loc_fin", [
    (None, (0.0, 0.0)),
    ((0.0, 0.0), None),
])
def test_a_single_earth_location_is_rejected(loc_ini, loc_fin):
    """A chord needs both ends. Accepting one would leave the other end
    silently defaulted, which is a wrong baseline rather than an error."""
    with pytest.raises(ValueError, match="loc_ini|loc_fin"):
        op.validate_input_osc_prob_earth('t', loc_ini=loc_ini, loc_fin=loc_fin)


@pytest.mark.parametrize("bad", [(1.0, 2.0, 3.0), 42])
def test_earth_locations_of_the_wrong_shape_are_rejected_clearly(bad):
    """Regression test. These two guards caught KeyError, which unpacking
    never raises: a three-entry tuple escaped as `too many values to unpack
    (expected 2)`, and a non-iterable as TypeError -- breaking the
    package-wide convention that bad input raises ValueError naming the
    parameter at fault.

    The message is asserted, not just the exception type. Checking only for
    ValueError would let the three-entry case pass against the old code,
    since `too many values to unpack` is itself a ValueError -- the test
    would then be green for a reason that has nothing to do with the guard
    working."""
    with pytest.raises(ValueError, match="loc_ini"):
        op.validate_input_osc_prob_earth('t', loc_ini=bad, loc_fin=(0.0, 0.0))

    with pytest.raises(ValueError, match="loc_fin"):
        op.validate_input_osc_prob_earth('t', loc_ini=((0, 0, 0), (0, 0, 0)), loc_fin=bad)


def test_a_single_earth_location_given_as_none_takes_the_costhz_path():
    """Both locations None is not an error -- it is the costhz route -- so
    it must not be swept up by the shape check above."""
    costhz, L = op.validate_input_osc_prob_earth('t', loc_ini=None, loc_fin=None,
                                                 costhz=-0.8, L=BASELINE)
    assert costhz == -0.8
    assert L == BASELINE


def test_earth_without_locations_requires_costhz():
    with pytest.raises(ValueError, match="costhz"):
        op.validate_input_osc_prob_earth('t', loc_ini=None, loc_fin=None, costhz=None,
                                         L=BASELINE)


def test_earth_with_costhz_requires_a_baseline():
    """costhz fixes the direction, not the distance."""
    with pytest.raises(ValueError, match="L|baseline"):
        op.validate_input_osc_prob_earth('t', loc_ini=None, loc_fin=None, costhz=-0.8, L=None)


def test_earth_with_two_locations_computes_the_chord():
    """The success path of the same function: two locations override costhz
    and define the baseline themselves. Coordinates are (degrees, minutes,
    seconds) per axis, the same form as earth.loc_coords_dms."""
    costhz, L = op.validate_input_osc_prob_earth(
        't', loc_ini=((51, 45, 54), (104, 24, 54)), loc_fin=((46, 14, 1.8), (6, 3, 11.4)),
        costhz=0.5, L=None)
    assert -1.0 <= costhz <= 1.0
    assert L > 0.0


# ----------------------------------------------------------------------
# validate_input_battery: the shared guard every wrapper delegates to
# ----------------------------------------------------------------------

def test_energy_must_be_a_number_or_a_flat_sequence():
    with pytest.raises(ValueError, match="energy"):
        op.validate_input_battery('t', energy='1 GeV', L=BASELINE)
    with pytest.raises(ValueError, match="energy"):
        op.validate_input_battery('t', energy=[[1.0, 2.0], [3.0, 4.0]], L=BASELINE)
    with pytest.raises(ValueError, match="energy"):
        op.validate_input_battery('t', energy=[1.0, 'two'], L=BASELINE)


def test_baseline_must_be_a_number_or_a_flat_sequence():
    with pytest.raises(ValueError, match="L"):
        op.validate_input_battery('t', energy=ENERGY, L='1000 km')
    with pytest.raises(ValueError, match="L"):
        op.validate_input_battery('t', energy=ENERGY, L=[[1.0], [2.0]])
    with pytest.raises(ValueError, match="L"):
        op.validate_input_battery('t', energy=ENERGY, L=[1.0, None])


def test_energy_and_baseline_arrays_must_be_the_same_length():
    """Paired arrays are consumed point by point, so mismatched lengths are
    a silent truncation rather than a broadcast."""
    with pytest.raises(ValueError):
        op.validate_input_battery('t', energy=[1.0, 2.0, 3.0], L=[1.0, 2.0])


def test_flavor_indices_must_be_given_together():
    with pytest.raises(ValueError, match="nu_i|nu_f"):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, nu_i=0, nu_f=None)
    with pytest.raises(ValueError, match="nu_i|nu_f"):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, nu_i=None, nu_f=1)


@pytest.mark.parametrize("nu_i, nu_f", [(-1, 0), (0, 99)])
def test_flavor_indices_must_be_in_range(nu_i, nu_f):
    with pytest.raises(ValueError, match="nu_i|nu_f"):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, num_flavors=3,
                                  nu_i=nu_i, nu_f=nu_f)


def test_oscillation_parameters_must_be_numbers():
    with pytest.raises(ValueError):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, osc_params=[0.3, 'big'])


def test_initial_position_must_be_a_number():
    """Only checked when the caller asks for it: the wrappers that take an
    L0 set validate_initial_position, the ones that do not, do not."""
    with pytest.raises(ValueError, match="L0"):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, osc_params=[0.3, 2.5e-3],
                                  L0='surface', validate_initial_position=True)


def test_matter_composition_ratios_must_be_non_negative():
    """Both are number ratios of constituents, so a negative value is not a
    physical composition; they are only checked under validate_density."""
    common = dict(energy=ENERGY, L=BASELINE, osc_params=[0.3, 2.5e-3], L0=0.0,
                  rho_func=3.0, validate_density=True)
    with pytest.raises(ValueError, match="ratio_number_neutrons_to_protons"):
        op.validate_input_battery('t', ratio_number_neutrons_to_protons=-1.0, **common)
    with pytest.raises(ValueError, match="electron_fraction"):
        op.validate_input_battery('t', electron_fraction=-0.5, **common)


def test_a_negative_matter_density_is_rejected():
    with pytest.raises(ValueError, match="rho|density"):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, osc_params=[0.3, 2.5e-3],
                                  L0=0.0, rho_func=-3.0, validate_density=True)


def test_a_density_profile_of_more_than_one_argument_is_rejected():
    """rho_func is called as rho_func(l) internally, so a two-argument
    callable would fail later, inside the integration, with no indication
    that the profile was the problem."""
    with pytest.raises(ValueError, match="rho_func"):
        op.validate_input_battery('t', energy=ENERGY, L=BASELINE, osc_params=[0.3, 2.5e-3],
                                  L0=0.0, rho_func=lambda l, extra: 3.0, validate_density=True)


def test_a_valid_battery_call_raises_nothing_and_returns_none():
    """The contract changed in 1.0.0: this used to return 1 on failure, and
    now raises instead, so success must be a plain None rather than a code
    the caller is expected to compare against."""
    assert op.validate_input_battery('t', energy=ENERGY, L=BASELINE, num_flavors=3,
                                     nu_i=0, nu_f=1, osc_params=[0.5, 0.6, 0.1, 0.0, 7.5e-5,
                                                                 2.5e-3]) is None


# ----------------------------------------------------------------------
# osc_prob's own guards
# ----------------------------------------------------------------------

def _flat_H(t):
    return np.array([[1.0e-12, 2.0e-13], [2.0e-13, -1.0e-12]], dtype=complex)


def test_osc_prob_rejects_a_backwards_interval():
    with pytest.raises(ValueError, match="t_fin"):
        op.osc_prob(_flat_H, t_ini=BASELINE, t_fin=0.0)


@pytest.mark.parametrize("order", [0, -2])
def test_osc_prob_rejects_an_out_of_range_expansion_order(order):
    """Zero and negative are meaningless, and silently clamping would return
    a different calculation from the one asked for."""
    with pytest.raises(ValueError, match="magnus_exp_order"):
        op.osc_prob(_flat_H, t_ini=0.0, t_fin=BASELINE, magnus_exp_order=order)


def test_osc_prob_rejects_an_order_above_the_implemented_maximum():
    """Above MAGNUS_EXP_ORDER_MAX there is no implementation at all. Checked
    with a quadrature method, because 'gl' refuses earlier and for its own
    reason: its schemes are separately derived integrators that stop at
    order 6 rather than continuing with the Magnus recursion."""
    with pytest.raises(ValueError, match="magnus_exp_order|order"):
        op.osc_prob(_flat_H, t_ini=0.0, t_fin=BASELINE, integration_method='trapezoid',
                    magnus_exp_order=gd.MAGNUS_EXP_ORDER_MAX + 1)


# ----------------------------------------------------------------------
# expansionterms' guards
# ----------------------------------------------------------------------

@pytest.mark.parametrize("order", [0, -1])
def test_omega_terms_rejects_orders_below_one(order):
    import magnus.expansionterms as et
    with pytest.raises(ValueError, match="order"):
        et.omega_terms(order)


@pytest.mark.parametrize("max_order", [0, -1])
def test_magnus_terms_rejects_orders_below_one(max_order):
    import magnus.expansionterms as et
    with pytest.raises(ValueError, match="max_order"):
        et.magnus_terms(max_order)


# ----------------------------------------------------------------------
# Predefined parameter sets carry labels, not just parameters
# ----------------------------------------------------------------------

@pytest.mark.parametrize("extra", [{}, {'average': True}])
def test_splatting_a_predefined_parameter_set_is_rejected_with_the_remedy(extra):
    """`**OSC_PARAMS_PREDEFINED['...']` is the natural thing to write and it
    does not work: the entries carry 'name' and 'description' strings
    alongside the six mixing parameters.

    Left unchecked, those two travel down the shared **kwargs chain until
    the Magnus core rejects them, naming the one function in the chain that
    has nothing to do with the mistake. Both the ordinary and the averaged
    path are checked here, because they diverge before the core is reached:
    the averaged one returns early, so a guard placed further down would
    catch the strings on one path and silently ignore them on the other --
    which is exactly what happened first time round."""
    predefined = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']

    with pytest.raises(ValueError, match="load_nufit_params|label"):
        op.osc_prob_3nu_vacuum(ENERGY, 1.0e8*gd.UNIT_KM, **extra, **predefined)


@pytest.mark.parametrize("extra", [{}, {'average': True}])
def test_load_nufit_params_is_the_form_that_works(extra):
    """The documented remedy has to actually run, on both paths."""
    osc = gd.load_nufit_params('NuFIT 6.1')

    P = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, 1.0e8*gd.UNIT_KM, **extra, **osc))

    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-9)


def test_only_the_labelling_keys_are_rejected():
    """The guard names two specific keys. An unrelated unexpected keyword is
    still the Magnus core's business to reject, and must not be swept into
    this diagnostic."""
    osc = gd.load_nufit_params('NuFIT 6.1')

    with pytest.raises(ValueError, match="name|description"):
        op.osc_prob_3nu_vacuum(ENERGY, BASELINE, name='mine', **osc)

    with pytest.raises((TypeError, ValueError)) as excinfo:
        op.osc_prob_3nu_vacuum(ENERGY, BASELINE, not_a_real_option=1, **osc)
    assert 'load_nufit_params' not in str(excinfo.value)


@pytest.mark.parametrize("bad", [
    {'max_n_slabs': 0},
    {'rtol': -1.0, 'atol': -1.0},
    {'max_num_loops': 0},
])
def test_a_constant_density_does_not_exempt_the_refinement_bounds(bad):
    """The batched constant engine answers before osc_prob is reached, and so before
    its validation runs.

    A bound osc_prob would have rejected therefore has to be rejected in the
    dispatcher too, or whether a bad parameter is reported at all depends on whether
    the caller's density happens to be constant -- the same request raising on a PREM
    profile and returning quietly on a uniform one.  The engine declines rather than
    raising, so the message the caller sees is still osc_prob's own.
    """
    osc = gd.load_nufit_params('NuFIT 6.1', 'NO')
    energy = np.linspace(1.0, 10.0, 6)*gd.UNIT_GEV
    with pytest.raises(ValueError):
        op.osc_prob_matter_std_potential(
            3, 2.848, energy, BASELINE, osc, L0=0.0,
            density_matter_is_in_g_per_cm3=True, **bad)


@pytest.mark.parametrize("bad, named", [
    ({'max_n_slabs': 1}, "max_n_slabs"),
    ({'max_n_tpts_per_slab': 2}, "max_n_tpts_per_slab"),
    ({'max_n_tpts_per_slab': 0}, "max_n_tpts_per_slab"),
])
def test_each_refinement_ceiling_is_checked_against_its_own_floor(bad, named):
    """`min_n_slabs` defaults to 1 and `min_n_tpts_per_slab` to 2, so each ceiling
    has to clear its own floor -- which is what the two messages say.

    One of the two conditions named the wrong variable: it tested `max_n_slabs`
    while reporting `max_n_tpts_per_slab`, so that parameter was never validated
    at all and `max_n_slabs` was bounded at > 2 while its own message promised
    > 1.  The error a caller got named a parameter they had not passed.
    """
    with pytest.raises(ValueError, match=named):
        op.osc_prob(np.diag([0.0, 1.0e-13, 2.0e-13]).astype(complex), 0.0, BASELINE,
                    rtol=1.0e-6, atol=1.0e-6, **bad)


def test_the_slab_ceiling_of_two_is_now_accepted():
    """The loosening half of the same fix: `max_n_slabs=2` clears a floor of 1 and
    always should have, but the mislabelled condition refused it."""
    P = op.osc_prob(np.diag([0.0, 1.0e-13, 2.0e-13]).astype(complex), 0.0, BASELINE,
                    rtol=1.0e-6, atol=1.0e-6, max_n_slabs=2)
    assert np.all(np.isfinite(np.asarray(P)))


@pytest.mark.parametrize('builder, args, bad_name', [
    ('hamiltonian_2nu_vacuum_energy_independent', (1.5, 2.5e-3), 'sth'),
    ('hamiltonian_3nu_vacuum_energy_independent',
     (1.5, 0.7, 0.15, 0.0, 7.4e-5, 2.5e-3), 's12'),
    ('hamiltonian_3nu_vacuum_energy_independent',
     (0.55, 0.7, -1.2, 0.0, 7.4e-5, 2.5e-3), 's13'),
    ('hamiltonian_4nu_vacuum_energy_independent',
     (0.55, 0.76, 0.15, 3.79, 2.0, 0.0, 0.22, 0.0, 0.0, 7.4e-5, 2.5e-3, 1.0), 's14'),
    ('hamiltonian_5nu_vacuum_energy_independent',
     (0.55, 0.76, 0.15, 3.79, 0.32, 0.0, 0.0, 0.0, 0.22, 0.0, 3.79, 0.0, 0.0, 0.0,
      7.4e-5, 2.5e-3, 1.0, 2.0), 's25'),
])
def test_a_sine_outside_the_unit_interval_is_rejected(builder, args, bad_name):
    """A cosine is built as sqrt(1 - s^2), and NumPy answers a bad sine with `nan` and
    a RuntimeWarning rather than an exception.

    The builder then returned a Hamiltonian full of `nan`, `osc_prob` propagated it into
    `nan` probabilities, and nothing raised -- a result-shaped object that is not a
    result, which is the silent-wrong-answer shape this package keeps being caught by.

    The mistake being caught is a *slot* error rather than a physics one: these
    signatures interleave each angle with its CP phase (`s14, d14, s15, d15, s24, d24,
    s25, ...`), so grouping the angles together as any reader would expect puts a phase
    -- 3.79 rad here -- into a sine slot.  That is how it was found.
    """
    import magnus.hamiltonians as hams
    with pytest.raises(ValueError, match=bad_name):
        getattr(hams, builder)(*args)


def test_the_guard_does_not_reject_legitimate_angles():
    """Including the boundary: maximal mixing is s = 1 exactly, and must be allowed."""
    import magnus.hamiltonians as hams
    hams.hamiltonian_2nu_vacuum_energy_independent(1.0, 2.5e-3)
    hams.hamiltonian_2nu_vacuum_energy_independent(-1.0, 2.5e-3)
    hams.hamiltonian_3nu_vacuum_energy_independent(0.55, 0.76, 0.15, 3.79,
                                                   7.4e-5, 2.5e-3)


# ----------------------------------------------------------------------
# The guards added by the pre-publish audit (items A1, B3, B4, B5)
#
# Each of these closed a path where a wrong argument produced an answer
# rather than a complaint, so what matters is that they fire at all -- an
# untested guard is indistinguishable from a guard that cannot fire, which
# is what the module docstring above records finding once already.
# ----------------------------------------------------------------------

@pytest.mark.parametrize("bad_energy", [-ENERGY, 0.0])
def test_non_positive_energy_is_refused(bad_energy):
    """E < 0 used to return the ANTINEUTRINO probability: unitary, in range,
    and an answer to a different question.  E = 0 returned NaN."""
    with pytest.raises(ValueError, match="energy"):
        op.osc_prob_3nu_vacuum(bad_energy, BASELINE,
                               **gd.load_nufit_params('NuFIT 6.1', 'NO'))


def test_a_negative_energy_hidden_in_an_array_is_refused():
    """The scalar case is the obvious one; a sign error usually arrives
    inside a scan, where one element of many has gone negative."""
    energies = np.array([1.0, -2.0, 3.0])*gd.UNIT_GEV
    with pytest.raises(ValueError, match="energy"):
        op.osc_prob_3nu_vacuum(energies, np.full(3, BASELINE),
                               **gd.load_nufit_params('NuFIT 6.1', 'NO'))


def test_a_positive_energy_and_the_antineutrino_flag_both_still_work():
    """The guard must not have closed the legitimate route to the answer
    that a negative energy used to return by accident."""
    osc = gd.load_nufit_params('NuFIT 6.1', 'NO')
    nu = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, BASELINE, **osc))
    nubar = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, BASELINE, nubar=True, **osc))
    np.testing.assert_allclose(nu.sum(axis=1), 1.0, atol=1.0e-12)
    np.testing.assert_allclose(nubar.sum(axis=1), 1.0, atol=1.0e-12)
    assert not np.allclose(nu, nubar)      # they are different physics


@pytest.mark.parametrize("kwargs", [
    {'n_slabs': 0}, {'n_slabs': -5}, {'min_n_slabs': 0}, {'n_tpts_per_slab': 1},
    {'min_n_tpts_per_slab': 1},
])
def test_non_positive_slab_and_sample_counts_are_refused(kwargs):
    """These were accepted and then ignored, so a typo looked like a setting
    that had been honoured."""
    key = list(kwargs)[0]
    with pytest.raises(ValueError, match=key):
        op.osc_prob(lambda t: np.zeros((3, 3), dtype=complex), 0.0, BASELINE, **kwargs)


@pytest.mark.parametrize("kwargs", [
    {'min_n_slabs': 100, 'max_n_slabs': 5},
    {'min_n_tpts_per_slab': 400, 'max_n_tpts_per_slab': 100},
])
def test_a_floor_above_its_own_ceiling_is_refused(kwargs):
    """A contradictory pair used to be answered, and which of the two bounds
    the ladder obeyed was an implementation detail."""
    with pytest.raises(ValueError, match="must be <="):
        op.osc_prob(lambda t: np.zeros((3, 3), dtype=complex), 0.0, BASELINE, **kwargs)


@pytest.mark.parametrize("bad_density", [float('nan'), float('inf')])
def test_a_non_finite_density_is_refused_and_not_called_a_unit_mistake(bad_density):
    """A NaN density used to be reported as a DensityUnitWarning saying it was
    "far too small to be in natural units" -- a confident diagnosis of the
    wrong problem, reached because `nan == 0.0` and `nan >= threshold` are
    both False and the guard fell through to its warning."""
    with pytest.raises(ValueError, match="finite"):
        op.osc_prob_3nu_matter_constant_density(
            ENERGY, BASELINE, bad_density,
            density_matter_is_in_g_per_cm3=True,
            **gd.load_nufit_params('NuFIT 6.1', 'NO'))


# ----------------------------------------------------------------------
# Zero NSI must be no NSI, at every flavour count
# ----------------------------------------------------------------------

@pytest.mark.parametrize("num_flavors", [3, 4, 5])
def test_nsi_with_zero_couplings_reproduces_the_standard_matter_case(num_flavors):
    """The NSI route built its standard piece as a literal `[1, 0, 0, 0]`
    diagonal, giving the sterile states zero where they carry
    -V_NC = (r/2) V_CC.  With every eps set to zero the two routes must be
    the same calculation, and they differed by 5.2e-02 at four flavours and
    5.1e-02 at five.  Three flavours is the control: no sterile state, so it
    agreed all along, which is why nothing noticed.

    This is the same omission that was found in the standard-potential path
    -- worth 0.29 in probability on a PREM chord -- surviving in the BSM
    one.  A fix applied per-call-site rather than at the shared definition
    leaves exactly this kind of remainder."""
    osc = dict(gd.load_nufit_params('NuFIT 6.1', 'NO'))
    eps_names = {
        3: ('eps_ee', 'eps_em', 'eps_et', 'eps_mm', 'eps_mt', 'eps_tt'),
        4: ('eps_ee', 'eps_em', 'eps_et', 'eps_es', 'eps_mm', 'eps_mt',
            'eps_ms', 'eps_tt', 'eps_ts', 'eps_ss'),
        5: ('eps_ee', 'eps_em', 'eps_et', 'eps_es1', 'eps_es2', 'eps_mm', 'eps_mt',
            'eps_ms1', 'eps_ms2', 'eps_tt', 'eps_ts1', 'eps_ts2',
            'eps_s1s1', 'eps_s1s2', 'eps_s2s2'),
    }[num_flavors]
    if num_flavors >= 4:
        osc.update(s14=0.3, s24=0.3, s34=0.0, d14=0.0, d24=0.0, D41=1.0)
    if num_flavors == 5:
        osc.update(s15=0.2, s25=0.0, s35=0.0, d15=0.0, d35=0.0, D51=2.0)

    std_fn = getattr(op, 'osc_prob_%dnu_matter_constant_density' % num_flavors)
    nsi_fn = getattr(op, 'osc_prob_%dnu_matter_nsi_constant_density' % num_flavors)
    common = dict(density_matter_is_in_g_per_cm3=True)

    std = np.asarray(std_fn(ENERGY, BASELINE, 3.0, **osc, **common))
    nsi = np.asarray(nsi_fn(ENERGY, BASELINE, 3.0, **osc,
                            **{k: 0.0 for k in eps_names}, **common))
    np.testing.assert_allclose(nsi, std, rtol=0.0, atol=1.0e-14)


def test_a_non_zero_nsi_coupling_still_changes_the_answer():
    """The guard above is satisfied by an NSI route that ignores eps
    entirely, so it needs this next to it."""
    osc = dict(gd.load_nufit_params('NuFIT 6.1', 'NO'))
    osc.update(s14=0.3, s24=0.3, s34=0.0, d14=0.0, d24=0.0, D41=1.0)
    eps = {k: 0.0 for k in ('eps_ee', 'eps_em', 'eps_et', 'eps_es', 'eps_mm',
                            'eps_mt', 'eps_ms', 'eps_tt', 'eps_ts', 'eps_ss')}
    common = dict(density_matter_is_in_g_per_cm3=True)
    std = np.asarray(op.osc_prob_4nu_matter_constant_density(
        ENERGY, BASELINE, 3.0, **osc, **common))
    nsi = np.asarray(op.osc_prob_4nu_matter_nsi_constant_density(
        ENERGY, BASELINE, 3.0, **osc, **{**eps, 'eps_ee': 0.1}, **common))
    assert np.max(np.abs(nsi - std)) > 1.0e-4
