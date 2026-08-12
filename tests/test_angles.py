# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""The ``angles`` convention: 'sin', 'sin2', 'rad', 'deg'.

Four ways of saying the same thing, so the tests that matter are the ones that check they
*are* the same thing -- at the builders, through the wrappers, and through the parameter
sources that fill in what the caller omitted.

The failure mode this whole feature has to avoid is the one this package keeps finding: a
plausible number in a slot that means something else, producing a converged, unitary,
entirely wrong probability rather than an error.  Two of those are pinned here explicitly,
with the size of the answer they would otherwise return.
"""

import warnings

import numpy as np
import pytest

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.oscprob as op
import magnus.oscprobstd as ops
from magnus.hamiltonians import _angles


P = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
TH = {k: float(np.arcsin(P[k])) for k in ('s12', 's23', 's13')}
D3 = dict(D21=P['D21'], D31=P['D31'])


def as_convention(theta, conv):
    """One angle, stated in `conv`."""
    if conv == 'sin':
        return float(np.sin(theta))
    if conv == 'sin2':
        return float(np.sin(theta))**2
    if conv == 'rad':
        return float(theta)
    return float(np.degrees(theta))


def phase(value, conv):
    """A CP phase follows the convention only for 'deg'."""
    return float(np.degrees(value)) if conv == 'deg' else float(value)


# ----------------------------------------------------------------------
# The four conventions describe one physical configuration
# ----------------------------------------------------------------------

@pytest.mark.parametrize('conv', gd.ANGLE_CONVENTIONS)
def test_the_three_flavour_hamiltonian_is_the_same_in_every_convention(conv):
    """The whole point.  'sin' is the reference because it is what the tables store."""
    ref = np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(
        P['s12'], P['s23'], P['s13'], P['dCP'], **D3))
    got = np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(
        *[as_convention(TH[k], conv) for k in ('s12', 's23', 's13')],
        phase(P['dCP'], conv), **D3, angles=conv))
    assert np.max(np.abs(got - ref)) < 1e-15, conv


@pytest.mark.parametrize('conv', gd.ANGLE_CONVENTIONS)
def test_a_probability_is_the_same_in_every_convention(conv):
    """Through the wrapper and its dispatcher, not just the builder.

    A builder that converts correctly proves nothing about whether the value survives the
    wrapper, the dispatcher and the defaults helper on the way down.
    """
    costhz = -0.8
    L = earth.distance_traveled_inside_earth(costhz)*gd.CONV_KM_TO_INV_EV
    kw = dict(energy=1.0*gd.UNIT_GEV, costhz=costhz, L=L, **D3)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ref = np.asarray(op.osc_prob_3nu_earth(
            s12=P['s12'], s23=P['s23'], s13=P['s13'], dCP=P['dCP'], **kw), dtype=float)
        got = np.asarray(op.osc_prob_3nu_earth(
            s12=as_convention(TH['s12'], conv), s23=as_convention(TH['s23'], conv),
            s13=as_convention(TH['s13'], conv), dCP=phase(P['dCP'], conv),
            angles=conv, **kw), dtype=float)
    # 'deg' cannot round-trip exactly -- degrees<->radians goes through pi -- so this is
    # not machine epsilon.  It is still far below any tolerance the package is asked for.
    assert np.max(np.abs(got - ref)) < 1e-9, conv


@pytest.mark.parametrize('conv', gd.ANGLE_CONVENTIONS)
def test_the_closed_form_two_flavour_probability_agrees_too(conv):
    """``oscprobstd`` consumes ``sth`` directly rather than through a builder."""
    ref = np.asarray(ops.osc_prob_2nu_vacuum_std(
        P['s12'], P['D21'], 1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM))
    got = np.asarray(ops.osc_prob_2nu_vacuum_std(
        as_convention(TH['s12'], conv), P['D21'], 1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM,
        angles=conv))
    assert np.max(np.abs(got - ref)) < 1e-15, conv


def test_the_default_path_is_untouched():
    """``angles='sin'`` must be the same object the caller passed, not a converted copy.

    The default is the overwhelmingly common path, and it should cost nothing and change
    nothing; this pins that it is a pass-through rather than a round trip.
    """
    out, _ = _angles.resolve('t', 'sin', {'s12': P['s12']})
    assert out['s12'] is P['s12']


# ----------------------------------------------------------------------
# The guards
# ----------------------------------------------------------------------

def test_an_unknown_convention_is_refused_and_the_message_lists_the_four():
    with pytest.raises(ValueError, match="angles must be one of"):
        hams.hamiltonian_3nu_vacuum_energy_independent(
            P['s12'], P['s23'], P['s13'], P['dCP'], **D3, angles='degrees')


@pytest.mark.parametrize('bad', [1.4, -1.4])
def test_a_sine_outside_the_unit_interval_is_refused(bad):
    with pytest.raises(ValueError, match='s12'):
        hams.hamiltonian_3nu_vacuum_energy_independent(
            bad, P['s23'], P['s13'], P['dCP'], **D3)


def test_a_negative_sin2_is_refused_because_it_means_the_caller_is_in_sin():
    """``sin2`` cannot express a negative sine; a negative value is a convention error."""
    with pytest.raises(ValueError, match='SQUARE'):
        hams.hamiltonian_3nu_vacuum_energy_independent(
            -0.3088, 0.47, 0.02248, P['dCP'], **D3, angles='sin2')


def test_degrees_in_a_radians_slot_are_refused():
    """33.4 is a perfectly ordinary angle in degrees and impossible in radians."""
    with pytest.raises(ValueError, match=r'radians'):
        hams.hamiltonian_3nu_vacuum_energy_independent(
            33.4, 49.0, 8.5, P['dCP'], **D3, angles='rad')


def test_sines_in_a_degrees_slot_warn_because_no_bound_can_catch_them():
    """Every angle under one degree is not a small-angle study, it is sines.

    theta_13, the smallest angle anyone measures, is about 8.5 degrees.  A warning rather
    than an error, following ``matter.DensityUnitWarning``: the threshold describes the
    mixing people currently study, not a law.
    """
    with pytest.warns(gd.MixingAngleConventionWarning, match='fifty times too small'):
        hams.hamiltonian_3nu_vacuum_energy_independent(
            P['s12'], P['s23'], P['s13'], P['dCP'], **D3, angles='deg')


def test_ordinary_degrees_do_not_warn():
    """The guard must not fire on the case it exists to protect."""
    with warnings.catch_warnings():
        warnings.simplefilter('error', gd.MixingAngleConventionWarning)
        hams.hamiltonian_3nu_vacuum_energy_independent(
            33.76, 43.28, 8.62, 212.0, **D3, angles='deg')


# ----------------------------------------------------------------------
# The parameter sources, which supply what the caller did not
# ----------------------------------------------------------------------

def test_load_nufit_params_reproduces_the_published_table():
    """NuFit quotes sin^2 and degrees; the loader must give back exactly those.

    An independent check: these numbers come from the publication, not from inverting
    Magnus's own stored sines with Magnus's own arithmetic.
    """
    sin2 = gd.load_nufit_params('NuFIT 6.1', 'NO', angles='sin2')
    assert sin2['s12'] == pytest.approx(0.3088, abs=5e-5)
    assert sin2['s13'] == pytest.approx(0.02248, abs=5e-6)
    deg = gd.load_nufit_params('NuFIT 6.1', 'NO', angles='deg')
    assert deg['s12'] == pytest.approx(33.76, abs=5e-3)
    assert deg['s13'] == pytest.approx(8.62, abs=5e-3)
    assert deg['dCP'] == pytest.approx(212.0, abs=5e-3)


@pytest.mark.parametrize('conv', gd.ANGLE_CONVENTIONS)
def test_the_loader_round_trips_through_every_convention(conv):
    """Its output, handed straight back with the same ``angles``, must be the same physics."""
    ref = np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(
        **gd.load_nufit_params('NuFIT 6.1', 'NO')))
    pars = gd.load_nufit_params('NuFIT 6.1', 'NO', angles=conv)
    got = np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(**pars, angles=conv))
    assert np.max(np.abs(got - ref)) < 1e-15, conv


def test_omitting_an_angle_does_not_mix_conventions_within_one_call():
    """The predefined sets are stored as SINES, and a caller may not be in sines.

    Filling one omitted angle from the stored value while the others arrived as degrees
    would hand the builder a parameter set in two conventions at once -- and it converts
    them all alike, so 0.1499 would be read as 0.1499 degrees.  Wrong only for the
    parameters the caller happened to leave out.

    The size of the answer that would otherwise come back is pinned below, because "it is
    subtly different" and "it is 14% different" are different bugs.
    """
    costhz = -0.8
    L = earth.distance_traveled_inside_earth(costhz)*gd.CONV_KM_TO_INV_EV
    deg = {k: float(np.degrees(TH[k])) for k in TH}
    kw = dict(energy=1.0*gd.UNIT_GEV, costhz=costhz, L=L, **D3,
              dCP=float(np.degrees(P['dCP'])), angles='deg')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        explicit = np.asarray(op.osc_prob_3nu_earth(
            s12=deg['s12'], s23=deg['s23'], s13=deg['s13'], **kw), dtype=float)
        omitted = np.asarray(op.osc_prob_3nu_earth(
            s12=deg['s12'], s23=deg['s23'], **kw), dtype=float)
        # what the un-fixed behaviour produced: the stored SINE left in a degrees call
        mixed = np.asarray(op.osc_prob_3nu_earth(
            s12=deg['s12'], s23=deg['s23'], s13=P['s13'], **kw), dtype=float)

    assert np.max(np.abs(omitted - explicit)) == 0.0
    assert np.max(np.abs(mixed - explicit)) > 0.1


# ----------------------------------------------------------------------
# Structural agreement: signature, docstring and body
# ----------------------------------------------------------------------

def test_every_public_function_taking_angles_documents_and_uses_it():
    """Three texts that have to agree, across ~95 functions.

    Adding this parameter went wrong in exactly these two ways during the sweep: one
    function got the argument but not the Parameters entry, and four accepted it and
    silently ignored it -- which no test, docs build or linter could see.
    """
    import ast
    import inspect
    import pathlib

    missing_doc, inert = [], []
    for mod in (hams, op, ops, gd):
        for name in dir(mod):
            if name.startswith('_'):
                continue
            fn = getattr(mod, name)
            if not callable(fn) or isinstance(fn, type):
                continue
            try:
                sig = inspect.signature(fn)
            except (TypeError, ValueError):
                continue
            if 'angles' in sig.parameters and 'angles : ' not in (inspect.getdoc(fn) or ''):
                missing_doc.append('%s.%s' % (mod.__name__, name))

    root = pathlib.Path(__file__).resolve().parent.parent/'src'/'magnus'
    for path in sorted(root.rglob('*.py')):
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.FunctionDef):
                continue
            if 'angles' not in [a.arg for a in node.args.args + node.args.kwonlyargs]:
                continue
            reads = any(isinstance(sub, ast.Name) and sub.id == 'angles'
                        for sub in ast.walk(node))
            if not reads:
                inert.append('%s:%s' % (path.name, node.name))

    assert not missing_doc, 'take `angles` but do not document it: %s' % missing_doc
    assert not inert, 'accept `angles` and never use it: %s' % inert
