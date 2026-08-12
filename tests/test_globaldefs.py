# -*- coding: utf-8 -*-
"""Tests for the historical NuFit data in magnus.globaldefs.

Covers NUFIT_GLOBAL_FITS structure/coverage and the load_nufit_params()
loader: version/ordering/category validation, default category selection,
mass-ordering sign conventions, and direct usability with oscprob wrappers.
"""

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as oscprob


EXPECTED_VERSIONS = [
    'NuFIT 1.0', 'NuFIT 1.1', 'NuFIT 1.2', 'NuFIT 1.3',
    'NuFIT 2.0', 'NuFIT 2.1', 'NuFIT 2.2',
    'NuFIT 3.0', 'NuFIT 3.1', 'NuFIT 3.2',
    'NuFIT 4.0', 'NuFIT 4.1',
    'NuFIT 5.0', 'NuFIT 5.1', 'NuFIT 5.2', 'NuFIT 5.3',
    'NuFIT 6.0', 'NuFIT 6.1',
]

PARAM_KEYS = {'s12', 's23', 's13', 'dCP', 'D21', 'D31'}


def test_all_expected_versions_present():
    assert set(gd.NUFIT_GLOBAL_FITS.keys()) == set(EXPECTED_VERSIONS)


@pytest.mark.parametrize('version', EXPECTED_VERSIONS)
def test_every_version_has_no_and_io_for_every_category(version):
    entry = gd.NUFIT_GLOBAL_FITS[version]
    assert entry['categories'], f'{version} has no categories'
    for category, params_by_ordering in entry['categories'].items():
        assert set(params_by_ordering.keys()) == {'NO', 'IO'}
        for ordering in ('NO', 'IO'):
            params = params_by_ordering[ordering]
            assert set(params.keys()) == PARAM_KEYS
            for key in PARAM_KEYS:
                assert np.isfinite(params[key])
            assert 0.0 < params['s12'] < 1.0
            assert 0.0 < params['s23'] < 1.0
            assert 0.0 < params['s13'] < 1.0
            assert params['D21'] > 0.0


@pytest.mark.parametrize('version', EXPECTED_VERSIONS)
def test_d31_sign_matches_ordering(version):
    for category, params_by_ordering in gd.NUFIT_GLOBAL_FITS[version]['categories'].items():
        assert params_by_ordering['NO']['D31'] > 0.0
        assert params_by_ordering['IO']['D31'] < 0.0


def test_nufit_6_0_with_sk_matches_existing_module_constants():
    # OSC_PARAMS_NU_FIT_6_0_SK_{NO,IO} predate this dict and are the values
    # already relied upon elsewhere (e.g. OSC_PARAMS_DEFAULT); the new
    # historical dict must reproduce them exactly for NuFIT 6.0/with_SK.
    no = gd.load_nufit_params('NuFIT 6.0', ordering='NO', category='with_SK')
    io = gd.load_nufit_params('NuFIT 6.0', ordering='IO', category='with_SK')
    for key in PARAM_KEYS:
        assert no[key] == pytest.approx(gd.OSC_PARAMS_NU_FIT_6_0_SK_NO[key])
        assert io[key] == pytest.approx(gd.OSC_PARAMS_NU_FIT_6_0_SK_IO[key])


def test_load_nufit_params_default_ordering_and_category():
    params = gd.load_nufit_params('NuFIT 6.1')
    assert set(params.keys()) == PARAM_KEYS
    assert params['D31'] > 0.0  # default ordering is 'NO'
    assert params == gd.load_nufit_params('NuFIT 6.1', ordering='NO', category='with_SK')


def test_load_nufit_params_category_specific_defaults():
    # Version families with a with_SK/without_SK split default to 'with_SK'.
    assert (gd.load_nufit_params('NuFIT 4.0', ordering='NO')
            == gd.load_nufit_params('NuFIT 4.0', ordering='NO', category='with_SK'))
    # v1.x families default to 'free_fluxes_rsbl'.
    assert (gd.load_nufit_params('NuFIT 1.0', ordering='NO')
            == gd.load_nufit_params('NuFIT 1.0', ordering='NO', category='free_fluxes_rsbl'))
    # v2.1 defaults to 'LEM'.
    assert (gd.load_nufit_params('NuFIT 2.1', ordering='NO')
            == gd.load_nufit_params('NuFIT 2.1', ordering='NO', category='LEM'))


def test_load_nufit_params_legacy_versions_share_angles_across_ordering():
    # v1.0-v1.3 predate separate per-ordering fits for the mixing angles.
    no = gd.load_nufit_params('NuFIT 1.2', ordering='NO')
    io = gd.load_nufit_params('NuFIT 1.2', ordering='IO')
    for key in ('s12', 's23', 's13', 'dCP'):
        assert no[key] == io[key]
    assert no['D31'] != io['D31']


def test_load_nufit_params_rejects_unknown_version():
    with pytest.raises(ValueError):
        gd.load_nufit_params('NuFIT 99.9')


def test_load_nufit_params_rejects_bad_ordering():
    with pytest.raises(ValueError):
        gd.load_nufit_params('NuFIT 6.1', ordering='NX')


def test_load_nufit_params_rejects_bad_category():
    with pytest.raises(ValueError):
        gd.load_nufit_params('NuFIT 6.1', category='not_a_real_category')


def test_load_nufit_params_output_feeds_directly_into_oscprob():
    for version, ordering in [('NuFIT 6.1', 'NO'), ('NuFIT 6.1', 'IO'),
                               ('NuFIT 1.0', 'NO'), ('NuFIT 3.0', 'IO')]:
        params = gd.load_nufit_params(version, ordering=ordering)
        P = oscprob.osc_prob_3nu_vacuum(1.0 * gd.UNIT_GEV, 1300.0 * gd.UNIT_KM,
                                         **params)
        assert np.all(np.isfinite(P))
        assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-8)


def test_load_nufit_params_returns_independent_copy():
    p1 = gd.load_nufit_params('NuFIT 6.1')
    p1['s12'] = -999.0
    p2 = gd.load_nufit_params('NuFIT 6.1')
    assert p2['s12'] != -999.0


# ----------------------------------------------------------------------
# set_color_output: the supported way to suppress ANSI escapes.
#
# Every call site in the package references the *_IN_COLOR names, so this
# function rebinding them is the only mechanism there is for turning
# colour off. It had no test and no caller, while the docs build worked
# around its absence by assigning gd.WARNING_MSG_IN_COLOR by hand to stop
# raw escape codes appearing in the rendered HTML.
# ----------------------------------------------------------------------

def test_set_color_output_round_trips():
    """Disabling must give exactly the plain-text constants, and enabling
    again must restore escapes -- a one-way switch would leave any process
    that turned colour off unable to get it back."""
    try:
        gd.set_color_output(False)
        assert gd.WARNING_MSG_IN_COLOR == gd.WARNING_MSG_NO_COLOR
        assert gd.ERROR_MSG_IN_COLOR == gd.ERROR_MSG_NO_COLOR
        assert gd.TOL_MSG_IN_COLOR == gd.TOL_MSG_NO_COLOR
        for msg in [gd.WARNING_MSG_IN_COLOR, gd.ERROR_MSG_IN_COLOR, gd.TOL_MSG_IN_COLOR]:
            assert '\x1b[' not in msg

        gd.set_color_output(True)
        for msg in [gd.WARNING_MSG_IN_COLOR, gd.ERROR_MSG_IN_COLOR, gd.TOL_MSG_IN_COLOR]:
            assert '\x1b[' in msg
    finally:
        # Whatever this test does, the rest of the session sees the default.
        gd.set_color_output(True)


def test_warnings_respect_set_color_output(recwarn):
    """The point of the switch is the message text actually emitted, not
    just the constant: a call site that had captured the coloured string at
    import time would keep printing escapes regardless."""
    import warnings

    import magnus.oscprob as op

    try:
        gd.set_color_output(False)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM, verbose=0)
            # A warning is not guaranteed here; what must never happen is an
            # escape code in one that is emitted.
            for w in caught:
                assert '\x1b[' not in str(w.message)
    finally:
        gd.set_color_output(True)


# ----------------------------------------------------------------------
# The default oscillation parameters are a contract
# ----------------------------------------------------------------------

def test_the_implicit_default_is_the_release_the_loader_returns():
    """There are two ways to get "the default parameters", and they used to
    disagree.  Omitting them fell back to constants named `*_BF_NUFIT_6_0`,
    while `load_nufit_params()` with no arguments returned 6.1 -- so the same
    script got different answers depending on which door it came through, by
    4.0e-03 in probability at 1 GeV over 1300 km.  Both were documented, which
    is why neither read as a mistake.

    Nothing pinned it: changing the fallback release broke no test in the
    suite, which is how the two drifted apart in the first place."""
    default = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    loaded = gd.load_nufit_params()          # no arguments: the documented default
    for key, value in loaded.items():
        assert default[key] == value, (
            'OSC_PARAMS_DEFAULT[%r] is %r but load_nufit_params() gives %r'
            % (key, default[key], value))


def test_the_default_and_the_loader_give_the_same_probability():
    """The contract that matters to a caller is not the dictionary, it is the
    number that comes out."""
    import magnus.oscprob as op
    energy, baseline = 1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM
    implicit = np.asarray(op.osc_prob_3nu_vacuum(energy, baseline))
    explicit = np.asarray(op.osc_prob_3nu_vacuum(energy, baseline,
                                                 **gd.load_nufit_params()))
    np.testing.assert_allclose(implicit, explicit, rtol=0.0, atol=0.0)


def test_the_superseded_release_is_still_reachable_by_name():
    """Changing which release is the default must not remove the old one: a
    caller pinned to 6.0 for reproducibility asks for it explicitly."""
    for name in ('OSC_PARAMS_NU_FIT_6_0_SK_NO', 'OSC_PARAMS_NU_FIT_6_0_SK_IO',
                 'OSC_PARAMS_NU_FIT_6_1_SK_NO', 'OSC_PARAMS_NU_FIT_6_1_SK_IO'):
        assert name in gd.OSC_PARAMS_PREDEFINED
    old = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_NU_FIT_6_0_SK_NO']
    assert old['D31'] == pytest.approx(2.513e-3)
