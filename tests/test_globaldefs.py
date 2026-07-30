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
