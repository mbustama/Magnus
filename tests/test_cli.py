# -*- coding: utf-8 -*-
"""Smoke tests for the magnus command-line calculator (magnus.cli).

These call magnus.cli.main() directly (not via subprocess/the installed
console script), so they run against the same in-repo source as the rest
of the suite without requiring `pip install -e .` first.
"""

import json

import numpy as np
import pytest

import magnus.cli as cli


def run(argv, capsys):
    exit_code = cli.main(argv)
    captured = capsys.readouterr()
    return exit_code, captured.out


def test_vacuum_3nu_prints_unitary_matrix(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'vacuum',
                      '--energy', '1', '--energy-unit', 'GeV',
                      '--baseline', '1300', '--baseline-unit', 'km'], capsys)
    assert code == 0
    assert 'osc_prob_3nu_vacuum' in out


def test_vacuum_3nu_json_matches_direct_call(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'vacuum',
                      '--energy', '1', '--energy-unit', 'GeV',
                      '--baseline', '1300', '--baseline-unit', 'km', '--json'], capsys)
    assert code == 0
    payload = json.loads(out)
    P = np.asarray(payload['probability'])
    assert np.allclose(np.sum(P, axis=-1), 1.0, atol=1e-8)

    import magnus.oscprob as oscprob
    import magnus.globaldefs as gd
    P_direct = oscprob.osc_prob_3nu_vacuum(1.0 * gd.UNIT_GEV, 1300.0 * gd.UNIT_KM)
    assert np.allclose(P, P_direct)


def test_matter_constant_density_2nu_requires_sth_dm2(capsys):
    with pytest.raises(SystemExit):
        cli.main(['prob', '--flavors', '2', '--environment', 'vacuum',
                  '--energy', '1', '--baseline', '1300'])


def test_matter_requires_rho(capsys):
    with pytest.raises(SystemExit):
        cli.main(['prob', '--flavors', '3', '--environment', 'matter',
                  '--energy', '1', '--baseline', '1000'])


def test_vacuum_nsi_is_rejected(capsys):
    with pytest.raises(SystemExit):
        cli.main(['prob', '--flavors', '3', '--environment', 'vacuum', '--scenario', 'nsi',
                  '--energy', '1', '--baseline', '1300'])


def test_earth_costhz_without_baseline_is_rejected(capsys):
    with pytest.raises(SystemExit):
        cli.main(['prob', '--flavors', '3', '--environment', 'earth',
                  '--energy', '1', '--costhz', '-0.8'])


def test_earth_with_locations_computes_baseline_automatically(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'earth',
                      '--energy', '1', '--energy-unit', 'GeV',
                      '--loc-ini', 'fermilab', '--loc-fin', 'homestake'], capsys)
    assert code == 0
    assert 'osc_prob_3nu_earth' in out


def test_matter_nsi_constant_density_dispatch(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'matter', '--scenario', 'nsi',
                      '--rho', '2.7', '--eps-ee', '0.06', '--eps-em', '-0.06',
                      '--energy', '1', '--energy-unit', 'GeV',
                      '--baseline', '1000', '--baseline-unit', 'km'], capsys)
    assert code == 0
    assert 'osc_prob_3nu_matter_nsi_constant_density' in out


def test_vacuum_liv_dispatch_and_effect(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'vacuum', '--scenario', 'liv',
                      '--sxi12', '0.3', '--b1', '6e-13', '--b2', '1.2e-12',
                      '--liv-lambda', '1e9', '--n-liv', '1',
                      '--energy', '1', '--energy-unit', 'GeV',
                      '--baseline', '1300', '--baseline-unit', 'km', '--json'], capsys)
    assert code == 0
    payload = json.loads(out)
    P_liv = np.asarray(payload['probability'])

    import magnus.oscprob as oscprob
    import magnus.globaldefs as gd
    P_vac = oscprob.osc_prob_3nu_vacuum(1.0 * gd.UNIT_GEV, 1300.0 * gd.UNIT_KM)
    assert not np.allclose(P_liv, P_vac)
    assert np.allclose(np.sum(P_liv, axis=-1), 1.0, atol=1e-8)


def test_channel_selection_returns_scalar(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'vacuum',
                      '--energy', '1', '--energy-unit', 'GeV',
                      '--baseline', '1300', '--baseline-unit', 'km',
                      '--nu-i', 'e', '--nu-f', 'mu'], capsys)
    assert code == 0
    assert 'P =' in out


def test_5nu_earth_liv_sterile_states_decouple_when_zero(capsys):
    code, out = run(['prob', '--flavors', '5', '--environment', 'earth', '--scenario', 'liv',
                      '--costhz', '-0.8', '--baseline', '10193.6',
                      '--sxi12', '0.2', '--b1', '1e-13', '--liv-lambda', '1e9',
                      '--energy', '1', '--energy-unit', 'GeV', '--json'], capsys)
    assert code == 0
    payload = json.loads(out)
    P = np.asarray(payload['probability'])
    # s14=s15=s24=s25=s34=s35=0.0 (defaults): the two sterile states must stay
    # perfectly decoupled from the three active flavors and from each other.
    assert np.allclose(P[3:, :3], 0.0)
    assert np.allclose(P[:3, 3:], 0.0)
    assert np.allclose(P[3, 3], 1.0)
    assert np.allclose(P[4, 4], 1.0)


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc_info:
        cli.main(['--version'])
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    from magnus.version import __version__
    assert __version__ in out
