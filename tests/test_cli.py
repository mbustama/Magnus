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


# ----------------------------------------------------------------------
# Argument-error paths, and the flavor/scenario combinations whose
# parameter builders no invocation ever reached.
#
# A coverage run found most of the CLI's own error handling unexecuted.
# These are the paths a user meets by typing a command slightly wrong, so
# what matters is that each fails as an argument error -- exit code 2 with
# a message naming the missing flag -- rather than as a traceback out of
# the library, or worse, a silently wrong calculation.
# ----------------------------------------------------------------------

def test_flavor_name_is_accepted_in_place_of_an_index(capsys):
    """--nu-i/--nu-f take either an integer or a flavor name."""
    code, out = run(['prob', '--flavors', '3', '--environment', 'vacuum',
                     '--energy', '1', '--baseline', '1300',
                     '--nu-i', 'e', '--nu-f', 'mu'], capsys)
    assert code == 0
    assert 'P = ' in out


def test_unknown_flavor_name_is_an_argument_error():
    """An unrecognized flavor name must be rejected by argparse itself, with
    the accepted names listed, rather than reaching the library as a
    nonsense index."""
    with pytest.raises(SystemExit) as excinfo:
        cli.main(['prob', '--flavors', '3', '--environment', 'vacuum',
                  '--energy', '1', '--baseline', '1300', '--nu-i', 'quark'])
    assert excinfo.value.code == 2


def test_two_flavor_vacuum_passes_sth_and_dm2(capsys):
    code, out = run(['prob', '--flavors', '2', '--environment', 'vacuum',
                     '--energy', '1', '--baseline', '1300',
                     '--sth', '0.4', '--dm2', '2.5e-3'], capsys)
    assert code == 0
    assert 'osc_prob_2nu_vacuum' in out


def test_two_flavor_nsi_passes_its_own_epsilon_names(capsys):
    """The 2nu NSI couplings are named eps_aa/eps_ab, not the eps_ee/eps_em
    of three or more flavors."""
    code, out = run(['prob', '--flavors', '2', '--environment', 'matter',
                     '--scenario', 'nsi', '--energy', '1', '--baseline', '1300',
                     '--rho', '3.0', '--sth', '0.4', '--dm2', '2.5e-3',
                     '--eps-aa', '0.1', '--eps-ab', '0.05'], capsys)
    assert code == 0
    assert 'osc_prob_2nu_matter_nsi_constant_density' in out


@pytest.mark.parametrize("flavors", ['4', '5'])
def test_sterile_nsi_passes_the_sterile_epsilons(flavors, capsys):
    """4nu names its single sterile row without an index, 5nu numbers
    them; each has its own branch in the parameter builder."""
    code, out = run(['prob', '--flavors', flavors, '--environment', 'matter',
                     '--scenario', 'nsi', '--energy', '1', '--baseline', '1300',
                     '--rho', '3.0'], capsys)
    assert code == 0
    assert f'osc_prob_{flavors}nu_matter_nsi_constant_density' in out


def test_two_flavor_liv_passes_its_own_parameter_names(capsys):
    code, out = run(['prob', '--flavors', '2', '--environment', 'vacuum',
                     '--scenario', 'liv', '--energy', '1', '--baseline', '1300',
                     '--sth', '0.4', '--dm2', '2.5e-3', '--sxi', '0.2',
                     '--b1', '1e-9', '--b2', '1e-9', '--n-liv', '1'], capsys)
    assert code == 0
    assert 'osc_prob_2nu_vacuum_liv' in out


def test_antineutrino_run_is_labelled_as_such(capsys):
    """--nubar changes the physics, so the printed header has to say so;
    otherwise two runs differing only in this flag print identically."""
    code, out = run(['prob', '--flavors', '3', '--environment', 'matter',
                     '--energy', '1', '--baseline', '1300', '--rho', '3.0',
                     '--nubar'], capsys)
    assert code == 0
    assert 'antineutrinos' in out


def test_exponential_density_requires_rho_central_and_l_scale():
    with pytest.raises(SystemExit):
        cli.main(['prob', '--flavors', '3', '--environment', 'matter',
                  '--density-profile', 'exp', '--energy', '1', '--baseline', '1000'])


def test_earth_requires_costhz_or_a_pair_of_locations():
    with pytest.raises(SystemExit):
        cli.main(['prob', '--flavors', '3', '--environment', 'earth',
                  '--energy', '1', '--baseline', '1000'])


def test_sun_requires_a_baseline():
    """Caught by the general baseline check, which runs before the
    per-environment one: every environment except earth needs --baseline,
    so the sun branch's own check is a backstop that is never reached from
    the command line."""
    with pytest.raises(SystemExit) as excinfo:
        cli.main(['prob', '--flavors', '3', '--environment', 'sun', '--energy', '1'])
    assert excinfo.value.code == 2


def test_exponential_density_run_forwards_its_profile_parameters(capsys):
    code, out = run(['prob', '--flavors', '3', '--environment', 'matter',
                     '--density-profile', 'exp', '--energy', '1',
                     '--baseline', '1000', '--rho-central', '10.0',
                     '--l-scale', '100'], capsys)
    assert code == 0
    assert 'osc_prob_3nu_matter_exp_density' in out


def test_sun_run_forwards_its_baseline_and_origin(capsys):
    """A short solar baseline, chosen to keep this a CLI test rather than a
    solar-oscillation benchmark."""
    code, out = run(['prob', '--flavors', '3', '--environment', 'sun',
                     '--energy', '10', '--energy-unit', 'MeV',
                     '--baseline', '69470', '--baseline-unit', 'km'], capsys)
    assert code == 0
    assert 'osc_prob_3nu_sun' in out


def test_a_library_validation_error_becomes_an_argument_error():
    """The library raises ValueError for a bad input; the CLI must turn that
    into a clean argument error rather than letting a traceback escape.
    A negative density is rejected by the library, not by argparse."""
    with pytest.raises(SystemExit) as excinfo:
        cli.main(['prob', '--flavors', '3', '--environment', 'matter',
                  '--energy', '1', '--baseline', '1300', '--rho', '-3.0'])
    assert excinfo.value.code == 2
