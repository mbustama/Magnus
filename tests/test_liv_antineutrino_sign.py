# -*- coding: utf-8 -*-
"""The LIV term for antineutrinos (issue #162).

For antineutrinos the Lorentz-violating operator is transposed in flavor space -- conjugated,
being Hermitian -- and its CPT-odd part changes sign (Kostelecky & Mewes, Phys. Rev. D 85,
096005 (2012), Eqs. 77-78).  An operator of dimension d is CPT-odd for odd d and CPT-even for
even d; with n_liv = d - 3, the antineutrino term is (-1)**(n_liv + 1) times the conjugate of the
neutrino one.  On 1.1 it was only conjugated, whatever n_liv: right at odd n_liv, of the wrong
sign at even n_liv, and at two flavors it was not even told it was computing antineutrinos.
"""

import numpy as np
import pytest
from scipy.linalg import expm

import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.oscprob as op

LIV3 = dict(sxi12=0.3, sxi23=0.5, sxi13=0.2, dxiCP=1.1, b1=1e-13, b2=-2e-13, b3=3e-13, Lambda=1.0)
LIV4 = dict(LIV3, sxi14=0.1, dxi14=0.4, sxi24=0.2, dxi24=-0.7, sxi34=0.3, b4=-1e-13)
LIV5 = dict(LIV4, sxi15=0.1, dxi15=0.9, sxi25=0.2, sxi35=0.1, dxi35=-0.3, b5=2e-13)
BUILDERS = [(ham.hamiltonian_3nu_liv_energy_independent, ham.hamiltonian_3nu_liv, LIV3),
            (ham.hamiltonian_4nu_liv_energy_independent, ham.hamiltonian_4nu_liv, LIV4),
            (ham.hamiltonian_5nu_liv_energy_independent, ham.hamiltonian_5nu_liv, LIV5)]


@pytest.mark.parametrize('n_liv', [0, 1, 2, 3])
@pytest.mark.parametrize('indep, full, liv', BUILDERS)
def test_the_antineutrino_term_is_conjugated_and_its_cpt_odd_part_flipped(indep, full, liv, n_liv):
    sign = (-1)**(n_liv + 1)
    H = indep(**liv, n_liv=n_liv)
    assert np.max(np.abs(indep(**liv, n_liv=n_liv, nubar=True) - sign*np.conj(H))) == 0.0
    E = 3.0e9
    H = full(E, **liv, n_liv=n_liv)
    assert np.max(np.abs(full(E, **liv, n_liv=n_liv, nubar=True) - sign*np.conj(H))) == 0.0


@pytest.mark.parametrize('n_liv', [0, 1, 2])
def test_two_flavors_flip_the_cpt_odd_term_too(n_liv):
    liv = dict(sxi=0.4, b1=1e-23, b2=-2e-23, Lambda=1.0, n_liv=n_liv)
    sign = (-1)**(n_liv + 1)
    H = ham.hamiltonian_2nu_liv_energy_independent(**liv)
    assert np.array_equal(ham.hamiltonian_2nu_liv_energy_independent(**liv, nubar=True), sign*H)
    assert np.array_equal(ham.hamiltonian_2nu_liv(3.0e9, **liv, nubar=True),
                          sign*ham.hamiltonian_2nu_liv(3.0e9, **liv))


def _reference(H_vac_nu, H_liv_nu, n_liv, L):
    """The antineutrino probability built by hand from the neutrino Hamiltonians, as P[a, b] =
    P(a -> b) = |U[b, a]|**2."""
    H = np.conj(H_vac_nu) + (-1)**(n_liv + 1)*np.conj(H_liv_nu)
    U = expm(-1j*H*L)
    return (np.abs(U)**2).T


@pytest.mark.parametrize('n_liv', [0, 1])
def test_three_flavor_vacuum_antineutrino_probability(n_liv):
    osc = {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k]
           for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
    E, L = 1.0e9, 1000.0*gd.UNIT_KM
    liv = dict(LIV3, Lambda=1.0e9) if n_liv else LIV3
    P = np.asarray(op.osc_prob_3nu_vacuum_liv(E, L, **osc, **liv, n_liv=n_liv, nubar=True))
    ref = _reference(ham.hamiltonian_3nu_vacuum(E, **osc),
                     ham.hamiltonian_3nu_liv(E, **liv, n_liv=n_liv), n_liv, L)
    assert np.max(np.abs(P - ref)) < 1e-9
    # The effect is not negligible: the other sign is far away.
    wrong = _reference(ham.hamiltonian_3nu_vacuum(E, **osc),
                       -ham.hamiltonian_3nu_liv(E, **liv, n_liv=n_liv), n_liv, L)
    assert np.max(np.abs(P - wrong)) > 1e-2


def test_two_flavor_vacuum_antineutrino_probability():
    E, L = 1.0e9, 1000.0*gd.UNIT_KM
    vac = dict(sth=0.5, Dm2=2.5e-3)
    liv = dict(sxi=0.8, b1=0.0, b2=3e-13, Lambda=1.0)
    P = np.asarray(op.osc_prob_2nu_vacuum_liv(E, L, **vac, **liv, n_liv=0, nubar=True))
    ref = _reference(ham.hamiltonian_2nu_vacuum(E, **vac),
                     ham.hamiltonian_2nu_liv(E, **liv, n_liv=0), 0, L)
    assert np.max(np.abs(P - ref)) < 1e-9
    P_nu = np.asarray(op.osc_prob_2nu_vacuum_liv(E, L, **vac, **liv, n_liv=0))
    assert np.max(np.abs(P - P_nu)) > 1e-2


def test_the_cli_nubar_reaches_the_two_flavor_vacuum_liv_wrapper(capsys):
    """On 1.1, osc_prob_2nu_vacuum_liv did not declare nubar, so the CLI dropped --nubar."""
    import json
    from magnus import cli
    argv = ['prob', '--flavors', '2', '--environment', 'vacuum', '--scenario', 'liv',
            '--energy', '1', '--energy-unit', 'GeV', '--baseline', '1000', '--baseline-unit', 'km',
            '--sth', '0.5', '--dm2', '2.5e-3', '--sxi', '0.8', '--b2', '3e-13', '--json']
    P = []
    for extra in ([], ['--nubar']):
        assert cli.main(argv + extra) == 0
        P.append(np.array(json.loads(capsys.readouterr().out)['probability']))
    E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
    liv = dict(sxi=0.8, b1=0.0, b2=3e-13, Lambda=1.0, n_liv=0)
    assert np.max(np.abs(P[1] - np.asarray(op.osc_prob_2nu_vacuum_liv(
        E, L, 0.5, 2.5e-3, **liv, nubar=True)))) < 1e-12
    assert np.max(np.abs(P[1] - P[0])) > 1e-2
