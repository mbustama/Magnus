# -*- coding: utf-8 -*-
"""Tests of the Hamiltonian builders (magnus.hamiltonians).

This module had no dedicated regression coverage before the audit that
found the bugs these tests guard against -- notably two mixing-matrix
formula bugs (4nu and 5nu) that invalidated every sterile-neutrino
calculation in the package. Every test here was derived and verified
independently of the fix (via unitarity, an independent SymPy
re-derivation, or a Hermiticity/reference-formula check), not merely by
asserting the code equals itself.
"""

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.hamiltonians.hamiltonians2nu as h2
import magnus.hamiltonians.hamiltonians3nu as h3
import magnus.hamiltonians.hamiltonians4nu as h4
import magnus.hamiltonians.hamiltonians5nu as h5

RNG = np.random.default_rng(11)


def rand_angles(n, lo=0.05, hi=0.5):
    return RNG.uniform(lo, hi, n)


def rand_phases(n, lo=0.1, hi=6.0):
    return RNG.uniform(lo, hi, n)


def maxabs(x):
    return np.max(np.abs(x))


def is_unitary(U, tol=1e-10):
    d = U.shape[0]
    return maxabs(U @ U.conj().T - np.eye(d)) < tol


# ----------------------------------------------------------------------
# A1: mixing_matrix_4x4
# ----------------------------------------------------------------------

@pytest.mark.parametrize("trial", range(8))
def test_mixing_matrix_4x4_unitary_and_paths_agree(trial):
    s12, s23, s13, s14, s24, s34 = rand_angles(6)
    d13, d14, d24 = rand_phases(3)
    args = (s12, s23, s13, d13, s14, d14, s24, d24, s34)
    U_fast = h4.mixing_matrix_4x4(*args, compute_matrix_multiplication=False)
    U_slow = h4.mixing_matrix_4x4(*args, compute_matrix_multiplication=True)
    assert is_unitary(U_fast)
    assert is_unitary(U_slow)
    assert maxabs(U_fast - U_slow) < 1e-10


def test_mixing_matrix_4x4_reduces_to_3nu_pmns():
    s12, s23, s13, dCP = 0.3, 0.4, 0.15, 1.2
    U4 = h4.mixing_matrix_4x4(s12, s23, s13, dCP, 0.0, 0.0, 0.0, 0.0, 0.0,
                              compute_matrix_multiplication=False)
    U3 = h3.pmns_mixing_matrix(s12, s23, s13, dCP)
    assert maxabs(U4[:3, :3] - U3) < 1e-12


def test_hamiltonian_4nu_vacuum_reduces_to_3nu():
    """With sterile mixing off, the 4nu vacuum Hamiltonian's active-flavor
    block must equal the 3nu vacuum Hamiltonian exactly."""
    s12, s23, s13, dCP, D21, D31 = 0.3, 0.4, 0.15, 1.2, 7.5e-5, 2.5e-3
    H4 = h4.hamiltonian_4nu_vacuum_energy_independent(
        s12, s23, s13, dCP, 0.0, 0.0, 0.0, 0.0, 0.0, D21, D31, D41=0.0)
    H3 = h3.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)
    assert maxabs(H4[:3, :3] - H3) < 1e-12


# ----------------------------------------------------------------------
# A2: mixing_matrix_5x5
# ----------------------------------------------------------------------

@pytest.mark.parametrize("trial", range(8))
def test_mixing_matrix_5x5_unitary_and_paths_agree(trial):
    s12, s23, s13, s14, s15, s24, s25, s34, s35 = rand_angles(9)
    d13, d14, d15, d24, d35 = rand_phases(5)
    args = (s12, s23, s13, d13, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35)
    U_fast = h5.mixing_matrix_5x5(*args, compute_matrix_multiplication=False)
    U_slow = h5.mixing_matrix_5x5(*args, compute_matrix_multiplication=True)
    assert is_unitary(U_fast)
    assert is_unitary(U_slow)
    assert maxabs(U_fast - U_slow) < 1e-10


def test_mixing_matrix_5x5_reduces_to_3nu_pmns():
    s12, s23, s13, dCP = 0.3, 0.4, 0.15, 1.2
    zeros10 = [0.0]*10
    U5 = h5.mixing_matrix_5x5(s12, s23, s13, dCP, *zeros10,
                              compute_matrix_multiplication=False)
    U3 = h3.pmns_mixing_matrix(s12, s23, s13, dCP)
    assert maxabs(U5[:3, :3] - U3) < 1e-12


def test_hamiltonian_5nu_vacuum_reduces_to_3nu():
    s12, s23, s13, dCP, D21, D31 = 0.3, 0.4, 0.15, 1.2, 7.5e-5, 2.5e-3
    H5 = h5.hamiltonian_5nu_vacuum_energy_independent(
        s12, s23, s13, dCP, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        D21, D31, D41=0.0, D51=0.0)
    H3 = h3.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)
    assert maxabs(H5[:3, :3] - H3) < 1e-12


# ----------------------------------------------------------------------
# A3: hamiltonian_2nu_nsi eps_aa
# ----------------------------------------------------------------------

def test_hamiltonian_2nu_nsi_eps_aa_has_real_effect():
    """Regression test: eps_aa used to sit on both diagonal entries,
    making it a pure multiple of the identity (an unobservable global
    phase, zero effect on any probability). The fix moves it to a single
    entry, so H(eps_aa=5) - H(eps_aa=0) must be nonzero on exactly one
    diagonal entry and exactly zero on the other -- checked by exact
    equality (both sides are exact arithmetic on floats, not accumulated
    numerical noise), not by a tolerance-based comparison whose default
    absolute tolerance could hide the effect at small VCC."""
    VCC = 1e-13
    H0 = h2.hamiltonian_2nu_nsi(VCC, eps_aa=0.0, eps_ab=0.2 + 0.1j)
    H5 = h2.hamiltonian_2nu_nsi(VCC, eps_aa=5.0, eps_ab=0.2 + 0.1j)
    diff = H5 - H0
    assert diff[0, 0] == 5.0 * VCC
    assert diff[1, 1] == 0.0
    assert diff[0, 0] != diff[1, 1]


def test_hamiltonian_2nu_nsi_stays_hermitian():
    H = h2.hamiltonian_2nu_nsi(1e-13, eps_aa=3.0, eps_ab=0.2 + 0.1j)
    assert maxabs(H - H.conj().T) < 1e-15


# ----------------------------------------------------------------------
# B5: hamiltonian_2nu_liv_energy_independent sign convention
# ----------------------------------------------------------------------

def test_hamiltonian_2nu_liv_matches_similarity_transform_convention():
    """H must equal R . diag(b1, b2) . R^T, with R = mixing_matrix_2nu(sxi)
    -- the same convention used by every sibling vacuum/LIV Hamiltonian
    (2nu vacuum's slow path, and 3/4/5nu LIV)."""
    sxi, b1, b2, Lambda, n_liv = 0.3, 1e-9, 2e-9, 1e12, 1
    H = h2.hamiltonian_2nu_liv_energy_independent(sxi, b1, b2, Lambda, n_liv)
    R = h2.mixing_matrix_2nu(sxi)
    H_ref = (1.0 / Lambda)**n_liv * R @ np.diag([b1, b2]) @ R.T
    assert maxabs(H - H_ref) < 1e-25


# ----------------------------------------------------------------------
# C1: the _nsi_td convenience functions
# ----------------------------------------------------------------------

def test_nsi_td_functions_match_their_non_td_counterparts():
    """Regression test: these used to crash with TypeError because they
    forwarded a single `eps` list to a function that takes separate
    named eps_* arguments."""
    l = 2.0
    VCC_func = lambda l: 1e-13 * l

    H2 = h2.hamiltonian_2nu_nsi_td(l, VCC_func, 0.1, 0.2 + 0.1j)
    assert maxabs(H2 - h2.hamiltonian_2nu_nsi(VCC_func(l), 0.1, 0.2 + 0.1j)) == 0.0

    H3 = h3.hamiltonian_3nu_nsi_td(l, VCC_func, 0.1, 0.2j, 0.0, 0.05, 0.0, 0.0)
    assert maxabs(H3 - h3.hamiltonian_3nu_nsi(VCC_func(l), 0.1, 0.2j, 0.0, 0.05,
                                              0.0, 0.0)) == 0.0

    args4 = (0.1, 0.2j, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0)
    H4 = h4.hamiltonian_4nu_nsi_td(l, VCC_func, *args4)
    assert maxabs(H4 - h4.hamiltonian_4nu_nsi(VCC_func(l), *args4)) == 0.0


# ----------------------------------------------------------------------
# General: no name collisions across the wildcard-imported split modules
# ----------------------------------------------------------------------

def test_hamiltonians_package_exposes_all_flavor_counts():
    for name in ['mixing_matrix_2nu', 'hamiltonian_2nu_vacuum_energy_independent',
                 'pmns_mixing_matrix', 'hamiltonian_3nu_vacuum_energy_independent',
                 'mixing_matrix_4x4', 'hamiltonian_4nu_vacuum_energy_independent',
                 'mixing_matrix_5x5', 'hamiltonian_5nu_vacuum_energy_independent']:
        assert hasattr(hams, name), name


# ----------------------------------------------------------------------
# Every exported builder, called once and checked for Hermiticity.
#
# A coverage run found 25 of the builders re-exported in this package's
# __all__ were executed by nothing: not by the tests, and not by the
# library either, since oscprob reaches for the *_energy_independent and
# *_nsi variants instead. They are public, documented API with no caller
# anywhere, which is the state in which a wrong sign or a transposed index
# survives indefinitely.
#
# Hermiticity is the one invariant every Hamiltonian here must satisfy
# regardless of flavor count or scenario, and it is a genuine physical
# constraint rather than an assertion that the code equals itself: it fails
# if a term is added to one triangle and not conjugated into the other,
# which is exactly the shape of the mixing-matrix bugs found in the audit.
# ----------------------------------------------------------------------

# One value per parameter *name*; each builder takes the subset its own
# signature asks for. Angles and phases are deliberately non-trivial, so a
# term dropped from the off-diagonal cannot pass by being zero anyway.
BUILDER_ARGS = {
    'l': 100.0*gd.UNIT_KM, 'energy': 1.0*gd.UNIT_GEV,
    'VCC': 1.0e-13, 'VCC_func': lambda l: 1.0e-13,
    'sth': 0.4, 'Dm2': 2.5e-3,
    's12': 0.55, 's23': 0.68, 's13': 0.15, 'dCP': 3.7, 'd13': 3.7,
    's14': 0.10, 'd14': 1.1, 's24': 0.05, 'd24': 2.2, 's34': 0.02,
    's15': 0.05, 'd15': 0.7, 's25': 0.02, 's35': 0.01, 'd35': 1.9,
    'D21': 7.5e-5, 'D31': 2.5e-3, 'D41': 1.0, 'D51': 2.0,
    'sxi': 0.2, 'sxi12': 0.2, 'sxi23': 0.1, 'sxi13': 0.05,
    'dxiCP': 1.3, 'dxi13': 1.3, 'sxi14': 0.05, 'dxi14': 0.4,
    'sxi24': 0.03, 'dxi24': 0.9, 'sxi34': 0.02,
    'sxi15': 0.03, 'dxi15': 0.6, 'sxi25': 0.02, 'sxi35': 0.01, 'dxi35': 1.5,
    'b1': gd.B1, 'b2': gd.B2, 'b3': gd.B3, 'b4': 3.0e-9, 'b5': 4.0e-9,
    'Lambda': gd.LAMBDA, 'n_liv': 1,
    # NSI: diagonal entries must be real for H to be Hermitian at all, but the
    # off-diagonal ones are deliberately complex -- a term placed in the upper
    # triangle without its conjugate in the lower one is invisible if every
    # epsilon is real, which is precisely the bug this check is here to catch.
    'eps_aa': 0.1, 'eps_ab': 0.05 + 0.02j,
    'eps_ee': 0.2, 'eps_em': 0.1 + 0.03j, 'eps_et': 0.05 - 0.01j,
    'eps_mm': 0.05, 'eps_mt': 0.02 + 0.01j, 'eps_tt': 0.01,
    # The 4nu builders name their single sterile row without an index
    # (eps_es), the 5nu ones number theirs (eps_es1, eps_es2).
    'eps_es': 0.02 + 0.005j, 'eps_ms': 0.01 - 0.003j, 'eps_ts': 0.005 + 0.002j,
    'eps_ss': 0.01,
    'eps_es1': 0.02 + 0.005j, 'eps_ms1': 0.01 - 0.003j, 'eps_ts1': 0.005 + 0.002j,
    'eps_es2': 0.01 - 0.004j, 'eps_ms2': 0.005 + 0.001j, 'eps_ts2': 0.002 - 0.001j,
    'eps_s1s1': 0.01, 'eps_s1s2': 0.005 + 0.001j, 'eps_s2s2': 0.01,
}


def exported_hamiltonian_builders():
    """Every `hamiltonian_*` name in the package's own __all__.

    Read from __all__ rather than hand-listed: __all__ is the package's
    statement about what is public, so anything added there is swept here
    without a second edit, and the test cannot drift from the API it is
    supposed to be guarding.
    """
    return sorted(n for n in hams.__all__ if n.startswith('hamiltonian_'))


@pytest.mark.parametrize("name", exported_hamiltonian_builders())
def test_every_exported_hamiltonian_builder_is_hermitian(name):
    import inspect
    fn = getattr(hams, name)
    params = inspect.signature(fn).parameters
    # Required parameters only; the optional ones (e.g. the
    # compute_matrix_multiplication fast/slow-path flag) are left at their
    # defaults, which is how the library itself calls these.
    required = [p for p, v in params.items() if v.default is inspect.Parameter.empty]
    missing = [p for p in required if p not in BUILDER_ARGS]
    assert not missing, f"{name}: no test value for {missing} -- add it to BUILDER_ARGS"

    H = np.asarray(fn(**{p: BUILDER_ARGS[p] for p in required}))
    dim = int(name[len('hamiltonian_')])
    assert H.shape == (dim, dim), f"{name}: expected {dim}x{dim}, got {H.shape}"
    assert np.all(np.isfinite(H)), f"{name}: non-finite entries"

    # Relative to the size of H itself, with no absolute floor. A floor of 1.0
    # would make this vacuous for exactly the builders most worth checking: the
    # matter and NSI Hamiltonians are proportional to VCC ~ 1e-13 eV, so an
    # absolute 1e-12 tolerance exceeds every entry in the matrix and the
    # assertion can never fail. Verified by deleting a conjugation in
    # hamiltonian_3nu_nsi: with the floor the test still passed, without it the
    # test fails as it should.
    scale = maxabs(H)
    assert scale > 0.0, f"{name}: returned an all-zero matrix"
    assert maxabs(H - H.conj().T) <= 1e-12*scale, f"{name}: not Hermitian"


def test_mixing_matrix_3x3_is_unitary_and_matches_pmns():
    """mixing_matrix_3x3 is exported alongside pmns_mixing_matrix and had no
    caller of its own. If the two ever disagree, one of them is building a
    different convention than the wrappers assume."""
    args = (0.55, 0.68, 0.15, 3.7)
    U = hams.mixing_matrix_3x3(*args)
    assert is_unitary(U)
    assert maxabs(U - hams.pmns_mixing_matrix(*args)) == 0.0
