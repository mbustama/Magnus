# -*- coding: utf-8 -*-
"""Tests of the pseudo-Dirac Hamiltonians (magnus.hamiltonians).

The physical content of the pseudo-Dirac case is a separation of scales: the
pair splittings are far below the standard ones, so over an astrophysical
baseline the standard phases average away while each pair stays mutually
coherent.  That is the regime the coherent-block averaging form exists for, and
the regime in which the naive sum over eigenstates is wrong by the number of
states sharing a block.

The central check here is an *invariant* rather than a stored constant: in the
fully coherent limit the block form has to collapse back onto the ordinary Dirac
answer, because each pair's two columns carry the parent state's mixing split
evenly between them.  Written that way the test does not depend on which global
fit is loaded.  The stored numbers below are recorded for NuFIT 6.1 NO and are
checked against the invariant, not trusted ahead of it.
"""

import itertools
import warnings

import numpy as np
import pytest

import magnus.avgprob as ap
import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.oscprob as op


MPC_IN_KM = 3.0856775814913673e19

# The configuration the handover brief specifies: 100 TeV, 100 Mpc, NuFIT 6.1
# normal ordering, every mass state paired.
ENERGY = 100.0*gd.UNIT_TEV
BASELINE = 100.0*MPC_IN_KM*gd.UNIT_KM
COHERENT_SPLITTING = 1.0e-18


def _osc():
    return gd.load_nufit_params('NuFIT 6.1', 'NO')


def _pmns():
    o = _osc()
    return hams.pmns_mixing_matrix(o['s12'], o['s23'], o['s13'], o['dCP'])


def _mass_squared():
    o = _osc()
    return [0.0, o['D21'], o['D31']]


def _block_form(W, blocks, a, b):
    """Eq. (17): the coherent-block averaged probability."""
    return float(sum(abs(sum(np.conj(W[a, i])*W[b, i] for i in block))**2
                     for block in blocks).real)


def _naive_sum(W, a, b):
    """Eq. (16): one term per eigenstate, valid only when all pairs decohere."""
    return float(sum(abs(np.conj(W[a, i])*W[b, i])**2
                     for i in range(W.shape[0])).real)


# ----------------------------------------------------------------------
# The Dirac limit
# ----------------------------------------------------------------------

def test_empty_pairing_reproduces_the_dirac_hamiltonian():
    """No pairs must reduce to the ordinary three-flavor Hamiltonian.

    To machine precision, not merely to a tolerance: with no pairs the
    construction is the same product of the same matrices, so any disagreement
    beyond round-off would mean the extended builder is not doing what the
    three-flavor one does.
    """
    o = _osc()
    H_pd = hams.hamiltonian_pseudo_dirac_vacuum_energy_independent(
        _pmns(), _mass_squared(), {})
    H_3nu = hams.hamiltonian_3nu_vacuum_energy_independent(
        o['s12'], o['s23'], o['s13'], o['dCP'], o['D21'], o['D31'])
    assert H_pd.shape == (3, 3)
    assert np.max(np.abs(H_pd - H_3nu)) < 1.0e-16*max(1.0, np.max(np.abs(H_3nu)))


def test_empty_pairing_leaves_the_mixing_matrix_untouched():
    U = _pmns()
    assert np.array_equal(hams.pseudo_dirac_mixing_matrix(U, {}), U)


def test_empty_pairing_leaves_the_masses_untouched():
    m2 = _mass_squared()
    assert np.allclose(hams.pseudo_dirac_mass_squared(m2, {}), m2, atol=0.0)


# ----------------------------------------------------------------------
# Structure, for every pairing pattern
# ----------------------------------------------------------------------

@pytest.mark.parametrize('paired', [c for r in range(4)
                                    for c in itertools.combinations(range(3), r)])
def test_dimension_and_unitarity_for_every_pairing_pattern(paired):
    """n = n_active + len(pairs), and the extension stays unitary.

    Partial patterns are the point: the interface is required to pair
    individual mass states, not to toggle all of them at once.
    """
    pairs = {j: COHERENT_SPLITTING*(j + 1) for j in paired}
    W = hams.pseudo_dirac_mixing_matrix(_pmns(), pairs)
    M2 = hams.pseudo_dirac_mass_squared(_mass_squared(), pairs)
    n = 3 + len(pairs)
    assert W.shape == (n, n)
    assert M2.shape == (n,)
    assert np.allclose(W @ np.conj(W.T), np.eye(n), atol=1.0e-14)


@pytest.mark.parametrize('paired', [(0,), (1,), (2,), (0, 2), (0, 1, 2)])
def test_paired_states_split_the_masses_and_unpaired_ones_do_not(paired):
    pairs = {j: 1.0e-18*(j + 1) for j in paired}
    m2 = _mass_squared()
    out = list(hams.pseudo_dirac_mass_squared(m2, pairs))
    expected = []
    for j in range(3):
        expected.append(m2[j])
        if j in pairs:
            expected.append(m2[j] + pairs[j])
    assert np.allclose(out, expected, atol=0.0)


def test_a_partially_paired_spectrum_is_the_mixed_case():
    """Two of three paired: two blocks of two, one singleton."""
    pairs = {0: COHERENT_SPLITTING, 2: COHERENT_SPLITTING}
    W = hams.pseudo_dirac_mixing_matrix(_pmns(), pairs)
    M2 = hams.pseudo_dirac_mass_squared(_mass_squared(), pairs)
    assert W.shape == (5, 5)
    blocks = ap.coherence_blocks(M2/(2.0*ENERGY), BASELINE)
    assert [len(b) for b in blocks] == [2, 1, 2]


# ----------------------------------------------------------------------
# The physics: the block form against the naive sum
# ----------------------------------------------------------------------

def test_coherent_pairs_are_grouped_and_wide_ones_are_not():
    """coherence_blocks recovers the pairing, and loses it when it should."""
    pairs = {j: COHERENT_SPLITTING for j in range(3)}
    M2 = hams.pseudo_dirac_mass_squared(_mass_squared(), pairs)
    assert ap.coherence_blocks(M2/(2.0*ENERGY), BASELINE) == [[0, 1], [2, 3], [4, 5]]

    # A splitting large enough to decohere the pair leaves six singletons.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', hams.PseudoDiracSplittingWarning)
        M2_wide = hams.pseudo_dirac_mass_squared(
            _mass_squared(), {j: 1.0e-10 for j in range(3)})
    assert ap.coherence_blocks(M2_wide/(2.0*ENERGY), BASELINE) == [[i] for i in range(6)]


def test_the_block_form_collapses_onto_the_dirac_answer():
    """The invariant, and the reason the naive sum is wrong.

    Each pair splits its parent state's mixing evenly between two columns, so
    summing the amplitudes *within* a block rebuilds |U_aj|^2 exactly and the
    block form returns the Dirac probability.  Summing the probabilities
    instead -- one term per eigenstate -- loses a factor equal to the number of
    states sharing the block, which is two.
    """
    U = _pmns()
    pairs = {j: COHERENT_SPLITTING for j in range(3)}
    W = hams.pseudo_dirac_mixing_matrix(U, pairs)
    M2 = hams.pseudo_dirac_mass_squared(_mass_squared(), pairs)
    blocks = ap.coherence_blocks(M2/(2.0*ENERGY), BASELINE)

    dirac = float(np.sum(np.abs(U[0, :])**4))
    block = _block_form(W, blocks, 0, 0)
    naive = _naive_sum(W, 0, 0)

    assert block == pytest.approx(dirac, abs=1.0e-12)
    assert naive == pytest.approx(0.5*dirac, abs=1.0e-12)
    # Recorded for NuFIT 6.1 NO.  The assertions above are the ones that carry
    # the physics; these pin the parameter set the numbers came from.
    assert block == pytest.approx(0.54814, abs=5.0e-5)
    assert naive == pytest.approx(0.27407, abs=5.0e-5)


def test_every_active_channel_collapses_onto_its_dirac_value():
    """Not just the survival channel: the identity holds for all of them."""
    U = _pmns()
    pairs = {j: COHERENT_SPLITTING for j in range(3)}
    W = hams.pseudo_dirac_mixing_matrix(U, pairs)
    M2 = hams.pseudo_dirac_mass_squared(_mass_squared(), pairs)
    blocks = ap.coherence_blocks(M2/(2.0*ENERGY), BASELINE)
    for a in range(3):
        for b in range(3):
            dirac = float(np.sum(np.abs(U[a, :])**2 * np.abs(U[b, :])**2))
            assert _block_form(W, blocks, a, b) == pytest.approx(dirac, abs=1.0e-12)


def test_the_six_state_hamiltonian_propagates_and_stays_unitary():
    """osc_prob accepts the 6x6 through the generic callable interface."""
    pairs = {j: COHERENT_SPLITTING for j in range(3)}
    H = hams.hamiltonian_pseudo_dirac_vacuum_energy_independent(
        _pmns(), _mass_squared(), pairs)/ENERGY
    P = np.asarray(op.osc_prob(H, 0.0, 1.0e-3*MPC_IN_KM*gd.UNIT_KM))
    assert P.shape == (6, 6)
    assert np.allclose(P.sum(axis=-1), 1.0, atol=1.0e-12)
    assert (P >= -1.0e-15).all()


# ----------------------------------------------------------------------
# Matter
# ----------------------------------------------------------------------

def test_the_matter_projector_gives_the_steriles_their_neutral_current_entry():
    """diag(1, 0, 0, r/2, ...): the sterile partners feel neither current."""
    H = hams.hamiltonian_pseudo_dirac_matter(
        1.0, 3, {0: 1.0e-18, 1: 1.0e-18, 2: 1.0e-18},
        ratio_number_neutrons_to_protons=1.2)
    assert np.allclose(np.diag(H).real, [1.0, 0.0, 0.0, 0.6, 0.6, 0.6])
    assert np.allclose(H - np.diag(np.diag(H)), 0.0)


def test_the_matter_projector_follows_a_partial_pairing():
    H = hams.hamiltonian_pseudo_dirac_matter(
        1.0, 3, {0: 1.0e-18, 2: 1.0e-18}, ratio_number_neutrons_to_protons=1.2)
    assert np.allclose(np.diag(H).real, [1.0, 0.0, 0.0, 0.6, 0.6])


def test_the_matter_term_broadcasts_over_positions():
    """An array of potentials returns a stack, position axis leading."""
    H = hams.hamiltonian_pseudo_dirac_matter(
        np.array([1.0, 2.0]), 3, {0: 1.0e-18})
    assert H.shape == (2, 4, 4)
    assert H[1][0][0] == pytest.approx(2.0)


def test_the_matter_term_refuses_a_flavor_count_the_projector_cannot_serve():
    """Refusing beats a silently mislabelled sterile state.

    magnus.matter.matter_potential_projector treats everything past the third
    state as sterile.  At two active flavors that convention would put a
    charged-current zero where a partner's r/2 belongs, and the answer would be
    wrong without saying so.
    """
    with pytest.raises(ValueError, match='three active flavors'):
        hams.hamiltonian_pseudo_dirac_matter(1.0, 2, {0: 1.0e-18})


# ----------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------

def test_a_pairing_index_outside_the_spectrum_raises():
    with pytest.raises(ValueError, match='outside the range'):
        hams.pseudo_dirac_mixing_matrix(_pmns(), {3: 1.0e-18})


def test_a_negative_or_zero_splitting_raises():
    for bad in (-1.0e-18, 0.0):
        with pytest.raises(ValueError, match='must be positive'):
            hams.pseudo_dirac_mass_squared(_mass_squared(), {0: bad})


def test_a_non_finite_splitting_raises():
    with pytest.raises(ValueError, match='not finite'):
        hams.pseudo_dirac_mass_squared(_mass_squared(), {0: np.inf})


def test_a_non_integer_pairing_key_raises():
    with pytest.raises(TypeError, match='integer mass-state indices'):
        hams.pseudo_dirac_mixing_matrix(_pmns(), {'0': 1.0e-18})


def test_a_non_mapping_pairing_raises():
    with pytest.raises(TypeError, match='must be a mapping'):
        hams.pseudo_dirac_mixing_matrix(_pmns(), [0, 2])


def test_a_non_square_mixing_matrix_raises():
    with pytest.raises(ValueError, match='must be square'):
        hams.pseudo_dirac_mixing_matrix(np.zeros((3, 2)), {})


def test_a_splitting_comparable_with_the_standard_ones_warns():
    """The separation of scales is the physics; losing it is worth saying."""
    with pytest.warns(hams.PseudoDiracSplittingWarning, match='not small'):
        hams.pseudo_dirac_mass_squared(_mass_squared(), {0: 1.0e-4})


def test_a_genuinely_small_splitting_does_not_warn():
    with warnings.catch_warnings():
        warnings.simplefilter('error', hams.PseudoDiracSplittingWarning)
        hams.pseudo_dirac_mass_squared(_mass_squared(), {0: 1.0e-18})
