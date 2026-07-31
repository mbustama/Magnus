# -*- coding: utf-8 -*-
"""Tests of the phase-averaged probabilities (magnus.avgprob).

The averaged probability is the exact L/E -> infinity limit, so the
ground truth is the package's own oscillation probability averaged
numerically over a window: the closed form has to be what that converges
to, in the regime where averaging applies, and must not be quoted in the
regime where it does not.

Both halves are tested here.  The second is not a formality: the first
attempt at the cross-check used a 1000 km baseline and disagreed by 0.4,
because at that distance the solar pair has accumulated only ~0.2 rad of
relative phase and has not decohered at all.  Nothing about the
arithmetic was wrong -- the comparison was being made in a regime where
no averaged expression describes the physics, which is exactly what
coherence_report exists to detect.
"""

import warnings

import numpy as np
import pytest

import magnus.avgprob as ap
import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.oscprob as op

RNG = np.random.default_rng(19)

S12, S23, S13, DCP = 0.55, 0.68, 0.15, 3.7
D21, D31 = 7.5e-5, 2.5e-3
OSC_PARAMS = dict(s12=S12, s23=S23, s13=S13, dCP=DCP, D21=D21, D31=D31)


def maxabs(x):
    return np.max(np.abs(np.asarray(x)))


def pmns():
    return hams.pmns_mixing_matrix(S12, S23, S13, DCP)


def random_unitary(dim, rng=RNG):
    """A Haar-ish random unitary, via QR of a complex Gaussian matrix."""
    X = rng.standard_normal((dim, dim)) + 1j*rng.standard_normal((dim, dim))
    Q, R = np.linalg.qr(X)
    return Q*(np.diag(R)/np.abs(np.diag(R)))


# ----------------------------------------------------------------------
# The invariants every averaged probability matrix must satisfy
# ----------------------------------------------------------------------

@pytest.mark.parametrize("dim", [2, 3, 4, 5])
def test_averaged_probabilities_are_doubly_stochastic_and_symmetric(dim):
    """Rows sum to one because the neutrino has to arrive as something;
    columns too, because the average of a unitary evolution is
    doubly stochastic; and the matrix is symmetric, so the averaged
    probability is the same in both directions."""
    P = ap.averaged_probabilities_from_eigenbasis(random_unitary(dim))

    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-12), "rows do not sum to one"
    assert np.allclose(P.sum(axis=-2), 1.0, atol=1e-12), "columns do not sum to one"
    assert maxabs(P - P.T) < 1e-12, "not symmetric"
    assert np.all(P >= -1e-15)


def test_averaged_probabilities_are_identical_for_antineutrinos():
    """Conjugating the mixing matrix leaves every |V|^2 unchanged, so CP
    violation does not survive the average even though dCP still enters
    through the magnitudes.  A caller expecting an asymmetry here is
    expecting something the physics does not contain."""
    U = pmns()
    P = ap.averaged_probabilities_from_eigenbasis(U)
    P_bar = ap.averaged_probabilities_from_eigenbasis(np.conj(U))

    assert maxabs(P - P_bar) == 0.0


def test_batched_input_matches_one_at_a_time():
    """A leading axis is broadcast over, so a scan over energies costs one
    contraction rather than a Python loop."""
    batch = np.stack([random_unitary(3) for _ in range(4)])
    P_batch = ap.averaged_probabilities_from_eigenbasis(batch)

    assert P_batch.shape == (4, 3, 3)
    for i in range(4):
        assert maxabs(P_batch[i] - ap.averaged_probabilities_from_eigenbasis(batch[i])) < 1e-14


# ----------------------------------------------------------------------
# Coherence blocks: the part the naive formula gets wrong
# ----------------------------------------------------------------------

def test_a_fully_degenerate_spectrum_does_not_oscillate():
    """With every eigenvalue equal there is no relative phase at all, so
    the evolution is the identity and no flavor conversion happens.  The
    naive incoherent sum would instead report a spurious mixture -- this
    is the case that motivates treating coherence blockwise."""
    U = pmns()
    blocks = ap.coherence_blocks([0.0, 0.0, 0.0], phase_scale=1.0e6)

    assert blocks == [[0, 1, 2]]
    P = ap.averaged_probabilities_from_eigenbasis(U, blocks=blocks)
    assert maxabs(P - np.eye(3)) < 1e-12

    naive = ap.averaged_probabilities_from_eigenbasis(U)
    assert maxabs(naive - np.eye(3)) > 0.3, \
        "the naive sum should differ sharply here; if it does not, the test proves nothing"


def test_one_coherent_pair_keeps_its_cross_term():
    """A spectrum with two close eigenvalues and one far away: the close
    pair keeps its interference, the distant one contributes
    incoherently.  Checked against the formula written out by hand for
    that specific block structure."""
    U = pmns()
    blocks = [[0, 1], [2]]

    P = ap.averaged_probabilities_from_eigenbasis(U, blocks=blocks)

    expected = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            coherent = sum(np.conj(U[a, i])*U[b, i] for i in (0, 1))
            expected[a, b] = abs(coherent)**2 + abs(np.conj(U[a, 2])*U[b, 2])**2
    assert maxabs(P - expected) < 1e-14

    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-12), \
        "partial coherence must still conserve probability"


def test_blocks_split_only_once_the_phase_exceeds_a_cycle():
    """The split is decided by accumulated phase, not by the eigenvalue
    spacing on its own: the same spectrum is coherent at a short baseline
    and decohered at a long one."""
    lam = [0.0, 1.0e-12]

    short = 1.0/(1.0e-12)*(0.5*ap.DECOHERENCE_PHASE_THRESHOLD)
    long = 1.0/(1.0e-12)*(2.0*ap.DECOHERENCE_PHASE_THRESHOLD)

    assert ap.coherence_blocks(lam, phase_scale=short) == [[0, 1]]
    assert ap.coherence_blocks(lam, phase_scale=long) == [[0], [1]]


def test_coherence_blocks_of_an_empty_spectrum():
    assert ap.coherence_blocks([], phase_scale=1.0) == []


def test_coherence_report_flags_pairs_in_neither_limit():
    """A pair whose phase is neither much larger than a cycle nor much
    smaller than a radian is in no valid limit, and has to be named rather
    than quietly assigned to one side."""
    lam = np.array([0.0, 1.0])
    scale = 1.0  # phase of exactly 1 rad: between the two thresholds

    blocks, undecided = ap.coherence_report(lam, scale)

    assert blocks == [[0, 1]]
    assert len(undecided) == 1
    i, j, phase = undecided[0]
    assert (i, j) == (0, 1)
    assert phase == pytest.approx(1.0)


def test_coherence_report_is_silent_when_the_limit_is_clean():
    """Well-separated eigenvalues over a long baseline: every pair is
    decohered, nothing to report."""
    blocks, undecided = ap.coherence_report([0.0, 1.0, 2.0], phase_scale=1.0e3)

    assert blocks == [[0], [1], [2]]
    assert undecided == []


# ----------------------------------------------------------------------
# The constant-Hamiltonian entry point
# ----------------------------------------------------------------------

def test_constant_hamiltonian_path_matches_the_mixing_matrix_path():
    """Diagonalizing the vacuum Hamiltonian must reproduce what the PMNS
    matrix gives directly -- the eigenvectors are the mixing matrix, up to
    ordering and phases that |V|^2 does not see."""
    H = hams.hamiltonian_3nu_vacuum_energy_independent(S12, S23, S13, DCP, D21, D31)

    P_from_H = ap.averaged_probabilities_constant_hamiltonian(H)
    P_from_U = ap.averaged_probabilities_from_eigenbasis(pmns())

    assert maxabs(P_from_H - P_from_U) < 1e-12


def test_constant_hamiltonian_accepts_a_baseline_and_uses_it():
    """With a baseline given, the coherence structure is computed from it,
    so a short baseline returns something closer to the identity than the
    fully decohered limit."""
    H = np.diag([0.0, 1.0e-14, 1.0e-12]).astype(complex)
    U = random_unitary(3)
    H = U @ H @ U.conj().T

    P_short = ap.averaged_probabilities_constant_hamiltonian(H, baseline=1.0e11)
    P_long = ap.averaged_probabilities_constant_hamiltonian(H, baseline=1.0e16)

    assert np.allclose(P_short.sum(axis=-1), 1.0, atol=1e-12)
    assert np.allclose(P_long.sum(axis=-1), 1.0, atol=1e-12)
    assert maxabs(P_short - np.eye(3)) < maxabs(P_long - np.eye(3)), \
        "the shorter baseline should retain more coherence, hence be closer to the identity"


def test_constant_hamiltonian_refuses_a_baseline_for_a_batch():
    """The coherence structure can differ from one Hamiltonian to the
    next, so a single set of blocks cannot describe a batch."""
    H = np.stack([np.eye(3, dtype=complex), 2.0*np.eye(3, dtype=complex)])
    with pytest.raises(ValueError, match="batch"):
        ap.averaged_probabilities_constant_hamiltonian(H, baseline=1.0e12)


# ----------------------------------------------------------------------
# Input validation
# ----------------------------------------------------------------------

@pytest.mark.parametrize("bad", [np.zeros((3, 2), dtype=complex), np.zeros(3, dtype=complex)])
def test_eigenvectors_must_be_square(bad):
    with pytest.raises(ValueError, match="square"):
        ap.averaged_probabilities_from_eigenbasis(bad)


def test_hamiltonian_must_be_square():
    with pytest.raises(ValueError, match="square"):
        ap.averaged_probabilities_constant_hamiltonian(np.zeros((3, 2), dtype=complex))


@pytest.mark.parametrize("blocks", [[[0, 1]], [[0, 1, 2], [2]], [[0, 1, 3]]])
def test_blocks_must_partition_the_indices(blocks):
    """A block list that misses an eigenvalue, or counts one twice, would
    silently produce a matrix whose rows do not sum to one."""
    with pytest.raises(ValueError, match="partition"):
        ap.averaged_probabilities_from_eigenbasis(pmns(), blocks=blocks)


# ----------------------------------------------------------------------
# Ground truth: the engine, averaged numerically
# ----------------------------------------------------------------------

def _numerically_averaged_3nu_vacuum(baseline, n_points=2001):
    """Averages the package's own probability over a window sampled
    uniformly in 1/E, in which the oscillation phase is linear."""
    inv_energy = np.linspace(1.0/(10.0*gd.UNIT_GEV), 1.0/(1.0*gd.UNIT_GEV), n_points)
    P = np.asarray(op.osc_prob_3nu_vacuum(1.0/inv_energy, baseline, **OSC_PARAMS))
    return P.mean(axis=0)


def test_closed_form_matches_the_engine_averaged_numerically():
    """The whole claim of this module, checked against the engine it
    replaces: at a baseline long enough for every pair to decohere, the
    closed form is what numerical averaging converges to."""
    P_closed = ap.averaged_probabilities_from_eigenbasis(pmns())
    P_numeric = _numerically_averaged_3nu_vacuum(1.0e8*gd.UNIT_KM)

    assert maxabs(P_numeric - P_closed) < 2.0e-2, \
        "the closed form is not what averaging the engine converges to"


def test_the_averaged_limit_does_not_apply_at_a_terrestrial_baseline():
    """The other half, and the reason coherence_report exists: at 1000 km
    the solar pair has accumulated a fraction of a radian, so the averaged
    expression is simply the wrong description -- it disagrees with the
    engine by tenths, and the diagnostic says so in advance."""
    baseline = 1000.0*gd.UNIT_KM
    energy = 1.0*gd.UNIT_GEV
    eigenvalues = np.array([0.0, D21, D31])/(2.0*energy)

    _, undecided = ap.coherence_report(eigenvalues, baseline)
    assert undecided, "the diagnostic failed to flag a regime where averaging does not apply"

    P_closed = ap.averaged_probabilities_from_eigenbasis(pmns())
    P_numeric = _numerically_averaged_3nu_vacuum(baseline)
    assert maxabs(P_numeric - P_closed) > 0.1, \
        "if these agreed, the diagnostic above would be flagging a regime that is in fact fine"


# ----------------------------------------------------------------------
# The astrophysical case this module exists for
# ----------------------------------------------------------------------

def test_pion_decay_flavor_ratios_arrive_near_equipartition():
    """1:2:0 at the source is the canonical pion-decay composition, and
    the standard result is that it arrives close to 1:1:1 -- the single
    most quoted consequence of averaged astrophysical oscillations."""
    P = ap.averaged_probabilities_from_eigenbasis(pmns())
    at_source = np.array([1.0, 2.0, 0.0])/3.0

    at_earth = at_source @ P

    assert at_earth.sum() == pytest.approx(1.0), "flux is not conserved"
    assert maxabs(at_earth - 1.0/3.0) < 0.05, \
        f"expected near-equipartition, got {at_earth*3.0} in units of 1/3"


# ----------------------------------------------------------------------
# average=True, reaching every routine through the existing kwargs chain
# ----------------------------------------------------------------------

FAR = 1.0e8*gd.UNIT_KM     # long enough that every pair has decohered
ENERGY = 1.0*gd.UNIT_GEV


def test_average_kwarg_reaches_a_wrapper_and_returns_the_closed_form():
    """The keyword is declared on the middle layer only; every wrapper
    inherits it through **kwargs, which is what makes this apply to any
    probability routine rather than to a hand-picked list."""
    P = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, FAR, average=True, **OSC_PARAMS))

    assert maxabs(P - ap.averaged_probabilities_from_eigenbasis(pmns())) < 1e-12


def test_average_defaults_to_off_and_changes_nothing():
    """Omitting the keyword must be indistinguishable from passing False,
    so an existing call cannot change meaning."""
    baseline = 1000.0*gd.UNIT_KM
    without = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, baseline, **OSC_PARAMS))
    explicit = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, baseline, average=False, **OSC_PARAMS))

    assert maxabs(without - explicit) == 0.0


def test_average_warns_where_the_oscillation_has_not_averaged():
    """1000 km: the solar pair has barely turned. The matrix returned is
    still a valid probability matrix, so this warns rather than raising --
    but it must not stay silent."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        P = np.asarray(op.osc_prob_3nu_vacuum(ENERGY, 1000.0*gd.UNIT_KM, average=True,
                                              **OSC_PARAMS))

    assert any(issubclass(w.category, op.PhaseAveragingWarning) for w in caught), \
        "no warning at a baseline where the averaged limit does not apply"
    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-12)


def test_average_is_silent_at_an_astrophysical_baseline():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        op.osc_prob_3nu_vacuum(ENERGY, FAR, average=True, **OSC_PARAMS)

    assert not [w for w in caught if issubclass(w.category, op.PhaseAveragingWarning)]


@pytest.mark.parametrize("num_flavors", [2, 3, 4, 5])
def test_average_works_for_constant_density_matter(num_flavors):
    """Matter of constant density is still a position-independent
    Hamiltonian, so the closed form applies -- with the matter eigenbasis
    rather than the vacuum one."""
    fn = getattr(op, f'osc_prob_{num_flavors}nu_matter_constant_density')
    kwargs = dict(validate_input=False)
    if num_flavors == 2:
        kwargs.update(sth=0.3, Dm2=2.5e-3)

    P = np.asarray(fn(ENERGY, FAR, 3.0*gd.UNIT_G_PER_CM3, average=True, **kwargs))

    assert P.shape == (num_flavors, num_flavors)
    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-12)
    assert maxabs(P - P.T) < 1e-12


def test_average_accepts_an_energy_array_and_a_single_flavor_pair():
    energies = np.array([1.0, 5.0, 10.0])*gd.UNIT_GEV

    P = np.asarray(op.osc_prob_3nu_vacuum(energies, FAR, average=True, **OSC_PARAMS))
    assert P.shape == (3, 3, 3)

    P_e_mu = np.asarray(op.osc_prob_3nu_vacuum(energies, FAR, nu_i=0, nu_f=1, average=True,
                                               **OSC_PARAMS))
    assert P_e_mu.shape == (3,)
    assert np.allclose(P_e_mu, P[:, 0, 1])


def test_vacuum_averaged_probability_does_not_depend_on_energy_or_baseline():
    """Scaling the vacuum Hamiltonian by 1/E leaves its eigenvectors
    alone, so once everything has decohered the answer is one matrix for
    the whole flux -- which is the property that makes this cheap."""
    P_low = np.asarray(op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, FAR, average=True, **OSC_PARAMS))
    P_high = np.asarray(op.osc_prob_3nu_vacuum(1.0e3*gd.UNIT_GEV, 10.0*FAR, average=True,
                                               **OSC_PARAMS))

    assert maxabs(P_low - P_high) < 1e-12


# ----------------------------------------------------------------------
# Position-dependent Hamiltonians: decohere at production, transport along
# the levels, read out at detection
# ----------------------------------------------------------------------

L_SCALE = 0.1*gd.SUN_RADIUS*gd.UNIT_KM
RHO_CENTRAL = 100.0*gd.UNIT_G_PER_CM3
SOLAR_ENERGY = 10.0*gd.UNIT_MEV


def _exponential_H(energy):
    """A solar-like exponential profile, whose central potential sits on the
    1-2 MSW resonance at SOLAR_ENERGY -- i.e. the case this path exists for.

    Built from the same helpers the wrappers use, and with the density
    already in internal units, so no conversion is applied twice.  Getting
    that wrong makes the matter term dominate by seventeen orders of
    magnitude, which looks exactly like a broken formula.
    """
    import magnus.matter as matter
    h_vac = hams.hamiltonian_3nu_vacuum_energy_independent(S12, S23, S13, DCP, D21, D31)
    h_matt = np.diag([1.0, 0.0, 0.0]).astype(complex)
    vcc = matter.vcc_func_from_rho_func(matter.exp_density_profile(RHO_CENTRAL, L_SCALE))

    def H_of_l(l):
        return (1.0/energy)*h_vac + np.asarray(vcc(l))*h_matt

    return H_of_l


def test_adiabatic_averaged_probabilities_are_a_valid_probability_matrix():
    P, report = ap.averaged_probabilities_adiabatic(_exponential_H(SOLAR_ENERGY), 0.0, 5.0*L_SCALE)

    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-10), "rows do not sum to one"
    assert np.all(P >= -1e-15)
    assert report['patches_converged']


def test_adiabatic_averaged_probabilities_need_not_be_symmetric():
    """Unlike the constant-Hamiltonian case, production and detection happen
    in *different* eigenbases, so there is no reason for the matrix to be
    symmetric -- and asserting that it were would be asserting the matter
    effect away."""
    P, _ = ap.averaged_probabilities_adiabatic(_exponential_H(SOLAR_ENERGY), 0.0, 5.0*L_SCALE)

    assert maxabs(P - P.T) > 1e-3


def test_adiabatic_result_matches_the_engine_averaged_numerically():
    """The ground truth, with its own uncertainty accounted for.

    The accumulated phase here is ~1e8 rad, so sampled probabilities are
    effectively independent draws and the reference has a standard error of
    its own.  The comparison is therefore against that error, not against a
    tolerance picked by hand: an earlier 61-point reference disagreed by
    0.12 and was simply undersampled (its own s.e.m. was 0.036).
    """
    energy, l1 = SOLAR_ENERGY, 5.0*L_SCALE
    P_adia, _ = ap.averaged_probabilities_adiabatic(_exponential_H(energy), 0.0, l1)

    n = 241
    inv_energy = np.linspace(1.0/(12.0*gd.UNIT_MEV), 1.0/(8.0*gd.UNIT_MEV), n)
    samples = np.array([
        np.asarray(op.osc_prob_3nu_matter_exp_density(E, l1, 0.0, RHO_CENTRAL, L_SCALE,
                                                      validate_input=False, rtol=1e-4, atol=1e-4,
                                                      **OSC_PARAMS))
        for E in 1.0/inv_energy])
    mean = samples.mean(axis=0)
    sem = samples.std(axis=0)/np.sqrt(n)

    assert maxabs(P_adia - mean) < 3.0*sem.max(), (
        f"closed form differs from the averaged engine by {maxabs(P_adia - mean):.2e}, "
        f"more than three times the reference's own standard error of {sem.max():.2e}")


def test_level_crossing_matrix_is_the_identity_when_evolution_is_adiabatic():
    """No resonance sharp enough to break adiabaticity: the neutrino stays
    on the level it was born on, so the crossing matrix has nothing to do."""
    crossing, windows, converged = ap.level_crossing_matrix(_exponential_H(SOLAR_ENERGY),
                                                            0.0, 5.0*L_SCALE)

    assert windows == []
    assert converged
    assert maxabs(crossing - np.eye(3)) == 0.0


def test_level_crossing_matrix_is_doubly_stochastic():
    """Whatever the crossings do, they redistribute probability among the
    levels without creating or destroying it."""
    crossing, _, _ = ap.level_crossing_matrix(_exponential_H(SOLAR_ENERGY), 0.0, 5.0*L_SCALE)

    assert np.allclose(crossing.sum(axis=-1), 1.0, atol=1e-10)
    assert np.allclose(crossing.sum(axis=-2), 1.0, atol=1e-10)


def test_adiabatic_phase_differences_are_antisymmetric_and_scale_with_length():
    H_of_l = _exponential_H(SOLAR_ENERGY)

    short = ap.adiabatic_phase_differences(H_of_l, 0.0, L_SCALE)
    long = ap.adiabatic_phase_differences(H_of_l, 0.0, 5.0*L_SCALE)

    assert maxabs(short + short.T) < 1e-6*maxabs(short), "not antisymmetric"
    assert maxabs(long) > maxabs(short), "a longer trajectory must accumulate more phase"


def test_average_kwarg_reaches_a_position_dependent_wrapper():
    """The keyword now covers the exponential-density and solar families,
    which previously refused."""
    P = np.asarray(op.osc_prob_3nu_matter_exp_density(
        SOLAR_ENERGY, 5.0*L_SCALE, 0.0, RHO_CENTRAL, L_SCALE, average=True,
        validate_input=False, **OSC_PARAMS))

    assert P.shape == (3, 3)
    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-10)


# ----------------------------------------------------------------------
# Genuine level crossings: the part that is not the identity
# ----------------------------------------------------------------------

CROSSING_L_SCALE = 1.0e5
CROSSING_V0, CROSSING_V1 = 0.3, 0.1


def _two_level_crossing_H(eps):
    """A two-level Hamiltonian with an exponentially falling diagonal that
    crosses a constant one, coupled by ``eps``.

    Note the direction: at a two-level crossing the gap *is* 2*eps, so the
    adiabaticity parameter goes like 1/eps and a **smaller** coupling is
    the more non-adiabatic one -- the opposite of the naive reading, and
    the same construction tests/test_adiabatic.py uses.
    """
    def H_func(l):
        v = CROSSING_V0*np.exp(-l/CROSSING_L_SCALE)
        return np.array([[v, eps], [eps, CROSSING_V1]], dtype=complex)

    return H_func, 0.0, 4.0*CROSSING_L_SCALE


def _landau_zener_hop_probability(eps):
    """exp(-2 pi eps^2 / |d(Delta)/dl|), the analytic transition
    probability for a linear crossing.  The sweep rate at the crossing is
    V1/l_scale, since the exponential passes through V1 there."""
    sweep_rate = CROSSING_V1/CROSSING_L_SCALE
    return np.exp(-2.0*np.pi*eps**2/sweep_rate)


def test_no_window_and_no_hopping_when_the_coupling_is_strong():
    H_func, l0, l1 = _two_level_crossing_H(1.0e-2)
    crossing, windows, converged = ap.level_crossing_matrix(H_func, l0, l1)

    assert windows == [], "a strongly coupled crossing should stay adiabatic"
    assert converged
    assert maxabs(crossing - np.eye(2)) == 0.0


@pytest.mark.parametrize("eps", [1.0e-4, 3.0e-4])
def test_crossing_probabilities_match_landau_zener(eps):
    """The crossing matrix is computed from the package's own
    convergence-checked Magnus patch, with no Landau-Zener assumption
    anywhere in it -- so agreement with the analytic formula is an
    independent check on the machinery, not a tautology."""
    H_func, l0, l1 = _two_level_crossing_H(eps)
    crossing, windows, converged = ap.level_crossing_matrix(H_func, l0, l1)

    assert len(windows) == 1, "expected exactly one non-adiabatic window"
    assert converged

    expected_hop = _landau_zener_hop_probability(eps)
    assert crossing[0, 1] == pytest.approx(expected_hop, rel=0.02), \
        f"hop probability {crossing[0, 1]:.4f} vs Landau-Zener {expected_hop:.4f}"
    assert crossing[0, 0] == pytest.approx(1.0 - expected_hop, rel=0.05)


@pytest.mark.parametrize("eps", [1.0e-4, 3.0e-4, 1.0e-3])
def test_crossing_matrix_stays_doubly_stochastic_through_a_crossing(eps):
    """However much probability moves between levels, none is created or
    destroyed: the patch is unitary, so its modulus-squared in the level
    basis is doubly stochastic."""
    H_func, l0, l1 = _two_level_crossing_H(eps)
    crossing, _, _ = ap.level_crossing_matrix(H_func, l0, l1)

    assert np.allclose(crossing.sum(axis=-1), 1.0, atol=1e-9)
    assert np.allclose(crossing.sum(axis=-2), 1.0, atol=1e-9)
    assert np.all(crossing >= -1e-12)


def test_averaged_probabilities_across_a_crossing_conserve_probability():
    """The full expression with a non-identity crossing matrix still has to
    return a probability matrix."""
    H_func, l0, l1 = _two_level_crossing_H(1.0e-4)
    P, report = ap.averaged_probabilities_adiabatic(H_func, l0, l1)

    assert len(report['windows']) == 1
    assert report['patches_converged']
    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-9)
    assert np.all(P >= -1e-12)


def test_hopping_changes_the_averaged_probability():
    """If the crossing matrix were quietly left as the identity, the
    diabatic and adiabatic cases would give the same answer.  They must
    not."""
    P_adiabatic, _ = ap.averaged_probabilities_adiabatic(*_two_level_crossing_H(1.0e-2))
    P_diabatic, _ = ap.averaged_probabilities_adiabatic(*_two_level_crossing_H(1.0e-4))

    assert maxabs(P_adiabatic - P_diabatic) > 0.3, \
        "hopping across the crossing left the averaged probability unchanged"


def test_a_short_trajectory_reports_pairs_that_have_not_decohered():
    """A stretch over which the levels dephase by about a radian is in
    neither limit: too much to call coherent, far too little to call
    averaged.  The report has to name it rather than answer anyway.

    Note the length is chosen so the phase lands *between* the two
    thresholds.  Shorter is not "more undecided" -- a trajectory short
    enough to accumulate 2e-3 rad is firmly in the coherent limit, which is
    a different statement and not what this checks."""
    H_func, l0, _ = _two_level_crossing_H(1.0e-2)
    l1 = l0 + 5.0                                  # about 1 rad of relative phase

    dphi = abs(ap.adiabatic_phase_differences(H_func, l0, l1)[0, 1])
    assert ap.COHERENCE_PHASE_THRESHOLD <= dphi <= ap.DECOHERENCE_PHASE_THRESHOLD, \
        f"fixture drifted out of the undecided band: {dphi:.3f} rad"

    _, report = ap.averaged_probabilities_adiabatic(H_func, l0, l1)

    assert report['undecided'], "a pair in neither limit was not reported"


def _two_window_H():
    """Two independent 2-level sectors with well-separated crossings, so the
    trajectory passes through two distinct non-adiabatic windows.  Same
    construction tests/test_adiabatic.py uses to check that far-apart
    resonances are kept separate."""
    l_scale = 1.0e5

    def H_func(l):
        v = np.exp(-l/l_scale)
        eps_a, eps_b = 3.0e-4, 1.0e-4
        v_res_a, v_res_b = 0.535, 0.0535
        base_b = 1.3
        return np.array([
            [v, eps_a, 0.0, 0.0],
            [eps_a, v_res_a, 0.0, 0.0],
            [0.0, 0.0, base_b + v, eps_b],
            [0.0, 0.0, eps_b, base_b + v_res_b],
        ], dtype=complex)

    return H_func, 0.0, 4.0*l_scale


def test_two_crossings_compose_into_one_doubly_stochastic_matrix():
    """Two windows in sequence: the level-transfer matrices multiply, and
    the product is still doubly stochastic."""
    H_func, l0, l1 = _two_window_H()
    crossing, windows, converged = ap.level_crossing_matrix(H_func, l0, l1)

    assert len(windows) == 2, f"expected two windows, got {windows}"
    assert windows[0][1] < windows[1][0], "windows should be disjoint and ordered"
    assert converged
    assert np.allclose(crossing.sum(axis=-1), 1.0, atol=1e-9)
    assert np.allclose(crossing.sum(axis=-2), 1.0, atol=1e-9)
    assert maxabs(crossing - np.eye(4)) > 0.1, "neither crossing did anything"


def test_between_crossing_coherence_is_checked_not_assumed():
    """Composing two crossings as *probabilities* rather than as amplitudes
    throws away any interference surviving between them, which is only
    legitimate if the levels dephase in the stretch separating the two
    windows.  That stretch is checked, and the check has to be reachable --
    it is the assumption the composition rests on."""
    H_func, l0, l1 = _two_window_H()
    _, report = ap.averaged_probabilities_adiabatic(H_func, l0, l1)

    assert len(report['windows']) == 2
    assert 'undecided_between_crossings' in report

    # Here the two windows are far apart and the levels dephase enormously
    # between them, so the composition is sound and nothing is reported.
    assert report['undecided_between_crossings'] == [], \
        "these windows are far apart; the levels between them have dephased"


def test_crossings_too_close_in_phase_are_reported():
    """The same two crossings, with the whole Hamiltonian scaled down so the
    levels barely dephase between them.

    Composing crossings as probabilities discards the interference that
    survives from one to the next, which is only legitimate once they have
    dephased.  Here they have not, and the report must say so rather than
    silently multiply the two matrices together anyway.  (Scaling H down
    leaves the crossings in place -- the adiabaticity parameter goes like
    1/scale, so they get sharper -- while shrinking the accumulated phase
    in proportion.)"""
    l_scale = 1.0e5
    scale = 1.0e-4

    def H_func(l):
        v = np.exp(-l/l_scale)
        return scale*np.array([
            [v, 3.0e-4, 0.0, 0.0],
            [3.0e-4, 0.535, 0.0, 0.0],
            [0.0, 0.0, 1.3 + v, 1.0e-4],
            [0.0, 0.0, 1.0e-4, 1.3 + 0.0535],
        ], dtype=complex)

    _, report = ap.averaged_probabilities_adiabatic(H_func, 0.0, 4.0*l_scale)

    assert len(report['windows']) == 2, "fixture no longer produces two separate windows"
    assert report['undecided_between_crossings'], \
        "levels that have not dephased between two crossings were composed without comment"

    # Every entry names the stretch and the pair it refers to, so the caller
    # can see which crossings and which levels are at issue.
    l_start, l_end, i, j, phase = report['undecided_between_crossings'][0]
    assert l_start < l_end
    assert 0 <= i < j < 4
    assert phase <= ap.DECOHERENCE_PHASE_THRESHOLD


# ----------------------------------------------------------------------
# The sampled fallback: no closed form, and honest about it
# ----------------------------------------------------------------------

def test_numerical_averaging_of_a_constant_is_exact():
    """A probability that does not vary across the window averages to
    itself, with no scatter -- the degenerate case that catches a sampling
    or weighting mistake."""
    target = np.array([[0.7, 0.3], [0.3, 0.7]])
    mean, sem = ap.averaged_probabilities_numerically(lambda e: target, 1.0*gd.UNIT_GEV)

    assert maxabs(mean - target) < 1e-15
    assert sem == pytest.approx(0.0, abs=1e-15)


def test_numerical_averaging_reports_its_own_uncertainty():
    """The sampled average has a standard error the closed forms do not,
    and it is returned rather than left for the caller to guess."""
    rng = np.random.default_rng(3)

    def noisy(energy):
        p = 0.5 + 0.4*np.sin(energy/1.0e3) + 0.01*rng.standard_normal()
        p = min(max(p, 0.0), 1.0)
        return np.array([[p, 1.0 - p], [1.0 - p, p]])

    mean, sem = ap.averaged_probabilities_numerically(noisy, 1.0*gd.UNIT_GEV)

    assert sem > 0.0
    assert np.allclose(mean.sum(axis=-1), 1.0, atol=1e-12)


def test_numerical_averaging_depends_on_the_window():
    """The point that justifies warning about the default: this is the
    average over a particular window, not a limit, so a different window is
    a different number."""
    def oscillating(energy):
        p = 0.5*(1.0 + np.cos(energy/1.0e5))
        return np.array([[p, 1.0 - p], [1.0 - p, p]])

    narrow, _ = ap.averaged_probabilities_numerically(oscillating, 1.0e6, relative_spread=0.01)
    wide, _ = ap.averaged_probabilities_numerically(oscillating, 1.0e6, relative_spread=0.5)

    assert maxabs(narrow - wide) > 1e-3, \
        "if the window did not matter, defaulting it would need no warning"


@pytest.mark.parametrize("spread", [0.0, -0.1, 1.0, 2.0])
def test_numerical_averaging_rejects_an_impossible_window(spread):
    with pytest.raises(ValueError, match="relative_spread"):
        ap.averaged_probabilities_numerically(lambda e: np.eye(2), 1.0, relative_spread=spread)


@pytest.mark.parametrize("n_samples", [1, 0, -5])
def test_numerical_averaging_needs_at_least_two_samples(n_samples):
    with pytest.raises(ValueError, match="n_samples"):
        ap.averaged_probabilities_numerically(lambda e: np.eye(2), 1.0, n_samples=n_samples)


def test_earth_averaging_falls_back_to_sampling_and_says_so():
    """PREM steps through layer boundaries, so there is no instantaneous
    eigenbasis and no closed form.  A number is still returned -- averaged
    over a default window -- and the warning has to name the window, say
    this is not the L/E limit, and quote the resulting uncertainty, so the
    figure is never silently dependent on a constant nobody chose."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        P = np.asarray(op.osc_prob_3nu_earth(1.0*gd.UNIT_GEV, costhz=-0.8,
                                             L=2*6371.0*0.8*gd.UNIT_KM, average=True,
                                             validate_input=False, **OSC_PARAMS))

    assert np.allclose(P.sum(axis=-1), 1.0, atol=1e-9)
    assert np.all(P >= -1e-12)

    messages = [str(w.message) for w in caught if issubclass(w.category, op.PhaseAveragingWarning)]
    assert messages, "sampled a window without saying so"
    text = messages[0]
    assert str(100.0*ap.AVG_DEFAULT_ENERGY_SPREAD) in text, "the window width is not named"
    assert str(ap.AVG_DEFAULT_N_SAMPLES) in text, "the sample count is not named"
    assert 'standard error' in text, "the uncertainty of the result is not quoted"
