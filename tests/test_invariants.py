# -*- coding: utf-8 -*-
"""Oracle-free invariants: the answer must not depend on which door you came in.

Every defect in ``docs/dev/FINDINGS_ADVERSARIAL_VALIDATION.md`` lived at a dispatch seam, and
these are the properties a seam cannot break without being noticed.  None of them needs a
reference solution, which is what makes them cheap enough to run here rather than only in a
battery.

**Two of them are exact and the rest are not, and the difference is the point.**  Reordering
baselines or distributing points over workers must not change a single bit -- there is no
numerical reason for it to.  ``strategy='auto'`` and ``strategy='magnus'`` *are* different
methods and will differ; the work is choosing how much, per invariant, from a measured
distribution rather than from a guess.  The bounds below come from
``docs/dev/adversarial_batteries/invariants.py``, which sweeps 60 configurations (5 profile
families x d in {2,3} x 2 energies x N in {1, 8, 30}); each test's docstring records what that
sweep measured.
"""

import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.oscprob as op

L_SCALE = gd.L_SCALE_SUN
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL
L1 = 1.0*L_SCALE
PARAMS = {'sth': 0.55, 'Dm2': 7.5e-5}
PARAMS_3NU = {'s12': 0.55, 's23': 0.68, 's13': 0.15, 'dCP': 3.7,
              'D21': 7.5e-5, 'D31': 2.5e-3}


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x))))


def solar_ne():
    return matter.exp_density_profile(NE0, L_SCALE)


def multi_resonance_ne():
    def ne(l):
        x = np.asarray(l, dtype=float)
        out = NE0*np.exp(-x/L_SCALE)*(1.0 + 0.9*np.sin(2.0*np.pi*6.0*x/L1))
        return out[()] if out.ndim == 0 else out
    return ne


PROFILES = {'solar': solar_ne, 'multi-resonance': multi_resonance_ne}


def call(profile, d, energy, L, **kw):
    params = PARAMS if d == 2 else PARAMS_3NU
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = op.osc_prob_matter_std_potential(
            d, PROFILES[profile](), energy, L, params, L0=0.0,
            density_is_of_number_of_electrons=True, **kw)
    return np.asarray(P).reshape(np.size(L), d, d)


MATRIX = [('solar', 2, 10.0e6), ('solar', 3, 50.0e6),
          ('multi-resonance', 2, 50.0e6), ('multi-resonance', 3, 50.0e6)]


# ----------------------------------------------------------------------
# exact invariants -- no tolerance, because there is no reason for one
# ----------------------------------------------------------------------

@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_shuffled_baselines_give_identical_results(profile, d, energy):
    """Measured over 60 configurations: 0.00e+00 every time, so this is asserted exactly."""
    Ls = np.linspace(0.05*L1, L1, 12)
    perm = np.random.default_rng(4).permutation(len(Ls))
    P = call(profile, d, energy, Ls)
    P_shuffled = call(profile, d, energy, Ls[perm])
    unshuffled = np.empty_like(P_shuffled)
    unshuffled[perm] = P_shuffled
    assert maxabs(P - unshuffled) == 0.0


@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_parallel_points_give_identical_results(profile, d, energy):
    """Measured over 60 configurations: 0.00e+00 every time.  The warm start carries state
    from point to point, so this is the invariant that would break first if it were carried
    wrongly across workers."""
    Ls = np.linspace(0.05*L1, L1, 12)
    assert maxabs(call(profile, d, energy, Ls) - call(profile, d, energy, Ls, n_jobs=2)) == 0.0


@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_repeated_calls_are_bit_identical(profile, d, energy):
    Ls = np.linspace(0.05*L1, L1, 8)
    assert maxabs(call(profile, d, energy, Ls) - call(profile, d, energy, Ls)) == 0.0


# ----------------------------------------------------------------------
# near-exact: unitarity and composition
# ----------------------------------------------------------------------

@pytest.mark.parametrize('profile,d,energy', MATRIX)
@pytest.mark.parametrize('N', [1, 30])
def test_probability_rows_and_columns_sum_to_one(profile, d, energy, N):
    """Unitarity read off the output the caller receives, not off an internal operator.

    Measured worst case over the sweep: 2.7e-11, and over d = 2...5 and N = 25...1e5 in the
    adversarial validation: 1.6e-11.  Four decades of N cost half a decade of unitarity, so
    1e-9 is a bound with room rather than a target being scraped."""
    Ls = np.linspace(0.05*L1, L1, N) if N > 1 else np.array([L1])
    P = call(profile, d, energy, Ls)
    assert maxabs(P.sum(axis=1) - 1.0) < 1e-9
    assert maxabs(P.sum(axis=2) - 1.0) < 1e-9


@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_hybrid_propagator_composes(profile, d, energy):
    """:math:`U(0 \\to L) = U(L/2 \\to L)\\,U(0 \\to L/2)` is exact for the true evolution, so
    any disagreement is the two sides' own approximation error, each budgeted at the default
    1e-3.  Measured worst case over the sweep: 1.3e-03; the bound is 5e-3."""
    import magnus.adiabatic as ad

    params = PARAMS if d == 2 else PARAMS_3NU
    h_vac = _h_vac(d, params)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0
    vcc = matter.vcc_func_from_rho_func(
        PROFILES[profile](), 0.0, 1.0, 0.5, nubar=False,
        density_matter_is_in_g_per_cm3=False, density_is_of_number_of_electrons=True)

    def H_func(l):
        v = np.asarray(vcc(l))
        return (1.0/energy)*h_vac + v[..., None, None]*proj

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        U_full, _, _ = ad.hybrid_propagator(H_func, 0.0, L1)
        U_a, _, _ = ad.hybrid_propagator(H_func, 0.0, 0.5*L1)
        U_b, _, _ = ad.hybrid_propagator(H_func, 0.5*L1, L1)
    assert maxabs(U_full - U_b @ U_a) < 5e-3


def _h_vac(d, params):
    import magnus.hamiltonians as hams
    if d == 2:
        return np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(
            params['sth'], params['Dm2']), dtype=complex)
    return np.asarray(hams.hamiltonian_3nu_vacuum_energy_independent(
        params['s12'], params['s23'], params['s13'], params['dCP'],
        params['D21'], params['D31']), dtype=complex)


# ----------------------------------------------------------------------
# cross-entry-point agreement -- these are allowed to differ, within a measured bound
# ----------------------------------------------------------------------

@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_strategies_agree_within_a_measured_bound(profile, d, energy):
    """'auto' and 'magnus' are different methods and 'auto' is usually the better one, so the
    assertion is "within tolerance", never "equal".

    Measured over the sweep: median 2.0e-04, p90 1.4e-03, worst 6.5e-03 -- the worst on a noisy
    profile at a single point, where both answers are inside their own requested 1e-3 of the
    truth and differ because they resolve it differently.  The bound is 2e-2: loose enough that
    an ordinary method difference does not fail, tight enough that the 4.6e-01 an unmarked
    discontinuity produces would."""
    for N in (1, 12):
        Ls = np.linspace(0.05*L1, L1, N) if N > 1 else np.array([L1])
        P_auto = call(profile, d, energy, Ls)
        P_magnus = call(profile, d, energy, Ls, strategy='magnus')
        assert maxabs(P_auto - P_magnus) < 2e-2, 'N=%d' % N


@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_cumulative_settings_agree_within_a_measured_bound(profile, d, energy):
    """cumulative True / False / 'auto' build different grids, so they differ; measured over
    the sweep at median 1.2e-04, p90 1.2e-03, worst 3.1e-03.  'auto' must equal one of the two
    exactly, since it resolves to one of them rather than being a third thing."""
    Ls = np.linspace(0.05*L1, L1, 12)
    P_auto = call(profile, d, energy, Ls, cumulative='auto')
    P_true = call(profile, d, energy, Ls, cumulative=True)
    P_false = call(profile, d, energy, Ls, cumulative=False)
    assert maxabs(P_true - P_false) < 2e-2
    assert min(maxabs(P_auto - P_true), maxabs(P_auto - P_false)) == 0.0


@pytest.mark.parametrize('profile,d,energy', MATRIX)
def test_a_scan_agrees_with_the_same_points_computed_singly(profile, d, energy):
    """A 12-point scan takes the cumulative engine; the same baselines one at a time take the
    hybrid one -- two different families answering the same question.

    Measured over the sweep: median 1.9e-09, p90 2.6e-04, worst 6.2e-04, i.e. inside the
    requested 1e-3 everywhere.  It was not always: the same sweep run before the resolution
    test was repaired reported a worst case of **1.1e-02**, which was this invariant doing its
    job -- the hybrid path was declining a smooth Gaussian bump as though it were
    discontinuous.  Bound 2e-2, for the same reason as the strategy invariant above."""
    Ls = np.linspace(0.05*L1, L1, 12)
    scan = call(profile, d, energy, Ls)
    singly = np.concatenate([call(profile, d, energy, np.array([L])) for L in Ls])
    assert maxabs(scan - singly) < 2e-2
