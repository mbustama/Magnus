# -*- coding: utf-8 -*-
"""Tests of the engine-reporting instrument: ``strategy_info`` and ``cross_check_strategies``.

Both are built on the same trace (``magnus.oscprob._engine_probe``), so they are tested
together.  The point of the cross-check is stated in its own docstring and reproduced by
``docs/dev/adversarial_batteries/crosscheck_acceptance.py``: every silently-wrong result the
adversarial validation found came from a method comparing itself with itself, and two genuinely
different engines disagreeing needs no oracle at all.
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
PARAMS_2NU = {'sth': 0.55, 'Dm2': 7.5e-5}
PARAMS_3NU = {'s12': 0.55, 's23': 0.68, 's13': 0.15, 'dCP': 3.7,
              'D21': 7.5e-5, 'D31': 2.5e-3}


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x))))


def solar_ne():
    return matter.exp_density_profile(NE0, L_SCALE)


def multi_resonance_ne(n_cycles=6.0):
    def ne(l):
        x = np.asarray(l, dtype=float)
        base = NE0*np.exp(-x/L_SCALE)
        out = base*(1.0 + 0.9*np.sin(2.0*np.pi*n_cycles*x/L1))
        return out[()] if out.ndim == 0 else out
    return ne


def step_ne():
    """The FINDINGS §3.1 construction: a density jump the caller does not declare."""
    def ne(l):
        x = np.asarray(l, dtype=float)
        out = np.where(x < 0.5*L1, 0.02*NE0, 0.30*NE0)
        return out[()] if out.ndim == 0 else out
    return ne


def call(ne, energy, L, params=PARAMS_2NU, d=2, **kw):
    return op.osc_prob_matter_std_potential(
        d, ne, energy, L, params, L0=0.0, density_is_of_number_of_electrons=True, **kw)


# ----------------------------------------------------------------------
# strategy_info
# ----------------------------------------------------------------------

def test_strategy_info_names_the_engine_that_answered():
    """The three engines a solar request can reach, each reported by name."""
    seen = {}
    for tag, kw in (('hybrid', dict(strategy='hybrid')),
                    ('magnus-strategy', dict(strategy='magnus')),
                    ('scan', dict())):
        info = {}
        L = np.linspace(0.05*L1, L1, 30) if tag == 'scan' else L1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            call(solar_ne(), 10.0e6, L, strategy_info=info, **kw)
        seen[tag] = info['engine']
        assert info['family'] == op.ENGINE_FAMILIES[info['engine']]
    assert seen['hybrid'] == 'hybrid'
    assert seen['scan'] == 'cumulative'
    # strategy='magnus' must reach a magnus-family or interaction-picture engine, never hybrid
    assert seen['magnus-strategy'] != 'hybrid'


def test_strategy_info_reports_a_silent_hybrid_fallback():
    """Under strategy='auto' an uncertified hybrid result falls back silently -- which is the
    right default and leaves a user debugging a moved result with nothing to look at.  This is
    the way to see it, and the unmarked step is the case that produces it.

    The reason must be the *specific* one, not merely "did not certify".  ``certified=False``
    on its own does not say which of two things went wrong, and their cures are opposite: an
    unresolved profile wants ``t_breakpoints``, an exhausted refinement wants a looser
    tolerance."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(step_ne(), 50.0e6, L1, strategy_info=info)
    assert info['engine'] != 'hybrid'
    assert [eng for eng, _ in info['declined']] == ['hybrid']
    assert 'not resolved' in info['declined'][0][1]

    # ... and the two reasons are genuinely distinguishable at the level where they are
    # decided.  Starved of refinement iterations, a perfectly resolved profile fails to
    # certify -- and reports `resolved` True, so a caller reading it is not sent after
    # t_breakpoints it does not need.
    import magnus.adiabatic as ad
    import magnus.hamiltonians as hams

    h_vac = np.asarray(hams.hamiltonian_2nu_vacuum_energy_independent(
        PARAMS_2NU['sth'], PARAMS_2NU['Dm2']), dtype=complex)
    proj = np.diag([1.0, 0.0]).astype(complex)
    vcc = matter.vcc_func_from_rho_func(
        multi_resonance_ne(), 0.0, 1.0, 0.5, nubar=False,
        density_matter_is_in_g_per_cm3=False, density_is_of_number_of_electrons=True)

    def H_func(l):
        return (1.0/50.0e6)*h_vac + np.asarray(vcc(l))[..., None, None]*proj

    hinfo = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, _, certified = ad.hybrid_propagator(H_func, 0.0, L1, rtol=1e-13, atol=1e-13,
                                               max_iters=1, info=hinfo)
    assert certified is False
    assert hinfo['resolved'] is True


def test_strategy_info_reports_hybrid_certification():
    """certified is True where the hybrid strategy answered under 'auto' (it cannot be
    otherwise, since 'auto' declines an uncertified result) and can be False under 'hybrid'."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, L1, strategy_info=info)
    assert info['engine'] == 'hybrid' and info['certified'] is True

    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(step_ne(), 50.0e6, L1, strategy='hybrid', strategy_info=info)
    assert info['engine'] == 'hybrid' and info['certified'] is False


def test_strategy_info_is_free_and_optional():
    """Omitting it must change nothing about the answer -- the instrument is not allowed to
    move the measurement."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P_plain = np.asarray(call(solar_ne(), 10.0e6, L1))
        P_watched = np.asarray(call(solar_ne(), 10.0e6, L1, strategy_info={}))
    assert maxabs(P_plain - P_watched) == 0.0


def test_strategy_info_trace_carries_no_private_keys():
    """The trace is public output; the Hamiltonian payload the cross-check uses is not."""
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, L1, strategy='magnus', strategy_info=info)
    assert info['trace']
    for entry in info['trace']:
        assert not any(k.startswith('_') for k in entry)


# ----------------------------------------------------------------------
# an explicit `cumulative` must configure the scan, not disable three engines
# ----------------------------------------------------------------------

def test_passing_cumulative_does_not_disable_the_other_engines():
    """Regression: ``cumulative`` reached the dispatchers inside ``**kwargs``, where any
    unrecognized key makes them decline.  Passing the documented default ``'auto'`` therefore
    changed which engine answered (hybrid -> general ladder) and moved a 10 MeV solar single
    point by 9.3e-06."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        info_plain, info_auto = {}, {}
        P_plain = np.asarray(call(solar_ne(), 10.0e6, L1, strategy_info=info_plain))
        P_auto = np.asarray(call(solar_ne(), 10.0e6, L1, cumulative='auto',
                                 strategy_info=info_auto))
    assert info_plain['engine'] == info_auto['engine'] == 'hybrid'
    assert maxabs(P_plain - P_auto) == 0.0


def test_explicit_cumulative_true_still_reaches_the_cumulative_scan():
    """The other half of the same fix: cumulative=True names one engine and is documented to
    raise rather than be quietly substituted, so the engines tried before it stand aside."""
    info = {}
    Ls = np.linspace(0.05*L1, L1, 8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, Ls, cumulative=True, strategy_info=info)
    assert info['engine'] == 'cumulative'
    with pytest.raises(ValueError):
        call(solar_ne(), np.linspace(5.0e6, 50.0e6, 8), L1, cumulative=True)


# ----------------------------------------------------------------------
# cross_check_strategies
# ----------------------------------------------------------------------

def test_cross_check_runs_several_engines_and_reports_which():
    out = op.cross_check_strategies(
        op.osc_prob_matter_std_potential, 2, solar_ne(), 10.0e6, L1, PARAMS_2NU,
        L0=0.0, density_is_of_number_of_electrons=True)
    assert 'hybrid' in out['ran']
    assert len(out['ran']) >= 3
    assert set(out['families']) == set(out['ran'])
    for a, b in out['spread']:
        assert a in out['ran'] and b in out['ran']
    assert out['max_spread'] >= out['max_spread_independent']


def test_cross_check_declines_are_reported_not_hidden():
    """Most engines decline most requests.  A decline must be visible with a reason, since a
    cross-check that quietly compared two engines and called it three would be worthless."""
    out = op.cross_check_strategies(
        op.osc_prob_matter_std_potential, 2, solar_ne(), 10.0e6, L1, PARAMS_2NU,
        L0=0.0, density_is_of_number_of_electrons=True)
    assert out['declined']
    for label, reason in out['declined'].items():
        assert label not in out['ran']
        assert isinstance(reason, str) and reason


def test_cross_check_uses_expm_only_where_it_is_exact():
    """A constant density makes expm the exact answer; a varying one makes it a first-order
    approximation, and it must decline rather than supply a bad reference."""
    out = op.cross_check_strategies(
        op.osc_prob_matter_std_potential, 2, 0.05*NE0, 10.0e6, L1, PARAMS_2NU,
        L0=0.0, density_is_of_number_of_electrons=True)
    assert 'expm' in out['ran']
    assert out['max_spread'] < 1e-9

    out = op.cross_check_strategies(
        op.osc_prob_matter_std_potential, 2, solar_ne(), 10.0e6, L1, PARAMS_2NU,
        L0=0.0, density_is_of_number_of_electrons=True)
    assert 'expm' in out['declined']
    assert 'varies with position' in out['declined']['expm']


def test_cross_check_sees_the_unmarked_step_the_hybrid_path_gets_wrong():
    """The acceptance case, run against today's code: with ``strategy='hybrid'`` forced, the
    hybrid engine still returns the 0.46-wrong answer for the FINDINGS §3.1 step profile (it
    warns, and under 'auto' it is declined) -- and the cross-check must show the disagreement
    without any oracle.  Measured on the pre-fix package, where this answer was returned
    *silently* by the default entry point, the same comparison reports 5.4e-01."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = op.cross_check_strategies(
            op.osc_prob_matter_std_potential, 2, step_ne(), 50.0e6, L1, PARAMS_2NU,
            L0=0.0, density_is_of_number_of_electrons=True)
    assert {'hybrid', 'magnus'} <= set(out['ran'])
    assert out['max_spread_independent'] > 1e-2
    a, b = out['max_spread_independent_pair']
    assert op.ENGINE_FAMILIES[a] != op.ENGINE_FAMILIES[b]


def test_cross_check_rejects_an_unknown_engine_label():
    with pytest.raises(ValueError):
        op.cross_check_strategies(
            op.osc_prob_matter_std_potential, 2, solar_ne(), 10.0e6, L1, PARAMS_2NU,
            L0=0.0, density_is_of_number_of_electrons=True, engines=('hybrid', 'nonesuch'))


def test_cross_check_leaves_no_engine_disabled_after_a_raising_call():
    """cumulative=True raises by design on a request it cannot serve.  If that escaped the
    context manager it would leave dispatchers disabled for the rest of the session -- a
    diagnostic silently changing every later result."""
    op.cross_check_strategies(
        op.osc_prob_matter_std_potential, 2, solar_ne(),
        np.linspace(5.0e6, 50.0e6, 4), L1, PARAMS_2NU,
        L0=0.0, density_is_of_number_of_electrons=True)
    assert op._ENGINES_DISABLED == frozenset()
    assert op._ENGINE_TRACE is None
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, L1, strategy_info=info)
    assert info['engine'] == 'hybrid'


def test_a_short_baseline_scan_reaches_the_cumulative_engine():
    """The N = 25 seam was lowered to 8, and this pins why rather than only that.

    ``HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`` used to be 25 on the grounds that yielding
    earlier cost "several times slower ... to buy accuracy that was already two orders inside
    what the caller asked for".  Measured over 42 workloads, both halves were wrong outside
    solar: the cumulative scan is *cheaper* on median at every size (0.25x at N = 8) and its
    worst error over the 28 workloads it serves is 1.13e-07 against the hybrid path's 1.68e-03
    -- including two baseline scans at N = 8 that were outside the requested tolerance with no
    warning at all.

    Asserted as routing rather than accuracy so it needs no oracle: below the seam the hybrid
    path keeps the request, at and above it the cumulative engine takes it."""
    Ls_below = np.linspace(0.05*L1, L1, op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS - 1)
    Ls_at = np.linspace(0.05*L1, L1, op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS)

    info_below, info_at = {}, {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, Ls_below, strategy_info=info_below)
        call(solar_ne(), 10.0e6, Ls_at, strategy_info=info_at)
    assert info_below['engine'] == 'hybrid'
    assert info_at['engine'] == 'cumulative'

    # And the fall-through stays safe: the hybrid threshold must not drop below the one that
    # makes cumulative='auto' engage, or a scan could be declined by both and land on the
    # general per-point path, which is slower AND less accurate than either.
    assert (op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS
            >= op.CUMULATIVE_AUTO_MIN_POINTS)


@pytest.mark.parametrize('profile,d,params,energy', [
    ('solar', 2, PARAMS_2NU, 10.0e6),
    ('solar', 3, PARAMS_3NU, 50.0e6),
    ('multi-resonance', 2, PARAMS_2NU, 50.0e6),
    ('multi-resonance', 3, PARAMS_3NU, 50.0e6),
])
def test_engines_agree_across_a_small_profile_matrix(profile, d, params, energy):
    """The CI wiring for the cross-check: on the profiles the package is meant to serve, every
    pair of engines from *different* families must agree to 5e-3.

    The bound is measured rather than aspirational.  Over the matrix in
    ``docs/dev/adversarial_batteries/invariants.py`` the largest cross-family disagreement on
    these families is 3.8e-03, between the hybrid strategy and the energy-batched scan; both
    are inside their own requested 1e-3 of the truth, and they differ because they are
    different methods, not because either is wrong.  5e-3 leaves room for that without leaving
    room for the 4.6e-01 the same comparison reports on an unmarked discontinuity.
    """
    ne = solar_ne() if profile == 'solar' else multi_resonance_ne()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = op.cross_check_strategies(
            op.osc_prob_matter_std_potential, d, ne, energy, L1, params,
            L0=0.0, density_is_of_number_of_electrons=True)
    assert len(out['ran']) >= 2, out['declined']
    assert out['max_spread_independent'] < 5e-3, (
        out['max_spread_independent_pair'], out['spread'])
