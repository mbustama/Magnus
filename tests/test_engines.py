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


def test_strategy_info_reports_sampling_only_when_asked():
    """The sampling report is opt-in, because paying for it by default would be wrong.

    A Nyquist criterion is objectively correct and fires on 98 % of realistic scans -- 4400
    points would be needed on a solar trajectory and 73 000 on a supernova ray -- so it is
    reported rather than warned about (``adversarial_batteries/alias_fp.py``).  Reporting still
    costs eigenvalues, 5.5 % of the cheapest scan measured, so it runs only when the caller
    passed ``strategy_info``.  This pins both halves: the key is present when asked for, and the
    machinery is not invoked when it is not.
    """
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, np.linspace(0.2*L1, L1, 8), strategy_info=info)
    rep = info.get('sampling')
    assert rep, 'strategy_info did not carry a sampling report'
    for key in ('oscillation_length', 'cycles_over_trajectory', 'nyquist_points',
                'spacing', 'cycles_per_step', 'aliased'):
        assert key in rep, 'sampling report is missing %r' % key
    # A solar trajectory is thousands of oscillations long, so an 8-point scan is aliased and
    # Nyquist would want far more points than anyone would ask for.
    assert rep['aliased'] is True
    assert rep['nyquist_points'] > 1000

    # And it is not computed when nobody asked: the report is the only consumer, so if it were
    # running unconditionally the guard would be the thing that broke.
    called = []
    orig = op._sampling_report

    def spy(*a, **k):
        called.append(1)
        return orig(*a, **k)

    op._sampling_report = spy
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            call(solar_ne(), 10.0e6, np.linspace(0.2*L1, L1, 8))
        assert not called, 'the sampling report ran without strategy_info being requested'
    finally:
        op._sampling_report = orig


def test_hybrid_does_not_stand_aside_for_a_disabled_engine():
    """The hybrid path must not yield to the cumulative scan when the caller switched it off.

    ``_cumulative_scan_would_serve``'s docstring argues the fall-through is safe: the hybrid
    dispatcher's threshold is the larger of the two, so whenever it declines on point count,
    ``cumulative='auto'`` is guaranteed to accept and "a scan can never be declined by both and
    land on the general per-point path -- slower AND less accurate than either, and silently
    so."

    That argument holds only while the cumulative scan is *available*.  ``cumulative=False``
    switches it off without touching the hybrid path, so above the seam the hybrid dispatcher
    stood aside for an engine that could not run; ip_exp needs every baseline equal and the
    separable engine a single shared baseline, so both declined too, and the request landed on
    exactly the general ladder the docstring rules out.

    Measured before the fix, on a tagged exponential at d = 2 and 10 MeV: adding one baseline to
    a seven-point scan moved the engine from hybrid to magnus and the error from 1.157e-05 to
    2.966e-03 -- a factor of 256, and outside the requested 1e-3.

    Asserted as routing rather than accuracy so it needs no oracle.
    """
    seam = op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS
    Ls_below = np.linspace(0.2*L1, L1, max(seam - 1, 1))
    Ls_at = np.linspace(0.2*L1, L1, seam)
    Ls_above = np.linspace(0.2*L1, L1, seam + 8)

    for Ls in (Ls_below, Ls_at, Ls_above):
        info = {}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            call(solar_ne(), 10.0e6, Ls, cumulative=False, strategy_info=info)
        assert info['engine'] == 'hybrid', (
            'with cumulative=False the hybrid path stood aside at N=%d and %r answered instead; '
            'the cumulative scan it yielded to was disabled by the caller'
            % (len(Ls), info['engine']))

    # The default path must be untouched: 'auto' still yields at and above the seam, which is
    # the behaviour HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS exists to produce.
    info_below, info_at = {}, {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, Ls_below, strategy_info=info_below)
        call(solar_ne(), 10.0e6, Ls_at, strategy_info=info_at)
    assert info_below['engine'] == 'hybrid'
    assert info_at['engine'] == 'cumulative'

    # And an explicit cumulative=True still names one engine and is never substituted for.
    info_true = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        call(solar_ne(), 10.0e6, Ls_at, cumulative=True, strategy_info=info_true)
    assert info_true['engine'] == 'cumulative'


# ----------------------------------------------------------------------
# cross_check_strategies: a spread of zero that means nothing was compared
# ----------------------------------------------------------------------

def test_cross_check_warns_when_no_engine_ran():
    """``osc_prob`` has no ``strategy`` parameter, so every engine declines -- and it is the
    entry point a reader is most likely to reach for, since four of the notebooks call it
    directly.  The result then carries max_spread = 0.0, which is exactly what perfect
    agreement looks like.  The numbers cannot distinguish the two cases, so the warning has
    to."""
    import magnus.hamiltonians as hamiltonians

    h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**PARAMS_3NU))
    vcc = matter.vcc_func_from_rho_func(solar_ne(), density_is_of_number_of_electrons=True)

    def H_func(l):
        return h_vac/10.0e6 + np.diag([vcc(l), 0.0, 0.0])

    with pytest.warns(op.CrossCheckInconclusiveWarning, match="no engine ran"):
        out = op.cross_check_strategies(op.osc_prob, H_func, 0.0, L1)

    assert out['ran'] == ()
    assert out['max_spread'] == 0.0
    assert out['declined']


def test_cross_check_warns_when_only_one_family_ran():
    """Two engines of the same family agreeing is the self-certification the cross-check
    exists to avoid relying on.  ``max_spread_independent`` is 0.0 there because no
    cross-family pair exists, not because independent methods agreed."""
    with pytest.warns(op.CrossCheckInconclusiveWarning, match="family"):
        out = op.cross_check_strategies(
            op.osc_prob_matter_std_potential, 2, solar_ne(), 10.0e6, L1, PARAMS_2NU,
            L0=0.0, density_is_of_number_of_electrons=True,
            engines=('magnus', 'cumulative'))

    assert set(out['ran']) == {'magnus', 'cumulative'}
    assert out['max_spread_independent'] == 0.0
    assert out['max_spread_independent_pair'] is None
    assert len({op.ENGINE_FAMILIES[lab] for lab in out['ran']}) == 1


def test_cross_check_is_quiet_when_it_actually_compared_something():
    """The ordinary case must not warn, or the warning is noise: several engines, more than
    one family, and a spread that means what it says."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", op.CrossCheckInconclusiveWarning)
        out = op.cross_check_strategies(
            op.osc_prob_matter_std_potential, 2, solar_ne(), 10.0e6, L1, PARAMS_2NU,
            L0=0.0, density_is_of_number_of_electrons=True)

    assert len(out['ran']) >= 2
    assert len({op.ENGINE_FAMILIES[lab] for lab in out['ran']}) >= 2
    assert out['max_spread_independent_pair'] is not None


# ----------------------------------------------------------------------
# The 'constant' engine: a position-independent Hamiltonian answered in one
# batched exponential instead of one osc_prob call per point.
#
# It must agree with the per-point route it replaces, *for every scenario
# wrapper*, and the reason is a trap this engine walked straight into.  Two of
# the three dispatch call sites rebound ``h_matt`` to ``VCC_func*h_matt``
# before calling the dispatcher, while the dispatcher is documented to take
# ``h_matt`` as "the constant matrix multiplying VCC_func(l)" -- and it goes on
# to multiply by VCC itself.  The separable engine never noticed, because the
# rebinding only happens on the constant branch it used to decline outright.
#
# The resulting error was invisible in exactly the ways this package fears:
# VCC^2 is ~1e-25 rather than ~1e-13, so the matter term all but vanished and
# the answer stayed a perfectly plausible, unitary, nearly-vacuum probability;
# and because VCC^2 has no sign, the neutrino and antineutrino results became
# *identical*.  Constant-density standard matter was unaffected (that call site
# passes the bare projector), so a test of the headline case alone passed.
# ----------------------------------------------------------------------

CONST_RHO = 10.0*gd.UNIT_G_PER_CM3
CONST_E = 1.0*gd.UNIT_GEV
CONST_L = 1000.0*gd.UNIT_KM
CONST_ESCAN = np.linspace(0.6, 8.0, 12)*gd.UNIT_GEV


def _constant_density_wrappers():
    r"""One call per scenario wrapper that can reach the constant engine."""
    return {
        'matter std 2nu': lambda E, **kw: op.osc_prob_2nu_matter_constant_density(
            E, CONST_L, CONST_RHO, validate_input=False, **PARAMS_2NU, **kw),
        'matter std 3nu': lambda E, **kw: op.osc_prob_3nu_matter_constant_density(
            E, CONST_L, CONST_RHO, validate_input=False, **PARAMS_3NU, **kw),
        'nsi 3nu': lambda E, **kw: op.osc_prob_3nu_matter_nsi_constant_density(
            E, CONST_L, CONST_RHO, eps_ee=0.3, eps_em=0.1j, eps_et=0.0,
            eps_mm=0.0, eps_mt=0.0, eps_tt=0.0, validate_input=False,
            **PARAMS_3NU, **kw),
        'liv 2nu': lambda E, **kw: op.osc_prob_2nu_matter_liv_constant_density(
            E, CONST_L, CONST_RHO, sxi=0.2, b1=gd.B1, b2=gd.B2,
            Lambda=gd.LAMBDA, n_liv=1, validate_input=False, **PARAMS_2NU, **kw),
        'liv 3nu': lambda E, **kw: op.osc_prob_3nu_matter_liv_constant_density(
            E, CONST_L, CONST_RHO, b1=gd.B1, b2=gd.B2, b3=gd.B3,
            Lambda=gd.LAMBDA, n_liv=1, validate_input=False, **PARAMS_3NU, **kw),
        'vacuum 3nu': lambda E, **kw: op.osc_prob_3nu_vacuum(
            E, CONST_L, validate_input=False, **PARAMS_3NU, **kw),
    }


@pytest.mark.parametrize('name', list(_constant_density_wrappers()))
@pytest.mark.parametrize('nubar', [False, True])
def test_constant_engine_agrees_with_the_refinement_ladder(name, nubar):
    r"""Every wrapper, both signs, single point and scan, against the *ladder*.

    Disabling ``'constant'`` alone is not enough to reach the ladder, and the
    first version of this test did exactly that.  ``osc_prob`` has its own
    constant-Hamiltonian shortcut, gated on ``not kwargs``, which then answers
    instead -- so both sides of the comparison were one batched exponential
    through the same ``magnus._expm_stack``, agreed to 12 digits by construction,
    and the tolerance could never fire.  Instrumented, the slab machinery was
    entered **zero** times on either side.

    The lever that does defeat it is a keyword ``osc_prob`` does not name but the
    Magnus core accepts, ``expm_backend`` -- ``n_slabs`` is a *named* parameter
    and so never reaches the ``kwargs`` the guard tests.  With it the reference is
    a genuine slab-quadrature product, which for a constant Hamiltonian is the
    same answer reached by a different route, and that is what makes it a useful
    check.  The test asserts the slab path was entered rather than trusting it.
    """
    fn = _constant_density_wrappers()[name]
    for E in (CONST_E, CONST_ESCAN):
        slab_calls = []
        real = op.compute_evolution_operator_multiple_slabs

        def counted(*a, **k):
            slab_calls.append(1)
            return real(*a, **k)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with op._engine_probe() as trace:
                new = np.asarray(fn(E, nubar=nubar), dtype=float)
            op.compute_evolution_operator_multiple_slabs = counted
            try:
                with op._engine_probe(disabled=('constant',)):
                    old = np.asarray(fn(E, nubar=nubar, expm_backend='numba'),
                                     dtype=float)
            finally:
                op.compute_evolution_operator_multiple_slabs = real
        assert any(t['engine'] == 'constant' for t in trace), (
            "the constant engine declined for %s, so this comparison is void" % name)
        assert slab_calls, (
            "%s: the reference route never entered the slab machinery, so it is "
            "another constant shortcut rather than the ladder" % name)
        assert maxabs(new - old) < 1e-12, (name, nubar, maxabs(new - old))
        assert maxabs(np.sum(new, axis=-1) - 1.0) < 1e-12


@pytest.mark.parametrize('name', ['matter std 2nu', 'matter std 3nu', 'nsi 3nu',
                                  'liv 2nu', 'liv 3nu'])
def test_constant_engine_keeps_the_antineutrino_sign(name):
    r"""nubar must still change the answer once matter is present.

    The matter potential carries the neutrino/antineutrino sign, so this is the
    assertion that fails if a call site's ``h_matt`` has had ``VCC`` folded into
    it and the engine multiplies by ``VCC`` a second time: ``VCC**2`` is
    sign-blind, and the two answers come back bit-identical.
    """
    fn = _constant_density_wrappers()[name]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # The engine guard is the point.  Without it this test passed at full
        # strength while exercising the route the constant engine replaced --
        # measured identical to 17 digits with the engine disabled -- so the one
        # test written to catch the VCC-folded-twice bug was not testing the code
        # that had the bug.  See the module docstring of test_expm_backend.py.
        with op._engine_probe() as trace:
            P_nu = np.asarray(fn(CONST_E, nubar=False), dtype=float)
            P_nubar = np.asarray(fn(CONST_E, nubar=True), dtype=float)
    assert any(t['engine'] == 'constant' for t in trace), (
        "%s: the constant engine did not answer, so this test is not looking at "
        "the code it is meant to protect" % name)
    assert maxabs(P_nu - P_nubar) > 1e-6, (
        "%s: nubar had no effect, which is what folding VCC into h_matt twice "
        "looks like" % name)


def test_constant_engine_does_not_take_a_position_dependent_profile():
    r"""PREM and an exponential profile must keep their own engines.

    A constant-H engine that captured a varying profile would propagate the whole
    trajectory with one exponential of one Hamiltonian -- wrong by O(1) while
    still unitary.
    """
    import magnus.earth as earth
    costhz = -0.85
    Lc = earth.distance_traveled_inside_earth(costhz)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for label, call in (
            ('PREM single point', lambda info: op.osc_prob_3nu_earth(
                CONST_E, costhz=costhz, L=Lc, validate_input=False,
                strategy_info=info)),
            ('PREM scan', lambda info: op.osc_prob_3nu_earth(
                CONST_ESCAN, costhz=costhz, L=Lc, validate_input=False,
                strategy_info=info)),
            ('exponential profile', lambda info: op.osc_prob_3nu_sun(
                10.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
                validate_input=False, strategy_info=info)),
        ):
            info = {}
            call(info)
            assert info.get('engine') != 'constant', (
                "%s was answered by the constant engine" % label)


def test_constant_appears_in_the_cross_check_tables():
    r"""A new engine that is missing from these tables is invisible to the cross-check.

    ``_CROSS_CHECK_FORCING`` must list it in *every other* row's forbid set as
    well as its own: it answers before ``osc_prob_energy_baseline``, which is
    what records the payload the independent ``expm`` reference is built from, so
    an unlisted constant engine silently removes the only non-Magnus oracle in
    the table.
    """
    assert 'constant' in op.ENGINE_FAMILIES
    assert 'constant' in op._CROSS_CHECK_FORCING
    for label, (_, _, forbid) in op._CROSS_CHECK_FORCING.items():
        if label != 'constant':
            assert 'constant' in forbid, (
                "%r does not forbid 'constant', so forcing it cannot work" % label)
