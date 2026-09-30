"""Issue #160: every public argument is checked once per call, by one set of rules.

One case per item of the issue's checklist.  Each bad value must be refused with a
ValueError or TypeError whose message names the argument and the public function the caller
called; each valid edge value must be accepted.  Refusals return before any propagation, so
the whole module runs in a few seconds.
"""

import warnings

import numpy as np
import pytest

import magnus.adiabatic as ad
import magnus.avgprob as ap
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.magnus as mg
import magnus.matter as matter
import magnus.oscprob as op
import magnus.oscprobstd as std
from magnus._validate import InputTypeError

KM = gd.UNIT_KM
E = 1.0*gd.UNIT_GEV
L = 1300.0*KM
RHO = 2.8
OSC = dict(s12=0.55, s23=0.75, s13=0.15, dCP=1.2, D21=7.4e-5, D31=2.5e-3)
H3 = np.diag([0.0, 1e-13, 3e-13]).astype(complex)


def vac(**kw):
    return op.osc_prob_3nu_vacuum(E, L, **kw)


def const(**kw):
    kw.setdefault('density_matter_is_in_g_per_cm3', True)
    return op.osc_prob_3nu_matter_constant_density(E, L, RHO, **kw)


def expd(**kw):
    kw.setdefault('density_matter_is_in_g_per_cm3', True)
    return op.osc_prob_3nu_matter_exp_density(E, L, 0.0, RHO, 1000.0*KM, **kw)


def hyb2(**kw):
    # A request the adiabatic engine answers before the general ladder is reached.
    return op.osc_prob_2nu_matter_exp_density(
        3*gd.UNIT_MEV, 6e5*KM, 0.0, rho_central=100*gd.UNIT_G_PER_CM3, l_scale=7e4*KM,
        sth=0.55, Dm2=7.5e-5, nu_i=0, nu_f=0, **kw)


REFUSED = {
    # §1 physics arguments
    'energy nan entry': (lambda: op.osc_prob_3nu_vacuum(np.array([E, np.nan]), L),
                         'energy', 'osc_prob_3nu_vacuum'),
    'energy inf entry': (lambda: op.osc_prob_3nu_vacuum(np.array([E, np.inf]), L),
                         'energy', 'osc_prob_3nu_vacuum'),
    'L nan entry': (lambda: op.osc_prob_3nu_vacuum(E, np.array([L, np.nan])), 'L', 'vacuum'),
    'L below L0': (lambda: op.osc_prob_3nu_vacuum(E, -L), 'L must be >= L0', 'vacuum'),
    'L True': (lambda: op.osc_prob_3nu_vacuum(E, True), 'L', 'vacuum'),
    'masked L': (lambda: op.osc_prob_3nu_vacuum(E, np.ma.masked_array([L, L], [0, 1])),
                 'masked', 'vacuum'),
    'dCP nan': (lambda: vac(dCP=np.nan), 'dCP', 'osc_prob_3nu_vacuum'),
    'D31 True': (lambda: vac(D31=True), 'D31', 'osc_prob_3nu_vacuum'),
    'rho negative': (lambda: op.osc_prob_3nu_matter_constant_density(E, L, -1.0),
                     'rho must', 'constant_density'),
    'electron_fraction 1.5': (lambda: const(electron_fraction=1.5), 'electron_fraction',
                              'constant_density'),
    'Sun average L0 > L': (lambda: op.osc_prob_3nu_sun(
        5*gd.UNIT_MEV, 0.5*gd.SUN_RADIUS*KM, 0.8*gd.SUN_RADIUS*KM, nu_i=0, nu_f=0,
        average=True), 'L must be >= L0', 'osc_prob_3nu_sun'),
    # §2 parameter dicts
    'osc_params typo': (lambda: op.osc_prob_matter_std_potential(
        3, RHO, E, L, dict(OSC, dcp=0.5), density_matter_is_in_g_per_cm3=True),
        "did you mean 'dCP'", 'osc_prob_matter_std_potential'),
    'osc_params list': (lambda: op.osc_prob_matter_std_potential(
        3, RHO, E, L, list(OSC.values()), density_matter_is_in_g_per_cm3=True),
        'osc_params must be a dict', 'osc_prob_matter_std_potential'),
    'complex eps_ee': (lambda: op.osc_prob_3nu_matter_nsi_constant_density(
        E, L, RHO*gd.UNIT_G_PER_CM3, eps_ee=0.1+0.2j), 'eps_ee', 'nsi_constant_density'),
    'complex eps_s1s1': (lambda: op.osc_prob_5nu_matter_nsi_constant_density(
        E, L, RHO*gd.UNIT_G_PER_CM3, eps_s1s1=0.2j), 'eps_s1s1', 'osc_prob_5nu'),
    'complex b1': (lambda: op.osc_prob_3nu_vacuum_liv(E, L, b1=1e-13+1e-13j, Lambda=1e9),
                   'b1', 'osc_prob_3nu_vacuum_liv'),
    'Lambda nan': (lambda: op.osc_prob_3nu_vacuum_liv(E, L, b1=1e-13, Lambda=np.nan),
                   'Lambda', 'osc_prob_3nu_vacuum_liv'),
    'n_liv -1': (lambda: op.osc_prob_3nu_vacuum_liv(E, L, b1=1e-13, Lambda=1e9, n_liv=-1),
                 'n_liv', 'osc_prob_3nu_vacuum_liv'),
    # §3 Earth
    'costhz -1.5 with L': (lambda: op.osc_prob_3nu_earth(E, costhz=-1.5, L=100*KM),
                           'costhz', 'osc_prob_3nu_earth'),
    'L beyond chord': (lambda: op.osc_prob_3nu_earth(E, -0.8, L=12000*KM), 'chord',
                       'osc_prob_3nu_earth'),
    'latitude 146': (lambda: op.osc_prob_3nu_earth(
        E, loc_ini=((146, 14, 0), (6, 3, 0)), loc_fin=((42, 25, 0), (13, 30, 0))),
        'latitude', 'osc_prob_3nu_earth'),
    'helper costhz': (lambda: earth.distance_traveled_inside_earth(-1.5), 'costhz',
                      'distance_traveled_inside_earth'),
    # §4 bools, ints, strings
    'average str': (lambda: vac(nu_i=1, nu_f=1, average='False'), 'average', 'vacuum'),
    'nubar str': (lambda: vac(nubar='yes'), 'nubar', 'vacuum'),
    'nu_i float': (lambda: vac(nu_i=1.0, nu_f=0), 'nu_i', 'vacuum'),
    'cumulative 2': (lambda: op.osc_prob_energy_baseline(H3, E, L, cumulative=2),
                     'cumulative', 'osc_prob_energy_baseline'),
    'integration_method': (lambda: vac(integration_method='foo'), 'integration_method',
                           'vacuum'),
    'strategy_info list': (lambda: vac(strategy_info=[]), 'strategy_info', 'vacuum'),
    # §5 refinement: independent of the engine that answers
    'max_n_slabs -1 (adiabatic)': (lambda: hyb2(max_n_slabs=-1), 'max_n_slabs',
                                   'osc_prob_2nu_matter_exp_density'),
    'rtol 0': (lambda: expd(rtol=0.0), 'rtol', 'exp_density'),
    'rtol True': (lambda: expd(rtol=True), 'rtol', 'exp_density'),
    'rtol nan': (lambda: expd(rtol=np.nan), 'rtol', 'exp_density'),
    'n_slabs 2.5': (lambda: expd(n_slabs=2.5), 'n_slabs', 'exp_density'),
    'growth 1.0': (lambda: expd(growth_factor_n_slabs=1.0), 'growth_factor_n_slabs',
                   'exp_density'),
    'magnus_exp_order True': (lambda: expd(magnus_exp_order=True), 'magnus_exp_order',
                              'exp_density'),
    'n_jobs 0': (lambda: vac(n_jobs=0), 'n_jobs', 'vacuum'),
    # §6 slab edges and breakpoints
    'edges gap': (lambda: expd(t_slab_edges=[[0, L/3], [L/2, L]]), 'gap', 'exp_density'),
    'edges short': (lambda: expd(t_slab_edges=[[0, L/2]]), 'span the whole path',
                    'exp_density'),
    'breakpoints nan': (lambda: expd(t_breakpoints=[np.nan]), 't_breakpoints',
                        'exp_density'),
    # §7 Hamiltonians and profiles
    'H None': (lambda: op.osc_prob_energy_baseline(None, E, L), 'H_func',
               'osc_prob_energy_baseline'),
    'H not Hermitian': (lambda: op.osc_prob_energy_baseline(np.triu(np.ones((3, 3)))*1e-13,
                        E, L), 'Hermitian', 'osc_prob_energy_baseline'),
    'H returns list': (lambda: op.osc_prob_energy_baseline(
        lambda e, l: (H3/e).tolist(), E, L), 'NumPy array', 'osc_prob_energy_baseline'),
    'rho_func nan partway': (lambda: op.osc_prob_matter_std_potential(
        3, lambda l: np.nan if l > L/2 else 2.7, E, L, OSC,
        density_matter_is_in_g_per_cm3=True), 'rho_func', 'osc_prob_matter_std_potential'),
    # §8 avgprob
    'avg constant H baseline': (lambda: ap.averaged_probabilities_constant_hamiltonian(
        H3, -100.0), 'baseline', 'averaged_probabilities_constant_hamiltonian'),
    'avg numerically energy': (lambda: ap.averaged_probabilities_numerically(
        lambda e: 0.5, -1.0, 0.1, 5), 'energy', 'averaged_probabilities_numerically'),
    'coherence phase_scale': (lambda: ap.coherence_blocks([0.0, 1.0], 0.0), 'phase_scale',
                              'coherence_blocks'),
    # §9 adiabatic
    'hybrid l0 > l1': (lambda: ad.hybrid_propagator(
        lambda l: np.array([[l-5, .3], [.3, 5-l]], dtype=complex), 10.0, 0.0), 'l1',
        'hybrid_propagator'),
    'fd_step_frac 0': (lambda: ad.find_resonance_candidates(
        lambda l: np.array([[l-5, .3], [.3, 5-l]], dtype=complex), 0.0, 10.0,
        fd_step_frac=0.0), 'fd_step_frac', 'find_resonance_candidates'),
    # §10 core
    'gl_nodes 0': (lambda: mg.gl_nodes(0), 'order', 'gl_nodes'),
    'multislab flat edges': (lambda: mg.magnus_expansion_multislab(
        lambda t: 1j*np.eye(2), [0, 1, 2]), 't_slab_edges', 'magnus_expansion_multislab'),
    # §11 builders, helpers, references
    'builder energy 0': (lambda: hams.hamiltonian_3nu_vacuum(0.0, **OSC), 'energy',
                         'hamiltonian_3nu_vacuum'),
    'builder complex eps_ee': (lambda: hams.hamiltonian_3nu_nsi(1e-13, 0.1j, 0, 0, 0, 0, 0),
                               'eps_ee', 'hamiltonian_3nu_nsi'),
    'exp profile l_scale 0': (lambda: matter.exp_density_profile(3.0, 0.0), 'l_scale',
                              'exp_density_profile'),
    'reference energy negative': (lambda: std.osc_prob_2nu_vacuum_std(0.5, 7e-5, -E, L),
                                  'energy', 'osc_prob_2nu_vacuum_std'),
    # §14 cross_check_strategies
    'cross-check engines []': (lambda: op.cross_check_strategies(
        op.osc_prob_3nu_vacuum, E, L, engines=[]), 'engines', 'cross_check_strategies'),
}


@pytest.mark.parametrize('case', sorted(REFUSED))
def test_bad_value_is_refused_by_name(case):
    call, argument, function = REFUSED[case]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises((ValueError, TypeError)) as info:
            call()
    message = str(info.value)
    assert argument in message, message
    assert function in message, message


def test_wrong_type_is_also_a_value_error():
    # Code written against the ValueError that every check raised before 1.2.0 keeps working.
    with pytest.raises(ValueError):
        vac(average='False')
    with pytest.raises(InputTypeError):
        vac(average='False')


ACCEPTED = {
    'np.float32 energy': lambda: op.osc_prob_3nu_vacuum(np.float32(1e9), L),
    'np.int64 phase': lambda: vac(dCP=np.int64(1)),
    '0-d array density': lambda: op.osc_prob_3nu_matter_constant_density(
        E, L, np.array(RHO), density_matter_is_in_g_per_cm3=True),
    'np.where profile': lambda: op.osc_prob_matter_std_potential(
        3, lambda l: np.where(l < L/2, 3.0, 8.0), E, L, OSC,
        density_matter_is_in_g_per_cm3=True),
    'real part of a complex coupling': lambda: op.osc_prob_3nu_matter_nsi_constant_density(
        E, L, RHO*gd.UNIT_G_PER_CM3, eps_ee=0.1+0j),
    'complex off-diagonal coupling': lambda: op.osc_prob_3nu_matter_nsi_constant_density(
        E, L, RHO*gd.UNIT_G_PER_CM3, eps_em=0.1+0.2j),
    'partial Earth path': lambda: op.osc_prob_3nu_earth(E, -0.8, L=5000*KM),
    'NuFIT any case': lambda: gd.load_nufit_params('nufit 6.1'),
    'zero average_spread': lambda: vac(average=True, average_spread=0.0),
    'reference Dm2 = 0': lambda: std.osc_prob_2nu_matter_std(0.55, 0.0, 1e-13, E, L),
}


@pytest.mark.parametrize('case', sorted(ACCEPTED))
def test_valid_edge_value_is_accepted(case):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = ACCEPTED[case]()
    assert np.all(np.isfinite(np.asarray(out if not isinstance(out, dict) else
                                         list(out.values()), dtype=float)))


def test_numpy_scalars_give_the_float_answer():
    assert np.array_equal(np.asarray(op.osc_prob_3nu_vacuum(np.float64(1e9), L)),
                          np.asarray(vac()))
    assert np.array_equal(np.asarray(vac(dCP=np.int64(1))), np.asarray(vac(dCP=1.0)))


def test_complex_off_diagonal_nsi_stays_unitary():
    P = op.osc_prob_3nu_matter_nsi_constant_density(E, L, RHO*gd.UNIT_G_PER_CM3,
                                                    eps_em=0.1+0.2j, eps_mt=0.05-0.1j)
    assert np.max(np.abs(np.asarray(P).sum(axis=1) - 1.0)) < 1e-12


def test_energy_unit_warning():
    with pytest.warns(gd.EnergyUnitWarning):
        op.osc_prob_3nu_vacuum(10.0, L)
