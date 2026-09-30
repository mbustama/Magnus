r"""Keywords that were accepted and then ignored now either work or raise (issues #110, #112, #114).

* #110: ``default_osc_params_set_name`` on a two-flavor call, where no parameter set applies.
* #112: ``density_matter_is_in_g_per_cm3`` and ``density_is_of_number_of_electrons`` together,
  which give ``rho`` two different units; the electron flag used to win in silence.
* #114: ``strategy_info`` on ``osc_prob_energy_baseline``, rejected on every route but
  ``average=True``, which left it empty and also ignored any misspelled keyword.
"""

import inspect
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.oscprob as op

E = 1.0*gd.UNIT_GEV
L = 1000.0*gd.UNIT_KM
KM = gd.UNIT_KM
OSC = gd.load_nufit_params('NuFIT 6.1')


def _call_args(name):
    """Positional and keyword arguments for a valid call of wrapper ``name``."""
    kw = {}
    if 'earth' in name:
        args = [E]
        kw.update(costhz=-0.5, L=6371.0*KM)
    elif 'sun' in name:
        args = [E, 0.5*gd.SUN_RADIUS*KM, 0.0]
    elif 'exp' in name:
        args = [E, L, 0.0, 3.0, 500.0*KM]
    else:
        args = [E, L]
    if 'constant' in name:
        kw['rho'] = 3.0
    if name.startswith('osc_prob_2nu_'):
        kw.update(sth=0.5, Dm2=2.5e-3)
    return args, kw


TWO_FLAVOR = sorted(n for n in dir(op) if n.startswith('osc_prob_2nu_')
                    and 'default_osc_params_set_name' not in inspect.signature(
                        getattr(op, n)).parameters
                    and 'kwargs' in inspect.signature(getattr(op, n)).parameters)


def test_every_two_flavor_wrapper_is_covered():
    assert len(TWO_FLAVOR) == 14


@pytest.mark.parametrize('name', TWO_FLAVOR)
def test_two_flavor_rejects_parameter_set_name(name):
    args, kw = _call_args(name)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        getattr(op, name)(*args, **kw)
        with pytest.raises(ValueError, match='no parameter set applies at 2 flavors'):
            getattr(op, name)(*args, **kw,
                              default_osc_params_set_name='OSC_PARAMS_NU_FIT_5_2_SK_IO')


def test_three_flavor_still_takes_parameter_set_name():
    P = op.osc_prob_3nu_vacuum(E, L, default_osc_params_set_name='OSC_PARAMS_NU_FIT_5_2_SK_IO')
    assert P.shape == (3, 3)


BOTH_UNITS = sorted(n for n in dir(op) if n.startswith('osc_prob') and callable(getattr(op, n))
                    and 'density_is_of_number_of_electrons' in inspect.signature(
                        getattr(op, n)).parameters)
SCENARIO = {
    'osc_prob_matter_std_potential': {},
    'osc_prob_matter_nsi': {'nsi_params': dict(eps_ee=0.0, eps_em=0.05, eps_et=0.0,
                                               eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)},
    'osc_prob_liv': {'liv_params': dict(sxi12=0.0, sxi23=0.0, sxi13=0.0, dxiCP=0.0, b1=0.0,
                                        b2=0.0, b3=1.0e-23, Lambda=1.0, n_liv=1)},
}


def test_every_density_unit_function_is_covered():
    assert len(BOTH_UNITS) == 27


@pytest.mark.parametrize('name', BOTH_UNITS)
def test_both_density_units_rejected(name):
    if name in SCENARIO:
        args, kw = [], dict(num_flavors=3, rho_func=lambda l: 3.0, energy=E, L=L,
                            osc_params=OSC, **SCENARIO[name])
    else:
        args, kw = _call_args(name)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises(ValueError, match='are both True'):
            getattr(op, name)(*args, **kw, density_matter_is_in_g_per_cm3=True,
                              density_is_of_number_of_electrons=True)


def test_each_density_unit_alone_unchanged():
    # The values quoted in #112: the g/cm^3 reading, and rho = 3 eV^3 of electrons (vacuum).
    P_g = op.osc_prob_3nu_matter_constant_density(E, L, rho=3.0,
                                                  density_matter_is_in_g_per_cm3=True)
    P_e = op.osc_prob_3nu_matter_constant_density(E, L, rho=3.0,
                                                  density_is_of_number_of_electrons=True)
    assert abs(P_g[1][0] - 0.013475) < 5.0e-7
    assert abs(P_e[1][0] - 0.003773) < 5.0e-7


H0 = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
    **{k: OSC[k] for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')})
ROUTES = {
    'point': (lambda **k: op.osc_prob_energy_baseline(H0/E, E, L, **k), 'magnus'),
    'average': (lambda **k: op.osc_prob_energy_baseline(H0/E, E, 1.0e8*KM, average=True, **k),
                'average'),
    'operator': (lambda **k: op.osc_prob_energy_baseline(H0/E, E, L,
                                                         return_evolution_operator=True, **k),
                 'magnus'),
}


@pytest.mark.parametrize('route', sorted(ROUTES))
def test_energy_baseline_fills_strategy_info(route):
    call, engine = ROUTES[route]
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with_info = call(strategy_info=info)
        without = call()
    assert info['engine'] == engine
    for a, b in zip(with_info if isinstance(with_info, tuple) else (with_info,),
                    without if isinstance(without, tuple) else (without,)):
        assert np.array_equal(np.asarray(a), np.asarray(b))


def test_energy_baseline_average_rejects_unknown_keyword():
    with pytest.raises(ValueError, match="unrecognized keyword argument.*'rtoll'"):
        op.osc_prob_energy_baseline(H0/E, E, 1.0e8*KM, average=True, rtoll=1.0e-6)


@pytest.mark.parametrize('call, bad, listed', [
    (lambda **kw: op.osc_prob_3nu_vacuum(E, L, **kw), 'theta12', 's12, s23, s13, dCP'),
    (lambda **kw: op.osc_prob_3nu_matter_nsi_constant_density(E, L, 3.0, density_matter_is_in_g_per_cm3=True, **kw),
     'eps_mue', 'eps_ee, eps_em, eps_et'),
])
def test_unknown_keyword_lists_wrapper_physics_keywords(call, bad, listed):
    """#152: a misspelled physics keyword names the wrapper's own keywords, not only the engine's."""
    with pytest.raises(ValueError) as err:
        call(**{bad: 0.1})
    msg = str(err.value)
    assert "'" + bad + "'" in msg
    assert listed in msg
    assert 'validate_input' not in msg and 'n_slabs, n_tpts' not in msg.split('engine keywords')[0]


def test_small_parallel_scan_runs_in_the_calling_process(monkeypatch):
    """#155 §3: n_jobs=2 on a millisecond 5-point scan cost about 1.4 s starting workers."""
    H0v = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC)
    H = lambda En, l: H0v/En + (1e-13*np.exp(-np.asarray(l, float)/(300*KM)))[..., None, None]*np.diag([1., 0, 0])
    E5 = np.linspace(1.0, 5.0, 5)*gd.UNIT_GEV
    started = []
    real = op.Parallel
    monkeypatch.setattr(op, 'Parallel', lambda *a, **k: started.append(k) or real(*a, **k))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        serial = op.osc_prob_energy_baseline(H, E5, L, n_jobs=1)
        par = op.osc_prob_energy_baseline(H, E5, L, n_jobs=2)
    assert started == []
    assert np.array_equal(serial, par)


def test_heavy_parallel_scan_still_uses_the_workers(monkeypatch):
    import time
    H0v = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC)
    def H(En, l):
        time.sleep(0.05)
        return H0v/En + 0.0*np.asarray(l, float)[..., None, None]
    started = []
    real = op.Parallel
    monkeypatch.setattr(op, 'Parallel', lambda *a, **k: started.append(k) or real(*a, **k))
    monkeypatch.setattr(op, 'N_JOBS_MIN_PARALLEL_WORK_S', 0.01)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        op.osc_prob_energy_baseline(H, np.linspace(1.0, 2.0, 3)*gd.UNIT_GEV, L, n_jobs=2)
    assert len(started) == 1


def test_parallel_and_batched_scans_agree_within_the_tolerance_as_documented():
    """#166 §1: n_jobs > 1 takes the per-point path, one ladder per point; it agrees with the
    batched scan to within the tolerance, not bit for bit, and the docstring says so."""
    import magnus.earth as earth
    E = np.geomspace(0.5, 20.0, 12)*gd.UNIT_GEV
    Lc = earth.distance_traveled_inside_earth(-0.7)*gd.UNIT_KM
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = np.asarray(op.osc_prob_3nu_earth(E, costhz=-0.7, L=Lc, n_jobs=1))
        b = np.asarray(op.osc_prob_3nu_earth(E, costhz=-0.7, L=Lc, n_jobs=2))
    assert np.max(np.abs(a - b)) < 1e-3
    assert 'not bit for bit' in ' '.join(op.osc_prob_energy_baseline.__doc__.split())
