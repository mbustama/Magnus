"""engines.rst's engine table, tight-tolerance column: the engine strategy_info reports, per request shape."""
import time, warnings, numpy as np
import magnus.oscprob as op, magnus.globaldefs as gd
osc = gd.load_nufit_params('NuFIT 6.1')
def exp3(E, L, h, **kw):
    return op.osc_prob_3nu_matter_exp_density(np.asarray(E, float)*gd.UNIT_GEV, L=np.asarray(L, float)*gd.UNIT_KM, L0=0.0, rho_central=3e3,
        l_scale=h*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **osc, **kw)
cases = [
    ('one point, small phase (25 km, 2 MeV)',        lambda **k: exp3([0.002], 25, 10, **k)),
    ('one point, large phase (2500 km, 2 MeV)',      lambda **k: exp3([0.002], 2500, 1000, **k)),
    ('12 energies, small phase (25 km)',             lambda **k: exp3(np.logspace(-2.7, -0.7, 12), 25, 10, **k)),
    ('12 energies, large phase (2500 km)',           lambda **k: exp3(np.logspace(-2.7, -0.7, 12), 2500, 1000, **k)),
    ('300 energies 3-100 MeV at 200 km (h 100 km)',  lambda **k: exp3(np.logspace(np.log10(3e-3), np.log10(0.1), 300), 200, 100, **k)),
    ('3 baselines, small phase (5-25 km, 20 MeV)',   lambda **k: exp3([0.02]*3, [5, 10, 25], 10, **k)),
    ('3 baselines, large phase (500-2500 km, 2 MeV)', lambda **k: exp3([0.002]*3, [500, 1000, 2500], 1000, **k)),
    ('Sun, 2nu 10 MeV 0.9 R_sun',                    lambda **k: op.osc_prob_2nu_sun(np.array([0.01])*gd.UNIT_GEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM, 0.0, sth=osc['s12'], Dm2=osc['D21'], nu_i=gd.NUE, nu_f=gd.NUE, **k)),
]
print('from', op.__file__[-45:])
for tol in (1e-3, 1e-8):
    for name, f in cases:
        info = {}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            t = time.perf_counter(); f(rtol=tol, atol=tol, strategy_info=info); dt = time.perf_counter() - t
        note = [x for x in info.get('trace', []) if 'estimated_phase' in x]
        ph = ('%.0f' % note[0]['estimated_phase']) if note else '-'
        print('%.0e  %-46s engine %-10s phase %-6s %6.2f s  %s' % (tol, name, info.get('engine'), ph, dt,
              sorted({x.category.__name__ for x in w} & {'ToleranceNotAchievedWarning', 'HybridCertificationWarning'})))
