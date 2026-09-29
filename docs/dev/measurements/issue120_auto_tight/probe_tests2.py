import time, warnings, numpy as np
import magnus.oscprob as op, magnus.globaldefs as gd
print('from', op.__file__[-40:])
OSC = dict(s12=0.55, s23=0.68, s13=0.15, dCP=3.7, D21=7.5e-5, D31=2.5e-3)
PROFILE = dict(L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.e3, l_scale=10.0*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE)
ENERGIES = np.logspace(np.log10(0.002), np.log10(0.2), 12)*gd.UNIT_GEV
p3 = lambda E, **kw: np.asarray(op.osc_prob_3nu_matter_exp_density(E, **PROFILE, **OSC, **kw))
def run(label, f):
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        t = time.perf_counter()
        try:
            f(info); err = None
        except Exception as e:
            err = '%s: %s' % (type(e).__name__, str(e)[:70])
        dt = time.perf_counter() - t
    print('%-34s engine %-9s declined %-26s warnings %s %s %.2f s' % (label, info.get('engine'), dict(info.get('declined', [])).get('hybrid'), sorted({x.category.__name__ for x in w}), err or '', dt))
for m in ('simpson', 'trapezoid'):
    run('%s 1e-9' % m, lambda i, m=m: p3(ENERGIES[:2], rtol=1e-9, atol=1e-9, integration_method=m, strategy_info=i))
run('full Sun 2nu 10 MeV 1e-9', lambda i: op.osc_prob_2nu_sun(np.array([10e-3])*gd.UNIT_GEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM, 0.0, sth=0.55, Dm2=7.5e-5, nu_i=gd.NUE, nu_f=gd.NUE, rtol=1e-9, atol=1e-9, strategy_info=i))
for p in (0, 10, 'x'):
    run('order %r 1e-9' % (p,), lambda i, p=p: p3(ENERGIES[:2], rtol=1e-9, atol=1e-11, magnus_exp_order=p, strategy_info=i))
for s in ('hybrid', 'magnus'):
    run('strategy %s 1e-12' % s, lambda i, s=s: p3(ENERGIES[:2], rtol=1e-12, atol=1e-12, strategy=s, strategy_info=i))
run('auto 2 energies 1e-9', lambda i: p3(ENERGIES[:2], rtol=1e-9, atol=1e-9, strategy_info=i))
run('auto 1 energy 1e-12', lambda i: p3(ENERGIES[:1], rtol=1e-12, atol=1e-12, strategy_info=i))
