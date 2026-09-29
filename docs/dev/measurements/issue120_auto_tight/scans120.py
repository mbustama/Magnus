"""#120: batched scans the rule routes beyond Listing 1 -- exp250-3nu and exp25-3nu over their DOP853 energies."""
import sys, json, time, warnings, numpy as np
sys.argv = ['x', 'unused']
sys.path.insert(0, '/home/mbustamante/Research/magnus/resources/handover_auto_tight'); sys.path.insert(0, '/home/mbustamante/Research/magnus/resources/handover_auto_tight/hp')
import tight_measure as tm, magnus.oscprob as op, magnus.globaldefs as gd
print('magnus from', op.__file__)
R = tm.refs(); osc = gd.load_nufit_params('NuFIT 6.1')
HONEST = tm.HONEST
for name, L, h in (('exp250-3nu', 250, 100), ('exp25-3nu', 25, 10)):
    Eg = np.array([0.002, 0.02, 0.2])
    f = lambda E, **k: op.osc_prob_3nu_matter_exp_density(E, L=L*gd.UNIT_KM, L0=0.0, rho_central=3e3, l_scale=h*gd.UNIT_KM,
            density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **osc, **k)
    for p in (4, 8):
        for rt in (1e-7, 1e-9, 1e-12):
            kw = dict(rtol=rt, atol=rt*1e-2, magnus_exp_order=p)
            res = {}
            for s in ('auto', 'hybrid'):
                info = {}
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    P = np.asarray(f(Eg*gd.UNIT_GEV, strategy=s, strategy_info=info, **kw), float).ravel()
                warned = sorted({x.category.__name__ for x in w} & HONEST)
                errs = []
                for Pi, E in zip(P, Eg):
                    r12, r13, _ = R[(name, E)]
                    errs.append(abs(Pi - r13)/max(rt*1e-2 + rt*abs(r13), abs(r12 - r13)))
                best = np.inf
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    for _ in range(5):
                        t0 = time.perf_counter(); f(Eg*gd.UNIT_GEV, strategy=s, **kw); best = min(best, time.perf_counter() - t0)
                res[s] = (info.get('engine'), max(errs), warned, best)
            a, hy = res['auto'], res['hybrid']
            print('%-10s p=%d rtol %.0e | rule auto -> %-9s err/lim %.2f %-28s %7.1f ms | hybrid: err/lim %.2f %-28s %7.1f ms | ratio %.2f'
                  % (name, p, rt, a[0], a[1], a[2] or '-', 1e3*a[3], hy[1], hy[2] or '-', 1e3*hy[3], a[3]/hy[3]))
