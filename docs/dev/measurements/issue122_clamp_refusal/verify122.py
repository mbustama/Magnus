"""#122: probabilities, warnings and times of a set of batched scans, for one source tree."""
import sys, time, json, warnings, numpy as np
import magnus.oscprob as op, magnus.globaldefs as gd
out = sys.argv[1]
osc = gd.load_nufit_params('NuFIT 6.1')
ster = {2: dict(sth=osc['s12'], Dm2=osc['D21']), 3: dict(osc), 4: dict(osc, s14=np.sqrt(0.1), s24=np.sqrt(0.1), D41=1.0),
        5: dict(osc, s14=np.sqrt(0.1), s24=np.sqrt(0.1), s15=np.sqrt(0.06), s25=np.sqrt(0.06), D41=1.0, D51=1.7)}
fns = {2: op.osc_prob_2nu_matter_exp_density, 3: op.osc_prob_3nu_matter_exp_density, 4: op.osc_prob_4nu_matter_exp_density, 5: op.osc_prob_5nu_matter_exp_density}
win = {2: (0.0005, 0.05), 3: (0.002, 0.2), 4: (2.0, 20.0), 5: (2.0, 20.0)}
jobs = []
for d in (2, 3, 4, 5):                                    # Listing 1's curves, 26 energies each
    E = np.logspace(*np.log10(win[d]), 26)*gd.UNIT_GEV
    base = dict(L=25*gd.UNIT_KM, L0=0.0, rho_central=3e3, l_scale=10*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **ster[d])
    for rtol in (1e-3, 1e-6, 1e-9, 1e-12, 5e-13):
        for floor in ((None, 13, 19) if rtol < 1e-11 else (None,)):
            extra = {} if floor is None else dict(n_slabs=floor)
            jobs.append(('L1-%dnu rtol %g floor %s' % (d, rtol, floor), fns[d], E, dict(base, strategy='magnus', magnus_exp_order=4, rtol=rtol, atol=rtol*1e-2, **extra)))
for L, h in ((250, 100), (2500, 1000)):                   # longer exponential profiles, 3nu
    E = np.logspace(np.log10(0.002), np.log10(0.2), 12)*gd.UNIT_GEV
    for rtol in (1e-3, 1e-9):
        jobs.append(('exp%d-3nu rtol %g' % (L, rtol), fns[3], E, dict(L=L*gd.UNIT_KM, L0=0.0, rho_central=3e3, l_scale=h*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **ster[3], strategy='magnus', rtol=rtol, atol=rtol*1e-2)))
E = np.logspace(np.log10(2.0), np.log10(20.0), 26)*gd.UNIT_GEV       # a breakpoint grid, 5nu
jobs.append(('bp-5nu rtol 5e-13 floor 13', fns[5], E, dict(L=25*gd.UNIT_KM, L0=0.0, rho_central=3e3, l_scale=10*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **ster[5], strategy='magnus', magnus_exp_order=4, rtol=5e-13, atol=5e-15, n_slabs=13, t_breakpoints=np.array([12.5*gd.UNIT_KM]))))
E = np.linspace(5e-3, 15e-3, 12)*gd.UNIT_GEV                          # a solar scan, 2nu, 0.2 R_sun
jobs.append(('sun2-0.2R rtol 1e-9', lambda E, **k: op.osc_prob_2nu_sun(E, 0.2*gd.SUN_RADIUS*gd.UNIT_KM, 0.0, sth=osc['s12'], Dm2=osc['D21'], nu_i=gd.NUE, nu_f=gd.NUE, **k), E, dict(strategy='magnus', rtol=1e-9, atol=1e-11)))
res = {}
for name, f, E, kw in jobs:
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = np.asarray(f(E, strategy_info=info, **kw), float)
    best = np.inf
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(3):
            t = time.perf_counter(); f(E, **kw); best = min(best, time.perf_counter() - t)
    res[name] = dict(P=P.tolist(), t=best, engine=info.get('engine'), warnings=sorted({x.category.__name__ for x in w}))
json.dump(res, open(out, 'w'))
print('magnus from', op.__file__[-45:], '| scans:', len(res))
