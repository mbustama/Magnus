"""Broad A/B benchmark of the public API: outputs (bitwise), warning categories, timings."""
import sys, time, pickle, warnings, numpy as np
src = sys.argv[1]; out = sys.argv[2]; mode = sys.argv[3]   # mode: 'values' or 'timing'
sys.path.insert(0, src)
import magnus.oscprob as o, magnus.globaldefs as gd, magnus.hamiltonians as hm, magnus.earth as earth
KM = gd.UNIT_KM; GEV = gd.UNIT_GEV; MEV = gd.UNIT_MEV
p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
std = {k: p[k] for k in ('s12','s23','s13','dCP','D21','D31')}
hv3 = np.asarray(hm.hamiltonian_3nu_vacuum_energy_independent(**std)); e00 = np.diag([1., 0, 0])
Vf = lambda rho: np.sqrt(2)*gd.GF*rho*gd.UNIT_G_PER_CM3/((gd.MASS_PROTON + gd.MASS_NEUTRON)/2)*0.5
R = gd.SUN_RADIUS*KM
Es = np.geomspace(0.5, 20, 12)*GEV; Esun = np.geomspace(1, 15, 8)*MEV
Lc = lambda cz: earth.distance_traveled_inside_earth(cz)*KM
C = {}
TWO = dict(sth=0.55, Dm2=2.5e-3)
def case(name, f, fl=None):
    if fl == 2:
        g = f
        f = lambda g=g: g(**TWO)
    C[name] = f
for fl, fv in ((2, o.osc_prob_2nu_vacuum), (3, o.osc_prob_3nu_vacuum), (4, o.osc_prob_4nu_vacuum), (5, o.osc_prob_5nu_vacuum)):
    case('vac%d single' % fl, lambda fv=fv, **kw: fv(1*GEV, 1300*KM, **kw), fl)
    case('vac%d scan' % fl, lambda fv=fv, **kw: fv(Es, 1300*KM, **kw), fl)
    case('vac%d avg' % fl, lambda fv=fv, **kw: fv(Es, 1e5*KM, average=True, **kw), fl)
for fl, fc in ((2, o.osc_prob_2nu_matter_constant_density), (3, o.osc_prob_3nu_matter_constant_density), (4, o.osc_prob_4nu_matter_constant_density), (5, o.osc_prob_5nu_matter_constant_density)):
    case('const%d single' % fl, lambda fc=fc, **kw: fc(3*GEV, 2000*KM, rho=3.0, **kw), fl)
    case('const%d scan' % fl, lambda fc=fc, **kw: fc(Es, 2000*KM, rho=3.0, **kw), fl)
for fl, fe in ((2, o.osc_prob_2nu_matter_exp_density), (3, o.osc_prob_3nu_matter_exp_density), (4, o.osc_prob_4nu_matter_exp_density), (5, o.osc_prob_5nu_matter_exp_density)):
    case('exp%d single' % fl, lambda fe=fe, **kw: fe(3*GEV, 3000*KM, 0.0, 10.0, 1000*KM, **kw), fl)
    case('exp%d escan' % fl, lambda fe=fe, **kw: fe(Es, 3000*KM, 0.0, 10.0, 1000*KM, **kw), fl)
    case('exp%d lscan' % fl, lambda fe=fe, **kw: fe(3*GEV, np.linspace(500, 5000, 10)*KM, 0.0, 10.0, 1000*KM, **kw), fl)
case('exp3 tight', lambda: o.osc_prob_3nu_matter_exp_density(3*GEV, 3000*KM, 0.0, 10.0, 1000*KM, rtol=1e-8, atol=1e-8))
case('exp3 magnus', lambda: o.osc_prob_3nu_matter_exp_density(3*GEV, 3000*KM, 0.0, 10.0, 1000*KM, strategy='magnus'))
case('exp3 hybrid', lambda: o.osc_prob_3nu_matter_exp_density(3*GEV, 3000*KM, 0.0, 10.0, 1000*KM, strategy='hybrid'))
case('exp3 avg', lambda: o.osc_prob_3nu_matter_exp_density(Es, 3000*KM, 0.0, 10.0, 1000*KM, average=True))
case('exp3 operator', lambda: o.osc_prob_3nu_matter_exp_density(3*GEV, 3000*KM, 0.0, 10.0, 1000*KM, return_evolution_operator=True)[0])
for fl, fe in ((2, o.osc_prob_2nu_earth), (3, o.osc_prob_3nu_earth), (4, o.osc_prob_4nu_earth), (5, o.osc_prob_5nu_earth)):
    for cz in (-1.0, -0.6, -0.2):
        case('earth%d cz%.1f' % (fl, cz), lambda fe=fe, cz=cz, **kw: fe(3*GEV, costhz=cz, L=Lc(cz), **kw), fl)
    case('earth%d escan' % fl, lambda fe=fe, **kw: fe(Es, costhz=-0.7, L=Lc(-0.7), **kw), fl)
case('earth3 tight', lambda: o.osc_prob_3nu_earth(3*GEV, costhz=-0.8, L=Lc(-0.8), rtol=1e-8, atol=1e-8))
case('earth3 magnus', lambda: o.osc_prob_3nu_earth(3*GEV, costhz=-0.8, L=Lc(-0.8), strategy='magnus'))
case('earth3 avg', lambda: o.osc_prob_3nu_earth(Es, costhz=-0.8, L=Lc(-0.8), average=True))
case('earth3 nubar', lambda: o.osc_prob_3nu_earth(Es, costhz=-0.8, L=Lc(-0.8), nubar=True))
case('earth3 nsi', lambda: o.osc_prob_3nu_earth_nsi(3*GEV, costhz=-0.8, L=Lc(-0.8), eps_ee=0.2, eps_mt=0.05))
case('earth3 liv', lambda: o.osc_prob_3nu_earth_liv(3*GEV, costhz=-0.8, L=Lc(-0.8), b1=1e-23))
case('earth3 nsi escan', lambda: o.osc_prob_3nu_earth_nsi(Es, costhz=-0.8, L=Lc(-0.8), eps_ee=0.2))
for fl, fs in ((2, o.osc_prob_2nu_sun), (3, o.osc_prob_3nu_sun), (4, o.osc_prob_4nu_sun), (5, o.osc_prob_5nu_sun)):
    case('sun%d single' % fl, lambda fs=fs, **kw: fs(8*MEV, R, 0.0, **kw), fl)
    if fl != 5:  # the 5-flavor averaged Sun at defaults is issue #148 (hangs on main)
        case('sun%d avg' % fl, lambda fs=fs, **kw: fs(Esun, R, 0.0, average=True, **kw), fl)
case('sun3 escan', lambda: o.osc_prob_3nu_sun(Esun, R, 0.0))
case('sun3 BS05', lambda: o.osc_prob_3nu_sun(8*MEV, R, 0.0, density_profile='BS05-AGS-OP'))
case('sun3 B16 avg', lambda: o.osc_prob_3nu_sun(Esun, R, 0.0, density_profile='B16-GS98', average=True))
case('sun3 magnus', lambda: o.osc_prob_3nu_sun(8*MEV, 0.5*R, 0.0, strategy='magnus'))
case('sun3 hybrid', lambda: o.osc_prob_3nu_sun(8*MEV, R, 0.0, strategy='hybrid'))
case('sun3 nsi', lambda: o.osc_prob_3nu_sun_nsi(8*MEV, R, 0.0, eps_ee=0.1))
profiles = {
    'linear': (lambda x: 2 + 8*np.asarray(x, float)/(5000*KM), 5000*KM),
    'sine': (lambda x: 5 + 3*np.sin(np.asarray(x, float)/(800*KM)), 6000*KM),
    'bump200': (lambda x: 4 + 6*np.exp(-((np.asarray(x, float) - 3000*KM)/(200*KM))**2), 6000*KM),
    'twoexp': (lambda x: 3 + 9*np.exp(-np.abs(np.asarray(x, float) - 4000*KM)/(500*KM)), 8000*KM),
    'poly': (lambda x: 3 + 4*(np.asarray(x, float)/(7000*KM))**2 - 2*(np.asarray(x, float)/(7000*KM))**3, 7000*KM),
    'step': (lambda x: np.where(np.asarray(x, float) < 2000*KM, 3.0, 8.0), 5000*KM),
    'spike5': (lambda x: 3 + 50*np.exp(-((np.asarray(x, float) - 1000*KM)/(5*KM))**2), 4000*KM),
    'spike50': (lambda x: 3 + 50*np.exp(-((np.asarray(x, float) - 1000*KM)/(50*KM))**2), 4000*KM),
}
for nm, (rho, L) in profiles.items():
    case('std %s single' % nm, lambda rho=rho, L=L: o.osc_prob_matter_std_potential(3, rho, 3*GEV, L, std, density_matter_is_in_g_per_cm3=True))
    case('std %s escan' % nm, lambda rho=rho, L=L: o.osc_prob_matter_std_potential(3, rho, Es[:6], L, std, density_matter_is_in_g_per_cm3=True))
    case('std %s lscan' % nm, lambda rho=rho, L=L: o.osc_prob_matter_std_potential(3, rho, 3*GEV, np.linspace(0.2, 1, 6)*L, std, density_matter_is_in_g_per_cm3=True))
    case('std %s magnus' % nm, lambda rho=rho, L=L: o.osc_prob_matter_std_potential(3, rho, 3*GEV, L, std, density_matter_is_in_g_per_cm3=True, strategy='magnus'))
case('std sine avg', lambda: o.osc_prob_matter_std_potential(3, profiles['sine'][0], Es[:6], 6000*KM, std, density_matter_is_in_g_per_cm3=True, average=True))
case('std sine 5nu', lambda: o.osc_prob_matter_std_potential(5, profiles['sine'][0], 3*GEV, 6000*KM, dict(std, s14=0.1, d14=0.0, s15=0.05, d15=0.0, s24=0.1, d24=0.0, s25=0.05, s34=0.1, s35=0.05, d35=0.0, D41=1.0, D51=2.0), density_matter_is_in_g_per_cm3=True))
case('nsi sine', lambda: o.osc_prob_matter_nsi(3, profiles['sine'][0], 3*GEV, 6000*KM, std, dict(eps_ee=0.1, eps_em=0.0, eps_et=0.0, eps_mm=0.0, eps_mt=0.02, eps_tt=0.0), density_matter_is_in_g_per_cm3=True))
rawH = {
    'smooth': lambda e, l: hv3/e + Vf(4 + 2*np.sin(l/(700*KM)))*e00,
    'castle': lambda e, l: hv3/e + Vf(3.3 if (l < 3000*KM or l > 9000*KM) else 11.5)*e00,
    'const': lambda e, l: hv3/e + Vf(4.0)*e00,
}
for nm, H in rawH.items():
    case('raw %s single' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, 6200*KM))
    case('raw %s escan' % nm, lambda H=H: o.osc_prob_energy_baseline(H, Es[:5], 6200*KM))
    case('raw %s lscan' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, np.linspace(1000, 11000, 6)*KM))
    case('raw %s tight' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, 6200*KM, rtol=1e-7, atol=1e-7))
    case('raw %s nocumul' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, np.linspace(1000, 11000, 4)*KM, cumulative=False))
    case('raw %s bp' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, 6200*KM, t_breakpoints=[3000*KM]))
    case('raw %s op' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, 6200*KM, return_evolution_operator=True)[1])
    case('raw %s strict' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, 6200*KM, strict_convergence=True))
    case('raw %s notol' % nm, lambda H=H: o.osc_prob_energy_baseline(H, 3*GEV, 6200*KM, rtol=None, atol=None, n_slabs=50))
    case('osc_prob %s' % nm, lambda H=H: o.osc_prob(lambda l, H=H: H(3*GEV, l), 0.0, 6200*KM, rtol=1e-3, atol=1e-3))
case('pd vac avg', lambda: o.osc_prob_pseudo_dirac_vacuum(100e12, 1e30, {0: 1e-17}, average=True, average_spread=0.01))
case('cross_check', lambda: np.array([v['P'] if isinstance(v, dict) and 'P' in v else np.nan for v in o.cross_check_strategies(o.osc_prob_3nu_matter_exp_density, 3*GEV, 3000*KM, 0.0, 10.0, 1000*KM).get('results', {}).values()], dtype=object) if False else 0)

def run(f):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        r = f()
    return np.asarray(r), sorted({x.category.__name__ for x in w})

ONLY = sys.argv[5].split('|') if len(sys.argv) > 5 else None
if ONLY: C = {k: v for k, v in C.items() if k in ONLY}
if mode == 'values':
    res = {}
    for k, f in C.items():
        try:
            res[k] = run(f)
        except Exception as e:
            res[k] = ('ERROR', repr(e))
    pickle.dump(res, open(out, 'wb'))
else:
    warnings.simplefilter('ignore')
    reps = int(sys.argv[4])
    t = {}
    for k, f in C.items():
        f()
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter(); f(); ts.append(time.perf_counter() - t0)
        t[k] = min(ts)
    pickle.dump(t, open(out, 'wb'))
print(len(C), 'cases')
