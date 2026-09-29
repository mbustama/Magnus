import time, warnings, numpy as np
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.oscprob as op, magnus.adiabatic as ad
osc = gd.load_nufit_params('NuFIT 6.1')
R = gd.SUN_RADIUS*gd.UNIT_KM
real_once = ad._hybrid_propagator_once
real_leo = ad._local_evolution_operator
MODE = {'fix': False, 'tol': None, 'log': []}
def leo(*a, **k):
    if MODE['fix'] and 'patch_atol' not in k:
        k['patch_atol'] = min(1e-7, MODE['tol']/10)
    U, ok = real_leo(*a, **k); MODE['log'].append(ok); return U, ok
ad._local_evolution_operator = leo
cases = [('2nu sun 5 MeV', lambda **k: op.osc_prob_2nu_sun(5e6, R, 0.0, sth=osc['s12'], Dm2=osc['D21'], **k)),
         ('3nu sun 5 MeV', lambda **k: op.osc_prob_3nu_sun(5e6, R, 0.0, **osc, **k)),
         ('3nu sun 15 MeV B16', lambda **k: op.osc_prob_3nu_sun(15e6, R, 0.0, **osc, density_profile='B16-GS98', **k))]
for name, f in cases:
    for tol in (1e-7, 1e-9):
        row = []
        for fix in (False, True):
            MODE.update(fix=fix, tol=2*tol, log=[])
            info = {}
            t = time.perf_counter()
            P = float(np.asarray(f(rtol=tol, atol=tol, strategy='hybrid', strategy_info=info)).ravel()[0])
            dt = time.perf_counter() - t
            row.append((P, info.get('certified'), dt, len(MODE['log']), all(MODE['log'])))
        (p0, c0, t0, n0, ok0), (p1, c1, t1, n1, ok1) = row
        print('%-20s tol=%g | today: certified=%s %.2fs patches=%d ok=%s | fixed: certified=%s %.2fs patches=%d ok=%s | |dP|=%.1e'
              % (name, tol, c0, t0, n0, ok0, c1, t1, n1, ok1, abs(p1 - p0)), flush=True)
