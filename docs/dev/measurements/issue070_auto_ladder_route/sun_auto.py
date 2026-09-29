import time, warnings, numpy as np
import magnus.globaldefs as gd, magnus.oscprob as op
osc = gd.load_nufit_params('NuFIT 6.1'); R = gd.SUN_RADIUS*gd.UNIT_KM
cases = [('2nu sun 5 MeV', lambda **k: op.osc_prob_2nu_sun(5e6, R, 0.0, sth=osc['s12'], Dm2=osc['D21'], **k)),
         ('3nu sun 5 MeV', lambda **k: op.osc_prob_3nu_sun(5e6, R, 0.0, **osc, **k))]
for name, f in cases:
    for tol in (1e-7, 1e-9):
        info = {}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            t = time.perf_counter(); f(rtol=tol, atol=tol, strategy_info=info); dt = time.perf_counter() - t
        print('%-14s tol=%g auto -> engine=%s declined=%s %.1fs warnings=%s' % (name, tol, info.get('engine'),
              [d[0] + ': ' + d[1][:40] for d in info.get('declined', [])], dt, sorted({x.category.__name__ for x in w})), flush=True)
