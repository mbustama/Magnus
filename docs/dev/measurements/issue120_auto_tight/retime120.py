"""#120: the four routed cases ever slower than the hybrid in B1, re-timed alone and interleaved.
Rule's auto (ladder, /1) against strategy='hybrid' (what main's auto runs below 1e-6)."""
import sys, time, warnings, numpy as np
sys.path.insert(0, '/home/mbustamante/Research/magnus/resources/handover_auto_tight/hp')
warnings.simplefilter('ignore')
sys.argv = ['x', 'unused']
import verify_b as vb, magnus.oscprob as op, magnus.globaldefs as gd
print('magnus from', op.__file__)
want = [('L1-3nu', 0.1663527542205342, 1e-12, 4), ('L1-2nu', 0.0045600541967795475, 1e-9, 4),
        ('L1-2nu', 0.0045600541967795475, 1e-12, 4), ('L1-3nu', 0.13836619418378734, 1e-9, 4)]
cases = {(n, E): s for n, E, s in vb.cases()}
for name, E_gev, tol, p in want:
    f, _ = vb.call_for(cases[(name, E_gev)]); E = E_gev*gd.UNIT_GEV
    kw = dict(rtol=tol, atol=tol*1e-2, magnus_exp_order=p)
    eng = {}
    for s in ('auto', 'hybrid'):
        info = {}; f(E, strategy=s, strategy_info=info, **kw); eng[s] = info.get('engine')
    ratios, best = [], {'auto': np.inf, 'hybrid': np.inf}
    for rnd in range(15):
        t = {}
        for s in (('auto', 'hybrid') if rnd % 2 == 0 else ('hybrid', 'auto')):
            b = np.inf
            for _ in range(3):
                t0 = time.perf_counter(); f(E, strategy=s, **kw); b = min(b, time.perf_counter() - t0)
            t[s] = b; best[s] = min(best[s], b)
        ratios.append(t['auto']/t['hybrid'])
    print('%-7s E=%.4g GeV rtol %.0e p=%d: auto->%s %.2f ms, hybrid->%s %.2f ms; ratio median %.2f, of bests %.2f, max %.2f'
          % (name, E_gev, tol, p, eng['auto'], 1e3*best['auto'], eng['hybrid'], 1e3*best['hybrid'],
             np.median(ratios), best['auto']/best['hybrid'], max(ratios)))
