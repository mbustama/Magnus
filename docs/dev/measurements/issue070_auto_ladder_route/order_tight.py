import sys, time, warnings, numpy as np
warnings.simplefilter('ignore')
exec(open('check_l1.py').read().split("mode = sys.argv")[0])
ref = np.load('l1_ref.npz')
tol = dict(rtol=1e-12, atol=1e-14)
for strat in ('magnus', 'hybrid'):
    for order in (4, 6, 8):
        kw = dict(tol, strategy=strat, magnus_exp_order=order)
        for k, (f, E, p) in curves.items():
            f(E[:2], **base, **p, **kw)
        T, eng, err = {}, set(), 0.0
        for k, (f, E, p) in curves.items():
            info = {}
            t = time.perf_counter(); P = np.asarray(f(E, **base, **p, **kw, strategy_info=info)); T[k] = time.perf_counter() - t
            eng.add(info.get('engine')); err = max(err, float(np.max(np.abs(P - ref[k]))))
        print('strategy=%-7s order=%d engines=%-12s total=%6.3f s  max|dP| vs ref %.1e  (%s)' % (strat, order, ','.join(sorted(eng)), sum(T.values()), err, ' '.join('%s=%.3f' % kv for kv in T.items())), flush=True)
