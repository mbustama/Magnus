import sys, time, warnings, numpy as np
warnings.simplefilter('ignore')
exec(open('check_l1.py').read().split("mode = sys.argv")[0])
ref = np.load('l1_ref.npz')
tols = {'1e-12': dict(rtol=1e-12, atol=1e-14), '1e-3': dict(rtol=1e-3, atol=1e-3)}
for tlabel, strat in (('1e-12', 'auto'), ('1e-12', 'magnus'), ('1e-3', 'auto'), ('1e-3', 'hybrid'), ('1e-3', 'magnus')):
    tol = tols[tlabel]
    for k, (f, E, p) in curves.items():
        f(E[:2], **base, **p, **tol, strategy=strat)
    T, eng, err = {}, set(), 0.0
    for k, (f, E, p) in curves.items():
        info = {}
        t = time.perf_counter(); P = np.asarray(f(E, **base, **p, **tol, strategy=strat, strategy_info=info)); T[k] = time.perf_counter() - t
        eng.add(info.get('engine')); err = max(err, float(np.max(np.abs(P - ref[k]))))
    print('tol %-6s strategy=%-7s engines=%-22s total=%7.3f s  max err vs 1e-12 hybrid %.1e   (%s)' % (tlabel, strat, ','.join(sorted(eng)), sum(T.values()), err, ' '.join('%s=%.3f' % kv for kv in T.items())), flush=True)
