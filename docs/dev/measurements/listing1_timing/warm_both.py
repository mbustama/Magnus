import sys, time, warnings, numpy as np
warnings.simplefilter('ignore')
exec(open('check_l1.py').read().split("mode = sys.argv")[0])
for label, tol in (('tight', dict(rtol=1e-12, atol=1e-14)), ('default', dict(rtol=1e-3, atol=1e-3))):
    for k, (f, E, p) in curves.items():
        f(E[:2], **base, **p, **tol)            # warm
    T = {}
    for k, (f, E, p) in curves.items():
        t = time.perf_counter(); f(E, **base, **p, **tol); T[k] = time.perf_counter() - t
    print('%-8s' % label, ' '.join('%s=%.3f' % kv for kv in T.items()), 'total=%.3f' % sum(T.values()))
