import sys, time, warnings, numpy as np
exec(open('check_l1.py').read().split("mode = sys.argv")[0])
import magnus.magnus as mg
tol = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
ref = np.load('l1_ref.npz')
orig = mg.evolution_operators_from_samples
log = []
def spy(At, widths, *a, **k):
    with mg._deferred_slab_norm() as sink:
        s0 = len(sink); U = orig(At, widths, *a, **k)
        log.append((len(widths), At.shape[0], max(sink[s0:], default=0.0)))
    return U
for k, (f, E, p) in curves.items():
    mg.evolution_operators_from_samples = spy
    log.clear(); info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = np.asarray(f(E, **base, **p, rtol=tol, atol=tol, strategy_info=info))
    d = {}
    for ns, ne, nm in log:
        a = d.setdefault(ns, [0, 0.0]); a[0] += ne; a[1] = max(a[1], nm)
    mg.evolution_operators_from_samples = orig
    f(E, **base, **p, rtol=tol, atol=tol)
    t = time.perf_counter(); f(E, **base, **p, rtol=tol, atol=tol); dt = time.perf_counter() - t
    err = float(np.max(np.abs(P - ref[k])))
    print('%s t=%.3fs err=%.1e warns=%d route=%s | %s' % (k, dt, err, len(w), info.get('route'),
          ' '.join('n=%d:E%d:%.2f' % (ns, v[0], v[1]) for ns, v in d.items())))
