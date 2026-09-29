"""Evidence for the root-cause claim: Listing 1's worst point (3nu, 2 MeV, 25 km exponential
profile), hybrid at rtol=1e-12/atol=1e-14, order 4, with only patch_atol changed.
DOP853 at rtol=1e-13/atol=1e-15 gives P_ee = 0.702322921066483 (see worst_dop.py)."""
import time, warnings, numpy as np
warnings.simplefilter('ignore')
exec(open('check_l1.py').read().split("mode = sys.argv")[0])
import magnus.oscprob as op, magnus.adiabatic as ad
CAP = {}
orig = op._hybrid_propagator_scan
def spy(H, *a, **k):
    CAP['H'] = H; return orig(H, *a, **k)
op._hybrid_propagator_scan = spy
f, E, p = curves['3nu']
P_wrapper = float(np.asarray(f(E[:1], **base, **p, rtol=1e-12, atol=1e-14, strategy='hybrid')).ravel()[0])
Hf = CAP['H'](E[0]); ref = 0.702322921066483
print('through the wrapper (strategy=hybrid, order 4): err=%.1e' % abs(P_wrapper - ref))
real = ad._local_evolution_operator
for order in (4, 6, 8):
    for patch_atol in (None, 1e-9, 1e-11, 1e-13):
        seen = []
        def wrapped(*a, **k):
            if patch_atol is not None:
                k['patch_atol'] = patch_atol
            U, ok = real(*a, **k); seen.append(ok); return U, ok
        ad._local_evolution_operator = wrapped
        info = {}; t = time.perf_counter()
        U, w, cert = ad.hybrid_propagator(Hf, 0.0, base['L'], rtol=1e-12, atol=1e-14,
                                          magnus_exp_order=order, info=info)
        P = (np.abs(U)**2).T[0, 0]
        print('order=%d patch_atol=%-8s err=%.1e certified=%s patch calls ok=%s windows=%s time=%.3fs'
              % (order, 'default' if patch_atol is None else patch_atol, abs(P - ref), cert, seen,
                 [(round(a/base['L'], 3), round(b/base['L'], 3)) for a, b in w], time.perf_counter() - t))
ad._local_evolution_operator = real
