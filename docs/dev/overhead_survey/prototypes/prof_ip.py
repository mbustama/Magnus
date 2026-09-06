# Profile the interaction-picture core directly: 2-level, exponential profile,
# fixed n_slabs, one pass.
import sys, warnings, cProfile, pstats, io
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
import gen_profile_benchmarks as gpb

warnings.simplefilter('ignore')
KM = gd.CONV_KM_TO_INV_EV if hasattr(gd, 'CONV_KM_TO_INV_EV') else None
l_scale = 3000.0*gpb.gd.UNIT_KM if hasattr(gpb.gd, 'UNIT_KM') else 1.0e12
L = 10.0*l_scale
rng = np.random.default_rng(0)
nE = 64
X = rng.normal(size=(nE, 2, 2)) + 1j*rng.normal(size=(nE, 2, 2))
H_E = ((X + np.conj(np.swapaxes(X, -1, -2)))/2)/l_scale     # phases O(10) over L
h_matt = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
v0 = 0.05/l_scale
VCC = lambda l: v0*np.exp(-np.asarray(l, dtype=float)/l_scale)

n_slabs = int(sys.argv[1]); reps = int(sys.argv[2])
call = lambda: oscprob._osc_prob_ip_exp_core(H_E, l_scale, VCC, h_matt, 0.0, L,
                                             None, None, 2.0, 25, 1, 10**6, n_slabs)
import time
call()
t0 = time.perf_counter(); call(); t1 = time.perf_counter()
print("one pass: %.1f ms at n_slabs=%d nE=%d" % (1e3*(t1-t0), n_slabs, nE))
pr = cProfile.Profile(); pr.enable()
for _ in range(reps):
    call()
pr.disable()
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(16)
print(s.getvalue())
