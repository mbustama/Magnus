import sys, time
import numpy as np
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
from jacobi_proto import _jacobi_expm_core
import gen_profile_benchmarks as gpb
timed = gpb.timed
c0 = timed(gpb.control)

rng = np.random.default_rng(1)
for d in (2, 3, 4, 5):
    nB = 1104
    X = rng.normal(size=(nB, d, d)) + 1j*rng.normal(size=(nB, d, d))
    K = (X + np.conj(np.swapaxes(X, -1, -2)))/2
    # scale to a realistic Magnus slab norm (|| Omega || < pi)
    K *= 0.5/np.max(np.abs(np.linalg.eigvalsh(K)))
    out = np.empty_like(K); lam = np.empty((nB, d))
    _jacobi_expm_core(K, out, lam)

    def t_jacobi():
        _jacobi_expm_core(K, out, lam)
        return out
    def t_eigh_path():
        lam_e, V = np.linalg.eigh(K)
        Vh = np.conj(np.swapaxes(V, -1, -2))
        return (V*np.exp(-1j*lam_e)[..., None, :]) @ Vh
    tj = timed(t_jacobi); te = timed(t_eigh_path)
    print("d=%d nB=%d  jacobi %8.1f us (%6.1f ns/mat)   eigh-path %8.1f us (%7.1f ns/mat)  ratio %5.2fx"
          % (d, nB, 1e6*tj, 1e9*tj/nB, 1e6*te, 1e9*te/nB, te/tj))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
