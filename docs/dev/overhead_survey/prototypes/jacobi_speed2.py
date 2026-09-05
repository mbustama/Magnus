import sys
import numpy as np
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
from jacobi_proto2 import _jacobi_expm_core2
import gen_profile_benchmarks as gpb
timed = gpb.timed
c0 = timed(gpb.control)

rng = np.random.default_rng(1)
for d in (3, 4, 5):
    nB = 1104
    X = rng.normal(size=(nB, d, d)) + 1j*rng.normal(size=(nB, d, d))
    K = (X + np.conj(np.swapaxes(X, -1, -2)))/2
    K *= 0.5/np.max(np.abs(np.linalg.eigvalsh(K)))
    out = np.empty_like(K); lam = np.empty((nB, d))
    sw = _jacobi_expm_core2(K, out, lam)
    lam_e, V = np.linalg.eigh(K)
    U_e = (V*np.exp(-1j*lam_e)[..., None, :]) @ np.conj(np.swapaxes(V, -1, -2))
    err = np.max(np.abs(out - U_e))
    def t_jacobi():
        return _jacobi_expm_core2(K, out, lam)
    def t_eigh_path():
        lam_e, V = np.linalg.eigh(K)
        Vh = np.conj(np.swapaxes(V, -1, -2))
        return (V*np.exp(-1j*lam_e)[..., None, :]) @ Vh
    tj = timed(t_jacobi); te = timed(t_eigh_path)
    print("d=%d  jacobi %7.1f us (%6.1f ns/mat, %.1f sweeps)  eigh-path %7.1f us  ratio %5.2fx  |dU|=%.1e"
          % (d, 1e6*tj, 1e9*tj/nB, sw, 1e6*te, te/tj, err))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
