import sys
import numpy as np
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
from jacobi_proto2 import _jacobi_expm_core2
from jacobi_proto4 import _jacobi_expm_warm_mgs
import gen_profile_benchmarks as gpb
timed = gpb.timed
c0 = timed(gpb.control)

for d in (4, 5):
    Om = np.load('om_d%d.npy' % d)
    K = np.ascontiguousarray(1j*Om)
    nB = K.shape[0]
    out = np.empty_like(K); lam = np.empty((nB, d))
    out_w = np.empty_like(K); lam_w = np.empty((nB, d))
    sw_c = _jacobi_expm_core2(K, out, lam)
    sw_w = _jacobi_expm_warm_mgs(K, out_w, lam_w)
    lam_e, V = np.linalg.eigh(K)
    U_e = (V*np.exp(-1j*lam_e)[..., None, :]) @ np.conj(np.swapaxes(V, -1, -2))
    ec = np.max(np.abs(out - U_e)); ew = np.max(np.abs(out_w - U_e))
    I = np.eye(d)
    uw = np.max(np.abs(np.conj(np.swapaxes(out_w, -1, -2)) @ out_w - I))
    def t_cold(): return _jacobi_expm_core2(K, out, lam)
    def t_warm(): return _jacobi_expm_warm_mgs(K, out_w, lam_w)
    def t_eigh_path():
        lam_e, V = np.linalg.eigh(K)
        Vh = np.conj(np.swapaxes(V, -1, -2))
        return (V*np.exp(-1j*lam_e)[..., None, :]) @ Vh
    tc = timed(t_cold); tw = timed(t_warm); te = timed(t_eigh_path)
    print("d=%d  cold %6.1f ns/mat (%.2f sw, dU %.1e) | warm+mgs %6.1f ns/mat (%.2f sw, dU %.1e, unit %.1e) | eigh %6.1f ns/mat -> warm %4.2fx"
          % (d, 1e9*tc/nB, sw_c, ec, 1e9*tw/nB, sw_w, ew, uw, 1e9*te/nB, te/tw))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
