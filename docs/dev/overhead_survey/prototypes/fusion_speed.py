import sys
import numpy as np
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import magnus.magnus as mg
import gen_profile_benchmarks as gpb
from fusion_proto import _antiherm_scale_dev, _gl4_omega
timed = gpb.timed
c0 = timed(gpb.control)

rng = np.random.default_rng(2)
for d in (2, 3, 4, 5):
    # realistic chunked shape: (chunk*1104, d, d)
    chunk = {2: 7, 3: 3, 4: 1, 5: 1}[d]
    nB = chunk*1104
    X = rng.normal(size=(nB, d, d)) + 1j*rng.normal(size=(nB, d, d))
    H = (X + np.conj(np.swapaxes(X, -1, -2)))/2
    Om = np.ascontiguousarray(-1j*H*0.1)
    A1 = np.ascontiguousarray(-1j*H*0.3)
    A2 = np.ascontiguousarray(np.roll(A1, 1, axis=0))
    h = np.abs(rng.normal(size=nB))*0.01
    out = np.empty_like(A1)
    _antiherm_scale_dev(Om); _gl4_omega(A1, A2, h, out)

    def t_frame_old():
        K = 1j*Om
        Kh = np.conj(np.swapaxes(K, -1, -2))
        scale = np.max(np.abs(K))
        return np.max(np.abs(K - Kh)) <= 1e-12*scale
    def t_frame_new():
        scale, dev = _antiherm_scale_dev(Om)
        return dev <= 1e-12*scale
    def t_gl4_old():
        C = mg._commutator_batched(A2, A1)
        hh = h[:, None, None]
        return 0.5*hh*(A1 + A2) + (np.sqrt(3.0)/12.0)*hh*hh*C
    def t_gl4_new():
        return _gl4_omega(A1, A2, h, out)
    # equivalence
    K = 1j*Om; Kh = np.conj(np.swapaxes(K, -1, -2))
    s_old, d_old = np.max(np.abs(K)), np.max(np.abs(K - Kh))
    s_new, d_new = _antiherm_scale_dev(Om)
    om_old = t_gl4_old(); om_new = t_gl4_new()
    fo, fn = timed(t_frame_old), timed(t_frame_new)
    go, gn = timed(t_gl4_old), timed(t_gl4_new)
    print("d=%d nB=%5d  frame %7.1f -> %6.1f us (%4.1fx; ds=%.1e dd=%.1e)  gl4 %7.1f -> %6.1f us (%4.1fx; dOm=%.2e)"
          % (d, nB, 1e6*fo, 1e6*fn, fo/fn, abs(s_old-s_new), abs(d_old-d_new),
             1e6*go, 1e6*gn, go/gn, np.max(np.abs(om_old-om_new))))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
