import sys
import numpy as np
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import magnus.magnus as mg
import gen_profile_benchmarks as gpb
from fusion_proto2 import _antiherm_scale_dev2, _gl4_omega_from_At
timed = gpb.timed
c0 = timed(gpb.control)

rng = np.random.default_rng(2)
for d in (2, 3, 4, 5):
    chunk = {2: 7, 3: 3, 4: 1, 5: 1}[d]
    nE, ns = chunk, 1104
    X = rng.normal(size=(nE, ns, 2, d, d)) + 1j*rng.normal(size=(nE, ns, 2, d, d))
    At = np.ascontiguousarray(-1j*(X + np.conj(np.swapaxes(X, -1, -2)))/2*0.1)
    widths = np.abs(rng.normal(size=ns))*0.01
    h4 = widths[..., None, None]
    Om = mg._magnus_gl(At, widths, 4)
    OmC = np.ascontiguousarray(Om.reshape(-1, d, d))
    nB = OmC.shape[0]
    Atf = np.ascontiguousarray(At.reshape(-1, 2, d, d))
    hf = np.ascontiguousarray(np.broadcast_to(widths, (nE, ns)).reshape(-1))
    out = np.empty((nB, d, d), dtype=complex)
    _antiherm_scale_dev2(OmC); _gl4_omega_from_At(Atf, hf, out)

    # current sequence, exactly as _magnus_gl runs it (strided views + wrappers)
    def t_gl4_old():
        return mg._magnus_gl(At, widths, 4)
    def t_gl4_new():
        _gl4_omega_from_At(Atf, hf, out)
        return out
    def t_frame_old():
        K = 1j*Om
        Kh = np.conj(np.swapaxes(K, -1, -2))
        scale = np.max(np.abs(K))
        return np.max(np.abs(K - Kh)) <= 1e-12*scale
    def t_frame_new():
        scale, dev = _antiherm_scale_dev2(OmC)
        return dev <= 1e-12*scale
    om_old = t_gl4_old(); om_new = t_gl4_new()
    K = 1j*Om; Kh = np.conj(np.swapaxes(K, -1, -2))
    s_old, d_old = np.max(np.abs(K)), np.max(np.abs(K - Kh))
    s_new, d_new = _antiherm_scale_dev2(OmC)
    fo, fn = timed(t_frame_old), timed(t_frame_new)
    go, gn = timed(t_gl4_old), timed(t_gl4_new)
    print("d=%d  frame %7.1f -> %6.1f us (%4.1fx; ds=%.1e dd=%.1e)  gl4(real path) %7.1f -> %6.1f us (%4.1fx; dOm=%.2e)"
          % (d, 1e6*fo, 1e6*fn, fo/fn, abs(s_old-s_new), abs(d_old-d_new),
             1e6*go, 1e6*gn, go/gn, np.max(np.abs(om_old.reshape(-1,d,d)-om_new))))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
