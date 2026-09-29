import sys, warnings, numpy as np
from scipy.integrate import solve_ivp
warnings.simplefilter('ignore')
exec(open('check_l1.py').read().split("mode = sys.argv")[0])
import magnus.oscprob as op
ref = np.load('l1_ref.npz')
CAP = {}
orig = op._hybrid_propagator_scan
def spy(H, *a, **k):
    CAP['H'] = H; return orig(H, *a, **k)
op._hybrid_propagator_scan = spy
tol = dict(rtol=1e-12, atol=1e-14)
for k, (f, E, p) in curves.items():
    P8 = np.asarray(f(E, **base, **p, **tol, strategy='magnus', magnus_exp_order=8))
    d = np.abs(P8 - ref[k]).reshape(len(E), -1).max(axis=1)
    i = int(np.argmax(d))
    f(E[i:i+1], **base, **p, **tol, strategy='hybrid')
    Hf = CAP['H'](E[i]); dim = np.asarray(Hf(0.0)).shape[-1]
    def rhs(x, y):
        U = y.reshape(dim, dim, 2); dU = -1j*np.asarray(Hf(x), dtype=complex) @ (U[..., 0] + 1j*U[..., 1])
        return np.stack([dU.real, dU.imag], -1).ravel()
    L = base['L']
    s = solve_ivp(rhs, (0.0, L), np.stack([np.eye(dim), np.zeros((dim, dim))], -1).ravel(), method='DOP853', rtol=1e-13, atol=1e-15)
    U = s.y[:, -1].reshape(dim, dim, 2); U = U[..., 0] + 1j*U[..., 1]
    Pd = (np.abs(U)**2).T
    # the listing's channel: nu_e -> nu_e, the (0, 0) entry
    h4, l8, dd = float(np.asarray(ref[k])[i].ravel()[0]), float(P8[i].ravel()[0]), float(Pd[0, 0])
    print('%s E=%.4g GeV  hybrid-o4=%.15f  ladder-o8=%.15f  DOP853=%.15f  |h4-dop|=%.1e  |l8-dop|=%.1e' % (k, E[i]/gd.UNIT_GEV, h4, l8, dd, abs(h4-dd), abs(l8-dd)))
