"""Proposed Listing 1 (strategy='magnus', order 8, rtol 1e-12, atol 1e-14) against DOP853 at the 26
reference energies of Fig. 1, per curve; one process per curve."""
import sys, json, time, warnings, numpy as np
from scipy.integrate import solve_ivp
from multiprocessing import Pool
warnings.simplefilter('ignore')

def run(k):
    exec(open('check_l1.py').read().split("mode = sys.argv")[0], globals())
    import magnus.oscprob as op
    CAP = {}
    orig = op._hybrid_propagator_scan
    def spy(H, *a, **kw):
        CAP['H'] = H; return orig(H, *a, **kw)
    op._hybrid_propagator_scan = spy
    f, E_plot, p = curves[k]
    lo, hi = E_plot[0], E_plot[-1]
    E_ref = np.logspace(np.log10(lo), np.log10(hi), 26)
    t = time.perf_counter()
    P = np.asarray(f(E_ref, **base, **p, rtol=1e-12, atol=1e-14, strategy='magnus', magnus_exp_order=8))
    t_lad = time.perf_counter() - t
    rows = []
    for i, E in enumerate(E_ref):
        f(E, **base, **p, rtol=1e-3, atol=1e-3, strategy='hybrid')     # only to capture H(E, l)
        Hf = CAP['H'](E); d = np.asarray(Hf(0.0)).shape[-1]
        def rhs(x, y):
            U = y.reshape(d, d, 2); dU = -1j*np.asarray(Hf(x), dtype=complex) @ (U[..., 0] + 1j*U[..., 1])
            return np.stack([dU.real, dU.imag], -1).ravel()
        def dop(rt, at):
            s = solve_ivp(rhs, (0.0, base['L']), np.stack([np.eye(d), np.zeros((d, d))], -1).ravel(), method='DOP853', rtol=rt, atol=at)
            U = s.y[:, -1].reshape(d, d, 2); U = U[..., 0] + 1j*U[..., 1]
            return (np.abs(U)**2).T
        r12, r13 = dop(1e-12, 1e-14), dop(1e-13, 1e-15)
        Pm = np.asarray(P[i]).reshape(r13.shape) if np.asarray(P[i]).size == r13.size else None
        pee = float(np.asarray(P[i]).ravel()[0])
        rows.append(dict(E=float(E), err_ee=abs(pee - r13[0, 0]), floor=float(np.max(np.abs(r12 - r13))),
                         err_full=(float(np.max(np.abs(Pm - r13))) if Pm is not None else None)))
    return k, t_lad, rows

if __name__ == '__main__':
    with Pool(4) as pool:
        for k, t_lad, rows in pool.imap_unordered(run, ['2nu', '3nu', '4nu', '5nu']):
            e = [r['err_ee'] for r in rows]; fl = [r['floor'] for r in rows]
            print('%s  ladder(26 E)=%.3f s  max|P_ee - DOP853(1e-13)| = %.1e (median %.1e)  oracle floor max %.1e  worst at %.4g GeV'
                  % (k, t_lad, max(e), float(np.median(e)), max(fl), rows[int(np.argmax(e))]['E']/1e9), flush=True)
            json.dump(rows, open('ladder_check_%s.json' % k, 'w'))
