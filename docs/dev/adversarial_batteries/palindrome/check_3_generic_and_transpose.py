"""(a) Sign rule on a generic random A(t) and at orders 7-8 (generated recursion).
(b) Handover 2.2(B): the transpose identity and its dependence on H^T = H.
"""
import numpy as np
import common
import magnus.magnus as mg

rng = np.random.default_rng(20260805)

print("=" * 78)
print("A.  Sign rule on a GENERIC random smooth A(t) -- no PREM structure")
print("=" * 78)

d = 4
c0, c1, c2 = (rng.normal(size=(3, d, d)) + 1j*rng.normal(size=(3, d, d)))


def A_rand(t):
    t = np.asarray(t, dtype=float)[..., None, None]
    return c0 + c1*np.sin(3.0*t) + c2*np.cos(2.0*t)*t


edges = np.array([[0.0, 0.37], [0.37, 0.81], [0.81, 1.25]])
widths = edges[:, 1] - edges[:, 0]

for order in (2, 4, 6, 7, 8):
    print(f"\n  order {order}, trapezoid")
    print(f"    {'n_tpts':>7} | {'max|Om_mirror - signrule|':>26} | {'relative':>10} | slope")
    print("    " + "-" * 60)
    prev = None
    for n_tpts in (21, 41, 81, 161, 321):
        s = np.linspace(0.0, 1.0, n_tpts)
        tg = edges[:, :1] + widths[:, None]*s
        At = A_rand(tg.ravel()).reshape(len(edges), n_tpts, d, d)
        Bt = widths[:, None, None, None]*At
        T_f = mg._magnus_terms_quadrature(Bt, order, 'trapezoid')
        T_m = mg._magnus_terms_quadrature(Bt[:, ::-1], order, 'trapezoid')
        sg = np.array([(-1.0)**(k + 1) for k in range(1, order + 1)])
        e = np.max(np.abs(T_m.sum(0) - (sg[:, None, None, None]*T_f).sum(0)))
        sc = np.max(np.abs(T_m.sum(0)))
        slope = "" if prev is None else f"{prev/e:6.2f}x"
        print(f"    {n_tpts:>7} | {e:>26.4e} | {e/sc:>10.2e} | {slope}")
        prev = e

print()
print("=" * 78)
print("B.  2.2(B) transpose identity  U(L,0) = U(L/2,0)^T @ U(L/2,0)")
print("=" * 78)

L, A, H, vcc = common.chord_setup()
h_vac = H(0.0) - vcc(0.0)*np.diag([1.0, 0, 0])
proj = np.zeros((3, 3))
proj[0][0] = 1.0
V0 = 4.0e-13


def make_A(dCP):
    Ls, Ad, Hd, vd = common.chord_setup(dCP=dCP)
    hv = Hd(0.0) - vd(0.0)*np.diag([1.0, 0, 0])

    def Af(l):
        l = np.asarray(l, dtype=float)
        v = V0*(1.0 + 0.5*np.cos(2.0*np.pi*l/L))
        return -1j*(hv + v[..., None, None]*proj), hv
    return Af


def total_U(Afun, t0, t1, n=64, order=4):
    e = np.linspace(t0, t1, n + 1)
    ed = np.stack([e[:-1], e[1:]], axis=1)
    U = mg.magnus_expansion_multislab(lambda x: Afun(x)[0], ed, order=order,
                                      integration_method='gl',
                                      validate_input=False)
    out = np.eye(3, dtype=complex)
    for k in range(n):
        out = U[k] @ out
    return out


print(f"  {'dCP':>8} | {'max|H - H^T|':>14} | {'max|U - F^T F|':>16}")
print("  " + "-" * 46)
for dCP in (0.0, common.DCP):
    Af = make_A(dCP)
    hv = Af(0.0)[1]
    Hs = hv + V0*proj
    asym = np.max(np.abs(Hs - Hs.T))
    U = total_U(Af, 0.0, L)
    F = total_U(Af, 0.0, L/2)
    err = np.max(np.abs(U - F.T @ F))
    print(f"  {dCP:>8.4f} | {asym:>14.3e} | {err:>16.3e}")

print()
print("  Generalisation  U(dCP) = F(-dCP)^T @ F(dCP):")
Ap = make_A(common.DCP)
Am = make_A(-common.DCP)
U = total_U(Ap, 0.0, L)
Fp = total_U(Ap, 0.0, L/2)
Fm = total_U(Am, 0.0, L/2)
print(f"    max|U(dCP) - F(-dCP)^T F(dCP)| = {np.max(np.abs(U - Fm.T @ Fp)):.3e}")
