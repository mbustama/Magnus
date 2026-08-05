"""Handover 2.2(B), done with MATCHED slab widths so the identity is what is
under test rather than a discretisation mismatch between U and F."""
import numpy as np
import common
import magnus.magnus as mg

L, A, H, vcc = common.chord_setup()
proj = np.zeros((3, 3)); proj[0][0] = 1.0
V0 = 4.0e-13
NH = 64          # slabs on the half; the full gets 2*NH, so widths match exactly


def make_A(dCP):
    _, _, Hd, vd = common.chord_setup(dCP=dCP)
    hv = Hd(0.0) - vd(0.0)*np.diag([1.0, 0, 0])

    def Af(l):
        l = np.asarray(l, dtype=float)
        v = V0*(1.0 + 0.5*np.cos(2.0*np.pi*l/L))     # symmetric about L/2
        return -1j*(hv + v[..., None, None]*proj)
    return Af, hv


def prod_U(Afun, t0, t1, n, order=4):
    e = np.linspace(t0, t1, n + 1)
    ed = np.stack([e[:-1], e[1:]], axis=1)
    U = mg.magnus_expansion_multislab(Afun, ed, order=order,
                                      integration_method='gl',
                                      validate_input=False)
    out = np.eye(3, dtype=complex)
    for k in range(n):
        out = U[k] @ out
    return out


print("=" * 74)
print("Transpose identity, matched slab widths (half: %d slabs, full: %d)" % (NH, 2*NH))
print("=" * 74)
print(f"  {'dCP':>8} | {'max|H-H^T|':>12} | {'rel to |H|':>11} | {'max|U - F^T F|':>16}")
print("  " + "-" * 60)
for dCP in (0.0, common.DCP):
    Af, hv = make_A(dCP)
    Hs = hv + V0*proj
    asym = np.max(np.abs(Hs - Hs.T))
    U = prod_U(Af, 0.0, L, 2*NH)
    F = prod_U(Af, 0.0, L/2, NH)
    print(f"  {dCP:>8.4f} | {asym:>12.3e} | {asym/np.max(np.abs(Hs)):>11.2e} "
          f"| {np.max(np.abs(U - F.T @ F)):>16.3e}")

print()
Ap, _ = make_A(common.DCP)
Am, _ = make_A(-common.DCP)
U = prod_U(Ap, 0.0, L, 2*NH)
Fp = prod_U(Ap, 0.0, L/2, NH)
Fm = prod_U(Am, 0.0, L/2, NH)
print(f"  Generalisation  max|U(dCP) - F(-dCP)^T F(dCP)| = "
      f"{np.max(np.abs(U - Fm.T @ Fp)):.3e}")
print(f"  (unitarity check max|U^dag U - I|            = "
      f"{np.max(np.abs(U.conj().T @ U - np.eye(3))):.3e})")
