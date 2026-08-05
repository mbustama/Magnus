"""Re-verify handover section 2.6 (grid palindromy) and 2.1 (U_j does not mirror)."""
import numpy as np
import common
import magnus.magnus as mg

L, A, H, vcc = common.chord_setup()

print("=" * 74)
print("2.6  Is the Earth chord's uniform grid palindromic?")
print("=" * 74)
for n in (64, 128):
    edges = common.uniform_edges(L, n)
    w = edges[:, 1] - edges[:, 0]
    mid = 0.5*(edges[:, 0] + edges[:, 1])
    dw = np.max(np.abs(w - w[::-1]))
    dm = np.max(np.abs(mid - (L - mid)[::-1]))
    print(f"  n={n:4d}  slab width      = {w[0]:.4e}")
    print(f"          max|w - w[::-1]|= {dw:.4e}   RELATIVE {dw/w[0]:.2e}")
    print(f"          max|mid-mirror| = {dm:.4e}   relative {dm/L:.2e}")
    print(f"          np.array_equal(w, w[::-1]) = {np.array_equal(w, w[::-1])}")
    edges_s, ws = common.symmetrised_edges(L, n)
    print(f"          after w=(w+w[::-1])/2: array_equal = "
          f"{np.array_equal(ws, ws[::-1])}")

print()
print("=" * 74)
print("2.1  max_j |U_j - U_{n-1-j}| on a PERFECTLY symmetric profile")
print("     (12 slabs, n_tpts_per_slab=41, trapezoid)")
print("=" * 74)

# A synthetic profile that is symmetric about L/2 by construction.
h_vac = H(0.0) - vcc(0.0)*np.diag([1.0, 0, 0])
proj = np.zeros((3, 3))
proj[0][0] = 1.0
V0 = 4.0e-13


def A_sym(l):
    l = np.asarray(l, dtype=float)
    v = V0*(1.0 + 0.5*np.cos(2.0*np.pi*l/L))          # v(l) == v(L-l)
    return -1j*(h_vac + v[..., None, None]*proj)


# confirm the profile really is symmetric to machine precision
probe = np.linspace(0.0, L, 1001)
asym = np.max(np.abs(A_sym(probe) - A_sym(L - probe)))
scale = np.max(np.abs(A_sym(probe)))
print(f"  profile asymmetry max|A(l)-A(L-l)| = {asym:.3e}  "
      f"(relative {asym/scale:.2e})")
print()

n_slabs = 12
edges = common.uniform_edges(L, n_slabs)
print(f"  {'order':>6} | {'max_j |U_j - U_{n-1-j}|':>26} | verdict")
print("  " + "-" * 62)
for order in (1, 2, 4, 6):
    U = mg.magnus_expansion_multislab(A_sym, edges, n_tpts_per_slab=41,
                                      order=order, integration_method='trapezoid')
    d = np.max(np.abs(U - U[::-1]))
    verdict = "mirrors" if d < 1e-10 else "BROKEN"
    print(f"  {order:>6} | {d:>26.3e} | {verdict}")

print()
print("  Same, with integration_method='gl':")
for order in (1, 2, 4, 6):
    U = mg.magnus_expansion_multislab(A_sym, edges, order=order,
                                      integration_method='gl')
    d = np.max(np.abs(U - U[::-1]))
    verdict = "mirrors" if d < 1e-10 else "BROKEN"
    print(f"  {order:>6} | {d:>26.3e} | {verdict}")
