# mpmath-referenced accuracy of old (BLAS chain) vs new (compiled) fold on the
# engine's own U_slab operators at n_slabs=32768, the battery's worst case.
import sys
import warnings
import numpy as np
import mpmath as mp

sys.path.insert(0, '/home/mbustamante/Research/magnus/src')
import magnus.magnus as mg

warnings.simplefilter('ignore')
mp.mp.dps = 40

l_scale = 1.0e12
L = 10.0*l_scale
nE, n_slabs = 4, 32768
rng = np.random.default_rng(0)
X = rng.normal(size=(nE, 2, 2)) + 1j*rng.normal(size=(nE, 2, 2))
H_E = ((X + np.conj(np.swapaxes(X, -1, -2)))/2)/l_scale
h_matt = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
v0 = 0.05/l_scale
VCC = lambda l: v0*np.exp(-np.asarray(l, dtype=float)/l_scale)

# U_slab exactly as the engine builds it (elementwise, so tiling-invariant)
Lam, W = np.linalg.eigh(H_E)
Wd = np.conj(np.swapaxes(W, -1, -2))
Mt = Wd @ h_matt[None] @ W
Delta = Lam[:, :, None] - Lam[:, None, :]
denom = 1j*Delta - 1.0/l_scale
grid = np.linspace(0.0, L, n_slabs + 1)
edges0 = grid[:-1]
widths = grid[1:] - edges0
V0 = np.asarray(VCC(edges0), dtype=complex)
arg = denom[:, None, :, :]*widths[None, :, None, None]
I = (np.exp(arg) - 1.0)/denom[:, None, :, :]
Omega_t = -1j*Mt[:, None, :, :]*V0[None, :, None, None]*I
U_free_diag = np.exp(-1j*Lam[:, None, :]*widths[None, :, None])
U_slab = U_free_diag[..., :, None]*mg._expm_stack(Omega_t, warn_wide=False)

# old fold: numpy @ chain, k descending
acc_old = None
for k in range(n_slabs - 1, -1, -1):
    acc_old = U_slab[:, k] if acc_old is None else acc_old @ U_slab[:, k]

# new fold: compiled kernel
acc_new = np.ascontiguousarray(U_slab[:, -1])
mg._ordered_product_into(acc_new, U_slab[:, :-1])

# exact fold in mpmath
for e in range(nE):
    a = [[mp.mpc(U_slab[e, -1, i, j]) for j in range(2)] for i in range(2)]
    for k in range(n_slabs - 2, -1, -1):
        b = U_slab[e, k]
        a = [[a[i][0]*mp.mpc(b[0, j]) + a[i][1]*mp.mpc(b[1, j]) for j in range(2)]
             for i in range(2)]
    exact = np.array([[complex(a[i][j]) for j in range(2)] for i in range(2)])
    e_old = float(np.max(np.abs(acc_old[e] - exact)))
    e_new = float(np.max(np.abs(acc_new[e] - exact)))
    d = float(np.max(np.abs(acc_old[e] - acc_new[e])))
    print("energy %d: |old-exact| = %.3e  |new-exact| = %.3e  |old-new| = %.3e" %
          (e, e_old, e_new, d))
