# Stage-level timing of the IP core body, replicated at the real tile shape.
import sys, warnings, time
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.magnus as mg
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb

warnings.simplefilter('ignore')
timed = gpb.timed
c0 = timed(gpb.control)

l_scale = 1.0e12
rng = np.random.default_rng(0)
nE, dim = 64, 2
X = rng.normal(size=(nE, dim, dim)) + 1j*rng.normal(size=(nE, dim, dim))
H_E = ((X + np.conj(np.swapaxes(X, -1, -2)))/2)/l_scale
h_matt = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
v0 = 0.05/l_scale
VCC = lambda l: v0*np.exp(-np.asarray(l, dtype=float)/l_scale)
L0, L_val = 0.0, 10.0*l_scale
n_slabs = 32768

Lam, W = np.linalg.eigh(H_E)
Wd = np.conj(np.swapaxes(W, -1, -2))
Mt = Wd @ h_matt[None] @ W
Delta = Lam[:, :, None] - Lam[:, None, :]
denom = 1j*Delta - 1.0/l_scale
e_chunk, blk = oscprob._tile_for_working_set(nE, n_slabs, dim*dim, live_arrays=8)
print("e_chunk=%d blk=%d blocks=%d" % (e_chunk, blk, int(np.ceil(n_slabs/blk))))
grid = np.linspace(L0, L_val, n_slabs + 1)
esel = slice(0, e_chunk)
b1 = n_slabs; b0 = b1 - blk
edges0 = grid[b0:b1]; widths = grid[b0+1:b1+1] - edges0
nblocks = int(np.ceil(n_slabs/blk))

V0 = np.asarray(VCC(edges0), dtype=complex)
arg = denom[esel, None, :, :]*widths[None, :, None, None]
I = (np.exp(arg) - 1.0)/denom[esel, None, :, :]
Omega_t = -1j*Mt[esel, None, :, :]*V0[None, :, None, None]*I
U_free_diag = np.exp(-1j*Lam[esel, None, :]*widths[None, :, None])
U_slab = U_free_diag[..., :, None]*mg._expm_stack(Omega_t, warn_wide=False)
acc0 = np.ascontiguousarray(U_slab[:, -1])

def t_v0():    return np.asarray(VCC(edges0), dtype=complex)
def t_argI():
    arg = denom[esel, None, :, :]*widths[None, :, None, None]
    return (np.exp(arg) - 1.0)/denom[esel, None, :, :]
def t_omega():
    Om = -1j*Mt[esel, None, :, :]*V0[None, :, None, None]*I
    return Om, float(np.max(np.abs(Om)))
def t_ufree():  return np.exp(-1j*Lam[esel, None, :]*widths[None, :, None])
def t_expm():   return mg._expm_stack(Omega_t, warn_wide=False)
def t_mult():   return U_free_diag[..., :, None]*mg._expm_stack(Omega_t, warn_wide=False)
def t_fold():
    acc = acc0
    for k in range(U_slab.shape[1] - 1, -1, -1):
        acc = acc @ U_slab[:, k]
    return acc
tot = 0.0
for name, f in (('V0 eval', t_v0), ('arg+exp+I', t_argI), ('Omega build+max', t_omega),
                ('U_free exp', t_ufree), ('expm_stack', t_expm),
                ('mult (incl expm again)', t_mult), ('fold (python @ loop)', t_fold)):
    t = timed(f)
    if name != 'mult (incl expm again)':
        tot += t
    print("%-24s %9.1f us/block  -> %7.1f ms per pass (%5.3f us/slab/E-chunk)"
          % (name, 1e6*t, 1e3*t*nblocks, 1e6*t*nblocks/n_slabs))
print("sum of stages: %.1f ms per pass (measured pass earlier: ~757 ms)" % (1e3*tot*nblocks))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
