# Stage-level timing of the separable-scan pipeline on the PREM chord,
# order 4 GL, fixed slabs, replicating _osc_prob_scan_separable's internals.
import sys, json, warnings, time
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.magnus as mg
import magnus.expmkernels as ek
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc

warnings.simplefilter('ignore')
refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
ch = pcc.chord()
E = np.asarray(refs['energy_ev'], dtype=float)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
timed = gpb.timed
c_first = timed(gpb.control)

d = int(sys.argv[1])
n_slabs = 1088
order = 4

# Build the same inputs the scan builds
import magnus.hamiltonians as hams
params = gpb.osc_params(d)
# get H_E and h_matt through the same dispatch the scan uses: call once and intercept
H_E_holder = {}
orig = oscprob._osc_prob_scan_separable
def spy(H_E, VCC_func, h_matt, *a, **k):
    H_E_holder['H_E'] = H_E; H_E_holder['h_matt'] = h_matt
    H_E_holder['VCC_func'] = VCC_func
    return orig(H_E, VCC_func, h_matt, *a, **k)
oscprob._osc_prob_scan_separable = spy
np.asarray(oscprob.osc_prob_matter_std_potential(
    d, lambda x: ch['vcc'](x)/per_ne, E, ch['baseline'], params,
    L0=0.0, density_is_of_number_of_electrons=True, rtol=None, atol=None,
    n_slabs=n_slabs, magnus_exp_order=order, strategy='magnus',
    t_breakpoints=ch['edges'][1:-1], validate_input=False))
oscprob._osc_prob_scan_separable = orig
H_E, h_matt, VCC_func = H_E_holder['H_E'], H_E_holder['h_matt'], H_E_holder['VCC_func']
L0, L_val = 0.0, ch['baseline']
dim = d; nE = len(E)

# replicate the scan body
grid = np.linspace(L0, L_val, n_slabs + 1)
bp = np.atleast_1d(np.asarray(ch['edges'][1:-1], dtype=float))
bp = bp[(bp > L0) & (bp < L_val)]
grid = np.unique(np.concatenate([grid, bp]))
edges = np.column_stack([grid[:-1], grid[1:]])
widths = edges[:, 1] - edges[:, 0]
s_nodes = mg.gl_nodes(order)
tgrid = edges[:, :1] + widths[:, None]*s_nodes
V = np.asarray(VCC_func(tgrid.ravel())).reshape(tgrid.shape)
mA = -1j*h_matt.astype(complex)
HE_c = -1j*H_E.astype(complex)
Vmat = V[:, :, None, None]*mA
ns = len(widths)
print("d=%d n_slabs(actual)=%d nE=%d" % (d, ns, nE))
chunk, _ = oscprob._tile_for_working_set(nE, 1, tgrid.size*dim*dim)
print("chunk =", chunk)
sel = np.arange(min(chunk, nE))
nchunks = int(np.ceil(nE/len(sel)))

# Stage arrays (one chunk)
At = HE_c[sel][:, None, None, :, :] + Vmat[None, :, :, :, :]
h = widths[..., None, None]
A1 = At[..., 0, :, :]; A2 = At[..., 1, :, :]
Om = mg._magnus_gl(At, widths, order)
K = 1j*Om
U = mg._expm_stack(Om)
C = mg._commutator_batched(A2, A1)

def t_V():        return np.asarray(VCC_func(tgrid.ravel())).reshape(tgrid.shape)
def t_Vmat():     return V[:, :, None, None]*mA
def t_At():       return HE_c[sel][:, None, None, :, :] + Vmat[None, :, :, :, :]
def t_sident():   return mg._samples_identical(A1, A2)
def t_comm():     return mg._commutator_batched(A2, A1)
def t_lincomb():  return 0.5*h*(A1 + A2) + (np.sqrt(3.0)/12.0)*h*h*C
def t_gl():       return mg._magnus_gl(At, widths, order)
def t_frame():
    K = 1j*Om
    Kh = np.conj(np.swapaxes(K, -1, -2))
    s = np.max(np.abs(K))
    return np.max(np.abs(K - Kh)) <= 1e-12*s
def t_kernel():
    if ek.HAVE_NUMBA and ek.supports_dim(dim):
        return ek.expm_herm_stack(K)
    lam, Vv = np.linalg.eigh(K)
    Vh = np.conj(np.swapaxes(Vv, -1, -2))
    return (Vv*np.exp(-1j*lam)[..., None, :]) @ Vh
def t_eigh():     return np.linalg.eigh(K)
def t_expm():     return mg._expm_stack(Om, warn_wide=True, A_is_const=False)
def t_evops():    return mg.evolution_operators_from_samples(At, widths, order, 'gl', validate_input=False)
def t_ordprod():  return mg._ordered_product_batched(U)
def t_prob():     return np.swapaxes(U[:, 0].real**2 + U[:, 0].imag**2, -1, -2)

stages = [('V eval (per level, shared)', t_V, 1),
          ('Vmat build (shared)', t_Vmat, 1),
          ('At build', t_At, nchunks),
          ('  samples_identical', t_sident, nchunks),
          ('  commutator', t_comm, nchunks),
          ('  lin.comb (Omega)', t_lincomb, nchunks),
          ('magnus_gl total', t_gl, nchunks),
          ('  expm framing', t_frame, nchunks),
          ('  expm kernel/eigh', t_kernel, nchunks),
          ('  (eigh alone)', t_eigh, nchunks),
          ('expm_stack total', t_expm, nchunks),
          ('evops total', t_evops, nchunks),
          ('ordered product', t_ordprod, nchunks),
          ('prob extract', t_prob, nchunks)]
res = {}
for name, f, mult in stages:
    t = timed(f)
    res[name] = (t, mult)
    print("%-28s %8.3f us x%d = %8.1f us/call  (%6.4f us/slab/E)"
          % (name, 1e6*t, mult, 1e6*t*mult, 1e6*t*mult/ns/nE))
c_last = timed(gpb.control)
print("control ratio: %.3f" % (c_last/c_first))
