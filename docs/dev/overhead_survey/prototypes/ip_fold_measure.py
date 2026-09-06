# Old-vs-new timing at three levels: the fold alone, the engine, and a public
# osc_prob caller.  Arms interleaved, best-of blocks, control workload timed
# before/between/after so machine drift is visible.
import subprocess
import sys
import ast
import warnings
import time
import numpy as np

sys.path.insert(0, '/home/mbustamante/Research/magnus/src')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import magnus.oscprob as op
import magnus.magnus as mg
import magnus.matter as matter
import magnus.globaldefs as gd
import gen_profile_benchmarks as gpb

warnings.simplefilter('ignore')
timed = gpb.timed
controls = []


def snap_control(tag):
    controls.append((tag, timed(gpb.control)))
    print("  [control %-12s %.4f ms]" % (tag, 1e3*controls[-1][1]))


# --- old engine core, exactly as committed at HEAD -------------------------
src = subprocess.run(
    ['git', '-C', '/home/mbustamante/Research/magnus', 'show', 'HEAD:src/magnus/oscprob.py'],
    capture_output=True, text=True, check=True).stdout
fn = next(n for n in ast.walk(ast.parse(src))
          if isinstance(n, ast.FunctionDef) and n.name == '_osc_prob_ip_exp_core')
old_ns = dict(op.__dict__)
exec(compile(ast.Module(body=[fn], type_ignores=[]), '<old>', 'exec'), old_ns)
old_core = old_ns['_osc_prob_ip_exp_core']
new_core = op._osc_prob_ip_exp_core

snap_control('start')

# ==========================================================================
# Level 1: the fold alone, at the real tile shape (e_chunk=64, blk=32).
# ==========================================================================
l_scale = 1.0e12
rng = np.random.default_rng(0)
nE, dim, n_slabs = 64, 2, 32768
X = rng.normal(size=(nE, dim, dim)) + 1j*rng.normal(size=(nE, dim, dim))
H_E = ((X + np.conj(np.swapaxes(X, -1, -2)))/2)/l_scale
h_matt = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
v0 = 0.05/l_scale
VCC = lambda l: v0*np.exp(-np.asarray(l, dtype=float)/l_scale)
L0, L_val = 0.0, 10.0*l_scale

Lam, W = np.linalg.eigh(H_E)
Wd = np.conj(np.swapaxes(W, -1, -2))
Mt = Wd @ h_matt[None] @ W
Delta = Lam[:, :, None] - Lam[:, None, :]
denom = 1j*Delta - 1.0/l_scale
grid = np.linspace(L0, L_val, n_slabs + 1)
blk = 32
b1 = n_slabs; b0 = b1 - blk
edges0 = grid[b0:b1]; widths = grid[b0+1:b1+1] - edges0
V0 = np.asarray(VCC(edges0), dtype=complex)
arg = denom[:, None, :, :]*widths[None, :, None, None]
I = (np.exp(arg) - 1.0)/denom[:, None, :, :]
Omega_t = -1j*Mt[:, None, :, :]*V0[None, :, None, None]*I
U_free_diag = np.exp(-1j*Lam[:, None, :]*widths[None, :, None])
U_slab = U_free_diag[..., :, None]*mg._expm_stack(Omega_t, warn_wide=False)
acc0 = np.ascontiguousarray(U_slab[:, -1])
nblocks = n_slabs//blk


def fold_old():
    acc = acc0
    for k in range(U_slab.shape[1] - 1, -1, -1):
        acc = acc @ U_slab[:, k]
    return acc


acc_buf = acc0.copy()


def fold_new():
    acc_buf[...] = acc0
    mg._ordered_product_into(acc_buf, U_slab)
    return acc_buf


print("\n== Level 1: fold alone, one (64, 32, 2, 2) block, x%d blocks/pass ==" % nblocks)
t_old, t_new = [], []
for arm in range(3):
    t_old.append(timed(fold_old))
    t_new.append(timed(fold_new))
to, tn = min(t_old), min(t_new)
print("old fold: %7.1f us/block -> %6.1f ms/pass   (arms: %s)"
      % (1e6*to, 1e3*to*nblocks, ["%.1f" % (1e6*t) for t in t_old]))
print("new fold: %7.1f us/block -> %6.1f ms/pass   (arms: %s)"
      % (1e6*tn, 1e3*tn*nblocks, ["%.1f" % (1e6*t) for t in t_new]))
print("fold gain: %.2fx" % (to/tn))
snap_control('after-fold')

# ==========================================================================
# Level 2: the engine, one fixed pass at n_slabs=32768, nE=64.
# ==========================================================================
print("\n== Level 2: engine, one pass, n_slabs=32768, nE=64 ==")


def eng_old():
    return old_core(H_E, l_scale, VCC, h_matt, L0, L_val, None, None, 2.0, 25, 1, 10**6, 32768)


def eng_new():
    return new_core(H_E, l_scale, VCC, h_matt, L0, L_val, None, None, 2.0, 25, 1, 10**6, 32768)


t_old, t_new = [], []
for arm in range(3):
    t_old.append(timed(eng_old, repeat=3))
    t_new.append(timed(eng_new, repeat=3))
to, tn = min(t_old), min(t_new)
print("old engine: %7.1f ms/pass   (arms: %s)" % (1e3*to, ["%.0f" % (1e3*t) for t in t_old]))
print("new engine: %7.1f ms/pass   (arms: %s)" % (1e3*tn, ["%.0f" % (1e3*t) for t in t_new]))
print("engine gain: %.2fx" % (to/tn))
snap_control('after-engine')

# ==========================================================================
# Level 3: a public osc_prob caller (osc_prob_matter_std_potential).
# ==========================================================================
L_SCALE = gd.L_SCALE_SUN
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL
L1 = 1.0*L_SCALE
PARAMS_2NU = {'sth': 0.55, 'Dm2': 7.5e-5}
ne = matter.exp_density_profile(NE0, L_SCALE)
E64 = np.linspace(5.0e6, 15.0e6, 64)


def caller(**kw):
    return op.osc_prob_matter_std_potential(
        2, ne, E64, L1, PARAMS_2NU, L0=0.0,
        density_is_of_number_of_electrons=True, strategy='magnus', **kw)


# routing check: both arms must land on ip_exp
for tag, core in (('new', new_core), ('old', old_core)):
    op._osc_prob_ip_exp_core = core
    info = {}
    caller(rtol=None, atol=None, n_slabs=256, strategy_info=info)
    print("caller fixed-slab arm=%s -> engine %s" % (tag, info['engine']))
    assert info['engine'] == 'ip_exp'
    info = {}
    caller(rtol=1e-7, atol=1e-7, strategy_info=info)
    print("caller ladder     arm=%s -> engine %s" % (tag, info['engine']))
    assert info['engine'] == 'ip_exp'
op._osc_prob_ip_exp_core = new_core

print("\n== Level 3a: caller, fixed n_slabs=32768, 64 energies ==")
t_old, t_new = [], []
for arm in range(3):
    op._osc_prob_ip_exp_core = old_core
    t_old.append(timed(lambda: caller(rtol=None, atol=None, n_slabs=32768), repeat=3))
    op._osc_prob_ip_exp_core = new_core
    t_new.append(timed(lambda: caller(rtol=None, atol=None, n_slabs=32768), repeat=3))
to, tn = min(t_old), min(t_new)
print("old caller: %7.1f ms/call   (arms: %s)" % (1e3*to, ["%.0f" % (1e3*t) for t in t_old]))
print("new caller: %7.1f ms/call   (arms: %s)" % (1e3*tn, ["%.0f" % (1e3*t) for t in t_new]))
print("caller fixed-slab gain: %.2fx" % (to/tn))
snap_control('after-3a')

print("\n== Level 3b: caller, tolerance ladder rtol=atol=1e-7, 64 energies ==")
t_old, t_new = [], []
for arm in range(3):
    op._osc_prob_ip_exp_core = old_core
    t_old.append(timed(lambda: caller(rtol=1e-7, atol=1e-7), repeat=3))
    op._osc_prob_ip_exp_core = new_core
    t_new.append(timed(lambda: caller(rtol=1e-7, atol=1e-7), repeat=3))
to, tn = min(t_old), min(t_new)
print("old caller: %7.1f ms/call   (arms: %s)" % (1e3*to, ["%.0f" % (1e3*t) for t in t_old]))
print("new caller: %7.1f ms/call   (arms: %s)" % (1e3*tn, ["%.0f" % (1e3*t) for t in t_new]))
print("caller ladder gain: %.2fx" % (to/tn))
snap_control('end')

base = controls[0][1]
print("\ncontrol drift across the run: " +
      ", ".join("%s %.3f" % (tag, t/base) for tag, t in controls))
