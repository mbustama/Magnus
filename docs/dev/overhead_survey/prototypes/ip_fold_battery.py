# Old-vs-new battery for the compiled IP-engine fold.
# "Old" is _osc_prob_ip_exp_core exactly as committed at HEAD, exec'd into the
# live module's namespace so both versions share every helper and constant.
import subprocess
import sys
import ast
import warnings
import numpy as np

sys.path.insert(0, '/home/mbustamante/Research/magnus/src')
import magnus.oscprob as op
import magnus.magnus as mg

warnings.simplefilter('ignore')

src = subprocess.run(
    ['git', '-C', '/home/mbustamante/Research/magnus', 'show', 'HEAD:src/magnus/oscprob.py'],
    capture_output=True, text=True, check=True).stdout
tree = ast.parse(src)
fn = next(n for n in ast.walk(tree)
          if isinstance(n, ast.FunctionDef) and n.name == '_osc_prob_ip_exp_core')
old_ns = dict(op.__dict__)
exec(compile(ast.Module(body=[fn], type_ignores=[]), '<old>', 'exec'), old_ns)
old_core = old_ns['_osc_prob_ip_exp_core']
new_core = op._osc_prob_ip_exp_core

# --- engine inputs: 2-level, exponential profile (prof_ip.py's setup) ------
l_scale = 1.0e12
L = 10.0*l_scale
rng = np.random.default_rng(0)
h_matt = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
v0 = 0.05/l_scale
VCC = lambda l: v0*np.exp(-np.asarray(l, dtype=float)/l_scale)


def H_E_for(nE, seed=0):
    r = np.random.default_rng(seed)
    X = r.normal(size=(nE, 2, 2)) + 1j*r.normal(size=(nE, 2, 2))
    return ((X + np.conj(np.swapaxes(X, -1, -2)))/2)/l_scale


def run(core, nE, rtol, atol, n_slabs, seed=0):
    return core(H_E_for(nE, seed), l_scale, VCC, h_matt, 0.0, L, rtol, atol,
                2.0, 25, 1, 10**6, n_slabs)


worst = 0.0
worst_cfg = None
ncfg = 0
budget0 = op.BATCH_WORKING_ENTRIES

# fixed slab counts x energy counts x tiling budgets
for n_slabs in (8, 512, 4096, 32768):
    for nE in (4, 64):
        for budget in (budget0, 4096):
            op.BATCH_WORKING_ENTRIES = budget
            P_old, c_old = run(old_core, nE, None, None, n_slabs)
            P_new, c_new = run(new_core, nE, None, None, n_slabs)
            assert c_old == c_new
            d = float(np.max(np.abs(P_new - P_old)))
            ncfg += 1
            if d > worst:
                worst, worst_cfg = d, ('fixed', n_slabs, nE, budget)
            print("fixed n_slabs=%5d nE=%2d budget=%8d  max|dP| = %.3e" %
                  (n_slabs, nE, budget, d))
op.BATCH_WORKING_ENTRIES = budget0

# refinement ladders: tolerance requested, certification decisions must not move
for tol in (1e-3, 1e-6, 1e-9):
    for nE in (4, 64):
        for seed in (0, 1, 2):
            P_old, c_old = run(old_core, nE, tol, tol, 8, seed=seed)
            P_new, c_new = run(new_core, nE, tol, tol, 8, seed=seed)
            assert c_old == c_new, ("certification decision moved", tol, nE, seed)
            d = float(np.max(np.abs(P_new - P_old)))
            ncfg += 1
            if d > worst:
                worst, worst_cfg = d, ('ladder', tol, nE, seed)
            print("ladder tol=%.0e nE=%2d seed=%d  converged=%s  max|dP| = %.3e" %
                  (tol, nE, seed, c_new, d))

print("\n%d configurations, worst max|dP| = %.3e at %s" % (ncfg, worst, worst_cfg))

# --- numba-less fallback arm: must equal the old code exactly --------------
saved = mg._ordered_product_into_kernel
mg._ordered_product_into_kernel = None
worst_fb = 0.0
for n_slabs in (8, 512, 4096):
    for nE in (4, 64):
        P_old, c_old = run(old_core, nE, None, None, n_slabs)
        P_fb, c_fb = run(new_core, nE, None, None, n_slabs)
        assert c_old == c_fb
        eq = np.array_equal(P_fb, P_old)
        d = float(np.max(np.abs(P_fb - P_old)))
        worst_fb = max(worst_fb, d)
        print("fallback n_slabs=%5d nE=%2d  bit-identical=%s  max|dP| = %.3e" %
              (n_slabs, nE, eq, d))
mg._ordered_product_into_kernel = saved
print("fallback worst max|dP| = %.3e (must be exactly 0.0)" % worst_fb)
