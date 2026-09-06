# Level 3 only: the public osc_prob caller, fixed-slab and ladder arms.
# Ladder moved to 40-100 MeV at rtol=atol=1e-3, where ip_exp certifies.
import subprocess
import sys
import ast
import warnings
import numpy as np

sys.path.insert(0, '/home/mbustamante/Research/magnus/src')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import magnus.oscprob as op
import magnus.matter as matter
import magnus.globaldefs as gd
import gen_profile_benchmarks as gpb

warnings.simplefilter('ignore')
timed = gpb.timed
controls = []


def snap_control(tag):
    controls.append((tag, timed(gpb.control)))
    print("  [control %-12s %.4f ms]" % (tag, 1e3*controls[-1][1]))


src = subprocess.run(
    ['git', '-C', '/home/mbustamante/Research/magnus', 'show', 'HEAD:src/magnus/oscprob.py'],
    capture_output=True, text=True, check=True).stdout
fn = next(n for n in ast.walk(ast.parse(src))
          if isinstance(n, ast.FunctionDef) and n.name == '_osc_prob_ip_exp_core')
old_ns = dict(op.__dict__)
exec(compile(ast.Module(body=[fn], type_ignores=[]), '<old>', 'exec'), old_ns)
old_core = old_ns['_osc_prob_ip_exp_core']
new_core = op._osc_prob_ip_exp_core

L_SCALE = gd.L_SCALE_SUN
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL
L1 = 1.0*L_SCALE
PARAMS_2NU = {'sth': 0.55, 'Dm2': 7.5e-5}
ne = matter.exp_density_profile(NE0, L_SCALE)
E_fix = np.linspace(5.0e6, 15.0e6, 64)
E_lad = np.linspace(40.0e6, 100.0e6, 64)


def caller(E, **kw):
    return op.osc_prob_matter_std_potential(
        2, ne, E, L1, PARAMS_2NU, L0=0.0,
        density_is_of_number_of_electrons=True, strategy='magnus', **kw)


# warm both cores, then routing checks for both arms
for tag, core in (('new', new_core), ('old', old_core)):
    op._osc_prob_ip_exp_core = core
    info = {}
    caller(E_fix, rtol=None, atol=None, n_slabs=256, strategy_info=info)
    print("caller fixed-slab arm=%s -> engine %s" % (tag, info['engine']))
    assert info['engine'] == 'ip_exp', info['engine']
    info = {}
    caller(E_lad, rtol=1e-3, atol=1e-3, strategy_info=info)
    print("caller ladder     arm=%s -> engine %s" % (tag, info['engine']))
    assert info['engine'] == 'ip_exp', info['engine']
op._osc_prob_ip_exp_core = new_core
snap_control('start')

print("\n== Level 3a: caller, fixed n_slabs=32768, 64 energies (5-15 MeV) ==")
t_old, t_new = [], []
for arm in range(3):
    op._osc_prob_ip_exp_core = old_core
    t_old.append(timed(lambda: caller(E_fix, rtol=None, atol=None, n_slabs=32768), repeat=3))
    op._osc_prob_ip_exp_core = new_core
    t_new.append(timed(lambda: caller(E_fix, rtol=None, atol=None, n_slabs=32768), repeat=3))
to, tn = min(t_old), min(t_new)
print("old caller: %7.1f ms/call   (arms: %s)" % (1e3*to, ["%.0f" % (1e3*t) for t in t_old]))
print("new caller: %7.1f ms/call   (arms: %s)" % (1e3*tn, ["%.0f" % (1e3*t) for t in t_new]))
print("caller fixed-slab gain: %.2fx" % (to/tn))
snap_control('after-3a')

print("\n== Level 3b: caller, ladder rtol=atol=1e-3, 64 energies (40-100 MeV) ==")
t_old, t_new = [], []
for arm in range(3):
    op._osc_prob_ip_exp_core = old_core
    t_old.append(timed(lambda: caller(E_lad, rtol=1e-3, atol=1e-3), repeat=3))
    op._osc_prob_ip_exp_core = new_core
    t_new.append(timed(lambda: caller(E_lad, rtol=1e-3, atol=1e-3), repeat=3))
to, tn = min(t_old), min(t_new)
print("old caller: %7.1f ms/call   (arms: %s)" % (1e3*to, ["%.0f" % (1e3*t) for t in t_old]))
print("new caller: %7.1f ms/call   (arms: %s)" % (1e3*tn, ["%.0f" % (1e3*t) for t in t_new]))
print("caller ladder gain: %.2fx" % (to/tn))
snap_control('end')

base = controls[0][1]
print("\ncontrol drift: " + ", ".join("%s %.3f" % (tag, t/base) for tag, t in controls))
