# Accuracy of the Jacobi backend against the mpmath references (numu->numu),
# plus unitarity, at fixed slab counts.
import sys, json, warnings
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.globaldefs as gd
import magnus.expmkernels as ek
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc
from jacobi_proto5 import _jacobi_expm_warm_mgs2 as _jacobi_expm_warm_mgs

warnings.simplefilter('ignore')
prem_refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
mp_refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/mp_reference_profile.json').read())
ch = pcc.chord()
E = np.asarray(prem_refs['energy_ev'], dtype=float)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
prof_exp = gpb.exponential_profile()

orig_supports = ek.supports_dim
orig_herm = ek.expm_herm_stack
def expm_herm_stack_j(K):
    d = K.shape[-1]
    if d in (2, 3):
        return orig_herm(K)
    flat = np.ascontiguousarray(K, dtype=complex).reshape(-1, d, d)
    out = np.empty_like(flat)
    lam = np.empty((flat.shape[0], d), dtype=float)
    _jacobi_expm_warm_mgs(flat, out, lam)
    return out.reshape(K.shape), lam.reshape(K.shape[:-1]), 0.0
def set_mode(on):
    ek.supports_dim = (lambda d: d in (2,3,4,5)) if on else orig_supports
    ek.expm_herm_stack = expm_herm_stack_j if on else orig_herm

def prem_call(d, n_slabs, rtol=None, order=4, im='gl'):
    kw = dict(rtol=None, atol=None) if rtol is None else dict(rtol=rtol, atol=rtol*1e-2)
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        d, lambda x: ch['vcc'](x)/per_ne, E, ch['baseline'], gpb.osc_params(d),
        L0=0.0, density_is_of_number_of_electrons=True,
        n_slabs=n_slabs, magnus_exp_order=order, integration_method=im,
        strategy='magnus', t_breakpoints=ch['edges'][1:-1],
        validate_input=False, **kw))
def exp_call(d, n_slabs, rtol=None, order=4, im='gl'):
    kw = dict(rtol=None, atol=None) if rtol is None else dict(rtol=rtol, atol=rtol*1e-2)
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        d, lambda x: prof_exp['vcc'](x)/per_ne, prof_exp['energies'],
        prof_exp['baseline'], gpb.osc_params(d), L0=0.0,
        density_is_of_number_of_electrons=True,
        n_slabs=n_slabs, magnus_exp_order=order, integration_method=im,
        strategy='magnus', validate_input=False, **kw))

def unit(P):
    return np.max(np.abs(np.sum(P, axis=-1) - 1.0))

worst = 0.0
print("profile  d  n_slabs   err_eigh     err_jacobi   shift(fullP)  unit_e     unit_j")
for name, callf, refs in (('PREM', prem_call, {c['flavours']: np.array(c['reference']) for c in prem_refs['cases']}),
                          ('expo', exp_call, {c['flavours']: np.array(c['reference']) for c in mp_refs['cases']})):
    for d in (4, 5):
        for n in (1088, 4352):
            set_mode(False); old = callf(d, n)
            set_mode(True);  new = callf(d, n)
            set_mode(False)
            ref = refs[d]
            eo = np.max(np.abs(old[:, gd.NUMU, gd.NUMU] - ref))
            en = np.max(np.abs(new[:, gd.NUMU, gd.NUMU] - ref))
            sh = np.max(np.abs(new - old)); worst = max(worst, sh)
            print("%-6s  %d  %6d   %.4e   %.4e   %.4e   %.2e   %.2e"
                  % (name, d, n, eo, en, sh, unit(old), unit(new)))
print("worst full-matrix shift:", worst)
