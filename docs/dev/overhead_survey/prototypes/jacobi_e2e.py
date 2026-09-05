# End-to-end: PREM chord scan with the Jacobi kernel patched in for d=4,5.
# A/B alternated, control interleaved; probabilities scored against the
# mpmath reference for both arms.
import sys, json, warnings
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.magnus as mg
import magnus.expmkernels as ek
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc
from jacobi_proto5 import _jacobi_expm_warm_mgs2 as _jacobi_expm_warm_mgs

warnings.simplefilter('ignore')
refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
ch = pcc.chord()
E = np.asarray(refs['energy_ev'], dtype=float)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
timed = gpb.timed
c0 = timed(gpb.control)

orig_supports = ek.supports_dim
orig_herm = ek.expm_herm_stack
def supports_dim_j(d):
    return d in (2, 3, 4, 5)
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
    ek.supports_dim = supports_dim_j if on else orig_supports
    ek.expm_herm_stack = expm_herm_stack_j if on else orig_herm

def prem_call(d, n_slabs):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        d, lambda x: ch['vcc'](x)/per_ne, E, ch['baseline'], gpb.osc_params(d),
        L0=0.0, density_is_of_number_of_electrons=True, rtol=None, atol=None,
        n_slabs=n_slabs, magnus_exp_order=4, strategy='magnus',
        t_breakpoints=ch['edges'][1:-1], validate_input=False))

cases = {c['flavours']: np.asarray(c['reference']) for c in refs['cases']}
NS = (544, 1088, 2176)
for d in (4, 5):
    ref = cases[d]
    res = {}
    for n in NS:
        for mode in ('eigh', 'jacobi'):
            set_mode(mode == 'jacobi')
            P = prem_call(d, n)
            # P_ee across the scan (reference stores one probability per energy)
            pe = P[:, 0, 0] if P.ndim == 3 else P
            t = timed(lambda: prem_call(d, n))
            res[(n, mode)] = (t, pe)
    set_mode(False)
    for n in NS:
        te, pe = res[(n, 'eigh')]; tj, pj = res[(n, 'jacobi')]
        print("d=%d n=%4d  eigh %7.2f ms  jacobi %7.2f ms  speedup %4.2fx  |dP(two arms)| %.2e  |P-ref| eigh %.2e jac %.2e"
              % (d, n, 1e3*te, 1e3*tj, te/tj, np.max(np.abs(pe - pj)),
                 np.max(np.abs(pe - ref)), np.max(np.abs(pj - ref))))
    n1, n2 = NS[-2] + 16, NS[-1] + 16
    me = (res[(NS[-1],'eigh')][0] - res[(NS[-2],'eigh')][0])/(n2 - n1)/12
    mj = (res[(NS[-1],'jacobi')][0] - res[(NS[-2],'jacobi')][0])/(n2 - n1)/12
    print("   marginal us/slab/E: eigh %6.3f  jacobi %6.3f  (%4.2fx)" % (1e6*me, 1e6*mj, me/mj))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
