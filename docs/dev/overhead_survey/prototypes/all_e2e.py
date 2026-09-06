# End-to-end: all three prototype changes together vs current code.
import sys, json, warnings
sys.path.insert(0, '.')
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.magnus as mg
import magnus.expmkernels as ek
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc
from fusion_proto2 import _antiherm_scale_dev2, _gl4_omega_from_At
from jacobi_proto5 import _jacobi_expm_warm_mgs2

warnings.simplefilter('ignore')
refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
ch = pcc.chord()
E = np.asarray(refs['energy_ev'], dtype=float)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
timed = gpb.timed
c0 = timed(gpb.control)

orig_gl = mg._magnus_gl
orig_expm = mg._expm_stack

def gl_fused(An, widths, order):
    if (order in (3, 4) and An.ndim >= 3 and An.shape[-3] == 2
            and An.dtype == np.complex128 and An.flags.c_contiguous):
        d = An.shape[-1]
        lead = An.shape[:-3]
        w = np.asarray(widths, dtype=float)
        hb = np.ascontiguousarray(np.broadcast_to(w, lead)).reshape(-1)
        out = np.empty(lead + (d, d), dtype=complex)
        _gl4_omega_from_At(An.reshape((-1, 2, d, d)), hb, out.reshape((-1, d, d)))
        return out
    return orig_gl(An, widths, order)

def expm_fused(Om, warn_wide=False, A_is_const=False, expm_backend=None):
    backend = mg._resolve_expm_backend(expm_backend)
    Om = np.asarray(Om)
    if Om.dtype == np.complex128 and Om.ndim >= 3 and Om.flags.c_contiguous:
        d = Om.shape[-1]
        scale, dev = _antiherm_scale_dev2(Om.reshape((-1, d, d)))
        if scale == 0.0:
            return np.broadcast_to(np.eye(d, dtype=complex), Om.shape).copy()
        if dev <= 1.e-12*scale:
            K = 1j*Om
            if backend != 'eigh' and ek.HAVE_NUMBA and d in (2, 3):
                U, lam, sev = ek.expm_herm_stack(K)
                if sev > ek.SEV_TOL:
                    lam, V = np.linalg.eigh(K)
                    Vh = np.conj(np.swapaxes(V, -1, -2))
                    U = (V*np.exp(-1j*lam)[..., None, :]) @ Vh
            elif backend != 'eigh' and ek.HAVE_NUMBA and d in (4, 5):
                flat = np.ascontiguousarray(K).reshape(-1, d, d)
                U = np.empty_like(flat)
                lam = np.empty((flat.shape[0], d), dtype=float)
                _jacobi_expm_warm_mgs2(flat, U, lam)
                U = U.reshape(K.shape); lam = lam.reshape(K.shape[:-1])
            else:
                lam, V = np.linalg.eigh(K)
                Vh = np.conj(np.swapaxes(V, -1, -2))
                U = (V*np.exp(-1j*lam)[..., None, :]) @ Vh
            if warn_wide and not A_is_const:
                mg._warn_slab_norm(np.max(np.abs(lam)))
            return U
    return orig_expm(Om, warn_wide=warn_wide, A_is_const=A_is_const,
                     expm_backend=expm_backend)

def set_mode(on):
    mg._magnus_gl = gl_fused if on else orig_gl
    mg._expm_stack = expm_fused if on else orig_expm

def prem_call(d, n_slabs):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        d, lambda x: ch['vcc'](x)/per_ne, E, ch['baseline'], gpb.osc_params(d),
        L0=0.0, density_is_of_number_of_electrons=True, rtol=None, atol=None,
        n_slabs=n_slabs, magnus_exp_order=4, strategy='magnus',
        t_breakpoints=ch['edges'][1:-1], validate_input=False))

NPE = {2: 0.049, 3: 0.130, 4: 1.70}
NS = (544, 1088, 2176)
for d in (2, 3, 4, 5):
    res = {}
    for n in NS:
        for mode in ('old', 'new'):
            set_mode(mode == 'new')
            P = prem_call(d, n)
            t = timed(lambda: prem_call(d, n))
            res[(n, mode)] = (t, P)
    set_mode(False)
    n1, n2 = NS[-2] + 16, NS[-1] + 16
    mo = (res[(NS[-1],'old')][0] - res[(NS[-2],'old')][0])/(n2 - n1)/12
    mn = (res[(NS[-1],'new')][0] - res[(NS[-2],'new')][0])/(n2 - n1)/12
    dP = max(np.max(np.abs(res[(n,'old')][1] - res[(n,'new')][1])) for n in NS)
    npe = NPE.get(d)
    tail = ("  [NPE %.3f: %0.2fx -> %0.2fx]" % (npe, 1e6*mo/npe, 1e6*mn/npe)) if npe else ""
    print("d=%d  marginal us/slab/E %6.3f -> %6.3f (%4.2fx)  |dP| %.1e%s"
          % (d, 1e6*mo, 1e6*mn, mo/mn, dP, tail))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
