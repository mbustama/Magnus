# 2x2 factorial: composition kernel x commutator kernel, PREM chord, order 4,
# fixed 544->1088, marginal us/slab/E, control interleaved.
import sys, json, warnings
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.magnus as mg
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc

warnings.simplefilter('ignore')
KER_P = mg._ordered_product_batched_kernel
KER_C = mg._commutator_batched_kernel
assert KER_P is not None and KER_C is not None
refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
ch = pcc.chord()
E = np.asarray(refs['energy_ev'], dtype=float)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
timed = gpb.timed
c0 = timed(gpb.control)

def prem_call(d, n_slabs):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        d, lambda x: ch['vcc'](x)/per_ne, E, ch['baseline'], gpb.osc_params(d),
        L0=0.0, density_is_of_number_of_electrons=True, rtol=None, atol=None,
        n_slabs=n_slabs, magnus_exp_order=4, strategy='magnus',
        t_breakpoints=ch['edges'][1:-1], validate_input=False))

NS = (544, 1088)
for d in (2, 3):
    res = {}
    for p_on in (0, 1):
        for c_on in (0, 1):
            mg._ordered_product_batched_kernel = KER_P if p_on else None
            mg._commutator_batched_kernel = KER_C if c_on else None
            for n in NS:
                res[(p_on, c_on, n)] = timed(lambda: prem_call(d, n))
    mg._ordered_product_batched_kernel = KER_P
    mg._commutator_batched_kernel = KER_C
    marg = {}
    for p_on in (0, 1):
        for c_on in (0, 1):
            marg[(p_on, c_on)] = (res[(p_on, c_on, 1088)] - res[(p_on, c_on, 544)])/(1104 - 560)/12
    base = marg[(0, 0)]
    print("d=%d marginal us/slab/E: none %.3f | compose-only %.3f (%.2fx) | comm-only %.3f (%.2fx) | both %.3f (%.2fx)"
          % (d, 1e6*base, 1e6*marg[(1,0)], base/marg[(1,0)],
             1e6*marg[(0,1)], base/marg[(0,1)], 1e6*marg[(1,1)], base/marg[(1,1)]))
    print("     additive check: c=%.3f k=%.3f  c+k+rest=%.3f (rest=%.3f)  product-of-speedups %.2fx vs actual %.2fx"
          % (1e6*(base - marg[(1,0)]), 1e6*(base - marg[(0,1)]),
             1e6*base, 1e6*(marg[(1,1)]),
             (base/marg[(1,0)])*(base/marg[(0,1)]), base/marg[(1,1)]))
print("control ratio: %.3f" % (timed(gpb.control)/c0))
