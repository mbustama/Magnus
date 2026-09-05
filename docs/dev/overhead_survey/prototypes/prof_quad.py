import sys, json, warnings, cProfile, pstats, io
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.oscprob as oscprob
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc

warnings.simplefilter('ignore')
refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
ch = pcc.chord()
E = np.asarray(refs['energy_ev'], dtype=float)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)

def prem_call(d, n_slabs, order, im, ntps):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        d, lambda x: ch['vcc'](x)/per_ne, E, ch['baseline'], gpb.osc_params(d),
        L0=0.0, density_is_of_number_of_electrons=True, rtol=None, atol=None,
        n_slabs=n_slabs, magnus_exp_order=order, integration_method=im,
        n_tpts_per_slab=ntps, strategy='magnus',
        t_breakpoints=ch['edges'][1:-1], validate_input=False))

d, n, order, im, ntps, reps = (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
                               sys.argv[4], int(sys.argv[5]), int(sys.argv[6]))
prem_call(d, n, order, im, ntps); prem_call(d, n, order, im, ntps)
pr = cProfile.Profile(); pr.enable()
for _ in range(reps):
    prem_call(d, n, order, im, ntps)
pr.disable()
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(28)
print(s.getvalue())
