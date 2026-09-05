# Profile the cumulative baseline scan: one energy, many baselines, PREM chord.
import sys, json, warnings, cProfile, pstats, io, time
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import numpy as np
import magnus.oscprob as oscprob
import magnus.hamiltonians.hamiltonians3nu as h3
import gen_profile_benchmarks as gpb
import prem_chord_common as pcc
import magnus.globaldefs as gd

warnings.simplefilter('ignore')
refs = json.loads(open('/home/mbustamante/Research/magnus/notebooks/prem_chord_reference.json').read())
ch = pcc.chord()
E0 = 5.0e9
params = gpb.osc_params(3)
per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
import magnus.matter as matter
h_vac = np.asarray(h3.hamiltonian_3nu_vacuum_energy_independent(
    params['s12'], params['s23'], params['s13'], params['dCP'],
    params['D21'], params['D31']))

L_out = np.linspace(ch['baseline']*0.01, ch['baseline'], 400)
def call():
    return np.asarray(oscprob.osc_prob_energy_baseline(
        lambda E: h_vac/(2.0*E), E0, L_out, L0=0.0,
        H_func_is_function_only_of_energy=True,
        magnus_exp_order=4, integration_method='gl',
        rtol=1e-6, atol=1e-8, cumulative=True, validate_input=False))
t0=time.perf_counter(); r = call(); t1=time.perf_counter()
print("call: %.1f ms for %d baselines (shape %s)" % (1e3*(t1-t0), len(L_out), np.shape(r)))
pr = cProfile.Profile(); pr.enable()
call(); call()
pr.disable()
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(24)
print(s.getvalue())
