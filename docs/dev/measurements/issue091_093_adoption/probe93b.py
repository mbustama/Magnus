import json, time, warnings, numpy as np
warnings.simplefilter('ignore')
import magnus.oscprob as op, magnus.globaldefs as gd, magnus.earth as earth
R = json.load(open('/home/mbustamante/Research/magnus/notebooks/magnus_own_reference.json')); p = R['oscillation_parameters']
L = earth.distance_traveled_inside_earth(-0.9)*gd.CONV_KM_TO_INV_EV
c = dict(costhz=-0.9, L=L, s12=np.sqrt(p['s12sq']), s23=np.sqrt(p['s23sq']), s13=np.sqrt(p['s13sq']), D21=p['dmsq21_ev2'], D31=p['dmsq31_ev2'], electron_fraction=0.5, ratio_number_neutrons_to_protons=1.0, nu_i=gd.NUMU, strategy='magnus')
E = np.array(R['three_flavor']['energy_gev'])*gd.UNIT_GEV
dcp0 = np.radians(p['dcp_deg'])
try:
    full = np.asarray(op.osc_prob_3nu_earth(E, dCP=dcp0, nu_f=None, rtol=1e-3, atol=1e-5, **c)); print('nu_f=None ->', full.shape)
except Exception as exc:
    print('nu_f=None ->', str(exc)[:120])
def scan(fn, steps=25):
    t = time.perf_counter()
    for k in range(steps): fn(dcp0 + 0.2*k/steps)
    return (time.perf_counter() - t)/(steps*len(E))*1e6
for rtol, stored in ((1e-3, 226), (1e-8, 2950)):
    kw = dict(rtol=rtol, atol=rtol*1e-2)
    one = lambda d: op.osc_prob_3nu_earth(E, dCP=d, nu_f=gd.NUMU, **kw, **c)
    three = lambda d: [op.osc_prob_3nu_earth(E, dCP=d, nu_f=f, **kw, **c) for f in (gd.NUE, gd.NUMU, gd.NUTAU)]
    one(dcp0); three(dcp0)
    print('rtol %g: stored %d us/prob | one channel %.0f | three channels %.0f' % (rtol, stored, min(scan(one) for _ in range(2)), min(scan(three) for _ in range(2))), flush=True)
