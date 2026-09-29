import json, time, warnings, numpy as np
warnings.simplefilter('ignore')
import magnus.oscprob as op, magnus.globaldefs as gd, magnus.earth as earth
R = json.load(open('/home/mbustamante/Research/magnus/notebooks/magnus_own_reference.json')); p = R['oscillation_parameters']
L = earth.distance_traveled_inside_earth(-0.9)*gd.CONV_KM_TO_INV_EV
c = dict(costhz=-0.9, L=L, s12=np.sqrt(p['s12sq']), s23=np.sqrt(p['s23sq']), s13=np.sqrt(p['s13sq']), D21=p['dmsq21_ev2'], D31=p['dmsq31_ev2'], electron_fraction=0.5, ratio_number_neutrons_to_protons=1.0, strategy='magnus')
st = dict(s14=np.sqrt(p['sinsq_th14']), s24=np.sqrt(p['sinsq_th24']), s34=np.sin(p['th34']), D41=p['dmsq41_ev2'])
dcp0 = np.radians(p['dcp_deg'])
E3 = np.array(R['three_flavor']['energy_gev'])*gd.UNIT_GEV; row = np.array(R['three_flavor']['numu_row'])
E4 = np.array(R['sterile_3plus1']['energy_gev'])*gd.UNIT_GEV; mm4 = np.array(R['sterile_3plus1']['p_numu_numu'])
def scan(fn, n, steps=25):
    fn(dcp0); best = np.inf
    for _ in range(2):
        t = time.perf_counter()
        for k in range(steps): fn(dcp0 + 0.2*k/steps)
        best = min(best, (time.perf_counter() - t)/(steps*n)*1e6)
    return best
full3 = lambda d, kw: np.asarray(op.osc_prob_3nu_earth(E3, dCP=d, nu_i=None, nu_f=None, **kw, **c), float)
full4 = lambda d, kw: np.asarray(op.osc_prob_4nu_earth(E4, dCP=d, nu_i=None, nu_f=None, **kw, **c, **st), float)
one4 = lambda d, kw: np.asarray(op.osc_prob_4nu_earth(E4, dCP=d, nu_i=gd.NUMU, nu_f=gd.NUMU, **kw, **c, **st), float)
P = full3(dcp0, dict(rtol=1e-3, atol=1e-5)); print('full-matrix shape', P.shape)
for rtol, err, us in ((1e-3, 3.111e-05, 226), (1e-8, 3.183e-10, 2950)):
    kw = dict(rtol=rtol, atol=rtol*1e-2); P = full3(dcp0, kw)
    r = P[:, gd.NUMU, :] if P.shape[-2:] == (3, 3) else None
    e = np.max(np.abs(r - row)) if r is not None else float('nan')
    print('3nu rtol %g: full-matrix row err %.3e (stored %.3e) | %.0f us/prob (stored %d)' % (rtol, e, err, scan(lambda d: full3(d, kw), len(E3)), us), flush=True)
for rtol, err, us in ((1e-4, 1.879e-06, 611), (1e-10, 7.937e-12, 14003)):
    kw = dict(rtol=rtol, atol=rtol*1e-2); P = full4(dcp0, kw)
    e = np.max(np.abs(P[:, gd.NUMU, gd.NUMU] - mm4))
    print('3+1 rtol %g: full err %.3e, one-channel err %.3e (stored %.3e) | full %.0f, one %.0f us/prob (stored %d)' % (rtol, e, np.max(np.abs(one4(dcp0, kw).ravel() - mm4)), err, scan(lambda d: full4(d, kw), len(E4)), scan(lambda d: one4(d, kw), len(E4)), us), flush=True)
