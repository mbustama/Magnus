import json, sys, warnings, numpy as np
warnings.simplefilter('ignore')
import magnus.oscprob as op, magnus.globaldefs as gd
print('magnus from', op.__file__)
R = json.load(open('/home/mbustamante/Research/magnus/notebooks/magnus_own_reference.json'))
p = R['oscillation_parameters']
import magnus.earth as earth
LCH = earth.distance_traveled_inside_earth(-0.9)*gd.CONV_KM_TO_INV_EV
base = dict(costhz=-0.9, L=LCH, s12=np.sqrt(p['s12sq']), s23=np.sqrt(p['s23sq']), s13=np.sqrt(p['s13sq']), dCP=np.radians(p['dcp_deg']),
            D21=p['dmsq21_ev2'], D31=p['dmsq31_ev2'], electron_fraction=0.5, ratio_number_neutrons_to_protons=1.0, nu_i=gd.NUMU)
ster = dict(s14=np.sqrt(p['sinsq_th14']), s24=np.sqrt(p['sinsq_th24']), s34=np.sin(p['th34']), D41=p['dmsq41_ev2'])
t3, t4 = R['three_flavor'], R['sterile_3plus1']
E3, E4 = np.array(t3['energy_gev'])*gd.UNIT_GEV, np.array(t4['energy_gev'])*gd.UNIT_GEV
row3 = np.array(t3['numu_row']); ref4 = np.array(t4['p_numu_numu'])
print('ref shapes', row3.shape, ref4.shape)
def err3(**k):
    P = np.array([np.asarray(op.osc_prob_3nu_earth(E3, **base, nu_f=f, **k), float).ravel() for f in (gd.NUE, gd.NUMU, gd.NUTAU)]).T
    rr = row3 if row3.shape == P.shape else row3.T
    return np.max(np.abs(P - rr)), np.max(np.abs(P[:, 1] - rr[:, 1]))
def err4(**k):
    P = np.asarray(op.osc_prob_4nu_earth(E4, **base, **ster, nu_f=gd.NUMU, **k), float).ravel()
    return np.max(np.abs(P - ref4))
targets3 = {('tol', 1e-3): 3.111e-05, ('tol', 1e-8): 3.183e-10, ('slab', 16): 1.845e-05}
targets4 = {('tol', 1e-4): 1.879e-06, ('slab', 128): 6.780e-07}
for strat in ('magnus', 'auto'):
    for aname, af in (('atol=rtol', 1.0), ('atol=rtol/100', 1e-2), ('atol=0', 0.0)):
        out = []
        for (kind, v), tgt in targets3.items():
            k = dict(rtol=v, atol=v*af) if kind == 'tol' else dict(n_slabs=v, rtol=None, atol=None)
            if kind == 'slab' and aname != 'atol=rtol': continue
            try: row, mm = err3(strategy=strat, **k); out.append('3nu %s %g: row %.3e mm %.3e (target %.3e)' % (kind, v, row, mm, tgt))
            except Exception as exc: out.append('3nu %s %g: %s' % (kind, v, str(exc)[:80]))
        for (kind, v), tgt in targets4.items():
            k = dict(rtol=v, atol=v*af) if kind == 'tol' else dict(n_slabs=v, rtol=None, atol=None)
            if kind == 'slab' and aname != 'atol=rtol': continue
            try: out.append('3+1 %s %g: %.3e (target %.3e)' % (kind, v, err4(strategy=strat, **k), tgt))
            except Exception as exc: out.append('3+1 %s %g: %s' % (kind, v, str(exc)[:80]))
        print('--', strat, aname); print('\n'.join('   ' + o for o in out), flush=True)
