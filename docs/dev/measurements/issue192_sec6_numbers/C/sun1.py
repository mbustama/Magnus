import os
import numpy as np, warnings
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
import magnus.solarmodels as solarmodels
model = 'BS05-AGS-OP'
r = solarmodels.load_solar_model(model)['r_over_r_sun']
ne_sun = solarmodels.electron_density_profile(model)
ne = ne_sun(r*gd.SUN_RADIUS*gd.UNIT_KM)
E = np.logspace(-1.0, np.log10(20.0), 90)*gd.UNIT_MEV
osc = gd.load_nufit_params('NuFIT 6.1')
R = solarmodels.table_edge(model)
run = dict(nu_i=gd.NUE, nu_f=gd.NUE, average=True, density_profile=model)
P3 = oscprob.osc_prob_3nu_sun(E, R, 0.0, **osc, **run)
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
P3_nsi = oscprob.osc_prob_3nu_sun_nsi(E, R, 0.0, **osc, **eps, **run)
ster = dict(s14=0.10**0.5, s24=0.10**0.5, D41=1.0)
P4 = oscprob.osc_prob_4nu_sun(E, R, 0.0, **osc, **ster, **run)
P5 = oscprob.osc_prob_5nu_sun(E, R, 0.0, **osc, **ster, s15=0.06**0.5, s25=0.06**0.5, D51=1.7, **run)
P3=np.asarray(P3)
print('P3 at 0.1, 20:', P3[0], P3[-1])
print('P3 at 1 MeV interp:', np.interp(0.0, np.log(E/gd.UNIT_MEV), P3))
print('P3 shape', P3.shape, 'P4 last', np.asarray(P4)[[0,-1]], 'P5', np.asarray(P5)[[0,-1]], 'nsi', np.asarray(P3_nsi)[[0,-1]])
print('min over E: P4<P3?', np.all(np.asarray(P4)<P3), 'P5<P4?', np.all(np.asarray(P5)<np.asarray(P4)))
# |Ue2|^2 and vacuum
import magnus.hamiltonians as ham
U = ham.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP'])
print('|Ue2|^2', abs(U[0,1])**2, 'vac Pee', sum(abs(U[0,i])**4 for i in range(3)), 'c13^4/2', (1-osc['s13']**2)**2/2)
for e in (1.0,):
    print('P3 at 1 MeV direct', oscprob.osc_prob_3nu_sun(e*gd.UNIT_MEV, R, 0.0, **osc, **run))
np.save('P3.npy', P3); np.save('E.npy', E)
# pysnip: instantaneous
dL = np.array([0, 50, 100])*gd.UNIT_KM
info={}
P = oscprob.osc_prob_3nu_sun(5.0*gd.UNIT_MEV, R - dL, 0.0, **osc, nu_i=gd.NUE, nu_f=gd.NUE, density_profile=model, strategy_info=info)
print('inst P', np.asarray(P))
s = info.get('sampling')
print('sampling', s if s is None else {k: s[k] for k in ('cycles_per_step','aliased','nyquist_points')})
# production point
P = oscprob.osc_prob_3nu_sun(5.0*gd.UNIT_MEV, R, 0.05*gd.SUN_RADIUS*gd.UNIT_KM, **osc, **run)
P0 = oscprob.osc_prob_3nu_sun(5.0*gd.UNIT_MEV, R, 0.0, **osc, **run)
print('5MeV r0=0.05', P, 'center', P0, 'diff', P-P0)
P = oscprob.osc_prob_3nu_sun(20.0*gd.UNIT_MEV, R, 0.05*gd.SUN_RADIUS*gd.UNIT_KM, **osc, **run)
P0 = oscprob.osc_prob_3nu_sun(20.0*gd.UNIT_MEV, R, 0.0, **osc, **run)
print('20MeV r0=0.05', P, 'center', P0)
# P2 approx
c13sq = 1.0 - osc['s13']**2
P2 = oscprob.osc_prob_matter_std_potential(2, lambda l: c13sq*ne_sun(l), E, R, dict(sth=osc['s12'], Dm2=osc['D21']), L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, average=True, density_is_of_number_of_electrons=True)
P_approx = osc['s13']**4 + c13sq**2*np.asarray(P2)
rel = (P_approx - P3)/P3
np.save('rel_bs05.npy', rel)
print('rel err 0.1,1,20:', rel[0], np.interp(0.0, np.log(E/gd.UNIT_MEV), rel), rel[-1], 'all >0', np.all(rel>0))
