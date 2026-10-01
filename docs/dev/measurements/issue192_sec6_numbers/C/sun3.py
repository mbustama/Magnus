import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.solarmodels as sm
import magnus.matter as matter, magnus.hamiltonians as ham
osc = gd.load_nufit_params('NuFIT 6.1')
RS = gd.SUN_RADIUS*gd.UNIT_KM
E = np.logspace(-1.0, np.log10(20.0), 90)*gd.UNIT_MEV
Em = E/gd.UNIT_MEV
ne_bs = sm.electron_density_profile('BS05-AGS-OP'); ne_b16 = sm.electron_density_profile('B16-GS98')
Rbs = sm.table_edge('BS05-AGS-OP'); Rb16 = sm.table_edge('B16-GS98')
fit = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
# density ratios over tabulated BS05 rows
rb = sm.load_solar_model('B16-GS98')['r_over_r_sun']; rb = rb[rb <= Rbs/RS]
ratio = ne_b16(rb*RS)/ne_bs(rb*RS)
print('B16/BS05 inner half max|1-ratio| %.4f, outer half %.4f' % (np.abs(ratio[rb<=0.5]-1).max(), np.abs(ratio[rb>0.5]-1).max()))
print('exp/BS05 at center %.3f; at 0.98: %.2f; at 0.9 %.2f' % (fit(0.0)/ne_bs(0.0), fit(Rbs)/ne_bs(Rbs), fit(0.9*RS)/ne_bs(0.9*RS)))
kw = dict(nu_i=gd.NUE, nu_f=gd.NUE, average=True)
P = {}
P['BS05'] = np.asarray(oscprob.osc_prob_3nu_sun(E, Rbs, 0.0, **osc, density_profile='BS05-AGS-OP', **kw))
P['B16'] = np.asarray(oscprob.osc_prob_3nu_sun(E, Rb16, 0.0, **osc, density_profile='B16-GS98', **kw))
P['exp'] = np.asarray(oscprob.osc_prob_3nu_sun(E, RS, 0.0, **osc, **kw))
for k in ('B16','exp'):
    d = P[k]-P['BS05']; i = np.argmax(np.abs(d)); print('%s - BS05 largest %+.5f at %.2f MeV' % (k, d[i], Em[i]))
c13 = 1-osc['s13']**2
REL={}
for k, ne, Rm in (('BS05',ne_bs,Rbs),('B16',ne_b16,Rb16),('exp',fit,RS)):
    P2 = np.asarray(oscprob.osc_prob_matter_std_potential(2, lambda l, ne=ne: c13*ne(l), E, Rm, dict(sth=osc['s12'], Dm2=osc['D21']), L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, average=True, density_is_of_number_of_electrons=True))
    # textbook 2nu formula
    PER = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
    th = np.arcsin(osc['s12'])
    def c2m(l):
        x = 2*E*PER*c13*float(ne(l))/osc['D21']; return np.cos(np.arctan2(np.sin(2*th), np.cos(2*th)-x))
    tb = 0.5+0.5*c2m(0.0)*c2m(Rm)
    REL[k] = (osc['s13']**4 + c13**2*P2 - P[k])/P[k]
    r = REL[k]
    print('%s: P2 vs textbook max %.1e; rel err 0.1 %.2e, 1 %.2e, 20 %.4f, all>0 %s; first E where >1e-2: %s' % (k, np.abs(P2-tb).max(), r[0], r[np.argmin(abs(Em-1))], r[-1], np.all(r>0), Em[np.argmax(r>1e-2)] if np.any(r>1e-2) else None))
for k in ('BS05','B16'):
    q = REL['exp']/REL[k]; print('exp/%s ratio of rel errs %.2f .. %.2f' % (k, q.min(), q.max()))
# adiabatic limit reference for P3 (flavor start), table-edge readout
PER = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
U = np.asarray(ham.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP']))
HV = U@np.diag([0, osc['D21'], osc['D31']])@U.conj().T
ref=[]
for e in E:
    _,u0=np.linalg.eigh(HV/(2*e)+np.diag([PER*float(ne_bs(0.0)),0,0]))
    _,u1=np.linalg.eigh(HV/(2*e)+np.diag([PER*float(ne_bs(Rbs)),0,0]))
    ref.append(np.sum(abs(u0[0])**2*abs(u1[0])**2))
print('P3 vs eq averaged_varying max diff %.1e' % np.abs(P['BS05']-np.array(ref)).max())
# 12 models center to surface, 40 energies (notebook 13 setup) and also 90 energies to table edge
OSCNO = gd.load_nufit_params('NuFIT 6.1','NO')
EB = np.logspace(np.log10(0.1), np.log10(20.0), 40)*gd.UNIT_MEV
tables = sm.available_solar_models(); print(tables)
import sys; sys.exit() if 0 else None
A = np.array([np.asarray(oscprob.osc_prob_3nu_sun(EB, RS, 0.0, **OSCNO, density_profile=m, **kw)) for m in tables])
sp = A.max(0)-A.min(0); i=np.argmax(sp)
print('12-model spread (to surface, 40 E): %.2e at %.2f MeV' % (sp[i], EB[i]/gd.UNIT_MEV))
b23 = np.array([A[j] for j,m in enumerate(tables) if m.startswith('B23-')])
print('B23 spread %.2e' % (b23.max(0)-b23.min(0)).max())
A2 = np.array([np.asarray(oscprob.osc_prob_3nu_sun(E, sm.table_edge(m), 0.0, **osc, density_profile=m, **kw)) for m in tables])
sp2 = A2.max(0)-A2.min(0); print('12-model spread (to table edges, 90 E): %.2e' % sp2.max())
for k in ('BS05','B16'):
    q = REL['exp']/REL[k]; print(k, 'min at %.3f MeV, max at %.3f MeV' % (Em[np.argmin(q)], Em[np.argmax(q)])); print(np.round(q[::6],3))
np.save('REL.npy', np.array([REL['BS05'],REL['B16'],REL['exp']]))
