import numpy as np, warnings
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
E = 1.0*gd.UNIT_GEV
L = 1300.0*gd.UNIT_KM
P = oscprob.osc_prob_3nu_vacuum(E, L)
print('1411 Pme = %.5f, Pmm = %.5f, Pmt = %.5f' % (P[gd.NUMU][gd.NUE], P[gd.NUMU][gd.NUMU], P[gd.NUMU][gd.NUTAU]), repr(P[1,0]),repr(P[1,1]),repr(P[1,2]))
P, U = oscprob.osc_prob_3nu_vacuum(E, L, return_evolution_operator=True)
print('1447', abs(U[gd.NUE][gd.NUMU])**2, '1448', abs(U[gd.NUMU][gd.NUE])**2)
P = oscprob.osc_prob_3nu_vacuum(E, L, default_osc_params_set_name='OSC_PARAMS_NU_FIT_5_2_SK_IO')
print('1503', P[1,0])
print('1506 count predefined', len(gd.OSC_PARAMS_PREDEFINED), [k for k in gd.OSC_PARAMS_PREDEFINED if 'DEFAULT' in k])
osc = gd.load_nufit_params('NuFIT 5.2', ordering='IO', category='without_SK')
P = oscprob.osc_prob_3nu_vacuum(E, L, **osc); print('1510 (no number)', P[1,0])
for c in ('LID','LEM'):
    osc = gd.load_nufit_params('NuFIT 2.1', category=c)
    P = oscprob.osc_prob_3nu_vacuum(E, L, **osc); print('1519', c, P[1,0])
for c in ('huber_fluxes_no_rsbl','free_fluxes_rsbl'):
    osc = gd.load_nufit_params('NuFIT 1.3', category=c)
    P = oscprob.osc_prob_3nu_vacuum(E, L, **osc); print('1525', c, P[1,0])
sets = gd.OSC_PARAMS_PREDEFINED
lines=[]
for name in sorted(sets):
    lines.append((name, sets[name]['description']))
print('1535 nlines', len(lines), lines[0], [l for l in lines if l[0]=='OSC_PARAMS_NU_FIT_1_0_IO'])
fits = gd.NUFIT_GLOBAL_FITS
n=0
for v in fits:
    n+=1
    if v=='NuFIT 2.1': print('1543', v, list(fits[v]['categories']))
print('1543 nlines', n)
P = oscprob.osc_prob_3nu_vacuum(E, L, dCP=0.0); print('1550', P[1,0])
P = oscprob.osc_prob_3nu_vacuum(E, L, default_osc_params_set_name='OSC_PARAMS_NU_FIT_6_1_SK_IO', dCP=0.0); print('1560', P[1,0])
P = oscprob.osc_prob_3nu_vacuum(E, L, default_osc_params_set_name='OSC_PARAMS_NU_FIT_6_1_SK_IO'); print('1561', P[1,0])
for conv in ('sin', 'sin2', 'rad', 'deg'):
    osc = gd.load_nufit_params('NuFIT 6.1', angles=conv)
    P = oscprob.osc_prob_3nu_vacuum(E, L, angles=conv, **osc)
    print('1577', conv, osc['s12'], P[1,0])
P = oscprob.osc_prob_2nu_vacuum(E, L, sth=0.5557, Dm2=7.537e-5); print('1591', P[0][1])
P = oscprob.osc_prob_2nu_vacuum(E, L, sth=33.759, Dm2=7.537e-5, angles='deg'); print('1595', P[0][1])
P = oscprob.osc_prob_4nu_vacuum(E, L, s12=0.3088, s23=0.4700, s13=0.02248, dCP=3.700, D21=7.537e-5, D31=2.511e-3,
 s14=0.10, s24=0.10, s34=0.0, d14=0.0, d24=0.0, D41=1.0, angles='sin2')
print('1612', P[1,0], P[1,1], P.shape)
P = oscprob.osc_prob_5nu_vacuum(E, L, s12=0.3088, s23=0.4700, s13=0.02248, dCP=3.700, D21=7.537e-5, D31=2.511e-3,
 s14=0.10, s24=0.10, s34=0.0, s15=0.06, s25=0.06, s35=0.0, d14=0.0, d15=0.0, d24=0.0, d35=0.0, D41=1.0, D51=1.7, angles='sin2')
print('1631', P[1,0], P[1,1], P.shape)
d = gd.OSC_PARAMS_DEFAULT if hasattr(gd,'OSC_PARAMS_DEFAULT') else None
print('default', d)
o=gd.load_nufit_params('NuFIT 6.1', angles='sin2'); print('sin2 defaults', o)
# batched
Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
print('1650', oscprob.osc_prob_3nu_vacuum(Es, L).shape)
Ls = np.linspace(0, 1300, 200)*gd.UNIT_KM
P = oscprob.osc_prob_3nu_vacuum(E, Ls); print('1656', P.shape, np.array_equal(P[-1], oscprob.osc_prob_3nu_vacuum(E, L)), np.abs(P[-1]-oscprob.osc_prob_3nu_vacuum(E, L)).max())
Es = np.array([1.0, 2.0, 5.0])*gd.UNIT_GEV
Ls = np.array([500, 1000, 1300])*gd.UNIT_KM
P = oscprob.osc_prob_3nu_vacuum(Es, Ls); print('1664', P.shape, P[0,1,0], P[1,1,0])
try:
    oscprob.osc_prob_3nu_vacuum(Es, Ls[:2]); print('no error')
except Exception as e: print('mismatch ->', type(e).__name__)
