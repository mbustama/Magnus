import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.solarmodels as sm
import magnus.matter as matter, magnus.hamiltonians as ham
from scipy.integrate import quad
model='BS05-AGS-OP'
osc = gd.load_nufit_params('NuFIT 6.1')
R = sm.table_edge(model); RS = gd.SUN_RADIUS*gd.UNIT_KM
ne_sun = sm.electron_density_profile(model)
PER_NE = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
tab = sm.load_solar_model(model)
rr = tab['r_over_r_sun']; ne_t = ne_sun(rr*RS)
print('R/RS', R/RS, 'ne fall center/edge', ne_t[0]/ne_t[-1])
U = np.asarray(ham.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP']))
HV = U@np.diag([0, osc['D21'], osc['D31']])@U.conj().T
def lam(l, E):
    H = HV/(2*E) + np.diag([PER_NE*float(ne_sun(l)),0,0])
    return np.linalg.eigvalsh(H)
E5=5*gd.UNIT_MEV
ls = np.linspace(0, R, 20001)
L = np.array([lam(l,E5) for l in ls])
ph21 = np.trapezoid(L[:,1]-L[:,0], ls); ph31 = np.trapezoid(L[:,2]-L[:,0], ls)
print('phase21 %.3e phase31 %.3e' % (ph21, ph31), '2pi/ph21', 2*np.pi/ph21)
print('local osc length at surface (km): 21 %.1f, 31 %.2f' % (2*np.pi/(L[-1,1]-L[-1,0])/gd.UNIT_KM, 2*np.pi/(L[-1,2]-L[-1,0])/gd.UNIT_KM))
print('vacuum Losc 21 at 5 MeV %.1f km, 31: %.2f km' % (4*np.pi*E5/osc['D21']/gd.UNIT_KM, 4*np.pi*E5/osc['D31']/gd.UNIT_KM))
print('local osc length at center 21: %.1f km' % (2*np.pi/(L[0,1]-L[0,0])/gd.UNIT_KM))
print('total phase 31 over 2pi -> nyquist ~', 2*ph31/(2*np.pi))
c2 = 1-2*osc['s12']**2
def nres(Emev): return osc['D21']*c2/(2*Emev*gd.UNIT_MEV)/PER_NE
rg = np.linspace(0, R, 200001); ng = ne_sun(rg)
for Em in (1,5,20):
    nr = nres(Em); i = np.argmax(ng < nr)
    print('E=%g res radius %.4f RS, central/res %.3f' % (Em, rg[i]/RS if ng[0]>nr else np.nan, ng[0]/nr))
# 8B band
frac=[]
for line in open('/home/user/Magnus/docs/dev/adversarial_batteries/bs2005agsopflux.csv'):
    f=line.split()
    if len(f)==13:
        try: frac.append([float(x) for x in f])
        except ValueError: pass
frac=np.array(frac); cum=np.cumsum(frac[:,6])/frac[:,6].sum()
band=[float(np.interp(q,cum,frac[:,0])) for q in (0.05,0.95)]
print('8B band', band)
ne_out = float(ne_sun(band[1]*RS))
print('E threshold where n_res = ne(outer band): %.3f MeV' % (osc['D21']*c2/(2*PER_NE*ne_out)/gd.UNIT_MEV))
# production curves
run = dict(nu_i=gd.NUE, nu_f=gd.NUE, average=True, density_profile=model)
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
R0 = np.linspace(0,0.8,41)*RS
for Em in (1.,5.,20.):
    E=Em*gd.UNIT_MEV
    std=np.array([float(oscprob.osc_prob_3nu_sun(E,R,float(l0),**osc,**run)) for l0 in R0])
    nsi=np.array([float(oscprob.osc_prob_3nu_sun_nsi(E,R,float(l0),**osc,**eps,**run)) for l0 in R0])
    print('E=%g: center %.4f, 0.5R %.4f, 0.8R %.4f; max|nsi-std| %.4f' % (Em, std[0], np.interp(0.5,R0/RS,std), std[-1], np.abs(nsi-std).max()))
    if Em==20.:
        nr=nres(20); rres=rg[np.argmax(ng<nr)]
        print('  20 MeV at resonance radius %.4f RS: %.4f' % (rres/RS, float(oscprob.osc_prob_3nu_sun(E,R,float(rres),**osc,**run))))
    if Em==5.:
        nr=nres(5); rres=rg[np.argmax(ng<nr)]
        print('  5 MeV at resonance radius %.4f RS: %.4f' % (rres/RS, float(oscprob.osc_prob_3nu_sun(E,R,float(rres),**osc,**run))))
# 2E VCC/D31 at center 20 MeV
print('2EV/D31 center 20MeV', 2*20*gd.UNIT_MEV*PER_NE*float(ne_sun(0))/osc['D31'], 's13^2', osc['s13']**2)
