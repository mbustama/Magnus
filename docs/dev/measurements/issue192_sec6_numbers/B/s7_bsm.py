import numpy as np, warnings, time
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
costhz = -0.9
L = earth.distance_traveled_inside_earth(costhz) *gd.UNIT_KM
kw = dict(costhz=costhz, L=L, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=1.0e-8, atol=1.0e-10)
E = np.logspace(0, np.log10(40), 260)*gd.UNIT_GEV
E_TEV = np.logspace(0, np.log10(30), 200)*gd.UNIT_TEV
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
liv = dict(b1=0.0, b2=0.0, b3=np.pi/L, Lambda=1.0, sxi12=0.0, sxi23=1.0/np.sqrt(2.0), sxi13=0.0, dxiCP=0.0, n_liv=0)
print('b3 eV', np.pi/L)
s4 = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, D41=1.0)
s5 = dict(**s4, s15=np.sqrt(0.06), s25=np.sqrt(0.06), s35=0.0, D51=1.7)
t=time.time()
P_std = oscprob.osc_prob_3nu_earth(E, **kw)
P_nsi = oscprob.osc_prob_3nu_earth_nsi(E, **eps, **kw)
P_liv = oscprob.osc_prob_3nu_earth_liv(E, **liv, **kw)
P_std_tev = oscprob.osc_prob_3nu_earth(E_TEV, **kw)
P_3p1 = oscprob.osc_prob_4nu_earth(E_TEV, **s4, **kw)
P_3p2 = oscprob.osc_prob_5nu_earth(E_TEV, **s5, **kw)
print('bsm listing ran', time.time()-t, np.shape(P_liv))
# timing one sterile scan point
E5 = 5.0*gd.UNIT_TEV
t=time.time()
for i in range(20):
    oscprob.osc_prob_4nu_earth(E5, s14=0.3, s24=np.sqrt(0.1), s34=0.0, D41=1.0, **kw)
print('per call', (time.time()-t)/20)
# oscillogram facts
print('chord -0.05, -1', earth.distance_traveled_inside_earth(-0.05), earth.distance_traveled_inside_earth(-1.0))
print('graze cos', -np.sqrt(1-(3480/6371)**2))
