import numpy as np
import magnus.globaldefs as gd, magnus.matter as matter, magnus.hamiltonians as ham
from shockdef import *
osc = gd.load_nufit_params('NuFIT 6.1')
PER = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
ne = sn_shock_ne(1e-6)
rr = np.linspace(R0_KM, R1_KM, 700001)
n = ne(rr*KM)
def nres31(Em): return osc['D31']*(1-2*osc['s13']**2)/(2*Em*gd.UNIT_MEV)/PER
def nres21(Em): return osc['D21']*(1-2*osc['s12']**2)/(2*Em*gd.UNIT_MEV)/PER
print('21 res at 15 MeV / min ne on ray: %.2f (min ne %.3e at r=%.0f)' % (n.min()/nres21(15), n.min(), rr[n.argmin()]))
print('  with cos2th13 factor: %.2f' % (n.min()/(nres21(15)/(1-osc['s13']**2))))
eps=1.0
fs_in, fs_out = float(ne((R_FORWARD_KM-eps)*KM)), float(ne((R_FORWARD_KM+eps)*KM))
cd_in, cd_out = float(ne((R_CONTACT_KM-eps)*KM)), float(ne((R_CONTACT_KM+eps)*KM))
print('FS jump %.2f, CD jump %.2f' % (fs_in/fs_out, cd_in/cd_out))
E_of = lambda nn: osc['D31']*(1-2*osc['s13']**2)/(2*PER*nn)/gd.UNIT_MEV
print('FS: resonance between sides for E in %.2f .. %.2f MeV' % (E_of(fs_in), E_of(fs_out)))
print('CD: resonance between sides for E in %.2f .. %.2f MeV' % (E_of(cd_in), E_of(cd_out)))
# also with cos^2 theta13 not; plus no-shock at 80000
print('15 MeV res density between FS sides?', fs_out < nres31(15) < fs_in)
s2 = np.sin(2*np.arcsin(osc['s13']))
for Em in (6, 15, 18):
    print('Losc at H resonance %g MeV: %.1f km' % (Em, 4*np.pi*Em*gd.UNIT_MEV/(osc['D31']*s2)/KM))
# full 3nu eigenvalue gap at resonance-ish: min over r of lambda3-lambda2 gap in matter at 15MeV
U = np.asarray(ham.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP']))
HV = U@np.diag([0, osc['D21'], osc['D31']])@U.conj().T
def lam(E, nn): return np.linalg.eigvalsh(HV/(2*E)+np.diag([PER*nn,0,0]))
for Em in (6., 15., 18.):
    E=Em*gd.UNIT_MEV; ng = np.logspace(np.log10(nres31(Em))-0.3, np.log10(nres31(Em))+0.3, 2001)
    gap = np.array([np.diff(lam(E,x))[1] for x in ng]); 
    print('   3nu min gap osc length %g MeV: %.1f km' % (Em, 2*np.pi/gap.min()/KM))
# phase between fronts and after FS
def phase(E, a, b, N=200001):
    r = np.linspace(a, b, N); L = np.array([lam(E, x) for x in ne(r*KM)])
    return np.trapezoid(L[:,2]-L[:,0], r*KM), np.trapezoid(L[:,2]-L[:,1], r*KM), np.trapezoid(L[:,1]-L[:,0], r*KM)
E=15*gd.UNIT_MEV; dE=1e-3*E
for name,(a,b) in (('between',(R_CONTACT_KM+1, R_FORWARD_KM-1)), ('after FS',(R_FORWARD_KM+1, R1_KM))):
    p = np.array(phase(E,a,b,40001)); p2=np.array(phase(E+dE,a,b,40001))
    print(name, 'phases /2pi', p/(2*np.pi), 'energy period keV', 2*np.pi/np.abs((p2-p)/dE)/gd.UNIT_MEV*1e3)
