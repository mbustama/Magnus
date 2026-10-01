import numpy as np, magnus.globaldefs as gd, magnus.solarmodels as solarmodels, magnus.matter as matter, magnus.hamiltonians as ham
osc = gd.load_nufit_params('NuFIT 6.1')
R = gd.SUN_RADIUS*gd.UNIT_KM
ne_sun = solarmodels.electron_density_profile('B16-GS98')
PER_NE = matter.VCC_func(0.0, lambda l: 1.0)
def chord(br):
    b=br*R; hl=np.sqrt(max(R**2-b**2,0)); return (lambda l: ne_sun(np.sqrt((np.asarray(l,dtype=float)-hl)**2+b**2))), hl
def matter_phase(br, n=4001):
    ne_b, hl = chord(br); l=np.linspace(0,2*hl,n); v=PER_NE*ne_b(l); return float(np.sum(0.5*(v[1:]+v[:-1])*np.diff(l)))
print('phase diameter %.1f' % matter_phase(0.0))
bb=np.linspace(0,0.5,4001); ph=np.array([matter_phase(b) for b in bb])
ring=2*np.pi/np.maximum(np.abs(np.gradient(ph,bb)),1e-300)
for b in (0.11,0.3,0.5): print('ring at %.2f: %.3g R' % (b, np.interp(b,bb,ring)))
step=np.minimum(0.0005, ring/5)
count=np.concatenate([[0.0],np.cumsum(0.5*(1/step[1:]+1/step[:-1])*np.diff(bb))])
core=np.interp(np.arange(0.0,count[-1]),count,bb)
grid=np.linspace(0,1,2001)
n=len(core)+np.sum(grid>=0.5); print('impact params', n, 'inside 0.5:', len(core))
sunarc=959.6  # arcsec, solar angular radius approx
print('ring arcsec: %.2f .. %.2f' % (ring[np.argmin(abs(bb-0.11))]*sunarc, ring[-1]*sunarc), 'min over grid', ring[bb>0.01].min()*sunarc)
# nu1-nu2 vacuum phase across Sun at 100 GeV
print('D21 phase 100 GeV across diameter', osc['D21']*2*R/(2*100*gd.UNIT_GEV), ' D21 L/(4E):', osc['D21']*2*R/(4*100*gd.UNIT_GEV))
U=np.linalg.eigh(np.array(ham.hamiltonian_3nu_vacuum(10*gd.UNIT_MEV, **osc),dtype=complex))[1]
print('sum|Ue|^4', np.sum(abs(U[0])**4))
