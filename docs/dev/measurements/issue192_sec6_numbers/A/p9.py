import numpy as np, inspect
import magnus.globaldefs as gd, magnus.hamiltonians as ham, magnus.solarmodels as sm, magnus.matter as matter
osc = gd.load_nufit_params('NuFIT 6.1')
R = gd.SUN_RADIUS*gd.UNIT_KM; E=5*gd.UNIT_MEV
hv = ham.hamiltonian_3nu_vacuum(E, **osc)
print(inspect.signature(sm.electron_density_profile))
for prof in ('BS2005-AGS,OP','B16-GS98'):
    try: ne=sm.electron_density_profile(prof)
    except Exception as e: print(prof, e); continue
    ls=np.linspace(0,R,200001); ph=0; ph21=0
    lam=[]
    for l in ls[::100]:
        h=hv.copy(); h[0,0]+=matter.VCC_func(l, ne); w=np.linalg.eigvalsh(h); lam.append(w)
    lam=np.array(lam); x=ls[::100]
    P=np.trapezoid(lam[:,2]-lam[:,0], x); P2=np.trapezoid(lam[:,1]-lam[:,0], x)
    print(prof, 'Phi %.3g rad, cycles %.3g; slow pair %.3g rad'%(P, P/2/np.pi, P2))
