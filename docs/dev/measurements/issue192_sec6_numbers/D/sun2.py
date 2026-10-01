import numpy as np, magnus.adiabatic as adiabatic, magnus.hamiltonians as ham, magnus.matter as matter, magnus.globaldefs as gd, magnus.solarmodels as solarmodels
osc = gd.load_nufit_params('NuFIT 6.1')
R = gd.SUN_RADIUS*gd.UNIT_KM
ne_sun = solarmodels.electron_density_profile('B16-GS98')
for bfrac in (0.0, 0.3):
  b = bfrac*R; half = np.sqrt(R**2-b**2)
  for Egev in (100.0, 1e3, 1e5, 1e6):
    E = Egev*gd.UNIT_GEV
    hv = ham.hamiltonian_3nu_vacuum(E, **osc)
    def ne(l): return ne_sun(np.sqrt((l-half)**2+b**2))
    def H(l):
        h = np.array(hv, dtype=complex); h[0,0] += matter.VCC_func(l, ne); return h
    info={}; probe={}
    adiabatic.find_nonadiabatic_windows(H, 0.0, 2*half, n_probe=20000, info=info, _probe_out=probe)
    ls, lam, W, dH = probe['ls'], probe['lam'], probe['W'], probe['dH']
    best=(0,None)
    for j in range(3):
        for k in range(j+1,3):
            num = np.abs(np.einsum('ni,nij,nj->n', W[:,:,j].conj(), dH, W[:,:,k]))
            g = num/(lam[:,k]-lam[:,j])**2
            i = g.argmax()
            if g[i]>best[0]: best=(g[i], (j,k,ls[i]))
    g,(j,k,l)=best
    r = np.sqrt((l-half)**2+b**2)/R
    print('b=%.1f E=%g GeV gamma_max(pkg)=%.4g  grid max %.4g pair %d-%d at r=%.4f R' % (bfrac, Egev, info['gamma_max'], g, j+1,k+1, r), flush=True)
