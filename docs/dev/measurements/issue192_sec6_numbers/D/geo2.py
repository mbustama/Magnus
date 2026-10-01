import numpy as np, magnus.globaldefs as gd, magnus.matter as matter, magnus.earth as earth, magnus.hamiltonians as ham
osc = gd.load_nufit_params('NuFIT 6.1'); print({k:osc[k] for k in osc})
dinv = (1/1.8-1/3.3)/gd.UNIT_MEV
for name,L in [('local',100.21),('far',3395.4),('mantle',4314.9),('core',7295.1)]:
    Lk=L*gd.UNIT_KM
    print(name, 'slow cycles %.2f fast cycles %.1f' % (osc['D21']*Lk*dinv/(4*np.pi), osc['D31']*Lk*dinv/(4*np.pi)))
E=2.5*gd.UNIT_MEV
print('Losc slow %.2f km fast %.3f km' % (4*np.pi*E/osc['D21']/gd.UNIT_KM, 4*np.pi*E/osc['D31']/gd.UNIT_KM))
# potential vs vacuum slow splitting at ends; densities
def V(rho, ye):
    ne = matter.num_density_e_func(0.0, lambda _: rho, electron_fraction=ye, density_matter_is_in_g_per_cm3=True)
    return matter.VCC_func(0.0, lambda _: ne)
R=gd.EARTH_RADIUS
for lab, r in [('det 1.4km', R-1.4), ('src 10km', R-10),('src 20km',R-20),('src 1000', R-1000),('src 2800',R-2800)]:
    rho = earth.density_matter_func_prem(r, density_matter_ocean=2.65)
    for Emev in (1.8, 3.3):
        print(lab, 'rho %.3f'%rho, Emev, 'V/(D21/2E) = %.4f' % (V(rho, 0.5)/(osc['D21']/(2*Emev*gd.UNIT_MEV))))
# CMB jump effect on mixing (2-flavour 12 sector, antineutrino: V -> -V)
def s2m(rho, Emev, ye=0.5):
    v = -V(rho, ye); d = osc['D21']/(2*Emev*gd.UNIT_MEV); s12=osc['s12']; c2=1-2*s12**2; s2=2*s12*np.sqrt(1-s12**2)
    c2m = (d*c2 - v)/np.hypot(d*c2 - v, d*s2); return (1-c2m)/2
for Emev in (1.8,2.5,3.3):
    a=s2m(5.566,Emev); b=s2m(9.903,Emev)
    print('CMB', Emev, a, b, 'rel change %.4f' % ((b-a)/a), 'angle change %.4f' %((np.arcsin(np.sqrt(b))-np.arcsin(np.sqrt(a)))/np.arcsin(np.sqrt(a))))
print(earth.density_matter_func_prem(3480.0-1e-6), earth.density_matter_func_prem(3480.0+1e-6), earth.Y_E_MANTLE_PREM)
