import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
fnal = ((41, 49, 55), (-88, 15, 26))
hs = ((44, 21, 5.76), (-103, 45, 4.68))
P = oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, loc_ini=fnal, loc_fin=hs, nu_i=gd.NUMU, nu_f=gd.NUE)
print('coords', repr(P))
P = oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, loc_ini='fermilab', loc_fin='homestake', nu_i=gd.NUMU, nu_f=gd.NUE)
print('names', repr(P))
E = np.logspace(np.log10(0.3), 1.0, 400)*gd.UNIT_GEV
sites = ('snolab', 'homestake', 'cern', 'south_pole')
P = {site: oscprob.osc_prob_3nu_earth(E, loc_ini='fermilab', loc_fin=site, nu_i=gd.NUMU, nu_f=gd.NUE, rtol=1.0e-8, atol=1.0e-10) for site in sites}
a = earth.loc_coords_dms['fermilab']
b = earth.loc_coords_dms['homestake']
L_km = earth.chord_length_inside_earth(a['lat'], a['lon'], b['lat'], b['lon'])
print('%.0f km' % L_km, repr(L_km))
print('table count', len(earth.loc_coords_dms), sorted(earth.loc_coords_dms))
R=6371.0
Eg=E/gd.UNIT_GEV
for s in sites:
    b = earth.loc_coords_dms[s]
    Lk = earth.chord_length_inside_earth(a['lat'], a['lon'], b['lat'], b['lon'])
    depth = R - np.sqrt(R**2 - (Lk/2)**2)
    p=np.asarray(P[s]); m=Eg>0.8; k=np.argmax(np.where(m,p,-1))
    print(s, 'L=%.1f km depth=%.1f km  max P(E>0.8)=%.4f at %.2f GeV; overall max %.4f'%(Lk, depth, p[k], Eg[k], p.max()))
for s in sites:
    p=np.asarray(P[s]); idx=[i for i in range(1,len(p)-1) if p[i]>p[i-1] and p[i]>p[i+1]]
    print(s, 'local maxima (highest E last):', [(round(Eg[i],2), round(p[i],3)) for i in idx][-4:])
