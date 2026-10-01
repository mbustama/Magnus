import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
E, costhz = 10.0*gd.UNIT_GEV, -0.8
kw = dict(nu_i=gd.NUMU, nu_f=gd.NUMU)
L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
a = oscprob.osc_prob_3nu_earth(E, costhz=costhz, L=L, **kw)
b = oscprob.osc_prob_3nu_earth(E, costhz=costhz, detector_depth=2.0*gd.UNIT_KM, **kw)
print('-0.8 surface %.6f  2km %.6f diff %.2e'%(a,b,b-a))
E, costhz = 10.0*gd.UNIT_GEV, 0.0
L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
a = oscprob.osc_prob_3nu_earth(E, costhz=costhz, L=L, **kw)
b = oscprob.osc_prob_3nu_earth(E, costhz=costhz, detector_depth=2.0*gd.UNIT_KM, **kw)
print('0 surface %.6f  2km %.6f diff %.2e'%(a,b,a-b))
R=6371.0; print('horizontal chord at 2 km depth', 2*np.sqrt(R**2-(R-2)**2))
L = earth.distance_traveled_inside_earth(-0.02)*gd.UNIT_KM
kw = dict(costhz=-0.02, L=L, nu_i=gd.NUMU, nu_f=gd.NUE)
E = 1.0*gd.UNIT_GEV
p0 = oscprob.osc_prob_3nu_earth(E, **kw)
p1 = oscprob.osc_prob_3nu_earth(E, **kw, density_matter_ocean=0.92)
p2 = oscprob.osc_prob_3nu_earth(E, **kw, density_matter_ocean=2.65, electron_fraction_ocean=0.4952)
print('ocean %.6f ice %.6f rock %.6f; ice rel %.3f%% rock rel %.2f%%'%(p0,p1,p2,100*(p1-p0)/p0,100*(p2-p0)/p0))
print('chord -0.02 km', L/gd.UNIT_KM, 'max depth', R-np.sqrt(R**2-(L/gd.UNIT_KM/2)**2))
