import numpy as np; import magnus.oscprob as oscprob; import magnus.globaldefs as gd
from magnus.earth import (
 distance_traveled_inside_earth as chord,
 earth_radial_distance_from_depth as radius,
 density_matter_func_prem as prem,
 Y_E_MANTLE_PREM)

osc = gd.load_nufit_params('NuFIT 6.1')
c, depth = -0.552, 1000.0   # mantle chord
geo = dict(source_depth=depth,
           detector_depth=1.4)          # km
L = chord(c, **geo)                     # km
def rho(l):   # PREM along the chord, g/cm^3
    r = radius(c, l/gd.UNIT_KM, **geo)
    return prem(r,
                density_matter_ocean=2.65)
E = 2.5*gd.UNIT_MEV
P_m = oscprob.osc_prob_matter_std_potential(
 3, rho, E, L*gd.UNIT_KM, osc, average=True,
 nubar=True, nu_i=gd.NUE, nu_f=gd.NUE,
 electron_fraction=Y_E_MANTLE_PREM,
 density_matter_is_in_g_per_cm3=True)
# 0.551
P_v = oscprob.osc_prob_3nu_vacuum(E,
 L*gd.UNIT_KM, average=True, nubar=True,
 nu_i=gd.NUE, nu_f=gd.NUE)
# 0.548
print('L', L, 'P_m', P_m, 'P_v', P_v)
import magnus.earth as earth
for d,c in [(20.0,-0.272),(1000.0,-0.552),(2800.0,-0.872)]:
    Lc = earth.distance_traveled_inside_earth(c, source_depth=d, detector_depth=1.4)
    print(d,c,Lc, len(earth.prem_layer_edges_along_chord(c, source_depth=d, detector_depth=1.4)))
print('local', earth.distance_traveled_inside_earth(0.078, source_depth=1.4, detector_depth=10.0), len(earth.prem_layer_edges_along_chord(0.078, source_depth=1.4, detector_depth=10.0)))
