import numpy as np, warnings
import magnus.earth as earth
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
costhz = -0.9
L = earth.distance_traveled_inside_earth(costhz)
print('L', repr(L))
print('closest approach', 6371*np.sqrt(1-costhz**2))
edges = earth.prem_layer_edges_along_chord(costhz)
print('len edges', len(edges), edges)
P = oscprob.osc_prob_3nu_earth(10.0*gd.UNIT_GEV, costhz=costhz, L=L*gd.UNIT_KM)
print('P numu nue 10GeV', repr(P[gd.NUMU][gd.NUE]))
kw = dict(costhz=costhz, L=L*gd.UNIT_KM, nu_i=gd.NUMU, nu_f=gd.NUE)
E = 5.0*gd.UNIT_GEV
a = oscprob.osc_prob_3nu_earth(E, **kw); b = oscprob.osc_prob_3nu_earth(E, **kw, electron_fraction=0.5)
print('5GeV layered', repr(a), 'Ye=0.5', repr(b), 'diff', a-b)
# PREM facts
print('Ye consts', earth.Y_E_CORE_PREM, earth.Y_E_MANTLE_PREM, earth.Y_E_CRUST_PREM, earth.Y_E_OCEAN_PREM)
print([n for n in dir(earth) if 'prem' in n.lower() or 'radius' in n.lower() or 'RADIUS' in n])
