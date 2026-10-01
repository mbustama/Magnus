import numpy as np
import magnus.oscprob as oscprob
import magnus.globaldefs as gd

# Twelve energies per cycle of the pair split by
# Dm31^2 on the longest chord, 7,300 km
E = 1/np.linspace(1/1.8, 1/3.3, 22500)*gd.UNIT_MEV

# Borexino, 1.4 km underground.  A production point
# enters as its depth and the cosine of its zenith
# angle at the detector; the call declares the PREM
# boundaries on the chord by itself.  The oscillation
# parameters are the defaults, NuFIT 6.1.
points = {'far crust': (20.0, -0.272),    # 3,400 km
          'mantle': (1000.0, -0.552),     # 4,300 km
          'core': (2800.0, -0.872)}       # 7,300 km
P = {}
for name, (depth, costhz) in points.items():
    P[name] = oscprob.osc_prob_3nu_earth(E,
     costhz=costhz, source_depth=depth*gd.UNIT_KM,
     detector_depth=1.4*gd.UNIT_KM, nubar=True,
     nu_i=gd.NUE, nu_f=gd.NUE,
     density_matter_ocean=2.65)

# The local point, 10 km deep and 100 km away, sits
# at the near end of its chord; the call runs every
# chord from the far end.  The survival probability
# is the same along a path and along its reverse, so
# run the segment from the detector, with the depths
# swapped and the zenith angle taken at the
# production point.
P['local'] = oscprob.osc_prob_3nu_earth(E,
 costhz=0.078, source_depth=1.4*gd.UNIT_KM,
 detector_depth=10.0*gd.UNIT_KM, nubar=True,
 nu_i=gd.NUE, nu_f=gd.NUE,
 density_matter_ocean=2.65)
import json
json.dump({k: np.asarray(v).tolist() for k,v in P.items()}, open('geo_energy.json','w'))
osc = gd.load_nufit_params('NuFIT 6.1')
Pv = float(oscprob.osc_prob_3nu_vacuum(E[0], 1.0e8*gd.UNIT_KM, average=True, nubar=True, nu_i=gd.NUE, nu_f=gd.NUE, **osc))
print('vac avg', Pv)
for k,v in P.items(): v=np.asarray(v); print(k, v.min(), v.max(), 'plain mean', v.mean())
