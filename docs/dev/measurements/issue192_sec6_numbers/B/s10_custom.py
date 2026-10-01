import numpy as np, warnings
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.matter as matter
import magnus.globaldefs as gd
EDGES = [1221.5, 3480.0, 6346.6]       # km
RHO = [12.894, 10.901, 4.476, 2.520]   # g/cm^3
YE = [0.4656, 0.4656, 0.4957, 0.4952]  # PREM's own
costhz = -0.9
r_of = earth.earth_radial_distance_from_depth
L = earth.distance_traveled_inside_earth(costhz) *gd.UNIT_KM
osc = gd.load_nufit_params('NuFIT 6.1')
def ne_coarse(l):
    r = r_of(costhz, l/gd.UNIT_KM)
    i = int(np.searchsorted(EDGES, r))
    return matter.num_density_e_func(
        r, lambda _: RHO[i], electron_fraction=YE[i],
        ratio_number_neutrons_to_protons=
            (1.0 - YE[i])/YE[i],
        density_matter_is_in_g_per_cm3=True)
cross = earth.prem_layer_edges_along_chord(costhz)
keep = cross[np.isclose(
    r_of(costhz, cross)[:, None], EDGES).any(axis=1)]
print('keep', len(keep), keep)
E = 10.0*gd.UNIT_GEV
kw = dict(nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=1.0e-8, atol=1.0e-10)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    P = oscprob.osc_prob_matter_std_potential(
        3, ne_coarse, E, L, osc_params=osc, L0=0.0,
        t_breakpoints=keep*gd.UNIT_KM,
        density_is_of_number_of_electrons=True, **kw)
    print('warnings', [str(x.message)[:120] for x in w])
print('coarse P', repr(P))
YE2=[0.5]*4
def ne_half(l):
    r = r_of(costhz, l/gd.UNIT_KM); i = int(np.searchsorted(EDGES, r))
    return matter.num_density_e_func(r, lambda _: RHO[i], electron_fraction=0.5, ratio_number_neutrons_to_protons=1.0, density_matter_is_in_g_per_cm3=True)
warnings.simplefilter('ignore')
P5 = oscprob.osc_prob_matter_std_potential(3, ne_half, E, L, osc_params=osc, L0=0.0, t_breakpoints=keep*gd.UNIT_KM, density_is_of_number_of_electrons=True, **kw)
print('coarse Ye=0.5', repr(P5))
Pp = oscprob.osc_prob_3nu_earth(E, costhz=costhz, L=L, **kw)
print('PREM', repr(Pp))
Pp1 = oscprob.osc_prob_3nu_earth(E, costhz=costhz, L=L, nu_i=gd.NUMU, nu_f=gd.NUMU)
print('PREM default tol', repr(Pp1))
# mass-weighted means
from scipy.integrate import quad
f = earth.density_matter_func_prem
B = list(earth.PREM_BOUNDARIES)+[6371.0]
def mean(r0, r1):
    pts=[0.0]+B
    num=0; den=0
    for a,b in zip(pts[:-1],pts[1:]):
        lo=max(a,r0); hi=min(b,r1)
        if hi<=lo: continue
        num+=quad(lambda r: float(f(r))*r*r, lo, hi, epsabs=0, epsrel=1e-12)[0]; den+=(hi**3-lo**3)/3
    return num/den
print('means', [round(mean(a,b),3) for a,b in ((0,1221.5),(1221.5,3480),(3480,6346.6),(6346.6,6371))])
