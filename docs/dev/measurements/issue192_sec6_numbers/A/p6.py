import numpy as np, warnings, time, sys
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.solarmodels as solarmodels
osc = gd.load_nufit_params('NuFIT 6.1')
R = gd.SUN_RADIUS*gd.UNIT_KM
ne_sun = solarmodels.electron_density_profile('B16-GS98')
def chord(b):
    half = np.sqrt(R**2 - b**2)
    return (lambda l: ne_sun(np.sqrt((l - half)**2 + b**2))), half
ne, half = chord(0.3*R)
E = 100.0*gd.UNIT_GEV
kw = dict(osc_params=osc, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, average=True, density_is_of_number_of_electrons=True)
f = oscprob.osc_prob_matter_std_potential
r={}
for start in ('flavor','decohered'):
    t=time.time(); r[start]=f(3, ne, E, 2*half, **kw, average_initial_state=start); print(start, r[start], time.time()-t, flush=True)
print('2021 diff', r['flavor']-r['decohered'], flush=True)
if len(sys.argv)>1:
    rng = np.random.default_rng(0)
    bs = np.linspace(0.0, 0.2, 5)*R
    Es = np.logspace(np.log10(30), np.log10(3000), 4)*gd.UNIT_GEV
    mx=0
    for b in bs:
        ne, half = chord(b)
        for e in Es:
            for start in ('decohered','flavor'):
                p0 = f(3, ne, e, 2*half, **kw, average_initial_state=start)
                p1 = f(3, ne, e, 2*half, **kw, average_initial_state=start, rtol=1e-5, atol=1e-5)
                d=abs(p0-p1); mx=max(mx,d)
                print('b=%.2f E=%.0f GeV %s %.8f %.8f %.2e'%(b/R, e/gd.UNIT_GEV, start, p0, p1, d), flush=True)
    print('max', mx)
