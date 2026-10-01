import numpy as np, warnings, time
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
E = 5.0*gd.UNIT_TEV
D41 = np.logspace(np.log10(0.1), np.log10(3.0), 80)
S14 = np.logspace(np.log10(0.02), 0.0, 80)
s24 = np.sqrt(0.10)
def kwf(costhz, nubar):
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
    return dict(costhz=costhz, L=L, nu_i=gd.NUMU, nu_f=gd.NUMU, nubar=nubar, rtol=1.0e-8, atol=1.0e-10)
def point(cz, nb, i, j, std):
    return oscprob.osc_prob_4nu_earth(E, s14=S14[j], s24=s24, s34=0.0, D41=D41[i], **kwf(cz,nb)) - std
res = {}
for cz in (-1.0, -0.5):
    for nb in (False, True):
        kw = kwf(cz, nb); std = oscprob.osc_prob_3nu_earth(E, **kw)
        I = range(0,80,2); J = range(0,80,2)
        g = np.array([[point(cz,nb,i,j,std) for j in J] for i in I])
        # refine around extreme
        a = np.unravel_index(np.argmax(np.abs(g)), g.shape); i0, j0 = 2*a[0], 2*a[1]
        best = (abs(g[a]), g[a], i0, j0)
        for i in range(max(0,i0-3), min(80,i0+4)):
            for j in range(max(0,j0-3), min(80,j0+4)):
                v = point(cz,nb,i,j,std)
                if abs(v) > best[0]: best = (abs(v), v, i, j)
        print('cz', cz, 'nubar', nb, 'coarse min %.4f max %.4f'%(g.min(), g.max()), 'refined extreme %.4f at D41=%.3f sin2=%.4f'%(best[1], D41[best[2]], S14[best[3]]**2), flush=True)
        np.save('ster_%s_%s.npy'%(cz,nb), g)
