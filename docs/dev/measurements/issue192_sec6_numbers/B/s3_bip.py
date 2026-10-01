import numpy as np, warnings, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
import magnus.plotting as plotting
E, L, rho = 2.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM, 3.0
kw = dict(nu_i=gd.NUMU, nu_f=gd.NUE, density_matter_is_in_g_per_cm3=True)
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
std = oscprob.osc_prob_3nu_matter_constant_density
nsi = oscprob.osc_prob_3nu_matter_nsi_constant_density
NPT = 181
def locus(wrapper, **couplings):
    out = []
    for dcp in np.linspace(-np.pi, np.pi, NPT):
        P = [wrapper(E, L, rho, nubar=nubar, dCP=dcp, **couplings, **kw) for nubar in (False, True)]
        out.append(P)
    return np.array(out)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    P_std = locus(std); P_nsi = locus(nsi, **eps)
    print('warnings', len(w), set(str(x.category.__name__) for x in w))
fig, ax = plotting.plot_biprobability([P_std[:, 0], P_nsi[:, 0]], [P_std[:, 1], P_nsi[:, 1]], labels=['Standard', 'NSI'])
fig, ax = plt.subplots()
for P, label in ((P_std, 'Standard'), (P_nsi, 'NSI')):
    ax.plot(P[:, 0], P[:, 1], label=label)
ax.legend()
print('listing ran; shapes', P_std.shape)
# crossings with fine grid
NPT = 4001
d = np.linspace(-np.pi, np.pi, NPT)
A = locus(std); B = locus(nsi, **eps)
def seg_int(p1,p2,q1,q2):
    r = p2-p1; s = q2-q1; den = r[0]*s[1]-r[1]*s[0]
    if den == 0: return None
    t = ((q1-p1)[0]*s[1]-(q1-p1)[1]*s[0])/den; u = ((q1-p1)[0]*r[1]-(q1-p1)[1]*r[0])/den
    if 0<=t<=1 and 0<=u<=1: return t,u
    return None
# coarse filtering by bounding boxes
res=[]
for i in range(NPT-1):
    p1,p2=A[i],A[i+1]
    lo=np.minimum(p1,p2); hi=np.maximum(p1,p2)
    bl=np.minimum(B[:-1],B[1:]); bh=np.maximum(B[:-1],B[1:])
    cand=np.where((bl[:,0]<=hi[0])&(bh[:,0]>=lo[0])&(bl[:,1]<=hi[1])&(bh[:,1]>=lo[1]))[0]
    for j in cand:
        r=seg_int(p1,p2,B[j],B[j+1])
        if r:
            t,u=r; x=p1+t*(p2-p1)
            res.append((x[0],x[1],(d[i]+t*(d[1]-d[0]))/np.pi,(d[j]+u*(d[1]-d[0]))/np.pi))
for r in res: print('cross P=%.4f Pbar=%.4f dCP_std=%.4f pi dCP_nsi=%.4f pi'%r)
print('bestfit dCP/pi', gd.load_nufit_params('NuFIT 6.1')['dCP']/np.pi)
