import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
osc = gd.load_nufit_params('NuFIT 6.1')
print({k: osc[k] for k in osc})
n, width = 24, 250.0*gd.UNIT_KM
edges = np.arange(n + 1)*width
def mk(r):
    def f(l):
        k = np.searchsorted(edges, l, side='right') - 1
        return r[np.clip(k, 0, n - 1)]
    return f
castle = np.where(np.arange(n) % 2 == 0, 2.0, 8.0)
uni = np.full(n,5.0)
def run(r, E, nb):
    return np.asarray(oscprob.osc_prob_matter_std_potential(3, mk(r), E*gd.UNIT_GEV, n*width, osc,
        t_breakpoints=edges[1:-1], nubar=nb, nu_i=gd.NUMU, nu_f=gd.NUE,
        density_matter_is_in_g_per_cm3=True, rtol=1e-8, atol=1e-10))
for nb, lo, hi in ((False, 0.450, 0.460),):
    E = np.linspace(lo, hi, 201)
    p = run(castle, E, nb); u = run(uni, E, nb)
    k = np.argmax(p)
    print('nubar', nb, 'castle peak %.5f at %.4f GeV; uniform there %.5f; ratio %.3f' % (p[k], E[k], u[k], p[k]/u[k]))
E = np.linspace(5.0, 6.0, 401)
for name, r in ():
    pass # name!='random' else np.linspace(5.3,6.3,401), False); k=np.argmax(p)
    pass # peak %.4f at %.3f'%(p[k], (E if name!='random' else np.linspace(5.3,6.3,401))[k]))
# matter oscillation lengths at 0.46 GeV
s12,s13,s23,dcp,D21,D31 = osc['s12'],osc['s13'],osc['s23'],osc['dCP'],osc['D21'],osc['D31']
c12,c13,c23=[np.sqrt(1-x**2) for x in (s12,s13,s23)]
e=np.exp(1j*dcp)
U = np.array([[1,0,0],[0,c23,s23],[0,-s23,c23]])@np.array([[c13,0,s13/e],[0,1,0],[-s13*e,0,c13]])@np.array([[c12,s12,0],[-s12,c12,0],[0,0,1]])
for nb in (False,True):
  for Eg in (0.46, 0.52):
    for rho in (2.0, 8.0):
        En = Eg*gd.UNIT_GEV
        Uu = U.conj() if nb else U
        Hv = Uu@np.diag([0,D21,D31])@Uu.conj().T/(2*En)
        ne = rho*gd.N_AV*0.5/gd.CONV_CM_TO_INV_EV**3
        V = np.sqrt(2)*gd.GF*ne*(-1 if nb else 1)
        w = np.linalg.eigvalsh(Hv+np.diag([V,0,0]))
        dl = [2*np.pi/abs(w[j]-w[i])/gd.UNIT_KM for i,j in ((0,1),(1,2),(0,2))]
        print('nubar',nb,'E',Eg,'rho',rho,'osc lengths km', ['%.0f'%x for x in dl])
