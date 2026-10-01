import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.oscprob as op
from scipy.stats import unitary_group
E=1.0*gd.UNIT_GEV; L=1000*gd.UNIT_KM
V0=gd.VCC_EARTH_CRUST
for n in (5,6,8):
    Ur=unitary_group.rvs(n, random_state=1)
    m2=np.r_[0, 7.5e-5, 2.5e-3, np.linspace(0.5,2,n-3)]
    Hv=Ur@np.diag(m2)@Ur.conj().T/(2*E)
    P0=np.zeros(n); P0[0]=1
    def H(l, Hv=Hv, n=n):
        h=Hv.copy(); v=V0*np.exp(-l/L)
        h[0,0]+=v
        for k in range(3,n): h[k,k]+=v/2
        return h
    P=np.asarray(op.osc_prob(H,0.0,L))
    print(n, np.abs(P.sum(0)-1).max(), np.abs(P.sum(1)-1).max())
