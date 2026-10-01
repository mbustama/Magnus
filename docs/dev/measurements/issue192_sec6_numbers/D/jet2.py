import time; t0=time.time()
import numpy as np
import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.oscprob as oscprob

OSC = gd.load_nufit_params('NuFIT 6.1')
CM = 1.0e-5*gd.UNIT_KM               # One cm in eV^-1
RSTAR, R0, RHE = 3.0e12, 6.3e10, 1.0e11      # cm
E = np.geomspace(0.1, 1.0e4, 161)*gd.UNIT_TEV

# The three envelopes, g/cm3, r in cm
def smooth(r):
    return 3.3e-6*(RSTAR/r - 1.0)**3

rng = np.random.default_rng(7)        
lam = np.geomspace(1.0e12, 1.0e10, 40)  # Modes, cm
k = 2*np.pi/lam
amp = k**(-1/3)                         # Kolmogorov
amp *= 0.20/np.sqrt(0.5*np.sum(amp**2)) # 20% rms
phi = rng.uniform(0.0, 2*np.pi, 40)
def turbulent(r):
    d = np.sum(amp*np.cos(k*r + phi))
    return smooth(r)*(1.0 + d)

def stepped(r):            # Mena et al., model C
    x = RSTAR/r - 1.0
    return 6.3e-6*(20.0*x**2.1 if r < RHE else x**2.5)

# Production points spread over one oscillation
# length at the jet head
r0s = R0 + 4.9e9*np.arange(8)/8.0

def p_earth(rho, **extra):
    """P at Earth, Eq. (26), averaged over the
    production points."""
    P_out = np.zeros((len(E), 3, 3))
    for r0 in r0s:
        P, U = oscprob.osc_prob_matter_std_potential(
            3, lambda l: rho(l/CM), E, RSTAR*CM, OSC,
            L0=r0*CM, electron_fraction=1.0,
            density_matter_is_in_g_per_cm3=True,
            rtol=1.0e-6, atol=1.0e-6,
            return_evolution_operator=True, **extra)
        for i, e in enumerate(E):
            hv = np.array(ham.hamiltonian_3nu_vacuum(
                e, **OSC), dtype=complex)
            w, R = np.linalg.eigh(hv)
            content = abs(R.conj().T @ U[i])**2
            P_out[i] += abs(R)**2 @ content/len(r0s)
    return P_out

P_smooth = p_earth(smooth)
P_turb = p_earth(turbulent)
P_step = p_earth(stepped, t_breakpoints=[RHE*CM])
# Drawn: P_smooth[:, gd.NUE, gd.NUE], and the others
print('t', time.time()-t0)
np.save('jet.npy', np.array([E/gd.UNIT_TEV, P_smooth[:,0,0], P_turb[:,0,0], P_step[:,0,0]]))
Et = E/gd.UNIT_TEV; ps, pt, pst = P_smooth[:,0,0], P_turb[:,0,0], P_step[:,0,0]
vac = 0.5481426092150248
print('smooth at 0.1 TeV %.4f, at 100 TeV %.4f, 1e4 TeV %.4f' % (ps[0], ps[np.argmin(abs(Et-100))], ps[-1]))
for x in (10, 30, 50, 100, 300, 1000):
    print('  smooth at %g TeV: %.4f (vac-P %.4f)' % (x, ps[np.argmin(abs(Et-x))], vac-ps[np.argmin(abs(Et-x))]))
d = pt-ps; i=abs(d).argmax(); print('turb max |dP| %.4f at %.3g TeV; max above 5 TeV %.4f' % (abs(d).max(), Et[i], abs(d[Et>5]).max()))
d = pst-ps; big = abs(d) > 0.02; print('stepped-smooth > 0.02 at E range', Et[big].min() if big.any() else None, Et[big].max() if big.any() else None, 'max', abs(d).max())
for x in Et[(Et>0.09)&(Et<1.5)][::3]: print('   %.3f  smooth %.3f step %.3f turb %.3f' % (x, ps[Et==x][0], pst[Et==x][0], pt[Et==x][0]))
