import sys, time, numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.solarmodels as solarmodels, magnus.avgprob as avgprob, magnus.hamiltonians as ham, magnus.matter as matter
osc = gd.load_nufit_params('NuFIT 6.1')
R = gd.SUN_RADIUS*gd.UNIT_KM
ne_sun = solarmodels.electron_density_profile('B16-GS98')
def chord(br):
    b=br*R; half=np.sqrt(R**2-b**2)
    return (lambda l: ne_sun(np.sqrt((np.asarray(l,dtype=float)-half)**2+b**2))), half
def P(Egev, br):
    ne, half = chord(br)
    return float(oscprob.osc_prob_matter_std_potential(3, ne, Egev*gd.UNIT_GEV, 2*half, average=True,
        average_initial_state='decohered', osc_params=osc, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True))
def Pvary(Egev, br):
    ne, half = chord(br); E=Egev*gd.UNIT_GEV
    hv = np.array(ham.hamiltonian_3nu_vacuum(E, **osc), dtype=complex)
    def H(l):
        h=hv.copy(); h[0,0]+=matter.VCC_func(l, ne); return h
    Pm, rep = avgprob.averaged_probabilities_adiabatic(H, 0.0, 2*half)
    return float(Pm[0,0])
if __name__=='__main__':
    t=time.time(); print(P(300.0, 0.81), time.time()-t); t=time.time(); print(Pvary(300.0,0.81), time.time()-t)
