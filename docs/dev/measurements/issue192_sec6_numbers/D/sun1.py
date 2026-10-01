import magnus.oscprob as oscprob; import magnus.globaldefs as gd; osc = gd.load_nufit_params('NuFIT 6.1'); import time; t0=time.time()
import numpy as np
import magnus.adiabatic as adiabatic
import magnus.hamiltonians as ham
import magnus.matter as matter
import magnus.globaldefs as gd
import magnus.solarmodels as solarmodels

R = gd.SUN_RADIUS*gd.UNIT_KM
ne_sun = solarmodels.electron_density_profile(
    'B16-GS98')
b = 0.3*R                # impact parameter
half = np.sqrt(R**2 - b**2)
E = 100.0*gd.UNIT_GEV
hv = ham.hamiltonian_3nu_vacuum(E, **osc)

def ne(l):
    """Electron density on the chord."""
    r = np.sqrt((l - half)**2 + b**2)
    return ne_sun(r)

def H(l):
    """Hamiltonian at l along the chord."""
    h = np.array(hv, dtype=complex)
    h[0, 0] += matter.VCC_func(l, ne)
    return h

info = {}
adiabatic.find_nonadiabatic_windows(
    H, 0.0, 2*half, n_probe=20000,
    info=info)
print("gamma_max", repr(info["gamma_max"]), info.keys())
print({k:v for k,v in info.items() if k!='gamma_max' and not hasattr(v,'__len__')})
import magnus.oscprob as oscprob

# The chord above: b = 0.3 R_sun
P = oscprob.osc_prob_matter_std_potential(
    3, ne, E, 2*half, average=True,
    average_initial_state='decohered',
    osc_params=osc, L0=0.0, nu_i=gd.NUE,
    nu_f=gd.NUE,
    density_is_of_number_of_electrons=True)
print("P", repr(P), time.time()-t0)
kw = dict(
    osc_params=osc, L0=0.0, nu_i=gd.NUE,
    nu_f=gd.NUE, average=True,
    density_is_of_number_of_electrons=True)
f = oscprob.osc_prob_matter_std_potential
for start in ('flavor', 'decohered'):
    P = f(3, ne, E, 2*half, **kw,
          average_initial_state=start)
    print(start, P)  # 0.535944, 0.283618
