import numpy as np
import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.oscprob as oscprob

OSC = gd.load_nufit_params('NuFIT 6.1')
CM = 1.0e-5*gd.UNIT_KM   # One cm in eV^-1
RSTAR = 3.0e12           # Star's radius, cm
R0 = 6.3e10              # Jet head, cm

def rho(l):
    """Density, g/cm3, at l in eV^-1."""
    return 3.3e-6*(RSTAR/(l/CM) - 1.0)**3
E = 1.0*gd.UNIT_TEV
out = oscprob.osc_prob_matter_std_potential(
    3, rho, E, RSTAR*CM, OSC, L0=R0*CM,
    electron_fraction=1.0,
    rtol=1.0e-6, atol=1.0e-6,
    density_matter_is_in_g_per_cm3=True,
    return_evolution_operator=True)
P, U = out
# Vacuum mass states; the content of each in
# the state that leaves the star, by flavor
hv = np.array(ham.hamiltonian_3nu_vacuum(
    E, **OSC), dtype=complex)
w, R = np.linalg.eigh(hv)
content = abs(R.conj().T @ U)**2

# Phases average on the way to Earth
P_earth = abs(R)**2 @ content
print("P_earth", P_earth[gd.NUE, gd.NUE])
import magnus.matter as matter
print('rho head', rho(R0*CM))
V1 = matter.VCC_func(0.0, lambda l: matter.num_density_e_func(0.0, lambda _: 1.0, electron_fraction=1.0, density_matter_is_in_g_per_cm3=True))
print('Losc head cm %.4g' % (2*np.pi/(V1*rho(R0*CM))/CM))
# without matter: vacuum averaged
print('no matter', (abs(R)**2 @ abs(R.conj().T)**2)[0,0])
# resonance density for 1-3 at E: V = D31 cos2th13/(2E)
for Etev in (0.1, 100.0, 1e3):
    Ee=Etev*gd.UNIT_TEV; Vres = OSC['D31']*(1-2*OSC['s13']**2)/(2*Ee); print('E %g TeV resonant rho %.3g' % (Etev, Vres/V1))
# jet density at 1e13-1e15 cm: profile is for r<RSTAR; (paper says ~1e-11) ; phase
