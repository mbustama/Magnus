import numpy as np
import magnus.oscprob as oscprob
import magnus.matter as matter
import magnus.globaldefs as gd

L0 = 1500.0                     # km, source to detector
RHO_CRUST, YE_CRUST = 3.3, 0.5  # Near PREM's crust
E = np.linspace(25.0, 150.0, 2500)*gd.UNIT_MEV
osc = gd.load_nufit_params('NuFIT 6.1')
kw = dict(osc_params=osc, L0=0.0, nu_i=gd.NUE,
          nu_f=gd.NUE, nubar=True, rtol=1.0e-8,
          atol=1.0e-10,
          density_is_of_number_of_electrons=True)

def n_e(rho, ye):
    """Electron density [eV^3] of uniform matter."""
    return matter.num_density_e_func(
        0.0, lambda _: rho, electron_fraction=ye,
        ratio_number_neutrons_to_protons=(1.0-ye)/ye,
        density_matter_is_in_g_per_cm3=True)

NE_CRUST = n_e(RHO_CRUST, YE_CRUST)

def cavity(rho, ye, w):
    """Crust holding a cavity of width w, centered."""
    d = (L0 - w)/2.0
    ne_in = n_e(rho, ye)

    def profile(l):
        x = np.asarray(l, dtype=float)/gd.UNIT_KM
        return NE_CRUST + (ne_in - NE_CRUST) \
            *((x >= d) & (x <= d + w))

    return profile, np.array([d, d + w])*gd.UNIT_KM
import numpy as np

BODY_R, BODY_D0 = 125.0, 750.0   # km: size, position
NE_BODY = n_e(10.0, 0.5)         # a mineral deposit
E = np.linspace(25.0, 150.0, 400)*gd.UNIT_MEV
ALPHA = np.linspace(-15.0, 15.0, 220)

def crossing(alpha_deg):
    """Entry and exit points in the body [km]."""
    a = np.radians(alpha_deg)
    miss = abs(BODY_D0*np.sin(a))
    if miss >= BODY_R:
        return None              # the beam misses it
    half = np.sqrt(BODY_R**2 - miss**2)
    mid = BODY_D0*np.cos(a)
    return mid - half, mid + half

P0 = oscprob.osc_prob_matter_std_potential(
    3, NE_CRUST, E, L0*gd.UNIT_KM, **kw)

dP = np.zeros((len(E), len(ALPHA)))
for j, alpha in enumerate(ALPHA):
    seg = crossing(alpha)
    if seg is None:
        continue                 # no body, no change
    lo, hi = seg

    def body(l, lo=lo, hi=hi):
        x = np.asarray(l, dtype=float)/gd.UNIT_KM
        return NE_CRUST + (NE_BODY - NE_CRUST) \
            *((x >= lo) & (x <= hi))

    # The walls move with the angle, so every call
    # declares its own pair
    dP[:, j] = oscprob.osc_prob_matter_std_potential(
        3, body, E, L0*gd.UNIT_KM,
        t_breakpoints=np.array(seg)*gd.UNIT_KM,
        **kw) - P0
print('sweep ok', dP.min(), dP.max(), E[P0.argmax()]/gd.UNIT_MEV)
