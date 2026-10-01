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

# An empty crust needs no profile, only a number
P0 = oscprob.osc_prob_matter_std_potential(
    3, NE_CRUST, E, L0*gd.UNIT_KM, **kw)

# Water, iron, a mineral deposit, a zone of faults:
# (rho, Y_e, w) for each
CASES = [(1.0, 0.555, 250.0), (5.0, 0.5, 250.0),
         (10.0, 0.5, 100.0), (25.0, 0.5, 50.0)]
dP = []
for rho, ye, w in CASES:
    profile, walls = cavity(rho, ye, w)
    # The two walls are density jumps: declare them
    P = oscprob.osc_prob_matter_std_potential(
        3, profile, E, L0*gd.UNIT_KM,
        t_breakpoints=walls, **kw)
    dP.append(P - P0)

Em = E/gd.UNIT_MEV
i = P0.argmax()
print('peak E %.4f MeV  P %.6f' % (Em[i], P0[i]))
for (rho, ye, w), d in zip(CASES, dP):
    s = np.sign(d)
    idx = np.where(np.diff(s) != 0)[0]
    # zero crossing nearest 49 MeV
    xs = [Em[k] - d[k]*(Em[k+1]-Em[k])/(d[k+1]-d[k]) for k in idx]
    xs = np.array(xs)
    j = np.argmin(abs(xs-49.2))
    print(rho, ye, w, 'crossing near peak %.3f' % xs[j], 'maxshift %+.3f' % d[np.abs(d).argmax()])
print('alpha', np.degrees(np.arcsin(125/750)))
