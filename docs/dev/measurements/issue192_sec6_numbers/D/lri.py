import numpy as np
import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.matter as matter
import magnus.oscprob as oscprob
import magnus.solarmodels as solarmodels

osc = gd.load_nufit_params('NuFIT 6.1')

# The BS2005-AGS,OP model that ships with Magnus:
# the electron density at its tabulated radii
model = 'BS05-AGS-OP'
tab = solarmodels.load_solar_model(model)
r = tab['r_over_r_sun']*gd.SUN_RADIUS*gd.UNIT_KM
R_SUN = r[-1]                       # Table's last row
n_e = solarmodels.electron_density_profile(model)(r)
vcc = matter.VCC_func(0.0, lambda l: 1.0)*n_e

def vcc_sun(l):
    """V_CC at l, log-linear between table rows."""
    return np.exp(np.interp(l, r, np.log(vcc)))

def cumulative(y, x):
    """Running trapezoidal integral, zero at x[0]."""
    steps = 0.5*(y[1:] + y[:-1])*np.diff(x)
    return np.concatenate([[0.0], np.cumsum(steps)])

def shc(x):
    """sinh(x)/x, equal to 1 at the origin."""
    x = np.asarray(x, dtype=float)
    return np.divide(np.sinh(x), x, out=np.ones_like(x),
                     where=(x != 0.0))

def long_range_potential(r, n_e, m):
    """Eq. (35): V_emu(r)/g'^2, spherical n_e."""
    I_in = cumulative(r**2*n_e*shc(m*r), r)
    I_out = cumulative(r*n_e*np.exp(-m*r), r)
    I_out = I_out[-1] - I_out      # From r outward
    r_safe = np.where(r > 0.0, r, 1.0e-30)
    return np.exp(-m*r)*I_in/r_safe + shc(m*r)*I_out

# Two mediator ranges.  g'^2 is fixed so that V_emu is
# a tenth of V_CC at the center in both cases.
v_emu = {}
for frac in (1.0, 0.1):            # 1/m in R_sun
    v = long_range_potential(r, n_e, 1.0/(frac*R_SUN))
    v_emu[frac] = 0.1*vcc[0]/v[0]*v

# The standard part from the shipped builders; the
# new term is the only line written here
h_vac = ham.\
    hamiltonian_3nu_vacuum_energy_independent(**osc)
q = np.diag([1.0, -1.0, 0.0])     # L_e-L_mu charges

def H_lri(v):
    """Eq. (33) with V_emu = v on r, as H(E, l)."""
    def H(E, l):
        h_matt = ham.hamiltonian_3nu_matter_td(
            l, vcc_sun)
        return (h_vac/E + h_matt
                + np.interp(l, r, v)[..., None, None]*q)
    return H
Es = np.logspace(-1.0, np.log10(20.0), 70)*gd.UNIT_MEV

def averaged(v):
    """<P_ee> from the center to the surface, per E."""
    return oscprob.osc_prob_energy_baseline(
        H_lri(v), Es, R_SUN, 0.0, nu_i=gd.NUE,
        nu_f=gd.NUE, average=True)

P_std = averaged(0.0*vcc)
P_lri = {frac: averaged(v_emu[frac]) for frac in v_emu}
half = np.searchsorted(r, 0.5*R_SUN)
print("2.49 ->", v_emu[1.0][half]/vcc[half]    )
print("0.38 ->", v_emu[0.1][half]/vcc[half]    )
i = np.argmin(abs(Es - 10.0*gd.UNIT_MEV))
print("0.318 ->", P_std[i])
print("0.322, 0.311 ->", P_lri[1.0][i], P_lri[0.1][i])
print('m', [1.0/(f*R_SUN) for f in (1.0, 0.1)], 'R_SUN km', R_SUN/gd.UNIT_KM)
print('surface ratios', v_emu[1.0][-1]/vcc[-1], v_emu[0.1][-1]/vcc[-1])
Em = Es/gd.UNIT_MEV
for f in (1.0,0.1):
    d = P_lri[f]-P_std; print(f, 'max |dP| %.4f' % abs(d).max(), 'min %.4f max %.4f' % (d.min(), d.max()))
    s=np.sign(d); idx=np.where(np.diff(s)!=0)[0]; print('  sign changes at', Em[idx], Em[idx+1])
for f in (1.0,):
    d = P_lri[f]-P_std; k=np.where(np.diff(np.sign(d))!=0)[0][0]
    print('interp crossing MeV', Em[k] - d[k]*(Em[k+1]-Em[k])/(d[k+1]-d[k]), 'log-interp', np.exp(np.log(Em[k]) - d[k]*(np.log(Em[k+1])-np.log(Em[k]))/(d[k+1]-d[k])))
