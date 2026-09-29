import numpy as np
import magnus.globaldefs as gd
import magnus.oscprob as oscprob

# Standard oscillation parameters
osc = gd.load_nufit_params('NuFIT 6.1')
# Active-sterile mixing
s14 = s24 = np.sqrt(0.10)
s15 = s25 = np.sqrt(0.06)

def energies(lo, hi):
    return np.logspace(np.log10(lo), np.log10(hi),
                       140)*gd.UNIT_GEV

# Exponentially falling density: 3e3 g/cm^3 at the
# origin, scale height 10 km, 25 km baseline.
# Every mixing parameter left out defaults to zero.
common = dict(L=25.0*gd.UNIT_KM, L0=0.0,
              rho_central=3.e3,
              l_scale=10.0*gd.UNIT_KM,
              density_matter_is_in_g_per_cm3=True,
              nu_i=gd.NUE, nu_f=gd.NUE,
              rtol=1e-12, atol=1e-14,
              magnus_exp_order=8, strategy='hybrid')

P_2nu = oscprob.osc_prob_2nu_matter_exp_density(
    energies(0.0005, 0.05), **common,
    sth=osc['s12'], Dm2=osc['D21'])

P_3nu = oscprob.osc_prob_3nu_matter_exp_density(
    energies(0.002, 0.2), **common, **osc)

P_4nu = oscprob.osc_prob_4nu_matter_exp_density(
    energies(2.0, 20.0), **common, **osc,
    s14=s14, s24=s24, D41=1.0)

P_5nu = oscprob.osc_prob_5nu_matter_exp_density(
    energies(2.0, 20.0), **common, **osc,
    s14=s14, s15=s15, s24=s24, s25=s25,
    D41=1.0, D51=1.7)
