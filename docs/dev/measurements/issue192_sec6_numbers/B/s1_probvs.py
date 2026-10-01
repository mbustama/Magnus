import numpy as np
import magnus.oscprob as oscprob
import magnus.globaldefs as gd

osc = gd.load_nufit_params('NuFIT 6.1')
KM, MEV = gd.UNIT_KM, gd.UNIT_MEV
NA_CM3 = gd.N_AV/gd.CONV_CM_TO_INV_EV**3
ne_constant = 10.0*NA_CM3
def ne_exponential(l):
    return 10.0*NA_CM3*np.exp(-l/(100.0*KM))
def ne_gaussian(l):
    return 8.0*NA_CM3*np.exp(
     -(l - 300.0*KM)**2/(2*(100.0*KM)**2))
KW = dict(density_is_of_number_of_electrons=True,
          rtol=1e-6, atol=1e-6)
PROFILES = (('constant', ne_constant),
            ('exponential', ne_exponential),
            ('gaussian', ne_gaussian))
L = np.linspace(20.0, 500.0, 5000)*KM
P_vs_L = {name: oscprob.osc_prob_matter_std_potential(
 3, ne, 10.0*MEV, L, osc, L0=0.0, **KW)
 for name, ne in PROFILES}
E = np.logspace(np.log10(3.0), 2.0, 3000)*MEV
P_vs_E = {name: oscprob.osc_prob_matter_std_potential(
 3, ne, E, 200.0*KM, osc, L0=0.0, **KW)
 for name, ne in PROFILES}
P_vac_L = oscprob.osc_prob_3nu_vacuum(10.0*MEV, L, **osc)
P_vac_E = oscprob.osc_prob_3nu_vacuum(E, 200.0*KM, **osc)
for k,v in P_vs_L.items(): print('L', k, np.shape(v), repr(np.asarray(v)[-1][gd.NUE][gd.NUE]))
for k,v in P_vs_E.items(): print('E', k, np.shape(v))

P = oscprob.osc_prob_3nu_matter_exp_density(
 10.0*MEV, 500.0*KM, L0=0.0,
 rho_central=10.0*NA_CM3,
 l_scale=100.0*KM, **osc,
 nu_i=gd.NUE, nu_f=gd.NUE, **KW)
print('exp wrapper', repr(P))
P = oscprob.\
 osc_prob_3nu_matter_constant_density(
 10.0*MEV, 500.0*KM, 10.0*NA_CM3, **osc,
 nu_i=gd.NUE, nu_f=gd.NUE, **KW)
print('const wrapper', repr(P))

info = {}
P = oscprob.osc_prob_matter_std_potential(
 3, ne_exponential, E, 200.0*KM, osc,
 L0=0.0, strategy_info=info, **KW)
print('engine', info['engine'])

# tight tolerances
KW2 = dict(density_is_of_number_of_electrons=True, rtol=1e-8, atol=1e-10)
mx = 0
for name, ne in PROFILES:
    a = oscprob.osc_prob_matter_std_potential(3, ne, 10.0*MEV, L, osc, L0=0.0, **KW2)
    b = oscprob.osc_prob_matter_std_potential(3, ne, E, 200.0*KM, osc, L0=0.0, **KW2)
    dL = np.max(np.abs(np.asarray(a)-np.asarray(P_vs_L[name]))); dE = np.max(np.abs(np.asarray(b)-np.asarray(P_vs_E[name])))
    print('tight diff', name, dL, dE)
