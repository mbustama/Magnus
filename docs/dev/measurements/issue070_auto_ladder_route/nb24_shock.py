import warnings, numpy as np
import magnus.globaldefs as gd, magnus.oscprob as oscprob
OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
MEAN_NUCLEON_D = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
R0_D, R1_D, W_D, RF_D = 1.0e4, 8.0e4, 1.0e-3, 3.0e4
KM_D = gd.CONV_KM_TO_INV_EV
def ne_shock_d(l):
    w_km = W_D*(R1_D - R0_D)
    r = np.asarray(l, dtype=float)/KM_D
    u = np.clip((RF_D + 0.5*w_km - r)/w_km, 0.0, 1.0)
    out = (1.0e14*r**(-2.4)*(1.0 + (u*u*(3.0 - 2.0*u))*9.0)*gd.UNIT_G_PER_CM3/MEAN_NUCLEON_D*0.5)
    return out[()] if np.ndim(out) == 0 else out
for rtol in np.logspace(-2, -10, 9):
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p = float(np.asarray(oscprob.osc_prob_matter_std_potential(3, ne_shock_d, 15.0*gd.UNIT_MEV, R1_D*KM_D, OSC, L0=R0_D*KM_D,
            density_is_of_number_of_electrons=True, nu_i=gd.NUE, nu_f=gd.NUE, rtol=rtol, atol=rtol*1.0e-2, strategy_info=info)))
    tr = [t for t in info.get('trace', []) if t.get('reason') == 'auto prefers the ladder']
    print('%9.0e %.6f %-10s %s %s' % (rtol, p, info.get('engine'), (info.get('declined') or [['', '--']])[0][1], ('phase=%.0f' % tr[0]['estimated_phase']) if tr else ''))
