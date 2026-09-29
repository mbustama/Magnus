"""Step 4a: every hybrid-path result at atol + rtol >= 1e-6, saved for a bitwise comparison
between the HEAD export and the fix.  Run with PYTHONPATH pointing at either src."""
import sys, warnings, numpy as np
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.oscprob as op, magnus.matter as matter, magnus.hamiltonians as hams
out = {}
osc = gd.load_nufit_params('NuFIT 6.1')
s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
def E(lo, hi, n=140): return np.logspace(np.log10(lo), np.log10(hi), n)*gd.UNIT_GEV
base = dict(L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.e3, l_scale=10.0*gd.UNIT_KM,
            density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE)
for tol in (1e-3, 1e-6):
    kw = dict(base, rtol=tol, atol=tol, strategy='hybrid')
    out['l1_2nu_%g' % tol] = op.osc_prob_2nu_matter_exp_density(E(0.0005, 0.05), **kw, sth=osc['s12'], Dm2=osc['D21'])
    out['l1_3nu_%g' % tol] = op.osc_prob_3nu_matter_exp_density(E(0.002, 0.2), **kw, **osc)
    out['l1_4nu_%g' % tol] = op.osc_prob_4nu_matter_exp_density(E(2.0, 20.0), **kw, **osc, s14=s14, s24=s24, D41=1.0)
    out['l1_5nu_%g' % tol] = op.osc_prob_5nu_matter_exp_density(E(2.0, 20.0), **kw, **osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7)
    for L_km, h_km, rho0 in ((250, 100, 3e3), (2500, 1000, 3e3)):
        k2 = dict(L=L_km*gd.UNIT_KM, L0=0.0, rho_central=rho0, l_scale=h_km*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, rtol=tol, atol=tol, strategy='hybrid')
        out['exp%d_3nu_%g' % (L_km, tol)] = op.osc_prob_3nu_matter_exp_density(E(0.002, 0.2, 7), **k2, **osc)
        out['exp%d_nsi_%g' % (L_km, tol)] = op.osc_prob_3nu_matter_nsi_exp_density(E(0.02, 0.2, 5), **k2, **osc, eps_ee=0.15, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
        out['exp%d_liv_%g' % (L_km, tol)] = op.osc_prob_3nu_matter_liv_exp_density(E(0.02, 0.2, 5), **k2, **osc, sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0, b1=gd.B1, b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1)
    R = gd.SUN_RADIUS*gd.UNIT_KM
    out['sun_2nu_%g' % tol] = op.osc_prob_2nu_sun(np.array([5e6, 10e6]), R, 0.0, sth=osc['s12'], Dm2=osc['D21'], rtol=tol, atol=tol, strategy='hybrid')
    out['sun_3nu_b16_%g' % tol] = op.osc_prob_3nu_sun(np.array([5e6, 15e6]), R, 0.0, **osc, density_profile='B16-GS98', rtol=tol, atol=tol, strategy='hybrid')
    NE0, LS = gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN
    def ne_multi(l):
        x = np.asarray(l, dtype=float); o = NE0*np.exp(-x/LS)*(1.0 + 0.9*np.sin(2.0*np.pi*6.0*x/LS)); return o[()] if o.ndim == 0 else o
    out['multi_%g' % tol] = op.osc_prob_matter_std_potential(2, ne_multi, np.array([10e6, 50e6]), LS, {'sth': 0.55, 'Dm2': 7.5e-5}, L0=0.0, density_is_of_number_of_electrons=True, rtol=tol, atol=tol, strategy='hybrid')
    hvac = hams.hamiltonian_2nu_vacuum_energy_independent(np.sqrt(0.308), 7.5e-5); e00 = np.diag([1.0, 0.0])
    def H(energy, l, VCC): return (1.0/energy)*hvac + np.asarray(VCC)[..., None, None]*e00
    out['generic_sun_%g' % tol] = op.osc_prob_sun(H, 10e6, 0.3*R, rtol=tol, atol=tol, strategy='hybrid', validate_input=False)
out['average_1e-3'] = op.osc_prob_3nu_sun(np.array([5e6, 10e6]), gd.SUN_RADIUS*gd.UNIT_KM, 0.0, **osc, average=True)
np.savez(sys.argv[1], **{k: np.asarray(v) for k, v in out.items()})
print(op.__file__, len(out), 'results')
