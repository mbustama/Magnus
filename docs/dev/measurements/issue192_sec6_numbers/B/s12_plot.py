import numpy as np, warnings, os
warnings.simplefilter('ignore')
import matplotlib; matplotlib.use('Agg')
import magnus.globaldefs as gd, magnus.hamiltonians as ham, magnus.matter as matter, magnus.oscprob as oscprob
osc = gd.load_nufit_params('NuFIT 6.1')
E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
H_vac = ham.hamiltonian_3nu_vacuum_energy_independent(s12=osc['s12'], s23=osc['s23'], s13=osc['s13'], dCP=osc['dCP'], D21=osc['D21'], D31=osc['D31'])
proj = matter.matter_potential_projector(3)
for which in ('wrapper','energy_baseline'):
    if which=='wrapper':
        P_scan = oscprob.osc_prob_3nu_matter_constant_density(Es, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc)
    else:
        P_scan = oscprob.osc_prob_energy_baseline(lambda e: H_vac/e + gd.VCC_EARTH_CRUST*proj, Es, L, H_func_is_function_only_of_energy=True)
    P_scan = np.asarray(P_scan); print(which, P_scan.shape)
    import magnus.plotting as plotting
    P_vac = oscprob.osc_prob_3nu_vacuum(Es, L)
    fig, ax = plotting.\
        plot_probability_vs_energy(
        Es/gd.UNIT_GEV,
        [dict(y=P_vac[:, gd.NUMU, gd.NUE],
              label='Vacuum'),
         dict(y=P_scan[:, gd.NUMU, gd.NUE],
              label='Matter')],
        nu_i=gd.NUMU, nu_f=gd.NUE)
    fig.set_size_inches(6, 3); ax.grid(True); fig.savefig('numu_to_nue_%s.pdf'%which)
    print('ok', which, os.path.getsize('numu_to_nue_%s.pdf'%which))
# oscillogram listing, reduced CZ
import matplotlib.pyplot as plt
import magnus.earth as earth
import magnus.plotting as plotting
chord = earth.distance_traveled_inside_earth
CZ = np.linspace(-1.0, -0.05, 170)[::42]
E_GEV = np.logspace(np.log10(2.0), np.log10(60.0), 200)*gd.UNIT_GEV
E_TEV = np.logspace(0.0, np.log10(30.0), 200)*gd.UNIT_TEV
kw = dict(nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=1.0e-8, atol=1.0e-10)
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
s4 = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, D41=1.0)
s5 = dict(**s4, s15=np.sqrt(0.06), s25=np.sqrt(0.06), s35=0.0, D51=1.7)
def oscillogram(wrapper, E, **params):
    return np.array([wrapper(E, costhz=c, L=chord(c)*gd.UNIT_KM, **params, **kw) for c in CZ]).T
PANELS = [(oscprob.osc_prob_3nu_earth, E_GEV, {}), (oscprob.osc_prob_3nu_earth_nsi, E_GEV, eps), (oscprob.osc_prob_4nu_earth, E_TEV, s4), (oscprob.osc_prob_5nu_earth, E_TEV, s5)]
P = [oscillogram(w, E, **p) for w, E, p in PANELS]
for (_, E, _), grid in zip(PANELS, P):
    fig, ax = plotting.plot_oscillogram(CZ, np.log10(E/gd.UNIT_GEV), grid, nu_i=gd.NUMU, nu_f=gd.NUMU)
for (_, E, _), grid in zip(PANELS, P):
    fig, ax = plt.subplots()
    im = ax.pcolormesh(CZ, np.log10(E/gd.UNIT_GEV), grid, shading='gouraud')
    fig.colorbar(im, ax=ax)
print('oscillogram listing ok (CZ subset of', len(CZ), ')', [g.shape for g in P])
