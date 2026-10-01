import numpy as np, warnings, inspect
import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.matter as matter
import magnus.oscprob as oscprob
print('nufit-named', len([k for k in gd.OSC_PARAMS_PREDEFINED if k.startswith('OSC_PARAMS_NU_FIT')]))
for f in ('osc_prob_4nu_vacuum','osc_prob_5nu_vacuum'):
    print(f, [p for p in inspect.signature(getattr(oscprob,f)).parameters])
osc = gd.load_nufit_params('NuFIT 6.1')
E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
for rho in (0.0, 3.0, 8.0, 13.0):
    P = oscprob.osc_prob_3nu_matter_constant_density(E, L, rho=rho, density_matter_is_in_g_per_cm3=True, **osc)
    print('%5.1f g/cm3  Pme = %.6f' % (rho, P[gd.NUMU][gd.NUE]), P[1,0])
P = oscprob.osc_prob_2nu_matter_constant_density(E, L, rho=3.0, sth=osc['s13'], Dm2=osc['D31'], density_matter_is_in_g_per_cm3=True)
print('1706', P[0][1])
Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
P = oscprob.osc_prob_3nu_matter_constant_density(Es, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc); print('1713', P.shape)
P = oscprob.osc_prob_3nu_matter_constant_density(E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, electron_fraction=0.2, **osc)
print('1724 Ye=0.2', P[1,0])
print('iron', 26/55.845)
st=dict(s14=np.sqrt(0.1), s24=np.sqrt(0.1), D41=1.0)
P0 = oscprob.osc_prob_4nu_matter_constant_density(E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc, **st)
P1 = oscprob.osc_prob_4nu_matter_constant_density(E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, ratio_number_neutrons_to_protons=1.5, **osc, **st)
print('1728', P0[1,1], P1[1,1], (P0[1,1]-P1[1,1])/P0[1,1], (P0[1,1]-P1[1,1])/P1[1,1])
# Listing constant
P_vac = oscprob.osc_prob_3nu_vacuum(E, L, **osc)
P_mat = oscprob.osc_prob_3nu_matter_constant_density(E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc)
print('1750 %.8f %.8f'%(P_vac[1,0], P_mat[1,0]))
P_scanw = oscprob.osc_prob_3nu_matter_constant_density(Es, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc)
print('1755', P_scanw.shape)
rho_func = lambda l: 3.0*gd.UNIT_G_PER_CM3
P_mat2 = oscprob.osc_prob_matter_std_potential(3, rho_func, E, L, osc_params=osc)
print('1761 diff', np.abs(P_mat2-P_mat).max(), abs(P_mat2[1,0]-P_mat[1,0]))
H_vac = ham.hamiltonian_3nu_vacuum_energy_independent(s12=osc['s12'], s23=osc['s23'], s13=osc['s13'], dCP=osc['dCP'], D21=osc['D21'], D31=osc['D31'])
proj = matter.matter_potential_projector(3)
P_vac3 = oscprob.osc_prob(H_vac/E, 0.0, L)
P_mat3 = oscprob.osc_prob(H_vac/E + gd.VCC_EARTH_CRUST*proj, 0.0, L)
print('1772 vac diff', np.abs(P_vac3-P_vac).max(), 'mat diff', np.abs(P_mat3-P_mat).max(), P_vac3[1,0], P_mat3[1,0])
P_scan3 = oscprob.osc_prob_energy_baseline(lambda e: H_vac/e + gd.VCC_EARTH_CRUST*proj, Es, L, H_func_is_function_only_of_energy=True)
print('1777', np.shape(P_scan3), np.abs(np.asarray(P_scan3)-P_scanw).max())
P_scan2 = oscprob.osc_prob_matter_std_potential(3, rho_func, Es, L, osc_params=osc)
print('scan scenario diff', np.abs(P_scan2-P_scanw).max())
print('VCC_EARTH_CRUST', gd.VCC_EARTH_CRUST, 'computed', gd.UNIT_G_PER_CM3*3.0 if False else '')
# NSI wrapper
P = oscprob.osc_prob_3nu_matter_nsi_constant_density(E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc, eps_em=0.05)
Pz = oscprob.osc_prob_3nu_matter_nsi_constant_density(E, L, rho=3.0, density_matter_is_in_g_per_cm3=True, **osc)
print('1841', P[1,0], 'no eps', Pz[1,0], 'bitwise eq std', np.array_equal(Pz, P_mat))
nsi = dict(eps_ee=0.0, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
Ps = oscprob.osc_prob_matter_nsi(3, rho_func, E, L, osc_params=osc, nsi_params=nsi)
print('1886', Ps[1,0])
H_vacE = ham.hamiltonian_3nu_vacuum(E, s12=osc['s12'], s23=osc['s23'], s13=osc['s13'], dCP=osc['dCP'], D21=osc['D21'], D31=osc['D31'])
print('1907 proj', proj)
V = gd.VCC_EARTH_CRUST; print('1910 V', V)
H_nsi = ham.hamiltonian_3nu_nsi(V, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0)
H = H_vacE + V*proj + H_nsi
P = oscprob.osc_prob(H, 0.0, L); print('1918', P[1,0])
P = oscprob.osc_prob(H_vacE + H_nsi, 0.0, L); print('1920 no Vproj', P[1,0], P.sum(0), P.sum(1))
