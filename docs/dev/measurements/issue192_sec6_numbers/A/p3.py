import numpy as np, warnings
import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.matter as matter
import magnus.oscprob as oscprob
osc = gd.load_nufit_params('NuFIT 6.1')
E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
H_vac = ham.hamiltonian_3nu_vacuum_energy_independent(s12=osc['s12'], s23=osc['s23'], s13=osc['s13'], dCP=osc['dCP'], D21=osc['D21'], D31=osc['D31'])
proj = matter.matter_potential_projector(3)
P, U = oscprob.osc_prob_3nu_vacuum(E, L, return_evolution_operator=True, **osc)
print('1954', np.allclose(P, abs(U.T)**2), np.allclose(U.conj().T @ U, np.eye(3)))
H1 = H_vac/E + gd.VCC_EARTH_CRUST*proj
H2 = H_vac/E + 1.5*gd.VCC_EARTH_CRUST*proj
_, U1 = oscprob.osc_prob(H1, 0.0, L, return_evolution_operator=True)
_, U2 = oscprob.osc_prob(H2, 0.0, L, return_evolution_operator=True)
Pc = abs((U2@U1).T)**2
Hf = lambda l: H1 if l < L else H2
Pb = oscprob.osc_prob(Hf, 0.0, 2*L, t_breakpoints=[L])
print('1972 compose vs breakpoints', np.abs(Pc-Pb).max())
FAR = 1.0e8*gd.UNIT_KM
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    P = oscprob.osc_prob_3nu_vacuum(E, FAR, average=True, **osc)
    print('1992', P[0,0], P[1,0], [str(x.message)[:80] for x in w])
R = ham.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP'])
print('   sum|Uei|^2|Umi|^2', np.sum(abs(R[0])**2*abs(R[1])**2), np.sum(abs(R[0])**4))
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    P = oscprob.osc_prob_3nu_vacuum(E, L, average=True, **osc)
    print('2003', P[0,0], P[1,1], [type(x.message).__name__+': '+str(x.message)[:200] for x in w])
P = oscprob.osc_prob_3nu_vacuum(E, L, **osc); print('2004', P[0,0], P[1,1])
print('2001 atm phase rad', osc['D31']*L/(2*E))
for start in ('flavor', 'decohered'):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = oscprob.osc_prob_3nu_sun(10.0*gd.UNIT_MEV, gd.SUN_RADIUS*gd.UNIT_KM, 0.0, nu_i=gd.NUE, nu_f=gd.NUE, average=True, average_initial_state=start)
    print('2019', start, repr(P), [type(x.message).__name__ for x in w])
H = H_vac/E + gd.VCC_EARTH_CRUST*proj
P = oscprob.osc_prob_energy_baseline(H, E, FAR, average=True)
print('2038', P[0,0])
