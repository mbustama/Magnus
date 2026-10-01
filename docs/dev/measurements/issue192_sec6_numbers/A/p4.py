import numpy as np, warnings, inspect
warnings.simplefilter('ignore')
import magnus.globaldefs as gd
import magnus.hamiltonians as ham
import magnus.matter as matter
import magnus.oscprob as oscprob
import magnus.earth as earth
osc = gd.load_nufit_params('NuFIT 6.1')
E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
H = ham.hamiltonian_3nu_vacuum(E, **osc); print('2105', H.shape)
st4 = dict(s14=0.1, s24=0.1, s34=0.1, d14=0.0, d24=0.0, D41=1.0)
print('2125', ham.hamiltonian_4nu_vacuum(E, **osc, **st4).shape)
st5 = dict(st4, s15=0.1, s25=0.1, s35=0.1, d15=0.0, d35=0.0, D51=1.7)
print('2131', ham.hamiltonian_5nu_vacuum(E, **osc, **st5).shape)
H1 = ham.hamiltonian_3nu_vacuum(E, **osc)
H0 = ham.hamiltonian_3nu_vacuum_energy_independent(**osc)
print('2160', H1 - H0/E, np.array_equal(H1, H0/E))
V = gd.VCC_EARTH_CRUST; print('2173', V)
print('2175', ham.hamiltonian_3nu_matter(V).shape, np.diag(ham.hamiltonian_3nu_matter(V))/V)
print('2185 2nu', np.diag(ham.hamiltonian_2nu_matter(V))/V)
Ye = earth.Y_E_CORE_PREM; r=(1.0-Ye)/Ye; print('2191', Ye, r)
print('2195', ham.hamiltonian_4nu_matter(V, r).shape, np.diag(ham.hamiltonian_4nu_matter(V, r)).real/V)
print('2198', ham.hamiltonian_5nu_matter(V, r).shape, np.diag(ham.hamiltonian_5nu_matter(V, r)).real/V)
print('4nu matter default r', inspect.signature(ham.hamiltonian_4nu_matter))
f=lambda l: V*(1+l/L)
print('2214', np.array_equal(ham.hamiltonian_3nu_matter_td(0.3*L, f), ham.hamiltonian_3nu_matter(f(0.3*L))))
eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = 0.1, 0.05+0.02j, 0.03, 0.0, 0.04, 0.02
print('2223', ham.hamiltonian_3nu_nsi(V, eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt).shape)
print('2nu nsi', inspect.signature(ham.hamiltonian_2nu_nsi), ham.hamiltonian_2nu_nsi(V, 0.1, 0.2j)/V)
eps_es=eps_ms=eps_ts=eps_ss=0.01
print('2252', ham.hamiltonian_4nu_nsi(V, eps_ee, eps_em, eps_et, eps_es, eps_mm, eps_mt, eps_ms, eps_tt, eps_ts, eps_ss).shape)
eps_es1=eps_es2=eps_ms1=eps_ms2=eps_ts1=eps_ts2=eps_s1s1=eps_s1s2=eps_s2s2=0.01
print('2259', ham.hamiltonian_5nu_nsi(V, eps_ee, eps_em, eps_et, eps_es1, eps_es2, eps_mm, eps_mt, eps_ms1, eps_ms2, eps_tt, eps_ts1, eps_ts2, eps_s1s1, eps_s1s2, eps_s2s2).shape)
for n in ('hamiltonian_3nu_liv','hamiltonian_4nu_liv','hamiltonian_5nu_liv','hamiltonian_2nu_liv'):
    print(n, list(inspect.signature(getattr(ham,n)).parameters))
sxi12=sxi23=sxi13=dxiCP=sxi14=dxi14=sxi24=dxi24=sxi34=sxi15=dxi15=sxi25=sxi35=dxi35=0.1
b1,b2,b3,b4,b5=1e-23,2e-23,3e-23,4e-23,5e-23; Lambda=1.0; n_liv=0
print('2279', ham.hamiltonian_3nu_liv(E, sxi12, sxi23, sxi13, dxiCP, b1, b2, b3, Lambda, n_liv).shape)
print('2301', ham.hamiltonian_4nu_liv(E, sxi12, sxi23, sxi13, dxiCP, sxi14, dxi14, sxi24, dxi24, sxi34, b1, b2, b3, b4, Lambda, n_liv).shape)
print('2308', ham.hamiltonian_5nu_liv(E, sxi12, sxi23, sxi13, dxiCP, sxi14, dxi14, sxi15, dxi15, sxi24, dxi24, sxi25, sxi34, sxi35, dxi35, b1, b2, b3, b4, b5, Lambda, n_liv).shape)
R = ham.pmns_mixing_matrix(osc['s12'], osc['s23'], osc['s13'], osc['dCP'])
H = ham.hamiltonian_pseudo_dirac_vacuum(E, R, [0.0, osc['D21'], osc['D31']], {1: 1.0e-18}); print('2341', H.shape)
H = ham.hamiltonian_pseudo_dirac_vacuum(E, R, [0.0, osc['D21'], osc['D31']], {})
print('2343 empty vs 3nu', np.abs(H - ham.hamiltonian_3nu_vacuum(E, **osc)).max(), np.abs(ham.hamiltonian_3nu_vacuum(E, **osc)).max())
print('PD matter', inspect.signature(ham.hamiltonian_pseudo_dirac_matter))
# building new
H_vac = ham.hamiltonian_3nu_vacuum(E, s12=osc['s12'], s23=osc['s23'], s13=osc['s13'], dCP=osc['dCP'], D21=osc['D21'], D31=osc['D31'])
VCC_func = lambda l: V*(1.0 + l/L)
def H_func(l):
    H_mat = ham.hamiltonian_3nu_matter_td(l, VCC_func)
    return H_vac + H_mat
P = oscprob.osc_prob(H_func, 0.0, L, rtol=1e-8, atol=1e-8); print('2402', P[1,0])
def H_own(l):
    H = H_vac.copy()
    H[0, 2] += 1.0e-13*(l/L)
    H[2, 0] += 1.0e-13*(l/L)
    return H
P = oscprob.osc_prob(H_own, 0.0, L, rtol=1e-8, atol=1e-8); print('2415', P[1,0])
