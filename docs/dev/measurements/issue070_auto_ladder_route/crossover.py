"""#70: per-point cost of strategy='auto' (hybrid) against strategy='magnus' (ladder) across the
accumulated phase, at the default tolerance, with each answer scored against a tight reference."""
import time, warnings, sys, numpy as np
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.oscprob as op, magnus.magnus as mm, magnus.hamiltonians as hm, magnus.matter as matter
osc = gd.load_nufit_params('NuFIT 6.1')

def phase_estimate(Hf, L):
    """What the dispatcher could afford to compute: suggest_n_slabs's norm of the traceless Omega_1."""
    return 2*np.pi*mm.suggest_n_slabs(lambda t: -1j*np.asarray(Hf(t)), 0.0, L)

def exp_case(L_km, rho0, h_km, E_gev):
    L = L_km*gd.UNIT_KM
    kw = dict(L=L, L0=0.0, rho_central=rho0, l_scale=h_km*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE)
    f = lambda E, **k: op.osc_prob_3nu_matter_exp_density(E, **kw, **osc, **k)
    hv = np.asarray(hm.hamiltonian_3nu_vacuum(E_gev*gd.UNIT_GEV, **osc), dtype=complex)
    V0 = matter.VCC_func(0.0, lambda l: 1.0)*rho0*gd.UNIT_G_PER_CM3/(gd.MASS_PROTON)*0.5 if False else None
    return f, L

def timed(f, E, **k):
    f(E, **k)
    best = np.inf
    for _ in range(3):
        t = time.perf_counter(); P = np.asarray(f(E, **k)); best = min(best, time.perf_counter() - t)
    return best, P

rows = []
cases = [(25, 3e3, 10, 0.02), (25, 3e3, 10, 0.2), (250, 3e3, 100, 0.2), (2500, 3e3, 1000, 0.2),
         (2500, 3e3, 1000, 0.02), (25000, 3e3, 10000, 0.02), (25000, 100, 10000, 0.005), (250000, 100, 60000, 0.005)]
for L_km, rho0, h_km, E_gev in cases:
    f, L = exp_case(L_km, rho0, h_km, E_gev)
    E1 = E_gev*gd.UNIT_GEV
    Escan = np.logspace(np.log10(E_gev*0.5), np.log10(E_gev*2), 40)*gd.UNIT_GEV
    info = {}
    f(E1, strategy_info=info)
    eng = info.get('engine')
    # phase estimate from the same Hamiltonian the wrapper builds
    Hf = lambda l: (np.asarray(hm.hamiltonian_3nu_vacuum(E1, **osc), dtype=complex)
                    + op.matter.VCC_func(l, lambda x: rho0*gd.UNIT_G_PER_CM3*np.exp(-x/(h_km*gd.UNIT_KM))/(0.5*(gd.MASS_PROTON+gd.MASS_NEUTRON))*0.5)*np.diag([1, 0, 0]))
    phi = phase_estimate(Hf, L)
    ta1, Pa1 = timed(f, E1); tm1, Pm1 = timed(f, E1, strategy='magnus')
    tas, Pas = timed(f, Escan); tms, Pms = timed(f, Escan, strategy='magnus')
    Pref = np.asarray(f(Escan, strategy='magnus', rtol=1e-9, atol=1e-9))
    rows.append((L_km, E_gev, phi, eng, ta1, tm1, tas/len(Escan), tms/len(Escan), np.max(np.abs(Pas - Pref)), np.max(np.abs(Pms - Pref))))
    print('L=%7g km E=%6g GeV phase~%9.0f rad engine(auto)=%-8s | 1 point: auto %8.4f s magnus %8.4f s | per point in a 40-scan: auto %8.5f s magnus %8.5f s | err auto %.1e magnus %.1e'
          % rows[-1], flush=True)
