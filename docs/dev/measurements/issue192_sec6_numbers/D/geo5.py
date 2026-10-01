import json, time, numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.earth as earth
trapz = np.trapezoid
OSC = gd.load_nufit_params('NuFIT 6.1')
D_DET, D_LOCAL = 1.4, 10.0
KW_GEO = dict(nubar=True, nu_i=gd.NUE, nu_f=gd.NUE, density_matter_ocean=2.65)
R_E = gd.EARTH_RADIUS
t0=time.time()
# matter average along three chords
E_GEO_AVG = 1.0/np.linspace(1.0/1.8, 1.0/3.3, 301)*gd.UNIT_MEV
Pv = float(oscprob.osc_prob_3nu_vacuum(E_GEO_AVG[0], 1.0e8*gd.UNIT_KM, average=True, nubar=True, nu_i=gd.NUE, nu_f=gd.NUE, **OSC))
for depth, costhz in [(20.0,-0.272),(1000.0,-0.552),(2800.0,-0.872)]:
    geo = dict(source_depth=depth, detector_depth=D_DET)
    L = earth.distance_traveled_inside_earth(costhz, **geo)
    def rho(l):
        r = earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM, **geo)
        return earth.density_matter_func_prem(r, density_matter_ocean=2.65)
    v = np.asarray(oscprob.osc_prob_matter_std_potential(3, rho, E_GEO_AVG, L*gd.UNIT_KM, OSC, average=True, nubar=True, nu_i=gd.NUE,
        nu_f=gd.NUE, electron_fraction=earth.Y_E_MANTLE_PREM, density_matter_is_in_g_per_cm3=True))
    print(depth, 'raise %+.3f%% .. %+.3f%%' % (100*(v.min()/Pv-1), 100*(v.max()/Pv-1)), flush=True)
print('t', time.time()-t0, flush=True)
_RS, _RD = R_E - D_DET, R_E - D_LOCAL
def local_L_of_c(c): return np.sqrt(_RS**2 - (_RD*np.sqrt(1.0 - c*c))**2) - _RD*c
def local_c_of_L(L):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5*(lo + hi)
        if local_L_of_c(mid) > L: lo = mid
        else: hi = mid
    return 0.5*(lo + hi)
L_LOCAL = np.arange(local_L_of_c(1.0), local_L_of_c(0.0)*0.999, 0.25)
P_LOCAL = np.array([float(oscprob.osc_prob_3nu_earth(2.5*gd.UNIT_MEV, costhz=local_c_of_L(L), source_depth=D_DET*gd.UNIT_KM,
                     detector_depth=D_LOCAL*gd.UNIT_KM, **KW_GEO)) for L in L_LOCAL])
np.save('P_LOCAL.npy', np.array([L_LOCAL, P_LOCAL]))
print('local scan done', len(L_LOCAL), time.time()-t0, flush=True)
N_AV, YR = 6.02214076e23, 3.15576e7
LAM_U, LAM_TH = np.log(2.0)/(4.468e9*YR), np.log(2.0)/(1.405e10*YR)
NU_U, NU_TH, A_U, A_TH = 0.392, 0.147, 238.0, 232.0
CRUST_KM, MANTLE_KM = 35.0, 2891.0
ABUND = dict(crust=(0.67e-6, 3.0e-6), mantle=(0.0127e-6, 0.0446e-6))
def geo_emissivity(r_km):
    d = R_E - r_km
    U = np.where(d <= CRUST_KM, ABUND['crust'][0], np.where(d <= MANTLE_KM, ABUND['mantle'][0], 0.0))
    Th = np.where(d <= CRUST_KM, ABUND['crust'][1], np.where(d <= MANTLE_KM, ABUND['mantle'][1], 0.0))
    return earth.density_matter_func_prem(r_km)*(U/A_U*NU_U*LAM_U + Th/A_TH*NU_TH*LAM_TH)*N_AV
L_FLUX = np.logspace(0.0, np.log10(2*R_E), 900)
psis = np.linspace(0.0, np.pi, 40001)
DPHI, DPHI_CRUST = np.zeros_like(L_FLUX), np.zeros_like(L_FLUX)
for i, L in enumerate(L_FLUX):
    r = np.sqrt(R_E**2 + L**2 - 2*R_E*L*np.cos(psis)); inside = r <= R_E
    if not inside.any(): continue
    em = geo_emissivity(r[inside]); crust = (R_E - r[inside]) <= CRUST_KM
    w = 0.5*np.sin(psis[inside])
    DPHI[i] = trapz(em*w, psis[inside])*1.0e5
    DPHI_CRUST[i] = trapz(np.where(crust, em, 0.0)*w, psis[inside])*1.0e5
print('flux done', time.time()-t0, flush=True)
L_HORIZON = local_L_of_c(0.0)
def geo_window_average(L):
    cycles = OSC['D31']*L*gd.UNIT_KM/(4*np.pi)*(1/(1.8*gd.UNIT_MEV) - 1/(3.3*gd.UNIT_MEV))
    E = 1.0/np.linspace(1.0/1.8, 1.0/3.3, max(240, int(12*cycles)))*gd.UNIT_MEV
    if local_L_of_c(1.0) <= L <= L_HORIZON:
        P = np.asarray(oscprob.osc_prob_3nu_earth(E, costhz=local_c_of_L(L), source_depth=D_DET*gd.UNIT_KM, detector_depth=D_LOCAL*gd.UNIT_KM, **KW_GEO))
    else:
        P = np.asarray(oscprob.osc_prob_3nu_vacuum(E, L*gd.UNIT_KM, nubar=True, nu_i=gd.NUE, nu_f=gd.NUE, **OSC))
    return float(trapz(P, E)/(E[-1] - E[0]))
P_WINDOW = np.array([geo_window_average(L) for L in L_FLUX])
np.save('flux.npy', np.array([L_FLUX, DPHI, DPHI_CRUST, P_WINDOW]))
DPHI_OSC = DPHI*P_WINDOW
TOTAL, TOTAL_OSC = trapz(DPHI, L_FLUX), trapz(DPHI_OSC, L_FLUX)
CUM = np.array([trapz(DPHI[:i + 1], L_FLUX[:i + 1]) for i in range(len(L_FLUX))])/TOTAL
_at = lambda x: int(np.searchsorted(L_FLUX, x))
print('crust %.2f%% within100 %.2f%% 350 %.2f%% 1000 %.2f%% osc leave %.2f%%' % (100*trapz(DPHI_CRUST, L_FLUX)/TOTAL, 100*CUM[_at(100)],
      100*CUM[_at(350)], 100*CUM[_at(1000)], 100*TOTAL_OSC/TOTAL))
print('horizon %.2f km  %.4f crust, %.4f vac' % (L_HORIZON, P_WINDOW[_at(L_HORIZON) - 1], P_WINDOW[_at(L_HORIZON)]))
print('t', time.time()-t0)
