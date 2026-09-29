"""#70 campaign: hybrid (strategy='auto' today) against the ladder at a tenth of the tolerance
(strategy='magnus', rtol/10, atol/10), by the phase estimate the dispatcher would compute."""
import json, sys, time, warnings, numpy as np
from scipy.integrate import solve_ivp
warnings.simplefilter('ignore')  # per call below
import magnus.globaldefs as gd, magnus.oscprob as op, magnus.magnus as mm
osc = gd.load_nufit_params('NuFIT 6.1')
OUT = sys.argv[1]

# Capture the H(E, l) the wrapper hands the hybrid engine, without changing what it returns.
CAP = {}
_orig = op._hybrid_propagator_scan
def _spy(H_at_energy, energy_arr, L_arr, L0, *a, **k):
    CAP['H'] = H_at_energy; CAP['L0'] = L0
    return _orig(H_at_energy, energy_arr, L_arr, L0, *a, **k)
op._hybrid_propagator_scan = _spy

def phase_hat(H_at_energy, E, L0, L1, n_probe=17):
    """The estimate the dispatcher can afford: norm of the traceless integral of H over the path."""
    ls = np.linspace(L0, L1, n_probe)
    Hs = np.array([np.asarray(H_at_energy(E)(l), dtype=complex) for l in ls])
    M = np.trapz(Hs, ls, axis=0); d = M.shape[-1]; M = M - np.trace(M)/d*np.eye(d)
    return float(np.max(np.linalg.svd(M, compute_uv=False)))

def dop853(H_at_energy, E, L0, L1):
    Hf = H_at_energy(E); d = np.asarray(Hf(L0)).shape[-1]
    def rhs(x, y):
        U = y.reshape(d, d, 2); dU = -1j*np.asarray(Hf(x), dtype=complex) @ (U[..., 0] + 1j*U[..., 1])
        return np.stack([dU.real, dU.imag], -1).ravel()
    s = solve_ivp(rhs, (L0, L1), np.stack([np.eye(d), np.zeros((d, d))], -1).ravel(), method='DOP853', rtol=1e-12, atol=1e-12)
    U = s.y[:, -1].reshape(d, d, 2); U = U[..., 0] + 1j*U[..., 1]; return np.abs(U.T)**2

def best_time(f, n=2):
    f(); b = np.inf
    for _ in range(n):
        t = time.perf_counter(); r = f(); b = min(b, time.perf_counter() - t)
    return b, r

S14 = dict(s14=np.sqrt(0.1), s24=np.sqrt(0.1), D41=1.0)
S15 = dict(s15=np.sqrt(0.06), s25=np.sqrt(0.06), D51=1.7)
NSI = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
LIV = dict(sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0, b1=gd.B1, b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1)

def exp_call(d, L_km, h_km, rho0, extra=None, fam='std'):
    kw = dict(L=L_km*gd.UNIT_KM, L0=0.0, rho_central=rho0, l_scale=h_km*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True)
    fn = {('std', 2): op.osc_prob_2nu_matter_exp_density, ('std', 3): op.osc_prob_3nu_matter_exp_density,
          ('std', 4): op.osc_prob_4nu_matter_exp_density, ('std', 5): op.osc_prob_5nu_matter_exp_density,
          ('nsi', 3): op.osc_prob_3nu_matter_nsi_exp_density, ('liv', 3): op.osc_prob_3nu_matter_liv_exp_density}[(fam, d)]
    pars = dict(sth=osc['s12'], Dm2=osc['D21']) if d == 2 else dict(osc)
    if d >= 4: pars.update(S14)
    if d == 5: pars.update(S15)
    if fam == 'nsi': pars.update(NSI)
    if fam == 'liv': pars.update(LIV)
    return (lambda E, **k: np.asarray(fn(E, **kw, **pars, **k))), 0.0, L_km*gd.UNIT_KM

def sun_call(profile):
    R = gd.SUN_RADIUS*gd.UNIT_KM
    extra = {} if profile == 'exp' else dict(density_profile=profile)
    return (lambda E, **k: np.asarray(op.osc_prob_3nu_sun(E, R, 0.0, **osc, **extra, **k))), 0.0, R

CASES = []
for L_km, h_km, rho0 in ((25, 10, 3e3), (250, 100, 3e3), (2500, 1000, 3e3), (25000, 10000, 3e3), (250000, 60000, 100)):
    for E_gev in (0.002, 0.02, 0.2):
        CASES.append(('exp L=%g d=3' % L_km, lambda L_km=L_km, h_km=h_km, rho0=rho0: exp_call(3, L_km, h_km, rho0), E_gev))
for d in (2, 4, 5):
    CASES.append(('exp L=25 d=%d' % d, lambda d=d: exp_call(d, 25, 10, 3e3), 0.02 if d == 2 else 5.0))
    CASES.append(('exp L=2500 d=%d' % d, lambda d=d: exp_call(d, 2500, 1000, 3e3), 0.02 if d == 2 else 5.0))
for fam in ('nsi', 'liv'):
    CASES.append(('exp L=250 %s' % fam, lambda fam=fam: exp_call(3, 250, 100, 3e3, fam=fam), 0.2))
for prof in ('exp', 'B16-GS98'):
    for E_mev in (1.0, 5.0, 15.0):
        CASES.append(('sun %s' % prof, lambda prof=prof: sun_call(prof), E_mev*1e-3))


with open(OUT, 'a') as out:
    for name, make, E_gev in CASES:
        f, L0, L1 = make()
        E1 = E_gev*gd.UNIT_GEV
        Es = np.logspace(np.log10(E_gev*0.7), np.log10(E_gev*1.4), 40)*gd.UNIT_GEV
        CAP.clear(); f(E1, strategy='hybrid', rtol=1e-3, atol=1e-3)
        H = CAP.get('H')
        ph = max(phase_hat(H, E, L0, L1) for E in Es[::13]) if H is not None else None
        use_dop = (H is not None) and (ph < 3e4)
        pick = [E1] + list(Es[::10])
        ref = {}
        for tol in (1e-3, 1e-6):
            rec = dict(case=name, E_gev=E_gev, phase=ph, tol=tol)
            try:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    info1 = {}; P1 = f(E1, rtol=tol, atol=tol, strategy_info=info1)
                    infos = {}; Ps = f(Es, rtol=tol, atol=tol, strategy_info=infos)
                wn = sorted(set(x.category.__name__ for x in w))
                t1, _ = best_time(lambda: f(E1, rtol=tol, atol=tol))
                ts, _ = best_time(lambda: f(Es, rtol=tol, atol=tol), 1)
                prefer = any(tr.get('reason') == 'auto prefers the ladder' for tr in infos.get('trace', []))
                if not ref:
                    for E in pick:
                        ref[float(E)] = dop853(H, E, L0, L1) if use_dop else f(E, strategy='hybrid', rtol=1e-9, atol=1e-9)
                e = [np.max(np.abs(P1 - ref[float(E1)]))]
                e += [np.max(np.abs(Ps[i] - ref[float(Es[i])])) for i in range(0, 40, 10)]
                rec.update(engine_1=info1.get('engine'), engine_scan=infos.get('engine'), prefer=prefer,
                           warnings=wn, t_1=t1, t_scan=ts/40, err=float(max(e)),
                           ref='dop853' if use_dop else 'hybrid-1e-9')
            except Exception as exc:
                rec.update(error=str(exc)[:200])
            out.write(json.dumps(rec) + '\n'); out.flush()
            print(json.dumps(rec), flush=True)
