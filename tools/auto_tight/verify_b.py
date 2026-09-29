"""Step 4b.  For each case and tolerance: the hybrid (strategy='hybrid') and strategy='auto', with
P_ee, certification, warnings; and (with --refs) DOP853 at 1e-12 and 1e-13 on the Hamiltonian the
wrapper hands the hybrid.  Run once against the fix (with refs) and once against the HEAD export.
Issue #120's tight_measure.py imports its case list and callables from here; b_fix.jsonl is the
output of the run with refs."""
import sys, json, time, warnings, numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp

TOLS = (1e-7, 1e-9, 1e-12)
WANT_REFS = '--refs' in sys.argv
OUT = sys.argv[1] if len(sys.argv) > 1 else None

def cases():
    import magnus.globaldefs as gd, magnus.oscprob as op
    osc = gd.load_nufit_params('NuFIT 6.1')
    s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
    c = []
    def expw(fn, L, h, rho, extra):
        return lambda E, **k: op.__dict__[fn](E, L=L*gd.UNIT_KM, L0=0.0, rho_central=rho, l_scale=h*gd.UNIT_KM,
                                              density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **extra, **k)
    two = dict(sth=osc['s12'], Dm2=osc['D21'])
    s4 = dict(osc, s14=s14, s24=s24, D41=1.0); s5 = dict(osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7)
    for name, fn, lo, hi, key in (('L1-2nu', 'osc_prob_2nu_matter_exp_density', 0.0005, 0.05, 'two'),
                                  ('L1-3nu', 'osc_prob_3nu_matter_exp_density', 0.002, 0.2, 'osc'),
                                  ('L1-4nu', 'osc_prob_4nu_matter_exp_density', 2.0, 20.0, 's4'),
                                  ('L1-5nu', 'osc_prob_5nu_matter_exp_density', 2.0, 20.0, 's5')):
        for E in np.logspace(np.log10(lo), np.log10(hi), 26):
            c.append((name, float(E), (fn, 25, 10, 3e3, key)))
    for L, h in ((25, 10), (250, 100), (2500, 1000), (25000, 10000)):
        for E in (0.002, 0.02, 0.2):
            c.append(('exp%d-3nu' % L, E, ('osc_prob_3nu_matter_exp_density', L, h, 3e3, 'osc')))
    for L, h in ((25, 10), (2500, 1000)):
        c.append(('exp%d-2nu' % L, 0.02, ('osc_prob_2nu_matter_exp_density', L, h, 3e3, 'two')))
        c.append(('exp%d-4nu' % L, 5.0, ('osc_prob_4nu_matter_exp_density', L, h, 3e3, 's4')))
        c.append(('exp%d-5nu' % L, 5.0, ('osc_prob_5nu_matter_exp_density', L, h, 3e3, 's5')))
    c.append(('exp250-nsi', 0.2, ('osc_prob_3nu_matter_nsi_exp_density', 250, 100, 3e3, 'nsi')))
    c.append(('exp250-liv', 0.2, ('osc_prob_3nu_matter_liv_exp_density', 250, 100, 3e3, 'liv')))
    for E in (0.01, 0.05):
        c.append(('multi-res', E, ('multi',)))
    for prof in ('exp', 'B16-GS98'):
        for frac in (0.1, 0.2, 0.3):
            for E in (0.005, 0.01):
                c.append(('sun2-%s-%.1fR' % (prof, frac), E, ('sun2', prof, frac)))
        c.append(('sun3-%s-0.1R' % prof, 0.01, ('sun3', prof, 0.1)))
    return c

def call_for(spec):
    import magnus.globaldefs as gd, magnus.oscprob as op
    osc = gd.load_nufit_params('NuFIT 6.1')
    s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
    extras = dict(osc=osc, two=dict(sth=osc['s12'], Dm2=osc['D21']), s4=dict(osc, s14=s14, s24=s24, D41=1.0),
                  s5=dict(osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7),
                  nsi=dict(osc, eps_ee=0.15, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0),
                  liv=dict(osc, sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0, b1=gd.B1, b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1))
    if spec[0] == 'multi':
        NE0, LS = gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN
        def ne(l):
            x = np.asarray(l, dtype=float); o = NE0*np.exp(-x/LS)*(1.0 + 0.9*np.sin(2.0*np.pi*6.0*x/LS)); return o[()] if o.ndim == 0 else o
        return (lambda E, **k: op.osc_prob_matter_std_potential(2, ne, np.array([E]), LS, {'sth': 0.55, 'Dm2': 7.5e-5}, L0=0.0,
                density_is_of_number_of_electrons=True, nu_i=gd.NUE, nu_f=gd.NUE, **k)), LS
    if spec[0] in ('sun2', 'sun3'):
        _, prof, frac = spec; R = frac*gd.SUN_RADIUS*gd.UNIT_KM
        ex = {} if prof == 'exp' else dict(density_profile=prof)
        if spec[0] == 'sun2':
            return (lambda E, **k: op.osc_prob_2nu_sun(np.array([E]), R, 0.0, sth=osc['s12'], Dm2=osc['D21'], nu_i=gd.NUE, nu_f=gd.NUE, **ex, **k)), R
        return (lambda E, **k: op.osc_prob_3nu_sun(np.array([E]), R, 0.0, **osc, nu_i=gd.NUE, nu_f=gd.NUE, **ex, **k)), R
    fn, L, h, rho, key = spec
    ex = extras[key]
    return (lambda E, **k: op.__dict__[fn](np.array([E]), L=L*gd.UNIT_KM, L0=0.0, rho_central=rho, l_scale=h*gd.UNIT_KM,
            density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, **ex, **k)), L*gd.UNIT_KM

def run(args):
    name, E_gev, spec = args
    warnings.simplefilter('ignore')
    import magnus.globaldefs as gd, magnus.adiabatic as ad
    f, L = call_for(spec)
    E = E_gev*gd.UNIT_GEV
    cap = {}
    real = ad.hybrid_propagator
    def hp(H_func, *a, **k):
        cap['H'] = H_func; return real(H_func, *a, **k)
    ad.hybrid_propagator = hp
    rec = dict(name=name, E_gev=E_gev, rows=[])
    for tol in TOLS:
        row = dict(tol=tol)
        for strat in ('hybrid', 'auto'):
            info = {}
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                t = time.perf_counter()
                try:
                    P = float(np.asarray(f(E, rtol=tol, atol=tol*1e-2, strategy=strat, strategy_info=info)).ravel()[0])
                except Exception as exc:
                    P = float('nan'); info['error'] = str(exc)[:120]
                dt = time.perf_counter() - t
            row[strat] = dict(P=P, t=dt, engine=info.get('engine'), certified=info.get('certified'),
                              warnings=sorted({x.category.__name__ for x in w}), error=info.get('error'))
        rec['rows'].append(row)
    if WANT_REFS and 'H' in cap:
        Hf = cap['H']; d = np.asarray(Hf(0.0)).shape[-1]
        ls = np.linspace(0.0, L, 257)
        ev = np.linalg.eigvalsh(np.array([np.asarray(Hf(l), dtype=complex) for l in ls]))
        phase = float(np.trapz(ev[:, -1] - ev[:, 0], ls)); rec['phase'] = phase
        if phase <= 3.0e4:
            def rhs(x, y):
                U = y.reshape(d, d, 2); dU = -1j*np.asarray(Hf(x), dtype=complex) @ (U[..., 0] + 1j*U[..., 1])
                return np.stack([dU.real, dU.imag], -1).ravel()
            refs = []
            for rt, at in ((1e-12, 1e-14), (1e-13, 1e-15)):
                s = solve_ivp(rhs, (0.0, L), np.stack([np.eye(d), np.zeros((d, d))], -1).ravel(), method='DOP853', rtol=rt, atol=at)
                U = s.y[:, -1].reshape(d, d, 2); U = U[..., 0] + 1j*U[..., 1]
                refs.append(float(abs(U[0, 0])**2))
            rec['ref12'], rec['ref13'] = refs
    return rec

if __name__ == '__main__':
    import magnus.oscprob as op
    print('magnus from', op.__file__, flush=True)
    with Pool(int(sys.argv[2])) as pool, open(OUT, 'w') as out:
        for rec in pool.imap_unordered(run, cases()):
            out.write(json.dumps(rec) + '\n'); out.flush()
    print('done', flush=True)
