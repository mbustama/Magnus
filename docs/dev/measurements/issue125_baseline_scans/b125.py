"""Issue #125: strategy='auto' at tight tolerances on baseline scans (hybrid vs cumulative) and
on energy scans (auto vs the energy-batched engine).  Scored against DOP853.

  python b125.py run OUT.jsonl NPROC [--smoke]
  python b125.py score OUT.jsonl
"""
import sys, json, time, warnings
import numpy as np
from scipy.integrate import solve_ivp

TOLS = (1e-7, 1e-9, 1e-12)            # rtol; atol = rtol/100 (as #120's B1)
ORDERS = (4, 8)
NPTS = (2, 4, 7)
HONEST = {'HybridCertificationWarning', 'ToleranceNotAchievedWarning', 'UnmarkedDiscontinuityWarning'}


def extras():
    import magnus.globaldefs as gd
    osc = gd.load_nufit_params('NuFIT 6.1')
    s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
    return dict(osc=osc, two=dict(sth=osc['s12'], Dm2=osc['D21']),
                s4=dict(osc, s14=s14, s24=s24, D41=1.0),
                s5=dict(osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7),
                nsi=dict(osc, eps_ee=0.15, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0),
                liv=dict(osc, sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0, b1=gd.B1, b2=gd.B2, b3=gd.B3,
                         Lambda=gd.LAMBDA, n_liv=1))


def caller(spec):
    """f(E, L, **kw) for the workload, and its full length [eV^-1]."""
    import magnus.globaldefs as gd, magnus.oscprob as op
    ex = extras(); osc = ex['osc']
    if spec[0] == 'multi':
        NE0, LS = gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN
        def ne(l):
            x = np.asarray(l, dtype=float)
            o = NE0*np.exp(-x/LS)*(1.0 + 0.9*np.sin(2.0*np.pi*6.0*x/LS))
            return o[()] if o.ndim == 0 else o
        return (lambda E, L, **k: op.osc_prob_matter_std_potential(
            2, ne, E, L, {'sth': 0.55, 'Dm2': 7.5e-5}, L0=0.0, density_is_of_number_of_electrons=True,
            nu_i=gd.NUE, nu_f=gd.NUE, **k)), LS
    if spec[0] in ('sun2', 'sun3'):
        _, prof, frac = spec
        R = frac*gd.SUN_RADIUS*gd.UNIT_KM
        e2 = {} if prof == 'exp' else dict(density_profile=prof)
        if spec[0] == 'sun2':
            return (lambda E, L, **k: op.osc_prob_2nu_sun(E, L, 0.0, sth=osc['s12'], Dm2=osc['D21'],
                                                          nu_i=gd.NUE, nu_f=gd.NUE, **e2, **k)), R
        return (lambda E, L, **k: op.osc_prob_3nu_sun(E, L, 0.0, **osc, nu_i=gd.NUE, nu_f=gd.NUE, **e2, **k)), R
    fn, L, h, rho, key = spec
    e = ex[key]
    return (lambda E, LL, **k: getattr(op, fn)(E, L=LL, L0=0.0, rho_central=rho, l_scale=h*gd.UNIT_KM,
                                               density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE,
                                               **e, **k)), L*gd.UNIT_KM


def workloads():
    """(kind, name, spec, energies_gev): kind 'bl' = one energy, baselines; 'en' = energy scan at L."""
    w = []
    def ex_(fn, L, h, key): return (fn, L, h, 3e3, key)
    for L, h in ((25, 10), (250, 100), (2500, 1000), (25000, 10000)):
        for E in (0.002, 0.02, 0.2):
            w.append(('bl', 'exp%d-3nu' % L, ex_('osc_prob_3nu_matter_exp_density', L, h, 'osc'), [E]))
    for L, h in ((25, 10), (2500, 1000)):
        w.append(('bl', 'exp%d-2nu' % L, ex_('osc_prob_2nu_matter_exp_density', L, h, 'two'), [0.02]))
        w.append(('bl', 'exp%d-4nu' % L, ex_('osc_prob_4nu_matter_exp_density', L, h, 's4'), [5.0]))
        w.append(('bl', 'exp%d-5nu' % L, ex_('osc_prob_5nu_matter_exp_density', L, h, 's5'), [5.0]))
    w.append(('bl', 'exp250-nsi', ex_('osc_prob_3nu_matter_nsi_exp_density', 250, 100, 'nsi'), [0.2]))
    w.append(('bl', 'exp250-liv', ex_('osc_prob_3nu_matter_liv_exp_density', 250, 100, 'liv'), [0.2]))
    for E in (0.01, 0.05):
        w.append(('bl', 'multi-res', ('multi',), [E]))
    for prof in ('exp', 'B16-GS98'):
        for frac in (0.1, 0.3, 1.0):
            for E in (0.005, 0.01):
                w.append(('bl', 'sun2-%s-%.1fR' % (prof, frac), ('sun2', prof, frac), [E]))
                w.append(('bl', 'sun3-%s-%.1fR' % (prof, frac), ('sun3', prof, frac), [E]))
    # energy scans, 16 energies each
    def en(name, spec, lo, hi): w.append(('en', name, spec, list(np.logspace(np.log10(lo), np.log10(hi), 16))))
    en('L1-3nu', ex_('osc_prob_3nu_matter_exp_density', 25, 10, 'osc'), 0.002, 0.2)
    en('L1-4nu', ex_('osc_prob_4nu_matter_exp_density', 25, 10, 's4'), 2.0, 20.0)
    en('exp250-3nu', ex_('osc_prob_3nu_matter_exp_density', 250, 100, 'osc'), 0.002, 0.2)
    en('exp2500-3nu', ex_('osc_prob_3nu_matter_exp_density', 2500, 1000, 'osc'), 0.002, 0.2)
    en('exp25000-3nu', ex_('osc_prob_3nu_matter_exp_density', 25000, 10000, 'osc'), 0.002, 0.2)
    en('exp2500-4nu', ex_('osc_prob_4nu_matter_exp_density', 2500, 1000, 's4'), 2.0, 20.0)
    en('multi-res', ('multi',), 0.01, 0.05)
    for prof in ('exp', 'B16-GS98'):
        for frac in (0.1, 0.3):
            en('sun2-%s-%.1fR' % (prof, frac), ('sun2', prof, frac), 0.005, 0.01)
            en('sun3-%s-%.1fR' % (prof, frac), ('sun3', prof, frac), 0.005, 0.01)
    return w


def capture_H(f, E, L):
    """The Hamiltonian the wrapper hands the hybrid at one energy, from a cheap hybrid call."""
    import magnus.adiabatic as ad
    cap = {}; real = ad.hybrid_propagator
    def hp(H, *a, **k):
        cap['H'] = H; return real(H, *a, **k)
    ad.hybrid_propagator = hp
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            f(np.array([E]), L, rtol=1e-3, atol=1e-3, strategy='hybrid')
    finally:
        ad.hybrid_propagator = real
    return cap.get('H')


def dop(H, Ls):
    d = np.asarray(H(0.0)).shape[-1]
    def rhs(x, y):
        U = y.reshape(d, d, 2); dU = -1j*np.asarray(H(x), dtype=complex) @ (U[..., 0] + 1j*U[..., 1])
        return np.stack([dU.real, dU.imag], -1).ravel()
    out = []
    for rt, at in ((1e-12, 1e-14), (1e-13, 1e-15)):
        s = solve_ivp(rhs, (0.0, float(np.max(Ls))), np.stack([np.eye(d), np.zeros((d, d))], -1).ravel(),
                      method='DOP853', rtol=rt, atol=at, t_eval=np.sort(Ls))
        U = s.y.T.reshape(-1, d, d, 2); U = U[..., 0] + 1j*U[..., 1]
        order = np.argsort(np.argsort(Ls))
        out.append([float(abs(U[i, 0, 0])**2) for i in order])
    return out


def timed(f, E, L, kw, reps=2):
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        t = time.perf_counter()
        try:
            P = np.asarray(f(E, L, strategy_info=info, **kw), dtype=float).ravel()
        except Exception as exc:
            return dict(error='%s: %s' % (type(exc).__name__, str(exc)[:150]))
        cold = time.perf_counter() - t
    best = cold
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(reps):
            if best > 3.0:
                break
            t = time.perf_counter(); f(E, L, **kw); best = min(best, time.perf_counter() - t)
    return dict(P=P.tolist(), t=best, engine=info.get('engine'),
                warnings=sorted({x.category.__name__ for x in w}))


def run_one(args):
    kind, name, spec, Es, n = args
    warnings.simplefilter('ignore')
    import magnus.oscprob as op, magnus.globaldefs as gd
    f, Lfull = caller(spec)
    E = np.asarray(Es)*gd.UNIT_GEV
    if kind == 'bl':
        Ls = Lfull*np.arange(1, n + 1)/n
        H = capture_H(f, E[0], Lfull)
        if H is None:
            return dict(kind=kind, name=name, E_gev=Es, n=n, error='no H captured')
        ref = dop(H, Ls)
        Larg = Ls
    else:
        ref12, ref13 = [], []
        for e in E:
            H = capture_H(f, e, Lfull)
            r = dop(H, np.array([Lfull])); ref12.append(r[0][0]); ref13.append(r[1][0])
        ref = [ref12, ref13]; Larg = Lfull
    rec = dict(kind=kind, name=name, E_gev=Es, n=n, ref12=ref[0], ref13=ref[1], rows=[])
    real_min = op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS
    for tol in TOLS:
        for order in ORDERS:
            row = dict(tol=tol, order=order)
            kw = dict(rtol=tol, atol=tol*1e-2, magnus_exp_order=order)
            if kind == 'bl':
                op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = real_min
                row['auto'] = timed(f, E, Larg, dict(kw, strategy='auto'))
                op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 2       # #125: baseline scans -> cumulative
                row['prop'] = timed(f, E, Larg, dict(kw, strategy='auto'))
                op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = real_min
            else:
                row['auto'] = timed(f, E, Larg, dict(kw, strategy='auto'))
                row['prop'] = timed(f, E, Larg, dict(kw, strategy='magnus'))   # energy-batched
            rec['rows'].append(row)
    return rec


def score(path):
    recs = [json.loads(l) for l in open(path)]
    bad = [r for r in recs if 'error' in r]
    print('%d records (%d without reference)' % (len(recs), len(bad)))
    def el(r, x, tol):
        worst = (0.0, 1.0)
        for P, a, b in zip(x['P'], r['ref12'], r['ref13']):
            lim = max(tol*1e-2 + tol*abs(b), abs(a - b)); e = abs(P - b)
            if e/lim > worst[0]/worst[1]: worst = (e, lim)
        return worst
    for kind, lab in (('bl', 'BASELINE SCANS: auto(main)=hybrid vs prop=cumulative'),
                      ('en', 'ENERGY SCANS: auto(main) vs prop=energy-batched (strategy=magnus)')):
        print('\n==== %s' % lab)
        combos = sorted({(x['tol'], x['order']) for r in recs if r.get('kind') == kind and 'rows' in r for x in r['rows']}, reverse=True)
        for tol, order in combos:
            if True:
                stats = {}
                ratios = []; slow = []
                for r in recs:
                    if r.get('kind') != kind or 'error' in r: continue
                    row = next((x for x in r['rows'] if x['tol'] == tol and x['order'] == order), None)
                    if row is None: continue
                    for c in ('auto', 'prop'):
                        x = row[c]; s = stats.setdefault(c, dict(n=0, silent=0, warned=0, err=0, worst=[] , eng={}))
                        if 'error' in x: s['err'] += 1; continue
                        s['n'] += 1; s['eng'][x['engine']] = s['eng'].get(x['engine'], 0) + 1
                        e, lim = el(r, x, tol)
                        if HONEST & set(x['warnings']):
                            s['warned'] += 1
                            if e <= lim: s['w_in'] = s.get('w_in', 0) + 1
                            if c == 'prop' and not (HONEST & set(row['auto'].get('warnings', []))):
                                s['new_w'] = s.get('new_w', 0) + 1
                                s.setdefault('new_list', []).append((e/lim, r['name'], r['E_gev'][0], r['n'], sorted(HONEST & set(x['warnings']))))
                        elif e > lim: s['silent'] += 1; s['worst'].append((e/lim, r['name'], r['E_gev'][0], r['n']))
                    if 'error' not in row['auto'] and 'error' not in row['prop']:
                        q = row['prop']['t']/row['auto']['t']; ratios.append(q)
                        if q > 1.2: slow.append((q, r['name'], r['E_gev'][0], r['n'], row['auto']['t'], row['prop']['t']))
                rs = np.array(ratios)
                print('-- rtol %g order %d: t(prop)/t(auto) median %.2f  max %.2f  n>1.2: %d of %d'
                      % (tol, order, np.median(rs), rs.max(), int((rs > 1.2).sum()), rs.size))
                for c in ('auto', 'prop'):
                    s = stats[c]
                    print('   %-4s n=%d silent=%d warned=%d (of which within tol %d; new vs auto %d) errors=%d engines=%s' % (c, s['n'], s['silent'], s['warned'], s.get('w_in', 0), s.get('new_w', 0), s['err'], s['eng']))
                    for q in sorted(s.get('new_list', []), reverse=True)[:3]:
                        print('        new warning, err/limit %.2g: %s E=%g n=%s %s' % q)
                    for q in sorted(s['worst'], reverse=True)[:3]:
                        print('        silent %.2fx outside: %s E=%g n=%s' % q)
                for q in sorted(slow, reverse=True)[:4]:
                    print('   slower %.2fx: %s E=%g n=%s (%.3fs -> %.3fs)' % q)


if __name__ == '__main__':
    mode, out = sys.argv[1], sys.argv[2]
    if mode == 'score':
        score(out); sys.exit()
    from multiprocessing import Pool
    import subprocess
    import magnus.oscprob as op
    print('magnus from', op.__file__, subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
          capture_output=True, text=True).stdout.strip(), flush=True)
    jobs = []
    for kind, name, spec, Es in workloads():
        if kind == 'bl':
            jobs += [(kind, name, spec, Es, n) for n in NPTS]
        else:
            jobs.append((kind, name, spec, Es, len(Es)))
    if '--smoke' in sys.argv:
        TOLS = (1e-9,); ORDERS = (4,)
        jobs = [j for j in jobs if j[1] in ('exp250-3nu', 'sun3-exp-0.1R')][:3]
    if '--no-sun' in sys.argv:
        jobs = [j for j in jobs if not j[1].startswith('sun')]
    if '--only-en' in sys.argv:
        jobs = [j for j in jobs if j[0] == 'en']
    mode_w = 'w'
    if '--resume' in sys.argv:
        import os
        seen = set()
        for path in [a.split('=', 1)[1] for a in sys.argv if a.startswith('--skip=')] + ([out] if os.path.exists(out) else []):
            for l in open(path):
                r = json.loads(l); seen.add((r['kind'], r['name'], round(r['E_gev'][0], 12), r.get('n')))
        jobs = [j for j in jobs if (j[0], j[1], round(j[3][0], 12), j[4]) not in seen]
        mode_w = 'a'
    print(len(jobs), 'jobs', flush=True)
    with Pool(int(sys.argv[3])) as pool, open(out, mode_w) as fh:
        for rec in pool.imap_unordered(run_one, jobs):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print('done', rec['kind'], rec['name'], rec['E_gev'][0], rec.get('n'), rec.get('error', ''), flush=True)
    print('all done', flush=True)
