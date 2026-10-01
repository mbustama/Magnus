"""#125 on the physical population (docs/dev/adversarial_batteries/physical_profiles.py): baseline
scans of 2/4/7 points at tight tolerances, auto (main) vs auto with HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS=2.
Records use b125.py's schema; score with `b125.py score`.
  python b125p.py OUT.jsonl NPROC"""
import sys, os, json, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'adversarial_batteries'))
import b125
PHASE_MAX = 6.0e4 if "--sn-high" in sys.argv else 3.0e4


def jobs():
    import physical_profiles as pp
    out = []
    keep = ('tabulated linear N=50', 'tabulated linear N=5000', 'BS05(AGS,OP) linear', 'SN shock w=1e-02',
            'SN shock w=1e-06', 'SN turbulence C*=0.1', 'Earth crust costhz=-1.0')
    for i, f in enumerate(pp.families()):
        if f['label'] not in keep: continue
        if '--sn-high' in sys.argv:
            if not f['kind'].startswith('sn_'): continue
            E = max(f['energies'])
        else:
            E = min(f['energies'])
        for n in (2, 7):
            for tol in (1e-7, 1e-12):
                out.append((i, f['label'], E, n, tol))
    return out


def run_one(args):
    i, label, E, n, tol_only = args
    b125.TOLS = (tol_only,); b125.ORDERS = (4,)
    warnings.simplefilter('ignore')
    import physical_profiles as pp, magnus.oscprob as op, magnus.globaldefs as gd
    f = pp.families()[i]
    osc = {k: v for k, v in gd.load_nufit_params('NuFIT 6.1').items() if k in ('s12','s23','s13','dCP','D21','D31')}
    ne, l0, l1 = f['ne'], f['l0'], f['l1']
    call = lambda EE, LL, **k: op.osc_prob_matter_std_potential(3, ne, EE, LL, osc, L0=l0,
                         density_is_of_number_of_electrons=True, nu_i=gd.NUE, nu_f=gd.NUE, **k)
    Ls = l0 + (l1 - l0)*np.arange(1, n + 1)/n
    rec = dict(kind='bl', name=label, E_gev=[E/gd.UNIT_GEV], n=n, rows=[])
    H = b125.capture_H(call, E, l1)
    if H is None:
        rec['error'] = 'no H captured'; return rec
    ls = np.linspace(l0, l1, 513)
    ev = np.linalg.eigvalsh(np.array([np.asarray(H(l), dtype=complex) for l in ls]))
    phase = float(np.trapezoid(ev[:, -1] - ev[:, 0], ls)); rec['phase'] = phase
    if phase > PHASE_MAX:
        rec['error'] = 'phase %.3g above %g' % (phase, PHASE_MAX); return rec
    # reference from l0: shift the variable
    Hs = lambda x: H(l0 + x)
    ref = b125.dop(Hs, Ls - l0)
    rec['ref12'], rec['ref13'] = ref
    real_min = op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS
    for tol in b125.TOLS:
        for order in b125.ORDERS:
            kw = dict(rtol=tol, atol=tol*1e-2, magnus_exp_order=order, strategy='auto')
            row = dict(tol=tol, order=order)
            op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = real_min
            row['auto'] = b125.timed(call, E, Ls, kw)
            op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 2
            row['prop'] = b125.timed(call, E, Ls, kw)
            op.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = real_min
            rec['rows'].append(row)
    return rec


if __name__ == '__main__':
    from multiprocessing import Pool
    b125.TOLS = (1e-7, 1e-12); b125.ORDERS = (4,)
    out = sys.argv[1]; J = jobs()
    if '--smoke' in sys.argv:
        b125.TOLS = (1e-9,); b125.ORDERS = (4,); J = [j for j in J if j[3] == 2][::12]
    print(len(J), 'jobs', flush=True)
    with Pool(int(sys.argv[2])) as pool, open(out, 'w') as fh:
        for rec in pool.imap_unordered(run_one, J):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print('done', rec['name'], rec['E_gev'][0], rec['n'], [x['tol'] for x in rec.get('rows', [])], rec.get('error', ''), flush=True)
    print('all done', flush=True)
