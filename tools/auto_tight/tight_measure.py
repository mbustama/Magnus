"""Issue #120's measurement (step B1): strategy='auto' at tight tolerances, ladder against hybrid.

For every step-4b case: at a tight tolerance, is strategy='auto' faster if it hands the request to
the Magnus ladder (the #70 route with its 1e-6 cut-off lifted) than on the fixed hybrid, and does the
ladder stay inside the tolerance?  Orders 4 and 8; tolerance margins /10 and /1.  Scored against the
DOP853 references step 4b saved (b_fix.jsonl: ref13 at rtol 1e-13/atol 1e-15, ref12 at
1e-12/1e-14, spread = |ref12 - ref13|; verify_b.py holds the case list and the code that computed
them).  Run from the repo root:

  python tools/auto_tight/tight_measure.py run OUT.jsonl NPROC [--smoke]
  python tools/auto_tight/tight_measure.py score OUT.jsonl

H4 and H8 are 'auto' as it was before issue #120, the hybrid at that order: they close the route
the issue opened below 1e-6 (AUTO_LADDER_TIGHT_MAX_PHASE = -1), so a run on later code measures
the same comparison.  L* force the #70 route open at any tolerance with the margin given.

Writes only OUT.jsonl.  Timings are warm (one untimed call first), best of two, under NPROC
parallel workers: ratios inside one record are fair, absolute times are not; re-time borderline
cases alone on an idle machine before quoting them.
"""
import sys, os, json, time, warnings
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
TOLS = (1e-7, 1e-9, 1e-12)    # rtol; atol = rtol/100 as in step 4b, so 1e-12 is Listing 1's 1e-12/1e-14
CONFIGS = (('H4', 4, None), ('H8', 8, None),        # auto as on the branch: the hybrid at that order
           ('L4/10', 4, 10.0), ('L4/1', 4, 1.0),    # the route forced open, ladder margin /10 or /1
           ('L8/10', 8, 10.0), ('L8/1', 8, 1.0))
HONEST = {'HybridCertificationWarning', 'ToleranceNotAchievedWarning', 'UnmarkedDiscontinuityWarning'}


def refs():
    out = {}
    for line in open(os.path.join(HERE, 'b_fix.jsonl')):
        r = json.loads(line)
        if 'ref13' in r:
            out[(r['name'], r['E_gev'])] = (r['ref12'], r['ref13'], r.get('phase'))
    return out


def route_of(info):
    """The #70 route's trace entry if auto took it (only that entry carries estimated_phase)."""
    for v in info.values():
        for d in (v if isinstance(v, list) else ()):
            if isinstance(d, dict) and 'estimated_phase' in d:
                return dict(phase=d['estimated_phase'], n_floor=d.get('min_n_slabs'),
                            reason=d.get('reason'))
    return None


def scan_cases():
    import verify_b as vb
    groups = {}
    for name, E, spec in vb.cases():
        if name.startswith('L1-'):
            groups.setdefault(name, (spec, []))[1].append(E)
    return [(name + '-scan', Es, spec) for name, (spec, Es) in groups.items()]


def scan_f(spec):
    import magnus.oscprob as op, magnus.globaldefs as gd
    osc = gd.load_nufit_params('NuFIT 6.1')
    s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
    ex = dict(osc=osc, two=dict(sth=osc['s12'], Dm2=osc['D21']),
              s4=dict(osc, s14=s14, s24=s24, D41=1.0),
              s5=dict(osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7))[spec[4]]
    fn, L, h, rho = spec[:4]
    return lambda E, **k: getattr(op, fn)(np.asarray(E, dtype=float), L=L*gd.UNIT_KM, L0=0.0,
                                          rho_central=rho, l_scale=h*gd.UNIT_KM,
                                          density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE,
                                          nu_f=gd.NUE, **ex, **k)


def call(f, E, tol, order, margin, reps=2):
    import magnus.oscprob as op
    op.AUTO_LADDER_MIN_TOLERANCE = 1.0e-6 if margin is None else 0.0
    op.AUTO_LADDER_TOLERANCE_MARGIN = 10.0 if margin is None else margin
    op.AUTO_LADDER_TIGHT_MAX_PHASE = -1.0 if margin is None else 2000.0
    kw = dict(rtol=tol, atol=tol*1e-2, magnus_exp_order=order, strategy='auto')
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        t = time.perf_counter()
        try:
            P = np.asarray(f(E, strategy_info=info, **kw), dtype=float).ravel()
        except Exception as exc:
            return dict(error='%s: %s' % (type(exc).__name__, str(exc)[:150]))
        cold = time.perf_counter() - t
    best = cold
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(reps):
            if best > 3.0:
                break
            t = time.perf_counter(); f(E, **kw); best = min(best, time.perf_counter() - t)
    return dict(P=P.tolist(), t=best, t_cold=cold, engine=info.get('engine'),
                certified=info.get('certified'), warnings=sorted({x.category.__name__ for x in w}),
                route=route_of(info))


def run_case(args):
    kind, name, E_gev, spec = args
    warnings.simplefilter('ignore')
    import verify_b as vb, magnus.oscprob as op, magnus.globaldefs as gd
    if kind == 'point':
        f, _ = vb.call_for(spec); E = E_gev*gd.UNIT_GEV
    else:
        f = scan_f(spec); E = np.asarray(E_gev)*gd.UNIT_GEV
    rec = dict(kind=kind, name=name, E_gev=E_gev, rows=[])
    # Eligibility: the #70 route's conditions (phase, starting slab count) do not depend on the
    # tolerance, so "auto takes the route at 1e-3" is "eligible once the cut-off is lifted".
    op.AUTO_LADDER_MIN_TOLERANCE = 1.0e-6; op.AUTO_LADDER_TOLERANCE_MARGIN = 10.0
    info = {}
    try:
        f(E, rtol=1e-3, atol=1e-3, strategy='auto', strategy_info=info)
        rec['route_1e-3'] = route_of(info)
    except Exception as exc:
        rec['route_1e-3'] = None; rec['error_1e-3'] = str(exc)[:150]
    rec['eligible'] = rec['route_1e-3'] is not None
    for tol in TOLS:
        row = dict(tol=tol)
        for label, order, margin in CONFIGS:
            if margin is None or rec['eligible']:
                row[label] = call(f, E, tol, order, margin)
        rec['rows'].append(row)
    return rec


def score(path):
    R = refs()
    recs = [json.loads(line) for line in open(path)]
    print('%d records, %d eligible (the #70 route taken at 1e-3)'
          % (len(recs), sum(r['eligible'] for r in recs)))

    def err_lim(r, x, tol):
        if r['kind'] == 'point':
            keys = [(r['name'], r['E_gev'])]
        else:
            keys = [(r['name'][:-len('-scan')], E) for E in r['E_gev']]
        out = []
        for P, key in zip(x['P'], keys):
            if key in R:
                r12, r13, _ = R[key]
                out.append((abs(P - r13), max(tol*1e-2 + tol*abs(r13), abs(r12 - r13))))
        return max(out, key=lambda e: e[0]/e[1]) if out else None

    for tol in TOLS:
        print('\n=== rtol %g, atol %g   (limit = max(atol + rtol*P_ref, DOP853 spread))' % (tol, tol*1e-2))
        print('%-6s %5s %6s %6s %5s   t/t(H same order): median   max  n>1.2  n<1' %
              ('config', 'n', 'silent', 'warned', 'err'))
        for label, order, margin in CONFIGS:
            n = silent = warned = errors = 0; ratios = []; worst = []
            for r in recs:
                row = {x['tol']: x for x in r['rows']}.get(tol, {})
                x = row.get(label)
                if x is None:
                    continue
                if 'error' in x:
                    errors += 1; continue
                el = err_lim(r, x, tol)
                if el is None:
                    continue
                n += 1
                if HONEST & set(x['warnings']):
                    warned += 1
                elif el[0] > el[1]:
                    silent += 1; worst.append((el[0]/el[1], r['name'], r['E_gev']))
                base = row.get('H%d' % order)
                if margin is not None and base and 'error' not in base:
                    ratios.append((x['t']/base['t'], r['name'], r['E_gev'],
                                   (r['route_1e-3'] or {}).get('phase'),
                                   (r['route_1e-3'] or {}).get('n_floor')))
            rs = np.array([q[0] for q in ratios]) if ratios else np.array([np.nan])
            print('%-6s %5d %6d %6d %5d   %28.2f %6.2f %6d %4d' % (
                label, n, silent, warned, errors, np.nanmedian(rs), np.nanmax(rs),
                int(np.sum(rs > 1.2)), int(np.sum(rs < 1.0))))
            for q in sorted(worst, reverse=True)[:3]:
                print('        silent miss %.1fx outside: %s E=%s' % q)
            for q in sorted(ratios, key=lambda q: -q[0])[:4]:
                if q[0] > 1.2:
                    print('        slower %.2fx: %s E=%s  route phase %.3g  start slabs %s' % q)
    print('\n=== Listing 1 scans (26 reference energies per curve): time, worst err/limit, engine')
    for r in sorted((r for r in recs if r['kind'] == 'scan'), key=lambda r: r['name']):
        for row in r['rows']:
            cells = []
            for label, _, _ in CONFIGS:
                x = row.get(label)
                if x and 'error' not in x:
                    el = err_lim(r, x, row['tol'])
                    cells.append('%s %.3gs %.2g %s' % (label, x['t'], el[0]/el[1] if el else np.nan,
                                                       x['engine']))
            print('%-12s %g: %s' % (r['name'], row['tol'], ' | '.join(cells)))


if __name__ == '__main__':
    mode, out = sys.argv[1], sys.argv[2]
    if mode == 'score':
        score(out)
        sys.exit()
    import subprocess
    from multiprocessing import Pool
    import verify_b as vb, magnus.oscprob as op
    print('magnus from', op.__file__, subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
          capture_output=True, text=True).stdout.strip(), flush=True)
    jobs = ([('scan', n, Es, s) for n, Es, s in scan_cases()]
            + [('point', n, E, s) for n, E, s in vb.cases()])
    if '--smoke' in sys.argv:
        TOLS = (1e-12,)
        jobs = [jobs[1]] + [j for j in jobs if j[1] in ('L1-3nu', 'sun2-exp-0.1R')][:1] \
            + [j for j in jobs if j[1] == 'sun2-exp-0.1R'][:1]
    with Pool(int(sys.argv[3])) as pool, open(out, 'w') as fh:
        for rec in pool.imap_unordered(run_case, jobs):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print('done', rec['kind'], rec['name'], rec['E_gev'] if rec['kind'] == 'point' else '',
                  'eligible' if rec['eligible'] else '', flush=True)
    print('all done', flush=True)
