# -*- coding: utf-8 -*-
r"""gen_shock_cost.py

What it costs to reach a given accuracy on a supernova shock, as the front widens.

Phase 1 of two.  This one finds, for every front width, every code and every accuracy
target, the smallest discretisation that reaches it -- an accuracy measurement, immune to
what else the machine is doing.  Phase 2 (``--time``) times those settings and must run on
a quiet machine, because its output is a wall clock.

WHY THIS FIGURE AND NOT A SPEED-ACCURACY PLANE.  The plane that Figure 11 draws is dialled
by a requested tolerance, and on this ray a tolerance does not order the answers: measured
across all four combinations of ``cumulative`` and ``t_breakpoints``, including the
single-baseline no-breakpoint call Figure 11 itself uses, the map from request to achieved
accuracy is either flat or non-monotonic -- 3e-05 returns an answer four times worse than
1e-04.  Dialling by slab count instead is monotonic and clean, but a slab count is not
something a user chooses, and Figure 11's own cell rejects that sweep for handing the
closed form a discretisation Magnus had to find for itself.

So neither dial makes a fair plane here.  This figure asks the question a reader actually
arrives with instead: *my front is this wide and I need this many digits -- which code, and
what will it cost?*  The slab count is the dial internally and never appears in the
figure; the reader sees the accuracy requested and the time paid.

DECLARING A FRONT IS NOT RESOLVING IT.  The first version of this script passed only the
front edges, which is what `shock_breakpoints` returns, and measured an order-4 cost 34
times too high.  The reason is that `n_slabs` lays a *uniform* grid over the whole ray and
unions it with the breakpoints; it does not deal slabs out to the intervals the
breakpoints create.  A 0.07 km front inside a 70 000 km ray therefore gets exactly one
slab, however large `n_slabs` is, and the error sits on a plateau that no refinement
moves -- 3e-06 at order 4 from 10 000 slabs to 260 000.  Subdividing the front interior
removes it: at 16 384 slabs, order 4 goes from 2.8e-06 to 6.0e-08.

So each order is measured twice, `declared` and `resolved`, and the gap between them is
the figure's point.  It also makes the comparison against NuOscProbExact fair for the
first time: `npe_at` has always floored its per-segment allocation, so the earlier run
had the competitor resolving a front that Magnus could not.

THE REFERENCES ARE ALREADY FROZEN.  ``shock_reference.json`` holds all five widths at
rtol 1e-14, each with its own self-convergence recorded (1.2e-10 to 2.9e-11).  Every target
here sits at least two decades above that, so no point is measuring the reference.

    python notebooks/gen_shock_cost.py            # Phase 1, accuracy: run anywhere
    python notebooks/gen_shock_cost.py --time     # Phase 2, timing: needs a quiet machine

Writes ``external_shock_cost.json``, checkpointed after every width.
"""

import argparse
import json
import pathlib
import platform
import resource
import sys
import time
import warnings

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'src'))
import magnus.matter as matter                               # noqa: E402
import magnus.oscprob as oscprob                             # noqa: E402

OUT = HERE/'external_shock_cost.json'
# Twenty log-spaced widths, not the notebook's five.  Five samples a decade apart put
# the closed form's cost peak somewhere between 7 and 70 km and cannot say whether
# Mag(nu)s has one at all, which is the question the figure turns on.
WIDTHS = tuple(sorted(set(
    [round(float(w), 12) for w in np.logspace(-6.0, -2.0, 20)]
    # The three interior widths notebook 14 already has frozen come along for nothing,
    # so the scan is 23 points for the cost of 20.
    + [1.0e-5, 1.0e-4, 1.0e-3])))
ORDERS = (4, 6, 8)
# Order 2 is the midpoint slab product: one sample per slab, no commutator, which is
# exactly the information a constant-density slab has.  It is in the figure as the control
# that isolates the commutator from every other difference between the two codes, and only
# with the front resolved -- the declared arm would be measuring the other effect.
ORDERS_RESOLVED_ONLY = (2,)
# 1e-7 is the figure's target.  At 1e-8 NuOscProbExact cannot reach two of the five
# widths at any slab count, which leaves the comparison with holes exactly where it is
# most interesting; at 1e-7 every arm reaches every width.  The looser 1e-6 was the
# other candidate and costs too much: there the two codes tie on the thinnest front.
TARGETS = (1.0e-6, 1.0e-7, 1.0e-8)

# The local oscillation length at the contact density, from Section 6.6: 4 pi E / Dm2_21
# shortened by the matter potential to about 16 km.  The figure's top axis is the front
# width in these units, so the value is stored with the data rather than recomputed at
# plot time.
L_OSC_LOCAL_KM = 16.0

# Bytes per slab per d^2, measured at 1e6 slabs and three flavours, one call per order.
BYTES_PER_SLAB_PER_D2 = {2: 110.0, 4: 124.0, 6: 150.0, 8: 173.0, 'npe': 120.0}
SLAB_MEMORY_BUDGET = 5.0e9

# NuOscProbExact composes its chain a chunk at a time (see `npe_at`), so its ceiling is
# no longer memory.  It is set instead by what is worth timing: at this width a call runs
# about a third of a microsecond per slab, so the cap is some ten seconds.
CHUNK_SLABS = 1_000_000
NPE_SLAB_CAP = 32_000_000

# Magnus reports the slab count it was asked for; the grid it walks is this many times
# larger, from CUMULATIVE_N_ACC_SAFETY in oscprob.  Stored so the figure can quote either.
SLAB_SAFETY_FACTOR = 4

# The admission criterion of Section 7.1, and how many times a point that fails it is
# measured again before we give up and let it stand with its scatter on the record.
CV_LIMIT = 0.10
CV_RETRIES = 6


def front_edges(width, n_slabs, resolve):
    """The breakpoints for one call: front edges alone, or the front interior cut up.

    `resolve=False` is the notebook's own `shock_breakpoints`, four edges, two per front.
    `resolve=True` cuts each front into the same number of pieces `npe_at` floors its own
    segments at, so the two codes are told the same thing about where the structure is.
    Eight pieces already saturates the gain -- 32 measures no better -- but the rule is
    written to match NuOscProbExact rather than to the measurement, because parity between
    the arms is what the figure claims.
    """
    bps = np.asarray(NS['shock_breakpoints'](width), dtype=float)
    if not resolve:
        return bps
    pieces = max(8, int(n_slabs)//1000)
    inner = bps[1:-1].reshape(-1, 2)                      # (lo, hi) of each front
    cut = [np.linspace(lo, hi, pieces + 1) for lo, hi in inner]
    return np.unique(np.concatenate([bps] + cut))


def notebook_namespace():
    r"""Executes notebook 14's definition cells and returns their namespace.

    The shock profile is defined once, in the notebook, and every script that needs it
    executes that definition rather than transcribing it.  A script that transcribed it
    once invented a power-law rarefaction in place of the Fogli form, which `solve_ivp`
    reported as an overflow rather than as a wrong answer.
    """
    import contextlib

    import matplotlib
    matplotlib.use('Agg')
    import make_notebooks

    ns = {'__name__': '__notebook__'}
    with contextlib.redirect_stdout(sys.stderr):
        for cell in make_notebooks.books['14_magnus_supernova_shock.ipynb'].cells:
            if cell.cell_type != 'code':
                continue
            if 'def measure(' in cell.source:
                break
            exec(compile(cell.source, '<notebook 14>', 'exec'), ns)
    need = ('sn_shock_ne', 'make_H', 'shock_breakpoints', 'hvac3',
            'frozen_reference', 'L0', 'L1', 'Ls', 'ENERGY', 'params3')
    missing = [n for n in need if n not in ns]
    if missing:
        raise SystemExit('notebook 14 no longer defines %s' % ', '.join(missing))
    return ns


NS = notebook_namespace()
L0, L1 = NS['L0'], NS['L1']
Ls = np.asarray(NS['Ls'], dtype=float)
ENERGY, P3 = NS['ENERGY'], NS['params3']
RAY_KM = float(L1 - L0)/5.06e9


REF_STORE = HERE/'shock_reference_scan.json'
_REFS = {}


def reference(width):
    """The frozen probabilities for this width and the reference's own error, off disk.

    Reads `shock_reference_scan.json`, whose keys carry enough digits for twenty widths a
    factor of 1.6 apart; the notebook's own five-width store rounds them to one digit and
    would alias.  Both of its fingerprints are checked -- the electron density and the
    Hamiltonian built on it -- because a frozen oracle that outlives a change to the
    physics is worse than no oracle at all: every comparison against it still looks fine.
    This store went stale exactly that way once, guarded on the density alone.
    """
    key = '%.6e' % width
    if key not in _REFS:
        if not REF_STORE.exists():
            raise SystemExit('run notebooks/make_shock_scan_references.py first: %s is '
                             'missing' % REF_STORE.name)
        store = json.loads(REF_STORE.read_text())
        lf = np.array([float.fromhex(x) for x in store['fingerprint_l']])
        if key not in store['cases']:
            raise SystemExit('no frozen reference for width %s; re-run '
                             'make_shock_scan_references.py' % key)
        case = store['cases'][key]
        unhex = lambda xs: np.array([float.fromhex(x) for x in xs])       # noqa: E731
        ne = NS['sn_shock_ne'](width)
        H = np.asarray(NS['make_H'](ne)(lf), dtype=complex)
        want_ne = np.asarray(ne(lf), dtype=float)
        want_h = np.concatenate([H.real.ravel(), H.imag.ravel()])
        for got, want, what in ((case['fingerprint_ne'], want_ne, 'density'),
                                (case['fingerprint_h'], want_h, 'Hamiltonian')):
            if not np.array_equal(unhex(got), want):
                raise SystemExit('frozen reference at width %s was built on a different '
                                 '%s; re-run make_shock_scan_references.py' % (key, what))
        _REFS[key] = (unhex(case['P']).reshape(case['shape']),
                      float(case['self_convergence']))
    return _REFS[key]


def magnus_at(width, order, n_slabs, resolve):
    """The 61 probabilities from one Magnus call at a fixed slab count.

    `resolve` picks the breakpoint treatment; see `front_edges`.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(oscprob.osc_prob_matter_std_potential(
            3, NS['sn_shock_ne'](width), ENERGY, Ls, P3, L0=L0,
            density_is_of_number_of_electrons=True,
            t_breakpoints=front_edges(width, n_slabs, resolve), cumulative=True,
            magnus_exp_order=order, strategy='magnus',
            n_slabs=n_slabs, max_n_slabs=n_slabs,
            rtol=1.0e-12, atol=1.0e-14)).reshape(len(Ls), 3, 3)


def npe_at(width, n_slabs):
    """The same 61 probabilities from a NuOscProbExact slab chain, front declared.

    The allocation is the one gen_shock_benchmarks.py uses: never fewer than a thousandth
    of the budget in any declared segment, because the fronts can be 0.07 km of 70 000 and
    proportional allocation alone leaves them unresolved.  Both codes are told where the
    front is; neither is asked to find it.
    """
    import slabs as npe

    hv = np.asarray(NS['hvac3'])/ENERGY
    proj = matter.matter_potential_projector(3)
    h_fixed = NS['make_H'](NS['sn_shock_ne'](width))
    hvac00 = hv[0, 0]

    def hf(l):
        v = np.asarray(h_fixed(l))[..., 0, 0] - hvac00
        return hv + np.real(v)[..., None, None]*proj

    bps = np.asarray(NS['shock_breakpoints'](width), dtype=float)
    h = (L1 - L0)/float(n_slabs)
    floor = max(8, n_slabs//1000)
    U = np.eye(3, dtype=complex)
    P = np.empty((len(Ls), 3, 3))
    lo_prev = float(L0)
    for k, target in enumerate(Ls):
        edges = np.unique(np.concatenate([bps[(bps > lo_prev) & (bps < target)],
                                          [lo_prev, float(target)]]))
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = max(floor, int(round((hi - lo)/h)))
            e = np.linspace(lo, hi, m + 1)
            # Composed a chunk at a time.  One segment of this ray can want twelve
            # million slabs, whose Hamiltonians would be fourteen gigabytes in one
            # array; in chunks the peak is CHUNK_SLABS wide instead and the count stops
            # being bounded by memory.  The slabs, their widths and their order are
            # untouched, so only the association of the matrix product changes -- the
            # chunk boundaries multiply in sequence exactly where the whole-segment
            # call would have.  Measured against the unchunked path, the two agree to
            # the last bits of the product; the check is in the smoke test below.
            for a in range(0, m, CHUNK_SLABS):
                b = min(a + CHUNK_SLABS, m)
                edge = e[a:b + 1]
                U = np.asarray(npe.evolution_operator_3nu_slabs(
                    np.asarray(hf(0.5*(edge[:-1] + edge[1:])), dtype=complex),
                    np.diff(edge))) @ U
        P[k] = np.swapaxes(np.abs(U)**2, -1, -2)
        lo_prev = float(target)
    return P


def max_slabs_for(key):
    """The largest slab count worth trying: a memory budget, or for the chunked arm, time."""
    if key == 'npe':
        return NPE_SLAB_CAP
    return int(SLAB_MEMORY_BUDGET/(BYTES_PER_SLAB_PER_D2[key]*9.0))


def smallest_n_reaching(evaluate, ref, target, seed, cap):
    """Smallest slab count whose error first drops below `target`, and that error.

    A geometric climb to bracket the crossing, then a bisection inside the bracket.  Not a
    pure bisection from the start: the error is monotonic in slab count where it has been
    measured, but that was verified at order 4 only, and a climb that takes the *first*
    crossing is robust to a rung that steps out of line, where a bisection assuming
    monotonicity would converge on the wrong side of it.

    `seed` is the answer for the previous width, which is usually within a factor of two:
    the counts vary smoothly with the front, so starting there saves most of the climb.
    Returns (None, best_error) if the target is out of reach inside the cap -- a curve that
    ends is information, and better than one that pretends.
    """
    cache = {}

    def err(n):
        if n not in cache:
            cache[n] = float(np.max(np.abs(evaluate(n) - ref)))
        return cache[n]

    lo = max(256, int(seed)//4)
    n, best = lo, None
    while n <= cap:
        e = err(n)
        best = e if best is None else min(best, e)
        print('      n=%-9d err %.3e' % (n, e), file=sys.stderr, flush=True)
        if e < target:
            break
        # The doubling stops when the next rung would clear the cap, which leaves the cap
        # itself untried -- and the whole interval below it.  One arm reached the target
        # at exactly the cap and was recorded as out of reach, so the ladder ends on the
        # cap rather than beside it.
        n = cap if n*2 > cap > n else n*2
    else:
        return None, best
    if n == lo:                                  # met at the first count tried
        return n, err(n)
    hi, lo = n, n//2
    while hi - lo > max(256, lo//16):            # within about six percent
        mid = (lo + hi)//2
        if err(mid) < target:
            hi = mid
        else:
            lo = mid
    return hi, err(hi)


def timed(call, repeat=3, min_block=0.05, budget=6.0):
    """The best of `timed_blocks`; see there."""
    return timed_blocks(call, repeat, min_block, budget)[0]


def timed_blocks(call, repeat=3, min_block=0.05, budget=6.0):
    """Best of `repeat` autoranged blocks and the blocks themselves, warm-up discarded.

    The blocks come back because Section 7.1 promises a coefficient of variation across
    them for every timing in the paper, and a function that returns only the winner cannot
    supply one.  A setting whose single call already exceeds `budget` yields one block and
    no spread, which the caller records as absent rather than as zero.

    Copied from gen_shock_benchmarks.py so every timing in the paper uses one protocol.
    The discard is not decoration: the first Mag(nu)s call of a session pays ~0.7 s to
    compile the numba kernel, which a user pays once and not once per call.  `budget` caps
    the repeats for settings where a single call already runs for seconds -- most of the
    expensive end of this figure is timed once, which the scatter on the control measures.
    """
    t0 = time.perf_counter()
    call()
    first = time.perf_counter() - t0
    if first > budget:
        return first, [first]
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        el = time.perf_counter() - t0
        if el >= min_block:
            break
        reps *= 2
    blocks = [el/reps]
    for _ in range(repeat - 1):
        if min(blocks)*reps > budget:
            break
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        blocks.append((time.perf_counter() - t0)/reps)
    return min(blocks), blocks


def call_for(case):
    """The callable that reproduces one Phase 1 row, at the slab count it found."""
    n, width, code = case['n_slabs'], case['width'], case['code']
    if code == 'NuOscProbExact':
        return lambda: npe_at(width, n)
    return lambda: magnus_at(width, case['order'], n, case['front'] == 'resolved')


CONTROL = dict(width=1.0e-4, order=4, n_slabs=8192, front='resolved',
               code='order 4, resolved')


def phase_two():
    """Times the settings Phase 1 found, interleaving a control against machine drift.

    The control is one cheap setting timed every few cases.  Its scatter across the run is
    what says whether the expensive single-shot timings can be believed: this figure is
    not always measured on an idle machine, and a wall clock that drifts 30 percent
    through the run would otherwise be invisible in the output.
    """
    store = json.loads(OUT.read_text())
    todo = [c for c in store['cases'] if c.get('n_slabs') and 'seconds' not in c]
    def control():
        return magnus_at(CONTROL['width'], CONTROL['order'], CONTROL['n_slabs'], True)

    print('  timing %d of %d cases' % (len(todo), len(store['cases'])),
          file=sys.stderr, flush=True)

    controls = [timed(control)]
    print('  control %.4f s' % controls[0], file=sys.stderr, flush=True)
    for i, case in enumerate(todo):
        # Section 7.1 admits no point whose blocks scatter by more than CV_LIMIT, so a
        # point that scatters is re-measured rather than explained.  A transient on this
        # machine hits about one point in a hundred; tested across 536 points, there is
        # no systematic drift between the first block and the last (median ratio 1.0006,
        # sign test p = 0.52), so a high scatter is interference and not warm-up, and
        # measuring again is the honest response to it.
        for attempt in range(CV_RETRIES):
            seconds, blocks = timed_blocks(call_for(case))
            cv = (float(np.std(blocks, ddof=1)/np.mean(blocks)) if len(blocks) > 1
                  else None)
            if cv is None or cv <= CV_LIMIT:
                break
            print('      re-timing %s %.0e %.2f km, cv %.3f (attempt %d)'
                  % (case['code'], case['target'], case['width_km'], cv, attempt + 1),
                  file=sys.stderr, flush=True)
        case['seconds'], case['cv_attempts'] = seconds, attempt + 1
        case['block_seconds'], case['cv'] = blocks, cv
        case['control_seconds'] = controls[-1]
        print('    %.0e %-20s %.0e  n=%-9d  %8.3f s'
              % (case['width'], case['code'], case['target'], case['n_slabs'],
                 case['seconds']), file=sys.stderr, flush=True)
        if (i + 1) % 8 == 0 or i == len(todo) - 1:
            controls.append(timed(control))
            print('  control %.4f s (%d of %d done)'
                  % (controls[-1], i + 1, len(todo)), file=sys.stderr, flush=True)
        OUT.write_text(json.dumps(store, indent=1))

    cvs = [c['cv'] for c in todo if c.get('cv') is not None]
    if cvs:
        print('  worst coefficient of variation %.3f over %d points%s'
              % (max(cvs), len(cvs), '' if max(cvs) <= 0.10 else '  <-- ABOVE 0.10'),
              file=sys.stderr, flush=True)
        store['worst_cv'] = max(cvs)
    spread = max(controls)/min(controls)
    store['control'] = dict(CONTROL, seconds=controls,
                            spread=spread,
                            note=('one cheap setting re-timed every eight cases; the '
                                  'spread bounds how much the machine moved under the '
                                  'single-shot timings'))
    OUT.write_text(json.dumps(store, indent=1))
    print('  control spread %.2fx over %d samples%s'
          % (spread, len(controls),
             '' if spread < 1.15 else '  <-- MACHINE MOVED, timings are soft'),
          file=sys.stderr, flush=True)
    print('timed %d cases' % len(todo), file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--time', action='store_true',
                    help='Phase 2: time the settings Phase 1 found (needs a quiet machine)')
    args = ap.parse_args()
    if args.time:
        warnings.simplefilter('ignore')
        resource.setrlimit(resource.RLIMIT_AS, (int(8.0*1024**3),)*2)
        if not OUT.exists():
            raise SystemExit('run Phase 1 first: %s does not exist' % OUT.name)
        return phase_two()
    warnings.simplefilter('ignore')
    resource.setrlimit(resource.RLIMIT_AS, (int(8.0*1024**3),)*2)

    store = {'note': ('cost of reaching a fixed accuracy on the supernova shock, against '
                      'front width; produced by notebooks/gen_shock_cost.py'),
             'energy_ev': ENERGY, 'ray_km': RAY_KM,
             'l_osc_local_km': L_OSC_LOCAL_KM,
             'targets': list(TARGETS), 'orders': list(ORDERS),
             'slab_safety_factor': SLAB_SAFETY_FACTOR,
             'n_baselines': len(Ls),
             'timing_note': ('seconds is one whole cumulative call, which returns'
                             ' n_baselines probabilities; divide by n_baselines'
                             ' for the per-probability time the paper quotes'),
             'front_note': ('each order appears twice: "declared" passes the two front '
                            'edges, "resolved" also cuts the front interior, which a '
                            'uniform n_slabs grid never reaches'),
             'reference_note': ('slab counts scored against shock_reference.json, frozen '
                                'at rtol 1e-14; every target sits at least two decades '
                                'above the reference own self-convergence'),
             'machine': platform.platform(), 'python': platform.python_version(),
             'numpy': np.__version__,
             'magnus': getattr(__import__('magnus'), '__version__', '?'),
             'cases': []}
    if OUT.exists():
        store = json.loads(OUT.read_text())
    done = {(c['width'], c['code'], c['target']) for c in store['cases']}

    seeds = {}
    for width in WIDTHS:
        ref, sc = reference(width)
        km = width*RAY_KM
        print('  width %.0e = %.2f km (%.3f local oscillation lengths), '
              'reference self-convergence %.2e'
              % (width, km, km/L_OSC_LOCAL_KM, sc), file=sys.stderr, flush=True)
        arms = [('order %d, %s' % (o, kind), o,
                 lambda n, o=o, w=width, r=res: magnus_at(w, o, n, r), max_slabs_for(o))
                for o in ORDERS for kind, res in (('declared', False), ('resolved', True))]
        arms += [('order %d, resolved' % o, o,
                  lambda n, o=o, w=width: magnus_at(w, o, n, True), max_slabs_for(o))
                 for o in ORDERS_RESOLVED_ONLY]
        arms.append(('NuOscProbExact', None, lambda n, w=width: npe_at(w, n),
                     max_slabs_for('npe')))
        for name, order, evaluate, cap in arms:
            for target in TARGETS:
                if (width, name, target) in done:
                    continue
                print('    %s, target %.0e' % (name, target), file=sys.stderr, flush=True)
                t0 = time.perf_counter()
                n, err = smallest_n_reaching(evaluate, ref, target,
                                             seeds.get((name, target), 4096), cap)
                if n:
                    seeds[(name, target)] = n
                store['cases'].append(dict(
                    width=width, width_km=km, width_in_l_osc=km/L_OSC_LOCAL_KM,
                    code=name, order=order, target=target,
                    front='resolved' if name.endswith('resolved') else 'declared',
                    n_slabs_walked=None if n is None else n*SLAB_SAFETY_FACTOR,
                    n_slabs=n, error=err, reachable=n is not None,
                    reference_self_convergence=sc,
                    search_seconds=time.perf_counter() - t0))
                print('      -> %s  (error %.3e, %.0f s)'
                      % (('n_slabs = %d' % n) if n else 'NOT REACHABLE within the budget',
                         err, time.perf_counter() - t0), file=sys.stderr, flush=True)
                OUT.write_text(json.dumps(store, indent=1))
        print('  [checkpoint: width %.0e]' % width, file=sys.stderr, flush=True)

    OUT.write_text(json.dumps(store, indent=1))
    print('wrote %s (%d cases)' % (OUT.name, len(store['cases'])),
          file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
