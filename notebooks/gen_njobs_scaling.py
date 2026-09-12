# -*- coding: utf-8 -*-
r"""What ``n_jobs`` costs and buys, against the number of workers.

``n_jobs`` reads like a performance knob and is not one.  On the common request -- many
energies sharing one baseline -- the batched engine answers at ``n_jobs = 1`` and
**declines** at anything else, at ``src/magnus/oscprob.py:4726``, so raising the worker
count does not add parallelism to the batched path: it trades that path away for process
parallelism over single points.  The docstring has said the default "is usually fastest"
for some time with no measurement behind it.  This measures it.

Two arms, both on an Earth chord at ``cos(theta_z) = -0.9``, which is the scan the
question actually arises for.

* ``shared baseline``  -- what a caller gets by default.  Every energy shares one
  baseline, so the batched engine is available at ``n_jobs = 1`` and is given up above it.
* ``own baselines``    -- the advanced case.  Each energy carries its own baseline, which
  the batched engine never accepts, so ``n_jobs`` here is parallelism added to a path that
  was serial anyway, with nothing traded for it.

Read together they separate the two effects: the second arm is what process parallelism
is worth, and the first is what it costs to reach it.

WHAT THIS FILE HAD TO GET RIGHT, each of which was wrong at least once while writing it.

* **Threads are pinned to one per process, before NumPy is imported.**  Unpinned, MKL
  takes ten threads and OpenMP twelve on this machine, so ``n_jobs = 1`` is not one core
  and ``n_jobs = 12`` is twelve processes each entitled to ten threads -- a hundred and
  twenty threads on twelve cores, measuring oversubscription rather than scaling.
* **The entry point matters.**  ``osc_prob_energy_baseline`` with a hand-built Hamiltonian
  never reaches the batched engine at all; only the standard-potential wrappers do.
* **The profile matters.**  A genuine exponential profile is claimed by the
  interaction-picture engine, and a smooth polynomial ramp by the hybrid engine, both of
  which run *before* the batched path and ignore ``n_jobs`` entirely.  On those profiles
  the honest answer is that the knob does nothing, which is worth knowing and is not what
  this figure is about.
* **joblib reuses its worker pool.**  The first call at a given worker count pays start-up
  and later ones do not, so each configuration is warmed before it is timed; otherwise the
  first configuration measured carries the pool for all the others.
* **Which engine answered is recorded, not assumed.**  The dispatchers are instrumented
  and the answer stored beside each point, so the figure can say what ran rather than what
  was expected to run.

    python notebooks/gen_njobs_scaling.py --smoke   # a few points, minutes
    python notebooks/gen_njobs_scaling.py           # the figure, needs a quiet machine

Writes ``external_njobs_scaling.json``, checkpointed after every worker count.
"""

import os

# Before NumPy, and before anything that imports it.  See the note above.
for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_var] = '1'

import argparse                                              # noqa: E402
import random                                                # noqa: E402
import json                                                  # noqa: E402
import pathlib                                               # noqa: E402
import platform                                              # noqa: E402
import sys                                                   # noqa: E402
import time                                                  # noqa: E402
import warnings                                              # noqa: E402

import numpy as np                                           # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'src'))

import joblib                                                # noqa: E402
import magnus.earth as earth                                 # noqa: E402
import magnus.globaldefs as gd                               # noqa: E402
import magnus.oscprob as oscprob                             # noqa: E402

OUT = HERE/'external_njobs_scaling.json'
N_JOBS = tuple(range(1, 13))                                 # one worker to one per core
COSTHZ = -0.9
ENERGY_RANGE_GEV = (0.5, 20.0)
RTOL, ATOL = 1.0e-3, 1.0e-3                                  # the defaults a caller gets

# This machine: cpu0/1 and cpu10/11 are the two performance cores, hyperthreaded; cpu2
# through cpu9 are eight efficiency cores with one thread each.  Ten physical cores, twelve
# logical.  Pinning to the eight identical ones is a DIAGNOSTIC -- no caller restricts
# magnus to a cpu mask -- and appears only in the appendix figure, where it shows how much
# of the small-scan wobble the fast/slow mixture accounts for.
E_CORES = frozenset(range(2, 10))
N_PHYSICAL_CORES = 10

# Passes over the whole sweep.  More for the short scans, which need them: a 500-point
# scan runs in 0.3 s and a 1000-point one in 0.6 s.
ROUNDS_BY_SIZE = {500: 30, 1000: 30, 5000: 14, 20000: 8}
DEFAULT_ROUNDS = 8

# Every arm the two figures need.  `order` is the one that matters and was wrong until
# now: a sweep that runs worker counts 1, 2, ... 12 in that order every round is not
# interleaved, because anything that ramps within a round gives count 1 the cold slot and
# count 12 the hot one, identically each time.  Averaging cannot remove that; randomising
# the order can.  Measured: it moved the 500-point roughness by a quarter.
ARMS = (
    # label,                          points, baselines, order,        pinned
    #
    # The three large arms are NOT here: they are already measured and stored, and were
    # carried into this file's schema rather than re-run.  They carry order='fixed',
    # which is honest about how they were taken.
    #
    # The paper figure needs this one, which the appendix shares:
    ('1000 points',                     1000, 'own',     'randomised', False),
    # and the appendix needs the same two scans measured three ways:
    ('500 points, fixed order',          500, 'own',     'fixed',      False),
    ('500 points, randomised',           500, 'own',     'randomised', False),
    ('500 points, pinned',               500, 'own',     'randomised', True),
    ('1000 points, fixed order',        1000, 'own',     'fixed',      False),
    ('1000 points, pinned',             1000, 'own',     'randomised', True),
)


OSC = gd.load_nufit_params('NuFIT 6.1')
PARAMS = {k: OSC[k] for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
# Every arm builds its own energy grid; see `main`.
RANDOM_SEED = 20260912       # so a re-run shuffles the same way
L_CHORD = earth.distance_traveled_inside_earth(COSTHZ)*gd.CONV_KM_TO_INV_EV

# Which engine answered, observed rather than assumed.  Each dispatcher returns
# NotImplemented when the request does not fit it, so the one that returns anything else
# is the one that ran.
DISPATCHERS = ('_osc_prob_hybrid_dispatch', '_osc_prob_ip_exp_dispatch',
               '_osc_prob_scan_separable_dispatch')
OBSERVED = {}


def _watch(name):
    real = getattr(oscprob, name)

    def wrapper(*args, **kwargs):
        out = real(*args, **kwargs)
        if out is not NotImplemented:
            OBSERVED[name] = OBSERVED.get(name, 0) + 1
        return out
    return wrapper


for _name in DISPATCHERS:
    setattr(oscprob, _name, _watch(_name))

_REAL_PARALLEL = joblib.Parallel


class _WatchedParallel(_REAL_PARALLEL):
    """Counts fan-outs, so 'n_jobs did nothing' is a measurement and not an inference."""

    def __call__(self, *args, **kwargs):
        OBSERVED['fanout'] = OBSERVED.get('fanout', 0) + 1
        return super().__call__(*args, **kwargs)


oscprob.Parallel = _WatchedParallel
joblib.Parallel = _WatchedParallel


def engine_and_fanouts():
    """What ran during the last call, as a short label and a fan-out count."""
    answered = [n.replace('_osc_prob_', '').replace('_dispatch', '')
                for n in DISPATCHERS if OBSERVED.get(n)]
    return (answered[0] if answered else 'per-point'), OBSERVED.get('fanout', 0)


def call(energies, baselines, n_jobs):
    """One whole scan: every energy, at the default strategy a caller would get."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return oscprob.osc_prob_3nu_earth(energies, costhz=COSTHZ, L=baselines,
                                          **PARAMS, rtol=RTOL, atol=ATOL, n_jobs=n_jobs)


def timed_blocks(run, repeat=3, min_block=0.05, budget=6.0):
    """Best of `repeat` autoranged blocks and the blocks themselves, warm-up discarded.

    The same protocol as the rest of the paper (Section 7.1); see
    notebooks/gen_shock_cost.py, from which this is copied so that the two figures are
    measured the same way.
    """
    t0 = time.perf_counter()
    run()
    first = time.perf_counter() - t0
    if first > budget:
        return first, [first]
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            run()
        elapsed = time.perf_counter() - t0
        if elapsed >= min_block:
            break
        reps *= 2
    blocks = [elapsed/reps]
    for _ in range(repeat - 1):
        if min(blocks)*reps > budget:
            break
        t0 = time.perf_counter()
        for _ in range(reps):
            run()
        blocks.append((time.perf_counter() - t0)/reps)
    return min(blocks), blocks


def one_block(energies, baselines, n_jobs):
    """One timed pass at one worker count, on a warm pool, and what answered it.

    The discarded call is not optional.  joblib keeps ONE reusable executor and rebuilds it
    when the worker count changes, so a sweep that changes the count every call -- which
    randomised order does by design -- would otherwise pay pool creation on every timing.
    Measured on a short scan that was 1.1 s against 0.03 s of work.
    """
    call(energies, baselines, n_jobs)
    OBSERVED.clear()
    t0 = time.perf_counter()
    call(energies, baselines, n_jobs)
    seconds = time.perf_counter() - t0
    engine, fanouts = engine_and_fanouts()
    return seconds, engine, fanouts


def sweep(energies, baselines, n_jobs_values, rounds, order):
    """Every worker count, `rounds` times, in fixed or randomised order.

    `order='fixed'` reproduces the flawed sweep, kept because the appendix figure exists to
    show what it costs.  `order='randomised'` shuffles within every round, which is what
    decorrelates a ramp from the worker count.
    """
    rng = random.Random(RANDOM_SEED)
    passes = {n: [] for n in n_jobs_values}
    engines = {}
    for index in range(rounds):
        sequence = list(n_jobs_values)
        if order == 'randomised':
            rng.shuffle(sequence)
        for n_jobs in sequence:
            seconds, engine, fanouts = one_block(energies, baselines, n_jobs)
            passes[n_jobs].append(seconds)
            engines[n_jobs] = (engine, fanouts)
        if (index + 1) % 5 == 0 or index == rounds - 1:
            print('      round %d of %d' % (index + 1, rounds), file=sys.stderr,
                  flush=True)

    rows = []
    for n_jobs in n_jobs_values:
        times = passes[n_jobs]
        engine, fanouts = engines[n_jobs]
        rows.append(dict(
            n_jobs=n_jobs, seconds=min(times), mean_seconds=float(np.mean(times)),
            rounds=len(times), pass_seconds=times, engine=engine, fanouts=fanouts,
            # Scatter across rounds separated in time, not across blocks taken back to
            # back: the first says the number is reproducible, the second only that the
            # clock held still for a few seconds.
            cv=float(np.std(times, ddof=1)/np.mean(times)) if len(times) > 1 else None))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true',
                        help='two arms and three worker counts, to check the machinery')
    args = parser.parse_args()

    arms, n_jobs_values = ARMS, N_JOBS
    if args.smoke:
        arms = (ARMS[0], ARMS[-1])
        n_jobs_values = (1, 2, 8)

    store = {'note': ('mean speed-up against n_jobs on an Earth chord, at several scan '
                      'sizes and three measurement protocols; produced by '
                      'notebooks/gen_njobs_scaling.py'),
             'costhz': COSTHZ, 'energy_range_gev': list(ENERGY_RANGE_GEV),
             'rtol': RTOL, 'atol': ATOL, 'threads_per_process': 1,
             'n_physical_cores': N_PHYSICAL_CORES,
             'e_cores': sorted(E_CORES),
             'protocol': ('each timed call runs on a pool warmed inside its own round; '
                          'arms marked randomised shuffle the worker order every round, '
                          'which arms marked fixed do not; pinned arms run on the eight '
                          'identical efficiency cores and stop at eight workers; '
                          '"seconds" is the best pass, "mean_seconds" the mean over the '
                          'rounds, and "speedup" the ratio of MEANS, which is what the '
                          'figure plots'),
             'machine': platform.platform(), 'python': platform.python_version(),
             'numpy': np.__version__,
             'magnus': getattr(__import__('magnus'), '__version__', '?'),
             'arms': {}}
    if OUT.exists() and not args.smoke:
        store = json.loads(OUT.read_text())

    print('  chord at cos(theta_z) = %.2f, one thread per process, rtol %.0e'
          % (COSTHZ, RTOL), file=sys.stderr, flush=True)

    for label, n_points, baseline_kind, order, pinned in arms:
        if label in store['arms']:
            continue
        rounds = ROUNDS_BY_SIZE.get(n_points, DEFAULT_ROUNDS)
        counts = tuple(n for n in n_jobs_values if not pinned or n <= len(E_CORES))
        energies = np.logspace(*np.log10(ENERGY_RANGE_GEV), n_points)*gd.UNIT_GEV
        # A scalar for the shared arm, not a list of identical values: handed N copies the
        # batched engine sees N baselines and declines, which is the opposite of what that
        # arm exists to show.
        baselines = (L_CHORD if baseline_kind == 'shared'
                     else list(np.linspace(0.5*L_CHORD, L_CHORD, n_points)))

        previous = os.sched_getaffinity(0)
        if pinned:
            os.sched_setaffinity(0, E_CORES)     # the workers inherit the mask
        print('    %s: %d worker counts, %s order, %d rounds%s'
              % (label, len(counts), order, rounds, ', pinned' if pinned else ''),
              file=sys.stderr, flush=True)
        try:
            rows = sweep(energies, baselines, counts, rounds, order)
        finally:
            os.sched_setaffinity(0, previous)

        one = [r for r in rows if r['n_jobs'] == 1][0]['mean_seconds']
        for row in rows:
            row['speedup'] = one/row['mean_seconds']
        store['arms'][label] = dict(
            n_points=n_points, baselines=baseline_kind, order=order, pinned=pinned,
            rounds=rounds, rows=rows,
            worst_cv=max(r['cv'] for r in rows if r['cv'] is not None))
        print('      best %.2fx at n_jobs=%d, worst cv %.3f'
              % (max(r['speedup'] for r in rows),
                 max(rows, key=lambda r: r['speedup'])['n_jobs'],
                 store['arms'][label]['worst_cv']), file=sys.stderr, flush=True)
        if not args.smoke:
            OUT.write_text(json.dumps(store, indent=1))

    if args.smoke:
        print('  smoke test only; nothing written', file=sys.stderr, flush=True)
        return
    OUT.write_text(json.dumps(store, indent=1))
    print('wrote %s (%d arms)' % (OUT.name, len(store['arms'])), file=sys.stderr,
          flush=True)


if __name__ == '__main__':
    main()
