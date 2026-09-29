"""#91: re-time the n_jobs = 1 rows of the shared-baseline arms of external_njobs_scaling.json."""
import json, os, pathlib, sys
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import gen_njobs_scaling as gn                       # sets one thread per process before numpy
import numpy as np
import magnus.globaldefs as gd
stored = json.loads(gn.OUT.read_text())
out = {}
for label, arm in stored['arms'].items():
    if arm.get('baselines') != 'shared':
        continue
    n = arm['n_points']
    rounds = arm.get('rounds') or gn.ROUNDS_BY_SIZE.get(n, gn.DEFAULT_ROUNDS)
    energies = np.logspace(*np.log10(gn.ENERGY_RANGE_GEV), n)*gd.UNIT_GEV
    previous = os.sched_getaffinity(0)
    if arm.get('pinned'):
        os.sched_setaffinity(0, gn.E_CORES)
    try:
        new = gn.sweep(energies, gn.L_CHORD, (1,), rounds, arm.get('order', 'randomised'))[0]
    finally:
        os.sched_setaffinity(0, previous)
    old = [r for r in arm['rows'] if r['n_jobs'] == 1][0]
    out[label] = dict(n_points=n, rounds=rounds, stored=old, new=new)
    print('%-32s %5d energies: n_jobs=1 mean %.4f -> %.4f s, best %.4f -> %.4f s, engine %s -> %s'
          % (label, n, old['mean_seconds'], new['mean_seconds'], old['seconds'], new['seconds'],
             old.get('engine'), new['engine']), flush=True)
pathlib.Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
