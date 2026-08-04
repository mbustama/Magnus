# -*- coding: utf-8 -*-
"""Alternating timing, against `main` and against the branch before these fixes.

The handover is emphatic about method, and every rule it gives is applied here:

  * **Alternate.** This machine moved `test_adiabatic.py` 124 s -> 69 s between two runs of
    identical code, and one whole-suite comparison (392.9 s vs 331 s) was pure artefact.  Trees
    are therefore interleaved round-robin, never run back to back.
  * **Carry a control.** Two workloads the branch cannot touch (a vacuum scan and a
    constant-density scan, both excluded from the cumulative scan by construction) are timed in
    the same rounds.  Any effect on the workloads under test that is not large compared with the
    control's own scatter is noise.
  * **Report the spread, not one number.**  min and median of every round are shown.

Usage:
    python timing.py <label>=<path-to-src> [<label>=<path-to-src> ...] [--rounds N]

The child process imports magnus from the given src directory; the parent only alternates and
tabulates, so no two trees are ever live in one interpreter.
"""

import json
import subprocess
import sys

WORKLOADS = r'''
import json, sys, time, warnings
import numpy as np
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.matter as matter, magnus.oscprob as oscprob

L_SCALE, R_SUN, NE0 = gd.L_SCALE_SUN, gd.SUN_RADIUS*gd.UNIT_KM, gd.NUM_DENSITY_E_SUN_CENTRAL
_p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
P2 = {'sth': _p['s12'], 'Dm2': _p['D21']}
P3 = {k: _p[k] for k in ('s12','s23','s13','dCP','D21','D31')}
solar = matter.exp_density_profile(NE0, L_SCALE)
E, LM = 50.0e6, 0.5*R_SUN

def modulated(l):
    x = np.asarray(l, dtype=float)
    y = NE0*np.exp(-x/L_SCALE)*(1.0 + 0.9*np.sin(2.0*np.pi*6.0*x/LM))
    a = np.asarray(y)
    return a[()] if a.ndim == 0 else a

def std(*a, **k):
    return oscprob.osc_prob_matter_std_potential(
        *a, L0=0.0, density_is_of_number_of_electrons=True, **k)

CASES = {
    'single point, solar (hybrid)':       lambda: std(2, solar, E, LM, P2),
    'single point, 3nu solar':            lambda: std(3, solar, E, LM, P3),
    'single point, multi-resonance':      lambda: std(2, modulated, E, LM, P2),
    'solar scan N=400 (cumulative)':      lambda: std(2, solar, E, np.linspace(0.05*LM, LM, 400), P2),
    'solar scan N=8 (hybrid)':            lambda: std(2, solar, E, np.linspace(0.05*LM, LM, 8), P2),
    'CONTROL vacuum scan N=300':          lambda: std(3, 0.0, E, np.linspace(0.05*LM, LM, 300), P3),
    'CONTROL const-density scan N=300':   lambda: std(3, 0.05*NE0, E, np.linspace(0.05*LM, LM, 300), P3),
}

out = {}
for name, fn in CASES.items():
    fn()                       # warm any caches so the first round is not special
    t0 = time.time(); fn(); out[name] = time.time() - t0
print(json.dumps(out))
'''


def one_round(src):
    r = subprocess.run([sys.executable, '-c', WORKLOADS], capture_output=True, text=True,
                       env={'PYTHONPATH': src, 'PATH': '/usr/bin:/bin',
                            'HOME': '/home/mbustamante'}, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


def main():
    argv = list(sys.argv[1:])
    rounds = 5
    if '--rounds' in argv:
        i = argv.index('--rounds')
        rounds = int(argv[i + 1])
        del argv[i:i + 2]
    trees = [a.split('=', 1) for a in argv]
    if len(trees) < 2:
        print(__doc__)
        return

    acc = {label: {} for label, _ in trees}
    for rnd in range(rounds):
        for label, src in trees:                      # round-robin: never back to back
            t = one_round(src)
            for k, v in t.items():
                acc[label].setdefault(k, []).append(v)
        print('  round %d/%d done' % (rnd + 1, rounds), flush=True)

    import statistics as st
    names = list(acc[trees[0][0]])
    w = max(len(n) for n in names) + 2
    hdr = ' '*w + ''.join('%22s' % lab for lab, _ in trees)
    print('\n' + hdr)
    print(' '*w + ''.join('%22s' % '(min / median, ms)' for _ in trees))
    print('-'*len(hdr))
    base = trees[0][0]
    for n in names:
        row = '%-*s' % (w, n)
        for lab, _ in trees:
            v = sorted(acc[lab][n])
            row += '%22s' % ('%.1f / %.1f' % (1e3*v[0], 1e3*st.median(v)))
        print(row)

    print('\nratio of medians, relative to %r  (criterion: no configuration >1.10x)' % base)
    for n in names:
        b = st.median(acc[base][n])
        parts = []
        for lab, _ in trees[1:]:
            parts.append('%s %.2fx' % (lab, st.median(acc[lab][n])/b))
        print('  %-*s %s' % (w, n, '  '.join(parts)))

    print('\nmachine scatter on the CONTROL workloads (max/min within one tree):')
    for lab, _ in trees:
        for n in names:
            if n.startswith('CONTROL'):
                v = sorted(acc[lab][n])
                print('  %-10s %-34s %.2fx' % (lab, n, v[-1]/v[0]))


if __name__ == '__main__':
    main()
