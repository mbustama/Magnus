# -*- coding: utf-8 -*-
"""Is the closed form's cost peak the commutator term, or something else?

RESULT, measured 2026-09-11 on the 23-width scan: it is the commutator term.  Across the
scan the slab count rises and falls by 25.9x at order 2 and by 17.3x for the closed form,
both peaking at the same width, 38 km; at order 4 the same ratio is 6.0x.  Order 2 carries
exactly the information per slab that the closed form does and no commutator, so this
separates the commutator from every other difference between the two codes.

The claim is that a constant-density slab cannot represent the matter eigenbasis rotating
inside it, that the rotation is what the commutator term carries, and that Magnus order 4
therefore pays a much smaller price at the adiabatic crossover.

Magnus order 2 IS the midpoint slab product: same information per slab as the closed form,
no commutator.  So if the claim is right, order 2's slab count should peak at the crossover
the way the closed form's does, and order 4's should not.  Slab counts, not times -- the
count is machine-independent, and it is the count the claim is about.
"""
import json
import pathlib
import resource
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'src'))
resource.setrlimit(resource.RLIMIT_AS, (8_000_000_000,)*2)

import gen_shock_cost as G                                   # noqa: E402

TARGET = 1.0e-7
d = json.loads((HERE/'external_shock_cost.json').read_text())
rows = {c['width']: c for c in d['cases']
        if c['code'] == 'order 4, resolved' and c['target'] == TARGET}
npe = {c['width']: c for c in d['cases']
       if c['code'] == 'NuOscProbExact' and c['target'] == TARGET}

print('  %9s %10s %10s %10s' % ('width_km', 'order 2', 'order 4', 'npe'), flush=True)
out = {}
for w in sorted(rows):
    ref, _ = G.reference(w)
    n, err = G.smallest_n_reaching(
        lambda k, w=w: G.magnus_at(w, 2, k, True), ref, TARGET,
        rows[w]['n_slabs'], G.max_slabs_for(4))
    out[w] = n
    print('  %9.3f %10s %10d %10d'
          % (rows[w]['width_km'], n if n else 'X', rows[w]['n_slabs'], npe[w]['n_slabs']),
          flush=True)

good = {w: n for w, n in out.items() if n}
if good:
    print('\npeak/floor ratio across the scan:', flush=True)
    print('  order 2 %.1fx   order 4 %.1fx   npe %.1fx'
          % (max(good.values())/min(good.values()),
             max(c['n_slabs'] for c in rows.values())/min(c['n_slabs'] for c in rows.values()),
             max(c['n_slabs'] for c in npe.values())/min(c['n_slabs'] for c in npe.values())),
          flush=True)
    pk = max(good, key=good.get)
    print('  order 2 peaks at %.1f km; npe peaks at %.1f km'
          % (rows[pk]['width_km'],
             rows[max(npe, key=lambda w: npe[w]['n_slabs'])]['width_km']), flush=True)
