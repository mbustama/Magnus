"""Adopt the re-measured #91 series and the #93 tolerance series into the stored files."""
import json, pathlib
NB = pathlib.Path('/home/mbustamante/Research/magnus/notebooks')
S = pathlib.Path('/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad')
WHY = ('the breakpoint-grid fix of the batched engine (#71: fe385c9, #94, #105) moved these '
       'points; re-timed on the paper\'s laptop')


def indent_of(path):
    lines = path.read_text().splitlines()
    return len(lines[1]) - len(lines[1].lstrip(' ')) if len(lines) > 1 else None


def write(path, data):
    path.write_text(json.dumps(data, indent=indent_of(path)) + ('\n' if path.read_text().endswith('\n') else ''))


# --- #91: PREM chord, 2nu and 3nu Magnus series ----------------------------------------------
p = NB/'external_prem_chord_benchmarks.json'
stored, new = json.loads(p.read_text()), json.loads((S/'prem_chord_rerun2.json').read_text())
names = ('Magnus', 'Magnus, order 6', 'Magnus, order 8')
for case in stored['cases']:
    if case['flavours'] not in (2, 3):
        continue
    fresh = {s['name']: s for c in new['cases'] if c['flavours'] == case['flavours'] for s in c['series']}
    case['series'] = [fresh[s['name']] if s['name'] in names else s for s in case['series']]
stored['magnus_rerun_71'] = dict(
    date='2026-09-28', why=WHY + '; at four and five flavors no point moves, so those series '
    'stay as measured on 2026-09-05', series='Magnus, orders 4, 6 and 8, at 2 and 3 flavors',
    control_ratio_this_run=0.995, magnus_commit='bd9f9d7',
    how='notebooks/gen_prem_benchmarks.py, unmodified, run against a copy without these series')
write(p, stored)
print('#91 adopted into', p.name)

# --- #93: the tolerance series of the two Earth planes -----------------------------------------
t93 = json.loads((S/'timed93_main.json').read_text())
worst_cv = max(q['block_cv'] for pan in t93['panels'].values() for q in pan['Magnus (tolerance)'])
for panel, fname, key in (('3nu', 'external_earth_plane.json', None),
                          ('3+1', 'external_prem_speed_accuracy_new.json', 'sterile_3plus1')):
    p = NB/fname
    data = json.loads(p.read_text())
    holder = data[key] if key else data
    floor = data['_magnus_provenance']['reference_confirmed_to']['disjoint_high_precision_extrapolation']
    old = next(s for s in holder['series'] if s['name'] == 'Magnus (tolerance)')
    print(panel, 'stored flags:', [(q['label'], q['at_reference_floor'], q['best_at_this_accuracy']) for q in old['points']])
    pts = [dict(at_reference_floor=q['max_abs_error'] <= floor, best_at_this_accuracy=True,
                block_cv=round(q['block_cv'], 4), knob=q['knob'], label=q['label'],
                max_abs_error=float('%.5g' % q['max_abs_error']), n_slabs=None, rtol=q['rtol'],
                us_per_probability=round(q['us_per_probability'], 4), us_sd=round(q['us_sd'], 4))
           for q in t93['panels'][panel]['Magnus (tolerance)']]
    best = {}                              # emit_figures.py's rule: the fastest at each error
    for q in pts:
        k = '%.6e' % q['max_abs_error']
        if k not in best or q['us_per_probability'] < best[k]['us_per_probability']:
            best[k] = q
    for q in pts:
        q['best_at_this_accuracy'] = q is best['%.6e' % q['max_abs_error']]
    old['points'] = pts
    data['_magnus_provenance']['retimed_tolerance_series'] = dict(
        date='2026-09-28', why=WHY, magnus_commit=t93['repository_head'],
        driver=('notebooks/gen_prem_plane_magnus.py, which reproduces all 16 stored errors of '
                'the tolerance series at 6751a8c'),
        protocol='manifest.json AMORTIZED (25-step delta_CP scan), 0.25 s blocks, 30 to 100 blocks',
        worst_block_cv=round(worst_cv, 4), load_average_at_start=t93['load_average_at_start'],
        untouched='the slab-count series, which the fix cannot move, as measured on 2026-09-02')
    write(p, data)
    print(panel, 'new flags:  ', [(q['label'], q['at_reference_floor'], q['best_at_this_accuracy']) for q in pts])
