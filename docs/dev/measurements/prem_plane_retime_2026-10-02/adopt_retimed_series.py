"""Adopt the re-timed Magnus series into the stored Earth-plane files (2026-10-02, a978203).

Run once from the repository root, inline, as

    python - <<'EOF' ... EOF

on the output of `python notebooks/gen_prem_plane_magnus.py timed /tmp/magnus_timed.json`,
which is `magnus_timed.json` in this folder.  It replaces the two Magnus series of both
Earth panels, with the fields, rounding and best_at_this_accuracy rule of the #93
adoption (`../issue091_093_adoption/adopt_91_93.py`), and writes the record
`_magnus_provenance.retimed_series`.  Kept as it was run.
"""
import json, pathlib
NB = pathlib.Path('notebooks')
t = json.loads(pathlib.Path('/tmp/magnus_timed.json').read_text())

def indent_of(path):
    lines = path.read_text().splitlines()
    return len(lines[1]) - len(lines[1].lstrip(' ')) if len(lines) > 1 else None

def write(path, data):
    path.write_text(json.dumps(data, indent=indent_of(path)) + ('\n' if path.read_text().endswith('\n') else ''))

worst_cv = max(q['block_cv'] for pan in t['panels'].values() for s in pan.values() for q in s)
for panel, fname, key in (('3nu', 'external_earth_plane.json', None),
                          ('3+1', 'external_prem_speed_accuracy_new.json', 'sterile_3plus1')):
    p = NB/fname
    data = json.loads(p.read_text())
    holder = data[key] if key else data
    floor = data['_magnus_provenance']['reference_confirmed_to']['disjoint_high_precision_extrapolation']
    for name in ('Magnus', 'Magnus (tolerance)'):
        old = next(s for s in holder['series'] if s['name'] == name)
        fields = list(old['points'][0])
        pts = []
        for q in t['panels'][panel][name]:
            d = dict(at_reference_floor=q['max_abs_error'] <= floor, best_at_this_accuracy=True,
                     block_cv=round(q['block_cv'], 4), knob=q['knob'], label=q['label'],
                     max_abs_error=float('%.5g' % q['max_abs_error']), n_slabs=q['n_slabs'])
            if 'rtol' in fields:
                d['rtol'] = q['rtol']
            d.update(us_per_probability=round(q['us_per_probability'], 4), us_sd=round(q['us_sd'], 4))
            assert list(d) == fields, (list(d), fields)
            pts.append(d)
        best = {}                          # emit_figures.py's rule: the fastest at each error
        for q in pts:
            k = '%.6e' % q['max_abs_error']
            if k not in best or q['us_per_probability'] < best[k]['us_per_probability']:
                best[k] = q
        for q in pts:
            q['best_at_this_accuracy'] = q is best['%.6e' % q['max_abs_error']]
        print(panel, name, 'flags old', [(q['label'], q['at_reference_floor'], q['best_at_this_accuracy']) for q in old['points']])
        print(panel, name, 'flags new', [(q['label'], q['at_reference_floor'], q['best_at_this_accuracy']) for q in pts])
        old['points'] = pts
    prov = data['_magnus_provenance']
    prov['retimed_tolerance_series']['superseded_by'] = 'retimed_series (2026-10-02)'
    prov['retimed_series'] = dict(
        date=t['date'],
        why=('the slab-count series, measured on 2026-09-02 at 6751a8c, was stale: Magnus became '
             '1.6 to 3.9 times faster on it since, and its 1-slab errors changed; both series were '
             're-timed together on the paper\'s laptop'),
        magnus_commit=t['repository_head'],
        machine='Intel Core i5-1334U, 12 logical CPUs, governor powersave (intel_pstate), on mains power',
        driver='notebooks/gen_prem_plane_magnus.py timed, after PR #210 (n/p ratio passed at 3+1 only)',
        protocol='manifest.json AMORTIZED (25-step delta_CP scan), 0.25 s blocks, 30 to 100 blocks',
        worst_block_cv=round(worst_cv, 4), load_average_at_start=[round(x, 2) for x in t['load_average_at_start']],
        series='both Magnus series, slab count and tolerance, every point')
    write(p, data)
    print('wrote', fname)
