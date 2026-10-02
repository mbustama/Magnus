"""Add the tolerance points rtol = 1e-5, 1e-7 and 1e-9 to the stored Earth-plane series.

Reads the output of

    python notebooks/gen_prem_plane_magnus.py timed OUT.json --series tolerance --knobs=-4,-5,-6,-7,-8,-9

in which -4, -6 and -8 are controls: they are already stored, are not written, and are only
compared with the stored time and error.  Every stored point stays as it is.  Nothing is
written unless every control is within MAX_CONTROL_DRIFT of its stored time, and then only with
--write.
"""
import json
import pathlib
import sys

NB = pathlib.Path('/home/mbustamante/Research/magnus/notebooks')
NAME = 'Magnus (tolerance)'
MAX_CONTROL_DRIFT = 0.10

run = json.loads(pathlib.Path(sys.argv[1]).read_text())
write = sys.argv[2:] == ['--write']
files, ratios, added, worst_cv = [], {}, [], 0.0
for panel, fname, key in (('3nu', 'external_earth_plane.json', None),
                          ('3+1', 'external_prem_speed_accuracy_new.json', 'sterile_3plus1')):
    path = NB/fname
    data = json.loads(path.read_text())
    holder = data[key] if key else data
    floor = data['_magnus_provenance']['reference_confirmed_to']['disjoint_high_precision_extrapolation']
    series = next(s for s in holder['series'] if s['name'] == NAME)
    stored = {q['knob']: q for q in series['points']}
    flags_before = {q['knob']: q['best_at_this_accuracy'] for q in series['points']}
    for q in run['panels'][panel][NAME]:
        if q['knob'] in stored:
            s = stored[q['knob']]
            r = q['us_per_probability']/s['us_per_probability']
            ratios['%s %s' % (panel, q['label'])] = round(r, 3)
            same = float('%.5g' % q['max_abs_error']) == s['max_abs_error']
            print('%-4s control %s  time x%.3f  error %s' % (panel, q['label'], r,
                                                             'same' if same else 'DIFFERS'))
            if not same:
                raise SystemExit('a control error differs from the stored one: nothing written')
            continue
        worst_cv = max(worst_cv, q['block_cv'])
        series['points'].append(dict(
            at_reference_floor=q['max_abs_error'] <= floor, best_at_this_accuracy=True,
            block_cv=round(q['block_cv'], 4), knob=q['knob'], label=q['label'],
            max_abs_error=float('%.5g' % q['max_abs_error']), n_slabs=None, rtol=q['rtol'],
            us_per_probability=round(q['us_per_probability'], 4), us_sd=round(q['us_sd'], 4)))
        added.append('%s %s' % (panel, q['label']))
        print('%-4s added   %s  err %.4g  %.1f us/prob (cv %.3f)' % (
            panel, q['label'], q['max_abs_error'], q['us_per_probability'], q['block_cv']))
    series['points'].sort(key=lambda q: -q['knob'])
    best = {}                              # emit_figures.py's rule: the fastest at each error
    for q in series['points']:
        k = '%.6e' % q['max_abs_error']
        if k not in best or q['us_per_probability'] < best[k]['us_per_probability']:
            best[k] = q
    for q in series['points']:
        q['best_at_this_accuracy'] = q is best['%.6e' % q['max_abs_error']]
        if q['knob'] in flags_before and q['best_at_this_accuracy'] != flags_before[q['knob']]:
            raise SystemExit('the new points change the flag of stored point %s' % q['label'])
    print('%-4s series now: %s' % (panel, [(q['label'], q['best_at_this_accuracy'])
                                           for q in series['points']]))
    files.append((path, data))

drift = max(abs(r - 1.0) for r in ratios.values())
print('controls (time this run / time stored):', ratios, ' worst drift %.3f' % drift)
if drift > MAX_CONTROL_DRIFT:
    raise SystemExit('a control moved by more than %.0f%%: nothing written' % (100*MAX_CONTROL_DRIFT))
print('added:', added, ' worst block cv %.4f' % worst_cv)
for path, data in files:
    data['_magnus_provenance']['added_tolerance_points'] = dict(
        date=run['date'],
        why=('rtol 1e-5, 1e-7 and 1e-9 fill the gaps of the decade sweep; the stored points '
             'were not re-timed'),
        magnus_commit='35dcdd5 (src/ unchanged at %s, where the run was made)'
                      % run['repository_head'],
        machine='Intel Core i5-1334U, 12 logical CPUs, governor powersave (intel_pstate), '
                'on mains power',
        driver=('notebooks/gen_prem_plane_magnus.py timed --series tolerance '
                '--knobs=-4,-5,-6,-7,-8,-9'),
        protocol='manifest.json AMORTIZED (25-step delta_CP scan), 0.25 s blocks, 30 to 100 blocks',
        worst_block_cv=round(worst_cv, 4),
        load_average_at_start=[round(x, 2) for x in run['load_average_at_start']],
        controls=dict(note=('rtol 1e-4, 1e-6 and 1e-8, timed in the same run, interleaved, '
                            'and not written: time this run over time stored'), ratios=ratios),
        points='rtol 1e-5, 1e-7 and 1e-9, both panels')
    if write:
        path.write_text(json.dumps(data, indent=1) + '\n')
print('written' if write else 'dry run: pass --write to write')
