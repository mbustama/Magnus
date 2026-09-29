"""Merge each copy's recomputed sections into one cache (scratch), for the diff against the committed one."""
import json, pathlib
REPO = pathlib.Path('/home/mbustamante/Research/magnus')
S = pathlib.Path('/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad')
merged = json.loads((REPO/'notebooks/paper_figure_cache.json').read_text())
for d in ['nb28_w0', 'nb28_w1', 'nb28_w2', 'nb28_w3', 'nb28_solo']:
    got = json.loads((S/('par_' + d)/'notebooks'/'paper_figure_cache.json').read_text())
    for s in (S/('par_' + d)/'notebooks'/'redo_only.txt').read_text().split(','):
        merged[s] = got[s]
(S/'paper_figure_cache_recomputed.json').write_text(json.dumps(merged, indent=1))
print('merged', len(merged), 'sections')
