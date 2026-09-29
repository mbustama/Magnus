"""#92 in parallel: four copies recompute disjoint sets of value-only sections; the timed
prob_vs sections go to one solo copy run afterwards on a quiet machine.  Patches the copies only."""
import json, pathlib, shutil
REPO = pathlib.Path('/home/mbustamante/Research/magnus')
S = pathlib.Path('/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad')
cache = json.loads((REPO/'notebooks/paper_figure_cache.json').read_text())
timed = {k for k in cache if k.startswith('prob_vs')}
skip = {'what', 'timings', 'prem_timings'} | timed
jobs = sorted(((v.get('seconds') or 0.0, k) for k, v in cache.items() if isinstance(v, dict) and k not in skip), reverse=True)
load, sets = [0.0]*4, [[] for _ in range(4)]
for secs, k in jobs:                               # longest first, onto the least loaded copy
    w = load.index(min(load)); load[w] += secs; sets[w].append(k)
OLD_T = "if stored.get('fingerprint') == key and not os.environ.get('MAGNUS_PAPER_RETIME'):"
NEW_T = ("if (stored.get('fingerprint') == key or (os.environ.get('MAGNUS_PAPER_NO_TIMING') and 'rows' in stored "
         "and print('  TIMING-SKIPPED: timings configuration moved; stored rows kept for a quiet pass') is None)) "
         "and not os.environ.get('MAGNUS_PAPER_RETIME'):")
ANCHOR = "    key = fingerprint(*key_parts)\n"
def make(name, sections):
    top = S/('par_' + name)          # the copy sits beside a docs link, as notebooks/ sits in the repo
    if top.exists(): shutil.rmtree(top)
    d = top/'notebooks'
    shutil.copytree(REPO/'notebooks', d); (d/'figs').mkdir(); (top/'docs').symlink_to(REPO/'docs')
    p = d/'28_magnus_paper_figures.ipynb'; nb = json.loads(p.read_text())
    n_cached = n_timing = 0
    for c in nb['cells']:
        if c['cell_type'] != 'code': continue
        src = ''.join(c['source'])
        if "def cached(section, key_parts, compute, what=''):" in src:
            assert src.count(ANCHOR) == 1 and src.count("not os.environ.get('MAGNUS_PAPER_REDO')") == 2
            src = src.replace(ANCHOR, ANCHOR + "    _REDO = bool(os.environ.get('MAGNUS_PAPER_REDO')) and section in os.environ.get('MAGNUS_PAPER_REDO_ONLY', '').split(',')\n")
            src = src.replace("not os.environ.get('MAGNUS_PAPER_REDO')", 'not _REDO'); n_cached += 1
        if OLD_T in src:
            src = src.replace(OLD_T, NEW_T); n_timing += 1
        c['source'] = src.splitlines(keepends=True)
    assert n_cached == 1 and n_timing == 1, (n_cached, n_timing)
    p.write_text(json.dumps(nb, indent=1))
    (d/'redo_only.txt').write_text(','.join(sections))
for w in range(4):
    make('nb28_w%d' % w, sets[w]); print('copy %d: %2d sections, %.2f h recorded' % (w, len(sets[w]), load[w]/3600))
make('nb28_solo', sorted(timed)); print('solo : %s' % sorted(timed))
