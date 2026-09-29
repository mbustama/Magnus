"""#91: re-time the 2nu and 3nu Magnus series of external_prem_chord_benchmarks.json into a copy."""
import json, pathlib, sys
sys.path.insert(0, '/home/mbustamante/Research/magnus/notebooks')
import gen_prem_benchmarks as g                      # the generator itself, unmodified
out = pathlib.Path(sys.argv[1])
store = json.loads(g.OUT.read_text())
for case in store['cases']:
    if case['flavours'] in (2, 3):
        case['series'] = [s for s in case['series']
                          if s['name'] not in ('Magnus', 'Magnus, order 6', 'Magnus, order 8')]
out.write_text(json.dumps(store, indent=1))
g.OUT = out                                          # it skips every series still present
g.main()
