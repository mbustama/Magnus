import os, sys, time, json, pathlib
import nbformat as nbf
from nbclient import NotebookClient
M = pathlib.Path(sys.argv[1]).resolve(); LOG = str(pathlib.Path(sys.argv[2]).resolve())
os.environ['ROUTE_LOG'] = LOG
os.environ['IPYTHONDIR'] = str(M.parent/'ipy')
os.environ['PYTHONPATH'] = str(M/'src')
os.environ['MAGNUS_PAPER_CACHE_ONLY'] = '1'
ORDER = [p for p in sys.argv[3].split(',') if p]
for path in [q for pre in ORDER for q in sorted((M/'notebooks').glob(pre + '*.ipynb')) if '.audited' not in q.name]:
    os.environ['NB_NAME'] = path.name
    nb = nbf.read(path, as_version=4)
    t = time.perf_counter()
    try:
        NotebookClient(nb, timeout=5400, kernel_name='python3',
                       resources={'metadata': {'path': str(path.parent)}}).execute()
        status = 'ok'
    except Exception as exc:
        status = 'FAILED %s' % str(exc).splitlines()[-1][:160] if str(exc) else 'FAILED'
    first = [o for c in nb.cells if c.cell_type == 'code' for o in c.get('outputs', []) if o.get('output_type') == 'stream'][:1]
    hook = ''.join(first[0].get('text', '')).splitlines()[0] if first else ''
    nbf.write(nb, M.parent/'audited'/path.name)
    print('%-50s %-8s %7.1f s  %s' % (path.name, status[:60], time.perf_counter() - t, hook[:90]), flush=True)
