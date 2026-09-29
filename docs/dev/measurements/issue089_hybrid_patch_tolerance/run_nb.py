import os, sys, time, pathlib
import nbformat as nbf
from nbclient import NotebookClient
S = pathlib.Path(__file__).resolve().parent
M = S/'mirror'
os.environ.update(HP_LOG=str(S/'hp_notebooks.log'), HP_DIR=str(S), IPYTHONDIR=str(S/'ipy'),
                  PYTHONPATH=str(M/'src'), MAGNUS_PAPER_CACHE_ONLY='1')
order = [p for p in sys.argv[1].split(',') if p]
paths = [q for pre in order for q in sorted((M/'notebooks').glob(pre + '*.ipynb'))]
for path in paths:
    os.environ['NB_NAME'] = path.name
    nb = nbf.read(path, as_version=4)
    t = time.perf_counter()
    try:
        NotebookClient(nb, timeout=5400, kernel_name='python3',
                       resources={'metadata': {'path': str(path.parent)}}).execute()
        status = 'ok'
    except Exception as exc:
        status = 'FAILED ' + (str(exc).strip().splitlines() or [''])[-1][:150]
    nbf.write(nb, S/'audited'/path.name)
    print('%-48s %s %.1f s' % (path.name, status, time.perf_counter() - t), flush=True)
