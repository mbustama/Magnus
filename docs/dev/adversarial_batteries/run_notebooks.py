# -*- coding: utf-8 -*-
"""Re-execute notebooks 02 and 03, hash every embedded PNG, and report what moved.

This is Battery 1.4, the gap `FINDINGS_ADVERSARIAL_VALIDATION.md` 7 recorded as not covered.

Two traps from the handover are handled here:

  * `nbclient` dispatches through `async_execute_cell`; overriding the *sync* `execute_cell`
    records **nothing** while the run still reports success.  Per-cell timing therefore hooks
    the async method.  (Getting this wrong cost a full 14-minute pass.)
  * A figure changing is not by itself a regression, and a finer grid is not automatically
    better -- error is not monotone in `n_slabs`.  This script only *identifies* which figures
    moved; each one still has to be scored against `solve_ivp` separately.

The notebooks are executed against a COPY, so the committed .ipynb files are never modified.

Usage:  python run_notebooks.py [02] [03]
"""

import hashlib
import json
import os
import pathlib
import sys
import tempfile
import time

import nbformat
from nbclient import NotebookClient

REPO = pathlib.Path(__file__).resolve().parents[3]
NB_DIR = REPO/'notebooks'
NAMES = {'02': '02_magnus_2nu_vacuum_matter', '03': '03_magnus_3nu_vacuum_matter'}


def png_hashes(nb):
    """{cell_index: [sha256 of each embedded PNG]} for a notebook node or dict."""
    cells = nb['cells'] if isinstance(nb, dict) else nb.cells
    out = {}
    for i, c in enumerate(cells):
        hs = []
        for o in (c.get('outputs', []) if isinstance(c, dict) else getattr(c, 'outputs', [])):
            data = o.get('data', {}) if hasattr(o, 'get') else {}
            if 'image/png' in data:
                payload = data['image/png']
                if isinstance(payload, list):
                    payload = ''.join(payload)
                hs.append(hashlib.sha256(payload.encode()).hexdigest()[:16])
        if hs:
            out[i] = hs
    return out


def run(tag):
    name = NAMES[tag]
    src = NB_DIR/(name + '.ipynb')
    print('=== %s ===' % name, flush=True)

    stored = json.loads(src.read_text())
    stored_h = png_hashes(stored)

    nb = nbformat.read(str(src), as_version=4)
    client = NotebookClient(nb, timeout=1800, kernel_name='python3',
                            resources={'metadata': {'path': str(NB_DIR)}},
                            allow_errors=False)

    # Per-cell timing must hook the ASYNC method; the sync one is never called.
    timings = {}
    original = client.async_execute_cell

    async def timed(cell, index, *a, **k):
        t0 = time.time()
        res = await original(cell, index, *a, **k)
        timings[index] = time.time() - t0
        return res

    client.async_execute_cell = timed

    t0 = time.time()
    try:
        client.execute()
    except Exception as ex:                    # noqa: BLE001
        print('  EXECUTION FAILED: %s: %s' % (type(ex).__name__, str(ex)[:300]))
        return
    total = time.time() - t0
    print('  executed cleanly in %.1f s (%d cells)' % (total, len(nb.cells)))

    # Outside the repo: an executed notebook carries megabytes of embedded PNGs, and it is a
    # byproduct of the check rather than a deliverable.  Override with $NB_RERUN_OUT.
    out_dir = pathlib.Path(os.environ.get('NB_RERUN_OUT', tempfile.gettempdir()))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir/('nb_%s_rerun.ipynb' % tag)
    nbformat.write(nb, str(out_path))

    new_h = png_hashes(nb)
    print('  embedded PNGs: stored %d cells / %d images, rerun %d cells / %d images'
          % (len(stored_h), sum(len(v) for v in stored_h.values()),
             len(new_h), sum(len(v) for v in new_h.values())))

    changed, added, removed = [], [], []
    for i in sorted(set(stored_h) | set(new_h)):
        a, b = stored_h.get(i), new_h.get(i)
        if a is None:
            added.append(i)
        elif b is None:
            removed.append(i)
        elif a != b:
            changed.append(i)
    print('  figures CHANGED: %s' % (changed or 'none'))
    print('  figures added  : %s' % (added or 'none'))
    print('  figures removed: %s' % (removed or 'none'))

    slow = sorted(timings.items(), key=lambda kv: -kv[1])[:8]
    print('  slowest cells  : %s'
          % ', '.join('#%d %.1fs' % (i, t) for i, t in slow))

    # Warning output, which is what first indicated "the machine, not the change" last time.
    warn_cells = []
    for i, c in enumerate(nb.cells):
        for o in getattr(c, 'outputs', []):
            txt = o.get('text', '') if hasattr(o, 'get') else ''
            if isinstance(txt, list):
                txt = ''.join(txt)
            if 'Warning' in str(txt):
                warn_cells.append(i)
                break
    print('  cells emitting warnings: %s' % (warn_cells or 'none'))
    print('  wrote %s' % out_path, flush=True)


if __name__ == '__main__':
    for tag in (sys.argv[1:] or ['02', '03']):
        run(tag)
