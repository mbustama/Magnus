"""After the #88 rebuild: what moved against HEAD, by pixels for figures and by text for outputs."""
import base64, io, json, os, subprocess, tempfile
import numpy as np
from PIL import Image

os.chdir('/home/mbustamante/Research/magnus')


def git_show(path):
    return subprocess.run(['git', 'show', 'HEAD:' + path], capture_output=True, check=True).stdout


def raster_pdf(data, dpi=100):
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'a.pdf')
        open(p, 'wb').write(data)
        subprocess.run(['pdftoppm', '-r', str(dpi), '-png', p, os.path.join(d, 'r')], check=True)
        return [np.asarray(Image.open(os.path.join(d, f)).convert('RGB'))
                for f in sorted(os.listdir(d)) if f.startswith('r')]


def pix(b64):
    return np.asarray(Image.open(io.BytesIO(base64.b64decode(b64))).convert('RGB'))


def frac(x, y):
    return 100*np.mean(np.any(x != y, axis=-1))


lines = subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True).stdout.splitlines()
paths = [l[3:] for l in lines if not l.startswith('??')]
print('modified: %d files; untracked: %s' % (len(paths), [l[3:] for l in lines if l.startswith('??')]))
restore = []
for p in paths:
    if p.endswith('.pdf'):
        old, new = raster_pdf(git_show(p)), raster_pdf(open(p, 'rb').read())
        if len(old) == len(new) and all(x.shape == y.shape for x, y in zip(old, new)):
            f = max(frac(x, y) for x, y in zip(old, new))
            print('%-60s %s' % (p, 'pixel-identical' if f == 0 else 'PIXELS DIFFER (%.4f%%)' % f))
            if f == 0:
                restore.append(p)
        else:
            print('%-60s PAGE SIZE OR COUNT DIFFERS' % p)
    elif p.endswith('.png'):
        x = np.asarray(Image.open(io.BytesIO(git_show(p))).convert('RGB'))
        y = np.asarray(Image.open(p).convert('RGB'))
        same = x.shape == y.shape and np.array_equal(x, y)
        print('%-60s %s' % (p, 'pixel-identical' if same else 'PIXELS DIFFER'
                            + ('' if x.shape != y.shape else ' (%.4f%%)' % frac(x, y))))
        if same:
            restore.append(p)
    elif not p.endswith('.ipynb'):
        print('%-60s (other file)' % p)

nbp = 'notebooks/28_magnus_paper_figures.ipynb'
if nbp in paths:
    a, b = json.loads(git_show(nbp)), json.load(open(nbp))
    print('\n%s: %d -> %d cells' % (nbp, len(a['cells']), len(b['cells'])))

    def text(c):
        out = []
        for o in c.get('outputs', []):
            if o.get('output_type') == 'stream':
                out.append(''.join(o.get('text', '')))
            elif 'data' in o and 'image/png' not in o['data'] and 'text/plain' in o['data']:
                out.append(''.join(o['data']['text/plain']))
        return ''.join(out)

    def images(c):
        return [o['data']['image/png'] for o in c.get('outputs', [])
                if 'data' in o and 'image/png' in o['data']]

    for i, (ca, cb) in enumerate(zip(a['cells'], b['cells'])):
        if ca['cell_type'] != 'code':
            continue
        if ''.join(ca['source']) != ''.join(cb['source']):
            print('cell %d: SOURCE DIFFERS' % i)
        if not ca.get('outputs'):
            if cb.get('outputs'):
                print('cell %d: outputs restored (%d text lines, %d images)'
                      % (i, len(text(cb).splitlines()), len(images(cb))))
            continue
        ta, tb = text(ca).splitlines(), text(cb).splitlines()
        if ta != tb:
            print('cell %d: text changed' % i)
            for x, y in [(x, y) for x, y in zip(ta, tb) if x != y][:4]:
                print('    - %s\n    + %s' % (x[:160], y[:160]))
            if len(ta) != len(tb):
                print('    line count %d -> %d' % (len(ta), len(tb)))
        ia, ib = images(ca), images(cb)
        if len(ia) != len(ib):
            print('cell %d: image count %d -> %d' % (i, len(ia), len(ib)))
            continue
        for k, (x, y) in enumerate(zip(ia, ib)):
            if x != y:
                X, Y = pix(x), pix(y)
                print('cell %d image %d: %s' % (i, k, 'bytes differ, pixels identical'
                      if X.shape == Y.shape and np.array_equal(X, Y) else 'PIXELS DIFFER'
                      + ('' if X.shape != Y.shape else ' (%.4f%%)' % frac(X, Y))))
print('\npixel-identical files (candidates to restore from HEAD):', restore)
