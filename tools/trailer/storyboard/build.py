r"""Builds the trailer's storyboard, one panel per shot of ``trailer.json``, as a single HTML page.

A shot with ``"scene"`` shows that scene's clip (or its still, if no clip was rendered) from
``build/``; a shot with ``"still_from"`` shows the paper figure or notebook animation it will be
remade from, dimmed and tagged as such; a card shows its words.  Run ``data.py`` and
``render.py`` first.  PDF thumbnails need PyMuPDF::

    python tools/trailer/storyboard/build.py        # -> tools/trailer/build/storyboard.html
"""
import base64
import html
import io
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import BUILD, REPO  # noqa: E402

T = json.loads((HERE.parent / 'trailer.json').read_text())
SPB = 60.0 / T['bpm']
ACTS = {1: ('Two journeys', 'One neutrino, vacuum against matter', 'A low drone; a single tone per track', 'a'),
        2: ('The problem', 'A new solver for every setup and theory, and the question', 'Pulse enters, one hit per cell', 'a'),
        3: ('The reveal', 'The name and the promise', 'First full chord, then the groove', 'a'),
        4: ('How it works', 'The Hamiltonian, the code, and how Magnus picks its solver',
            'Rising arpeggio as the curve redraws', 'b'),
        5: ('What you can do', 'From two flavors to cosmic neutrinos', 'Full groove, building', 'b'),
        6: ('From textbook to frontier', "The paper's research examples",
            'Climax; the melody plays a computed probability', 'b'),
        7: ('Accurate, fast, flexible', 'Three words, each over its evidence', 'Three hits, resolving', 'c'),
        8: ('Get it', 'The call to action', 'Final chord, long fade', 'c')}
e = html.escape


def tc(s):
    return '%d:%04.1f' % (int(s // 60), s % 60)


def data_uri(data, mime):
    return 'data:%s;base64,%s' % (mime, base64.b64encode(data).decode())


def jpeg(im, size=1280, quality=80):
    im = im.convert('RGB')
    im.thumbnail((size, size))
    buf = io.BytesIO()
    im.save(buf, 'JPEG', quality=quality, optimize=True)
    return data_uri(buf.getvalue(), 'image/jpeg')


_paper = {}


def paper_still(path):
    """First page of a paper figure, or a frame 60 % into a notebook GIF, as a JPEG data URI."""
    if path not in _paper:
        from PIL import Image, ImageSequence
        if path.endswith('.pdf'):
            try:
                import pymupdf
            except ImportError:                     # PyMuPDF before 1.24 is imported as fitz
                import fitz as pymupdf
            page = pymupdf.open(REPO / path)[0]
            pix = page.get_pixmap(matrix=pymupdf.Matrix(900 / page.rect.width, 900 / page.rect.width), alpha=False)
            im = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
        else:
            frames = [f.convert('RGB') for f in ImageSequence.Iterator(Image.open(REPO / path))]
            im = frames[int(0.6 * (len(frames) - 1))]
        _paper[path] = jpeg(im, 720, 74)
    return _paper[path]


def scene_media(name, label):
    """The scene's clip, with its still as the poster; or the still alone."""
    from PIL import Image
    still = BUILD / 'stills' / (name + '.png')
    poster = jpeg(Image.open(still)) if still.exists() else ''
    for ext, mime in (('.mp4', 'video/mp4'), ('.webm', 'video/webm')):
        clip = BUILD / 'clips' / (name + ext)
        if clip.exists():
            return ('<video class="own" src="%s" poster="%s" autoplay muted loop playsinline aria-label="%s"></video>'
                    % (data_uri(clip.read_bytes(), mime), poster, e(label)))
    return '<img class="own" src="%s" alt="%s">' % (poster, e(label)) if poster else ''


def frame(s):
    words, kind = s['words'], s.get('kind')
    if kind == 'card':
        inner = ''.join('<div class="cw%s">%s</div>' % (' big' if i == 0 else '', e(w)) for i, w in enumerate(words))
        return '<div class="frame card">%s</div>' % inner
    tag, pic = '', ''
    if 'scene' in s:
        pic = scene_media(s['scene'], 'Scene for %s' % s['id'])
    elif 'still_from' in s:
        pic = '<img src="%s" alt="Source still for %s">' % (paper_still(s['still_from']), e(s['id']))
        tag = '<span class="tag">paper still, to be remade</span>'
    if 'code' in s:
        lines = ''.join('<span class="%s">%s</span>\n' % ('hl' if i in s.get('highlight', []) else '', e(ln) or ' ')
                        for i, ln in enumerate(s['code']))
        return ('<div class="frame split"><pre class="code">%s</pre><div class="pic">%s</div><div class="w">%s</div></div>'
                % (lines, pic, e(words[0])))
    if kind == 'pillar':
        return ('<div class="frame pillar">%s%s<div class="pw"><b>%s</b><span>%s</span></div></div>'
                % (pic, tag, e(words[0]), e(s['sub'])))
    if kind == 'schematic':
        return '<div class="frame">%s<div class="w top">%s</div></div>' % (pic, '<br>'.join(e(w) for w in words))
    return '<div class="frame">%s%s<div class="w">%s</div></div>' % (pic, tag, e(' '.join(words)))


def build():
    t, shots = 0.0, []
    for s in T['shots']:
        d = s['beats'] * SPB
        shots.append(dict(s, t0=t, t1=t + d))
        t += d
    total = t
    segs = ''.join('<div class="seg g%s" style="flex:%d" title="%s  %s"></div>'
                   % (ACTS[s['act']][3], s['beats'], e(s['id']), tc(s['t0'])) for s in shots)
    ticks = ''.join('<span style="left:%.3f%%">%s</span>' % (100 * x / total, tc(x)[:-2])
                    for x in range(0, int(total) + 1, 10))
    actlabels, acts = '', ''
    for a, (name, beat, music, g) in sorted(ACTS.items()):
        ss = [s for s in shots if s['act'] == a]
        actlabels += ('<div class="alab g%s" style="flex:%d"><b>%d</b><span>%s</span></div>'
                      % (g, sum(s['beats'] for s in ss), a, e(name)))
        cards = ''
        for s in ss:
            wide = ' full' if s.get('kind') == 'schematic' else (' wide' if (len(ss) == 1 or 'code' in s) else '')
            cards += ('<article class="shot%s">%s<div class="meta"><div class="tc"><span>%s &rarr; %s</span>'
                      '<span>%d beats &middot; %.1f s</span></div><p>%s</p><p class="src">%s</p></div></article>'
                      % (wide, frame(s), tc(s['t0']), tc(s['t1']), s['beats'], s['beats'] * SPB,
                         e(s.get('visual', '')), e(s.get('source', 'text card'))))
        acts += ('<section class="act g%s"><header><div class="num">Act %d</div><h2>%s</h2><div class="sub">'
                 '<span>%s &ndash; %s</span><span>%s</span><span class="mus">&#9834; %s</span></div></header>'
                 '<div class="grid">%s</div></section>'
                 % (g, a, e(name), tc(ss[0]['t0']), tc(ss[-1]['t1']), e(beat), e(music), cards))
    decisions = T.get('decisions', [])
    page = (HERE / 'template.html').read_text()
    for key, val in (('TOTAL', tc(total)[:-2]), ('NSHOTS', str(len(shots))), ('BPM', str(T['bpm'])),
                     ('SEGS', segs), ('TICKS', ticks), ('ACTLABELS', actlabels), ('ACTS', acts),
                     ('DECISIONS', ''.join('<li>%s</li>' % e(d.replace('{length}', tc(total)[:-2])) for d in decisions))):
        page = page.replace('{{%s}}' % key, val)
    out = BUILD / 'storyboard.html'
    out.write_text(page)
    print('%s: %d KB, %d shots, %s' % (out.relative_to(REPO), len(page) // 1024, len(shots), tc(total)))


if __name__ == '__main__':
    build()
