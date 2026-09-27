r"""Assembles the whole trailer, shot by shot from ``trailer.json``, as 1920x1080 frames, and muxes
them with the music (``music.py``) into ``build/magnus_trailer.mp4``.

Every frame is drawn fresh: a shot with a ``scene`` draws that scene (``scenes.py``) at the shot's
progress; the code moment types its code and redraws the opening's matter curve; cards and the
"Flexible." pillar are drawn here.  Every scene is drawn below a title band at the top of the frame,
where its words go, one phrase after another (clear of a video player's controls).  Cuts are hard,
on the beat.
Run ``data.py`` first.  Uses every core it is given::

    nice -n 19 python tools/trailer/cut.py [--jobs 4] [--width 1920] [--frames-only]
"""
import multiprocessing as mp
import os
import subprocess
import sys

import numpy as np

from common import (BUILD, REPO, HERE, BG, INK, BLUE, AMBER, TEAL, VIOLET, MONO, DISP, CODE, FPS, SPB,
                    seg, setup_matplotlib)

T = __import__('json').loads((HERE / 'trailer.json').read_text())
HOLD = 0.12
ACT_COLOR = {1: BLUE, 2: BLUE, 3: BLUE, 4: AMBER, 5: AMBER, 6: AMBER, 7: TEAL, 8: TEAL}
WIDTH = 1920


def timeline():
    """[(shot, first frame, number of frames)], the frame counts rounded so that shot boundaries
    stay on the beat grid overall."""
    out, t = [], 0.0
    for s in T['shots']:
        f0, t = int(round(t * FPS)), t + s['beats'] * SPB
        out.append((s, f0, int(round(t * FPS)) - f0))
    return out


# ------------------------------------------------------------------ drawing helpers
def text(ax, x, y, s, size, col=INK, fam=MONO, alpha=1.0, **kw):
    if alpha > 0:
        return ax.text(x, y, s, fontsize=size, color=col, fontfamily=fam, alpha=alpha, **kw)


def logotype(fig, ax, x, y, size, alpha=1.0, ha='center'):
    """'Magνs' with the ν in amber, laid out piece by piece."""
    if alpha <= 0:
        return
    r = fig.canvas.get_renderer()
    fam = DISP
    pieces = [('Mag', INK), ('ν', AMBER), ('s', INK)]
    arts = [ax.text(0, y, p, fontsize=size, color=c, fontfamily=fam, weight=700, alpha=alpha, va='center')
            for p, c in pieces]
    inv = ax.transData.inverted()
    widths = [inv.transform((a.get_window_extent(r).width, 0))[0] - inv.transform((0, 0))[0] for a in arts]
    x0 = x - sum(widths) / 2 if ha == 'center' else x
    for a, w in zip(arts, widths):
        a.set_x(x0)
        x0 += w


CONTENT_FRAC = 0.85                         # scenes are drawn at 0.85 of the frame, below the band
BAND = 1.35                                 # the title band: the top 1.35 of the 9 units


def caption(ax, phrases, u, n_frames, act):
    """The scene's words in a band at the top of the frame, clear of a player's controls: one phrase
    after another, each fading in, with the act's color as an accent."""
    if not phrases:
        return
    from matplotlib.patches import Rectangle
    k = len(phrases)
    j = min(k - 1, int(u * k))
    fade_in = 0.3 / max(n_frames / FPS / k, 0.3)            # 0.3 s per phrase
    a = seg(u * k - j, 0, fade_in)
    ax.add_patch(Rectangle((0.8, 9 - BAND / 2 - 0.3), 0.08, 0.6, fc=ACT_COLOR[act], ec='none', zorder=21))
    text(ax, 1.1, 9 - BAND / 2, phrases[j], 28, INK, DISP, alpha=a, weight=700, va='center',
         zorder=21)


# ------------------------------------------------------------------ shots drawn here
def _glow():
    """The cards' ground: a soft radial glow from #16203a at the center to the background."""
    y, x = np.mgrid[0:90, 0:160]
    r = np.hypot((x - 80) / 95.0, (y - 40) / 60.0)
    w = np.clip(1 - r, 0, 1)[..., None] ** 1.6
    c0, c1 = np.array([0x06, 0x08, 0x0d]) / 255, np.array([0x16, 0x20, 0x3a]) / 255
    return c0 + (c1 - c0) * w


GLOW = _glow()


def card(fig, ax, s, u, n):
    words, sid = s['words'], s['id']
    beats = s['beats']
    at = lambda b: seg(u * beats, b, b + 0.6)             # noqa: E731  fade in over 0.6 beat from beat b
    ax.imshow(GLOW, extent=(0, 16, 0, 9), zorder=0, interpolation='bilinear')
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    if sid == 'reveal':
        logotype(fig, ax, 8, 5.9, 100, alpha=at(0))      # a fade only: a growing logotype jitters
        for i, w in enumerate(words[1:]):
            last = i == len(words) - 2
            text(ax, 8, 4.4 - 0.68 * i, w, 25 if last else 24, TEAL if last else INK, DISP,
                 alpha=at(3 + 2 * i), ha='center', va='center', weight=700 if last else 500)
    elif sid == 'cta':
        typed = words[0][:int(len(words[0]) * seg(u * beats, 0.3, 3))]
        text(ax, 8, 6.9, '$ ' + typed, 40, TEAL, MONO, ha='center', va='center')
        for i, w in enumerate(words[1:6]):
            text(ax, 8, 5.75 - 0.58 * i, w, 21, INK, MONO, alpha=at(3.5 + 1.2 * i), ha='center', va='center')
        logotype(fig, ax, 8, 2.35, 70, alpha=at(9))
        text(ax, 8, 1.25, 'Accurate  ·  Fast  ·  Flexible', 24, AMBER, DISP, alpha=at(10),
             ha='center', va='center')
    else:
        for i, w in enumerate(words):
            text(ax, 8, 5.4 - 1.15 * i + 0.575 * (len(words) - 1), w, 46,
                 INK, DISP, alpha=at(2 * i), ha='center', va='center', weight=700)


CODE_PT = 16


def code_moment(fig, ax, s, u, n):
    import scenes
    from matplotlib.patches import FancyBboxPatch
    lines = s['code']
    total = sum(len(ln) for ln in lines)
    shown = int(total * seg(u, 0.02, 0.5))
    ax.add_patch(FancyBboxPatch((0.6, 2.0), 8.4, 6.2, boxstyle='round,pad=0,rounding_size=0.15', fc='#0e131e',
                                ec='#232c3f', lw=1.5))
    left = shown
    for i, ln in enumerate(lines):
        vis, left = ln[:max(0, left)], left - len(ln)
        y = 7.5 - 0.62 * i
        if i in s.get('highlight', []) and u > 0.52:
            ax.add_patch(FancyBboxPatch((0.75, y - 0.26), 8.1, 0.52, boxstyle='round,pad=0,rounding_size=0.05',
                                        fc=AMBER, ec='none', alpha=0.18 * seg(u, 0.52, 0.6)))
        t = text(ax, 0.95, y, vis if vis else ' ', CODE_PT, '#cfd6e4', CODE, va='center')
        if 0 < left + len(ln) <= len(ln) and u < 0.5 and int(u * 60) % 2 == 0:
            x1 = ax.transData.inverted().transform(t.get_window_extent(fig.canvas.get_renderer()))[1, 0]
            text(ax, x1 + 0.05, y, '|', CODE_PT, TEAL, CODE, va='center', weight=700)
    d = scenes.load('opening.npz')
    L, Pm = d['L'], d['Pm']
    g = seg(u, 0.55, 0.95)
    a = fig.add_axes([0.6, 0.26, 0.36, 0.62])
    scenes.plain(a)
    a.set_xlim(0, L[-1])
    top = 1.15 * Pm.max()
    a.set_ylim(-0.02, top)
    a.plot(L, Pm, color=AMBER, lw=1.2, alpha=0.12)
    k = int(g * (len(L) - 1))
    if k > 0:
        a.plot(L[:k + 1], Pm[:k + 1], color=AMBER, lw=3)
        a.fill_between(L[:k + 1], 0, Pm[:k + 1], color=AMBER, alpha=0.1, lw=0)
    a.set_xticks([0, 5000, 10000])
    a.set_xticklabels(['0', '5000', '10,000 km'])
    scenes.monoticks(a)
    a.text(0, top, r'$P(\nu_\mu \to \nu_e)$ from this call', color=INK, fontfamily=MONO, fontsize=17, va='top',
           alpha=max(0.35, g))


def flexible(fig, ax, s, u, n):
    from matplotlib.patches import FancyBboxPatch
    chips = [('Vacuum', BLUE), ('Matter', AMBER), ('Constant or not', AMBER), ('Standard Model or beyond', VIOLET),
             ('2 to 5 flavors, ready-made', TEAL)]
    for i, (c, col) in enumerate(chips):
        p = seg(u, 0.05 + 0.1 * i, 0.25 + 0.1 * i)
        y = 7.3 - 0.95 * i
        x = 8.8 + 6 * (1 - p)                          # slide in from the right and lock
        w = 0.19 * len(c) + 0.9
        ax.add_patch(FancyBboxPatch((x, y - 0.33), w, 0.66, boxstyle='round,pad=0,rounding_size=0.33',
                                    fc='#101626', ec=col, lw=2, alpha=p))
        text(ax, x + w / 2, y, c, 19, col, MONO, alpha=p, ha='center', va='center')
    text(ax, 0.8, 4.6, r'$H(E,\,x)$', 64, INK, alpha=seg(u, 0, 0.2))


def pillar(fig, ax, s, u, n):
    a = seg(u, 0, 0.18)
    text(ax, 0.8, 9 - BAND / 2, s['words'][0], 44, TEAL, DISP, alpha=a, weight=700, va='center')
    import textwrap
    text(ax, 5.2 if len(s['words'][0]) < 10 else 6.4, 9 - BAND / 2, textwrap.fill(s['sub'], 52), 18, INK, MONO,
         alpha=seg(u, 0.1, 0.3), va='center', linespacing=1.5)


# ------------------------------------------------------------------ one frame of the cut
_plt = None


# Scenes are drawn below the title band: a box of the frame's own aspect, 0.85 of its height.
CONTENT = ((1 - CONTENT_FRAC) / 2, 0.0, CONTENT_FRAC, CONTENT_FRAC)


def in_content(fig, draw, dy=0.0):
    """Runs draw(fig, ax) with ax covering CONTENT (raised by dy, a fraction of the frame) in the
    scene's 16 x 9 units, and every axes the scene adds (given in fractions of that box) mapped
    into it too."""
    add = fig.add_axes
    cx, cy, cw, ch = CONTENT
    cy += dy

    def mapped(rect, **kw):
        x, y, w, h = rect
        return add([cx + x * cw, cy + y * ch, w * cw, h * ch], **kw)

    ax = add([cx, cy, cw, ch])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis('off')
    fig.add_axes = mapped
    try:
        draw(fig, ax)
    finally:
        fig.add_axes = add


# The scenes were laid out with titles of their own at the top, now moved into the band, which left
# them sitting low.  Each is measured once, on its finished picture (u = 1): drawn at low
# resolution, the first and last rows of pixels that differ from the background are its extent.
# (Axes' bounding boxes cannot be used for this: an empty axes reports its full extent.)  It is
# then raised so that its top sits TOP_GAP below the band; one offset per shot, the same for every
# frame, so nothing moves from frame to frame.
BAND_BOTTOM = 1 - BAND / 9                    # where the band ends, in fractions of the frame
TOP_GAP = 0.03
_offsets = {}


def content_offset(key, draw):
    if key not in _offsets:
        fig = _plt.figure(figsize=(16, 9), dpi=40, facecolor=BG)
        in_content(fig, draw)
        fig.canvas.draw()
        img = np.asarray(fig.canvas.buffer_rgba())[..., :3].astype(int)
        _plt.close(fig)
        bg = np.array([int(BG[i:i + 2], 16) for i in (1, 3, 5)])
        rows = np.where((np.abs(img - bg).sum(axis=2) > 24).any(axis=1))[0]
        if rows.size == 0:
            _offsets[key] = 0.0
        else:
            h = img.shape[0]
            top, bottom = 1 - rows[0] / h, 1 - (rows[-1] + 1) / h     # fractions of the frame, from below
            dy = (BAND_BOTTOM - TOP_GAP) - top                          # top just under the band
            _offsets[key] = float(max(dy, 0.01 - bottom)) if top - bottom < BAND_BOTTOM - TOP_GAP else 0.0
    return _offsets[key]


def draw_frame(job):
    global _plt
    i, path = job
    if _plt is None:
        _plt = setup_matplotlib()
    import scenes
    for s, f0, n in TL:
        if f0 <= i < f0 + n:
            break
    j = i - f0
    u_scene = min(1.0, j / max(1, (n - 1) * (1 - HOLD)))
    u = j / max(1, n - 1)
    fig = _plt.figure(figsize=(16, 9), dpi=WIDTH / 16, facecolor=BG)
    kind = s.get('kind')
    if kind == 'card':
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 16)
        ax.set_ylim(0, 9)
        ax.axis('off')
        card(fig, ax, s, u, n)
    else:
        if kind == 'code':
            final, now = (lambda f, a: code_moment(f, a, s, 1.0, n)), (lambda f, a: code_moment(f, a, s, u, n))
        elif s['id'] == 't_flexible':
            final, now = (lambda f, a: flexible(f, a, s, 1.0, n)), (lambda f, a: flexible(f, a, s, u, n))
        else:
            draw = scenes.all_scenes()[s['scene']][0]
            final, now = (lambda f, a: draw(f, a, 1.0)), (lambda f, a: draw(f, a, u_scene))
        in_content(fig, now, content_offset(s['id'], final))
        ax = fig.add_axes([0, 0, 1, 1])                    # the title band, over everything
        ax.set_xlim(0, 16)
        ax.set_ylim(0, 9)
        ax.axis('off')
        if kind == 'pillar':
            pillar(fig, ax, s, u, n)
        else:
            caption(ax, s['words'], u, n, s['act'])
    # fade in from black at the start, out to black at the end
    total = TL[-1][1] + TL[-1][2]
    fade = max(1 - i / 12, (i - (total - 45)) / 45, 0)
    if fade > 0:
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((0, 0), 16, 9, fc='#000000', ec='none', alpha=min(1, fade), zorder=50))
    fig.savefig(path, facecolor=BG)
    _plt.close(fig)
    return i


TL = timeline()

if __name__ == '__main__':
    args = sys.argv[1:]
    jobs = int(args[args.index('--jobs') + 1]) if '--jobs' in args else os.cpu_count()
    if '--width' in args:
        WIDTH = int(args[args.index('--width') + 1])
    out = BUILD / 'cut'
    out.mkdir(parents=True, exist_ok=True)
    total = TL[-1][1] + TL[-1][2]
    todo = [(i, str(out / ('f%05d.png' % i))) for i in range(total) if not (out / ('f%05d.png' % i)).exists()]
    print('%d frames (%.2f s), %d to draw on %d processes' % (total, total / FPS, len(todo), jobs), flush=True)
    with mp.get_context('fork').Pool(jobs) as pool:
        for k, _ in enumerate(pool.imap_unordered(draw_frame, todo, chunksize=8), 1):
            if k % 250 == 0:
                print(k, flush=True)
    if '--frames-only' in args:
        sys.exit()
    import music
    wav = music.write(total / FPS)
    import render
    exe = render.ffmpeg()
    dest = BUILD / 'magnus_trailer.mp4'
    subprocess.run([exe, '-y', '-loglevel', 'error', '-framerate', str(FPS), '-i', str(out / 'f%05d.png'),
                    '-i', str(wav), '-c:v', 'libx264', '-crf', '18', '-preset', 'slow', '-pix_fmt', 'yuv420p',
                    '-c:a', 'aac', '-b:a', '192k', '-movflags', '+faststart', '-shortest', str(dest)], check=True)
    print('%s: %.1f MB' % (dest.relative_to(REPO), dest.stat().st_size / 1e6))
