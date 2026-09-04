# -*- coding: utf-8 -*-
r"""Drag Figure 11's rtol labels, and print the offsets to paste back.

Run it, drag any label with the mouse, close the window.  It prints the
``RTOL_LABEL_OFFSETS`` dict with whatever you left them at, in exactly the form
``notebooks/make_notebooks.py`` holds it, and writes the same text to
``fig11_offsets.txt`` beside this file.

    python3 docs/dev/nudge_fig11_labels.py

It does NOT duplicate the figure.  It lifts the figure cell out of the generator and
executes it, so what you drag is what the paper prints, and a change to the generator
cannot leave this tool drawing something else.  Only the offsets are yours to move.

WHY BY HAND AT ALL.  On the exponential profile the order-4 curve is nearly vertical --
five points spanning eight decades in error and barely one in time -- so consecutive
markers are a few points apart on screen and no placement rule keeps their labels from
landing on each other.  Alternating sides was tried and merely moved the collision.

Nothing is written back automatically: the printed dict is meant to be read before it is
pasted, because two labels can be dragged into agreement on screen and still collide at
the paper's own scale.
"""

import os
import pathlib
import re
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
GEN = ROOT/'notebooks'/'make_notebooks.py'
NB = ROOT/'notebooks'
OUT = pathlib.Path(__file__).resolve().parent/'fig11_offsets.txt'

START = "# ------------------------------------------------------------ smooth profile"
END = "save(fig, 'smooth_reach.pdf')"


def cell_source():
    text = GEN.read_text()
    i = text.index(START)
    j = text.index(END, i)
    return text[i:j]


def main():
    sys.path.insert(0, str(ROOT/'src'))
    import matplotlib
    interactive = bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'))
    if not interactive:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.text import Annotation
    from matplotlib.ticker import FuncFormatter, LogLocator

    from magnus import earth, globaldefs as gd, matter

    gen = GEN.read_text()

    def grab(name):
        m = re.search(r'^def %s\(.*?(?=\n(?:def |[A-Za-z_]+ = |#))' % name, gen,
                      re.S | re.M)
        return m.group(0)

    def const(name):
        m = re.search(r'^%s = .*$' % name, gen, re.M)
        return m.group(0)

    ns = {'np': np, 'plt': plt, 'json': __import__('json'), 'os': os, 'pe': pe,
          'HERE': NB, 'earth': earth, 'gd': gd, 'matter': matter,
          'COL': 3.487, 'WIDE': 7.224,
          'INK': '#333333', 'GRID': '#cccccc', 'RED': '#c01c28',
          'BLUE': '#1c71d8', 'GREEN': '#2ec27e',
          'LogLocator': LogLocator, 'FuncFormatter': FuncFormatter,
          'save': lambda *a, **k: None}
    for c in ('FLAVOR_LABEL',):
        exec(const(c), ns)
    for fn in ('logx', 'logy', 'corner', 'stamp', '_plain', 'snug'):
        try:
            exec(grab(fn), ns)
        except AttributeError:
            pass
    exec(cell_source(), ns)

    # Both columns carry labels, on their own first flavour row, with their own
    # offsets.  Dragging one must not write the other's, so each panel is collected
    # under its own key and printed as its own block.
    panels = {}
    for key, bench_name in (('exp', 'BENCH'), ('earth', 'PREM_BENCH')):
        bench = ns.get(bench_name)
        if bench is None:
            continue
        ax = ns['cols'][0 if key == 'exp' else 1][0]
        anns = [t for t in ax.texts if isinstance(t, Annotation) and t.get_text()]
        if not anns:
            continue
        keys = [p['label'] for p in
                {s['name']: s for s in bench['cases'][0]['series']}['Magnus']['points']]
        if len(keys) != len(anns):
            print('%s: mismatch, %d rtol values against %d annotations drawn'
                  % (key, len(keys), len(anns)))
        panels[key] = (keys, anns)
    anns = [a for _k, aa in panels.values() for a in aa]
    if not anns:
        raise SystemExit('no annotations drawn in either column')
    for a in anns:
        a.draggable(True)

    if interactive:
        # Bigger on screen, same figure in inches.  Resizing the figure would stretch the
        # curve without stretching the labels, so offsets set here would not be the ones
        # that print; raising the dpi changes only how many screen pixels a point is
        # drawn across, which is exactly the precision that was missing.
        ns['fig'].set_dpi(3.0*ns['fig'].get_dpi())

        # Dragging still quantises to whole screen pixels.  Arrow keys do not: click a
        # label to select it, then nudge in fifths of a point, or whole points with shift.
        state = {'sel': None}

        def on_pick(ev):
            if isinstance(ev.artist, Annotation):
                state['sel'] = ev.artist
                ns['fig'].canvas.manager.set_window_title(
                    'selected: %s   (arrows nudge 0.2 pt, shift+arrows 1 pt)'
                    % ev.artist.get_text())

        def on_key(ev):
            a = state['sel']
            if a is None or ev.key is None:
                return
            step = 1.0 if 'shift' in ev.key else 0.2
            k = ev.key.split('+')[-1]
            dx, dy = a.xyann
            if k == 'left':
                dx -= step
            elif k == 'right':
                dx += step
            elif k == 'up':
                dy += step
            elif k == 'down':
                dy -= step
            else:
                return
            a.xyann = (dx, dy)
            print('  %-22s -> (%.1f, %.1f)' % (a.get_text(), dx, dy), flush=True)
            ns['fig'].canvas.draw_idle()

        for a in anns:
            a.set_picker(True)
        ns['fig'].canvas.mpl_connect('pick_event', on_pick)
        ns['fig'].canvas.mpl_connect('key_press_event', on_key)

        print('Drag a label, or click it and use the arrow keys: 0.2 pt a press,')
        print('1 pt with shift.  Zoom with the toolbar to aim.  Do NOT resize the')
        print('window.  Close it when done.')
        plt.show()

    lines = ['RTOL_LABEL_OFFSETS = {']
    for col_key, (keys, aa) in panels.items():
        lines.append('    %r: {' % col_key)
        for key, a in zip(keys, aa):
            dx, dy = a.xyann
            # The alignment the annotation actually carries, not one inferred from the
            # sign of the offset.  A drag moves xyann and leaves ha/va alone, so
            # recomputing them shifts the text by its own width or height on top of the
            # offset -- which is exactly how the first round of drags came back wrong.
            lines.append("        %r: (%.1f, %.1f, %r, %r),"
                         % (key, dx, dy, a.get_ha(), a.get_va()))
        lines.append('    },')
    lines.append('}')
    if OUT.exists():                       # never clobber the previous round of drags
        OUT.with_suffix('.txt.bak').write_text(OUT.read_text())
    OUT.write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print('\nwritten to %s' % OUT)


if __name__ == '__main__':
    main()
