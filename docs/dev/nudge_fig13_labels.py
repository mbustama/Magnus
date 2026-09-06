# -*- coding: utf-8 -*-
r"""Drag Figure 13's knob-rung labels, and print the offsets to paste back.

Run it, drag any label with the mouse, close the window.  It prints each
panel's annotation list with the offsets you left them at, in exactly the form
``notebooks/make_notebooks.py`` holds them, and writes the same text to
``fig13_offsets.txt`` beside this file.

    python3 docs/dev/nudge_fig13_labels.py

It does NOT duplicate the figure.  It lifts the figure cell out of the
generator and executes it, so what you drag is what the paper prints, and a
change to the generator cannot leave this tool drawing something else.  Only
the offsets are yours to move; everything else the cell decides.

Nothing is written back automatically: the printed lists are meant to be read
before they are pasted, because two labels can be dragged into agreement on
screen and still collide at the paper's own scale.
"""

import ast
import os
import pathlib
import re
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
GEN = ROOT/'notebooks'/'make_notebooks.py'
NB = ROOT/'notebooks'
OUT = pathlib.Path(__file__).resolve().parent/'fig13_offsets.txt'

START = "# ------------------------------------------------ three speed-accuracy planes"
END = "save(fig, 'speed_accuracy_combined.pdf')"


def cell_source():
    r"""The figure cell, lifted from the generator and stripped of its save()."""
    text = GEN.read_text()
    i = text.index(START)
    j = text.index(END, i)
    return text[i:j]


def item_lists(src):
    r"""The three annotate() item lists, in panel order, as Python objects."""
    tree = ast.parse(src)
    found = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and getattr(node.func, 'id', '') == 'annotate'):
            axes_arg = node.args[0]
            panel = axes_arg.slice.value if hasattr(axes_arg, 'slice') else None
            found.append((panel, ast.literal_eval(node.args[2])))
    return [items for _, items in sorted(found, key=lambda p: p[0])]


def main():
    sys.path.insert(0, str(ROOT/'src'))
    import matplotlib
    interactive = bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'))
    if not interactive:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.text import Annotation
    from matplotlib.ticker import FuncFormatter, LogLocator

    # The generator's own helpers, imported rather than re-implemented.
    gen = GEN.read_text()

    def grab(name):
        m = re.search(r'^def %s\(.*?(?=\n(?:def |[A-Za-z_]+ = |#))' % name, gen,
                      re.S | re.M)
        return m.group(0)

    import matplotlib.patheffects as pe
    ns = {'np': np, 'plt': plt, 'json': __import__('json'), 'os': os, 'pe': pe,
          'HERE': NB, 'WIDE': 7.224, 'INK': '#333333', 'GRID': '#cccccc',
          'LogLocator': LogLocator, 'FuncFormatter': FuncFormatter,
          'save': lambda *a, **k: None}
    for fn in ('logx', 'corner', 'stamp', '_plain'):
        exec(grab(fn), ns)
    exec(cell_source(), ns)

    fig = ns['fig']
    axes = ns['axes']
    anchors = [ns['a0'], ns['a1'], ns['a2']]
    lists = item_lists(cell_source())

    drawn = []
    for ax, anch, items in zip(axes, anchors, lists):
        kept = [it for it in items if tuple(it[0]) in anch]
        # Empty-text annotations are leader lines, not labels: the two that join
        # the shared "N_layers = 1" to both NuFast-Earth curves are drawn that way.
        anns = [t for t in ax.texts if isinstance(t, Annotation) and t.get_text()]
        if len(kept) != len(anns):
            print('panel mismatch: %d items kept, %d annotations drawn'
                  % (len(kept), len(anns)))
        for it, a in zip(kept, anns):
            a.draggable(True)
        drawn.append((items, kept, anns))

    if interactive:
        fig.set_size_inches(13.0, 13.0)      # room to aim; offsets are in points
        print('Drag any label.  Close the window when done.')
        plt.show()

    out = []
    for k, (items, kept, anns) in enumerate(drawn):
        moved = {id(it): tuple(a.xyann) for it, a in zip(kept, anns)}
        out.append('# --- panel %d ---' % k)
        for it in items:
            dx, dy = moved.get(id(it), (it[2], it[3]))
            key, text, ha = it[0], it[1], it[4]
            tail = ', True' if len(it) > 5 and it[5] else ''
            out.append("    ((%r, %r), %r, %.1f, %.1f, %r%s),"
                       % (key[0], key[1], text, dx, dy, ha, tail))
    OUT.write_text('\n'.join(out) + '\n')
    print('\n'.join(out))
    print('\nwritten to %s' % OUT)


if __name__ == '__main__':
    main()
