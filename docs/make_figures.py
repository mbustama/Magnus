#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Regenerate the data-driven figure used by docs/source/adiabatic_strategy.rst.

``docs/source/_static/adiabatic_speedup.svg`` plots measured speedups, and was originally
produced by an ad-hoc script that was never committed.  With no way to re-derive it, it
drifted out of step with the validation table it sits directly beneath: the chart carried
a 25,800x bar for the case that table reports as ~30x (25,800x came from a *different*
measurement -- a standard-3nu run at 2 MeV that does not appear in the table at all).

The numbers now live in exactly one place in this repository, ``VALIDATION_GRID`` below,
which mirrors the table in ``adiabatic_strategy.rst``.  Change a measurement there, rerun
this script, and the chart follows::

    python3 docs/make_figures.py

The other two figures on that page (``adiabatic_avoided_crossing.svg`` and
``adiabatic_segmentation.svg``) are deliberately *not* generated here.  They are
hand-authored schematics that illustrate concepts rather than measurements, so they carry
no data that can go stale, and reproducing them from a script would mean redrawing them
less carefully than they were drawn.
"""

import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    sys.exit("matplotlib is required to regenerate the figures: pip install matplotlib")

# Make the SVG byte-reproducible: matplotlib otherwise stamps the current time into the
# file's metadata and derives glyph ids from a random salt, so regenerating an unchanged
# figure would still show up as a diff.
plt.rcParams['svg.hashsalt'] = 'magnus-docs'
SVG_METADATA = {'Date': None}

STATIC = Path(__file__).resolve().parent / 'source' / '_static'

# Palette matching the hand-authored figures on the same page.
GREEN = '#2ec27e'   # purely adiabatic cases (no Magnus patch needed)
RED = '#c01c28'     # cases needing one or more Magnus patches
INK = '#333333'
GRID = '#cccccc'

# (label, number of non-adiabatic windows, measured speedup vs. tight-tolerance solve_ivp).
# Must stay in step with the validation table in docs/source/adiabatic_strategy.rst.
VALIDATION_GRID = [
    ('Standard 3ν', 0, 3600),
    ('Standard 4ν (3+1)', 0, 4670),
    ('Standard 5ν (3+2)', 0, 4800),
    ('NSI 3ν, resonance', 1, 88),
    ('NSI 4ν, resonance', 1, 60),
    ('NSI 5ν, resonance', 1, 52),
    ('Synthetic, 2 nearby', 1, 91),
    ('Synthetic, 2 separate', 2, 30),
]


def speedup_chart(path):
    labels = [row[0] for row in VALIDATION_GRID]
    windows = [row[1] for row in VALIDATION_GRID]
    speedups = [row[2] for row in VALIDATION_GRID]
    colors = [GREEN if w == 0 else RED for w in windows]

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    bars = ax.bar(range(len(labels)), speedups, color=colors, edgecolor=INK, linewidth=0.6)
    ax.set_yscale('log')
    ax.set_ylim(10, 12000)
    ax.set_title('Measured speedup vs. a tight-tolerance solve_ivp ground truth',
                 color=INK, fontsize=11)
    ax.set_ylabel('Speedup (log scale)', color=INK)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
    ax.yaxis.grid(True, which='both', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(INK)
    ax.tick_params(colors=INK)

    for bar, value in zip(bars, speedups):
        ax.annotate(f'{value:,}x',
                    (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    textcoords='offset points', xytext=(0, 4),
                    ha='center', fontsize=9, color=INK)

    handles = [plt.Rectangle((0, 0), 1, 1, color=GREEN, ec=INK, lw=0.6),
               plt.Rectangle((0, 0), 1, 1, color=RED, ec=INK, lw=0.6)]
    ax.legend(handles, ['Purely adiabatic (0 windows)', 'One or more Magnus patches'],
              loc='upper right', frameon=False, fontsize=9)

    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.28)
    fig.savefig(path, metadata=SVG_METADATA)
    plt.close(fig)
    print(f'wrote {path}')


def main():
    STATIC.mkdir(parents=True, exist_ok=True)
    speedup_chart(STATIC / 'adiabatic_speedup.svg')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
