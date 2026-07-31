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
BLUE = '#1c71d8'    # the oscillating probability
ORANGE = '#e66100'  # its phase-averaged value
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


def averaging_regimes(path):
    """The three regimes a pair of eigenvalues can be in, and where averaging applies.

    Data-driven rather than hand-drawn: the curve is a real two-flavour vacuum
    probability, so the figure cannot drift away from the formula it illustrates.

    The phase axis is logarithmic because the three regimes span decades -- on a linear
    axis the coherent one (below ~1e-2 rad) is an invisible sliver against the many
    cycles of the decohered one.
    """
    import numpy as np

    # P = 1 - sin^2(2 theta) sin^2(Delta), with Delta the accumulated relative phase.
    sin2_2theta = 0.85
    phase = np.logspace(-3.0, np.log10(40.0), 6000)
    prob = 1.0 - sin2_2theta*np.sin(phase)**2
    averaged = 1.0 - 0.5*sin2_2theta

    coherent_edge, decohered_edge = 1.0e-2, 2.0*np.pi

    fig, ax = plt.subplots(figsize=(7.8, 3.6))

    ax.axvspan(1.0e-3, coherent_edge, color=GREEN, alpha=0.13)
    ax.axvspan(coherent_edge, decohered_edge, color=RED, alpha=0.11)
    ax.axvspan(decohered_edge, 40.0, color=GRID, alpha=0.45)

    ax.plot(phase, prob, color=BLUE, lw=1.0, label='Oscillation probability')
    ax.axhline(averaged, color=ORANGE, lw=1.6, ls='--', label='Phase-averaged value')

    ax.axvline(coherent_edge, color=INK, lw=0.8, ls=':')
    ax.axvline(decohered_edge, color=INK, lw=0.8, ls=':')

    # Annotations sit below the oscillation's minimum (1 - sin^2 2theta = 0.15), so they
    # never overlap the curve however many cycles are drawn.
    for x, text in [(3.0e-3, 'coherent:\nnothing to average'),
                    (0.22, 'neither limit:\nno averaged expression\ndescribes this'),
                    (13.0, 'decohered:\nthe average is exact')]:
        ax.annotate(text, (x, 0.015), fontsize=8.5, color=INK, ha='center', va='bottom')

    ax.set_xlabel(r'Accumulated relative phase, $(\lambda_i-\lambda_j)\,L$  [rad]', fontsize=9)
    ax.set_ylabel('Probability', fontsize=9)
    ax.set_xscale('log')
    ax.set_xlim(1.0e-3, 40.0)
    ax.set_ylim(0.0, 1.05)
    ax.tick_params(colors=INK, labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(INK)

    # Above the axes: inside, it sits on top of the oscillation at large phase.
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), frameon=False, fontsize=9, ncol=2)

    fig.subplots_adjust(left=0.09, right=0.98, top=0.87, bottom=0.19)
    fig.savefig(path, metadata=SVG_METADATA)
    plt.close(fig)
    print(f'wrote {path}')


def main():
    STATIC.mkdir(parents=True, exist_ok=True)
    speedup_chart(STATIC / 'adiabatic_speedup.svg')
    averaging_regimes(STATIC / 'averaging_regimes.svg')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
