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


# ---------------------------------------------------------------------------
# Schematics for docs/source/architecture.rst
#
# These replace three Mermaid blocks.  Mermaid sizes a node by the width of its
# label, and those labels carried whole function signatures, so the boxes grew
# wider than the canvas and the text collided.  Here the geometry is explicit:
# every box is placed on a fixed grid and every arrow is clipped to the box
# borders, so nothing can overlap however the labels change.
#
# The rule the old diagrams broke, and these keep: **the figure carries the
# shape, the surrounding table carries the words.**  Boxes hold module or
# function names and nothing else.
# ---------------------------------------------------------------------------

PANEL = '#eef3f9'    # module box fill
PANEL_EDGE = '#5a6b7d'
BAND = '#dce7f5'     # the wide band for oscprob / layer groups

# The internal import graph, top-level imports only, as read from the source.
# Re-derive with:
#     python3 -c "import ast,pathlib;[print(p.stem, sorted({(n.module or '').split('.')[1]
#       for n in ast.walk(ast.parse(p.read_text()))
#       if isinstance(n,(ast.Import,ast.ImportFrom))}))
#       for p in pathlib.Path('src/magnus').glob('*.py')]"
# `main()` checks it against the real package and fails rather than drawing a
# picture that has quietly stopped being true -- which is how the diagram this
# replaces came to claim that `hamiltonians` does not import `matter`.
MODULE_TIERS = [
    ['magnus', 'expmkernels', 'expansionterms'],
    ['globaldefs', 'adiabatic'],
    ['earth', 'matter', 'avgprob'],
    ['hamiltonians'],
    ['oscprobstd'],
]
MODULE_EDGES = [
    ('magnus', 'globaldefs'), ('magnus', 'adiabatic'),
    ('globaldefs', 'earth'), ('globaldefs', 'matter'),
    ('adiabatic', 'avgprob'),
    ('matter', 'hamiltonians'), ('globaldefs', 'hamiltonians'),
    ('hamiltonians', 'oscprobstd'),
]


def _rounded(ax, x, y, w, h, label, face=PANEL, edge=PANEL_EDGE, fontsize=9.5,
             weight='normal', color=INK, zorder=2):
    """A centred rounded box.  Returns its rectangle for arrow clipping."""
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch(
        (x - w/2.0, y - h/2.0), w, h,
        boxstyle='round,pad=0.0,rounding_size=0.10',
        linewidth=1.1, facecolor=face, edgecolor=edge, zorder=zorder))
    ax.text(x, y, label, ha='center', va='center', fontsize=fontsize,
            color=color, fontweight=weight, zorder=zorder + 1)
    return (x, y, w, h)


def _border(rect, toward):
    """Where the segment from a box's centre to `toward` leaves the box.

    Clipping every arrow to the border is what keeps arrowheads off the text,
    whatever the boxes are labelled.
    """
    x, y, w, h = rect
    dx, dy = toward[0] - x, toward[1] - y
    if dx == 0 and dy == 0:
        return (x, y)
    scale = min((w/2.0)/abs(dx) if dx else float('inf'),
                (h/2.0)/abs(dy) if dy else float('inf'))
    return (x + dx*scale, y + dy*scale)


def _connect(ax, a, b, color=PANEL_EDGE, dashed=False, lw=1.1, zorder=1):
    from matplotlib.patches import FancyArrowPatch
    start = _border(a, (b[0], b[1]))
    end = _border(b, (a[0], a[1]))
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle='-|>', mutation_scale=11, linewidth=lw,
        color=color, linestyle=(0, (4, 2)) if dashed else 'solid',
        shrinkA=0, shrinkB=0, zorder=zorder))


def module_layout(path):
    """The internal import graph, by dependency tier.

    Every arrow points from a module to one that imports it, so the picture
    reads bottom-up: nothing at the bottom knows anything above it.
    """
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.45, 6.45)
    ax.axis('off')

    # The dependency graph lives in x < 7.5; `plotting` gets its own lane to the
    # right, because it belongs to no tier -- it depends on one module and
    # nothing depends on it.
    left, width = 0.35, 7.05
    bw, bh, dy = 1.62, 0.52, 1.05
    pos, rects = {}, {}
    for tier, names in enumerate(MODULE_TIERS):
        y = 0.55 + tier*dy
        span = width/(len(names) + 1)
        for i, name in enumerate(names):
            x = left + span*(i + 1)
            pos[name] = (x, y)
            rects[name] = _rounded(ax, x, y, bw, bh, name)

    for src, dst in MODULE_EDGES:
        deferred = (src, dst) == ('globaldefs', 'hamiltonians')
        _connect(ax, rects[src], rects[dst],
                 color=RED if deferred else PANEL_EDGE, dashed=deferred)

    # Top-left is the only region no box or arrow occupies; the label is the
    # sole red element on the figure, so it needs no leader line to the edge.
    ax.text(1.35, 4.28, 'deferred import —\nwithout it, this graph\nhas a cycle',
            fontsize=7.4, color=RED, ha='center', va='center', linespacing=1.4)

    # oscprob spans the top rather than receiving eight crossing arrows: that it
    # imports everything below IS the architectural point, and one band says it
    # more clearly than the spaghetti would.
    top = 0.55 + len(MODULE_TIERS)*dy + 0.12
    _rounded(ax, left + width/2.0, top, width + 0.30, 0.58,
             'magnus.oscprob  —  imports every module below',
             face=BAND, edge=PANEL_EDGE, fontsize=9.4, weight='bold')
    ax.plot([left, left + width], [top - 0.45, top - 0.45], color=PANEL_EDGE,
            lw=0.9, ls=(0, (3, 3)), zorder=0)

    ax.plot([7.85, 7.85], [0.15, top + 0.30], color='#c8d0d8', lw=0.9,
            ls=(0, (3, 3)), zorder=0)
    # No arrow across the divider: it would have to cross the whole graph to
    # reach globaldefs, and the caption says the same thing without the line.
    _rounded(ax, 9.00, 1.60, 1.85, 0.52, 'plotting',
             face='#f7f7f7', edge='#9aa5b1', color='#555555')
    ax.text(9.00, 0.92, 'imports globaldefs\nonly; nothing\nimports it',
            ha='center', va='center', fontsize=7.4, color='#666666',
            linespacing=1.4)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(path, metadata=SVG_METADATA)
    plt.close(fig)
    print(f'wrote {path}')


def api_layers(path):
    """The three layers of ``magnus.oscprob``, names only.

    The Mermaid version put full signatures in the boxes; that is what made it
    overlap.  The signatures belong in the table beside it.
    """
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.6)
    ax.axis('off')

    from matplotlib.patches import FancyBboxPatch

    # Rows per band: four long names do not fit across one row at a legible
    # size, and shrinking the font to make them fit is how a figure becomes
    # unreadable at print size.  They wrap instead.
    layers = [
        ('Layer 3 — wrappers', '~60 functions',
         [['osc_prob_3nu_earth', 'osc_prob_2nu_matter_nsi…', 'osc_prob_5nu_vacuum_liv']]),
        ('Layer 2 — scenario', '4 functions, generic in num_flavors',
         [['osc_prob_vacuum', 'osc_prob_matter_std_potential'],
          ['osc_prob_matter_nsi', 'osc_prob_liv']]),
        ('Layer 1 — primordial', '2 functions',
         [['osc_prob_energy_baseline', 'osc_prob']]),
    ]

    bw_band, x0 = 7.5, 0.35
    row_h, head_h, gap = 0.50, 0.46, 0.34
    y_cursor = 5.35
    band_rects = []
    for title, count, rows in layers:
        height = head_h + row_h*len(rows) + 0.14
        y_top = y_cursor
        y_bot = y_top - height
        ax.add_patch(FancyBboxPatch(
            (x0, y_bot), bw_band, height,
            boxstyle='round,pad=0.0,rounding_size=0.10',
            linewidth=1.1, facecolor=BAND, edgecolor=PANEL_EDGE, zorder=1))
        ax.text(x0 + 0.27, y_top - 0.25, title, fontsize=9.6, color=INK,
                fontweight='bold', va='center', zorder=3)
        ax.text(x0 + bw_band - 0.27, y_top - 0.25, count, fontsize=8.0,
                color='#5a6b7d', va='center', ha='right', zorder=3)
        for r, names in enumerate(rows):
            span = bw_band/len(names)
            yr = y_top - head_h - row_h*(r + 0.5)
            for i, name in enumerate(names):
                _rounded(ax, x0 + span*(i + 0.5), yr, span - 0.26, 0.42,
                         name, face='#ffffff', fontsize=8.0, zorder=2)
        band_rects.append((x0 + bw_band/2.0, (y_top + y_bot)/2.0, bw_band, height))
        y_cursor = y_bot - gap

    for upper, lower in zip(band_rects, band_rects[1:]):
        _connect(ax, upper, lower, lw=1.3)

    core = _rounded(ax, x0 + bw_band/2.0, y_cursor - 0.28, 5.6, 0.50,
                    'magnus.magnus — the Magnus expansion core',
                    face=PANEL, fontsize=8.8, weight='bold')
    _connect(ax, band_rects[-1], core, lw=1.3)

    # The escape hatch keeps to the right of every band, so "bypasses" is what
    # the geometry shows as well as what the caption says.
    _rounded(ax, 9.05, band_rects[0][1], 1.62, 0.80,
             'your own\nH_func(l)', face='#fff4e8', edge=ORANGE,
             fontsize=8.2, color='#8a4b00')
    ax.annotate('', xy=(x0 + bw_band, band_rects[-1][1]),
                xytext=(9.05, band_rects[0][1] - 0.44),
                arrowprops=dict(arrowstyle='-|>', mutation_scale=11,
                                color=ORANGE, linewidth=1.2,
                                linestyle=(0, (4, 2)), shrinkA=2, shrinkB=0,
                                connectionstyle='angle,angleA=-90,angleB=0,rad=6'))
    # Right-aligned so it ends before the dashed lane rather than sitting on it.
    ax.text(8.88, (band_rects[0][1] + band_rects[-1][1])/2.0 + 0.12,
            'bypasses layers\n3 and 2 entirely', ha='right', va='center',
            fontsize=7.4, color=ORANGE, linespacing=1.4)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(path, metadata=SVG_METADATA)
    plt.close(fig)
    print(f'wrote {path}')


def call_sequence(path):
    """What is built, in what order, on one call.

    The point this figure exists to make is that the potential and the
    Hamiltonian are built **once per call**, not once per slab, and stay plain
    callables all the way down -- which is what lets `probe_eval_mode` decide
    vectorization once.  A static dependency graph cannot show that; an ordered
    one can.
    """
    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.4)
    ax.axis('off')

    actors = ['wrapper', 'scenario', 'matter', 'hamiltonians', 'baseline', 'osc_prob', 'core']
    # Spacing derived from the count so the last lifeline stays inside the
    # canvas; hard-coding a stride put `core` off the right edge.
    margin, span = 0.72, 8.56
    xs = [margin + i*span/(len(actors) - 1) for i in range(len(actors))]
    top, bottom = 4.95, 0.45

    for x, name in zip(xs, actors):
        _rounded(ax, x, top, 1.30, 0.44, name, face=PANEL, fontsize=8.2)
        ax.plot([x, x], [top - 0.24, bottom], color='#9aa5b1', lw=0.8,
                ls=(0, (2, 3)), zorder=0)

    steps = [
        (0, 1, 4.35, 'osc_params, rho_func'),
        (1, 2, 3.90, 'build VCC_func(l)'),
        (1, 3, 3.45, 'wrap into H(l)'),
        (1, 4, 3.00, 'H_func, energy[], L[]'),
        (4, 5, 2.20, 'one (energy, L) point'),
        (5, 6, 1.75, 'refine slabs'),
        (6, 5, 1.30, 'unitary U'),
        (5, 4, 0.85, 'P = |U|²'),
    ]
    for a, b, y, label in steps:
        back = xs[b] < xs[a]
        ax.annotate('', xy=(xs[b], y), xytext=(xs[a], y),
                    arrowprops=dict(arrowstyle='-|>', mutation_scale=10,
                                    color=ORANGE if back else PANEL_EDGE,
                                    linewidth=1.0,
                                    linestyle=(0, (4, 2)) if back else 'solid',
                                    shrinkA=0, shrinkB=0))
        ax.text((xs[a] + xs[b])/2.0, y + 0.10, label, ha='center', va='bottom',
                fontsize=7.4, color=ORANGE if back else INK)

    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch(
        (xs[4] - 0.62, 0.60), xs[6] - xs[4] + 1.24, 1.92,
        boxstyle='round,pad=0.0,rounding_size=0.06', linewidth=0.9,
        facecolor='none', edgecolor='#9aa5b1', linestyle=(0, (5, 3)), zorder=0))
    ax.text(xs[4] - 0.52, 2.63, 'for every (energy, L) point — in parallel if n_jobs ≠ 1',
            fontsize=7.6, color='#5a6b7d', va='bottom')

    ax.text(0.15, 0.12,
            'Everything passed downward is a plain callable, built once per call '
            'rather than once per slab.',
            fontsize=7.8, color=INK, va='bottom')

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(path, metadata=SVG_METADATA)
    plt.close(fig)
    print(f'wrote {path}')


def _check_module_graph():
    """Fail rather than draw a module diagram that has stopped being true."""
    src = Path(__file__).resolve().parent.parent / 'src' / 'magnus'
    drawn = {m for tier in MODULE_TIERS for m in tier} | {'oscprob', 'plotting'}
    actual = {p.stem for p in src.glob('*.py')
              if p.stem not in ('__init__', '__main__', 'version', 'authors', 'cli')}
    actual.add('hamiltonians')
    missing, extra = actual - drawn, drawn - actual
    if missing or extra:
        raise SystemExit(
            'module_layout is out of step with src/magnus: missing %s, stale %s'
            % (sorted(missing) or 'none', sorted(extra) or 'none'))


def main():
    STATIC.mkdir(parents=True, exist_ok=True)
    speedup_chart(STATIC / 'adiabatic_speedup.svg')
    averaging_regimes(STATIC / 'averaging_regimes.svg')
    _check_module_graph()
    module_layout(STATIC / 'module_layout.svg')
    api_layers(STATIC / 'api_layers.svg')
    call_sequence(STATIC / 'call_sequence.svg')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
