# -*- coding: utf-8 -*-
r"""Keep the validation table in the docs and the figure's data in one story.

``docs/source/adiabatic_strategy.rst`` prints a table of measured speed-ups, and
``docs/make_figures.py`` draws the same numbers as a bar chart from
``VALIDATION_GRID``.  The figure caption tells the reader the two cannot drift apart.
Nothing enforced that until this test: the script writes an SVG and never touches the
page, so the two were a hand-maintained pair with a comment asking someone to remember.

The check is on the data, not on the prose.  The page's row labels carry detail the
grid does not ("18 MeV, 0.9 R_sun-scale baseline") and the two order their rows
differently, so what is compared is the multiset of (windows, speed-up) pairs -- which
is exactly what the caption claims the two share.

``VALIDATION_GRID`` is read with ``ast`` rather than imported, so this costs no
Matplotlib import and runs no drawing code.
"""

import ast
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RST = REPO / 'docs' / 'source' / 'adiabatic_strategy.rst'
FIGURES = REPO / 'docs' / 'make_figures.py'


def _validation_grid():
    """``VALIDATION_GRID`` from docs/make_figures.py, without importing it."""
    tree = ast.parse(FIGURES.read_text(encoding='utf-8'), filename=str(FIGURES))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == 'VALIDATION_GRID'
                for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError("no VALIDATION_GRID assignment in %s" % FIGURES)


def _documented_rows():
    """(windows, speed-up) for every row of the validation table in the page."""
    lines = RST.read_text(encoding='utf-8').splitlines()
    try:
        start = lines.index('   * - Case')
    except ValueError:
        raise AssertionError("the validation table's header row moved in %s" % RST)

    rows, cells = [], None
    for line in lines[start:]:
        if line.startswith('   * - '):
            if cells is not None:
                rows.append(cells)
            cells = [line[len('   * - '):].strip()]
        elif line.startswith('     - ') and cells is not None:
            cells.append(line[len('     - '):].strip())
        elif not line.strip():
            continue
        else:
            break
    if cells is not None:
        rows.append(cells)

    out = []
    for cells in rows[1:]:                       # rows[0] is the header
        assert len(cells) == 5, "unexpected column count in row %r" % (cells,)
        windows = int(cells[1])
        speedup = int(re.fullmatch(r'~?([\d,]+)x', cells[2]).group(1).replace(',', ''))
        out.append((windows, speedup))
    return out


def test_validation_table_matches_the_figures_data():
    """The page's table and the figure's VALIDATION_GRID carry the same numbers."""
    documented = sorted(_documented_rows())
    drawn = sorted((windows, speedup) for _, windows, speedup in _validation_grid())
    assert documented == drawn, (
        "docs/source/adiabatic_strategy.rst and docs/make_figures.py have drifted.\n"
        "  page:          %s\n"
        "  VALIDATION_GRID: %s\n"
        "Update both, then rerun docs/make_figures.py to redraw the chart."
        % (documented, drawn))
