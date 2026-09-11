# -*- coding: utf-8 -*-
r"""Every figure the paper includes is a file the repository actually carries.

`.gitignore` carries a blanket ``*.pdf``, and the paper's figures are tracked only because
each was force-added.  So ``git add -A`` does not pick up a new one, says nothing, and the
omission surfaces as a LaTeX error for whoever next builds from a clean checkout -- not for
whoever introduced it, whose working tree still has the file.

That happened: ``shock_cost.pdf`` was committed in neither the commit that added the figure
nor the one that rewrote the section around it, while ``main.tex`` had already been changed
to include it.  Nothing in the suite noticed, because nothing was reading `main.tex`.
"""

import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent
PAPER = ROOT/'resources'/'paper'/'main.tex'
FIGS = ROOT/'resources'/'paper'/'figs'


def included_figures():
    """Every file named by an \\includegraphics in the paper."""
    return sorted(set(re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',
                                 PAPER.read_text())))


def tracked():
    """The repository's tracked files under the paper's figure directory."""
    out = subprocess.run(['git', 'ls-files', str(FIGS.relative_to(ROOT))],
                         cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return {pathlib.Path(line).name for line in out.split()}


def test_the_paper_includes_some_figures():
    """The pattern above still matches, so the real test cannot pass by finding nothing."""
    assert included_figures(), (
        'no \\includegraphics found in %s; either the paper changed shape or this test no '
        'longer recognises its figures' % PAPER.relative_to(ROOT))


def test_every_included_figure_is_tracked():
    """A clean checkout can build the paper."""
    have = tracked()
    missing = [name for name in included_figures() if pathlib.Path(name).name not in have]
    assert not missing, (
        'main.tex includes %s, which git does not track. .gitignore carries a blanket '
        '*.pdf, so these need `git add -f`; without it the paper does not build from a '
        'clean checkout.' % ', '.join(missing))


def test_no_tracked_figure_has_been_orphaned():
    """And nothing is carried that the paper stopped including.

    Not a correctness problem, but a figure nobody includes is a figure nobody regenerates,
    and it will eventually be read as current by someone.
    """
    included = {pathlib.Path(name).name for name in included_figures()}
    orphans = sorted(name for name in tracked()
                     if name.endswith('.pdf') and name not in included)
    assert not orphans, (
        'tracked but included nowhere in main.tex: %s' % ', '.join(orphans))
