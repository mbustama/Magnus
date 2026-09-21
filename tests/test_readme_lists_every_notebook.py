# -*- coding: utf-8 -*-
r"""`notebooks/README.md` describes every notebook the generator builds.

The README is the folder's table of contents.  Nothing forces it to grow when a notebook
is added, so it drifts silently: it listed twelve notebooks while the generator built
twenty-nine, and every notebook added after the twelfth was invisible to a reader browsing
the folder.

A count alone would not catch a rename, so this compares the names themselves, in both
directions: a notebook the README omits, and a notebook the README names that no longer
exists.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
NOTEBOOKS = ROOT/'notebooks'
README = NOTEBOOKS/'README.md'

sys.path.insert(0, str(NOTEBOOKS))

_BUILT = {}

# A link to a notebook, as the README writes it: `(01_magnus_introduction.ipynb)`.
_LINK = re.compile(r'\((\d\d_magnus_[a-z0-9_]+\.ipynb)\)')


def built_names():
    """The notebook file names `make_notebooks.py` writes, without executing any."""
    if not _BUILT:
        import matplotlib
        matplotlib.use('Agg')               # importing the module touches pyplot
        import make_notebooks
        _BUILT.update(make_notebooks.books)
    return set(_BUILT)


def described_names():
    """The notebook file names the README links to."""
    return set(_LINK.findall(README.read_text(encoding='utf-8')))


def test_the_readme_exists_and_links_to_notebooks():
    """So that an empty or restructured README cannot pass the comparison below."""
    assert README.exists(), 'notebooks/README.md is missing'
    assert described_names(), (
        'notebooks/README.md links to no notebook.  This test finds them by the pattern '
        '(NN_magnus_name.ipynb); if the README now links to them some other way, update '
        '_LINK here rather than dropping the guard.')


def test_the_readme_describes_every_notebook():
    """Every notebook is in the README, and every notebook in the README exists."""
    built, described = built_names(), described_names()
    undescribed = sorted(built - described)
    stale = sorted(described - built)
    assert not undescribed, (
        'make_notebooks.py builds these notebooks, and notebooks/README.md does not '
        'describe them: %s.  A reader browsing the folder cannot find them.'
        % ', '.join(undescribed))
    assert not stale, (
        'notebooks/README.md describes these, and make_notebooks.py no longer builds '
        'them: %s.  Their links are dead.' % ', '.join(stale))
