# -*- coding: utf-8 -*-
r"""The committed notebooks are the ones `make_notebooks.py` builds.

Continuous integration executes the ``.ipynb`` files in ``notebooks/``, not the generator
that writes them.  So an edit to ``make_notebooks.py`` changes nothing that any gate runs
until the notebooks are rebuilt and committed, and until then the repository holds a
generator and an artifact that disagree, with every check reading the artifact.

That is not hypothetical.  A missing ``import matplotlib as mpl`` was fixed in the
generator, verified locally by executing ``make_notebooks.books[...]`` -- the freshly built
in-memory notebook, which had the fix -- and then failed on the runner forty-six minutes
later against the committed file, which did not.  The same commit's new Figure 12 cell was
absent from the artifact for the same reason.

`make_notebooks.py` already refuses to finish while any notebook on disk is stale.  That
check only runs when somebody rebuilds, which is exactly what an author who forgot to
rebuild has not done.  Here it runs on every push, and it executes nothing: only the cell
sources are compared.
"""

import pathlib
import sys

import nbformat

ROOT = pathlib.Path(__file__).resolve().parent.parent
NOTEBOOKS = ROOT/'notebooks'
sys.path.insert(0, str(NOTEBOOKS))


_BUILT = {}


def built():
    """The notebooks `make_notebooks.py` writes, without executing any of them.

    `add_footers` comes first because `build` calls it first: it appends a cell to every
    notebook, so `books` on its own is one cell short of what reaches disk.  Comparing
    against the unfootered dictionary reports all twenty-nine as stale, which is a test
    that cannot be acted on.
    """
    if not _BUILT:
        import matplotlib
        matplotlib.use('Agg')                   # importing the module touches pyplot
        import make_notebooks
        # Once, and memoised: `add_footers` appends rather than replaces, so calling it
        # a second time leaves every notebook a cell longer than the one on disk and the
        # comparison below fails for a reason that has nothing to do with the repository.
        make_notebooks.add_footers()
        _BUILT.update(make_notebooks.books)
    return _BUILT


def test_the_generator_describes_some_notebooks():
    """So that a generator which stopped building anything cannot pass this file."""
    assert built(), 'make_notebooks.py builds no notebooks'


def test_every_committed_notebook_matches_the_generator():
    """What is on disk is what the generator writes, cell for cell."""
    stale, missing = [], []
    for name, nb in built().items():
        path = NOTEBOOKS/name
        if not path.exists():
            missing.append(name)
            continue
        on_disk = nbformat.read(path, as_version=4)
        if [c.source for c in on_disk.cells] != [c.source for c in nb.cells]:
            stale.append(name)
    assert not missing, ('the generator builds notebooks that are not committed: %s'
                         % ', '.join(missing))
    assert not stale, (
        'these notebooks on disk are not what make_notebooks.py builds: %s.  Continuous '
        'integration executes the committed files, so until they are rebuilt no gate sees '
        'the change.  Run:  MAGNUS_PAPER_CACHE_ONLY=1 python notebooks/make_notebooks.py'
        % ', '.join(stale))
