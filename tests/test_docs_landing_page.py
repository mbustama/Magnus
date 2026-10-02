r"""The documentation's landing page does not repeat the sidebar or link the changelog.

Every toctree on ``docs/source/index.rst`` is ``:hidden:``: Sphinx then builds the sidebar
navigation from it, without also printing the full list of sections in the page body, which
repeated the sidebar.  The changelog stays in that navigation but is not linked from the
landing page's own text.
"""
import pathlib
import re

import pytest

INDEX = pathlib.Path(__file__).resolve().parents[1]/'docs'/'source'/'index.rst'


def _toctrees(text):
    r"""The option lines of each toctree directive."""
    blocks = re.split(r'^\.\. toctree::\n', text, flags=re.M)[1:]
    return [re.findall(r'^   (:[a-z_]+:.*)$', b.split('\n\n')[0], flags=re.M) for b in blocks]


@pytest.mark.checkout_only     # reads docs/
def test_every_toctree_is_hidden():
    toctrees = _toctrees(INDEX.read_text())
    assert toctrees, 'no toctree found on the landing page'
    assert all(':hidden:' in opts for opts in toctrees)


@pytest.mark.checkout_only     # reads docs/
def test_landing_page_does_not_link_the_changelog():
    body = re.sub(r'^\.\. toctree::\n(?:   .*\n|\n)*', '', INDEX.read_text(), flags=re.M)
    assert ':doc:`changelog`' not in body
