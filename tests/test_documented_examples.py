# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
"""The code in the README and the quickstart has to run.

Those two pages are the most-read documentation in the project and were the
least verified: the docs build executes ``.. jupyter-execute::`` blocks, but
``.. code-block:: python`` is only syntax-highlighted, and nothing at all
reads the README.  A pre-publish audit found the quickstart's second example
raising ``NameError`` (it used ``energy`` and ``L`` fifteen lines before the
page defined them) and its last one raising ``ValueError`` from the library's
own validation (the Hamiltonian was an ``...`` stub, so it returned None).
The README had the same stub problem.

Both are the first thing a new user copies, so an example that cannot run is
worse than an undocumented feature -- it looks like the library is broken.

The blocks are executed in order, sharing one namespace, because that is how
a reader follows the page: a name bound in one block is expected to still be
there in the next.
"""

import pathlib
import re
import textwrap
import warnings

import matplotlib
import pytest

matplotlib.use('Agg')

ROOT = pathlib.Path(__file__).resolve().parent.parent


def _readme_blocks():
    md = (ROOT/'README.md').read_text()
    return re.findall(r'```(?:python|py)\n(.*?)```', md, re.S)


def _rst_blocks(path):
    lines = path.read_text().splitlines(True)
    out, i = [], 0
    while i < len(lines):
        if re.match(r'\s*\.\. (?:code-block:: python|jupyter-execute::)', lines[i]):
            j = i + 1
            while j < len(lines) and (lines[j].strip() == ''
                                      or re.match(r'\s+:\w+:', lines[j])):
                j += 1
            base = len(lines[j]) - len(lines[j].lstrip()) if j < len(lines) else 0
            body = []
            while j < len(lines) and (lines[j].strip() == ''
                                      or (len(lines[j]) - len(lines[j].lstrip())) >= base):
                body.append(lines[j])
                j += 1
            block = textwrap.dedent(''.join(body))
            if block.strip():
                out.append(block)
            i = j
        else:
            i += 1
    return out


CASES = [
    ('README.md', _readme_blocks),
    ('quickstart.rst', lambda: _rst_blocks(ROOT/'docs'/'source'/'quickstart.rst')),
]


@pytest.mark.parametrize("name,loader", CASES, ids=[c[0] for c in CASES])
def test_documented_examples_run(name, loader):
    blocks = loader()
    assert blocks, "no code blocks found in %s -- the extractor has drifted" % name
    namespace = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for k, block in enumerate(blocks, start=1):
            try:
                exec(compile(block, '<%s block %d>' % (name, k), 'exec'), namespace)
            except Exception as exc:                       # noqa: BLE001 - report, don't mask
                pytest.fail('%s block %d raised %s: %s\n\n%s'
                            % (name, k, type(exc).__name__, exc, block))
