#!/usr/bin/env python3
r"""Check the committed notebooks for names that are used but never defined.

Why this exists
---------------
``ruff check`` already runs over ``notebooks/make_notebooks.py``, but the generator
holds every notebook cell as a string literal, so the linter never looks inside them.
That blind spot is not theoretical.  ``28_magnus_paper_figures.ipynb`` called an
undefined ``TURB4_KW`` for weeks: the cell reads its value from
``paper_figure_cache.json`` and reaches the call only on a cache miss, so the notebook
ran cleanly until the configuration moved, and then continuous integration stopped with
a ``NameError``.

The ``notebooks/`` directory is not handed to ``ruff`` wholesale, because a notebook
re-imports in later cells so that a reader can start anywhere, and ruff reads that as an
F811 redefinition.  This script concatenates the code cells of one notebook into a
single module and asks for F821 alone, which that pattern does not trigger.

Continuous integration runs the committed ``.ipynb`` files rather than the generator, so
the committed files are what this checks.

Usage
-----
    python tools/lint_notebook_cells.py [notebook ...]

With no arguments it checks every notebook in ``notebooks/``.  It prints one line per
finding, naming the notebook, the undefined name and the line that uses it, and exits
non-zero if there is any.
"""
import json
import pathlib
import re
import subprocess
import sys
import tempfile

# A magic or shell escape, anchored at column zero.  The anchor matters: a continuation
# line may legitimately begin with the ``%`` formatting operator, and commenting one of
# those out breaks the parse of the whole file and hides every real finding behind it.
MAGIC = re.compile(r'^(%{1,2}[A-Za-z_]|![A-Za-z_])')

# Names a notebook gets from IPython rather than from an import.
BUILTINS = ['display', 'get_ipython']

FINDING = re.compile(r':(\d+):\d+: F821 Undefined name `([^`]+)`')

ROOT = pathlib.Path(__file__).resolve().parents[1]


def cells_as_module(path):
    """The code cells of one notebook, concatenated into the lines of one module."""
    nb = json.loads(path.read_text())
    lines = []
    for cell in nb['cells']:
        if cell['cell_type'] != 'code':
            continue
        for line in cell['source']:
            line = line.rstrip('\n')
            lines.append('# ' + line if MAGIC.match(line) else line)
        lines.append('')
    return lines


def check(path):
    """Findings for one notebook, as (message, is_error) pairs."""
    lines = cells_as_module(path)
    with tempfile.TemporaryDirectory() as tmp:
        module = pathlib.Path(tmp) / (path.stem + '.py')
        module.write_text('\n'.join(lines) + '\n')
        # --isolated: the repository's own ruff settings must not widen or narrow the
        # one rule this asks for.
        result = subprocess.run(
            ['ruff', 'check', '--isolated', '--no-cache', '--select', 'F821',
             '--config', 'builtins=%s' % json.dumps(BUILTINS),
             '--output-format', 'concise', str(module)],
            capture_output=True, text=True)
    if result.returncode not in (0, 1):
        return [('%s: ruff could not run: %s'
                 % (path.name, result.stderr.strip().splitlines()[0]), True)]
    out = []
    for line in result.stdout.splitlines():
        m = FINDING.search(line)
        if m:
            n, name = int(m.group(1)), m.group(2)
            out.append(('%s: undefined name %r in: %s'
                        % (path.name, name, lines[n - 1].strip()), False))
    return out


def main(argv):
    paths = ([pathlib.Path(a) for a in argv[1:]] if len(argv) > 1
             else sorted((ROOT / 'notebooks').glob('*.ipynb')))
    results = [r for p in paths for r in check(p)]
    for message, _ in results:
        print(message)
    if results:
        print('\nA cell that reads a cached value hides a name like this until the '
              'cache misses, and then the notebook stops with a NameError.')
        return 1
    print('%d notebooks checked, every name defined.' % len(paths))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
