#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Executes every ``.. jupyter-execute::`` block in docs/source/*.rst.

The documentation's code examples are run for real at build time by
jupyter_sphinx, so a broken one fails CI -- but only in the *full* build, which
takes minutes and needs a Jupyter kernel.  The fast build used while writing
docs stubs those directives out precisely so it does not need one, which means
it validates the prose and the cross-references while saying nothing at all
about the code.

That gap is not hypothetical.  ``averaged_probability.rst`` shipped a snippet
that splatted a whole ``OSC_PARAMS_PREDEFINED`` entry into a probability
function; the fast build reported success, and only the full build (or running
the snippet by hand) would have found it.

This script closes the gap in about a second::

    python3 docs/check_doc_snippets.py            # every page
    python3 docs/check_doc_snippets.py methodology averaged_probability

Each page's blocks run in order, in one shared namespace, so a later block can
use a name an earlier one defined -- the same way jupyter_sphinx executes them.
Blocks are *not* shared between pages.

Two directive options are honoured, because ignoring them would report failures
that the real build does not: ``:raises:`` marks a block that is *expected* to
raise, and ``:stderr:`` merely routes stderr into the page.
"""

import ast
import re
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / 'docs' / 'source'
PACKAGE = ROOT / 'src' / 'magnus'

DIRECTIVE = re.compile(r'^(?P<indent>[ \t]*)\.\.[ \t]+jupyter-execute::[ \t]*$')
OPTION = re.compile(r'^[ \t]*:(?P<name>[a-z-]+):(?P<value>.*)$')


def extract_blocks(path):
    """[(line_number, source, options)] for every jupyter-execute block in `path`."""
    return extract_blocks_from_text(path.read_text(encoding='utf-8'))


def extract_blocks_from_text(text):
    """The directive parser itself, shared by the RST and docstring routes."""
    lines = text.split('\n')
    blocks = []
    i = 0
    while i < len(lines):
        match = DIRECTIVE.match(lines[i])
        if match is None:
            i += 1
            continue

        start = i + 1
        indent = len(match.group('indent'))
        i += 1

        # Directive options come first, then a blank line, then the content.
        options = {}
        while i < len(lines):
            if not lines[i].strip():
                i += 1
                continue
            if len(lines[i]) - len(lines[i].lstrip()) <= indent:
                break
            opt = OPTION.match(lines[i])
            if opt is None:
                break
            options[opt.group('name')] = opt.group('value').strip()
            i += 1

        body = []
        while i < len(lines):
            line = lines[i]
            if not line.strip():
                body.append('')
                i += 1
                continue
            if len(line) - len(line.lstrip()) <= indent:
                break
            body.append(line)
            i += 1

        if body:
            # Strip the common indentation the directive imposes.
            widths = [len(l) - len(l.lstrip()) for l in body if l.strip()]
            cut = min(widths) if widths else 0
            blocks.append((start + 1, '\n'.join(l[cut:] for l in body).strip('\n'), options))

    return blocks


def extract_blocks_from_docstrings(path):
    """[(line_number, source, options)] for every jupyter-execute block in the
    docstrings of a Python module.

    Most of this package's executable examples live in docstrings, not in the
    RST pages: autoapi renders them into the API reference, where jupyter_sphinx
    runs them exactly as it runs the ones written on a page.  A checker that
    looked only at docs/source/ would therefore miss the large majority of them.
    """
    tree = ast.parse(path.read_text(encoding='utf-8'))
    blocks = []

    # Every documentation string in a module -- a function's, a class's, the module's own,
    # and the bare strings this package writes under its module-level constants -- is a
    # string expression standing alone as a statement.  Collecting them by that one rule
    # covers all four uniformly; using ast.get_docstring as well would count the first
    # three twice, since a docstring is also the first statement of its own body.
    for node in ast.walk(tree):
        body = getattr(node, 'body', None)
        # Lambda and IfExp also have a `body`, but it is a single node, not a list.
        for child in (body if isinstance(body, list) else []):
            if (isinstance(child, ast.Expr) and isinstance(child.value, ast.Constant)
                    and isinstance(child.value.value, str)
                    and 'jupyter-execute' in child.value.value):
                for lineno, source, options in extract_blocks_from_text(child.value.value):
                    blocks.append((child.lineno + lineno - 1, source, options))

    return blocks


def check_page(path, verbose=False):
    """Runs every block in one page or module. Returns the number of failures."""
    blocks = (extract_blocks_from_docstrings(path) if path.suffix == '.py'
              else extract_blocks(path))
    if not blocks:
        return 0

    namespace = {'__name__': '__doc_snippet__'}
    failures = 0
    for lineno, source, options in blocks:
        expects_raise = 'raises' in options
        try:
            exec(compile(source, f'{path.name}:{lineno}', 'exec'), namespace)
        except Exception:
            if expects_raise:
                continue
            failures += 1
            print(f'\n{path.name}:{lineno}: jupyter-execute block failed\n')
            print('\n'.join('    ' + l for l in source.split('\n')))
            print()
            print('\n'.join('    ' + l for l in traceback.format_exc().rstrip().split('\n')))
        else:
            if expects_raise:
                failures += 1
                print(f'\n{path.name}:{lineno}: block is marked :raises: but did not raise\n')
                print('\n'.join('    ' + l for l in source.split('\n')))

    if verbose and failures == 0:
        print(f'{path.name}: {len(blocks)} block(s) OK')

    return failures


def main(argv):
    wanted = [a for a in argv[1:] if not a.startswith('--')]
    rst_only = '--rst-only' in argv[1:]

    pages = sorted(SOURCE.glob('*.rst'))
    if not rst_only:
        pages += sorted(PACKAGE.rglob('*.py'))

    if wanted:
        stems = {w.removesuffix('.rst').removesuffix('.py') for w in wanted}
        pages = [p for p in pages if p.stem in stems]
        missing = stems - {p.stem for p in pages}
        if missing:
            sys.exit(f'no such page or module: {sorted(missing)}')

    total_blocks = sum(len(extract_blocks_from_docstrings(p) if p.suffix == '.py'
                           else extract_blocks(p)) for p in pages)
    failures = sum(check_page(p, verbose=True) for p in pages)

    print(f'\n{total_blocks} jupyter-execute block(s) in {len(pages)} file(s): '
          f'{"all OK" if failures == 0 else str(failures) + " FAILED"}')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
