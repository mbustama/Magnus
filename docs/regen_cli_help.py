#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Regenerate the verbatim ``--help`` blocks embedded in docs/source/cli.rst.

``cli.rst`` quotes the output of ``magnus --help`` and ``magnus prob --help`` so that
readers get the full flag reference without running anything.  Hand-copied output goes
stale silently the moment a flag is added (which is how ``--version`` came to be missing
from that page), so this script rewrites those two blocks straight from the live
``argparse`` parser instead.

Run it from anywhere after touching ``src/magnus/cli.py``::

    python3 docs/regen_cli_help.py

It edits ``docs/source/cli.rst`` in place and reports whether anything changed, so it is
safe to run in a CI check (``--check`` exits non-zero if the page is out of date).
"""

import argparse
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src'))

CLI_RST = REPO / 'docs' / 'source' / 'cli.rst'

# argparse wraps to the terminal width; pin it so the committed output is reproducible
# regardless of the terminal the script happens to run in.
os.environ['COLUMNS'] = '90'

from magnus.cli import build_parser  # noqa: E402


def _help_text(which: str) -> str:
    parser = build_parser()
    if which == 'top':
        return parser.format_help()
    # The 'prob' subparser is reachable only through the subparsers action.
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices['prob'].format_help()
    raise AssertionError("no subparsers found on the top-level parser")


def _indent(text: str) -> str:
    return ''.join('   ' + line if line.strip() else line
                   for line in text.splitlines(keepends=True)).rstrip() + '\n'


def _replace_block(rst: str, first_line: str, new_body: str) -> str:
    """Replace the ``.. code-block:: text`` whose first content line starts with first_line.

    Works line-wise rather than with a regex: the block runs from the directive to the first
    line that is neither indented nor blank, which is unambiguous and, unlike a lazy regex over
    optional blank lines, gives byte-identical output when run twice (so ``--check`` is stable).
    """
    lines = rst.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.strip() != '.. code-block:: text':
            continue
        # First non-blank line after the directive decides whether this is the block we want.
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines) or not lines[j].lstrip().startswith(first_line):
            continue
        # The block ends at the first line that is neither blank nor indented.
        end = j
        while end < len(lines) and (not lines[end].strip() or lines[end].startswith('   ')):
            end += 1
        return ''.join(lines[:i + 1]) + '\n' + new_body + '\n' + ''.join(lines[end:])
    raise SystemExit(f"could not find a code block starting with {first_line!r} in {CLI_RST}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', action='store_true',
                    help='Do not write; exit 1 if cli.rst is out of date.')
    args = ap.parse_args()

    original = CLI_RST.read_text(encoding='utf-8')
    # Only the `magnus prob --help` output is embedded verbatim; the top-level help appears
    # on the page solely inside a short error-handling example, which is prose, not a dump.
    updated = _replace_block(original, 'usage: magnus prob', _indent(_help_text('prob')))

    if updated == original:
        print(f"{CLI_RST.relative_to(REPO)} is up to date.")
        return 0
    if args.check:
        print(f"{CLI_RST.relative_to(REPO)} is OUT OF DATE; run: python3 docs/regen_cli_help.py")
        return 1
    CLI_RST.write_text(updated, encoding='utf-8')
    print(f"{CLI_RST.relative_to(REPO)} updated.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
