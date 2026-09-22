# -*- coding: utf-8 -*-
r"""The worked examples in docs/source/cli.rst print what the CLI prints.

``cli.rst`` says its examples are "output captured from this version", and quotes
six of them verbatim.  Nothing checked that, so all of them went stale: the banner
read ``Magnus 1.0.0rc1`` against a package at 1.1.0, and every probability matrix
predated the NuFIT 6.1 parameter move -- the Earth example's nu_mu -> nu_tau entry
by nine percent.  ``docs/regen_cli_help.py`` already keeps the ``--help`` dump
honest; this does the same for the examples, at about 3.6 s for all of them.

Blocks are skipped rather than checked when they cannot be compared verbatim: the
JSON example elides its array with ``...``, and the error-handling example runs two
commands in one block to show both messages.

To adopt new output after a deliberate change::

    python tests/test_cli_examples_match.py --write
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CLI_RST = REPO / 'docs' / 'source' / 'cli.rst'
BLOCK = '.. code-block:: text'
INDENT = '   '


def _blocks(lines):
    """Yield (first_line_index, end_index, body_lines) for each text code block."""
    for i, line in enumerate(lines):
        if line.rstrip() != BLOCK:
            continue
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        end = j
        while end < len(lines) and (not lines[end].strip()
                                    or lines[end].startswith(INDENT)):
            end += 1
        while end > j and not lines[end - 1].strip():
            end -= 1
        yield j, end, lines[j:end]


def _examples():
    """Every verbatim ``$ magnus prob`` example: (command, expected stdout, span)."""
    lines = CLI_RST.read_text(encoding='utf-8').splitlines()
    out = []
    for start, end, body in _blocks(lines):
        text = [b[len(INDENT):] if b.startswith(INDENT) else b for b in body]
        if not text or not text[0].startswith('$ magnus prob'):
            continue
        if sum(1 for t in text if t.startswith('$ ')) > 1:
            continue                      # the two-command error example
        argv, k = [], 0
        while k < len(text):
            piece = text[k]
            argv.append(piece[2:] if k == 0 else piece)
            k += 1
            if not piece.rstrip().endswith('\\'):
                break
        command = ' '.join(a.rstrip().rstrip('\\').strip() for a in argv)
        expected = text[k:]
        if any(t.strip() == '...' or '...' in t for t in expected):
            continue                      # the JSON example elides its array
        out.append((command, expected, (start + k, end)))
    return out


def _run(command):
    assert command.startswith('magnus prob ')
    env = dict(os.environ, COLUMNS='90')
    r = subprocess.run([sys.executable, '-m', 'magnus', 'prob']
                       + command[len('magnus prob '):].split(),
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, "%r exited %d: %s" % (command, r.returncode, r.stderr[-400:])
    return [ln.rstrip() for ln in r.stdout.rstrip('\n').splitlines()]


EXAMPLES = _examples()


def test_the_page_quotes_examples_at_all():
    """A parser that silently matches nothing would make every check below vacuous."""
    assert len(EXAMPLES) >= 4, "found %d verbatim examples in cli.rst" % len(EXAMPLES)


@pytest.mark.parametrize('command,expected,_span',
                         EXAMPLES, ids=[e[0][:48] for e in EXAMPLES])
def test_documented_cli_example_still_prints_what_the_page_says(command, expected, _span):
    actual = _run(command)
    assert actual == [e.rstrip() for e in expected], (
        "docs/source/cli.rst is out of date for:\n  %s\n"
        "expected:\n%s\nactual:\n%s\n"
        "Run: python tests/test_cli_examples_match.py --write"
        % (command, '\n'.join(expected), '\n'.join(actual)))


def _write():
    lines = CLI_RST.read_text(encoding='utf-8').splitlines(keepends=True)
    changed = 0
    for command, expected, (start, end) in reversed(_examples()):
        actual = _run(command)
        if actual == [e.rstrip() for e in expected]:
            continue
        lines[start:end] = [(INDENT + ln).rstrip() + '\n' if ln else '\n'
                            for ln in actual]
        changed += 1
    if changed:
        CLI_RST.write_text(''.join(lines), encoding='utf-8')
    print("%s: %d example(s) rewritten." % (CLI_RST.relative_to(REPO), changed))
    return 0


if __name__ == '__main__':
    raise SystemExit(_write() if '--write' in sys.argv else
                     pytest.main([__file__, '-q']))
