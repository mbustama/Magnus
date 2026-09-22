# -*- coding: utf-8 -*-
r"""The worked examples in docs/source/cli.rst and README.md print what the CLI prints.

``cli.rst`` says its examples are "output captured from this version", and quotes
six of them verbatim.  Nothing checked that, so all of them went stale: the banner
read ``Magnus 1.0.0rc1`` against a package at 1.1.0, and every probability matrix
predated the NuFIT 6.1 parameter move -- the Earth example's nu_mu -> nu_tau entry
by nine percent.  ``docs/regen_cli_help.py`` already keeps the ``--help`` dump
honest; this does the same for the examples, at about 3.6 s for all of them.

The README quotes one more, in a fenced block, and went stale the same way: its banner
still read 1.1.0 when the package moved to 1.1.1, with nothing to say so.  Both files are
read here, so a version bump that forgets either fails the build.

Blocks are skipped rather than checked when they cannot be compared verbatim: the
JSON example elides its array with ``...``, and the error-handling example runs two
commands in one block to show both messages.

To adopt new output after a deliberate change::

    python tests/test_cli_examples_match.py --write
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CLI_RST = REPO / 'docs' / 'source' / 'cli.rst'
README = REPO / 'README.md'
BLOCK = '.. code-block:: text'
FENCE = '```'
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


def _fenced(lines):
    """Yield (first_line_index, end_index, body_lines) for each fenced block with a language."""
    i = 0
    while i < len(lines):
        if lines[i].startswith(FENCE) and lines[i].strip() != FENCE:
            end = i + 1
            while end < len(lines) and not lines[end].startswith(FENCE):
                end += 1
            yield i + 1, end, lines[i + 1:end]
            i = end + 1
        else:
            i += 1


# Where each file keeps its examples, and the indentation its blocks carry.
SOURCES = ((CLI_RST, _blocks, INDENT), (README, _fenced, ''))


def _parse(text):
    """(command, expected stdout, lines consumed by the command) for one quoted call, or None."""
    if not text or not text[0].startswith('$ magnus prob'):
        return None
    if sum(1 for t in text if t.startswith('$ ')) > 1:
        return None                       # the two-command error example
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
        return None                       # the JSON example elides its array
    return command, expected, k


def _examples():
    """Every verbatim ``$ magnus prob`` example: (file, command, expected stdout, span)."""
    out = []
    for path, blocks, indent in SOURCES:
        lines = path.read_text(encoding='utf-8').splitlines()
        for start, end, body in blocks(lines):
            text = [b[len(indent):] if indent and b.startswith(indent) else b for b in body]
            parsed = _parse(text)
            if parsed is not None:
                command, expected, k = parsed
                out.append((path, command, expected, (start + k, end)))
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


def test_each_file_quotes_examples_at_all():
    """A parser that silently matches nothing would make every check below vacuous."""
    found = {path: sum(1 for e in EXAMPLES if e[0] == path) for path, _, _ in SOURCES}
    assert found[CLI_RST] >= 4, "found %d verbatim examples in cli.rst" % found[CLI_RST]
    assert found[README] >= 1, "found no verbatim example in README.md"


@pytest.mark.parametrize('path,command,expected,_span', EXAMPLES,
                         ids=['%s: %s' % (e[0].name, e[1][:40]) for e in EXAMPLES])
def test_documented_cli_example_still_prints_what_the_page_says(path, command, expected,
                                                                _span):
    actual = _run(command)
    assert actual == [e.rstrip() for e in expected], (
        "%s is out of date for:\n  %s\n"
        "expected:\n%s\nactual:\n%s\n"
        "Run: python tests/test_cli_examples_match.py --write"
        % (path.relative_to(REPO), command, '\n'.join(expected), '\n'.join(actual)))


def _write():
    for path, _, indent in SOURCES:
        lines = path.read_text(encoding='utf-8').splitlines(keepends=True)
        changed = 0
        for _, command, expected, (start, end) in reversed(
                [e for e in _examples() if e[0] == path]):
            actual = _run(command)
            if actual == [e.rstrip() for e in expected]:
                continue
            lines[start:end] = [(indent + ln).rstrip() + '\n' if ln else '\n'
                                for ln in actual]
            changed += 1
        if changed:
            path.write_text(''.join(lines), encoding='utf-8')
        print("%s: %d example(s) rewritten." % (path.relative_to(REPO), changed))
    return 0


if __name__ == '__main__':
    raise SystemExit(_write() if '--write' in sys.argv else
                     pytest.main([__file__, '-q']))
