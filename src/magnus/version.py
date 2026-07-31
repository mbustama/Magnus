# -*- coding: utf-8 -*-
r"""version.py

Resolves Magnus's version number.

There is exactly one place where the version is written down: the
``version`` field of ``pyproject.toml`` at the repository root.  This
module reads it from there rather than repeating it, so the number can
never drift between the packaging metadata, the ``magnus --version``
output, and the docs.

Resolution order:

#. ``importlib.metadata``, when Magnus is installed (``pip install .`` /
   ``pip install -e .``).  This is the normal case and costs one cheap
   lookup of the installed distribution metadata.
#. ``pyproject.toml``, read directly, when Magnus is *not* installed and
   is merely being imported off ``src/`` on ``sys.path`` (the layout used
   by ``tests/conftest.py``).
#. ``'0.0.0+unknown'``, if neither is available (e.g., the package was
   vendored without its metadata).  A recognizable placeholder is better
   than a stale hard-coded number that silently disagrees with the real
   release.

Routine listings
----------------

    (none; only resolves the __version__ string)
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import re
from pathlib import Path


def _version_from_metadata() -> str:
    """Version of the installed ``magnus`` distribution, or '' if not installed."""
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:                                       # pragma: no cover
        return ''
    try:
        return version('magnus')
    except PackageNotFoundError:
        return ''


def _version_from_pyproject() -> str:
    """Version field of the repository's pyproject.toml, or '' if unavailable.

    Parsed with a regex rather than ``tomllib``, which is only in the standard
    library from Python 3.11 on, while Magnus supports 3.10.
    """
    # src/magnus/version.py -> src/magnus -> src -> repository root
    pyproject = Path(__file__).resolve().parents[2] / 'pyproject.toml'
    try:
        text = pyproject.read_text(encoding='utf-8')
    except OSError:
        return ''
    match = re.search(r'^\s*version\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
    return match.group(1) if match else ''


__version__ = _version_from_metadata() or _version_from_pyproject() or '0.0.0+unknown'
