# -*- coding: utf-8 -*-
r"""version.py

Package version string, the single source of truth for Magnus's
version number.  Re-exported by magnus/__init__.py and printed by
oscprob.print_banner() and oscprob.print_run_parameters().  Must be kept
in sync with the ``version`` field in pyproject.toml and with
``release`` in docs/source/conf.py.

Routine listings
----------------

    (none; only defines the __version__ string)

Created: 2025/01/26 21:37
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


__version__ = '0.10.0'
