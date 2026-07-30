# -*- coding: utf-8 -*-
r"""__init__.py

Subpackage initializer for magnus.magnus.  Re-exports everything from
magnus.py (the Magnus-expansion numerical core: term recursion,
Gauss-Legendre integrators, and the batched matrix-exponential kernel).

Routine listings
----------------

    (none; only re-exports magnus.py's public names)

Created: 2024/12/30 02:00
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


name = 'magnus'

from .magnus import *

__all__ = [s for s in dir() if not s.startswith('_')]
