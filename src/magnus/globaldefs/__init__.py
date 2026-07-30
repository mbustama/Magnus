# -*- coding: utf-8 -*-
r"""__init__.py

Subpackage initializer for magnus.globaldefs.  Re-exports everything
from globaldefs.py: physical constants, unit-conversion factors, and
predefined oscillation-parameter sets (e.g. NuFit 6.0).

Routine listings
----------------

    (none; only re-exports globaldefs.py's public names)

Created: 2024/12/30 02:00
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


name = 'globaldefs'

from .globaldefs import *

__all__ = [s for s in dir() if not s.startswith('_')]
