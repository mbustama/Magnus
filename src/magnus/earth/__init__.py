# -*- coding: utf-8 -*-
r"""__init__.py

Subpackage initializer for magnus.earth.  Re-exports everything from
earth.py: the PREM matter-density profile and Earth chord/zenith-angle
geometry.

Routine listings
----------------

    (none; only re-exports earth.py's public names)

Created: 2024/12/31 18:14
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


name = 'earth'

from .earth import *

__all__ = [s for s in dir() if not s.startswith('_')]
