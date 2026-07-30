# -*- coding: utf-8 -*-
r"""__init__.py

Subpackage initializer for magnus.matter.  Re-exports everything from
matter.py: generic (non-Earth, non-Sun) matter-density profiles,
electron number density, and the coherent forward-scattering potential.

Routine listings
----------------

    (none; only re-exports matter.py's public names)

Created: 2024/12/30 02:00
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


name = 'matter'

from .matter import *

__all__ = [s for s in dir() if not s.startswith('_')]
