# -*- coding: utf-8 -*-
r"""__init__.py

Subpackage initializer for magnus.oscprob.  Re-exports everything from
oscprob.py (the Magnus-based oscillation-probability API: osc_prob and
every physics-scenario wrapper) and oscprobstd.py (closed-form
probabilities used to validate oscprob.py).

Routine listings
----------------

    (none; only re-exports oscprob.py's and oscprobstd.py's public names)

Created: 2024/12/30 02:00
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


name = 'oscprob'

from .oscprob import *
from .oscprobstd import *

__all__ = [s for s in dir() if not s.startswith('_')]
