# -*- coding: utf-8 -*-
r"""__init__.py

Subpackage initializer for magnus.hamiltonians.  Re-exports everything
from hamiltonians2nu.py, hamiltonians3nu.py, hamiltonians4nu.py, and
hamiltonians5nu.py: the mixing matrices and vacuum/matter/NSI/LIV
Hamiltonians for 2, 3, 4, and 5 neutrino flavors.

Routine listings
----------------

    (none; only re-exports the four hamiltonians{2,3,4,5}nu.py modules)

Created: 2024/12/30 02:00
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


name = 'hamiltonians'

from .hamiltonians2nu import *
from .hamiltonians3nu import *
from .hamiltonians4nu import *
from .hamiltonians5nu import *

__all__ = [s for s in dir() if not s.startswith('_')]
