# -*- coding: utf-8 -*-
r"""__init__.py

Top-level package initializer for magnus. Exposes the package version and
explicitly imports the seven modules that make up Magnus's public API:
earth, globaldefs, hamiltonians, magnus, matter, oscprob, and oscprobstd.
See :doc:`/architecture` for how they fit together.

``authors`` and ``version`` are internal metadata modules (used by
:func:`magnus.oscprob.print_banner` and the ``magnus`` command-line
calculator), not part of the public API, and are intentionally not
listed in ``submodules``/``__all__`` or documented in the API reference.
``__version__`` itself is still exposed directly, as is standard Python
packaging convention.

Routine listings
----------------

    (none; only re-exports __version__ and defines submodules/__all__)
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from .version import __version__

from . import earth
from . import globaldefs
from . import hamiltonians
from . import magnus
from . import matter
from . import oscprob
from . import oscprobstd

submodules = [
    'earth',
    'globaldefs',
    'hamiltonians',
    'magnus',
    'matter',
    'oscprob',
    'oscprobstd',
]

__all__ = submodules
