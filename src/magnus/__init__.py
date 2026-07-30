# -*- coding: utf-8 -*-
r"""__init__.py

Top-level package initializer for magnus.  Exposes the package version
and lists the six subpackages that make up Magnus: earth, globaldefs,
hamiltonians, magnus, matter, and oscprob.  See
:doc:`/architecture` for how they fit together.

Routine listings
----------------

    (none; only re-exports __version__ and defines submodules/__all__)

Created: 2024/12/30 02:00
Last modified: 2026/07/30
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# __path__ = __import__('pkgutil').extend_path(__path__, __name__)

from .version import __version__

submodules = [
    'earth',
    'globaldefs',
    'hamiltonians',
    'magnus',
    'matter',
    'oscprob'
]

__all__ = submodules
