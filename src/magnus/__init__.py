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