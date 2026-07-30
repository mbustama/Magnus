# -*- coding: utf-8 -*-
r"""__main__.py

Entry point for ``python -m magnus``. Delegates to :func:`magnus.cli.main`.

Routine listings
----------------

    (none; only calls magnus.cli.main() when run as a script)
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import sys

from magnus.cli import main

if __name__ == "__main__":
    sys.exit(main())
