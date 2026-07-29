# -*- coding: utf-8 -*-
"""Pytest configuration: make the magnus package importable without installation."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The inner path (src/magnus) is needed because oscprob.py imports the
# top-level modules `version` and `authors`; the outer path (src) exposes the
# `magnus` package itself. Order matters: src must come first so that the
# package `magnus` wins over the inner `magnus/magnus` subdirectory.
sys.path.insert(0, os.path.join(ROOT, 'src', 'magnus'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
