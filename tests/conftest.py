# -*- coding: utf-8 -*-
"""Pytest configuration: make the magnus package importable without installation."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, os.path.join(ROOT, 'src'))
