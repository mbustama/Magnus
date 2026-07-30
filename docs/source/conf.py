# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Magnus'
copyright = '2026, Mauricio Bustamante'
author = 'Mauricio Bustamante'
release = '0.10.0'

# -- General configuration ----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# The package (src/magnus/) is not pip-installed in the doc-build environment,
# and several of its internal modules resolve their own imports relative to
# the current working directory rather than the package location, which makes
# them unsafe to actually *import* from a Sphinx worker (as sphinx.ext.autodoc
# would need to).  sphinx-autoapi instead parses the source with static
# analysis, so it needs only a path to the source tree, never an import.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))

extensions = [
    'sphinx.ext.napoleon',      # NumPy/Google-style docstrings
    'sphinx.ext.viewcode',      # Link API entries to highlighted source
    'sphinx.ext.mathjax',       # Render the LaTeX math in docstrings
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',        # Copy-to-clipboard button on code blocks
    'sphinxcontrib.bibtex',     # References page (refs.bib)
    'sphinxcontrib.mermaid',    # Architecture diagrams (architecture.rst)
    'myst_parser',              # Lets changelog.rst .. include:: the root CHANGELOG.md
    'jupyter_sphinx',           # Executes .. jupyter-execute:: blocks in docstring Examples
    'autoapi.extension',        # API reference, generated from src/magnus/
]

# jupyter-sphinx runs every `.. jupyter-execute::` block (used in the "Examples" section of
# several oscprob.py docstrings) through a real Jupyter kernel at build time, so the output
# shown in the docs is always the actual current output, not hand-typed text that can drift.
# Requires a "python3" kernel with `magnus` and its dependencies importable -- see the
# "Install magnus (editable)" / "Register the Jupyter kernel" steps in .github/workflows/pages.yml.
# For a local build: `pip install -e .` from the repo root, then
# `python -m ipykernel install --user --name python3` (harmless if a "python3" kernel already
# exists and points to an environment where `magnus` is importable).
jupyter_execute_default_kernel = 'python3'

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# -- API reference (sphinx-autoapi) -------------------------------------------

autoapi_type = 'python'
autoapi_dirs = [str(Path(__file__).resolve().parents[2] / 'src' / 'magnus')]
autoapi_root = 'api'
autoapi_own_page_level = 'module'
autoapi_add_toctree_entry = False  # Added explicitly in index.rst instead
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
]
# Skip the packaging-only version/authors modules and the 'old' prototypes
autoapi_ignore = ['*/old/*']

napoleon_google_docstring = False
napoleon_numpy_docstring = True

bibtex_bibfiles = ['refs.bib']

templates_path = ['_templates']
exclude_patterns = ['_build', 'sandbox']

# The master toctree document.
master_doc = 'index'
source_suffix = '.rst'

# -- Options for HTML output --------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/magnus_logo.png'

html_theme_options = {
    'logo_only': True,
    'navigation_depth': 4,
    'vcs_pageview_mode': 'edit',
}

html_context = {
    'display_github': True,
    'github_user': 'mbustama',
    'github_repo': 'Magnus',
    'github_version': 'main',
    'conf_py_path': '/docs/source/',
}
