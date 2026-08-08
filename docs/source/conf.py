# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Magnus'
copyright = '2026, Mauricio Bustamante'
author = 'Mauricio Bustamante'

# -- General configuration ----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# The package (src/magnus/) is not pip-installed in the doc-build environment,
# and several of its internal modules resolve their own imports relative to
# the current working directory rather than the package location, which makes
# them unsafe to actually *import* from a Sphinx worker (as sphinx.ext.autodoc
# would need to).  sphinx-autoapi instead parses the source with static
# analysis, so it needs only a path to the source tree, never an import.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))

# Read the version from the same single source everything else uses -- the
# ``version`` field of pyproject.toml, resolved by magnus/version.py -- rather
# than repeating the number here, where it would silently go stale.  This is a
# leaf module with no third-party imports, so it is safe to import even though
# the rest of the package is not (see the note above).
from magnus.version import __version__ as release  # noqa: E402,F401

# sphinx_rtd_theme renders `version`, not `release`, under the logo; setting
# only the latter is why the sidebar showed no version at all.
version = release

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

# Without this Sphinx waits with `connect timeout=None`, so an unreachable
# inventory host stalls the build rather than failing it; ten seconds bounds
# that.  A fetch that fails is still a warning and so still fatal under -W,
# which is intended: a mapping listed here is one whose links are meant to
# resolve.
intersphinx_timeout = 10

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
# Source order, not alphabetical: the modules are written so that reading them
# top to bottom follows the method, and sorting the members destroys that.
autoapi_member_order = 'bysource'

autoapi_ignore = ['*/old/*', '*/authors.py', '*/version.py']

# Excluding version.py/authors.py above is deliberate (they are internal
# metadata modules, not part of the public API), but other modules still
# import from them (e.g. magnus/__init__.py, oscprob.py, cli.py), and
# autoapi warns that it cannot resolve those imports since it never scans
# the excluded modules. That warning is expected here, not a real problem.
# 'myst.header': changelog.rst includes CHANGELOG.md from its second line,
# dropping the file's own "# Changelog" heading so the sidebar does not carry a
# second, identical entry nested under the page title.  myst then reports that
# the included document starts at H2, which is the intent rather than a mistake.
# 'toc.not_included': sphinx-autoapi generates a package index page titled after
# the package, which put a redundant "magnus" level between API Reference and the
# modules.  api_reference.rst lists the module pages directly instead, leaving
# that one page in no toctree on purpose.
suppress_warnings = ['autoapi.python_import_resolution', 'myst.header',
                     'toc.not_included']

napoleon_google_docstring = False
napoleon_numpy_docstring = True

bibtex_bibfiles = ['refs.bib']

templates_path = ['_templates']
exclude_patterns = ['_build']

# The master toctree document.
master_doc = 'index'
source_suffix = '.rst'

# -- Options for HTML output --------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/magnus_logo.png'

# `release` is read from pyproject.toml at the top of this file, and putting it
# in the title is the only place a reader sees which version these pages
# describe -- which matters while the version is a release candidate.
html_title = 'Magnus %s' % release

html_theme_options = {
    # False, not True: `logo_only` hides `html_title`, and with it the version.
    'logo_only': False,
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
