Installation & Requirements
============================

Requirements
------------

Magνs requires **Python 3.10+**.  Dependencies:

* ``numpy``
* ``scipy >= 1.9`` (needed for stacked-input ``scipy.linalg.expm`` and
  ``scipy.integrate.cumulative_trapezoid``/``cumulative_simpson``)
* ``joblib`` (used to parallelize probability scans over energy/baseline
  points; a single-core install works fine with ``n_jobs=1``, the default)

See :download:`src/requirements.txt <../../src/requirements.txt>`.

Installation
------------

Magνs is not yet published on PyPI.  Clone the repository and add
``src/`` (and ``src/magnus/``, needed by a few modules that resolve
sibling imports directly) to your Python path:

.. code-block:: bash

   git clone https://github.com/mbustama/Magnus.git
   cd Magnus
   pip install -r src/requirements.txt

.. code-block:: python

   import sys
   sys.path.extend(['src', 'src/magnus'])

   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

Alternatively, an editable ``pip`` install (using the ``pyproject.toml`` at
the repository root) avoids the ``sys.path`` step and also installs the
``magnus`` command-line calculator (see :doc:`cli`):

.. code-block:: bash

   pip install -e .

   magnus prob --flavors 3 --environment vacuum --energy 1 --energy-unit GeV \
       --baseline 1300 --baseline-unit km

Verifying the Installation
---------------------------

After installing the dependencies, run the test suite to confirm everything
is configured correctly for your system:

.. code-block:: bash

   pip install pytest
   pytest tests/ -v

The same suite runs in CI on Python 3.10-3.12 on every push; see the badge
on the :doc:`index` page.

File Tree
---------

.. code-block:: text

   Magnus/
   ├── .github/
   │   └── workflows/
   │       ├── lint.yml                 # Ruff linter/formatter (informational; does not block CI)
   │       ├── pages.yml                # GitHub Pages deployment for this documentation
   │       ├── publish.yml              # PyPI (OIDC) automated publishing workflow, on GitHub Release
   │       └── tests.yml                # GitHub Actions CI testing pipeline (Python 3.10-3.12)
   ├── docs/                            # Sphinx documentation configuration and source
   │   ├── source/
   │   │   ├── conf.py                  # Sphinx build configuration (autoapi + napoleon + bibtex + mermaid + myst)
   │   │   ├── index.rst                # Master documentation page: overview, features, when Magnus wins
   │   │   ├── installation.rst         # This page
   │   │   ├── quickstart.rst           # Worked Python-API code examples for every entry point
   │   │   ├── cli.rst                  # Command-line calculator: flag reference and examples
   │   │   ├── functions.rst            # Full osc_prob_{2,3,4,5}nu_* listing, by environment/scenario
   │   │   ├── architecture.rst         # The wrapper/middle/primordial layering, with diagrams
   │   │   ├── methodology.rst          # The Magnus expansion, integrators, and performance engineering
   │   │   ├── tutorials.rst            # Guide to the numbered example notebooks in notebooks/
   │   │   ├── references.rst           # Bibliography page rendering
   │   │   ├── refs.bib                 # BibTeX citations for the Magnus-expansion and PREM literature
   │   │   └── changelog.rst            # Renders the root CHANGELOG.md via myst-parser
   │   ├── requirements.txt             # Sphinx + theme + extensions needed to build the docs
   │   ├── Makefile                     # Build commands for Unix
   │   └── make.bat                     # Build commands for Windows
   ├── fig/                             # Plots produced by the example notebooks
   ├── notebooks/                       # Numbered Jupyter notebooks -- see the Tutorials page
   │   ├── 01_magnus_introduction.ipynb
   │   ├── 02_magnus_2nu_vacuum_matter.ipynb
   │   ├── 03_magnus_3nu_vacuum_matter.ipynb
   │   ├── 04_magnus_long_baseline.ipynb
   │   ├── 05_magnus_biprobability.ipynb
   │   ├── 06_magnus_oscillograms.ipynb
   │   ├── 07_magnus_bsm_sterile_nu.ipynb
   │   ├── 08_magnus_bsm_nsi.ipynb
   │   ├── 09_magnus_bsm_liv.ipynb
   │   └── 10_magnus_matrix_exponential.ipynb
   ├── src/
   │   ├── magnus/                      # Main Python package
   │   │   ├── magnus/                  # Magnus-expansion numerical core
   │   │   │   └── magnus.py            # Term recursion, Gauss-Legendre integrators, batched kernel
   │   │   ├── oscprob/                 # Oscillation probabilities: main API
   │   │   │   ├── oscprob.py           # osc_prob and every physics-scenario wrapper
   │   │   │   └── oscprobstd.py        # Closed-form 2nu/3nu probabilities (used to validate oscprob.py)
   │   │   ├── hamiltonians/            # 2nu-5nu Hamiltonians: vacuum, matter, NSI, LIV
   │   │   ├── earth/                   # PREM density profile, chord/zenith-angle geometry
   │   │   ├── matter/                  # Density profiles, electron number density, CC potential
   │   │   ├── globaldefs/              # Units, physical constants, NuFit parameter sets
   │   │   ├── cli.py                   # `magnus` command-line calculator (also `python -m magnus`)
   │   │   └── __main__.py              # Entry point for `python -m magnus`
   │   └── requirements.txt
   ├── tests/                           # Test suite (pytest; runs in CI)
   │   ├── conftest.py                  # Path setup so magnus is importable without installation
   │   ├── test_magnus_expansion.py     # Magnus-core correctness (terms, orders, GL rates, unitarity)
   │   ├── test_oscprob.py              # Oscillation-probability engine, closed-form and ODE cross-checks
   │   ├── test_earth_matter.py         # PREM profile, chord geometry, electron density
   │   ├── test_hamiltonians.py         # Hamiltonian/mixing-matrix builders
   │   └── test_cli.py                  # magnus command-line calculator
   ├── .gitignore
   ├── CHANGELOG.md                     # Version history (Keep a Changelog format)
   ├── pyproject.toml                   # Build system, dependencies, and the `magnus` console-script entry point
   └── README.md
