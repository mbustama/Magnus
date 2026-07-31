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

Install the current release candidate from PyPI:

.. code-block:: bash

   pip install --pre magnus

The ``--pre`` flag is needed while the current release is a release
candidate (see :doc:`changelog`); it can be dropped once 1.0.0 final is
out.  This installs the dependencies and the ``magnus`` command-line
calculator (see :doc:`cli`).

To work from a checkout instead -- to follow development, or to modify the
code -- install it in editable mode using the ``pyproject.toml`` at the
repository root:

.. code-block:: bash

   git clone https://github.com/mbustama/Magnus.git
   cd Magnus
   pip install -e .

.. code-block:: bash

   magnus prob --flavors 3 --environment vacuum --energy 1 --energy-unit GeV \
       --baseline 1300 --baseline-unit km

If you would rather not install anything, put ``src/`` on your Python
path instead.  That is all that is needed -- every module imports through
the ``magnus`` package (``import magnus.globaldefs``), so the package
directory itself does not belong on the path:

.. code-block:: bash

   pip install -r src/requirements.txt

.. code-block:: python

   import sys
   sys.path.insert(0, 'src')

   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

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
   │       ├── lint.yml                 # Ruff lint (blocking) + CLI-reference drift check
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
   │   │   ├── adiabatic_strategy.rst   # The adiabatic + Magnus hybrid strategy: derivation, diagrams, validation
   │   │   ├── tutorials.rst            # Guide to the numbered example notebooks in notebooks/
   │   │   ├── references.rst           # Bibliography page rendering
   │   │   ├── refs.bib                 # BibTeX citations for the Magnus-expansion and PREM literature
   │   │   └── changelog.rst            # Renders the root CHANGELOG.md via myst-parser
   │   ├── requirements.txt             # Sphinx + theme + extensions needed to build the docs
   │   ├── make_figures.py              # Regenerates the data-driven SVG in source/_static/
   │   ├── regen_cli_help.py            # Regenerates the --help block quoted in source/cli.rst
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
   │   ├── 10_magnus_matrix_exponential.ipynb
   │   ├── 11_magnus_adiabatic_hybrid_strategy.ipynb
   │   ├── matplotlibrc                 # Shared plot styling for the notebooks
   │   └── README.md                    # Per-notebook description and suggested reading order
   ├── src/
   │   ├── magnus/                      # Main Python package
   │   │   ├── __init__.py              # Explicitly imports/exposes the 8 public modules below
   │   │   ├── magnus.py                # Magnus-expansion numerical core: term recursion, GL integrators, batched kernel
   │   │   ├── adiabatic.py             # Adiabatic transport + Magnus-patch hybrid strategy (strategy='hybrid'/'auto')
   │   │   ├── oscprob.py               # osc_prob and every physics-scenario wrapper (main API)
   │   │   ├── oscprobstd.py            # Closed-form 2nu/3nu probabilities (used to validate the wrapper API)
   │   │   ├── hamiltonians/            # 2nu-5nu Hamiltonians: vacuum, matter, NSI, LIV (the one true subpackage)
   │   │   │   ├── __init__.py          # Explicit named imports from the four hamiltonians{2,3,4,5}nu.py modules
   │   │   │   ├── hamiltonians2nu.py
   │   │   │   ├── hamiltonians3nu.py
   │   │   │   ├── hamiltonians4nu.py
   │   │   │   └── hamiltonians5nu.py
   │   │   ├── earth.py                 # PREM density profile, chord/zenith-angle geometry
   │   │   ├── matter.py                # Density profiles, electron number density, CC potential
   │   │   ├── globaldefs.py            # Units, physical constants, NuFit parameter sets
   │   │   ├── cli.py                   # `magnus` command-line calculator (also `python -m magnus`)
   │   │   ├── __main__.py              # Entry point for `python -m magnus`
   │   │   ├── authors.py               # Package author string (internal; not part of the public API)
   │   │   └── version.py               # Resolves the version from pyproject.toml (internal)
   │   └── requirements.txt
   ├── tests/                           # Test suite (pytest; runs in CI)
   │   ├── conftest.py                  # Path setup so magnus is importable without installation
   │   ├── test_magnus_expansion.py     # Magnus-core correctness (terms, orders, GL rates, unitarity)
   │   ├── test_adiabatic.py            # Adiabatic propagator, resonance detection, hybrid certification
   │   ├── test_oscprob.py              # Oscillation-probability engine, closed-form and ODE cross-checks
   │   ├── test_earth_matter.py         # PREM profile, chord geometry, electron density
   │   ├── test_hamiltonians.py         # Hamiltonian/mixing-matrix builders
   │   ├── test_cli.py                  # magnus command-line calculator
   │   └── test_globaldefs.py           # NuFit historical parameter dict/loader
   ├── .gitignore
   ├── CHANGELOG.md                     # Version history (Keep a Changelog format)
   ├── pyproject.toml                   # Build system, dependencies, and the `magnus` console-script entry point
   └── README.md
