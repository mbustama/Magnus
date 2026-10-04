Installation and requirements
=============================

Requirements
------------

Magνs requires **Python 3.10+**.  Dependencies:

* ``numpy >= 1.22``
* ``scipy >= 1.9`` (stacked-input ``scipy.linalg.expm`` and
  ``scipy.integrate.cumulative_trapezoid``; ``cumulative_simpson`` is used
  when available, from SciPy 1.12)
* ``joblib >= 1.2`` (used to parallelize probability scans over energy/baseline
  points; a single-core install works fine with ``n_jobs=1``, the default)
* ``matplotlib >= 3.5`` (for :mod:`magnus.plotting`, imported lazily)
* ``numba >= 0.59`` (the compiled matrix-exponential kernels of
  :mod:`magnus.expmkernels`)

See :download:`src/requirements.txt <../../src/requirements.txt>`.

Magνs is licensed under the GNU General Public License v3.0 only
(``GPL-3.0-only``); see :doc:`index` for a summary, and the ``LICENSE``
file in the repository root, which is also shipped inside the installed
distribution, for the full text.

Installation
------------

Install from PyPI:

.. code-block:: bash

   pip install magnuspy

This installs the dependencies and the ``magnus`` command-line
calculator (see :doc:`cli`).

.. note::
   The distribution is published as **magnuspy**, because ``magnus`` was
   already taken on PyPI.  The import package and the command-line tool are
   both **magnus**: ``pip install magnuspy``, then ``import magnus``.

For most users that is the whole installation.  A dedicated virtual environment
(``conda`` or ``venv``) is supported but not required.

The PyPI package does not include the notebooks or the test suite.  The notebooks have
saved output and can be read in the repository; to run or modify them, or to get changes
not yet released, clone the repository and install from it:

.. code-block:: bash

   git clone https://github.com/mbustama/Magnus.git
   cd Magnus
   pip install -e .                  # the package
   pip install -e ".[notebooks]"     # plus what the notebooks need
   pip install -e ".[test]"          # plus what the tests need

The notebooks are in ``notebooks/``; opening them interactively also needs a Jupyter front
end, such as JupyterLab.  Only contributors typically need the test suite.  Wherever these
pages cite a file under ``tests/``, ``docs/`` or ``notebooks/``, it is in such a checkout (or
on GitHub), not in the installed package.

Either way, one command confirms it worked:

.. code-block:: bash

   magnus --flavors 3 --environment vacuum --energy 1 --energy-unit GeV \
       --baseline 1300 --baseline-unit km

If you would rather not install the package at all, put ``src/`` on your Python
path instead.  You still need its five runtime dependencies:

.. code-block:: bash

   pip install -r src/requirements.txt

Then point at ``src/`` itself.  Every module imports through the ``magnus``
package (``import magnus.globaldefs``), so the package directory does not
belong on the path:

.. code-block:: python

   import sys
   sys.path.insert(0, 'src')

   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

Verifying the installation
---------------------------

From a clone, the test suite confirms that everything is configured correctly for your
system:

.. code-block:: bash

   pip install -e '.[test]'
   pytest tests/ -n auto

The suite has about 2400 tests and takes a few minutes spread over all cores.
``-n auto`` comes from ``pytest-xdist``, which the ``test`` extra installs.  CI runs
the same suite on Python 3.10 to 3.13 on every push, and the badge on :doc:`index`
shows the result.

**What passing means.**  The suite checks the package against closed forms, an
independent ODE solver, an independently generated expansion, and properties that must
hold exactly; :ref:`what-accuracy-means` lists each check and its result.  It also
fixes the conventions that would give self-consistent but wrong answers if changed:
the slab ordering, the sign of the antineutrino potential, the mass ordering and the
channel indexing (:ref:`conventions`).  Finally,
``tests/test_documented_examples.py`` runs the code blocks in ``README.md`` and the
:doc:`quickstart`.  The ``jupyter-execute`` examples in these pages are run by the
documentation build rather than by pytest, and the notebooks by their own CI job.

Some tests are skipped, and that is expected.  Tests that need an optional tool skip
when it is absent.  In an unpacked sdist, the tests of the repository rather than the
library (the documentation, the notebooks, the paper's assets and the CI workflows) skip
too, because the sdist does not ship those files; each states this in its skip reason.

Measuring test coverage
~~~~~~~~~~~~~~~~~~~~~~~~

The ``test`` extra also installs ``pytest-cov``, so the same suite can report
which lines and branches of the package it exercises:

.. code-block:: bash

   pytest tests/ --cov --cov-report=term-missing

The settings, including branch coverage, are in ``[tool.coverage.run]`` in
``pyproject.toml``, so a bare ``--cov`` measures what CI measures.
``--cov-report=html`` writes a browsable ``htmlcov/`` tree, which shows which branch
of a function is untested.  Branch coverage is on because ``oscprob.py`` consists
mostly of thin wrappers that one parametrized test covers, so line coverage alone
would overstate how well the dispatch, the refinement caps and the warnings are
tested.

The run fails below **90%**, a few points under what the suite measures.  The floor
applies to every coverage run, so measuring a single test file fails it; pass
``--cov-fail-under=0`` to measure part of the suite:

.. code-block:: bash

   pytest tests/test_cli.py --cov --cov-report=term-missing --cov-fail-under=0

Coverage roughly doubles the run time (188 s to 394 s on one machine), because
most of the time is spent in a per-slab Python loop.  Run it when you need the
number, not on every change.

File tree
---------

The top level of the repository, and the package under ``src/``.  The complete
listing, with a comment on every tracked file, is ``TREE`` in
``tests/test_file_tree.py``, which generates this block.

.. code-block:: text

   Magnus/
   ├── .github/                        # GitHub Actions workflows: tests, lint, notebooks, docs, publishing
   ├── .gitignore                      # Build, cache and generated-output artifacts
   ├── CHANGELOG.md                    # Version history (Keep a Changelog format)
   ├── CITATION.cff                    # Machine-readable citation metadata; drives GitHub's "Cite this repository"
   ├── LICENSE                         # GNU GPL v3 (GPL-3.0-only), the full license text
   ├── MANIFEST.in                     # Adds tests/conftest.py to the sdist, so that tests needing a source checkout are skipped there
   ├── README.md                       # Project overview; also the PyPI project description
   ├── docs/                           # Sphinx documentation configuration and source
   ├── fig/                            # Plots produced by the example notebooks
   ├── img/                            # Figures used by the documentation
   ├── notebooks/                      # Numbered Jupyter notebooks -- see docs/source/tutorials.rst
   ├── pyproject.toml                  # Build system, dependencies, and the `magnus` console-script entry point
   ├── resources/                      # Paper sources and benchmark drivers, kept in the repository but not distributed
   ├── tools/                          # Standalone utilities that are not part of the package
   ├── src/                            # The package itself -- the only thing a `pip install` delivers
   │   ├── magnus/                     # Main Python package
   │   │   ├── __init__.py             # Imports the twelve public modules and expmkernels, and exposes __version__
   │   │   ├── __main__.py             # Entry point for `python -m magnus`
   │   │   ├── _validate.py            # The argument checks every public entry point applies once per call
   │   │   ├── adiabatic.py            # Adiabatic transport + Magnus-patch hybrid strategy (strategy='hybrid'/'auto')
   │   │   ├── authors.py              # Package author string (internal; not part of the public API)
   │   │   ├── avgprob.py              # The phase average over an energy spread, and the decohered limit
   │   │   ├── cli.py                  # `magnus` command-line calculator (also `python -m magnus`)
   │   │   ├── data/                   # Package data, installed with the code
   │   │   │   └── solar_models/       # Twelve standard solar models: three columns each, with provenance
   │   │   ├── earth.py                # PREM density profile, chord/zenith-angle geometry
   │   │   ├── expansionterms.py       # Generates the Omega_k terms symbolically, to any order
   │   │   ├── expmkernels.py          # Compiled Cayley-Hamilton matrix exponential for 2x2/3x3 (the numba backend)
   │   │   ├── globaldefs.py           # Units, physical constants, NuFIT parameter sets
   │   │   ├── hamiltonians/           # 2nu-5nu Hamiltonians: vacuum, matter, NSI, LIV (the only subpackage)
   │   │   │   ├── __init__.py         # Explicit named imports from the four hamiltonians{2,3,4,5}nu.py modules
   │   │   │   ├── _angles.py          # Interprets the four angle conventions; rejects an out-of-range sine
   │   │   │   ├── _broadcast.py       # Lets every builder take an array of its one scalar argument
   │   │   │   ├── hamiltonians2nu.py
   │   │   │   ├── hamiltonians3nu.py
   │   │   │   ├── hamiltonians4nu.py
   │   │   │   ├── hamiltonians5nu.py
   │   │   │   └── hamiltonians_pseudodirac.py  # Pseudo-Dirac spectra: per-mass-state pairing, and the sterile partners
   │   │   ├── magnus.py               # Magnus-expansion numerical core: term recursion, GL integrators, batched kernel
   │   │   ├── matter.py               # Density profiles, electron number density, CC potential
   │   │   ├── oscprob.py              # osc_prob and every physics-scenario wrapper (main API)
   │   │   ├── oscprobstd.py           # Closed-form 2nu/3nu probabilities (used to validate the wrapper API)
   │   │   ├── plotting.py             # Pre-packaged plotting tools: one call per standard figure
   │   │   ├── py.typed                # PEP 561 marker: tells type checkers the package is typed
   │   │   ├── solarmodels.py          # Tabulated standard solar models, as profiles for the Sun wrappers
   │   │   └── version.py              # Resolves the version from pyproject.toml (internal)
   │   └── requirements.txt            # The five runtime dependencies: numpy, scipy, joblib, matplotlib, numba
   └── tests/                          # Test suite (pytest; runs in CI)
