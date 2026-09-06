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
   The distribution is published as **magnuspy**, but the import package is
   **magnus** -- so you ``pip install magnuspy`` and then ``import magnus``.
   The two names are independent in Python packaging, and they differ here
   only because ``magnus`` was already taken on PyPI by an unrelated project.
   The command-line tool is ``magnus`` as well.

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

   pip install -e '.[test]'
   pytest tests/ -v

It is about 1200 tests and takes some fifteen minutes with ``-n auto``; the same suite runs in CI on
Python 3.10-3.13 on every push, so the badge on the :doc:`index` page tells
you whether it passes there.

**What passing means.**  The suite is not only a smoke test, so it is worth
knowing what it establishes:

* **Against closed forms.**  Two- and three-flavor vacuum probabilities, and
  two-flavor constant-density matter, for neutrinos and antineutrinos, to
  machine precision.
* **Against an independent integrator.**  Asymmetric profiles with complex
  Hamiltonians and full PREM Earth crossings, scored against
  ``scipy.integrate.solve_ivp``/DOP853 at ``rtol=1e-12``.
* **Against an independently coded recursion.**  The Magnus terms at orders
  1--6, and the Gauss--Legendre convergence rates under slab halving (error
  ratios 4, 16, 64).
* **Properties that must hold exactly.**  Unitarity, and two *bit-identity*
  assertions rather than tolerances: the energy-batched scan against the
  per-point path, and ``n_jobs > 1`` against serial.  An optimization that
  changed an answer fails those rather than passing quietly.
* **Conventions.**  Slab ordering, the antineutrino potential sign, the mass
  ordering and the channel indexing -- each of which has been wrong here at
  some point, and each of which is self-consistent when wrong.  See
  :ref:`conventions`.
* **The documentation.**  Every ``jupyter-execute`` block in the docstrings is
  run when the docs are built, so an example that no longer works fails the
  build rather than misleading a reader; the notebooks are executed by their
  own CI job for the same reason.

Skips are expected rather than a sign of trouble: tests that need an optional
tool stand down when it is absent.

Measuring test coverage
~~~~~~~~~~~~~~~~~~~~~~~~

The ``test`` extra also installs ``pytest-cov``, so the same suite can report
which lines and branches of the package it exercises:

.. code-block:: bash

   pytest tests/ --cov --cov-report=term-missing

What to measure -- the source tree, the omitted files, and branch coverage --
is configured once in ``[tool.coverage.run]`` in ``pyproject.toml``, so a bare
``--cov`` here measures exactly what CI measures.  ``--cov-report=html`` writes
a browsable ``htmlcov/`` tree instead, which is the more useful form when the
question is *which* branch of a particular function is untested.

Branch coverage is on deliberately.  A plain line-coverage figure overstates
how well this package is tested: ``oscprob.py`` is dominated by thin wrappers
that one parametrized test sweeps in a single pass, so the number to read is
whether the dispatch chain, the refinement caps and the warning paths are each
taken in *both* directions.

The run fails below **90%**, which is a floor rather than a target: the suite
measures 93%, and the three points of headroom keep the check from tripping on
the fraction of a percent that moves between interpreters while still catching
a module added without tests or a test file deleted.  The floor is in
``pyproject.toml``, so it applies to *every* coverage run -- measuring a single
test file therefore reports far below 90 and exits non-zero.  That is expected;
pass ``--cov-fail-under=0`` when deliberately measuring part of the suite:

.. code-block:: bash

   pytest tests/test_cli.py --cov --cov-report=term-missing --cov-fail-under=0

Instrumentation is expensive for this suite: measured on one machine, the run
goes from 188 s to 394 s, a factor of 2.1.  That is more than the usual
coverage overhead, and it is what one would expect here, since the cost is
dominated by a per-slab Python loop rather than by time spent inside numpy.
Run it when you want the number, not on every iteration.

File Tree
---------

.. code-block:: text

   Magnus/
   ├── .github/                        # GitHub Actions workflows: tests, lint, notebooks, docs, publishing
   │   └── workflows/
   │       ├── lint.yml                # Ruff lint (blocking) + CLI-reference drift check
   │       ├── notebooks.yml           # Executes every notebook; paths-filtered, so docs-only changes skip it
   │       ├── pages.yml               # GitHub Pages deployment for the Sphinx documentation
   │       ├── publish.yml             # PyPI (OIDC) automated publishing workflow, on GitHub Release
   │       └── tests.yml               # GitHub Actions CI testing pipeline (Python 3.10-3.13) + coverage
   ├── .gitignore                      # Build, cache and generated-output artifacts
   ├── CHANGELOG.md                    # Version history (Keep a Changelog format)
   ├── CITATION.cff                    # Machine-readable citation metadata; drives GitHub's "Cite this repository"
   ├── LICENSE                         # GNU GPL v3 (GPL-3.0-only), the full license text
   ├── README.md                       # This file
   ├── docs/                           # Sphinx documentation configuration and source
   │   ├── Makefile                    # Build commands for Unix
   │   ├── check_doc_snippets.py       # Checks the code snippets quoted in the prose pages still run
   │   ├── dev/
   │   ├── make.bat                    # Build commands for Windows
   │   ├── make_figures.py             # Regenerates the data-driven SVG in source/_static/
   │   ├── regen_cli_help.py           # Regenerates the --help block quoted in source/cli.rst
   │   ├── requirements.txt            # Sphinx + theme + extensions needed to build the docs
   │   └── source/
   │       ├── _static/
   │       │   ├── adiabatic_avoided_crossing.svg  # Hand-authored: adiabatic against diabatic at a crossing
   │       │   ├── api_layers.svg      # Generated by docs/make_figures.py: the three layers of oscprob
   │       │   ├── adiabatic_segmentation.svg  # Hand-authored: adiabatic / patch / adiabatic along the ray
   │       │   ├── adiabatic_speedup.svg  # Generated by docs/make_figures.py from the measured grid
   │       │   ├── averaging_regimes.svg  # Hand-authored: when averaging removes an error and when it does not
   │       │   ├── call_sequence.svg   # Generated by docs/make_figures.py: what is built, in what order, per call
   │       │   ├── magnus_logo.png     # Sidebar logo
   │       │   └── module_layout.svg   # Generated by docs/make_figures.py: the real internal import graph
   │       ├── adiabatic_strategy.rst  # The adiabatic + Magnus hybrid strategy: derivation, diagrams, validation
   │       ├── api_reference.rst       # Wraps the autoapi-generated module pages
   │       ├── architecture.rst        # The wrapper/middle/primordial layering, with diagrams
   │       ├── averaged_probability.rst  # Phase-averaged probabilities: derivation, diagram, validation
   │       ├── changelog.rst           # Renders the root CHANGELOG.md via myst-parser
   │       ├── citing.rst              # How to cite the software, and what to state in the text
   │       ├── cli.rst                 # Command-line calculator: flag reference and examples
   │       ├── comparison.rst          # Against NuOscProbExact and nuSQuIDS: where each method wins, from notebook 25
   │       ├── conf.py                 # Sphinx build configuration (autoapi + napoleon + bibtex + mermaid + myst)
   │       ├── diagnostics.rst         # What rtol really controls, what each safeguard cannot do, every warning
   │       ├── engines.rst             # Which engine answers a call, and how the dispatch order is decided
   │       ├── expansion_terms.rst     # The Omega_k terms to any order, and how they are generated
   │       ├── functions.rst           # Full osc_prob_{2,3,4,5}nu_* listing, grouped by environment/scenario
   │       ├── index.rst               # Master documentation page: overview, features, when Magnus wins
   │       ├── installation.rst        # Requirements, install instructions, file tree
   │       ├── methodology.rst         # The Magnus expansion, integrators, and performance engineering
   │       ├── performance.rst         # Where the time goes, and the population every tuned constant was measured on
   │       ├── plotting.rst            # The pre-packaged plotting tools
   │       ├── quickstart.rst          # Worked Python-API code examples for every entry point
   │       ├── recipes.rst             # What Magnus can compute, with the code -- executed at build time
   │       ├── references.rst          # Bibliography page rendering
   │       ├── refs.bib                # BibTeX citations for the Magnus-expansion and PREM literature
   │       └── tutorials.rst           # Guide to the numbered example notebooks in notebooks/
   ├── fig/                            # Plots produced by the example notebooks
   ├── img/                            # Figures used by the documentation
   │   ├── anim_cp.gif                 # Animated: the CP phase running through 2 pi
   │   ├── anim_earth.gif              # Animated: a chord swinging to a detector at the South Pole
   │   ├── anim_shock.gif              # Animated: a supernova shock front sweeping outward
   │   ├── anim_slabs.gif              # Animated: a profile cut into more and more slabs
   │   ├── anim_solar_nsi.gif          # Animated: the Sun, with a non-standard interaction dialed up
   │   ├── anim_sterile.gif            # Animated: a sterile state as its mass splitting grows
   │   ├── anim_wave.gif               # Animated: a density crest traveling along the baseline
   │   └── gallery/                    # Figures lifted from the executed notebooks, embedded in the docs
   ├── notebooks/                      # Numbered Jupyter notebooks -- see docs/source/tutorials.rst
   │   ├── 01_magnus_introduction.ipynb  # The shortest path to a probability
   │   ├── 02_magnus_2nu_vacuum_matter.ipynb  # Two flavors, across seven matter profiles
   │   ├── 03_magnus_3nu_vacuum_matter.ipynb  # The same, with three flavors and a CP phase
   │   ├── 04_magnus_long_baseline.ipynb  # Between two points on the surface
   │   ├── 05_magnus_biprobability.ipynb  # The CP ellipse
   │   ├── 06_magnus_oscillograms.ipynb  # Zenith angle against energy, in one call
   │   ├── 07_magnus_bsm_sterile_nu.ipynb  # Four and five flavors
   │   ├── 08_magnus_bsm_nsi.ipynb     # Non-standard interactions
   │   ├── 09_magnus_bsm_liv.ipynb     # Lorentz-invariance violation
   │   ├── 10_magnus_averaged_probability.ipynb  # What survives when the phase is unresolvable
   │   ├── 11_magnus_matrix_exponential.ipynb  # How exp(Omega) is actually built
   │   ├── 12_magnus_adiabatic_hybrid_strategy.ipynb  # 'auto' against 'magnus', timed against solve_ivp
   │   ├── 13_magnus_tabulated_solar_model.ipynb  # A real BS05 profile: an error that is a phase
   │   ├── 14_magnus_supernova_shock.ipynb  # A shock front: an error that is an envelope
   │   ├── 15_magnus_antineutrinos.ipynb  # Conjugate and flip, and two ways to get it half right
   │   ├── 16_magnus_exact_vs_approximations.ipynb  # Where the textbook formulas are exact, and where the substitution breaks
   │   ├── 17_magnus_ordering_and_octant.ipynb  # The sign of D31, and how large the two open questions are
   │   ├── 18_magnus_unusual_density_profiles.ipynb  # Arrangement beats the mean, except for one exact symmetry
   │   ├── 19_magnus_custom_hamiltonian.ipynb  # The H_func contract, and the vectorization trick
   │   ├── 20_magnus_numerical_edge_cases.ipynb  # Degeneracies that return numbers, and the nine warnings
   │   ├── 21_magnus_what_tolerance_means.ipynb  # rtol is a stopping criterion, not an error bound
   │   ├── 22_magnus_which_engine_answered.ipynb  # strategy_info, and an error bar with no oracle
   │   ├── 23_magnus_when_averaging_helps.ipynb  # Phase error falls away, envelope error does not
   │   ├── 24_magnus_performance.ipynb  # What is worth doing, and when each trick is worth nothing
   │   ├── 25_magnus_against_other_codes.ipynb  # Where a closed form wins, and a conventions trap that looks like accuracy
   │   ├── 26_magnus_nufit_evolution.ipynb  # How the NuFIT likelihood, not just the best fit, moves the probability
   │   ├── 27_magnus_animations.ipynb  # Ten sweeps as filmstrips; RENDER = True writes them as GIFs
   │   ├── 28_magnus_paper_figures.ipynb  # Every figure in the CPC article, in one run
   │   ├── 29_magnus_pseudo_dirac.ipynb  # Tiny splittings, coherent blocks, and where the effect is invisible
   │   ├── README.md                   # This file
   │   ├── make_notebooks.py           # BUILDS the notebooks above -- edit this, not the .ipynb
   │   ├── external_speed_accuracy.json  # Five external codes' speed and accuracy (NuOscProbExact project)
   │   ├── external_prem_speed_accuracy.json  # Notebook 25 section 5: the same, on a PREM chord, both codes batched
   │   ├── external_speed_accuracy_const.json  # Figure 12, top panel: constant density, seven codes plus Magnus
   │   ├── external_earth_plane.json   # Figure 12, middle panel: a PREM chord at three flavors
   │   ├── external_prem_speed_accuracy_new.json  # Figure 12, bottom panel: the same chord at 3+1
   │   ├── magnus_own_reference.json   # Magnus's own 50-digit references, in its own conventions, on those three grids
   │   ├── external_profile_benchmarks.json  # Notebook 25 section 9: smooth-profile speed/accuracy, all codes on one machine
   │   ├── external_shock_benchmarks.json  # Notebook 25 section 11: the supernova shock, both front widths
   │   ├── external_shock_4nu.json     # Notebook 25 section 12: the same shock at 3+1
   │   ├── external_shock_nsi.json     # Notebook 25 section 13: the same shock with NSI
   │   ├── external_solar_nusquids.json  # Notebook 25 section 10: nuSQuIDS's energy-averaged solar survival probability
   │   ├── gen_profile_benchmarks.py   # GENERATES external_profile_benchmarks.json -- needs the external codes
   │   ├── gen_mp_reference.py         # GENERATES mp_reference_profile.json -- the mpmath referee for Figure 11
   │   ├── mp_reference_profile.json   # Triple-Richardson mpmath reference, exponential profile, 2-5 flavors
   │   ├── rescore_against_mp_reference.py  # RE-SCORES external_profile_benchmarks.json against it; timings untouched
   │   ├── append_order_series.py      # ADDS the order-6 and order-8 Magnus series to that file
   │   ├── probe_commensurability.py   # Asks whether a timing taken today is comparable with the stored ones
   │   ├── gen_solar_average_cost.py   # GENERATES external_solar_average_cost.json -- cost per configuration
   │   ├── external_solar_average_cost.json  # Averaged-probability cost on BS2005-AGS,OP, eight configurations
   │   ├── sterile_projector_check.py  # Reproduces the sterile projector defect and its fix, three arms, one command
   │   ├── retime_magnus_series.py     # RE-TIMES both codes in Figure 11; references and grids untouched
   │   ├── prem_chord_common.py        # The PREM chord at cos(theta_z) = -0.9, shared by the two scripts below
   │   ├── gen_prem_reference.py       # GENERATES prem_chord_reference.json -- segment-aligned, layer edges respected
   │   ├── prem_chord_reference.json   # That reference; PARTIAL, 4nu stops at 6 of 12 energies and 5nu is unstarted
   │   ├── gen_prem_benchmarks.py      # GENERATES external_prem_chord_benchmarks.json -- the Earth analogue of Fig. 11
   │   ├── external_prem_chord_benchmarks.json  # That file: both codes on one requested tolerance, Earth chord, 2-5 flavors
   │   ├── append_npe_rtol_series.py   # ADDS a tolerance-dialled NuOscProbExact series to the smooth-profile file
   │   ├── append_npe_rtol_prem.py     # The same for the Earth chord, via earth_slabs and the librarys own refinement
   │   ├── gen_shock_benchmarks.py     # GENERATES external_shock_benchmarks.json -- runs notebook 14s own cells
   │   ├── gen_shock_4nu.py            # GENERATES external_shock_4nu.json -- the shock at 3+1, own DOP853 referee
   │   ├── gen_shock_nsi.py            # GENERATES external_shock_nsi.json -- the shock with NSI, own DOP853 referee
   │   ├── gen_solar_nusquids.py       # GENERATES external_solar_nusquids.json -- needs nuSQuIDS
   │   ├── make_nufit_chi2.py          # Extracts notebook 26's NuFIT chi^2 profiles
   │   ├── make_shock_reference.py     # Freezes notebook 14's solve_ivp oracle
   │   ├── matplotlibrc                # Shared plot styling for the notebooks
   │   ├── paper_figure_cache.json     # Every paper-figure input that depends on the configuration and not on the run: reference probabilities, order curves, and timings
   │   ├── nufit_chi2.json             # Those profiles, v2.0-v6.1 (NuFIT collaboration)
   │   └── shock_reference.json        # That oracle, as exact hex floats
   ├── pyproject.toml                  # Build system, dependencies, and the `magnus` console-script entry point
   ├── resources/                      # Travels with the code; reaches neither the wheel nor the sdist
   │   ├── benchmarks/                 # The cross-code benchmark harness and its frozen artifacts, copied from NuOscProbExact so its measurements can be reproduced here
   │   └── paper/                      # The Computer Physics Communications article documenting this package
   │       ├── README.md               # How to build the paper, and where each of its numbers comes from
   │       ├── HANDOVER-pseudodirac.md  # Brief for adding pseudo-Dirac neutrinos to the library
   │       ├── API-pseudodirac.md      # The pseudo-Dirac API, summarized for the session writing the panel
   │       ├── audit-criteria.md       # What the manuscript audit checks
   │       ├── pending-edits.md        # Edits and re-runs the manuscript still owes, with what each one moves
   │       ├── review-crossread.md     # A cross-read of the manuscript against the code
   │       ├── main.tex                # The paper -- ordinary LaTeX; a revision diff is mechanical
   │       ├── refs.bib                # NuOscProbExact's bibliography, with the Magnus entries appended below a separator
   │       ├── elsarticle.cls          # Bundled, so the folder compiles without the Elsevier bundle
   │       ├── elsarticle-num.bst
   │       └── figs/                   # Its eight figures, written by notebook 28
   ├── tools/                          # Standalone utilities that are not part of the package
   │   └── make_demo_video.py          # Joins and shrinks notebook 27's clips; shared with NuOscProbExact
   ├── src/                            # The package itself -- the only thing a `pip install` delivers
   │   ├── magnus/                     # Main Python package
   │   │   ├── __init__.py             # Explicit named imports from the four hamiltonians{2,3,4,5}nu.py modules
   │   │   ├── __main__.py             # Entry point for `python -m magnus`
   │   │   ├── adiabatic.py            # Adiabatic transport + Magnus-patch hybrid strategy (strategy='hybrid'/'auto')
   │   │   ├── authors.py              # Package author string (internal; not part of the public API)
   │   │   ├── avgprob.py              # Phase-averaged (decohered) probabilities
   │   │   ├── cli.py                  # `magnus` command-line calculator (also `python -m magnus`)
   │   │   ├── earth.py                # PREM density profile, chord/zenith-angle geometry
   │   │   ├── expansionterms.py       # Generates the Omega_k terms symbolically, to any order
   │   │   ├── expmkernels.py          # Compiled Cayley-Hamilton matrix exponential for 2x2/3x3 (the numba backend)
   │   │   ├── globaldefs.py           # Units, physical constants, NuFit parameter sets
   │   │   ├── hamiltonians/           # 2nu-5nu Hamiltonians: vacuum, matter, NSI, LIV (the one true subpackage)
   │   │   │   ├── __init__.py         # Explicit named imports from the four hamiltonians{2,3,4,5}nu.py modules
   │   │   │   ├── _angles.py          # Interprets the four angles conventions; rejects an out-of-range sine
   │   │   │   ├── hamiltonians2nu.py
   │   │   │   ├── hamiltonians3nu.py
   │   │   │   ├── hamiltonians4nu.py
   │   │   │   ├── hamiltonians5nu.py
   │   │   │   └── hamiltonians_pseudodirac.py  # Pseudo-Dirac spectra: per-mass-state pairing, and the sterile partners
   │   │   ├── magnus.py               # Magnus-expansion numerical core: term recursion, GL integrators, batched kernel
   │   │   ├── matter.py               # Density profiles, electron number density, CC potential
   │   │   ├── oscprob.py              # osc_prob and every physics-scenario wrapper (main API)
   │   │   ├── oscprobstd.py           # Closed-form 2nu/3nu probabilities (used to validate the wrapper API)
   │   │   ├── plotting.py             # Pre-packaged plotting tools: one call instead of thirty lines
   │   │   ├── py.typed                # PEP 561 marker: tells type checkers the annotations are real
   │   │   └── version.py              # Resolves the version from pyproject.toml (internal)
   │   └── requirements.txt            # Sphinx + theme + extensions needed to build the docs
   └── tests/                          # Test suite (pytest; runs in CI)
       ├── conftest.py                 # Path setup so magnus is importable without installation
       ├── test_adiabatic.py           # Adiabatic + Magnus hybrid strategy: detection, merging, ODE cross-checks
       ├── test_angles.py              # The four `angles` conventions and the guards between them
       ├── test_avgprob.py             # Phase-averaged probabilities
       ├── test_cli.py                 # magnus command-line calculator
       ├── test_pseudodirac.py         # Pseudo-Dirac Hamiltonians: the Dirac limit, blocks, and the factor of two
       ├── test_documented_examples.py  # Runs the code blocks in README.md and quickstart.rst
       ├── test_earth_matter.py        # PREM profile, chord geometry, electron density
       ├── test_engines.py             # Which engine answers, and the cross-checks between them
       ├── test_expansionterms.py      # The symbolic term generator against the hand-written orders
       ├── test_expm_backend.py        # The two matrix-exponential backends, their switch, and degeneracies
       ├── test_fuzz_statistics.py     # Randomized profiles, scored in bulk
       ├── test_file_tree.py           # This file: generates the tree above and checks it against git
       ├── test_globaldefs.py          # NuFit historical parameter dict/loader
       ├── test_hamiltonians.py        # Hamiltonian/mixing-matrix builders
       ├── test_invariants.py          # Properties that must hold across the whole engine matrix
       ├── test_magnus_expansion.py    # Magnus-core correctness (terms, orders, GL rates, unitarity)
       ├── test_oscprob.py             # Oscillation-probability engine, closed-form and ODE cross-checks
       ├── test_palindrome.py          # The palindromic-profile optimization and its gate
       ├── test_plotting.py            # Pre-packaged plotting tools: house-style defaults, layouts
       ├── test_routine_listings.py    # Each module's Routine listings names every public function it defines
       ├── test_tolerance.py           # What rtol/atol promise, and the effective-refinement gate
       ├── test_validation.py          # Input-validation guards and their error messages
       └── test_version.py             # Version resolution from pyproject.toml / installed metadata
