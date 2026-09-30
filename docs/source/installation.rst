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
   The distribution is published as **magnuspy**, but the import package is
   **magnus** -- so you ``pip install magnuspy`` and then ``import magnus``.
   The two names are independent in Python packaging, and they differ here
   only because ``magnus`` was already taken on PyPI by an unrelated project.
   The command-line tool is ``magnus`` as well.

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

It is about 1900 tests and takes a few minutes spread over all cores (``-n auto``,
from ``pytest-xdist``, which the ``test`` extra installs); the same suite
runs in CI on Python 3.10-3.13 on every push, so the badge on the :doc:`index`
page tells you whether it passes there.

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
* **Properties that must hold exactly.**  Unitarity, and *bit-identity*
  assertions rather than tolerances: repeated calls, and a baseline scan given in
  shuffled order, return the same answer to the last bit, so an optimization that
  changed an answer fails rather than passing quietly.  A parallel scan is held to
  the tolerance against a serial one, not to the last bit.  The
  energy-batched scan is held to 1e-12 against the per-point path, with the
  grid and tolerances pinned so that the two are arithmetically the same
  problem.
* **Conventions.**  Slab ordering, the antineutrino potential sign, the mass
  ordering and the channel indexing -- each of which has been wrong here at
  some point, and each of which is self-consistent when wrong.  See
  :ref:`conventions`.
* **The documentation.**  Every ``jupyter-execute`` block in the docstrings is
  run when the docs are built, so an example that no longer works fails the
  build rather than misleading a reader; the notebooks are executed by their
  own CI job for the same reason.

Skips are expected rather than a sign of trouble: tests that need an optional
tool stand down when it is absent.  Run from an unpacked sdist, the tests that check
the repository rather than the library -- the documentation, the notebooks, the
paper's assets, the CI workflows -- skip as well, since the sdist does not ship those
files; each says so in its skip reason.

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

File tree
---------

.. code-block:: text

   Magnus/
   ├── .github/                        # GitHub Actions workflows: tests, lint, notebooks, docs, publishing
   │   └── workflows/
   │       ├── lint.yml                # Ruff lint (blocking) + CLI-reference drift check
   │       ├── notebooks.yml           # Executes every notebook; paths-filtered, so docs-only changes skip it
   │       ├── pages.yml               # GitHub Pages deployment for the Sphinx documentation
   │       ├── publish.yml             # PyPI (OIDC) automated publishing workflow, on GitHub Release
   │       └── tests.yml               # GitHub Actions CI testing pipeline (Python 3.10-3.13), coverage, and the sdist's own tests
   ├── .gitignore                      # Build, cache and generated-output artifacts
   ├── CHANGELOG.md                    # Version history (Keep a Changelog format)
   ├── CITATION.cff                    # Machine-readable citation metadata; drives GitHub's "Cite this repository"
   ├── LICENSE                         # GNU GPL v3 (GPL-3.0-only), the full license text
   ├── MANIFEST.in                     # Adds tests/conftest.py to the sdist, which skips the checkout-only tests there
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
   │       │   ├── adiabatic_segmentation.svg  # Hand-authored: adiabatic / patch / adiabatic along the ray
   │       │   ├── adiabatic_speedup.svg  # Generated by docs/make_figures.py from the measured grid
   │       │   ├── averaging_regimes.svg  # Hand-authored: when averaging removes an error and when it does not
   │       │   ├── call_sequence.svg   # Generated by docs/make_figures.py: what is built, in what order, per call
   │       │   ├── magnus.css          # Places equation numbers beside display math
   │       │   ├── magnus_logo.png     # Sidebar logo
   │       │   └── module_layout.svg   # Generated by docs/make_figures.py: the real internal import graph
   │       ├── adiabatic_strategy.rst  # The adiabatic + Magnus hybrid strategy: derivation, diagrams, validation
   │       ├── api_reference.rst       # Wraps the autoapi-generated module pages
   │       ├── architecture.rst        # Modules, the path of a request, and the four layers of oscprob
   │       ├── averaged_probability.rst  # Phase-averaged probabilities: derivation, diagram, validation
   │       ├── changelog.rst           # Renders the root CHANGELOG.md via myst-parser
   │       ├── citing.rst              # How to cite the software, and what to state in the text
   │       ├── cli.rst                 # Command-line calculator: flag reference and examples
   │       ├── comparison.rst          # Against other codes: the paper's eight-code comparison, then NuOscProbExact and nuSQuIDS in detail
   │       ├── conf.py                 # Sphinx build configuration (autoapi + napoleon + bibtex + mermaid + myst)
   │       ├── diagnostics.rst         # What rtol really controls, what each safeguard cannot do, every warning
   │       ├── engines.rst             # Which engine answers a call, and how the dispatch order is decided
   │       ├── examples.rst            # Usage and examples: index of the pages in examples/
   │       ├── examples/               # Usage and examples, one page per section of Sec. 6 of the paper
   │       ├── expansion_terms.rst     # The Omega_k terms to any order, and how they are generated
   │       ├── conventions.rst         # Flavor order, signs, mass ordering, parameters and units, in one place
   │       ├── functions.rst           # Full osc_prob_{2,3,4,5}nu_* listing, grouped by environment/scenario
   │       ├── index.rst               # Master documentation page: overview, features, when Magnus wins
   │       ├── installation.rst        # Requirements, install instructions, file tree
   │       ├── methodology.rst         # The Magnus expansion, integrators, and performance engineering
   │       ├── performance.rst         # Where the time goes, and the population every tuned constant was measured on
   │       ├── plotting.rst            # The pre-packaged plotting tools
   │       ├── quickstart.rst          # Worked Python-API code examples for every entry point
   │       ├── recipes.rst             # What Magnus can compute, with the code -- executed at build time
   │       ├── references.rst          # Bibliography page rendering
   │       ├── refs.bib                # BibTeX citations for the Magnus-expansion, PREM and solar-model literature
   │       ├── solar_models.rst        # The twelve tabulated standard solar models, and how the Sun wrappers use them
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
   │   ├── gallery/                    # Figures lifted from the executed notebooks, embedded in the docs
   │   └── paper/                      # Figures from the paper, as SVG, embedded in the docs
   ├── notebooks/                      # Numbered Jupyter notebooks -- see docs/source/tutorials.rst
   │   ├── 01_magnus_introduction.ipynb  # The shortest path to a probability
   │   ├── 02_magnus_2nu_vacuum_matter.ipynb  # Two flavors, across seven matter profiles
   │   ├── 03_magnus_3nu_vacuum_matter.ipynb  # The same, with three flavors and a CP phase
   │   ├── 04_magnus_long_baseline.ipynb  # Between two points on the surface
   │   ├── 05_magnus_biprobability.ipynb  # The CP ellipse
   │   ├── 06_magnus_oscillograms.ipynb  # Zenith angle against energy, one energy scan per angle
   │   ├── 07_magnus_bsm_sterile_nu.ipynb  # Four and five flavors
   │   ├── 08_magnus_bsm_nsi.ipynb     # Non-standard interactions
   │   ├── 09_magnus_bsm_liv.ipynb     # Lorentz-invariance violation
   │   ├── 10_magnus_averaged_probability.ipynb  # What survives when the phase is unresolvable
   │   ├── 11_magnus_matrix_exponential.ipynb  # How exp(Omega) is actually built
   │   ├── 12_magnus_adiabatic_hybrid_strategy.ipynb  # 'auto' against 'magnus', timed against solve_ivp
   │   ├── 13_magnus_tabulated_solar_model.ipynb  # The twelve solar models by name, and the observable an experiment measures
   │   ├── 14_magnus_supernova_shock.ipynb  # A shock front: an error that is an envelope
   │   ├── 15_magnus_antineutrinos.ipynb  # Conjugate and flip, and two ways to get it half right
   │   ├── 16_magnus_exact_vs_approximations.ipynb  # Where the textbook formulas are exact, and where the substitution breaks
   │   ├── 17_magnus_ordering_and_octant.ipynb  # The sign of D31, and how large the two open questions are
   │   ├── 18_magnus_unusual_density_profiles.ipynb  # Arrangement beats the mean, except for one exact symmetry
   │   ├── 19_magnus_custom_hamiltonian.ipynb  # The H_func contract, and the vectorization trick
   │   ├── 20_magnus_numerical_edge_cases.ipynb  # Degeneracies that return numbers, and the fifteen warnings
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
   │   ├── external_speed_accuracy_const.json  # Figure 15, top panel: constant density, seven codes plus Magnus
   │   ├── external_earth_plane.json   # Figure 15, middle panel: a PREM chord at three flavors
   │   ├── external_prem_speed_accuracy_new.json  # Figure 15, bottom panel: the same chord at 3+1
   │   ├── magnus_own_reference.json   # Magnus's own 50-digit references, in its own conventions, on those three grids
   │   ├── external_profile_benchmarks.json  # Notebook 25 section 9: smooth-profile speed/accuracy, all codes on one machine
   │   ├── external_shock_benchmarks.json  # Notebook 25 section 11: the supernova shock, both front widths
   │   ├── external_shock_4nu.json     # Notebook 25 section 12: the same shock at 3+1
   │   ├── external_shock_nsi.json     # Notebook 25 section 13: the same shock with NSI
   │   ├── external_solar_nusquids.json  # Notebook 25 section 10: nuSQuIDS's energy-averaged solar survival probability
   │   ├── gen_profile_benchmarks.py   # GENERATES external_profile_benchmarks.json -- needs the external codes
   │   ├── gen_mp_reference.py         # GENERATES mp_reference_profile.json -- the mpmath referee for Figure 12
   │   ├── mp_reference_profile.json   # Triple-Richardson mpmath reference, exponential profile, 2-5 flavors
   │   ├── rescore_against_mp_reference.py  # RE-SCORES external_profile_benchmarks.json against it; timings untouched
   │   ├── append_order_series.py      # ADDS the order-6 and order-8 Magnus series to that file
   │   ├── probe_commensurability.py   # Asks whether a timing taken today is comparable with the stored ones
   │   ├── gen_solar_average_cost.py   # GENERATES external_solar_average_cost.json -- cost per configuration
   │   ├── external_solar_average_cost.json  # Averaged-probability cost on BS2005-AGS,OP, eight configurations
   │   ├── gen_shock_cost.py           # GENERATES external_shock_cost.json -- Figure 13, in two phases
   │   ├── external_shock_cost.json    # Cost of a fixed accuracy on the shock, 23 front widths, five arms
   │   ├── make_shock_scan_references.py  # FREEZES shock_reference_scan.json -- one DOP853 oracle per front width
   │   ├── shock_reference_scan.json   # Figure 13's frozen references, keyed at full precision, with both fingerprints
   │   ├── check_shock_adiabaticity.py  # Tests whether Figure 13's cost peak sits at the adiabatic crossover
   │   ├── check_shock_commutator.py   # Tests whether that peak is the commutator term, using order 2 as the control
   │   ├── gen_njobs_scaling.py        # GENERATES external_njobs_scaling.json -- Figure 11, one thread per process
   │   ├── external_njobs_scaling.json  # Wall clock against n_jobs on an Earth chord, nine arms at four scan sizes
   │   ├── sterile_projector_check.py  # Reproduces the sterile projector defect and its fix, three arms, one command
   │   ├── retime_magnus_series.py     # RE-TIMES both codes in Figure 12; references and grids untouched
   │   ├── prem_chord_common.py        # The PREM chord at cos(theta_z) = -0.9, shared by the two scripts below
   │   ├── gen_prem_reference.py       # GENERATES prem_chord_reference.json -- segment-aligned, layer edges respected
   │   ├── prem_chord_reference.json   # That reference; PARTIAL, 4nu stops at 6 of 12 energies and 5nu is unstarted
   │   ├── gen_prem_benchmarks.py      # GENERATES external_prem_chord_benchmarks.json -- the Earth analogue of Fig. 11
   │   ├── external_prem_chord_benchmarks.json  # That file: both codes on one requested tolerance, Earth chord, 2-5 flavors
   │   ├── gen_prem_plane_magnus.py    # RE-TIMES the Magnus series of the two Earth speed-accuracy planes through the harness
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
   │   ├── shock_reference.json        # That oracle, as exact hex floats
   │   └── solar_models_cache.json     # Notebook 13's comparison of the twelve solar models, keyed on its inputs
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
   │       ├── HANDOVER-audit.md       # Handover for the manuscript audit
   │       ├── HANDOVER-nuoscprobexact-batched-tolerance.md  # Handover: giving NuOscProbExact a tolerance dial, so both codes answer one request
   │       ├── audit-report.md         # What the manuscript audit found
   │       ├── AUDIT_2026-09-13.md     # The deep audit of 2026-09-13: 98 numbered items, ticked as they are done
   │       ├── AUDIT_2026-09-13_replacements.md  # The replacement text the deep audit proposed, item by item
   │       ├── PLAN_fig12_revamp.md    # Scoping for rebuilding Fig. 12 in Fig. 11 shape -- not started
   │       ├── PLAN_fig13_solar_average.md  # Scoping for Fig. 13, the cost of one averaged solar probability
   │       ├── main.tex                # The paper -- ordinary LaTeX; a revision diff is mechanical
   │       ├── refs.bib                # NuOscProbExact's bibliography, with the Magnus entries appended below a separator
   │       ├── elsarticle.cls          # Bundled, so the folder compiles without the Elsevier bundle
   │       ├── elsarticle-num.bst
   │       └── figs/                   # Its twenty-three figures, written by notebook 28
   ├── tools/                          # Standalone utilities that are not part of the package
   │   ├── auto_tight/                 # Issue #120's measurement of 'auto' at tight tolerances, with its DOP853 references
   │   ├── build_solar_model_tables.py  # Trims the authors' solar-model files to the shipped tables, checking each hash
   │   ├── lint_notebook_cells.py      # Finds names the notebooks use but never define; run by the lint workflow
   │   ├── make_demo_video.py          # Joins and shrinks notebook 27's clips; shared with NuOscProbExact
   │   └── trailer/                    # The trailer (issue #99): its script, and the scripts that compute and draw its scenes
   ├── src/                            # The package itself -- the only thing a `pip install` delivers
   │   ├── magnus/                     # Main Python package
   │   │   ├── __init__.py             # Explicit named imports from the four hamiltonians{2,3,4,5}nu.py modules
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
   │   │   ├── hamiltonians/           # 2nu-5nu Hamiltonians: vacuum, matter, NSI, LIV (the one true subpackage)
   │   │   │   ├── __init__.py         # Explicit named imports from the four hamiltonians{2,3,4,5}nu.py modules
   │   │   │   ├── _angles.py          # Interprets the four angles conventions; rejects an out-of-range sine
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
   │   │   ├── plotting.py             # Pre-packaged plotting tools: one call instead of thirty lines
   │   │   ├── py.typed                # PEP 561 marker: tells type checkers the annotations are real
   │   │   ├── solarmodels.py          # Tabulated standard solar models, as profiles for the Sun wrappers
   │   │   └── version.py              # Resolves the version from pyproject.toml (internal)
   │   └── requirements.txt            # The five runtime dependencies: numpy, scipy, joblib, matplotlib, numba
   └── tests/                          # Test suite (pytest; runs in CI)
       ├── test_paper_cache_only.py    # MAGNUS_PAPER_CACHE_ONLY stops notebook 28 on a cache miss instead of recomputing
       ├── test_paper_cache_key_is_portable.py  # The figure cache's key survives a change of machine: a ULP must not move it
       ├── test_ci_honours_the_docs.py  # Every MAGNUS_* variable the docs tell CI to set, a workflow actually sets
       ├── test_notebooks_match_their_generator.py  # The committed .ipynb files are the ones make_notebooks.py builds
       ├── test_paper_assets_are_tracked.py  # Every figure main.tex includes is tracked, which .gitignore's *.pdf defeats
       ├── test_readme_lists_every_notebook.py  # notebooks/README.md describes every notebook make_notebooks.py builds
       ├── test_adiabatic_validation_table.py  # adiabatic_strategy.rst's speed-up table and make_figures.py's chart agree
       ├── test_cli_examples_match.py  # cli.rst's worked examples still print what the CLI prints
       ├── test_diagnostics_documents_every_warning.py  # diagnostics.rst's catalogue covers every warning class the package defines
       ├── conftest.py                 # Path setup so magnus is importable without installation
       ├── test_adiabatic.py           # Adiabatic + Magnus hybrid strategy: detection, merging, ODE cross-checks
       ├── test_array_arguments.py     # Arguments other than energy and L refuse arrays, naming them (issue #116)
       ├── test_angles.py              # The four `angles` conventions and the guards between them
       ├── test_avgprob.py             # Phase-averaged probabilities
       ├── test_builders_broadcast.py  # Every Hamiltonian builder broadcasts over its scalar argument (issue #155)
       ├── test_phase_average.py       # The phase average over an energy spread (issue #64)
       ├── test_phase_groups.py        # Batched scans split into groups of similar phase (issue #111)
       ├── test_silently_ignored_keywords.py  # Keywords that were accepted and ignored now work or raise (issues #110, #112, #114)
       ├── test_scan_potential_forms.py  # Energy scans with scalar-valued or scalar-only density functions (issue #113)
       ├── test_cli.py                 # magnus command-line calculator
       ├── test_pseudodirac.py         # Pseudo-Dirac Hamiltonians: the Dirac limit, blocks, and the factor of two
       ├── test_documented_examples.py  # Runs the code blocks in README.md and quickstart.rst
       ├── test_earth_matter.py        # PREM profile, chord geometry, electron density
       ├── test_engines.py             # Which engine answers, and the cross-checks between them
       ├── test_auto_ladder.py         # When strategy='auto' takes the ladder (issue #70)
       ├── test_expansionterms.py      # The symbolic term generator against the hand-written orders
       ├── test_expm_backend.py        # The two matrix-exponential backends, their switch, and degeneracies
       ├── test_fuzz_statistics.py     # Randomized profiles, scored in bulk
       ├── test_file_tree.py           # This file: generates the tree above and checks it against git
       ├── test_globaldefs.py          # NuFIT historical parameter dict/loader
       ├── test_hamiltonians.py        # Hamiltonian/mixing-matrix builders
       ├── test_input_fuzz.py          # Issue #160: each bad argument refused by name, each valid edge accepted
       ├── test_invariants.py          # Properties that must hold across the whole engine matrix
       ├── test_magnus_expansion.py    # Magnus-core correctness (terms, orders, GL rates, unitarity)
       ├── test_one_sided_breakpoints.py  # 'trapezoid'/'simpson' sample each side of a declared breakpoint with its own values
       ├── test_oscprob.py             # Oscillation-probability engine, closed-form and ODE cross-checks
       ├── test_palindrome.py          # The palindromic-profile optimization and its gate
       ├── test_plotting.py            # Pre-packaged plotting tools: house-style defaults, layouts
       ├── test_routine_listings.py    # Each module's Routine listings names every public function it defines
       ├── test_separable_breakpoints.py  # The energy-batched engine on grids with breakpoints: real refinement only (issue #71)
       ├── test_separable_gl_gate.py   # The energy-batched 'gl' ladder refuses an agreement while an energy's slabs are wide (issue #71)
       ├── test_separable_slab_cap.py  # The energy-batched engine at its slab cap: an agreement there warns (issue #71)
       ├── test_solarmodels.py         # Solar-model tables, their profiles, and the Sun wrappers that use them
       ├── test_thread_safety.py       # Concurrent calls from several threads match the serial answer (issue #153)
       ├── test_tolerance.py           # What rtol/atol promise, and the effective-refinement gate
       ├── test_undeclared_features.py  # Undeclared jumps and narrow spikes do not end the refinement ladder on a wrong answer (issue #161)
       ├── test_validation.py          # Input-validation guards and their error messages
       ├── test_validation_160_reopened.py  # The #160 sub-cases found open after the first validation pass, one test each
       ├── test_validate_helpers.py    # The shared argument checks: reals, integers, bools, slab edges, Hamiltonians
       └── test_version.py             # Version resolution from pyproject.toml / installed metadata
