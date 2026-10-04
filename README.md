# Magnus

[![CI Tests](https://github.com/mbustama/Magnus/actions/workflows/tests.yml/badge.svg)](https://github.com/mbustama/Magnus/actions/workflows/tests.yml)
[![Code Quality](https://github.com/mbustama/Magnus/actions/workflows/lint.yml/badge.svg)](https://github.com/mbustama/Magnus/actions/workflows/lint.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://mbustama.github.io/Magnus/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![codecov](https://codecov.io/gh/mbustama/Magnus/branch/main/graph/badge.svg)](https://codecov.io/gh/mbustama/Magnus)
[![PyPI](https://img.shields.io/pypi/v/magnuspy.svg)](https://pypi.org/project/magnuspy/)
[![Downloads](https://pepy.tech/badge/magnuspy)](https://pepy.tech/project/magnuspy)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)


**Magνs** computes neutrino oscillation probabilities for two to five flavors, or
for any Hermitian Hamiltonian, in vacuum, in matter, through the Earth and through
the Sun.  Its evolution operator is unitary to round-off at any accuracy setting,
so every probability is non-negative and every row sums to one.

## Installation

```shell
pip install magnuspy
```

Python 3.10 or newer.  The distribution is **magnuspy** on PyPI, because plain
`magnus` was taken; the import package is `magnus`.  From a checkout, use
`pip install -e .`, or `pip install -e '.[test]'` to run the tests.

## Your first probability

```python
import magnus.oscprob as oscprob
import magnus.globaldefs as gd

# 3 flavors, 1 GeV, 1300 km of vacuum
P = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM)
print(P[gd.NUMU][gd.NUE])        # P(nu_mu -> nu_e) = 0.0313
```

`P[i][f]` is the probability for $\nu_i \to \nu_f$, with `gd.NUE`, `gd.NUMU`,
`gd.NUTAU` = 0, 1, 2.  **Energies and distances are in natural units**: multiply
by `gd.UNIT_GEV` and `gd.UNIT_KM` on the way in.  Oscillation parameters default
to the NuFIT 6.1 best fit, normal ordering.

The [quick start guide](https://mbustama.github.io/Magnus/quickstart.html) takes
you from here through matter, the Earth, the Sun, antineutrinos, new physics and
a Hamiltonian of your own, with every example runnable as written.

## More in a few lines

```python
import numpy as np
import magnus.earth as earth

osc = gd.load_nufit_params('NuFIT 6.1')          # best-fit parameters, as sines
# Angles come back as sines (s12 = 0.556 is sin(theta12), so sin^2 = 0.309),
# phases in radians.  Leaving **osc out gives the same numbers: unset
# parameters default to this same fit.

# 200 energies through the Earth in one call: cos(theta_z) = -0.8 fixes the
# direction of the chord, and magnus.earth gives its length
E = np.linspace(0.5, 10.0, 200)*gd.UNIT_GEV
L = earth.distance_traveled_inside_earth(-0.8)*gd.UNIT_KM
P = oscprob.osc_prob_3nu_earth(E, costhz=-0.8, L=L, **osc)          # shape (200, 3, 3)

# Antineutrinos, one channel
P_bar = oscprob.osc_prob_3nu_earth(E, costhz=-0.8, L=L, nubar=True,
                                   nu_i=gd.NUMU, nu_f=gd.NUE, **osc)   # shape (200,)

# Solar neutrinos from the center of a standard solar model, phase-averaged
# over a 10% energy resolution, as a solar experiment measures them
E_sun = np.logspace(-1, np.log10(20.0), 50)*gd.UNIT_MEV
P = oscprob.osc_prob_3nu_sun(E_sun, gd.SUN_RADIUS*gd.UNIT_KM, 0.0, **osc,
                             density_profile='B16-GS98', average=True)

# The evolution operator alongside the probabilities
P, U = oscprob.osc_prob_3nu_earth(1.0*gd.UNIT_GEV, costhz=-0.8, L=L,
                                  return_evolution_operator=True, **osc)
```

A phase-averaged result (`average=True`) depends on the energy spread,
`average_spread`, wherever some oscillation phase is small, and a
`PhaseAveragingWarning` flags those cases.  The
[numerical recipes](https://mbustama.github.io/Magnus/recipes.html) page has
a runnable snippet for each common task.

### Or from the command line

```bash
$ magnus --flavors 3 --environment vacuum \
    --energy 1 --energy-unit GeV --baseline 1300 --baseline-unit km
Magνs 1.2.0 -- osc_prob_3nu_vacuum
E = 1 GeV, L = 1300 km

            nu_e   nu_mu  nu_tau
nu_e      0.9289  0.0085  0.0625
nu_mu     0.0313  0.3923  0.5764
nu_tau    0.0398  0.5992  0.3611
```

Add `--scenario nsi`, `--flavors 5`, `--nu-i e --nu-f mu` for one channel, or
`--json` to pipe the result into another program.  The command computes one
probability at a time; scans, averages and custom Hamiltonians need Python.
Full reference: [CLI](https://mbustama.github.io/Magnus/cli.html).

## Why Magνs

**Flexible.**  The Hamiltonian is an argument, not an assumption.  Standard
oscillations, non-standard interactions, Lorentz-invariance violation, sterile
states, pseudo-Dirac pairs and a model of your own all go through the same call.
Two to five flavors ship ready-made; the generic entry points take any
dimension and any profile, given as a function of position.

**Fast.**  An energy scan is one batched call rather than a loop, which makes each
probability one to two orders of magnitude cheaper.  The cost of a probability
follows how fast the density varies, not how many times the neutrino oscillates.
[Performance](#performance) has the timings.

**Accurate.**  Magνs propagates the evolution operator with the **Magnus
expansion**: it exponentiates truncated integrals of the Hamiltonian over a chain
of position slabs.  At a tight tolerance, it
agrees with an independent integration to a few parts in **10¹²** at two to five
flavors.  Where it cannot certify its own answer, it says so.

The [Against other codes](https://mbustama.github.io/Magnus/comparison.html#when-to-use-magnus)
page measures when to use Magνs, and when another code is cheaper.

> **How do I say that?** Like the name **Magnus**: the Greek letter **ν** (nu),
> the neutrino's symbol, stands in for the "nu" syllable.

## What it can compute

- Oscillations through a **varying profile**: the Earth's PREM layers, any of
  twelve tabulated standard solar models, a supernova shock front, or any
  density you supply.
- The **phase-averaged** probability a solar or astrophysical experiment
  measures: each interference term damped by the spread of its phase across
  the energy resolution, 10% by default, without resolving the oscillation.
- The **evolution operator** itself, alongside the probabilities, for
  observables built from amplitudes.
- The same probabilities **from a shell**, with no Python, through the
  `magnus` command.

## Examples

Each of these is one call with a different Hamiltonian, profile or observable.

- **Beam experiments** — appearance probabilities along the DUNE, T2K,
  Hyper-K and ESS chords, from two named sites
  ([notebook 04](https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb)).
- **Atmospheric oscillograms** — probability over zenith angle and energy, one
  batched energy scan per zenith angle
  ([notebook 06](https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb)).
- **Solar neutrinos** — twelve standard solar models, taken by name, and the
  averaged probability an experiment sees
  ([notebook 13](https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb)).
- **Supernova shock fronts** — where a traveling discontinuity changes the
  conversion probability itself
  ([notebook 14](https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb)).
- **Astrophysical flavor composition**, including pseudo-Dirac pairs that stay
  coherent after everything else has averaged
  ([notebook 29](https://github.com/mbustama/Magnus/blob/main/notebooks/29_magnus_pseudo_dirac.ipynb)).
- **A Hamiltonian of your own** — a long-range $L_e - L_\mu$ interaction
  sourced by the Sun's electrons, a cavity in the Earth's crust, geoneutrinos,
  or a jet inside a collapsing star
  ([notebook 19](https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb)).

Every figure below is taken from the executed output of a notebook in
[`notebooks/`](https://github.com/mbustama/Magnus/tree/main/notebooks/).

| | |
|:--:|:--:|
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_3nu_vacuum.png" width="380"/><br/>**Oscillation probabilities** against baseline or energy, for two to five flavors, in vacuum and in matter.<br/>[notebook 03](https://github.com/mbustama/Magnus/blob/main/notebooks/03_magnus_3nu_vacuum_matter.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_long_baseline.png" width="380"/><br/>**Between two points on the Earth's surface** — Fermilab to SNOLAB, Homestake, CERN and the South Pole, through PREM.<br/>[notebook 04](https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_oscillogram.png" width="380"/><br/>**Oscillograms** across zenith angle and energy, one batched energy scan per zenith angle.<br/>[notebook 06](https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_biprobability.png" width="380"/><br/>**CP violation**, as bi-probability ellipses traced by the CP phase.<br/>[notebook 05](https://github.com/mbustama/Magnus/blob/main/notebooks/05_magnus_biprobability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_sterile_3plus2.png" width="380"/><br/>**Five flavors: a 3+2 sterile spectrum**, its fast oscillation filling the three-flavor envelope.<br/>[notebook 07](https://github.com/mbustama/Magnus/blob/main/notebooks/07_magnus_bsm_sterile_nu.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_custom_h.png" width="380"/><br/>**A Hamiltonian of your own** — here a long-range $L_e - L_\mu$ force through the Earth, against the standard curve.<br/>[notebook 19](https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_density_arrangement.png" width="380"/><br/>**Arrangement matters**: the same average density over the same path length, ordered differently, gives different probabilities.<br/>[notebook 18](https://github.com/mbustama/Magnus/blob/main/notebooks/18_magnus_unusual_density_profiles.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_averaged.png" width="380"/><br/>**Phase-averaged probabilities**: what remains when the oscillation is too fast for a detector to resolve.<br/>[notebook 10](https://github.com/mbustama/Magnus/blob/main/notebooks/10_magnus_averaged_probability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_solar_averaged.png" width="380"/><br/>**The averaged solar survival probability**, returned directly in about 0.7 s. The green trace is the *instantaneous* probability that another code returns, which oscillates between 0.15 and 0.9.<br/>[notebook 25](https://github.com/mbustama/Magnus/blob/main/notebooks/25_magnus_against_other_codes.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_solar_bsm.png" width="380"/><br/>**BSM against the standard curve**: NSI and a sterile state on a real BS2005 solar model, with the departure below.<br/>[notebook 13](https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_shock_bsm.png" width="380"/><br/>**The same two scenarios on a supernova shock**, where the same $\varepsilon$ moves the probability by up to 0.14, against 0.014 in the Sun.<br/>[notebook 14](https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb) | |

## When is Magνs not the right tool?

Magνs solves the Schrödinger equation for a Hermitian Hamiltonian fixed before
the propagation.  Decoherence, coupling to a bath and decay to invisible states
need a non-unitary evolution; [nuSQuIDS](https://github.com/arguelles/nuSQuIDS)
is the tool for those.  Collective oscillations need a Hamiltonian that depends
on the solution, and Magνs does not ship the self-consistent iteration.  A
feature narrower than every sampling grid goes unseen by all the engines alike;
Magνs scans the profile for such features and warns, but
the scan does not catch every one.  Magνs is not a flux, cross-section or detector
code, a fitting framework, or an event generator: it computes oscillation
probabilities and stops there.

## Performance

A single three-flavor probability through the Earth takes about 2 ms at the
default tolerance of 10⁻³.  Over 164 Earth, solar, vacuum and constant-density
configurations, the median call takes 2 ms and the slowest under a second.  These
times are per call, on one laptop, without the first call of a session, which
also loads the compiled kernels: 0.1 to 0.3 s, or about 2 s the first time on a
machine, when they compile.  The timing harness is
[`docs/dev/adversarial_batteries/timing.py`](https://github.com/mbustama/Magnus/blob/main/docs/dev/adversarial_batteries/timing.py).

**Pass arrays instead of looping.**  Every wrapper takes an array of energies,
of baselines or of both, and shares work across the points: worth one to two
orders of magnitude, at every number of flavors from two to five.

**Write your `H_func` to accept an array of positions.**  It is then called once
per refinement stage rather than once per quadrature node, which is several times
faster, with identical output; the
[numerical recipes](https://mbustama.github.io/Magnus/recipes.html#write-h-func-vectorized)
show how.  A Hamiltonian that ignores its argument is detected and broadcast
automatically.

**`n_jobs` helps only where no batched engine applies**: there, ten workers make
a scan 2 to 3 times faster.  Where a batched engine applies, a single process is
faster.  See [performance](https://mbustama.github.io/Magnus/performance.html).

## What "accurate" means here

Unitarity holds at any order and tolerance (see above).  Accuracy is measured
against closed forms, an independent ODE solver, an independently generated
expansion and populations of random profiles; the [table of checks and
results](https://mbustama.github.io/Magnus/index.html#what-accuracy-means) is in
the documentation.  `rtol` and `atol` are a **stopping rule**, not a guarantee;
see [what they
control](https://mbustama.github.io/Magnus/diagnostics.html#what-rtol-and-atol-control).

## Salient features

- **Any number of flavors, any Hamiltonian**: validated wrappers for 2ν, 3ν,
  4ν (3+1) and 5ν (3+2), named by environment and scenario
  ([functions](https://mbustama.github.io/Magnus/functions.html)), plus a
  generic `osc_prob` that takes a Hermitian Hamiltonian of any dimension.
- **Vacuum, matter, Earth and Sun**: constant and exponentially falling
  density, the Earth through PREM including chords between named sites, the
  Sun on an exponential fit or any of twelve standard solar models
  ([solar models](https://mbustama.github.io/Magnus/solar_models.html)), or
  any profile you supply.
- **The Magnus expansion to order 10**, with the terms of orders 1 to 10
  checked against an independently coded recursion.  The default Gauss–Legendre
  integrators reach orders 2, 4, 6 and 8 from 1, 2, 3 and 4 evaluations of the
  Hamiltonian per slab; cumulative trapezoid and Simpson quadrature reach
  order 10.
- **Adaptive refinement** to a requested tolerance, slab edges on density
  discontinuities, and an energy-batched engine for scans; seven engines in
  all, chosen from the request
  ([engines](https://mbustama.github.io/Magnus/engines.html)).
- **Plotting**: the figures the notebooks use, one call each, from curves and
  small multiples to bi-probability ellipses and oscillograms
  ([plotting](https://mbustama.github.io/Magnus/plotting.html)).

## Documentation

Full documentation: **[mbustama.github.io/Magnus](https://mbustama.github.io/Magnus/)**.

| | |
|---|---|
| [Numerical recipes](https://mbustama.github.io/Magnus/recipes.html) | What it can compute, with code |
| [Tutorials](https://mbustama.github.io/Magnus/tutorials.html) | All 29 notebooks, with what each is for |
| [Phase-averaged probabilities](https://mbustama.github.io/Magnus/averaged_probability.html) | What `average=True` returns, and at what cost |
| [Standard solar models](https://mbustama.github.io/Magnus/solar_models.html) | The twelve tables that ship with the package |
| [Mathematical method](https://mbustama.github.io/Magnus/methodology.html) | The expansion derived term by term, and why truncation is unitary |
| [Engines and dispatch](https://mbustama.github.io/Magnus/engines.html) | Which of the seven engines answers a call, and why |
| [Accuracy and diagnostics](https://mbustama.github.io/Magnus/diagnostics.html) | What each safeguard cannot catch, and every warning explained |
| [Against other codes](https://mbustama.github.io/Magnus/comparison.html) | The full cross-code comparison, with the measurements behind it |
| [API reference](https://mbustama.github.io/Magnus/api_reference.html) | Every public function, generated from the source |

## Repository layout

The top level of the repository.  The
[documentation](https://mbustama.github.io/Magnus/installation.html#file-tree) also
shows the package under `src/`.  The complete listing, with a comment on every file,
is `TREE` in
[`tests/test_file_tree.py`](https://github.com/mbustama/Magnus/blob/main/tests/test_file_tree.py),
which checks it against the files in the repository.

```text
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
└── tests/                          # Test suite (pytest; runs in CI)
```

## Continuous integration

Every push runs the test suite on Python 3.10 through 3.13 and lints with Ruff;
the documentation is built with warnings as errors.  The 29 notebooks
execute too, across four parallel shards, with a notebook served from cache
when neither the package nor that notebook has changed.  The badges at the top
report those runs.

## Requirements

`numpy (>= 1.22)`, `scipy (>= 1.9)`, `joblib (>= 1.2)`, `matplotlib (>= 3.5)` (so that
`magnus.plotting` works in any installation) and `numba (>= 0.59)` (for the
fast matrix exponential).  Tests:

```bash
pip install -e '.[test]' && pytest tests/
```

Add `--cov --cov-report=term-missing` for coverage, configured in
`pyproject.toml` so it measures the same thing locally and in CI.

## Changelog

See [CHANGELOG.md](https://github.com/mbustama/Magnus/blob/main/CHANGELOG.md),
also [rendered in the docs](https://mbustama.github.io/Magnus/changelog.html).

## How to cite

If Magnus contributed to work you are publishing, please cite it and say which
version you used, since results can depend on it.  The
[citing](https://mbustama.github.io/Magnus/citing.html) page has the BibTeX
entry.

Mauricio Bustamante (2026).  *Magνs: neutrino oscillation probabilities for any
Hermitian Hamiltonian, any number of flavors, and any matter profile*.  GitHub repository:
https://github.com/mbustama/Magnus.

**Methodology references:**
* Sergio Blanes, Fernando Casas, José A. Oteo & José Ros (2009). The Magnus
  expansion and some of its applications. *Physics Reports, 470*(5–6),
  151–238. [doi:10.1016/j.physrep.2008.11.001](https://doi.org/10.1016/j.physrep.2008.11.001).
* Sergio Blanes, Fernando Casas & Javier Ros (2000). Improved high order
  integrators based on the Magnus expansion. *BIT Numerical Mathematics,
  40*(3), 434–450. [doi:10.1023/A:1022311628317](https://doi.org/10.1023/A:1022311628317).

## License

GNU General Public License v3.0 only (`GPL-3.0-only`); the full text is in
[LICENSE](https://github.com/mbustama/Magnus/blob/main/LICENSE).  You may use,
study, modify and redistribute it, including commercially, provided derivative
works carry the same license with source available.

## Author

Mauricio Bustamante (mbustamante@gmail.com).  Bug reports and questions are
best raised as [GitHub issues](https://github.com/mbustama/Magnus/issues), which
leave a public record others can find.
