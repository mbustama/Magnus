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
for any Hamiltonian you write, in vacuum, in matter, through the Earth and through
the Sun.  Its evolution operator is exactly unitary by construction, so every
probability is non-negative and every row sums to one, at any accuracy setting.

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

The [numerical recipes](https://mbustama.github.io/Magnus/recipes.html) page has
a runnable snippet for each common task.

### Or from the command line

```bash
$ magnus --flavors 3 --environment vacuum \
    --energy 1 --energy-unit GeV --baseline 1300 --baseline-unit km
Magνs 1.1.1 -- osc_prob_3nu_vacuum
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

**Fast.**  An energy scan is one batched call rather than a loop, worth about an
order of magnitude per probability; an oscillogram is one such call per zenith
angle.  The median call over 164 Earth and solar configurations is
**2 ms**; [Performance](#performance) has the rest.

**Accurate.**  Internally, Magνs propagates the evolution operator with the
**Magnus expansion**: it exponentiates truncated integrals of the Hamiltonian over
a chain of position slabs, and every truncation is exactly unitary.  Asked for it, it agrees with an
independent integration to a few parts in **10¹²** at two to five
flavors.  Where it cannot certify its own answer, it says so.

> **How do I say that?** Just like the name **Magnus** — the Greek letter
> **ν** (nu), the neutrino's symbol, stands in for the "nu" syllable.  (Most of
> this package was written in Denmark, so
> [the Danish way](https://translate.google.com/?sl=da&tl=en&text=Magnus&op=translate)
> is welcome too.)

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

## What it has been used for

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
- **Supernova shock fronts** — where a travelling discontinuity changes the
  conversion probability itself
  ([notebook 14](https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb)).
- **Astrophysical flavor composition**, including pseudo-Dirac pairs that stay
  coherent after everything else has averaged
  ([notebook 29](https://github.com/mbustama/Magnus/blob/main/notebooks/29_magnus_pseudo_dirac.ipynb)).
- **A Hamiltonian of your own** — a long-range $L_e - L_\mu$ interaction
  sourced by the Sun's electrons, a cavity in the Earth's crust, geoneutrinos,
  or a jet inside a collapsing star
  ([notebook 19](https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb)).

Every figure below comes from a notebook in
[`notebooks/`](https://github.com/mbustama/Magnus/tree/main/notebooks/), lifted
out of the executed file, so what you see is what that notebook produced.

| | |
|:--:|:--:|
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_3nu_vacuum.png" width="380"/><br/>**Oscillation probabilities** against baseline or energy, for two to five flavors, in vacuum and in matter.<br/>[notebook 03](https://github.com/mbustama/Magnus/blob/main/notebooks/03_magnus_3nu_vacuum_matter.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_long_baseline.png" width="380"/><br/>**Between two points on the Earth's surface** — Fermilab to SNOLAB, Homestake, CERN and the South Pole, through PREM.<br/>[notebook 04](https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_oscillogram.png" width="380"/><br/>**Oscillograms** across zenith angle and energy, one batched energy scan per zenith angle.<br/>[notebook 06](https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_biprobability.png" width="380"/><br/>**CP violation**, as bi-probability ellipses traced by the CP phase.<br/>[notebook 05](https://github.com/mbustama/Magnus/blob/main/notebooks/05_magnus_biprobability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_sterile_3plus2.png" width="380"/><br/>**Five flavors: a 3+2 sterile spectrum**, its fast oscillation filling the three-flavor envelope.<br/>[notebook 07](https://github.com/mbustama/Magnus/blob/main/notebooks/07_magnus_bsm_sterile_nu.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_custom_h.png" width="380"/><br/>**A Hamiltonian of your own** — here a long-range $L_e - L_\mu$ force through the Earth, against the standard curve.<br/>[notebook 19](https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_density_arrangement.png" width="380"/><br/>**Arrangement beats the mean**: the same average density and the same path length, ordered differently, give different probabilities.<br/>[notebook 18](https://github.com/mbustama/Magnus/blob/main/notebooks/18_magnus_unusual_density_profiles.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_averaged.png" width="380"/><br/>**Phase-averaged probabilities** — what survives when the oscillation is faster than anything can resolve.<br/>[notebook 10](https://github.com/mbustama/Magnus/blob/main/notebooks/10_magnus_averaged_probability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_solar_averaged.png" width="380"/><br/>**The averaged solar survival probability**, returned directly in about 0.7 s. The green trace is the *instantaneous* probability another code returns, thrashing between 0.15 and 0.9.<br/>[notebook 25](https://github.com/mbustama/Magnus/blob/main/notebooks/25_magnus_against_other_codes.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_solar_bsm.png" width="380"/><br/>**BSM against the standard curve**: NSI and a sterile state on a real BS2005 solar model, with the departure below.<br/>[notebook 13](https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_shock_bsm.png" width="380"/><br/>**The same two scenarios on a supernova shock**, where the identical $\varepsilon$ moves the answer thirty times further.<br/>[notebook 14](https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb) | |

## When is Magνs the right tool?

Three advantages follow from the method, whatever the Hamiltonian contains.
Magνs is **robust**: the evolution operator is exactly unitary.  It is **fast**
without giving up accuracy: its cost follows the profile, not the phase, where
an ODE solver pays for every radian it resolves, and a scan over energy or
direction is one batched call on compiled kernels.  And it is **flexible**: the
Hamiltonian is a callable returning a Hermitian matrix of any size.

Reach for it when the density varies fast against the oscillation length, as in
the Sun; when the profile has structure at a known place, declared with
`t_breakpoints`; for more than three flavors; for a Hamiltonian with no closed
form; for an averaged observable; or for an accuracy below where a composition
of constant-density slabs floors.  Elsewhere another code may be cheaper: at
constant density NuFast-LBL is about sixty times faster, and on an Earth chord
resolved to no better than about 10⁻⁶
[NuOscProbExact](https://github.com/mbustama/NuOscProbExact) is cheaper.  The
measurements are on the
[Against other codes](https://mbustama.github.io/Magnus/comparison.html) page.

## When is it not the right tool?

Magνs solves the Schrödinger equation for a Hermitian Hamiltonian fixed before
the propagation.  Decoherence, coupling to a bath and decay to invisible states
need a non-unitary evolution; [nuSQuIDS](https://github.com/arguelles/nuSQuIDS)
is the tool for those.  Collective oscillations need a Hamiltonian that depends
on the solution, and Magνs does not ship the self-consistent iteration.  And a
feature narrower than every sampling grid is missed by every engine together;
the matter scenario functions scan the profile for such features and warn, but
the scan catches most, not all.  Magνs is not a flux, cross-section or detector
code, a fitting framework, or an event generator: it computes oscillation
probabilities and stops there.

## Performance

A single three-flavor Earth probability takes about 2 ms at the default
tolerance of 10⁻³; across 164 Earth and solar configurations the median call
takes 2 ms and the slowest under a second.

**Pass arrays instead of looping.**  Every wrapper takes an array of energies,
of baselines or of both, and shares work across the points: worth about an order
of magnitude at two and three flavors, and several-fold at four and five.

**Write your `H_func` to accept an array of positions.**  It is then called once
per refinement stage rather than once per quadrature node: several times faster,
with identical output.  The trick is broadcasting the potential into a stack of
matrices:

```python
def H_func(l):
    l = np.asarray(l, dtype=float)
    VCC = VCC_central*np.exp(-(l/gd.UNIT_KM)/l_scale)   # an array
    return (1.0/energy)*h_vac + VCC[..., None, None]*e00
```

A Hamiltonian that ignores its argument is detected and broadcast already.  The
fallback warns once per session, naming the fix.

**`n_jobs` helps only where no batched engine applies**, 2 to 3 times with ten
workers; where one applies, a single process is faster.  More:
[performance](https://mbustama.github.io/Magnus/performance.html).

## What "accurate" means here

Every truncation of the Magnus series is anti-Hermitian, so the evolution is
unitary at any order and tolerance.  The rest is measured:

| Checked against | Result |
|---|---|
| Every expansion term, Ω₁ to Ω₁₀, against an independent recursion | relative 10⁻¹¹ |
| 2ν and 3ν vacuum, 2ν constant-density matter: the closed forms | 10⁻¹² |
| Halving the slab width at orders 2 / 4 / 6, against DOP853 | error ÷ 4 / 16 / 64 |
| Repeated calls, and a baseline scan in any order | identical, bit for bit |
| A parallel run against a serial one | to the requested tolerance |
| 40 random smooth profiles, default tolerance | median 10⁻⁸, one silent miss |
| 120 random step profiles, edges undeclared | 19 outside tolerance, all warned |

A *silent miss* is an answer outside the tolerance with no warning.  `rtol` and
`atol` are a **stopping rule**, not a guarantee: at the default 10⁻³ an Earth
probability is usually far more accurate (median 10⁻⁶ from the same call at
10⁻⁷, largest about 10⁻³).  See
[what they actually control](https://mbustama.github.io/Magnus/diagnostics.html#what-rtol-and-atol-actually-control).

## Salient features

- **Two ways to use it**: as a Python module, or as a `magnus` command-line
  calculator for a single probability.
- **Any number of flavors, any Hamiltonian**: validated wrappers for 2ν, 3ν,
  4ν (3+1) and 5ν (3+2), named by environment and scenario
  ([functions](https://mbustama.github.io/Magnus/functions.html)), plus a
  generic `osc_prob` that takes a Hermitian Hamiltonian of any dimension.
- **Vacuum, matter, Earth and Sun**: constant and exponentially falling
  density, the Earth through PREM including chords between named sites, the
  Sun on an exponential fit or any of twelve standard solar models
  ([solar models](https://mbustama.github.io/Magnus/solar_models.html)), or
  any profile you supply.
- **Beyond the Standard Model**: non-standard interactions, Lorentz-invariance
  violation, sterile states and pseudo-Dirac pairs.
- **The Magnus expansion to order 10**, with the terms of orders 1 to 6
  checked against an independently coded recursion.  The default Gauss–Legendre integrators reach orders
  2, 4, 6 and 8 from 1, 2, 3 and 4 evaluations of the Hamiltonian per slab;
  cumulative trapezoid and Simpson quadrature reach order 10.
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

The top level only; the [complete listing with a comment on every
file](https://mbustama.github.io/Magnus/installation.html) is in the
documentation.  Both are generated from `git ls-files`, so neither can drift.

```text
Magnus/
├── .github/                        # GitHub Actions workflows: tests, lint, notebooks, docs, publishing
├── .gitignore                      # Build, cache and generated-output artifacts
├── CHANGELOG.md                    # Version history (Keep a Changelog format)
├── CITATION.cff                    # Machine-readable citation metadata; drives GitHub's "Cite this repository"
├── LICENSE                         # GNU GPL v3 (GPL-3.0-only), the full license text
├── README.md                       # This file
├── docs/                           # Sphinx documentation configuration and source
├── fig/                            # Plots produced by the example notebooks
├── img/                            # Figures used by the documentation
├── notebooks/                      # Numbered Jupyter notebooks -- see docs/source/tutorials.rst
├── pyproject.toml                  # Build system, dependencies, and the `magnus` console-script entry point
├── resources/                      # Travels with the code; reaches neither the wheel nor the sdist
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

`numpy (>= 1.22)`, `scipy (>= 1.9)`, `joblib`, `matplotlib` (so that
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

Mauricio Bustamante (2026).  *Magnus: neutrino oscillation probabilities via
the Magnus expansion*.  GitHub repository:
https://github.com/mbustama/Magnus.

**Methodology references:**
* Sergio Blanes, Fernando Casas, José A. Oteo & José Ros (2009). The Magnus
  expansion and some of its applications. *Physics Reports, 470*(5-6),
  151-238. [doi:10.1016/j.physrep.2008.11.001](https://doi.org/10.1016/j.physrep.2008.11.001).
* Sergio Blanes, Fernando Casas & Javier Ros (2000). Improved high order
  integrators based on the Magnus expansion. *BIT Numerical Mathematics,
  40*(3), 434-450. [doi:10.1023/A:1022311628317](https://doi.org/10.1023/A:1022311628317).

## License

GNU General Public License v3.0 only (`GPL-3.0-only`); the full text is in
[LICENSE](https://github.com/mbustama/Magnus/blob/main/LICENSE).  You may use,
study, modify and redistribute it, including commercially, provided derivative
works carry the same license with source available.

## Author

Mauricio Bustamante (mbustamante@gmail.com).  Bug reports and questions are
best raised as [GitHub issues](https://github.com/mbustama/Magnus/issues), which
leave a public record others can find.
