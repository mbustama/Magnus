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

Neutrino oscillation probabilities for **any Hermitian Hamiltonian, any number
of flavors, and any matter profile**.

> **How do I say that?** Just like the name **Magnus** — the Greek letter
> **ν** (nu), the neutrino's symbol, stands in for the "nu" syllable.  (Most of
> this package was written in Denmark, so
> [the Danish way](https://translate.google.com/?sl=da&tl=en&text=Magnus&op=translate)
> is welcome too.)

Magnus computes the evolution operator with the **Magnus expansion**: over each
slab of the trajectory the operator is one exponential, built from time-ordered
integrals of the Hamiltonian.  Two consequences matter in practice.

- **Exactly unitary, by construction.**  Every truncation stays in the Lie
  algebra, so probabilities sum to one at machine precision, at any order and
  any tolerance.
- **The cost follows the profile, not the phase.**  A slab is exact for a
  constant Hamiltonian however many oscillations it spans, so the Sun, with its
  hundreds of thousands of cycles, is as tractable as the Earth.

## Installation

```shell
pip install magnuspy
```

Python 3.10 or newer.  The distribution is **magnuspy** on PyPI, because plain
`magnus` was taken; the import package is `magnus`.  From a checkout, use
`pip install -e .`, or `pip install -e '.[test]'` to run the tests.

## What you can compute

Every figure comes from a notebook in
[`notebooks/`](https://github.com/mbustama/Magnus/tree/main/notebooks/), lifted
out of the executed file, so what you see is what that notebook produced.

| | |
|:--:|:--:|
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_3nu_vacuum.png" width="380"/><br/>**Oscillation probabilities** against baseline or energy, for two to five flavors, in vacuum and in matter.<br/>[notebook 03](https://github.com/mbustama/Magnus/blob/main/notebooks/03_magnus_3nu_vacuum_matter.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_long_baseline.png" width="380"/><br/>**Between two points on the Earth's surface** — Fermilab to SNOLAB, Homestake, CERN and the South Pole, through PREM.<br/>[notebook 04](https://github.com/mbustama/Magnus/blob/main/notebooks/04_magnus_long_baseline.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_oscillogram.png" width="380"/><br/>**Oscillograms** across zenith angle and energy, in a single call.<br/>[notebook 06](https://github.com/mbustama/Magnus/blob/main/notebooks/06_magnus_oscillograms.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_biprobability.png" width="380"/><br/>**CP violation**, as bi-probability ellipses traced by the CP phase.<br/>[notebook 05](https://github.com/mbustama/Magnus/blob/main/notebooks/05_magnus_biprobability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_sterile_3plus2.png" width="380"/><br/>**Five flavors: a 3+2 sterile spectrum**, its fast oscillation filling the three-flavor envelope.<br/>[notebook 07](https://github.com/mbustama/Magnus/blob/main/notebooks/07_magnus_bsm_sterile_nu.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_custom_h.png" width="380"/><br/>**A Hamiltonian of your own** — here a long-range $L_e - L_\mu$ force through the Earth, against the standard curve.<br/>[notebook 19](https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_density_arrangement.png" width="380"/><br/>**Arrangement beats the mean**: the same average density and the same path length, ordered differently, give different probabilities.<br/>[notebook 18](https://github.com/mbustama/Magnus/blob/main/notebooks/18_magnus_unusual_density_profiles.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_averaged.png" width="380"/><br/>**Phase-averaged probabilities** — what survives when the oscillation is faster than anything can resolve.<br/>[notebook 10](https://github.com/mbustama/Magnus/blob/main/notebooks/10_magnus_averaged_probability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_solar_averaged.png" width="380"/><br/>**The averaged solar survival probability**, returned directly in about 0.7 s. The green trace is the *instantaneous* probability another code returns, thrashing between 0.15 and 0.9.<br/>[notebook 25](https://github.com/mbustama/Magnus/blob/main/notebooks/25_magnus_against_other_codes.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_solar_bsm.png" width="380"/><br/>**BSM against the standard curve**: NSI and a sterile state on a real BS2005 solar model, with the departure below.<br/>[notebook 13](https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/Magnus/main/img/gallery/gallery_shock_bsm.png" width="380"/><br/>**The same two scenarios on a supernova shock**, where the identical $\varepsilon$ moves the answer thirty times further.<br/>[notebook 14](https://github.com/mbustama/Magnus/blob/main/notebooks/14_magnus_supernova_shock.ipynb) | |

## Quick start

```python
import numpy as np
import magnus.oscprob as oscprob
import magnus.globaldefs as gd

osc = gd.load_nufit_params('NuFIT 6.1')          # best-fit parameters

# One probability: 3 flavors, 1 GeV, 1300 km of vacuum
P = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM, **osc)

# A scan: 200 energies through the Earth, in one call.  The zenith angle
# fixes the direction of the chord; magnus.earth gives its length
import magnus.earth as earth
E = np.linspace(0.5, 10.0, 200)*gd.UNIT_GEV
L = earth.distance_traveled_inside_earth(-0.8)*gd.UNIT_KM
P = oscprob.osc_prob_3nu_earth(E, costhz=-0.8, L=L, **osc)

# The phase-averaged limit an astrophysical measurement sees
P = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_TEV, 1.0e8*gd.UNIT_KM,
                                average=True, **osc)
```

Pass an array where a single energy or baseline would go, and one probability
comes back per entry.  Pass `nubar=True` for antineutrinos, `nu_i`/`nu_f` for a
single channel.  The
[numerical recipes](https://mbustama.github.io/Magnus/recipes.html) page has a
runnable snippet for each common task.

## Or from the command line

```bash
$ magnus prob --flavors 3 --environment vacuum \
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
probability at a time: scans, averages and custom Hamiltonians belong in
Python.  Full reference: [CLI](https://mbustama.github.io/Magnus/cli.html).

## When to reach for it

| Use Magnus when | Use something else when |
|---|---|
| The density varies along the path — the Sun, a supernova, a tabulated profile | The density is constant or piecewise constant, where a closed form is exact and about 20× cheaper |
| You need accuracy below 10⁻¹¹: composing constant slabs floors at 2.5 × 10⁻¹¹ and then worsens, while Magnus reaches 2.9 × 10⁻¹³ | A per-cent answer is enough |
| You have more than four flavors, or a Hamiltonian nobody has diagonalized | Your problem is standard three-flavor and speed is everything |
| You want the phase-averaged probability an experiment measures, returned rather than reconstructed | Your problem is not unitary |

That last row is the hard boundary.  Magnus solves the Schrödinger equation for
a Hermitian Hamiltonian, so decoherence, decay and collective oscillations lie
outside it: each needs a density matrix with a non-unitary term, which no
truncation of this expansion can carry.

Numbers behind each row, and a case-by-case table, are on the
[Against other codes](https://mbustama.github.io/Magnus/comparison.html) page;
[notebook 25](notebooks/25_magnus_against_other_codes.ipynb) runs the
comparison in full, every code timed in one process and refereed by a method
that is neither code's.

## What it covers

- **Flavors:** 2, 3, 4 (3+1) and 5 (3+2) through named wrappers; any number
  through the generic `osc_prob`.
- **Environments:** vacuum, constant density, exponentially falling density,
  the Earth through PREM including chords between named sites, the Sun, or any
  profile you supply.
- **Beyond the Standard Model:** non-standard interactions, Lorentz-invariance
  violation, sterile states, pseudo-Dirac pairs, or your own Hamiltonian.
- **Either kind of output:** the full probability matrix or one channel,
  single points or arrays, neutrinos or antineutrinos, instantaneous or
  phase-averaged.

Each environment and scenario has an explicitly named function, for `N` in
2, 3, 4, 5:

| Environment | Scenario | Function pattern (`{N}` = 2, 3, 4, 5) |
|---|---|---|
| Vacuum | Standard | `osc_prob_{N}nu_vacuum` |
| Vacuum | LIV | `osc_prob_{N}nu_vacuum_liv` |
| Matter, constant density | Standard | `osc_prob_{N}nu_matter_constant_density` |
| Matter, constant density | NSI | `osc_prob_{N}nu_matter_nsi_constant_density` |
| Matter, constant density | LIV | `osc_prob_{N}nu_matter_liv_constant_density` |
| Matter, exponential density | Standard | `osc_prob_{N}nu_matter_exp_density` |
| Matter, exponential density | NSI | `osc_prob_{N}nu_matter_nsi_exp_density` |
| Matter, exponential density | LIV | `osc_prob_{N}nu_matter_liv_exp_density` |
| Earth (PREM) | Standard | `osc_prob_{N}nu_earth` |
| Earth (PREM) | NSI | `osc_prob_{N}nu_earth_nsi` |
| Earth (PREM) | LIV | `osc_prob_{N}nu_earth_liv` |
| Sun | Standard | `osc_prob_{N}nu_sun` |
| Sun | NSI | `osc_prob_{N}nu_sun_nsi` |
| Sun | LIV | `osc_prob_{N}nu_sun_liv` |

Anything the table does not anticipate goes through `osc_prob`,
`osc_prob_earth` or `osc_prob_sun` directly.  Full signatures:
[functions](https://mbustama.github.io/Magnus/functions.html).

## Two features worth knowing about

**A large accumulated phase costs nothing extra.**  Where a profile varies
slowly, Magnus transports the state along the instantaneous eigenstates and
reserves the expansion for the level crossings, choosing the hand-over itself.
This is what makes a solar neutrino, with some 10⁴ radians of phase, a
sub-second calculation.  The `strategy` keyword (`'auto'`, `'hybrid'`,
`'magnus'`) overrides the choice; see
[adiabatic strategy](https://mbustama.github.io/Magnus/adiabatic_strategy.html).

**`average=True` returns the decohered limit in closed form.**  An
astrophysical neutrino arrives with a phase of order 10¹⁵, which nothing
resolves, so the measurement sees

```math
P(\nu_\alpha \to \nu_\beta) = \sum_i |V_{\alpha i}|^2 |V_{\beta i}|^2
```

This is the exact limit rather than an approximation, and it costs one matrix
product.  Magnus checks that the limit applies instead of assuming it: a pair
of eigenvalues in neither limit warns rather than returning a number the
physics does not support.  See
[phase-averaged probabilities](https://mbustama.github.io/Magnus/averaged_probability.html)
and [notebook 10](notebooks/10_magnus_averaged_probability.ipynb).

## Performance

Measured on a laptop, three flavors through PREM, default tolerance 10⁻³:

| Workload | Time |
|---|---|
| One probability, 1 GeV, cos θz = −0.8 | ~2 ms |
| 200-energy scan, one direction | 76 ms |
| 100 × 100 oscillogram | ~2 s |
| `solve_ivp` DOP853 reference, one probability | ~360 ms |

**One thing is under your control.**  If you pass your own Hamiltonian, write
it to accept an array of positions: the engine tries a single vectorized call
and falls back to a Python loop, which is 4.6× slower for a bit-identical
answer.  The trick is broadcasting the potential into a stack of matrices:

```python
def H_func(l):
    l = np.asarray(l, dtype=float)
    VCC = VCC_central*np.exp(-(l/gd.UNIT_KM)/l_scale)   # an array
    return (1.0/energy)*h_vac + VCC[..., None, None]*e00
```

A Hamiltonian that ignores its argument is detected and broadcast already.  The
fallback warns once per session, naming the fix.  More:
[performance](https://mbustama.github.io/Magnus/performance.html).

## Plotting

Magnus ships the figures its own notebooks use, so a plot is one call.
Matplotlib is a dependency, so nothing extra is needed.

| Function | Shape it draws |
|---|---|
| `plot_curves` | curves against any swept variable, with an optional relative-error subpanel |
| `plot_probability_vs_energy`, `..._vs_baseline` | presets over it, with labels, scales and tick spacings fixed |
| `plot_curves_stacked` | small multiples -- one panel per configuration down a shared abscissa |
| `plot_probability_with_profile` | a probability above the matter profile that produced it |
| `plot_probability_with_average` | instantaneous against phase-averaged |
| `plot_biprobability` | the CP ellipse, neutrino against antineutrino |
| `plot_oscillogram` | the two-dimensional map across zenith angle and energy |

Every function returns the `fig` and `ax`, so anything can be overridden.  The
point is consistency: panels that must share limits, scales and ticks are what
drifts when each figure is built by hand.  More:
[plotting](https://mbustama.github.io/Magnus/plotting.html).

## Accuracy and validation

1409 tests, 92% coverage, on Python 3.10 through 3.13, on every push.  The
suite checks the mathematics against independent sources — a separately coded
Bernoulli recursion, closed-form expressions, `solve_ivp` at `rtol=1e-12` — and
asserts identities rather than tolerances where an identity holds, so an
optimization that changed an answer fails rather than passing quietly.

One caveat worth stating plainly: `rtol` and `atol` are a **stopping
criterion**, not an error bound.  They compare two of the code's own
approximations, not the distance to the truth.  See
[what they actually control](https://mbustama.github.io/Magnus/diagnostics.html#what-rtol-and-atol-actually-control).

## Documentation

Full documentation: **[mbustama.github.io/Magnus](https://mbustama.github.io/Magnus/)**.

| | |
|---|---|
| [Mathematical method](https://mbustama.github.io/Magnus/methodology.html) | The expansion derived term by term, and why truncation is unitary |
| [Expansion terms](https://mbustama.github.io/Magnus/expansion_terms.html) | The $\Omega_k$ at any order, and how they are generated |
| [Architecture](https://mbustama.github.io/Magnus/architecture.html) | How the modules fit together, and which layer to call |
| [Engines and dispatch](https://mbustama.github.io/Magnus/engines.html) | Which of the seven engines answers a call, and why |
| [Accuracy and diagnostics](https://mbustama.github.io/Magnus/diagnostics.html) | What each safeguard cannot catch, and every warning explained |
| [Against other codes](https://mbustama.github.io/Magnus/comparison.html) | The full cross-code comparison, with the measurements behind it |
| [Tutorials](https://mbustama.github.io/Magnus/tutorials.html) | All 29 notebooks, with what each is for |
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

Every push runs the test suite on Python 3.10 through 3.13, builds the
documentation with warnings as errors, and lints with Ruff.  The 29 notebooks
execute too, across four parallel shards, with a notebook served from cache
when neither the package nor that notebook has changed.  The badges at the top
report those runs.

## Requirements

`numpy`, `scipy (>= 1.9)`, `joblib` and `matplotlib` — the last so that
`magnus.plotting` works in any installation.  Tests:

```bash
pip install -e '.[test]' && pytest tests/
```

Add `--cov --cov-report=term-missing` for coverage, configured in
`pyproject.toml` so it measures the same thing locally and in CI.

## Changelog

See [CHANGELOG.md](CHANGELOG.md), also
[rendered in the docs](https://mbustama.github.io/Magnus/changelog.html).

## How to cite

Mauricio Bustamante (2026).  *Magnus: neutrino oscillation probabilities via
the Magnus expansion*.  GitHub repository:
https://github.com/mbustama/Magnus.

**Methodology References:**
* Sergio Blanes, Fernando Casas, José A. Oteo & José Ros (2009). The Magnus
  expansion and some of its applications. *Physics Reports, 470*(5-6),
  151-238. [doi:10.1016/j.physrep.2008.11.001](https://doi.org/10.1016/j.physrep.2008.11.001).
* Sergio Blanes, Fernando Casas & Javier Ros (2000). Improved high order
  integrators based on the Magnus expansion. *BIT Numerical Mathematics,
  40*(3), 434-450. [doi:10.1023/A:1022311628317](https://doi.org/10.1023/A:1022311628317).

## License

GNU General Public License v3.0 only (`GPL-3.0-only`); the full text is in
[LICENSE](LICENSE).  You may use, study, modify and redistribute it, including
commercially, provided derivative works carry the same license with source
available.

## Author

Mauricio Bustamante (mbustamante@gmail.com)
