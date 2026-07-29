# Mag$`\nu`$s

[![tests](https://github.com/mbustama/Magnus/actions/workflows/tests.yml/badge.svg)](https://github.com/mbustama/Magnus/actions/workflows/tests.yml)

Code to compute neutrino oscillation probabilities between an arbitrary number
of flavors, for any given Hamiltonian, time-dependent or -independent.

Mag$`\nu`$s computes the neutrino evolution operator via the **Magnus
expansion**: instead of integrating the Schrödinger equation step by step, it
exponentiates truncated time-ordered integrals of the Hamiltonian over a chain
of position slabs.  Any truncation of the Magnus series lives in the Lie
algebra, so the resulting evolution operator is **exactly unitary by
construction** — probabilities are non-negative and sum to one at machine
precision, at any accuracy setting.

## When is Mag$`\nu`$s a win?

Compared to solving the propagation ODE directly (e.g., with an adaptive
Runge–Kutta solver), Mag$`\nu`$s wins when one or more of these apply:

1. **The matter profile varies slowly compared to the oscillation length.**
   A Magnus slab is *exact* for a constant Hamiltonian no matter how many
   oscillation cycles it spans, so the slab size is set by how fast the
   *profile* changes, not by how fast the phase winds.  An ODE solver must
   resolve every oscillation.  For a 1-GeV neutrino crossing the Earth
   (PREM profile), Mag$`\nu`$s needs ~10 slabs plus the ~16 layer crossings,
   versus thousands of right-hand-side evaluations for `solve_ivp` — measured:
   **~2 ms vs ~360–700 ms per probability at comparable accuracy**.

2. **You scan over energy and/or direction** (spectra, oscillograms,
   sensitivity studies).  The Magnus kernel is built from fixed, data-independent
   matrix operations, so slabs — and, for the standard/NSI/LIV Hamiltonians,
   the *entire energy axis* — evaluate as batched NumPy/BLAS calls.  Adaptive
   ODE integration is inherently sequential and cannot share steps across
   energies.  Measured: a 200-energy Earth-crossing scan takes **76 ms**
   (0.4 ms per energy); a 100×100 oscillogram takes **~2 s**.

3. **Unitarity matters more than raw local error** — long baselines, small
   probabilities, CP/T asymmetries.  Runge–Kutta iterates drift off the
   unitary manifold (probability leaks of ~10⁻⁶ at typical tolerances, growing
   with baseline); the Magnus route has no leakage to leak, ever
   (rows sum to 1 to ~10⁻¹⁴).

4. **You want arbitrary physics with no per-model work**: any number of
   flavors, any Hermitian Hamiltonian — sterile neutrinos, NSI,
   Lorentz-invariance violation, or your own matrix function of energy and
   position.

When is it *not* the best tool?  For a single probability at a single energy,
any method is fast enough.  For **extreme accumulated phases** — e.g.,
~10-MeV neutrinos crossing most of the Sun (~10⁴ rad of matter-dominated
phase) — the required slab count can exceed the default caps; Mag$`\nu`$s then
warns (`ToleranceNotAchievedWarning`) instead of failing silently, and you
should raise `max_n_slabs` (or use an adiabatic approximation, which is the
natural method in that regime).  And a tight-tolerance ODE solver remains the
best *reference* for validation — Mag$`\nu`$s's own test suite uses
`scipy.integrate.solve_ivp` at `rtol=1e-12` as ground truth.

## Quick start

```python
import sys
sys.path.extend(['src', 'src/magnus'])  # until pip packaging lands

import numpy as np
import magnus.oscprob as oscprob
import magnus.globaldefs as gd

# --- Three-flavor vacuum probability at 1 GeV over 1000 km ---
energy = 1.0*gd.UNIT_GEV      # [eV]
L = 1000.0*gd.UNIT_KM         # [eV^-1]
P = oscprob.osc_prob_3nu_vacuum(energy, L)   # 3x3 matrix, P[i][j] = P(nu_i -> nu_j)

# --- Energy scan through the Earth (PREM), nu_e -> nu_mu ---
energies = np.logspace(-0.3, 1.3, 200)*gd.UNIT_GEV
P_scan = oscprob.osc_prob_3nu_earth(
    energies, costhz=-0.8, L=2.0*6371.0*0.8*gd.UNIT_KM,
    nu_i=gd.NUE, nu_f=gd.NUMU,
    integration_method='gl')  # Gauss-Legendre: fastest method

# --- Your own Hamiltonian through the Earth ---
# H(energy, l, VCC): VCC is the PREM charged-current potential at position l
h_vac = ...  # your energy-independent vacuum term [eV^2]
def H(energy, l, VCC):
    return (1.0/energy)*h_vac + np.asarray(VCC)[..., None, None]*np.diag([1.0, 0.0, 0.0])

P = oscprob.osc_prob_earth(H, energy, loc_ini='fermilab', loc_fin='homestake')

# --- Fully generic: any Hamiltonian function of position ---
P = oscprob.osc_prob(H_func, t_ini=0.0, t_fin=L)   # H_func(l) -> (d, d) array
```

Oscillation parameters default to the NuFit 6.0 best fit (normal ordering);
pass `s12`, `D31`, `dCP`, ..., or `nubar=True`, to change them.  Find many
worked examples — vacuum, matter, Earth, Sun, oscillograms, biprobability
plots, steriles, NSI, LIV — in the [Jupyter notebooks](notebooks/).

## What Mag$`\nu`$s computes

- **Flavors:** 2ν, 3ν, 4ν (3+1), 5ν (3+2) via dedicated wrappers; any number
  of flavors via the generic `osc_prob`.
- **Environments:** vacuum, constant-density matter, exponentially falling
  density, the Earth (PREM, including chords between named detector sites),
  the Sun; or any density profile you supply.
- **Beyond the Standard Model:** non-standard interactions (NSI),
  CPT-odd Lorentz-invariance violation (LIV), and — via `osc_prob`,
  `osc_prob_earth`, `osc_prob_sun` — arbitrary user Hamiltonians.
- **Neutrinos and antineutrinos**, single energies or arrays, full probability
  matrices or single channels.

## Numerical engine

- Magnus expansion to **order 6**, with the term recursion verified
  term-by-term against Blanes, Casas, Oteo & Ros,
  [Phys. Rep. 470, 151 (2009)](https://doi.org/10.1016/j.physrep.2008.11.001).
- Three integration methods: `'trapezoid'`, `'simpson'` (cumulative
  quadrature), and `'gl'` — **Gauss–Legendre commutator-free integrators**
  of orders 2/4/6 that need only 1/2/3 Hamiltonian evaluations per slab
  (Blanes, Casas & Ros, BIT 40, 434 (2000)).  `'gl'` is the fastest and the
  recommended choice for smooth-per-slab profiles.
- **Exactly unitary exponentials** from the eigendecomposition of the
  (anti-Hermitian) Magnus operator, batched over slabs and energies.
- **Adaptive refinement** to a requested tolerance (`rtol`, `atol`), with a
  phase-based starting slab count, warm starts across scan points, and an
  always-on warning if the refinement caps are hit before convergence.
- **Slab edges aligned with the PREM layer boundaries**, so the high-order
  quadrature never integrates across a density discontinuity.
- **Silent vectorization:** Hamiltonian and density-profile functions that
  accept position arrays are detected and used automatically (with a safe
  scalar fallback); profile evaluations are cached across repeated grids; for
  standard/NSI/LIV scans, the whole energy axis is batched through the kernel.

## Performance

Measured on a laptop, 3ν through the Earth (PREM), default tolerance 10⁻³:

| Workload | Time |
|---|---|
| Single probability (1 GeV, cos θz = −0.8) | ~2 ms |
| 200-energy scan, one direction | 76 ms (0.4 ms/energy) |
| 100 × 100 oscillogram (energy × direction) | ~2 s |
| Reference: `solve_ivp` DOP853, single probability, rtol 10⁻⁶ | ~360 ms |

## Accuracy and validation

The [test suite](tests/) (running in CI on Python 3.10–3.12) validates:

- Magnus terms against an independently coded Bernoulli-number recursion
  (orders 1–6, machine precision) and Gauss–Legendre convergence rates
  (measured error ratios 4/16/64 under slab halving);
- probabilities against closed-form expressions (2ν/3ν vacuum, 2ν
  constant-density matter, for neutrinos and antineutrinos) and against
  high-accuracy `solve_ivp` integrations (asymmetric profiles with complex
  Hamiltonians, full PREM Earth crossings);
- unitarity, slab-ordering, channel conventions, vectorized-vs-scalar
  consistency, and the energy-batched scan against the per-point path.

Requested tolerances are targets for the difference between successive
refinements — the standard adaptive heuristic — not strict global error
bounds; in practice the default 10⁻³ setting delivers ~5 × 10⁻⁴ on Earth
crossings (verified against 10⁻⁷-tolerance references).

## Recent changes (2026-07, `dev`)

The Magnus core and the oscillation-probability chain were audited against
the literature and against exact solutions, then overhauled:

- **Correctness:** fixed dropped coefficients in Ω₄–Ω₆ and a crash at order
  6; fixed the time-ordering of the slab product (affected asymmetric
  profiles with CP violation, e.g., the Sun); fixed the antineutrino sign of
  the matter potential (it was applied twice); fixed the 2ν mass-ordering
  convention (the MSW resonance sat in the wrong channel); implemented
  `osc_prob_earth` / `osc_prob_sun`, which previously were silent stubs;
  removed hand-tuned caps in the Sun wrappers that silently limited solar
  probabilities to ~10⁻² accuracy; NumPy 2 compatibility.
  **Results obtained with earlier versions for antineutrinos in matter, for
  2ν in matter, and for the Sun should be recomputed.**
- **Performance:** ~100× on single probabilities and ~180× on scans
  (see table above), from the Gauss–Legendre integrators, batched slab and
  energy axes, vectorized PREM, probe/profile caching, warm-started
  refinement, and layer-aligned slabs.
- **Testing:** new pytest suite (59 tests) and GitHub Actions workflow.

## Requirements

`numpy`, `scipy (>= 1.9)`, `joblib` — see
[src/requirements.txt](src/requirements.txt).  Run the tests with:

```bash
pip install -r src/requirements.txt pytest && pytest tests/
```

## Repository layout

| Path | Contents |
|---|---|
| [src/magnus/magnus/](src/magnus/magnus/) | Magnus-expansion numerical core |
| [src/magnus/oscprob/](src/magnus/oscprob/) | Oscillation probabilities (main API) + closed forms |
| [src/magnus/hamiltonians/](src/magnus/hamiltonians/) | 2ν–5ν Hamiltonians: vacuum, matter, NSI, LIV |
| [src/magnus/earth/](src/magnus/earth/), [src/magnus/matter/](src/magnus/matter/) | PREM, chord geometry, density profiles, potentials |
| [src/magnus/globaldefs/](src/magnus/globaldefs/) | Units, constants, NuFit parameter sets |
| [notebooks/](notebooks/) | Worked examples and plots |
| [tests/](tests/) | Test suite (pytest; runs in CI) |

## Author

Mauricio Bustamante
