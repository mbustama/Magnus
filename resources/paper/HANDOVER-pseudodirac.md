# Handover: pseudo-Dirac neutrinos in Mag(nu)s

Self-contained brief for a fresh session. Repository: `~/Research/magnus`, branch
`fig2-2x2-and-audit` (or a branch off it). Nothing here has been implemented; the
verification below was done from the outside, against the shipped library.

---

## 1. Why this is wanted

The CPC paper (`resources/paper/main.tex`) derives two averaging formulas and
illustrates neither properly:

- **Eq. (17)**, `equ:averaged_blocks`, the *coherent-block* form
  `<P_ab> = sum_b | sum_{i in b} V*_ai V_bi |^2`, which is what you must use when some
  eigenvalue pairs are still mutually coherent while others have averaged away. It is
  referenced exactly **once** in the whole manuscript, in a sentence saying what it is
  *not*.
- **Eq. (18)**, `equ:averaged_varying`, is referenced four times, but always for the Sun,
  which is adiabatic — so its crossing matrix is the identity and the general form never
  does any work.

Pseudo-Dirac neutrinos are the natural physical case for Eq. (17): tiny splittings that
leave each pair coherent long after the standard splittings have averaged. The paper will
add a panel built on whatever this work produces. Eq. (18) needs a different example (a
supernova shock) and is **out of scope here**.

## 2. The physics to implement

Each Dirac neutrino may in fact be two Majorana states separated by a tiny
`delta_m2`. With `n_active` active flavors and a chosen subset of mass eigenstates
carrying a partner:

- a **paired** mass state `j` becomes two eigenstates `(|nu_j> +- |s_j>)/sqrt(2)`, with
  masses `m_j^2` and `m_j^2 + delta_m2_j`;
- an **unpaired** mass state `j` stays a single eigenstate `|nu_j>`;
- total dimension `n = n_active + (number of pairs)`.

**The user requirement is that pairing is selectable per mass state.** Not every
eigenstate need have a partner: 3 active states with pairs on states 1 and 3 only is a
5-dimensional problem, and must work. Do not hard-code "three pairs".

For active-to-active channels the amplitude is
`A_ab = sum_j U_aj U*_bj exp(-i m_j^2 L / 2E) cos(delta_m2_j L / 4E)`,
and averaging the fast standard phases gives
`<P_ab> = sum_j |U_aj|^2 |U_bj|^2 cos^2(delta_m2_j L / 4E)`.

**That last expression is partial averaging, which this package deliberately does not
provide** (see `avgprob`'s three regimes). Do not implement it. The library's contract is:
fully coherent pair -> block form; fully decohered pair -> naive sum; in between -> refuse
to average and tell the caller to use the un-averaged probability. Build the Hamiltonian
and let the existing engines do the rest.

## 3. What has already been verified (do not redo)

Built by hand, outside the library, at `E = 100` TeV and `L = 100` Mpc, three PMNS states
each paired, NuFIT 6.1 NO parameters:

| check | result |
|---|---|
| 6x6 mixing matrix unitary | yes, to machine precision |
| `osc_prob` on a 6x6 Hamiltonian | runs; rows sum to 1.0 exactly |
| `avgprob.coherence_blocks` at `delta_m2 = 1e-18` | `[[0,1],[2,3],[4,5]]` — correct |
| block form `P_ee` (coherent limit) | **0.56462**, equal to the Dirac 3nu answer |
| naive Eq. (16) on the same six states | **0.28231**, wrong by exactly a factor of two |
| `coherence_report` at `1e-17` and `3e-17` | flags 3 pairs in neither limit |
| thresholds | `COHERENCE_PHASE_THRESHOLD = 0.01`, `DECOHERENCE_PHASE_THRESHOLD = 2*pi` |

Use `0.56462` and `0.28231` as regression targets.

Scale, computed in package units: the pair phase is O(1) at
`delta_m2 ~ 2.6e-17 eV^2` for 100 Mpc / 100 TeV, and `2.6e-16` for 10 Mpc / 100 TeV,
against `2.9e12` and `9.8e13` radians for the standard `Dm21` and `Dm31` phases at the
same point. That separation of scales is the whole physical content.

## 4. What to build

### 4.1 Library functions

Put these in `src/magnus/hamiltonians.py`, following the existing naming
(`hamiltonian_3nu_matter`, `pmns_mixing_matrix`, ...). Suggested surface, adjust if the
module's conventions push another way:

- `pseudo_dirac_mixing_matrix(mixing_matrix, pairs)` — takes an `n_active x n_active`
  mixing matrix and the pairing specification, returns the `n x n` matrix in the
  (active flavors, then sterile partners) x (mass eigenstates) basis.
- `pseudo_dirac_mass_squared(mass_squared, pairs)` — returns the `n` mass-squared values
  with each paired state split.
- `hamiltonian_pseudo_dirac_vacuum(energy, ..., pairs)` — the vacuum Hamiltonian.
- `hamiltonian_pseudo_dirac_matter(...)` — vacuum plus `V_CC(l) * P`.

**Pairing specification.** Make it explicit and order-independent; a mapping from
mass-state index to its splitting is the clearest thing, e.g.
`pairs={0: 1e-18, 2: 4e-18}` meaning states 0 and 2 are paired and state 1 is not.
Validate: indices in range, no duplicates, splittings positive and (warn) small compared
with the standard splittings. Document what happens for an empty mapping — it must
reduce **exactly** to the ordinary Dirac case, and there should be a test asserting
bit-level or machine-precision agreement with `hamiltonian_3nu_vacuum`.

**Matter is the part to get right.** The sterile partners feel no matter potential, so
the projector is not the 3nu `diag(1,0,0)`. With sterile states present the neutral-current
term no longer cancels among the active flavors, exactly as in the 3+N case of Sec. 3 of the
paper: the projector carries `r/2` entries for the steriles, `r` being the
neutron-to-proton ratio. Reuse `matter.matter_potential_projector` rather than writing a
second copy — Sec. 3 records that writing it out by hand a second time is what once
produced a wrong answer. Note also the shipped limitation that the projector takes one
scalar per trajectory while the density varies per layer; it applies here too and should
be mentioned in the docstrings rather than solved.

Check that `H = H_E(E) + V_CC(l) * P` separability still holds for the pseudo-Dirac
Hamiltonian, since the energy-batched engine depends on it. If it does, batched scans get
the speedup for free; if it does not, say so in the docstring.

### 4.2 Docstrings

House style, and it is enforced: numpydoc sections, `.. versionadded:: 1.0.0` (never
`versionchanged` — 1.0.0 has not shipped), and **a runnable Examples block in every new
public function**. Look at `avgprob.coherence_blocks` and `hamiltonians.hamiltonian_3nu_matter`
for the pattern. The examples must actually run: the docs build treats warnings as errors
and CI executes them.

### 4.3 Notebook

**Notebooks are generated.** Edit `notebooks/make_notebooks.py` and never the `.ipynb`.
Add `29_magnus_pseudo_dirac.ipynb`: register it in the `books` dict alongside
`28_magnus_paper_figures.ipynb` (around line 12872), and add its title and one-line
description to the index list near line 15350. Build with

```
python make_notebooks.py --only 29
```

Never run it bare — that regenerates and executes all notebooks and takes the better part
of an hour, and killing it mid-run leaves the others output-stripped.

The user asked for **both regimes**, and the notebook should make the contrast its spine:

1. **Oscillatory.** A baseline and energy where the pair phase is O(1) and the pseudo-Dirac
   oscillation is a real, visible modulation of `P_ee`. Use the un-averaged `osc_prob`;
   this is the regime where the averaged expressions are invalid by construction. Choose
   the configuration so the modulation is resolvable — pseudo-Dirac splittings are small,
   so this is a long-baseline or astrophysical-distance case with a moderate `delta_m2`,
   not a reactor one. Show the Dirac limit (`pairs={}`) on the same axes.
2. **Non-oscillatory.** The high-energy astrophysical case: everything averaged except
   the pairs. Show the block form against the naive sum and recover the factor of two
   above. Then sweep `delta_m2` upward through the three regimes — coherent, undecided,
   decohered — and show the `PhaseAveragingWarning` firing in the middle band. That
   demonstrates a piece of the library's honesty machinery that nothing else exercises.

Also show a case with **only some states paired**, since that is the requirement that
makes the interface more than a toggle.

Explicitly not confined to astrophysical: include at least one terrestrial or solar
configuration, even if only to show that the pseudo-Dirac effect is invisible there
because the pair phase is far below threshold. That is a useful negative result and it
guards the reader against assuming the feature matters everywhere.

### 4.4 Tests

`tests/`, mirroring the existing files. At minimum:

- empty pairing reduces exactly to the Dirac Hamiltonian;
- the mixing matrix is unitary for every pairing pattern, including partial ones;
- dimension is `n_active + len(pairs)`;
- the coherent-limit averaged probability equals the Dirac answer (`0.56462` for the
  configuration above) and the naive sum does not (`0.28231`);
- `coherence_blocks` recovers the intended pairing at a small splitting and splits into
  singletons at a large one;
- a partial pairing (say two of three) behaves as the corresponding mixed case;
- input validation raises rather than silently accepting a bad index or a negative
  splitting.

Note `tests/test_file_tree.py` enforces the repository layout — a new module or notebook
must be registered there or the suite fails.

### 4.5 Docs

`docs/source/` is Sphinx and the build runs with warnings as errors. New public functions
need to appear in `api_reference.rst` / `functions.rst`. Consider a short prose page if the
model needs more explanation than a docstring holds; `averaged_probability.rst` is the
natural neighbour and should probably gain a paragraph pointing at the new capability,
since that page owns the three coherence regimes.

## 5. Traps, all of them paid for already

- **Units.** Distances are in `eV^-1` internally. `earth.distance_traveled_inside_earth`
  returns **km**; multiply by `gd.UNIT_KM`. Passing km where `eV^-1` is wanted returns a
  converged, unitary, wrong answer, and a self-convergence study of exactly `0.0` is the tell.
- **`density_is_of_number_of_electrons=True`** must be passed when the profile callable
  returns an electron *number* density. Omitting it reads the profile as a mass density
  and returns a survival probability that is flat in energy. This cost me a run today.
- **The averaged path warns, one helper does not.** `oscprob` (around line 3606) calls
  `avgprob.coherence_report` and raises `PhaseAveragingWarning`;
  `avgprob.averaged_probabilities_constant_hamiltonian(H, baseline=...)` calls
  `coherence_blocks` instead and answers silently in the undecided band. A separate task is
  already filed for that. Do not build on the silent helper, and if the two are reconciled
  while this work is in flight, prefer the warning path.
- **The paper claims a flavor ceiling.** Sec. 2 says the package ships two, three, four and
  five flavors. Pseudo-Dirac with three pairs is a six-state problem reached through the
  generic callable interface, not a named wrapper. **Do not add a six-flavor wrapper
  family**: it would be ~14 entry points needing fifteen angles and ten phases, and its
  vocabulary (`s14`, `s25`, ...) does not map onto the physical pseudo-Dirac parameters,
  which are the splittings. If any paper sentence needs adjusting because six states are now
  reachable through a named builder, flag it rather than editing the manuscript — the paper
  is being audited in a separate session.
- **Do not regenerate `notebooks/paper_figure_cache.json`.** It holds ~190 s of physics and
  is fingerprint-keyed. Nothing in this work should touch it.

## 6. Definition of done

- New functions with numpydoc docstrings carrying runnable examples.
- `29_magnus_pseudo_dirac.ipynb` generated from `make_notebooks.py`, covering the
  oscillatory case, the astrophysical averaged case, a partially paired spectrum, and one
  configuration where the effect is negligible.
- Tests passing, including the regression targets above.
- Docs building clean with warnings as errors.
- A one-paragraph summary of the final API, so the paper session can write the panel
  against it without reading the source.
