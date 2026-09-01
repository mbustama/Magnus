# The pseudo-Dirac API, for the paper panel

*Written by the session that implemented it, so the panel can be written
without reading the source. Branch `pseudodirac`, module
`src/magnus/hamiltonians/hamiltonians_pseudodirac.py`, notebook 29.*

---

## The paragraph

Pseudo-Dirac spectra are built by four functions in
`magnus.hamiltonians`, all of which take the pairing as a single
`pairs` mapping from mass-state index to that state's splitting
`delta_m2_j` in eV^2 — `pairs={0: 1e-18, 2: 4e-18}` pairs the first and
third mass states and leaves the second alone, giving a five-dimensional
problem, and `pairs={}` (the default) reduces every one of them exactly
to the ordinary Dirac case. `pseudo_dirac_mixing_matrix(U, pairs)` and
`pseudo_dirac_mass_squared(m2, pairs)` return the extended `n x n`
mixing matrix and the `n` mass-squared values, `n = n_active +
len(pairs)`, in a basis whose rows are the `n_active` active flavors
followed by one sterile partner per pair and whose columns run over mass
states in ascending order, each paired state emitting its `+` and `-`
eigenstates consecutively. `hamiltonian_pseudo_dirac_vacuum(energy, U,
m2, pairs, nubar=False)` returns the vacuum Hamiltonian (with
`hamiltonian_pseudo_dirac_vacuum_energy_independent` giving the same
thing without the `1/E`), and `hamiltonian_pseudo_dirac_matter(VCC,
n_active, pairs, ratio_number_neutrons_to_protons=1.0)` returns
`V_CC * P` with the projector taken from `magnus.matter`, so the sterile
partners carry `r/2` and the actives `diag(1,0,0)`; the two are added by
the caller, which keeps `H = H_E(E) + V_CC(l) * P` separable and so
leaves the energy-batched engine applicable unchanged. There is no new
`osc_prob_*` wrapper family and no six-flavor vocabulary: the resulting
`n x n` Hamiltonian is handed to the existing generic callable
interface, and the existing `avgprob` machinery groups it into one
coherent block per pair without being told to.

## Signatures

```python
pseudo_dirac_mixing_matrix(mixing_matrix, pairs=None)
pseudo_dirac_mass_squared(mass_squared, pairs=None)
hamiltonian_pseudo_dirac_vacuum_energy_independent(
    mixing_matrix, mass_squared, pairs=None, nubar=False)
hamiltonian_pseudo_dirac_vacuum(
    energy, mixing_matrix, mass_squared, pairs=None, nubar=False)
hamiltonian_pseudo_dirac_matter(
    VCC, n_active, pairs=None, ratio_number_neutrons_to_protons=1.0)
PseudoDiracSplittingWarning   # a splitting not small vs the standard ones
```

## Three things the panel should not get wrong

1. **The brief's regression targets are wrong.** It gives `P_ee = 0.56462`
   in the coherent-block limit and `0.28231` for the naive sum. The
   correct values for NuFIT 6.1 NO are **0.54814** and **0.27407**. The
   block form must equal `sum_j |U_ej|^4`, which is what the brief itself
   says it should equal, and that sum is 0.54814. The factor of two
   between the two numbers — the actual physical point — is unaffected.

2. **`n_active` must be 3 for the matter term.** `matter_potential_projector`
   hard-codes the charged-current entry on the first of three flavors, so
   `hamiltonian_pseudo_dirac_matter` raises `ValueError` for any other
   value rather than returning a silently mislabeled potential. Vacuum has
   no such restriction.

3. **A visible pair modulation and resolvable standard oscillations cannot
   coexist at physical splittings** — the ratio `delta_m2 / Dm2_31` fixes
   that. Notebook 29's oscillatory panel therefore uses a deliberately
   exaggerated splitting and says so on the figure; the physical-splitting
   case is shown separately, as a null result at a terrestrial baseline.
