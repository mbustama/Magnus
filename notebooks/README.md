# Mag$`\nu`$s notebooks

Twenty-nine notebooks. Each one is written to stand on its own, so you can open the one
that matches your problem without reading the others. If you are new here, start with
[`01_magnus_introduction.ipynb`](01_magnus_introduction.ipynb), then take whichever of
notebooks 2 to 6 covers the setting you care about.

The notebooks carry their outputs, so they can be read on GitHub without running anything.

## Running them

```bash
pip install -e '.[notebooks]'      # adds nbformat, nbclient, ipykernel, scipy, mpmath, python-ternary
jupyter lab notebooks/
```

## They are generated, so do not edit them by hand

Every notebook in this folder is written by
[`make_notebooks.py`](make_notebooks.py). Edit the generator, then rebuild:

```bash
python notebooks/make_notebooks.py            # all of them
python notebooks/make_notebooks.py --only 28  # one of them
```

Rebuilding executes each notebook, which is why the committed files carry their outputs.
A test, `tests/test_notebooks_match_their_generator.py`, fails if a committed notebook
differs from what the generator writes, because continuous integration runs the committed
files rather than the generator.

Notebook 28 reads its expensive inputs from `paper_figure_cache.json` and recomputes only
what the configuration has changed. Setting `MAGNUS_PAPER_CACHE_ONLY=1` forbids it from
recomputing anything, which is how continuous integration runs it. Notebook 13 keeps its
comparison of the twelve solar models the same way, in `solar_models_cache.json`, and honors
the same variable.

---

## Start here

**1. Introduction** [[`01_magnus_introduction.ipynb`](01_magnus_introduction.ipynb)]
The shortest path to a probability, and the conventions the other notebooks assume.

## Probabilities, standard oscillations

**2. Two-neutrino probabilities** [[`02_magnus_2nu_vacuum_matter.ipynb`](02_magnus_2nu_vacuum_matter.ipynb)]
Two flavors against energy and against direction, in seven settings: vacuum, constant
density, exponential and Gaussian profiles, a periodic castle wall, a noisy profile, then
the Earth and the Sun. Each is checked against a closed form where one exists, which makes
this the notebook to read before trusting any of the others.

**3. Three-neutrino probabilities** [[`03_magnus_3nu_vacuum_matter.ipynb`](03_magnus_3nu_vacuum_matter.ipynb)]
The same seven settings with three flavors and a CP-violating phase. The method does not
change; the Hamiltonian is a 3×3 matrix in the same slot.

**4. Long baselines** [[`04_magnus_long_baseline.ipynb`](04_magnus_long_baseline.ipynb)]
Probabilities between two points on the Earth's surface, the geometry of DUNE, T2K,
Hyper-K and ESS. Give the two sets of coordinates; the chord and its density profile
follow.

**5. Biprobability plots** [[`05_magnus_biprobability.ipynb`](05_magnus_biprobability.ipynb)]
The neutrino probability against the antineutrino one, traced as the CP phase runs over
its range. The area enclosed is the CP violation an experiment tries to measure.

**6. Oscillograms** [[`06_magnus_oscillograms.ipynb`](06_magnus_oscillograms.ipynb)]
Probability across zenith angle and energy at once, the map an atmospheric-neutrino
detector sees. This is the workload that most rewards passing arrays instead of looping.

## New physics in the same slot

**7. Sterile neutrinos** [[`07_magnus_bsm_sterile_nu.ipynb`](07_magnus_bsm_sterile_nu.ipynb)]
Four- and five-flavor systems. Only the dimension of the Hamiltonian and the number of
angles and phases grow.

**8. Non-standard interactions** [[`08_magnus_bsm_nsi.ipynb`](08_magnus_bsm_nsi.ipynb)]
A new matter potential in the slot the standard one occupies, with off-diagonal couplings
the Standard Model does not have.

**9. Lorentz-invariance violation** [[`09_magnus_bsm_liv.ipynb`](09_magnus_bsm_liv.ipynb)]
A term that grows as the vacuum term falls, so standard oscillations switch off at high
energy while this one switches on. The operator dimension sets where that happens.

**29. Pseudo-Dirac neutrinos** [[`29_magnus_pseudo_dirac.ipynb`](29_magnus_pseudo_dirac.ipynb)]
Each mass state split into a nearly degenerate pair. Over an astrophysical baseline the
standard splittings have averaged away while each pair is still part-way through its first
cycle, which is the regime the coherent-block averaging form exists for.

## Real profiles

**13. Tabulated solar models** [[`13_magnus_tabulated_solar_model.ipynb`](13_magnus_tabulated_solar_model.ipynb)]
The twelve standard solar models that ship with the package, taken by name through
`density_profile`. On BS2005-AGS,OP: the difference between the instantaneous probability and
the observable. Averaging a scan does not converge on a solar trajectory; `average=True` gives
the limit in closed form, matching the textbook adiabatic MSW expression to 3e-16. Then all
twelve compared on the averaged observable: they agree to 2.4e-3, and the exponential fit is
off by 0.1.

**14. A supernova shock front** [[`14_magnus_supernova_shock.ipynb`](14_magnus_supernova_shock.ipynb)]
The opposite case. A shock front changes the adiabaticity of the level crossing, so it
moves the conversion probability itself. The averaged observable is wrong by 0.21 and no
averaging repairs it. Mag$`\nu`$s says so every time.

**18. Unusual density profiles** [[`18_magnus_unusual_density_profiles.ipynb`](18_magnus_unusual_density_profiles.ipynb)]
Five profiles with the same mean density and the same length, giving different
probabilities: a neutrino responds to how matter is arranged, not only to how much there
is. Reversing a profile end to end is the one rearrangement that changes nothing, and only
under a condition the notebook makes exact.

## The observable: phase-averaged probabilities

**10. `average=True`** [[`10_magnus_averaged_probability.ipynb`](10_magnus_averaged_probability.ipynb)]
The exact limit an astrophysical measurement reaches, for two through five flavors and for
a custom Hamiltonian, plotted against the oscillation it settles around. Ends with the
pion-decay composition at Earth and the check that refuses the limit where it does not
apply.

**23. When averaging rescues you** [[`23_magnus_when_averaging_helps.ipynb`](23_magnus_when_averaging_helps.ipynb)]
Why notebooks 13 and 14 end in opposite places. An error in the accumulated phase cancels
when integrated over several oscillations; an error in the envelope does not. Both kinds
are injected deliberately on a cheap vacuum probability.

## Controlling the calculation

**11. The matrix exponential** [[`11_magnus_matrix_exponential.ipynb`](11_magnus_matrix_exponential.ipynb)]
How exp(Ω) is computed and why the choice matters: the truncated series is anti-Hermitian,
so its exponential is unitary only if the exponential preserves that. Ends on the
expansion coefficients themselves, derived in exact rational arithmetic at any order.

**12. The `strategy` keyword** [[`12_magnus_adiabatic_hybrid_strategy.ipynb`](12_magnus_adiabatic_hybrid_strategy.ipynb)]
A live comparison of `'auto'`, `'hybrid'` and `'magnus'` for two through five flavors,
each checked against a tight-tolerance `solve_ivp` solution in both runtime and accuracy.

**19. Bring your own Hamiltonian** [[`19_magnus_custom_hamiltonian.ipynb`](19_magnus_custom_hamiltonian.ipynb)]
The actual interface is smaller than the wrapper list suggests: a callable returning a
Hermitian matrix. This covers that contract, the one way of writing it that is worth real
time, and what the Earth entry point declares on your behalf.

**21. What `rtol` and `atol` promise** [[`21_magnus_what_tolerance_means.ipynb`](21_magnus_what_tolerance_means.ipynb)]
They are a stopping criterion, not an error bound: the quantity compared is the difference
between two of the code's own approximations. Both the conservative case and the case
where two grids agree while both are wrong, measured against an independent solution.

**22. Which engine answered** [[`22_magnus_which_engine_answered.ipynb`](22_magnus_which_engine_answered.ipynb)]
Six algorithms in five families, and how `strategy='auto'` picks between them. The second
half is about why a method certifying itself cannot find its own blind spot, and why
running two genuinely different engines needs no oracle.

## Conventions and traps

**15. Antineutrinos** [[`15_magnus_antineutrinos.ipynb`](15_magnus_antineutrinos.ipynb)]
Two things change, the conjugated mixing matrix and the sign of the matter potential.
Apply one and not the other and nothing complains: the result is still a valid probability
matrix, answering a different question. `nubar=True` does both.

**16. Exact versus the textbook approximations** [[`16_magnus_exact_vs_approximations.ipynb`](16_magnus_exact_vs_approximations.ipynb)]
The familiar formulas are exact for the problem they were derived for, and Mag$`\nu`$s
reproduces each to machine precision. What breaks is the substitution: using a
constant-density formula on a profile that is not constant.

**17. Mass ordering and the θ₂₃ octant** [[`17_magnus_ordering_and_octant.ipynb`](17_magnus_ordering_and_octant.ipynb)]
The ordering is carried entirely by the sign of `D31`. The shipped best-fit sets for the
two orderings differ in the ordering, the octant and the CP phase at once, so comparing
them conflates three effects; this isolates one at a time.

**20. Numerical edge cases and the warnings** [[`20_magnus_numerical_edge_cases.ipynb`](20_magnus_numerical_edge_cases.ipynb)]
Which degenerate or empty inputs return a number rather than a `NaN`, and what each of the
nine warning classes means. Some report a bad input, some an expensive choice, some a
condition that was not met and may not matter.

## Cost, and against other codes

**24. Performance** [[`24_magnus_performance.ipynb`](24_magnus_performance.ipynb)]
Three things that make a scan substantially faster without changing any answer, measured
live. Two of the three are worth nothing in the wrong circumstances, which is the more
useful half of the result.

**25. Against other codes** [[`25_magnus_against_other_codes.ipynb`](25_magnus_against_other_codes.ipynb)]
Measured comparisons against whichever other codes are installed, with the two things that
make such a comparison misleading handled explicitly: what problem each code solves, and
what conventions each assumes. Nothing fails if a code is missing.

## Parameters, animation, and the paper

**26. Fourteen years of NuFIT** [[`26_magnus_nufit_evolution.ipynb`](26_magnus_nufit_evolution.ipynb)]
How much of a probability is the parameters. Samples from the published Δχ² profiles of
eighteen NuFIT releases and pushes each sample through Mag$`\nu`$s. The data are the
NuFIT collaboration's; cite them and <http://www.nu-fit.org/> if you use them.

**27. Animated scenes** [[`27_magnus_animations.ipynb`](27_magnus_animations.ipynb)]
Nine scenes, each sweeping one parameter. Four match the ones NuOscProbExact draws, so the
two can be read side by side; the other five animate something a closed-form slab code
does not have. Drawn as filmstrips by default, with the animation behind an opt-in switch.

**28. The paper's figures** [[`28_magnus_paper_figures.ipynb`](28_magnus_paper_figures.ipynb)]
Every figure in `resources/paper/`, from one run. Mag$`\nu`$s's own numbers are computed as
it runs; every other code's are read from `external_*.json`, so none of them has to be
installed. Figures are written to `resources/paper/figs/` as PDF.
