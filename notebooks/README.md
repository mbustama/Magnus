# Magnus notebooks

## 1. Introduction [`1_magnus_introduction.ipynb`]

   Contains examples of basic usage of Mag$`\nu`$s.  Open this if you are looking for a quick introduction, but you will likely have to also look at least at the notebooks `2_magnus_2nu_vacuum_matter.ipynb` and `3_magnus_3nu_vacuum_matter.ipynb` to do something useful.

## 2. Two-neutrino oscillation probabilties [`2_magnus_2nu_vacuum_matter.ipynb`]

   Contains example calculations and plots of oscillation probabilities in a two-neutrino system in vacuum and matter.  The notebook contains the following sections:

   2.1 Probabilities $`2\nu`$: in vacuum

   2.2 Probabilities $`2\nu`$: in matter with constant density

   2.3 Probabilities $`2\nu`$: in matter with varying density

   2.4 Probabilities $`2\nu`$: in matter with castle-wall density potential

   2.5 Probabilities $`2\nu`$: in matter with noisy density potential

   2.6 Probabilities $`2\nu`$: in the Earth

   2.7 Probabilities $`2\nu`$: in the Sun

## 3. Three-neutrino oscillation probabilties [`3_magnus_3nu_vacuum_matter.ipynb`]

   In analogy to the previous notebook, this one contains example calculations and plots of oscillation probabilities in a three-neutrino system in vacuum and matter.  The notebook contains the following sections:

   3.1 Probabilities $`3\nu`$: in vacuum

   3.2 Probabilities $`3\nu`$: in matter with constant density

   3.3 Probabilities $`3\nu`$: in matter with varying density

   3.4 Probabilities $`3\nu`$: in matter with castle-wall density potential

   3.5 Probabilities $`3\nu`$: in matter with noisy density potential

   3.6 Probabilities $`3\nu`$: in the Earth

   3.7 Probabilities $`3\nu`$: in the Sun

## 4. Long-baseline oscillation probabilities [`4_magnus_long_baseline.ipynb`]

   Contains examples of how to compute oscillation probabiltiies between two points on the surface of the Earth (e.g., one neutrino source and one neutrino detector).  This is especially useful to study oscillations in long-baseline neutrino experiments, like DUNE, Super-K, Hyper-K, T2K, and ESS.

## 5. Biprobability plots [`5_magnus_biprobability.ipynb`]

   Contains examples of how to generate a biprobability plot, i.e., a plot of the $3\nu$ oscillation probability ($`\nu_\alpha \nu_\beta`$) vs.~the corresponding anti-neutrino oscillation probability ($`\bar{\nu}_\alpha \bar{\nu}_\beta`$), for different values of the CP violation parameter.

## 6. Oscillograms in Earth [`6_magnus_oscillograms.ipynb`]

   Contains examples of how to generate oscillograms of neutrinos propagating inside Earth, i.e., plots of probability vs.~neutrino direction (expressed as the zenith angle measured from the point of neutrino detection) vs.~neutrino energy. 

## 7. BSM: sterile neutrinos [`7_magnus_bsm_sterile_nu.ipynb`]

   Contains examples of how to compute oscillation probabilities in systems of more than three neutrinos, i.e., containing one or more sterile neutrino.  Specifically, we show examples for 3+1 and 3+2 systems.

## 8. BSM: non-standard interactions [`8_magnus_bsm_nsi.ipynb`]

   Contains examples of how to compute $2\nu$ and $3\nu$ oscillation probabilities in matter if the neutrinos undergo non-standard neutral-current interactions with the medium, conventionally parametrized by the $\epsilon$ parameters.

## 9. BSM: Lorentz-invariance violation [`9_magnus_bsm_liv.ipynb`]

   Contains examples of how to compute $2\nu$ and $3\nu$ oscillation probabilities including an additional effective, energy-dependent Hamiltonian that represents the effect of Lorentz-invariance violation.

## 10. Using Mag$`\nu`$s to expand time-dependent matrix exponentials [`10_magnus_matrix_exponential.ipynb`]

   Contains example of how to use Mag$`\nu`$s to compute the Magnus expansion of the matrix exponential of the exponential of a time-dependent matrix, $A(t)$, such as the ones that appear in the calculation of the time-evolution operator of a time-dependent Hamiltonian, i.e., $\exp(\int_{t_i}^{t_f} A(t))$.