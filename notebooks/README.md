# Mag$`\nu`$s notebooks

## 1. Introduction [[`01_magnus_introduction.ipynb`](01_magnus_introduction.ipynb)]
   Contains examples of basic usage of Mag$`\nu`$s.  Open this if you are looking for a quick introduction, but you will likely have to also look at least at the notebooks `02_magnus_2nu_vacuum_matter.ipynb` and `03_magnus_3nu_vacuum_matter.ipynb` to do something useful.

## 2. Two-neutrino oscillation probabilties [[`02_magnus_2nu_vacuum_matter.ipynb`](02_magnus_2nu_vacuum_matter.ipynb)]
   Contains example calculations and plots of oscillation probabilities in a two-neutrino system in vacuum and matter.  The notebook contains the following sections:

   ### 2.1 Probabilities $`2\nu`$: in vacuum
   Oscillation proabilities in vacuum vs. energy and vs. direction, validated against the standard probability expression

   ### 2.2 Probabilities $`2\nu`$: in matter with constant density
   Oscillation proabilities in constant-density matter vs. energy and vs. direction, validated against the standard probability expression

   ### 2.3 Probabilities $`2\nu`$: in matter with varying density
   Oscillation proabilities in varying-density matter vs. energy and vs. direction for a matter density profile that falls exponentially with distance and a Gaussian density profile.

   ### 2.4 Probabilities $`2\nu`$: in matter with castle-wall density potential
   Oscillation proabilities in matter vs. energy and vs. direction for periodic, castle-wall density potential with different wall widths.

   ### 2.5 Probabilities $`2\nu`$: in matter with noisy density potential
   Oscillation proabilities in matter vs. energy and vs. direction for density profiles that are not smooth, but rather noisy around a central value.

   ### 2.6 Probabilities $`2\nu`$: in the Earth
   Oscillation proabilities vs. energy and vs. direction for neutrinos propagating inside the Earth.  The example use matter density profile from the popular Preliminary Reference Earth Model (PREM), but others can be used.

   ### 2.7 Probabilities $`2\nu`$: in the Sun
   Oscillation proabilities vs. energy and vs. direction for neutrinos propagating inside the Sun, showing the effect of the MSW resonance.

## 3. Three-neutrino oscillation probabilties [[`03_magnus_3nu_vacuum_matter.ipynb`](03_magnus_3nu_vacuum_matter.ipynb)]

   In analogy to the previous notebook, this one contains example calculations and plots of oscillation probabilities in a three-neutrino system in vacuum and matter.  The notebook contains the following sections (the same descriptions as in the $2\nu$ examples above apply below):

   ### 3.1 Probabilities $`3\nu`$: in vacuum

   ### 3.2 Probabilities $`3\nu`$: in matter with constant density

   ### 3.3 Probabilities $`3\nu`$: in matter with varying density

   ### 3.4 Probabilities $`3\nu`$: in matter with castle-wall density potential

   ### 3.5 Probabilities $`3\nu`$: in matter with noisy density potential

   ### 3.6 Probabilities $`3\nu`$: in the Earth

   ### 3.7 Probabilities $`3\nu`$: in the Sun

## 4. Long-baseline oscillation probabilities [[`04_magnus_long_baseline.ipynb`](04_magnus_long_baseline.ipynb)]

   Contains examples of how to compute oscillation probabiltiies between two points on the surface of the Earth (e.g., one neutrino source and one neutrino detector).  This is especially useful to study oscillations in long-baseline neutrino experiments, like DUNE, Super-K, Hyper-K, T2K, and ESS.

## 5. Biprobability plots [[`05_magnus_biprobability.ipynb`](05_magnus_biprobability.ipynb)]

   Contains examples of how to generate a biprobability plot, i.e., a plot of the $3\nu$ oscillation probability ($`\nu_\alpha \nu_\beta`$) vs.~the corresponding anti-neutrino oscillation probability ($`\bar{\nu}_\alpha \bar{\nu}_\beta`$), for different values of the CP violation parameter.

## 6. Oscillograms in Earth [[`06_magnus_oscillograms.ipynb`](06_magnus_oscillograms.ipynb)]

   Contains examples of how to generate oscillograms of neutrinos propagating inside Earth, i.e., plots of probability vs.~neutrino direction (expressed as the zenith angle measured from the point of neutrino detection) vs.~neutrino energy. 

## 7. BSM: sterile neutrinos [[`07_magnus_bsm_sterile_nu.ipynb`](07_magnus_bsm_sterile_nu.ipynb)]

   Contains examples of how to compute oscillation probabilities in systems of more than three neutrinos, i.e., containing one or more sterile neutrino.  Specifically, we show examples for 3+1 and 3+2 systems.

## 8. BSM: non-standard interactions [[`08_magnus_bsm_nsi.ipynb`](08_magnus_bsm_nsi.ipynb)]

   Contains examples of how to compute $2\nu$ and $3\nu$ oscillation probabilities in matter if the neutrinos undergo non-standard neutral-current interactions with the medium, conventionally parametrized by the $\epsilon$ parameters.

## 9. BSM: Lorentz-invariance violation [[`09_magnus_bsm_liv.ipynb`](09_magnus_bsm_liv.ipynb)]

   Contains examples of how to compute $2\nu$ and $3\nu$ oscillation probabilities including an additional effective, energy-dependent Hamiltonian that represents the effect of Lorentz-invariance violation.

## 10. Using Mag$`\nu`$s to expand time-dependent matrix exponentials [[`10_magnus_matrix_exponential.ipynb`](10_magnus_matrix_exponential.ipynb)]

   Contains example of how to use Mag$`\nu`$s to compute the Magnus expansion of the matrix exponential of the exponential of a time-dependent matrix, $A(t)$, such as the ones that appear in the calculation of the time-evolution operator of a time-dependent Hamiltonian, i.e., $\exp\left(\int_{t_i}^{t_f} A(t)\right)$.