# The 8B flux at Earth: spectrum and production profile

Inputs to notebook 28's Figure 5c (`solar_8b_flux.pdf`), the flux of nu_e from 8B decay at
Earth under the long-range L_e - L_mu potential of Figure 5b.

| File | What it is |
|---|---|
| `winter2006_b8_spectrum.csv` | The 8B neutrino spectrum, Table IV of Winter et al., Phys. Rev. C 73, 025503 (2006), arXiv:nucl-ex/0406019, transcribed from the arXiv PDF.  Normalized to 1000 over E in MeV; it integrates to 999.98 and peaks at 6.5 MeV. |
| `b16_gs98_struct_columns.csv` | The columns of `struct_b16_gs98.dat` (Vinyoles et al. 2017) that the production rate needs: radius, temperature, density, and the 1H, 4He and 3He mass fractions, copied as written. |
| `b8_production.py` | Rebuilds the 8B production profile from those columns. |
| `b16_gs98_b8_production.csv` | Its output: 8B production per unit radius in B16-GS98, normalized over r in R_sun. |

## Why the profile is computed

The B16 release published its neutrino production distributions on the authors' page
(Vinyoles et al. 2017, footnote 1).  That page is offline, and the Internet Archive copy could
not be reached from the session that built this.  The structure file was supplied by the
author; it matches the SHA-256 recorded in `src/magnus/data/solar_models/b16_gs98.dat`.

`b8_production.py` computes the rate per volume of 7Be(p,gamma)8B, with 7Be in equilibrium
between 3He(4He,gamma)7Be and its destruction by electron and proton capture: Gamow-peak
rates with the SFII S-factors (Adelberger et al. 2011), Salpeter weak screening, and the 7Be
electron-capture rate of the same review.  Only the shape is used; the total flux is the
model's, 5.46e6 cm^-2 s^-1 (Vinyoles et al. 2017, Table 6).

## Check

On BS2005-AGS,OP, whose 8B distribution Bahcall published
(`../../adversarial_batteries/bs2005agsopflux.csv`), the same calculation gives

| | peak | 10% | median | 90% |
|---|---|---|---|---|
| published | 0.0443 | 0.0238 | 0.0484 | 0.0803 |
| computed | 0.0434 | 0.0242 | 0.0488 | 0.0803 |

in R_sun, and the two curves differ by at most 1.1% of the peak.  The computed 7Be abundance is
12% above Bahcall's at every radius, a normalization that drops out.  For B16-GS98 the profile
peaks at 0.0455 R_sun, with median 0.0490 R_sun and 80% between 0.0245 and 0.0805 R_sun.
`tests/test_paper_long_range_readout.py` reruns the BS05 check.

Run `python b8_production.py` from this folder to reproduce both.
