Assorted long examples
----------------------

.. contents:: On this page
   :local:
   :depth: 1

.. _ex-sec-cavity:

A cavity in the Earth’s crust
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _ex-fig-cavity:

.. figure:: ../../../img/paper/cavity.svg
   :width: 95%
   :alt: A cavity in the Earth’s crust

   **A cavity in the Earth’s crust.** Effect of a cavity in the Earth's crust on the :math:`\bar{\nu}_e` survival probability over a baseline of :math:`1\,500` km, computed with Magνs, after Ref.  :cite:p:`Arguelles:2012nw`. *Top*: the crust alone, uniform at 3.3 g cm\ :math:`^{-3}` with :math:`Y_e = 0.5`. *Bottom*: the change that four cavities make to it, each centered on the baseline, with the density and width its label gives. The water cavity carries :math:`Y_e = 0.555`, the other three 0.5. The dashed line marks 49 MeV, where the crust curve peaks. The beam of Ref.  :cite:p:`Arguelles:2012nw` spans 5–150 MeV; below 25 MeV, the oscillation is too rapid to draw. Listing :ref:`A cavity in the Earth's crust <ex-lst-cavity>` computes every curve. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

.. _ex-lst-cavity:

**A cavity in the Earth's crust.** A cavity search in the Earth's crust, after Ref.  :cite:p:`Arguelles:2012nw`: the :math:`\bar{\nu}_e` survival probability over :math:`1\,500` km of Earth's crust, with and without a cavity of different density centered on the baseline. A cavity is a density profile with two walls, whose positions are declared as breakpoints so that the refinement ladder runs inside the cavity. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

.. code-block:: python

   import numpy as np
   import magnus.oscprob as oscprob
   import magnus.matter as matter
   import magnus.globaldefs as gd

   L0 = 1500.0                     # km, source to detector
   RHO_CRUST, YE_CRUST = 3.3, 0.5  # Near PREM's crust
   E = np.linspace(25.0, 150.0, 2500)*gd.UNIT_MEV
   osc = gd.load_nufit_params('NuFIT 6.1')
   kw = dict(osc_params=osc, L0=0.0, nu_i=gd.NUE,
             nu_f=gd.NUE, nubar=True, rtol=1.0e-8,
             atol=1.0e-10,
             density_is_of_number_of_electrons=True)

   def n_e(rho, ye):
       """Electron density [eV^3] of uniform matter."""
       return matter.num_density_e_func(
           0.0, lambda _: rho, electron_fraction=ye,
           ratio_number_neutrons_to_protons=(1.0-ye)/ye,
           density_matter_is_in_g_per_cm3=True)

   NE_CRUST = n_e(RHO_CRUST, YE_CRUST)

   def cavity(rho, ye, w):
       """Crust holding a cavity of width w, centered."""
       d = (L0 - w)/2.0
       ne_in = n_e(rho, ye)

       def profile(l):
           x = np.asarray(l, dtype=float)/gd.UNIT_KM
           return NE_CRUST + (ne_in - NE_CRUST) \
               *((x >= d) & (x <= d + w))

       return profile, np.array([d, d + w])*gd.UNIT_KM

   # An empty crust needs no profile, only a number
   P0 = oscprob.osc_prob_matter_std_potential(
       3, NE_CRUST, E, L0*gd.UNIT_KM, **kw)

   # Water, iron, a mineral deposit, a zone of faults:
   # (rho, Y_e, w) for each
   CASES = [(1.0, 0.555, 250.0), (5.0, 0.5, 250.0),
            (10.0, 0.5, 100.0), (25.0, 0.5, 50.0)]
   dP = []
   for rho, ye, w in CASES:
       profile, walls = cavity(rho, ye, w)
       # The two walls are density jumps: declare them
       P = oscprob.osc_prob_matter_std_potential(
           3, profile, E, L0*gd.UNIT_KM,
           t_breakpoints=walls, **kw)
       dP.append(P - P0)


A region of anomalous density along the path of a neutrino through an otherwise uniform medium changes the oscillation probability measured at the far end. Using neutrinos to explore the interior of the Earth was proposed decades ago :cite:p:`DeRujula:1983ya,Ermilova:1986ph`; oscillations in matter offer one way to do it (see, e.g., Ref. :cite:p:`Winter:2006vg`). Reference :cite:p:`Arguelles:2012nw` applied this to a search for cavities in the Earth's crust. We use Magνs to implement the setup of Ref.  :cite:p:`Arguelles:2012nw`: a low-energy :math:`\bar{\nu}_e` beam crossing :math:`1\,500` km of crust, with the survival probability measured from 5 to 150 MeV. The example is illustrative; we do not claim that such an experiment is feasible.

Listing :ref:`A cavity in the Earth's crust <ex-lst-cavity>` sets this up. The crust is uniform at 3.3 g cm\ :math:`^{-3}` with :math:`Y_e = 0.5`, a round value close to the 0.4952 of the PREM crust. A cavity is a slab of different density, centered on the baseline. We show four cavities, with the densities and labels of Ref.  :cite:p:`Arguelles:2012nw`: water at 1 g cm\ :math:`^{-3}` with :math:`Y_e = 0.555`, an iron-banded formation at 5 g cm\ :math:`^{-3}`, a mineral deposit at 10 g cm\ :math:`^{-3}`, and a zone of seismic faults at 25 g cm\ :math:`^{-3}`; the last two are denser than any natural rock, and serve to show how the effect grows with density. Their widths along the path range from 250 to 50 km, the densest being the narrowest. Each wall of a cavity is a density jump, so we pass the positions of both walls to the scenario function via ``t_breakpoints``; this places a slab edge on each wall, so that no slab straddles a jump (:doc:`/recipes`).

Figure :ref:`A cavity in the Earth’s crust <ex-fig-cavity>` shows the resulting probabilities. Without a cavity, the survival probability peaks at 0.9997 near 49 MeV, where the oscillation driven by :math:`\Delta m^2_{21}` completes half a cycle (the peak would reach exactly 1 if :math:`\theta_{13}` were zero). A cavity shifts this peak along the energy axis, leaving its shape almost unchanged. A denser cavity adds column density and pushes the peak to higher energy; a lighter one removes column density and pulls it to lower energy. Therefore, just below the peak of the crust curve, a denser cavity lowers the probability, and just above it, raises it; a lighter cavity does the opposite. Every curve of the change crosses zero between the peak of the crust curve and that of the cavity curve, between 48.2 and 49.5 MeV.

Figure :ref:`A beam swept across a buried body <ex-fig-cavity-sweep>` sweeps the beam across a buried body. The source stays fixed and the baseline keeps its length of :math:`1\,500` km; only the direction of the beam changes, by an angle :math:`\alpha`. The body is a sphere of radius 125 km at the density of the mineral deposit, 10 g cm\ :math:`^{-3}`, centered 750 km from the source along the baseline at :math:`\alpha = 0`. (Ref.  :cite:p:`Arguelles:2012nw` uses elliptical cavities, but a sphere makes the same point.) The beam crosses the body only for :math:`\lvert \alpha \rvert < 9.59^\circ`; the width it crosses is largest, 250 km, at :math:`\alpha = 0`, and falls to zero at the edges of that band.

The change in probability follows that width: it is largest near :math:`\alpha = 0` and fades toward the edges of the band. Outside the band, the beam misses the body, the profile is the uniform crust, and the change is exactly zero. The sweep is what separates a large, light body from a small, dense one. Along a single direction, the probability depends mainly on the excess column density, i.e., the density contrast times the width, so density and width cannot be told apart; the angular width of the band fixes the size of the body independently. Listing :ref:`A beam swept across a buried body <ex-lst-cavity-sweep>` computes the map.

.. _ex-fig-cavity-sweep:

.. figure:: ../../../img/paper/cavity_sweep.svg
   :width: 95%
   :alt: A beam swept across a buried body

   **A beam swept across a buried body.** *Top*: the geometry. The beam leaves the source at an angle :math:`\alpha` and reaches a detector on the dashed arc, :math:`1\,500` km away. The two dashed straight lines are the beams tangent to the body, at :math:`\alpha = \pm 9.59^\circ`. The beam crosses a width :math:`w` of the body. *Bottom left*: width against angle, i.e., the cavity's silhouette. *Bottom right*: change in probability over energy and angle. The dashed vertical line marks the :math:`49` MeV of Figure :ref:`A cavity in the Earth’s crust <ex-fig-cavity>`. Listing :ref:`A beam swept across a buried body <ex-lst-cavity-sweep>` computes the map. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

.. _ex-lst-cavity-sweep:

**A beam swept across a buried body.** The map of Figure :ref:`A beam swept across a buried body <ex-fig-cavity-sweep>`, which reuses ``n_e``, ``NE_CRUST``, ``L0``, and ``kw`` from Listing :ref:`A cavity in the Earth's crust <ex-lst-cavity>`. For each beam angle, ``crossing`` finds where the chord enters and leaves the body. Each call uses these two points twice: as the edges of the body in the density profile, and as ``t_breakpoints``. Angles at which the beam misses the body are skipped, since there the probability equals the crust-only ``P0``. A cavity is a density profile with two walls, whose positions are declared as breakpoints, so that no slab straddles a wall. See :ref:`ex-sec-cavity` for details.

.. code-block:: python

   import numpy as np

   BODY_R, BODY_D0 = 125.0, 750.0   # km: size, position
   NE_BODY = n_e(10.0, 0.5)         # a mineral deposit
   E = np.linspace(25.0, 150.0, 400)*gd.UNIT_MEV
   ALPHA = np.linspace(-15.0, 15.0, 220)

   def crossing(alpha_deg):
       """Entry and exit points in the body [km]."""
       a = np.radians(alpha_deg)
       miss = abs(BODY_D0*np.sin(a))
       if miss >= BODY_R:
           return None              # the beam misses it
       half = np.sqrt(BODY_R**2 - miss**2)
       mid = BODY_D0*np.cos(a)
       return mid - half, mid + half

   P0 = oscprob.osc_prob_matter_std_potential(
       3, NE_CRUST, E, L0*gd.UNIT_KM, **kw)

   dP = np.zeros((len(E), len(ALPHA)))
   for j, alpha in enumerate(ALPHA):
       seg = crossing(alpha)
       if seg is None:
           continue                 # no body, no change
       lo, hi = seg

       def body(l, lo=lo, hi=hi):
           x = np.asarray(l, dtype=float)/gd.UNIT_KM
           return NE_CRUST + (NE_BODY - NE_CRUST) \
               *((x >= lo) & (x <= hi))

       # The walls move with the angle, so every call
       # declares its own pair
       dP[:, j] = oscprob.osc_prob_matter_std_potential(
           3, body, E, L0*gd.UNIT_KM,
           t_breakpoints=np.array(seg)*gd.UNIT_KM,
           **kw) - P0


.. _ex-sec-geoneutrinos:

Geoneutrinos
~~~~~~~~~~~~

.. _ex-fig-geoneutrinos:

.. figure:: ../../../img/paper/geoneutrinos.svg
   :width: 95%
   :alt: Geoneutrino geometry at Borexino

   **Geoneutrino geometry at Borexino.** The geometry of a geoneutrino measurement at Borexino, at Gran Sasso, 1.4 km underground. A quarter of the Earth is cut away to expose the PREM layers of Figure :ref:`The Earth’s density profile <ex-fig-prem>`. Three chords reach production points beyond the local crust: in the far crust, 20 km deep and 3 400 km away; in the mantle, 1 000 km deep and 4 300 km away; and at the base of the mantle, 2 800 km deep and 7 300 km away, on a path through the outer core. The inset shows the local crust, to scale in distance and stretched in depth, with the crust layers of PREM and a production point 10 km deep and 100 km away. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

Geoneutrinos are the :math:`\bar{\nu}_e` emitted in the decay chains of uranium-238 and thorium-232 inside the Earth  :cite:p:`Bellini:2013wsa`. Inverse beta decay can detect only those above 1.8 MeV. The most energetic geoneutrinos, from the uranium chain, reach 3.3 MeV, so the detectable spectrum spans 1.8–3.3 MeV. The sources are the crust and the mantle. Uranium and thorium are lithophile elements: they concentrate in silicate rock, and are expected to be absent from the metallic core. The core still enters the calculation, though, as the medium that neutrinos from the far side of the mantle cross. For a detector at or near the surface, the production points lie from a few km to a full Earth diameter away. As an example, we take Borexino :cite:p:`Borexino:2008gab,Borexino:2019gps` as the detector, at Gran Sasso, 1.4 km underground.

Figure :ref:`Geoneutrino geometry at Borexino <ex-fig-geoneutrinos>` shows the setup: four production points illustrate the kinds of path a geoneutrino takes to Borexino. One lies in the local crust, 10 km deep and 100 km away. One lies in the far crust, 20 km deep and 3 400 km away, on a path that dips into the upper mantle. One lies in the mantle, 1 000 km deep and 4 300 km away. The last lies at the base of the mantle, 2 800 km deep and 7 300 km away, on a path through the outer core.

.. _ex-fig-geoneutrino-energy:

.. figure:: ../../../img/paper/geoneutrino_energy.svg
   :width: 95%
   :alt: Geoneutrino survival against energy

   **Geoneutrino survival against energy.** Survival probability of :math:`\bar{\nu}_e` geoneutrinos reaching Borexino, from the four production points of Figure :ref:`Geoneutrino geometry at Borexino <ex-fig-geoneutrinos>`. In every panel, the fast oscillation of the pair split by :math:`\Delta m^2_{31}` is resolved, riding as a fine ripple on the slower oscillation of the pair split by :math:`\Delta m^2_{21}`. This is visible clearly only in the near-production case (*top panel*). In the three distant-production cases (*bottom three panels*), the line marks the phase average in vacuum, 0.548; matter raises it along these chords by 0.3–0.8%, less than the width of the line. Listing :ref:`Geoneutrino survival against energy <ex-lst-geoneutrinos>` computes the four panels. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

Figure :ref:`Geoneutrino survival against energy <ex-fig-geoneutrino-energy>` shows the survival probability against energy for these production points. Two oscillations appear: a slow one, from the pair of mass eigenstates split by :math:`\Delta m^2_{21}`, and a fast one, from the pairs split by :math:`\Delta m^2_{31}`. From the local crust, the slow oscillation completes less than one cycle across the window. From the three distant points, both complete tens to thousands of cycles, far more than a detector resolves in energy, so a measurement of a distant reservoir sees only the average probability. Within a few hundred kilometers of the detector, the slow oscillation does not average out across the window. For such nearby production points, the full probability must be used instead of the average; Figure :ref:`Where the geoneutrino flux comes from <ex-fig-geoneutrino-flux>` uses it out to about 330 km.

For the three distant production points, Figure :ref:`Geoneutrino survival against energy <ex-fig-geoneutrino-energy>` marks the phase-averaged probability in vacuum, 0.548. Matter changes it little. The passage through the Earth is adiabatic: even the largest jump in electron density, at the core-mantle boundary, changes the mixing angle in matter by only about a percent at these energies, so the neutrino stays in the eigenstate in which it was produced. The average then depends only on the eigenstates at the two ends of the path [:ref:`averaged limit on a varying profile <avg-varying>`, with :math:`P^{\rm cross}` the identity; :doc:`/averaged_probability`]. At both ends, the matter potential is at most 2% of :math:`\Delta m^2_{21}/2E`, so those eigenstates are nearly the vacuum ones. As a result, matter raises the average by only 0.3–0.8%, less than the width of the line in Figure :ref:`Geoneutrino survival against energy <ex-fig-geoneutrino-energy>`, and the vacuum value serves our calculation. The snippet below checks this along the mantle chord, with the density of PREM and the electron fraction of the mantle, using ``average=True`` (:ref:`ex-sec-average-keyword`):

.. code-block:: python

   from magnus.earth import (
    distance_traveled_inside_earth as chord,
    earth_radial_distance_from_depth as radius,
    density_matter_func_prem as prem,
    Y_E_MANTLE_PREM)

   osc = gd.load_nufit_params('NuFIT 6.1')

   # The mantle chord: cosine of the zenith
   # angle at the detector, and the depths of
   # the production point and of Borexino [km]
   c = -0.552
   geo = dict(source_depth=1000.0,
              detector_depth=1.4)

   # Length of the chord [km]
   L = chord(c, **geo)

   # Density along the chord [g/cm^3]: PREM,
   # with its ocean layer replaced by rock
   def rho(l):
       # Radius at distance l along the chord
       r = radius(c, l/gd.UNIT_KM, **geo)
       return prem(r,
                   density_matter_ocean=2.65)

   E = 2.5*gd.UNIT_MEV

   # Phase-averaged survival probability in
   # matter, at electron fraction of the mantle
   P_m = oscprob.osc_prob_matter_std_potential(
    3, rho, E, L*gd.UNIT_KM, osc,
    average=True, nubar=True,
    nu_i=gd.NUE, nu_f=gd.NUE,
    electron_fraction=Y_E_MANTLE_PREM,
    density_matter_is_in_g_per_cm3=True)
   # 0.551

   # The same, in vacuum
   P_v = oscprob.osc_prob_3nu_vacuum(E,
    L*gd.UNIT_KM, average=True, nubar=True,
    nu_i=gd.NUE, nu_f=gd.NUE)
   # 0.548


Listing :ref:`Geoneutrino survival against energy <ex-lst-geoneutrinos>` computes the panels of Figure :ref:`Geoneutrino survival against energy <ex-fig-geoneutrino-energy>` with the Earth wrapper, ``osc_prob_3nu_earth``. Each production point enters through its depth and the cosine of its zenith angle at the detector. Every chord crosses PREM boundaries—six on the far-crust chord, nine on the chord through the core—and the wrapper internally places a slab edge on each. Since Gran Sasso is a continental site, ``density_matter_ocean`` replaces the ocean layer of PREM with rock.

(When a source is buried, the wrapper automatically places it on the far side of the chord, where the chord first reaches the source depth. However, the local-crust production point, 10 km deep and 100 km away, lies instead on the near side: along its chord, a depth of 10 km is first reached about 1 100 km from the detector. To overcome this, the last call in Listing :ref:`Geoneutrino survival against energy <ex-lst-geoneutrinos>` swaps the two ends of the path: it treats the detector as the source, 1.4 km deep, and the production point as the detector, 10 km deep, with the zenith angle measured at the production point. The result is unchanged, because the survival probability is the same along a path and along its reverse.)

.. _ex-lst-geoneutrinos:

**Geoneutrino survival against energy.** The four panels of Figure :ref:`Geoneutrino survival against energy <ex-fig-geoneutrino-energy>`: the survival probability of a geoneutrino reaching Borexino from the far crust, the mantle, the base of the mantle, and the local crust, across the detectable window, on one energy grid fine enough to resolve the pair split by :math:`\Delta m^2_{31}` on the longest chord. The call declares the PREM boundaries on each chord by itself. The ocean layer of PREM is replaced by rock. The local production point is run in reverse, from the detector. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

.. code-block:: python

   import numpy as np
   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

   # Twelve energies per cycle of the pair split by
   # Dm31^2 on the longest chord, 7,300 km
   E = 1/np.linspace(1/1.8, 1/3.3, 22500)*gd.UNIT_MEV

   # Borexino, 1.4 km underground.  A production point
   # enters as its depth and the cosine of its zenith
   # angle at the detector; the call declares the PREM
   # boundaries on the chord by itself.  The oscillation
   # parameters are the defaults, NuFIT 6.1.
   points = {'far crust': (20.0, -0.272),    # 3,400 km
             'mantle': (1000.0, -0.552),     # 4,300 km
             'core': (2800.0, -0.872)}       # 7,300 km
   P = {}
   for name, (depth, costhz) in points.items():
       P[name] = oscprob.osc_prob_3nu_earth(E,
        costhz=costhz, source_depth=depth*gd.UNIT_KM,
        detector_depth=1.4*gd.UNIT_KM, nubar=True,
        nu_i=gd.NUE, nu_f=gd.NUE,
        density_matter_ocean=2.65)

   # The local point, 10 km deep and 100 km away,
   # run in reverse, from the detector (see text)
   P['local'] = oscprob.osc_prob_3nu_earth(E,
    costhz=0.078, source_depth=1.4*gd.UNIT_KM,
    detector_depth=10.0*gd.UNIT_KM, nubar=True,
    nu_i=gd.NUE, nu_f=gd.NUE,
    density_matter_ocean=2.65)


.. _ex-fig-geoneutrino-flux:

.. figure:: ../../../img/paper/geoneutrino_flux.svg
   :width: 95%
   :alt: Where the geoneutrino flux comes from

   **Where the geoneutrino flux comes from.** Where the detectable geoneutrino flux at Borexino comes from and what oscillations leave of it. *Top*: the survival probability at 2.5 MeV against the distance to a production point 10 km deep, across the local crust, with the pair split by :math:`\Delta m^2_{31}` resolved; the horizontal line is the phase average. *Bottom*: the distribution of the flux of Eq. :eq:`ex-equ-geoflux` in the distance to the production point, per unit :math:`\log_{10} L` and as a fraction of the total, for a spherically symmetric Earth with a 35-km crust holding 7 TW of radiogenic power, a mantle with the abundances of the geochemical model of Ref. :cite:p:`Bellini:2013wsa`, and no uranium or thorium in the core. The unoscillated curve has :math:`\langle P_{\bar\nu_e \to \bar\nu_e} \rangle = 1`. The oscillated curve folds in the survival probability averaged over the detectable window with equal weight in energy; its cumulative fraction is on the right axis. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

Figure :ref:`Where the geoneutrino flux comes from <ex-fig-geoneutrino-flux>` shows how the detectable flux is distributed over the distance to the production point and how much of it survives the oscillations. A flux measurement weighs the curves of Figure :ref:`Geoneutrino survival against energy <ex-fig-geoneutrino-energy>` by that distribution. The Earth model chosen for the figure is the simplest one that carries both reservoirs: a spherically symmetric crust 35-km thick holding 7 TW of radiogenic power at a thorium-to-uranium ratio of 4.5, a mantle with the abundances of the geochemical model of Ref. :cite:p:`Bellini:2013wsa`, no uranium or thorium in the core, and the density of PREM throughout. The flux at the detector, located at :math:`\mathbf{r}_{\rm det}`, is

.. math::
   :label: ex-equ-geoflux

   F = \int_\oplus d^3r \, \frac{\varepsilon(\mathbf{r})}{4 \pi L^2} \, \langle P_{\bar\nu_e \to \bar\nu_e} (L) \rangle \,,

where :math:`L \equiv \lvert \mathbf{r} - \mathbf{r}_{\rm det} \rvert` and :math:`\varepsilon(\mathbf{r})` is the number of detectable :math:`\bar{\nu}_e` emitted per unit volume and time at position :math:`\mathbf{r}`, fixed by the density of PREM and the abundances above. The survival probability from that point, :math:`\langle P_{\bar\nu_e \to \bar\nu_e}(L) \rangle`, is the mean over the detectable window, 1.8–3.3 MeV. In spherical coordinates centered on the detector, :math:`d^3r = L^2 \, dL \, d\Omega` and the :math:`L^2` cancels: the production points between :math:`L` and :math:`L + dL` contribute :math:`dL / 4\pi` times the emissivity integrated over the directions in which that shell lies inside the Earth. The lower panel of Figure :ref:`Where the geoneutrino flux comes from <ex-fig-geoneutrino-flux>` shows this distribution per unit :math:`\log_{10} L`, normalized to the total flux without oscillations, so that the area under a curve between two distances is the share of the flux from between them. Without oscillations, the crust supplies 61% of the detectable flux, and the region within 350 km of the detector supplies 30%, all of it close enough that the full probability must be used instead of the average.

The oscillated curve is the unoscillated one multiplied by :math:`\langle P_{\bar\nu_e \to \bar\nu_e}(L) \rangle`, the survival probability averaged over the detectable window. Up to 331 km, we compute the probability in full, at each distance on the grid, with the last call of Listing :ref:`Geoneutrino survival against energy <ex-lst-geoneutrinos>`, and then average it over the window. Beyond 331 km, we use the average in vacuum instead, which differs from the full result by 0.01 at 331 km and by less farther out. The upper panel of Figure :ref:`Where the geoneutrino flux comes from <ex-fig-geoneutrino-flux>` shows the full probability at 2.5 MeV: the slow oscillation has a period of 82 km, and the fast one, with a period of 2.5 km, appears as a ripple on it. Close to the detector, the window holds too few cycles of the slow oscillation to average it out, so the oscillated curve in the lower panel ripples out to a few hundred kilometers; farther out, it is smooth. Over the whole Earth, oscillations leave 54% of the detectable flux.

Both curves in the lower panel of Figure :ref:`Where the geoneutrino flux comes from <ex-fig-geoneutrino-flux>` flatten at large distances, for a geometric reason. The production points at a distance :math:`L` from the detector lie on a sphere of radius :math:`L` centered on it. The part of this sphere at a distance between :math:`r` and :math:`r + dr` from the center of the Earth has area :math:`2 \pi L r \, dr / R_\oplus`. In Eq. :eq:`ex-equ-geoflux`, this area is divided by :math:`4 \pi L^2`, and on the logarithmic axis of the panel it is multiplied by :math:`L`, so the result does not depend on :math:`L`. Hence, a layer of the Earth contributes the same flux per decade of distance as long as the sphere reaches every depth of that layer. For the crust, this holds from 35 km out to a full Earth diameter; for the mantle, from about 2 900 to 9 900 km. The total stays flat over that range because the core, which holds no uranium or thorium in this model, adds nothing to it.

.. _ex-sec-solar-tomography:

Neutrino tomography of the Sun
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _ex-fig-solar-adiabaticity:

.. figure:: ../../../img/paper/solar_adiabaticity.svg
   :width: 95%
   :alt: Adiabaticity along a solar chord

   **Adiabaticity along a solar chord.** *Top*: the setup. An isotropic flux of :math:`\nu_e` crosses the Sun and continues to Earth, where it may oscillate. The impact parameter :math:`b` is marked on one ray. *Bottom:* Adiabaticity of a neutrino trajectory through the Sun against neutrino energy, computed with Magνs on the B16-GS98 solar density profile  :cite:p:`Vinyoles:2016djt`, tabulated up to the surface, for four impact parameters, each defining a chord through the Sun. The parameter :math:`\gamma_{\rm max}` is the adiabaticity parameter of :doc:`/adiabatic_strategy`, evaluated at its maximum value along the chord and over the three pairs of eigenvalue levels. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

Neutrinos made outside the Solar System arrive from all directions, so a share of them crosses the Sun on its way to Earth. At around 10 MeV, they come mainly from the diffuse supernova neutrino background; from TeV to PeV, they are part of the high-energy astrophysical flux. Below, we compute what that does to their oscillations, using the density profile of the B16-GS98 solar model  :cite:p:`Vinyoles:2016djt` of :ref:`ex-tab-solar-models`.

At the energies of this section, up to 1 PeV, a neutrino from a source outside the Solar System has traveled far enough for its oscillations to average out on the way (Figure :ref:`The averaged flavor composition settles once the source is far enough. <ex-fig-astro-baseline>`), so it reaches the Sun as an incoherent mixture of :math:`\nu_1`, :math:`\nu_2`, and :math:`\nu_3`. Therefore, the quantity to compute is the phase-averaged probability, propagated through the varying profile of the Sun as described in :doc:`/averaged_probability`. Where the passage is adiabatic, the neutrino stays in the instantaneous eigenstates of :math:`\mathbb{H}`; where it is not, Magνs crosses the region with a Magnus patch. If no region is non-adiabatic, the result is :ref:`averaged limit on a varying profile <avg-varying>`, with :math:`P^{\rm cross}` the level-crossing matrix between the eigenbases at entry and exit. Both ends of the path lie in vacuum, where the eigenvectors of :math:`\mathbb{H}` form the vacuum mixing matrix, so :math:`\mathbb{V}(l_0) = \mathbb{V}(l_1) = \mathbb{U}` in :ref:`averaged limit on a varying profile <avg-varying>`. If the passage is adiabatic, :math:`P^{\rm cross}` is the identity, and :ref:`averaged limit on a varying profile <avg-varying>` reduces to the vacuum :ref:`averaged limit <avg-limit>`: the Sun leaves no trace. Whether the passage is adiabatic depends on the energy, through the adiabaticity parameter of :doc:`/adiabatic_strategy`.

This may seem at odds with the MSW effect, which changes solar neutrinos strongly (:ref:`ex-sec-sun`). The difference lies in where the path starts. A solar neutrino is born deep inside the Sun, at high density, as an eigenstate of :math:`\mathbb{H}` in matter; an adiabatic passage outward turns it into a single vacuum mass eigenstate, which is a different state. A neutrino crossing the Sun from outside starts and ends in vacuum, so an adiabatic passage returns it to the same mass eigenstate it entered as.

Figure :ref:`Adiabaticity along a solar chord <ex-fig-solar-adiabaticity>` shows how adiabatic the passage is, from MeV to PeV, at three flavors, for four impact parameters. The adiabaticity parameter :math:`\gamma_{jk}(l)` of :doc:`/adiabatic_strategy` is defined for each pair of eigenvalues of :math:`\mathbb{H}`, 1-2, 1-3, and 2-3, at each point :math:`l` of the path; it is large where the two eigenvalues come close. The figure plots its largest value along the chord and over the three pairs, :math:`\gamma_{\rm max}`. Where :math:`\gamma_{\rm max} > 1`, the passage is not adiabatic: somewhere along the chord, the neutrino moves between eigenstates of :math:`\mathbb{H}`. The ``adiabatic`` module of Magνs computes :math:`\gamma_{\rm max}`:

.. code-block:: python

   import numpy as np
   import magnus.adiabatic as adiabatic
   import magnus.hamiltonians as ham
   import magnus.matter as matter
   import magnus.globaldefs as gd
   import magnus.solarmodels as solarmodels

   osc = gd.load_nufit_params('NuFIT 6.1')
   R = gd.SUN_RADIUS*gd.UNIT_KM
   prof = solarmodels.electron_density_profile
   ne_sun = prof('B16-GS98')
   b = 0.3*R                # impact parameter
   half = np.sqrt(R**2 - b**2)
   E = 100.0*gd.UNIT_GEV
   hv = ham.hamiltonian_3nu_vacuum(E, **osc)

   def ne(l):
       """Electron density on the chord."""
       r = np.sqrt((l - half)**2 + b**2)
       return ne_sun(r)

   def H(l):
       """Hamiltonian at l along the chord."""
       h = np.array(hv, dtype=complex)
       h[0, 0] += matter.VCC_func(l, ne)
       return h

   info = {}
   adiabatic.find_nonadiabatic_windows(
       H, 0.0, 2*half, n_probe=20000,
       info=info)
   info['gamma_max']        # 13.7580


Here, ``ne_sun`` is the B16-GS98 table, interpolated in its logarithm, ``ne`` evaluates it along the chord, and ``hv`` is the vacuum Hamiltonian at this energy. The value shown is for a :math:`100`-GeV neutrino crossing at an impact parameter of :math:`b = 0.3\,R_\odot`.

The parameter :math:`\gamma_{\rm max}` grows with energy because the vacuum splitting :math:`\Delta m^2/2E` shrinks while the matter potential stays fixed. The largest value comes from the 1-2 pair, near the point where its two eigenvalues come closest: the 1-2 resonance, where the matter potential equals :math:`\Delta m^2_{21}\cos 2\theta_{12}/2E`. This resonance occurs at a single radius in the Sun, set by the energy alone, so every chord reaches its :math:`\gamma_{\rm max}` at that radius. As the energy grows, that density is reached further out: :math:`\gamma_{\rm max}` sits at :math:`0.97\,R_\odot` at :math:`100` GeV, at :math:`0.990\,R_\odot` at 1 TeV, and at :math:`0.998\,R_\odot` at 100 TeV. By 1 PeV, it reaches :math:`5\times10^6` on the diameter.

.. _ex-fig-solar-tomography:

.. figure:: ../../../img/paper/solar_tomography.svg
   :width: 95%
   :alt: The Sun in the electron-neutrino channel

   **The Sun in the electron-neutrino channel.** The Sun seen face-on in the :math:`\nu_e` survival channel, computed with Magνs on the B16-GS98 solar density profile  :cite:p:`Vinyoles:2016djt`. The line of sight runs into the page, so each point of the disk is an impact parameter and the neutrino crosses the whole Sun along it. The color is the phase-averaged :math:`P_{\nu_e \to \nu_e}` at exit, for a relative energy spread of :math:`10%`, averaged over the area of each pixel. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

Figure :ref:`The Sun in the electron-neutrino channel <ex-fig-solar-tomography>` shows the averaged survival probability across the face of the Sun, at selected energies. The line of sight runs into the page, so each point in the disk fixes the impact parameter of one trajectory, along which the neutrino crosses the whole Sun. For one impact parameter, the probability is computed via

.. code-block:: python

   import magnus.oscprob as oscprob

   # The chord above: b = 0.3 R_sun
   P = oscprob.osc_prob_matter_std_potential(
       3, ne, E, 2*half, average=True,
       average_initial_state='decohered',
       osc_params=osc, L0=0.0, nu_i=gd.NUE,
       nu_f=gd.NUE,
       density_is_of_number_of_electrons=True)
   P                        # 0.283618


Here, ``average=True`` returns the phase average along the chord, and ``average_initial_state='decohered'`` tells Magνs that the neutrino reaches the Sun decohered, i.e., as :math:`\nu_1`, :math:`\nu_2`, or :math:`\nu_3` (:ref:`ex-sec-average-keyword`). The default, ``'flavor'``, would instead treat it as a :math:`\nu_e` produced at the edge of the Sun, which keeps the interference between :math:`\nu_1` and :math:`\nu_2`; their phase across the Sun is only about 1 rad at 100 GeV. The two choices give very different results on this chord:

.. code-block:: python

   kw = dict(
       osc_params=osc, L0=0.0, nu_i=gd.NUE,
       nu_f=gd.NUE, average=True,
       density_is_of_number_of_electrons=True)
   f = oscprob.osc_prob_matter_std_potential
   for start in ('flavor', 'decohered'):
       P = f(3, ne, E, 2*half, **kw,
             average_initial_state=start)
       print(start, P)  # 0.535944, 0.283618


A neutrino from a distant source arrives without that interference, which is why Figure :ref:`The Sun in the electron-neutrino channel <ex-fig-solar-tomography>` uses ``'decohered'``. Magνs computes the phase-averaged probability over the default relative energy spread of 10%, without sampling energies. Averaging a scan of probabilities over energy instead does not, in general, give the same result (:ref:`ex-sec-averaging-is-not-estimating`). On this chord, a little interference survives the 10% spread, so the phase-averaged probability depends slightly on the spread, by about 0.01 per :math:`e`-fold; Magνs issues a ``PhaseAveragingWarning`` saying so.

One phase survives the average: the matter phase, :math:`\int V_{\rm CC}\,dl`, which is :math:`8\,699` rad across the solar diameter. Because it does not depend on energy, no energy spread averages it out. Each chord, however, crosses a different amount of matter, so the matter phase varies with the impact parameter. As :math:`b` changes, the phase cycles through :math:`2\pi` over and over, and :math:`P` oscillates with it, drawing fine rings in the core of the disk. The rings are :math:`1.8\times10^{-4}\,R_\odot` apart at :math:`b = 0.11\,R_\odot` and :math:`5.1\times10^{-4}\,R_\odot` apart at :math:`b = 0.3\,R_\odot`, and at 1 TeV they change :math:`P` by up to about :math:`\pm 0.03`. The probability computed in the snippets above, :math:`P = 0.2836` at :math:`b = 0.3\,R_\odot`, lies on one of these rings: as :math:`b` moves by one ring spacing, :math:`5.1\times10^{-4}\,R_\odot`, around that value, :math:`P` changes appreciably, from 0.282 to 0.308. These rings are far finer than any realistic experimental angular resolution, so an observation would average over many of them. Therefore, each pixel of Figure :ref:`The Sun in the electron-neutrino channel <ex-fig-solar-tomography>` is averaged over its area, which also keeps the rings from aliasing into false patterns. Each panel from 30 GeV to 3 TeV uses :math:`7\,791` impact parameters. Most of them, :math:`6\,790`, lie inside :math:`b = 0.5\,R_\odot`, where the rings are densest.

At :math:`10` MeV, the disk is uniform at :math:`\sum_i |\mathbb{U}_{ei}|^4 \approx 0.55`: the Sun is transparent, as :math:`P^{\rm cross} = \mathbb{1}` requires. At intermediate energies, where Figure :ref:`Adiabaticity along a solar chord <ex-fig-solar-adiabaticity>` shows the passage to be non-adiabatic, rings appear: there, the level-crossing probabilities depart from the identity. The phase average also keeps part of the interference between eigenstates that :ref:`averaged limit on a varying profile <avg-varying>` drops, and this changes the depth of the rings: at 300 GeV, the deepest ring reaches :math:`P = 0.03` at :math:`b = 0.81\,R_\odot`, against :math:`0.22` from :ref:`averaged limit on a varying profile <avg-varying>`. At higher energies, the deepest ring moves outward, to :math:`b = 0.88\,R_\odot` at 1 TeV and :math:`0.91\,R_\odot` at 10 and 50 TeV. At 50 TeV, the disk is uniform again, but for a different reason than at 10 MeV: inside the Sun, the matter potential dominates, so a :math:`\nu_e` that enters the Sun leaves it as a :math:`\nu_e`, at every impact parameter, having only gained a phase, which does not affect the survival probability.

No neutrino telescope can resolve this structure. Super-Kamiokande has imaged the Sun in MeV neutrinos :cite:p:`Super-Kamiokande:2005wtt`, but the image is blurred to about :math:`25^\circ` at 10 MeV, because the recoil electron in neutrino-electron scattering follows the neutrino direction only loosely; the solar disk spans half a degree. At TeV energies, neutrino telescopes reconstruct muon tracks to within about :math:`1^\circ`, comparable to the size of the disk, and particle showers only to within :math:`10^\circ`–:math:`15^\circ` :cite:p:`IceCube:2019cia`. The rings of Figure :ref:`The Sun in the electron-neutrino channel <ex-fig-solar-tomography>`, by contrast, are a few to a few tens of arcseconds wide. Moreover, at GeV energies, where the rings first appear, no astrophysical neutrino flux has been identified above the atmospheric background. Finally, both Figure :ref:`Adiabaticity along a solar chord <ex-fig-solar-adiabaticity>` and Figure :ref:`The Sun in the electron-neutrino channel <ex-fig-solar-tomography>` include only oscillations: absorption in the Sun, which becomes important above about 100 GeV, is outside what Magνs models (:ref:`limitations <what-magnus-is-not>`).

.. _ex-sec-jet:

Astrophysical relativistic jet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gamma-ray bursts are the relativistic jets that some massive stars launch as they collapse :cite:p:`Woosley:1993wj,MacFadyen:1998vz,Piran:2004ba,Meszaros:2006rc`.  Shells of plasma ejected by the central engine collide with one another along the jet :cite:p:`Rees:1994nw`; protons accelerated at those collisions can make TeV–PeV neutrinos :cite:p:`Paczynski:1994uv,Waxman:1997ti,Guetta:2003wi,Murase:2005hy,Hummer:2011ms,Bustamante:2014oka`. The jet itself does not change their flavor composition: at the radii where the collisions happen, :math:`10^{13}`–:math:`10^{15}` cm, its density is around :math:`10^{-11}` g cm\ :math:`^{-3}`, nine orders of magnitude below the density at which a TeV neutrino resonates; the phase that the jet matter adds to the evolution is below :math:`10^{-4}`. Matter effects arise only while the jet is still inside the star, in the choked jets that never break out and in the precursor phase of the jets that do  :cite:p:`Mena:2006eq,Razzaque:2009kq,Sahu:2010ap,Varela:2014mma,Xiao:2015gea`. The neutrinos then cross the stellar envelope. Its density falls from about 0.3 g cm\ :math:`^{-3}` at the head of the jet to zero at the surface, so every neutrino above about 50 GeV meets its 1-3 resonance somewhere along the way. At low energies, the crossing is adiabatic, and the neutrino leaves the star in a mass eigenstate; at high energies, it is not, and the neutrino keeps its flavor content. Therefore, the flavor ratios at Earth depend on the neutrino energy, in a way that reflects the density profile of the star  :cite:p:`Razzaque:2009kq,Xiao:2015gea`.

As example, we take a jet whose head is at :math:`r_0 = 6.3 \cdot 10^{10}` cm, inside a blue supergiant of radius :math:`R_\star = 3 \cdot 10^{12}` cm, with the hydrogen envelope of Refs.  :cite:p:`Mena:2006eq,Razzaque:2009kq`, :math:`\rho(r) = 3.3 \cdot 10^{-6}\,(R_\star/r - 1)^3` g cm\ :math:`^{-3}`, and an electron fraction of :math:`Y_e = 1`:

.. code-block:: python

   import numpy as np
   import magnus.globaldefs as gd
   import magnus.hamiltonians as ham
   import magnus.oscprob as oscprob

   OSC = gd.load_nufit_params('NuFIT 6.1')
   CM = 1.0e-5*gd.UNIT_KM   # One cm in eV^-1
   RSTAR = 3.0e12           # Star's radius, cm
   R0 = 6.3e10              # Jet head, cm

   def rho(l):
       """Density, g/cm3, at l in eV^-1."""
       return 3.3e-6*(RSTAR/(l/CM) - 1.0)**3


The density at the head is :math:`0.33` g cm\ :math:`^{-3}`. The oscillation length there is set by the potential, :math:`2\pi/V_{\rm CC} \approx 5 \cdot 10^{9} {\rm cm} \approx 0.1 r_0`.

A neutrino born as :math:`\nu_\alpha` at the jet head leaves the star in the state :math:`\mathbb{U}(R_\star, r_0)\,\lvert\nu_\alpha\rangle`, where :math:`\mathbb{U}(R_\star, r_0)` is the evolution operator across the envelope. Its weight in the mass eigenstate :math:`\nu_i` is

.. math::

   w_{\alpha i} = \bigl\lvert [\mathbb{V}^\dagger(R_\star)\,\mathbb{U}(R_\star, r_0)]_{i\alpha} \bigr\rvert^2 \;,

where :math:`\mathbb{V}(R_\star)` is the vacuum mixing matrix. On its way to Earth, over a cosmological distance, the mass eigenstates lose coherence  :cite:p:`Razzaque:2009kq`. The probability at Earth is then

.. math::

   \langle P_{\nu_\alpha \to \nu_\beta} \rangle = \sum_i w_{\alpha i}\, \lvert \mathbb{V}_{\beta i}(R_\star) \rvert^2 \;.

This is :ref:`averaged limit on a varying profile <avg-varying>`, with the contribution of the star computed directly as :math:`w_{\alpha i}` and the vacuum eigenbasis at detection.

The weights :math:`w_{\alpha i}` depend on the relative phases of the flavor amplitudes at the surface of the star. The flavor probabilities, :math:`\lvert \mathbb{U}_{\beta\alpha}(R_\star, r_0) \rvert^2`, discard those phases, so the weights must be computed from the evolution operator. The scenario functions return the operator alongside the probabilities when called with ``return_evolution_operator=True``. As usual, the refinement ladder chooses the slabs:

.. code-block:: python

   E = 1.0*gd.UNIT_TEV  # 1 TeV
   # Evolve from the jet head (L0) to the
   # surface (L); also return the operator
   out = oscprob.osc_prob_matter_std_potential(
       3, rho, E, RSTAR*CM, OSC, L0=R0*CM,
       electron_fraction=1.0,  # hydrogen
       rtol=1.0e-6, atol=1.0e-6,
       density_matter_is_in_g_per_cm3=True,
       return_evolution_operator=True)
   # Probabilities and evolution operator
   P, U = out


Here, ``P`` holds the flavor probabilities at the surface of the star, and ``U`` is :math:`\mathbb{U}(R_\star, r_0)` in the flavor basis.

.. _ex-fig-jet:

.. figure:: ../../../img/paper/jet.svg
   :width: 95%
   :alt: Jet neutrinos in a stellar envelope

   **Jet neutrinos in a stellar envelope.** Neutrinos from a relativistic jet inside a collapsing star. *Top:* the setup. Neutrinos are made at the head of the jet, at :math:`r_0 = 6.3 \cdot 10^{10}` cm from the center, and leave the star at :math:`R_\star = 3 \cdot 10^{12}` cm. *Middle:* matter density along their path, for a smooth envelope, the same envelope with turbulence overlaid, and the same envelope with a drop in density at the edge of the helium core. The right axis gives the energy whose 1-3 resonance lies at each density. *Bottom:* probability that a :math:`\nu_e` produced at the jet head is detected as :math:`\nu_e` at Earth, against its energy, for the three envelopes. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

From ``U``, the weights and the probability at Earth follow:

.. code-block:: python

   # Vacuum mixing matrix, V(R_star)
   hv = np.array(ham.hamiltonian_3nu_vacuum(
       E, **OSC), dtype=complex)
   _, R = np.linalg.eigh(hv)
   # Weights, w[i, alpha]
   w = abs(R.conj().T @ U)**2
   # Probability at Earth, P[beta, alpha]
   P_earth = abs(R)**2 @ w
   P_earth[gd.NUE, gd.NUE]    # 0.315


Without matter effects, the same calculation gives :math:`0.55`.

Figure :ref:`Jet neutrinos in a stellar envelope <ex-fig-jet>` shows this probability from 0.1 TeV to 10 PeV for three density envelopes. The first envelope is the smooth profile above. The second adds turbulence to it, as one realization of forty Kolmogorov modes with wavelengths from :math:`10^{10}` to :math:`10^{12}` cm and a root-mean-square amplitude of 20%. The third is model C of Refs.  :cite:p:`Mena:2006eq`, in which the density drops by a factor of five at the edge of the helium core, at :math:`10^{11}` cm. The code receives this edge as a breakpoint. The neutrinos are made across a region about one oscillation length wide, so each curve averages over eight production points spread across that length. At the lowest energies, the resonance lies only a few oscillation lengths from the head. A single production point would then show an interference pattern, which the average removes.

Listing :ref:`Jet neutrinos in a stellar envelope <ex-lst-jet>` contains the calculation behind the figure. Each envelope takes well under a minute, over 161 energies and eight production points. At the requested tolerance of :math:`10^{-6}`, the ladder uses a few hundred slabs per energy at PeV energies and a few thousand at 0.1 TeV.

.. _ex-lst-jet:

**Jet neutrinos in a stellar envelope.** Computing the data in Figure :ref:`Jet neutrinos in a stellar envelope <ex-fig-jet>`: the averaged :math:`\nu_e` survival probability at Earth for neutrinos from a jet head at :math:`r_0`, through the three envelopes. Each call returns the evolution operator from the refinement ladder alongside the probabilities; the projection onto the vacuum mass states and the sum over them are :ref:`averaged limit on a varying profile <avg-varying>`, with the eigenbasis at detection the vacuum one. The eight production points span one oscillation length at the jet head, :math:`4.9 \cdot 10^{9}` cm. The turbulent envelope is one realization of forty Kolmogorov modes between :math:`10^{10}` and :math:`10^{12}` cm at 20% rms; the stepped one is model C of Refs.  :cite:p:`Mena:2006eq`, with the drop at the helium-core edge declared through ``t_breakpoints``. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

.. code-block:: python

   import numpy as np
   import magnus.globaldefs as gd
   import magnus.hamiltonians as ham
   import magnus.oscprob as oscprob

   OSC = gd.load_nufit_params('NuFIT 6.1')
   CM = 1.0e-5*gd.UNIT_KM               # One cm in eV^-1
   RSTAR, R0, RHE = 3.0e12, 6.3e10, 1.0e11      # cm
   E = np.geomspace(0.1, 1.0e4, 161)*gd.UNIT_TEV

   # The three envelopes, g/cm3, r in cm
   def smooth(r):
       return 3.3e-6*(RSTAR/r - 1.0)**3

   rng = np.random.default_rng(7)        
   lam = np.geomspace(1.0e12, 1.0e10, 40)  # Modes, cm
   k = 2*np.pi/lam
   amp = k**(-1/3)                         # Kolmogorov
   amp *= 0.20/np.sqrt(0.5*np.sum(amp**2)) # 20% rms
   phi = rng.uniform(0.0, 2*np.pi, 40)
   def turbulent(r):
       d = np.sum(amp*np.cos(k*r + phi))
       return smooth(r)*(1.0 + d)

   def stepped(r):            # Mena et al., model C
       x = RSTAR/r - 1.0
       return 6.3e-6*(20.0*x**2.1 if r < RHE else x**2.5)

   # Production points spread over one oscillation
   # length at the jet head
   r0s = R0 + 4.9e9*np.arange(8)/8.0

   def p_earth(rho, **extra):
       """P at Earth, phase-averaged, and averaged
       over the production points."""
       P_out = np.zeros((len(E), 3, 3))
       for r0 in r0s:
           P, U = oscprob.osc_prob_matter_std_potential(
               3, lambda l: rho(l/CM), E, RSTAR*CM, OSC,
               L0=r0*CM, electron_fraction=1.0,
               density_matter_is_in_g_per_cm3=True,
               rtol=1.0e-6, atol=1.0e-6,
               return_evolution_operator=True, **extra)
           for i, e in enumerate(E):
               hv = np.array(ham.hamiltonian_3nu_vacuum(
                   e, **OSC), dtype=complex)
               w, R = np.linalg.eigh(hv)
               content = abs(R.conj().T @ U[i])**2
               P_out[i] += abs(R)**2 @ content/len(r0s)
       return P_out

   P_smooth = p_earth(smooth)
   P_turb = p_earth(turbulent)
   P_step = p_earth(stepped, t_breakpoints=[RHE*CM])
   # Drawn: P_smooth[:, gd.NUE, gd.NUE], and the others


For the smooth envelope, the :math:`\nu_e` survival probability at Earth rises with energy. At 0.1 TeV, a :math:`\nu_e` leaves the star mostly as :math:`\nu_3`, the mass eigenstate with the smallest electron content, so it is rarely detected as :math:`\nu_e`. At 100 TeV, the probability reaches its vacuum value. The probability rises because the resonance crossing becomes less adiabatic as the energy grows. The resonance density falls as :math:`1/E`. Deep in the envelope, :math:`\rho \propto r^{-3}`, so the resonance moves outward as :math:`r_{\rm res} \propto E^{1/3}`. There, the density changes over a distance :math:`\rho/\lvert d\rho/dr \rvert = r_{\rm res}/3 \propto E^{1/3}`, while the oscillation length at the resonance grows as :math:`E`. Their ratio, which sets how adiabatic the crossing is, therefore falls as :math:`E^{-2/3}`. Thus, at high energies, the crossing is non-adiabatic, and the neutrino keeps its flavor. At the PeV energies that IceCube detects, the star leaves the flavor content unchanged.

Turbulence shifts the probability by up to 0.05 below a few TeV. At these energies, the turbulent modes near the resonance have wavelengths comparable to the oscillation length (as in :ref:`ex-sec-turbulence`). At the drop, the density falls by a factor of five within a small fraction of an oscillation length. Neutrinos between 0.1 and 0.55 TeV have their 1-3 resonance density inside the drop. They skip over their resonance, rather than crossing it gradually, so the crossing is non-adiabatic. This produces the plateau at about 0.36 in Figure :ref:`Jet neutrinos in a stellar envelope <ex-fig-jet>`, which extends to about 1 TeV. Above about 1.5 TeV, the resonance lies well outside the helium core, and the drop no longer matters: the stepped and smooth envelopes give the same probability to within 0.01.

.. _ex-sec-lri-sun:

Long-range interactions in the Sun
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _ex-fig-solar-lri:

.. figure:: ../../../img/paper/solar_long_range.svg
   :width: 95%
   :alt: A long-range force in the Sun

   **A long-range force in the Sun.** Averaged survival probability of :math:`\nu_e` solar neutrinos under a long-range :math:`L_e - L_\mu` interaction sourced by the electrons of the Sun. *Top:* the electrons that the neutrino feels at each mediator range. A :math:`\nu_e` born at :math:`0.05\,R_\odot` leaves along the ray drawn; the discs enclose the electrons within a distance of :math:`1/m` of it along the way. The charged-current potential is local, so it has no disc. At :math:`1/m = R_\odot/10`, the neutrino feels a narrow tube of electrons. At :math:`1/m = R_\odot`, it feels most of the Sun at once, even beyond the surface. *Middle:* :math:`\langle P_{\nu_e \to \nu_e}\rangle` against energy, in the B16-GS98 solar model :cite:p:`Vinyoles:2016djt`, without the new potential and with it at the two ranges. Both ranges share one coupling, fixed so that the longer-ranged potential is a tenth of :math:`V_{\rm CC}` where the neutrino is born. The probability is read out at :math:`20\,R_\odot`, where the new potential has faded. *Bottom:* the difference from the standard case. See notebooks `#19 <https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb>`__ and `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.


This example adds a new term to the Hamiltonian, one that Magνs does not ship. We write it as described in :ref:`ex-sec-building-new-hamiltonian`, and average the probability as in the standard solar case of :ref:`ex-sec-sun`. The new term comes from gauging :math:`L_e - L_\mu`, the difference between the electron and muon lepton numbers, which is anomaly-free. Its gauge boson, :math:`Z'`, is a new light mediator :cite:p:`He:1990pn,Foot:1994vd`. The charge :math:`L_e - L_\mu` is :math:`+1` for the electron and :math:`\nu_e`, :math:`-1` for the muon and :math:`\nu_\mu`, and zero for every other field of the Standard Model. Therefore, the electrons of the Sun source a potential that :math:`\nu_e` and :math:`\nu_\mu` feel with opposite signs, and that :math:`\nu_\tau` does not feel. The Hamiltonian is

.. math::
   :label: ex-equ-h3-lri

   \mathbb{H}^{\rm LRI}(E, l)
   =
   \frac{\mathbb{A}}{E} + V_{\rm CC}(l)\,\mathbb{P}
   + V_{e\mu}(l)\,{\rm diag}\left(1, -1, 0\right) \;,

where the first two terms are the standard vacuum and charged-current terms of :doc:`/conventions`, and the third is the new one.

Unlike :math:`V_{\rm CC}`, the new potential is non-local. Each electron sources a Yukawa potential around it, with a range set by the mass :math:`m` of the :math:`Z'`. A neutrino at :math:`\mathbf{r}` feels the electrons within about :math:`1/m` of it:

.. math::
   :label: ex-equ-yukawa

   V_{e\mu}(\mathbf{r})
   =
   \frac{g^{\prime 2}}{4\pi}
   \int d^3r^\prime \, n_e(\mathbf{r}^\prime) \,
   \frac{e^{-m \lvert \mathbf{r} - \mathbf{r}^\prime \rvert}}
   {\lvert \mathbf{r} - \mathbf{r}^\prime \rvert} \;,

where :math:`g^\prime` is the new gauge coupling.

When the interaction range, :math:`1/m`, is larger than the body, the potential depends mainly on how many electrons the body holds, and only weakly on where they are :cite:p:`Wise:2018rnb`. When :math:`1/m` is much shorter than the distance over which the density changes, only nearby electrons contribute, and :math:`V_{e\mu} \to g^{\prime 2} n_e(\mathbf{r})/m^2`. In this limit, the potential is local and proportional to :math:`V_{\rm CC}`, like a non-standard interaction with constant couplings :cite:p:`Coloma:2020gfv`. Between the two limits, the potential depends on the shape of the density profile. For a spherical body, like the Sun, the angular integral can be done analytically. The potential then splits into the contributions of the electrons inside and outside radius :math:`r`:

.. math::
   :label: ex-equ-yukawa-1d

   \begin{split}
   V_{e\mu}(r)
   = {} &
   \frac{g^{\prime 2}}{r}\, e^{-m r}
   \int_0^{r} dr^\prime\, r^{\prime 2}\, n_e(r^\prime)\, {\rm shc}(m r^\prime) \\
   & + g^{\prime 2}\, {\rm shc}(m r)
   \int_r^{R_\odot} dr^\prime\, r^\prime\, n_e(r^\prime)\, e^{-m r^\prime} \;,
   \end{split}

where :math:`{\rm shc}(x) \equiv \sinh(x)/x`. Both are running integrals over the density profile, so one pass over the table of a solar model gives the potential at every radius. For a massless mediator, :math:`{\rm shc}` and the exponentials tend to 1, and Eq. :eq:`ex-equ-yukawa-1d` becomes the Coulomb potential of the electrons. Notebook `#19 <https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb>`__ derives Eq. :eq:`ex-equ-yukawa-1d` and checks it against the closed form for a uniform ball.

Figure :ref:`A long-range force in the Sun <ex-fig-solar-lri>` shows the averaged :math:`\nu_e` survival probability from 0.1 to 20 MeV. The density profile is from the B16-GS98 solar model :cite:p:`Vinyoles:2016djt`. The neutrino is born at :math:`0.05\,R_\odot`, near where :math:`^8`\ B production peaks in that model. The figure compares the standard case with two mediator ranges, :math:`1/m = R_\odot` and :math:`R_\odot/10`, corresponding to :math:`m \approx 3 \cdot 10^{-16}` and :math:`3 \cdot 10^{-15}` eV, respectively. Both ranges share one coupling. We fix it so that the longer-ranged potential is a tenth of :math:`V_{\rm CC}` where the neutrino is born. The longer-ranged potential lowers the probability by up to 0.018, near 4 MeV. The shorter-ranged one lowers it by at most 0.005. Listing :ref:`A long-range force in the Sun <ex-lst-arbitrary>` contains the calculation behind the figure.

At the same coupling, a longer range lets the neutrino feel more of the Sun's electrons. Where the neutrino is born, the longer-ranged potential is :math:`0.10\,V_{\rm CC}`, about four times the shorter-ranged one, :math:`0.027\,V_{\rm CC}`. This difference at birth sets the difference between the two curves. The neutrino leaves the Sun adiabatically, so its averaged probability depends only on the Hamiltonian where it is born and where it is detected.

Along the path, both new potentials fall more slowly than :math:`V_{\rm CC}`. The electron density of the Sun falls by about five orders of magnitude from the center to :math:`0.98\,R_\odot`, and :math:`V_{\rm CC}` falls with it. In contrast, :math:`V_{e\mu}` is an integral over the electrons within about :math:`1/m`, which for the longer range include those of the dense core. Using the variables defined in Listing :ref:`A long-range force in the Sun <ex-lst-arbitrary>`, the ratio :math:`V_{e\mu}/V_{\rm CC}` at half a solar radius is

.. code-block:: python

   half = np.searchsorted(r, 0.5*R_SUN)
   v_emu[1.0][half]/vcc[half]     # 2.39
   v_emu[0.1][half]/vcc[half]     # 0.098

At :math:`0.98\,R_\odot`, the two ratios are about 700 and 1.2. The top panel of Figure :ref:`A long-range force in the Sun <ex-fig-solar-lri>` illustrates this. With :math:`1/m = R_\odot`, the neutrino feels most of the Sun wherever it is, and the potential extends beyond the surface. We read out the probability at :math:`20\,R_\odot` away from the solar center, where the new potential has faded.

In Listing :ref:`A long-range force in the Sun <ex-lst-arbitrary>`, the Hamiltonian accepts an array of positions, so the engine can evaluate it at many positions in one call (:doc:`/performance`). The averaged probability comes from ``osc_prob_energy_baseline`` with ``average=True``. The standard case uses the same call, with :math:`V_{e\mu} = 0`. At 10 MeV,

.. code-block:: python

   i = np.argmin(abs(Es - 10.0*gd.UNIT_MEV))
   P_std[i]                      # 0.325
   P_lri[1.0][i], P_lri[0.1][i]  # 0.317, 0.323

.. _ex-fig-solar-8b-flux:

.. figure:: ../../../img/paper/solar_8b_flux.svg
   :width: 95%
   :alt: The 8B neutrino flux at Earth under a long-range force

   **The** :math:`^8`\ **B neutrino flux at Earth under a long-range force.** *Top:* flux of :math:`\nu_e` from :math:`^8`\ B decay at Earth, against energy, for the standard case and for the two mediator ranges of Figure :ref:`A long-range force in the Sun <ex-fig-solar-lri>`, with the same coupling. The dotted curve has :math:`1/m = R_\odot` and a coupling :math:`g^{\prime 2}` three times larger. The fluxes are the :math:`^8`\ B spectrum of Ref. :cite:p:`Winter:2004kf` times the total flux of the B16-GS98 model, :math:`5.46 \cdot 10^6` cm\ :sup:`-2` s\ :sup:`-1` :cite:p:`Vinyoles:2016djt`, times the survival probability averaged over where :math:`^8`\ B neutrinos are made in the model. *Bottom:* each case divided by the standard one. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

Figure :ref:`The 8B neutrino flux at Earth under a long-range force <ex-fig-solar-8b-flux>` shows the flux of :math:`\nu_e` from :math:`^8`\ B decay at Earth, from 1 to 15 MeV. Since :math:`^8`\ B neutrinos are produced over a range of radii, 80% of them between :math:`0.025\,R_\odot` and :math:`0.081\,R_\odot` in the B16-GS98 model, we average the survival probability over the production distribution:

.. math::
   :label: ex-equ-b8-flux

   \frac{d\Phi_{\nu_e}}{dE}
   =
   \Phi_{^8{\rm B}}\, \lambda(E)
   \int_0^{R_\odot} dr\, f(r)\,
   \langle P_{\nu_e \to \nu_e}(E; r) \rangle \;,

where :math:`\Phi_{^8{\rm B}} = 5.46 \cdot 10^6` cm\ :sup:`-2` s\ :sup:`-1` is the total flux of the model :cite:p:`Vinyoles:2016djt`, :math:`\lambda(E)` is the :math:`^8`\ B spectrum normalized to unity :cite:p:`Winter:2004kf`, :math:`f(r)` is the radial distribution of production, normalized to unity, and :math:`\langle P_{\nu_e \to \nu_e}(E; r) \rangle` is the averaged survival probability of a neutrino born at radius :math:`r`. We evaluate the integral on 28 radii between :math:`0.005\,R_\odot` and :math:`0.14\,R_\odot`, which contain 99.9% of production.

(We compute :math:`f(r)` from the temperature, density, and composition of B16-GS98, as the rate of :math:`^7`\ Be(p, γ) :math:`^8`\ B. In this rate, :math:`^7`\ Be is in equilibrium between its production by :math:`^3`\ He(α, γ) :math:`^7`\ Be and its destruction by electron and proton capture, with the reaction rates of Ref. :cite:p:`Adelberger:2010qa`. To test this calculation, we apply it to an older solar model, BS2005-AGS,OP, whose :math:`^8`\ B production distribution is published :cite:p:`Bahcall:2004pz`. The two distributions agree within 1% of the peak. Notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__ contains the calculations.)

Extending Listing :ref:`A long-range force in the Sun <ex-lst-arbitrary>`, Eq. :eq:`ex-equ-b8-flux` for :math:`1/m = R_\odot` is

.. code-block:: python

   # Birth radii, weighted by where 8B is made
   r_b = np.linspace(0.005, 0.14, 28)*R_SUN
   w = f(r_b)*np.gradient(r_b)
   w = w/w.sum()
   # Average over them, then scale by the 8B
   # spectrum and total flux, cm^-2 s^-1 MeV^-1
   P = sum(wk*averaged(v_emu[1.0], rk)
           for rk, wk in zip(r_b, w))
   flux = 5.46e6*spectrum(Es)*P

where ``f`` and ``spectrum`` interpolate :math:`f(r)` and :math:`\lambda(E)`. The longer-ranged potential lowers the flux by up to 4.5%, near 5 MeV, and by 3.6% in total above 1 MeV. The shorter-ranged one lowers it by at most 1.3%. Tripling :math:`g^{\prime 2}` for :math:`1/m = R_\odot` lowers it by up to 11%, near 4.5 MeV.


.. _ex-lst-arbitrary:

**A long-range force in the Sun.** Computing the data in Figure :ref:`A long-range force in the Sun <ex-fig-solar-lri>`. The electron density comes from the B16-GS98 table :cite:p:`Vinyoles:2016djt`, continued past the surface. The potential :math:`V_{e\mu}` is Eq. :eq:`ex-equ-yukawa-1d` on that profile. The vacuum and matter terms come from the shipped builders; only the new term is written here. The neutrino is born at :math:`0.05\,R_\odot`, and the :ref:`averaged probability <avg-varying>` is read out at :math:`20\,R_\odot`. See notebooks `#19 <https://github.com/mbustama/Magnus/blob/main/notebooks/19_magnus_custom_hamiltonian.ipynb>`__ and `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.

.. code-block:: python

   import numpy as np
   import magnus.globaldefs as gd
   import magnus.hamiltonians as ham
   import magnus.matter as matter
   import magnus.oscprob as oscprob
   import magnus.solarmodels as solarmodels

   osc = gd.load_nufit_params('NuFIT 6.1')

   # The B16-GS98 model that ships with Magnus,
   # on its tabulated radii and on to 20 R_sun
   model = 'B16-GS98'
   tab = solarmodels.load_solar_model(model)
   r_tab = tab['r_over_r_sun']*gd.SUN_RADIUS*gd.UNIT_KM
   R_SUN = r_tab[-1]                   # The surface
   R0, R_END = 0.05*R_SUN, 20.0*R_SUN  # Birth, readout
   r = np.concatenate([r_tab, np.geomspace(
       1.002*R_SUN, R_END, 800)])
   ne_func = solarmodels.electron_density_profile(model)
   n_e = ne_func(r)                    # Fades past R_SUN
   per_ne = matter.VCC_func(0.0, lambda l: 1.0)
   vcc = per_ne*n_e

   def vcc_sun(l):
       """V_CC at l."""
       return per_ne*ne_func(l)

   def cumulative(y, x):
       """Running trapezoidal integral, zero at x[0]."""
       steps = 0.5*(y[1:] + y[:-1])*np.diff(x)
       return np.concatenate([[0.0], np.cumsum(steps)])

   def shc(x):
       """sinh(x)/x, equal to 1 at the origin."""
       x = np.asarray(x, dtype=float)
       return np.divide(np.sinh(x), x, out=np.ones_like(x),
                        where=(x != 0.0))

   def long_range_potential(r, n_e, m):
       """V_emu(r)/g'^2 for a spherical n_e."""
       I_in = cumulative(r**2*n_e*shc(m*r), r)
       I_out = cumulative(r*n_e*np.exp(-m*r), r)
       I_out = I_out[-1] - I_out      # From r outward
       r_safe = np.where(r > 0.0, r, 1.0e-30)
       return np.exp(-m*r)*I_in/r_safe + shc(m*r)*I_out

   # Two mediator ranges, one coupling: g'^2 makes V_emu
   # a tenth of V_CC at R0 when 1/m = R_sun
   v_emu = {}
   for frac in (1.0, 0.1):            # 1/m in R_sun
       v_emu[frac] = long_range_potential(
           r, n_e, 1.0/(frac*R_SUN))
   g2 = 0.1*vcc_sun(R0)/np.interp(R0, r, v_emu[1.0])
   v_emu = {frac: g2*v for frac, v in v_emu.items()}

   # The standard part from the shipped builders; the
   # new term is the only line written here
   h_vac = ham.\
       hamiltonian_3nu_vacuum_energy_independent(**osc)
   q = np.diag([1.0, -1.0, 0.0])     # L_e-L_mu charges

   def H_lri(v):
       """H(E, l) of the text, with V_emu = v."""
       def H(E, l):
           h_matt = ham.hamiltonian_3nu_matter_td(
               l, vcc_sun)
           return (h_vac/E + h_matt
                   + np.interp(l, r, v)[..., None, None]*q)
       return H

   Es = np.logspace(-1.0, np.log10(20.0), 70)*gd.UNIT_MEV

   def averaged(v, r0=R0):
       """<P_ee> from r0 to R_END, per E."""
       return oscprob.osc_prob_energy_baseline(
           H_lri(v), Es, R_END, r0, nu_i=gd.NUE,
           nu_f=gd.NUE, average=True)

   P_std = averaged(0.0*vcc)
   P_lri = {frac: averaged(v_emu[frac]) for frac in v_emu}

.. _ex-sec-turbulence:

Turbulent matter profile
~~~~~~~~~~~~~~~~~~~~~~~~

.. _ex-fig-turbulence-rabi:

.. figure:: ../../../img/paper/turbulence_rabi.svg
   :width: 95%
   :alt: A density mode against its closed form

   **A density mode against its closed form.** Probability of transition between two matter levels coupled by a single density mode, against the wavenumber of the Fourier mode of the matter density in the medium, for a 5-GeV neutrino in a region of mean density :math:`\rho_0 = 4` g cm\ :math:`^{-3}` carrying one mode of amplitude :math:`C = 0.03`. The panels are region lengths of :math:`1`, :math:`3`, :math:`10`, and :math:`30\,L_{\rm osc}`, with :math:`L_{\rm osc} = 2\pi/\Delta_{32} = 11\,030` km the oscillation length driven by the 2-3 sector. The closed form of Refs.  :cite:p:`Patton:2013dba,Patton:2014lza`, Eq. :eq:`ex-equ-turb-rabi`, is evaluated for the two-flavor reduction of the 1-3 sector; Magνs is run at two, three, and four flavors, the last with :math:`\sin^2\theta_{14} = \sin^2\theta_{24} = 0.10` and :math:`\Delta m^2_{41} = 1` eV\ :math:`^2`. The levels are those of the mean density, between which no transition occurs without the mode. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.


.. _ex-fig-turbulence:

.. figure:: ../../../img/paper/turbulence.svg
   :width: 95%
   :alt: The same mode in the flavor channel

   **The same mode in the flavor channel.** Flavor transition probability :math:`P_{\nu_\mu \to \nu_e}` for the same neutrino, medium, mode, and region lengths as Figure :ref:`A density mode against its closed form <ex-fig-turbulence-rabi>`, at two, three, and four flavors. The line labeled :math:`q = \Delta_{32}` marks the mode that matches the three-flavor gap; the two- and four-flavor gaps lie a few percent away from it. The curves sit at different heights because the three flavor counts give different probabilities over the same region even without the mode. The closed-form curve is Eq. :eq:`ex-equ-turb-rw-flavor`, the flavor probability built from the propagator behind Eq. :eq:`ex-equ-turb-rabi`, for the same two-flavor system as in Figure :ref:`A density mode against its closed form <ex-fig-turbulence-rabi>`. See notebook `#28 <https://github.com/mbustama/Magnus/blob/main/notebooks/28_magnus_paper_figures.ipynb>`__.


Matter behind a supernova shock does not settle: convection keeps stirring it. Its density fluctuates around its mean value by tens of percent, over a range of length scales spanning from the size of the star to far below it  :cite:p:`Loreti:1995ae,Schirato:2002tg,Fogli:2006xy,Patton:2013dba,Patton:2014lza`. Therefore, an escaping neutrino meets a profile with a density that has a regular component—decreasing with distance from the center—and a random component overlaid on it. In the Fourier decomposition of that random component, a mode of wavenumber :math:`q` moves neutrinos between a pair of eigenstates of :math:`\mathbb{H}` when :math:`q` coincides with the eigenvalue gap :math:`\Delta_{jk} = \lambda_j - \lambda_k`, the random density mode playing the role that an oscillating field plays in parametric resonance (:doc:`/methodology`).

References  :cite:p:`Patton:2013dba,Patton:2014lza` solved that problem in closed form. Treating the random component as a perturbation on the regular one and keeping the resonant term, they obtained a Rabi formula for the transition between the two levels,

.. math::
   :label: ex-equ-turb-rabi

   P_{jk}
   =
   \frac{\kappa^2}{p^2 + \kappa^2}\,
   \sin^2\!\left(\sqrt{p^2 + \kappa^2}\; L\right) ,
   \qquad
   p = \tfrac{1}{2}\left(\Delta_{jk} - q\right) ,

where :math:`L` is the length of the turbulent region and :math:`\kappa \equiv \tfrac{1}{2} C V_{\rm CC} \lvert U^m_{ej} U^m_{ek}\rvert` is the coupling induced between the two levels by a density mode :math:`\rho(l) = \rho_0 \left[1 + C \cos(q l)\right]`, with :math:`\rho_0` the mean density and :math:`U^m` the mixing matrix in matter. On resonance, where :math:`q = \Delta_{jk}`, the probability is :math:`\sin^2(\kappa L)`: the amplitude :math:`C` sets how long the region must be for complete conversion to occur.

As example, we consider a region of mean density :math:`\rho_0 = 4` g cm\ :math:`^{-3}` carrying one mode of amplitude :math:`C = 0.03`, crossed by a 5-GeV neutrino:

.. code-block:: python

   import numpy as np
   import magnus.globaldefs as gd
   import magnus.matter as matter
   import magnus.hamiltonians as ham

   OSC = gd.load_nufit_params('NuFIT 6.1')
   E = 5.0*gd.UNIT_GEV    # Neutrino energy
   RHO0 = 4.0             # Mean density, g/cm3
   C = 0.03               # Mode amplitude

   def one_mode(q):
       """Density along the ray: the mean, plus
       one mode of wavenumber q."""
       def rho(l):
           x = np.asarray(l, dtype=float)
           return RHO0*(1.0 + C*np.cos(q*x))
       return rho


The gaps of the Hamiltonian at :math:`\rho_0` fix the range of :math:`q` to scan:

.. code-block:: python

   # n_e and V_CC of the smooth medium
   ne0 = matter.num_density_e_func(
       0.0, lambda l: RHO0,
       density_matter_is_in_g_per_cm3=True)
   V0 = matter.VCC_func(0.0, lambda l: ne0)

   # H at the mean density: vacuum plus matter
   H0 = np.array(ham.hamiltonian_3nu_vacuum(
       E, **OSC), dtype=complex)
   H0 += ham.hamiltonian_3nu_matter(V0)

   # Eigenvalues w, matter eigenvectors Vm, and
   # the gap the mode has to match
   w, Vm = np.linalg.eigh(H0)
   d32 = abs(w[2] - w[1])
   LOSC = 2*np.pi/d32     # 11 030 km


Here, :math:`\Delta_{32}` is matched by a mode of wavelength :math:`L_{\rm osc} \equiv 2\pi/\Delta_{32} = 11\,030` km, the oscillation length of that pair of levels.

Eq. :eq:`ex-equ-turb-rabi` gives the transition probability between two eigenstates of :math:`\mathbb{H}`\ :cite:p:`Patton:2013dba,Patton:2014lza`. The same quantity is contained in the evolution operator :math:`\mathbb{U}` (:ref:`conv-hamiltonian`), once :math:`\mathbb{U}` is written in the basis of those eigenstates: :math:`\lvert [\mathbb{U}]_{jk} \rvert^2` is the transition probability between levels :math:`k` and :math:`j`. If the density were constant, :math:`\mathbb{H}` would be the same at every point of the region; its eigenstates would then evolve independently of one another, so :math:`\lvert [\mathbb{U}]_{jk} \rvert^2` would vanish for :math:`j \neq k`. The Fourier mode makes the density vary along the region, so :math:`\mathbb{H}` varies with it. In the basis of the eigenstates at :math:`\rho_0`, that variation appears as the off-diagonal elements :math:`C V_{\rm CC} \cos(ql)\, U^{m*}_{ej} U^m_{ek}` of :math:`\mathbb{H}`. Those elements mix the levels as the neutrino advances, so :math:`\mathbb{U}` acquires off-diagonal elements too; :math:`\lvert [\mathbb{U}]_{jk} \rvert^2` is what that mixing amounts to over the whole region. Magνs returns :math:`\mathbb{U}`; rotating it into that basis is one line:

.. code-block:: python

   from magnus.oscprob import (
     compute_evolution_operator_multiple_slabs
     as chain)

   L = LOSC             # The region, one L_osc
   rho = one_mode(d32)  # Mode tuned to the gap

   def H(l):
       """H at position l, with the mode on."""
       h = np.array(ham.hamiltonian_3nu_vacuum(
           E, **OSC), dtype=complex)
       # V_CC follows the density
       return h + ham.hamiltonian_3nu_matter(
           V0*rho(l)/RHO0)

   # One operator per slab, then the ordered
   # product over the chain
   edges = np.linspace(0.0, L, 501)
   slabs = np.stack([edges[:-1], edges[1:]], 1)
   Us = chain(H, slabs, 9, 6)
   U = Us[0]
   for u in Us[1:]:     # Earliest slab first
       U = u @ U

   # Rotate to the matter basis, then read the
   # transition between levels 2 and 3
   P32 = abs((Vm.conj().T @ U @ Vm)[2, 1])**2
   P32                    # 0.001758


Figure :ref:`A density mode against its closed form <ex-fig-turbulence-rabi>` compares Eq. :eq:`ex-equ-turb-rabi` with Magνs at four lengths of the turbulent region; the closed form is evaluated for the two-flavor system it was derived for, the 1–3 sector, with :math:`\kappa` from the two-flavor :math:`\mathbb{U}^m` at :math:`\rho_0`. The closed form of Refs. :cite:p:`Patton:2013dba,Patton:2014lza` keeps only the part of the coupling that varies slowly along the trajectory; the part it discards averages out only over many of its own cycles. One :math:`L_{\rm osc}` (Figure :ref:`A density mode against its closed form <ex-fig-turbulence-rabi>`, top panel) contains too few of them, so the closed form and the two-flavor Magνs curve differ appreciably; the difference shrinks as the region lengthens (Figure :ref:`A density mode against its closed form <ex-fig-turbulence-rabi>`, bottom three panels). A third or fourth flavor changes the eigenvalue gap and moves the resonance to a different wavenumber, while the length of the region sets its width: the longer the region, the more closely the mode has to match the gap.

Figure :ref:`The same mode in the flavor channel <ex-fig-turbulence>` shows the associated oscillation probabilities :math:`P_{\nu_\mu \to \nu_e}`, computed via

.. code-block:: python

   import magnus.oscprob as oscprob

   kw = dict(
       nu_i=gd.NUMU, nu_f=gd.NUE,
       density_matter_is_in_g_per_cm3=True)

   # Two flavors: the 1-3 sector
   OSC2 = dict(sth=OSC['s13'], Dm2=OSC['D31'])
   P2 = oscprob.osc_prob_matter_std_potential(
       2, one_mode(d32), E, L, OSC2, **kw)

   # Three flavors
   P3 = oscprob.osc_prob_matter_std_potential(
       3, one_mode(d32), E, L, OSC, **kw)

   # The same call with a sterile state added
   OSC4 = dict(OSC, s14=np.sqrt(0.10), D41=1.0,
               s24=np.sqrt(0.10), s34=0.0,
               d14=0.0, d24=0.0)
   P4 = oscprob.osc_prob_matter_std_potential(
       4, one_mode(d32), E, L, OSC4, **kw)


In Figure :ref:`The same mode in the flavor channel <ex-fig-turbulence>`, the largest change of the probability sits around the wavenumber :math:`q = \Delta_{32}`. This is the resonance of Figure :ref:`A density mode against its closed form <ex-fig-turbulence-rabi>`, seen in the flavor channel; as there, it narrows as the region lengthens. Above :math:`\Delta_{32}` the curves are flat: a mode faster than every gap cancels its effect along the trajectory.

Figure :ref:`The same mode in the flavor channel <ex-fig-turbulence>` also shows the flavor probability that follows from the closed form. Near the resonance it follows the two-flavor Magνs curve, more closely the longer the region. At low wavenumber it does not: there it stays at the value of the smooth medium, while Magνs swings above it. That curve is built as follows. Eq. :eq:`ex-equ-turb-rabi` is the modulus squared of one element of the propagator of the two-level system in the rotating-wave approximation, the solution behind the closed form  :cite:p:`Patton:2013dba,Patton:2014lza`. In the matter basis of the mean density, that propagator is

.. math::
   :label: ex-equ-turb-rw-u

   \tilde{\mathbb{U}}(L)
   =
   {\rm diag}\!\left(e^{iqL/2}, e^{-iqL/2}\right)
   \exp\!\left[
   -i L
   \begin{pmatrix}
   \bar\lambda - p & \kappa_c \\
   \kappa_c^* & \bar\lambda + p
   \end{pmatrix}
   \right] ,

where :math:`\bar\lambda = \tfrac{1}{2}(\lambda_1 + \lambda_2)` and :math:`p = \tfrac{1}{2}(\lambda_2 - \lambda_1 - q)` are the mean of the two levels and the detuning, while :math:`\kappa_c = \tfrac{1}{2} C V_{\rm CC}\, U^{m*}_{e1} U^m_{e2}` is the coupling, with :math:`\lvert \kappa_c \rvert = \kappa`. The matrix in the exponent is the Hamiltonian in the frame that rotates with the mode, in which it is constant; the first factor undoes that rotation. The flavor probability reads the same propagator with flavor indices,

.. math::
   :label: ex-equ-turb-rw-flavor

   P_{\nu_\mu \to \nu_e}
   =
   \left\lvert \left[ U^m\, \tilde{\mathbb{U}}(L)\, U^{m\dagger} \right]_{e\mu} \right\rvert^2 .

The matrix in Eq. :eq:`ex-equ-turb-rw-u` couples the two levels but does not move them: the rotating-wave approximation drops the part of the mode that raises and lowers the levels along with the density. At low wavenumber that part is all there is. Such a mode changes the density slowly; the levels follow it, with no transition driven between them. Eq. :eq:`ex-equ-turb-rw-u` misses that effect, so its curve stays at the value of the smooth medium at low wavenumber.
