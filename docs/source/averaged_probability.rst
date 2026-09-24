Phase-Averaged Probabilities
============================

This page documents the ``average`` keyword, and the ``average_spread`` keyword that sets its
width, accepted by
:func:`magnus.oscprob.osc_prob_vacuum`,
:func:`magnus.oscprob.osc_prob_matter_std_potential`,
:func:`magnus.oscprob.osc_prob_matter_nsi` and
:func:`magnus.oscprob.osc_prob_liv` -- and therefore, through the shared
``**kwargs`` chain, by every ``osc_prob_{2,3,4,5}nu_*`` wrapper built on
them -- together with the module that implements them,
:mod:`magnus.avgprob`. See :doc:`adiabatic_strategy` for the
position-dependent machinery this reuses, and :doc:`methodology` for the
plain Magnus engine both sit alongside.

The problem: a phase nobody can resolve
------------------------------------------

A neutrino from an astrophysical source arrives with an oscillation phase

.. math::

   \Delta \phi = \frac{\Delta m^2 L}{2E}
   \simeq 1.27 \times \frac{\Delta m^2/\text{eV}^2 \times L/\text{km}}
   {E/\text{GeV}}

of order :math:`10^{15}` for a TeV neutrino from 100 Mpc away.  No
ingredient of that number is known to anything close to the precision the
phase would demand: not the source distance, not the size of the
production region, and not the detector's energy resolution.  Whatever
the true phase is, the measurement integrates over many complete cycles
of it.

Computing such a probability by propagation is therefore doubly
unattractive.  It is expensive -- resolving :math:`10^{15}` radians is
exactly the regime that defeats slab refinement -- and it is pointless,
because every oscillatory term is about to be averaged away by the
integration the measurement performs anyway.

The averaged limit
---------------------

Write the amplitude in the basis that diagonalizes the Hamiltonian,
:math:`H = V \,\mathrm{diag}(\lambda_i)\, V^\dagger`:

.. math::

   A(\nu_\alpha \to \nu_\beta) = \sum_i V^*_{\alpha i} V_{\beta i}\,
   e^{-i \lambda_i L} .

The probability :math:`|A|^2` contains a diagonal part and interference
terms carrying :math:`e^{-i(\lambda_i - \lambda_j)L}`.  Averaging over
the phase leaves only the terms whose phase does not vary:

.. math::

   \boxed{\;P(\nu_\alpha \to \nu_\beta) = \sum_i |V_{\alpha i}|^2\,
   |V_{\beta i}|^2\;}

This is the exact :math:`L/E \to \infty` limit, reached when every relative phase runs
through many cycles across what the measurement cannot resolve, and it costs one matrix
product rather than an integration.  Three properties of the limit are worth stating,
because each surprises someone eventually, and none survives away from it:

* The result is **symmetric** in :math:`\alpha \leftrightarrow \beta`, so
  the averaged probability is the same in both directions.
* It is **identical for neutrinos and antineutrinos**, since
  :math:`|V^*|^2 = |V|^2`.  CP violation does not survive the average,
  even though :math:`\delta_{\rm CP}` still enters through the
  magnitudes :math:`|V_{\alpha i}|`.
* For **vacuum** oscillations it does not depend on energy or baseline at
  all: scaling :math:`H` by :math:`1/E` leaves its eigenvectors
  untouched, so a single matrix serves an entire flux calculation.

The phase average
-------------------

A measurement averages over what it cannot resolve, and the limit assumes that this includes
every relative phase.  It need not.  A detector with an energy resolution of 10% averages a
phase over the range the phase covers across that resolution, and a phase that barely changes
with energy is not averaged at all: over 1000 km at 1 GeV the atmospheric phase is about
6 rad, and the limit misses the oscillation entirely.

``average=True`` therefore returns the **phase average**.  Every interference term keeps its
phase at the central energy and is weighted by the spread of that phase across a relative
energy spread :math:`\sigma`:

.. math::

   P(\nu_\alpha \to \nu_\beta) = \sum_{ij} V^*_{\alpha i} V_{\beta i} V_{\alpha j}
   V^*_{\beta j}\, e^{-i\phi_{ij}}\, e^{-\sigma^2 \phi_{ij}'^2/2} ,
   \qquad \phi_{ij} = (\lambda_i - \lambda_j) L , \quad
   \phi'_{ij} = \frac{d\phi_{ij}}{d\ln E} .

The factor :math:`e^{-\sigma^2\phi'^2/2}` is the average of :math:`e^{-i\phi}` over a
Gaussian spread of width :math:`\sigma` in :math:`\ln E`, with the phase taken to first
order in the offset.  In vacuum :math:`\phi' = -\phi`, and at :math:`\sigma = 10\%` the
fraction of an interference term that survives is

=============== ======= ======= ======= ======= ======= ======= =======
phase [rad]     1       2π      10      20      30      40      50
surviving       0.995   0.82    0.61    0.14    0.011   3e-4    4e-6
=============== ======= ======= ======= ======= ======= ======= =======

So the phase average is the boxed limit where every phase runs through many cycles, the
oscillation probability itself where none does, and a smooth weighting in between.  The
spread is ``average_spread``, 0.1 by default (:data:`magnus.avgprob.AVG_PHASE_SPREAD`), and
it should be the resolution of the measurement.

Only the phases are averaged.  Mixing, the crossing amplitudes of the next sections, and the
eigenbases at the two ends of the path stay at the central energy, so a probability without
interference is returned unchanged.  Averaging the whole probability over energy would smooth
those too: on the BS05 solar curve from 0.1 to 20 MeV it moves 26 of 40 energies by more than
1e-4, and 5.2e-4 at most, where the phase average moves none.

The slopes come from Hellmann-Feynman, :math:`d\lambda_i/d\ln E = \langle v_i|\, dH/d\ln
E\, |v_i\rangle`, with :math:`dH/d\ln E` a central difference of the Hamiltonian in
:math:`\ln E` with step 1e-3: two more evaluations of :math:`H` per energy, and no further
eigendecomposition.  A slope smaller than the round-off of that difference is replaced by
minus the phase of its pair, its value in vacuum.  Without that, a pseudo-Dirac pair split by
1e-21 eV² reads a slope of 0.2 rad at 100 TeV over 100 Mpc, against a true 8e-5 rad, and its
result depends on the spread although the pair is coherent.

A Hamiltonian that does not depend on energy -- a matrix, or a function of position alone,
passed on the direct route below -- has no slope for a spread to act on, and ``average=True``
returns the limit for it.  Away from the limit, the three properties above fail:
:math:`P_{\alpha\beta} \neq P_{\beta\alpha}` in general, and CP violation survives in the
terms that do.

Every point is computed as the limit first, and returned as such, bit for bit, wherever the
phase average agrees with it to 1e-4, so a result that was right before stays exactly what it
was.

Coherence, and where the spread matters
-----------------------------------------

Whether a pair of eigenvalues has decohered is a statement about that pair, not about the
spectrum as a whole.  The phase average decides it per pair through the weight; the functions
of :mod:`magnus.avgprob` that return the limit decide it with a threshold.
:func:`magnus.avgprob.averaged_probabilities_constant_hamiltonian` groups the spectrum into
blocks of mutually coherent eigenvalues, pairs whose phase is below
:data:`magnus.avgprob.DECOHERENCE_PHASE_THRESHOLD` (:math:`2\pi`), and sums coherently inside
each block,

.. math::

   P(\nu_\alpha \to \nu_\beta) = \sum_{b} \Big|
   \sum_{i \in b} V^*_{\alpha i} V_{\beta i} \Big|^2 ,

which reduces to the boxed expression when every block is a singleton.  Until 1.1.1 this
was what ``average=True`` returned for a constant Hamiltonian, with the phase inside a block
set to zero: a pair at 1 rad was kept as fully coherent, and one at 20 rad dropped as fully
decohered, although 14% of its interference survives a 10% spread.  The distinction between
the block form and the naive sum is not academic.  A sterile state with a small
:math:`\Delta m^2_{41}`, or any exactly degenerate spectrum, makes the
naive sum quietly wrong: with *all* eigenvalues equal the correct answer
is the identity -- nothing oscillates at all -- while the naive sum
returns a spurious mixture.

Pseudo-Dirac neutrinos are the case the block form was written for.
Each mass eigenstate that carries a sterile partner splits into two states
separated by a :math:`\delta m^2` many orders below :math:`\Delta m^2_{21}`,
so over an astrophysical baseline the standard phases have long since
averaged away while every pair is still coherent -- one block per pair, and
the naive sum is wrong by a factor of two.  The phase average needs no blocks for them: a pair
whose phase barely moves with energy keeps a weight near one.
:mod:`magnus.hamiltonians.hamiltonians_pseudodirac` builds those
Hamiltonians, with the pairing selectable per mass state; notebook 29 walks
the splitting up through the three regimes below with
:func:`magnus.avgprob.coherence_report`.

The same per-pair phase decides whether the limit applies at all.  A pair
is in one of three regimes:

.. image:: _static/averaging_regimes.svg
   :width: 100%
   :align: center
   :alt: The coherent, intermediate and decohered regimes of a pair of eigenvalues

|

* far below :data:`magnus.avgprob.COHERENCE_PHASE_THRESHOLD`, the
  relative phase has barely advanced: the pair is coherent and there is
  nothing to average;
* far above :data:`magnus.avgprob.DECOHERENCE_PHASE_THRESHOLD`, the cross
  term has averaged away and the boxed expression is exact;
* **in between, neither statement holds**, and the limit does not describe
  the result.

:func:`magnus.avgprob.coherence_report` names the pairs in that middle
band.  For the phase average it is the band where the answer depends on
:math:`\sigma`, and that is what :class:`magnus.oscprob.PhaseAveragingWarning`
reports: it fires where :math:`|\sigma\, \partial P/\partial\sigma|` exceeds
:data:`magnus.avgprob.PHASE_SPREAD_SENSITIVITY_THRESHOLD`, 1e-3, the default
tolerance.  The number is then the average over the spread asked for.  Asking
for the average at a 1000 km beamline does exactly this: at 1 GeV,
:math:`P_{\mu\mu}` is 0.91 at a 10% spread and 0.97 at 5%.

Position-dependent Hamiltonians
-----------------------------------

When the Hamiltonian varies along the trajectory there is no single
eigenbasis to decohere in.  A neutrino produced at :math:`l_0` decoheres
in the eigenbasis *there*, is carried along the levels of the
instantaneous Hamiltonian, and is detected in the eigenbasis at
:math:`l_1`:

.. math::

   P(\nu_\alpha \to \nu_\beta) = \sum_{ij} |V_{\alpha i}(l_0)|^2 \,
   P^{\rm cross}_{ij} \, |V_{\beta j}(l_1)|^2 ,

the standard MSW-plus-decoherence result, generalized here to any number
of levels and any number of crossings.  :math:`P^{\rm cross}` is the
probability of ending on level :math:`j` having started on level
:math:`i`.  Adiabatic evolution keeps a neutrino on its level, so
:math:`P^{\rm cross}` is the **identity** wherever the adiabatic
approximation holds, and departs from it only across a non-adiabatic
window.

Those windows are located with the Hellmann-Feynman diagnostic of
:mod:`magnus.adiabatic` (see :doc:`adiabatic_strategy`), and the transfer
across each one is computed with that module's own convergence-checked
Magnus patch -- *not* with a Landau-Zener formula.  The result is an
exact treatment of the crossing rather than an asymptotic approximation
to it.  As a check, the computed hop probability reproduces the analytic
Landau-Zener value :math:`\exp(-2\pi\epsilon^2/|d\Delta/dl|)` to a few
parts in a thousand for a linear crossing, with nothing in the
implementation assuming that formula:

.. list-table::
   :header-rows: 1
   :widths: 20 25 25 15

   * - Coupling :math:`\epsilon`
     - Computed hop probability
     - Landau-Zener
     - Difference
   * - :math:`10^{-4}`
     - 0.9382
     - 0.9391
     - 0.1%
   * - :math:`3\times10^{-4}`
     - 0.5665
     - 0.5680
     - 0.3%

The search for windows runs on a grid of 200 probes.  A feature narrower
than their spacing falls between two of them and is never examined: until
1.1.1 no window opened there, :math:`P^{\rm cross}` came out the identity,
and the fully adiabatic answer was returned without a warning -- on a
supernova shock ray, 0.04 where the averaged probability is 0.37 to 0.59.
The profile is now checked first for features that sharp *and* able to
move probability between levels -- an instantaneous change across them
would move more than
:data:`magnus.avgprob.SUDDEN_TRANSFER_THRESHOLD`, the default tolerance.
Where there is one, the windows are taken from
:func:`magnus.adiabatic.hybrid_propagator`, which refines its search until
it certifies.  Of 24 fronts on that ray, 0.07 to 2000 km wide, 16 then
come back within 0.01 of a decohered reference, and the other 8 warn:
where the feature is a discontinuity no
refinement resolves, the call raises
:class:`magnus.oscprob.UnmarkedDiscontinuityWarning`, and the cure is to
declare it with ``t_breakpoints``, which selects the energy-window average
described below.  Everywhere else -- every solar profile measured,
tabulated models included -- nothing is escalated and the result is what
it was, bit for bit.

The phase average on a profile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The expression above composes the crossings as probabilities and reads the result out in the
eigenbasis at :math:`l_1`, so it keeps no interference at all: between two crossings, or
between the last crossing and the end of the path, every phase is dropped however short it is.
On a chord through the edge of the Sun at 10 TeV the last window ends where the path does, the
phase that follows is zero, and the expression gives 0.40 where the neutrinos leave with 0.54.
:func:`magnus.avgprob.averaged_probabilities_adiabatic`, which returns it, reports the pairs
for which that composition is unsafe.

The phase average carries its definition onto the profile.  An energy offset
:math:`u = \delta\ln E` moves every instantaneous eigenvalue by
:math:`u\, d\lambda_i/d\ln E` and leaves the eigenvectors alone; the result is the Gaussian
average over :math:`u` of the evolution under :math:`H + u\, D_{\rm diag}`, with
:math:`D_{\rm diag}` the part of :math:`dH/d\ln E` diagonal in the instantaneous eigenbasis.
:func:`magnus.avgprob.phase_averaged_probabilities_adiabatic` computes it without sampling
energies along the path.  Each non-adiabatic window is an amplitude matrix, computed with the
same Magnus patch at a few Gauss-Hermite nodes in :math:`u` (at most 31, set by how fast the
phases inside the window run with energy, and a uniform grid beyond that).  Each adiabatic
stretch between windows is a diagonal phase with an exact slope in :math:`u`, carried with the
parallel-transport phase of its eigenvectors.  The density matrix is carried as terms labelled
by accumulated slope, whose Gaussian average is analytic, and a term is dropped only once no
later stretch can bring its slope back.

Because the definition reaches inside the windows, the answer does not depend on where they are
drawn: one window over a stretch, and two windows with the stretch between them, return the same
number.  Where there is no window, a decohered start carried adiabatically has no interference
to keep, and ``average=True`` returns the expression above, bit for bit; every solar MSW curve
is of that kind.  Transfer between levels outside the windows is neglected, as in any adiabatic
calculation; the limit drops the interference such a transfer carries and the phase average
keeps it, so the windows are searched at an adiabaticity threshold of 0.01
(:data:`magnus.avgprob.PHASE_AVERAGE_WINDOW_THRESHOLD`) rather than 0.1.  Against a brute-force
average of the same definition on five solar chords from 10 GeV to 10 TeV, the result is within
4.2e-5; on the same chords the limit is off by up to 0.14.

When there is no closed form
--------------------------------

A profile with discontinuities -- the PREM layer boundaries an
Earth-crossing trajectory steps through -- has no instantaneous
eigenbasis to decohere in, so neither construction above applies.  There,
``average=True`` propagates the probability for real across an energy
window and averages over it
(:func:`magnus.avgprob.averaged_probabilities_numerically`).

**This returns a different quantity from the other two paths.**  They
return the phase average, which weights phases and samples no energy;
this returns the probability averaged over one particular window -- a
top-hat of 41 energies, from a flavor state at the start -- and the
answer depends on its width.  The default width,
:data:`magnus.avgprob.AVG_DEFAULT_ENERGY_SPREAD`, is 10% -- the order of
a real detector's energy resolution -- and every use of it raises
:class:`magnus.oscprob.PhaseAveragingWarning` naming the width, the
number of samples, and the standard error of the resulting mean, so the
figure is never silently dependent on a constant the caller did not
choose.  Callers with a known resolution should pass their own.

Cost
-------

.. list-table::
   :header-rows: 1
   :widths: 40 32 28

   * - Case
     - Method
     - Cost
   * - Vacuum, constant density (and their NSI/LIV variants)
     - Closed form
     - 10 :math:`\mu`\ s per energy where every phase has decohered, 23 :math:`\mu`\ s
       where some has not
   * - Exponential density, Sun (and their NSI/LIV variants), no window
     - Adiabatic transport
     - 6 ms per energy
   * - A profile with windows whose phases survive: a solar chord, 10 GeV to 1 TeV
     - Windows at Gauss-Hermite nodes
     - 5.9 s (10 GeV) to 0.35 s (1 TeV) per point
   * - Earth (PREM)
     - Sampled over an energy window
     - ~0.1 s

Measured on an idle machine with NuFIT 6.1, the solar rows on BS05 and the chords on
B16-GS98 at :math:`b = 0.6\,R_\odot`.  For comparison, averaging the engine numerically over
2001 energies for the vacuum number takes about 0.25 s, and is itself an approximation.

Usage
--------

.. jupyter-execute::

    import numpy as np

    import magnus.oscprob as oscprob
    import magnus.globaldefs as gd

    # load_nufit_params returns just the six mixing parameters; the
    # OSC_PARAMS_PREDEFINED entries also carry 'name'/'description' strings,
    # which the propagation machinery would reject.
    osc = gd.load_nufit_params('NuFIT 6.1')

    # 1 TeV from 100 Mpc: every phase has decohered, and this is the limit
    P = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_TEV, 3.0857e21*gd.UNIT_KM,
                                    average=True, **osc)
    np.round(np.asarray(P), 4)

The most quoted consequence of averaged astrophysical oscillations
follows in one line: a source producing the pion-decay composition
:math:`(1:2:0)` delivers something close to equipartition at Earth.

.. jupyter-execute::

    at_source = np.array([1.0, 2.0, 0.0])/3.0
    at_earth = at_source @ np.asarray(P)
    np.round(at_earth*3.0, 3)

At a terrestrial baseline the phases have not decohered, and the answer depends on the spread:

.. jupyter-execute::

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')     # PhaseAveragingWarning: it depends on the spread
        P_mumu = [float(oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM,
                        nu_i=gd.NUMU, nu_f=gd.NUMU, average=True, average_spread=s, **osc))
                  for s in (0.05, 0.1, 0.3)]
    np.round(P_mumu, 4)

The keyword is not confined to the wrappers.  On the direct route, a
Hamiltonian of your own goes through ``osc_prob_energy_baseline`` (or through
``osc_prob_earth`` and ``osc_prob_sun``, which add the geometry), and
``average=True`` there takes the same three routes the wrappers take: the closed
form when the Hamiltonian does not depend on position, adiabatic transport when
it does and the profile is smooth, and an energy-window average when
``t_breakpoints`` or ``t_slab_edges`` declare discontinuities.  The Hamiltonian
is passed as ``H(E, l)``, or as ``H(E)`` with
``H_func_is_function_only_of_energy=True``, or as a matrix; a matrix, or a
function of position alone, does not depend on energy and returns the limit.

.. jupyter-execute::

    import magnus.hamiltonians as hamiltonians
    import magnus.matter as matter

    h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**osc)
    proj = matter.matter_potential_projector(3)
    vcc = matter.vcc_func_from_rho_func(
        lambda l: 5.0*np.exp(-l/(1000.0*gd.UNIT_KM)), 0.0, 1.0, 0.5,
        nubar=False, density_matter_is_in_g_per_cm3=True,
        density_is_of_number_of_electrons=False)

    def H(E, l):
        return h_vac/E + np.asarray(vcc(l))[..., None, None]*proj

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')     # this one, too, depends on the spread
        P = oscprob.osc_prob_energy_baseline(
            H, 0.1*gd.UNIT_GEV, 5000.0*gd.UNIT_KM, 0.0,
            nu_i=gd.NUE, nu_f=gd.NUE, average=True)
    round(float(P), 4)

``osc_prob`` itself, which computes one point, does not take the keyword and
says so if handed it.

Am I computing the wrong thing?  ``strategy_info['sampling']``
----------------------------------------------------------------

The hardest part of this page in practice is not the mathematics -- it is
noticing that it applies to you.  A scan of instantaneous probabilities
over a long trajectory returns perfectly correct numbers, and they can
still be the wrong quantity, because the observable is an average over a
phase nobody resolves.

Every ``osc_prob_*`` entry point that accepts ``strategy_info`` now
reports how coarsely the request samples the oscillation it is
computing::

    info = {}
    P = magnus.oscprob.osc_prob_3nu_sun(energy, L, info_kwargs..., strategy_info=info)
    info['sampling']
    # {'oscillation_length': 2.53e+10,   'cycles_over_trajectory': 1.32e+04,
    #  'spacing': 3.82e+13,              'cycles_per_step': 1.51e+03,
    #  'nyquist_points': 26446,          'aliased': True}

``cycles_per_step`` is the number to read.  Above about 0.5 the scan
takes less than two samples per oscillation, so the returned array
**cannot represent the oscillation** and must not be plotted or
interpolated as a curve -- the individual values are right, the curve
through them is an artifact.  ``nyquist_points`` says how many baselines
would be needed to sample it properly.

Those numbers are usually stark.  Measured over the physically-motivated
profile families in ``docs/dev/adversarial_batteries/``:

=========================== ============================ =========================
trajectory                  oscillations across it        baselines for Nyquist
=========================== ============================ =========================
Earth chord                 ~430                          861
Solar, one scale height     ~2200                         4 390
Supernova ray               ~37 000                       73 392
=========================== ============================ =========================

**This is reported and never warned about, deliberately.**  The Nyquist
criterion is objectively correct and would fire on 44 of 45 realistic
scan sizes -- a warning firing on 98 % of calls is noise however right
each firing is, and it would teach users to silence a category that also
carries genuine discontinuity warnings.  The measurement behind that
decision is ``adversarial_batteries/alias_fp.py``.

The report costs eigenvalues at eight points along the trajectory, so it
is computed **only when ``strategy_info`` is supplied**: callers who do
not ask pay nothing, and callers who do pay 5.5 % of the cheapest scan
measured and under 0.1 % of a substantial one.

When ``aliased`` is ``True``, the question worth asking is whether you
wanted the average all along.  If you did, ``average=True`` gives it
without resolving the oscillation, with the spread set by
``average_spread``, and :class:`magnus.oscprob.PhaseAveragingWarning` says
where the answer depends on that spread.

How much does the phase actually matter?
-------------------------------------------

It depends on the profile, and the difference is measurable rather than a
matter of taste.  Averaging an instantaneous scan over six oscillation
lengths and comparing against a ``solve_ivp`` reference
(``adversarial_batteries/avg_check.py`` and ``avg_check2.py``).  The solar
row is the log-linear interpolant of the BS05 table, which is the profile
notebook 13 works from; ``avg_check.py`` prints a cubic-spline variant of
the same ray beside it, and that one reads 8.889e-04 and 6.051e-04 for the
two columns --- a different profile, and the same verdict:

=============================== ================== ================== ====================
configuration                   instantaneous      averaged           averaged inside 1e-3
=============================== ================== ================== ====================
Solar model, d = 2, 5 MeV       6.000e-04          7.110e-04          yes
Supernova turbulence, 45 MeV    5.584e-03          **3.843e-04**      yes
Supernova shock, 70 km front    4.917e-04          **2.151e-04**      yes
Supernova shock, 0.07 km front  1.988e-01          **2.222e-01**      **no**, and warned
=============================== ================== ================== ====================

The last row is the one that matters, and it is the only one where the
*observable* is wrong.  A shock front changes the adiabaticity of the level
crossing, so it moves the conversion probability itself rather than the
phase at which it oscillates; that is an **envelope** error and no
averaging operation removes it.  Everywhere else the averaged answer lands
inside the target even where a single baseline does not, because the
instantaneous error is largely **phase** -- the profile perturbs *when* the
oscillation is, and no observable resolves that.

.. warning::

   **Do not read the ratio of these two columns as a diagnostic.**  Both are
   finite-window means, and such a mean is an estimator with a bias of its
   own.  On a profile whose density varies across the averaging window the
   bias does not shrink as the window widens, because a wider window also
   averages over different matter conditions: on the solar ray above, the
   window mean moves from 0.5924 to 0.6023 between six and forty-eight
   oscillation lengths, drifting away from rather than towards a limit.
   Notebook 13 prints that sweep.
   The reduction factor is meaningful only on a controlled comparison at
   fixed matter conditions, as in notebook 23.

   To obtain the averaged probability, ask for it rather than estimating
   it.  ``average=True`` evaluates it with no window to choose: the solar
   ray has no non-adiabatic window, so there the phase average is the
   decohered limit, and it reproduces the adiabatic MSW expression
   :math:`\langle P_{ee}\rangle = \tfrac12 + \tfrac12\cos2\theta_m(L_0)\cos2\theta_m(L_1)`
   to 3e-16 across 1--20 MeV.

Limitations and scope
-------------------------

* The phase average is a statement about a **measurement with an energy
  resolution**.  It is not a model of quantum decoherence: there is no
  dissipative term here.  See the "When is Magνs not the right tool?"
  section of :doc:`index`.
* A spread in geometry -- the size of the production region, the range of
  impact parameters a pixel covers -- is not modeled.  A phase that does not
  depend on energy, as the matter phases on a solar chord at TeV energies
  do, is kept whatever its size, and a figure has to average over its own
  geometry.
* The Earth/PREM path is a windowed average of the probability, as
  described above.

See :doc:`functions` and the API reference for the full listing of
:mod:`magnus.avgprob`.
