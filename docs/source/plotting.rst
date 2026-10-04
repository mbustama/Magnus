Pre-packaged plotting tools
===========================

This page documents :mod:`magnus.plotting`, a small set of functions that
draw the figures of the example notebooks, each in one call.  See
:doc:`tutorials` for the notebooks.

What it draws
---------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Function
     - Layout
   * - :func:`~magnus.plotting.plot_curves`
     - Curves against any swept variable, with an optional relative-error
       subpanel: the shape of most of the notebooks' figures.  Also serves the
       convergence and error studies.
   * - :func:`~magnus.plotting.plot_probability_vs_baseline`
     - The same, preset to a logarithmic baseline axis and a unit ordinate.
   * - :func:`~magnus.plotting.plot_probability_vs_energy`
     - The same, against neutrino energy.
   * - :func:`~magnus.plotting.plot_curves_stacked`
     - Small multiples: one panel per case, stacked over a shared abscissa,
       with the axes forced to match so the panels can be read against each
       other.
   * - :func:`~magnus.plotting.plot_probability_with_profile`
     - Matter-density panel above one or more probability panels sharing an
       abscissa.  Given trajectories instead of curves, it computes the
       probabilities through the Earth wrappers, against baseline or energy,
       and draws the electron density they integrate.
   * - :func:`~magnus.plotting.plot_probability_with_average`
     - An oscillating probability with its phase average overlaid; see
       :doc:`averaged_probability`.
   * - :func:`~magnus.plotting.plot_biprobability`
     - Neutrino against antineutrino appearance probability, as
       :math:`\delta_{\rm CP}` runs over its range.  Given configurations
       (a path, an energy, parameters) instead of probabilities, it computes
       them through the Earth wrappers, markers at named phases included.
   * - :func:`~magnus.plotting.plot_oscillogram`
     - Probability over zenith angle and energy, as a filled contour map.
       Given no probability, it computes one through the Earth wrappers, with
       the PREM layer boundaries declared and a layered electron fraction that
       can be overridden.

Installation
--------------

Nothing to do: Matplotlib is a dependency of Magνs, so ``pip install magnuspy``
brings it and :mod:`magnus.plotting` is available straight away.

.. code-block:: python

    from magnus import plotting

A first figure
----------------

.. jupyter-execute::

    import numpy as np
    import magnus.globaldefs as gd
    import magnus.oscprob as oscprob
    from magnus.plotting import plot_probability_vs_baseline

    osc = gd.load_nufit_params('NuFIT 6.1')
    energy = 1.0 * gd.UNIT_GEV
    distances = np.logspace(1.0, 4.0, 300)                      # [km]

    prob = np.asarray(oscprob.osc_prob_3nu_vacuum(
        np.full(distances.size, energy), distances*gd.CONV_KM_TO_INV_EV,
        **osc))[:, gd.NUMU, gd.NUMU]

    fig, ax = plot_probability_vs_baseline(
        distances,
        [dict(y=prob, label='Magnus expansion', color='C1')],
        nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=3,
        xlim=(distances[0], distances[-1]),
        title=r'$3\nu$ vacuum, $E_\nu = 1$ GeV',
        legend_title='Calculation method',
    )

Note that ``osc`` comes from :func:`magnus.globaldefs.load_nufit_params`,
which returns exactly the six mixing parameters.  Splatting
``gd.OSC_PARAMS_PREDEFINED[...]`` instead would also forward its ``name`` and
``description`` strings, which the probability functions reject.

Adding the error subpanel
---------------------------

Passing ``residual`` adds the short lower panel the notebooks use to compare a
Magnus result against a closed-form one.  The two panels then share their
abscissa limits and scale, and the upper panel's tick labels are suppressed,
so they read as a single figure:

.. jupyter-execute::

    import magnus.oscprobstd as oscprobstd
    from magnus.plotting import plot_curves

    sth, Dm2 = osc['s12'], osc['D21']
    energy = 10.0 * gd.UNIT_MEV
    L = np.logspace(1.0, 5.0, 400)
    L_nat = L * gd.CONV_KM_TO_INV_EV

    # One call each. The point axis lands last in the closed form and first in
    # osc_prob, which is the only difference between the two lines.
    exact = np.asarray(
        oscprobstd.osc_prob_2nu_vacuum_std(sth, Dm2, energy, L_nat))[0, 0]
    approx = np.asarray(oscprob.osc_prob_2nu_vacuum(
        np.full(L.size, energy), L_nat, sth, Dm2))[:, 0, 0]

    fig, ax = plot_curves(
        L,
        [dict(y=approx, label='Magnus expansion', color='C1'),
         dict(y=exact, label='Standard formula', color='k', ls='--')],
        xlabel=r'Baseline, $L$ [km]',
        ylabel=r'Two-neutrino probability, $P_{\nu_e \to \nu_e}$',
        xlim=(L[0], L[-1]), ylim=(0.0, 1.0), xscale='log',
        ymajor=0.10, yminor=0.02,
        residual=(approx - exact) / np.maximum(exact, 1.0e-300) / 1.0e-12,
        residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
        legend_title='Calculation method',
    )

The other layouts
-------------------

Three of the functions can compute the probabilities themselves, through the
Earth wrappers: :func:`~magnus.plotting.plot_probability_with_profile`,
:func:`~magnus.plotting.plot_biprobability` and
:func:`~magnus.plotting.plot_oscillogram`.  Leave out the probability, and give
``num_flavors``, the channel (``nu_i``, ``nu_f``), the energy where the plot does not
sweep it, and optionally ``osc_params`` and ``wrapper_kw`` (any other keyword of the
Earth wrapper, such as ``nubar``).  The docstrings list the keywords each one reserves
for itself.

**Probability against energy**, from a probability you computed:

.. jupyter-execute::

    from magnus.plotting import plot_probability_vs_energy

    E = np.linspace(0.5, 5.0, 300)                              # [GeV]
    P_mue = np.asarray(oscprob.osc_prob_3nu_matter_constant_density(
        E*gd.UNIT_GEV, 1300.0*gd.UNIT_KM, 2.848*gd.UNIT_G_PER_CM3,
        nu_i=gd.NUMU, nu_f=gd.NUE, **osc))

    fig, ax = plot_probability_vs_energy(
        E, [dict(y=P_mue, label='DUNE, 2.848 g/cm$^3$')],
        nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3, xscale='linear',
        xlim=(E[0], E[-1]))

``energy_unit`` sets only the label (GeV by default): pass the energies already in
that unit.

**Density profile over the probability**, computed along a chord through the core:

.. jupyter-execute::

    from magnus.plotting import plot_probability_with_profile

    L_km = np.linspace(100.0, 11000.0, 200)      # the chord at cos = -0.9 is 11 468 km
    fig, ax = plot_probability_with_profile(
        L_km, trajectories=[dict(costhz=-0.9, label=r'$\cos\theta_z = -0.9$')],
        energy=5.0*gd.UNIT_GEV, nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3,
        xscale='linear')

**CP violation as a bi-probability plot**, Fermilab to Homestake at 2 GeV, both
orderings:

.. jupyter-execute::

    from magnus.plotting import plot_biprobability

    where = dict(loc_ini='fermilab', loc_fin='homestake')
    fig, ax = plot_biprobability(
        configurations=[dict(where, **gd.load_nufit_params('NuFIT 6.1', 'NO')),
                        dict(where, **gd.load_nufit_params('NuFIT 6.1', 'IO'))],
        dcp=np.linspace(-np.pi, np.pi, 25), energy=2.0*gd.UNIT_GEV, num_flavors=3,
        labels=['NO', 'IO'], markers=[dict(dcp=0.0, marker='*', label='0')])

**An oscillogram**, over :math:`\cos\theta_z` and :math:`\log_{10}(E/{\rm GeV})`.
A probability you pass must have shape ``(len(log10_energy), len(costhz))``: one
row per energy.

.. jupyter-execute::

    import warnings
    from magnus.magnus import MagnusConvergenceWarning
    from magnus.plotting import plot_oscillogram

    cos_grid = np.linspace(-1.0, -0.1, 25)
    log10_E = np.linspace(0.0, 1.0, 20)
    with warnings.catch_warnings():      # expected on a few chords; see diagnostics
        warnings.simplefilter('ignore', MagnusConvergenceWarning)
        fig, ax = plot_oscillogram(cos_grid, log10_E, nu_i=gd.NUMU, nu_f=gd.NUE,
                                   num_flavors=3)

**Small multiples** (:func:`~magnus.plotting.plot_curves_stacked`), one panel per
case with matching axes, and **a probability under its phase average**
(:func:`~magnus.plotting.plot_probability_with_average`) take probabilities you
computed; their docstrings have an example each.

.. _plotting-api-conventions:

API conventions
-----------------

Named arguments, and no catch-all
""""""""""""""""""""""""""""""""""""

Each function takes **named arguments for the things every figure has** --
data, labels, limits, scales, tick spacings, title, legend placement, output
path -- and **explicit pass-through dictionaries for the long tail** of
Matplotlib settings: ``legend_kw``, ``grid_kw``, ``savefig_kw``,
``subplots_kw``, and, per curve, any
:class:`~matplotlib.lines.Line2D` keyword.

No function accepts unknown keywords, so a misspelled keyword raises an error.  Every
keyword either appears in a signature, where a misspelling raises a :class:`TypeError`,
or goes into a dictionary for one Matplotlib call, which then raises an error naming
the key.  The three functions that take extra keywords
(``plot_probability_vs_baseline``, ``plot_probability_vs_energy`` and
``plot_probability_with_average``) pass them to a function that applies the same
check:

.. jupyter-execute::

    try:
        plot_curves(L, [dict(y=approx)], ylabl='typo')
    except TypeError as error:
        print(type(error).__name__, '->', error)

Curves
"""""""""

``curves`` is a sequence, one entry per line.  An entry is either a bare
ordinate array or a dictionary carrying the ordinate under ``'y'`` plus any
Line2D keyword.  Entries without an explicit color take the ``'C0'``,
``'C1'``, ... cycle in order, matching the notebooks; reference curves are
conventionally given ``color='k', ls='--'``.

Returning ``(fig, ax)``
""""""""""""""""""""""""""

Every function returns both, so that the figure can be edited further.  Each
call creates its own figure: there is no
``ax=`` argument for drawing into existing axes.  ``subplots_kw`` cannot set
``nrows`` or ``ncols``, which the layout fixes, nor ``figsize``, which has its own
``figsize=`` argument.  ``ax`` is a single
:class:`~matplotlib.axes.Axes` for the single-panel layouts and an array for
the multi-panel ones -- with a residual subpanel, ``ax[0]`` is the main panel
and ``ax[1]`` the residual:

.. jupyter-execute::

    fig, ax = plot_curves(L, [dict(y=approx, label='m')], xscale='log')
    ax.axvline(1.0e3, color='0.6', ls=':', lw=1)
    ax.set_title('annotated after the fact', fontsize=20)

House style
"""""""""""""

The module sets the text and tick sizes it draws with
(:data:`~magnus.plotting.HOUSE_RC`: 25-point axis labels, 23-point tick labels, ticks
pointing in on all four sides), so a figure looks the same in a script, a notebook or
these pages; it reads no ``matplotlibrc``.  A size you have changed yourself, in
``rcParams`` or a style, is kept.  Font family and LaTeX rendering follow your
Matplotlib settings.  The other house values are exposed as
:data:`~magnus.plotting.HOUSE_FIGSIZE`,
:data:`~magnus.plotting.HOUSE_LEGEND_KW`,
:data:`~magnus.plotting.HOUSE_GRID_KW` and
:data:`~magnus.plotting.HOUSE_SAVEFIG_KW`, so a user can build on them
rather than restate them.

Labels
--------

:func:`~magnus.plotting.prob_label` builds the LaTeX for a probability from a
flavor pair, the sterile states included:

.. jupyter-execute::

    from magnus.plotting import prob_label

    print(prob_label(gd.NUMU, gd.NUE))
    print(prob_label(gd.NUMU, gd.NUE, nubar=True))
    print(prob_label(gd.NUE, gd.NUS1))

Saving
--------

Pass ``savefig`` to write the figure; ``savefig_kw`` is merged over
:data:`~magnus.plotting.HOUSE_SAVEFIG_KW`, which is ``dpi=200``.  The
notebooks write PDFs into ``../fig/``, whose contents are ignored by git.
