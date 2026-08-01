Pre-Packaged Figures
======================

This page documents :mod:`magnus.plotting`, a small set of functions that
produce the figures the example notebooks use, so that a plot costs one call
rather than thirty lines of Matplotlib.  It is the only part of Magνs that
needs Matplotlib, and the only one behind an optional dependency; see
:doc:`tutorials` for the notebooks it was extracted from.

Why it exists
---------------

Every figure across the twelve notebooks used to be built by hand.  A typical
one ran to twenty-five or forty lines: a ``gridspec_kw`` dictionary, a
``subplots_adjust`` call, the plotting itself, a nine-keyword ``legend``
invocation, four ``MultipleLocator`` assignments, axis limits and scales, and
a ``savefig``.  That block was copied from figure to figure and varied
slightly each time, which is exactly the arrangement in which a figure
quietly stops matching its neighbours.

Cataloguing the figures first turned out to matter more than writing the code.
Most of them are the *same figure* with different data.  Probability against
baseline, probability against energy, probability against a sterile mixing
angle, and the convergence studies of the matrix-exponential notebook all
share one shape:

   a set of curves against a swept variable, optionally over a short
   relative-error subpanel.

:func:`magnus.plotting.plot_curves` is that shape.
:func:`~magnus.plotting.plot_probability_vs_baseline` and
:func:`~magnus.plotting.plot_probability_vs_energy` are thin presets over it
that fix labels, scales and tick spacings.  Only three layouts are genuinely
distinct: the density profile stacked over one or more probability panels,
the bi-probability plane, and the oscillogram.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Function
     - Layout
   * - :func:`~magnus.plotting.plot_curves`
     - Curves against any swept variable, with an optional relative-error
       subpanel.  Also serves the convergence and error studies.
   * - :func:`~magnus.plotting.plot_probability_vs_baseline`
     - The same, preset to a logarithmic baseline axis and a unit ordinate.
   * - :func:`~magnus.plotting.plot_probability_vs_energy`
     - The same, against neutrino energy.
   * - :func:`~magnus.plotting.plot_probability_with_profile`
     - Matter-density panel above one or more probability panels sharing an
       abscissa.
   * - :func:`~magnus.plotting.plot_probability_with_average`
     - An oscillating probability with its phase-averaged limit overlaid; see
       :doc:`averaged_probability`.
   * - :func:`~magnus.plotting.plot_biprobability`
     - Neutrino against antineutrino appearance probability, as
       :math:`\delta_{\rm CP}` runs over its range.
   * - :func:`~magnus.plotting.plot_oscillogram`
     - Probability over zenith angle and energy, as a filled contour map.

Installation
--------------

Matplotlib is **not** a core dependency: the engine needs only NumPy, SciPy
and joblib, and a user computing probabilities for their own analysis code
should not have to install a plotting stack.  It is declared as the ``plot``
extra instead::

    pip install 'magnuspy[plot]'

or, from a source checkout::

    pip install -e '.[plot]'

``import magnus`` works either way.  :mod:`magnus.plotting` imports cleanly
without Matplotlib too -- it defers the import to the calls that actually
draw -- so only a plotting call raises, and it raises
:class:`~magnus.plotting.MatplotlibNotFoundError` (a subclass of
:class:`ImportError`) naming the command above.

A first figure
----------------

.. jupyter-execute::

    import matplotlib
    matplotlib.use('Agg')

    import numpy as np
    import magnus.globaldefs as gd
    import magnus.oscprob as oscprob
    from magnus.plotting import plot_probability_vs_baseline

    osc = gd.load_nufit_params('NuFIT 6.1')
    energy = 1.0 * gd.UNIT_GEV
    distances = np.logspace(1.0, 4.0, 300)                      # [km]

    prob = np.array([
        oscprob.osc_prob_3nu_vacuum(energy, L * gd.CONV_KM_TO_INV_EV,
                                    **osc)[gd.NUMU][gd.NUMU]
        for L in distances])

    fig, ax = plot_probability_vs_baseline(
        distances,
        [dict(y=prob, label='Magnus expansion', color='C1')],
        nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=3,
        xlim=(distances[0], distances[-1]),
        title=r'$3\nu$~vacuum, $E_\nu = 1$~GeV',
        legend_title='Calculation method',
    )
    print(ax.get_ylabel())

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

    sth, Dm2 = gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0
    energy = 10.0 * gd.UNIT_MEV
    L = np.logspace(1.0, 5.0, 400)

    exact = np.array([
        oscprobstd.osc_prob_2nu_vacuum_std(sth, Dm2, energy,
                                           l * gd.CONV_KM_TO_INV_EV)[0][0]
        for l in L])
    approx = np.array([
        oscprob.osc_prob_2nu_vacuum(energy, l * gd.CONV_KM_TO_INV_EV,
                                    sth, Dm2)[0][0]
        for l in L])

    fig, ax = plot_curves(
        L,
        [dict(y=approx, label='Magnus expansion', color='C1'),
         dict(y=exact, label='Standard formula', color='k', ls='--')],
        xlabel=r'Baseline, $L$ [km]',
        ylabel=r'Two-neutrino probability,~$P_{\nu_e \to \nu_e}$',
        xlim=(L[0], L[-1]), ylim=(0.0, 1.0), xscale='log',
        ymajor=0.10, yminor=0.02,
        residual=(approx - exact) / np.maximum(exact, 1.0e-300) / 1.0e-12,
        residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
        legend_title='Calculation method',
    )
    print(f'panels: {len(ax)}')

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

There is deliberately no bare ``**kwargs`` on any of them.  A catch-all
signature accepts a misspelled keyword in silence, and this project has
already paid for that: ``oscprob``'s keyword chain used to forward unknown
names down several layers before failing somewhere unrecognisable.  Here
every keyword either appears in the signature, so a typo is a
:class:`TypeError` at the call site, or lands in a dictionary destined for one
specific Matplotlib call, so a typo is an error from that call naming the
offending key.  Nothing is swallowed:

.. jupyter-execute::

    try:
        plot_curves(L, [dict(y=approx)], ylabl='typo')
    except TypeError as error:
        print(type(error).__name__, '->', error)

Curves
"""""""""

``curves`` is a sequence, one entry per line.  An entry is either a bare
ordinate array or a dictionary carrying the ordinate under ``'y'`` plus any
Line2D keyword.  Entries without an explicit colour take the ``'C0'``,
``'C1'``, ... cycle in order, matching the notebooks; reference curves are
conventionally given ``color='k', ls='--'``.

Returning ``(fig, ax)``
""""""""""""""""""""""""""

Every function returns both, so that a pre-packaged figure is a starting
point rather than a dead end.  ``ax`` is a single
:class:`~matplotlib.axes.Axes` for the single-panel layouts and an array for
the multi-panel ones -- with a residual subpanel, ``ax[0]`` is the main panel
and ``ax[1]`` the residual:

.. jupyter-execute::

    fig, ax = plot_curves(L, [dict(y=approx, label='m')], xscale='log')
    ax.axvline(1.0e3, color='0.6', ls=':', lw=1)
    ax.set_title('annotated after the fact', fontsize=20)
    print(ax.get_title())

What is *not* set here
""""""""""""""""""""""""

Global styling -- fonts, tick sizes and directions, LaTeX rendering -- lives
in ``notebooks/matplotlibrc`` and is picked up because the notebooks run from
that directory.  This module sets only what the notebooks were overriding per
figure: figure size, the legend keyword block, gridspec ratios, tick
spacings.  The house values are exposed as
:data:`~magnus.plotting.HOUSE_FIGSIZE`,
:data:`~magnus.plotting.HOUSE_LEGEND_KW`,
:data:`~magnus.plotting.HOUSE_GRID_KW` and
:data:`~magnus.plotting.HOUSE_SAVEFIG_KW`, so a caller can build on them
rather than restate them.

Labels
--------

:func:`~magnus.plotting.prob_label` builds the LaTeX for a probability from a
flavour pair.  A helper of this name was defined separately in several
notebooks, each covering only the three active flavours; this one also covers
the sterile states, so the sterile-neutrino notebook can use it:

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
