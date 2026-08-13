# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""Pre-packaged figures for Mag :math:`\nu` s.

Every figure in ``notebooks/`` used to be built by hand: roughly 25--40 lines of
``gridspec_kw``, tick locators, legend keywords and ``savefig`` per plot,
copy-pasted and lightly varied. This module collapses that into one call per
figure while reproducing the same output, so that switching an existing figure
over to it leaves the figure unchanged.

Taking stock of the fifty-odd figures showed that most of them are the *same*
figure with different data. Curves against baseline, curves against energy,
curves against a mixing angle, and the convergence studies of the
matrix-exponential notebook all share one shape: a set of curves plotted
against a swept variable, optionally over a short relative-error subpanel.
:func:`plot_curves` is that shape, and the ``plot_probability_vs_*`` helpers
are thin wrappers that only preset labels and limits.
:func:`plot_curves_stacked` is its small-multiples form: the same plot repeated
once per case down a shared abscissa, where the comparison is between panels.
The genuinely distinct layouts are the profile-plus-probability stack, the
bi-probability plane, and the oscillogram.

API conventions
---------------
The functions take **named arguments for the quantities every figure has**
(data, labels, limits, scales, ticks, title, legend placement, output path) and
**explicit pass-through dictionaries for the long tail** of Matplotlib
settings: ``legend_kw``, ``grid_kw``, ``savefig_kw``, ``subplots_kw``, and per
curve any :class:`~matplotlib.lines.Line2D` keyword.

There is deliberately **no bare** ``**kwargs`` **on any of these functions**. A
catch-all signature accepts a misspelled keyword in silence, and this project
has already paid for that once: ``oscprob``'s keyword chain forwarded unknown
names down several layers before failing somewhere unrecognizable, which is why
:func:`magnus.oscprob.osc_prob` now raises on stray keys. Here every keyword is
either named in the signature -- so a typo is a :class:`TypeError` at the call
site -- or lands in a dictionary destined for one specific Matplotlib call --
so a typo is an error from that call, naming the offending key. Nothing is
swallowed.

Styling that is global (fonts, tick sizes and directions, LaTeX rendering)
belongs to ``notebooks/matplotlibrc`` and is deliberately **not** set here; the
defaults below cover only what the notebooks were overriding per figure.

All functions return ``(fig, ax)`` so that the caller can keep customizing:
``fig`` for figure-level work and saving, ``ax`` for anything Matplotlib
exposes on an axes.

Requirements
------------
Matplotlib ships with Magνs, so this module is available in any installation
and needs nothing extra.

.. versionadded:: 1.0.0
"""

from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    'MatplotlibNotFoundError',
    'HOUSE_FIGSIZE',
    'HOUSE_LEGEND_KW',
    'HOUSE_GRID_KW',
    'HOUSE_SAVEFIG_KW',
    'HOUSE_RESIDUAL_HEIGHT',
    'prob_label',
    'plot_curves',
    'plot_curves_stacked',
    'plot_probability_vs_baseline',
    'plot_probability_vs_energy',
    'plot_probability_with_profile',
    'plot_probability_with_average',
    'plot_biprobability',
    'plot_oscillogram',
]


class MatplotlibNotFoundError(ImportError):
    r"""Raised when :mod:`magnus.plotting` is used without Matplotlib installed.

    .. versionadded:: 1.0.0
    """


_MPL_HINT = (
    "magnus.plotting requires Matplotlib, which ships with Magnus, so this "
    "means it has been removed from the environment. Reinstall it "
    "with:\n\n    pip install matplotlib"
)


def _mpl():
    r"""Import Matplotlib on demand, with an actionable error if it is absent.

    Imported lazily rather than at module scope so that the error is raised by
    the call the user actually made, and so that merely listing the package's
    submodules does not require Matplotlib.

    Returns
    -------
    tuple
        The ``(matplotlib, matplotlib.pyplot)`` modules.

    Raises
    ------
    MatplotlibNotFoundError
        If Matplotlib cannot be imported.

    .. versionadded:: 1.0.0
    """
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    except ImportError as error:                      # pragma: no cover
        raise MatplotlibNotFoundError(_MPL_HINT) from error
    return mpl, plt


# --------------------------------------------------------------------------
# House style: exactly what the notebooks were overriding per figure.
# --------------------------------------------------------------------------

HOUSE_FIGSIZE: Tuple[float, float] = (18.0, 9.0)
r"""Default figure size, in inches, used throughout the notebooks.

.. versionadded:: 1.0.0
"""

HOUSE_RESIDUAL_HEIGHT: float = 0.3
r"""Height of the relative-error subpanel, relative to the main panel.

.. versionadded:: 1.0.0
"""

HOUSE_LEGEND_KW: Dict[str, Any] = {
    'fontsize': 17,
    'frameon': True,
    'handlelength': 1.2,
    'handleheight': 0.7,
    'borderpad': 0.8,
    'title_fontsize': 20,
    'edgecolor': 'k',
    'labelspacing': 0.7,
    'ncol': 1,
}
r"""Legend keywords repeated verbatim on essentially every notebook figure.

.. versionadded:: 1.0.0
"""

HOUSE_GRID_KW: Dict[str, Any] = {'visible': True, 'c': '0.8', 'which': 'both'}
r"""Grid keywords used by the notebooks.

.. versionadded:: 1.0.0
"""

HOUSE_SAVEFIG_KW: Dict[str, Any] = {'dpi': 200}
r"""Default :func:`~matplotlib.pyplot.savefig` keywords; figures go to ``../fig/`` as PDF.

.. versionadded:: 1.0.0
"""

_HSPACE = 0.05
_WSPACE = 0.05

# Flavor index -> LaTeX, covering the sterile states used by the 4nu/5nu notebooks.
_FLAVOR_TEX = {
    0: r'\nu_e',
    1: r'\nu_\mu',
    2: r'\nu_\tau',
    3: r'\nu_{s_1}',
    4: r'\nu_{s_2}',
}


def prob_label(nu_i: int, nu_f: int, nubar: Optional[bool] = False) -> str:
    r"""Return the LaTeX label for an oscillation probability.

    A ``prob_label`` helper was defined separately in several notebooks, with
    hand-written ``if``/``elif`` chains covering only the three active flavors.
    This version also covers the sterile states, so the sterile-neutrino
    notebook can use it too.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    nu_i : int
        Initial flavor, as one of the ``magnus.globaldefs`` constants
        ``NUE``, ``NUMU``, ``NUTAU``, ``NUS1``, ``NUS2``.
    nu_f : int
        Final flavor, same encoding.
    nubar : bool, optional
        If ``True``, label the antineutrino channel. Default is ``False``.

    Returns
    -------
    str
        A LaTeX string such as ``'$P_{\\nu_e \\to \\nu_\\mu}$'``.

    Raises
    ------
    ValueError
        If either flavor index is not one of the known values.

    Examples
    --------
    .. jupyter-execute::

        import magnus.globaldefs as gd
        from magnus.plotting import prob_label

        print(prob_label(gd.NUMU, gd.NUE))
        print(prob_label(gd.NUMU, gd.NUE, nubar=True))
    """
    for name, value in (('nu_i', nu_i), ('nu_f', nu_f)):
        if value not in _FLAVOR_TEX:
            raise ValueError(
                f'Error in magnus: plotting.prob_label: {name} must be one of '
                f'{sorted(_FLAVOR_TEX)} (the flavor constants in '
                f'magnus.globaldefs), not {value!r}'
            )
    ini, fin = _FLAVOR_TEX[nu_i], _FLAVOR_TEX[nu_f]
    if nubar:
        # The bar belongs over the nu alone, not over the subscript too:
        # \bar{\nu}_\mu, not \bar{\nu_\mu}.
        ini = ini.replace(r'\nu', r'\bar{\nu}', 1)
        fin = fin.replace(r'\nu', r'\bar{\nu}', 1)
    return r'$P_{%s \to %s}$' % (ini, fin)


def _as_curve_list(curves):
    r"""Normalize the ``curves`` argument into a list of ``(y, plot_kwargs)``.

    .. versionadded:: 1.0.0
    """
    out = []
    for i, c in enumerate(curves):
        if isinstance(c, dict):
            d = dict(c)
            if 'y' not in d:
                raise ValueError(
                    f"Error in magnus: plotting: curve {i} is a dict without a 'y' entry; each "
                    "curve must provide its ordinate as 'y', with any "
                    'remaining keys passed through to Axes.plot'
                )
            y = d.pop('y')
        else:
            y, d = c, {}
        out.append((np.asarray(y), d))
    return out


def _apply_locators(axis, major, minor):
    r"""Set major/minor :class:`~matplotlib.ticker.MultipleLocator` spacings.

    .. versionadded:: 1.0.0
    """
    mpl, _ = _mpl()
    if major is not None:
        axis.set_major_locator(mpl.ticker.MultipleLocator(base=major))
    if minor is not None:
        axis.set_minor_locator(mpl.ticker.MultipleLocator(base=minor))


def _finish(fig, savefig, savefig_kw, tight):
    r"""Apply the shared tail of every plotting function: layout and saving.

    .. versionadded:: 1.0.0
    """
    _, plt = _mpl()
    if tight:
        fig.tight_layout()
    if savefig is not None:
        kw = dict(HOUSE_SAVEFIG_KW)
        kw.update(savefig_kw or {})
        fig.savefig(savefig, **kw)
    return fig


def plot_curves(
    x: Sequence[float],
    curves: Sequence[Union[Sequence[float], Dict[str, Any]]],
    *,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    xscale: str = 'linear',
    yscale: str = 'linear',
    xmajor: Optional[float] = None,
    xminor: Optional[float] = None,
    ymajor: Optional[float] = None,
    yminor: Optional[float] = None,
    residual: Optional[Sequence[float]] = None,
    residual_label: Optional[str] = None,
    residual_ylim: Optional[Tuple[float, float]] = None,
    residual_ymajor: Optional[float] = None,
    residual_yminor: Optional[float] = None,
    residual_height: float = HOUSE_RESIDUAL_HEIGHT,
    residual_kw: Optional[Dict[str, Any]] = None,
    annotations: Optional[Sequence[Dict[str, Any]]] = None,
    legend: bool = True,
    legend_title: Optional[str] = None,
    legend_loc: Optional[str] = None,
    legend_kw: Optional[Dict[str, Any]] = None,
    grid: bool = False,
    grid_kw: Optional[Dict[str, Any]] = None,
    ylabel_labelpad: float = 25.0,
    title_fontsize: float = 20.0,
    figsize: Tuple[float, float] = HOUSE_FIGSIZE,
    subplots_kw: Optional[Dict[str, Any]] = None,
    savefig: Optional[str] = None,
    savefig_kw: Optional[Dict[str, Any]] = None,
    tight_layout: bool = True,
):
    r"""Plot a set of curves against a swept variable, with an optional error subpanel.

    This is the workhorse: most notebook figures are an instance of it. The
    ``plot_probability_vs_baseline`` and ``plot_probability_vs_energy``
    wrappers differ from it only in their preset labels and limits, and the
    convergence studies of the matrix-exponential notebook use it directly with
    a slab count or grid size on the abscissa.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    x : sequence of float
        Abscissa, shared by every curve and by the residual panel.
    curves : sequence
        One entry per curve. An entry is either a bare ordinate array, or a
        dict carrying the ordinate under ``'y'`` plus any
        :class:`~matplotlib.lines.Line2D` keyword (``label``, ``color``,
        ``ls``, ``lw``, ...). Entries without an explicit color take the
        ``'C0'``, ``'C1'``, ... cycle in order.
    xlabel, ylabel, title : str, optional
        Axis labels and title. ``ylabel`` goes on the main panel.
    xlim, ylim : tuple of float, optional
        Axis limits. ``xlim`` is applied to the residual panel too, so the two
        panels stay aligned.
    xscale, yscale : str, optional
        Matplotlib axis scales, e.g. ``'log'``. Default is ``'linear'``.
        ``xscale`` is applied to the residual panel as well.
    xmajor, xminor, ymajor, yminor : float, optional
        Major/minor tick spacings for the main panel.
    residual : sequence of float, optional
        If given, a short subpanel is added below the main panel and this is
        plotted in it -- typically a relative error against a reference curve.
        The main panel's tick labels are then suppressed, as in the notebooks.
    residual_label : str, optional
        Ordinate label for the residual subpanel.
    residual_ylim : tuple of float, optional
        Ordinate limits for the residual subpanel.
    residual_ymajor, residual_yminor : float, optional
        Tick spacings for the residual subpanel.
    residual_height : float, optional
        Height of the residual subpanel relative to the main panel. Default is
        :data:`HOUSE_RESIDUAL_HEIGHT`.
    residual_kw : dict, optional
        Extra :class:`~matplotlib.lines.Line2D` keywords for the residual
        curve. Defaults to a thin black solid line.
    annotations : sequence of dict, optional
        Text placed on the main panel. Each entry needs ``'text'`` and
        ``'xy'`` (axes fractions by default) and may carry any other
        :meth:`~matplotlib.axes.Axes.annotate` keyword. Used by the BSM
        notebooks to record the parameter values a figure was made with.
    legend : bool, optional
        Whether to draw a legend. Default is ``True``; it is drawn only if at
        least one curve carries a ``label``.
    legend_title : str, optional
        Legend title.
    legend_loc : str, optional
        Legend location.
    legend_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_LEGEND_KW` and forwarded to
        :meth:`~matplotlib.axes.Axes.legend`.
    grid : bool, optional
        Whether to draw a grid. Default is ``False``.
    grid_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_GRID_KW`.
    ylabel_labelpad : float, optional
        Padding of the main ordinate label. Default is ``25.0``.
    title_fontsize : float, optional
        Title font size. Default is ``20.0``.
    figsize : tuple of float, optional
        Figure size in inches. Default is :data:`HOUSE_FIGSIZE`.
    subplots_kw : dict, optional
        Extra keywords for :func:`~matplotlib.pyplot.subplots`.
    savefig : str, optional
        If given, the figure is written here.
    savefig_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_SAVEFIG_KW`.
    tight_layout : bool, optional
        Whether to call :meth:`~matplotlib.figure.Figure.tight_layout`.
        Default is ``True``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure, ready for further customization or saving.
    ax : matplotlib.axes.Axes or numpy.ndarray of Axes
        A single axes when there is no residual panel; an array of two
        (main, residual) when there is.

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        from magnus.plotting import plot_curves

        # starts away from zero: the reference appears in a denominator below,
        # and sin(0)**2 is exactly 0
        L = np.linspace(50.0, 1000.0, 200)
        exact = np.sin(L / 200.0) ** 2
        approx = exact + 1e-3 * np.cos(L / 50.0)

        fig, ax = plot_curves(
            L,
            [dict(y=approx, label='Magnus expansion', color='C1'),
             dict(y=exact, label='Standard formula', color='k', ls='--')],
            xlabel=r'Baseline, $L$ [km]', ylabel='Probability',
            ylim=(0, 1), residual=(approx - exact) / exact,
            residual_label=r'$\epsilon_{\rm rel}$', legend_title='Method',
        )
        print(len(ax), ax[0].get_ylim())
    """
    _, plt = _mpl()
    entries = _as_curve_list(curves)

    has_res = residual is not None
    skw = dict(subplots_kw or {})
    if has_res:
        gs_kw = skw.pop('gridspec_kw', None) or dict(
            height_ratios=[1.0, residual_height], width_ratios=[1.0])
        fig, ax = plt.subplots(ncols=1, nrows=2, gridspec_kw=gs_kw,
                               figsize=figsize, **skw)
        fig.subplots_adjust(hspace=_HSPACE, wspace=_WSPACE)
        main, res = ax[0], ax[1]
    else:
        fig, main = plt.subplots(ncols=1, nrows=1, figsize=figsize, **skw)
        ax, res = main, None

    for i, (y, kw) in enumerate(entries):
        kw.setdefault('lw', 1)
        kw.setdefault('color', f'C{i}')
        main.plot(x, y, **kw)

    if has_res:
        rkw = dict(lw=1, color='k', ls='-')
        rkw.update(residual_kw or {})
        res.plot(x, residual, **rkw)

    for a in (annotations or []):
        a = dict(a)
        text, xy = a.pop('text'), a.pop('xy')
        a.setdefault('xycoords', 'axes fraction')
        a.setdefault('ha', 'left')
        main.annotate(text, xy=xy, **a)

    if legend and any('label' in kw for _, kw in entries):
        lkw = dict(HOUSE_LEGEND_KW)
        if legend_title is not None:
            lkw['title'] = legend_title
        if legend_loc is not None:
            lkw['loc'] = legend_loc
        lkw.update(legend_kw or {})
        main.legend(**lkw)

    if ylabel is not None:
        main.set_ylabel(ylabel, labelpad=ylabel_labelpad)
    if title is not None:
        main.set_title(title, fontsize=title_fontsize, pad=10)
    main.set_xscale(xscale)
    main.set_yscale(yscale)
    if xlim is not None:
        main.set_xlim(*xlim)
    if ylim is not None:
        main.set_ylim(*ylim)
    _apply_locators(main.xaxis, xmajor, xminor)
    _apply_locators(main.yaxis, ymajor, yminor)

    bottom = main
    if has_res:
        main.xaxis.set_ticklabels([])
        res.set_xscale(xscale)
        if xlim is not None:
            res.set_xlim(*xlim)
        if residual_ylim is not None:
            res.set_ylim(*residual_ylim)
        if residual_label is not None:
            res.set_ylabel(residual_label, labelpad=7)
        _apply_locators(res.xaxis, xmajor, xminor)
        _apply_locators(res.yaxis, residual_ymajor, residual_yminor)
        bottom = res
    if xlabel is not None:
        bottom.set_xlabel(xlabel)

    if grid:
        gkw = dict(HOUSE_GRID_KW)
        gkw.update(grid_kw or {})
        for axx in ([main, res] if has_res else [main]):
            axx.grid(**gkw)

    _finish(fig, savefig, savefig_kw, tight_layout)
    return fig, ax


def plot_curves_stacked(
    x: Sequence[float],
    panels: Sequence[Sequence[Union[Sequence[float], Dict[str, Any]]]],
    *,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    xscale: str = 'linear',
    yscale: str = 'linear',
    xmajor: Optional[float] = None,
    xminor: Optional[float] = None,
    ymajor: Optional[float] = None,
    yminor: Optional[float] = None,
    panel_labels: Optional[Sequence[str]] = None,
    panel_label_xy: Tuple[float, float] = (0.02, 0.10),
    panel_label_kw: Optional[Dict[str, Any]] = None,
    annotations: Optional[Sequence[Dict[str, Any]]] = None,
    legend: bool = True,
    legend_panel: int = 0,
    legend_proxies: Optional[Sequence[Dict[str, Any]]] = None,
    legend_title: Optional[str] = None,
    legend_loc: Optional[str] = None,
    legend_kw: Optional[Dict[str, Any]] = None,
    grid: bool = False,
    grid_kw: Optional[Dict[str, Any]] = None,
    ylabel_kw: Optional[Dict[str, Any]] = None,
    title_fontsize: float = 23.0,
    figsize: Optional[Tuple[float, float]] = None,
    height_ratios: Optional[Sequence[float]] = None,
    subplots_kw: Optional[Dict[str, Any]] = None,
    savefig: Optional[str] = None,
    savefig_kw: Optional[Dict[str, Any]] = None,
    tight_layout: bool = True,
):
    r"""Plot small multiples: one panel per case, stacked over a shared abscissa.

    The layout for "the same quantity, once per configuration" -- one panel per
    detector, per baseline, per zenith angle -- where the comparison the reader
    makes is *between* panels, so every panel must share limits, scales and tick
    spacings exactly. Only the bottom panel keeps its tick labels and abscissa
    label, and the ordinate label is a single figure-level label spanning the
    stack.

    This differs from :func:`plot_probability_with_profile`, whose panels show
    *different* quantities (a density profile above a probability), and from
    :func:`plot_curves`, whose optional second panel is a relative error rather
    than another instance of the same plot.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    x : sequence of float
        Abscissa, shared by every panel.
    panels : sequence of sequence
        One entry per panel, each a sequence of curves in the form
        :func:`plot_curves` accepts: a bare ordinate array, or a dict carrying
        the ordinate under ``'y'`` plus any
        :class:`~matplotlib.lines.Line2D` keyword. Curves without an explicit
        color take the ``'C0'``, ``'C1'``, ... cycle *within* their panel, so
        the n-th curve of every panel matches by default.
    xlabel : str, optional
        Abscissa label, placed on the bottom panel only.
    ylabel : str, optional
        Ordinate label. Drawn once for the whole stack with
        :meth:`~matplotlib.figure.Figure.supylabel`, since every panel shows
        the same quantity. Being figure-level, it takes no ``labelpad``; use
        ``ylabel_kw`` for its placement.
    title : str, optional
        Title, placed above the top panel.
    xlim, ylim : tuple of float, optional
        Axis limits, applied to every panel.
    xscale, yscale : str, optional
        Matplotlib axis scales, applied to every panel. Default ``'linear'``.
    xmajor, xminor, ymajor, yminor : float, optional
        Major/minor tick spacings, applied to every panel.
    panel_labels : sequence of str, optional
        One caption per panel, annotated inside it -- the usual way of saying
        which case a panel is. Must match the number of panels.
    panel_label_xy : tuple of float, optional
        Position of those captions, in axes fractions. Default ``(0.02, 0.10)``.
    panel_label_kw : dict, optional
        Extra :meth:`~matplotlib.axes.Axes.annotate` keywords for them.
    annotations : sequence of dict, optional
        Free-form text. Each entry needs ``'text'`` and ``'xy'``, may name a
        ``'panel'`` (index, default 0), and may carry any other
        :meth:`~matplotlib.axes.Axes.annotate` keyword.
    legend : bool, optional
        Whether to draw a legend. Default ``True``; drawn only if there is
        something to put in it.
    legend_panel : int, optional
        Which panel carries the legend. Default ``0``.
    legend_proxies : sequence of dict, optional
        Legend entries that describe a *style* shared across panels rather than
        any one curve -- e.g. "solid: 3+1, dashed: standard" when the color
        varies from panel to panel. Each entry is a set of
        :class:`~matplotlib.lines.Line2D` keywords including ``label``, drawn
        as an empty proxy artist. When given, these replace the labels picked
        up from the curves themselves. This exists because the alternative, and
        what the notebooks did, is plotting dummy points outside the axis
        limits to manufacture legend handles.
    legend_title, legend_loc : str, optional
        Legend title and location.
    legend_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_LEGEND_KW`.
    grid : bool, optional
        Whether to draw a grid on every panel. Default ``False``.
    grid_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_GRID_KW`.
    ylabel_kw : dict, optional
        Extra keywords for :meth:`~matplotlib.figure.Figure.supylabel`.
    title_fontsize : float, optional
        Title font size. Default ``23.0``.
    figsize : tuple of float, optional
        Figure size in inches. Defaults to :data:`HOUSE_FIGSIZE`'s width and
        half its height per panel, which reproduces the notebooks' proportions.
    height_ratios : sequence of float, optional
        Relative panel heights. Default: equal.
    subplots_kw : dict, optional
        Extra keywords for :func:`~matplotlib.pyplot.subplots`.
    savefig : str, optional
        If given, the figure is written here.
    savefig_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_SAVEFIG_KW`.
    tight_layout : bool, optional
        Whether to call :meth:`~matplotlib.figure.Figure.tight_layout`.
        Default ``True``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure, ready for further customization or saving.
    ax : numpy.ndarray of Axes
        One axes per panel, top to bottom. Always an array, including for a
        single panel, so that indexing does not depend on the panel count.

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        from magnus.plotting import plot_curves_stacked

        E = np.linspace(1.0, 40.0, 200)
        cases = [0.5, 1.0, 2.0]
        panels = [
            [dict(y=np.sin(k*E/8.0)**2, color=f'C{i}'),
             dict(y=np.sin(k*E/8.0)**2*0.8, color='0.7', ls='--')]
            for i, k in enumerate(cases)
        ]

        fig, ax = plot_curves_stacked(
            E, panels,
            xlabel=r'Neutrino energy, $E_\nu$ [GeV]', ylabel='Probability',
            ylim=(0, 1), xlim=(1.0, 40.0),
            panel_labels=[f'baseline {k:.1f} kton-yr' for k in cases],
            legend_proxies=[dict(label='3+1', color='k', ls='-'),
                            dict(label=r'standard', color='k', ls='--')],
        )
        print(ax.shape, ax[0].get_xticklabels()[0].get_text() == '')
    """
    mpl, plt = _mpl()

    n = len(panels)
    if n == 0:
        raise ValueError(
            'Error in magnus: plotting.plot_curves_stacked: panels is empty; it needs at least '
            'one panel, each a sequence of curves'
        )
    if panel_labels is not None and len(panel_labels) != n:
        raise ValueError(
            f'Error in magnus: plotting.plot_curves_stacked: got {len(panel_labels)} panel_labels '
            f'for {n} panels; there must be exactly one label per panel'
        )
    if not (-n <= legend_panel < n):
        raise ValueError(
            f'Error in magnus: plotting.plot_curves_stacked: legend_panel={legend_panel} is out of '
            f'range for {n} panels'
        )

    if figsize is None:
        figsize = (HOUSE_FIGSIZE[0], 0.5*HOUSE_FIGSIZE[1]*n)
    skw = dict(subplots_kw or {})
    gs_kw = skw.pop('gridspec_kw', None) or dict(
        height_ratios=list(height_ratios) if height_ratios is not None else [1.0]*n,
        width_ratios=[1.0])
    # squeeze=False so a one-panel stack still indexes like every other one.
    fig, ax = plt.subplots(ncols=1, nrows=n, gridspec_kw=gs_kw, figsize=figsize,
                           squeeze=False, **skw)
    ax = ax[:, 0]
    fig.subplots_adjust(hspace=_HSPACE, wspace=_WSPACE)

    for panel, axx in zip(panels, ax):
        for i, (y, kw) in enumerate(_as_curve_list(panel)):
            kw.setdefault('lw', 1)
            kw.setdefault('color', f'C{i}')
            axx.plot(x, y, **kw)

    for i, axx in enumerate(ax):
        axx.set_xscale(xscale)
        axx.set_yscale(yscale)
        if xlim is not None:
            axx.set_xlim(*xlim)
        if ylim is not None:
            axx.set_ylim(*ylim)
        _apply_locators(axx.xaxis, xmajor, xminor)
        _apply_locators(axx.yaxis, ymajor, yminor)
        if grid:
            gkw = dict(HOUSE_GRID_KW)
            gkw.update(grid_kw or {})
            axx.grid(**gkw)
        if i != n - 1:
            axx.xaxis.set_ticklabels([])

    for axx, text in zip(ax, panel_labels or []):
        akw = dict(xycoords='axes fraction', ha='left', fontsize=20)
        akw.update(panel_label_kw or {})
        axx.annotate(text, xy=panel_label_xy, **akw)

    for a in (annotations or []):
        a = dict(a)
        text, xy = a.pop('text'), a.pop('xy')
        panel = a.pop('panel', 0)
        a.setdefault('xycoords', 'axes fraction')
        a.setdefault('ha', 'left')
        ax[panel].annotate(text, xy=xy, **a)

    if legend:
        lkw = dict(HOUSE_LEGEND_KW)
        if legend_title is not None:
            lkw['title'] = legend_title
        if legend_loc is not None:
            lkw['loc'] = legend_loc
        lkw.update(legend_kw or {})
        if legend_proxies:
            handles = [mpl.lines.Line2D([], [], **dict(kw)) for kw in legend_proxies]
            ax[legend_panel].legend(handles=handles, **lkw)
        elif any('label' in kw for _, kw in _as_curve_list(panels[legend_panel])):
            ax[legend_panel].legend(**lkw)

    if title is not None:
        ax[0].set_title(title, fontsize=title_fontsize, pad=10)
    if xlabel is not None:
        ax[-1].set_xlabel(xlabel)
    if ylabel is not None:
        # Match the abscissa label rather than Matplotlib's figure-label default:
        # supylabel takes its size from rcParams['figure.labelsize'] ('large'),
        # while every axis label in these figures takes rcParams['axes.labelsize']
        # (the notebooks' matplotlibrc sets it to 25). Left alone, the shared
        # ordinate label comes out visibly smaller than the abscissa label beneath
        # it, which is not what the hand-built version did.
        ykw = {'fontsize': plt.rcParams['axes.labelsize']}
        ykw.update(ylabel_kw or {})
        fig.supylabel(ylabel, **ykw)

    _finish(fig, savefig, savefig_kw, tight_layout)
    return fig, ax


def plot_probability_vs_baseline(
    distances: Sequence[float],
    curves: Sequence[Union[Sequence[float], Dict[str, Any]]],
    *,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    num_flavors: Optional[int] = None,
    xlabel: str = r'Baseline, $L$ [km]',
    ylabel: Optional[str] = None,
    ylim: Tuple[float, float] = (0.0, 1.0),
    xscale: str = 'log',
    ymajor: Optional[float] = 0.10,
    yminor: Optional[float] = 0.02,
    **_forbidden: Any,
):
    r"""Plot oscillation probabilities against baseline.

    A thin preset over :func:`plot_curves`: log abscissa, ordinate on
    :math:`[0, 1]` with the notebooks' tick spacings, and an ordinate label
    built from the flavor pair.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    distances : sequence of float
        Baselines [km].
    curves : sequence
        As in :func:`plot_curves`.
    nu_i, nu_f : int, optional
        Flavor pair, used to build the ordinate label via :func:`prob_label`
        when ``ylabel`` is not given.
    num_flavors : int, optional
        If given, prefixes the ordinate label with ``'Two-'``, ``'Three-'``,
        ``'Four-'`` or ``'Five-neutrino probability'``.
    xlabel : str, optional
        Abscissa label.
    ylabel : str, optional
        Ordinate label; overrides the one built from the flavor pair.
    ylim : tuple of float, optional
        Ordinate limits. Default is ``(0.0, 1.0)``.
    xscale : str, optional
        Abscissa scale. Default is ``'log'``.
    ymajor, yminor : float, optional
        Ordinate tick spacings.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes or numpy.ndarray of Axes

    Other Parameters
    ----------------
    **_forbidden
        Every remaining keyword of :func:`plot_curves` is accepted and
        forwarded unchanged; unknown names raise :class:`TypeError` there.

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        import magnus.globaldefs as gd
        from magnus.plotting import plot_probability_vs_baseline

        L = np.logspace(1, 5, 200)
        P = np.sin(L / 3000.0) ** 2
        fig, ax = plot_probability_vs_baseline(
            L, [dict(y=P, label='Magnus expansion')],
            nu_i=gd.NUE, nu_f=gd.NUE, num_flavors=2, xlim=(L[0], L[-1]),
        )
        print(ax.get_xlabel())
    """
    if ylabel is None and nu_i is not None and nu_f is not None:
        ylabel = _probability_ylabel(nu_i, nu_f, num_flavors)
    return plot_curves(
        distances, curves, xlabel=xlabel, ylabel=ylabel, ylim=ylim,
        xscale=xscale, ymajor=ymajor, yminor=yminor, **_forbidden)


def plot_probability_vs_energy(
    energies: Sequence[float],
    curves: Sequence[Union[Sequence[float], Dict[str, Any]]],
    *,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    num_flavors: Optional[int] = None,
    energy_unit: str = 'GeV',
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    ylim: Tuple[float, float] = (0.0, 1.0),
    xscale: str = 'log',
    ymajor: Optional[float] = 0.10,
    yminor: Optional[float] = 0.02,
    **_forbidden: Any,
):
    r"""Plot oscillation probabilities against neutrino energy.

    The energy counterpart of :func:`plot_probability_vs_baseline`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    energies : sequence of float
        Neutrino energies, in the unit named by ``energy_unit``.
    curves : sequence
        As in :func:`plot_curves`.
    nu_i, nu_f : int, optional
        Flavor pair for the ordinate label.
    num_flavors : int, optional
        Flavor count, for the ordinate label prefix.
    energy_unit : str, optional
        Unit shown in the abscissa label. Default is ``'GeV'``.
    xlabel : str, optional
        Abscissa label; overrides the one built from ``energy_unit``.
    ylabel : str, optional
        Ordinate label.
    ylim : tuple of float, optional
        Ordinate limits. Default is ``(0.0, 1.0)``.
    xscale : str, optional
        Abscissa scale. Default is ``'log'``.
    ymajor, yminor : float, optional
        Ordinate tick spacings.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes or numpy.ndarray of Axes

    Other Parameters
    ----------------
    **_forbidden
        Forwarded to :func:`plot_curves`.

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        import magnus.globaldefs as gd
        from magnus.plotting import plot_probability_vs_energy

        E = np.logspace(-1, 1, 200)
        P = np.cos(1.0 / E) ** 2
        fig, ax = plot_probability_vs_energy(
            E, [dict(y=P, label='Magnus expansion')],
            nu_i=gd.NUMU, nu_f=gd.NUE, xlim=(E[0], E[-1]),
        )
        print(ax.get_xlabel())
    """
    if xlabel is None:
        xlabel = r'Neutrino energy, $E_\nu$ [%s]' % energy_unit
    if ylabel is None and nu_i is not None and nu_f is not None:
        ylabel = _probability_ylabel(nu_i, nu_f, num_flavors)
    return plot_curves(
        energies, curves, xlabel=xlabel, ylabel=ylabel, ylim=ylim,
        xscale=xscale, ymajor=ymajor, yminor=yminor, **_forbidden)


_FLAVOR_WORD = {2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five'}


def _probability_ylabel(nu_i, nu_f, num_flavors):
    r"""Build the notebooks' ordinate label for a probability panel.

    .. versionadded:: 1.0.0
    """
    label = prob_label(nu_i, nu_f)
    if num_flavors is None:
        return 'Probability,~' + label
    if num_flavors not in _FLAVOR_WORD:
        raise ValueError(
            'Error in magnus: plotting: num_flavors must be one of '
            f'{sorted(_FLAVOR_WORD)}, not {num_flavors!r}'
        )
    return f'{_FLAVOR_WORD[num_flavors]}-neutrino probability,~' + label


def plot_probability_with_profile(
    x: Sequence[float],
    profiles: Optional[Sequence[Union[Sequence[float], Dict[str, Any]]]],
    panels: Sequence[Sequence[Union[Sequence[float], Dict[str, Any]]]],
    *,
    xlabel: str = r'Baseline, $L$~[km]',
    profile_ylabel: str = r'$\frac{N_e}{N_{\rm Av}}$~[cm$^{-3}$]',
    panel_ylabels: Optional[Sequence[Optional[str]]] = None,
    panel_annotations: Optional[Sequence[Optional[str]]] = None,
    panel_annotation_xy: Tuple[float, float] = (0.02, 0.88),
    panel_annotation_fontsize: float = 23.0,
    shared_ylabel: Optional[str] = None,
    shared_ylabel_labelpad: float = 20.0,
    title: Optional[str] = None,
    title_fontsize: float = 23.0,
    xlim: Optional[Tuple[float, float]] = None,
    xscale: str = 'log',
    xmajor: Optional[float] = None,
    xminor: Optional[float] = None,
    profile_ylim: Optional[Tuple[float, float]] = None,
    profile_ymajor: Optional[float] = None,
    profile_yminor: Optional[float] = None,
    profile_height: float = 0.4,
    panel_ylim: Optional[Tuple[float, float]] = (0.0, 1.0),
    panel_yscale: str = 'linear',
    panel_ymajor: Optional[float] = 0.10,
    panel_yminor: Optional[float] = 0.02,
    legend: bool = True,
    legend_title: Optional[str] = None,
    legend_loc: Optional[str] = None,
    legend_kw: Optional[Dict[str, Any]] = None,
    legend_on_panel: int = 0,
    grid: bool = True,
    grid_kw: Optional[Dict[str, Any]] = None,
    ylabel_labelpad: float = 25.0,
    figsize: Optional[Tuple[float, float]] = None,
    subplots_kw: Optional[Dict[str, Any]] = None,
    savefig: Optional[str] = None,
    savefig_kw: Optional[Dict[str, Any]] = None,
    tight_layout: bool = False,
):
    r"""Stack a matter-density panel above one or more probability panels.

    This is the layout of the long-baseline notebook: the electron-density
    profile along the trajectory on top, then one probability panel per
    detector or per profile, sharing the abscissa. With a single probability
    panel it is the profile-plus-probability figure of the introduction and the
    two-flavor notebooks.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    x : sequence of float
        Shared abscissa, or, when the panels have different abscissae, the one
        used by the profile panel. Individual curves may carry their own ``x``.
    profiles : sequence or None
        Curves for the density panel, in the form :func:`plot_curves` takes.
        A curve may add its own abscissa under ``'x'``. Pass ``None`` (or an
        empty sequence) to omit the density panel entirely and get a plain
        stack of probability panels sharing an abscissa -- the layout the
        long-baseline notebook uses for a probability above its
        energy-smoothed version.
    panels : sequence of sequence
        One entry per probability panel; each entry is a sequence of curves.
    xlabel : str, optional
        Abscissa label, placed under the bottom panel.
    profile_ylabel : str, optional
        Ordinate label of the density panel.
    panel_ylabels : sequence of str, optional
        Ordinate labels for the probability panels. Entries may be ``None``.
    panel_annotations : sequence, optional
        Text placed inside each probability panel, one entry per panel, at
        ``panel_annotation_xy`` in axes coordinates. An entry is a string, or a
        dict with ``'text'`` plus any other
        :meth:`~matplotlib.axes.Axes.annotate` keyword -- a ``bbox``, for
        instance, when the text would otherwise sit over dense curves. Entries
        may be ``None``. The three-flavor notebook uses this to name the
        channel each panel shows, rather than repeating it in the ordinate
        label.
    panel_annotation_xy : tuple of float, optional
        Position of those annotations, in axes fractions. Default
        ``(0.02, 0.88)``.
    panel_annotation_fontsize : float, optional
        Their font size. Default is ``23.0``.
    shared_ylabel : str, optional
        A single ordinate label spanning the whole stack, drawn on a frameless
        overlay axes. Use it instead of ``panel_ylabels`` when every panel
        shows the same quantity.
    shared_ylabel_labelpad : float, optional
        Padding of that shared label. Default is ``20.0``.
    title : str, optional
        Title, placed above the density panel.
    title_fontsize : float, optional
        Title font size. Default is ``23.0``.
    xlim : tuple of float, optional
        Shared abscissa limits.
    xscale : str, optional
        Shared abscissa scale. Default is ``'log'``.
    xmajor, xminor : float, optional
        Major/minor tick spacings on the shared abscissa. Only meaningful on a
        linear scale.
    profile_ylim : tuple of float, optional
        Ordinate limits of the density panel.
    profile_ymajor, profile_yminor : float, optional
        Tick spacings for the density panel.
    profile_height : float, optional
        Height of the density panel relative to a probability panel. Default
        is ``0.4``.
    panel_ylim : tuple of float, optional
        Ordinate limits shared by the probability panels. ``None`` autoscales.
    panel_yscale : str, optional
        Ordinate scale shared by the panels, e.g. ``'log'`` when they carry
        something other than a probability. Default is ``'linear'``.
    panel_ymajor, panel_yminor : float, optional
        Tick spacings for the probability panels.
    legend : bool, optional
        Whether to draw a legend.
    legend_title : str, optional
        Legend title.
    legend_loc : str, optional
        Legend location.
    legend_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_LEGEND_KW`.
    legend_on_panel : int, optional
        Index of the probability panel carrying the legend, or ``-1`` to give
        every panel its own. Default is ``0``.
    grid : bool, optional
        Whether to draw grids. Default is ``True``.
    grid_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_GRID_KW`.
    ylabel_labelpad : float, optional
        Padding of the ordinate labels.
    figsize : tuple of float, optional
        Figure size. Defaults to ``(18, 9)`` for one probability panel, growing
        by 4.5 inches per extra panel.
    subplots_kw : dict, optional
        Extra keywords for :func:`~matplotlib.pyplot.subplots`.
    savefig : str, optional
        If given, the figure is written here.
    savefig_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_SAVEFIG_KW`.
    tight_layout : bool, optional
        Whether to call ``tight_layout``. Default is ``False``, matching the
        notebooks, whose explicit ``subplots_adjust`` this would override.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : numpy.ndarray of Axes
        Length ``1 + len(panels)`` with a density panel, which comes first;
        length ``len(panels)`` without one.

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        from magnus.plotting import plot_probability_with_profile

        L = np.logspace(2, 4, 300)
        n_e = 5.0 * np.exp(-L / 5000.0)
        P = np.sin(L / 900.0) ** 2

        fig, ax = plot_probability_with_profile(
            L, [dict(y=n_e, color='C0')], [[dict(y=P, label='PREM')]],
            xlim=(L[0], L[-1]), profile_ylim=(0, 6),
        )
        print(len(ax))
    """
    _, plt = _mpl()
    n_panels = len(panels)
    if n_panels == 0:
        raise ValueError('Error in magnus: plotting.plot_probability_with_profile: panels is '
                         'empty; at least one probability panel is required')
    has_profile = bool(profiles)
    n_rows = n_panels + (1 if has_profile else 0)
    if figsize is None:
        figsize = (HOUSE_FIGSIZE[0], HOUSE_FIGSIZE[1] + 4.5 * (n_panels - 1))

    skw = dict(subplots_kw or {})
    ratios = ([profile_height] if has_profile else []) + [1.0] * n_panels
    gs_kw = skw.pop('gridspec_kw', None) or dict(height_ratios=ratios,
                                                 width_ratios=[1.0])
    fig, ax = plt.subplots(ncols=1, nrows=n_rows, gridspec_kw=gs_kw,
                           figsize=figsize, squeeze=False, **skw)
    ax = ax[:, 0]
    fig.subplots_adjust(hspace=0.1, wspace=0.1)

    def _draw(axx, entries):
        for i, (y, kw) in enumerate(_as_curve_list(entries)):
            kw.setdefault('lw', 1)
            kw.setdefault('color', f'C{i}')
            xi = kw.pop('x', x)
            axx.plot(xi, y, **kw)

    off = 1 if has_profile else 0
    if has_profile:
        _draw(ax[0], profiles)
    for j, panel in enumerate(panels):
        _draw(ax[j + off], panel)

    if legend:
        lkw = dict(HOUSE_LEGEND_KW)
        if legend_title is not None:
            lkw['title'] = legend_title
        if legend_loc is not None:
            lkw['loc'] = legend_loc
        lkw.update(legend_kw or {})
        targets = range(n_panels) if legend_on_panel == -1 else [legend_on_panel]
        for j in targets:
            if ax[j + off].get_legend_handles_labels()[1]:
                ax[j + off].legend(**lkw)

    for i, axx in enumerate(ax):
        axx.set_xscale(xscale)
        if xlim is not None:
            axx.set_xlim(*xlim)
        _apply_locators(axx.xaxis, xmajor, xminor)
        if grid:
            gkw = dict(HOUSE_GRID_KW)
            gkw.update(grid_kw or {})
            axx.grid(**gkw)
        if i != len(ax) - 1:
            axx.xaxis.set_ticklabels([])

    if has_profile:
        ax[0].set_ylabel(profile_ylabel, labelpad=ylabel_labelpad)
        if profile_ylim is not None:
            ax[0].set_ylim(*profile_ylim)
        _apply_locators(ax[0].yaxis, profile_ymajor, profile_yminor)
    if title is not None:
        ax[0].set_title(title, fontsize=title_fontsize, pad=10)

    for j in range(n_panels):
        axx = ax[j + off]
        axx.set_yscale(panel_yscale)
        if panel_ylim is not None:
            axx.set_ylim(*panel_ylim)
        _apply_locators(axx.yaxis, panel_ymajor, panel_yminor)
        if panel_ylabels is not None and j < len(panel_ylabels):
            if panel_ylabels[j] is not None:
                axx.set_ylabel(panel_ylabels[j], labelpad=15)
        if panel_annotations is not None and j < len(panel_annotations):
            entry = panel_annotations[j]
            if entry is not None:
                # a bare string, or a dict carrying extra annotate keywords
                # (a white bbox, say, so the text stays readable over curves)
                akw = dict(xy=panel_annotation_xy, xycoords='axes fraction',
                           ha='left', fontsize=panel_annotation_fontsize)
                if isinstance(entry, dict):
                    entry = dict(entry)
                    text = entry.pop('text')
                    akw.update(entry)
                else:
                    text = entry
                axx.annotate(text, **akw)

    ax[-1].set_xlabel(xlabel)

    if shared_ylabel is not None:
        # A frameless axes spanning the figure carries one label for the whole
        # stack. Its tick *marks* are switched off but its tick *labels* are
        # only made invisible, not removed: they still reserve the width that
        # pushes the shared label clear of the panels' own tick labels. Calling
        # set_yticks([]) here instead would drop that reservation and the label
        # would land on top of the numbers.
        overlay = fig.add_subplot(111, frameon=False)
        overlay.tick_params(labelcolor='none', top=False, bottom=False,
                            left=False, right=False)
        overlay.set_ylabel(shared_ylabel, labelpad=shared_ylabel_labelpad)

    _finish(fig, savefig, savefig_kw, tight_layout)
    return fig, ax


def plot_probability_with_average(
    x: Sequence[float],
    probabilities: Union[Sequence[float], Sequence[Sequence[float]]],
    averages: Union[float, Sequence[float], Sequence[Sequence[float]]],
    *,
    labels: Optional[Sequence[str]] = None,
    colors: Optional[Sequence[str]] = None,
    average_label: str = 'Phase-averaged',
    oscillating_kw: Optional[Dict[str, Any]] = None,
    average_kw: Optional[Dict[str, Any]] = None,
    **_forbidden: Any,
):
    r"""Overlay phase-averaged probabilities on the oscillating ones.

    The figure of the averaged-probability notebook: rapidly oscillating
    curves, each with its decohered limit drawn through it as a dashed line of
    the same color -- the value :func:`magnus.oscprob.osc_prob` returns with
    ``average=True``.

    Several channels are usually shown at once, so the legend carries one
    entry per channel plus a single entry explaining the dashed style, rather
    than repeating "averaged" once per curve.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    x : sequence of float
        Abscissa, typically baseline [km].
    probabilities : sequence of float or sequence of sequence of float
        One oscillating probability, or several.
    averages : float or sequence
        The corresponding phase-averaged values: a scalar per curve (broadcast
        across ``x``), or a full curve each. Must match ``probabilities`` in
        number.
    labels : sequence of str, optional
        Legend label per channel.
    colors : sequence of str, optional
        Color per channel. Defaults to the ``'C0'``, ``'C1'``, ... cycle;
        each average takes its channel's color.
    average_label : str, optional
        Text of the single legend entry explaining the dashed lines.
    oscillating_kw, average_kw : dict, optional
        Extra :class:`~matplotlib.lines.Line2D` keywords applied to every
        oscillating or every averaged curve.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes or numpy.ndarray of Axes

    Raises
    ------
    ValueError
        If the number of averages does not match the number of probabilities.

    Other Parameters
    ----------------
    **_forbidden
        Forwarded to :func:`plot_probability_vs_baseline`.

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        from magnus.plotting import plot_probability_with_average

        L = np.linspace(1.0, 1000.0, 500)
        P = np.sin(L / 7.0) ** 2
        fig, ax = plot_probability_with_average(L, P, 0.5, xscale='linear')
        print(len(ax.get_lines()))
    """
    _, plt = _mpl()
    import matplotlib.lines as mlines

    probs = np.asarray(probabilities, dtype=float)
    if probs.ndim == 1:
        probs = probs[None, :]
    avgs = np.asarray(averages, dtype=float)
    if avgs.ndim == 0:
        avgs = avgs[None]
    if len(avgs) != len(probs):
        raise ValueError(
            'Error in magnus: plotting.plot_probability_with_average: got '
            f'{len(probs)} probability curve(s) but {len(avgs)} average(s); '
            'they must correspond one to one'
        )

    n = len(np.asarray(x))
    colors = list(colors) if colors else [f'C{i}' for i in range(len(probs))]
    curves = []
    for i, p in enumerate(probs):
        okw = dict(lw=1, color=colors[i], ls='-')
        okw.update(oscillating_kw or {})
        if labels is not None and i < len(labels):
            okw['label'] = labels[i]
        curves.append(dict(y=p, **okw))

        a = avgs[i]
        avg_curve = np.full(n, float(a)) if np.ndim(a) == 0 else np.asarray(a)
        akw = dict(lw=1.5, color=colors[i], ls='--')
        akw.update(average_kw or {})
        curves.append(dict(y=avg_curve, **akw))

    fig, ax = plot_probability_vs_baseline(x, curves, legend=False, **_forbidden)

    main = ax[0] if isinstance(ax, np.ndarray) else ax
    handles, labs = main.get_legend_handles_labels()
    style = dict(lw=1.5, color='k', ls='--')
    style.update({k: v for k, v in (average_kw or {}).items() if k != 'color'})
    handles.append(mlines.Line2D([], [], **style))
    labs.append(average_label)
    lkw = dict(HOUSE_LEGEND_KW)
    lkw.update(_forbidden.get('legend_kw') or {})
    for key, name in (('legend_title', 'title'), ('legend_loc', 'loc')):
        if _forbidden.get(key) is not None:
            lkw[name] = _forbidden[key]
    main.legend(handles, labs, **lkw)
    return fig, ax


def plot_biprobability(
    prob_nu: Sequence[Sequence[float]],
    prob_nubar: Sequence[Sequence[float]],
    *,
    labels: Optional[Sequence[str]] = None,
    curve_kw: Optional[Sequence[Dict[str, Any]]] = None,
    markers: Optional[Sequence[Dict[str, Any]]] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    title_fontsize: float = 20.0,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    xmajor: Optional[float] = None,
    xminor: Optional[float] = None,
    ymajor: Optional[float] = None,
    yminor: Optional[float] = None,
    annotations: Optional[Sequence[Dict[str, Any]]] = None,
    legend: bool = True,
    legend_title: str = r'$\delta_{\rm CP}$',
    legend_loc: Optional[str] = None,
    legend_kw: Optional[Dict[str, Any]] = None,
    figsize: Tuple[float, float] = (9.0, 9.0),
    subplots_kw: Optional[Dict[str, Any]] = None,
    savefig: Optional[str] = None,
    savefig_kw: Optional[Dict[str, Any]] = None,
    tight_layout: bool = False,
):
    r"""Plot neutrino against antineutrino appearance probability.

    The bi-probability plane: for each configuration, the locus traced out as
    :math:`\delta_{\rm CP}` runs over :math:`[-\pi, \pi]`, with optional
    markers at selected phases.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    prob_nu, prob_nubar : sequence of sequence of float
        One entry per curve, each a sequence of probabilities over the same
        grid of :math:`\delta_{\rm CP}` values.
    labels : sequence of str, optional
        Legend label per curve.
    curve_kw : sequence of dict, optional
        Per-curve :class:`~matplotlib.lines.Line2D` keywords.
    markers : sequence of dict, optional
        Markers at selected phases. Each entry gives its position either as
        ``'index'`` (a position along the curve) or as ``'xy'`` (an explicit
        coordinate pair, which is what you have when the marked phases were
        computed separately from the curve). Optionally ``'marker'``,
        ``'label'``, ``'filled'`` and ``'curve'`` (which curve it belongs to,
        default all).
    xlabel, ylabel : str, optional
        Axis labels. Default to the :math:`\nu_\mu \to \nu_e` pair.
    title : str, optional
        Title.
    title_fontsize : float, optional
        Title font size.
    xlim, ylim : tuple of float, optional
        Axis limits.
    xmajor, xminor, ymajor, yminor : float, optional
        Tick spacings.
    annotations : sequence of dict, optional
        Passed to :meth:`~matplotlib.axes.Axes.annotate`; each entry needs
        ``'text'`` and ``'xy'``, and may carry any other annotate keyword.
        Coordinates are axes fractions.
    legend : bool, optional
        Whether to draw a legend.
    legend_title : str, optional
        Legend title. Default is ``'$\\delta_{\\rm CP}$'``.
    legend_loc : str, optional
        Legend location.
    legend_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_LEGEND_KW`.
    figsize : tuple of float, optional
        Figure size. Default is ``(9.0, 9.0)``, the square panel this plot uses.
    subplots_kw : dict, optional
        Extra keywords for :func:`~matplotlib.pyplot.subplots`.
    savefig : str, optional
        If given, the figure is written here.
    savefig_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_SAVEFIG_KW`.
    tight_layout : bool, optional
        Whether to call ``tight_layout``. Default is ``False``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        from magnus.plotting import plot_biprobability

        d = np.linspace(-np.pi, np.pi, 100)
        P_nu = 0.05 + 0.02 * np.sin(d)
        P_nubar = 0.04 + 0.02 * np.sin(d + 0.4)

        fig, ax = plot_biprobability([P_nu], [P_nubar], labels=['NO'])
        print(ax.get_xlabel())
    """
    import magnus.globaldefs as gd

    _, plt = _mpl()
    if len(prob_nu) != len(prob_nubar):
        raise ValueError(
            'Error in magnus: plotting.plot_biprobability: prob_nu and prob_nubar must have '
            f'the same number of curves, got {len(prob_nu)} and '
            f'{len(prob_nubar)}'
        )
    if xlabel is None:
        xlabel = prob_label(gd.NUMU, gd.NUE)
    if ylabel is None:
        ylabel = prob_label(gd.NUMU, gd.NUE, nubar=True)

    skw = dict(subplots_kw or {})
    gs_kw = skw.pop('gridspec_kw', None) or dict(height_ratios=[1.0],
                                                 width_ratios=[1.0])
    fig, ax = plt.subplots(ncols=1, nrows=1, gridspec_kw=gs_kw,
                           figsize=figsize, **skw)
    fig.subplots_adjust(hspace=_HSPACE, wspace=_WSPACE)

    colors = []
    for i, (yn, yb) in enumerate(zip(prob_nu, prob_nubar)):
        kw = dict(lw=1, color=f'C{i}', ls='-')
        if curve_kw is not None and i < len(curve_kw):
            kw.update(curve_kw[i] or {})
        if labels is not None and i < len(labels):
            kw.setdefault('label', labels[i])
        colors.append(kw['color'])
        ax.plot(np.asarray(yn), np.asarray(yb), **kw)

    for m in (markers or []):
        if 'index' not in m and 'xy' not in m:
            raise ValueError(
                "Error in magnus: plotting.plot_biprobability: each marker needs either an "
                "'index' along the curve or an explicit 'xy' coordinate pair; "
                f'got keys {sorted(m)}'
            )
        filled = m.get('filled', True)
        which = m.get('curve')
        targets = range(len(prob_nu)) if which is None else [which]
        for i in targets:
            c = colors[i]
            if 'xy' in m:
                x, y = m['xy']
            else:
                x = np.asarray(prob_nu[i])[m['index']]
                y = np.asarray(prob_nubar[i])[m['index']]
            ax.scatter(x, y, marker=m.get('marker', 'o'), ls='-', edgecolors=c,
                       s=70, c=c if filled else 'none')

    if legend:
        # Proxy handles, drawn off-axis, so the marker legend does not depend
        # on which curve happened to be plotted last.
        labelled = [m for m in (markers or []) if m.get('label')]
        for m in labelled:
            ax.scatter(np.nan, np.nan, marker=m.get('marker', 'o'), ls='-',
                       edgecolors='k', s=70, label=m['label'],
                       c='k' if m.get('filled', True) else 'none')
        if ax.get_legend_handles_labels()[1]:
            lkw = dict(HOUSE_LEGEND_KW)
            lkw['title_fontsize'] = 18
            lkw['title'] = legend_title
            if legend_loc is not None:
                lkw['loc'] = legend_loc
            lkw.update(legend_kw or {})
            ax.legend(**lkw)

    for a in (annotations or []):
        a = dict(a)
        text, xy = a.pop('text'), a.pop('xy')
        a.setdefault('xycoords', 'axes fraction')
        ax.annotate(text, xy=xy, **a)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title, fontsize=title_fontsize, pad=10)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    _apply_locators(ax.xaxis, xmajor, xminor)
    _apply_locators(ax.yaxis, ymajor, yminor)

    _finish(fig, savefig, savefig_kw, tight_layout)
    return fig, ax


def plot_oscillogram(
    costhz: Sequence[float],
    log10_energy: Sequence[float],
    probability: Sequence[Sequence[float]],
    *,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    levels: int = 120,
    cmap: str = 'plasma',
    xlabel: str = r'Zenith angle, $\cos(\theta_z)$',
    ylabel: str = r'Neutrino energy, $\log_{10}(E_\nu/{\rm GeV})$',
    cbar_label: Optional[str] = None,
    cbar_label_prefix: str = '',
    cbar_fontsize: float = 25.0,
    cbar_labelsize: float = 25.0,
    annotation: Optional[str] = None,
    annotation_fontsize: float = 23.0,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    xmajor: Optional[float] = 0.2,
    xminor: Optional[float] = 0.02,
    ymajor: Optional[float] = 0.1,
    yminor: Optional[float] = 0.02,
    figsize: Tuple[float, float] = (9.0, 9.0),
    contourf_kw: Optional[Dict[str, Any]] = None,
    subplots_kw: Optional[Dict[str, Any]] = None,
    savefig: Optional[str] = None,
    savefig_kw: Optional[Dict[str, Any]] = None,
    tight_layout: bool = False,
):
    r"""Plot an oscillogram: probability over zenith angle and energy.

    A filled contour map of the oscillation probability in the plane of
    :math:`\cos\theta_z` (equivalently, baseline through the Earth) and
    :math:`\log_{10} E_\nu`, with a color bar and the channel annotated in the
    corner over a white stroke so it stays legible against the color map.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    costhz : sequence of float
        Zenith-angle cosines, the abscissa.
    log10_energy : sequence of float
        :math:`\log_{10}` of the energy in GeV, the ordinate.
    probability : sequence of sequence of float
        Probability with shape ``(len(log10_energy), len(costhz))``.
    nu_i, nu_f : int, optional
        Flavor pair, used for the color-bar label and the annotation when
        those are not given explicitly.
    levels : int, optional
        Number of filled contour levels. Default is ``120``.
    cmap : str, optional
        Color map. Default is ``'plasma'``.
    xlabel, ylabel : str, optional
        Axis labels.
    cbar_label : str, optional
        Color-bar label; overrides the one built from the flavor pair.
    cbar_label_prefix : str, optional
        Text placed before the probability label on the color bar.
    cbar_fontsize, cbar_labelsize : float, optional
        Color-bar label and tick-label sizes.
    annotation : str, optional
        Corner annotation. Defaults to the probability label when the flavor
        pair is given; pass ``''`` to suppress it.
    annotation_fontsize : float, optional
        Corner annotation size.
    xlim, ylim : tuple of float, optional
        Axis limits. Default to the data range.
    xmajor, xminor, ymajor, yminor : float, optional
        Tick spacings.
    figsize : tuple of float, optional
        Figure size. Default is ``(9.0, 9.0)``.
    contourf_kw : dict, optional
        Extra keywords for :meth:`~matplotlib.axes.Axes.contourf`.
    subplots_kw : dict, optional
        Extra keywords for :func:`~matplotlib.pyplot.subplots`.
    savefig : str, optional
        If given, the figure is written here.
    savefig_kw : dict, optional
        Extra keywords merged over :data:`HOUSE_SAVEFIG_KW`.
    tight_layout : bool, optional
        Whether to call ``tight_layout``. Default is ``False``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes

    Examples
    --------
    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        import magnus.globaldefs as gd
        from magnus.plotting import plot_oscillogram

        c = np.linspace(-1.0, 0.0, 40)
        lE = np.linspace(-1.0, 1.0, 30)
        P = np.sin(np.outer(10 ** lE, 1.0 + c)) ** 2

        fig, ax = plot_oscillogram(c, lE, P, nu_i=gd.NUMU, nu_f=gd.NUMU)
        print(ax.get_xlabel())
    """
    mpl, plt = _mpl()
    import matplotlib.patheffects as path_effects

    prob = np.asarray(probability)
    expected = (len(log10_energy), len(costhz))
    if prob.shape != expected:
        raise ValueError(
            'Error in magnus: plotting.plot_oscillogram: probability must have shape '
            f'(len(log10_energy), len(costhz)) = {expected}, got {prob.shape}'
        )

    skw = dict(subplots_kw or {})
    gs_kw = skw.pop('gridspec_kw', None) or dict(height_ratios=[1.0],
                                                 width_ratios=[1.0])
    fig, ax = plt.subplots(ncols=1, nrows=1, gridspec_kw=gs_kw,
                           figsize=figsize, **skw)
    fig.subplots_adjust(hspace=_HSPACE, wspace=_WSPACE)

    ckw = dict(levels=levels, cmap=plt.get_cmap(cmap))
    ckw.update(contourf_kw or {})
    cs = ax.contourf(costhz, log10_energy, prob, **ckw)

    label = cbar_label
    if label is None and nu_i is not None and nu_f is not None:
        label = cbar_label_prefix + prob_label(nu_i, nu_f)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.ax.tick_params(labelsize=cbar_labelsize)
    if label is not None:
        cbar.set_label(label=label, fontsize=cbar_fontsize)

    if annotation is None and nu_i is not None and nu_f is not None:
        annotation = prob_label(nu_i, nu_f)
    if annotation:
        text = ax.text(0.96, 0.95, annotation, ha='right', va='center',
                       size=annotation_fontsize, color='k', rotation=0.0,
                       transform=ax.transAxes)
        text.set_path_effects([path_effects.Stroke(linewidth=12,
                                                   foreground='white'),
                               path_effects.Normal()])

    ax.set_xlim(*(xlim if xlim is not None else (min(costhz), max(costhz))))
    ax.set_ylim(*(ylim if ylim is not None
                  else (min(log10_energy), max(log10_energy))))
    _apply_locators(ax.xaxis, xmajor, xminor)
    _apply_locators(ax.yaxis, ymajor, yminor)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    _finish(fig, savefig, savefig_kw, tight_layout)
    return fig, ax
