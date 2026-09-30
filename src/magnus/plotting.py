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

There is deliberately **no bare** ``**kwargs`` **on any of these functions** except
the three presets, which forward verbatim to :func:`plot_curves` under the name
``**_forbidden``, so a typo still raises and names the offending key. A
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
exposes on an axes.  With ``return_probability=True``, the functions that
compute through the Earth wrappers also return what they drew.

Requirements
------------
Matplotlib ships with Magνs, so this module is available in any installation
and needs nothing extra.

.. versionadded:: 1.0.0
"""

from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np

from magnus import _validate as _v

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


@_v.validated(dict(nu_i=_v.r_int(lo=0), nu_f=_v.r_int(lo=0), nubar=_v.r_bool))
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


# ---------------------------------------------------------------------------------------------
# Argument rules (issue #160 §12).  Applied once per figure, to calls from outside the package;
# the presets forward to plot_curves with arguments already checked.  Each names the routine
# the user called and the argument at fault, where Matplotlib would otherwise fail several
# calls later with a message about its own internals -- or draw something meaningless.
# ---------------------------------------------------------------------------------------------

_SCALES = ('linear', 'log', 'symlog', 'logit')
_ENERGY_UNITS = ('eV', 'keV', 'MeV', 'GeV', 'TeV', 'PeV', 'EeV')
# Set by the routines themselves, so a value in subplots_kw would collide with theirs.
_RESERVED_SUBPLOTS_KW = ('nrows', 'ncols', 'figsize')


def _r_pair(positive=False):
    def rule(name, x, where, a):
        if x is None:
            return
        try:
            lo, hi = x
        except (TypeError, ValueError):
            raise ValueError(_v._msg(where, name + " must be a pair of numbers; got " +
                                     repr(x) + ".")) from None
        _v.check_real(name + "[0]", lo, where, positive=positive)
        _v.check_real(name + "[1]", hi, where, positive=positive)
        if not positive and lo == hi:
            raise ValueError(_v._msg(where, name + " must span a range; both ends are " +
                                     repr(lo) + "."))
    return rule


_TICK = _v.r_real(positive=True, allow_none=True)
_FONT = _v.r_real(positive=True)
_OPEN_UNIT = _v.r_real(lo=0.0, hi=1.0, lo_open=True, hi_open=True, what="in (0, 1)")
_BOOL_OR_NONE = lambda name, x, where, a: _v.check_bool(name, x, where, allow_none=True)  # noqa: E731


def _r_subplots_kw(name, x, where, a):
    _v.check_dict(name, x, where)
    for key in _RESERVED_SUBPLOTS_KW:
        if x and key in x:
            raise ValueError(_v._msg(where, name + " must not set " + repr(key) + ", which "
                                     "this routine sets itself" + (" (use figsize=)"
                                     if key == 'figsize' else "") + "."))


def _r_abscissa(name, x, where, a):
    _v.check_real_array(name, x, where, allow_scalar=False)
    scale = a.get('xscale')
    if scale == 'log' and np.min(np.asarray(x, dtype=float)) <= 0.0:
        raise ValueError(_v._msg(where, name + " has entries <= 0, which a log axis cannot "
                                 "show and would drop without a word; pass xscale='linear'."))


def _r_curves(xname, probability=False):
    r"""Each curve 1-D, finite, as long as the abscissa; in [0, 1] if a probability.

    .. versionchanged:: 1.2.0
       A curve with no point above 0 is refused on a log y axis (issue #160 §12).
    """
    def rule(name, x, where, a):
        n = len(np.atleast_1d(a[xname]))
        for i, (y, _) in enumerate(_as_curve_list(x)):
            _v.check_real_array(name + "[" + str(i) + "]", y, where, allow_scalar=False)
            if len(y) != n:
                raise ValueError(_v._msg(where, name + "[" + str(i) + "] has " + str(len(y)) +
                                         " points and " + xname + " has " + str(n) + "."))
            if probability:
                _check_probability(name + "[" + str(i) + "]", y, where)
            # A curve with no positive point has nothing a log axis can show: it was drawn as
            # an empty axis, without a word (issue #160 §12).
            if a.get('yscale') == 'log' and not np.any(np.asarray(y, dtype=float) > 0.0):
                raise ValueError(_v._msg(where, name + "[" + str(i) + "] has no entry above 0, "
                                         "so a log y axis cannot show any of it; pass "
                                         "yscale='linear'."))
    return rule


def _r_panels(xname):
    r"""Each panel of :func:`plot_curves_stacked` checked as :func:`plot_curves` checks its curves.

    .. versionadded:: 1.2.0
    """
    curves = _r_curves(xname)

    def rule(name, x, where, a):
        for j, panel in enumerate(x):
            curves(name + "[" + str(j) + "]", panel, where, a)
    return rule


def _check_probability(name, y, where):
    y = np.asarray(y, dtype=float)
    if np.any(y < -1e-9) or np.any(y > 1.0 + 1e-9):
        raise ValueError(_v._msg(where, name + " is a probability, so it must lie in [0, 1]; "
                                 "its range is [" + format(float(np.min(y)), '.4g') + ", " +
                                 format(float(np.max(y)), '.4g') + "]."))


def _r_probability_array(name, x, where, a):
    if x is None:
        return
    arr = np.asarray(x)
    if arr.dtype.kind not in 'iuf':
        raise _v.InputTypeError(_v._msg(where, name + " must hold real numbers."))
    if not np.all(np.isfinite(arr)):
        raise ValueError(_v._msg(where, name + " must be finite."))
    _check_probability(name, arr, where)


def _r_channel(name, x, where, a):
    if x is None:
        return
    x = _v.check_int(name, x, where, lo=0)
    n = a.get('num_flavors')
    if n is not None and x >= n:
        raise ValueError(_v._msg(where, name + " = " + str(x) + " is not a flavor at " +
                                 str(n) + " flavors (indices run from 0 to " + str(n - 1) +
                                 ")."))


def _r_grid(name, x, where, a):
    _v.check_real_array(name, x, where, allow_scalar=False)
    if len(x) < 2:
        raise ValueError(_v._msg(where, name + " must have at least two entries to draw a "
                                 "contour; got " + str(len(x)) + "."))
    if name == 'costhz' and (np.min(x) < -1.0 or np.max(x) > 1.0):
        raise ValueError(_v._msg(where, "costhz must lie in [-1, 1]."))


def _r_oscillogram_probability(name, x, where, a):
    if x is None:
        return
    shape = (len(a['log10_energy']), len(a['costhz']))
    if np.shape(x) != shape:
        raise ValueError(_v._msg(where, name + " must have shape (len(log10_energy), "
                                 "len(costhz)) = " + str(shape) + "; got " +
                                 str(np.shape(x)) + "."))
    # Finite reals, not [0, 1]: a precomputed map may be any quantity with a matching
    # cbar_label, such as a difference of two probabilities (notebook 06).
    arr = np.asarray(x)
    if arr.dtype.kind not in 'iuf':
        raise _v.InputTypeError(_v._msg(where, name + " must hold real numbers."))
    if not np.all(np.isfinite(arr)):
        raise ValueError(_v._msg(where, name + " must be finite."))


# The rules shared by every routine, by argument name; each routine adds its data rules.
_COMMON_RULES = dict(
    xlim=_r_pair(), ylim=_r_pair(), residual_ylim=_r_pair(), profile_ylim=_r_pair(),
    panel_ylim=_r_pair(), figsize=_r_pair(positive=True), panel_label_xy=_r_pair(),
    panel_annotation_xy=_r_pair(),
    xscale=_v.r_choice(_SCALES), yscale=_v.r_choice(_SCALES), panel_yscale=_v.r_choice(_SCALES),
    xmajor=_TICK, xminor=_TICK, ymajor=_TICK, yminor=_TICK, residual_ymajor=_TICK,
    residual_yminor=_TICK, profile_ymajor=_TICK, profile_yminor=_TICK, panel_ymajor=_TICK,
    panel_yminor=_TICK,
    title_fontsize=_FONT, cbar_fontsize=_FONT, cbar_labelsize=_FONT,
    annotation_fontsize=_FONT, panel_annotation_fontsize=_FONT,
    ylabel_labelpad=_v.r_real(), shared_ylabel_labelpad=_v.r_real(),
    residual_height=_OPEN_UNIT, profile_height=_OPEN_UNIT,
    legend=_v.r_bool, grid=_v.r_bool, tight_layout=_v.r_bool, return_probability=_v.r_bool,
    show_profile=_BOOL_OR_NONE, panel_per_trajectory=_BOOL_OR_NONE, nubar=_BOOL_OR_NONE,
    legend_panel=_v.r_int(lo=0), legend_on_panel=_v.r_int(lo=-1), levels=_v.r_int(lo=1),
    num_flavors=_v.r_int(allow_none=True), nu_i=_r_channel, nu_f=_r_channel,
    energy_unit=_v.r_choice(_ENERGY_UNITS),
    x_unit=_v.r_real(positive=True, allow_none=True),
    energy=_v.r_real(positive=True, allow_none=True),
    electron_fraction=_v.r_real(lo=0.0, hi=1.0, lo_open=True, allow_none=True),
    electron_fraction_core=_v.r_real(lo=0.0, hi=1.0, lo_open=True, allow_none=True),
    electron_fraction_mantle=_v.r_real(lo=0.0, hi=1.0, lo_open=True, allow_none=True),
    electron_fraction_crust=_v.r_real(lo=0.0, hi=1.0, lo_open=True, allow_none=True),
    electron_fraction_ocean=_v.r_real(lo=0.0, hi=1.0, lo_open=True, allow_none=True),
    ratio_number_neutrons_to_protons=_v.r_real(nonnegative=True, allow_none=True),
    dcp=_v.r_real_array(), subplots_kw=_r_subplots_kw, legend_kw=_v.r_dict,
    grid_kw=_v.r_dict, savefig_kw=_v.r_dict, residual_kw=_v.r_dict, contourf_kw=_v.r_dict,
    wrapper_kw=_v.r_dict, osc_params=_v.r_dict,
)


def _rules(**data_rules):
    return dict(_COMMON_RULES, **data_rules)


def _check_forwarded(where: str, kw: dict) -> None:
    r"""Refuse, naming the preset, a keyword :func:`plot_curves` does not take (#146 §2).

    Without this the TypeError came from plot_curves, a function the caller never called.
    """
    import inspect
    known = inspect.signature(plot_curves).parameters
    unknown = [k for k in kw if k not in known]
    if unknown:
        raise TypeError("Error in magnus: plotting." + where + "() got unexpected keyword "
                        "argument" + ("s " if len(unknown) > 1 else " ") +
                        ", ".join(repr(k) for k in unknown) + (".  Each call makes its own "
                        "figure, so there is no ax= to draw into." if 'ax' in unknown else "."))


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


@_v.validated(_rules(x=_r_abscissa, curves=_r_curves('x'), residual=_v.r_real_array(allow_scalar=False)))
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


@_v.validated(_rules(x=_r_abscissa, panels=_r_panels('x')))
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


@_v.validated(_rules(distances=_r_abscissa, curves=_r_curves('distances', probability=True)))
def plot_probability_vs_baseline(
    distances: Sequence[float],
    curves: Sequence[Union[Sequence[float], Dict[str, Any]]],
    *,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    nubar: Optional[bool] = None,
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

    .. versionchanged:: 1.2.0
       Takes nubar, for the default ordinate label (issue #145 §2).

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
    nubar : bool, optional
        True if the curves are antineutrino probabilities, for the default ordinate label.
        Default: None, read as False.
    xlabel : str, optional
        Abscissa label.
    ylabel : str, optional
        Ordinate label; overrides the one built from the flavor pair.
    ylim : tuple of float, optional
        Ordinate limits. Default is ``(0.0, 1.0)``.
    xscale : str, optional
        Abscissa scale. Default is ``'log'``.
    ymajor, yminor : float, optional
        Ordinate tick spacings.  Pass None to hand the axis back to Matplotlib's own
        locator. Default: 0.10 and 0.02.

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
    _check_forwarded('plot_probability_vs_baseline', _forbidden)
    if ylabel is None and nu_i is not None and nu_f is not None:
        ylabel = _probability_ylabel(nu_i, nu_f, num_flavors, nubar=bool(nubar))
    return plot_curves(
        distances, curves, xlabel=xlabel, ylabel=ylabel, ylim=ylim,
        xscale=xscale, ymajor=ymajor, yminor=yminor, **_forbidden)


@_v.validated(_rules(energies=_r_abscissa, curves=_r_curves('energies', probability=True)))
def plot_probability_vs_energy(
    energies: Sequence[float],
    curves: Sequence[Union[Sequence[float], Dict[str, Any]]],
    *,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    nubar: Optional[bool] = None,
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

    .. versionchanged:: 1.2.0
       Takes nubar, for the default ordinate label (issue #145 §2); energies that look like eV
       under a larger unit raise EnergyUnitWarning (issue #160 §12).

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
    nubar : bool, optional
        True if the curves are antineutrino probabilities, for the default ordinate label.
        Default: None, read as False.
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
        Ordinate tick spacings.  Pass None to hand the axis back to Matplotlib's own
        locator. Default: 0.10 and 0.02.

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
    _check_forwarded('plot_probability_vs_energy', _forbidden)
    # Energies left in eV, the package's own unit, under the default 'GeV' label: axis values
    # a billion times too large, labelled as if right (issue #160 §12).  Read as eV-sized when
    # the numbers reach 1e6 and would put the axis above 100 PeV in the unit named.
    _e_max = float(np.max(np.asarray(energies, dtype=float)))
    _scale = {'eV': 1.0, 'keV': 1e3, 'MeV': 1e6, 'GeV': 1e9, 'TeV': 1e12, 'PeV': 1e15,
              'EeV': 1e18}[energy_unit]
    if energy_unit != 'eV' and _e_max >= 1.0e6 and _e_max*_scale >= 1.0e17:
        import warnings
        import magnus.globaldefs as gd
        warnings.warn(gd.WARNING_MSG_NO_COLOR + " plotting.plot_probability_vs_energy: the "
            "energies reach " + format(_e_max, '.4g') + ", labelled as " + energy_unit + ", "
            "which is " + format(_e_max*_scale, '.3g') + " eV.  They look like energies in eV, "
            "the package's own unit, left unconverted: divide by gd.UNIT_" + energy_unit.upper() +
            " or pass energy_unit='eV'.", gd.EnergyUnitWarning, stacklevel=2)
    if xlabel is None:
        xlabel = r'Neutrino energy, $E_\nu$ [%s]' % energy_unit
    if ylabel is None and nu_i is not None and nu_f is not None:
        ylabel = _probability_ylabel(nu_i, nu_f, num_flavors, nubar=bool(nubar))
    return plot_curves(
        energies, curves, xlabel=xlabel, ylabel=ylabel, ylim=ylim,
        xscale=xscale, ymajor=ymajor, yminor=yminor, **_forbidden)


_FLAVOR_WORD = {2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five'}


def _resolve_nubar(caller, nubar, wrapper_kw, computing):
    r"""The ``nubar`` a figure is drawn and labelled for, and the ``wrapper_kw`` to compute with.

    ``nubar=None`` means not given.  In compute mode a given ``nubar`` is also passed to the
    wrapper, so that the grid and its label cannot disagree; one that contradicts
    ``wrapper_kw['nubar']`` is refused (issue #145 §2).

    .. versionadded:: 1.2.0
    """
    in_kw = (wrapper_kw or {}).get('nubar')
    if nubar is not None and in_kw is not None and bool(in_kw) != bool(nubar):
        raise ValueError('Error in magnus: plotting.%s: nubar=%r contradicts '
                         'wrapper_kw[\'nubar\']=%r; give it once.' % (caller, nubar, in_kw))
    if nubar is not None and computing:
        wrapper_kw = dict(wrapper_kw or {}, nubar=bool(nubar))
    effective = bool(nubar) if nubar is not None else bool(in_kw)
    return effective, wrapper_kw


def _probability_ylabel(nu_i, nu_f, num_flavors, nubar=False):
    r"""Build the notebooks' ordinate label for a probability panel.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.2.0
       Takes nubar and labels an antineutrino probability as such (issue #145).
    """
    label = prob_label(nu_i, nu_f, nubar=bool(nubar))
    if num_flavors is None:
        return 'Probability, ' + label
    if num_flavors not in _FLAVOR_WORD:
        raise ValueError(
            'Error in magnus: plotting: num_flavors must be one of '
            f'{sorted(_FLAVOR_WORD)}, not {num_flavors!r}'
        )
    return f'{_FLAVOR_WORD[num_flavors]}-neutrino probability, ' + label


@_v.validated(_rules(x=_r_abscissa))
def plot_probability_with_profile(
    x: Sequence[float],
    profiles: Optional[Sequence[Union[Sequence[float], Dict[str, Any]]]] = None,
    panels: Optional[Sequence[Sequence[Union[Sequence[float], Dict[str, Any]]]]] = None,
    *,
    trajectories: Optional[Sequence[Dict[str, Any]]] = None,
    x_axis: Optional[str] = None,
    x_unit: Optional[float] = None,
    energy: Optional[float] = None,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    nubar: Optional[bool] = None,
    num_flavors: Optional[int] = None,
    osc_params: Optional[Dict[str, float]] = None,
    electron_fraction: Optional[float] = None,
    electron_fraction_core: Optional[float] = None,
    electron_fraction_mantle: Optional[float] = None,
    electron_fraction_crust: Optional[float] = None,
    electron_fraction_ocean: Optional[float] = None,
    ratio_number_neutrons_to_protons: Optional[float] = None,
    wrapper_kw: Optional[Dict[str, Any]] = None,
    show_profile: Optional[bool] = None,
    panel_per_trajectory: Optional[bool] = None,
    return_probability: bool = False,
    xlabel: Optional[str] = None,
    profile_ylabel: str = r'$\frac{N_e}{N_{\rm Av}}$ [cm$^{-3}$]',
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

    It draws curves you computed, or, given ``trajectories`` instead of
    ``profiles`` and ``panels``, computes them through the Earth wrappers
    (:func:`magnus.oscprob.osc_prob_3nu_earth` and its two-, four- and
    five-flavor siblings), which declare the PREM layer boundaries as slab
    edges and use a layered electron fraction unless told otherwise.  The
    density panel then shows the electron density those wrappers integrate
    along each trajectory.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.1
       Computes the probabilities through the Earth wrappers when given
       ``trajectories``; ``profiles`` and ``panels`` default to None.

    .. versionchanged:: 1.2.0
       Takes nubar (issue #145 §2), labels the computed probability panel (issue #146 §2), and
       computes NSI and LIV through the Earth wrappers (issue #146 §1).

    Parameters
    ----------
    x : sequence of float
        Shared abscissa, or, when the panels have different abscissae, the one
        used by the profile panel. Individual curves may carry their own ``x``.
        When computing, baselines measured from where the neutrino enters the
        Earth (``x_axis='baseline'``) or energies (``x_axis='energy'``), in
        units of ``x_unit``.
    profiles : sequence or None
        Curves for the density panel, in the form :func:`plot_curves` takes.
        A curve may add its own abscissa under ``'x'``. Pass ``None`` (or an
        empty sequence) to omit the density panel entirely and get a plain
        stack of probability panels sharing an abscissa -- the layout the
        long-baseline notebook uses for a probability above its
        energy-smoothed version.  Leave it None when computing.
    panels : sequence of sequence
        One entry per probability panel; each entry is a sequence of curves.
        Required unless computing, and then left None.
    trajectories : sequence of dict, optional
        Compute instead of draw: one curve per entry, each a path through the
        Earth given as ``'costhz'`` (the cosine of the zenith angle) or as
        ``'loc_ini'`` and ``'loc_fin'`` (two named locations of
        :data:`magnus.earth.loc_coords_dms`, or their coordinates, joined by
        the chord between them).  An entry may carry its own abscissa under
        ``'x'`` -- a baseline grid ending at its own chord length, say -- and
        any :class:`~matplotlib.lines.Line2D` keyword (``color``, ``ls``,
        ``label``...) for its curve.
    x_axis : str, optional
        What the abscissa is when computing, ``'baseline'`` or ``'energy'``.
        ``'baseline'`` (the default): one
        wrapper call per trajectory over all its baselines, at ``energy``; a
        baseline longer than the trajectory's chord is an error.
        ``'energy'``: one call per trajectory over all its energies, at the
        full chord.
    x_unit : float, optional
        What one unit of the abscissa is, in the wrappers' units: a length in
        :math:`\text{eV}^{-1}` or an energy in eV.  Default
        :data:`magnus.globaldefs.UNIT_KM` for baselines (``x`` in km) and
        :data:`magnus.globaldefs.UNIT_GEV` for energies (``x`` in GeV).
    energy : float, optional
        Neutrino energy [eV], required on the baseline axis and refused on the
        energy axis.
    nu_i, nu_f : int, optional
        Initial and final flavors when computing, e.g.
        :data:`magnus.globaldefs.NUMU` and :data:`magnus.globaldefs.NUE`.
    num_flavors : int, optional
        2, 3, 4 or 5: which Earth wrapper computes the probabilities.
    nubar : bool, optional
        True for antineutrinos.  In compute mode it is also passed to the wrapper, so the
        curves and their label agree; a value contradicting ``wrapper_kw['nubar']`` is refused.
        Default: None, which reads ``wrapper_kw['nubar']`` in compute mode and False otherwise.
    osc_params : dict, optional
        Mixing parameters, passed to the wrapper as keywords: ``sth`` and
        ``Dm2`` (required) at two flavors; at three to five flavors any of the
        wrapper's own (``s12``, ..., ``D41``...), the rest taking the wrapper's
        defaults.  :func:`magnus.globaldefs.load_nufit_params` returns such a
        dict.
    electron_fraction, electron_fraction_core, electron_fraction_mantle, electron_fraction_crust, electron_fraction_ocean : float, optional
        The Earth's composition, as the wrappers take it: one :math:`Y_e` for
        the whole Earth, or per PREM layer.  Default: None, the layered values.
        The density panel follows them too.
    ratio_number_neutrons_to_protons : float, optional
        Passed to the wrapper; matters from four flavors up.  Default: None,
        derived from the same :math:`Y_e`.
    wrapper_kw : dict, optional
        Any other wrapper keyword, e.g. ``rtol`` or ``nubar``.
    show_profile : bool, optional
        Whether to draw the density panel when computing.  Default: yes on the
        baseline axis, no on the energy axis, where the abscissa is not a
        position along the trajectory (so ``True`` there is an error).
    panel_per_trajectory : bool, optional
        When computing, one probability panel per trajectory (the default) or
        all curves in one panel.
    return_probability : bool, optional
        When computing, also return the probabilities drawn.  Default False.
    xlabel : str, optional
        Abscissa label, placed under the bottom panel.  Default
        ``'Baseline, $L$ [km]'``, or, when computing on the energy axis, a
        neutrino-energy label in GeV or MeV (for those two ``x_unit`` values).
    profile_ylabel : str, optional
        Ordinate label of the density panel.
    panel_ylabels : sequence of str, optional
        Ordinate labels for the probability panels. Entries may be ``None``.  In compute
        mode, when neither this nor ``shared_ylabel`` is given, every probability panel is
        labelled with the channel of ``nu_i`` and ``nu_f`` (and ``nubar``).
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
        Title, placed above the top panel, which is the density panel only when
        ``profiles`` is given.
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
        Padding of the density panel's ordinate label.  The probability panels use a
        fixed padding of 15 and do not read this. Default: 25.0.
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
    probability : list of numpy.ndarray
        Only with ``return_probability=True``: one array per trajectory, over
        its abscissa.

    Raises
    ------
    ValueError
        If ``panels`` is missing or empty when drawing; if ``trajectories`` is
        given together with ``profiles`` or ``panels``, or any computing
        argument without ``trajectories``; if a trajectory names neither or
        both of ``costhz`` and the two locations; if ``x_axis`` is unknown,
        ``energy`` is missing on the baseline axis or given on the energy axis,
        or ``show_profile=True`` on the energy axis; if ``trajectories`` is
        empty, or ``wrapper_kw`` holds ``source_depth`` or ``detector_depth``;
        and for the channel, flavor-count and duplicate-keyword checks of
        :func:`plot_oscillogram`.

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

    Computed through the Earth wrappers instead, along a chord that crosses
    the core, with the electron density it samples on top:

    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        import magnus.globaldefs as gd
        from magnus.plotting import plot_probability_with_profile

        L = np.linspace(100.0, 11000.0, 200)    # [km]; the chord is 11 467 km long
        fig, ax, P = plot_probability_with_profile(
            L, trajectories=[dict(costhz=-0.9, label=r'$\cos\theta_z = -0.9$')],
            energy=5.0*gd.UNIT_GEV, nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3,
            xscale='linear', return_probability=True)
        print(len(ax), P[0].shape)
    """
    nubar_label, wrapper_kw = _resolve_nubar('plot_probability_with_profile', nubar,
                                             wrapper_kw, trajectories is not None)
    computing = dict(x_axis=x_axis, x_unit=x_unit, energy=energy, nu_i=nu_i, nu_f=nu_f,
                     num_flavors=num_flavors, osc_params=osc_params,
                     electron_fraction=electron_fraction,
                     electron_fraction_core=electron_fraction_core,
                     electron_fraction_mantle=electron_fraction_mantle,
                     electron_fraction_crust=electron_fraction_crust,
                     electron_fraction_ocean=electron_fraction_ocean,
                     ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
                     wrapper_kw=wrapper_kw, show_profile=show_profile,
                     panel_per_trajectory=panel_per_trajectory)
    probability = None
    if trajectories is not None:
        if profiles is not None or panels is not None:
            raise ValueError(
                'Error in magnus: plotting.plot_probability_with_profile: give either curves '
                'to draw (profiles, panels) or trajectories to compute them from, not both.')
        profiles, panels, probability, computed_xlabel = _profile_through_earth_wrappers(
            x, trajectories, **computing)
        if xlabel is None:
            xlabel = computed_xlabel
        # The channel is known here, so the probability panels get the label the curve
        # routines give (issue #146 §2); it used to be left empty.
        if (panel_ylabels is None and shared_ylabel is None and nu_i is not None
                and nu_f is not None):
            panel_ylabels = [_probability_ylabel(nu_i, nu_f, num_flavors,
                                                 nubar=nubar_label)]*len(panels)
    else:
        given = [name for name, value in computing.items() if value is not None]
        if return_probability:
            given.append('return_probability')
        if given:
            raise ValueError(
                'Error in magnus: plotting.plot_probability_with_profile: '
                + ', '.join(given) + ' compute the probabilities, which needs trajectories; '
                'give either curves to draw or trajectories to compute them from, not both.')
        if panels is None:
            raise ValueError('Error in magnus: plotting.plot_probability_with_profile: give '
                             'panels to draw, or trajectories to compute them from.')
    if xlabel is None:
        xlabel = r'Baseline, $L$ [km]'

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
    if return_probability:
        return fig, ax, probability
    return fig, ax


_TRAJECTORY_GEOMETRY = ('costhz', 'loc_ini', 'loc_fin')


def _costhz_of_trajectory(caller, entry):
    r"""The cosine of the zenith angle of one trajectory: given, or of the chord between two
    surface locations (named, or as coordinates), resolved as the Earth wrappers do."""
    locations = [k for k in ('loc_ini', 'loc_fin') if k in entry]
    if ('costhz' in entry) == bool(locations) or len(locations) == 1:
        raise ValueError(
            'Error in magnus: plotting.%s: each trajectory needs either costhz or both '
            'loc_ini and loc_fin, and not both; got keys %s.' % (caller, sorted(entry)))
    if 'costhz' in entry:
        return float(entry['costhz'])
    from magnus import earth
    ends = []
    for key in ('loc_ini', 'loc_fin'):
        loc = entry[key]
        if isinstance(loc, str):
            loc = earth.coordinates_of_named_location(caller, loc_name=loc)
        ends.append(loc)
    (lat1, lon1), (lat2, lon2) = ends
    return float(earth.costhz_between_points_on_surface(lat1, lon1, lat2, lon2))


def _profile_through_earth_wrappers(x, trajectories, x_axis, x_unit, energy, nu_i, nu_f,
                                    num_flavors, osc_params, wrapper_kw, show_profile,
                                    panel_per_trajectory, **composition):
    r"""Density curves, probability panels, the probabilities and a default abscissa label.

    One wrapper call per trajectory: over its baselines at one energy (the wrapper's
    baseline-array path), or over its energies at the full chord (the energy-batched
    engine).  The density drawn is the one the wrapper integrates: the same PREM
    profile, electron fraction and ocean density, from the same resolution.
    """
    caller = 'plot_probability_with_profile'
    where = 'Error in magnus: plotting.%s: ' % caller
    x_axis = 'baseline' if x_axis is None else x_axis
    if x_axis not in ('baseline', 'energy'):
        raise ValueError(where + "x_axis is 'baseline' or 'energy', not %r." % (x_axis,))
    on_baseline = x_axis == 'baseline'
    if on_baseline and energy is None:
        raise ValueError(where + 'on the baseline axis, computing needs energy [eV].')
    if not on_baseline and energy is not None:
        raise ValueError(where + 'on the energy axis the energies are the abscissa; energy '
                         'cannot also be given.')
    if show_profile is None:
        show_profile = on_baseline
    elif show_profile and not on_baseline:
        raise ValueError(where + 'show_profile=True needs the baseline axis: on the energy '
                         'axis the abscissa is not a position along the trajectory.')
    if len(trajectories) == 0:
        raise ValueError(where + 'trajectories is empty; at least one is required.')
    depths = sorted({'source_depth', 'detector_depth'} & set(wrapper_kw or {}))
    if depths:
        raise ValueError(where + 'the trajectories here run from surface to surface; '
                         + ', '.join(depths) + ' cannot be given.')
    fn, kwargs = _earth_wrapper_and_arguments(
        caller, nu_i, nu_f, num_flavors, osc_params, wrapper_kw, composition,
        reserved={'energy', 'costhz', 'L', 'loc_ini', 'loc_fin'})

    from magnus import earth, globaldefs as gd, oscprob
    if x_unit is None:
        x_unit = gd.UNIT_KM if on_baseline else gd.UNIT_GEV
    if on_baseline:
        xlabel = r'Baseline, $L$ [km]' if x_unit == gd.UNIT_KM else r'Baseline, $L$'
    else:
        unit = {gd.UNIT_GEV: ' [GeV]', gd.UNIT_MEV: ' [MeV]'}.get(x_unit, '')
        xlabel = r'Neutrino energy, $E_\nu$' + unit
    # n_e in units of N_Av per cm^3, which reads as rho*Y_e in g cm^-3
    per_n_av = gd.N_AV/gd.CONV_CM_TO_INV_EV**3

    profiles, curves, probability = [], [], []
    for entry in trajectories:
        entry = dict(entry)
        costhz = _costhz_of_trajectory(caller, entry)
        for key in _TRAJECTORY_GEOMETRY:
            entry.pop(key, None)
        xi = np.asarray(entry.pop('x', x), dtype=float)
        if on_baseline:
            prob = fn(energy, costhz=costhz, L=xi*x_unit, **kwargs)
        else:
            chord = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
            prob = fn(xi*x_unit, costhz=costhz, L=chord, **kwargs)
        prob = np.asarray(prob, dtype=float).reshape(len(xi))
        probability.append(prob)
        curves.append(dict(entry, x=xi, y=prob))
        if show_profile:
            # num_flavors=3: the wrapper call above already raised any composition
            # warning, and the electron density does not depend on the neutron ratio
            rho_func, _ = oscprob._earth_composition(
                costhz, composition.get('electron_fraction'), None,
                composition.get('electron_fraction_core'),
                composition.get('electron_fraction_mantle'),
                composition.get('electron_fraction_crust'),
                composition.get('electron_fraction_ocean'),
                caller, num_flavors=3,
                density_matter_ocean=(wrapper_kw or {}).get('density_matter_ocean'))
            style = {k: v for k, v in entry.items() if k != 'label'}
            profiles.append(dict(style, x=xi, y=rho_func(xi*x_unit)/per_n_av))

    if panel_per_trajectory is None or panel_per_trajectory:
        panels = [[curve] for curve in curves]
    else:
        panels = [curves]
    return profiles, panels, probability, xlabel


@_v.validated(_rules(x=_r_abscissa, probabilities=_r_probability_array, averages=_r_probability_array))
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
    curves, each with its phase-averaged value drawn through it as a dashed line of
    the same color -- for instance what :func:`magnus.oscprob.osc_prob` returns with
    ``average=True`` (the average over ``average_spread``).

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
    _check_forwarded('plot_probability_with_average', _forbidden)
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

    # Three keywords are handled here rather than forwarded.  `legend` collided with the
    # legend=False below, so a caller who forwarded it got "multiple values" and could not
    # suppress the legend at all.  `savefig` would have written the file inside the call
    # below, before the averaged entry is added to the legend, so the saved figure was
    # missing the very thing the plot exists to show.
    _legend = _forbidden.pop('legend', True)
    _savefig = _forbidden.pop('savefig', None)
    _savefig_kw = _forbidden.pop('savefig_kw', None)

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
    if _legend:
        main.legend(handles, labs, **lkw)
    # Saved last, so the file carries the averaged entry.
    _finish(fig, _savefig, _savefig_kw, False)
    return fig, ax


@_v.validated(_rules(prob_nu=_r_probability_array, prob_nubar=_r_probability_array))
def plot_biprobability(
    prob_nu: Optional[Sequence[Sequence[float]]] = None,
    prob_nubar: Optional[Sequence[Sequence[float]]] = None,
    *,
    configurations: Optional[Sequence[Dict[str, Any]]] = None,
    dcp: Optional[Sequence[float]] = None,
    energy: Optional[float] = None,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    num_flavors: Optional[int] = None,
    osc_params: Optional[Dict[str, float]] = None,
    electron_fraction: Optional[float] = None,
    electron_fraction_core: Optional[float] = None,
    electron_fraction_mantle: Optional[float] = None,
    electron_fraction_crust: Optional[float] = None,
    electron_fraction_ocean: Optional[float] = None,
    ratio_number_neutrons_to_protons: Optional[float] = None,
    wrapper_kw: Optional[Dict[str, Any]] = None,
    return_probability: bool = False,
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

    It draws probabilities you computed, or, given ``configurations`` instead,
    computes them through the Earth wrappers
    (:func:`magnus.oscprob.osc_prob_3nu_earth` and its four- and five-flavor
    siblings), which declare the PREM layer boundaries as slab edges and use a
    layered electron fraction unless told otherwise.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.1
       Computes the probabilities through the Earth wrappers when given
       ``configurations``; markers may name a phase (``'dcp'``); ``nu_i`` and
       ``nu_f`` set the default axis labels.

    Parameters
    ----------
    prob_nu, prob_nubar : sequence of sequence of float
        One entry per curve, each a sequence of probabilities over the same
        grid of :math:`\delta_{\rm CP}` values.  Required unless computing,
        and then left None.
    configurations : sequence of dict, optional
        Compute instead of draw: one curve per entry, each a path through the
        Earth -- ``'costhz'``, with ``'L'`` [:math:`\text{eV}^{-1}`] or, by
        default, the full chord; or ``'loc_ini'`` and ``'loc_fin'``, two
        locations joined by the chord between them -- plus, optionally, its
        own ``'energy'`` and any other wrapper keyword for that curve alone
        (the inverted-ordering parameters for an inverted-ordering curve, say).
        A configuration's entry overrides the shared ``energy``,
        ``osc_params``, composition or ``wrapper_kw`` value for its curve.
        Style the curves with ``labels`` and ``curve_kw``, as when drawing.
    dcp : sequence of float, optional
        The phases the curves run over, in the wrapper's convention (radians,
        or degrees under ``angles='deg'``).  Default: 100 points from
        :math:`-\pi` to :math:`\pi`.
    energy : float, optional
        Neutrino energy [eV] shared by the configurations that give none.
    nu_i, nu_f : int, optional
        The channel: the default axis labels, and, when computing, what is
        computed.  Default :math:`\nu_\mu \to \nu_e`.  Give both or
        neither.
    num_flavors : int, optional
        3, 4 or 5: which Earth wrapper computes the probabilities (two flavors
        have no CP phase).
    osc_params : dict, optional
        Mixing parameters shared by all curves, as the wrapper takes them.  A
        ``dCP`` in it, or in a configuration, is replaced by the phases in
        ``dcp`` (and a marker's), so a set from
        :func:`magnus.globaldefs.load_nufit_params` can be passed whole.
    electron_fraction, electron_fraction_core, electron_fraction_mantle, electron_fraction_crust, electron_fraction_ocean : float, optional
        The Earth's composition, as the wrappers take it.  Default: None, the
        layered values.
    ratio_number_neutrons_to_protons : float, optional
        Passed to the wrapper; matters from four flavors up.
    wrapper_kw : dict, optional
        Any other wrapper keyword, e.g. ``rtol``.
    return_probability : bool, optional
        Also return the probabilities drawn.  Default False.
    labels : sequence of str, optional
        Legend label per curve.
    curve_kw : sequence of dict, optional
        Per-curve :class:`~matplotlib.lines.Line2D` keywords.
    markers : sequence of dict, optional
        Markers at selected phases. Each entry gives its position either as
        ``'index'`` (a position along the curve) or as ``'xy'`` (an explicit
        coordinate pair, which is what you have when the marked phases were
        computed separately from the curve) or, when computing, as ``'dcp'``
        (a phase, computed for each curve it marks). Optionally ``'marker'``,
        ``'label'``, ``'filled'`` and ``'curve'`` (which curve it belongs to,
        default all).
    xlabel, ylabel : str, optional
        Axis labels. Default to the ``nu_i``, ``nu_f`` pair.
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
    prob_nu, prob_nubar : numpy.ndarray or sequence
        Only with ``return_probability=True``: the probabilities drawn, arrays
        of shape ``(len(configurations), len(dcp))`` when computed, the inputs
        as given otherwise.

    Raises
    ------
    ValueError
        If the two probabilities have different numbers of curves, or one is
        missing when drawing; if ``configurations`` is given together with
        them, or a computing argument without ``configurations``; if a marker
        gives none or more than one of ``index``, ``xy`` and ``dcp``, or
        ``dcp`` when drawing; if a configuration has no energy, no geometry or
        both geometries, ``L`` with two locations, or a per-point keyword
        (``nubar``, ``nu_i``, ``nu_f``); if only one of ``nu_i``,
        ``nu_f`` is given; and for the flavor-count and duplicate-keyword
        checks of :func:`plot_oscillogram`, here with ``num_flavors`` 3, 4 or 5.

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

    Computed through the Earth wrappers instead, from Fermilab to Homestake at
    2 GeV, both orderings, with the phase :math:`\delta_{\rm CP} = 0` marked:

    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        import magnus.globaldefs as gd
        from magnus.plotting import plot_biprobability

        where = dict(loc_ini='fermilab', loc_fin='homestake')
        fig, ax, P_nu, P_nubar = plot_biprobability(
            configurations=[dict(where, **gd.load_nufit_params('NuFIT 6.1', 'NO')),
                            dict(where, **gd.load_nufit_params('NuFIT 6.1', 'IO'))],
            dcp=np.linspace(-np.pi, np.pi, 25), energy=2.0*gd.UNIT_GEV, num_flavors=3,
            labels=['NO', 'IO'], markers=[dict(dcp=0.0, marker='*', label='0')],
            return_probability=True)
        print(P_nu.shape)
    """
    import magnus.globaldefs as gd

    where = 'Error in magnus: plotting.plot_biprobability: '
    if (nu_i is None) != (nu_f is None):
        raise ValueError(where + 'give both nu_i and nu_f, or neither.')
    channel = (gd.NUMU, gd.NUE) if nu_i is None else (nu_i, nu_f)
    computing = dict(dcp=dcp, energy=energy, num_flavors=num_flavors, osc_params=osc_params,
                     electron_fraction=electron_fraction,
                     electron_fraction_core=electron_fraction_core,
                     electron_fraction_mantle=electron_fraction_mantle,
                     electron_fraction_crust=electron_fraction_crust,
                     electron_fraction_ocean=electron_fraction_ocean,
                     ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
                     wrapper_kw=wrapper_kw)
    points = {}
    if configurations is not None:
        if prob_nu is not None or prob_nubar is not None:
            raise ValueError(where + 'give either probabilities to draw or configurations to '
                             'compute them from, not both.')
        prob_nu, prob_nubar, points = _biprobability_through_earth_wrappers(
            configurations, markers, channel, **computing)
    else:
        given = [name for name, value in computing.items() if value is not None]
        if given:
            raise ValueError(where + ', '.join(given) + ' compute the probabilities, which '
                             'needs configurations; give either probabilities to draw or '
                             'configurations to compute them from, not both.')
        if prob_nu is None or prob_nubar is None:
            raise ValueError(where + 'give prob_nu and prob_nubar to draw, or configurations '
                             'to compute them from.')
        if any('dcp' in m for m in (markers or [])):
            raise ValueError(where + "a marker's 'dcp' is a phase to compute, which needs "
                             'configurations; mark a drawn curve by its index or xy.')

    _, plt = _mpl()
    if len(prob_nu) != len(prob_nubar):
        raise ValueError(
            'Error in magnus: plotting.plot_biprobability: prob_nu and prob_nubar must have '
            f'the same number of curves, got {len(prob_nu)} and '
            f'{len(prob_nubar)}'
        )
    if xlabel is None:
        xlabel = prob_label(*channel)
    if ylabel is None:
        ylabel = prob_label(*channel, nubar=True)

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

    for j, m in enumerate(markers or []):
        if sum(key in m for key in ('index', 'xy', 'dcp')) != 1:
            raise ValueError(
                "Error in magnus: plotting.plot_biprobability: each marker needs one of an "
                "'index' along the curve, an explicit 'xy' coordinate pair, or, when "
                f"computing, a phase 'dcp'; got keys {sorted(m)}"
            )
        filled = m.get('filled', True)
        which = m.get('curve')
        targets = range(len(prob_nu)) if which is None else [which]
        for i in targets:
            c = colors[i]
            if 'dcp' in m:
                x, y = points[(j, i)]
            elif 'xy' in m:
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
    if return_probability:
        return fig, ax, prob_nu, prob_nubar
    return fig, ax


def _biprobability_through_earth_wrappers(configurations, markers, channel, dcp, energy,
                                          num_flavors, osc_params, wrapper_kw, **composition):
    r"""Neutrino and antineutrino probabilities, ``(len(configurations), len(dcp))`` each,
    and the marked points, keyed by (marker index, curve index).

    Two calls per phase and curve, one per ``nubar``: the phase changes the Hamiltonian,
    so nothing is batched across it.
    """
    caller = 'plot_biprobability'
    where = 'Error in magnus: plotting.%s: ' % caller
    fn, shared = _earth_wrapper_and_arguments(
        caller, channel[0], channel[1], num_flavors, osc_params, wrapper_kw, composition,
        reserved={'energy', 'costhz', 'L', 'loc_ini', 'loc_fin', 'nubar'},
        flavors=(3, 4, 5))
    # A parameter set such as load_nufit_params returns carries its own dCP; the phases
    # scanned take its place, which is what this plot is (documented under osc_params).
    shared.pop('dCP', None)
    if len(configurations) == 0:
        raise ValueError(where + 'configurations is empty; at least one is required.')

    from magnus import earth, globaldefs as gd
    dcp = np.linspace(-np.pi, np.pi, 100) if dcp is None else np.asarray(dcp, dtype=float)
    marked = [(j, m) for j, m in enumerate(markers or []) if 'dcp' in m]
    prob_nu = np.empty((len(configurations), len(dcp)))
    prob_nubar = np.empty_like(prob_nu)
    points = {}
    for i, conf in enumerate(configurations):
        conf = dict(conf)
        conf.pop('dCP', None)
        per_point = sorted({'nubar', 'nu_i', 'nu_f'} & set(conf))
        if per_point:
            raise ValueError(where + 'these are set by plot_biprobability itself and cannot be '
                             'given in a configuration: ' + ', '.join(per_point) + '.')
        E = conf.pop('energy', energy)
        if E is None:
            raise ValueError(where + 'configuration %d gives no energy, and no shared energy '
                             'is given.' % i)
        locations = [k for k in ('loc_ini', 'loc_fin') if k in conf]
        if ('costhz' in conf) == bool(locations) or len(locations) == 1:
            raise ValueError(where + 'each configuration needs either costhz or both loc_ini '
                             'and loc_fin, and not both; got keys %s.' % sorted(conf))
        if locations and 'L' in conf:
            raise ValueError(where + 'L goes with costhz; two locations fix the baseline '
                             'themselves (configuration %d).' % i)
        if 'costhz' in conf and 'L' not in conf:
            conf['L'] = earth.distance_traveled_inside_earth(float(conf['costhz']))*gd.UNIT_KM
        kwargs = {**shared, **conf}

        def prob(phase, nubar):
            return np.asarray(fn(E, dCP=phase, nubar=nubar, **kwargs), dtype=float).item()

        prob_nu[i] = [prob(d, False) for d in dcp]
        prob_nubar[i] = [prob(d, True) for d in dcp]
        for j, m in marked:
            if m.get('curve') in (None, i):
                points[(j, i)] = (prob(m['dcp'], False), prob(m['dcp'], True))
    return prob_nu, prob_nubar, points


@_v.validated(_rules(costhz=_r_grid, log10_energy=_r_grid, probability=_r_oscillogram_probability))
def plot_oscillogram(
    costhz: Sequence[float],
    log10_energy: Sequence[float],
    probability: Optional[Sequence[Sequence[float]]] = None,
    *,
    nu_i: Optional[int] = None,
    nu_f: Optional[int] = None,
    nubar: Optional[bool] = None,
    num_flavors: Optional[int] = None,
    osc_params: Optional[Dict[str, float]] = None,
    electron_fraction: Optional[float] = None,
    electron_fraction_core: Optional[float] = None,
    electron_fraction_mantle: Optional[float] = None,
    electron_fraction_crust: Optional[float] = None,
    electron_fraction_ocean: Optional[float] = None,
    ratio_number_neutrons_to_protons: Optional[float] = None,
    wrapper_kw: Optional[Dict[str, Any]] = None,
    return_probability: bool = False,
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

    Given a ``probability`` array, it draws it.  Given none, it computes it first
    through the Earth wrappers (:func:`magnus.oscprob.osc_prob_3nu_earth` and its
    two-, four- and five-flavor siblings), one call per zenith angle over the whole
    energy array, so each is an energy scan the energy-batched engine answers.  The
    wrappers declare the PREM layer boundaries as slab edges and use a layered
    electron fraction unless it is overridden here.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.1
       ``probability`` is optional: without it, the oscillogram is computed through
       the Earth wrappers.  Adds ``num_flavors``, ``osc_params``, the
       electron-fraction keywords, ``wrapper_kw`` and ``return_probability``.

    .. versionchanged:: 1.2.0
       Takes nubar and labels an antineutrino oscillogram as such (issue #145), computes NSI and
       LIV through the Earth wrappers (issue #146 §1), and refuses log10_energy above 19 in
       compute mode (issue #160 §12).

    Parameters
    ----------
    costhz : sequence of float
        Zenith-angle cosines, the abscissa.  The physical range is
        :math:`[-1, 0]`; at :math:`\cos\theta_z \geq 0` the path inside the Earth
        has zero length and the computed probability is the identity.
    log10_energy : sequence of float
        :math:`\log_{10}` of the energy in GeV, the ordinate.
    probability : sequence of sequence of float, optional
        Probability with shape ``(len(log10_energy), len(costhz))``.  If None
        (default), it is computed through the Earth wrappers, which needs
        ``nu_i``, ``nu_f`` and ``num_flavors``.  Give either this or the
        arguments that compute it, not both.
    nu_i, nu_f : int, optional
        Flavor pair, used for the color-bar label and the annotation when
        those are not given explicitly, and as the channel when the
        probability is computed.  For two flavors they index the wrapper's
        :math:`2 \times 2` matrix (0 or 1).
    num_flavors : int, optional
        2, 3, 4 or 5: which Earth wrapper computes the probability.
    nubar : bool, optional
        True for antineutrinos.  In compute mode it is also passed to the wrapper, so the
        grid and its label agree; a value contradicting ``wrapper_kw['nubar']`` is refused.
        Default: None, which reads ``wrapper_kw['nubar']`` in compute mode and False otherwise.
    osc_params : dict, optional
        Mixing parameters, passed to the wrapper by name (``sth`` and ``Dm2``
        at two flavors, where they are required; ``s12``, ..., ``D31`` and the
        sterile ones at three to five flavors, where the wrapper's default set
        applies to any left out).
    electron_fraction, electron_fraction_core, electron_fraction_mantle, electron_fraction_crust, electron_fraction_ocean : float, optional
        Electron fraction, everywhere or per PREM layer.  None (default) keeps
        the wrappers' layered values.  Passed to the wrapper unchanged, with
        the wrapper's meaning and validation.
    ratio_number_neutrons_to_protons : float, optional
        Passed to the wrapper unchanged; None keeps its default.
    wrapper_kw : dict, optional
        Any other keyword the Earth wrapper accepts (``rtol``, ``atol``,
        ``nubar``, ``integration_method``, ``magnus_exp_order``, ...).  The
        channel, ``energy``, ``costhz``, ``L`` and the composition keywords
        above cannot be given here too.  ``loc_ini``, ``loc_fin`` and the
        depths are not refused but do not belong here: each column of the
        oscillogram is one zenith angle, from surface to surface.
    return_probability : bool, optional
        If True, also return the probability drawn.  Default is False.
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
        Tick spacings.  Pass None to hand an axis back to Matplotlib's own locator.
        Default: 0.2, 0.02, 0.1 and 0.02.
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
    probability : np.ndarray
        Only with ``return_probability=True``: the probability drawn, shape
        ``(len(log10_energy), len(costhz))``.

    Raises
    ------
    ValueError
        If the probability has the wrong shape; if it is given together with
        any argument that would compute it; if it is to be computed without
        ``nu_i``, ``nu_f`` or a valid ``num_flavors``, or at two flavors
        without ``osc_params``; or if a keyword arrives both explicitly and in
        ``wrapper_kw``.

    Examples
    --------
    A precomputed probability is drawn as given:

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

    Without one, it is computed through the Earth wrappers, here at three flavors
    with the default mixing parameters and the layered electron fraction:

    .. jupyter-execute::

        import matplotlib
        matplotlib.use('Agg')
        import numpy as np
        import magnus.globaldefs as gd
        from magnus.plotting import plot_oscillogram

        import warnings
        from magnus.magnus import MagnusConvergenceWarning

        c = np.linspace(-1.0, -0.1, 25)
        lE = np.linspace(0.0, 1.0, 20)
        # A few of these chords need the adaptive refinement to widen its first grids;
        # this is the expected, informational MagnusConvergenceWarning discussed in the
        # package README.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', MagnusConvergenceWarning)
            fig, ax, P = plot_oscillogram(c, lE, nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3,
                                          return_probability=True)
        print(P.shape, round(float(P.max()), 3))
    """
    mpl, plt = _mpl()
    import matplotlib.patheffects as path_effects

    computing = dict(num_flavors=num_flavors, osc_params=osc_params,
                     electron_fraction=electron_fraction,
                     electron_fraction_core=electron_fraction_core,
                     electron_fraction_mantle=electron_fraction_mantle,
                     electron_fraction_crust=electron_fraction_crust,
                     electron_fraction_ocean=electron_fraction_ocean,
                     ratio_number_neutrons_to_protons=ratio_number_neutrons_to_protons,
                     wrapper_kw=wrapper_kw)
    nubar_label, wrapper_kw = _resolve_nubar('plot_oscillogram', nubar, wrapper_kw,
                                             probability is None)
    computing['wrapper_kw'] = wrapper_kw
    if probability is None:
        # 10^19 GeV is the Planck scale: above it no oscillation wrapper is meaningful, and a
        # grid there was computed as if it were (issue #160 §12).
        _lg = np.asarray(log10_energy, dtype=float)
        if np.max(_lg) > 19.0:
            raise ValueError(_v._msg('plotting.plot_oscillogram', "log10_energy is log10 of the "
                                     "energy in GeV, so it must stay at or below 19, the Planck "
                                     "scale; its largest entry is " + format(float(np.max(_lg)),
                                                                              'g') + "."))
        probability = _oscillogram_through_earth_wrappers(costhz, log10_energy, nu_i, nu_f,
                                                          **computing)
    else:
        given = [name for name, value in computing.items() if value is not None]
        if given:
            raise ValueError(
                'Error in magnus: plotting.plot_oscillogram: give either a probability to draw '
                'or the arguments that compute it, not both; got a probability and '
                + ', '.join(given) + '.')

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

    # A grid computed with nubar=True is the antineutrino probability, and is labelled as one
    # (issue #145 §1): the labels used to say nu whatever the grid held.
    label = cbar_label
    if label is None and nu_i is not None and nu_f is not None:
        label = cbar_label_prefix + prob_label(nu_i, nu_f, nubar=nubar_label)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.ax.tick_params(labelsize=cbar_labelsize)
    if label is not None:
        cbar.set_label(label=label, fontsize=cbar_fontsize)

    if annotation is None and nu_i is not None and nu_f is not None:
        annotation = prob_label(nu_i, nu_f, nubar=nubar_label)
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
    if return_probability:
        return fig, ax, prob
    return fig, ax


_EARTH_COMPOSITION_KEYS = ('electron_fraction', 'electron_fraction_core',
                           'electron_fraction_mantle', 'electron_fraction_crust',
                           'electron_fraction_ocean', 'ratio_number_neutrons_to_protons')


def _oscillogram_through_earth_wrappers(costhz, log10_energy, nu_i, nu_f, num_flavors,
                                        osc_params, wrapper_kw, **composition):
    r"""Probability of shape ``(len(log10_energy), len(costhz))`` from the Earth wrappers.

    One call per zenith angle over the whole energy array: an energy scan at one
    baseline, which the energy-batched engine answers.  The wrapper declares the
    PREM layer boundaries itself; the composition keywords left at None are not
    passed, so the wrapper's own defaults (a layered electron fraction) apply.
    """
    fn, kwargs = _earth_wrapper_and_arguments(
        'plot_oscillogram', nu_i, nu_f, num_flavors, osc_params, wrapper_kw, composition,
        reserved={'energy', 'costhz', 'L'})

    from magnus import earth, globaldefs as gd
    energy = 10.0**np.asarray(log10_energy, dtype=float)*gd.UNIT_GEV
    probability = np.empty((len(energy), len(costhz)))
    for j, cz in enumerate(costhz):
        L = earth.distance_traveled_inside_earth(float(cz))*gd.UNIT_KM
        probability[:, j] = np.asarray(fn(energy, costhz=float(cz), L=L, **kwargs),
                                       dtype=float).reshape(len(energy))
    return probability


def _earth_wrapper_and_arguments(caller, nu_i, nu_f, num_flavors, osc_params, wrapper_kw,
                                 composition, reserved, flavors=(2, 3, 4, 5)):
    r"""The Earth wrapper for ``num_flavors`` and the keywords every call to it shares.

    The checks the plotting functions that compute through the Earth wrappers have in
    common: a channel, a flavor count the wrapper exists for, the two-flavor parameters
    that have no default, and no keyword arriving twice -- ``reserved`` names what the
    caller sets per call, which ``wrapper_kw`` and ``osc_params`` may not also set.  The
    composition keywords left at None are dropped, so the wrapper's own defaults (a
    layered electron fraction) apply.

    .. versionchanged:: 1.2.0
       Routes NSI and LIV parameters to the Earth wrapper that takes them (issue #146 §1), and
       the wrapper it returns names the plotting routine in its errors (issue #160 §12).
    """
    where = 'Error in magnus: plotting.%s: ' % caller
    if nu_i is None or nu_f is None:
        raise ValueError(where + 'computing the probability needs nu_i and nu_f.')
    if num_flavors not in flavors:
        allowed = ', '.join(str(n) for n in flavors[:-1]) + ' or %d' % flavors[-1]
        raise ValueError(where + 'computing the probability needs num_flavors = %s, not %r.'
                         % (allowed, num_flavors))
    osc_params = dict(osc_params or {})
    if num_flavors == 2 and not {'sth', 'Dm2'} <= set(osc_params):
        raise ValueError(where + 'at two flavors, osc_params must give sth and Dm2.')
    wrapper_kw = dict(wrapper_kw or {})
    reserved = set(reserved) | {'nu_i', 'nu_f'} | set(_EARTH_COMPOSITION_KEYS)
    clash = sorted(reserved & (set(wrapper_kw) | set(osc_params)))
    if clash:
        raise ValueError(where + 'these belong to %s itself and cannot be given in '
                         'wrapper_kw or osc_params: ' % caller + ', '.join(clash) + '.')
    composition = {k: v for k, v in composition.items() if v is not None}

    # NSI and LIV through the Earth (issue #146 §1): the wrapper follows the parameters given.
    # Before, the standard wrapper was always chosen, and refused the eps_* or b* keys.
    import re
    given = set(wrapper_kw) | set(osc_params)
    nsi = sorted(k for k in given if k.startswith('eps_'))
    liv = sorted(k for k in given if re.fullmatch(r'b\d|sxi\d*|dxi\w*|Lambda|n_liv', k))
    if nsi and liv:
        raise ValueError(where + 'NSI (%s) and LIV (%s) parameters together: no Earth wrapper '
                         'takes both.  Compute the grid yourself and pass it as probability.'
                         % (', '.join(nsi), ', '.join(liv)))
    suffix = '_nsi' if nsi else ('_liv' if liv else '')

    from magnus import oscprob
    fn = getattr(oscprob, 'osc_prob_%dnu_earth%s' % (num_flavors, suffix), None)
    if fn is None:
        raise ValueError(where + 'there is no %d-flavor Earth wrapper for %s parameters.  '
                         'Compute the grid yourself and pass it as probability.'
                         % (num_flavors, 'NSI' if nsi else 'LIV'))
    return _named_by(caller, fn), dict(nu_i=nu_i, nu_f=nu_f, **osc_params, **composition,
                                       **wrapper_kw)


def _named_by(caller: str, fn):
    r"""``fn``, with its argument errors prefixed by the plotting routine the caller called.

    In compute mode the physics arguments reach an ``oscprob`` wrapper, and its refusals named
    only that wrapper -- a function the caller never called (issue #160 §12).  The wrapper's
    own message is kept whole after the plotting routine's name.

    .. versionadded:: 1.2.0
    """
    import functools

    @functools.wraps(fn)
    def named(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (TypeError, ValueError) as error:
            text = str(error)
            if 'Error in magnus' not in text:
                raise
            body = text.split('Error in magnus: ', 1)[1]
            raise type(error)('Error in magnus: plotting.' + caller + ': computing the '
                              'probability with ' + body) from error
    return named
