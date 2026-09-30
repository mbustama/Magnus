# -*- coding: utf-8 -*-
"""Tests of the pre-packaged figures (magnus.plotting).

These assert on the objects the functions return -- line counts, labels, limits,
tick locators -- rather than on rendered pixels, which are fragile and say
little. The recurring question each test answers is whether the house style
actually survived the move into the module: a figure switched over from a
hand-built block must come out the same, so the defaults are pinned here
against the values the notebooks used.

Every test closes its figures. Matplotlib is forced onto the Agg backend at
import, before pyplot is first imported anywhere, so nothing tries to reach a
display.
"""

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import magnus.globaldefs as gd  # noqa: E402
import magnus.plotting as mp  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures():
    """Close every figure a test leaves behind, so the suite cannot leak them."""
    yield
    plt.close('all')


@pytest.fixture
def sample():
    """A baseline grid, a probability, and a slightly perturbed reference."""
    L = np.logspace(1.0, 5.0, 128)
    exact = np.sin(L / 3000.0) ** 2
    approx = exact + 1.0e-6 * np.cos(L / 700.0)
    return L, exact, approx


# ----------------------------------------------------------------------
# prob_label
# ----------------------------------------------------------------------

@pytest.mark.parametrize('nu_i, nu_f, expected', [
    (gd.NUE, gd.NUE, r'$P_{\nu_e \to \nu_e}$'),
    (gd.NUMU, gd.NUE, r'$P_{\nu_\mu \to \nu_e}$'),
    (gd.NUTAU, gd.NUMU, r'$P_{\nu_\tau \to \nu_\mu}$'),
    (gd.NUE, gd.NUS1, r'$P_{\nu_e \to \nu_{s_1}}$'),
    (gd.NUE, gd.NUS2, r'$P_{\nu_e \to \nu_{s_2}}$'),
])
def test_prob_label_covers_active_and_sterile_flavors(nu_i, nu_f, expected):
    """The notebooks' hand-written chains only covered the three active flavours."""
    assert mp.prob_label(nu_i, nu_f) == expected


def test_prob_label_bars_the_nu_not_the_subscript():
    """\\bar{\\nu}_\\mu, not \\bar{\\nu_\\mu}: the bar belongs over the nu alone."""
    assert mp.prob_label(gd.NUMU, gd.NUE, nubar=True) == \
        r'$P_{\bar{\nu}_\mu \to \bar{\nu}_e}$'


def test_prob_label_rejects_an_unknown_flavor():
    """A bad index must name the offending parameter, per the package convention."""
    with pytest.raises(ValueError, match='nu_f'):
        mp.prob_label(gd.NUE, 99)
    with pytest.raises(ValueError, match='nu_i'):
        mp.prob_label(-1, gd.NUE)


# ----------------------------------------------------------------------
# plot_curves -- the workhorse
# ----------------------------------------------------------------------

def test_plot_curves_returns_one_axes_without_a_residual(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [exact])
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    assert len(ax.get_lines()) == 1


def test_plot_curves_adds_a_residual_panel_and_mutes_the_main_ticks(sample):
    """With a residual, the main panel loses its tick labels -- the notebooks'
    ax[0].xaxis.set_ticklabels([]) -- so the two panels read as one figure."""
    L, exact, approx = sample
    fig, ax = mp.plot_curves(
        L, [dict(y=approx, label='Magnus expansion'),
            dict(y=exact, label='Standard formula', color='k', ls='--')],
        residual=(approx - exact) / exact, residual_label=r'$\epsilon_{\rm rel}$')
    assert len(ax) == 2
    assert len(ax[0].get_lines()) == 2
    assert len(ax[1].get_lines()) == 1
    assert all(t.get_text() == '' for t in ax[0].get_xticklabels())
    assert ax[1].get_ylabel() == r'$\epsilon_{\rm rel}$'


def test_plot_curves_applies_labels_limits_and_scales(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(
        L, [exact], xlabel='X', ylabel='Y', title='T',
        xlim=(10.0, 1.0e5), ylim=(0.0, 1.0), xscale='log')
    assert (ax.get_xlabel(), ax.get_ylabel(), ax.get_title()) == ('X', 'Y', 'T')
    assert ax.get_xlim() == (10.0, 1.0e5)
    assert ax.get_ylim() == (0.0, 1.0)
    assert ax.get_xscale() == 'log'


def test_plot_curves_shares_xlim_and_xscale_with_the_residual_panel(sample):
    """Misaligned panels are the classic failure of a hand-built two-panel figure."""
    L, exact, approx = sample
    fig, ax = mp.plot_curves(L, [exact], residual=approx - exact,
                             xlim=(10.0, 1.0e5), xscale='log')
    assert ax[0].get_xlim() == ax[1].get_xlim()
    assert ax[0].get_xscale() == ax[1].get_xscale() == 'log'


def test_plot_curves_xlabel_goes_under_the_bottom_panel(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_curves(L, [exact], residual=approx - exact, xlabel='L')
    assert ax[1].get_xlabel() == 'L'
    assert ax[0].get_xlabel() == ''


def test_plot_curves_assigns_the_default_color_cycle_in_order(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [exact, exact + 0.1, exact + 0.2])
    assert [ln.get_color() for ln in ax.get_lines()] == ['C0', 'C1', 'C2']


def test_plot_curves_honours_an_explicit_color(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [dict(y=exact, color='k', ls='--', lw=3)])
    line, = ax.get_lines()
    assert (line.get_color(), line.get_linestyle(), line.get_linewidth()) == \
        ('k', '--', 3)


def test_plot_curves_draws_a_legend_only_when_a_curve_is_labelled(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [exact])
    assert ax.get_legend() is None
    fig, ax = mp.plot_curves(L, [dict(y=exact, label='Magnus expansion')])
    assert ax.get_legend() is not None


def test_plot_curves_legend_defaults_match_the_house_style(sample):
    """These exact values are repeated on essentially every notebook figure."""
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [dict(y=exact, label='m')],
                             legend_title='Calculation method')
    leg = ax.get_legend()
    assert leg.get_title().get_text() == 'Calculation method'
    assert leg.get_frame_on()
    assert HOUSE_LEGEND_FONTSIZE == mp.HOUSE_LEGEND_KW['fontsize'] == 17
    assert mp.HOUSE_LEGEND_KW['title_fontsize'] == 20
    assert mp.HOUSE_LEGEND_KW['edgecolor'] == 'k'


HOUSE_LEGEND_FONTSIZE = 17


def test_plot_curves_legend_kw_overrides_the_house_default(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [dict(y=exact, label='m')],
                             legend_kw=dict(ncol=3))
    assert ax.get_legend()._ncols == 3


def test_plot_curves_sets_multiple_locators_from_the_tick_spacings(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [exact], ymajor=0.10, yminor=0.02)
    major = ax.yaxis.get_major_locator()
    minor = ax.yaxis.get_minor_locator()
    assert isinstance(major, matplotlib.ticker.MultipleLocator)
    assert isinstance(minor, matplotlib.ticker.MultipleLocator)
    # MultipleLocator keeps its spacing on a private _edge; the supported way to
    # read it back is the spacing between consecutive ticks it produces.
    assert np.diff(major.tick_values(0.0, 1.0))[0] == pytest.approx(0.10)
    assert np.diff(minor.tick_values(0.0, 1.0))[0] == pytest.approx(0.02)


def test_plot_curves_default_figsize_is_the_house_one(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [exact])
    assert tuple(fig.get_size_inches()) == mp.HOUSE_FIGSIZE == (18.0, 9.0)


def test_plot_curves_rejects_a_curve_dict_without_an_ordinate(sample):
    """A silent skip here would drop a curve from the figure without a word."""
    L, _, _ = sample
    with pytest.raises(ValueError, match="'y'"):
        mp.plot_curves(L, [dict(label='no ordinate')])


def test_plot_curves_rejects_an_unknown_keyword(sample):
    """No bare **kwargs anywhere: a typo must fail at the call site."""
    L, exact, _ = sample
    with pytest.raises(TypeError):
        mp.plot_curves(L, [exact], ylabl='typo')


def test_plot_curves_forwards_a_typo_in_a_curve_dict_to_matplotlib(sample):
    """Curve keywords land in Axes.plot, so a bad one raises there by name."""
    L, exact, _ = sample
    with pytest.raises(AttributeError):
        mp.plot_curves(L, [dict(y=exact, colour='k')])


def test_plot_curves_annotates_the_main_panel(sample):
    """The BSM notebooks record the parameter values on the figure itself."""
    L, exact, approx = sample
    fig, ax = mp.plot_curves(
        L, [exact], residual=approx - exact,
        annotations=[dict(text=r'$\epsilon_{ee} = 0.06$', xy=(0.02, 0.03)),
                     dict(text='second', xy=(0.02, 0.10), fontsize=18)])
    assert [t.get_text() for t in ax[0].texts] == \
        [r'$\epsilon_{ee} = 0.06$', 'second']
    # they go on the main panel, not the residual one
    assert len(ax[1].texts) == 0


def test_plot_curves_saves_when_asked(sample, tmp_path):
    L, exact, _ = sample
    out = tmp_path / 'fig.pdf'
    mp.plot_curves(L, [exact], savefig=str(out))
    assert out.exists() and out.stat().st_size > 0


def test_plot_curves_places_the_legend_where_asked(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_curves(L, [dict(y=exact, label='m')],
                             legend_loc='lower left')
    assert ax.get_legend()._loc == matplotlib.legend.Legend.codes['lower left']


def test_plot_curves_sets_the_residual_ylim_and_ticks(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_curves(
        L, [exact], residual=(approx - exact) / exact,
        residual_ylim=(-0.5, 0.5), residual_ymajor=0.20, residual_yminor=0.05)
    assert ax[1].get_ylim() == (-0.5, 0.5)
    step = np.diff(ax[1].yaxis.get_major_locator().tick_values(-0.5, 0.5))[0]
    assert step == pytest.approx(0.20)


def test_plot_curves_draws_a_grid_on_both_panels(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_curves(L, [exact], residual=approx - exact, grid=True)
    assert ax[0].xaxis.get_gridlines()[0].get_visible()
    assert ax[1].xaxis.get_gridlines()[0].get_visible()


def test_plot_curves_applies_x_tick_spacings_to_both_panels(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_curves(L, [exact], residual=approx - exact,
                             xmajor=10000.0, xminor=2000.0)
    for axx in ax:
        step = np.diff(axx.xaxis.get_major_locator().tick_values(0.0, 1.0e5))[0]
        assert step == pytest.approx(10000.0)


# ----------------------------------------------------------------------
# plot_curves_stacked: small multiples
#
# The point of this layout is that the reader compares panels against each
# other, so the tests below are mostly about what every panel has in common --
# limits, scales, tick spacings -- and about which single panel gets the parts
# that must not repeat: the abscissa labels, the title, the legend.
# ----------------------------------------------------------------------

@pytest.fixture
def stack():
    """An abscissa and three panels of two curves each."""
    E = np.linspace(1.0, 40.0, 64)
    panels = [[np.sin(k*E/8.0)**2, 0.8*np.sin(k*E/8.0)**2]
              for k in (0.5, 1.0, 2.0)]
    return E, panels


def test_stacked_returns_one_axes_per_panel(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels)
    assert ax.shape == (3,)
    assert all(len(axx.get_lines()) == 2 for axx in ax)


def test_stacked_returns_an_array_even_for_one_panel(stack):
    """A one-panel stack must index like any other, so callers that loop or
    subscript do not need a special case for it."""
    E, _ = stack
    fig, ax = mp.plot_curves_stacked(E, [[np.sin(E)]])
    assert ax.shape == (1,)


def test_stacked_mutes_every_panel_but_the_bottom(stack):
    """Repeating the tick labels on every panel is what makes a hand-built
    stack read as separate figures rather than one."""
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels)
    assert all(t.get_text() == '' for axx in ax[:-1]
               for t in axx.get_xticklabels())
    assert any(t.get_text() != '' for t in ax[-1].get_xticklabels())


def test_stacked_shares_limits_scales_and_ticks_across_panels(stack):
    """The comparison is between panels, so any axis that differs is a bug."""
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(
        E, panels, xlim=(1.0, 40.0), ylim=(0.0, 1.0), xmajor=10.0, ymajor=0.1)
    assert len({axx.get_xlim() for axx in ax}) == 1
    assert len({axx.get_ylim() for axx in ax}) == 1
    for axx in ax:
        assert np.diff(axx.xaxis.get_major_locator().tick_values(
            0.0, 40.0))[0] == pytest.approx(10.0)


def test_stacked_puts_the_xlabel_and_title_on_one_panel_each(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels, xlabel='X', title='T')
    assert [axx.get_xlabel() for axx in ax] == ['', '', 'X']
    assert [axx.get_title() for axx in ax] == ['T', '', '']


def test_stacked_ylabel_is_one_figure_level_label(stack):
    """Every panel shows the same quantity, so the ordinate is labelled once
    for the stack -- the notebooks did this by adding a frameless full-figure
    subplot purely to hang a label on."""
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels, ylabel='P')
    assert fig.get_supylabel() == 'P'
    assert all(axx.get_ylabel() == '' for axx in ax)


def test_stacked_shared_ylabel_matches_the_axis_label_size(stack):
    """supylabel takes rcParams['figure.labelsize'], every axis label takes
    rcParams['axes.labelsize'], and the notebooks set the latter to 25. Left on
    the default the shared ordinate label renders visibly smaller than the
    abscissa label under it -- caught by looking at the rendered figure, not by
    the code running."""
    E, panels = stack
    with matplotlib.rc_context({'axes.labelsize': 25, 'figure.labelsize': 11}):
        fig, ax = mp.plot_curves_stacked(E, panels, ylabel='P', xlabel='X')
    assert fig._supylabel.get_fontsize() == 25


def test_stacked_shared_ylabel_size_can_be_overridden(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels, ylabel='P',
                                     ylabel_kw=dict(fontsize=8))
    assert fig._supylabel.get_fontsize() == 8


def test_stacked_colour_cycle_restarts_in_each_panel(stack):
    """The n-th curve of every panel must match, or the panels cannot be read
    against each other."""
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels)
    first = [axx.get_lines()[0].get_color() for axx in ax]
    second = [axx.get_lines()[1].get_color() for axx in ax]
    assert set(first) == {'C0'}
    assert set(second) == {'C1'}


def test_stacked_labels_each_panel(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels,
                                     panel_labels=['a', 'b', 'c'])
    assert [axx.texts[0].get_text() for axx in ax] == ['a', 'b', 'c']


def test_stacked_rejects_a_panel_label_count_mismatch(stack):
    E, panels = stack
    with pytest.raises(ValueError, match='one label per panel'):
        mp.plot_curves_stacked(E, panels, panel_labels=['a'])


def test_stacked_rejects_an_empty_stack():
    with pytest.raises(ValueError, match='at least'):
        mp.plot_curves_stacked(np.linspace(0.0, 1.0, 4), [])


def test_stacked_rejects_a_legend_panel_out_of_range(stack):
    E, panels = stack
    with pytest.raises(ValueError, match='legend_panel'):
        mp.plot_curves_stacked(E, panels, legend_panel=7)


def test_stacked_annotation_can_name_its_panel(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(
        E, panels, annotations=[dict(text='params', xy=(0.5, 0.5), panel=1)])
    assert ax[1].texts[0].get_text() == 'params'
    assert len(ax[0].texts) == 0


def test_stacked_legend_proxies_describe_a_style_not_a_curve(stack):
    """When colour varies panel to panel, the legend has to describe line style
    instead. The proxy handles replace the notebooks' trick of plotting dummy
    points outside the axis limits to manufacture legend entries -- so no curve
    carries a label here, and the legend is still correct."""
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(
        E, panels, legend_panel=0,
        legend_proxies=[dict(label='3+1', color='k', ls='-'),
                        dict(label='standard', color='k', ls='--')])
    leg = ax[0].get_legend()
    assert [t.get_text() for t in leg.get_texts()] == ['3+1', 'standard']
    assert all(axx.get_legend() is None for axx in ax[1:])
    assert all(len(axx.get_lines()) == 2 for axx in ax)


def test_stacked_legend_falls_back_to_curve_labels(stack):
    E, panels = stack
    labelled = [[dict(y=panels[0][0], label='first'), panels[0][1]]] + panels[1:]
    fig, ax = mp.plot_curves_stacked(E, labelled)
    assert [t.get_text() for t in ax[0].get_legend().get_texts()] == ['first']


def test_stacked_draws_no_legend_when_there_is_nothing_to_say(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels)
    assert all(axx.get_legend() is None for axx in ax)


def test_stacked_legend_defaults_match_the_house_style(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(
        E, panels, legend_proxies=[dict(label='a', color='k')])
    leg = ax[0].get_legend()
    assert leg.get_frame_on() is mp.HOUSE_LEGEND_KW['frameon']
    assert leg.get_texts()[0].get_fontsize() == mp.HOUSE_LEGEND_KW['fontsize']


def test_stacked_default_figsize_scales_with_the_panel_count(stack):
    E, panels = stack
    fig, ax = mp.plot_curves_stacked(E, panels)
    assert tuple(fig.get_size_inches()) == (
        mp.HOUSE_FIGSIZE[0], 0.5*mp.HOUSE_FIGSIZE[1]*3)


def test_stacked_rejects_an_unknown_keyword(stack):
    E, panels = stack
    with pytest.raises(TypeError):
        mp.plot_curves_stacked(E, panels, colour='red')


def test_stacked_saves_when_asked(stack, tmp_path):
    E, panels = stack
    out = tmp_path / 'stack.pdf'
    mp.plot_curves_stacked(E, panels, savefig=str(out))
    assert out.exists() and out.stat().st_size > 0


# ----------------------------------------------------------------------
# the probability presets
# ----------------------------------------------------------------------

def test_vs_baseline_presets_log_axis_unit_range_and_label(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_vs_baseline(
        L, [exact], nu_i=gd.NUE, nu_f=gd.NUE, num_flavors=2)
    assert ax.get_xscale() == 'log'
    assert ax.get_ylim() == (0.0, 1.0)
    assert ax.get_xlabel() == r'Baseline, $L$ [km]'
    assert ax.get_ylabel() == \
        r'Two-neutrino probability, $P_{\nu_e \to \nu_e}$'


def test_vs_energy_builds_its_abscissa_label_from_the_unit(sample):
    _, exact, _ = sample
    E = np.logspace(-1.0, 1.0, len(exact))
    fig, ax = mp.plot_probability_vs_energy(E, [exact], energy_unit='MeV')
    assert ax.get_xlabel() == r'Neutrino energy, $E_\nu$ [MeV]'


def test_vs_baseline_without_a_flavor_pair_leaves_the_ylabel_unset(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_vs_baseline(L, [exact])
    assert ax.get_ylabel() == ''


def test_probability_ylabel_rejects_an_impossible_flavor_count(sample):
    L, exact, _ = sample
    with pytest.raises(ValueError, match='num_flavors'):
        mp.plot_probability_vs_baseline(L, [exact], nu_i=gd.NUE, nu_f=gd.NUE,
                                        num_flavors=7)


def test_probability_ylabel_without_a_flavor_count_is_generic(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_vs_baseline(L, [exact], nu_i=gd.NUE,
                                              nu_f=gd.NUMU)
    assert ax.get_ylabel() == r'Probability, $P_{\nu_e \to \nu_\mu}$'


def test_vs_energy_builds_its_ylabel_from_the_flavor_pair(sample):
    _, exact, _ = sample
    E = np.logspace(-1.0, 1.0, len(exact))
    fig, ax = mp.plot_probability_vs_energy(E, [exact], nu_i=gd.NUMU,
                                           nu_f=gd.NUTAU, num_flavors=3)
    assert ax.get_ylabel() == \
        r'Three-neutrino probability, $P_{\nu_\mu \to \nu_\tau}$'


def test_no_default_label_has_a_tilde_outside_math():
    """Issue #124: ``~`` is a space only under usetex; under mathtext it prints.

    Every string literal in the module that is not a docstring is checked, with
    its ``$...$`` segments removed, since inside math ``~`` is a space in both.
    """
    import ast
    import inspect
    import re
    tree = ast.parse(inspect.getsource(mp))
    docstrings = {id(n.value) for n in ast.walk(tree)
                  if isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant)
                  and isinstance(n.value.value, str)}
    offenders = [(n.lineno, n.value) for n in ast.walk(tree)
                 if isinstance(n, ast.Constant) and isinstance(n.value, str)
                 and id(n) not in docstrings
                 and '~' in re.sub(r'\$[^$]*\$', '', n.value)]
    assert offenders == []


def test_presets_forward_the_rest_to_plot_curves(sample):
    """The wrappers must not quietly drop plot_curves' arguments."""
    L, exact, approx = sample
    fig, ax = mp.plot_probability_vs_baseline(
        L, [exact], residual=approx - exact, title='T', grid=True)
    assert len(ax) == 2
    assert ax[0].get_title() == 'T'


# ----------------------------------------------------------------------
# profile + probability panels
# ----------------------------------------------------------------------

def test_profile_figure_has_one_density_panel_plus_one_per_probability(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=np.exp(-L / 5.0e3))],
        [[dict(y=exact, label='a')], [dict(y=approx, label='b')]])
    assert len(ax) == 3
    assert len(ax[0].get_lines()) == 1


def test_profile_figure_labels_only_the_bottom_abscissa(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=np.exp(-L / 5.0e3))], [[dict(y=exact)]], xlabel='L')
    assert ax[-1].get_xlabel() == 'L'
    assert all(t.get_text() == '' for t in ax[0].get_xticklabels())


def test_profile_figure_lets_a_curve_carry_its_own_abscissa(sample):
    """The long-baseline notebook gives each detector a different baseline grid."""
    L, exact, _ = sample
    other = np.logspace(1.0, 4.0, 64)
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=np.exp(-L / 5.0e3))],
        [[dict(x=other, y=np.sin(other / 900.0) ** 2)]])
    line, = ax[1].get_lines()
    assert len(line.get_xdata()) == 64


def test_profile_figure_omits_the_density_panel_when_asked(sample):
    """The long-baseline notebook stacks a probability over its smoothed
    version, with no density panel at all."""
    L, exact, approx = sample
    for profiles in (None, []):
        fig, ax = mp.plot_probability_with_profile(
            L, profiles, [[dict(y=exact, label='raw')], [dict(y=approx)]],
            xlabel='L')
        assert len(ax) == 2
        assert len(ax[0].get_lines()) == 1
        assert ax[1].get_xlabel() == 'L'
        assert all(t.get_text() == '' for t in ax[0].get_xticklabels())
        plt.close(fig)


def test_profile_figure_titles_the_top_panel_without_a_profile(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=exact)]], title='T')
    assert ax[0].get_title() == 'T'


def test_profile_figure_annotates_each_panel(sample):
    """The three-flavour notebook names the channel inside each panel rather
    than in its ordinate label."""
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=exact)], [dict(y=approx)]],
        panel_annotations=[r'$P_{ee}$', r'$P_{e\mu}$'])
    assert [t.get_text() for t in ax[0].texts] == [r'$P_{ee}$']
    assert [t.get_text() for t in ax[1].texts] == [r'$P_{e\mu}$']


def test_profile_figure_annotation_can_carry_its_own_style(sample):
    """Over dense curves the text needs a background to stay readable."""
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=exact)]],
        panel_annotations=[dict(text='n = 0', fontsize=12,
                                bbox=dict(facecolor='white', edgecolor='none'))])
    ann, = ax[0].texts
    assert ann.get_text() == 'n = 0'
    assert ann.get_fontsize() == 12
    assert ann.get_bbox_patch() is not None


def test_profile_figure_skips_a_none_annotation(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=exact)], [dict(y=approx)]],
        panel_annotations=[r'$P_{ee}$', None])
    assert len(ax[0].texts) == 1
    assert len(ax[1].texts) == 0


def test_profile_figure_draws_one_shared_ylabel_for_the_stack(sample):
    """Four panels of the same quantity get one label, not four."""
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=exact)], [dict(y=approx)]],
        shared_ylabel='Three-neutrino probability')
    # the overlay is an extra axes beyond the two panels
    assert len(fig.axes) == 3
    overlay = fig.axes[-1]
    assert overlay.get_ylabel() == 'Three-neutrino probability'
    assert not overlay.get_frame_on()
    # Its tick labels must still be present but invisible: they reserve the
    # width that keeps the shared label off the panels' own numbers. Removing
    # them (set_yticks([])) makes the label collide.
    assert overlay.get_yticks().size > 0
    assert all(t.get_color() == 'none' for t in overlay.get_yticklabels())


def test_profile_figure_panels_can_take_a_log_ordinate(sample):
    """The matrix-exponential notebook stacks convergence curves spanning
    decades, which are not probabilities."""
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=np.abs(exact) + 1e-3)], [dict(y=np.abs(exact) + 1e-3)]],
        panel_yscale='log', panel_ylim=(1e-3, 10.0),
        panel_ymajor=None, panel_yminor=None)
    assert ax[0].get_yscale() == 'log'
    assert ax[1].get_yscale() == 'log'


def test_profile_figure_panel_ylim_can_autoscale(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_profile(
        L, None, [[dict(y=exact * 100.0)]], panel_ylim=None)
    assert ax[0].get_ylim()[1] > 1.0


def test_profile_figure_applies_x_tick_spacings_to_every_panel(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=L)], [[dict(y=exact)], [dict(y=approx)]],
        xscale='linear', xmajor=10000.0, xminor=2000.0)
    for axx in ax:
        step = np.diff(axx.xaxis.get_major_locator().tick_values(0.0, 1.0e5))[0]
        assert step == pytest.approx(10000.0)


def test_profile_figure_requires_at_least_one_probability_panel(sample):
    L, _, _ = sample
    with pytest.raises(ValueError, match='panels'):
        mp.plot_probability_with_profile(L, [dict(y=L)], [])


def test_profile_figure_puts_the_legend_on_the_requested_panel(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=L)], [[dict(y=exact, label='a')], [dict(y=approx, label='b')]],
        legend_on_panel=1)
    assert ax[1].get_legend() is None
    assert ax[2].get_legend() is not None


def test_profile_figure_applies_title_limits_and_panel_labels(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=np.exp(-L / 5.0e3))],
        [[dict(y=exact)], [dict(y=approx)]],
        title=r'$3\nu$ inside the Earth', profile_ylim=(0.0, 6.0),
        profile_ymajor=2.0, profile_yminor=1.0, xlim=(10.0, 1.0e5),
        panel_ylabels=[r'$P_{ee}$', None])
    assert ax[0].get_title() == r'$3\nu$ inside the Earth'
    assert ax[0].get_ylim() == (0.0, 6.0)
    assert ax[0].get_xlim() == ax[1].get_xlim() == (10.0, 1.0e5)
    assert ax[1].get_ylabel() == r'$P_{ee}$'
    assert ax[2].get_ylabel() == ''


def test_profile_figure_places_its_legend_where_asked(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=L)], [[dict(y=exact, label='a')]],
        legend_title='Matter profile', legend_loc='lower left')
    leg = ax[1].get_legend()
    assert leg.get_title().get_text() == 'Matter profile'
    assert leg._loc == matplotlib.legend.Legend.codes['lower left']


def test_profile_figure_can_legend_every_panel(sample):
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_profile(
        L, [dict(y=L)], [[dict(y=exact, label='a')], [dict(y=approx, label='b')]],
        legend_on_panel=-1)
    assert ax[1].get_legend() is not None
    assert ax[2].get_legend() is not None


# ----------------------------------------------------------------------
# averaged overlay
# ----------------------------------------------------------------------

def test_average_overlay_broadcasts_a_scalar_average(sample):
    L, exact, _ = sample
    fig, ax = mp.plot_probability_with_average(L, exact, 0.5)
    lines = ax.get_lines()
    assert len(lines) == 2
    assert np.allclose(lines[1].get_ydata(), 0.5)
    assert len(lines[1].get_ydata()) == len(L)


def test_average_overlay_accepts_a_varying_average(sample):
    L, exact, _ = sample
    avg = np.linspace(0.3, 0.7, len(L))
    fig, ax = mp.plot_probability_with_average(L, [exact], [avg])
    assert np.allclose(ax.get_lines()[1].get_ydata(), avg)


def test_average_overlay_draws_one_dashed_line_per_channel(sample):
    """The averaged-probability notebook shows several channels at once, each
    with its own averaged value in its own colour."""
    L, exact, approx = sample
    P = [exact, approx, 1.0 - exact]
    fig, ax = mp.plot_probability_with_average(
        L, P, [0.5, 0.42, 0.31], labels=['Pee', 'Pem', 'Pet'])
    lines = ax.get_lines()
    assert len(lines) == 6
    for i in range(3):
        solid, dashed = lines[2*i], lines[2*i + 1]
        assert solid.get_color() == dashed.get_color() == f'C{i}'
        assert (solid.get_linestyle(), dashed.get_linestyle()) == ('-', '--')


def test_average_overlay_legend_names_channels_once_plus_the_style(sample):
    """One entry per channel, and a single entry explaining the dashes --
    rather than repeating 'averaged' once per curve."""
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_average(
        L, [exact, approx], [0.5, 0.4], labels=['Pee', 'Pem'])
    assert [t.get_text() for t in ax.get_legend().get_texts()] == \
        ['Pee', 'Pem', 'Phase-averaged']


def test_average_overlay_honours_legend_placement_and_colors(sample):
    """The legend is rebuilt here rather than by plot_curves, so its title,
    location and the per-channel colours have to survive that."""
    L, exact, approx = sample
    fig, ax = mp.plot_probability_with_average(
        L, [exact, approx], [0.5, 0.4], labels=['a', 'b'],
        colors=['C3', 'C4'], legend_title='Channel', legend_loc='upper left')
    leg = ax.get_legend()
    assert leg.get_title().get_text() == 'Channel'
    assert leg._loc == matplotlib.legend.Legend.codes['upper left']
    assert [ln.get_color() for ln in ax.get_lines()] == \
        ['C3', 'C3', 'C4', 'C4']


def test_average_overlay_rejects_mismatched_counts(sample):
    L, exact, approx = sample
    with pytest.raises(ValueError, match='one to one'):
        mp.plot_probability_with_average(L, [exact, approx], [0.5])


# ----------------------------------------------------------------------
# bi-probability
# ----------------------------------------------------------------------

@pytest.fixture
def biprob():
    d = np.linspace(-np.pi, np.pi, 100)
    return d, 0.05 + 0.02 * np.sin(d), 0.04 + 0.02 * np.sin(d + 0.4)


def test_biprobability_plots_one_locus_per_configuration(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability([p_nu, p_nu * 0.9], [p_nubar, p_nubar * 1.1])
    assert len(ax.get_lines()) == 2


def test_biprobability_defaults_to_the_appearance_channel_labels(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability([p_nu], [p_nubar])
    assert ax.get_xlabel() == r'$P_{\nu_\mu \to \nu_e}$'
    assert ax.get_ylabel() == r'$P_{\bar{\nu}_\mu \to \bar{\nu}_e}$'


def test_biprobability_rejects_mismatched_curve_counts(biprob):
    _, p_nu, p_nubar = biprob
    with pytest.raises(ValueError, match='same number of curves'):
        mp.plot_biprobability([p_nu, p_nu], [p_nubar])


def test_biprobability_marks_selected_phases(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability(
        [p_nu], [p_nubar],
        markers=[dict(index=0, marker='o', label=r'$-\pi$'),
                 dict(index=50, marker='s', label=r'$0$', filled=False)])
    # one scatter per marker per curve, plus one proxy per labelled marker
    assert len(ax.collections) == 4
    assert ax.get_legend() is not None


def test_biprobability_marker_proxies_are_off_the_data(biprob):
    """The notebooks parked proxies at (-10, -10), which is inside a rescaled
    axes; NaN keeps them out of the autoscale entirely."""
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability(
        [p_nu], [p_nubar], markers=[dict(index=0, marker='o', label='p')])
    xs = np.concatenate([c.get_offsets().data[:, 0] for c in ax.collections])
    assert np.isnan(xs).any()
    assert np.nanmin(xs) >= 0.0


def test_biprobability_applies_labels_curve_styles_and_limits(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability(
        [p_nu], [p_nubar], labels=['Normal ordering'],
        curve_kw=[dict(lw=3, ls=':')],
        title='T', xlim=(0.0, 0.1), ylim=(0.0, 0.1),
        legend_loc='upper left')
    line, = ax.get_lines()
    assert (line.get_linewidth(), line.get_linestyle()) == (3, ':')
    assert line.get_label() == 'Normal ordering'
    assert ax.get_title() == 'T'
    assert ax.get_xlim() == (0.0, 0.1) and ax.get_ylim() == (0.0, 0.1)
    assert ax.get_legend()._loc == matplotlib.legend.Legend.codes['upper left']


def test_biprobability_markers_accept_explicit_coordinates(biprob):
    """The marked phases are often computed separately from the curve, so they
    arrive as coordinates rather than as positions along it."""
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability(
        [p_nu], [p_nubar],
        markers=[dict(xy=(0.055, 0.045), marker='o', label='bf')])
    pts = np.concatenate([c.get_offsets().data for c in ax.collections])
    assert any(np.allclose(row, (0.055, 0.045)) for row in pts)


def test_biprobability_marker_needs_a_position(biprob):
    _, p_nu, p_nubar = biprob
    with pytest.raises(ValueError, match="'index'.*'xy'"):
        mp.plot_biprobability([p_nu], [p_nubar], markers=[dict(marker='o')])


def test_biprobability_markers_can_target_one_curve(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability(
        [p_nu, p_nu * 0.9], [p_nubar, p_nubar * 1.1],
        markers=[dict(index=0, marker='o', curve=1)])
    assert len(ax.collections) == 1


def test_biprobability_annotates(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability(
        [p_nu], [p_nubar],
        annotations=[dict(text='NO', xy=(0.1, 0.9), fontsize=20)])
    assert [t.get_text() for t in ax.texts] == ['NO']


# ----------------------------------------------------------------------
# oscillogram
# ----------------------------------------------------------------------

@pytest.fixture
def oscillogram():
    c = np.linspace(-1.0, 0.0, 24)
    lE = np.linspace(-1.0, 1.0, 18)
    P = np.sin(np.outer(10.0 ** lE, 1.0 + c)) ** 2
    return c, lE, P


def test_oscillogram_draws_a_filled_map_with_a_colorbar(oscillogram):
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(c, lE, P, nu_i=gd.NUMU, nu_f=gd.NUMU)
    assert len(ax.collections) >= 1
    # the colour bar lives on its own axes
    assert len(fig.axes) == 2
    assert mp.prob_label(gd.NUMU, gd.NUMU) in fig.axes[1].get_ylabel()


def test_oscillogram_checks_the_probability_orientation(oscillogram):
    """Transposing the grid is the easy mistake, and contourf's own error does
    not say which way round the array should have been."""
    c, lE, P = oscillogram
    with pytest.raises(ValueError, match='len\\(log10_energy\\), len\\(costhz\\)'):
        mp.plot_oscillogram(c, lE, P.T)


def test_oscillogram_draws_a_precomputed_difference(oscillogram):
    """A precomputed map need not be a probability: notebook 06 draws a difference of two,
    with its own colour-bar label, so values below 0 are accepted (#160)."""
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(c, lE, P - 0.5, cmap='RdBu_r', cbar_label='difference')
    assert len(ax.collections) >= 1
    with pytest.raises(ValueError, match='finite'):
        mp.plot_oscillogram(c, lE, np.full_like(P, np.nan))


def test_oscillogram_annotates_the_channel_over_a_white_stroke(oscillogram):
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(c, lE, P, nu_i=gd.NUE, nu_f=gd.NUMU)
    assert [t.get_text() for t in ax.texts] == [mp.prob_label(gd.NUE, gd.NUMU)]


def test_oscillogram_annotation_can_be_suppressed(oscillogram):
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(c, lE, P, nu_i=gd.NUE, nu_f=gd.NUMU,
                                  annotation='')
    assert len(ax.texts) == 0


def test_oscillogram_honours_explicit_limits_and_labels(oscillogram):
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(
        c, lE, P, xlim=(-0.9, -0.1), ylim=(-0.5, 0.5),
        cbar_label='custom', contourf_kw=dict(levels=8))
    assert ax.get_xlim() == (-0.9, -0.1)
    assert ax.get_ylim() == (-0.5, 0.5)
    assert fig.axes[1].get_ylabel() == 'custom'


def test_oscillogram_prefixes_the_colorbar_label(oscillogram):
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(c, lE, P, nu_i=gd.NUMU, nu_f=gd.NUMU,
                                  cbar_label_prefix='Average~')
    assert fig.axes[1].get_ylabel().startswith('Average~')


def test_oscillogram_defaults_its_limits_to_the_data_range(oscillogram):
    c, lE, P = oscillogram
    fig, ax = mp.plot_oscillogram(c, lE, P)
    assert ax.get_xlim() == (float(c.min()), float(c.max()))
    assert ax.get_ylim() == (float(lE.min()), float(lE.max()))


# ----------------------------------------------------------------------
# plot_oscillogram computing through the Earth wrappers
# ----------------------------------------------------------------------

EARTH_C = np.linspace(-1.0, -0.2, 5)
EARTH_LE = np.linspace(0.0, 1.0, 4)


def _by_hand(fn, nu_i, nu_f, **kw):
    """The loop plot_oscillogram is meant to run, written out: one wrapper call per zenith."""
    import magnus.earth as earth
    energy = 10.0**EARTH_LE*gd.UNIT_GEV
    cols = [np.asarray(fn(energy, costhz=float(cz),
                          L=earth.distance_traveled_inside_earth(float(cz))*gd.UNIT_KM,
                          nu_i=nu_i, nu_f=nu_f, **kw), dtype=float)
            for cz in EARTH_C]
    return np.column_stack(cols)


@pytest.mark.parametrize('composition', [
    {},                                                   # the wrappers' layered default
    {'electron_fraction': 0.5},
    {'electron_fraction_core': 0.45},
])
def test_oscillogram_computed_is_the_earth_wrapper_loop(composition):
    """Without a probability, the oscillogram is the Earth wrapper's answer per zenith angle,
    bit for bit, with the composition passed through unchanged."""
    import magnus.oscprob as op
    fig, ax, P = mp.plot_oscillogram(EARTH_C, EARTH_LE, nu_i=gd.NUMU, nu_f=gd.NUE,
                                     num_flavors=3, return_probability=True, **composition)
    assert P.shape == (len(EARTH_LE), len(EARTH_C))
    assert np.array_equal(P, _by_hand(op.osc_prob_3nu_earth, gd.NUMU, gd.NUE, **composition))


def test_the_electron_fraction_override_reaches_the_wrapper():
    """A uniform Y_e = 0.5 is a different Earth from the layered default, so the maps differ."""
    kw = dict(nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3, return_probability=True)
    *_, P_layered = mp.plot_oscillogram(EARTH_C, EARTH_LE, **kw)
    *_, P_uniform = mp.plot_oscillogram(EARTH_C, EARTH_LE, electron_fraction=0.5, **kw)
    assert not np.array_equal(P_layered, P_uniform)


def test_oscillogram_computed_at_two_flavors():
    import magnus.oscprob as op
    params = dict(sth=0.55, Dm2=2.5e-3)
    *_, P = mp.plot_oscillogram(EARTH_C, EARTH_LE, nu_i=0, nu_f=0, num_flavors=2,
                                osc_params=params, return_probability=True)
    assert np.array_equal(P, _by_hand(op.osc_prob_2nu_earth, 0, 0, **params))


def test_oscillogram_without_return_probability_returns_the_figure_only():
    out = mp.plot_oscillogram(EARTH_C, EARTH_LE, nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=3)
    assert len(out) == 2


@pytest.mark.parametrize('kw, message', [
    (dict(num_flavors=3), 'needs nu_i and nu_f'),
    (dict(nu_i=gd.NUMU, nu_f=gd.NUMU), 'num_flavors = 2, 3, 4 or 5'),
    (dict(nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=6), 'num_flavors = 2, 3, 4 or 5'),
    (dict(nu_i=0, nu_f=0, num_flavors=2), 'must give sth and Dm2'),
    (dict(nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=3, electron_fraction=0.5,
          wrapper_kw=dict(electron_fraction=0.4)), 'cannot be given in wrapper_kw'),
    (dict(nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=3, wrapper_kw=dict(costhz=-0.5)),
     'cannot be given in wrapper_kw'),
])
def test_oscillogram_computation_refuses_what_it_cannot_honour(kw, message):
    with pytest.raises(ValueError, match=message):
        mp.plot_oscillogram(EARTH_C, EARTH_LE, **kw)


def test_oscillogram_refuses_a_probability_and_the_arguments_that_compute_it(oscillogram):
    """Drawing one array while ignoring the parameters given with it would be a silent lie."""
    c, lE, P = oscillogram
    with pytest.raises(ValueError, match='not both'):
        mp.plot_oscillogram(c, lE, P, num_flavors=3)
    with pytest.raises(ValueError, match='not both'):
        mp.plot_oscillogram(c, lE, P, electron_fraction=0.5)


# ----------------------------------------------------------------------
# plot_probability_with_profile computing through the Earth wrappers
# ----------------------------------------------------------------------

PROFILE_L = np.linspace(100.0, 9000.0, 8)          # [km], inside the costhz=-0.9 chord
PROFILE_E = np.logspace(0.0, 1.0, 6)                # [GeV]


def _profile(**kw):
    kw.setdefault('nu_i', gd.NUMU)
    kw.setdefault('nu_f', gd.NUE)
    kw.setdefault('num_flavors', 3)
    return mp.plot_probability_with_profile(PROFILE_L, return_probability=True, **kw)


@pytest.mark.parametrize('composition', [{}, {'electron_fraction': 0.5}])
def test_profile_computed_over_baseline_is_the_earth_wrapper_call(composition):
    """One wrapper call per trajectory over its baselines, bit for bit."""
    import magnus.oscprob as op
    fig, ax, P = _profile(trajectories=[dict(costhz=-0.9, color='C3')],
                          energy=3.0*gd.UNIT_GEV, **composition)
    expected = op.osc_prob_3nu_earth(3.0*gd.UNIT_GEV, costhz=-0.9, L=PROFILE_L*gd.UNIT_KM,
                                     nu_i=gd.NUMU, nu_f=gd.NUE, **composition)
    assert np.array_equal(P[0], np.asarray(expected, dtype=float))
    assert len(ax) == 2                                 # density panel + one probability
    assert np.array_equal(ax[1].get_lines()[0].get_ydata(), P[0])
    assert ax[0].get_lines()[0].get_color() == 'C3'


def test_profile_computed_over_energy_is_the_earth_wrapper_call_at_the_full_chord():
    import magnus.earth as earth
    import magnus.oscprob as op
    fig, ax, P = mp.plot_probability_with_profile(
        PROFILE_E, trajectories=[dict(loc_ini='fermilab', loc_fin='cern')], x_axis='energy',
        nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3, return_probability=True)
    expected = op.osc_prob_3nu_earth(PROFILE_E*gd.UNIT_GEV, loc_ini='fermilab', loc_fin='cern',
                                     nu_i=gd.NUMU, nu_f=gd.NUE)
    assert np.array_equal(P[0], np.asarray(expected, dtype=float))
    assert len(ax) == 1                                 # no density panel on this axis
    assert 'GeV' in ax[-1].get_xlabel()
    del earth


def test_profile_density_panel_is_the_electron_density_the_wrapper_integrates():
    """rho * Y_e per layer (to the nucleon-mass factor), following a Y_e override."""
    import magnus.earth as earth
    r = earth.earth_radial_distance_from_depth(-0.9, PROFILE_L)
    rho = earth.density_matter_func_prem(r)
    layered = _profile(trajectories=[dict(costhz=-0.9)], energy=gd.UNIT_GEV)[1][0]
    uniform = _profile(trajectories=[dict(costhz=-0.9)], energy=gd.UNIT_GEV,
                       electron_fraction=0.5)[1][0]
    y_layered = layered.get_lines()[0].get_ydata()
    y_uniform = uniform.get_lines()[0].get_ydata()
    assert np.allclose(y_layered, rho*earth.electron_fraction_func_prem(r), rtol=1e-2)
    assert np.allclose(y_uniform, rho*0.5, rtol=1e-2)
    assert not np.allclose(y_layered, y_uniform, rtol=1e-3)


def test_profile_trajectories_carry_their_own_abscissa_and_can_share_a_panel():
    other = np.linspace(100.0, 1000.0, 5)
    fig, ax, P = _profile(trajectories=[dict(costhz=-0.9), dict(costhz=-0.1, x=other)],
                          energy=gd.UNIT_GEV, panel_per_trajectory=False)
    assert len(ax) == 2 and len(ax[1].get_lines()) == 2
    assert [len(p) for p in P] == [len(PROFILE_L), 5]


@pytest.mark.parametrize('kw, message', [
    (dict(trajectories=[dict(costhz=-0.9)]), 'needs energy'),
    (dict(trajectories=[dict(costhz=-0.9)], x_axis='energy', energy=1e9), 'cannot also be given'),
    (dict(trajectories=[dict(costhz=-0.9)], x_axis='energy', show_profile=True),
     'needs the baseline axis'),
    (dict(trajectories=[dict(costhz=-0.9)], energy=1e9, x_axis='zenith'), "'baseline' or 'energy'"),
    (dict(trajectories=[dict(color='C0')], energy=1e9), 'either costhz or both'),
    (dict(trajectories=[dict(costhz=-0.9, loc_ini='cern', loc_fin='fermilab')], energy=1e9),
     'either costhz or both'),
    (dict(trajectories=[dict(costhz=-0.9)], energy=1e9, wrapper_kw=dict(L=1.0)),
     'cannot be given in wrapper_kw'),
    (dict(trajectories=[dict(costhz=-0.9)], energy=1e9, wrapper_kw=dict(detector_depth=1.0)),
     'surface to surface'),
    (dict(trajectories=[], energy=1e9), 'trajectories is empty'),
])
def test_profile_computation_refuses_what_it_cannot_honour(kw, message):
    with pytest.raises(ValueError, match=message):
        _profile(**kw)


def test_profile_refuses_curves_and_trajectories_together(sample):
    L, exact, _ = sample
    with pytest.raises(ValueError, match='not both'):
        mp.plot_probability_with_profile(L, None, [[dict(y=exact)]],
                                         trajectories=[dict(costhz=-0.5)], energy=1e9)
    with pytest.raises(ValueError, match='not both'):
        mp.plot_probability_with_profile(L, None, [[dict(y=exact)]], energy=1e9)
    with pytest.raises(ValueError, match='not both'):
        mp.plot_probability_with_profile(L, None, [[dict(y=exact)]], return_probability=True)
    with pytest.raises(ValueError, match='give panels'):
        mp.plot_probability_with_profile(L)


# ----------------------------------------------------------------------
# plot_biprobability computing through the Earth wrappers
# ----------------------------------------------------------------------

BIPROB_DCP = np.linspace(-np.pi, np.pi, 4)
DUNE = dict(loc_ini='fermilab', loc_fin='homestake')


def _biprob_by_hand(fn, energy, phases, **kw):
    one = lambda d, nubar: np.asarray(fn(energy, dCP=d, nubar=nubar, nu_i=gd.NUMU,   # noqa: E731
                                         nu_f=gd.NUE, **kw), dtype=float).item()
    return (np.array([one(d, False) for d in phases]),
            np.array([one(d, True) for d in phases]))


def test_biprobability_computed_is_the_earth_wrapper_loop():
    """Per curve and phase, the wrapper's neutrino and antineutrino answers, bit for bit;
    a configuration's own parameters override the shared ones for its curve only, and a
    dCP carried by a parameter set gives way to the phases scanned."""
    import magnus.oscprob as op
    no = gd.load_nufit_params('NuFIT 6.1', 'NO')
    io = gd.load_nufit_params('NuFIT 6.1', 'IO')
    fig, ax, P_nu, P_nubar = mp.plot_biprobability(
        configurations=[dict(DUNE), dict(DUNE, **io)], dcp=BIPROB_DCP,
        energy=2.0*gd.UNIT_GEV, num_flavors=3, osc_params=no, return_probability=True)
    for i, params in enumerate((no, io)):
        params = {k: v for k, v in params.items() if k != 'dCP'}
        nu, nubar = _biprob_by_hand(op.osc_prob_3nu_earth, 2.0*gd.UNIT_GEV, BIPROB_DCP,
                                    **DUNE, **params)
        assert np.array_equal(P_nu[i], nu) and np.array_equal(P_nubar[i], nubar)
        assert np.array_equal(ax.get_lines()[i].get_xdata(), nu)


def test_biprobability_costhz_configuration_defaults_to_the_full_chord():
    import magnus.earth as earth
    import magnus.oscprob as op
    *_, P_nu, P_nubar = mp.plot_biprobability(
        configurations=[dict(costhz=-0.3, energy=3.0*gd.UNIT_GEV)], dcp=BIPROB_DCP,
        num_flavors=3, electron_fraction=0.5, return_probability=True)
    L = earth.distance_traveled_inside_earth(-0.3)*gd.UNIT_KM
    nu, nubar = _biprob_by_hand(op.osc_prob_3nu_earth, 3.0*gd.UNIT_GEV, BIPROB_DCP,
                                costhz=-0.3, L=L, electron_fraction=0.5)
    assert np.array_equal(P_nu[0], nu) and np.array_equal(P_nubar[0], nubar)


def test_biprobability_marker_at_a_phase_is_computed():
    import magnus.oscprob as op
    fig, ax = mp.plot_biprobability(configurations=[dict(DUNE)], dcp=BIPROB_DCP,
                                    energy=2.0*gd.UNIT_GEV, num_flavors=3,
                                    markers=[dict(dcp=0.5, marker='*')], legend=False)
    nu, nubar = _biprob_by_hand(op.osc_prob_3nu_earth, 2.0*gd.UNIT_GEV, [0.5], **DUNE)
    point = ax.collections[0].get_offsets()[0]
    assert np.array_equal(np.asarray(point), [nu[0], nubar[0]])


def test_biprobability_channel_sets_the_default_labels(biprob):
    _, p_nu, p_nubar = biprob
    fig, ax = mp.plot_biprobability([p_nu], [p_nubar], nu_i=gd.NUE, nu_f=gd.NUMU)
    assert ax.get_xlabel() == mp.prob_label(gd.NUE, gd.NUMU)
    assert ax.get_ylabel() == mp.prob_label(gd.NUE, gd.NUMU, nubar=True)


@pytest.mark.parametrize('kw, message', [
    (dict(num_flavors=2, osc_params=dict(sth=0.5, Dm2=2e-3)), 'num_flavors = 3, 4 or 5'),
    (dict(num_flavors=3, wrapper_kw=dict(nubar=True)), 'cannot be given in wrapper_kw'),
    (dict(num_flavors=3, wrapper_kw=dict(costhz=-0.5)), 'cannot be given in wrapper_kw'),
    (dict(num_flavors=3, configurations=[dict(DUNE, nubar=True)]), 'cannot be given in a configuration'),
    (dict(num_flavors=3, configurations=[dict(costhz=-0.5)], energy=None), 'gives no energy'),
    (dict(num_flavors=3, configurations=[dict(DUNE, L=1.0)]), 'L goes with costhz'),
    (dict(num_flavors=3, configurations=[dict(loc_ini='cern')]), 'either costhz or both'),
    (dict(num_flavors=3, configurations=[]), 'configurations is empty'),
    (dict(num_flavors=3, nu_i=gd.NUMU), 'both nu_i and nu_f'),
    (dict(num_flavors=3, markers=[dict(dcp=0.0, index=1)]), "one of an 'index'"),
])
def test_biprobability_computation_refuses_what_it_cannot_honour(kw, message):
    kw = dict(dict(configurations=[dict(DUNE)], dcp=BIPROB_DCP, energy=2.0*gd.UNIT_GEV), **kw)
    with pytest.raises(ValueError, match=message):
        mp.plot_biprobability(**kw)


def test_biprobability_refuses_probabilities_and_the_arguments_that_compute_them(biprob):
    _, p_nu, p_nubar = biprob
    with pytest.raises(ValueError, match='not both'):
        mp.plot_biprobability([p_nu], [p_nubar], configurations=[dict(DUNE)])
    with pytest.raises(ValueError, match='not both'):
        mp.plot_biprobability([p_nu], [p_nubar], num_flavors=3)
    with pytest.raises(ValueError, match="'dcp' is a phase to compute"):
        mp.plot_biprobability([p_nu], [p_nubar], markers=[dict(dcp=0.0)])
    with pytest.raises(ValueError, match='give prob_nu and prob_nubar'):
        mp.plot_biprobability([p_nu])


# ----------------------------------------------------------------------
# the optional dependency
# ----------------------------------------------------------------------

def test_missing_matplotlib_names_the_fix(monkeypatch):
    """The error has to say what to do about it, not just that an import failed.

    Matplotlib ships with Magnus, so reaching this at all means it was removed
    from the environment rather than never installed -- which is why the message
    names the package itself rather than an extra.
    """
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith('matplotlib'):
            raise ImportError('No module named matplotlib')
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', fake_import)
    with pytest.raises(mp.MatplotlibNotFoundError, match=r"pip install matplotlib"):
        mp._mpl()


def test_matplotlib_not_found_error_is_an_import_error():
    """So that `except ImportError` around an optional-feature import works."""
    assert issubclass(mp.MatplotlibNotFoundError, ImportError)


def test_importing_magnus_does_not_require_matplotlib():
    """The core package must stay installable and usable without the extra."""
    import magnus
    assert 'plotting' in magnus.submodules
