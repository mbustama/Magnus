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
        r'Two-neutrino probability,~$P_{\nu_e \to \nu_e}$'


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
    assert ax.get_ylabel() == r'Probability,~$P_{\nu_e \to \nu_\mu}$'


def test_vs_energy_builds_its_ylabel_from_the_flavor_pair(sample):
    _, exact, _ = sample
    E = np.logspace(-1.0, 1.0, len(exact))
    fig, ax = mp.plot_probability_vs_energy(E, [exact], nu_i=gd.NUMU,
                                           nu_f=gd.NUTAU, num_flavors=3)
    assert ax.get_ylabel() == \
        r'Three-neutrino probability,~$P_{\nu_\mu \to \nu_\tau}$'


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
# the optional dependency
# ----------------------------------------------------------------------

def test_missing_matplotlib_names_the_extra_to_install(monkeypatch):
    """The error has to say what to do about it, not just that an import failed."""
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith('matplotlib'):
            raise ImportError('No module named matplotlib')
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', fake_import)
    with pytest.raises(mp.MatplotlibNotFoundError, match=r"magnuspy\[plot\]"):
        mp._mpl()


def test_matplotlib_not_found_error_is_an_import_error():
    """So that `except ImportError` around an optional-feature import works."""
    assert issubclass(mp.MatplotlibNotFoundError, ImportError)


def test_importing_magnus_does_not_require_matplotlib():
    """The core package must stay installable and usable without the extra."""
    import magnus
    assert 'plotting' in magnus.submodules
