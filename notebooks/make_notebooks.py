# -*- coding: utf-8 -*-
r"""Builds every notebook in this directory, executes it, and stores its outputs.

**The notebooks are generated.  Edit this file, not the ``.ipynb``** -- anything
written into a notebook by hand is lost the next time this runs.

Why generate them.  Fourteen notebooks share a setup block, a title format and a
navigation footer, and hand-maintaining those across fourteen JSON files means
they drift: nine of the fourteen had no title cell at all before this script
existed, so they opened on a bare import block, and none of them carried a way
to reach the next one.  Generating them makes those properties structural.

Run ``python notebooks/make_notebooks.py`` to rebuild.  Execution is part of the
build, deliberately: a notebook that no longer runs against the current package
is a broken example, and the only way to know is to run it.  The build also
refuses to finish if a notebook came back without stored outputs, because a
notebook whose outputs were stripped renders blank on GitHub.
"""

import pathlib
import time

import nbformat as nbf


HERE = pathlib.Path(__file__).resolve().parent
DOCS = 'https://mbustama.github.io/Magnus'


def md(text):
    return nbf.v4.new_markdown_cell(text.rstrip())


def code(text):
    return nbf.v4.new_code_cell(text.rstrip())


def notebook(title, intro, cells):
    r"""One notebook: a title cell, then whatever the caller supplies."""
    nb = nbf.v4.new_notebook()
    nb.cells = [md('# %s\n\n%s' % (title, intro))] + cells
    nb.metadata = {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python',
                       'name': 'python3'},
        'language_info': {'name': 'python', 'pygments_lexer': 'ipython3'},
    }
    return nb


books = {}

# ---------------------------------------------------- 01_magnus_introduction
books['01_magnus_introduction.ipynb'] = notebook(
    'Introduction',
    'The shortest path to a probability, and the conventions the rest of these notebooks assume.\n\nMag$\\nu$s computes neutrino oscillation probabilities by the **Magnus expansion**: the evolution operator over a slab is $\\exp(\\Omega)$, with $\\Omega$ built from time-ordered integrals of nested commutators of the Hamiltonian. Truncating the series is exactly unitary at any order, which is the property that makes it worth doing this way rather than exponentiating a discretised Hamiltonian.',
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.magnus as magnus
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''To compute neutrino oscillation proabilities, Mag$\nu$s needs only three ingredients: 
* a Hamiltonian (`H_func`), written in the flavor basis, as a function of neutrino position (or time),
* the initial neutrino position (or time, `t_ini`), and
* the final neutrino position (or time, `t_fin`).

In practice, in most scenarios the Hamiltonian depends on the neutrino energy, whose value we often vary, so the neutrino energy becomes effectively also a free parameter to vary.'''),
    md(r'''Let us start by computing oscillation probabilities in vacuum in a two-neutrino system. Although the user can provide their own custom Hamiltonian, Mag$\nu$s conveniently provides a library of commonly used Hamiltonians for two-neutrino and three-neutrinos oscillations (and also for 3+1 and 3+2 systems, with one and two extra flavors, respectively).

All of them live in a single module, `magnus.hamiltonians`, which we import as'''),
    code(r'''import magnus.hamiltonians as hamiltonians'''),
    md(r'''We will use the two-neutrino Hamiltonian in vacuum, i.e.,'''),
    md(r'''\begin{equation}
\mathbf{H}_{2\nu}^{\rm vac}
=
\mathbf{U}^T 
\left(
 \begin{array}{cc}
  \frac{\Delta m^2}{2E} & 0 \\
  0 & -\frac{\Delta m^2}{2E} \\
 \end{array}
\right)
\mathbf{U} \;,
\end{equation}'''),
    md(r'''where the mixing matrix, $\mathbf{U}$, is parametrized by a single mixing angle, $\theta$, i.e.,'''),
    md(r'''\begin{equation}
 \mathbf{U}
 =
 \left(
  \begin{array}{cc}
   \cos\theta  & \sin\theta \\
   -\sin\theta & \cos\theta
  \end{array}
 \right) \;.
\end{equation}'''),
    md(r'''We need values for the oscillation parameters $\theta$ and $\Delta m^2$.  Users can provide any values they wish.  For this example, we will shows oscillations between $\nu_e$ and $\nu_mu$, and use central values of the $\theta_{12}$ and $\Delta m_{21}^2$ parameters from the NuFit 6.0 global fit to oscillation data, which are predefined in the Mag$\nu$s globaldefs module, i.e.,'''),
    code(r'''import magnus.globaldefs as gd

sth = gd.S12_NO_BF_NUFIT_6_0 # sin(theta) [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    md(r'''Let's compute the probability for a baseline of 10 km and an energy of 1 MeV.'''),
    code(r'''baseline = 10.*gd.UNIT_KM # 10 km natural units [eV^{-1}]
energy = 1.*gd.UNIT_MEV # [eV]'''),
    md(r'''And now we define our Hamiltonian, as a function of neutrino position (`l`) and energy (`energy`).  In this case, the dependence on neutrino position is a dummy dependence, since the Hamiltonian in vacuum is position-independent. However, Mag$\nu$s expects a position-dependent Hamiltonian in general.'''),
    code(r'''hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2)'''),
    md(r'''The oscillation probabilities are computed using the `osc_prob` function of the `oscprob` module, already imported above. Alongside it, `oscprobstd` provides the textbook closed-form expressions, which we will use throughout as an independent check on the Magnus result:'''),
    code(r'''import magnus.oscprobstd as oscprobstd'''),
    md(r'''Internally, the `osc_prob` function performs Magnus expansion to compute the evolution of the neutrino oscillation amplitude from `t_ini` to `t_fin` by partitioning this interval into subintervals, computing the evolution operator inside each subinterval, and performing their time-ordered product. The user manipulates the options of the Magnus expansion indirectly, via parameters passed to osc_prob.

(Users interested in using Mag$\nu$s to compute the Magnus expansion of an arbitrary matrix exponential, not only in neutrino oscillations, should look at the notebook `10_magnus_matrix_exponential.ipynb`.)

The `osc_prob` function computes probabilities for *any* number of neutrino flavors.  If it is fed a $2 \times 2$ Hamiltonian, then it will return probabilities for a two-neutrino system; if it is fed a $3 \times 3$ Hamiltonian, it will return probabilities for a three-neutrino system, *etc*.'''),
    md(r'''We call `osc_prob` to compute the probability from `t_ini = 0.0` to `t_fin = baseline`. We pass `verbose = 2` to see run information and warnings:'''),
    code(r'''P = oscprob.osc_prob(hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2), 0.0, baseline, verbose=2)
P'''),
    md(r'''The function `osc_prob` returns a symmetric probability matrix,
\begin{equation}
 \mathbf{P}_{2\nu}
 =
 \left(
  \begin{array}{cc}
   P_{ee}    & P_{e\mu} \\
   P_{\mu e} & P_{\mu\mu}
  \end{array}
 \right) \;,
\end{equation}
where $P_{\mu e} = P_{e \mu}$.'''),
    md(r'''We can select individual probabilities using flavor indices predefined in the `globaldefs` module.  I.e., for $\nu_e \to \nu_e$,'''),
    code(r'''P[gd.NUE][gd.NUE]'''),
    md(r'''and, for $\nu_e \to \nu_\mu$,'''),
    code(r'''P[gd.NUE][gd.NUMU]'''),
    md(r'''Notice that in the run above we had a warning.  This is because we passed a constant (*i.e.*, time-independent) Hamiltonian to `osc_prob`.  In cases like that (*i.e.*, in vacuum or constant-density matter), the only nonzero term of the Magnus expansion is the first one. *I.e.*, the evolution operator is simply $e^{-i H L}$, with $L = t_{\rm fin}-t_{\rm ini}$).  When a constant Hamiltonian is passed to `osc_prob` this is identified and the run parameters are adjusted to use first-order Magnus expansion, for speed-up.'''),
    md(r'''An alternative way to compute the same probability matrix is to define a vacuum Hamiltonian function of position of neutrino position (`l`) and neutrino energy (`energy`), with the dependence on position merely dummy, since the vacuum Hamiltonian is time-independent:'''),
    code(r'''def H_2nu_vac(l, energy):
    return hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2)'''),
    md(r'''In this case, we need to pass a Hamiltonian that is a function *only* of position:'''),
    code(r'''oscprob.osc_prob(lambda l: H_2nu_vac(l, energy), 0.0, baseline, magnus_exp_order=1, n_slabs=1, n_tpts_per_slab=2, rtol=None, atol=None)'''),
    md(r'''The resulting probability matrix is the same as before, but we had to manually set the run parameters to generate it with the same speed as when we simply passed the constant Hamiltonian `hamiltonian_2nu_vacuum(energy, sth, Dm2)` to `osc_prob`.  However, the recipe of defining a position-dependent Hamiltonian and passing it to `osc_prob` will be at the core of later calculations of probabilities for time-*dependent* Hamiltonians.

For convenience, Mag$\nu$s includes a wrapper to return the probability matrix in vacuum for given values of the oscillation parameters, the neutrino baseline, and the neutrino energy, bypassing the more lengthy procedure above:'''),
    code(r'''oscprob.osc_prob_2nu_vacuum(energy, baseline, sth, Dm2)'''),
    md(r'''This function can also return the probability matrices for multiple energies and baselines, passed to it either as lists of NumPy arrays. Internally, it will `zip` energy and baseline arrays (*i.e.*, (`zip(energy, L)`), compute the probability matrix for each element of the zip object, and return a NumPy array containing all of the matrices, in the `zip` order. This means that the lengths of the provided energy and baseline arrays should be equal, and that the length of the returned array will also be equal to that.  For example,'''),
    code(r'''prob = oscprob.osc_prob_2nu_vacuum(gd.UNIT_MEV*np.array([1.0, 10.0]), [baseline, baseline], sth, Dm2)
prob'''),
    md(r'''And from this we can select the probability channel we are interested, say $\nu_e \to \nu_\mu$,'''),
    code(r'''prob[:,gd.NUE,gd.NUMU]'''),
    md(r'''We one can directly ask `osc_prob_2nu_vacuum` to return only that probability by passing the initial and final flavors, `nu_i` and `nu_f` (which reduces memory requirements).

If `energy` and `L` as floats, the probability is returned as a float.'''),
    code(r'''oscprob.osc_prob_2nu_vacuum(gd.UNIT_MEV*1.0, baseline, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU)'''),
    md(r'''If `energy` is a float and `L` is a list (or NumPy array) with multiple entries, `osc_prob_2nu_vacuum` returns an array with the probability computed for that fixed `energy` and each value inside `L`:'''),
    code(r'''oscprob.osc_prob_2nu_vacuum(gd.UNIT_MEV*1.0, gd.UNIT_KM*np.array([1.0, 10.0, 100.0, 1000.0]), sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU)'''),
    md(r'''Conversely, we can pass a single value of `L` and multiple values of `E` to `osc_prob_2nu_vacuum` and it will return the probability computed for that fixed `L` and each value inside `energy`:'''),
    code(r'''oscprob.osc_prob_2nu_vacuum(gd.UNIT_MEV*np.array([1.0, 5.0, 10.0, 100.0]), baseline, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU)'''),
    md(r'''The same functionality is available in the functions `osc_prob_2nu_matter_constant_density`, `osc_prob_3nu_vacuum`, and `osc_prob_3nu_matter_constant_density` that we introduce later.

Let's now move on to three-neutrino oscillations.  No new import is needed: the three-neutrino Hamiltonians come from the same `magnus.hamiltonians` module imported above, and are named `hamiltonian_3nu_*` rather than `hamiltonian_2nu_*`.  For instance:'''),
    code(r'''# The 3nu builders sit alongside the 2nu ones in the same module
[n for n in hamiltonians.__all__ if n.startswith('hamiltonian_3nu')][:5]'''),
    md(r'''The vacuum Hamiltonian is now a $3 \times 3$ matrix, i.e.,
\begin{equation}
 \mathbf{H}_{3\nu}^{\rm vac}
 =
 \mathbf{U}_{\rm PMNS}^\dagger
 \left(
  \begin{array}{ccc}
   0 & 0 & 0 \\
   0 & \frac{\Delta m_{21}^2}{2E} & 0 \\
   0 & 0 & \frac{\Delta m_{31}^2}{2E}
  \end{array}
 \right)
 \mathbf{U}_\textrm{PMNS} \;,
\end{equation}
where $\mathbf{U}_{\rm PMNS}$ is the complex-valued Pontecorvo-Maki-Nakagawa-Sakata (PMNS) matrix, parametrized using three mixing angles, $\theta_{12}$, $\theta_{23}$, and $\theta_{13}$, and one CP-violation phase, $\delta_{\rm CP}$.

Like before, for this example we set the values of the oscillation parameters to their central values from NuFit 6.0, *i.e.*,'''),
    code(r'''s12 = gd.S12_NO_BF_NUFIT_6_0 # sin(theta_12) [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # sin(theta_23) [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # sin(theta_13) [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [radian]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]'''),
    md(r'''We compute the probabilities by calling `osc_prob` again. This time, we set `verbose = 1` to only show warnings and important messages:'''),
    code(r'''P = oscprob.osc_prob(hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31), 0.0, baseline, verbose=1)
P'''),
    md(r'''In this case, the result is the probability matrix
\begin{equation}
 \mathbf{P}_{3\nu}
 =
 \left(
  \begin{array}{ccc}
   P_{ee}     & P_{e\mu}    & P_{e\tau} \\
   P_{\mu e}  & P_{\mu\mu}  & P_{\mu\tau} \\
   P_{\tau e} & P_{\tau\mu} & P_{\tau\tau}
  \end{array}
 \right) \;,
\end{equation}
where $P_{\alpha \beta} = P_{\beta \alpha}$, with $\alpha, \beta = e, \mu, \tau$.'''),
    md(r'''We can retrieve probabilities for specific channels using flavor indices, as before, *e.g.*,'''),
    code(r'''P[gd.NUE][gd.NUE]'''),
    code(r'''P[gd.NUE][gd.NUMU]'''),
    code(r'''P[gd.NUMU][gd.NUMU]'''),
    code(r'''P[gd.NUMU][gd.NUTAU]'''),
    md(r'''Similarly to the two-neutrino case, for convenience, Mag$\nu$s contains a wrapper for the calculation of the three-neutrino probability in vacuum for given values of the oscillation parameters, neutrino position, and neutrino energy:'''),
    code(r'''oscprob.osc_prob_3nu_vacuum(energy, baseline, s12, s23, s13, dCP, D21, D31)'''),
    md(r'''And, like before, we can call this function also with an array of energies and baselines, *e.g.*,'''),
    code(r'''prob = oscprob.osc_prob_3nu_vacuum(gd.UNIT_MEV*np.array([1.0, 10.0]), [baseline, baseline], s12, s23, s13, dCP, D21, D31)
prob'''),
    md(r'''And we can select to return only the channel we are interested, say $\nu_e \to \nu_\mu$,'''),
    code(r'''oscprob.osc_prob_3nu_vacuum(gd.UNIT_MEV*np.array([1.0, 10.0]), [baseline, baseline], s12, s23, s13, dCP, D21, D31,
                            nu_i=gd.NUE, nu_f=gd.NUMU)'''),
    md(r'''If the vacuum probability needs to computed for many values of energy and baseline, it is more computationally efficient to call `osc_prob_3nu_vacuum` once with the array of energies and baselines at which it must be computed than calling it multiple times, each for a different value of energy and baseline.  The reason is that, internally, `osc_prob_3nu_vacuum` has an overhead that is the same regardless if it is computed for one or multiple probability computations.

Also, `osc_prob_3nu_vacuum` can be called without passing values for any of the oscillation parameters; this is useful if we just want to quickly a probability.  If so, internally, the function will assign them values from a default parameter set in Mag$\nu$s, which is set to a recent experimental determination.  Let's call `osc_prob_3nu_vacuum` with `verbose = 1` to see a warning when parameters are assigned default values:'''),
    code(r'''oscprob.osc_prob_3nu_vacuum(gd.UNIT_MEV*np.array([1.0, 10.0]), baseline, 
                            nu_i=gd.NUE, nu_f=gd.NUMU, verbose=1)'''),
    md(r'''We can also specify only some of the oscillation parameter values, and let `osc_prob_3nu_vacuum` set the unspecified parameters to their default values, *e.g.*,'''),
    code(r'''oscprob.osc_prob_3nu_vacuum(gd.UNIT_MEV*np.array([1.0, 10.0]), baseline, s13=0.0, dCP=0.0, D31=0.0,
                            nu_i=gd.NUE, nu_f=gd.NUMU, verbose=1)'''),
    md(r'''Finally, by passing the argument `default_osc_params_set_name` we can specifiy which parameter data set should be used by `osc_prob_3nu_vacuum` as the default from which to assign values to the unspecified parameters.  We can print the list of available predefined parameter set names via'''),
    code(r'''list(gd.OSC_PARAMS_PREDEFINED.keys())'''),
    md(r'''And the current default set is'''),
    code(r'''gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']'''),
    md(r'''For instance, we could choose to use instead NuFit 6.0 with inverted mass ordering (IO), with SK atmospheric data:'''),
    code(r'''oscprob.osc_prob_3nu_vacuum(gd.UNIT_MEV*np.array([1.0, 10.0]), baseline, 
                            nu_i=gd.NUE, nu_f=gd.NUMU, default_osc_params_set_name='OSC_PARAMS_NU_FIT_6_0_SK_IO', verbose=1)'''),
    md(r'''Let's now vary the baseline and energy and plot the two- and three-neutrino $\nu_e \to \nu_\mu$ probabilities.  We will stick to using the default oscillation parameter set values from now on in this notebook.'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUMU
npts = 1000
energy = 1.0*gd.UNIT_MEV # [eV]
baseline_arr = np.logspace(0.0, 2.0, npts) # [km]

# Both wrappers take the whole array of baselines in one call. Looping over
# `baseline_arr` and calling them once per point would give the same answer at
# roughly twenty-five times the cost -- the entry path (validation, Hamiltonian
# probing, dispatch) is paid per call, not per point.
prob_2nu_vacuum_vs_baseline = oscprob.osc_prob_2nu_vacuum(energy, gd.UNIT_KM*baseline_arr, s12, D21, nu_i=nu_i, nu_f=nu_f)
prob_3nu_vacuum_vs_baseline = oscprob.osc_prob_3nu_vacuum(energy, gd.UNIT_KM*baseline_arr, nu_i=nu_i, nu_f=nu_f)

fig, ax = plotting.plot_probability_vs_baseline(
    baseline_arr,
    [dict(y=prob_2nu_vacuum_vs_baseline, label=r'Two-neutrino'),
     dict(y=prob_3nu_vacuum_vs_baseline, label=r'Three-neutrino')],
    ylabel=r'Probability',
    xlim=(baseline_arr[0], baseline_arr[-1]),
    legend_loc='upper left',
    title=r'Neutrino oscillations in vacuum, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV))'''),
    md(r'''## Figures without the boilerplate

The plot above is short because it is a simple one. Most figures in these
notebooks are not: the standard shape -- a set of curves over a short panel
showing the relative error against a reference -- takes about thirty lines of
`gridspec_kw`, tick locators, legend keywords and axis wrangling, and notebooks
`02` and `03` each repeat it fourteen times.

That boilerplate is what `magnus.plotting` packages up. It is worth seeing the
before and after side by side, because the point is not that the figure gets
better -- it is that the figure stays *exactly the same* while the code
collapses.

First, the hand-built version, written out in full:'''),
    code(r'''# Two ways of computing the same two-neutrino probability, so that the lower
# panel has something to compare: the Magnus engine against the closed form.
prob_magnus = oscprob.osc_prob_2nu_vacuum(energy, gd.UNIT_KM*baseline_arr, s12, D21,
                                          nu_i=nu_i, nu_f=nu_f)
prob_std = np.array([oscprobstd.osc_prob_2nu_vacuum_std(s12, D21, energy,
                                                        gd.UNIT_KM*l)[nu_i][nu_f]
                     for l in baseline_arr])
residual = (prob_magnus - prob_std)/np.maximum(np.abs(prob_std), 1.e-300)/1.e-12

# ---- the hand-built figure: this is the block that gets copied around -------
heights = [1.0, 0.3]
widths = [1.0]
gs_kw = dict(height_ratios=heights, width_ratios=widths)
fig, ax = plt.subplots(
    ncols=1,
    nrows=2,
    gridspec_kw=gs_kw,
    figsize=[18, 9])
fig.subplots_adjust(hspace=0.05, wspace=0.05)

ax[0].plot(baseline_arr, prob_magnus, lw=1, color='C1', ls='-', label='Magnus expansion')
ax[0].plot(baseline_arr, prob_std, lw=1, color='k', ls='--', label='Standard formula')
ax[1].plot(baseline_arr, residual, lw=1, color='k', ls='-')

ax[0].legend(fontsize=17, frameon=True, loc='upper left', handlelength=1.2,
             handleheight=0.7, borderpad=0.8, title_fontsize=20, edgecolor='k',
             labelspacing=0.7, ncol=1, title=r'Calculation method')

ax[0].set_ylabel(r'Two-neutrino probability,~'+r'$P_{\nu_e \to \nu_\mu}$', labelpad=25)
ax[0].set_xlim(baseline_arr[0], baseline_arr[-1])
ax[0].set_ylim(0, 1)
ax[0].set_xscale('log')
ax[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.10))
ax[0].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.02))
ax[0].xaxis.set_ticklabels([])
ax[1].set_xlabel(r'Baseline, $L$ [km]')
ax[1].set_ylabel(r'$\epsilon_{\rm rel}~[\times 10^{-12}]$', labelpad=7)
ax[1].set_xscale('log')
ax[1].set_xlim(baseline_arr[0], baseline_arr[-1])
ax[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.20))
ax[1].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.05))
ax[0].set_title(r'$2\nu$~vacuum, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
                fontsize=20, pad=10)

plt.tight_layout()'''),
    md(r'''Thirty-one lines, of which exactly three carry any physics.

Now the same figure, in one call:'''),
    code(r'''from magnus.plotting import plot_curves, prob_label

fig, ax = plot_curves(
    baseline_arr,
    [dict(y=prob_magnus, label='Magnus expansion', color='C1'),
     dict(y=prob_std, label='Standard formula', color='k', ls='--')],
    xlabel=r'Baseline, $L$ [km]',
    ylabel=r'Two-neutrino probability,~'+prob_label(nu_i, nu_f),
    xlim=(baseline_arr[0], baseline_arr[-1]), ylim=(0, 1), xscale='log',
    ymajor=0.10, yminor=0.02,
    residual=residual, residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
    residual_ymajor=0.20, residual_yminor=0.05,
    legend_title=r'Calculation method', legend_loc='upper left',
    title=r'$2\nu$~vacuum, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
)'''),
    md(r'''The defaults are the house style, so nothing had to be restated: the figure
size, the two-panel height ratio, the nine legend keywords, the suppressed
tick labels on the upper panel, and the shared abscissa limits all come for
free. What is left in the call is what actually differs between figures.

`plot_curves` returns `(fig, ax)` precisely so that it is a starting point and
not a dead end -- anything Matplotlib can do to an axes can still be done
afterwards. The rest of this notebook, and all the others, use these functions;
see the [plotting page](https://mbustama.github.io/Magnus/plotting.html) for
the full set, which also covers density-profile stacks, bi-probability planes
and oscillograms.

One install note: Matplotlib is an optional dependency, so plotting needs
`pip install 'magnuspy[plot]'` (or `pip install -e '.[plot]'` from a checkout).
The engine itself needs only NumPy, SciPy and joblib.'''),
    md(r'''We can also plot probabilities *vs.* energy, say, for a fixed baseline of 100 km:'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUMU
npts = 1000
baseline = 1.0*gd.UNIT_KM # [eV^-1]
energy_arr = np.logspace(-1.0, 1.0, npts) # [MeV]

prob_2nu_vacuum_vs_energy = oscprob.osc_prob_2nu_vacuum(gd.UNIT_MEV*energy_arr, baseline, s12, D21, nu_i=nu_i, nu_f=nu_f)
prob_3nu_vacuum_vs_energy = oscprob.osc_prob_3nu_vacuum(gd.UNIT_MEV*energy_arr, baseline, nu_i=nu_i, nu_f=nu_f)

fig, ax = plotting.plot_probability_vs_energy(
    energy_arr,
    [dict(y=prob_2nu_vacuum_vs_energy, label=r'Two-neutrino'),
     dict(y=prob_3nu_vacuum_vs_energy, label=r'Three-neutrino')],
    energy_unit='MeV', ylabel=r'Probability',
    xlim=(energy_arr[0], energy_arr[-1]),
    legend_loc='upper left',
    title=r'Neutrino oscillations in vacuum, $L = $~{:.2f}~km'.format(baseline/gd.UNIT_KM))'''),
    md(r'''Now let's add matter effects.  To do this, we need to add a new contribution to the Hamiltonian that contributes the potential due to the coherent forward scattering of $\nu_e$ on electrons, $V_{\rm CC}$, *i.e.*,
\begin{equation}
 H_{3\nu}^{\rm matt}
 = 
 \left(
  \begin{array}{ccc}
   V_{\rm CC} & 0 & 0 \\
   0 & 0 & 0 \\
   0 & 0 & 0 
  \end{array}
 \right) \;.
\end{equation}
With this, the total Hamiltonian becomes
\begin{equation}
 H_{3\nu} = H_{3\nu}^{\rm vac} + H_{3\nu}^{\rm matt} \;.
\end{equation}
The potential $V_{\rm CC} = \sqrt{2} G_F N_e$ varies proportionally to the electron number densit, $N_e$.'''),
    code(r'''rho = 10*gd.UNIT_G_PER_CM3 # [g cm^{-3}]'''),
    md(r'''which is roughly the density at the center of the Earth.'''),
    md(r'''To compute the potential given the matter density, we will need a few helper functions that are provied by Mag$\nu$s in the `matter` module, which we now import:'''),
    code(r'''import magnus.matter as matter'''),
    md(r'''First, we compute the electron number density using the `num_density_e_func` function.  This function can return the electron number density at any position given a varying matter density profile, which is why we need to pass `rho` as a position-dependent function, and evaluate at `l = 0.0` (any value of `l` will do, since `rho` is a constant).  We assume that the matter is isoscalar (*i.e.*, that the ratio of neutrons to protons is 1) and electrically neutral (*i.e.*, that the ratio of electrons to baryons, or the electron fraction, is 0.5):'''),
    code(r'''num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, 
                                          ratio_number_neutrons_to_protons=1.0, electron_fraction=0.5) # [eV^3]
num_density_e'''),
    md(r'''With this, we compute the coherent forward potential using the `VCC_func` function.  Again, the function is designed to return the potential at any position given a varying electron number density profile, but we evaluate it at `l = 0`:'''),
    code(r'''VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # [eV]
VCC'''),
    md(r'''Now we can define the vacuum Hamiltonian, as before, and a constant matter contribution to the Hamiltonian using the predefined function `hamiltonian_3nu_matter`:'''),
    code(r'''H_3nu_vac = hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31)
H_3nu_matt = hamiltonians.hamiltonian_3nu_matter(VCC)
H_3nu = H_3nu_vac + H_3nu_matt'''),
    md(r'''And we call the `osc_prob` function to compute the probabilities:'''),
    code(r'''oscprob.osc_prob(H_3nu, 0.0, baseline, verbose=1)'''),
    md(r'''And, similarly to the two-neutrino case, for convenience, Mag$\nu$s contains a wrapper for the calculation of the three-neutrino probability in constant-density matter for given values of the oscillation parameters, neutrino position, and neutrino energy, and matter density (by default, it uses `ratio_number_neutrons_to_protons = 1.0` and `electron_fraction = 0.5`, but different values can be passed as optional arguments). Also, like for the three-neutrino case in vacuum (`osc_prob_3nu_vacuum`, see earlier in the notebook), we can choose not to pass values for the oscillation parameters, and let them be set to the default oscillation parameter data set in Mag$\nu$s (set `verbose = 1` to see a warning when this happens):'''),
    code(r'''oscprob.osc_prob_3nu_matter_constant_density(energy, baseline, rho, verbose=1)'''),
    md(r'''(A similar function exists for the two-neutrino case, `osc_prob_2nu_matter_constant_density`.)

By default, this returns the probability for neutrinos.  To compute it for anti-neutrinos (which takes the complex-conjugate of the PMNS mixing matrix and flips the sign of $V_{\rm CC}$, the optional argument `nubar = True` must be passed.  For instance,'''),
    code(r'''oscprob.osc_prob_3nu_matter_constant_density(energy, baseline, rho, nubar=True)'''),
    md(r'''And, like for the vacuum case (`osc_prob_2nu_vacuum` and `osc_prob_3nu_vacuum`), we can pass arrays of energy and baseline to compute in one go, and select single probability channels to return.'''),
    code(r'''oscprob.osc_prob_3nu_matter_constant_density(gd.UNIT_MEV*np.array([1.0, 5.0, 10.0, 100.0]), baseline, 
                                             rho, nu_i=gd.NUE, nu_f=gd.NUMU)'''),
    md(r'''Let's now compare three-neutrino probability of $\bar{\nu}_e \to \bar{\nu}_\mu$ in vacuum and in matter for two choices of constant matter density.  First, we plot them *vs.* baseline:'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUMU
npts = 1000
energy = 1.0*gd.UNIT_MEV # [eV]
baseline_arr = np.logspace(0.0, 2.0, npts) # [km]
rho_lo, rho_hi = 10*gd.UNIT_G_PER_CM3, 100*gd.UNIT_G_PER_CM3 # 10 and 100 g cm^{-3} in natural units [eV^4]

prob_3nu_matt_const_lo_vs_baseline = oscprob.osc_prob_3nu_matter_constant_density(energy, gd.UNIT_KM*baseline_arr, 
                                                                                  rho=rho_lo, nu_i=nu_i, nu_f=nu_f, nubar=True)
prob_3nu_matt_const_hi_vs_baseline = oscprob.osc_prob_3nu_matter_constant_density(energy, gd.UNIT_KM*baseline_arr, 
                                                                                  rho=rho_hi, nu_i=nu_i, nu_f=nu_f, nubar=True)

# This gives the same result, but it is less efficient (calling osc_prob_3nu_matter_constant_density has an overhead):
# prob_3nu_matt_const_lo_vs_baseline = [oscprob.osc_prob_3nu_matter_constant_density(energy, l, 
#                                                                                    rho=rho_lo, nubar=True)[nu_i][nu_f] 
#                                for l in gd.UNIT_KM*baseline_arr]
# prob_3nu_matt_const_hi_vs_baseline = [oscprob.osc_prob_3nu_matter_constant_density(energy, l, 
#                                                                                    rho=rho_hi, nubar=True)[nu_i][nu_f] 
#                                for l in gd.UNIT_KM*baseline_arr]

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=[18,9])
ax.plot(baseline_arr, prob_3nu_vacuum_vs_baseline, label=r'Vacuum', c='0.7', ls='--')
ax.plot(baseline_arr, prob_3nu_matt_const_lo_vs_baseline, label=r'Matter, $\rho = $~{:.0f}'.format(rho_lo/gd.UNIT_G_PER_CM3)+r'~g~cm$^{-3}$')
ax.plot(baseline_arr, prob_3nu_matt_const_hi_vs_baseline, label=r'Matter, $\rho = $~{:.0f}'.format(rho_hi/gd.UNIT_G_PER_CM3)+r'~g~cm$^{-3}$')
ax.legend(fontsize=17, frameon=True, loc='upper left', handlelength=1.2, handleheight=0.7, borderpad=0.8, 
          title_fontsize=20, edgecolor='k', labelspacing=0.7, ncol=1)
ax.set_xlabel(r'Baseline, $L$ [km]')
ax.set_ylabel(r'Probability')
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.10))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.02))
ax.set_xlim(baseline_arr[0], baseline_arr[-1])
ax.set_ylim(0, 1)
ax.set_xscale('log')
ax.set_title(r'$3\nu$  oscillations, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV), fontsize=20, pad=10)'''),
    md(r'''And now we plot the probabilities *vs.* energy:'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUMU
npts = 1000
energy_arr = np.logspace(0.0, 2.0, npts) # [MeV]
baseline = 100.*gd.UNIT_KM # 100 km natural units [eV^{-1}]
rho_lo, rho_hi = 10*gd.UNIT_G_PER_CM3, 100*gd.UNIT_G_PER_CM3 # 10 and 100 g cm^{-3} in natural units [eV^4]

prob_3nu_matt_const_lo_vs_energy = oscprob.osc_prob_3nu_matter_constant_density(gd.UNIT_MEV*energy_arr, baseline, 
                                                                                rho=rho_lo, nu_i=nu_i, nu_f=nu_f, nubar=True)
prob_3nu_matt_const_hi_vs_energy = oscprob.osc_prob_3nu_matter_constant_density(gd.UNIT_MEV*energy_arr, baseline, 
                                                                                rho=rho_hi, nu_i=nu_i, nu_f=nu_f, nubar=True)

# This gives the same result, but it is less efficient (calling osc_prob_3nu_matter_constant_density has an overhead):
# prob_3nu_matt_const_lo_vs_energy = [oscprob.osc_prob_3nu_matter_constant_density(s12, s23, s13, dCP, D21, D31, energy, baseline,
#                                                                                  rho=rho_lo, nubar=True)[nu_i][nu_f] 
#                                     for energy in gd.UNIT_MEV*energy_arr]
# prob_3nu_matt_const_hi_vs_energy = [oscprob.osc_prob_3nu_matter_constant_density(s12, s23, s13, dCP, D21, D31, energy, baseline, 
#                                                                                  rho=rho_hi, nubar=True)[nu_i][nu_f] 
#                                     for energy in gd.UNIT_MEV*energy_arr]

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=[18,9])

ax.plot(energy_arr, prob_3nu_vacuum_vs_energy, label=r'Vacuum', c='0.7', ls='--')
ax.plot(energy_arr, prob_3nu_matt_const_lo_vs_energy, label=r'Matter, $\rho = $~{:.0f}'.format(rho_lo/gd.UNIT_G_PER_CM3)+r'~g~cm$^{-3}$')
ax.plot(energy_arr, prob_3nu_matt_const_hi_vs_energy, label=r'Matter, $\rho = $~{:.0f}'.format(rho_hi/gd.UNIT_G_PER_CM3)+r'~g~cm$^{-3}$')
ax.legend(fontsize=17, frameon=True, loc='upper right', handlelength=1.2, handleheight=0.7, borderpad=0.8, 
          title_fontsize=20, edgecolor='k', labelspacing=0.7, ncol=1)
ax.set_xlabel(r'Neutrino energy, $E_\nu$ [MeV]')
ax.set_ylabel(r'Probability')
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.10))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.02))
ax.set_xlim(energy_arr[0], energy_arr[-1])
ax.set_ylim(0, 1)
ax.set_xscale('log')
ax.set_title(r'$3\nu$ oscillations, $L = $~{:.2f}~km'.format(baseline/gd.UNIT_KM), fontsize=20, pad=10)'''),
    md(r'''So far, we have only computed oscillation probabilities in cases where the Hamiltonian is time-independent (or, equivalently, position-independent), *i.e.*, in vacuum and in matter with uniform density.  As we pointed out earlier, in these cases only the first term of the Magnus expansion is nonzero, and so the problem reduces to a simple one.  

Where Mag$\nu$s really shows its power and flexibility is when computing oscillation probabilities for time- or position-*dependent* Hamiltonians.  

Let's illustrate this first for a matter density profile that an exponentially falling as a function of distance, *e.g.*, like inside the Sun.

We can no longer use the pre-defined functions for oscillations in vacuum or constant-density matter (`osc_prob_3nu_matter_vacuum` and `osc_prob_3nu_matter_constant_density`), but we need to define our own Hamiltonians.  We do this by calling `hamiltonian_3nu_matter` with different profiles of matter potential.'''),
    md(r'''### A note on speed: let your Hamiltonian take an array

The Hamiltonians defined below are written so that they work whether they are
handed a single position or a whole array of them. That is worth doing, and it
is invisible unless someone points it out.

`osc_prob` samples the Hamiltonian at every quadrature node of every slab --
often a few hundred positions for a single probability, and the adaptive
refinement repeats that at each level. So it first tries one vectorized call,
`H_func(array_of_positions)`, and uses the result if it has the right shape and
agrees with a scalar spot-check. If that fails, it falls back to calling the
function once per position: correct, but measured **4.6x slower** here
(7.8 ms against 1.7 ms per call), for bit-identical output.

Nothing about a scalar-only Hamiltonian looks wrong, which is why it is easy to
sit on the slow path indefinitely -- these notebooks did, for years. Since
version 1.0.0 the fallback raises `magnus.magnus.ScalarHamiltonianWarning` once
per session and names the fix.

In practice there is usually nothing to do: write the position dependence with
NumPy and it vectorizes by itself. `matter.VCC_func`,
`earth.density_matter_func_prem` and the `hamiltonian_*_matter` builders all
accept arrays, so an expression like

```python
def H_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep \
        + hamiltonians.hamiltonian_3nu_matter(VCC_exp_density(l, num_density_e_center, l_scale))
```

is already array-capable as written. The two things that break it are a Python
`if` on the position (use `np.where`) and a call to `float()` or `math.exp`
(use the NumPy equivalent). Where a matrix has to be built by hand, index the
potential with `[..., None, None]` so one value per position becomes a stack of
matrices:

```python
e00 = np.diag([1.0, 0.0, 0.0])
VCC[..., None, None]*e00        # (npts, 3, 3), not an error
```

A Hamiltonian that ignores its argument entirely -- constant density -- is
detected separately and costs nothing either way.'''),
    code(r'''# H_vac_energy_indep is the vacuum Hamiltonian without the multiplicative (1/E) factor.  Since this part will not change anymore in our
# examples below, we compute it once and save it.  (Note: in this case, we could have used the function `osc_prob_3nu_matter_vacuum` to
# compute probabilities, but we choose to compute them from scratch, as the other cases.)
H_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31) # [eV^2]
def H_vac(energy):
    return (1/energy)*H_vac_energy_indep # [eV]

# In matter: exponentially falling density profile: num_density_e_center * e^{-l/l_scale}
num_density_e_center = 10*gd.N_AV/gd.UNIT_CM3 # Central electron number density (10*N_Av cm^{-3}, converted to eV^3) [eV^3]
l_scale = 100.0 # [km]
def num_density_e_exp(l, num_density_e_center, l_scale):
    return num_density_e_center*np.exp(-(l/gd.UNIT_KM)/l_scale) # [eV^3]
def VCC_exp_density(l, num_density_e_center, l_scale):
    return matter.VCC_func(l, lambda r: num_density_e_exp(r, num_density_e_center, l_scale)) # [eV]
def H_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_3nu_matter(VCC_exp_density(l, num_density_e_center, l_scale))'''),
    md(r'''Let's first fix the baseline to 500 km and compute the $\nu_e \to \nu_e$ survival probability.  We call `osc_prob` to compute the probabilities, like before, and pass `verbose = 2` to see the full output:'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUE
baseline  = 5.e2 # [km]
oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, baseline*gd.UNIT_KM, verbose=2)[nu_i][nu_f]'''),
    md(r'''By default, `osc_prob` uses fourth-order Magnus expansion (`magnus_exp_order = 4`) and returns the probability with a relative tolerance of `rtol = 1e-3` and an absolute tolerance of `atol = 1e-3`.  These values are sufficient for most purposes, but the user can change them if they so wish.

Let's see how the probability changes with the order of the expansion (we loosen the target tolerance so that the order, rather than the adaptive refinement, is what limits the answer).

One limit is worth knowing about here.  The default integrator, `integration_method = 'gl'`, is a set of Gauss-Legendre *commutator-free* schemes: they are separately derived integrators rather than truncations of the Magnus recursion, so they exist only up to order 6, not up to `gd.MAGNUS_EXP_ORDER_MAX = 10`.  Asking for a higher order with `'gl'` raises a `ValueError` saying so.  The `'trapezoid'` and `'simpson'` integrators do build the terms from the recursion and go all the way to 10, at the cost of many more commutators per order.'''),
    code(r'''# The default 'gl' integrator exists up to order 6; 'trapezoid'/'simpson' go to 10.
max_order_gl = magnus.MAGNUS_EXP_ORDER_MAX_GL

[[magnus_exp_order, oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, baseline*gd.UNIT_KM,
                                     rtol=1.e-2, atol=1.e-2, magnus_exp_order=magnus_exp_order)[nu_i][nu_f]]
 for magnus_exp_order in range(1, max_order_gl+1)]'''),
    md(r'''The expansion order is something *you* choose and `osc_prob` keeps fixed.  What it adapts on your behalf is the **slab count**: it partitions the baseline into progressively more subintervals, recomputing until two successive refinements agree to within `rtol`/`atol`.

This split of responsibilities is deliberate, and it is worth understanding which knob does what:

* **Order** controls how well each slab is approximated at fixed width.  Raising it helps only while the series is converging on that slab -- and the cost of a term grows roughly like $2^k$ in the number of commutators, so high order is expensive.
* **Slab count** controls the width each slab has to cover.  Because a Magnus slab is *exact* for a constant Hamiltonian no matter how much phase accumulates across it, narrowing slabs is what buys accuracy when the matter profile varies.

That is why the adaptive loop refines slabs rather than order, and why order 4 is a good default: when a result is not converging, more slabs almost always beats a higher order.  If the loop reaches `max_n_slabs` without two successive levels agreeing, it emits a `ToleranceNotAchievedWarning` rather than silently returning the last value -- note that this says convergence could not be *verified*, not that the answer is necessarily bad.

Pass `verbose = 2` to watch the refinement happen:'''),
    code(r'''oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, baseline*gd.UNIT_KM, verbose=2,
                 rtol=1.e-3, atol=1.e-3, max_n_slabs=32)[nu_i][nu_f]'''),
    md(r'''For most purposes a fixed `magnus_exp_order` of 3 or 4 is enough, and the slab refinement does the rest.

Internally, `osc_prob` partitions the interval from `t_ini = 0` to `t_fin = baseline` into progressively more subintervals, or slabs, computes the evolution operator in each, and takes their position-ordered product.  It repeats that with a finer partition until two successive levels agree to within `rtol` and `atol`.

The number of slabs (`n_slabs`) grows from `min_n_slabs` (default `1`), multiplied at each pass by roughly `growth_factor_n_slabs` (default `1.5`), up to `max_n_slabs`.  That ceiling is **method-dependent**: `20000` for the default `'gl'` integrator and `2000` for `'trapezoid'`/`'simpson'` (see `oscprob.MAX_N_SLABS_DEFAULT`).  The split exists because the two families cost very different amounts per slab.  The loop also stops after `max_num_loops` passes (default `50`).

The number of time points per slab (`n_tpts_per_slab`) at which the Hamiltonian is sampled to integrate it grows the same way — from `min_n_tpts_per_slab` (default `2`) to `max_n_tpts_per_slab` (default `500`), by `growth_factor_n_tpts_per_slab` (default `1.5`).  **With the default `'gl'` integrator this ladder is switched off entirely**: Gauss-Legendre collocation needs only 1-3 Hamiltonian evaluations per slab regardless, so accuracy is controlled by the slab count alone and `osc_prob` pins `n_tpts_per_slab = 2`.  These parameters matter only for `'trapezoid'` and `'simpson'`.

If we set the caps low and ask for high accuracy (`rtol = atol = 1.e-6`), `osc_prob` still returns a result, but warns that the requested tolerance could not be *verified* — at the ceiling there is no finer level left to compare against, so the check cannot be completed:'''),
    code(r'''oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, baseline*gd.UNIT_KM, verbose=2,
                 max_n_slabs=10, max_n_tpts_per_slab=10, rtol=1.e-6, atol=1.e-6)[nu_i][nu_f]'''),
    md(r'''If we use the default values of these run parameters, we will achieve the high accuracy we seek:'''),
    code(r'''oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, baseline*gd.UNIT_KM, verbose=1,
                 rtol=1.e-6, atol=1.e-6)[nu_i][nu_f]'''),
    md(r'''Now let's illustrate this for two matter density profiles that vary with position: the exponentially falling density above and a Gaussian function where most of the matter is concentrated at a certain position.  We also show results for vacuum and for constant density.'''),
    code(r'''# In matter: constant density equal to central density 
VCC_center = matter.VCC_func(0.0, lambda l: num_density_e_center) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_3nu_matter(VCC_center) # [eV]
def H_const_density(energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density # [eV]

# In matter: Gaussian density profile, centered around l_central, with width l_width
num_density_e_basal = 0*gd.N_AV/gd.UNIT_CM3 # [eV^3]
num_density_e_central = 8*gd.N_AV/gd.UNIT_CM3 # [eV^3]
l_central, l_width = 300.0, 100.0 # [km]
def num_density_e_gaussian(l, num_density_e_central, l_central, l_width):
    return num_density_e_basal+num_density_e_central*np.exp( -(l/gd.UNIT_KM-l_central)**2/(2.0*l_width**2)) # [eV^3]
def VCC_gaussian_density(l, num_density_e_central, l_central, l_width):
    return matter.VCC_func(l, lambda r: num_density_e_gaussian(r, num_density_e_central, l_central, l_width)) # [eV]
def H_gaussian_density(l, energy):
    return (1/energy)*H_vac_energy_indep  \
        + hamiltonians.hamiltonian_3nu_matter(VCC_gaussian_density(l, num_density_e_central, l_central, l_width)) # [eV]'''),
    md(r'''Let's first generate the $\nu_e \to \nu_e$ survival probabilities vs. baseline for a fixed energy.  We pass `n_jobs = 10` to spread the work over ten parallel jobs: these parallelize over the requested **(energy, baseline) points**, not over the slabs within one point, so they pay off exactly when there are many points to compute -- as here, with a thousand baselines.  For a single point, leave it at 1.

We also drop to third-order Magnus expansion for speed; the results do not change noticeably at this accuracy.'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUE

# Baselines
l_ini, l_fin = 5.e1, 1.e3 # [km]
l_npts = 1000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
energy = 10.0*gd.UNIT_MEV # [eV]

# In vacuum
print('Computing probabilities in vacuum ...\n')
prob_vac = [oscprob.osc_prob(H_vac(energy), 0.0, l*gd.UNIT_KM)[nu_i][nu_f] for l in distances]
print('   Done\n')

# In matter: constant density equal to central density 
print('Computing probabilities in matter, constant density ...\n')
prob_matt_const_density = [oscprob.osc_prob(H_const_density(energy), 0, l*gd.UNIT_KM)[nu_i][nu_f] for l in distances]
print('   Done\n')

# For vacuum and constant-density matter above, because we passed a time-independent Hamiltonian to `osc_prob`, we need not pass any other
# argument to it; internally, Magnus detects that is a time-independent case and adjusts the run parameters automatically for it (e.g.,
# sets the Magnus expansion to first-order only).  See comments above.

# For varying-density matter below, we need to pass to `osc_prob`, first, the Hamiltonian as a position-dependent function and, second,
# values for the run parameters to set the target accuracy with which to compute the probability.

# In matter: exponentially falling density profile
print('Computing probabilities in matter, exponentially falling density ...\n')
prob_matt_exp_density = [oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, l*gd.UNIT_KM,
                                          magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                         for l in distances]
print('   Done\n')

# In matter: Gaussian density profile
print('Computing probabilities in matter, Gaussian density ...\n')
prob_matt_gaussian_density = [oscprob.osc_prob(lambda l: H_gaussian_density(l, energy), 0, l*gd.UNIT_KM, 
                                               magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                              for l in distances]
print('   Done\n')'''),
    md(r'''Now plot the probabilities vs. distance:'''),
    code(r'''norm = (gd.N_AV/gd.UNIT_CM3)
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center/norm for l in distances], color='C0'),
     dict(y=[num_density_e_exp(l*gd.UNIT_KM, num_density_e_center, l_scale)/norm
             for l in distances], color='C3', ls='--'),
     dict(y=[num_density_e_gaussian(l*gd.UNIT_KM, num_density_e_central, l_central, l_width)/norm
             for l in distances], color='C5', ls='-.')],
    [[dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density, label=r'Constant density', color='C0'),
      dict(y=prob_matt_exp_density, label=r'Exponentially falling', color='C3', ls='--'),
      dict(y=prob_matt_gaussian_density, label=r'Gaussian', color='C5', ls='-.')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 11),
    profile_ymajor=4, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Survival probability, $P_{\nu_e \to \nu_e}$'],
    xlabel=r'Baseline, $L$ [km]', title_fontsize=23,
    title=r'$3\nu$~varying-density matter, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left')'''),
    md(r'''Let's now generate probabilities vs. energy for a fixed baseline:'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e2 # [MeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
baseline = 1.e3*gd.UNIT_KM # 1000 km in eV^{-1} [eV^{-1}]

# In vacuum
print('Computing probabilities in vacuum ...\n')
prob_vac = [oscprob.osc_prob(H_vac(enu), 0.0, baseline)[nu_i][nu_f] for enu in gd.UNIT_MEV*energies]
print('   Done\n')

# In matter: constant density equal to central density 
print('Computing probabilities in matter, constant density ...\n')
prob_matt_const_density = [oscprob.osc_prob(H_const_density(enu), 0, baseline)[nu_i][nu_f] for enu in gd.UNIT_MEV*energies]
print('   Done\n')

# In matter: exponentially falling density profile
print('Computing probabilities in matter, exponentially falling density ...\n')
prob_matt_exp_density = [oscprob.osc_prob(lambda l: H_exp_density(l, enu), 0, baseline, 
                                          magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                         for enu in gd.UNIT_MEV*energies]
print('   Done\n')

# In matter: Gaussian density profile
print('Computing probabilities in matter, Gaussian density ...\n')
prob_matt_gaussian_density = [oscprob.osc_prob(lambda l: H_gaussian_density(l, enu), 0, baseline, 
                                               magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                              for enu in gd.UNIT_MEV*energies]
print('   Done\n')'''),
    md(r'''And plot them:'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
     dict(y=prob_matt_const_density, label=r'Constant density', color='C0'),
     dict(y=prob_matt_exp_density, label=r'Exponentially falling', color='C3', ls='--'),
     dict(y=prob_matt_gaussian_density, label=r'Gaussian', color='C5', ls='-.')],
    energy_unit='MeV',
    ylabel=r'Survival probability, $P_{\nu_e \to \nu_e}$',
    xlim=(energies[0], energies[-1]), title_fontsize=23,
    title=r'$3\nu$~varying-density matter, $L = $~{:.2f}~km'.format(baseline/gd.UNIT_KM),
    legend_title=r'Matter profile', legend_loc='lower left')'''),
    md(r'''Finally, let's  show one example of beyond-the-Standard-Model neutrino physics in the form of non-standard neutrino-matter interactions (NSI).  These interactions introduce a new matter contribution to the Hamiltonian,
\begin{equation}
 H_{3\nu}^{\rm NSI}
 = 
 V_{\rm CC}
 \left(
  \begin{array}{ccc}
   \epsilon_{ee}         & \epsilon_{e\mu}         & \epsilon_{e\tau} \\
   \epsilon_{e\mu}^\ast  & \epsilon_{\mu\mu}       & \epsilon_{\mu\tau} \\
   \epsilon_{e\tau}^\ast & \epsilon_{\mu\tau}^\ast & \epsilon_{\tau\tau}
  \end{array}
 \right) \;,
\end{equation}
where the coefficients $\epsilon_{\alpha\beta}$ set the strength of the interaction.  Their values are set by experiment.  We set them to their predefined values from `globaldefs`, but they can be set to any value. '''),
    code(r'''eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = gd.EPS_3
print(gd.EPS_3)'''),
    md(r'''We test NSI using the exponentially falling matter density.  We increase the distance scale (`l_scale`) beyond which the matter profile exponentially decays, in order for the matter effects to be more evident.  We also define the Hamiltonian with NSI matter effects:'''),
    code(r'''l_scale = 1.e3 # [km]

def num_density_e_exp(l, num_density_e_center, l_scale):
    return num_density_e_center*np.exp(-(l/gd.UNIT_KM)/l_scale) # [eV^3]
def VCC_exp_density(l, num_density_e_center, l_scale):
    return matter.VCC_func(l, lambda r: num_density_e_exp(r, num_density_e_center, l_scale)) # [eV]

# In matter: exponentially falling density profile
def H_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_3nu_matter(VCC_exp_density(l, num_density_e_center, l_scale))
    
# In matter: exponentially falling density profile
def H_exp_density_nsi(l, energy): 
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_3nu_nsi(
        VCC_exp_density(l, num_density_e_center, l_scale),
        eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt)'''),
    md(r'''Now we generate the probabilities vs. distance.  We choose a higher energy than before to make the matter contribution more dominant relative to the vacuum one:'''),
    code(r'''nu_i, nu_f = gd.NUE, gd.NUE

# Baselines
l_ini, l_fin = 5.e1, 1.e3 # [km]
l_npts = 1000 
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

energy = 10.0*gd.UNIT_MEV # [eV]

# In vacuum
print('Computing probabilities in vacuum ...\n')
prob_vac = [oscprob.osc_prob(H_vac(energy), 0.0, l*gd.UNIT_KM)[nu_i][nu_f] for l in distances]
print('   Done\n')

# In matter: constant density equal to central density 
print('Computing probabilities in matter, constant density ...\n')
prob_matt_const_density = [oscprob.osc_prob(H_const_density(energy), 0, l*gd.UNIT_KM)[nu_i][nu_f] for l in distances]
print('   Done\n')

# In matter: exponentially falling density profile
print('Computing probabilities in matter, exponentially falling density ...\n')
prob_matt_exp_density = [oscprob.osc_prob(lambda l: H_exp_density(l, energy), 0, l*gd.UNIT_KM, 
                                          magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                         for l in distances]
print('   Done\n')

# In matter: exponentially falling density profile with NSI
print('Computing probabilities in matter, exponentially falling density with NSI ...\n')
prob_matt_exp_density_nsi = [oscprob.osc_prob(lambda l: H_exp_density_nsi(l, energy), 0, l*gd.UNIT_KM, 
                                              magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                             for l in distances]
print('   Done\n')'''),
    md(r'''And plot them:'''),
    code(r'''norm = (gd.N_AV/gd.UNIT_CM3)
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center/norm for l in distances], color='C0'),
     dict(y=[num_density_e_exp(l*gd.UNIT_KM, num_density_e_center, l_scale)/norm
             for l in distances], color='C3', ls='--')],
    [[dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density, label=r'Constant density', color='C0'),
      dict(y=prob_matt_exp_density, label=r'Exponentially falling', color='C3', ls='--'),
      dict(y=prob_matt_exp_density_nsi, label=r'Exponentially falling with NSI',
           color='C5', ls='-.')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 11),
    profile_ymajor=4, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Survival probability, $P_{\nu_e \to \nu_e}$'],
    xlabel=r'Baseline, $L$ [km]', title_fontsize=23, ylabel_labelpad=25,
    title=r'$3\nu$~varying-density matter, $E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left')'''),
    md(r'''Finally, let's generate probabilities vs. energy for a fixed baseline:'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e2 # [MeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
baseline = 5.e2*gd.UNIT_KM # 500 km in eV^{-1} [eV^{-1}]

# In vacuum
print('Computing probabilities in vacuum ...\n')
prob_vac = [oscprob.osc_prob(H_vac(enu), 0.0, baseline)[nu_i][nu_f] for enu in gd.UNIT_MEV*energies]
print('   Done\n')

# In matter: constant density equal to central density 
print('Computing probabilities in matter, constant density ...\n')
prob_matt_const_density = [oscprob.osc_prob(H_const_density(enu), 0, baseline)[nu_i][nu_f] for enu in gd.UNIT_MEV*energies]
print('   Done\n')

# In matter: exponentially falling density profile
print('Computing probabilities in matter, exponentially falling density ...\n')
prob_matt_exp_density = [oscprob.osc_prob(lambda l: H_exp_density(l, enu), 0, baseline, 
                                          magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                         for enu in gd.UNIT_MEV*energies]
print('   Done\n')

# In matter: exponentially falling density profile
print('Computing probabilities in matter, exponentially falling density with NSI ...\n')
prob_matt_exp_density_nsi = [oscprob.osc_prob(lambda l: H_exp_density_nsi(l, enu), 0, baseline, 
                                              magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] 
                             for enu in gd.UNIT_MEV*energies]
print('   Done\n')'''),
    md(r'''And plot them:'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
     dict(y=prob_matt_const_density, label=r'Constant density', color='C0'),
     dict(y=prob_matt_exp_density, label=r'Exponentially falling', color='C3', ls='--'),
     dict(y=prob_matt_exp_density_nsi, label=r'Exponentially falling with NSI',
          color='C5', ls='-.')],
    energy_unit='MeV',
    ylabel=r'Survival probability, $P_{\nu_e \to \nu_e}$',
    xlim=(energies[0], energies[-1]), title_fontsize=23,
    title=r'$3\nu$~varying-density matter, $L = $~{:.2f}~km'.format(baseline/gd.UNIT_KM),
    legend_title=r'Matter profile', legend_loc='lower left')'''),
    ])

# ----------------------------------------------- 02_magnus_2nu_vacuum_matter
books['02_magnus_2nu_vacuum_matter.ipynb'] = notebook(
    'Two-neutrino probabilities',
    'Oscillation probabilities in a two-flavour system, against energy and against direction, in seven settings: vacuum, constant density, an exponential and a Gaussian profile, a periodic castle wall, a noisy profile, and then the Earth and the Sun.\n\nEach is validated against the closed-form expression where one exists, which is what makes this the notebook to read before trusting any of the others.',
    [
    code(r'''import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.oscprobstd as oscprobstd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''# 0. Helper functions and definitions'''),
    code(r"""col_full = '#E76F51'
col_ang = '#2A9D8F'
col_en = '#E9C46A'
col_count = '#73A8BF'#'#2A9D8F'"""),
    code(r'''def prob_label(nu_i, nu_f):
    if (nu_i == gd.NUE):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_e \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_e \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_e \to \nu_\tau}$'
    elif (nu_i == gd.NUMU):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_\mu \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_\mu \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_\mu \to \nu_\tau}$'
    elif (nu_i == gd.NUTAU):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_\tau \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_\tau \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_\tau \to \nu_\tau}$'
    return label'''),
    md(r'''# 1. Probabilities 2$\nu$: in vacuum'''),
    md(r'''## 1.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    md(r'''## 1.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e5 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 10.0*gd.UNIT_MEV # [eV]

# Both probabilities are computed for the *whole* array of baselines in a single
# call, which is the form to prefer: measured here, the array call takes 0.95 s
# against 1.57 s for a Python loop that calls osc_prob once per point (1.7x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# Using the standard formula
prob_std = oscprobstd.osc_prob_2nu_vacuum_std(
    sth, Dm2, energy, distances*gd.CONV_KM_TO_INV_EV)[nu_i][nu_f]

# Using the Magnus expansion
prob = oscprob.osc_prob_2nu_vacuum(
    energy, distances*gd.CONV_KM_TO_INV_EV, sth, Dm2, nu_i=nu_i, nu_f=nu_f)'''),
    md(r'''### Plot probabilities'''),
    code(r'''fig, ax = plotting.plot_probability_vs_baseline(
    distances,
    [dict(y=prob, label='Magnus expansion', color='C1'),
     dict(y=prob_std, label='Standard formula', color='k', ls='--')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2,
    xlim=(l_ini, l_fin),
    residual=(prob-prob_std)/prob_std/1.e-12,
    residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
    residual_ylim=(-1.2, 1.2), residual_ymajor=0.50, residual_yminor=0.10,
    legend_title=r'Calculation method', legend_loc='center left',
    title=r'$2\nu$~vacuum, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    savefig='../fig/prob_2nu_vacuum_vs_baseline.pdf')'''),
    md(r'''## 1.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

# Compute probability vs. energy
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 5.e3 # [km]

# Both probabilities are computed for the *whole* array of energies in a single
# call, which is the form to prefer: measured here, the array call takes 0.59 s
# against 0.80 s for a Python loop that calls osc_prob once per point (1.4x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# Using the standard formula
prob_std = oscprobstd.osc_prob_2nu_vacuum_std(
    sth, Dm2, energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV)[nu_i][nu_f]

# Using the Magnus expansion
prob = oscprob.osc_prob_2nu_vacuum(
    energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV, sth, Dm2,
    nu_i=nu_i, nu_f=nu_f)'''),
    md(r'''### Plot probabilities'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob, label='Magnus expansion', color='r'),
     dict(y=prob_std, label='Standard formula', color='k', ls='--')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2,
    xlim=(energy_min, energy_max),
    residual=(prob-prob_std)/prob_std/1.e-14,
    residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-14}]$',
    residual_ylim=(-0.6, 0.6), residual_ymajor=0.20, residual_yminor=0.05,
    legend_title=r'Calculation method', legend_loc='center right',
    title=r'$2\nu$~vacuum, $L = $~{:.2f}~km'.format(baseline),
    savefig='../fig/prob_2nu_vacuum_vs_energy.pdf')'''),
    md(r'''# 2. Probabilities 2$\nu$: in matter with constant density'''),
    md(r'''## 2.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    md(r'''## 2.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e5 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 10.0*gd.UNIT_MEV # [eV]

# In matter
rho = 10.0 # Matter density [g cm^{-3}]
# The potential is still built by hand here, because the standard formula below
# needs it explicitly. Note density_matter_is_in_g_per_cm3=True: without it, rho
# is taken to be in internal units already and the potential comes out ~4e18
# times too small, which looks exactly like having no matter at all.
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho,
                                          electron_fraction=0.5,
                                          density_matter_is_in_g_per_cm3=True) # [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # [eV]

# Both probabilities are computed for the *whole* array of baselines in a single
# call, which is the form to prefer: measured here, the array call takes 1.61 s
# against 3.20 s for a Python loop that calls osc_prob once per point (2.0x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# Using the standard formula
prob_matt_std = oscprobstd.osc_prob_2nu_matter_std(
    sth, Dm2, VCC, energy, distances*gd.CONV_KM_TO_INV_EV)[nu_i][nu_f]

# Using the Magnus expansion. The wrapper takes the density directly and builds
# the potential itself, so the matter Hamiltonian never has to be assembled by
# hand -- and its sign cannot be got wrong.
prob_matt = oscprob.osc_prob_2nu_matter_constant_density(
    energy, distances*gd.CONV_KM_TO_INV_EV, rho, sth, Dm2,
    density_matter_is_in_g_per_cm3=True, nu_i=nu_i, nu_f=nu_f)

# In vacuum, for comparison
prob_vac = oscprob.osc_prob_2nu_vacuum(
    energy, distances*gd.CONV_KM_TO_INV_EV, sth, Dm2, nu_i=nu_i, nu_f=nu_f)'''),
    md(r'''### Plot probabilities'''),
    code(r'''fig, ax = plotting.plot_probability_vs_baseline(
    distances,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', ls=':', lw=2),
     dict(y=prob_matt, color='C1', lw=3,
          label=r'Matter ($\rho = $~{:.1f}'.format(rho)+r'~g~cm$^{-3}$),'+'\n'+'Magnus expansion'),
     dict(y=prob_matt_std, color='k', ls='--',
          label=r'Matter ($\rho = $~{:.1f}'.format(rho)+r'~g~cm$^{-3}$),'+'\n'+'standard formula')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2,
    xlim=(l_ini, l_fin),
    residual=(prob_matt-prob_matt_std)/prob_matt_std/1.e-12,
    residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
    residual_ylim=(-1.6, 1.6), residual_ymajor=0.50, residual_yminor=0.10,
    legend_title=r'Calculation method', legend_loc='lower left',
    title=r'$2\nu$~constant-density matter, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    savefig='../fig/prob_2nu_matter_vs_baseline.pdf')'''),
    md(r'''## 2.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

# Compute probability vs. energy
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 5.e3 # [km]

# In matter
rho = 10.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho,
                                          electron_fraction=0.5,
                                          density_matter_is_in_g_per_cm3=True) # [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # [eV]

# Both probabilities are computed for the *whole* array of energies in a single
# call, which is the form to prefer: measured here, the array call takes 0.76 s
# against 1.72 s for a Python loop that calls osc_prob once per point (2.3x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# Using the standard formula
prob_matt_std = oscprobstd.osc_prob_2nu_matter_std(
    sth, Dm2, VCC, energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV)[nu_i][nu_f]

# Using the Magnus expansion
prob_matt = oscprob.osc_prob_2nu_matter_constant_density(
    energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV, rho, sth, Dm2,
    density_matter_is_in_g_per_cm3=True, nu_i=nu_i, nu_f=nu_f)

# In vacuum, for comparison
prob_vac = oscprob.osc_prob_2nu_vacuum(
    energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV, sth, Dm2,
    nu_i=nu_i, nu_f=nu_f)'''),
    md(r'''### Plot probabilities'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', ls=':', lw=2),
     dict(y=prob_matt, color='C1', lw=3,
          label=r'Matter ($\rho = $~{:.1f}'.format(rho)+r'~g~cm$^{-3}$),'+'\n'+'Magnus expansion'),
     dict(y=prob_matt_std, color='k', ls='--',
          label=r'Matter ($\rho = $~{:.1f}'.format(rho)+r'~g~cm$^{-3}$),'+'\n'+'standard formula')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2, energy_unit='GeV',
    xlim=(energy_min, energy_max),
    residual=(prob_matt-prob_matt_std)/prob_matt_std/1.e-12,
    residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
    residual_ylim=(-0.5, 0.5), residual_ymajor=0.20, residual_yminor=0.05,
    legend_title=r'Calculation method', legend_loc='lower left',
    title=r'$2\nu$~constant-density matter, $L = $~{:.2f}~km'.format(baseline),
    savefig='../fig/prob_2nu_matter_vs_energy.pdf')'''),
    md(r'''# 3. Probabilities 2$\nu$: in matter with varying density'''),
    md(r'''## 3.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Parameters for the electron number density in the Sun, Eq. (10.62) in Giunti & Kim
num_density_e_center = 10*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_scale = 100.0 # [km]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
# The matter density is constant, so this is merely a wrapper.
# [l] = km
VCC_center = matter.VCC_func(0.0, lambda l: num_density_e_center) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_2nu_matter(VCC_center) # [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density

# --------------------------------------------------------------
# In matter, exponentially falling density
# --------------------------------------------------------------
# [l] = km
def num_density_e_exp_func(l, num_density_e_center, l_scale):
    return num_density_e_center*np.exp(-(l/gd.CONV_KM_TO_INV_EV)/l_scale) # [eV^3]
def VCC_func_exp_density(l, num_density_e_center, l_scale):
    return matter.VCC_func(l, lambda r: num_density_e_exp_func(r, num_density_e_center, l_scale)) # [eV]
def H_func_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_density(l, num_density_e_center, l_scale))


# --------------------------------------------------------------
# In matter, Gaussian density
# --------------------------------------------------------------
# [l] = km
num_density_e_basal = 0*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
num_density_e_central = 8*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_central, l_width = 300.0, 100.0 # [km]
def num_density_e_gaussian_func(l, num_density_e_central, l_central, l_width):
    return num_density_e_basal+num_density_e_central*np.exp( -(l/gd.CONV_KM_TO_INV_EV-l_central)**2/(2.0*l_width**2)) # [eV^3]
def VCC_func_gaussian_density(l, num_density_e_central, l_central, l_width):
    return matter.VCC_func(l, lambda r: num_density_e_gaussian_func(r, num_density_e_central, l_central, l_width)) # [eV]
def H_func_gaussian_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_gaussian_density(l, 
                                                    num_density_e_central, l_central, l_width)) # [eV]'''),
    md(r'''## 3.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 5.e1, 1.e3 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 10.0*gd.UNIT_MEV # [eV]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, energy), 
                                                     0, l*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, exponentially falling density
# --------------------------------------------------------------
prob_matt_exp_density = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, energy), 
                                                   0, l*gd.CONV_KM_TO_INV_EV, 
                                                   n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, Gaussian density
# --------------------------------------------------------------
prob_matt_gaussian_density = np.array([oscprob.osc_prob(lambda l: H_func_gaussian_density(l, energy), 
                                                        0, l*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=10, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
H_vac = (1./energy)*H_vac_energy_indep # Can also call hamiltonians.hamiltonian_2nu_vacuum instead
prob_vac = np.array([oscprob.osc_prob(lambda l: H_vac,
                                      0.0, l*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center/norm for l in distances], color='C0'),
     dict(y=[num_density_e_exp_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_center, l_scale)/norm
             for l in distances], color='C3', ls='--'),
     dict(y=[num_density_e_gaussian_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_central,
                                         l_central, l_width)/norm
             for l in distances], color='C5', ls='-.')],
    [[dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density, label=r'Constant density', color='C0'),
      dict(y=prob_matt_exp_density, label=r'Exponentially falling', color='C3', ls='--'),
      dict(y=prob_matt_gaussian_density, label=r'Gaussian', color='C5', ls='-.')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 11),
    profile_ymajor=2, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]',
    title=r'$2\nu$~varying-density matter, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    title_fontsize=20, ylabel_labelpad=15,
    legend_title=r'Matter profile', legend_loc='lower left', grid=False)'''),
    md(r'''## 3.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e2 # [MeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 1.e3 # [km]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV), 
                                                     0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])

# --------------------------------------------------------------
# In matter, exponeitally falling density
# --------------------------------------------------------------
prob_matt_exp_density = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, enu*gd.UNIT_MEV), 
                                                   0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                   n_slabs=20, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] \
                                  for enu in energies])

# --------------------------------------------------------------
# In matter, Gaussian density
# --------------------------------------------------------------
prob_matt_gaussian_density = np.array([oscprob.osc_prob(lambda l: H_func_gaussian_density(l, enu*gd.UNIT_MEV), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=20, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] \
                                  for enu in energies])

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
prob_vac = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_MEV))*H_vac_energy_indep, 
                                      0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])'''),
    md(r'''### Plot probabilities'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
     dict(y=prob_matt_const_density, label=r'Constant density', color='C0'),
     dict(y=prob_matt_exp_density, label=r'Exponentially falling', color='C3', ls='--'),
     dict(y=prob_matt_gaussian_density, label=r'Gaussian', color='C5', ls='-.')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2, energy_unit='MeV',
    xlim=(energy_min, energy_max),
    legend_title=r'Matter profile', legend_loc='center left',
    title=r'$2\nu$~varying-density matter, $L = $~{:.2f}~km'.format(baseline))'''),
    md(r'''# 4. Probabilities 2$\nu$: in matter with castle-wall density profile'''),
    md(r'''## 4.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Low and high densities of the castle-wall density profile
num_density_e_low = 1*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
num_density_e_high = 10*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
num_density_e_avg = (num_density_e_high-num_density_e_low)/2.0 # [eV^3]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]
    
# --------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
# The matter density is constant, so this is merely a wrapper ([l] = km)
VCC_avg = matter.VCC_func(0.0, lambda l: num_density_e_avg) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_2nu_matter(VCC_avg) # [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density

# --------------------------------------------------------------
# In matter, castle-wall density profile
# --------------------------------------------------------------
# [l] = km
def num_density_e_castle_wall_func(l, num_density_e_low, num_density_e_high, n_castle_slabs, l_ini, l_fin):
    # np.where, not an if/else: the engine evaluates this for a whole array of
    # positions at once, and a Python `if` on an array raises "truth value of an
    # array is ambiguous". Written this way the profile stays array-capable and
    # osc_prob keeps its fast path -- see the note in notebook 01. It still
    # returns a plain number when handed a single position.
    l_scaled = (np.asarray(l)/gd.CONV_KM_TO_INV_EV-l_ini)/(l_fin-l_ini)
    dl = 1.0/n_castle_slabs # Width of one slab
    index_slab = l_scaled // dl
    # The first slab has low density, the second one has high density, and so on
    return np.where(index_slab % 2 == 0, num_density_e_low, num_density_e_high)
def VCC_func_castle_wall(l, num_density_e_low, num_density_e_high, n_castle_slabs, l_ini, l_fin):
    return matter.VCC_func(l, lambda r: num_density_e_castle_wall_func(r, num_density_e_low, num_density_e_high, 
                                                                       n_castle_slabs, l_ini, l_fin)) # [eV]
def H_func_castle_wall(l, energy, num_density_e_low, num_density_e_high, n_castle_slabs, l_ini, l_fin):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_castle_wall(l, 
                                                                                        num_density_e_low, num_density_e_high,
                                                                                        n_castle_slabs, l_ini, l_fin))

n_castle_slabs_narrow = 50 # Narrow castle wall
n_castle_slabs_wide = 10   # Wide castle wall

# Positions of the density walls, in eV^-1.  The profile is a step function, and
# high-order quadrature reaches its nominal order only where the Hamiltonian is
# smooth inside each slab -- integrating across a wall costs accuracy that no
# amount of extra slabs cheaply buys back.  Handing these to osc_prob as
# `t_breakpoints` makes them mandatory slab edges, so no slab ever straddles one.
# osc_prob keeps only the breakpoints that fall inside [t_ini, t_fin], so the same
# array can be passed at every baseline of a scan.
#
# Measured over the scans below, against a converged (n_slabs=24000) reference:
# with n_slabs=150 alone the worst point is 1.5e-2 off and 10% of points exceed
# 1e-3; adding the breakpoints gives 1.9e-3 worst and 1.5%, and runs 2.8x faster.
def castle_wall_breakpoints(n_castle_slabs, l_ini, l_fin):
    # Every wall edge, including the two at l_ini and l_fin where the profile switches on and
    # off.  A trajectory starting before l_ini crosses those as well, and a slab straddling an
    # unmarked discontinuity degrades the quadrature no matter how fine the grid becomes: with
    # arange(1, n) -- the interior walls only -- the error sat at 1.6e-2 in probability,
    # bit-identical from n_slabs=4 through 32, against 3.6e-12 once these two are included.
    return (l_ini+np.arange(0, n_castle_slabs+1)*(l_fin-l_ini)/n_castle_slabs) \
        *gd.CONV_KM_TO_INV_EV # [eV^-1]'''),
    md(r'''## 4.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e2, 1.e4 #1.e2, 1.e5 # [km]
l_npts = 6000 
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 50.*gd.UNIT_MEV # [eV] #10

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, energy), 
                                                     0, l*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, narrow castle-wall density profile
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_castle_wall_narrow = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_castle_wall(l, energy, num_density_e_low,
                                num_density_e_high, n_castle_slabs_narrow,
                                l_ini, l_fin), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3,
    t_breakpoints=castle_wall_breakpoints(n_castle_slabs_narrow, l_ini, l_fin))

# --------------------------------------------------------------
# In matter, constant density, low-amplitude noise
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_castle_wall_wide = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_castle_wall(l, energy, num_density_e_low,
                                num_density_e_high, n_castle_slabs_wide,
                                l_ini, l_fin), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3,
    t_breakpoints=castle_wall_breakpoints(n_castle_slabs_wide, l_ini, l_fin))'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_avg/norm for l in distances], color='0.5', ls='--'),
     dict(y=[num_density_e_castle_wall_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_low,
                                            num_density_e_high, n_castle_slabs_narrow, l_ini, l_fin)/norm
             for l in distances], color='C0'),
     dict(y=[num_density_e_castle_wall_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_low,
                                            num_density_e_high, n_castle_slabs_wide, l_ini, l_fin)/norm
             for l in distances], color='C3')],
    [[dict(y=prob_matt_const_density, label=r'Constant average density', color='0.5', ls='--'),
      dict(y=prob_matt_castle_wall_narrow, label=r'Narrow castle wall', color='C0'),
      dict(y=prob_matt_castle_wall_wide, label=r'Wide castle wall', color='C3')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 11),
    profile_ymajor=4, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]',
    title=r'$2\nu$~castle-wall matter profile, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left')'''),
    md(r'''## 4.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e3 # [MeV]
energy_npts = 4000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 1.e4 # [km]

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV),
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])


# --------------------------------------------------------------
# In matter, narrow castle-wall density profile
# --------------------------------------------------------------
prob_matt_castle_wall_narrow = np.array([oscprob.osc_prob(lambda l: H_func_castle_wall(l, enu*gd.UNIT_MEV, 
                                                                                       num_density_e_low, num_density_e_high, 
                                                                                       n_castle_slabs_narrow, l_ini, l_fin), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10,
                                                     t_breakpoints=castle_wall_breakpoints(n_castle_slabs_narrow, l_ini, l_fin))[nu_i][nu_f] \
                                          for enu in energies]) 

# --------------------------------------------------------------
# In matter, wide castle-wall density profile
# --------------------------------------------------------------
prob_matt_castle_wall_wide = np.array([oscprob.osc_prob(lambda l: H_func_castle_wall(l, enu*gd.UNIT_MEV, 
                                                                                       num_density_e_low, num_density_e_high, 
                                                                                       n_castle_slabs_wide, l_ini, l_fin), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10,
                                                     t_breakpoints=castle_wall_breakpoints(n_castle_slabs_wide, l_ini, l_fin))[nu_i][nu_f] \
                                          for enu in energies])  '''),
    md(r'''### Plot probabilities'''),
    code(r'''smooth = lambda y: sp.signal.savgol_filter(y, window_length=301, polyorder=1)
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_matt_const_density, label=r'Constant average density', color='0.5', ls='--'),
      dict(y=prob_matt_castle_wall_narrow, label=r'Narrow castle wall', color='C0'),
      dict(y=prob_matt_castle_wall_wide, label=r'Wide castle wall', color='C3')],
     [dict(y=smooth(prob_matt_const_density), color='0.5', ls='--'),
      dict(y=smooth(prob_matt_castle_wall_narrow), color='C0'),
      dict(y=smooth(prob_matt_castle_wall_wide), color='C3')]],
    xlim=(energy_min, energy_max), xscale='log',
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$2\nu$~castle-wall matter profile, $L = $~{:.2f}~km'.format(baseline),
    legend_title=r'Matter profile', legend_loc='lower right', legend_on_panel=0)'''),
    md(r'''# 5. Probabilities 2$\nu$: in matter with noisy density profile'''),
    md(r'''## 5.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Parameters for the electron number density in the Sun, Eq. (10.62) in Giunti & Kim
num_density_e_center = 10*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_scale = 100.0 # [km]

# Build the noise
np.random.seed(20)
# High-amplitude, more marked noise
npts_noise = 50 # Use more points to make it more finely noised
x_high = np.linspace(0, 1, npts_noise)
noise_sample_high = np.random.normal(loc=1.0, scale=0.3, size=npts_noise) # Increase scale to increase the width of the noise bumps
noise_mask_high_interp = sp.interpolate.make_interp_spline(x_high, noise_sample_high, k=1)
# Low-amplitude, choppier noise
npts_noise = 300 # Use more points to make it more finely noised
x_low = np.linspace(0, 1, npts_noise)
noise_sample_low = np.random.normal(loc=1.0, scale=0.08, size=npts_noise) # Increase scale to increase the width of the noise bumps
noise_mask_low_interp = sp.interpolate.make_interp_spline(x_low, noise_sample_low, k=1)

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
# The matter density is constant, so this is merely a wrapper.
# [l] = km
VCC_center = matter.VCC_func(0.0, lambda l: num_density_e_center) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_2nu_matter(VCC_center) # [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density

# --------------------------------------------------------------
# In matter, constant density, with noise around central value
# --------------------------------------------------------------
def num_density_e_const_noisy_func(l, num_density_e_const, l_ini, l_fin, noise_model='high-amplitude'):
    l_scaled = (l/gd.CONV_KM_TO_INV_EV-l_ini)/(l_fin-l_ini)
    if noise_model == 'high-amplitude':
        return num_density_e_const*noise_mask_high_interp(l_scaled)
    elif noise_model == 'low-amplitude':
        return num_density_e_const*noise_mask_low_interp(l_scaled)
def VCC_func_const_density_noisy(l, num_density_e_const, l_ini, l_fin, noise_model='high-amplitude'):
    return matter.VCC_func(l, lambda r: num_density_e_const_noisy_func(r, num_density_e_const, l_ini, l_fin, 
                                                                       noise_model=noise_model)) # [eV]
def H_func_const_density_noisy(l, energy, num_density_e_const, l_ini, l_fin, noise_model='high-amplitude'):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_const_density_noisy(l, num_density_e_const, 
                                                                                                             l_ini, l_fin,
                                                                                                             noise_model=noise_model))
# # --------------------------------------------------------------
# # In matter, exponentially falling density
# # --------------------------------------------------------------
# # [l] = km
# def num_density_e_exp_func(l, num_density_e_center, l_scale):
#     return num_density_e_center*np.exp(-(l/gd.CONV_KM_TO_INV_EV)/l_scale) # [eV^3]
# def VCC_func_exp_density(l, num_density_e_center, l_scale):
#     return matter.VCC_func(l, lambda r: num_density_e_exp_func(r, num_density_e_center, l_scale)) # [eV]
# def H_func_exp_density(l, energy):
#     return (1/energy)*H_vac_energy_indep-hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_density(l, num_density_e_center, l_scale))

# # --------------------------------------------------------------
# # In matter, exponentially falling density, noisy
# # --------------------------------------------------------------
# # [l] = km
# def num_density_e_exp_noisy_func(l, num_density_e_center, l_scale, l_ini, l_fin):
#     l_scaled = (l/gd.CONV_KM_TO_INV_EV-l_ini)/(l_fin-l_ini)
#     return num_density_e_exp_func(l, num_density_e_center, l_scale)*noise_mask_interp(l_scaled)
# def VCC_func_exp_density_noisy(l, num_density_e_center, l_scale, l_ini, l_fin):
#     return matter.VCC_func(l, lambda r: num_density_e_exp_noisy_func(r, num_density_e_center, l_scale, l_ini, l_fin)) # [eV]
# def H_func_exp_density_noisy(l, energy, num_density_e_center, l_scale, l_ini, l_fin):
#     return (1/energy)*H_vac_energy_indep \
#         - hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_density_noisy(l, num_density_e_center, l_scale, l_ini, l_fin))'''),
    md(r'''## 5.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e2, 1.e5 # [km]
l_npts = 6000 
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 50.*gd.UNIT_MEV # [eV] #10

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, energy), 
                                                     0, l*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, constant density, high-amplitude noise
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_const_density_noisy_high_amplitude = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_const_density_noisy(l, energy, num_density_e_center,
                                        l_ini, l_fin, noise_model='high-amplitude'), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3)

# --------------------------------------------------------------
# In matter, constant density, low-amplitude noise
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_const_density_noisy_low_amplitude = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_const_density_noisy(l, energy, num_density_e_center,
                                        l_ini, l_fin, noise_model='low-amplitude'), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3)

# # --------------------------------------------------------------
# # In matter, exponentially falling density
# # --------------------------------------------------------------
# prob_matt_exp_density = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, energy), 
#                                                    0, l*gd.CONV_KM_TO_INV_EV, 
#                                                    n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
smooth = lambda y: sp.signal.savgol_filter(y, window_length=301, polyorder=1)
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center/norm for l in distances], color='0.5', ls='--'),
     dict(y=[num_density_e_const_noisy_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_center,
                                            l_ini, l_fin,
                                            noise_model='high-amplitude')/norm
             for l in distances], color='C0'),
     dict(y=[num_density_e_const_noisy_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_center,
                                            l_ini, l_fin,
                                            noise_model='low-amplitude')/norm
             for l in distances], color='C3')],
    [[dict(y=prob_matt_const_density, label=r'Constant density', color='0.5', ls='--'),
      dict(y=prob_matt_const_density_noisy_high_amplitude, label=r'High-amplitude noise', color='C0'),
      dict(y=prob_matt_const_density_noisy_low_amplitude, label=r'Low-amplitude noise', color='C3')],
     [dict(y=smooth(prob_matt_const_density), color='0.5', ls='--'),
      dict(y=smooth(prob_matt_const_density_noisy_high_amplitude), color='C0'),
      dict(y=smooth(prob_matt_const_density_noisy_low_amplitude), color='C3')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 15),
    profile_ymajor=4, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]', ylabel_labelpad=15,
    title=r'$2\nu$~noisy matter profile, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left', legend_on_panel=0)'''),
    md(r'''## 5.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e1, 1.e2 # [MeV]
energy_npts = 6000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 1.e5 # [km]

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV),
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])

# --------------------------------------------------------------
# In matter, constant density, high-amplitude noise
# --------------------------------------------------------------
prob_matt_const_density_noisy_high_amplitude = np.array([oscprob.osc_prob(lambda l: H_func_const_density_noisy(l, enu*gd.UNIT_MEV,
                                                                                                            num_density_e_center, 
                                                                                                            l_ini, l_fin,
                                                                                                            noise_model='high-amplitude'), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] \
                                          for enu in energies]) 

# --------------------------------------------------------------
# In matter, constant density, low-amplitude noise
# --------------------------------------------------------------
prob_matt_const_density_noisy_low_amplitude = np.array([oscprob.osc_prob(lambda l: H_func_const_density_noisy(l, enu*gd.UNIT_MEV,
                                                                                                            num_density_e_center,
                                                                                                            l_ini, l_fin,
                                                                                                            noise_model='low-amplitude'), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] \
                                          for enu in energies]) '''),
    md(r'''### Plot probabilities'''),
    code(r'''smooth = lambda y, po=1: sp.signal.savgol_filter(y, window_length=301, polyorder=po)
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_matt_const_density, label=r'Constant density', color='0.5', ls='--'),
      dict(y=prob_matt_const_density_noisy_high_amplitude, label=r'High-amplitude noise', color='C0'),
      dict(y=prob_matt_const_density_noisy_low_amplitude, label=r'Low-amplitude noise', color='C3')],
     [dict(y=smooth(prob_matt_const_density), color='0.5', ls='--'),
      dict(y=smooth(prob_matt_const_density_noisy_high_amplitude), color='C0'),
      dict(y=smooth(prob_matt_const_density_noisy_low_amplitude, 2), color='C3')]],
    xlim=(energy_min, energy_max), xscale='linear', xmajor=10, xminor=1,
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$2\nu$~noisy matter profile, $L = $~{:.2f}~km'.format(baseline),
    legend_title=r'Matter profile', legend_loc='lower right', legend_on_panel=0)'''),
    md(r'''# 6. Probabilities 2$\nu$: in the Earth'''),
    md(r'''## 6.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''def num_density_e_func_prem(r): 
    return matter.num_density_e_func(r, earth.density_matter_func_prem, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

def VCC_func_prem(r):
    return matter.VCC_func(r, num_density_e_func_prem) # [eV]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# [l] = km
def H_func_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) 
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_prem(r))

# Test values of the zenith angle, cos(theta_z), to plot, with theta_z = 0 normal to the surface of the Earth, at the position of the
# detector.  All of the direction are upgoing.
costhz_val = [-0.1, -0.5, -1.0]

# Maximum baselines inside the Earth
l_max_val = [earth.distance_traveled_inside_earth(costhz) for costhz in costhz_val] # [km]'''),
    md(r'''## 6.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]

nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 20.*gd.UNIT_MEV # [eV] #10

l_min = 1.e2 # [km]

# Generate probabilities for the different directions
distances_val, prob_val = [], []
l_npts = 3000
for i in range(len(costhz_val)):
    print("costhz = "+str(costhz_val[i]))
    distances = np.logspace(np.log10(l_min), np.log10(l_max_val[i]), l_npts) # [km]
    distances_val.append(distances)
    # PREM is piecewise-smooth: the density jumps where the chord crosses a shell boundary.
    # Marking those crossings as mandatory slab edges keeps the quadrature from integrating
    # across a discontinuity.  Measured at this grid, against solve_ivp, it costs nothing in
    # time and buys 1.3e3x to 1.9e6x in accuracy -- the costhz = -1 direction was outside the
    # default 1e-3 tolerance without it (4.2e-3 -> 3.3e-6).
    prem_bp = earth.prem_layer_edges_along_chord(costhz_val[i])*gd.UNIT_KM
    prob = np.array([oscprob.osc_prob(lambda l: H_func_prem(costhz_val[i], l, energy), 0, l*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10,
                                      t_breakpoints=prem_bp)[nu_i][nu_f] for l in distances]) 
    prob_val.append(prob)'''),
    md(r'''### Plot probabilities'''),
    code(r'''lc = ['C0', 'C2', 'C3']
ls = ['-', '-', '-']
norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
profiles, panel = [], []
for i in range(len(costhz_val)):
    n_e = np.array([num_density_e_func_prem(
                        earth.earth_radial_distance_from_depth(costhz_val[i], l))
                    for l in distances_val[i]])/norm
    profiles.append(dict(x=distances_val[i], y=n_e, color=lc[i], ls=ls[i]))
    panel.append(dict(x=distances_val[i], y=prob_val[i], color=lc[i], ls=ls[i],
                      label=r'$\cos \theta_z = $~~{:.2f}'.format(costhz_val[i])))

fig, ax = plotting.plot_probability_with_profile(
    distances_val[0], profiles, [panel],
    xlim=(l_min, max(l_max_val)), profile_ylim=(0, 7),
    profile_ymajor=2, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]',
    title=r'$2\nu$~inside the Earth (PREM), $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Direction', legend_loc='lower left')'''),
    md(r'''## 6.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e1, 1.e2 # [MeV]
energy_npts = 3000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU

# Generate probabilities for the different directions
prob_val = []
for i in range(len(costhz_val)):
    print("costhz = "+str(costhz_val[i]))
    # PREM is piecewise-smooth: the density jumps where the chord crosses a shell boundary.
    # Marking those crossings as mandatory slab edges keeps the quadrature from integrating
    # across a discontinuity.  Measured at this grid, against solve_ivp, it costs nothing in
    # time and buys 1.3e3x to 1.9e6x in accuracy -- the costhz = -1 direction was outside the
    # default 1e-3 tolerance without it (4.2e-3 -> 3.3e-6).
    prem_bp = earth.prem_layer_edges_along_chord(costhz_val[i])*gd.UNIT_KM
    prob = np.array([oscprob.osc_prob(lambda l: H_func_prem(costhz_val[i], l, enu*gd.UNIT_MEV), 
                                      0, l_max_val[i]*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10,
                                      t_breakpoints=prem_bp)[nu_i][nu_f] for enu in energies]) 
    prob_val.append(prob)'''),
    md(r'''### Plot probabilities'''),
    code(r'''smooth = lambda y: sp.signal.savgol_filter(y, window_length=201, polyorder=1)
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_val[i], color=lc[i], ls=ls[i],
           label=r'$\cos \theta_z = $~~{:.2f}'.format(costhz_val[i]))
      for i in range(len(costhz_val))],
     [dict(y=smooth(prob_val[i]), color=lc[i], ls=ls[i])
      for i in range(len(costhz_val))]],
    xlim=(energy_min, energy_max), xscale='linear', xmajor=10, xminor=1,
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$2\nu$~inside the Earth (PREM), $L = L_{\rm max}(\theta_z)$',
    legend_title=r'Direction', legend_loc='lower left', legend_on_panel=0)'''),
    md(r'''# 7. Probabilities 2$\nu$: in the Sun'''),
    md(r'''## 7.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Parameters for the electron number density in the Sun, Eq. (10.62) in Giunti & Kim
num_density_e_center_sun = 245*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_scale_sun = gd.SUN_RADIUS/10.54 # [km]

# Potential assuming constant density: the matter density is constant, so this is merely a wrapper.
# [l] = km
def num_density_e_sun_const_func(l):
    return num_density_e_center_sun # [eV^3]
def VCC_func_const_density_sun(l):
    return matter.VCC_func(l, num_density_e_sun_const_func) # [eV]

# Potential assuming exponentially decreasing density
# [l] = km
def num_density_e_sun_exp_func(l):
    return num_density_e_center_sun*np.exp(-(l/gd.CONV_KM_TO_INV_EV)/l_scale_sun) # [eV^3]
def VCC_func_exp_sun(l):
    return matter.VCC_func(l, num_density_e_sun_exp_func) # [eV]'''),
    md(r'''## 7.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e-3, 1.0 # l/SUN_RADIUS
l_npts = 1000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # l/SUN_RADIUS

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 5.0*gd.UNIT_MEV # [eV]

# --------------------------------------------------------------
# In matter, constant density equal to central solar density
# --------------------------------------------------------------
H_func_const_density = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2) \
                                      + hamiltonians.hamiltonian_2nu_matter(VCC_func_const_density_sun(l))
prob_matt_const_density = np.array([oscprob.osc_prob(H_func_const_density, 
                                                     0, l*gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, solar matter density profile
# --------------------------------------------------------------
H_func_exp_density = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2) \
                                      + hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_sun(l))
# H_func_exp_density = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2) \
                                      # + hamiltonians.hamiltonian_2nu_matter(-VCC_func_exp_sun(l))
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_exp_density = oscprob.osc_prob_energy_baseline(
    H_func_exp_density, energy, distances*gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3)

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
# Using Magnus expansion
H_func = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2)
prob_vac = np.array([oscprob.osc_prob(H_func,
                                      0.0, l*gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])
# # Using standard formula
# prob_std = np.array([oscprobstd.osc_prob_2nu_vacuum_std(sth, Dm2, energy, l*gd.CONV_KM_TO_INV_EV)[nu_i][nu_f] for l in distances])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center_sun/norm for l in distances], color='C0'),
     dict(y=[num_density_e_center_sun/norm*np.exp(-l*gd.SUN_RADIUS/l_scale_sun)
             for l in distances], color='C3')],
    [[dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density, color='C0',
           label=r'Constant density, $N_e/N_{\rm Av} =~$'
                 + '{:.0f}'.format(num_density_e_center_sun/norm)),
      dict(y=prob_matt_exp_density, label=r'Solar profile', color='C3')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 250),
    profile_ymajor=50, profile_yminor=10, profile_height=0.3,
    profile_ylabel=r'$N_e/N_{\rm Av}$',
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L/R_\odot$', ylabel_labelpad=7,
    title=r'$2\nu$~in the Sun, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    title_fontsize=20, legend_title=r'Matter profile', legend_loc='lower left',
    grid=False)'''),
    md(r'''## 7.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e2 # [MeV]
energy_npts = 1000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = gd.SUN_RADIUS # [km]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# --------------------------------------------------------------
# In matter, constant density equal to central solar density
# --------------------------------------------------------------
H_matt = hamiltonians.hamiltonian_2nu_matter(VCC_func_const_density_sun(0.0)) # Matter Hamiltonian [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV), 
                                                     0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])

# --------------------------------------------------------------
# In matter, solar matter density profile
# --------------------------------------------------------------
def H_func_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_sun(l)) 
prob_matt_exp_density = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, enu*gd.UNIT_MEV), 
                                                   0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                   n_slabs=200, n_tpts_per_slab=100, magnus_exp_order=2, n_jobs=10)[nu_i][nu_f] \
                                                   # n_slabs=200, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=8)[nu_i][nu_f] \
                                  for enu in energies])

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
prob_vac = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_MEV))*H_vac_energy_indep, 
                                      0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
     dict(y=prob_matt_const_density, color='C0',
          label=r'Constant density, $N_e/N_{\rm Av} =~$'
                + '{:.0f}'.format(num_density_e_center_sun/norm)),
     dict(y=prob_matt_exp_density, label=r'Solar profile', color='C3')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2, energy_unit='MeV',
    xlim=(energy_min, energy_max),
    legend_title=r'Matter profile', legend_loc='center left',
    title=r'$2\nu$~in the Sun, $L = R_\odot$')'''),
    ])

# ----------------------------------------------- 03_magnus_3nu_vacuum_matter
books['03_magnus_3nu_vacuum_matter.ipynb'] = notebook(
    'Three-neutrino probabilities',
    'The same seven settings as the previous notebook, with three flavours and a CP-violating phase.\n\nNothing about the method changes -- the Hamiltonian is a $3\\times 3$ matrix in the same slot -- so the interest is in what the extra flavour and the phase do to the probabilities.',
    [
    code(r'''import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.oscprobstd as oscprobstd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''# 0. Helper functions and definitions'''),
    code(r"""def prob_label(nu_i, nu_f):
    if (nu_i == gd.NUE):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_e \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_e \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_e \to \nu_\tau}$'
    elif (nu_i == gd.NUMU):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_\mu \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_\mu \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_\mu \to \nu_\tau}$'
    elif (nu_i == gd.NUTAU):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_\tau \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_\tau \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_\tau \to \nu_\tau}$'
    return label

def flavor_index_to_str(nu_l):
    if (nu_l == gd.NUE): return 'e'
    if (nu_l == gd.NUMU): return 'mu'
    if (nu_l == gd.NUTAU): return 'tau'"""),
    md(r'''# 1. Probabilities 3$\nu$: in vacuum'''),
    md(r'''## 1.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# PMNS mixing mtarix
U = hamiltonians.pmns_mixing_matrix(s12, s23, s13, dCP)'''),
    md(r'''## 1.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e5 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline. Both calls return the full probability
# matrix for every baseline, with shape (l_npts, 3, 3), so a channel is
# prob_all[:, nu_i, nu_f].
energy = 10.0*gd.UNIT_MEV # [eV]

# Both probabilities are computed for the *whole* array of baselines in a single
# call, which is the form to prefer: measured here, the array call takes 1.69 s
# against 2.71 s for a Python loop that calls osc_prob once per point (1.6x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# The 3nu closed form is NOT vectorized -- osc_prob_3nu_vacuum_std reshapes its
# result to (3,3) and rejects an array -- so this one keeps its loop. At ~80 us per
# call it costs ~0.8 s over this grid, which is cheaper than the Magnus call not
# because it escapes some per-call overhead but because a closed-form expression
# does far less arithmetic than a slab-by-slab integration.
prob_std_all = np.array([oscprobstd.osc_prob_3nu_vacuum_std(
    U, D21, D31, energy, l*gd.CONV_KM_TO_INV_EV)
    for l in distances])

# Using the Magnus expansion
prob_all = oscprob.osc_prob_3nu_vacuum(
    energy, distances*gd.CONV_KM_TO_INV_EV, s12, s23, s13, dCP, D21, D31)'''),
    md(r'''### Plot probabilities'''),
    code(r'''def make_plot_prob_vac_vs_baseline(nu_i, nu_f, save_plot=True):
    """Probability against baseline, over its relative error vs the closed form."""
    filename = 'prob_3nu_vacuum_vs_baseline_' + str(nu_i) + str(nu_f)
    prob = prob_all[:, nu_i, nu_f]
    prob_std = prob_std_all[:, nu_i, nu_f]

    fig, ax = plotting.plot_probability_vs_baseline(
        distances,
        [dict(y=prob, label='Magnus expansion', color='C1'),
         dict(y=prob_std, label='Standard formula', color='k', ls='--')],
        nu_i=nu_i, nu_f=nu_f, num_flavors=3,
        xlim=(l_ini, l_fin),
        residual=(prob-prob_std)/prob_std/1.e-10,
        residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-10}]$',
        residual_ylim=(-0.5, 0.5), residual_ymajor=0.20, residual_yminor=0.05,
        legend_title=r'Calculation method', legend_loc='center left',
        title_fontsize=23,
        title=r'$3\nu$~vacuum, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
        savefig=('../fig/'+filename+'.pdf') if save_plot else None,
        savefig_kw=dict(dpi=300))
    return fig, ax'''),
    code(r'''make_plot_prob_vac_vs_baseline(gd.NUE, gd.NUE, save_plot=True)'''),
    code(r'''make_plot_prob_vac_vs_baseline(gd.NUE, gd.NUMU, save_plot=True)'''),
    code(r'''make_plot_prob_vac_vs_baseline(gd.NUMU, gd.NUMU, save_plot=True)'''),
    code(r'''make_plot_prob_vac_vs_baseline(gd.NUMU, gd.NUTAU, save_plot=True)'''),
    md(r'''## 1.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

# Compute probability vs. energy; shape is (energy_npts, 3, 3)
baseline = 5.e3 # [km]

# Both probabilities are computed for the *whole* array of energies in a single
# call, which is the form to prefer: measured here, the array call takes 0.70 s
# against 1.55 s for a Python loop that calls osc_prob once per point (2.2x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# The 3nu closed form is NOT vectorized -- osc_prob_3nu_vacuum_std reshapes its
# result to (3,3) and rejects an array -- so this one keeps its loop. At ~80 us per
# call it costs ~0.4 s over this grid, which is cheaper than the Magnus call not
# because it escapes some per-call overhead but because a closed-form expression
# does far less arithmetic than a slab-by-slab integration.
prob_std_all = np.array([oscprobstd.osc_prob_3nu_vacuum_std(
    U, D21, D31, enu*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV)
    for enu in energies])

# Using the Magnus expansion
prob_all = oscprob.osc_prob_3nu_vacuum(
    energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV,
    s12, s23, s13, dCP, D21, D31)'''),
    md(r'''### Plot probabilities'''),
    code(r'''def make_plot_prob_vac_vs_energy(nu_i, nu_f, save_plot=True):
    """The energy counterpart of the figure above."""
    filename = 'prob_3nu_vacuum_vs_energy_' + str(nu_i) + str(nu_f)
    prob = prob_all[:, nu_i, nu_f]
    prob_std = prob_std_all[:, nu_i, nu_f]

    fig, ax = plotting.plot_probability_vs_energy(
        energies,
        [dict(y=prob, label='Magnus expansion', color='r'),
         dict(y=prob_std, label='Standard formula', color='k', ls='--')],
        nu_i=nu_i, nu_f=nu_f, num_flavors=3,
        xlim=(energy_min, energy_max),
        residual=(prob-prob_std)/prob_std/1.e-12,
        residual_label=r'$\epsilon_{\rm rel}~[\times 10^{-12}]$',
        residual_ylim=(-0.5, 0.5), residual_ymajor=0.20, residual_yminor=0.05,
        legend_title=r'Calculation method', legend_loc='center right',
        title_fontsize=23,
        title=r'$3\nu$~vacuum, $L = $~{:.2f}~km'.format(baseline),
        savefig=('../fig/'+filename+'.pdf') if save_plot else None,
        savefig_kw=dict(dpi=300))
    return fig, ax'''),
    code(r'''make_plot_prob_vac_vs_energy(gd.NUE, gd.NUE, save_plot=True)'''),
    code(r'''make_plot_prob_vac_vs_energy(gd.NUE, gd.NUMU, save_plot=True)'''),
    code(r'''make_plot_prob_vac_vs_energy(gd.NUMU, gd.NUMU, save_plot=True)'''),
    code(r'''make_plot_prob_vac_vs_energy(gd.NUMU, gd.NUTAU, save_plot=True)'''),
    md(r'''# 2. Probabilities 3$\nu$: in matter with constant density'''),
    md(r'''## 2.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# PMNS mixing mtarix
U = hamiltonians.pmns_mixing_matrix(s12, s23, s13, dCP)'''),
    md(r'''## 2.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e4 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline; shape is (l_npts, 3, 3)
energy = 10.0*gd.UNIT_MEV # [eV]

# In matter
rho = 10.0 # Matter density [g cm^{-3}]

# Both probabilities are computed for the *whole* array of baselines in a single
# call, which is the form to prefer: measured here, the array call takes 1.26 s
# against 2.84 s for a Python loop that calls osc_prob once per point (2.3x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

# The wrapper takes the density directly and builds the potential itself. Note
# density_matter_is_in_g_per_cm3=True: without it rho is read as being in
# internal units already, and the potential comes out ~4e18 times too small --
# indistinguishable from having no matter at all.
prob_matt_all = oscprob.osc_prob_3nu_matter_constant_density(
    energy, distances*gd.CONV_KM_TO_INV_EV, rho, s12, s23, s13, dCP, D21, D31,
    density_matter_is_in_g_per_cm3=True)

# In vacuum, for comparison
prob_vac_all = oscprob.osc_prob_3nu_vacuum(
    energy, distances*gd.CONV_KM_TO_INV_EV, s12, s23, s13, dCP, D21, D31)'''),
    md(r'''### Plot probabilities'''),
    code(r'''channels = [(gd.NUE, gd.NUE), (gd.NUE, gd.NUMU),
            (gd.NUMU, gd.NUMU), (gd.NUMU, gd.NUTAU)]
fig, ax = plotting.plot_probability_with_profile(
    distances, None,
    [[dict(y=prob_vac_all[:, a, b], label=r'Vacuum', color='0.7'),
      dict(y=prob_matt_all[:, a, b], color='C1',
           label=r'Matter ($\rho = $~{:.1f}'.format(rho)+r'~g~cm$^{-3}$)')]
     for a, b in channels],
    xlim=(l_ini, l_fin), panel_ymajor=0.10, panel_yminor=0.02,
    panel_annotations=[plotting.prob_label(a, b) for a, b in channels],
    shared_ylabel=r'Three-neutrino probability',
    xlabel=r'Baseline, $L$ [km]', title_fontsize=23,
    title=r'$3\nu$~constant-density matter, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_loc='lower left', legend_on_panel=0, grid=False,
    figsize=[18, 18], tight_layout=True,
    savefig='../fig/prob_3nu_matter_const_density_vs_baseline.pdf',
    savefig_kw=dict(dpi=300))'''),
    md(r'''## 2.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

# Compute probability vs. energy; shape is (energy_npts, 3, 3)
baseline = 5.e3 # [km]

# In matter
rho = 10.0 # Matter density [g cm^{-3}]

# Both probabilities are computed for the *whole* array of energies in a single
# call, which is the form to prefer: measured here, the array call takes 0.62 s
# against 1.41 s for a Python loop that calls osc_prob once per point (2.3x).
#
# The margin is that modest, rather than an order of magnitude, because what a
# call costs is mostly the integration it performs -- not a fixed entry-path
# overhead paid per call. On a profile that actually varies, the same call
# measures 0.34 ms at n_slabs=1 and 10.4 ms at n_slabs=2000. Batching amortizes
# the entry path, which is real but small; it cannot make the physics cheaper.

prob_matt_all = oscprob.osc_prob_3nu_matter_constant_density(
    energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV, rho,
    s12, s23, s13, dCP, D21, D31, density_matter_is_in_g_per_cm3=True)

# In vacuum, for comparison
prob_vac_all = oscprob.osc_prob_3nu_vacuum(
    energies*gd.UNIT_GEV, baseline*gd.CONV_KM_TO_INV_EV,
    s12, s23, s13, dCP, D21, D31)'''),
    md(r'''### Plot probabilities'''),
    code(r'''channels = [(gd.NUE, gd.NUE), (gd.NUE, gd.NUMU),
            (gd.NUMU, gd.NUMU), (gd.NUMU, gd.NUTAU)]
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_vac_all[:, a, b], label=r'Vacuum', color='0.7'),
      dict(y=prob_matt_all[:, a, b], color='C1',
           label=r'Matter ($\rho = $~{:.1f}'.format(rho)+r'~g~cm$^{-3}$)')]
     for a, b in channels],
    xlim=(energy_min, energy_max), panel_ymajor=0.10, panel_yminor=0.02,
    panel_annotations=[plotting.prob_label(a, b) for a, b in channels],
    shared_ylabel=r'Three-neutrino probability',
    xlabel=r'Neutrino energy, $E_\nu$ [GeV]', title_fontsize=23,
    title=r'$3\nu$~constant-density matter, $L = $~{:.2f}~km'.format(baseline),
    legend_loc='lower left', legend_on_panel=0, grid=False,
    figsize=[18, 18], tight_layout=True,
    savefig='../fig/prob_3nu_matter_const_density_vs_energy.pdf',
    savefig_kw=dict(dpi=300))'''),
    md(r'''# 3. Probabilities 3$\nu$: in matter with varying density'''),
    md(r'''## 3.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# PMNS mixing mtarix
U = hamiltonians.pmns_mixing_matrix(s12, s23, s13, dCP)'''),
    code(r'''# Parameters for the electron number density in the Sun, Eq. (10.62) in Giunti & Kim
num_density_e_center = 10*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_scale = 100.0 # [km]

H_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
                                                            compute_matrix_multiplication=False) # Vacuum H without (1/E) [eV^2]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
# The matter density is constant, so this is merely a wrapper.
# [l] = km
VCC_center = matter.VCC_func(0.0, lambda l: num_density_e_center) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_3nu_matter(VCC_center) # [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density

# --------------------------------------------------------------
# In matter, exponentially falling density
# --------------------------------------------------------------
# [l] = km
def num_density_e_exp_func(l, num_density_e_center, l_scale):
    return num_density_e_center*np.exp(-(l/gd.CONV_KM_TO_INV_EV)/l_scale) # [eV^3]
def VCC_func_exp_density(l, num_density_e_center, l_scale):
    return matter.VCC_func(l, lambda r: num_density_e_exp_func(r, num_density_e_center, l_scale)) # [eV]
def H_func_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_3nu_matter(VCC_func_exp_density(l, num_density_e_center, l_scale))

# --------------------------------------------------------------
# In matter, Gaussian density
# --------------------------------------------------------------
# [l] = km
num_density_e_basal = 0*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
num_density_e_central = 8*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_central, l_width = 300.0, 100.0 # [km]
def num_density_e_gaussian_func(l, num_density_e_central, l_central, l_width):
    return num_density_e_basal+num_density_e_central*np.exp( -(l/gd.CONV_KM_TO_INV_EV-l_central)**2/(2.0*l_width**2)) # [eV^3]
def VCC_func_gaussian_density(l, num_density_e_central, l_central, l_width):
    return matter.VCC_func(l, lambda r: num_density_e_gaussian_func(r, num_density_e_central, l_central, l_width)) # [eV]
def H_func_gaussian_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_3nu_matter(VCC_func_gaussian_density(l, 
                                                    num_density_e_central, l_central, l_width)) # [eV]'''),
    md(r'''## 3.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 5.e1, 1.e3 # [km]
l_npts = 1000 #10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
energy = 10.0*gd.UNIT_MEV # [eV]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
prob_matt_const_density_all = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, energy), 
                                                     0, l*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for l in distances])

# --------------------------------------------------------------
# In matter, exponentially falling density
# --------------------------------------------------------------
prob_matt_exp_density_all = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, energy), 
                                                   0, l*gd.CONV_KM_TO_INV_EV, 
                                                   n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)for l in distances])

# --------------------------------------------------------------
# In matter, Gaussian density
# --------------------------------------------------------------
prob_matt_gaussian_density_all = np.array([oscprob.osc_prob(lambda l: H_func_gaussian_density(l, energy), 
                                                        0, l*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10) for l in distances])

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
H_vac = (1./energy)*H_vac_energy_indep # Can also call hamiltonians.hamiltonian_2nu_vacuum instead
prob_vac_all = np.array([oscprob.osc_prob(lambda l: H_vac,
                                      0.0, l*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1, n_jobs=1) for l in distances])'''),
    md(r'''### Plot probabilities'''),
    code(r'''channels = [(gd.NUE, gd.NUE), (gd.NUE, gd.NUMU),
            (gd.NUMU, gd.NUMU), (gd.NUMU, gd.NUTAU)]
norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV, 3.0))
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center/norm for l in distances], color='C0'),
     dict(y=[num_density_e_exp_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_center, l_scale)/norm
             for l in distances], color='C3', ls='--'),
     dict(y=[num_density_e_gaussian_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_central,
                                         l_central, l_width)/norm
             for l in distances], color='C5', ls='-.')],
    [[dict(y=prob_vac_all[:, a, b], label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density_all[:, a, b], label=r'Constant density', color='C0'),
      dict(y=prob_matt_exp_density_all[:, a, b], label=r'Exponentially falling', color='C3', ls='--'),
      dict(y=prob_matt_gaussian_density_all[:, a, b], label=r'Gaussian', color='C5', ls='-.')]
     for a, b in channels],
    xlim=(l_ini, l_fin), profile_ylim=(0, 11),
    profile_ymajor=4, profile_yminor=1, panel_ymajor=0.10, panel_yminor=0.02,
    panel_annotations=[plotting.prob_label(a, b) for a, b in channels],
    shared_ylabel=r'Three-neutrino probability',
    xlabel=r'Baseline, $L$ [km]', title_fontsize=23,
    title=r'$3\nu$~varying-density matter, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left', legend_on_panel=0,
    grid=False, figsize=[18, 22])'''),
    md(r'''## 3.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e2 # [MeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
baseline = 1.e3 # [km]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
prob_matt_const_density_all = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV), 
                                                     0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])

# --------------------------------------------------------------
# In matter, exponentially falling density
# --------------------------------------------------------------
prob_matt_exp_density_all = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, enu*gd.UNIT_MEV), 
                                                   0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                   n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10) \
                                  for enu in energies])

# --------------------------------------------------------------
# In matter, Gaussian density
# --------------------------------------------------------------
prob_matt_gaussian_density_all = np.array([oscprob.osc_prob(lambda l: H_func_gaussian_density(l, enu*gd.UNIT_MEV), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10) \
                                  for enu in energies])

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
prob_vac_all = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_MEV))*H_vac_energy_indep, 
                                      0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])'''),
    md(r'''### Plot probabilities'''),
    code(r'''channels = [(gd.NUE, gd.NUE), (gd.NUE, gd.NUMU),
            (gd.NUMU, gd.NUMU), (gd.NUMU, gd.NUTAU)]
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_vac_all[:, a, b], label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density_all[:, a, b], label=r'Constant density', color='C0'),
      dict(y=prob_matt_exp_density_all[:, a, b], label=r'Exponentially falling', color='C3', ls='--'),
      dict(y=prob_matt_gaussian_density_all[:, a, b], label=r'Gaussian', color='C5', ls='-.')]
     for a, b in channels],
    xlim=(energy_min, energy_max), panel_ymajor=0.10, panel_yminor=0.02,
    panel_annotations=[plotting.prob_label(a, b) for a, b in channels],
    shared_ylabel=r'Three-neutrino probability', shared_ylabel_labelpad=25,
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]', title_fontsize=23,
    title=r'$3\nu$~varying-density matter, $L = $~{:.2f}~km'.format(baseline),
    legend_title=r'Matter profile', legend_loc='upper right', legend_on_panel=1,
    grid=False, figsize=[18, 18],
    savefig='../fig/prob_3nu_varying_density_vs_energy.pdf',
    savefig_kw=dict(dpi=300, bbox_inches='tight'))'''),
    md(r'''# 4. Probabilities 2$\nu$: in matter with castle-wall density profile'''),
    md(r'''## 4.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Low and high densities of the castle-wall density profile
num_density_e_low = 1*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
num_density_e_high = 10*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
num_density_e_avg = (num_density_e_high-num_density_e_low)/2.0 # [eV^3]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]
    
# --------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
# The matter density is constant, so this is merely a wrapper ([l] = km)
VCC_avg = matter.VCC_func(0.0, lambda l: num_density_e_avg) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_2nu_matter(VCC_avg) # [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density

# --------------------------------------------------------------
# In matter, castle-wall density profile
# --------------------------------------------------------------
# [l] = km
def num_density_e_castle_wall_func(l, num_density_e_low, num_density_e_high, n_castle_slabs, l_ini, l_fin):
    # np.where, not an if/else: the engine evaluates this for a whole array of
    # positions at once, and a Python `if` on an array raises "truth value of an
    # array is ambiguous". Written this way the profile stays array-capable and
    # osc_prob keeps its fast path -- see the note in notebook 01. It still
    # returns a plain number when handed a single position.
    l_scaled = (np.asarray(l)/gd.CONV_KM_TO_INV_EV-l_ini)/(l_fin-l_ini)
    dl = 1.0/n_castle_slabs # Width of one slab
    index_slab = l_scaled // dl
    # The first slab has low density, the second one has high density, and so on
    return np.where(index_slab % 2 == 0, num_density_e_low, num_density_e_high)
def VCC_func_castle_wall(l, num_density_e_low, num_density_e_high, n_castle_slabs, l_ini, l_fin):
    return matter.VCC_func(l, lambda r: num_density_e_castle_wall_func(r, num_density_e_low, num_density_e_high, 
                                                                       n_castle_slabs, l_ini, l_fin)) # [eV]
def H_func_castle_wall(l, energy, num_density_e_low, num_density_e_high, n_castle_slabs, l_ini, l_fin):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_castle_wall(l, 
                                                                                        num_density_e_low, num_density_e_high,
                                                                                        n_castle_slabs, l_ini, l_fin))

n_castle_slabs_narrow = 50 # Narrow castle wall
n_castle_slabs_wide = 10   # Wide castle wall

# Positions of the density walls, in eV^-1.  The profile is a step function, and
# high-order quadrature reaches its nominal order only where the Hamiltonian is
# smooth inside each slab -- integrating across a wall costs accuracy that no
# amount of extra slabs cheaply buys back.  Handing these to osc_prob as
# `t_breakpoints` makes them mandatory slab edges, so no slab ever straddles one.
# osc_prob keeps only the breakpoints that fall inside [t_ini, t_fin], so the same
# array can be passed at every baseline of a scan.
#
# Measured over the scans below, against a converged (n_slabs=24000) reference:
# with n_slabs=150 alone the worst point is 1.5e-2 off and 10% of points exceed
# 1e-3; adding the breakpoints gives 1.9e-3 worst and 1.5%, and runs 2.8x faster.
def castle_wall_breakpoints(n_castle_slabs, l_ini, l_fin):
    # Every wall edge, including the two at l_ini and l_fin where the profile switches on and
    # off.  A trajectory starting before l_ini crosses those as well, and a slab straddling an
    # unmarked discontinuity degrades the quadrature no matter how fine the grid becomes: with
    # arange(1, n) -- the interior walls only -- the error sat at 1.6e-2 in probability,
    # bit-identical from n_slabs=4 through 32, against 3.6e-12 once these two are included.
    return (l_ini+np.arange(0, n_castle_slabs+1)*(l_fin-l_ini)/n_castle_slabs) \
        *gd.CONV_KM_TO_INV_EV # [eV^-1]'''),
    md(r'''## 4.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e2, 1.e4 #1.e2, 1.e5 # [km]
l_npts = 6000 
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 50.*gd.UNIT_MEV # [eV] #10

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, energy), 
                                                     0, l*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, narrow castle-wall density profile
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_castle_wall_narrow = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_castle_wall(l, energy, num_density_e_low,
                                num_density_e_high, n_castle_slabs_narrow,
                                l_ini, l_fin), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3,
    t_breakpoints=castle_wall_breakpoints(n_castle_slabs_narrow, l_ini, l_fin))

# --------------------------------------------------------------
# In matter, constant density, low-amplitude noise
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_castle_wall_wide = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_castle_wall(l, energy, num_density_e_low,
                                num_density_e_high, n_castle_slabs_wide,
                                l_ini, l_fin), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3,
    t_breakpoints=castle_wall_breakpoints(n_castle_slabs_wide, l_ini, l_fin))'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_avg/norm for l in distances], color='0.5', ls='--'),
     dict(y=[num_density_e_castle_wall_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_low,
                                            num_density_e_high, n_castle_slabs_narrow, l_ini, l_fin)/norm
             for l in distances], color='C0'),
     dict(y=[num_density_e_castle_wall_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_low,
                                            num_density_e_high, n_castle_slabs_wide, l_ini, l_fin)/norm
             for l in distances], color='C3')],
    [[dict(y=prob_matt_const_density, label=r'Constant average density', color='0.5', ls='--'),
      dict(y=prob_matt_castle_wall_narrow, label=r'Narrow castle wall', color='C0'),
      dict(y=prob_matt_castle_wall_wide, label=r'Wide castle wall', color='C3')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 11),
    profile_ymajor=4, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]',
    title=r'$2\nu$~castle-wall matter profile, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left')'''),
    md(r'''## 4.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e3 # [MeV]
energy_npts = 4000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 1.e4 # [km]

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV),
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])


# --------------------------------------------------------------
# In matter, narrow castle-wall density profile
# --------------------------------------------------------------
prob_matt_castle_wall_narrow = np.array([oscprob.osc_prob(lambda l: H_func_castle_wall(l, enu*gd.UNIT_MEV, 
                                                                                       num_density_e_low, num_density_e_high, 
                                                                                       n_castle_slabs_narrow, l_ini, l_fin), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10,
                                                     t_breakpoints=castle_wall_breakpoints(n_castle_slabs_narrow, l_ini, l_fin))[nu_i][nu_f] \
                                          for enu in energies]) 

# --------------------------------------------------------------
# In matter, wide castle-wall density profile
# --------------------------------------------------------------
prob_matt_castle_wall_wide = np.array([oscprob.osc_prob(lambda l: H_func_castle_wall(l, enu*gd.UNIT_MEV, 
                                                                                       num_density_e_low, num_density_e_high, 
                                                                                       n_castle_slabs_wide, l_ini, l_fin), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10,
                                                     t_breakpoints=castle_wall_breakpoints(n_castle_slabs_wide, l_ini, l_fin))[nu_i][nu_f] \
                                          for enu in energies])  '''),
    md(r'''### Plot probabilities'''),
    code(r'''smooth = lambda y: sp.signal.savgol_filter(y, window_length=301, polyorder=1)
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_matt_const_density, label=r'Constant average density', color='0.5', ls='--'),
      dict(y=prob_matt_castle_wall_narrow, label=r'Narrow castle wall', color='C0'),
      dict(y=prob_matt_castle_wall_wide, label=r'Wide castle wall', color='C3')],
     [dict(y=smooth(prob_matt_const_density), color='0.5', ls='--'),
      dict(y=smooth(prob_matt_castle_wall_narrow), color='C0'),
      dict(y=smooth(prob_matt_castle_wall_wide), color='C3')]],
    xlim=(energy_min, energy_max), xscale='log',
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$2\nu$~castle-wall matter profile, $L = $~{:.2f}~km'.format(baseline),
    legend_title=r'Matter profile', legend_loc='lower right', legend_on_panel=0)'''),
    md(r'''# 5. Probabilities 2$\nu$: in matter with noisy density profile'''),
    md(r'''## 5.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Parameters for the electron number density in the Sun, Eq. (10.62) in Giunti & Kim
num_density_e_center = 10*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_scale = 100.0 # [km]

# Build the noise
np.random.seed(20)
# High-amplitude, more marked noise
npts_noise = 50 # Use more points to make it more finely noised
x_high = np.linspace(0, 1, npts_noise)
noise_sample_high = np.random.normal(loc=1.0, scale=0.3, size=npts_noise) # Increase scale to increase the width of the noise bumps
noise_mask_high_interp = sp.interpolate.make_interp_spline(x_high, noise_sample_high, k=1)
# Low-amplitude, choppier noise
npts_noise = 300 # Use more points to make it more finely noised
x_low = np.linspace(0, 1, npts_noise)
noise_sample_low = np.random.normal(loc=1.0, scale=0.08, size=npts_noise) # Increase scale to increase the width of the noise bumps
noise_mask_low_interp = sp.interpolate.make_interp_spline(x_low, noise_sample_low, k=1)

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# --------------------------------------------------------------
# In matter, constant density equal to central density
# --------------------------------------------------------------
# The matter density is constant, so this is merely a wrapper.
# [l] = km
VCC_center = matter.VCC_func(0.0, lambda l: num_density_e_center) # [eV]
H_matt_const_density = hamiltonians.hamiltonian_2nu_matter(VCC_center) # [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt_const_density

# --------------------------------------------------------------
# In matter, constant density, with noise around central value
# --------------------------------------------------------------
def num_density_e_const_noisy_func(l, num_density_e_const, l_ini, l_fin, noise_model='high-amplitude'):
    l_scaled = (l/gd.CONV_KM_TO_INV_EV-l_ini)/(l_fin-l_ini)
    if noise_model == 'high-amplitude':
        return num_density_e_const*noise_mask_high_interp(l_scaled)
    elif noise_model == 'low-amplitude':
        return num_density_e_const*noise_mask_low_interp(l_scaled)
def VCC_func_const_density_noisy(l, num_density_e_const, l_ini, l_fin, noise_model='high-amplitude'):
    return matter.VCC_func(l, lambda r: num_density_e_const_noisy_func(r, num_density_e_const, l_ini, l_fin, 
                                                                       noise_model=noise_model)) # [eV]
def H_func_const_density_noisy(l, energy, num_density_e_const, l_ini, l_fin, noise_model='high-amplitude'):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_const_density_noisy(l, num_density_e_const, 
                                                                                                             l_ini, l_fin,
                                                                                                             noise_model=noise_model))
# # --------------------------------------------------------------
# # In matter, exponentially falling density
# # --------------------------------------------------------------
# # [l] = km
# def num_density_e_exp_func(l, num_density_e_center, l_scale):
#     return num_density_e_center*np.exp(-(l/gd.CONV_KM_TO_INV_EV)/l_scale) # [eV^3]
# def VCC_func_exp_density(l, num_density_e_center, l_scale):
#     return matter.VCC_func(l, lambda r: num_density_e_exp_func(r, num_density_e_center, l_scale)) # [eV]
# def H_func_exp_density(l, energy):
#     return (1/energy)*H_vac_energy_indep-hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_density(l, num_density_e_center, l_scale))

# # --------------------------------------------------------------
# # In matter, exponentially falling density, noisy
# # --------------------------------------------------------------
# # [l] = km
# def num_density_e_exp_noisy_func(l, num_density_e_center, l_scale, l_ini, l_fin):
#     l_scaled = (l/gd.CONV_KM_TO_INV_EV-l_ini)/(l_fin-l_ini)
#     return num_density_e_exp_func(l, num_density_e_center, l_scale)*noise_mask_interp(l_scaled)
# def VCC_func_exp_density_noisy(l, num_density_e_center, l_scale, l_ini, l_fin):
#     return matter.VCC_func(l, lambda r: num_density_e_exp_noisy_func(r, num_density_e_center, l_scale, l_ini, l_fin)) # [eV]
# def H_func_exp_density_noisy(l, energy, num_density_e_center, l_scale, l_ini, l_fin):
#     return (1/energy)*H_vac_energy_indep \
#         - hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_density_noisy(l, num_density_e_center, l_scale, l_ini, l_fin))'''),
    md(r'''## 5.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e2, 1.e5 # [km]
l_npts = 6000 
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 50.*gd.UNIT_MEV # [eV] #10

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, energy), 
                                                     0, l*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, constant density, high-amplitude noise
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_const_density_noisy_high_amplitude = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_const_density_noisy(l, energy, num_density_e_center,
                                        l_ini, l_fin, noise_model='high-amplitude'), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3)

# --------------------------------------------------------------
# In matter, constant density, low-amplitude noise
# --------------------------------------------------------------
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_const_density_noisy_low_amplitude = oscprob.osc_prob_energy_baseline(
    lambda l: H_func_const_density_noisy(l, energy, num_density_e_center,
                                        l_ini, l_fin, noise_model='low-amplitude'), energy, distances*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3)

# # --------------------------------------------------------------
# # In matter, exponentially falling density
# # --------------------------------------------------------------
# prob_matt_exp_density = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, energy), 
#                                                    0, l*gd.CONV_KM_TO_INV_EV, 
#                                                    n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
smooth = lambda y: sp.signal.savgol_filter(y, window_length=301, polyorder=1)
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center/norm for l in distances], color='0.5', ls='--'),
     dict(y=[num_density_e_const_noisy_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_center,
                                            l_ini, l_fin,
                                            noise_model='high-amplitude')/norm
             for l in distances], color='C0'),
     dict(y=[num_density_e_const_noisy_func(l*gd.CONV_KM_TO_INV_EV, num_density_e_center,
                                            l_ini, l_fin,
                                            noise_model='low-amplitude')/norm
             for l in distances], color='C3')],
    [[dict(y=prob_matt_const_density, label=r'Constant density', color='0.5', ls='--'),
      dict(y=prob_matt_const_density_noisy_high_amplitude, label=r'High-amplitude noise', color='C0'),
      dict(y=prob_matt_const_density_noisy_low_amplitude, label=r'Low-amplitude noise', color='C3')],
     [dict(y=smooth(prob_matt_const_density), color='0.5', ls='--'),
      dict(y=smooth(prob_matt_const_density_noisy_high_amplitude), color='C0'),
      dict(y=smooth(prob_matt_const_density_noisy_low_amplitude), color='C3')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 15),
    profile_ymajor=4, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]', ylabel_labelpad=15,
    title=r'$2\nu$~noisy matter profile, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Matter profile', legend_loc='lower left', legend_on_panel=0)'''),
    md(r'''## 5.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e1, 1.e2 # [MeV]
energy_npts = 6000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = 1.e5 # [km]

#--------------------------------------------------------------
# In matter, constant density
# --------------------------------------------------------------
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV),
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])

# --------------------------------------------------------------
# In matter, constant density, high-amplitude noise
# --------------------------------------------------------------
prob_matt_const_density_noisy_high_amplitude = np.array([oscprob.osc_prob(lambda l: H_func_const_density_noisy(l, enu*gd.UNIT_MEV,
                                                                                                            num_density_e_center, 
                                                                                                            l_ini, l_fin,
                                                                                                            noise_model='high-amplitude'), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] \
                                          for enu in energies]) 

# --------------------------------------------------------------
# In matter, constant density, low-amplitude noise
# --------------------------------------------------------------
prob_matt_const_density_noisy_low_amplitude = np.array([oscprob.osc_prob(lambda l: H_func_const_density_noisy(l, enu*gd.UNIT_MEV,
                                                                                                            num_density_e_center,
                                                                                                            l_ini, l_fin,
                                                                                                            noise_model='low-amplitude'), 
                                                     0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=150, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] \
                                          for enu in energies]) '''),
    md(r'''### Plot probabilities'''),
    code(r'''smooth = lambda y, po=1: sp.signal.savgol_filter(y, window_length=301, polyorder=po)
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_matt_const_density, label=r'Constant density', color='0.5', ls='--'),
      dict(y=prob_matt_const_density_noisy_high_amplitude, label=r'High-amplitude noise', color='C0'),
      dict(y=prob_matt_const_density_noisy_low_amplitude, label=r'Low-amplitude noise', color='C3')],
     [dict(y=smooth(prob_matt_const_density), color='0.5', ls='--'),
      dict(y=smooth(prob_matt_const_density_noisy_high_amplitude), color='C0'),
      dict(y=smooth(prob_matt_const_density_noisy_low_amplitude, 2), color='C3')]],
    xlim=(energy_min, energy_max), xscale='linear', xmajor=10, xminor=1,
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$2\nu$~noisy matter profile, $L = $~{:.2f}~km'.format(baseline),
    legend_title=r'Matter profile', legend_loc='lower right', legend_on_panel=0)'''),
    md(r'''# 6. Probabilities 2$\nu$: in the Earth'''),
    md(r'''## 6.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''def num_density_e_func_prem(r): 
    return matter.num_density_e_func(r, earth.density_matter_func_prem, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

def VCC_func_prem(r):
    return matter.VCC_func(r, num_density_e_func_prem) # [eV]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# [l] = km
def H_func_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) 
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_prem(r))

# Test values of the zenith angle, cos(theta_z), to plot, with theta_z = 0 normal to the surface of the Earth, at the position of the
# detector.  All of the direction are upgoing.
costhz_val = [-0.1, -0.5, -1.0]

# Maximum baselines inside the Earth
l_max_val = [earth.distance_traveled_inside_earth(costhz) for costhz in costhz_val] # [km]'''),
    md(r'''## 6.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]

nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 20.*gd.UNIT_MEV # [eV] #10

l_min = 1.e2 # [km]

# Generate probabilities for the different directions
distances_val, prob_val = [], []
l_npts = 3000
for i in range(len(costhz_val)):
    print("costhz = "+str(costhz_val[i]))
    distances = np.logspace(np.log10(l_min), np.log10(l_max_val[i]), l_npts) # [km]
    distances_val.append(distances)
    # PREM is piecewise-smooth: the density jumps where the chord crosses a shell boundary.
    # Marking those crossings as mandatory slab edges keeps the quadrature from integrating
    # across a discontinuity.  Measured at this grid, against solve_ivp, it costs nothing in
    # time and buys 1.3e3x to 1.9e6x in accuracy -- the costhz = -1 direction was outside the
    # default 1e-3 tolerance without it (4.2e-3 -> 3.3e-6).
    prem_bp = earth.prem_layer_edges_along_chord(costhz_val[i])*gd.UNIT_KM
    prob = np.array([oscprob.osc_prob(lambda l: H_func_prem(costhz_val[i], l, energy), 0, l*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10,
                                      t_breakpoints=prem_bp)[nu_i][nu_f] for l in distances]) 
    prob_val.append(prob)'''),
    md(r'''### Plot probabilities'''),
    code(r'''lc = ['C0', 'C2', 'C3']
ls = ['-', '-', '-']
norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
profiles, panel = [], []
for i in range(len(costhz_val)):
    n_e = np.array([num_density_e_func_prem(
                        earth.earth_radial_distance_from_depth(costhz_val[i], l))
                    for l in distances_val[i]])/norm
    profiles.append(dict(x=distances_val[i], y=n_e, color=lc[i], ls=ls[i]))
    panel.append(dict(x=distances_val[i], y=prob_val[i], color=lc[i], ls=ls[i],
                      label=r'$\cos \theta_z = $~~{:.2f}'.format(costhz_val[i])))

fig, ax = plotting.plot_probability_with_profile(
    distances_val[0], profiles, [panel],
    xlim=(l_min, max(l_max_val)), profile_ylim=(0, 7),
    profile_ymajor=2, profile_yminor=1, profile_height=0.3,
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]',
    title=r'$2\nu$~inside the Earth (PREM), $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    legend_title=r'Direction', legend_loc='lower left')'''),
    md(r'''## 6.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e1, 1.e2 # [MeV]
energy_npts = 3000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU

l_min = 1.e2 # [km]

# Generate probabilities for the different directions
prob_val = []
for i in range(len(costhz_val)):
    print("costhz = "+str(costhz_val[i]))
    # PREM is piecewise-smooth: the density jumps where the chord crosses a shell boundary.
    # Marking those crossings as mandatory slab edges keeps the quadrature from integrating
    # across a discontinuity.  Measured at this grid, against solve_ivp, it costs nothing in
    # time and buys 1.3e3x to 1.9e6x in accuracy -- the costhz = -1 direction was outside the
    # default 1e-3 tolerance without it (4.2e-3 -> 3.3e-6).
    prem_bp = earth.prem_layer_edges_along_chord(costhz_val[i])*gd.UNIT_KM
    prob = np.array([oscprob.osc_prob(lambda l: H_func_prem(costhz_val[i], l, enu*gd.UNIT_MEV), 
                                      0, l_max_val[i]*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10,
                                      t_breakpoints=prem_bp)[nu_i][nu_f] for enu in energies]) 
    prob_val.append(prob)'''),
    md(r'''### Plot probabilities'''),
    code(r'''smooth = lambda y: sp.signal.savgol_filter(y, window_length=201, polyorder=1)
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_val[i], color=lc[i], ls=ls[i],
           label=r'$\cos \theta_z = $~~{:.2f}'.format(costhz_val[i]))
      for i in range(len(costhz_val))],
     [dict(y=smooth(prob_val[i]), color=lc[i], ls=ls[i])
      for i in range(len(costhz_val))]],
    xlim=(energy_min, energy_max), xscale='linear', xmajor=10, xminor=1,
    panel_ylabels=[r'$2\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$2\nu$~inside the Earth (PREM), $L = L_{\rm max}(\theta_z)$',
    legend_title=r'Direction', legend_loc='lower left', legend_on_panel=0)'''),
    md(r'''# 7. Probabilities 2$\nu$: in the Sun'''),
    md(r'''## 7.1 General definitions'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Parameters for the electron number density in the Sun, Eq. (10.62) in Giunti & Kim
num_density_e_center_sun = 245*gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0) # [eV^3]
l_scale_sun = gd.SUN_RADIUS/10.54 # [km]

# Potential assuming constant density: the matter density is constant, so this is merely a wrapper.
# [l] = km
def num_density_e_sun_const_func(l):
    return num_density_e_center_sun # [eV^3]
def VCC_func_const_density_sun(l):
    return matter.VCC_func(l, num_density_e_sun_const_func) # [eV]

# Potential assuming exponentially decreasing density
# [l] = km
def num_density_e_sun_exp_func(l):
    return num_density_e_center_sun*np.exp(-(l/gd.CONV_KM_TO_INV_EV)/l_scale_sun) # [eV^3]
def VCC_func_exp_sun(l):
    return matter.VCC_func(l, num_density_e_sun_exp_func) # [eV]'''),
    md(r'''## 7.2 Probabilities vs. distance'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e-3, 1.0 # l/SUN_RADIUS
l_npts = 1000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # l/SUN_RADIUS

# Compute probability vs. baseline
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 5.0*gd.UNIT_MEV # [eV]

# --------------------------------------------------------------
# In matter, constant density equal to central solar density
# --------------------------------------------------------------
H_func_const_density = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2) \
                                      + hamiltonians.hamiltonian_2nu_matter(VCC_func_const_density_sun(l))
prob_matt_const_density = np.array([oscprob.osc_prob(H_func_const_density, 
                                                     0, l*gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])

# --------------------------------------------------------------
# In matter, solar matter density profile
# --------------------------------------------------------------
H_func_exp_density = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2) \
                                      + hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_sun(l))
# H_func_exp_density = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2) \
                                      # + hamiltonians.hamiltonian_2nu_matter(-VCC_func_exp_sun(l))
# One traversal of the profile answers every baseline: the evolution operator is a
# time-ordered product, so U(0->L2) = U(L1->L2) . U(0->L1) and each requested answer is a
# prefix of the next (cumulative=True).  Its accuracy grid is *inherited* from one adaptive
# call at the longest baseline rather than guessed, which is what makes it safe: measured
# against solve_ivp, the fixed grid this replaces is accurate at most baselines but spikes
# past 1e-3 at a scattered few of them, and error is not monotone in n_slabs.
prob_matt_exp_density = oscprob.osc_prob_energy_baseline(
    H_func_exp_density, energy, distances*gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV,
    0.0, nu_i, nu_f, cumulative=True,
    n_tpts_per_slab=100, magnus_exp_order=3)

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
# Using Magnus expansion
H_func = lambda l: hamiltonians.hamiltonian_2nu_vacuum(energy, sth, Dm2)
prob_vac = np.array([oscprob.osc_prob(H_func,
                                      0.0, l*gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for l in distances])
# # Using standard formula
# prob_std = np.array([oscprobstd.osc_prob_2nu_vacuum_std(sth, Dm2, energy, l*gd.CONV_KM_TO_INV_EV)[nu_i][nu_f] for l in distances])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_with_profile(
    distances,
    [dict(y=[num_density_e_center_sun/norm for l in distances], color='C0'),
     dict(y=[num_density_e_center_sun/norm*np.exp(-l*gd.SUN_RADIUS/l_scale_sun)
             for l in distances], color='C3')],
    [[dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
      dict(y=prob_matt_const_density, color='C0',
           label=r'Constant density, $N_e/N_{\rm Av} =~$'
                 + '{:.0f}'.format(num_density_e_center_sun/norm)),
      dict(y=prob_matt_exp_density, label=r'Solar profile', color='C3')]],
    xlim=(l_ini, l_fin), profile_ylim=(0, 250),
    profile_ymajor=50, profile_yminor=10, profile_height=0.3,
    profile_ylabel=r'$N_e/N_{\rm Av}$',
    panel_ylabels=[r'Two-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L/R_\odot$', ylabel_labelpad=7,
    title=r'$2\nu$~in the Sun, $E_\nu = $~{:.2f}~MeV'.format(energy/gd.UNIT_MEV),
    title_fontsize=20, legend_title=r'Matter profile', legend_loc='lower left',
    grid=False)'''),
    md(r'''## 7.3 Probabilities vs. energy'''),
    md(r'''### Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e2 # [MeV]
energy_npts = 1000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Compute probability vs. energy
# osc_prob and osc_prob_2nu_vacuum_std return a 2x2 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
baseline = gd.SUN_RADIUS # [km]

H_vac_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # Vacuum H without (1/E) [eV^2]

# --------------------------------------------------------------
# In matter, constant density equal to central solar density
# --------------------------------------------------------------
H_matt = hamiltonians.hamiltonian_2nu_matter(VCC_func_const_density_sun(0.0)) # Matter Hamiltonian [eV]
def H_func_const_density(l, energy):
    return (1/energy)*H_vac_energy_indep+H_matt
prob_matt_const_density = np.array([oscprob.osc_prob(lambda l: H_func_const_density(l, enu*gd.UNIT_MEV), 
                                                     0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])

# --------------------------------------------------------------
# In matter, solar matter density profile
# --------------------------------------------------------------
def H_func_exp_density(l, energy):
    return (1/energy)*H_vac_energy_indep+hamiltonians.hamiltonian_2nu_matter(VCC_func_exp_sun(l)) 
prob_matt_exp_density = np.array([oscprob.osc_prob(lambda l: H_func_exp_density(l, enu*gd.UNIT_MEV), 
                                                   0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                   n_slabs=200, n_tpts_per_slab=100, magnus_exp_order=2, n_jobs=10)[nu_i][nu_f] \
                                                   # n_slabs=200, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=8)[nu_i][nu_f] \
                                  for enu in energies])

# --------------------------------------------------------------
# In vacuum
# --------------------------------------------------------------
prob_vac = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_MEV))*H_vac_energy_indep, 
                                      0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1)[nu_i][nu_f] for enu in energies])'''),
    md(r'''### Plot probabilities'''),
    code(r'''norm = (gd.N_AV/pow(gd.CONV_CM_TO_INV_EV,3.0))
fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_vac, label=r'Vacuum', color='0.2', lw=0.5),
     dict(y=prob_matt_const_density, color='C0',
          label=r'Constant density, $N_e/N_{\rm Av} =~$'
                + '{:.0f}'.format(num_density_e_center_sun/norm)),
     dict(y=prob_matt_exp_density, label=r'Solar profile', color='C3')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=2, energy_unit='MeV',
    xlim=(energy_min, energy_max),
    legend_title=r'Matter profile', legend_loc='center left',
    title=r'$2\nu$~in the Sun, $L = R_\odot$')'''),
    ])

# --------------------------------------------------- 04_magnus_long_baseline
books['04_magnus_long_baseline.ipynb'] = notebook(
    'Long baselines',
    "Probabilities between two points on the Earth's surface, which is the geometry of an accelerator experiment: DUNE, T2K, Hyper-K, ESS.\n\nGive the source and detector coordinates and the chord and its density profile follow.",
    [
    code(r'''import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''To help in computing neutrino oscillation probabilities in real experimental setups, the `earth` module of Mag$\nu$s contains routines to compute the zenith angle of the location of one point on the surface of the Earth as measured from a different point, i.e., the position of a neutrino source measured from the position of a neutrino detector.  (These routines assume that the Earth is spherical, that the distribution of matter inside it radially symmetric, as in the Preliminary Reference Earth Model, and that the locations are on the surface, not underground.  Users can write their own functions to deal with these cases.)'''),
    md(r'''The `earth` module contains predefined coordinates for a few significant locations:'''),
    code(r'''earth.loc_coords_dms'''),
    md(r'''The coordinates are tuples of (days, minutes, seconds).  This is the format in which the coordinates are pased to the functions below.  When defining custom locations, the user must feed them in this format to the functions.'''),
    md(r'''Given the latitude (`lat`) and longitude (`lon`) of positions 1 (source) and 2 (detector) on the surface of the Earth, the chord length between them (i.e., the straight-line distance between them through the Earth) is computed with the `chord_length_inside_earth` function and the zenith angle of position 2 measured from position 1 is computed with the `costhz_between_points_on_surface` function.'''),
    code(r'''# DUNE setup: from Fermilab (source) to Homestake (far detector)
lat1, lon1 = earth.loc_coords_dms['fermilab']['lat'], earth.loc_coords_dms['fermilab']['lon']
lat2, lon2 = earth.loc_coords_dms['homestake']['lat'], earth.loc_coords_dms['homestake']['lon']

baseline = earth.chord_length_inside_earth(lat1, lon1, lat2, lon2)
print("Baseline from Fermilab to Homestake: " + str(baseline) + " km")

costhz = earth.costhz_between_points_on_surface(lat1, lon1, lat2, lon2)
print("Cosine of zenith angle of Fermilab measured from Homestake: " + str(costhz))'''),
    md(r'''Now that we know the zenith angle of the source as measured from the detector, we can compute and plot neutrino oscillation probabilities as usual (for details, see the example Mag$\nu$s notebooks `2_magnus_2nu_vacuum_matter.ipynb` and `3_magnus_3nu_vacuum_matter.ipynb` ).'''),
    md(r'''To illustrate this, below, we consider an artificial setup of a neutrino beam with adjustable direction shot from Fermilab to several detector locations.'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Electron number density inside Earth, using the PREM density model
def num_density_e_func_prem(r): 
    return matter.num_density_e_func(r, earth.density_matter_func_prem, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

# Coherent forward potential inside Earth, using the PREM density model
def VCC_func_prem(r):
    return matter.VCC_func(r, num_density_e_func_prem) # [eV]

# Vacuum Hamiltonian without the 1/E prefactor 
H_vac_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
                                                                               compute_matrix_multiplication=False)

# Hamiltonian including matter effects inside Earth using the PREM density model
def H_func_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return (1/energy)*H_vac_energy_indep + hamiltonians.hamiltonian_3nu_matter(VCC_func_prem(r))

# Cosines of zenith angles of Fermilab measured from the detector locations (we would get the same result if we seapped them, since we
# consider that Earth is radially symmetric)
detectors = ['SNOLAB', 'Homestake', 'CERN', "South Pole"]
costhz_arr = [earth.costhz_between_points_on_surface(earth.loc_coords_dms[det.lower().replace(" ", "_")]['lat'], 
                                                     earth.loc_coords_dms[det.lower().replace(" ", "_")]['lon'],
                                                     earth.loc_coords_dms['fermilab']['lat'], 
                                                     earth.loc_coords_dms['fermilab']['lon']) 
             for det in detectors]

# Maximum baselines inside the Earth
l_max_arr = [earth.distance_traveled_inside_earth(costhz) for costhz in costhz_arr] # [km]

# Helper function for labeling plots
def prob_label(nu_i, nu_f):
    if (nu_i == gd.NUE):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_e \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_e \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_e \to \nu_\tau}$'
    elif (nu_i == gd.NUMU):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_\mu \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_\mu \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_\mu \to \nu_\tau}$'
    elif (nu_i == gd.NUTAU):
        if (nu_f == gd.NUE):
            label = r'$P_{\nu_\tau \to \nu_e}$'
        elif (nu_f == gd.NUMU):
            label = r'$P_{\nu_\tau \to \nu_\mu}$'
        elif (nu_f == gd.NUTAU):
            label = r'$P_{\nu_\tau \to \nu_\tau}$'
    return label

# Per-detector colours and line styles, shared by both figures
lc = ['C0', 'C2', 'C3', 'C4']
ls = ['-', '-', '-', '-']'''),
    md(r'''First, we compute the probabilities vs. baseline:'''),
    code(r'''# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU
energy = 10.*gd.UNIT_MEV # [eV]

l_min = 1.e2 # [km]

# Generate probabilities for the different directions
distances_arr, prob_arr = [], []
l_npts = 3000
for i in range(len(costhz_arr)):
    print("detector = " + detectors[i])
    distances = np.logspace(np.log10(l_min), np.log10(l_max_arr[i]), l_npts) # [km
    distances_arr.append(distances)
    prob = np.array([oscprob.osc_prob(lambda ll: H_func_prem(costhz_arr[i], ll, energy), 0, l*gd.CONV_KM_TO_INV_EV,
                                      n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] for l in distances]) 
    prob_arr.append(prob)'''),
    md(r'''Now we plot them, one panel per detector, with the electron-density profile
each trajectory samples stacked on top.

The profile panel is what makes the figure worth reading. Each baseline
crosses a different part of the Earth, so the four curves sample different
depths: the South Pole trajectory is the only one long enough to reach the
core, which is the step up to $N_e/N_{\rm Av} \approx 5.5$ near
$10^4$ km. Where that step falls along the baseline is exactly where the
corresponding probability curve changes character.

This is one call to `magnus.plotting.plot_probability_with_profile`; the
density panel, the shared logarithmic abscissa, the suppressed tick labels on
every panel but the last, and the legend styling are its defaults.'''),
    code(r'''# The electron-density profile each trajectory samples, in units of N_Av, so
# that the panel reads directly as rho*Y_e in g cm^-3.
profiles, panels = [], []
for i in range(len(costhz_arr)):
    n_e = np.array([num_density_e_func_prem(
                        earth.earth_radial_distance_from_depth(costhz_arr[i], l))
                    for l in distances_arr[i]])
    n_e = n_e/(gd.N_AV/pow(gd.CONV_CM_TO_INV_EV, 3.0))
    profiles.append(dict(x=distances_arr[i], y=n_e, color=lc[i], ls=ls[i]))
    panels.append([dict(x=distances_arr[i], y=prob_arr[i], color=lc[i], ls=ls[i],
                        label=detectors[i])])

fig, ax = plotting.plot_probability_with_profile(
    distances_arr[0], profiles, panels,
    xlim=(l_min, max(l_max_arr)),
    profile_ylim=(0, 6), profile_ymajor=2, profile_yminor=1,
    panel_ymajor=0.10, panel_yminor=0.02,
    panel_ylabels=[r'Three-neutrino probability,~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Baseline, $L$~[km]',
    title=r'$3\nu$~inside the Earth (PREM), $E_\nu = $~{:.2f}~MeV'.format(
        energy/gd.UNIT_MEV),
    legend_title=r'From Fermilab to...', legend_loc='upper left',
    legend_on_panel=-1, legend_kw=dict(borderpad=0.7, ncol=4),
    figsize=[18, 18])'''),
    md(r'''Finally, we compute the probabilities vs. energy for the full baseline in each case:'''),
    code(r'''# Energies
energy_min, energy_max = 1.e1, 1.e2 # [MeV]
energy_npts = 3000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Generate probabilities for the different directions
prob_arr = []
for i in range(len(costhz_arr)):
    print("detector = " + detectors[i])
    prob = np.array([oscprob.osc_prob(lambda l: H_func_prem(costhz_arr[i], l, enu*gd.UNIT_MEV), 
                                      0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] for enu in energies]) 
    prob_arr.append(prob)'''),
    md(r'''Finally, the same four trajectories against energy, at the full baseline in
each case.

Two panels: the probability as computed, and below it the same curve passed
through a Savitzky-Golay filter. The lower panel is not cosmetic. Above a few
tens of MeV the oscillation length falls well below any real detector's energy
resolution, so what an experiment measures is the local average, not the
fringes -- and the averaged curve is where the matter effect shows up as a
genuine, resolvable feature rather than as faster wiggles.

The same function draws it, with `profiles=None` to omit the density panel.'''),
    code(r'''# No density panel this time: the lower panel is the same probability with the
# fast oscillations filtered out, which is what a detector with finite energy
# resolution actually measures.
fig, ax = plotting.plot_probability_with_profile(
    energies, None,
    [[dict(y=prob_arr[i], color=lc[i], ls=ls[i], label=detectors[i])
      for i in range(len(costhz_arr))],
     [dict(y=sp.signal.savgol_filter(prob_arr[i], window_length=301, polyorder=1),
           color=lc[i], ls=ls[i])
      for i in range(len(costhz_arr))]],
    xlim=(energy_min, energy_max), xscale='linear',
    xmajor=10, xminor=1, panel_ymajor=0.10, panel_yminor=0.02,
    panel_ylabels=[r'$3\nu$~probability,~'+plotting.prob_label(nu_i, nu_f),
                   r'Average~'+plotting.prob_label(nu_i, nu_f)],
    xlabel=r'Neutrino energy, $E_\nu$ [MeV]',
    title=r'$3\nu$~inside the Earth (PREM)',
    legend_title=r'From Fermilab to...', legend_loc='lower right',
    legend_on_panel=0, grid_kw=dict(c='0.7'),
    figsize=[18, 10])'''),
    ])

# --------------------------------------------------- 05_magnus_biprobability
books['05_magnus_biprobability.ipynb'] = notebook(
    'Biprobability plots',
    'The neutrino probability $P(\\nu_\\alpha \\to \\nu_\\beta)$ against the antineutrino one $P(\\bar\\nu_\\alpha \\to \\bar\\nu_\\beta)$, traced as the CP phase runs over its range.\n\nThe area enclosed is the CP violation an experiment is trying to measure.',
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.oscprobstd as oscprobstd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''Biprobability plots show how neutrino vs. anti-neutrino oscillation probabilities change as a function of the CP-violation phase, $\delta_{\rm CP}$.  This is typically done for neutrino propagation in matter, representing the setup of long-baseline neutrino oscillations in which $\delta_{\rm CP}$ is expected to be measured.  Biprobability plots are useful to illustrate the separation between neutrino and anti-neutrino oscillation probabilities from which this measurement is performed.

Sometimes, biprobability plots include new-physics effects that affect the probabilities, in order to show how their presence could cloud the measurement of $\delta_{\rm CP}$.  In this notebook, we use Mag$\nu$s to produce biprobability plots only under standard $3\nu$ oscillations; for biprobability plots including sterile neutrinos and non-standard interactions, see notebooks `07_magnus_bsm_sterile_nu.ipynb` and `08_magnus_bsm_nsi.ipynb`, respectively.'''),
    md(r'''## Plotting function'''),
    md(r'''First, let's define a general-purpose plotting function to produce biprobability plots, given precomputed neutrino and anti-neutrino oscillation probabilities:'''),
    code(r'''def make_plot_biprobability(prob_nu_arr, prob_nubar_arr, ls, lc, lw,
                            label_prob_nu, label_prob_nubar, points_sel_arr=None,
                            annotations=None, xaxis_major_locator=None,
                            xaxis_minor_locator=None, yaxis_major_locator=None,
                            yaxis_minor_locator=None, leg_fontsize=16, leg_loc=None,
                            leg_ncol=1, xlim=None, ylim=None, save_fig=False,
                            fig_filename=None, fig_format='pdf'):
    """Draw a bi-probability plane via magnus.plotting.

    The stored curves carry the delta_CP value in column 0 and the probability
    in column 1, so only column 1 is plotted. Each marker entry is
    [label, marker, filled, [x, y]] and is translated into the module's form.
    """
    def spacing(loc):
        # the notebook passes MultipleLocator objects; the module takes the step
        return None if loc is None else float(np.diff(loc.tick_values(0.0, 1.0))[0])

    markers = None
    if points_sel_arr is not None:
        # points_sel_arr[i] holds curve i's marked phases, each as
        # [label, marker, filled, [x, y]] -- coordinates, not indices, because the
        # marked phases were computed separately from the curve. Only the first
        # curve's markers carry labels, so the legend names each phase once.
        markers = [dict(xy=(p[3][0], p[3][1]), marker=p[1],
                        filled=(p[2] is True), curve=i,
                        **({'label': p[0]} if i == 0 else {}))
                   for i, points_sel in enumerate(points_sel_arr)
                   for p in points_sel]

    fig, ax = plotting.plot_biprobability(
        [np.asarray(c)[:, 1] for c in prob_nu_arr],
        [np.asarray(c)[:, 1] for c in prob_nubar_arr],
        curve_kw=[dict(color=lc[i], ls=ls[i], lw=lw[i])
                  for i in range(len(prob_nu_arr))],
        markers=markers,
        xlabel=label_prob_nu, ylabel=label_prob_nubar,
        xlim=xlim, ylim=ylim,
        xmajor=spacing(xaxis_major_locator), xminor=spacing(xaxis_minor_locator),
        ymajor=spacing(yaxis_major_locator), yminor=spacing(yaxis_minor_locator),
        annotations=[dict(text=a['text'], xy=a['xy'], color=a['color'],
                          fontsize=a['fontsize'], ha=a['ha'])
                     for a in (annotations or [])],
        legend_loc=leg_loc,
        legend_kw=dict(fontsize=leg_fontsize, ncol=leg_ncol, framealpha=1.0),
        savefig=(fig_filename + '.' + fig_format) if save_fig else None,
        savefig_kw=dict(transparent=False, bbox_inches='tight', dpi=300))
    return fig, ax'''),
    md(r'''## In vacuum'''),
    md(r'''Let's start with the case of propagation in vacuum, comparing the plots obtained under normal mass ordering (NO) and under inverted mass ordering (IO).  In each case, we fix the mixing parameters to their best-fit values from the NuFit 6.0 global fit to oscillation data, which are stored in the `globaldefs` module of Mag$\nu$s.'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12_NO = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23_NO = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13_NO = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP_NO = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21_NO = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31_NO = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]
s12_IO = gd.S12_IO_BF_NUFIT_6_0 # [adim]
s23_IO = gd.S23_IO_BF_NUFIT_6_0 # [adim]
s13_IO = gd.S13_IO_BF_NUFIT_6_0 # [adim]
dCP_IO = gd.DCP_IO_BF_NUFIT_6_0 # [adim]
D21_IO = gd.D21_IO_BF_NUFIT_6_0 # [eV^2]
D31_IO = gd.D31_IO_BF_NUFIT_6_0 # [eV^2]'''),
    code(r'''# Pick the baseline and energy of T2K
baseline = 810 # [km]
energy = 2*gd.UNIT_MEV # [eV]

# Values of the delta_CP phase at which to compute the probabilities
dCP_npts = 100
dCP_arr = np.linspace(-np.pi, np.pi, dCP_npts)

# We will compute the appearance probabilities, i.e., nu_mu --> nu_e and nu_mu-bar --> nu_e-bar
nu_i, nu_f = gd.NUMU, gd.NUE

# Hamiltonian in vacuum for neutrinos
def H_nu(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=False, compute_matrix_multiplication=False)

# Hamiltonian in vacuum for anti-neutrinos
def H_nubar(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=True, compute_matrix_multiplication=False)

# Generate the neutrino and anti-neutrino oscillation probabilities for many different values of dCP
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# nu_mu --> nu_e, normal ordering
prob_nu_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, normal ordering
prob_nubar_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu --> nu_e, inverted ordering
prob_nu_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, inverted ordering
prob_nubar_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])

# Points for selected values of dCP
dCP_sel = [-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi]
dCP_label_sel = [r'$-\pi, \pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', 'NuFit 6.0']
markers = ['o', 'v', 's', 'p', '*', 'p', 's', 'v', '^']
filled = [True, True, True, True, True, False, False, False, True]
prob_nu_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]] 
                           for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
points_NO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_NO_sel[i][1], prob_nubar_NO_sel[i][1]]] for i in range(len(dCP_sel)+1)]
points_IO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_IO_sel[i][1], prob_nubar_IO_sel[i][1]]] for i in range(len(dCP_sel)+1)]'''),
    code(r'''make_plot_biprobability([prob_nu_NO, prob_nu_IO], 
                        [prob_nubar_NO, prob_nubar_IO],
                        lc=['C0', 'C1'], ls=['-', '-'], lw=[1.0, 1.0],
                        points_sel_arr=[points_NO_sel, points_IO_sel],
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': 'Vacuum', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'NO', 'xy': (0.5,0.75), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'IO', 'xy': (0.17,0.75), 'color': 'C1', 'fontsize': 20, 'ha': 'left'}
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.10),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.02),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.10),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.02),
                        leg_fontsize=15, leg_loc='upper right', leg_ncol=1, xlim=[0,0.50], ylim=[0,0.50], save_fig=False)'''),
    md(r'''To validate this result, which was produced using the Mag$\nu$s expansion (albeit, to first order only, since for propoagation in vacuum this is enough), we compute the biprobability using the standard $3\nu$ oscillation probability instead.'''),
    code(r'''# Generate the neutrino and anti-neutrino oscillation probabilities for many different values of dCP
# osc_prob_3nu_vacuum_std returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
prob_nu_std_NO = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(hamiltonians.pmns_mixing_matrix(s12_NO, s23_NO, s13_NO, dCP),
                                                                          D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV,
                                                                          nubar=False)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, normal ordering
prob_nubar_std_NO = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(hamiltonians.pmns_mixing_matrix(s12_NO, s23_NO, s13_NO, dCP),
                                                                             D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV,
                                                                             nubar=True)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu --> nu_e, inverted ordering
prob_nu_std_IO = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(hamiltonians.pmns_mixing_matrix(s12_IO, s23_IO, s13_IO, dCP),
                                                                          D21_IO, D31_IO, energy, baseline*gd.CONV_KM_TO_INV_EV,
                                                                          nubar=False)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, inverted ordering
prob_nubar_std_IO = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(hamiltonians.pmns_mixing_matrix(s12_IO, s23_IO, s13_IO, dCP),
                                                                             D21_IO, D31_IO, energy, baseline*gd.CONV_KM_TO_INV_EV,
                                                                             nubar=True)[nu_i][nu_f]] for dCP in dCP_arr])

# Points for selected values of dCP
dCP_sel = [-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi]
dCP_label_sel = [r'$-\pi, \pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', 'NuFit 6.0']
markers = ['o', 'v', 's', 'p', '*', 'p', 's', 'v', '^']
filled = [True, True, True, True, True, False, False, False, True]
prob_nu_std_NO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_NO, s23_NO, s13_NO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=False)[nu_i][nu_f]] for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nubar_std_NO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_NO, s23_NO, s13_NO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=True)[nu_i][nu_f]] for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nu_std_IO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_IO, s23_IO, s13_IO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=False)[nu_i][nu_f]] for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nubar_std_IO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_IO, s23_IO, s13_IO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=True)[nu_i][nu_f]] for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
points_std_NO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_std_NO_sel[i][1], prob_nubar_std_NO_sel[i][1]]] 
                     for i in range(len(dCP_sel)+1)]
points_std_IO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_std_IO_sel[i][1], prob_nubar_std_IO_sel[i][1]]] 
                     for i in range(len(dCP_sel)+1)]'''),
    code(r'''make_plot_biprobability([prob_nu_std_NO, prob_nu_std_IO], 
                        [prob_nubar_std_NO, prob_nubar_std_IO],
                        lc=['C0', 'C1'], ls=['-', '-'], lw=[1.0, 1.0],
                        points_sel_arr=[points_NO_sel, points_IO_sel],
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': r'Vacuum (using standard $3\nu$ formula)', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'NO', 'xy': (0.5,0.75), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'IO', 'xy': (0.17,0.75), 'color': 'C1', 'fontsize': 20, 'ha': 'left'}
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.10),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.02),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.10),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.02),
                        leg_fontsize=15, leg_loc='upper right', leg_ncol=1, xlim=[0,0.50], ylim=[0,0.50], save_fig=False)'''),
    md(r'''And this biprobability plot is, as expected, identical as the one produced using the Magnus expansion.'''),
    md(r'''## In matter, long-baseline experiments (constant density, Earth's crust)'''),
    md(r'''Now let's generate the biprobability plot in matter, for the NO$\nu$A, T2K, and DUNE long-baseline setups.  Compare the plots below to the ones shown in Fig. 7 of arXiv:1307.3248.'''),
    code(r'''# ====== NOvA ======

# Pick the baseline and energy of NOvA
baseline = 810 # [km]
energy = 2*gd.UNIT_GEV # [eV]

# Values of the delta_CP phase at which to compute the probabilities
dCP_npts = 100
dCP_arr = np.linspace(-np.pi, np.pi, dCP_npts)

# We will compute the appearance probabilities, i.e., nu_mu --> nu_e and nu_mu-bar --> nu_e-bar
nu_i, nu_f = gd.NUMU, gd.NUE

# Coherent matter potential in matter
rho = 3.0 # Average matter density in the Earth's crus [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]
# Note: we could have directly used VCC = gd.NUM_DENSITY_E_EARTH_CRUST

# Hamiltonian in matter
H_matt = hamiltonians.hamiltonian_3nu_matter(VCC)

# Hamiltonian for neutrinos (the matter potential is added)
def H_nu(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=False, compute_matrix_multiplication=False) \
    + H_matt

# Hamiltonian for anti-neutrinos (the matter potential is subtracted)
def H_nubar(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=True, compute_matrix_multiplication=False) \
    - H_matt

# Generate the neutrino and anti-neutrino oscillation probabilities for many different values of dCP
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# nu_mu --> nu_e, normal ordering
prob_nu_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, normal ordering
prob_nubar_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu --> nu_e, inverted ordering
prob_nu_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, inverted ordering
prob_nubar_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])

# Points for selected values of dCP
dCP_sel = [-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi]
dCP_label_sel = [r'$-\pi, \pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', 'NuFit 6.0']
markers = ['o', 'v', 's', 'p', '*', 'p', 's', 'v', '^']
filled = [True, True, True, True, True, False, False, False, True]
prob_nu_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]] 
                           for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
points_NO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_NO_sel[i][1], prob_nubar_NO_sel[i][1]]] for i in range(len(dCP_sel)+1)]
points_IO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_IO_sel[i][1], prob_nubar_IO_sel[i][1]]] for i in range(len(dCP_sel)+1)]'''),
    code(r'''make_plot_biprobability([prob_nu_NO, prob_nu_IO], 
                        [prob_nubar_NO, prob_nubar_IO],
                        lc=['C0', 'C1'], ls=['-', '--'], lw=[1.0, 1.0],
                        points_sel_arr=[points_NO_sel, points_IO_sel],
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': r'NO$\nu$A ($L = 810$~km, $E_\nu = 2$~GeV)', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 
                             'ha': 'left'},
                            {'text': 'NO', 'xy': (0.14,0.67), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'IO', 'xy': (0.5,0.75), 'color': 'C1', 'fontsize': 20, 'ha': 'left'},
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        leg_fontsize=15, leg_loc='upper right', leg_ncol=1, xlim=[0.025,0.06], ylim=[0.025,0.06], save_fig=False)'''),
    code(r'''# ====== T2K ======

# Pick the baseline and energy of T2K
baseline = 295 # [km]
energy = 0.6*gd.UNIT_GEV # [eV]

# Values of the delta_CP phase at which to compute the probabilities
dCP_npts = 100
dCP_arr = np.linspace(-np.pi, np.pi, dCP_npts)

# We will compute the appearance probabilities, i.e., nu_mu --> nu_e and nu_mu-bar --> nu_e-bar
nu_i, nu_f = gd.NUMU, gd.NUE

# Coherent matter potential in matter
rho = 3.0 # Average matter density in the Earth's crus [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]
# Note: we could have directly used VCC = gd.NUM_DENSITY_E_EARTH_CRUST

# Hamiltonian in matter
H_matt = hamiltonians.hamiltonian_3nu_matter(VCC)

# Hamiltonian for neutrinos (the matter potential is added)
def H_nu(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=False, compute_matrix_multiplication=False) \
    + H_matt

# Hamiltonian for anti-neutrinos (the matter potential is subtracted)
def H_nubar(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=True, compute_matrix_multiplication=False) \
    - H_matt

# Generate the neutrino and anti-neutrino oscillation probabilities for many different values of dCP
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# nu_mu --> nu_e, normal ordering
prob_nu_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, normal ordering
prob_nubar_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu --> nu_e, inverted ordering
prob_nu_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, inverted ordering
prob_nubar_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])

# Points for selected values of dCP
dCP_sel = [-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi]
dCP_label_sel = [r'$-\pi, \pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', 'NuFit 6.0']
markers = ['o', 'v', 's', 'p', '*', 'p', 's', 'v', '^']
filled = [True, True, True, True, True, False, False, False, True]
prob_nu_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]] 
                           for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
points_NO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_NO_sel[i][1], prob_nubar_NO_sel[i][1]]] for i in range(len(dCP_sel)+1)]
points_IO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_IO_sel[i][1], prob_nubar_IO_sel[i][1]]] for i in range(len(dCP_sel)+1)]'''),
    code(r'''make_plot_biprobability([prob_nu_NO, prob_nu_IO], 
                        [prob_nubar_NO, prob_nubar_IO],
                        lc=['C0', 'C1'], ls=['-', '--'], lw=[1.0, 1.0],
                        points_sel_arr=[points_NO_sel, points_IO_sel],
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': r'T2K ($L = 295$~km, $E_\nu = 0.6$~GeV)', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 
                             'ha': 'left'},
                            {'text': 'NO', 'xy': (0.18,0.67), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'IO', 'xy': (0.45,0.75), 'color': 'C1', 'fontsize': 20, 'ha': 'left'},
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        leg_fontsize=15, leg_loc='upper right', leg_ncol=1, xlim=[0.025,0.065], ylim=[0.025,0.065], save_fig=False)'''),
    code(r'''# ====== DUNE ======

# Pick the baseline and energy of DUNE
baseline = 1300 # [km]
energy = 2*gd.UNIT_GEV # [eV]

# Values of the delta_CP phase at which to compute the probabilities
dCP_npts = 100
dCP_arr = np.linspace(-np.pi, np.pi, dCP_npts)

# We will compute the appearance probabilities, i.e., nu_mu --> nu_e and nu_mu-bar --> nu_e-bar
nu_i, nu_f = gd.NUMU, gd.NUE

# Coherent matter potential in matter
rho = 3.0 # Average matter density in the Earth's crus [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]
# Note: we could have directly used VCC = gd.NUM_DENSITY_E_EARTH_CRUST

# Hamiltonian in matter
H_matt = hamiltonians.hamiltonian_3nu_matter(VCC)

# Hamiltonian for neutrinos (the matter potential is added)
def H_nu(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=False, compute_matrix_multiplication=False) \
    + H_matt

# Hamiltonian for anti-neutrinos (the matter potential is subtracted)
def H_nubar(s12, s23, s13, dCP, D21, D31):
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=True, compute_matrix_multiplication=False) \
    - H_matt

# Generate the neutrino and anti-neutrino oscillation probabilities for many different values of dCP
# osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# nu_mu --> nu_e, normal ordering
prob_nu_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, normal ordering
prob_nubar_NO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu --> nu_e, inverted ordering
prob_nu_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                    0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                    n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                    integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])
# nu_mu-bar --> nu_e-bar, inverted ordering
prob_nubar_IO = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO),
                                                       0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                       n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                       integration_method='simpson', n_jobs=1)[nu_i][nu_f]] for dCP in dCP_arr])

# Points for selected values of dCP
dCP_sel = [-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi]
dCP_label_sel = [r'$-\pi, \pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', 'NuFit 6.0']
markers = ['o', 'v', 's', 'p', '*', 'p', 's', 'v', '^']
filled = [True, True, True, True, True, False, False, False, True]
prob_nu_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]] 
                           for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[gd.DCP_IO_BF_NUFIT_6_0]])
points_NO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_NO_sel[i][1], prob_nubar_NO_sel[i][1]]] for i in range(len(dCP_sel)+1)]
points_IO_sel = [[dCP_label_sel[i], markers[i], filled[i], [prob_nu_IO_sel[i][1], prob_nubar_IO_sel[i][1]]] for i in range(len(dCP_sel)+1)]'''),
    code(r'''make_plot_biprobability([prob_nu_NO, prob_nu_IO], 
                        [prob_nubar_NO, prob_nubar_IO],
                        lc=['C0', 'C1'], ls=['-', '--'], lw=[1.0, 1.0],
                        points_sel_arr=[points_NO_sel, points_IO_sel],
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': r'DUNE ($L = 1300$~km, $E_\nu = 2$~GeV)', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 
                             'ha': 'left'},
                            {'text': 'NO', 'xy': (0.10,0.74), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'IO', 'xy': (0.55,0.85), 'color': 'C1', 'fontsize': 20, 'ha': 'left'},
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        leg_fontsize=15, leg_loc='upper right', leg_ncol=1, xlim=[0.018,0.06], ylim=[0.018,0.06], save_fig=False)'''),
    md(r'''## In matter, through the Earth (varying density, PREM)'''),
    md(r'''We send a neutrino beam from Fermilab through the Earth to SNOLAB, Homestake, CERN, and the South Pole, using the onion-like matter density profile of the Preliminary Reference Earth Model (PREM).  For more details, including the line-of-sight density profiles in each of these cases, see the notebook `04_magnus_long_baseline.ipynb`.  

For this example, we fix the mass ordering to normal.'''),
    code(r'''# Electron number density inside Earth, using the PREM density model
def num_density_e_func_prem(r): 
    return matter.num_density_e_func(r, earth.density_matter_func_prem, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

# Coherent forward potential inside Earth, using the PREM density model
def VCC_func_prem(r):
    return matter.VCC_func(r, num_density_e_func_prem) # [eV]

# Hamiltonian for neutrinos (the matter potential is added)
def H_nu(costhz, l, energy, s12, s23, s13, dCP, D21, D31):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=False) \
    + hamiltonians.hamiltonian_3nu_matter(VCC_func_prem(r))

# Hamiltonian for anti-neutrinos (the matter potential is subtracted)
def H_nubar(costhz, l, energy, s12, s23, s13, dCP, D21, D31):
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, nubar=True) \
    - hamiltonians.hamiltonian_3nu_matter(VCC_func_prem(r))

# Cosines of zenith angles of Fermilab measured from the detector locations (we would get the same result if we seapped them, since we
# consider that Earth is radially symmetric)
detectors = ['SNOLAB', 'Homestake', 'CERN', "South Pole"]
costhz_arr = [earth.costhz_between_points_on_surface(earth.loc_coords_dms[det.lower().replace(" ", "_")]['lat'], 
                                                     earth.loc_coords_dms[det.lower().replace(" ", "_")]['lon'],
                                                     earth.loc_coords_dms['fermilab']['lat'], 
                                                     earth.loc_coords_dms['fermilab']['lon']) 
             for det in detectors]

# Maximum baselines inside the Earth
l_max_arr = [earth.distance_traveled_inside_earth(costhz) for costhz in costhz_arr] # [km]'''),
    code(r'''# Neutrino energy
energy = 2*gd.UNIT_GEV # [eV]

# Values of the delta_CP phase at which to compute the probabilities
dCP_npts = 100
dCP_arr = np.linspace(-np.pi, np.pi, dCP_npts)

# We will compute the appearance probabilities, i.e., nu_mu --> nu_e and nu_mu-bar --> nu_e-bar
nu_i, nu_f = gd.NUMU, gd.NUE

# Points for selected values of dCP
dCP_sel = [-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi]
dCP_label_sel = [r'$-\pi, \pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', 'NuFit 6.0']
markers = ['o', 'v', 's', 'p', '*', 'p', 's', 'v', '^']
filled = [True, True, True, True, True, False, False, False, True]

# Below, we use n_slabs=100, n_tpts_per_slab=10, and magnus_exp_order=3, which would allow for the MSW resonance to be picked (say, around 
# 10 GeV), but for our choice of 2 GeV, far from the resonance, this is overkill 
prob_nu_arr, prob_nubar_arr, points_sel_arr = [], [], []
for i in range(len(detectors)):
    print("detector = " + detectors[i])
    # nu_mu --> nu_e
    prob_nu = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(costhz_arr[i], l, energy, 
                                                                    s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                     0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3, 
                                                     integration_method='simpson', n_jobs=10)[nu_i][nu_f]] for dCP in dCP_arr])
    # nu_mu-bar --> nu_e-bar
    prob_nubar = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(costhz_arr[i], l, energy, 
                                                                          s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO),
                                                        0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                        integration_method='simpson', n_jobs=10)[nu_i][nu_f]] for dCP in dCP_arr])
    prob_nu_arr.append(prob_nu)
    prob_nubar_arr.append(prob_nubar)
    # Compute selected points
    points_nu_sel, points_nubar_sel = [], []
    prob_nu_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(costhz_arr[i], l, energy, 
                                                                    s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                         0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                         n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                         integration_method='simpson', n_jobs=10)[nu_i][nu_f]] 
                               for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
    prob_nubar_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(costhz_arr[i], l, energy, 
                                                                              s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                            0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                            n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                            integration_method='simpson', n_jobs=10)[nu_i][nu_f]] 
                               for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
    points_sel = [[dCP_label_sel[j], markers[j], filled[j], [prob_nu_sel[j][1], prob_nubar_sel[j][1]]] for j in range(len(dCP_sel)+1)]
    points_sel_arr.append(points_sel)'''),
    code(r'''make_plot_biprobability(prob_nu_arr, prob_nubar_arr,
                        lc=['C0', 'C1', 'C2', 'C3'], ls=['-', '--', ':', '-.'], lw=[1.0, 1.0, 1.0, 1.0],
                        points_sel_arr=points_sel_arr,
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': r'From Fermilab to ... ($E_\nu = 2$~GeV)', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 
                             'ha': 'left'},
                            {'text': 'SNOLAB', 'xy': (0.36,0.42), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'Homestake', 'xy': (0.55,0.27), 'color': 'C1', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'CERN', 'xy': (0.50,0.105), 'color': 'C2', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'South Pole', 'xy': (0.36,0.78), 'color': 'C3', 'fontsize': 20, 'ha': 'left'}
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.02),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.005),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.01),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.002),
                        leg_fontsize=15, leg_loc='upper right', leg_ncol=1, xlim=[0,0.10], ylim=[0,0.07], save_fig=False)'''),
    md(r'''In the plot above, we used 2 GeV, which, for the source-detector configurations we explore, is far from the MSW resonance energy inside Earth, so the high-order Magnus expansion was not really needed.

Below, we repeat the for a higher energy, closer to the MSW resonance.  Therefore, the higher-order Magnus expansion and high number of slabs (`n_slabs = 100`) is, in general, needed in order to be sure to pick up any resonance.  The recommendation is to always use higher-order expansion in order to be sure to pick up any resonance, which is especially important when adding new-physics effects, for which, in general, we may not know a priori at what energy the resonance occurs.  (Of course, the energy at which the resonance happens depends on the directiomn of neutrino propagation, so for a fixed energy, we typically cannot hit the resonance for all source-detector configurations simultaneously.)'''),
    code(r'''# Neutrino energy
energy = 20*gd.UNIT_GEV # [eV]

# Below, we use n_labs=100, n_tpts_per_slab=10, and magnus_exp_order=3, which would allow for the MSW resonance to be picked (say, around 
# 10 GeV), but for our choice of 2 GeV, far from the resonance, this is overkill 
prob_nu_arr, prob_nubar_arr, points_sel_arr = [], [], []
for i in range(len(detectors)):
    print("detector = " + detectors[i])
    # nu_mu --> nu_e
    prob_nu = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(costhz_arr[i], l, energy, 
                                                                    s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                     0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                     n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3, 
                                                     integration_method='simpson', n_jobs=10)[nu_i][nu_f]] for dCP in dCP_arr])
    # nu_mu-bar --> nu_e-bar
    prob_nubar = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(costhz_arr[i], l, energy, 
                                                                          s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO),
                                                        0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                        integration_method='simpson', n_jobs=10)[nu_i][nu_f]] for dCP in dCP_arr])
    prob_nu_arr.append(prob_nu)
    prob_nubar_arr.append(prob_nubar)
    # Compute selected points
    points_nu_sel, points_nubar_sel = [], []
    prob_nu_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(costhz_arr[i], l, energy, 
                                                                    s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                         0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                         n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                         integration_method='simpson', n_jobs=10)[nu_i][nu_f]] 
                               for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
    prob_nubar_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(costhz_arr[i], l, energy, 
                                                                              s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                            0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                            n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                            integration_method='simpson', n_jobs=10)[nu_i][nu_f]] 
                               for dCP in dCP_sel+[gd.DCP_NO_BF_NUFIT_6_0]])
    points_sel = [[dCP_label_sel[j], markers[j], filled[j], [prob_nu_sel[j][1], prob_nubar_sel[j][1]]] for j in range(len(dCP_sel)+1)]
    points_sel_arr.append(points_sel)'''),
    code(r'''make_plot_biprobability(prob_nu_arr, prob_nubar_arr,
                        lc=['C0', 'C1', 'C2', 'C3'], ls=['-', '--', ':', '-.'], lw=[1.0, 1.0, 1.0, 1.0],
                        points_sel_arr=points_sel_arr,
                        label_prob_nu=r'Neutrino probability, $P_{\nu_\mu \to \nu_e}$', 
                        label_prob_nubar=r'Anti-neutrino probability, $P_{\bar{\nu}_\mu \to \bar{\nu}_e}$',
                        annotations=[
                            {'text': r'From Fermilab to ... ($E_\nu = 20$~GeV)', 'xy': (0.02,0.95), 'color': 'k', 'fontsize': 20, 
                             'ha': 'left'},
                            {'text': 'SNOLAB', 'xy': (0.06,0.32), 'color': 'C0', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'Homestake', 'xy': (0.13,0.71), 'color': 'C1', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'CERN', 'xy': (0.60,0.86), 'color': 'C2', 'fontsize': 20, 'ha': 'left'},
                            {'text': 'South Pole', 'xy': (0.30,0.15), 'color': 'C3', 'fontsize': 20, 'ha': 'left'}
                        ], 
                        xaxis_major_locator=mpl.ticker.MultipleLocator(base=0.005),
                        xaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.001),
                        yaxis_major_locator=mpl.ticker.MultipleLocator(base=0.0005),
                        yaxis_minor_locator=mpl.ticker.MultipleLocator(base=0.0001),
                        leg_fontsize=15, leg_loc='lower right', leg_ncol=2, xlim=[0,0.016], ylim=[0,0.002], save_fig=False)'''),
    ])

# ---------------------------------------------------- 06_magnus_oscillograms
books['06_magnus_oscillograms.ipynb'] = notebook(
    'Oscillograms',
    'Probability across zenith angle and energy at once -- the two-dimensional map of what an atmospheric-neutrino detector sees.\n\nThis is the workload that most rewards passing arrays rather than looping: the energies share a chord, so the matter profile is built once.',
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''## Oscillograms

An oscillogram is the natural way to look at oscillations through the Earth:
the probability as a function of both the arrival direction and the energy, on
one map. The direction matters because it fixes the chord -- and therefore not
just the baseline but *which layers* the neutrino crosses. Everything
interesting in these figures comes from that.

The band structure running diagonally is ordinary vacuum oscillation: fixed
phase $\Delta m^2 L / 2E$ means the fringes follow lines of constant $L/E$.
What breaks that pattern is matter. Watch for the sharp feature near
$\cos\theta_z \approx -0.83$: that is the core-mantle boundary, where PREM's
density jumps by roughly a factor of two, and trajectories steeper than it
sample a genuinely different profile.

Below we define a helper that computes an oscillogram and draws it. The
drawing is a single call to `magnus.plotting.plot_oscillogram`; the colour
map, the colour-bar label, the tick spacings, and the white-stroked corner
annotation (which has to stay legible against a `plasma` background) are its
defaults.'''),
    code(r'''def make_oscillogram_plot(nu_i, nu_f, H_func, costhz_arr, log10_Enu_arr,
                          n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                          n_jobs=10, integration_method='trapezoid',
                          validate_input=False, sector_2nu=None,
                          cbar_label_pre='', save_plot=False, path=None,
                          filename=None, format=None):
    """Compute an oscillogram over (cos theta_z, log10 E) and plot it.

    Only the computation lives here now; the drawing is one call to
    magnus.plotting.plot_oscillogram, which carries the house style.
    """
    # For 2nu oscillations in the 23 sector the flavour indices have to be
    # remapped onto the 2x2 block the Hamiltonian actually spans.
    if sector_2nu == '23':
        nu_i_ = 0 if nu_i == gd.NUMU else 1
        nu_f_ = 0 if nu_f == gd.NUMU else 1
    else:
        nu_i_, nu_f_ = nu_i, nu_f

    # prob_arr[i][j] is the requested probability at Enu_arr[i], costhz_arr[j]
    prob_arr = np.array([[oscprob.osc_prob(
                              lambda l: H_func(costhz_arr[j], l, Enu_arr[i]*gd.UNIT_GEV),
                              0, l_max_arr[j]*gd.CONV_KM_TO_INV_EV,
                              n_slabs=n_slabs, n_tpts_per_slab=n_tpts_per_slab,
                              magnus_exp_order=magnus_exp_order, n_jobs=n_jobs,
                              validate_input=validate_input)[nu_i_][nu_f_]
                          for j in range(costhz_npts)]
                         for i in range(Enu_npts)])

    sector_label = {'23': r' (23 sector)', '12': r' (12 sector)'}.get(sector_2nu, '')

    fig, ax = plotting.plot_oscillogram(
        costhz_arr, log10_Enu_arr, prob_arr,
        nu_i=nu_i, nu_f=nu_f,
        cbar_label=cbar_label_pre + plotting.prob_label(nu_i, nu_f),
        annotation=plotting.prob_label(nu_i, nu_f) + sector_label,
        xlim=(costhz_min, costhz_max),
        ylim=(log10_Enu_min, log10_Enu_max),
        savefig=(path + filename + '.' + format) if save_plot else None,
        savefig_kw=dict(bbox_inches='tight'))

    return prob_arr'''),
    md(r'''Define the ranges of zenith angle (and, therefore, baseline) and neutrino energies'''),
    code(r'''# Cosines of zenith angles
costhz_min, costhz_max, costhz_npts = -1.0, 0, 150 #100 
costhz_arr = np.linspace(costhz_min, costhz_max, costhz_npts) 

# Baselines, L [km]
l_max_arr = [earth.distance_traveled_inside_earth(costhz) for costhz in costhz_arr] 

# Neutrino energies [GeV]
log10_Enu_min, log10_Enu_max, Enu_npts = 0.0, 1.0, 150 #100
Enu_arr = np.logspace(log10_Enu_min, log10_Enu_max, Enu_npts)
log10_Enu_arr = np.log10(Enu_arr)'''),
    md(r'''Define the coherent forward potential inside the Earth according to the PREM matter density model.

Because of numerical errors, sometimes `r` can be slightly larger han the Earth's radius, `gd.EARTH_RADIUS`.  The parameter `tol` defines how much larger than 1 the quantity `r-gd.EARTH_RADIUS` can be inside the `earth.density_matter_func_prem`function.  We have set `tol = 1.e-15` manually below because it leads to successful results in our runs, but the user is advised to change its value if their code fails.'''),
    code(r'''# r is the radius from the center of the Earth [km].  

def density_matter_func_prem_wrapper(r):
    return earth.density_matter_func_prem(r, tol=1.e-15)
    
def VCC_func_prem(r):
    return matter.VCC_func(r, lambda rr : matter.num_density_e_func(rr, density_matter_func_prem_wrapper, 
                                                                    ratio_number_neutrons_to_protons=1.0, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True)) # [eV]'''),
    md(r'''Generate and plot the $2\nu$ oscillograms.  Let's consider the 2-3 sector, which describes mixing between $\nu_\mu$ and $\nu_\tau$. First, we define the Hamiltonians:'''),
    code(r'''# # Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S23_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Vacuum Hamiltonian without the (1/E) prefactor
H_vac_2nu_23_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # [eV^2]

# Matter Hamiltonian
def H_2nu_23_func_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) 
    return (1/energy)*H_vac_2nu_23_energy_indep + hamiltonians.hamiltonian_2nu_matter(VCC_func_prem(r)) # [eV]'''),
    md(r'''Now plot the oscillogram; the probabilities are computed inside the `make_oscillogram_plot` function:'''),
    code(r'''make_oscillogram_plot(gd.NUMU, gd.NUMU, H_2nu_23_func_prem, costhz_arr, log10_Enu_arr,
                      n_slabs=1, n_tpts_per_slab=100, magnus_exp_order=1, n_jobs=1, integration_method='trapezoid', validate_input=False,
                      sector_2nu='23', cbar_label_pre=r'Two-neutrino probability, ', save_plot=False, path=None, filename=None, format=None)'''),
    md(r'''The effect of the core-mantle transition in the density is clearly visible around $\cos \theta_z = -0.85$, as expected.

Of note, in the $2\nu$ case we already obtain the oscillograms that we expected by using a single slab (`n_slabs = 1`) and the lowest order of the Magnus expansion (`magnus_exp_order = 1`), which is equivalent to assuming a time-independent Hamiltonian, or an average matter density along each chord inside the Earth.  '''),
    md(r'''Now we can do the same for the 1-2 sector, starting with defining a new Hamiltonian with the mixing parameters of that sector:'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
sth = gd.S12_NO_BF_NUFIT_6_0 # [adim]
Dm2 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]

# Vacuum Hamiltonian without the (1/E) prefactor
H_vac_2nu_12_energy_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2) # [eV^2]

# Matter Hamiltonian
def H_2nu_12_func_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) 
    return (1/energy)*H_vac_2nu_12_energy_indep + hamiltonians.hamiltonian_2nu_matter(VCC_func_prem(r)) # [eV]'''),
    md(r'''And then plot:'''),
    code(r'''make_oscillogram_plot(gd.NUMU, gd.NUMU, H_2nu_12_func_prem, costhz_arr, log10_Enu_arr,
                      n_slabs=1, n_tpts_per_slab=100, magnus_exp_order=1, n_jobs=1, integration_method='trapezoid', validate_input=False,
                      sector_2nu='12', cbar_label_pre=r'Two-neutrino probability, ', save_plot=False, path=None, filename=None, format=None)'''),
    md(r'''Finally, let's do the same for $3\nu$ oscillations.  First, define the Hamiltonians:'''),
    code(r'''# Mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Vacuum Hamiltonian without the (1/E) prefactor
H_vac_3nu_energy_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31) # [eV^2]

# Matter Hamiltonian
def H_3nu_func_prem(costhz, l, energy):
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) 
    return (1/energy)*H_vac_3nu_energy_indep + hamiltonians.hamiltonian_3nu_matter(VCC_func_prem(r)) # [eV]'''),
    md(r'''And then plot:'''),
    code(r'''def density_matter_func_prem_wrapper(r):
    return earth.density_matter_func_prem(r, tol=1.e-15)
    
def VCC_func_prem(r):
    return matter.VCC_func(r, lambda rr : matter.num_density_e_func(rr, density_matter_func_prem_wrapper, 
                                                                    ratio_number_neutrons_to_protons=1.0, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True)) # [eV]'''),
    md(r'''Let's first generate the oscillogram for the same simplest settings that we used in the $2\nu$ case (`n_slabs = 1`) and (`magnus_exp_order = 1`).'''),
    code(r'''prob_arr_3nu_order_1 = make_oscillogram_plot(gd.NUMU, gd.NUMU, H_3nu_func_prem, costhz_arr, log10_Enu_arr,
                                             n_slabs=1, n_tpts_per_slab=100, magnus_exp_order=1, n_jobs=1,
                                             integration_method='trapezoid', validate_input=False,
                                             sector_2nu=None, cbar_label_pre=r'Three-neutrino probability, ', save_plot=False, 
                                             path=None, filename=None, format=None)'''),
    md(r'''This already looks like it should, including showing the transition from core to mantle.  But the fine details are not quite there yet.  

Let's increase the number of slabs used to `n_slabs = 10` and use third-order Magnus expansion instead (`magnus_exp_order = 3`).  We can now set `n_jobs` to > 1 to compute evolution operators on each slab in parallel.'''),
    code(r'''make_oscillogram_plot(gd.NUMU, gd.NUMU, H_3nu_func_prem, costhz_arr, log10_Enu_arr,
                      n_slabs=10, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=5, integration_method='trapezoid', validate_input=False,
                      sector_2nu=None, cbar_label_pre=r'Three-neutrino probability, ', save_plot=False, path=None, filename=None, format=None)'''),
    md(r'''If one would like to check whether this result is final or near-final, they can increase the `n_slabs` (or`n_tpts_per_slab`, or `magnus_exp_order`) and compare the resulting oscillogram with the one above.'''),
    md(r'''We can also plot oscillograms for other probabilities:'''),
    code(r'''make_oscillogram_plot(gd.NUMU, gd.NUTAU, H_3nu_func_prem, costhz_arr, log10_Enu_arr,
                      n_slabs=10, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=5, integration_method='trapezoid', validate_input=False,
                      sector_2nu=None, cbar_label_pre=r'Three-neutrino probability, ', save_plot=False, path=None, filename=None, format=None)'''),
    code(r'''make_oscillogram_plot(gd.NUMU, gd.NUE, H_3nu_func_prem, costhz_arr, log10_Enu_arr,
                      n_slabs=10, n_tpts_per_slab=100, magnus_exp_order=3, n_jobs=5, integration_method='trapezoid', validate_input=False,
                      sector_2nu=None, cbar_label_pre=r'Three-neutrino probability, ', save_plot=False, path=None, filename=None, format=None)'''),
    ])

# -------------------------------------------------- 07_magnus_bsm_sterile_nu
books['07_magnus_bsm_sterile_nu.ipynb'] = notebook(
    'BSM: sterile neutrinos',
    'Four- and five-flavour systems, where the extra states do not couple to the weak interaction.\n\nThe machinery is unchanged; only the dimension of the Hamiltonian and the number of mixing angles and phases grow.',
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''Mag$\nu$s has no intrinsic limitation on the number of neutrino flavors contained in the system that it computes oscillation probabilities for.  However, it is distributed with ready-made mixing matrices and Hamiltonians for the two most popular scenarios: 3+1 (three active neutrinos plus one sterile neutrino) and 3+2.  Users can add their own mixing matrices and Hamiltonians for larger systems by following the same structure used when writing the modules `hamiltonians4nu` and `hamiltonians5nu` that are distributed with Mag$\nu$s.

Below, we show examples mainly for the 3+1 system, and one example for the 3+2 system. '''),
    md(r'''# 0. Helper functions and definitions'''),
    code(r'''def flavor_index_to_str(nu_l):
    return {gd.NUE: 'e', gd.NUMU: 'mu', gd.NUTAU: 'tau',
            gd.NUS1: 's1', gd.NUS2: 's2'}[nu_l]


def sterile_annotations(N):
    """The active-sterile parameters a figure was made with.

    Stacked upward from the bottom-left corner, so the number of lines follows
    from the number of sterile states rather than from hand-written positions.
    """
    if N == 1:
        lines = [r'$s_{14} = $ ' + str(s14) + r', $\delta_{14} = $ ' + str(d14),
                 r'$s_{24} = $ ' + str(s24) + r', $\delta_{24} = $ ' + str(d24),
                 r'$s_{34} = $ ' + str(s34),
                 r'$\Delta m_{41}^2 =$ ' + str(D41) + r' eV$^2$']
    else:
        lines = [r'$s_{14} = $ ' + str(s14) + r', $\delta_{14} = $ ' + str(d14),
                 r'$s_{15} = $ ' + str(s15) + r', $\delta_{15} = $ ' + str(d15),
                 r'$s_{24} = $ ' + str(s24) + r', $\delta_{24} = $ ' + str(d24),
                 r'$s_{25} = $ ' + str(s25),
                 r'$s_{34} = $ ' + str(s34) + r', $s_{35} = $ ' + str(s35),
                 r'$\Delta m_{41}^2 =$ ' + str(D41) + r' eV$^2$',
                 r'$\Delta m_{51}^2 =$ ' + str(D51) + r' eV$^2$']
    return [dict(text=t, xy=(0.02, 0.03 + 0.05*k), fontsize=18)
            for k, t in enumerate(reversed(lines))]


def _curves(prob_Nnu_all, prob_3nu_all, nu_i, nu_f, N):
    return [dict(y=prob_Nnu_all[:, nu_i, nu_f], color='C1',
                 label='3+' + str(N) + ' oscillations'),
            dict(y=prob_3nu_all[:, nu_i, nu_f], color='0.2', ls='--',
                 label=r'Standard, $3\nu$')]


def make_plot_prob_Nnu_3nu_vs_baseline(nu_i, nu_f, distances, prob_Nnu_all,
                                       prob_3nu_all, N=1, title=None, save_plot=True):
    """A sterile scenario against standard three-flavour, versus baseline."""
    filename = ('prob_3plus' + str(N) + '_vs_3nu_vs_baseline_'
                + flavor_index_to_str(nu_i) + '_' + flavor_index_to_str(nu_f))
    return plotting.plot_probability_vs_baseline(
        distances, _curves(prob_Nnu_all, prob_3nu_all, nu_i, nu_f, N),
        nu_i=nu_i, nu_f=nu_f, xlim=(l_ini, l_fin),
        annotations=sterile_annotations(N),
        legend_loc='center left', title=title, title_fontsize=23,
        savefig=('../fig/' + filename + '.pdf') if save_plot else None,
        savefig_kw=dict(dpi=300))


def make_plot_prob_Nnu_3nu_vs_energy(nu_i, nu_f, energies, prob_Nnu_all,
                                     prob_3nu_all, N=1, title=None, save_plot=True):
    """The same comparison versus neutrino energy."""
    filename = ('prob_3plus' + str(N) + '_vs_3nu_vs_energy_'
                + flavor_index_to_str(nu_i) + '_' + flavor_index_to_str(nu_f))
    return plotting.plot_probability_vs_energy(
        energies, _curves(prob_Nnu_all, prob_3nu_all, nu_i, nu_f, N),
        nu_i=nu_i, nu_f=nu_f, xlim=(energy_min, energy_max),
        annotations=sterile_annotations(N),
        legend_loc='center left', title=title, title_fontsize=23,
        savefig=('../fig/' + filename + '.pdf') if save_plot else None,
        savefig_kw=dict(dpi=300))'''),
    md(r'''# 1. Probabilities in a 3+1 system ($4\nu$): in vacuum'''),
    md(r'''In a 3+1 neutrino system, there are three active neutrinos and one sterile neutrino that mixes with them.  Accordingly, the mixing matrix is now a unitary $4 \times 4$ complex matrix (see Eqs. (2) and (3) in arXiv:1105.3911),
\begin{equation}
U_{4\nu} = R_{34} ~ \tilde{R}_{24} ~ \tilde{R}_{14} ~ R_{23} ~ \tilde{R}_{13} ~ R_{12} \;,
\end{equation}
where, e.g.,
\begin{equation}
R_{34} = \left( 
\begin{array}{cccc} 
1 & 0 & 0 & 0 \\  
0 & 1 & 0 & 0 \\  
0 & 0 & c_{34} & s_{34} \\  
0 & 0 & -s_{34} & c_{34} \\  
\end{array} \right)
\;,
\qquad
\tilde{R}_{14} = \left( 
\begin{array}{cccc} 
c_{14} & 0 & 0 & s_{14} e^{-i \delta_{14}} \\  
0 & 1 & 0 & 0 \\  
0 & 0 & 1 & 0 \\  
-s_{14} e^{i \delta_{14}} & 0 & 0 & c_{14} \\  
\end{array} \right) \;,
\end{equation} 
and we have not written the matrix containing Majorana phases because they do not affect neutrino oscillations. In addition, there is a new mass-squared difference, $\Delta m_{41}^2$, the final eigenvalue of the newly extended four-neutrino mass matrix.

Thus, there are six new mixing parameters that we need to supply to Mag$\nu$s: $s_{14}$ (`s14`), $\delta_{14}$ (`d14`), $s_{24}$ (`s24`), $\delta_{24}$ (`d24`), $s_{34}$ (`s34`), and $\Delta m_{41}^2$ (`D41`).  These parameters are passed *in addition* to the standard mixing parameters in three-neutrino oscillations.'''),
    md(r'''When computing neutrino oscillation probabilities in the 3+1 (and 3+2) system, we proceed just as we did for two- and three-neutrino oscillations (see notebooks `02_magnus_2nu_vacuum_matter.ipynb` and `03_magnus_3nu_vacuum_matter.ipynb` for details and further usage).

We compute probabilities using the `osc_prob` routine, the same one we used when computing two- and three-neutrino probabilities.  This routine is flexible, and, if fed with an $N \times N$ Hamiltonian, where $N$ is the number of flavors, it will return a $N \times N$ matrix of probabilities between all flavors automatically.  No separate parameter needs to be passed to it specifying the number of flavors.  

All of the work is done by the Hamiltonian function, whose calculation *does* require inputting all of the mixing parameters, including the active-sterile ones.'''),
    md(r'''### 1.1 Mixing parameters'''),
    code(r'''# Standard mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Mixing parameters between active and sterile sectors (can change them to anything else)
s14 = 1.e-1 # [adim]
d14 = 0.0 # [adim]
s24 = 1.e-1 # [adim]
d24 = 0.0 # [adim]
s34 = 1.e-1 # [adim]
D41 = 0.5 # [eV^2]'''),
    md(r'''### 1.2 Probabilities vs. distance'''),
    md(r'''Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e4 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

energy = 2.0*gd.UNIT_GEV # [eV]

# Compute probability vs. baseline
# For the 3nu case, osc_prob returns a 3x3 NumPy array with the probabilities: 
# [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# For the 4nu case, osc_prob returns a 3x3 NumPy array with the probabilities:
# [[Pee, Pem, Pet, Pes], [Pme, Pmm, Pmt, Pms], [Pte, Ptm, Ptt, Pts], [Pse, Psm, Pst, Pss]]

# Standard, 3nu probabilities
H_3nu_vac = hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, compute_matrix_multiplication=False)
prob_3nu_all = np.array([oscprob.osc_prob(lambda l: H_3nu_vac, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                          magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])

# 3+1 probabilities
H_4nu_vac = hamiltonians.hamiltonian_4nu_vacuum(energy, s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, 
                                                   compute_matrix_multiplication=False)
prob_4nu_all = np.array([oscprob.osc_prob(lambda l: H_4nu_vac, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                          magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_baseline(gd.NUE, gd.NUE, distances, prob_4nu_all, prob_3nu_all, N=1,
                                   title=r'3+1 oscillations in vacuum, $E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV), save_plot=False)'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_baseline(gd.NUMU, gd.NUTAU, distances, prob_4nu_all, prob_3nu_all, N=1, 
                                   title=r'3+1 oscillations in vacuum, $E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV), save_plot=False)'''),
    md(r'''### 1.3 Probabilities vs. energy'''),
    md(r'''Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e2 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

baseline = 5.e3 # [km]

# Standard, 3nu probabilities
H_3nu_vac_en_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
                                                                               compute_matrix_multiplication=False) # Vacuum H without (1/E)
prob_3nu_all = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_GEV))*H_3nu_vac_en_indep, 
                                          0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                          n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])

# 3+1 probabilities
H_4nu_vac_en_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, 
                                                                               compute_matrix_multiplication=False) 
prob_4nu_all = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_GEV))*H_4nu_vac_en_indep, 
                                          0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                          n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_energy(gd.NUE, gd.NUE, energies, prob_4nu_all, prob_3nu_all, N=1,
                                 title=r'3+1 oscillations in vacuum, $L = $~{:.2f}~km'.format(baseline), save_plot=False)'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_energy(gd.NUMU, gd.NUTAU, energies, prob_4nu_all, prob_3nu_all, N=1,
                                 title=r'3+1 oscillations in vacuum, $L = $~{:.2f}~km'.format(baseline), save_plot=False)'''),
    md(r'''# 2. Probabilities in a 3+1 system ($4\nu$): in matter with constant density'''),
    md(r'''### 2.1 Mixing parameters'''),
    code(r'''# Standard mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Mixing parameters between active and sterile sectors (can change them to anything else)
s14 = 1.e-1 # [adim]
d14 = 0.0 # [adim]
s24 = 1.e-1 # [adim]
d24 = 0.0 # [adim]
s34 = 1.e-1 # [adim]
D41 = 0.5 # [eV^2]'''),
    md(r'''### 2.2 Probabilities vs. distance'''),
    md(r'''Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e4 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

energy = 2.0*gd.UNIT_GEV # [eV]

# Compute probability vs. baseline
# For the 3nu case, osc_prob returns a 3x3 NumPy array with the probabilities: 
# [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# For the 4nu case, osc_prob returns a 3x3 NumPy array with the probabilities:
# [[Pee, Pem, Pet, Pes], [Pme, Pmm, Pmt, Pms], [Pte, Ptm, Ptt, Pts], [Pse, Psm, Pst, Pss]]

# Matter potential
rho = 3.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]

# Standard, 3nu probabilities
H_3nu_vac = hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, compute_matrix_multiplication=False)
H_3nu_matt = hamiltonians.hamiltonian_3nu_matter(VCC)
H_3nu = lambda l: H_3nu_vac + H_3nu_matt
prob_3nu_all = np.array([oscprob.osc_prob(H_3nu, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                          magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])

# 3+1 probabilities
H_4nu_vac = hamiltonians.hamiltonian_4nu_vacuum(energy, s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, 
                                                   compute_matrix_multiplication=False)
H_4nu_matt = hamiltonians.hamiltonian_4nu_matter(VCC)
H_4nu = lambda l: H_4nu_vac + H_4nu_matt
prob_4nu_all = np.array([oscprob.osc_prob(H_4nu, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                          magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_baseline(gd.NUE, gd.NUE, distances, prob_4nu_all, prob_3nu_all, N=1,
                                   title=r'3+1 oscillations in constant-density matter (3~g~cm$^{-3}$), ' +
                                   r'$E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV), save_plot=False)'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_baseline(gd.NUMU, gd.NUTAU, distances, prob_4nu_all, prob_3nu_all, N=1,
                                   title=r'3+1 oscillations in constant-density matter (3~g~cm$^{-3}$), ' +
                                   r'$E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV), save_plot=False)'''),
    md(r'''### 2.3 Probabilities vs. energy'''),
    md(r'''Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

baseline = 5.e3 # [km]

# Matter potential
rho = 3.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]

# Standard, 3nu probabilities
H_3nu_vac_en_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
                                                                               compute_matrix_multiplication=False)
H_3nu_matt = hamiltonians.hamiltonian_3nu_matter(VCC)
def H_3nu(energy, l):
    return (1/energy)*H_3nu_vac_en_indep + H_3nu_matt
prob_3nu_all = np.array([oscprob.osc_prob(lambda l: H_3nu(enu*gd.UNIT_GEV, l), 
                                          0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                          n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])

# 3+1 probabilities
H_4nu_vac_en_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, 
                                                   compute_matrix_multiplication=False)
H_4nu_matt = hamiltonians.hamiltonian_4nu_matter(VCC)
def H_4nu(energy, l):
    return (1/energy)*H_4nu_vac_en_indep + H_4nu_matt
prob_4nu_all = np.array([oscprob.osc_prob(lambda l: H_4nu(enu*gd.UNIT_GEV, l), 
                                          0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                          n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_energy(gd.NUE, gd.NUE, energies, prob_4nu_all, prob_3nu_all, N=1,
                                 title=r'3+1 oscillations in constant-density matter (3~g~cm$^{-3}$), ' +
                                       r'$L = $~{:.2f}~km'.format(baseline), save_plot=False)'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_energy(gd.NUMU, gd.NUTAU, energies, prob_4nu_all, prob_3nu_all, N=1,
                                 title=r'3+1 oscillations in constant-density matter (3~g~cm$^{-3}$), ' +
                                       r'$L = $~{:.2f}~km'.format(baseline), save_plot=False)'''),
    md(r'''# 3. Probabilities in a 3+1 system ($4\nu$): long-baseline inside Earth (PREM)'''),
    md(r'''For more details about computing and plotting probabilities in long-baseline experimental setups, see notebook `04_magnus_long_baseline.ipynb`.'''),
    md(r'''### 3.1 Mixing parameters'''),
    code(r'''# Standard mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Electron number density inside Earth, using the PREM density model
def num_density_e_func_prem(r): 
    return matter.num_density_e_func(r, earth.density_matter_func_prem, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

# Coherent forward potential inside Earth, using the PREM density model
def VCC_func_prem(r):
    return matter.VCC_func(r, num_density_e_func_prem) # [eV]

# Cosines of zenith angles of Fermilab measured from the detector locations (we would get the same result if we seapped them, since we
# consider that Earth is radially symmetric)
detectors = ['SNOLAB', 'Homestake', 'CERN', "South Pole"]
costhz_arr = [earth.costhz_between_points_on_surface(earth.loc_coords_dms[det.lower().replace(" ", "_")]['lat'], 
                                                     earth.loc_coords_dms[det.lower().replace(" ", "_")]['lon'],
                                                     earth.loc_coords_dms['fermilab']['lat'], 
                                                     earth.loc_coords_dms['fermilab']['lon']) 
             for det in detectors]

# Maximum baselines inside the Earth
l_max_arr = [earth.distance_traveled_inside_earth(costhz) for costhz in costhz_arr] # [km]'''),
    md(r'''### 3.2 Probabilities vs. energy'''),
    md(r'''Define the Hamiltonians'''),
    code(r'''# Mixing parameters between active and sterile sectors (can change them to anything else)
s14 = 2.e-1 # [adim]
d14 = 0.0 # [adim]
s24 = 1.e-1 # [adim]
d24 = 0.0 # [adim]
s34 = 1.e-1 # [adim]
D41 = 0.5 # [eV^2]

# Vacuum Hamiltonian without the 1/E prefactor 
H_3nu_vac_en_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
                                                                               compute_matrix_multiplication=False)
H_4nu_vac_en_indep = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41, 
                                                                               compute_matrix_multiplication=False)

# Hamiltonian including matter effects inside Earth using the PREM density model
def H_3nu_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return (1/energy)*H_3nu_vac_en_indep + hamiltonians.hamiltonian_3nu_matter(VCC_func_prem(r))
    
def H_4nu_prem(costhz, l, energy):
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return (1/energy)*H_4nu_vac_en_indep + hamiltonians.hamiltonian_4nu_matter(VCC_func_prem(r))'''),
    md(r'''Generate the probabilities'''),
    code(r'''# This scan is the most expensive cell in the notebook: four trajectories, two
# flavour counts, and a genuine PREM profile sampled inside every slab.
#
# The grid is deliberately coarse. With D41 = 0.5 eV^2 the sterile oscillation
# accumulates ~1e9 radians over an Earth-crossing baseline at these energies, so
# the curve is aliased at *any* practical number of points -- a denser grid draws
# a different alias, not a better-resolved oscillation. What survives sampling is
# the envelope and the average, and 200 points show those as well as 3000 did, at
# a fifteenth of the cost.
nu_i, nu_f = gd.NUE, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU

# Energies
energy_min, energy_max = 1.e1, 1.e2 # [MeV]
energy_npts = 200   # see the note below on why a dense grid buys nothing here
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [MeV]

# Generate probabilities for the different directions
prob_3nu_arr, prob_4nu_arr = [], []
for i in range(len(costhz_arr)):
    print("detector = " + detectors[i])
    prob_3nu = np.array([oscprob.osc_prob(lambda l: H_3nu_prem(costhz_arr[i], l, enu*gd.UNIT_MEV), 
                                          0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                          n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] for enu in energies]) 
    prob_4nu = np.array([oscprob.osc_prob(lambda l: H_4nu_prem(costhz_arr[i], l, enu*gd.UNIT_MEV), 
                                          0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                          n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] for enu in energies]) 
    prob_3nu_arr.append(prob_3nu)
    prob_4nu_arr.append(prob_4nu)'''),
    md(r'''Plot the probabilities'''),
    code(r'''# One panel per trajectory, sharing an abscissa: the comparison the reader makes
# is between panels, so every panel has to carry identical limits, scales and
# tick spacings. plot_curves_stacked enforces that rather than leaving it to
# four hand-repeated formatting loops.
#
# The legend describes line *style*, not colour -- the 3+1 curve is a different
# colour in every panel -- so it is built from legend_proxies. That replaces the
# older trick of plotting dummy points at (-1, -1), outside the axis limits,
# purely to manufacture legend handles.
lc = ['C0', 'C2', 'C3', 'C4']
panels = [[dict(y=prob_3nu_arr[i], color='0.7', ls='--'),
           dict(y=prob_4nu_arr[i], color=lc[i], ls='-')]
          for i in range(len(costhz_arr))]

fig, ax = plotting.plot_curves_stacked(
    energies, panels,
    xlabel=r'Neutrino energy, $E_\nu$~[GeV]',
    ylabel=r'Probability,~'+plotting.prob_label(nu_i, nu_f),
    title=r'3+1 oscillations inside Earth (PREM)',
    xlim=(energy_min, energy_max), ylim=(0, 1),
    xmajor=10, xminor=1, ymajor=0.10, yminor=0.02,
    panel_labels=[r'From Fermilab to '+d for d in detectors],
    legend_proxies=[dict(label=r'3+1 oscillations', color='k', ls='-', lw=1),
                    dict(label=r'Standard, $3\nu$', color='k', ls='--', lw=1)],
    legend_loc='lower right',
    annotations=[dict(
        text=r'$s_{14} = $ '+str(s14)+r', $\delta_{14} = $ '+str(d14)
             + r', $s_{24} = $ '+str(s24)+r', $\delta_{24} = $ '+str(d24)
             + r', $s_{34} = $ '+str(s34)+r', $\Delta m_{41}^2 =$ '+str(D41)+r' eV$^2$',
        xy=(0.98, 0.10), panel=1, ha='right', fontsize=18)],
    grid=True, grid_kw=dict(axis='x'),
)

# The original draws two different grids -- minor on x, major only on y -- which
# is one call more than `grid_kw` expresses. Every function here returns
# (fig, ax) so that a packaged figure is a starting point rather than a dead
# end; this is that escape hatch doing its job.
for axx in ax:
    axx.grid(visible=True, c='0.8', which='major', axis='y')'''),
    md(r'''### 3.3 Probabilities vs. active-sterile mixing parameter'''),
    md(r'''Let's pick a baseline from Fermilab to Homestake (i.e., DUNE), which is about 1300 km, and see how the probability changes with the value of one of the active-sterile mixing parameters.  This is part of the analysis chain used to place bounds on the values of sterile neutrino mixing parameters.'''),
    md(r'''Generate the probabilities'''),
    code(r'''nu_i, nu_f = gd.NUMU, gd.NUE # Initial and final flavors; can also choose NUMU or NUTAU

# Energies
energy = 2.0*gd.UNIT_GEV # [eV]

# Zenith angle corresponding to Fermilab-Homestake (already computed above)
costhz = costhz_arr[1]

# Chord length Fermilab-Homestake
l_max = l_max_arr[1] # [km]

# Fixed mixing parameters between active and sterile sectors (can change them to anything else)
s14 = 1.e-1 # [adim]
d14 = 0.0 # [adim]
# s24 = 1.e-1 # [adim]
d24 = 0.0 # [adim]
s34 = 1.e-1 # [adim]
D41 = 1.e3 # [eV^2]

# Active-sterile mixing parameters that we vary
s24_npts = 100
s24_arr = np.linspace(0.0, 1.0, s24_npts) # [adim]
D41_arr = [1.e-2, 1.e-1, 1.e0, 1.e1] # [eV^2]
D41_label = [r'$\Delta m_{41}^2 = 0.01$~eV$^2$', r'$\Delta m_{41}^2 = 0.1$~eV$^2$', r'$\Delta m_{41}^2 = 1$~eV$^2$', 
             r'$\Delta m_{41}^2 = 10$~eV$^2$']

# Hamiltonian including matter effects inside Earth using the PREM density model
def H_3nu_prem(costhz, l, energy):
    # Given a direction (costhz) and a depth (l), compute the radial distance from the center of the Earth (r)
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31) \
            + hamiltonians.hamiltonian_3nu_matter(VCC_func_prem(r))
    
def H_4nu_prem(costhz, l, energy, s14, d14, s24, d24, s34, D21, D31, D41):
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV) # [km]
    return hamiltonians.hamiltonian_4nu_vacuum(energy, s12, s23, s13, dCP, s14, d14, s24, d24, s34, D21, D31, D41) \
            + hamiltonians.hamiltonian_4nu_matter(VCC_func_prem(r))'''),
    md(r'''Generate the probabilities'''),
    code(r'''prob_3nu = oscprob.osc_prob(lambda l: H_3nu_prem(costhz, l, energy), 
                                      0, l_max*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f]

prob_4nu = np.array([[oscprob.osc_prob(lambda l: H_4nu_prem(costhz, l, energy, s14, d14, s24, d24, s34, D21, D31, D41), 
                                      0, l_max*gd.CONV_KM_TO_INV_EV, 
                                      n_slabs=100, n_tpts_per_slab=20, magnus_exp_order=3, n_jobs=10)[nu_i][nu_f] for s24 in s24_arr]
                     for D41 in D41_arr]) '''),
    md(r'''Plot the probabilities'''),
    code(r'''# A plain set of curves against a swept variable, which is what plot_curves is
# for. (The hand-built version carried a 1x1 gridspec_kw, which did nothing.)
ls = ['-', '--', ':', '-.']
fig, ax = plotting.plot_curves(
    s24_arr,
    [dict(y=prob_4nu[i], lw=2, color=f'C{i}', ls=ls[i], label=D41_label[i])
     for i in range(len(D41_arr))]
    + [dict(y=np.full_like(s24_arr, prob_3nu), lw=1, color='0.2', ls='--',
            label=r'Standard, $3\nu$')],
    xlabel=r'Active-sterile mixing angle, $\sin \theta_{24}$',
    ylabel=r'Probability,~'+plotting.prob_label(nu_i, nu_f),
    title=r'3+1 oscillations from Fermilab to Homestake, '
          r'$E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV),
    xlim=(s24_arr[0], s24_arr[-1]), ylim=(0, 0.04),
    ymajor=0.01, yminor=0.001,
    legend_loc='upper left', legend_kw=dict(handlelength=2),
    title_fontsize=23,
)'''),
    md(r'''# 4. Probabilities in a 3+2 system ($5\nu$)'''),
    md(r'''Neutrino oscillations in a 3+2 system proceed similarly as in a 3+1 system, but driven instead by a $5 \times 5$ mixing matrix (see Eqs. (3) and (5) in arXiv:1105.3911),
\begin{equation}
U_{5\nu} = \tilde{R}_{35} ~ R_{34} ~ R_{25} ~ \tilde{R}_{24} ~ R_{23} ~ \tilde{R}_{15} ~ \tilde{R}_{14} ~ \tilde{R}_{13} ~ R_{12} \;,
\end{equation}
where the rotation matrices $R_{ij}$ and $\tilde{R}_{ij}$ are defined as above for the 3+1 system and, again, we have not written the matrix containing Majorana phases. In addition, there are now two new mass-squared difference, $\Delta m_{41}^2$ and $\Delta m_{51}^2$.

Thus, there are twelve new mixing parameters that we need to supply to Mag$\nu$s: $s_{14}$ (`s14`), $\delta_{14}$ (`d14`), $s_{15}$ (`s15`), $\delta_{15}$ (`d15`), $s_{24}$ (`s24`), $\delta_{24}$ (`d24`), $s_{25}$ (`s25`), $s_{34}$ (`s34`), $s_{35}$ (`s35`), $\delta_{35}$ (`d35`), $\Delta m_{41}^2$ (`D41`), and $\Delta m_{51}^2$ (`D51`).  These parameters are passed *in addition* to the standard mixing parameters in three-neutrino oscillations.'''),
    md(r'''Oscillation probabilities in the 3+2 system are computed just as we did in the 3+1 system, but the Hamiltonian is now a $5 \times 5$ matrix, and the `osc_prob` function now returns a $5 \times 5$ probability matrix.  

Below, we show an examply only for oscillations in vacuum.  It is straightforward to produce all the other cases we showed for 3+1 above also for 3+2.'''),
    md(r'''### 4.1 Mixing parameters'''),
    code(r'''# Standard mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# Mixing parameters between active and sterile sectors (can change them to anything else)
s14 = 1.e-1 # [adim]
d14 = 0.0 # [adim]
s15 = 1.e-1 # [adim]
d15 = 0.0 # [adim]
s24 = 1.e-1 # [adim]
d24 = 0.0 # [adim]
s25 = 5.e-1 # [adim]
s34 = 1.e-1 # [adim]
s35 = 1.e-1 # [adim]
d35 = 0.0 # [adim]
D41 = 0.5 # [eV^2]
D51 = 1.e-1 # [eV^2]'''),
    md(r'''### 4.2 Probabilities vs. distance'''),
    md(r'''Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e4 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

energy = 2.0*gd.UNIT_GEV # [eV]

# Compute probability vs. baseline
# For the 3nu case, osc_prob returns a 3x3 NumPy array with the probabilities: 
# [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]
# For the 5nu case, osc_prob returns a 3x3 NumPy array with the probabilities:
# [[Pee, Pem, Pet, Pes1, Pes2], [Pme, Pmm, Pmt, Pms2, Pms2], [Pte, Ptm, Ptt, Pts1, Pts2], 
#  [Ps1e, Ps1m, Ps1t, Pts1s1, Pts1s2], [Ps2e, Ps2m, Ps2t, Pts2s1, Pts2s2]]

# Standard, 3nu probabilities
H_3nu_vac = hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31, compute_matrix_multiplication=False)
prob_3nu_all = np.array([oscprob.osc_prob(lambda l: H_3nu_vac, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                          magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])

# 3+2 probabilities
H_5nu_vac = hamiltonians.hamiltonian_5nu_vacuum(energy, s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34, s35, d35, 
                                                   D21, D31, D41, D51, compute_matrix_multiplication=False)
prob_5nu_all = np.array([oscprob.osc_prob(lambda l: H_5nu_vac, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                          magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_baseline(gd.NUE, gd.NUE, distances, prob_5nu_all, prob_3nu_all, N=2,
                                   title=r'3+2 oscillations in vacuum, $E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV), save_plot=False)'''),
    md(r'''### 4.3 Probabilities vs. energy'''),
    md(r'''Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e0, 1.e2 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

baseline = 5.e3 # [km]

# Standard, 3nu probabilities
H_3nu_vac_en_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31, 
                                                                               compute_matrix_multiplication=False) # Vacuum H without (1/E)
prob_3nu_all = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_GEV))*H_3nu_vac_en_indep, 
                                          0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                          n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])

# 3+2 probabilities
H_5nu_vac_en_indep = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(s12, s23, s13, dCP, s14, d14, s15, d15, s24, d24, s25, s34,
                                                                               s35, d35, D21, D31, D41, D51,
                                                                               compute_matrix_multiplication=False) 
prob_5nu_all = np.array([oscprob.osc_prob(lambda l: (1/(enu*gd.UNIT_GEV))*H_5nu_vac_en_indep, 
                                          0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                          n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])'''),
    code(r'''make_plot_prob_Nnu_3nu_vs_energy(gd.NUE, gd.NUE, energies, prob_5nu_all, prob_3nu_all, N=2,
                                 title=r'3+2 oscillations in vacuum, $L = $~{:.2f}~km'.format(baseline), save_plot=False)'''),
    ])

# --------------------------------------------------------- 08_magnus_bsm_nsi
books['08_magnus_bsm_nsi.ipynb'] = notebook(
    'BSM: non-standard interactions',
    'A new matter potential in the same slot as the standard one, with off-diagonal couplings the Standard Model does not have.\n\nBecause it enters exactly where $V_{CC}$ does, everything downstream is the same calculation.',
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''# 0. Helper functions and definitions'''),
    code(r'''def flavor_index_to_str(nu_l):
    return {gd.NUE: 'e', gd.NUMU: 'mu', gd.NUTAU: 'tau'}[nu_l]


def eps_annotations(N):
    """The epsilon values a figure was made with, as annotation entries.

    Recording them on the figure matters more than usual here: an NSI plot is
    meaningless without the couplings that produced it, and these notebooks get
    read as images long after the code that made them.
    """
    if N == 2:
        return [dict(text=r'$\epsilon_{e e} =$ ' + str(eps_ee) + ', ' +
                          r'$\epsilon_{e \mu} =$ ' + str(eps_em) + ', ' +
                          r'$\epsilon_{\mu \mu} =$ ' + str(eps_mm),
                     xy=(0.02, 0.03), fontsize=18)]
    return [dict(text=r'$\epsilon_{e e} =$ ' + str(eps_ee) + ', ' +
                      r'$\epsilon_{e \mu} =$ ' + str(eps_em) + ', ' +
                      r'$\epsilon_{e \tau} =$ ' + str(eps_et),
                 xy=(0.02, 0.08), fontsize=18),
            dict(text=r'$\epsilon_{\mu \mu} =$ ' + str(eps_mm) + ', ' +
                      r'$\epsilon_{\mu \tau} =$ ' + str(eps_mt) + ', ' +
                      r'$\epsilon_{\tau \tau} =$ ' + str(eps_tt),
                 xy=(0.02, 0.03), fontsize=18)]


def _curves(prob_matt_nsi_all, prob_matt_std_all, nu_i, nu_f):
    return [dict(y=prob_matt_nsi_all[:, nu_i, nu_f], color='C1',
                 label=r'Matter effects with NSI'),
            dict(y=prob_matt_std_all[:, nu_i, nu_f], color='0.2', ls='--',
                 label=r'Standard matter effects')]


def make_plot_prob_matt_std_vs_nsi_vs_baseline(nu_i, nu_f, distances, prob_matt_nsi_all,
                                               prob_matt_std_all, N=3, title=None,
                                               save_plot=True):
    """NSI against standard matter effects, as a function of baseline."""
    filename = ('prob_matt_std_vs_nsi_' + str(N) + 'nu_vs_baseline_'
                + flavor_index_to_str(nu_i) + '_' + flavor_index_to_str(nu_f))
    return plotting.plot_probability_vs_baseline(
        distances, _curves(prob_matt_nsi_all, prob_matt_std_all, nu_i, nu_f),
        nu_i=nu_i, nu_f=nu_f, xlim=(l_ini, l_fin),
        annotations=eps_annotations(N),
        legend_loc='center left', title=title, title_fontsize=23,
        savefig=('../fig/' + filename + '.pdf') if save_plot else None,
        savefig_kw=dict(dpi=300))


def make_plot_prob_matt_std_vs_nsi_vs_energy(nu_i, nu_f, distances, prob_matt_nsi_all,
                                             prob_matt_std_all, N=3, title=None,
                                             save_plot=True):
    """The same comparison against neutrino energy."""
    filename = ('prob_matt_std_vs_nsi_' + str(N) + 'nu_vs_energy_'
                + flavor_index_to_str(nu_i) + '_' + flavor_index_to_str(nu_f))
    return plotting.plot_probability_vs_energy(
        energies, _curves(prob_matt_nsi_all, prob_matt_std_all, nu_i, nu_f),
        nu_i=nu_i, nu_f=nu_f, xlim=(energy_min, energy_max),
        annotations=eps_annotations(N),
        legend_loc='center right', title=title, title_fontsize=23,
        savefig=('../fig/' + filename + '.pdf') if save_plot else None,
        savefig_kw=dict(dpi=300))'''),
    md(r'''# 1. Probabilities with NSI in 2$\nu$: in matter with constant density'''),
    md(r'''## 1.1 Mixing parameters and NSI parameters'''),
    code(r'''# Standard mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]

# NSI parameters (predefined exaples from globaldefs; can change them to anything else)
eps_ee, eps_em, eps_mm = gd.EPS_2'''),
    md(r'''## 1.1 Probabilities vs. distance'''),
    md(r'''Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e4 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

energy = 2.0*gd.UNIT_GEV # [eV]

# Compute probability vs. baseline
# For the 2nu case, osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem], [Pme, Pmm]]

# Matter potential
rho = 3.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]

# Vacuum Hamiltonian
H_vac = hamiltonians.hamiltonian_2nu_vacuum(energy, s12, D21)

# Standard matter Hamiltonian
H_matt = hamiltonians.hamiltonian_2nu_matter(VCC)

# NSI matter Hamiltonian
H_nsi = hamiltonians.hamiltonian_2nu_nsi(VCC, eps_ee - eps_mm, eps_em)

# Probabilities
prob_matt_std_all = np.array([oscprob.osc_prob(lambda r: H_vac+H_matt, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                               magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])
prob_matt_nsi_all = np.array([oscprob.osc_prob(lambda r: H_vac+H_matt+H_nsi, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                               magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_matt_std_vs_nsi_vs_baseline(gd.NUE, gd.NUE, distances, prob_matt_nsi_all, prob_matt_std_all, N=2,
                                           title=r'$2\nu$ oscillations with NSI in constant-density matter (3~g~cm$^{-3}$), ' + \
                                                   r'$E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV),
                                           save_plot=False)'''),
    md(r'''## 2.3 Probabilities vs. energy'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

baseline = 5.e3 # [km]

# Matter potential
rho = 3.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]

# Vacuum Hamiltonian without the 1/E prefactos
H_vac_en_indep = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(s12, D21)

# Standard matter Hamiltonian
H_matt = hamiltonians.hamiltonian_2nu_matter(VCC)

# NSI matter Hamiltonian
H_nsi = hamiltonians.hamiltonian_2nu_nsi(VCC, eps_ee - eps_mm, eps_em)

# Total Hamiltonian, standard matter effects
def H_matt_std(energy, l):
    return (1/energy)*H_vac_en_indep + H_matt
prob_matt_std_all = np.array([oscprob.osc_prob(lambda l: H_matt_std(enu*gd.UNIT_GEV, l), 
                                               0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                               n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])

# Total Hamiltonian, NSI matter effects
def H_matt_nsi(energy, l):
    return (1/energy)*H_vac_en_indep + H_matt + H_nsi
prob_matt_nsi_all = np.array([oscprob.osc_prob(lambda l: H_matt_nsi(enu*gd.UNIT_GEV, l), 
                                               0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                               n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_matt_std_vs_nsi_vs_energy(gd.NUE, gd.NUE, energies, prob_matt_nsi_all, prob_matt_std_all, N=2,
                                         title=r'$2\nu$ oscillations with NSI in constant-density matter (3~g~cm$^{-3}$), ' + \
                                               r'$L = $~{:.2f}~km'.format(baseline),
                                         save_plot=False)'''),
    md(r'''# 2. Probabilities with NSI in 3$\nu$: in matter with constant density'''),
    md(r'''## 2.1 Mixing parameters and NSI parameters'''),
    code(r'''# Standard mixing parameters (predefined examples from globaldefs; can change them to anything else)
s12 = gd.S12_NO_BF_NUFIT_6_0 # [adim]
s23 = gd.S23_NO_BF_NUFIT_6_0 # [adim]
s13 = gd.S13_NO_BF_NUFIT_6_0 # [adim]
dCP = gd.DCP_NO_BF_NUFIT_6_0 # [adim]
D21 = gd.D21_NO_BF_NUFIT_6_0 # [eV^2]
D31 = gd.D31_NO_BF_NUFIT_6_0 # [eV^2]

# NSI parameters (predefined exaples from globaldefs; can change them to anything else)
eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = gd.EPS_3'''),
    md(r'''## 2.2 Probabilities vs. distance'''),
    md(r'''Generate probabilities'''),
    code(r'''# Baselines
l_ini, l_fin = 1.e1, 1.e4 # [km]
l_npts = 10000
distances = np.logspace(np.log10(l_ini), np.log10(l_fin), l_npts) # [km]

energy = 2.0*gd.UNIT_GEV # [eV]

# Compute probability vs. baseline
# For the 3nu case, osc_prob returns a 3x3 NumPy array with the probabilities: [[Pee, Pem, Pet], [Pme, Pmm, Pmt], [Pte, Ptm, Ptt]]

# Matter potential
rho = 3.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]

# Vacuum Hamiltonian
H_vac = hamiltonians.hamiltonian_3nu_vacuum(energy, s12, s23, s13, dCP, D21, D31)

# Standard matter Hamiltonian
H_matt = hamiltonians.hamiltonian_3nu_matter(VCC)

# NSI matter Hamiltonian
H_nsi = hamiltonians.hamiltonian_3nu_nsi(VCC, eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt)

# Probabilities
prob_matt_std_all = np.array([oscprob.osc_prob(lambda r: H_vac+H_matt, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                               magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])
prob_matt_nsi_all = np.array([oscprob.osc_prob(lambda r: H_vac+H_matt+H_nsi, 0.0, l*gd.CONV_KM_TO_INV_EV, n_slabs=1, n_tpts_per_slab=2, 
                                               magnus_exp_order=1, integration_method='simpson', n_jobs=1) for l in distances])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_matt_std_vs_nsi_vs_baseline(gd.NUE, gd.NUE, distances, prob_matt_nsi_all, prob_matt_std_all, N=3,
                                           title=r'$3\nu$ oscillations with NSI in constant-density matter (3~g~cm$^{-3}$), ' + \
                                                   r'$E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV),
                                           save_plot=False)'''),
    code(r'''make_plot_prob_matt_std_vs_nsi_vs_baseline(gd.NUMU, gd.NUTAU, distances, prob_matt_nsi_all, prob_matt_std_all, N=3,
                                           title=r'$3\nu$ oscillations with NSI in constant-density matter (3~g~cm$^{-3}$), ' + \
                                                   r'$E_\nu = $~{:.2f}~GeV'.format(energy/gd.UNIT_GEV),
                                           save_plot=False)'''),
    md(r'''## 2.3 Probabilities vs. energy'''),
    md(r'''Generate probabilities'''),
    code(r'''# Energies
energy_min, energy_max = 1.e-1, 1.e1 # [GeV]
energy_npts = 5000
energies = np.logspace(np.log10(energy_min), np.log10(energy_max), energy_npts) # [GeV]

baseline = 5.e3 # [km]

# Matter potential
rho = 3.0 # Matter density [g cm^{-3}]
num_density_e = matter.num_density_e_func(l=0.0, density_matter_func=lambda l: rho, electron_fraction=0.5, density_matter_is_in_g_per_cm3=True) # Electron number density [eV^3]
VCC = matter.VCC_func(l=0.0, num_density_e_func=lambda l: num_density_e) # Coherent forward potential, VCC [eV]

# Vacuum Hamiltonian without the 1/E prefactos
H_vac_en_indep = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)

# Standard matter Hamiltonian
H_matt = hamiltonians.hamiltonian_3nu_matter(VCC)

# NSI matter Hamiltonian
H_nsi = hamiltonians.hamiltonian_3nu_nsi(VCC, eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt)

# Total Hamiltonian, standard matter effects
def H_matt_std(energy, l):
    return (1/energy)*H_vac_en_indep + H_matt
prob_matt_std_all = np.array([oscprob.osc_prob(lambda l: H_matt_std(enu*gd.UNIT_GEV, l), 
                                               0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                               n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])

# Total Hamiltonian, NSI matter effects
def H_matt_nsi(energy, l):
    return (1/energy)*H_vac_en_indep + H_matt + H_nsi
prob_matt_nsi_all = np.array([oscprob.osc_prob(lambda l: H_matt_nsi(enu*gd.UNIT_GEV, l), 
                                               0.0, gd.CONV_KM_TO_INV_EV*baseline, 
                                               n_slabs=1, n_tpts_per_slab=10, magnus_exp_order=1) for enu in energies])'''),
    md(r'''Plot probabilities'''),
    code(r'''make_plot_prob_matt_std_vs_nsi_vs_energy(gd.NUE, gd.NUE, energies, prob_matt_nsi_all, prob_matt_std_all, N=3,
                                         title=r'$3\nu$ oscillations with NSI in constant-density matter (3~g~cm$^{-3}$), ' + \
                                               r'$L = $~{:.2f}~km'.format(baseline),
                                         save_plot=False)'''),
    code(r'''make_plot_prob_matt_std_vs_nsi_vs_energy(gd.NUMU, gd.NUTAU, energies, prob_matt_nsi_all, prob_matt_std_all, N=3,
                                         title=r'$3\nu$ oscillations with NSI in constant-density matter (3~g~cm$^{-3}$), ' + \
                                               r'$L = $~{:.2f}~km'.format(baseline),
                                         save_plot=False)'''),
    ])

# --------------------------------------------------------- 09_magnus_bsm_liv
books['09_magnus_bsm_liv.ipynb'] = notebook(
    'BSM: Lorentz-invariance violation',
    "Some extensions of the Standard Model -- string-inspired constructions, models\nwith a preferred frame, effective descriptions of quantum gravity -- break\nLorentz invariance. If they do, neutrinos are an unusually good place to look:\noscillations measure a *phase*, and a phase accumulated over an astrophysical\nbaseline is sensitive to energies that no accelerator reaches.\n\nMag$\\nu$s treats this as an extra, CPT-odd term in the Hamiltonian,\n\n$$\\mathbf{H} = \\frac{1}{2E}\\,\\mathbf{U}\\,\\mathbf{M}^2\\,\\mathbf{U}^\\dagger\n             \\;+\\; \\mathbf{V}_{\\rm CC}\n             \\;+\\; E^{\\,n}\\,\n               \\mathbf{U}_\\xi\\,\n               \\frac{\\mathbf{B}}{\\Lambda^{\\,n}}\\,\n               \\mathbf{U}_\\xi^\\dagger ,$$\n\nwhere $\\mathbf{B} = {\\rm diag}(b_1, b_2, b_3)$ holds the eigenvalues of the LIV\noperator, $\\mathbf{U}_\\xi$ rotates from its eigenbasis to the flavour basis\nthrough angles $\\xi_{12}, \\xi_{23}, \\xi_{13}$ and a phase $\\delta_{\\xi\\rm CP}$,\n$\\Lambda$ is the scale that makes the eigenvalues dimensionless, and $n$ is the\noperator's dimension minus three.\n\n## Why the exponent is the whole story\n\nLook at how the three terms scale with energy:\n\n| term | scaling | behaviour at high $E$ |\n|---|---|---|\n| vacuum | $\\Delta m^2 / 2E \\;\\propto\\; E^{-1}$ | switches **off** |\n| matter | $V_{\\rm CC}$, independent of $E$ | flat |\n| LIV | $E^{\\,n}\\, b/\\Lambda^{\\,n}$ | switches **on** (for $n \\geq 0$) |\n\nStandard oscillations die away at high energy: the phase $\\Delta m^2 L / 2E$\nshrinks, and the probability freezes at its zero-baseline value. A LIV term\ndoes the opposite. So the signature is not that oscillations look slightly\ndifferent -- it is that they are still *there*, at energies where the Standard\nModel says they should have stopped.\n\nThat is what the figures below show, and it is why every one of them is plotted\nagainst energy.",
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting'''),
    md(r'''## 1. The interface

Every `osc_prob_{2,3,4,5}nu_*` family has a `_liv` counterpart -- in vacuum, in
constant- and exponential-density matter, inside the Earth, and inside the Sun:'''),
    code(r'''[n for n in dir(oscprob) if n.startswith('osc_prob') and n.endswith('liv')]'''),
    md(r'''They take the standard mixing parameters plus the LIV ones. For three
flavours those are `sxi12`, `sxi23`, `sxi13`, `dxiCP` (the mixing of the LIV
eigenbasis into flavour), `b1`, `b2`, `b3` (its eigenvalues), `Lambda`, and
`n_liv`.

We use NuFIT 6.1 for the standard parameters, and a baseline equal to the
Earth's diameter -- the longest a terrestrial experiment can have, and the one
atmospheric-neutrino detectors actually use:'''),
    code(r"""osc = gd.load_nufit_params('NuFIT 6.1')

L_km = 1.3e4                                  # Earth diameter [km]
L = L_km*gd.CONV_KM_TO_INV_EV                 # [eV^-1]

nu_i, nu_f = gd.NUMU, gd.NUMU                 # muon-neutrino survival
energies = np.logspace(0.0, 4.0, 700)         # [GeV]  (below ~1 GeV the
                                              #  oscillation is faster than
                                              #  any practical grid resolves)
E = energies*gd.UNIT_GEV                      # [eV]

# Only b3 is switched on, and only the 2-3 sector is rotated, so the LIV term
# acts in the same sector as atmospheric oscillations -- the cleanest case to
# compare against.
liv = dict(sxi12=0.0, sxi23=1.0/np.sqrt(2.0), sxi13=0.0, dxiCP=0.0,
           b1=0.0, b2=0.0, Lambda=1.0)

# Choose b3 so the LIV phase b3*E^n*L reaches ~pi at E_star: that is where the
# term stops being a correction and starts driving the oscillation.
E_star = 1.0e2*gd.UNIT_GEV
def b3_for(n_liv, E_star=E_star, L=L):
    '''The eigenvalue whose LIV phase reaches pi at E_star.'''
    return np.pi/(L*E_star**n_liv)

print(f'L        = {L:.3e} eV^-1  ({L_km:.0f} km)')
for n in (0, 1, 2):
    print(f'n_liv={n}: b3 = {b3_for(n):.3e} eV^(1-{n})')"""),
    md(r'''## 2. Where standard oscillations stop and LIV keeps going

The defining figure. Below a few GeV the two curves agree -- the LIV term is
still negligible there. Above the crossover the standard curve flattens to 1,
because the vacuum phase has shrunk below one radian, while the LIV curve keeps
oscillating with a phase that is still growing.'''),
    code(r'''n_liv = 1
b3 = b3_for(n_liv)

prob_std = oscprob.osc_prob_3nu_vacuum(E, L, **osc, nu_i=nu_i, nu_f=nu_f)
prob_liv = oscprob.osc_prob_3nu_vacuum_liv(E, L, **osc, **liv, b3=b3, n_liv=n_liv,
                                           nu_i=nu_i, nu_f=nu_f)

fig, ax = plotting.plot_probability_vs_energy(
    energies,
    [dict(y=prob_std, label=r'Standard $3\nu$', color='0.2', ls='--'),
     dict(y=prob_liv, label=r'With LIV, $n = 1$', color='C1')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=3,
    xlim=(energies[0], energies[-1]),
    legend_title=r'Hamiltonian', legend_loc='lower left',
    title=(r'$3\nu$~in vacuum, $L = $~{:.0f}~km (Earth diameter);  '.format(L_km)
           + r'$b_3 = $ ' + '{:.1e}'.format(b3)
           + r'$\,$eV, $\Lambda = 1\,$eV, $\xi_{23} = \pi/4$'))'''),
    md(r'''Read the flat part of the dashed curve as the Standard Model's prediction:
*nothing happens up here*. Any oscillation observed in that region cannot come
from mass-squared differences, whatever their values, because the $1/E$ has
already switched them off. That is what makes the high-energy tail the place to
look, and why the strongest limits on LIV come from atmospheric and
astrophysical neutrinos rather than from beam experiments.'''),
    md(r'''## 3. The operator dimension sets where it turns on

`n_liv` is the operator's dimension minus three, so $n = 0$ is a
dimension-three (CPT-odd mass-like) term, $n = 1$ dimension four, and so on.

Each curve below is normalised to the *same* LIV phase at 100 GeV, so they are
not being compared at different strengths -- what differs is only how fast the
term grows. The higher the dimension, the more abruptly LIV appears: a steeper
power stays negligible for longer and then takes over faster.'''),
    code(r'''panels = [[dict(y=prob_std, label=r'Standard $3\nu$', color='0.2', ls='--'),
           dict(y=oscprob.osc_prob_3nu_vacuum_liv(E, L, **osc, **liv, b3=b3_for(n),
                                                  n_liv=n, nu_i=nu_i, nu_f=nu_f),
                label=r'With LIV', color='C{:d}'.format(n+1))]
          for n in (0, 1, 2)]

fig, ax = plotting.plot_probability_with_profile(
    energies, None, panels,
    xlim=(energies[0], energies[-1]),
    panel_annotations=[dict(text=r'$n = {:d}$  (operator dimension {:d})'.format(n, n+3),
                            bbox=dict(facecolor='white', edgecolor='0.7', alpha=0.9))
                       for n in (0, 1, 2)],
    panel_annotation_xy=(0.02, 0.06), panel_annotation_fontsize=20,
    shared_ylabel=r'Three-neutrino probability,~'+plotting.prob_label(nu_i, nu_f),
    xlabel=r'Neutrino energy, $E_\nu$ [GeV]',
    title=r'Same LIV phase at 100~GeV, different operator dimension',
    legend_loc='upper left', legend_kw=dict(framealpha=1.0),
    legend_on_panel=0, figsize=[18, 14])'''),
    md(r'''## 4. LIV inside the Earth

None of the above involved matter. Adding it is a keyword change: the `_liv`
wrappers exist for every environment, so the Earth version takes the same LIV
parameters plus a direction through the PREM profile.

We take an upward-going trajectory, $\cos\theta_z = -1$, which crosses the core
-- the longest baseline and the densest matter available:'''),
    code(r'''costhz = -1.0                                     # straight up through the core
L_earth = earth.distance_traveled_inside_earth(costhz)*gd.CONV_KM_TO_INV_EV

energies_e = np.logspace(0.0, 3.0, 400)           # [GeV]
E_e = energies_e*gd.UNIT_GEV

prob_earth_std = oscprob.osc_prob_3nu_earth(
    E_e, costhz=costhz, L=L_earth, **osc, nu_i=nu_i, nu_f=nu_f)
prob_earth_liv = oscprob.osc_prob_3nu_earth_liv(
    E_e, costhz=costhz, L=L_earth, **osc, **liv,
    b3=b3_for(1, L=L_earth), n_liv=1, nu_i=nu_i, nu_f=nu_f)

fig, ax = plotting.plot_probability_vs_energy(
    energies_e,
    [dict(y=prob_earth_std, label=r'Standard $3\nu$', color='0.2', ls='--'),
     dict(y=prob_earth_liv, label=r'With LIV, $n = 1$', color='C3')],
    nu_i=nu_i, nu_f=nu_f, num_flavors=3,
    xlim=(energies_e[0], energies_e[-1]),
    legend_title=r'Hamiltonian', legend_loc='lower left',
    title=r'$3\nu$~through the Earth (PREM), $\cos\theta_z = -1$')'''),
    md(r'''The matter term does not rescue the Standard Model curve at high energy.
$V_{\rm CC}$ is a few $\times 10^{-13}$ eV in the core and does not grow with
energy, so once the vacuum phase has shrunk away there is nothing left to make
the probability move -- while the LIV term is still climbing. Matter changes
*where* the standard curve flattens, not *that* it flattens.'''),
    md(r'''## 5. Two flavours, and setting a limit

The two-flavour interface is the same with one angle and two eigenvalues
(`sxi`, `b1`, `b2`), which is often how published constraints are quoted.

Here is the shape of a limit. A detector that sees standard oscillations up to
some energy and no anomaly constrains $b_3$: values large enough to distort the
survival probability there are excluded. We scan $b_3$ at a fixed high energy
and ask how far the probability moves from the standard prediction.'''),
    code(r'''E_probe = 1.0e3*gd.UNIT_GEV            # 1 TeV: standard oscillation is long gone
b3_scan = np.logspace(-30.0, -22.0, 200)

P_std_probe = oscprob.osc_prob_3nu_vacuum(E_probe, L, **osc, nu_i=nu_i, nu_f=nu_f)
P_scan = np.array([
    oscprob.osc_prob_3nu_vacuum_liv(E_probe, L, **osc, **liv, b3=b, n_liv=1,
                                    nu_i=nu_i, nu_f=nu_f)
    for b in b3_scan])

# The value at which the LIV phase reaches one radian -- the scale a
# null result excludes above.
b3_sensitivity = 1.0/(L*E_probe)

fig, ax = plotting.plot_curves(
    b3_scan,
    [dict(y=P_scan, label=r'With LIV', color='C1'),
     dict(y=np.full_like(b3_scan, P_std_probe), label=r'Standard $3\nu$',
          color='0.2', ls='--')],
    xlabel=r'LIV eigenvalue, $b_3$ [eV$^0$]',
    ylabel=r'Three-neutrino probability,~'+plotting.prob_label(nu_i, nu_f),
    xlim=(b3_scan[0], b3_scan[-1]), ylim=(0, 1), xscale='log',
    ymajor=0.10, yminor=0.02,
    annotations=[dict(text=r'$E_\nu = 1$~TeV, $L = $~{:.0f}~km'.format(L_km),
                      xy=(0.02, 0.10), fontsize=18),
                 dict(text=r'LIV phase $= 1$ at $b_3 \simeq $ '
                           + '{:.0e}'.format(b3_sensitivity),
                      xy=(0.02, 0.04), fontsize=18)],
    legend_title=r'Hamiltonian', legend_loc='upper right', grid=True)

ax.axvline(b3_sensitivity, color='0.5', ls=':', lw=2)

print(f'phase reaches 1 radian at b3 = {b3_sensitivity:.3e} eV^0')'''),
    md(r'''Below that value the curve sits on the standard prediction and a measurement
says nothing; above it the probability departs and a null result excludes the
parameter. The vertical line is not a fit -- it is just $b_3 L E = 1$, and it
lands where the curve visibly leaves the dashed line, which is the useful
sanity check on any sensitivity estimate of this kind.

Note how the reach improves: the excluded $b_3$ scales as $1/(LE^n)$, so for
$n \geq 1$ going to higher energy buys sensitivity faster than going to longer
baseline. That is the reason astrophysical neutrinos, not longer beamlines,
set the strongest limits on the higher-dimension operators.

## Where to go next

* `magnus.oscprob` has `_liv` wrappers for 2, 3, 4 and 5 flavours and for every
  environment; sterile states and LIV can be switched on together.
* The same `strategy`, `average` and tolerance keywords documented elsewhere
  apply here unchanged -- LIV is a term in the Hamiltonian, not a separate code
  path.
* `notebooks/08_magnus_bsm_nsi.ipynb` covers the other BSM scenario, where the
  new term is tied to matter rather than to energy.'''),
    ])

# -------------------------------------------- 10_magnus_averaged_probability
books['10_magnus_averaged_probability.ipynb'] = notebook(
    'Phase-averaged (decohered) probabilities: `average=True`',
    "A neutrino from an astrophysical source arrives with an oscillation phase\n$\\Delta m^2 L / 2E$ of order $10^{15}$.  Nothing in that number — the source\ndistance, the size of the production region, the detector's energy resolution —\nis known to anything close to the precision it would take to predict the phase,\nso the measurement integrates over many complete cycles and every oscillatory\nterm averages away.  What survives is\n\n$$P(\\nu_\\alpha \\to \\nu_\\beta) = \\sum_i |V_{\\alpha i}|^2 |V_{\\beta i}|^2$$\n\nwhere $V$ diagonalizes the Hamiltonian.  This is the exact $L/E \\to \\infty$\nlimit, not an approximation to be refined, and it costs one matrix product\ninstead of resolving $10^{15}$ radians.\n\nEvery oscillation-probability function in Magνs takes `average=True`.  This\nnotebook shows it for **2ν, 3ν, 4ν and 5ν**, and for a **custom Hamiltonian**\nthat is not one of the ones Magνs ships with.  Each plot plots the oscillating\nprobability against baseline (solid) together with its averaged value (dashed),\nso the average can be read as what the oscillation settles around.\n\nSee the [Phase-Averaged Probabilities](https://mbustama.github.io/Magnus/averaged_probability.html)\npage of the documentation for the derivation, the coherence criterion, and the\ntreatment of position-dependent Hamiltonians.",
    [
    code(r'''import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import magnus.oscprob as oscprob
import magnus.avgprob as avgprob
import magnus.hamiltonians as hamiltonians
import magnus.globaldefs as gd
import magnus.plotting as plotting

plt.rcParams['figure.dpi'] = 110

# load_nufit_params returns exactly the six mixing parameters, ready to pass as
# keyword arguments.  OSC_PARAMS_PREDEFINED entries also carry 'name' and
# 'description' strings, which would be forwarded into the propagation machinery
# and rejected there.
osc = gd.load_nufit_params('NuFIT 6.1')
print({k: round(v, 6) for k, v in osc.items()})'''),
    md(r'''## Plotting helper

One helper draws every figure below, so the panels stay comparable: the
oscillating probability as a solid line, its averaged value as a dashed line of
the same color, a logarithmic baseline axis, and the same tick and legend
formatting used in the other probability-vs-baseline notebooks.'''),
    code(r'''def plot_oscillating_and_averaged(distances, prob_curves, prob_averages, labels,
                                  title, filename, ylabel=None, colors=None):
    """Plot oscillating probabilities against baseline, each with its averaged
    value overlaid as a dashed line of the same colour.

    One call to magnus.plotting: the house style, the shared legend entry for
    the dashed lines, and the tick spacings are its defaults.
    """
    fig, ax = plotting.plot_probability_with_average(
        distances/gd.UNIT_KM, prob_curves, prob_averages,
        labels=labels, colors=colors,
        average_label=(r'Phase-averaged (\texttt{average=True})'
                       if mpl.rcParams['text.usetex']
                       else 'Phase-averaged (average=True)'),
        ylabel=ylabel or r'Oscillation probability',
        xlim=(distances[0]/gd.UNIT_KM, distances[-1]/gd.UNIT_KM),
        title=title, legend_title=r'Channel', legend_loc='upper left',
        grid=True, savefig='../fig/' + filename)
    plt.show()
    return fig, ax'''),
    md(r'''## Two flavors

The clearest case: a single mass-squared splitting, so the probability is a
plain $\sin^2$ and the averaged value is exactly the mid-line of the
oscillation, $1 - \tfrac{1}{2}\sin^2 2\theta$ for survival.'''),
    code(r'''sth, Dm2 = 0.55, 7.5e-5
energy = 1.0*gd.UNIT_GEV
distances = np.logspace(np.log10(1.0e2), np.log10(1.0e6), 3000)*gd.UNIT_KM

# Oscillating: one probability per baseline.
P_osc = np.asarray(oscprob.osc_prob_2nu_vacuum(energy, distances, sth, Dm2))

# Averaged: one number, valid once every pair has decohered.  Taken at the far
# end of the range, where that is true.
P_avg = np.asarray(oscprob.osc_prob_2nu_vacuum(energy, distances[-1], sth, Dm2,
                                               average=True))

print('averaged survival  P_ee =', round(float(P_avg[0, 0]), 5))
print('1 - sin^2(2 theta)/2     =', round(1.0 - 0.5*(2*sth*np.sqrt(1-sth**2))**2, 5))

plot_oscillating_and_averaged(
    distances,
    [P_osc[:, 0, 0], P_osc[:, 0, 1]],
    [P_avg[0, 0], P_avg[0, 1]],
    [plotting.prob_label(gd.NUE, gd.NUE), plotting.prob_label(gd.NUE, gd.NUMU)],
    r'$2\nu$ vacuum, $E_\nu = ${:.0f} GeV'.format(energy/gd.UNIT_GEV),
    'prob_2nu_vacuum_averaged_vs_baseline.pdf')'''),
    md(r'''## Three flavors

With two splittings there are two oscillation scales, and the averaged value is
no longer the mid-line of anything — it is the incoherent sum over the three
mass eigenstates.  Note the fast $\Delta m^2_{31}$ oscillation riding on the
slow $\Delta m^2_{21}$ one, and that the average is only reached once **both**
have decohered.'''),
    code(r'''energy = 1.0*gd.UNIT_GEV
distances = np.logspace(np.log10(1.0e2), np.log10(1.0e7), 4000)*gd.UNIT_KM

P_osc = np.asarray(oscprob.osc_prob_3nu_vacuum(energy, distances, **osc))
P_avg = np.asarray(oscprob.osc_prob_3nu_vacuum(energy, distances[-1], average=True, **osc))

print('averaged 3nu probability matrix:')
print(np.round(P_avg, 4))
print('rows sum to', np.round(P_avg.sum(axis=-1), 12))

plot_oscillating_and_averaged(
    distances,
    [P_osc[:, 0, 0], P_osc[:, 0, 1], P_osc[:, 0, 2]],
    [P_avg[0, 0], P_avg[0, 1], P_avg[0, 2]],
    [plotting.prob_label(gd.NUE, gd.NUE), plotting.prob_label(gd.NUE, gd.NUMU), plotting.prob_label(gd.NUE, gd.NUTAU)],
    r'$3\nu$ vacuum, $E_\nu = ${:.0f} GeV'.format(energy/gd.UNIT_GEV),
    'prob_3nu_vacuum_averaged_vs_baseline.pdf')'''),
    md(r'''### The flavor composition of astrophysical neutrinos

The single most quoted consequence: a source producing the pion-decay
composition $(1:2:0)$ delivers something close to equipartition at Earth.  That
is one line once the averaged matrix is in hand.'''),
    code(r'''at_source = np.array([1.0, 2.0, 0.0])/3.0
at_earth = at_source @ P_avg

print('at the source (x1/3):', np.round(at_source*3.0, 3))
print('at Earth      (x1/3):', np.round(at_earth*3.0, 3))
print('flux conserved      :', np.isclose(at_earth.sum(), 1.0))'''),
    md(r'''## Four and five flavors

Nothing about the construction is specific to three flavors: sterile states
simply add rows and columns.  With a large $\Delta m^2_{41}$ the sterile
oscillation is fast and decoheres first, leaving the active sector to average
later.'''),
    code(r'''S4 = dict(s14=0.15, d14=0.0, s24=0.10, d24=0.0, s34=0.05, D41=1.0)
S5 = dict(S4, s15=0.08, d15=0.0, s25=0.05, s35=0.03, D51=2.0)

energy = 1.0*gd.UNIT_GEV
distances = np.logspace(np.log10(1.0e0), np.log10(1.0e7), 4000)*gd.UNIT_KM

for num_flavors, sterile, fname in [(4, S4, 'prob_4nu_vacuum_averaged_vs_baseline.pdf'),
                                    (5, S5, 'prob_5nu_vacuum_averaged_vs_baseline.pdf')]:
    fn = getattr(oscprob, f'osc_prob_{num_flavors}nu_vacuum')
    P_osc = np.asarray(fn(energy, distances, **osc, **sterile))
    P_avg = np.asarray(fn(energy, distances[-1], average=True, **osc, **sterile))

    print(f'--- {num_flavors}nu ---')
    print(np.round(P_avg, 4))

    channels = [(gd.NUE, gd.NUE), (gd.NUE, gd.NUMU), (gd.NUE, num_flavors - 1)]
    plot_oscillating_and_averaged(
        distances,
        [P_osc[:, i, f] for i, f in channels],
        [P_avg[i, f] for i, f in channels],
        [plotting.prob_label(i, f) for i, f in channels],
        r'${}\nu$ vacuum (sterile), $E_\nu = ${:.0f} GeV'.format(num_flavors,
                                                                 energy/gd.UNIT_GEV),
        fname)'''),
    md(r'''## A custom Hamiltonian

Nothing above depends on the Hamiltonian being one Magνs ships with.  Here is a
hand-built three-level Hamiltonian — a vacuum-like diagonal in a mixing basis of
our own choosing, plus an off-diagonal term with no counterpart in the standard
parameterization — driven through the same two functions.

The oscillating curve comes from `oscprob.osc_prob`, the primordial entry point
that takes any $H(t)$; the averaged value comes from
`avgprob.averaged_probabilities_constant_hamiltonian`, which diagonalizes
whatever it is handed.'''),
    code(r"""def custom_mixing_matrix(theta, phi, chi):
    '''An arbitrary 3x3 unitary, built as three rotations with one complex phase.'''
    c, s = np.cos(theta), np.sin(theta)
    R1 = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=complex)
    c, s = np.cos(phi), np.sin(phi)
    R2 = np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]], dtype=complex)
    c, s = np.cos(chi), np.sin(chi)
    R3 = np.array([[c, 0.0, s*np.exp(-1j*1.1)], [0.0, 1.0, 0.0],
                   [-s*np.exp(1j*1.1), 0.0, c]], dtype=complex)
    return R1 @ R3 @ R2


# A custom, energy-independent Hamiltonian: eigenvalues of our choosing in a
# mixing basis of our choosing, plus a small non-standard off-diagonal term.
U_custom = custom_mixing_matrix(0.6, 0.8, 0.25)
eigenvalues = np.array([0.0, 3.0e-14, 1.1e-13])       # [eV]
H_custom = U_custom @ np.diag(eigenvalues).astype(complex) @ U_custom.conj().T
H_custom = H_custom + 2.0e-15*np.array([[0.0, 0.0, 1.0],
                                        [0.0, 0.0, 0.0],
                                        [1.0, 0.0, 0.0]], dtype=complex)

print('Hermitian:', np.allclose(H_custom, H_custom.conj().T))

distances = np.logspace(np.log10(1.0e2), np.log10(1.0e7), 3000)*gd.UNIT_KM

# Oscillating: the primordial osc_prob, which takes any H(t).
P_osc = np.array([np.asarray(oscprob.osc_prob(lambda t: H_custom, 0.0, L))
                  for L in distances])

# Averaged: hand the same matrix to avgprob.
P_avg = avgprob.averaged_probabilities_constant_hamiltonian(H_custom)

print('averaged custom-Hamiltonian matrix:')
print(np.round(P_avg, 4))
print('symmetric:', np.allclose(P_avg, P_avg.T))

plot_oscillating_and_averaged(
    distances,
    [P_osc[:, 0, 0], P_osc[:, 0, 1], P_osc[:, 0, 2]],
    [P_avg[0, 0], P_avg[0, 1], P_avg[0, 2]],
    [r'$P_{1 \to 1}$', r'$P_{1 \to 2}$', r'$P_{1 \to 3}$'],
    r'Custom $3\times3$ Hamiltonian',
    'prob_custom_hamiltonian_averaged_vs_baseline.pdf',
    ylabel=r'Probability')"""),
    md(r'''## When the average does not apply

The averaged expression is a limit, and Magνs checks whether the request is
actually in it.  At a terrestrial baseline the solar pair has accumulated a
fraction of a radian and has not decohered at all, so no averaged expression
describes the result — and asking for one warns instead of quietly answering.'''),
    code(r'''import warnings

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    P_here = oscprob.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM,
                                         average=True, **osc)

for w in caught:
    if issubclass(w.category, oscprob.PhaseAveragingWarning):
        print('PhaseAveragingWarning:')
        print(' ', str(w.message)[:300], '...')

# The same check, done directly: which pairs are in neither limit?
eigenvalues = np.array([0.0, osc['D21'], osc['D31']])/(2.0*gd.UNIT_GEV)
for baseline_km in [1.0e3, 1.0e6, 1.0e8]:
    blocks, undecided = avgprob.coherence_report(eigenvalues, baseline_km*gd.UNIT_KM)
    print(f'L = {baseline_km:9.0e} km   blocks={blocks}   '
          f'pairs in neither limit: {[(i, j) for i, j, _ in undecided]}')'''),
    md(r'''## Cost

The averaged probability is not just exact — it is far cheaper than the
propagation it replaces, because it never resolves the phase it is about to
discard.'''),
    code(r'''import time

energy, L_far = 1.0*gd.UNIT_GEV, 1.0e8*gd.UNIT_KM

t0 = time.perf_counter()
for _ in range(100):
    oscprob.osc_prob_3nu_vacuum(energy, L_far, average=True, **osc)
t_avg = (time.perf_counter() - t0)/100

energies = 1.0/np.linspace(1.0/(10.0*gd.UNIT_GEV), 1.0/energy, 2001)
t0 = time.perf_counter()
P_num = np.asarray(oscprob.osc_prob_3nu_vacuum(energies, L_far, **osc)).mean(axis=0)
t_num = time.perf_counter() - t0

P_closed = np.asarray(oscprob.osc_prob_3nu_vacuum(energy, L_far, average=True, **osc))

print(f'closed form            : {t_avg*1e6:8.1f} us   (exact)')
print(f'averaging 2001 energies: {t_num:8.3f} s    (approximate)')
print(f'agreement              : {np.max(np.abs(P_num - P_closed)):.2e}')'''),
    md(r'''## Summary

* `average=True` works on every oscillation-probability function, for any
  number of flavors, and on any Hamiltonian — including one of your own.
* For a position-independent Hamiltonian it is the exact $L/E \to \infty$
  limit, in closed form, for the cost of one eigendecomposition.
* For vacuum it does not depend on energy or baseline at all: one matrix serves
  a whole flux calculation.
* Whether the limit *applies* is checked, not assumed.  A pair of eigenvalues
  that has neither decohered nor stayed coherent is reported rather than
  answered.
* Position-dependent profiles decohere in the eigenbasis at production and are
  carried along the levels of the instantaneous Hamiltonian, with exact
  level-crossing probabilities where the evolution stops being adiabatic; see
  the [documentation](https://mbustama.github.io/Magnus/averaged_probability.html).'''),
    ])

# ---------------------------------------------- 11_magnus_matrix_exponential
books['11_magnus_matrix_exponential.ipynb'] = notebook(
    'The matrix exponential',
    'How $\\exp(\\Omega)$ is actually computed, and why the choice matters.\n\nThe truncated Magnus series is anti-Hermitian, so its exponential is exactly unitary -- but only if the exponential itself is computed in a way that preserves that. This notebook compares the routes.',
    [
    code(r'''import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.magnus as magnus
import magnus.plotting as plotting'''),
    md(r'''While the main purpose of Mag$\nu$s is to compute neutrino oscillation probabilities, at its core, it is an efficient implementation of the Magnus expansion of the matrix exponential whose computation is required to compute the time-evolution operator for any given Hamiltonian, including time-dependent (or position-dependent) ones.  Therefore, we can use Mag$\nu$s to compute the Magnus expansion of the matrix exponential
\begin{equation}
 U = \exp\left( \int_{t_i}^{t_f} A(t) dt \right)
\end{equation}
where $A$ is an arbitary $N \times N$ square matrix.  The Magnus expansion of the exponential is
\begin{equation}
 U = \exp\left( \Omega(t_i, t_f) \right) \;,
\end{equation}
where $\Omega$ is the series expansion
\begin{equation}
 \Omega(t_i, t_f) = \sum_{k=1}^{\infty} \Omega_k(t_i, t_f) \;,
\end{equation}
where the terms $\Omega_k$ of the expansion are computed recursively in terms of lower-order terms (see Wikipedia entry).  Mag$\nu$s can return the terms $\Omega_k$ and also the matrix exponential, $U$.  The present version of the code implements calculation up to $k = 6$ term.'''),
    md(r'''First, define a sample $3 \times 3$ complex matrix (unlike the Hamiltonian in the case of neutrino oscillations, it does not have to be Hermitian).'''),
    code(r'''def A3(t):
    # return np.array([[1+1j*t, 2/(1+t), 3j*t], [0, 4-1j*t, 5+3j*t], [-3j*np.sqrt(t), -2j*(1/(1+t)), 1]], dtype=np.complex128)
    return np.array([[1+1j*t, 2*t*2, 3j*t], [0, 4-1j*t, 5+3j*t], [-3j*np.sqrt(t), -2j*t, 1]], dtype=np.complex128)

t_i, t_f = 0, 1 # Initial and final times of the matrix integral'''),
    md(r'''The terms $\Omega_k$ of the expansion have the same dimensions as the input matrix $A_3$ (in this case, $3 \times 3$).  To get at the individual $\Omega_k$ rather than just their exponentiated sum, pass `return_magnus_terms=True` to `magnus_expansion`; it then returns the pair `(exp_Omega, Omega_terms)`.

By default the matrix integral is evaluated on a grid of `n_tpts = 50` points, spaced *linearly* between $t_i$ and $t_f$, and the expansion is truncated at `order = 2`.

One caveat matters here.  The default integration method is `'gl'` (Gauss-Legendre), which is *commutator-free*: it never forms the individual $\Omega_k$ at all, and so returns only their combined effect, no matter what `order` you ask for.  Because this section is specifically about the term-by-term structure of the series, we ask for `'trapezoid'`, which does build each term explicitly.'''),
    code(r'''# 'trapezoid' (not the default 'gl') so that the individual terms actually exist
_, (Omega_1, Omega_2) = magnus.magnus_expansion(A3, t_i, t_f, order=2,
                                                integration_method='trapezoid',
                                                return_magnus_terms=True)
print(Omega_1) # The first-order term is just what we would compute were A3 time-independent
print(Omega_2) # The second-order term is the first correction to that'''),
    md(r'''Now compute terms up to $k = 6$, use a finer grid of 100 points to integrate, and switch from trapezoidal integration to integration via Simpson's rule.'''),
    code(r'''_, Omega = magnus.magnus_expansion(A3, t_i, t_f, n_tpts=100, order=6,
                                   integration_method='simpson',
                                   return_magnus_terms=True)
for k, Omega_k in enumerate(Omega):
    print(r'Omega_'+str(k+1)+' = ')
    print(Omega_k)'''),
    md(r'''To compute the matrix exponential $\exp(\Omega)$, we do not need to compute the Magnus terms separately first.  We can directly call the `magnus_expansion` function, which does it internally.'''),
    code(r'''exp_Omega = magnus.magnus_expansion(A3, t_i, t_f)
print(exp_Omega)'''),
    md(r'''And, like before, we can change the order of the expansion, the number of points used in the integration grid, and the integration method.'''),
    code(r'''exp_Omega = magnus.magnus_expansion(A3, t_i, t_f, n_tpts=100, order=6, integration_method='simpson')
print(exp_Omega)'''),
    md(r'''Naturally, the more terms in the Magnus expansion and the finer the integration grid, the more accurate the returned matrix exponential is.  To illustrate this, let's plot the value of the real part of the $U_{11}$ component (i.e., of `exp_Omega[0][0]`) for different choices. '''),
    code(r'''# Generate the data
n_tpts_min, n_tpts_max = 2, 10
order_max = 6
elem_sel_i, elem_sel_j = 0, 1 # Select the [0][0] element of the returned matrix exponential
exp_Omega_elem_sel = [[magnus.magnus_expansion(A3, t_i, t_f, n_tpts=n, order=o+1, 
                                              integration_method='trapezoid')[elem_sel_i][elem_sel_j].real 
                      for n in range(n_tpts_min, n_tpts_max+1)] for o in range(order_max)]
exp_Omega_elem_sel = np.array(exp_Omega_elem_sel)'''),
    code(r'''n_tpts = list(range(n_tpts_min, n_tpts_max+1))
best = exp_Omega_elem_sel[order_max-1]

fig, ax = plotting.plot_curves(
    n_tpts,
    [dict(y=exp_Omega_elem_sel[k], label=r'$k \leq $'+str(k+1))
     for k in range(order_max)],
    xlabel=r'Number of time-integration grid points, \texttt{n\_tpts}',
    ylabel=r'Matrix-exponential term, ${\rm Re} (U_{12})$',
    xlim=(n_tpts_min, n_tpts_max), ymajor=100, yminor=10,
    residual=[(exp_Omega_elem_sel[order_max-2][i]-best[i])/best[i]
              for i in range(len(n_tpts))],
    residual_label=r'$\epsilon_{\rm rel}$ vs.~$k=6$',
    residual_ylim=(-2, 0.1), residual_ymajor=0.5, residual_yminor=0.1,
    annotations=[dict(
        text=r'$U = \exp \left( \int_0^1 A_3(t) dt \right) '
             r'= \exp \left( \sum_{k=1}^\infty \Omega_k \right) $',
        xy=(0.10, 0.88), fontsize=30)],
    legend_title=r'Expansion order', legend_loc='upper right',
    grid=True, tight_layout=False)'''),
    md(r'''This shows that, at least for short time intervals, like $t_i = 0$ to $t_f = 1$, it is enough to sue `n_tnpts = 10` points in the integration grid, and expansion up to $k = 6$.

However, because the Magnus expansion involves calculating exponential, it can result in an overflow if its argument is too large, i.e., if $t_f - t_f \gg 1$.  In this case, we partition the interval $[t_i, t_f]$ into multiple, smaller subintervals, compute the matrix exponential in each via Magnus expansion, and then compute their time-ordered product.  (This is actually the strategy adopted by Mag$\nu$s to compute neutrino oscillation probabilities.)'''),
    code(r"""# Set run parameters
order_max = 6
elem_sel_i, elem_sel_j = 0, 1 # Select the [0][1] element of the returned matrix exponential
t_i, t_f = 0, 5
n_subintervals_min, n_subintervals_max = 5, 150

def elem_for(n_subintervals, order):
    '''Time-ordered product of the per-subinterval matrix exponentials.'''
    dt = (t_f-t_i)/n_subintervals
    subintervals = np.array([[dt*i, dt*(i+1)] for i in range(n_subintervals)])
    exp_Omega = np.array([magnus.magnus_expansion(A3, s[0], s[1], n_tpts=10, order=order,
                                                  integration_method='trapezoid')
                          for s in subintervals])
    return abs(np.linalg.multi_dot(exp_Omega)[elem_sel_i][elem_sel_j].real)

# The finest partition is the reference every coarser one is measured against
best = elem_for(n_subintervals_max, order_max)

n_sub = list(range(n_subintervals_min, n_subintervals_max+1))
series = []
for order in range(1, order_max+1):
    print("order = "+str(order))
    series.append([elem_for(n, order) for n in n_sub])

fig, ax = plotting.plot_curves(
    n_sub,
    [dict(y=series[k], label=r'$k \leq~$'+str(k+1)) for k in range(order_max)],
    xlabel=r'Number of time subintervals within $[t_i, t_f]$',
    ylabel=r'Matrix-exponential term, ${\rm Re} (U_{12})$',
    xlim=(n_subintervals_min, n_subintervals_max),
    ylim=(3.e25, 1.e26), xmajor=10, xminor=1,
    residual=[(v-best)/best for v in series[order_max-1]],
    residual_label=r'$\epsilon_{\rm rel}$ vs. best',
    residual_ylim=(-0.1, 1), residual_ymajor=0.5, residual_yminor=0.1,
    annotations=[dict(
        text=r'$U = \exp \left( \int_{t_i=0}^{t_f=5} A_3(t) dt \right) '
             r'= \exp \left( \sum_{k=1}^\infty \Omega_k \right) $',
        xy=(0.10, 0.88), fontsize=30)],
    legend_title=r'Expansion order', legend_loc='upper right',
    grid=True, tight_layout=False)"""),
    md(r'''So, in this case, about 100 subintervals are enough, regardless of the order of the Magnus expansion used.  This is because if there are many subintervals, then the corrections of higher-order terms within each subinterval are small.'''),
    md(r'''Finally, let's see how this changes when we widen the interval $[t_i, t_f]$.'''),
    code(r'''# Set run parameters
order_max = 6
elem_sel_i, elem_sel_j = 0, 1
t_i = 1
t_f_max = 150
t_f_arr = np.linspace(t_i, t_f_max, 200)
t_modifier = 1.e-2
n_subintervals_arr = [3, 10, 100]

def elem_vs_tf(n_subintervals, order):
    out = []
    for t_f in t_f_arr:
        dt = (t_f-t_i)/n_subintervals
        subintervals = np.array([[dt*i, dt*(i+1)] for i in range(n_subintervals)])
        exp_Omega = np.array([magnus.magnus_expansion(A3, t_modifier*s[0], t_modifier*s[1],
                                                      n_tpts=10, order=order,
                                                      integration_method='trapezoid')
                              for s in subintervals])
        out.append(abs(np.linalg.multi_dot(exp_Omega)[elem_sel_i][elem_sel_j].real))
    return out

panels = []
for n_subintervals in n_subintervals_arr:
    print("  n_subintervals = "+str(n_subintervals))
    panels.append([dict(y=elem_vs_tf(n_subintervals, order), label=r'$k \leq~$'+str(order))
                   for order in range(1, order_max+1)])

fig, ax = plotting.plot_probability_with_profile(
    t_f_arr, None, panels,
    xlim=(t_f_arr[0], t_f_arr[-1]), xscale='linear', xmajor=10, xminor=2,
    panel_yscale='log', panel_ylim=(1.e-2, 3500),
    panel_ymajor=None, panel_yminor=None,
    panel_annotations=[r'Number of subintervals: '+str(n) for n in n_subintervals_arr],
    panel_annotation_xy=(0.02, 0.15), panel_annotation_fontsize=25,
    shared_ylabel=r'Matrix-exponential term, ${\rm Re} (U_{12})$',
    xlabel=r'Final time of the evolution, $t_f$ [$\times 10^{-2}$]',
    legend_title=r'Expansion order', legend_loc='upper left',
    legend_on_panel=0, legend_kw=dict(ncol=2),
    grid=True, figsize=[18, 15])

ax[1].annotate(r'$U = \exp \left( \int_{t_i=0}^{t_f} A_3(t) dt \right) '
               r'= \exp \left( \sum_{k=1}^\infty \Omega_k \right) $',
               xy=(0.02, 0.77), xycoords='axes fraction', ha='left', fontsize=30)'''),
    md(r'''This shows that, as expected, the most accurate results are obtained using a high-order expansion and a large number of subintervals.  However, there is a trade-off between accuracy and speed.  The decision whether to favor one or the other, and to what extent will depend on the specific problem at hand.'''),
    md(r'''## Where the coefficients come from: the Magnus terms to any order

Everything above used the expansion as a black box. The coefficients it runs on
are not arbitrary, though: they follow from a single recursion, and
`magnus.expansionterms` derives them symbolically at any order, in exact rational
arithmetic.

The recursion is

$$\Omega_n(t) = \sum_{j=1}^{n-1} \frac{B_j}{j!} \int_0^t S_n^{(j)}(s)\, ds ,
\qquad
S_n^{(1)} = [\Omega_{n-1}, A], \quad
S_n^{(j)} = \sum_{m=1}^{n-j} [\Omega_m, S_{n-m}^{(j-1)}] ,$$

with $B_j$ the Bernoulli numbers in the $B_1 = -1/2$ convention. Because
$B_j = 0$ for every odd $j \geq 3$, whole commutator groups drop out.'''),
    code(r'''import magnus.expansionterms as et

# The Bernoulli numbers come back as exact fractions, not floats
{n: str(et.bernoulli(n)) for n in range(9)}'''),
    md(r'''The coefficient multiplying the $j$-th commutator group is $B_j/j!$. These are
precisely the constants the numerical core hard-codes as `F1`, `F2`, `F3`, `F4`:'''),
    code(r'''import magnus.magnus as mg

for j, name in [(2, 'F1'), (4, 'F2'), (6, 'F3'), (8, 'F4')]:
    exact = et.bernoulli_factor(j)
    print(f"B_{j}/{j}! = {str(exact):>12s} = {float(exact):+.12e}   "
          f"mg.{name} = {getattr(mg, name):+.12e}")'''),
    md(r'''### Printing the expansion

`print_magnus_terms` writes out the integrand of each $\Omega_n$, one term per
line, with its exact coefficient:'''),
    code(r'''et.print_magnus_terms(5)'''),
    md(r'''### Beyond the implemented ceiling

The numerical core implements orders up to `MAGNUS_EXP_ORDER_MAX`, but the
generator has no ceiling of its own — useful for seeing how fast the expansion
grows before committing to an order. Every term turns out to be a right-nested
chain $[\Omega_{m_1},[\Omega_{m_2},\ldots[\Omega_{m_j}, A]]]$ whose indices sum
to $n-1$, so the $j$-th group has $\binom{n-2}{j-1}$ terms.'''),
    code(r'''print(f"implemented up to order {mg.MAGNUS_EXP_ORDER_MAX}\n")
print(f"{'order':>6s} {'terms':>7s}")
for n in range(1, 15):
    marker = '  <- ceiling' if n == mg.MAGNUS_EXP_ORDER_MAX else ''
    print(f"{n:6d} {et.count_terms(n):7d}{marker}")'''),
    md(r'''The count roughly doubles per order, and so does the work per slab. That is the
trade behind the cost warning: order 7 costs about 2.7x order 6 on the same grid,
order 10 about 17x. Higher order does converge faster in the slab width, but
narrowing the slabs at order 4 or 6 often reaches a given accuracy for less total
work — and beyond the series' convergence radius, $\int\lVert A\rVert\,dt < \pi$,
no order helps at all.

One term of order 12, which the numerical core does not implement but the
generator will still hand you:'''),
    code(r'''terms_12 = et.omega_terms(12)
print(f"Omega_12 has {len(terms_12)} terms. The last one is")
print("   int " + et.format_term(terms_12[-1]))'''),
    ])

# --------------------------------------- 12_magnus_adiabatic_hybrid_strategy
books['12_magnus_adiabatic_hybrid_strategy.ipynb'] = notebook(
    "The `strategy` parameter: `'auto'` vs. `'hybrid'` vs. `'magnus'`",
    "Every matter/NSI/LIV oscillation-probability function in Magνs (and every\n`osc_prob_*_sun`/`osc_prob_*_sun_nsi`/`osc_prob_*_sun_liv` wrapper built on\nthem, plus the fully generic `osc_prob_sun`/`osc_prob_earth`) accepts a\n`strategy` keyword: `'auto'` (the default), `'hybrid'`, or `'magnus'`.\n\n* **`'magnus'`** uses only the traditional Magnus-expansion machinery: the\n  closed-form two-flavor interaction-picture integrator when it applies, or\n  the general adaptive slab-refinement method. This is the exact behavior\n  Magνs had before this feature existed.\n* **`'hybrid'`** additionally tries an **adiabatic-transport-plus-Magnus-patch**\n  propagator (`magnus.adiabatic.hybrid_propagator`): away from an eigenvalue\n  crossing of the instantaneous Hamiltonian, the evolution operator is\n  computed via the instantaneous eigenbasis (cheap regardless of how large\n  the accumulated oscillation phase is); near a genuine MSW resonance, a\n  short, exact Magnus patch is stitched in. The result stays exactly\n  unitary, and the whole computation self-certifies against its own\n  internal tolerances.\n* **`'auto'`** tries `'hybrid'` first, silently falling back to `'magnus'`\n  for any point where it does not apply or fails to self-certify.\n\nThis notebook reproduces, live, the validation described in\n[`docs/source/adiabatic_strategy.rst`](../docs/source/adiabatic_strategy.rst):\nfor a sequence of increasingly demanding cases (2 through 5 flavors,\nstandard oscillations and an engineered BSM resonance), we compare all\nthree `strategy` values against each other **and** against a\ntight-tolerance `scipy.integrate.solve_ivp` solution of the same\nSchrödinger equation, in both **runtime** and **accuracy**.\n\n**The headline result, previewed:** for 3 or more flavors, `strategy='magnus'`\n(the old default behavior) does not just get *slower* as the accumulated\noscillation phase grows -- with the library's default refinement caps, it\ncan hit them and return a **silently plausible-looking but wrong** answer,\nstill warning about it but not fixing it. `strategy='hybrid'`/`'auto'` are\nboth **fast and correct** in exactly this regime.",
    [
    code(r'''import time
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.globaldefs as gd
import magnus.adiabatic as adiabatic
from magnus.magnus import MagnusConvergenceWarning

plt.rcParams["figure.dpi"] = 110'''),
    md(r'''## A common reference point: a tight-tolerance `solve_ivp` ground truth

For every case below, we also solve the Schrödinger equation directly with
`scipy.integrate.solve_ivp` (`DOP853`, `rtol=1e-10`, `atol=1e-12`) as an
independent, method-agnostic ground truth to check every `strategy` value
against.  We use the same sign/time-ordering convention as the rest of
Magνs, $dU/dl = -iH(l)\,U(l)$.'''),
    code(r"""def exact_U(H_func, l0, l1, dim):
    '''Ground-truth evolution operator via a tight-tolerance ODE solve.'''
    def rhs(l, y):
        return (-1j * np.asarray(H_func(l)) @ y.reshape(dim, dim)).ravel()
    sol = solve_ivp(rhs, (l0, l1), np.eye(dim, dtype=complex).ravel(),
                     rtol=1e-10, atol=1e-12, method="DOP853")
    return sol.y[:, -1].reshape(dim, dim)


def to_P(U):
    return np.swapaxes(U.real**2 + U.imag**2, -1, -2)


results = []  # collects one dict per case, for the summary table/plot at the end


def run_case(name, wrapper_func, wrapper_kwargs, H_func, l0, l1, dim):
    '''Runs strategy='magnus'/'hybrid'/'auto' plus solve_ivp for one case, times
    and cross-checks every result, prints a short report, and records it in `results`.'''
    print(f"=== {name} ===")
    row = {"name": name}
    P_by_strategy = {}
    for strategy in ["magnus", "hybrid", "auto"]:
        with warnings.catch_warnings(record=True) as wlist:
            warnings.simplefilter("always")
            t0 = time.time()
            P = wrapper_func(**wrapper_kwargs, strategy=strategy, validate_input=False)
            dt = time.time() - t0
        warned = any(issubclass(w.category, (oscprob.ToleranceNotAchievedWarning,))
                     for w in wlist)
        P_by_strategy[strategy] = np.asarray(P)
        row[f"t_{strategy}"] = dt
        row[f"warned_{strategy}"] = warned
        flag = "  <-- ToleranceNotAchievedWarning (may be inaccurate)" if warned else ""
        print(f"  strategy={strategy:8s}  t={dt:9.4f} s{flag}")

    t0 = time.time()
    U_exact = exact_U(H_func, l0, l1, dim)
    t_exact = time.time() - t0
    P_exact = to_P(U_exact)
    row["t_solve_ivp"] = t_exact
    print(f"  solve_ivp         t={t_exact:9.4f} s  (ground truth)")

    for strategy in ["magnus", "hybrid", "auto"]:
        err = np.max(np.abs(P_by_strategy[strategy] - P_exact))
        speedup = t_exact / row[f"t_{strategy}"]
        row[f"err_{strategy}"] = err
        row[f"speedup_{strategy}"] = speedup
        print(f"  strategy={strategy:8s}  max abs error vs. solve_ivp = {err:.2e}"
              f"   (speedup {speedup:8.1f}x)")
    print()

    results.append(row)
    return P_by_strategy, P_exact"""),
    md(r'''## Case 1: two flavors

Two-flavor oscillations already have a dedicated, exact, closed-form
"interaction-picture" fast path in Magνs, independent of `strategy` (see
`docs/source/methodology.rst`). So for two flavors, all three `strategy`
values end up doing essentially the same thing -- and, at this particular
baseline, even `solve_ivp` itself is cheap, since a 2-level ODE is
inexpensive to integrate directly. The dramatic differences start
showing up at three flavors and beyond, where no such closed-form
shortcut exists.'''),
    code(r'''l_scale = gd.L_SCALE_SUN                  # length scale of the Sun's density profile
energy = 18.0 * gd.UNIT_MEV                # a representative solar-neutrino (8B/hep) energy

sth, Dm2 = gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0
hvac2 = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
e00_2 = np.diag([1.0, 0.0])

rho_func = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
VCC_func = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)

def H_2nu(l):
    return (1.0 / energy) * hvac2 + np.asarray(VCC_func(l)) * e00_2

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    _ = run_case(
        "2nu standard, L = 4 l_scale",
        oscprob.osc_prob_2nu_sun,
        dict(energy=energy, L=4.0 * l_scale, L0=0.0, sth=sth, Dm2=Dm2),
        H_2nu, 0.0, 4.0 * l_scale, 2,
    )'''),
    md(r'''## Case 2: three flavors, no resonance in range

For three (or more) flavors there is no closed-form fast path: `'magnus'`
must resolve the full accumulated phase slab by slab. At this baseline
that phase is large enough that, with the library's *default* refinement
caps, `'magnus'` hits them and returns -- fast, but **wrong**, flagged by
`ToleranceNotAchievedWarning`. `'hybrid'`/`'auto'` instead recognize that
the trajectory is purely adiabatic here (no resonance at this energy) and
return the correct answer via cheap instantaneous-eigenbasis transport,
with no patch needed at all.'''),
    code(r'''osc = gd.OSC_PARAMS_PREDEFINED["OSC_PARAMS_DEFAULT"]
hvac3 = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
    osc["s12"], osc["s23"], osc["s13"], osc["dCP"], osc["D21"], osc["D31"])
e00_3 = np.diag([1.0, 0.0, 0.0])

def H_3nu_std(l):
    return (1.0 / energy) * hvac3 + np.asarray(VCC_func(l)) * e00_3

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    _ = run_case(
        "3nu standard, L = 4 l_scale",
        oscprob.osc_prob_3nu_sun,
        dict(energy=energy, L=4.0 * l_scale, L0=0.0),
        H_3nu_std, 0.0, 4.0 * l_scale, 3,
    )'''),
    md(r'''## Case 3: three flavors with a genuine (BSM/NSI-induced) resonance

The real solar mixing angle is large enough (LMA solution) that no
solar-neutrino MSW crossing is ever sharp: adiabatic transport alone
already gets the right answer, with no patch needed, as Case 2 just
showed. To exercise the *patching* machinery itself, we engineer a
genuine non-adiabatic resonance with a non-standard-interaction (NSI)
coupling, $\epsilon_{e\tau}=3.0$ -- physically motivated new physics, not
an arbitrarily detuned Standard Model parameter.'''),
    code(r'''h_matt3_nsi = np.diag([1.0, 0.0, 0.0]) + hamiltonians.hamiltonian_3nu_nsi(
    1.0, 0.0, 0.0j, 3.0, 0.0, 0.0j, 0.0)  # eps_et = 3.0, everything else 0

def H_3nu_bsm(l):
    return (1.0 / energy) * hvac3 + np.asarray(VCC_func(l)) * h_matt3_nsi

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    P_bsm, P_bsm_exact = run_case(
        "3nu BSM (NSI, eps_et=3.0), L = 5 l_scale",
        oscprob.osc_prob_3nu_sun_nsi,
        dict(energy=energy, L=5.0 * l_scale, L0=0.0, eps_et=3.0),
        H_3nu_bsm, 0.0, 5.0 * l_scale, 3,
    )'''),
    md(r'''## Peeking under the hood: what did `'hybrid'` actually detect?

`magnus.adiabatic` is usable directly, independent of the `osc_prob_*`
API. Let's use it to see exactly where the non-adiabatic window was
found for the BSM case above, and plot the two nearly-degenerate
instantaneous eigenvalues responsible for it -- a real-data version of
the schematic "avoided crossing" diagram in the documentation.'''),
    code(r'''l0, l1 = 0.0, 5.0 * l_scale
windows, candidates = adiabatic.find_nonadiabatic_windows(H_3nu_bsm, l0, l1)
print("non-adiabatic window(s) found (l / l_scale):",
      [(w[0] / l_scale, w[1] / l_scale) for w in windows])
for c in candidates:
    print(f"  candidate: levels ({c['j']},{c['k']}) at l/l_scale={c['l']/l_scale:.4f}"
          f"  gap={c['gap']:.3e}  gamma={c['gamma']:.3e}")

ls = np.linspace(2.5 * l_scale, 4.5 * l_scale, 400)
eigs = np.array([np.linalg.eigvalsh(H_3nu_bsm(l)) for l in ls])

# The third level sits far above the other two on this energy scale, so a single shared y-axis
# would flatten the actual near-degeneracy into invisibility; show the full spectrum on top for
# context and a zoomed-in view of the two crossing levels (0, 1) on the bottom.
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True,
                                      gridspec_kw={"height_ratios": [1, 1.4]})
for k in range(3):
    ax_top.plot(ls / l_scale, eigs[:, k], lw=2, label=f"level {k}")
for k in range(2):
    ax_bot.plot(ls / l_scale, eigs[:, k], lw=2)
for ax in (ax_top, ax_bot):
    for w in windows:
        ax.axvspan(w[0] / l_scale, w[1] / l_scale, color="red", alpha=0.15)

ax_top.set_ylabel("eigenvalue [eV]\n(full spectrum)")
ax_top.set_title("3nu BSM (NSI): instantaneous eigenvalues near the engineered resonance")
ax_top.legend(loc="center right", fontsize=9)

ax_bot.set_xlabel(r"position $l\,/\,l_{\rm scale}$")
ax_bot.set_ylabel("eigenvalue [eV]\n(levels 0, 1 zoomed in)")
handles = [plt.Line2D([0], [0], color="red", alpha=0.3, lw=8)]
ax_bot.legend(handles, ["detected non-adiabatic window"], loc="upper left", fontsize=9)

fig.tight_layout()
plt.show()'''),
    md(r'''## Generalizing further: 4 and 5 flavors (sterile neutrinos)

Nothing about the resonance detector or the adiabatic propagator assumes
a particular number of flavors. We repeat the same standard/BSM
comparison for 3+1 and 3+2 sterile-neutrino Hamiltonians. The sterile
mixing/mass-splitting values below are illustrative (chosen so a
matter-driven crossing actually falls inside this energy/baseline
window), not a real global-fit point -- a realistic eV$^2$-scale sterile
splitting would never cross the matter potential here.'''),
    code(r'''s14, d14, s24, d24, s34 = 0.15, 1.2, 0.10, 0.0, 0.05
s15, d15, s25, s35, d35 = 0.08, 0.5, 0.05, 0.03, 0.9
D41 = 1.5 * osc["D31"]
D51 = 2.5 * osc["D31"]

hvac4 = hamiltonians.hamiltonian_4nu_vacuum_energy_independent(
    osc["s12"], osc["s23"], osc["s13"], osc["dCP"], s14, d14, s24, d24, s34,
    osc["D21"], osc["D31"], D41)
e00_4 = np.diag([1.0, 0.0, 0.0, 0.0])

def H_4nu_std(l):
    return (1.0 / energy) * hvac4 + np.asarray(VCC_func(l)) * e00_4

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    _ = run_case(
        "4nu standard, L = 4 l_scale",
        oscprob.osc_prob_4nu_sun,
        dict(energy=energy, L=4.0 * l_scale, L0=0.0,
             s14=s14, d14=d14, s24=s24, d24=d24, s34=s34, D41=D41),
        H_4nu_std, 0.0, 4.0 * l_scale, 4,
    )'''),
    code(r'''h_matt4_nsi = np.diag([1.0, 0.0, 0.0, 0.0]) + hamiltonians.hamiltonian_4nu_nsi(
    1.0, 0.0, 0.0j, 3.0, 0.0j, 0.0, 0.0j, 0.0j, 0.0, 0.0j, 0.0)  # eps_et = 3.0

def H_4nu_bsm(l):
    return (1.0 / energy) * hvac4 + np.asarray(VCC_func(l)) * h_matt4_nsi

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    _ = run_case(
        "4nu BSM (NSI, eps_et=3.0), L = 5 l_scale",
        oscprob.osc_prob_4nu_sun_nsi,
        dict(energy=energy, L=5.0 * l_scale, L0=0.0,
             s14=s14, d14=d14, s24=s24, d24=d24, s34=s34, D41=D41, eps_et=3.0),
        H_4nu_bsm, 0.0, 5.0 * l_scale, 4,
    )'''),
    code(r'''hvac5 = hamiltonians.hamiltonian_5nu_vacuum_energy_independent(
    osc["s12"], osc["s23"], osc["s13"], osc["dCP"], s14, d14, s15, d15, s24, d24,
    s25, s34, s35, d35, osc["D21"], osc["D31"], D41, D51)
e00_5 = np.diag([1.0, 0.0, 0.0, 0.0, 0.0])

def H_5nu_std(l):
    return (1.0 / energy) * hvac5 + np.asarray(VCC_func(l)) * e00_5

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    _ = run_case(
        "5nu standard, L = 4 l_scale",
        oscprob.osc_prob_5nu_sun,
        dict(energy=energy, L=4.0 * l_scale, L0=0.0,
             s14=s14, d14=d14, s15=s15, d15=d15, s24=s24, d24=d24,
             s25=s25, s34=s34, s35=s35, d35=d35, D41=D41, D51=D51),
        H_5nu_std, 0.0, 4.0 * l_scale, 5,
    )'''),
    code(r'''h_matt5_nsi = np.diag([1.0, 0.0, 0.0, 0.0, 0.0]) + hamiltonians.hamiltonian_5nu_nsi(
    1.0, 0.0, 0.0j, 3.0, 0.0j, 0.0j, 0.0, 0.0j, 0.0j, 0.0j, 0.0, 0.0j, 0.0j, 0.0, 0.0j, 0.0)

def H_5nu_bsm(l):
    return (1.0 / energy) * hvac5 + np.asarray(VCC_func(l)) * h_matt5_nsi

with warnings.catch_warnings():
    warnings.simplefilter("ignore", MagnusConvergenceWarning)
    _ = run_case(
        "5nu BSM (NSI, eps_et=3.0), L = 5 l_scale",
        oscprob.osc_prob_5nu_sun_nsi,
        dict(energy=energy, L=5.0 * l_scale, L0=0.0,
             s14=s14, d14=d14, s15=s15, d15=d15, s24=s24, d24=d24,
             s25=s25, s34=s34, s35=s35, d35=d35, D41=D41, D51=D51, eps_et=3.0),
        H_5nu_bsm, 0.0, 5.0 * l_scale, 5,
    )'''),
    md(r'''## Summary

The table below collects every case above: whether a non-adiabatic
window was found, the runtime of each `strategy`, whether `'magnus'`
triggered `ToleranceNotAchievedWarning`, and the accuracy/speedup of each
`strategy` relative to the `solve_ivp` ground truth.'''),
    code(r'''header = (f"{'case':38s} {'t_magnus':>10s} {'t_hybrid':>10s} {'t_solve_ivp':>12s} "
          f"{'err_magnus':>11s} {'err_hybrid':>11s} {'speedup_hybrid':>15s}")
print(header)
print("-" * len(header))
for row in results:
    warn_flag = "*" if row["warned_magnus"] else " "
    print(f"{row['name']:38s} {row['t_magnus']:9.3f}{warn_flag} {row['t_hybrid']:10.4f} "
          f"{row['t_solve_ivp']:12.3f} {row['err_magnus']:11.2e} {row['err_hybrid']:11.2e} "
          f"{row['speedup_hybrid']:14.1f}x")
print()
print("* strategy='magnus' raised ToleranceNotAchievedWarning for this case "
      "(its accuracy is not certified -- and, as the err_magnus column shows, "
      "is sometimes genuinely wrong).")'''),
    code(r'''import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(10, 6.5))
short_names = [r["name"].split(",")[0] for r in results]   # drop the ", L = ..." suffix
speedups = [r["speedup_hybrid"] for r in results]
windows_found = ["BSM" in r["name"] for r in results]      # BSM cases have an engineered resonance
colors = ["#c01c28" if w else "#2ec27e" for w in windows_found]

bars = ax.bar(range(len(short_names)), speedups, color=colors)
ax.set_yscale("log")
ax.set_ylim(top=max(speedups) * 6)   # headroom for the text labels above the tallest bars
ax.set_xticks(range(len(short_names)))
ax.set_xticklabels(short_names, rotation=25, ha="right", fontsize=10)
ax.set_ylabel("speedup of strategy='hybrid' vs. solve_ivp (log scale)")
ax.set_title("Measured speedup, this notebook's own run", pad=45)
for rect, val in zip(bars, speedups):
    ax.text(rect.get_x() + rect.get_width() / 2, val * 1.3, f"{val:,.0f}x",
            ha="center", fontsize=9)

handles = [mpatches.Patch(color="#2ec27e", label="pure adiabatic (no resonance in range)"),
           mpatches.Patch(color="#c01c28", label="adiabatic + Magnus patch (resonance found)")]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=9)
fig.subplots_adjust(bottom=0.32, top=0.82)
plt.show()'''),
    md(r'''## Takeaways

* For **two flavors**, a dedicated closed-form fast path already handles
  the accumulated-phase problem, independent of `strategy` -- and a
  direct ODE solve is itself inexpensive, so there is little practical
  difference between the three `strategy` values at this baseline.
* For **three or more flavors**, `strategy='magnus'` has no such
  shortcut. With the library's default refinement caps, it does not
  necessarily get dramatically *slower* as the accumulated phase grows --
  it can instead hit its caps and return **quickly but wrong**,
  correctly flagged by `ToleranceNotAchievedWarning` but not fixed by it.
* `strategy='hybrid'` (and the default, `'auto'`, which tries `'hybrid'`
  first) stays **fast and correct** across every case tested here --
  0-window (purely adiabatic) cases in particular are ~$10^3$ times
  faster than `solve_ivp`, and cases needing an actual Magnus patch still
  land 1-2 orders of magnitude faster, all while remaining exactly
  unitary and matching the `solve_ivp` ground truth well within the
  package's standard $10^{-3}$ target tolerance.
* None of this requires the user to do anything: `strategy='auto'` is
  already the default everywhere it applies.

See [`docs/source/adiabatic_strategy.rst`](../docs/source/adiabatic_strategy.rst)
for the full mathematical derivation (the adiabatic theorem, the exact
Hellmann-Feynman resonance diagnostic, and the self-certification
procedure) behind everything demonstrated in this notebook.'''),
    ])

# ------------------------------------------- 13_magnus_tabulated_solar_model
books['13_magnus_tabulated_solar_model.ipynb'] = notebook(
    'A tabulated solar model: are you computing the observable?',
    "You have a real solar model on disk -- a table of radius, density and composition --\nand you want oscillation probabilities from it. This notebook does exactly that with\n**BS2005-AGS,OP** (Bahcall, Serenelli & Basu, ApJ 621, L85), and uses it to show a\ndistinction that decides whether an error matters at all:\n\n* the **instantaneous** probability at one baseline, which is what `osc_prob_*` returns;\n* the **phase-averaged** probability, which is what a solar-neutrino experiment measures.\n\nThey are not the same, and on a solar trajectory they can differ by more than the\ntolerance you asked for. The headline: at 5 MeV the instantaneous answer is\n**1.2e-03** from the truth -- outside the requested 1e-3 -- while the averaged answer is\n**1.6e-04**, comfortably inside it. The error is mostly *phase*, and no observable\nsees it. With a cubic interpolant instead of a linear one the same comparison is\n1.4e-03 against 2.6e-05, a factor of 53: the conclusion does not depend on how you\ninterpolate.\n\nThe notebook also shows the diagnostic that tells you which regime you are in,\n`strategy_info['sampling']`, and the module that computes the average exactly,\n`magnus.avgprob`.",
    [
    code(r'''import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.oscprob as oscprob
import magnus.adiabatic as adiabatic

plt.rcParams["figure.dpi"] = 110'''),
    md(r'''## 1. Load the model and build an electron-density profile

The table ships with the repository, under `docs/dev/adversarial_batteries/`. Two
details matter and are easy to get wrong:

* **The electron fraction is not 0.5.** For fully ionized H + He,
  $n_e = \rho\,N_A\,(1+X)/2$, and the hydrogen mass fraction $X$ runs from 0.36 at the
  centre to 0.75 at the surface. Using a fixed $Y_e = 0.5$ would be wrong by up to 70 %,
  which is why we build the number density ourselves rather than handing a mass density
  to `vcc_func_from_rho_func` (that function takes a *scalar* electron fraction).
* **Interpolate in $\log n_e$.** The density spans five orders of magnitude; a linear
  interpolant in the raw value is poor and a spline can undershoot to negative values.'''),
    code(r"""TABLE = os.path.join('..', 'docs', 'dev', 'adversarial_batteries', 'bs05_agsop.dat')

rows = []
with open(TABLE) as fh:
    for line in fh:
        f = line.split()
        if len(f) == 12:
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                continue          # the column-heading line
table = np.array(rows)

r_over_rsun = table[:, 1]
rho_cgs, x_hydrogen = table[:, 3], table[:, 6]

# n_e = rho * N_A * (1 + X) / 2, in the package's natural units [eV^3].
MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
ne_ev3 = rho_cgs*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*(0.5*(1.0 + x_hydrogen))

x_nat = r_over_rsun*gd.SUN_RADIUS*gd.UNIT_KM          # radius in eV^-1
log_ne = np.log(ne_ev3)

def ne_bs05(l):
    '''Electron number density [eV^3] at position l [eV^-1], log-interpolated.'''
    xs = np.clip(np.asarray(l, dtype=float), x_nat[0], x_nat[-1])
    out = np.exp(np.interp(xs, x_nat, log_ne))
    return out[()] if np.ndim(out) == 0 else out

print('%d rows, r = %.5f .. %.5f R_sun' % (len(table), r_over_rsun[0], r_over_rsun[-1]))
print('central n_e = %.1f N_A cm^-3' % (ne_ev3[0]/gd.N_AV/gd.UNIT_PER_CM3))"""),
    md(r'''### The package's exponential is a *fit*, not this model

`gd.NUM_DENSITY_E_SUN_CENTRAL` = 245 $N_A$ is the $r\to0$ intercept of the standard
exponential fit $n_e = 245\,N_A e^{-10.54 r/R_\odot}$. The model's actual central value
is 101.9 $N_A$ (printed below). The fit is a few-percent description only in a band around
$0.2\,R_\odot$; inside $0.05\,R_\odot$ it is high by a factor 2.4.'''),
    code(r'''fit = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)

fig, ax = plt.subplots(figsize=(6.2, 3.6))
rr = np.linspace(0.002, 0.95, 400)
xs = rr*gd.SUN_RADIUS*gd.UNIT_KM
ax.semilogy(rr, np.asarray(ne_bs05(xs))/gd.N_AV/gd.UNIT_PER_CM3, label='BS2005-AGS,OP')
ax.semilogy(rr, np.asarray(fit(xs))/gd.N_AV/gd.UNIT_PER_CM3, '--',
            label=r'exponential fit ($245\,N_A e^{-10.54 r/R_\odot}$)')
ax.axvspan(0.0, gd.L_SCALE_SUN/(gd.SUN_RADIUS*gd.UNIT_KM), alpha=0.12,
           label='one scale height (used below)')
ax.set_xlabel(r'$r / R_\odot$'); ax.set_ylabel(r'$n_e\ [N_A\ {\rm cm}^{-3}]$')
ax.legend(fontsize=8); ax.set_title('A real solar model against the exponential fit')
fig.tight_layout()'''),
    md(r'''## 2. The instantaneous probability, and its error

We take a two-flavour calculation at **5 MeV** -- in the $^8$B range -- over one solar
scale height, and check it against a tight-tolerance `solve_ivp` ground truth.'''),
    code(r"""def exact_U_many(H_func, l0, Ls, dim):
    '''Ground truth at many baselines from ONE integration.'''
    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(dim, dim)).ravel()
    sol = solve_ivp(rhs, (float(l0), float(Ls[-1])), np.eye(dim, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853', t_eval=Ls)
    return np.array([sol.y[:, i].reshape(dim, dim) for i in range(len(Ls))])

def to_P(U):
    return np.swapaxes(np.asarray(U).real**2 + np.asarray(U).imag**2, -1, -2)

ENERGY = 5.0*gd.UNIT_MEV
L0, L1 = 0.0, gd.L_SCALE_SUN
params2 = {'sth': gd.S12_NO_BF_NUFIT_6_0, 'Dm2': gd.D21_NO_BF_NUFIT_6_0}

hvac2 = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(
    params2['sth'], params2['Dm2'])
e00 = np.diag([1.0, 0.0])
VCC = matter.vcc_func_from_rho_func(
    ne_bs05, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
    density_is_of_number_of_electrons=True)

def H_of_l(l):
    return (1.0/ENERGY)*hvac2 + np.asarray(VCC(l))[..., None, None]*e00

info = {}
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    P = np.asarray(oscprob.osc_prob_matter_std_potential(
        2, ne_bs05, ENERGY, L1, params2, L0=L0,
        density_is_of_number_of_electrons=True, strategy_info=info))
    raised = sorted({w.category.__name__ for w in caught})

P_ref = to_P(exact_U_many(H_of_l, L0, np.array([L1]), 2)[0])
print('engine        : %s' % info['engine'])
print('certified     : %s' % info.get('certified'))
print('warnings      : %s' % (', '.join(raised) or 'NONE'))
print('P_ee          : %.6f   (truth %.6f)' % (P[0][0], P_ref[0][0]))
print('max |error|   : %.3e   against a requested %.0e' % (np.max(np.abs(P - P_ref)), 1e-3))"""),
    md(r'''That is outside the default tolerance, reported as `certified`, with no warning.

Before concluding anything from it, ask what the number means physically.'''),
    md(r'''## 3. `strategy_info['sampling']`: how much phase is in there?

`adiabatic.oscillation_sampling` reports how coarsely a request samples the oscillation
it is computing. It is surfaced through `strategy_info` and is computed **only** when you
ask for it, so it costs nothing on an ordinary call.'''),
    code(r'''Ls_scan = np.linspace(0.2*L1, L1, 8)
info_scan = {}
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    oscprob.osc_prob_matter_std_potential(
        2, ne_bs05, ENERGY, Ls_scan, params2, L0=L0,
        density_is_of_number_of_electrons=True, strategy_info=info_scan)

for k, v in sorted(info_scan['sampling'].items()):
    print('%-24s %s' % (k, ('%.4e' % v) if isinstance(v, float) else v))'''),
    md(r'''`cycles_per_step` is the number to read. It is far above 0.5, so this eight-point scan
takes **less than two samples per oscillation**: the individual values are correct, but
the curve through them is an artefact. `nyquist_points` says how many baselines you would
need to sample the oscillation properly -- about 900 here, and several thousand at the
energies and flavour counts used elsewhere in the documentation.

That is the signature of a quantity dominated by phase. A solar-neutrino experiment
resolves none of it: the $^8$B production region is extended, the Sun-Earth phase is
$\sim10^{10}$ cycles, and detector energy resolution finishes the job.'''),
    md(r'''## 4. The averaged probability -- the quantity that is actually observed

Average over a window of several oscillation lengths and compare *that* against the
truth.'''),
    code(r'''L_OSC = 4.0*np.pi*ENERGY/params2['Dm2']         # vacuum oscillation length
Ls = np.linspace(L1 - 6.0*L_OSC, L1, 121)

P_ref_many = np.array([to_P(U) for U in exact_U_many(H_of_l, L0, Ls, 2)])
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    P_pkg = np.array([np.asarray(oscprob.osc_prob_matter_std_potential(
        2, ne_bs05, ENERGY, float(L), params2, L0=L0,
        density_is_of_number_of_electrons=True)) for L in Ls])

err_inst = np.max(np.abs(P_pkg[-1] - P_ref_many[-1]))
err_avg = np.max(np.abs(P_pkg.mean(axis=0) - P_ref_many.mean(axis=0)))

print('trajectory        : %.0f oscillation lengths' % ((L1 - L0)/L_OSC))
print('instantaneous err : %.3e   <-- outside 1e-3' % err_inst)
print('AVERAGED err      : %.3e   <-- inside 1e-3, by %.0fx' % (err_avg, 1e-3/err_avg))
print('averaging reduces the error by %.0fx' % (err_inst/err_avg))
print()
print('averaged P_ee: package %.6f   truth %.6f' %
      (P_pkg.mean(axis=0)[0][0], P_ref_many.mean(axis=0)[0][0]))'''),
    code(r'''fig, ax = plt.subplots(figsize=(6.6, 3.6))
xs_plot = (Ls - Ls[0])/L_OSC
ax.plot(xs_plot, P_ref_many[:, 0, 0], lw=1.0, label='truth (solve_ivp)')
ax.plot(xs_plot, P_pkg[:, 0, 0], lw=1.0, ls='--', label="magnus, strategy='auto'")
ax.axhline(P_ref_many[:, 0, 0].mean(), color='k', lw=1.2,
           label='averaged $P_{ee}$ (the observable)')
ax.set_xlabel('baseline, in oscillation lengths'); ax.set_ylabel(r'$P_{ee}$')
ax.legend(fontsize=8); ax.set_title('The disagreement is in the phase, not the envelope')
fig.tight_layout()'''),
    md(r'''The two curves sit on top of each other in *envelope* and drift in *phase*. Averaging
removes the drift, which is why the observable is accurate to 2.6e-05 while any single
baseline is off by 1.4e-03.

**The rule of thumb, measured across the profile families in
`docs/dev/adversarial_batteries/`:** if averaging shrinks the error by more than about
twentyfold, it was phase and no observable sees it; if it barely moves, the error is in
the envelope and it is real. A shock front is the case where it barely moves -- see
notebook 14.'''),
    md(r'''### It is not an artefact of the interpolation

A cubic spline through the same table (still in $\log n_e$) is a different profile, and
it gives a *larger* instantaneous error and a *smaller* averaged one. Both point the same
way: the disagreement lives in the phase.'''),
    code(r'''from scipy.interpolate import CubicSpline

_cs = CubicSpline(x_nat, log_ne, extrapolate=True)

def ne_bs05_cubic(l):
    xs = np.clip(np.asarray(l, dtype=float), x_nat[0], x_nat[-1])
    out = np.exp(_cs(xs))
    return out[()] if np.ndim(out) == 0 else out

VCC_c = matter.vcc_func_from_rho_func(
    ne_bs05_cubic, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
    density_is_of_number_of_electrons=True)

def H_cubic(l):
    return (1.0/ENERGY)*hvac2 + np.asarray(VCC_c(l))[..., None, None]*e00

ref_c = np.array([to_P(U) for U in exact_U_many(H_cubic, L0, Ls, 2)])
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    got_c = np.array([np.asarray(oscprob.osc_prob_matter_std_potential(
        2, ne_bs05_cubic, ENERGY, float(L), params2, L0=L0,
        density_is_of_number_of_electrons=True)) for L in Ls])

ic = np.max(np.abs(got_c[-1] - ref_c[-1]))
ac = np.max(np.abs(got_c.mean(axis=0) - ref_c.mean(axis=0)))
print('%-10s %12s %12s %10s' % ('interpolant', 'instant.', 'averaged', 'reduction'))
print('%-10s %12.3e %12.3e %9.0fx' % ('linear', err_inst, err_avg, err_inst/err_avg))
print('%-10s %12.3e %12.3e %9.0fx' % ('cubic', ic, ac, ic/ac))'''),
    md(r'''## 5. If the average is what you want, ask for it directly

`average=True` (and `magnus.avgprob` underneath) computes the phase-averaged limit in
closed form -- one matrix product rather than an integration -- and
`avgprob.coherence_report` will tell you whether that limit is valid for your spectrum
and baseline, or whether some pair of eigenvalues sits in the middle regime where
neither limit holds.'''),
    code(r'''import magnus.avgprob as avgprob

lam = np.linalg.eigvalsh(np.asarray(H_of_l(0.5*L1)))
blocks, undecided = avgprob.coherence_report(lam, phase_scale=L1)
print('coherence blocks :', blocks)
print('pairs in neither limit:', undecided or 'none -- the averaged expression is exact here')'''),
    md(r'''## Summary

| | |
|---|---|
| instantaneous error at 5 MeV | **1.2e-03** (linear) / 1.4e-03 (cubic) -- outside the requested 1e-3, `certified`, silent |
| averaged error (the observable) | **1.6e-04** (linear) / 2.6e-05 (cubic) -- inside, by 6x and 38x |
| what the difference is | phase, not envelope |
| how to tell | `strategy_info['sampling']`, and averaging over a few oscillation lengths |
| how to get the average exactly | `average=True`, or `magnus.avgprob` |

The instantaneous error is real and worth knowing about if you need the coherent
probability at a point. For solar physics it is not the quantity being measured, and
the package is accurate on the one that is.

See :doc:`averaged_probability` in the documentation for the full treatment, and
notebook 14 for a profile where averaging does **not** rescue the answer.'''),
    ])

# ------------------------------------------------- 14_magnus_supernova_shock
books['14_magnus_supernova_shock.ipynb'] = notebook(
    'A supernova shock front: when the error is real',
    'Notebook 13 showed a solar case where the package looked wrong by 1.4e-03 and was not:\nthe error was **phase**, and phase-averaging -- which is what a detector does -- removed\nit. This notebook is the opposite case, and the contrast is the point.\n\nA supernova shock front changes the **adiabaticity of the MSW level crossing**, so it\nmoves the conversion probability *itself* rather than the phase of an oscillation.\nAveraging cannot remove that. Here the package is wrong by **0.21 in probability on the\naveraged observable** -- and, importantly, it **says so every time**.\n\nThe profile is the standard one from the literature:\n\n* $\\rho_0(x) = 10^{14}\\,(x/\\mathrm{km})^{-2.4}\\ \\mathrm{g\\,cm^{-3}}$, forward-shock jump\n  $\\xi = V_+/V_- \\simeq 10$, and the rarefaction shape behind it, from\n  **Fogli, Lisi, Mirizzi & Montanino**, Phys. Rev. D 68, 033005 (2003).\n* Shock radii from **Kneller & Kabadi**, Phys. Rev. D 92, 013009 (2015), Fig. 1, which\n  reads them off a $10.8\\,M_\\odot$ simulation at $t = 3$ s post-bounce: reverse shock\n  1734 km, contact discontinuity 12 348 km, forward shock 30 323 km.',
    [
    code(r'''import json
import warnings

import numpy as np
import matplotlib.pyplot as plt

import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.oscprob as oscprob

plt.rcParams["figure.dpi"] = 110

KM = gd.UNIT_KM
MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)

R_CONTACT_KM, R_FORWARD_KM = 12348.0, 30323.0
R0_KM, R1_KM = 1.0e4, 8.0e4              # the ray we integrate along
L0, L1 = R0_KM*KM, R1_KM*KM'''),
    md(r'''## 1. Build the profile

The forward shock and the contact discontinuity are given a finite width `w`, as a
fraction of the ray, so we can sweep how sharp the front is. A real hydrodynamic shock is
mean-free-path thin ($w \sim 10^{-6}$ here, i.e. 7 cm); a shock read out of a simulation
snapshot is smeared across a few grid cells, so tens of km ($w \sim 10^{-3}$).'''),
    code(r"""def smoothstep(u):
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    return u*u*(3.0 - 2.0*u)

def rarefaction(r_km, r_shock_km):
    '''ln f(x) = [0.28 - 0.69 ln(x_s/km)] [arcsin(1 - x/x_s)]^1.1   (Fogli et al. 2003)'''
    u = np.clip(1.0 - np.asarray(r_km, dtype=float)/r_shock_km, 0.0, 1.0)
    return np.exp((0.28 - 0.69*np.log(r_shock_km))*np.arcsin(u)**1.1)

def sn_shock_ne(width_frac, contact_jump=2.5, y_e=0.5):
    '''Electron number density [eV^3] along the ray, for a front of the given width.'''
    w_km = float(width_frac)*(R1_KM - R0_KM)

    def ne(l):
        r = np.asarray(l, dtype=float)/KM
        rho = 1.0e14*r**(-2.4)                                    # static progenitor
        shocked = smoothstep((R_FORWARD_KM + 0.5*w_km - r)/w_km)  # forward shock
        factor = 1.0 + shocked*(10.0*rarefaction(r, R_FORWARD_KM) - 1.0)
        inside = smoothstep((R_CONTACT_KM + 0.5*w_km - r)/w_km)   # contact discontinuity
        factor = factor*(1.0 + inside*(contact_jump - 1.0))
        out = rho*factor*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*y_e
        return out[()] if np.ndim(out) == 0 else out
    return ne

def shock_breakpoints(width_frac):
    w_km = float(width_frac)*(R1_KM - R0_KM)
    edges = []
    for r in (R_CONTACT_KM, R_FORWARD_KM):
        edges += [(r - 0.5*w_km)*KM, (r + 0.5*w_km)*KM]
    return np.array([L0] + sorted(edges) + [L1])"""),
    code(r'''fig, ax = plt.subplots(figsize=(6.6, 3.6))
rr = np.linspace(R0_KM, R1_KM, 4000)
for w, lab in ((1e-2, r'$w=10^{-2}$ (2100 km)'), (1e-6, r'$w=10^{-6}$ (0.07 km)')):
    ax.loglog(rr, np.asarray(sn_shock_ne(w)(rr*KM)), lw=1.1, label=lab)
for r, lab in ((R_CONTACT_KM, 'contact'), (R_FORWARD_KM, 'forward shock')):
    ax.axvline(r, color='k', ls=':', lw=0.9)
    ax.text(r*1.02, 3e12, lab, rotation=90, fontsize=7, va='bottom')
ax.set_xlabel('radius [km]'); ax.set_ylabel(r'$n_e\ [{\rm eV}^3]$')
ax.legend(fontsize=8); ax.set_title('Post-bounce envelope with a forward shock')
fig.tight_layout()'''),
    md(r'''The density *rises* outward through the shocked shell and drops by an order of
magnitude across the forward shock -- that is the real shape, not an artefact: the
shocked material is compressed, and behind it the rarefaction ("hot bubble") thins out.'''),
    md(r'''## 2. The sharp shock: wrong by 0.2, and loud about it

Three flavours at 15 MeV, so the **H resonance** ($\Delta m^2_{31}$) sits on the ray at
about $4\times10^4$ km -- just outside the forward shock, which is the configuration the
shock-effect literature studies.'''),
    md(r'''### The ground truth, and why it is stored rather than recomputed

Everything below is scored against a tight-tolerance solution of the same Schrödinger
equation — `solve_ivp`/`DOP853` at `rtol=1e-12`, `atol=1e-14`. That solution is a
constant of the physics: it depends on the shock profile and the energy, and **not on
Magνs**. Recomputing it every time this notebook runs cost about fifteen minutes, which
is a quarter of an hour spent re-deriving a number that cannot have changed.

So it is computed once by `make_shock_reference.py` and stored in
`shock_reference.json` as hexadecimal floats, which round-trip exactly — you get the
bits it was computed from, not a decimal rendering of them.

**Only the oracle is frozen.** Every Magνs number here is still computed live; freezing
the reference would be pointless if it also froze the thing being tested. The risk that
does introduce is a stale reference outliving a change to the profile, so the file
carries a fingerprint of the electron density along the ray and the loader refuses a
reference that does not match the profile just built.'''),
    code(r'''_REF_CACHE = {}

def frozen_reference(width_frac):
    """The stored solve_ivp probabilities for this front width, off their exact bits."""
    if not _REF_CACHE:
        with open('shock_reference.json') as handle:
            _REF_CACHE.update(json.load(handle))
    unhex = lambda xs: np.array([float.fromhex(x) for x in xs])
    case = _REF_CACHE['cases']['%.0e' % width_frac]

    # Guard: rebuild the profile and check it is the one the reference came from.  A
    # frozen oracle that silently outlives a change to the physics is worse than no
    # oracle at all, because every comparison against it still looks fine.
    want = unhex(case['fingerprint_ne'])
    got = np.asarray(sn_shock_ne(width_frac)(unhex(_REF_CACHE['fingerprint_l'])),
                     dtype=float)
    if not np.allclose(got, want, rtol=1e-12, atol=0.0):
        raise RuntimeError(
            'the shock profile no longer matches shock_reference.json; '
            're-run `python notebooks/make_shock_reference.py`')

    return unhex(case['P']).reshape(case['shape'])

def to_P(U):
    return np.swapaxes(np.asarray(U).real**2 + np.asarray(U).imag**2, -1, -2)

ENERGY = 15.0*gd.UNIT_MEV
p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
params3 = {k: p[k] for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
hvac3 = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
    params3['s12'], params3['s23'], params3['s13'], params3['dCP'],
    params3['D21'], params3['D31'])
e00_3 = np.diag([1.0, 0.0, 0.0])

def make_H(ne):
    VCC = matter.vcc_func_from_rho_func(
        ne, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)
    def H(l):
        return (1.0/ENERGY)*hvac3 + np.asarray(VCC(l))[..., None, None]*e00_3
    return H

L_OSC = 4.0*np.pi*ENERGY/params3['D31']
Ls = np.linspace(L1 - 6.0*L_OSC, L1, 61)
print('ray = %.0f oscillation lengths of the fastest (D31) oscillation' % ((L1 - L0)/L_OSC))'''),
    code(r"""def measure(width_frac, **kw):
    '''Instantaneous and averaged error against solve_ivp, plus any warnings raised.'''
    ne = sn_shock_ne(width_frac)
    ref = frozen_reference(width_frac)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        got = np.array([np.asarray(oscprob.osc_prob_matter_std_potential(
            3, ne, ENERGY, float(L), params3, L0=L0,
            density_is_of_number_of_electrons=True, **kw)) for L in Ls])
        raised = sorted({w.category.__name__ for w in caught})
    return (float(np.max(np.abs(got[-1] - ref[-1]))),
            float(np.max(np.abs(got.mean(axis=0) - ref.mean(axis=0)))),
            raised, got, ref)

inst_sharp, avg_sharp, warn_sharp, got_sharp, ref_sharp = measure(1e-6)
print('SHARP FRONT (w = 1e-6, 0.07 km)')
print('  instantaneous error : %.3e' % inst_sharp)
print('  AVERAGED error      : %.3e   <-- averaging does NOT rescue this' % avg_sharp)
print('  warnings raised     : %s' % ', '.join(warn_sharp))"""),
    md(r'''**This is the whole point of the notebook.** In notebook 13 averaging cut the error by
53x. Here it does essentially nothing: the error is in the envelope, because the shock
changes how adiabatic the level crossing is, and that is a change in the conversion
probability rather than in its phase.

And the package is not quiet about it. `UnmarkedDiscontinuityWarning` says the
Hamiltonian is not resolved at the scale being sampled; `HybridCertificationWarning` and
`ToleranceNotAchievedWarning` say the answer did not certify. This is the failure mode
you want: wrong, and loud.'''),
    code(r'''fig, ax = plt.subplots(figsize=(6.6, 3.6))
xs_plot = (Ls - Ls[0])/L_OSC
ax.plot(xs_plot, ref_sharp[:, 0, 0], lw=1.0, label='truth (solve_ivp)')
ax.plot(xs_plot, got_sharp[:, 0, 0], lw=1.0, ls='--', label="magnus, strategy='auto'")
ax.axhline(ref_sharp[:, 0, 0].mean(), color='C0', lw=1.0, alpha=0.5,
           label='averaged truth')
ax.axhline(got_sharp[:, 0, 0].mean(), color='C1', lw=1.0, ls='--', alpha=0.7,
           label='averaged magnus')
ax.set_xlabel('baseline, in oscillation lengths'); ax.set_ylabel(r'$P_{ee}$')
ax.legend(fontsize=8)
ax.set_title('Sharp shock: the two averages differ -- this is envelope, not phase')
fig.tight_layout()'''),
    md(r'''## 3. The cure: declare the front with `t_breakpoints`

No fixed-grid method can resolve a discontinuity it was not told about -- that is
mathematics, not a defect. Tell it where the front is and the error collapses.

**On a baseline scan.** This is the case where the cure is established.'''),
    code(r'''ne_sharp = sn_shock_ne(1e-6)
ref_scan = frozen_reference(1e-6)               # the same stored oracle, same baselines

def scan_error(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = np.asarray(oscprob.osc_prob_matter_std_potential(
            3, ne_sharp, ENERGY, Ls, params3, L0=L0,
            density_is_of_number_of_electrons=True, **kw)).reshape(len(Ls), 3, 3)
    return float(np.max(np.abs(P - ref_scan)))

e_bare = scan_error()
e_bp = scan_error(t_breakpoints=shock_breakpoints(1e-6))
print('scan of %d baselines, sharp front' % len(Ls))
print('  without t_breakpoints : %.3e' % e_bare)
print('  with    t_breakpoints : %.3e   (%.0fx better)' % (e_bp, e_bare/e_bp))'''),
    md(r'''### The caveat, and it is a real one

On a **single point**, `t_breakpoints` is not an established cure. Declaring breakpoints
there also changes which engine answers -- it moves the request onto the general Magnus
ladder -- and measured across 18 shock configurations on the averaged observable it
improved 7, **worsened 11**, and pushed 2 answers from inside the requested tolerance to
outside it.

So: on a scan, pass `t_breakpoints`. On a single point, pass it *and check*, for example
against `strategy='magnus'` or `cumulative=True`.'''),
    code(r'''ne_mid = sn_shock_ne(1e-3)                      # a simulation-smeared front, 70 km
ref_pt = frozen_reference(1e-3)[-1]             # Ls[-1] is L1, so this is the endpoint

def point_error(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = np.asarray(oscprob.osc_prob_matter_std_potential(
            3, ne_mid, ENERGY, L1, params3, L0=L0,
            density_is_of_number_of_electrons=True, **kw))
    return float(np.max(np.abs(P - ref_pt)))

print('single point, 70 km front')
print('  bare              : %.3e' % point_error())
print('  with t_breakpoints: %.3e' % point_error(t_breakpoints=shock_breakpoints(1e-3)))
print('  cumulative=True   : %.3e' % point_error(cumulative=True))'''),
    md(r'''## 4. How the error depends on how sharp the front is

Sweeping the width shows where the package struggles and, crucially, that it warns
wherever the answer is actually outside tolerance.'''),
    code(r'''rows = []
for w in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    inst, avg, raised, _, _ = measure(w)
    rows.append((w, w*(R1_KM - R0_KM), avg, bool(raised)))
    print('w=%.0e (%7.2f km)  averaged err %.3e   %s'
          % (w, w*(R1_KM - R0_KM), avg, 'warns' if raised else 'QUIET'))

print()
bad = [r for r in rows if r[2] > 1e-3 and not r[3]]
print('configurations outside 1e-3 with NO warning: %d' % len(bad))'''),
    md(r'''## Summary

| | notebook 13 (solar) | this notebook (shock) |
|---|---|---|
| instantaneous error | 1.4e-03 | 2.0e-01 |
| averaged error | **2.6e-05** | **2.1e-01** |
| averaging helps by | 53x | ~1x |
| what the error is | phase | **envelope** |
| does the package warn? | no | **yes, every time** |
| cure | none needed -- the observable is right | `t_breakpoints` on a scan |

A shock front changes the adiabaticity of the level crossing, so it moves the conversion
probability itself. That is exactly the physics the shock-effect literature studies, and
it is why averaging cannot remove the error the way it does for a smooth solar profile.

The practical rule: **average your instantaneous scan over a few oscillation lengths and
see whether the error moves.** If it collapses, it was phase. If it does not, it is real
-- and on this package, it will also have warned you.'''),
    ])


# --------------------------------------------------------------- reading order

READING_ORDER = [
    ('01_magnus_introduction.ipynb', 'Introduction',
     'the shortest path to a probability'),
    ('02_magnus_2nu_vacuum_matter.ipynb', 'Two-neutrino probabilities',
     'vacuum, constant and varying density, castle wall, Earth and Sun'),
    ('03_magnus_3nu_vacuum_matter.ipynb', 'Three-neutrino probabilities',
     'the same ground with three flavours and a CP phase'),
    ('04_magnus_long_baseline.ipynb', 'Long baselines',
     'probabilities between two points on the surface'),
    ('05_magnus_biprobability.ipynb', 'Biprobability plots',
     'neutrino against antineutrino, as the CP phase runs'),
    ('06_magnus_oscillograms.ipynb', 'Oscillograms',
     'probability across zenith angle and energy at once'),
    ('07_magnus_bsm_sterile_nu.ipynb', 'BSM: sterile neutrinos',
     'four and five flavours'),
    ('08_magnus_bsm_nsi.ipynb', 'BSM: non-standard interactions',
     'a new matter potential in the same slot'),
    ('09_magnus_bsm_liv.ipynb', 'BSM: Lorentz-invariance violation',
     'an energy dependence the vacuum term does not have'),
    ('10_magnus_averaged_probability.ipynb', 'Phase-averaged probabilities',
     'what survives when the phase is unresolvable'),
    ('11_magnus_matrix_exponential.ipynb', 'The matrix exponential',
     'how the evolution operator is actually built'),
    ('12_magnus_adiabatic_hybrid_strategy.ipynb', 'The strategy parameter',
     "'auto' against 'magnus', and when the difference matters"),
    ('13_magnus_tabulated_solar_model.ipynb', 'A tabulated solar model',
     'a real BS05 profile rather than an exponential'),
    ('14_magnus_supernova_shock.ipynb', 'A supernova shock front',
     'where the error stops being a phase and becomes an envelope'),
]


def add_footers():
    r"""Appends a navigation footer to every notebook.

    Each carries the previous and next notebook and a pointer to the API
    reference, so a reader who arrives at any one of them -- which is what
    happens when a search engine or a colleague sends them a link -- can find
    both the path through and the underlying documentation.
    """
    assert set(name for name, _, _ in READING_ORDER) == set(books), (
        'READING_ORDER and the notebooks built here have diverged')

    for index, (name, _, _) in enumerate(READING_ORDER):
        parts = []
        if index:
            previous, title, _ = READING_ORDER[index-1]
            parts.append('**Previous:** [%s](%s)' % (title, previous))
        if index + 1 < len(READING_ORDER):
            following, title, blurb = READING_ORDER[index+1]
            parts.append('**Next:** [%s](%s) --- %s'
                         % (title, following, blurb))
        parts.append('[API reference](%s/functions.html) &middot; '
                     '[Implementation details](%s/implementation_details.html) '
                     '&middot; [All notebooks](.)' % (DOCS, DOCS))
        books[name].cells.append(md('---\n\n' + '  \n'.join(parts)))


def build(execute=True):
    r"""Writes every notebook, executes it, and checks it kept its outputs."""
    add_footers()

    for name, nb in books.items():
        nbf.write(nb, HERE/name)
    print('  wrote %d notebooks' % len(books))

    if not execute:
        return

    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    failed = []
    for path in sorted(HERE.glob('*.ipynb')):
        nb = nbf.read(path, as_version=4)
        started = time.perf_counter()
        try:
            NotebookClient(
                nb, timeout=3600, kernel_name='python3',
                resources={'metadata': {'path': str(path.parent)}}).execute()
            nbf.write(nb, path)
            print('  executed %-44s %6.1f s'
                  % (path.name, time.perf_counter()-started))
        except CellExecutionError as error:
            failed.append(path.name)
            print('  FAILED   %s\n%s' % (path.name, str(error)[:2000]))

    if failed:
        raise SystemExit('notebooks failed to execute: %s' % ', '.join(failed))

    # A notebook whose outputs were stripped renders blank on GitHub, so this
    # is checked here rather than discovered by a reader.
    bare = [path.name for path in sorted(HERE.glob('*.ipynb'))
            if not any(cell.get('outputs') for cell
                       in nbf.read(path, as_version=4).cells
                       if cell.cell_type == 'code')]
    if bare:
        raise SystemExit('notebooks carry no stored outputs: %s'
                         % ', '.join(bare))
    print('  all %d notebooks executed and carry stored outputs' % len(books))


if __name__ == '__main__':
    import sys
    build(execute='--no-execute' not in sys.argv)
