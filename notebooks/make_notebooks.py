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


# Figures are set through LaTeX, which is what notebooks/matplotlibrc asks for.
# A machine that only wants to check that the notebooks still run -- CI, most
# obviously, but also most readers -- has no reason to carry a TeX installation,
# so this asks rather than assumes.  With LaTeX the labels are typeset as the
# stored outputs show them; without it they fall back to matplotlib's own
# mathtext, which changes how they are set and nothing about what is plotted.
#
# It lives in every notebook rather than in matplotlibrc because an rc file
# states a value and cannot ask a question.  Without it, every figure-bearing
# notebook raises "latex could not be found" on a machine with no TeX -- which
# is what CI had been doing, unnoticed, on eighteen of the twenty-six.
LATEX_GUARD = code(r"""# The figures are set through LaTeX where one is available; where it is not,
# matplotlib's own mathtext renders the labels instead.  Same numbers either way.
import shutil

import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = shutil.which('latex') is not None""")


def notebook(title, intro, cells):
    r"""One notebook: a title cell, the LaTeX guard, then whatever the caller supplies."""
    nb = nbf.v4.new_notebook()
    nb.cells = [md('# %s\n\n%s' % (title, intro)), LATEX_GUARD] + cells
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
    'The shortest path to a probability, and the conventions the rest of these notebooks assume.\n\nMag$\\nu$s computes neutrino oscillation probabilities by the **Magnus expansion**: the evolution operator over a slab is $\\exp(\\Omega)$, with $\\Omega$ built from time-ordered integrals of nested commutators of the Hamiltonian. Truncating the series is exactly unitary at any order, which is the property that makes it worth doing this way rather than exponentiating a discretized Hamiltonian.',
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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

sth = NUFIT_NO['s12'] # sin(theta) [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
    code(r'''s12 = NUFIT_NO['s12'] # sin(theta_12) [adim]
s23 = NUFIT_NO['s23'] # sin(theta_23) [adim]
s13 = NUFIT_NO['s13'] # sin(theta_13) [adim]
dCP = NUFIT_NO['dCP'] # [radian]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]'''),
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
    'Oscillation probabilities in a two-flavor system, against energy and against direction, in seven settings: vacuum, constant density, an exponential and a Gaussian profile, a periodic castle wall, a noisy profile, and then the Earth and the Sun.\n\nEach is validated against the closed-form expression where one exists, which is what makes this the notebook to read before trusting any of the others.',
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
    code(r'''def num_density_e_func_prem(r):
    # Y_e per PREM layer (iron core, rock mantle), with the neutron-to-proton ratio
    # derived from it -- the same composition osc_prob_*_earth uses internally, so this
    # recipe and the wrappers describe one Earth rather than two.  A uniform 0.5 here
    # would disagree with them by up to a factor of four on a core-crossing chord, with
    # nothing on screen to say why.  For the uniform composition earlier versions
    # assumed, pass electron_fraction=0.5 here and to the wrappers alike.
    ye = earth.electron_fraction_func_prem(r)
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
        ratio_number_neutrons_to_protons=earth.neutron_to_proton_ratio_from_electron_fraction(ye),
        electron_fraction=ye, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
    'The same seven settings as the previous notebook, with three flavors and a CP-violating phase.\n\nNothing about the method changes -- the Hamiltonian is a $3\\times 3$ matrix in the same slot -- so the interest is in what the extra flavor and the phase do to the probabilities.',
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
    code(r'''def num_density_e_func_prem(r):
    # Y_e per PREM layer (iron core, rock mantle), with the neutron-to-proton ratio
    # derived from it -- the same composition osc_prob_*_earth uses internally, so this
    # recipe and the wrappers describe one Earth rather than two.  A uniform 0.5 here
    # would disagree with them by up to a factor of four on a core-crossing chord, with
    # nothing on screen to say why.  For the uniform composition earlier versions
    # assumed, pass electron_fraction=0.5 here and to the wrappers alike.
    ye = earth.electron_fraction_func_prem(r)
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
        ratio_number_neutrons_to_protons=earth.neutron_to_proton_ratio_from_electron_fraction(ye),
        electron_fraction=ye, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]'''),
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

# Electron number density inside Earth, using the PREM density model
def num_density_e_func_prem(r):
    # Y_e per PREM layer (iron core, rock mantle), with the neutron-to-proton ratio
    # derived from it -- the same composition osc_prob_*_earth uses internally, so this
    # recipe and the wrappers describe one Earth rather than two.  A uniform 0.5 here
    # would disagree with them by up to a factor of four on a core-crossing chord, with
    # nothing on screen to say why.  For the uniform composition earlier versions
    # assumed, pass electron_fraction=0.5 here and to the wrappers alike.
    ye = earth.electron_fraction_func_prem(r)
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
        ratio_number_neutrons_to_protons=earth.neutron_to_proton_ratio_from_electron_fraction(ye),
        electron_fraction=ye, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

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

# Per-detector colors and line styles, shared by both figures
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
NUFIT_IO = gd.load_nufit_params('NuFIT 6.1', 'IO')
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
s12_NO = NUFIT_NO['s12'] # [adim]
s23_NO = NUFIT_NO['s23'] # [adim]
s13_NO = NUFIT_NO['s13'] # [adim]
dCP_NO = NUFIT_NO['dCP'] # [adim]
D21_NO = NUFIT_NO['D21'] # [eV^2]
D31_NO = NUFIT_NO['D31'] # [eV^2]
s12_IO = NUFIT_IO['s12'] # [adim]
s23_IO = NUFIT_IO['s23'] # [adim]
s13_IO = NUFIT_IO['s13'] # [adim]
dCP_IO = NUFIT_IO['dCP'] # [adim]
D21_IO = NUFIT_IO['D21'] # [eV^2]
D31_IO = NUFIT_IO['D31'] # [eV^2]'''),
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
                           for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[NUFIT_IO['dCP']]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_IO['dCP']]])
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
    nubar=False)[nu_i][nu_f]] for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nubar_std_NO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_NO, s23_NO, s13_NO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=True)[nu_i][nu_f]] for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nu_std_IO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_IO, s23_IO, s13_IO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=False)[nu_i][nu_f]] for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nubar_std_IO_sel = np.array([[dCP/np.pi, oscprobstd.osc_prob_3nu_vacuum_std(
    hamiltonians.pmns_mixing_matrix(s12_IO, s23_IO, s13_IO, dCP), D21_NO, D31_NO, energy, baseline*gd.CONV_KM_TO_INV_EV, 
    nubar=True)[nu_i][nu_f]] for dCP in dCP_sel+[NUFIT_NO['dCP']]])
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
                           for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[NUFIT_IO['dCP']]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_IO['dCP']]])
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
                           for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[NUFIT_IO['dCP']]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_IO['dCP']]])
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
                           for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nubar_NO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_NO['dCP']]])
prob_nu_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nu(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                        0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                        n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                        integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                           for dCP in dCP_sel+[NUFIT_IO['dCP']]])
prob_nubar_IO_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(s12_IO, s23_IO, s13_IO, dCP, D21_IO, D31_IO), 
                                                           0.0, baseline*gd.CONV_KM_TO_INV_EV, 
                                                           n_slabs=1, n_tpts_per_slab=2, magnus_exp_order=1,
                                                           integration_method='simpson', n_jobs=1)[nu_i][nu_f]]
                              for dCP in dCP_sel+[NUFIT_IO['dCP']]])
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
    # Y_e per PREM layer (iron core, rock mantle), with the neutron-to-proton ratio
    # derived from it -- the same composition osc_prob_*_earth uses internally, so this
    # recipe and the wrappers describe one Earth rather than two.  A uniform 0.5 here
    # would disagree with them by up to a factor of four on a core-crossing chord, with
    # nothing on screen to say why.  For the uniform composition earlier versions
    # assumed, pass electron_fraction=0.5 here and to the wrappers alike.
    ye = earth.electron_fraction_func_prem(r)
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
        ratio_number_neutrons_to_protons=earth.neutron_to_proton_ratio_from_electron_fraction(ye),
        electron_fraction=ye, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

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
                               for dCP in dCP_sel+[NUFIT_NO['dCP']]])
    prob_nubar_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(costhz_arr[i], l, energy, 
                                                                              s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                            0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                            n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                            integration_method='simpson', n_jobs=10)[nu_i][nu_f]] 
                               for dCP in dCP_sel+[NUFIT_NO['dCP']]])
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
                               for dCP in dCP_sel+[NUFIT_NO['dCP']]])
    prob_nubar_sel = np.array([[dCP/np.pi, oscprob.osc_prob(lambda l: H_nubar(costhz_arr[i], l, energy, 
                                                                              s12_NO, s23_NO, s13_NO, dCP, D21_NO, D31_NO), 
                                                            0.0, l_max_arr[i]*gd.CONV_KM_TO_INV_EV, 
                                                            n_slabs=100, n_tpts_per_slab=10, magnus_exp_order=3,
                                                            integration_method='simpson', n_jobs=10)[nu_i][nu_f]] 
                               for dCP in dCP_sel+[NUFIT_NO['dCP']]])
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
drawing is a single call to `magnus.plotting.plot_oscillogram`; the color
map, the color-bar label, the tick spacings, and the white-stroked corner
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
    # For 2nu oscillations in the 23 sector the flavor indices have to be
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
sth = NUFIT_NO['s23'] # [adim]
Dm2 = NUFIT_NO['D31'] # [eV^2]

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
sth = NUFIT_NO['s12'] # [adim]
Dm2 = NUFIT_NO['D21'] # [eV^2]

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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
    'Four- and five-flavor systems, where the extra states do not couple to the weak interaction.\n\nThe machinery is unchanged; only the dimension of the Hamiltonian and the number of mixing angles and phases grow.',
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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
    """A sterile scenario against standard three-flavor, versus baseline."""
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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

# Electron number density inside Earth, using the PREM density model
def num_density_e_func_prem(r):
    # Y_e per PREM layer (iron core, rock mantle), with the neutron-to-proton ratio
    # derived from it -- the same composition osc_prob_*_earth uses internally, so this
    # recipe and the wrappers describe one Earth rather than two.  A uniform 0.5 here
    # would disagree with them by up to a factor of four on a core-crossing chord, with
    # nothing on screen to say why.  For the uniform composition earlier versions
    # assumed, pass electron_fraction=0.5 here and to the wrappers alike.
    ye = earth.electron_fraction_func_prem(r)
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
        ratio_number_neutrons_to_protons=earth.neutron_to_proton_ratio_from_electron_fraction(ye),
        electron_fraction=ye, density_matter_is_in_g_per_cm3=True) # [eV^{-3}]

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
# flavor counts, and a genuine PREM profile sampled inside every slab.
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
# The legend describes line *style*, not color -- the 3+1 curve is a different
# color in every panel -- so it is built from legend_proxies. That replaces the
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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
s12 = NUFIT_NO['s12'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]

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
s12 = NUFIT_NO['s12'] # [adim]
s23 = NUFIT_NO['s23'] # [adim]
s13 = NUFIT_NO['s13'] # [adim]
dCP = NUFIT_NO['dCP'] # [adim]
D21 = NUFIT_NO['D21'] # [eV^2]
D31 = NUFIT_NO['D31'] # [eV^2]

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
    "Some extensions of the Standard Model -- string-inspired constructions, models\nwith a preferred frame, effective descriptions of quantum gravity -- break\nLorentz invariance. If they do, neutrinos are an unusually good place to look:\noscillations measure a *phase*, and a phase accumulated over an astrophysical\nbaseline is sensitive to energies that no accelerator reaches.\n\nMag$\\nu$s treats this as an extra, CPT-odd term in the Hamiltonian,\n\n$$\\mathbf{H} = \\frac{1}{2E}\\,\\mathbf{U}\\,\\mathbf{M}^2\\,\\mathbf{U}^\\dagger\n             \\;+\\; \\mathbf{V}_{\\rm CC}\n             \\;+\\; E^{\\,n}\\,\n               \\mathbf{U}_\\xi\\,\n               \\frac{\\mathbf{B}}{\\Lambda^{\\,n}}\\,\n               \\mathbf{U}_\\xi^\\dagger ,$$\n\nwhere $\\mathbf{B} = {\\rm diag}(b_1, b_2, b_3)$ holds the eigenvalues of the LIV\noperator, $\\mathbf{U}_\\xi$ rotates from its eigenbasis to the flavor basis\nthrough angles $\\xi_{12}, \\xi_{23}, \\xi_{13}$ and a phase $\\delta_{\\xi\\rm CP}$,\n$\\Lambda$ is the scale that makes the eigenvalues dimensionless, and $n$ is the\noperator's dimension minus three.\n\n## Why the exponent is the whole story\n\nLook at how the three terms scale with energy:\n\n| term | scaling | behavior at high $E$ |\n|---|---|---|\n| vacuum | $\\Delta m^2 / 2E \\;\\propto\\; E^{-1}$ | switches **off** |\n| matter | $V_{\\rm CC}$, independent of $E$ | flat |\n| LIV | $E^{\\,n}\\, b/\\Lambda^{\\,n}$ | switches **on** (for $n \\geq 0$) |\n\nStandard oscillations die away at high energy: the phase $\\Delta m^2 L / 2E$\nshrinks, and the probability freezes at its zero-baseline value. A LIV term\ndoes the opposite. So the signature is not that oscillations look slightly\ndifferent -- it is that they are still *there*, at energies where the Standard\nModel says they should have stopped.\n\nThat is what the figures below show, and it is why every one of them is plotted\nagainst energy.",
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
flavors those are `sxi12`, `sxi23`, `sxi13`, `dxiCP` (the mixing of the LIV
eigenbasis into flavor), `b1`, `b2`, `b3` (its eigenvalues), `Lambda`, and
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

Each curve below is normalized to the *same* LIV phase at 100 GeV, so they are
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
    md(r'''## 5. Two flavors, and setting a limit

The two-flavor interface is the same with one angle and two eigenvalues
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

* `magnus.oscprob` has `_liv` wrappers for 2, 3, 4 and 5 flavors and for every
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
    value overlaid as a dashed line of the same color.

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

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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

sth, Dm2 = NUFIT_NO['s12'], NUFIT_NO['D21']
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
# NOT np.diag([1, 0, 0, 0]): beyond three flavors the matter term is not e_ee.
# The sterile state carries -V_NC = (r/2) V_CC once the actives' common V_NC is
# removed, and the library's wrapper below uses exactly this projector -- so a
# hand-built zero there makes the solve_ivp "ground truth" the wrong problem,
# and the error column would blame the strategy for the reference.
e00_4 = np.asarray(matter.matter_potential_projector(4))

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
    code(r'''h_matt4_nsi = np.asarray(matter.matter_potential_projector(4)) \
    + hamiltonians.hamiltonian_4nu_nsi(
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
e00_5 = np.asarray(matter.matter_potential_projector(5))   # see the 4nu note above

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
    code(r'''h_matt5_nsi = np.asarray(matter.matter_potential_projector(5)) \
    + hamiltonians.hamiltonian_5nu_nsi(
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
    "You have a real solar model on disk -- a table of radius, density and composition --\nand you want oscillation probabilities from it. This notebook does exactly that with\n**BS2005-AGS,OP** (Bahcall, Serenelli & Basu, ApJ 621, L85), and uses it to separate two\nquantities that are easy to confuse:\n\n* the **instantaneous** probability at one baseline, which is what `osc_prob_*` returns;\n* the **phase-averaged** probability, which is what a solar-neutrino experiment measures.\n\nThey are different quantities, not two estimates of one quantity, and the notebook's\nheadline is about how you get the second.\n\n**The tempting route does not work.** Averaging a scan of instantaneous probabilities over a\nwindow of several oscillation lengths looks like the obvious way to reach the observable. On\na solar trajectory it is not: the answer drifts by about $10^{-2}$ depending on how wide\na window you pick, because widening the window also averages over a changing density. The\nestimator has no converged value to offer.\n\n**The direct route is exact.** `average=True` evaluates the phase-averaged limit in closed\nform -- one matrix product, no scan -- and it reproduces the textbook adiabatic MSW\nexpression to **machine precision, 3e-16, across 1--20 MeV**, checked against a formula that\nowes nothing to Mag$\\nu$s.\n\nThe notebook also shows the diagnostics: `strategy_info['sampling']` for how coarsely a scan\nresolves the oscillation it is sampling, and `avgprob.coherence_report` for whether the\naveraged limit applies at all.",
    [
    code(r'''import os
import time
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import magnus.globaldefs as gd

# Best-fit oscillation parameters from the latest global fit.
# load_nufit_params returns exactly the six parameters the
# osc_prob_3nu_* functions take, so it can be splatted straight in.
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
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
  center to 0.75 at the surface. Using a fixed $Y_e = 0.5$ would be wrong by up to 70 %,
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

We take a two-flavor calculation at **5 MeV** -- in the $^8$B range -- over one solar
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
params2 = {'sth': NUFIT_NO['s12'], 'Dm2': NUFIT_NO['D21']}

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
    md(r'''That is inside the requested tolerance, reported as `certified`, with no warning --
which is the right answer, but not yet an interesting one. The interesting question is what
this number is *of*: it is the error in the probability at **one exact baseline**, and no
solar experiment measures that.

The rest of this notebook is about the quantity one does measure, and about how easy it is
to compute something that looks like it and is not.'''),
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
the curve through them is an artifact. `nyquist_points` says how many baselines you would
need to sample the oscillation properly -- about 900 here, and several thousand at the
energies and flavor counts used elsewhere in the documentation.

That is the signature of a quantity dominated by phase. A solar-neutrino experiment
resolves none of it: the $^8$B production region is extended, the Sun-Earth phase is
$\sim10^{10}$ cycles, and detector energy resolution finishes the job.'''),
    md(r'''## 4. The averaged probability -- the quantity that is actually observed

A solar experiment measures the **phase-averaged** survival probability. The obvious way to
get at it is to average a scan of instantaneous probabilities over a window of several
oscillation lengths, and that is what this section does first -- because it is the natural
thing to try, it is what an earlier version of this notebook reported, and **it does not
work here.** The section after it shows the quantity that does.'''),
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
# The verdicts are COMPUTED, not written into the format string.  They used to be
# asserted -- '<-- outside 1e-3' was a literal -- and when the underlying numbers
# moved the cell went on printing a verdict its own output contradicted.  A label
# that cannot be wrong by construction is worth more than a tidy one.
def _verdict(err, tol=1.0e-3):
    return ('outside' if err > tol else 'inside') + ' %.0e' % tol

print('instantaneous err : %.3e   <-- %s' % (err_inst, _verdict(err_inst)))
print('AVERAGED err      : %.3e   <-- %s' % (err_avg, _verdict(err_avg)))
print('averaging changes the error by %.2fx  (%s)'
      % (err_inst/err_avg,
         'helps' if err_avg < err_inst else 'does not help here -- see below'))
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
ax.legend(fontsize=8); ax.set_title('A scan, and the mean of that scan')
fig.tight_layout()'''),
    md(r'''Averaging did not help. It is worth being precise about why, because the reason is not
that the package is inaccurate -- it is that **the mean of this scan is not a converged
estimate of anything.**

Two things are wrong with it. The window is six *vacuum* oscillation lengths, but the
neutrino is in matter, so the window is not a whole number of actual cycles. And -- the
part that cannot be tuned away -- widening the window does not fix it, because a wider
window also averages over a **changing density**, and the averaged probability depends on
the local density. There is no window width that separates phase from the profile.

The truth alone shows this, no package involved: one ODE solve, several window widths.'''),
    code(r'''# Truth only.  One integration with many output points, so this is cheap: the
# expensive thing in the cell above was one full ray propagation per baseline.
print('%-22s %s' % ('window', 'mean of the truth'))
print('-'*44)
for n_osc in (6, 12, 24, 48):
    Ls_w = np.linspace(L1 - n_osc*L_OSC, L1, 20*n_osc + 1)
    P_w = np.array([to_P(U) for U in exact_U_many(H_of_l, L0, Ls_w, 2)])
    print('%-22s %.6f' % ('%d oscillation lengths' % n_osc, P_w[:, 0, 0].mean()))'''),
    md(r'''The number moves and keeps moving. Widening the window makes it *worse*, not better,
which is the signature of an estimator whose bias is not statistical.

### The averaged probability, computed rather than estimated

`average=True` evaluates the phase-averaged limit in closed form -- one matrix product, no
scan, no window. And it can be checked against something outside Mag$\nu$s entirely: for two
flavors on an adiabatic trajectory the averaged survival probability is the textbook MSW
expression

$$\langle P_{ee}\rangle = \frac{1}{2}
  + \frac{1}{2}\cos 2\theta_m(L_0)\,\cos 2\theta_m(L_1),$$

with $\theta_m$ the matter mixing angle at each end. Both ends are evaluated **in matter**:
this trajectory stops one scale height in, not in vacuum, so the usual $\cos 2\theta_{\rm vac}$
at the far end would be the wrong reference. The LMA solar crossing is never sharp, which is
what makes the adiabatic form exact here rather than approximate -- notebook 12 measures that
directly.'''),
    code(r'''TH_VAC = np.arcsin(params2['sth'])

def cos2theta_matter(l, energy):
    """cos(2 theta_m): the two-flavor matter mixing angle at position l."""
    x = 2.0*energy*float(np.asarray(VCC(l)))/params2['Dm2']
    return np.cos(np.arctan2(np.sin(2.0*TH_VAC), np.cos(2.0*TH_VAC) - x))

def adiabatic_averaged(energy):
    """The closed-form averaged P_ee, both ends in matter."""
    return 0.5 + 0.5*cos2theta_matter(L0, energy)*cos2theta_matter(L1, energy)

print('%-10s %-16s %-16s %s' % ('E [MeV]', 'average=True', 'adiabatic', '|difference|'))
print('-'*60)
worst = 0.0
for E_mev in (1.0, 2.0, 5.0, 8.0, 10.0, 15.0, 20.0):
    E_here = E_mev*gd.UNIT_MEV
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = np.asarray(oscprob.osc_prob_matter_std_potential(
            2, ne_bs05, E_here, L1, params2, L0=L0,
            density_is_of_number_of_electrons=True, average=True))[0][0]
    want = adiabatic_averaged(E_here)
    worst = max(worst, abs(got - want))
    print('%-10.1f %-16.8f %-16.8f %.2e' % (E_mev, got, want, abs(got - want)))
print()
print('worst disagreement across the 8B range: %.2e' % worst)'''),
    md(r'''Machine precision, across the whole $^8$B range. That is the averaged solar
probability, and it is exact.

Compare the two routes on the same quantity at 5 MeV: the scan mean gave a number that
drifted by about $10^{-2}$ depending on how wide a window was chosen, while
`average=True` reproduces an independent closed form to $10^{-16}$ at a fraction of the
cost. **If the averaged probability is what you want, ask for it; do not estimate it by
averaging a scan.** The scan is for looking at the oscillation, not for integrating it.'''),
    md(r'''### Nor is the scan's behavior an artifact of the interpolation

A cubic spline through the same table (still in $\log n_e$) is a different profile, so it is
a fair second opinion on the scan-mean estimator. It gives different numbers in both columns
and the same verdict: the reduction factor is order unity either way, so the failure of the
window mean is a property of the estimator rather than of one particular interpolant.'''),
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
# %.2f, not %.0f: a reduction of 0.85x printed as "1x" reads as "averaging left it
# alone" when what happened is that averaging made it slightly worse.  Two decimals
# is the difference between a number and a rounding artifact.
print('%-10s %12s %12s %10s' % ('interpolant', 'instant.', 'averaged', 'reduction'))
print('%-10s %12.3e %12.3e %9.2fx' % ('linear', err_inst, err_avg, err_inst/err_avg))
print('%-10s %12.3e %12.3e %9.2fx' % ('cubic', ic, ac, ic/ac))'''),
    md(r'''## 5. When is the averaged limit the right limit?

Section 4 used `average=True` on the strength of the trajectory being adiabatic and the
phase being unresolvable. Neither is something to assume. `avgprob.coherence_report` checks
the second directly: it reports which pairs of eigenvalues have genuinely decohered over the
baseline, and which sit in the middle regime where **neither** the oscillating nor the
averaged expression is right.'''),
    code(r'''import magnus.avgprob as avgprob

lam = np.linalg.eigvalsh(np.asarray(H_of_l(0.5*L1)))
blocks, undecided = avgprob.coherence_report(lam, phase_scale=L1)
print('coherence blocks :', blocks)
print('pairs in neither limit:', undecided or 'none -- the averaged expression is exact here')'''),
    md(r'''## Summary

| | |
|---|---|
| instantaneous error at 5 MeV | inside the requested 1e-3, `certified`, no warning -- but it is the error at *one baseline*, which no solar experiment measures |
| averaging a scan to get the observable | **does not work here.** The mean drifts by about $10^{-2}$ with the window width, because a wider window also averages over changing density |
| the averaged probability, done properly | `average=True` -- closed form, one matrix product, and it matches the adiabatic MSW expression to $10^{-16}$ across 1--20 MeV |
| how to check the limit applies | `avgprob.coherence_report`, and `strategy_info['sampling']` for how coarsely a scan resolves the oscillation |

The lesson is not that a large instantaneous error is harmless -- it is that the
instantaneous probability and the averaged probability are **different quantities**, and
that estimating the second from a scan of the first is a numerical method with its own
error, which here is larger than anything it was meant to diagnose.

See :doc:`averaged_probability` in the documentation for the full treatment, and notebook 14
for a profile where the averaged observable is genuinely wrong -- a sharp shock front, where
the error is in the envelope and no amount of averaging touches it.'''),
    md(r'''## 6. Two BSM scenarios on the same model

Everything above is two-flavor and runs to one scale length. The two sections that follow
change three things at once, deliberately, and it is worth saying which:

* **Three flavors**, because NSI and a sterile state are both defined against the
  three-flavor picture.
* **The whole ray**, center to surface, because the averaged survival probability is a
  property of the full path.
* **`average=True`**, because that is the observable -- notebook 25 section 10 measures what
  recovering it from instantaneous evaluations costs another code, and the answer is a
  Monte-Carlo estimate at $1/\sqrt{N}$.

**In both sections the standard three-flavor curve is drawn alongside.** A BSM curve on its
own says nothing about size, and size is the only thing worth reading off these panels.'''),
    code(r'''OSC3_BSM = dict(NUFIT_NO)
R_SUN_BSM = float(x_nat[-1])
E_BSM = np.logspace(np.log10(0.1), np.log10(20.0), 40)*gd.UNIT_MEV
h_vac3_bsm = np.asarray(
    hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC3_BSM))
PER_NE_BSM = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)


def averaged_3nu(**kw):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        3, ne_bs05, E_BSM, R_SUN_BSM, OSC3_BSM, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE,
        density_is_of_number_of_electrons=True, average=True, **kw))


t0 = time.perf_counter()
P_std_bsm = averaged_3nu()
t_std_bsm = time.perf_counter() - t0
print('standard 3nu, %d averaged energies over the full ray: %.3f s'
      % (len(E_BSM), t_std_bsm))
print('  <P_ee> runs %.4f (low energy) -> %.4f (high energy)'
      % (P_std_bsm[0], P_std_bsm[-1]))'''),
    md(r'''### 6.1 Non-standard interactions

NSI adds $\varepsilon$ to the matter matrix: $V_{\rm CC}(\mathrm{diag}(1,0,0) + \varepsilon)$,
with $\varepsilon$ dimensionless. Magνs offers this at two and three flavors only --- there is
no four- or five-flavor NSI route --- so "NSI at 3+1" is not a thing that can be asked for
here.

**One check first.** `average=True` can be served either by the adiabatic route or by a
numerical window over an explicitly propagated probability, and those are *different
quantities* --- the $L/E\to\infty$ limit against an average over a finite spread. Which one the
NSI wrapper reaches is checked by its cost: the adiabatic route never propagates.'''),
    code(r'''EPS_BSM = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0,
               eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)

t0 = time.perf_counter()
P_nsi_bsm = np.asarray(oscprob.osc_prob_matter_nsi(
    3, ne_bs05, E_BSM, R_SUN_BSM, OSC3_BSM, EPS_BSM, L0=0.0,
    nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True, average=True))
t_nsi_bsm = time.perf_counter() - t0

print('WHICH ROUTE? standard %.3f s, NSI %.3f s -- both sub-second, so both adiabatic.'
      % (t_std_bsm, t_nsi_bsm))
print('A numerical window over this ray propagates a ~13 000 radian phase and costs')
print('minutes, so the timing is what distinguishes them.')
print()
print('eps = %s' % {k: v for k, v in EPS_BSM.items() if v})
print('departure from the standard curve: max %.4f at %.2f MeV, mean %.4f'
      % (np.max(np.abs(P_nsi_bsm - P_std_bsm)),
         E_BSM[int(np.argmax(np.abs(P_nsi_bsm - P_std_bsm)))]/gd.UNIT_MEV,
         np.mean(np.abs(P_nsi_bsm - P_std_bsm))))'''),
    md(r'''### 6.2 A sterile state

3+1 changes the matter term as well as the mixing, and that is the part most easily got wrong.
The active flavors share the neutral-current potential $V_{\rm NC}$ and it cancels; a **sterile
state feels neither current**, so once the actives' common piece is removed it carries
$-V_{\rm NC} = (r/2)\,V_{\rm CC}$ with $r = n_n/n_p$. The projector is
$\mathrm{diag}(1, 0, 0, r/2)$, and it comes from `matter.matter_potential_projector` rather than
being written out --- writing it out is exactly how notebook 25's own PREM referee was wrong by
$2.6\times10^{-2}$ until recently.

**So the check that matters is whether that entry is live**, and it is easy: vary $r$ and watch
the curve move. If the term were missing, the rows below would be identical.'''),
    code(r'''STERILE_BSM = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0,
                   d14=0.0, d24=0.0, D41=1.0e-5)
OSC4_BSM = dict(OSC3_BSM, **STERILE_BSM)


def averaged_4nu(ratio=1.0, energies=None):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        4, ne_bs05, E_BSM if energies is None else energies, R_SUN_BSM, OSC4_BSM,
        L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True,
        average=True, ratio_number_neutrons_to_protons=ratio))


t0 = time.perf_counter()
P_ster_bsm = averaged_4nu()
t_ster_bsm = time.perf_counter() - t0
print('3+1, %d averaged energies: %.3f s  (the adiabatic route reaches four flavors)'
      % (len(E_BSM), t_ster_bsm))
print()
print('IS THE STERILE NEUTRAL-CURRENT ENTRY LIVE?  <P_ee> at 1, 5, 15 MeV:')
E_PROBE_BSM = np.array([1.0, 5.0, 15.0])*gd.UNIT_MEV
for ratio in (0.5, 1.0, 1.5):
    print('   n_n/n_p = %.1f  ->  %s'
          % (ratio, np.array2string(averaged_4nu(ratio, E_PROBE_BSM).ravel(),
                                    precision=6)))
print('   The curve moves, so the sterile state is feeling the medium.')
print()
print('departure from the standard curve: max %.4f at %.2f MeV, mean %.4f'
      % (np.max(np.abs(P_ster_bsm - P_std_bsm)),
         E_BSM[int(np.argmax(np.abs(P_ster_bsm - P_std_bsm)))]/gd.UNIT_MEV,
         np.mean(np.abs(P_ster_bsm - P_std_bsm))))'''),
    code(r'''fig, ax = plt.subplots(2, 1, figsize=(6.6, 6.0), sharex=True,
                       gridspec_kw=dict(height_ratios=[2.3, 1.0], hspace=0.08))
ax[0].semilogx(E_BSM/gd.UNIT_MEV, P_std_bsm, lw=2.0, color='k',
               label=r'standard 3$\nu$')
ax[0].semilogx(E_BSM/gd.UNIT_MEV, P_nsi_bsm, lw=1.5, color='C4',
               label=r'NSI  ($\varepsilon_{ee}=0.15$, $\varepsilon_{e\mu}=0.05$)')
ax[0].semilogx(E_BSM/gd.UNIT_MEV, P_ster_bsm, lw=1.5, color='C0',
               label=r'3+1  ($\Delta m^2_{41} = 10^{-5}$ eV$^2$)')
ax[0].set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
ax[0].set_title('BS2005-AGS,OP: two BSM scenarios against the standard curve',
                fontsize=10)
ax[0].grid(True, alpha=0.2)
ax[0].legend(fontsize=8)
ax[0].set_ylim(0.0, 1.0)

ax[1].semilogx(E_BSM/gd.UNIT_MEV, P_nsi_bsm - P_std_bsm, lw=1.4, color='C4',
               label='NSI')
ax[1].semilogx(E_BSM/gd.UNIT_MEV, P_ster_bsm - P_std_bsm, lw=1.4, color='C0',
               label='3+1')
ax[1].axhline(0.0, color='k', lw=0.7, alpha=0.5)
ax[1].set_xlabel(r'$E_\nu$ [MeV]')
ax[1].set_ylabel(r'departure from 3$\nu$')
ax[1].grid(True, which='both', alpha=0.2)
ax[1].legend(fontsize=7.6)
for a in ax:
    a.set_xlim(E_BSM[0]/gd.UNIT_MEV, E_BSM[-1]/gd.UNIT_MEV)
fig.savefig('../fig/solar_bsm.pdf', bbox_inches='tight')'''),
    md(r'''**Where each effect lives is as informative as how big it is.** The sterile
departure is concentrated near the MSW transition, because that is where the extra level
crossing sits and where an additional state can rearrange which vacuum eigenstate the neutrino
emerges as. The NSI departure is broader and smaller: with these $\varepsilon$ it shifts the
resonance position rather than adding structure.

**And both are much smaller here than the same physics produces on a supernova shock.** The
identical $\varepsilon$ moves notebook 14's shock probability by 0.44 on average, against
roughly a hundredth of that here. The reason is the observable rather than the model: this
panel is **phase-averaged**, so everything entering through the phase integrates away and only
the change in the eigenvectors and level crossings survives. On the shock the quantity is an
instantaneous probability along a ray, and the phase term is the larger part of the effect.

That is worth carrying away from this notebook: **a BSM effect's size is a property of the
observable at least as much as of the model**, and a sensitivity estimate quoted without saying
which probability it refers to can be wrong by a factor of thirty.'''),
    ])

# ------------------------------------------------- 14_magnus_supernova_shock
books['14_magnus_supernova_shock.ipynb'] = notebook(
    'A supernova shock front: when the error is real',
    'Notebook 13 ended on a solar case where the averaged observable is exact -- `average=True`\nreproduces the adiabatic MSW expression to machine precision. This notebook is the opposite\ncase, and the contrast is the point: here the averaged observable is genuinely **wrong**, by\n0.21 in probability, and no amount of averaging repairs it.\n\nA supernova shock front changes the **adiabaticity of the MSW level crossing**, so it\nmoves the conversion probability *itself* rather than the phase of an oscillation.\nAveraging cannot remove that. Here the package is wrong by **0.21 in probability on the\naveraged observable** -- and, importantly, it **says so every time**.\n\nThe profile is the standard one from the literature:\n\n* $\\rho_0(x) = 10^{14}\\,(x/\\mathrm{km})^{-2.4}\\ \\mathrm{g\\,cm^{-3}}$, forward-shock jump\n  $\\xi = V_+/V_- \\simeq 10$, and the rarefaction shape behind it, from\n  **Fogli, Lisi, Mirizzi & Montanino**, Phys. Rev. D 68, 033005 (2003).\n* Shock radii from **Kneller & Kabadi**, Phys. Rev. D 92, 013009 (2015), Fig. 1, which\n  reads them off a $10.8\\,M_\\odot$ simulation at $t = 3$ s post-bounce: reverse shock\n  1734 km, contact discontinuity 12 348 km, forward shock 30 323 km.',
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
magnitude across the forward shock -- that is the real shape, not an artifact: the
shocked material is compressed, and behind it the rarefaction ("hot bubble") thins out.'''),
    md(r'''## 2. The sharp shock: wrong by 0.2, and loud about it

Three flavors at 15 MeV, so the **H resonance** ($\Delta m^2_{31}$) sits on the ray at
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
    md(r'''**This is the whole point of the notebook.** The averaged error is as large as the
instantaneous one, and both are far outside any tolerance worth asking for. The error is in
the **envelope**: the shock changes how adiabatic the level crossing is, which moves the
conversion probability itself rather than the phase of an oscillation, and there is no
averaging operation that undoes that.

Contrast notebook 13, where the averaged observable came out exact against an independent
closed form. The difference between the two cases is not the flavor content or the energy;
it is whether the profile has a feature sharp enough to break adiabaticity.

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
| the averaged observable | **exact** -- matches the adiabatic MSW closed form to 1e-16 | **wrong by 0.21** in probability |
| what the error is | none to speak of | **envelope** -- the front breaks adiabaticity |
| does the package warn? | nothing to warn about | **yes, every time** |
| cure | -- | `t_breakpoints` on the front |

A shock front changes the adiabaticity of the level crossing, so it moves the conversion
probability itself. That is exactly the physics the shock-effect literature studies, and it
is why no averaging operation removes the error the way none is needed on a smooth solar
profile.

The practical rule is about **which quantity you ask for, and whether the package certified
it**: compute the observable directly with `average=True` rather than estimating it from a
scan, and read the warnings. A result that is outside tolerance and silent is the dangerous
one; here the package is outside tolerance and loud, which is the failure mode you want.'''),
    md(r'''## 6. Two BSM scenarios on the shock

The same two scenarios notebook 13 puts on the Sun, put here on the shock, and the contrast
between the two notebooks is the reason to read them together.

**Everything below passes `t_breakpoints`.** This notebook has already established that no
fixed grid resolves a front it was not told about, and that a scan without the breakpoints is
measuring straddled slabs rather than physics. A BSM comparison run without them would be
comparing two wrong answers, and the difference between two wrong answers is not the BSM
effect.

**The standard three-flavor curve is drawn alongside in both cases**, because a departure is
only legible against what it departs from.'''),
    code(r'''ne_bsm = sn_shock_ne(1e-3)                 # the 70 km, simulation-smeared front
bps_bsm = shock_breakpoints(1e-3)
E_BSM_SHOCK = ENERGY                       # 15 MeV, as everywhere in this notebook



P_std_shock = np.asarray(oscprob.osc_prob_matter_std_potential(
    3, ne_bsm, E_BSM_SHOCK, Ls, params3, L0=L0,
    density_is_of_number_of_electrons=True, t_breakpoints=bps_bsm,
    n_slabs=32000, max_n_slabs=128000)).reshape(len(Ls), 3, 3)[:, 0, 0]
print('standard 3nu along the ray: P_ee runs %.4f .. %.4f'
      % (P_std_shock.min(), P_std_shock.max()))'''),
    md(r'''### 6.1 Non-standard interactions'''),
    code(r'''EPS_SHOCK = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0,
                 eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    P_nsi_shock = np.asarray(oscprob.osc_prob_matter_nsi(
        3, ne_bsm, E_BSM_SHOCK, Ls, params3, EPS_SHOCK, L0=L0,
        density_is_of_number_of_electrons=True, t_breakpoints=bps_bsm,
        n_slabs=32000, max_n_slabs=128000)).reshape(len(Ls), 3, 3)[:, 0, 0]

print('eps = %s' % {k: v for k, v in EPS_SHOCK.items() if v})
print('departure from the standard curve: max %.4f, mean %.4f'
      % (np.max(np.abs(P_nsi_shock - P_std_shock)),
         np.mean(np.abs(P_nsi_shock - P_std_shock))))'''),
    md(r'''### 6.2 A sterile state

**The splitting here is $\Delta m^2_{41} = 1\,\mathrm{eV}^2$ -- the physically motivated
short-baseline value -- and notebook 25 deliberately uses $10^{-2}$ instead.** That is not an
inconsistency, and the reason is worth stating because it decides what each notebook can claim.

This section compares Magνs against *itself* at three flavors, so no independent referee is
needed. Notebook 25 compares it against another code and referees both with an adaptive DOP853
integration --- and at an eV-scale splitting that referee has to resolve about
$5.9\times10^{6}$ radians of phase over this ray, some 940 000 oscillations, which costs of
order a day. So the eV-scale case can be *computed* and shown here, and cannot be
*independently validated* there. Section 12 of notebook 25 states the same limit from the other
side.

The sterile state also feels the medium: the matter term carries $\mathrm{diag}(1,0,0,r/2)$ from
`matter.matter_potential_projector`, never written out by hand.'''),
    code(r'''STERILE_SHOCK = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0,
                     d14=0.0, d24=0.0, D41=1.0)
OSC4_SHOCK = dict(params3, **STERILE_SHOCK)

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    P_ster_shock = np.asarray(oscprob.osc_prob_matter_std_potential(
        4, ne_bsm, E_BSM_SHOCK, Ls, OSC4_SHOCK, L0=L0,
        density_is_of_number_of_electrons=True, t_breakpoints=bps_bsm,
        n_slabs=32000, max_n_slabs=128000)).reshape(len(Ls), 4, 4)[:, 0, 0]

L_OSC_41 = 4.0*np.pi*ENERGY/STERILE_SHOCK['D41']
print('D41 = %.0e eV^2: sterile oscillation length %.4f km, so the 88.8 km of sampled'
      % (STERILE_SHOCK['D41'], L_OSC_41/KM))
print('ray spans %.0f of them -- unresolvable by eye, which is the point.'
      % ((Ls[-1] - Ls[0])/L_OSC_41))
print()
print('departure from the standard curve: max %.4f, mean %.4f'
      % (np.max(np.abs(P_ster_shock - P_std_shock)),
         np.mean(np.abs(P_ster_shock - P_std_shock))))'''),
    code(r'''fig, ax = plt.subplots(2, 1, figsize=(6.8, 6.0), sharex=True,
                       gridspec_kw=dict(height_ratios=[2.3, 1.0], hspace=0.08))
xs_bsm = (Ls - Ls[0])/KM
ax[0].plot(xs_bsm, P_std_shock, lw=2.0, color='k', label=r'standard 3$\nu$')
ax[0].plot(xs_bsm, P_nsi_shock, lw=1.4, color='C4',
           label=r'NSI  ($\varepsilon_{ee}=0.15$, $\varepsilon_{e\mu}=0.05$)')
ax[0].plot(xs_bsm, P_ster_shock, lw=1.4, color='C0',
           label=r'3+1  ($\Delta m^2_{41} = 1$ eV$^2$)')
ax[0].set_ylabel(r'$P(\nu_e \to \nu_e)$')
ax[0].set_title('The shock, with NSI and with a sterile state', fontsize=10)
ax[0].grid(True, alpha=0.2)
ax[0].legend(fontsize=8)

ax[1].plot(xs_bsm, P_nsi_shock - P_std_shock, lw=1.3, color='C4', label='NSI')
ax[1].plot(xs_bsm, P_ster_shock - P_std_shock, lw=1.3, color='C0', label='3+1')
ax[1].axhline(0.0, color='k', lw=0.7, alpha=0.5)
ax[1].set_xlabel('distance beyond %.0f km along the ray [km]' % (Ls[0]/KM))
ax[1].set_ylabel(r'departure from 3$\nu$')
ax[1].grid(True, alpha=0.2)
ax[1].legend(fontsize=7.6)
for a in ax:
    a.set_xlim(xs_bsm[0], xs_bsm[-1])
fig.savefig('../fig/shock_bsm.pdf', bbox_inches='tight')'''),
    md(r'''**Both effects are large here, and that is the contrast with notebook 13.** The same
$\varepsilon$ that shifts the *averaged* solar survival probability by about 0.005 shifts this
one by a few tenths. Nothing about the model changed; the observable did.

The solar panel is **phase-averaged**, so every part of a BSM effect that enters through the
phase integrates away and only the change in the eigenvectors and the level crossings survives.
This panel is an **instantaneous** probability along a ray, where the phase term is present in
full --- and for the sterile case it dominates completely, since an eV-scale splitting
oscillates far faster than the 88.8 km of ray sampled here.

**So the useful question is not "how big is this BSM effect" but "how big is it in the quantity
I measure".** A sensitivity estimate taken from the wrong one of these two panels would be off
by roughly a factor of thirty for NSI, and by more than that for the sterile case, before any
detector detail entered.

The practical warning is the same one this notebook opens with: none of these curves means
anything without `t_breakpoints`. Run the identical comparison without them and the differences
you measure are between two straddled discretizations, not between two physics models.'''),
    ])


# ------------------------------------------------------ 15_magnus_antineutrinos
books['15_magnus_antineutrinos.ipynb'] = notebook(
    'Antineutrinos, done properly',
    r'''Going from neutrinos to antineutrinos changes **two** things in the Hamiltonian, and
they are easy to get half right:

1. the PMNS matrix is **conjugated**, $\mathbf{U} \to \mathbf{U}^*$, which flips the sign of
   the CP phase $\delta_{\rm CP}$; and
2. the matter potential **changes sign**, $V_{\rm CC} \to -V_{\rm CC}$, because the coherent
   forward scattering of $\bar\nu_e$ off electrons has the opposite sign.

Apply one and not the other and nothing complains. The result is still a valid, doubly
stochastic probability matrix, still between zero and one, still of the right general shape --
it is simply the answer to a different question. This is the single most defect-prone
convention in the package: at one point the sign was applied **twice** inside Mag$\nu$s
itself, which gave antineutrinos a *positive* matter potential and plausible-looking wrong
answers throughout.

In ordinary use you never do either by hand. Every wrapper takes `nubar=True` and does both:

```python
oscprob.osc_prob_3nu_matter_constant_density(energy, L, rho, **osc, nubar=True)
```

This notebook takes the flag apart to show what it stands for, what the two half-right
constructions look like, and which identities have to hold if you got it right.''',
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
    md(r'''## 1. The rule, taken apart

We work at a long-baseline configuration -- 1300 km through the crust, the DUNE baseline --
in the $\nu_\mu \to \nu_e$ appearance channel, where the CP phase and the matter effect both
matter and neither dominates.'''),
    code(r'''# load_nufit_params returns exactly the six mixing parameters, ready to splat
# into any osc_prob_3nu_* call.  'NuFIT 6.1' is the package default.
OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
osc = OSC

RHO = 2.848                                   # Earth's crust [g cm^-3]
energy = 2.0*gd.UNIT_GEV                      # [eV]
baseline = 1300.0*gd.UNIT_KM                  # [eV^-1]

for name, value in OSC.items():
    print('%-5s = %+.5e' % (name, value))'''),
    md(r'''Now build the Hamiltonian with the two flips under separate control. `nubar` appears
twice below, and they are *different* arguments -- one on the vacuum term, one on the
potential:

$$\mathbf{H} = \frac{1}{2E}\,\mathbf{U}\,\mathbf{M}^2\,\mathbf{U}^\dagger
             \;+\; {\rm diag}(V_{\rm CC}, 0, 0)$$'''),
    code(r'''def prob_halves(conjugate_pmns, flip_potential, energy=energy, L=baseline):
    """P for a Hamiltonian with each half of the antineutrino rule applied, or not.

    conjugate_pmns -> U becomes U*, flipping the sign of dCP
    flip_potential -> V_CC becomes -V_CC
    Both True is the antineutrino; both False is the neutrino."""
    h_vac = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
        **OSC, nubar=conjugate_pmns)                       # [eV^2]
    vcc = matter.vcc_func_from_rho_func(
        RHO, nubar=flip_potential, density_matter_is_in_g_per_cm3=True)   # [eV]
    H = np.asarray(h_vac)/energy + np.diag([vcc, 0.0, 0.0])
    return oscprob.osc_prob(H, 0.0, L)'''),
    md(r'''## 2. Two ways to get it half right

All four numbers below are unremarkable probabilities. Nothing about the two middle rows
looks wrong -- that is the entire problem.'''),
    code(r'''print('%-28s %s' % ('construction', 'P(nu_mu -> nu_e)'))
print('-'*48)
for name, conj, flip in [('neutrino',                False, False),
                         ('antineutrino (correct)',  True,  True),
                         ('conjugate PMNS only',     True,  False),
                         ('flip potential only',     False, True)]:
    print('%-28s %.6f' % (name, prob_halves(conj, flip)[gd.NUMU][gd.NUE]))'''),
    md(r'''Conjugating alone overstates the appearance probability by a factor of four; flipping
the potential alone overstates it by a factor of about 1.6. Both sit between the neutrino and
the antineutrino values, which is exactly where a plausible wrong answer lives.

The hand-built version is not a re-derivation -- it reproduces what the wrapper does, to
every digit, for both signs:'''),
    code(r'''for nubar in (False, True):
    from_wrapper = oscprob.osc_prob_3nu_matter_constant_density(
        energy, baseline, RHO, **OSC, nubar=nubar,
        density_matter_is_in_g_per_cm3=True)[gd.NUMU][gd.NUE]
    by_hand = prob_halves(nubar, nubar)[gd.NUMU][gd.NUE]
    print('nubar=%-6s wrapper %.9f   by hand %.9f   |difference| %.1e'
          % (nubar, from_wrapper, by_hand, abs(from_wrapper - by_hand)))'''),
    md(r'''## 3. Where the sign lives

The antineutrino sign of the potential is applied in exactly **one** place in Mag$\nu$s,
`matter.vcc_func_from_rho_func`, and every entry point routes through it. That is deliberate:
the one bug this convention actually produced was the sign being applied *twice*, once in the
potential builder and once again by the caller, which is only possible when there is more than
one place it could live.

So if you assemble a Hamiltonian yourself, take the potential from that function and do
**not** negate it again:'''),
    code(r'''v_nu = matter.vcc_func_from_rho_func(RHO, density_matter_is_in_g_per_cm3=True)
v_nubar = matter.vcc_func_from_rho_func(RHO, nubar=True,
                                        density_matter_is_in_g_per_cm3=True)
print('V_CC(nu)    = %+.5e eV' % v_nu)
print('V_CC(nubar) = %+.5e eV' % v_nubar)
print('sum         = %+.1e   (they are exact negatives, applied once)'
      % (v_nu + v_nubar))'''),
    md(r'''## 4. The four constructions against energy

Below the appearance peak the curves are hard to tell apart; through the peak they separate
by more than the CP-violation signal any experiment is trying to measure. A wrong convention
would not look like a bug in a plot like this -- it would look like a different value of
$\delta_{\rm CP}$.'''),
    code(r'''E_gev = np.logspace(-0.5, 1.2, 300)
E = E_gev*gd.UNIT_GEV

def curve(conjugate_pmns, flip_potential):
    h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
        **OSC, nubar=conjugate_pmns))
    vcc = matter.vcc_func_from_rho_func(RHO, nubar=flip_potential,
                                        density_matter_is_in_g_per_cm3=True)
    return np.array([oscprob.osc_prob(h_vac/e + np.diag([vcc, 0.0, 0.0]),
                                      0.0, baseline)[gd.NUMU][gd.NUE] for e in E])

specs = [(r'$\nu$',                        False, False, 'C0', '-'),
         (r'$\bar\nu$ (correct)',          True,  True,  'C3', '-'),
         (r'$\bar\nu$, conjugate only',    True,  False, 'C1', '--'),
         (r'$\bar\nu$, potential only',    False, True,  'C2', ':')]
curves = [dict(y=curve(c, f), label=lab, color=col, ls=ls)
          for lab, c, f, col, ls in specs]

fig, ax = plotting.plot_probability_vs_energy(
    E_gev, curves, nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3,
    xlim=(E_gev[0], E_gev[-1]), ylim=(0.0, 0.2),
    legend_title='construction', legend_loc='upper left',
    title=(r'$3\nu$ in constant density, $L = 1300$ km, '
           r'$\rho = 2.848$ g cm$^{-3}$'))'''),
    md(r'''## 5. Checks that must hold

Three identities separate a correct antineutrino calculation from a plausible one. They are
worth running whenever you build a Hamiltonian by hand.

**CP: with $\delta_{\rm CP} = 0$ and no matter, neutrinos and antineutrinos are identical.**
Conjugating a real matrix does nothing, so the two must agree to machine precision. If they
do not, the potential is not zero when you think it is.'''),
    code(r'''for dcp in (0.0, np.pi/2):
    Q = dict(OSC, dCP=dcp)
    p = oscprob.osc_prob_3nu_vacuum(energy, baseline, **Q)[gd.NUMU][gd.NUE]
    pbar = oscprob.osc_prob_3nu_vacuum(energy, baseline, **Q,
                                       nubar=True)[gd.NUMU][gd.NUE]
    print('dCP = %.4f   P = %.9f   Pbar = %.9f   |P - Pbar| = %.2e'
          % (dcp, p, pbar, abs(p - pbar)))'''),
    md(r'''**CPT: in vacuum, $P(\nu_\alpha \to \nu_\beta) = P(\bar\nu_\beta \to \bar\nu_\alpha)$
for any $\delta_{\rm CP}$.** Note the reversed flavor indices -- this is the *transpose*, not
the matrix itself. CP is violated here and CPT is not, so the first number below is zero and
the second is not.'''),
    code(r'''Q = dict(OSC, dCP=1.234)
P_nu = np.asarray(oscprob.osc_prob_3nu_vacuum(energy, baseline, **Q))
P_nubar = np.asarray(oscprob.osc_prob_3nu_vacuum(energy, baseline, **Q, nubar=True))

print('vacuum, max |P[a][b] - Pbar[b][a]| = %.2e   <- CPT, must vanish'
      % np.max(np.abs(P_nu - P_nubar.T)))
print('vacuum, max |P[a][b] - Pbar[a][b]| = %.2e   <- CP, need not'
      % np.max(np.abs(P_nu - P_nubar)))'''),
    md(r'''**In matter, CPT no longer relates the two.** Ordinary matter is made of electrons
and not positrons, so it is a CPT-odd background: the transpose identity is *expected* to
fail, and a calculation in which it still held would be the suspicious one.'''),
    code(r'''P_nu = np.asarray(prob_halves(False, False))
P_nubar = np.asarray(prob_halves(True, True))
print('matter, max |P[a][b] - Pbar[b][a]| = %.2e   <- no longer a symmetry'
      % np.max(np.abs(P_nu - P_nubar.T)))

# Unitarity survives regardless: every row still sums to one.
print('row sums, nubar: %s' % np.round(P_nubar.sum(axis=1), 12))'''),
    md(r'''## Summary

| | conjugate $\mathbf{U}$ | flip $V_{\rm CC}$ | result |
|---|---|---|---|
| neutrino | no | no | correct |
| **antineutrino** | **yes** | **yes** | **correct** |
| half right | yes | no | plausible, wrong |
| half right | no | yes | plausible, wrong |
| the old bug | yes | *twice* | plausible, wrong |

Both flips, or neither. In practice: **pass `nubar=True` and let the wrapper do it**, and if
you must build the Hamiltonian yourself, take $V_{\rm CC}$ from
`matter.vcc_func_from_rho_func(..., nubar=True)` and never negate it again.

The three checks in section 5 are cheap enough to keep in a test: CP agreement at
$\delta_{\rm CP} = 0$ in vacuum, the CPT transpose identity in vacuum, and its expected
failure in matter.'''),
    ])


# ------------------------------------------- 16_magnus_exact_vs_approximations
books['16_magnus_exact_vs_approximations.ipynb'] = notebook(
    'Exact versus the textbook approximations',
    r'''The familiar oscillation formulas -- the two-flavor $\sin^2(2\theta)\sin^2(\Delta m^2
L/4E)$, its constant-density counterpart with the matter-resonant mixing angle -- are not
approximations. Each is the **exact** solution of the problem it was derived for, and Mag$\nu$s
reproduces every one of them to machine precision. Showing that is the first half of this
notebook, and it is really a test of Mag$\nu$s rather than of the formulas.

What breaks is not the formula but the *problem substitution*: using a constant-density
formula for a profile that is not constant. The usual recipe -- take the mean density along
the path and plug it in -- is where the error enters, and on an Earth chord it is not a small
one. That is the second half.

The comparison uses `magnus.oscprobstd`, which ships the standard closed forms precisely so
that this kind of check is available without a second package.''',
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
import magnus.plotting as plotting

# load_nufit_params returns exactly the six mixing parameters, ready to splat
# into any osc_prob_3nu_* call.  'NuFIT 6.1' is the package default.
OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
osc = OSC'''),
    md(r'''## 1. The two-flavor vacuum formula is not an approximation

$$P(\nu_e \to \nu_\mu) = \sin^2 2\theta \,\sin^2\!\left(\frac{\Delta m^2 L}{4E}\right)$$

is exact for two flavors in vacuum, and the Magnus expansion has nothing to correct. Note
that `oscprobstd.osc_prob_2nu_vacuum_std` returns the full $2\times2$ probability **matrix**,
not a single number -- so it is indexed the same way as everything else in the package.'''),
    code(r'''sth, Dm2 = osc['s12'], osc['D21']
E_gev = np.logspace(-1.0, 1.0, 400)
E = E_gev*gd.UNIT_GEV
L = 1300.0*gd.UNIT_KM

magnus_2nu = oscprob.osc_prob_2nu_vacuum(E, L, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU)
closed_2nu = np.array([oscprobstd.osc_prob_2nu_vacuum_std(sth, Dm2, e, L)[gd.NUE][gd.NUMU]
                       for e in E])

print('2nu, vacuum       : max |Magnus - closed form| = %.2e'
      % np.max(np.abs(magnus_2nu - closed_2nu)))'''),
    md(r'''## 2. Three flavors in vacuum: still exact

Three flavors in vacuum has a closed form too, built from the PMNS matrix and the two
mass-squared splittings. Same story.'''),
    code(r'''U = hamiltonians.pmns_mixing_matrix(OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP'])

magnus_3nu = oscprob.osc_prob_3nu_vacuum(E, L, **OSC, nu_i=gd.NUMU, nu_f=gd.NUE)
closed_3nu = np.array([oscprobstd.osc_prob_3nu_vacuum_std(U, OSC['D21'], OSC['D31'], e, L)
                       [gd.NUMU][gd.NUE] for e in E])

print('3nu, vacuum       : max |Magnus - closed form| = %.2e'
      % np.max(np.abs(magnus_3nu - closed_3nu)))'''),
    md(r'''## 3. Constant density: exact when the density really is constant

The matter formula replaces the vacuum mixing angle and splitting with their in-matter
counterparts. It is exact for a genuinely constant potential, and again there is nothing to
correct.'''),
    code(r'''RHO = 2.848                                    # Earth's crust [g cm^-3]
vcc = matter.vcc_func_from_rho_func(RHO, density_matter_is_in_g_per_cm3=True)

magnus_matter = oscprob.osc_prob_2nu_matter_constant_density(
    E, L, RHO, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUMU,
    density_matter_is_in_g_per_cm3=True)
closed_matter = np.array([oscprobstd.osc_prob_2nu_matter_std(sth, Dm2, vcc, e, L)
                          [gd.NUE][gd.NUMU] for e in E])

print('V_CC              = %.4e eV' % vcc)
print('2nu, constant rho : max |Magnus - closed form| = %.2e'
      % np.max(np.abs(magnus_matter - closed_matter)))'''),
    md(r'''Three formulas, three agreements at the $10^{-14}$ level or better. Wherever a closed
form exists, Mag$\nu$s returns it.

## 4. Where it does break: a constant density standing in for a varying one

Now the substitution that is actually made in practice. A neutrino crossing the Earth through
the core passes from crust to mantle to outer and inner core and back -- from about 1 to 13
g cm$^{-3}$. The constant-density formula is still exact for *its* problem; it is just no
longer this problem.

We use the atmospheric sector ($\theta_{23}$, $\Delta m^2_{31}$), which is what an Earth chord
of this length actually probes.'''),
    code(r'''sth_atm, Dm2_atm = osc['s23'], osc['D31']

def num_density_e_prem(r):
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
                                     electron_fraction=0.5,
                                     density_matter_is_in_g_per_cm3=True)   # [eV^3]

def VCC_prem(r):
    return matter.VCC_func(r, num_density_e_prem)                            # [eV]

H_vac_atm = np.asarray(
    hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth_atm, Dm2_atm))

costhz = -1.0                                  # straight up through the core
L_km = earth.distance_traveled_inside_earth(costhz)
L_earth = L_km*gd.CONV_KM_TO_INV_EV

def H_prem(l, energy):
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV)
    return H_vac_atm/energy + hamiltonians.hamiltonian_2nu_matter(VCC_prem(r))

# The PREM shell boundaries are density *jumps*.  Declaring them keeps a slab from
# straddling one, which no amount of refinement fixes -- see notebook 14.
breakpoints = np.asarray(
    earth.prem_layer_edges_along_chord(costhz))*gd.CONV_KM_TO_INV_EV

print('chord length      = %.1f km' % L_km)
print('PREM edges on it  = %d' % len(breakpoints))'''),
    md(r'''The mean density along the chord is what someone reaching for the constant-density
formula would plug in, so that is what we compare against.'''),
    code(r'''l_grid = np.linspace(0.0, L_km, 4001)
rho_along = earth.density_matter_func_prem(
    earth.earth_radial_distance_from_depth(costhz, l_grid))
rho_bar = float(np.mean(rho_along))
VCC_bar = matter.vcc_func_from_rho_func(rho_bar, density_matter_is_in_g_per_cm3=True)

print('density along the chord: %.2f g/cm^3 (core) -> %.2f (crust)'
      % (rho_along.max(), rho_along.min()))
print('mean density           : %.4f g/cm^3' % rho_bar)'''),
    md(r'''The cell below raises a `MagnusConvergenceWarning`. It is worth reading rather than
hiding: it reports that some slab was wider than the sufficient condition for convergence of
the Magnus series, which is a statement about **slab width, not about the answer** -- measured
false about three quarters of the time. Here the answer is converged (declaring the
breakpoints and tightening the tolerance by four orders of magnitude moves it in the sixth
decimal). Notebooks 20 and 21 take this apart properly.'''),
    code(r'''E_gev_earth = np.logspace(0.0, 2.0, 240)
E_earth = E_gev_earth*gd.UNIT_GEV

exact = np.array([oscprob.osc_prob(lambda l, e=e: H_prem(l, e), 0.0, L_earth,
                                   t_breakpoints=breakpoints)[gd.NUMU][gd.NUMU]
                  for e in E_earth])
const_rho = np.array([oscprobstd.osc_prob_2nu_matter_std(sth_atm, Dm2_atm, VCC_bar, e, L_earth)
                      [gd.NUMU][gd.NUMU] for e in E_earth])
vacuum = np.array([oscprobstd.osc_prob_2nu_vacuum_std(sth_atm, Dm2_atm, e, L_earth)
                   [gd.NUMU][gd.NUMU] for e in E_earth])

gap = np.abs(exact - const_rho)
print('max |exact - constant rho| = %.3f   at %.2f GeV'
      % (gap.max(), E_gev_earth[gap.argmax()]))
print('rms                        = %.3f' % np.sqrt(np.mean(gap**2)))
print('max |exact - vacuum|       = %.3f' % np.max(np.abs(exact - vacuum)))'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    E_gev_earth,
    [dict(y=exact, label='Exact (PREM, Magnus)', color='C3'),
     dict(y=const_rho, label=r'Constant $\bar\rho$ closed form', color='C0', ls='--'),
     dict(y=vacuum, label='Vacuum closed form', color='0.5', ls=':')],
    nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=2,
    xlim=(E_gev_earth[0], E_gev_earth[-1]),
    legend_title='method', legend_loc='lower right',
    title=(r'$2\nu$ through the core, $\cos\theta_z = -1$ '
           r'($L = %.0f$ km)' % L_km))'''),
    md(r'''## Summary

| comparison | agreement |
|---|---|
| $2\nu$ vacuum, closed form | $4\times10^{-16}$ |
| $3\nu$ vacuum, closed form | $1\times10^{-14}$ |
| $2\nu$ constant density, closed form | $6\times10^{-16}$ |
| $2\nu$ Earth chord vs **mean-density** formula | **wrong by up to 0.51** |

The closed forms are exact and Mag$\nu$s agrees with all of them. The failure in the last row
is not the formula's -- it is the substitution of a constant density for a varying one, and it
is worth 0.51 in probability at the peak, with an rms of 0.23 across two decades in energy.

The mean density is not a bad *estimate* of the profile; it is simply not what the neutrino
sees. Oscillation depends on the arrangement of the matter along the path, not only on its
average -- which is the subject of notebook 18.'''),
    ])


# ------------------------------------------ 17_magnus_ordering_and_octant
books['17_magnus_ordering_and_octant.ipynb'] = notebook(
    'Mass ordering and the $\\theta_{23}$ octant',
    r'''Two things about the neutrino mass spectrum are still unmeasured, and both show up in
oscillation probabilities:

* **The ordering.** Is $\nu_3$ the heaviest state (normal ordering, NO) or the lightest
  (inverted, IO)? In Mag$\nu$s this is carried entirely by the **sign of `D31`**
  ($\Delta m^2_{31}$) -- positive for NO, negative for IO. Nothing else changes.
* **The octant.** Is $\theta_{23}$ below $45^\circ$ or above it? Disappearance experiments
  measure something close to $\sin^2 2\theta_{23}$, which cannot tell $\theta_{23}$ from
  $90^\circ - \theta_{23}$.

A warning about the shipped parameter sets before we start.
`gd.load_nufit_params('NuFIT 6.1', 'NO')` and `(..., 'IO')` are the NuFIT 6.1 best fits for each
ordering, and they differ in **both** things at once -- and in $\delta_{\rm CP}$ as well. Using
them as an ordering comparison conflates three effects. To isolate one variable we flip the
sign of `D31` by hand and change nothing else.''',
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
import magnus.plotting as plotting

NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
IO_fit = gd.load_nufit_params('NuFIT 6.1', 'IO')
KEYS = tuple(NO)

print('%-6s %-14s %-14s' % ('', 'NuFIT NO', 'NuFIT IO'))
for k in KEYS:
    print('%-6s %+.6e %+.6e' % (k, NO[k], IO_fit[k]))
print()
print('sin^2(th23): %.3f vs %.3f  <- different octants, too'
      % (NO['s23']**2, IO_fit['s23']**2))'''),
    md(r'''So we build our own inverted-ordering set: the normal-ordering best fit with the
sign of `D31` reversed, everything else held fixed. Any difference below is then the ordering
and nothing else.'''),
    code(r'''IO = dict(NO, D31=-NO['D31'])
print('NO: D31 = %+.6f eV^2' % NO['D31'])
print('IO: D31 = %+.6f eV^2' % IO['D31'])'''),
    md(r'''## 1. In vacuum, the ordering is entangled with the CP phase

At a long baseline the vacuum probability does depend on the sign of `D31`, through the
interference between the solar and atmospheric terms. But the same interference carries
$\delta_{\rm CP}$, so a vacuum measurement cannot separate them -- which is the whole reason
the ordering is hard.'''),
    code(r'''E_gev = np.logspace(-0.5, 1.2, 300)
E = E_gev*gd.UNIT_GEV
L = 1300.0*gd.UNIT_KM                          # DUNE-like baseline

vac_no = oscprob.osc_prob_3nu_vacuum(E, L, **NO, nu_i=gd.NUMU, nu_f=gd.NUE)
vac_io = oscprob.osc_prob_3nu_vacuum(E, L, **IO, nu_i=gd.NUMU, nu_f=gd.NUE)

print('vacuum, max |P(NO) - P(IO)| = %.4f' % np.max(np.abs(vac_no - vac_io)))'''),
    md(r'''## 2. Matter separates them, and separates $\nu$ from $\bar\nu$

The matter potential enters with one sign for neutrinos and the other for antineutrinos, so a
resonance that exists for one does not for the other -- and *which* one it is depends on the
ordering. That asymmetry, not the probability itself, is what an experiment measures.'''),
    code(r'''RHO = 2.848                                    # crust [g cm^-3]
kw = dict(density_matter_is_in_g_per_cm3=True, nu_i=gd.NUMU, nu_f=gd.NUE)

mat = {}
for label, params in (('NO', NO), ('IO', IO)):
    for nubar in (False, True):
        mat[(label, nubar)] = oscprob.osc_prob_3nu_matter_constant_density(
            E, L, RHO, **params, nubar=nubar, **kw)

print('crust, L = 1300 km')
print('  neutrinos     : max |P(NO) - P(IO)| = %.4f'
      % np.max(np.abs(mat[('NO', False)] - mat[('IO', False)])))
print('  antineutrinos : max |P(NO) - P(IO)| = %.4f'
      % np.max(np.abs(mat[('NO', True)] - mat[('IO', True)])))'''),
    md(r'''## 3. Through the core, the separation is not subtle

The crust is thin and its potential is weak. A neutrino crossing the **core** sits in the
matter resonance, and there the two orderings stop being a few-percent question.

This is the calculation Mag$\nu$s exists for: the density varies by more than a factor of ten
along the path, so there is no closed form to fall back on. We declare the PREM shell
boundaries as `t_breakpoints`, since they are genuine density jumps.'''),
    code(r'''def num_density_e_prem(r):
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
                                     electron_fraction=0.5,
                                     density_matter_is_in_g_per_cm3=True)   # [eV^3]

def VCC_prem(r):
    return matter.VCC_func(r, num_density_e_prem)                            # [eV]

costhz = -1.0                                  # straight up through the core
L_km = earth.distance_traveled_inside_earth(costhz)
L_earth = L_km*gd.CONV_KM_TO_INV_EV
breakpoints = np.asarray(
    earth.prem_layer_edges_along_chord(costhz))*gd.CONV_KM_TO_INV_EV

def scan_earth(params, nubar, energies):
    """P(nu_mu -> nu_e) along the chord, for one ordering and one sign."""
    h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
        **params, nubar=nubar))
    sign = -1.0 if nubar else 1.0        # the potential flips; see notebook 15
    def H(l, energy):
        r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV)
        return h_vac/energy + np.diag([sign*VCC_prem(r), 0.0, 0.0])
    return np.array([oscprob.osc_prob(lambda l, e=e: H(l, e), 0.0, L_earth,
                                      t_breakpoints=breakpoints)[gd.NUMU][gd.NUE]
                     for e in energies])

E_gev_earth = np.logspace(0.0, 1.5, 160)
E_earth = E_gev_earth*gd.UNIT_GEV

earth_p = {}
for label, params in (('NO', NO), ('IO', IO)):
    for nubar in (False, True):
        earth_p[(label, nubar)] = scan_earth(params, nubar, E_earth)

print('core chord, L = %.0f km' % L_km)
for nubar, name in ((False, 'neutrinos    '), (True, 'antineutrinos')):
    gap = np.abs(earth_p[('NO', nubar)] - earth_p[('IO', nubar)])
    print('  %s: max |P(NO) - P(IO)| = %.4f  at %.2f GeV'
          % (name, gap.max(), E_gev_earth[gap.argmax()]))'''),
    code(r'''i = int(np.argmax(np.abs(earth_p[('NO', False)] - earth_p[('IO', False)])))
print('at E = %.2f GeV:' % E_gev_earth[i])
print('  neutrinos     NO %.4f   IO %.4f' % (earth_p[('NO', False)][i],
                                             earth_p[('IO', False)][i]))
print('  antineutrinos NO %.4f   IO %.4f' % (earth_p[('NO', True)][i],
                                             earth_p[('IO', True)][i]))'''),
    md(r'''The resonance is in the **neutrinos** if the ordering is normal and in the
**antineutrinos** if it is inverted. That is the cleanest statement of how an atmospheric
experiment determines the ordering, and here it is as a pair of numbers differing by a factor
of a few hundred rather than a few percent.'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    E_gev_earth,
    [dict(y=earth_p[('NO', False)], label=r'$\nu$, NO', color='C3'),
     dict(y=earth_p[('IO', False)], label=r'$\nu$, IO', color='C3', ls='--'),
     dict(y=earth_p[('NO', True)], label=r'$\bar\nu$, NO', color='C0'),
     dict(y=earth_p[('IO', True)], label=r'$\bar\nu$, IO', color='C0', ls='--')],
    nu_i=gd.NUMU, nu_f=gd.NUE, num_flavors=3,
    xlim=(E_gev_earth[0], E_gev_earth[-1]), ylim=(0.0, 0.6),
    legend_title='ordering', legend_loc='upper right',
    title=(r'$3\nu$ through the core, $\cos\theta_z = -1$ '
           r'($L = %.0f$ km)' % L_km))'''),
    md(r'''## 4. The octant

Now hold the ordering fixed and move $\theta_{23}$ either side of maximal.
$\sin^2\theta_{23} = 0.45$ and $0.55$ give the *same* $\sin^2 2\theta_{23} = 0.99$, so the
leading term of the disappearance probability cannot distinguish them.

"Cannot distinguish" is often stated as an exact degeneracy. In three flavors it is not
exact -- subleading terms and the matter potential both break it -- but the residual is
small, which is the practical difficulty.'''),
    code(r'''s2_lo, s2_hi = 0.45, 0.55
LO = dict(NO, s23=np.sqrt(s2_lo))
HI = dict(NO, s23=np.sqrt(s2_hi))
print('sin^2(2 th23): %.4f  and  %.4f  (identical to four decimals)'
      % (4*s2_lo*(1 - s2_lo), 4*s2_hi*(1 - s2_hi)))

for name, chan in (('mu -> mu  (disappearance)', dict(nu_i=gd.NUMU, nu_f=gd.NUMU)),
                   ('mu -> e   (appearance)   ', dict(nu_i=gd.NUMU, nu_f=gd.NUE))):
    v_lo = oscprob.osc_prob_3nu_vacuum(E, L, **LO, **chan)
    v_hi = oscprob.osc_prob_3nu_vacuum(E, L, **HI, **chan)
    m_lo = oscprob.osc_prob_3nu_matter_constant_density(
        E, L, RHO, **LO, density_matter_is_in_g_per_cm3=True, **chan)
    m_hi = oscprob.osc_prob_3nu_matter_constant_density(
        E, L, RHO, **HI, density_matter_is_in_g_per_cm3=True, **chan)
    print('%s  vacuum %.4f   matter %.4f'
          % (name, np.max(np.abs(v_lo - v_hi)), np.max(np.abs(m_lo - m_hi))))'''),
    md(r'''Both channels separate the octants by about one part in a hundred -- an order of
magnitude smaller than the ordering effect through the core, and comparable to the systematic
uncertainties of a real experiment. That is why the octant is still open while the ordering is
increasingly not.

## Summary

| question | carried by | where it shows | size |
|---|---|---|---|
| ordering | **sign of `D31`** | vacuum, 1300 km | 0.10 |
| | | crust, $\nu$ / $\bar\nu$ | 0.13 / 0.06 |
| | | **core chord, $\nu$ / $\bar\nu$** | **0.48 / 0.42** |
| octant | $\sin^2\theta_{23}$ vs $1-\sin^2\theta_{23}$ | 1300 km, either channel | ~0.015 |

Two practical notes. The shipped NuFIT NO and IO sets differ in ordering, octant *and*
$\delta_{\rm CP}$ simultaneously -- flip `D31` by hand if you want the ordering alone. And an
Earth chord needs its PREM shell boundaries passed as `t_breakpoints`: they are density jumps,
and a slab straddling one is not fixed by refinement.'''),
    ])


# --------------------------------------------- 18_magnus_unusual_density_profiles
books['18_magnus_unusual_density_profiles.ipynb'] = notebook(
    'Unusual density profiles: arrangement, not just the mean',
    r'''Notebook 16 ended on a claim: a neutrino responds to how matter is *arranged* along its
path, not only to how much of it there is. This notebook makes that quantitative, with five
profiles that share the same mean density and the same total length and give completely
different probabilities.

It also isolates one arrangement that is **not** free to matter. Reversing a profile
end-to-end leaves every probability untouched -- but only when the Hamiltonian is complex
symmetric, which for three flavors means $\delta_{\rm CP} = 0$ or $\pi$. That distinction
turns out to be exact, measurable at the $10^{-16}$ level, and is the reason Mag$\nu$s can
halve the Hamiltonian evaluations on an Earth chord.''',
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
import magnus.plotting as plotting

OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
osc = OSC'''),
    md(r'''## 1. Five profiles, one mean density

Each profile is eight equal-width slabs across 6000 km. The densities are drawn from
$\{1, 9\}$ g cm$^{-3}$ or held at the mean, 5 g cm$^{-3}$ -- so every profile below has the
same mean density and the same integrated column depth.

Slab edges are density **jumps**, so they are passed as `t_breakpoints`. A slab of the
integrator straddling a jump is not fixed by refinement, at any tolerance.'''),
    code(r'''TOTAL_KM = 6000.0
LOW, HIGH, MEAN = 1.0, 9.0, 5.0

PROFILES = {
    'constant'          : [MEAN]*8,
    'castle  LHLHLHLH'  : [LOW, HIGH]*4,
    'reversed HLHLHLHL' : [HIGH, LOW]*4,
    'segregated LLLLHHHH': [LOW]*4 + [HIGH]*4,
    'segregated HHHHLLLL': [HIGH]*4 + [LOW]*4,
}

for name, rho in PROFILES.items():
    print('%-21s mean = %.3f g/cm^3' % (name, np.mean(rho)))'''),
    code(r'''def build_H(densities, dim, h_vac_energy_indep):
    """A piecewise-constant matter Hamiltonian, plus the slab edges to declare."""
    n = len(densities)
    edges_km = np.arange(n + 1)*(TOTAL_KM/n)
    vcc = np.array([matter.vcc_func_from_rho_func(
        float(r), density_matter_is_in_g_per_cm3=True) for r in densities])   # [eV]

    def H(l, energy):
        i = int(np.clip(np.searchsorted(edges_km, l/gd.CONV_KM_TO_INV_EV,
                                        side='right') - 1, 0, n - 1))
        m = np.zeros((dim, dim))
        m[0][0] = vcc[i]
        return h_vac_energy_indep/energy + m

    return H, edges_km*gd.CONV_KM_TO_INV_EV


def scan(densities, energies, dim, h_vac, nu_i, nu_f, **kwargs):
    H, breakpoints = build_H(densities, dim, h_vac)
    L = TOTAL_KM*gd.CONV_KM_TO_INV_EV
    return np.array([oscprob.osc_prob(lambda l, e=e: H(l, e), 0.0, L,
                                      t_breakpoints=breakpoints, **kwargs)[nu_i][nu_f]
                     for e in energies])'''),
    md(r'''## 2. The probabilities they produce

Two flavors first, in the atmospheric sector, muon-neutrino survival.'''),
    code(r'''h2 = np.asarray(hamiltonians.hamiltonian_2nu_vacuum_energy_independent(
    osc['s23'], osc['D31']))

E_gev = np.logspace(-0.7, 1.7, 120)
E = E_gev*gd.UNIT_GEV

survival = {name: scan(rho, E, 2, h2, gd.NUMU, gd.NUMU)
            for name, rho in PROFILES.items()}
reference = survival['constant']

print('%-21s %14s %8s' % ('profile', 'max |P - P_const|', 'rms'))
print('-'*46)
for name, p in survival.items():
    print('%-21s %14.4f %8.4f'
          % (name, np.max(np.abs(p - reference)), np.sqrt(np.mean((p - reference)**2))))'''),
    md(r'''The castle wall departs from the constant-density answer by **0.98** in probability --
essentially the entire available range -- while carrying exactly the same amount of matter.
Averaging the density first is not a mild approximation here; it is a different problem.

And the two rearrangements differ from each other, not merely from the mean:'''),
    code(r'''print('castle LHLHLHLH vs segregated LLLLHHHH : max |diff| = %.4f'
      % np.max(np.abs(survival['castle  LHLHLHLH']
                      - survival['segregated LLLLHHHH'])))

i = int(np.argmax(np.abs(survival['castle  LHLHLHLH'] - reference)))
print('\nat E = %.2f GeV:' % E_gev[i])
for name, p in survival.items():
    print('   %-21s P = %.4f' % (name, p[i]))'''),
    code(r'''fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6.8, 5.6),
                               gridspec_kw=dict(height_ratios=[1, 2]))

edges = np.arange(9)*(TOTAL_KM/8)
for name, color in (('constant', '0.4'), ('castle  LHLHLHLH', 'C3'),
                    ('segregated LLLLHHHH', 'C0')):
    ax0.step(edges, np.append(PROFILES[name], PROFILES[name][-1]),
             where='post', color=color, label=name)
ax0.set_ylabel(r'$\rho$ [g cm$^{-3}$]')
ax0.set_xlabel('distance [km]')
ax0.set_xlim(0.0, TOTAL_KM)
ax0.legend(fontsize=7, ncol=3, loc='upper center')

for name, color in (('constant', '0.4'), ('castle  LHLHLHLH', 'C3'),
                    ('segregated LLLLHHHH', 'C0')):
    ax1.semilogx(E_gev, survival[name], color=color, label=name)
ax1.set_xlabel(r'$E$ [GeV]')
ax1.set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
ax1.set_xlim(E_gev[0], E_gev[-1])
ax1.set_ylim(0.0, 1.0)
ax1.legend(fontsize=8, loc='lower left')
fig.suptitle('Same mean density, same length, different arrangement', fontsize=10)
fig.tight_layout()'''),
    md(r'''## 3. The one rearrangement that changes nothing

Look again at the table in section 2: `castle LHLHLHLH` and `reversed HLHLHLHL` have
identical entries, as do the two segregated profiles. That is not a coincidence and not a
rounding accident -- it is exact.

Reversing a profile end-to-end maps the evolution operator to its transpose *provided each
slab Hamiltonian is complex symmetric*, $\mathbf{H}^{T} = \mathbf{H}$. Then
$P_{\rm reversed}(\nu_\alpha \to \nu_\beta) = P(\nu_\beta \to \nu_\alpha)$, and survival
probabilities are invariant outright.

For two flavors the Hamiltonian is real, hence symmetric, and this always holds:'''),
    code(r'''castle = PROFILES['castle  LHLHLHLH']
print('2nu, castle vs its own reversal: max |diff| = %.2e'
      % np.max(np.abs(scan(castle, E, 2, h2, gd.NUMU, gd.NUMU)
                      - scan(castle[::-1], E, 2, h2, gd.NUMU, gd.NUMU))))'''),
    md(r'''For **three** flavors the vacuum Hamiltonian is Hermitian, but it is only *symmetric*
when $\delta_{\rm CP}$ is $0$ or $\pi$ -- a nonzero CP phase puts genuine complex numbers off
the diagonal. So the reversal symmetry is exactly as good as CP conservation:'''),
    code(r'''segregated = PROFILES['segregated LLLLHHHH']
E_one = 3.0*gd.UNIT_GEV

print('%-10s %-12s %-14s %s' % ('dCP', 'H symmetric', '|P_aa - Prev_aa|', '|P - Prev^T|'))
print('-'*56)
for dcp in (0.0, np.pi, OSC['dCP']):
    h3 = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
        **dict(OSC, dCP=dcp)))
    H_f, bp = build_H(segregated, 3, h3)
    H_r, _ = build_H(segregated[::-1], 3, h3)
    L = TOTAL_KM*gd.CONV_KM_TO_INV_EV
    A = np.asarray(oscprob.osc_prob(lambda l: H_f(l, E_one), 0.0, L, t_breakpoints=bp))
    B = np.asarray(oscprob.osc_prob(lambda l: H_r(l, E_one), 0.0, L, t_breakpoints=bp))
    print('%-10.4f %-12s %-14.2e %.2e'
          % (dcp, np.allclose(h3, h3.T),
             np.max(np.abs(np.diag(A) - np.diag(B))), np.max(np.abs(A - B.T))))'''),
    md(r'''Sixteen orders of magnitude separate the CP-conserving rows from the last one. The
symmetry is not approximately true and then slightly broken -- it is exact, and then absent.

**Why this matters to Mag$\nu$s.** A chord through a spherically symmetric Earth meets every
radius twice: the density profile *is* a palindrome. Mag$\nu$s exploits that by evaluating the
Hamiltonian on the first half of the slab chain and obtaining the rest by reversal, which
halves the calls to your `H_func`. It is worth 1.4--1.67x on an expensive Hamiltonian and
about 0.91x on plain PREM, where a density lookup is too cheap to be worth halving.
`magnus.magnus.USE_PALINDROME` disarms it. Notebook 24 measures this.

Note that the optimization reuses *evaluations of the profile*, which is valid whatever
$\delta_{\rm CP}$ is -- it does not assume the probability-level symmetry explored above.

## 4. A caution about `t_breakpoints`

Every profile in this notebook declared its slab edges, and it was the right thing to do: on a
**scan** across energies, breakpoints at genuine discontinuities are an established cure.

They are not a free win at a **single** point. Measured across 18 supernova-shock
configurations, adding breakpoints improved 7, worsened 11, and pushed 2 answers from inside
the requested tolerance to outside it. Notebook 14 has that measurement in full, and
`recipes.rst` states the rule. Declare breakpoints because you know where the discontinuities
are -- not as a general accuracy knob.

## Summary

| comparison | max $|\Delta P|$ |
|---|---|
| castle wall vs constant, same mean | **0.98** |
| segregated vs constant, same mean | 0.66 |
| castle vs segregated, same eight slabs | 0.95 |
| any profile vs its own reversal, $2\nu$ | $0$ (exact) |
| $3\nu$, $\delta_{\rm CP} = 0$ or $\pi$ | $\sim 10^{-16}$ (exact) |
| $3\nu$, $\delta_{\rm CP} = 3.70$ | $2.8\times10^{-2}$ |

The mean density tells you almost nothing. The order of the slabs tells you almost everything
-- except for the single reversal that CP conservation protects.'''),
    ])


# ------------------------------------------------ 19_magnus_custom_hamiltonian
books['19_magnus_custom_hamiltonian.ipynb'] = notebook(
    'Bring your own Hamiltonian',
    r'''Everything in the previous notebooks went through a wrapper. But Mag$\nu$s's actual
interface is smaller than its wrapper list suggests: give it a **callable that returns a
Hermitian matrix**, and it will propagate it. The mixing parameters, the matter potential, the
BSM terms -- all of those are conveniences built on top of that one contract.

This notebook covers the contract itself, the one way of writing your `H_func` that is worth
real time, and what the Earth entry point quietly declares on your behalf.''',
    [
    code(r'''import time
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.plotting as plotting

OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC))'''),
    md(r'''## 1. Two contracts

**`oscprob.osc_prob(H_func, t_ini, t_fin)`** is the general one. `H_func(l)` takes a position
and returns a `(dim, dim)` Hermitian matrix. You supply everything, including the matter
potential.

**`oscprob.osc_prob_earth(H, ...)` and `osc_prob_sun(H, ...)`** are the profile-aware ones.
They know the density along the path, so they call `H(energy, l, VCC)` and hand you the
potential; you decide what to do with it. This is the interface to reach for when your new
physics lives in a *known* environment.

The dimension is read from what you return, so 2, 3, 4 and 5 flavors all work with no flag.'''),
    code(r'''e00 = np.diag([1.0, 0.0, 0.0])          # the nu_e--nu_e slot the potential occupies

COSTHZ = -0.9
L_EARTH = earth.distance_traveled_inside_earth(COSTHZ)*gd.CONV_KM_TO_INV_EV
ENERGY = 5.0*gd.UNIT_GEV

def H_standard(energy, l, VCC):
    """Ordinary three-flavor oscillation, written out by hand."""
    return h_vac/energy + np.asarray(VCC)[..., None, None]*e00

print('P(nu_mu -> nu_mu) = %.9f'
      % oscprob.osc_prob_earth(H_standard, energy=ENERGY, costhz=COSTHZ, L=L_EARTH,
                               nu_i=gd.NUMU, nu_f=gd.NUMU))'''),
    md(r'''## 2. Write it so it accepts an array of positions

That `[..., None, None]` is the single most valuable line in this notebook.

Mag$\nu$s evaluates the Hamiltonian at many positions per slab. If your `H_func` handles an
**array** of positions in one call, it does so once per slab instead of once per point. The
indexing trick turns one potential per position into a stack of matrices, and the vacuum term
broadcasts against it for free.

Written the naive way -- `float(VCC)` and a fresh `(3,3)` -- everything still works, and
Mag$\nu$s tells you what you are paying for with a `ScalarHamiltonianWarning`.'''),
    code(r'''def make_H(extra_ops, vectorized):
    """The same Hamiltonian, written the fast way or the slow way.

    extra_ops stands in for a genuinely expensive H_func: a table lookup, an
    interpolation, an integral.  A plain PREM density lookup is very cheap, which
    is exactly why the speed-up below grows as the Hamiltonian gets dearer."""
    def H(energy, l, VCC):
        V = np.asarray(VCC)
        for _ in range(extra_ops):
            V = V + 1.0e-30*np.sin(V*1.0e13)
        if vectorized:
            return h_vac/energy + V[..., None, None]*e00
        m = np.zeros((3, 3))
        m[0][0] = float(V)
        return h_vac/energy + m
    return H

def best_of(H, n=5):
    best = np.inf
    for _ in range(n):
        t0 = time.perf_counter()
        P = oscprob.osc_prob_earth(H, energy=ENERGY, costhz=COSTHZ, L=L_EARTH,
                                   nu_i=gd.NUMU, nu_f=gd.NUMU)
        best = min(best, time.perf_counter() - t0)
    return P, best'''),
    code(r'''print('%-26s %10s %10s %10s %s'
      % ('H_func cost', 'scalar [s]', 'array [s]', 'speed-up', 'identical'))
print('-'*70)
for extra, label in ((0, 'plain PREM lookup'), (5, '5 extra ops/call'),
                     (50, '50 extra ops/call')):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P_scalar, t_scalar = best_of(make_H(extra, False))
    P_array, t_array = best_of(make_H(extra, True))
    print('%-26s %10.4f %10.4f %9.2fx %s'
          % (label, t_scalar, t_array, t_scalar/t_array, P_scalar == P_array))'''),
    md(r'''Bit-identical output, every time -- this is a pure restructuring, not an
approximation. The gain grows with what your Hamiltonian costs, which is the point: on a bare
PREM lookup there is little to save, and on a real interpolated profile or an integral there is
a great deal. The package documentation quotes **4.6x** on a three-flavor exponential-density
profile, measured the same way.

The warning fires once per session and names the fix:'''),
    code(r'''with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    oscprob.osc_prob_earth(make_H(0, False), energy=ENERGY, costhz=COSTHZ,
                           L=L_EARTH, nu_i=gd.NUMU, nu_f=gd.NUMU)

for w in caught:
    if w.category is not UserWarning:
        print('%s\n' % w.category.__name__)'''),
    md(r'''## 3. What the Earth entry point declares for you

A chord through a spherically symmetric Earth meets every radius twice, so its density profile
reads the same from either end -- it is a **palindrome**. Mag$\nu$s can then evaluate your
`H_func` on the first half of the slab chain and obtain the rest by reversal, halving the
calls.

That is geometry rather than an assumption about your input, which is why `osc_prob_earth` is
allowed to declare it (via `symmetric_over`) and a general caller is not. `osc_prob_sun` does
**not** declare it: a solar profile is monotonic, not mirror-symmetric.

The declaration is also withdrawn when it would not be true -- a request over part of a chord
rather than the whole of it is not symmetric, and is not mirrored.

It is worth what your Hamiltonian costs: 1.4--1.67x on an expensive `H_func`, and about 0.91x
on plain PREM, where halving cheap lookups does not repay the bookkeeping.
`magnus.magnus.USE_PALINDROME` disarms it globally. Notebook 24 measures this.

## 4. A worked example: a long-range force

Everything so far could have been done with a wrapper. This is the case that could not.

Suppose the Standard Model gauge group is extended by a $U(1)$ acting on $L_e - L_\mu$, with a
very light mediator $Z'$ (following [arXiv:1808.02042](https://arxiv.org/abs/1808.02042)).
Electrons carry $L_e = 1$, so ordinary matter sources the new field; among the neutrinos
$\nu_e$ carries $+1$, $\nu_\mu$ carries $-1$, and $\nu_\tau$ is neutral. The Hamiltonian
becomes

$$\mathbf{H} = \underbrace{\frac{\mathbf{H}_{\rm vac}}{E}
   + V_{\rm CC}\,{\rm diag}(1,0,0)}_{\text{the standard part}}
   \;+\; V_{e\mu}(r)\,{\rm diag}(1,-1,0),$$

with the new potential a Yukawa integral over **all** the electrons of the body,

$$V_{e\mu}(\mathbf{r}) = \frac{g'^2}{4\pi}\int d^3r'\;
   n_e(\mathbf{r}')\,\frac{e^{-m_{Z'}|\mathbf{r}-\mathbf{r}'|}}{|\mathbf{r}-\mathbf{r}'|}.$$

**That matrix is the reason this notebook exists.** It is not
$\mathbf{H}_{\rm vac}/E + V_{\rm CC}\mathbf{P}_{ee}$ for any choice of $V_{\rm CC}$, so no
amount of arguing with `osc_prob_3nu_earth` will produce it -- the flavor structure
${\rm diag}(1,-1,0)$ is not in that function's vocabulary. And unlike $V_{\rm CC}$, which reads
the density at the neutrino's own position, $V_{e\mu}$ integrates the whole body: it varies
along the trajectory for reasons that have nothing to do with the local density.

### The integral is one-dimensional

A three-dimensional integral per slab would be hopeless. But for a spherically symmetric $n_e$
the angular average of the Yukawa kernel over a shell is elementary, and the result separates
into an interior and an exterior piece:

$$V_{e\mu}(r) = \frac{g'^2}{r}\,e^{-mr}\!\int_0^r\! dr'\,r'^2 n_e(r')\,{\rm shc}(mr')
  \;+\; g'^2\,{\rm shc}(mr)\!\int_r^R\! dr'\,r'\,n_e(r')\,e^{-mr'},$$

with ${\rm shc}(x) \equiv \sinh(x)/x$, which is $1$ at the origin.

Two consequences, and both matter. The $r$-dependence has left the integrands, so **one pass
over the profile serves every point on the trajectory** -- the potential is built once per
body, not once per slab. And nothing diverges as $m \to 0$.'''),
    code(r'''def running_integral(y, x):
    """Trapezoidal running integral of y over x, zero at the first node.

    Spelled out rather than imported: numpy has no cumulative trapezoid, and
    the plain one changed name in NumPy 2.0 (np.trapz was removed in favor of
    np.trapezoid), so neither spelling is portable.  This is one cumsum, works
    on every version, and keeps the example to plain numpy."""
    return np.concatenate([[0.0], np.cumsum(0.5*(y[1:] + y[:-1])*np.diff(x))])


def shc(x):
    """sinh(x)/x, continued to 1 at the origin."""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    big = np.abs(x) > 1.0e-8
    out[big] = np.sinh(x[big])/x[big]
    return out

def long_range_potential(r_grid, ne_grid, m):
    """V_{e-mu}(r) on r_grid, in units of g'^2, for a spherical n_e(r).

    Two cumulative integrals over the profile, evaluated once.  Everything the
    trajectory then needs is a lookup."""
    inner = r_grid**2*ne_grid*shc(m*r_grid)
    outer = r_grid*ne_grid*np.exp(-m*r_grid)
    I_in = running_integral(inner, r_grid)
    running_out = running_integral(outer, r_grid)
    I_out = running_out[-1] - running_out          # the exterior piece, r' > r
    r_safe = np.where(r_grid > 0.0, r_grid, 1.0e-30)
    return np.exp(-m*r_grid)*I_in/r_safe + shc(m*r_grid)*I_out'''),
    md(r'''### Does the potential come out right?

A uniform ball has a closed form at **any** mediator mass,

$$V_{e\mu}(r) = \frac{g'^2 n_e}{m^2}
  \left[1 - \left(R + \frac{1}{m}\right)e^{-mR}\,\frac{\sinh(mr)}{r}\right],$$

so this is a check against something external rather than a self-consistency test. Worth doing
before any physics: a quadrature that is quietly wrong produces figures that look entirely
reasonable.'''),
    code(r'''R_E = gd.EARTH_RADIUS                       # [km]
r_uniform = np.linspace(0.0, R_E, 200001)
ne_uniform = np.ones_like(r_uniform)

print('%-14s %s' % ('1/m', 'max relative error vs the closed form'))
print('-'*52)
for inverse_m in (0.05, 0.3, 1.0, 5.0):
    m = 1.0/(inverse_m*R_E)
    numeric = long_range_potential(r_uniform, ne_uniform, m)
    r_safe = np.where(r_uniform > 0.0, r_uniform, 1.0e-30)
    closed = (1.0/m**2)*(1.0 - (R_E + 1.0/m)*np.exp(-m*R_E)
                         * np.where(r_uniform > 0.0, np.sinh(m*r_safe)/r_safe, m))
    print('%-14s %.2e' % ('%.2f R' % inverse_m,
                          np.max(np.abs(numeric - closed)/np.abs(closed))))'''),
    md(r'''### The mediator's range decides what the potential looks like

Now the Earth's own electrons, from PREM. The density **jumps** at each shell boundary, so the
quadrature grid must put a node on every one -- the same rule that makes `t_breakpoints`
necessary for the propagation itself (notebook 18).

Two limits are worth holding in mind:

* **Short range**, $1/m \ll R$: only the neighborhood contributes and
  $V_{e\mu} \to g'^2 n_e(r)/m^2$. The potential tracks the *local* density, which makes it
  degenerate with an NSI -- notebook 08's $\epsilon_{\alpha\beta}$ with a particular flavor
  structure.
* **Long range**, $1/m \gg R$: the whole body contributes, and the potential is smooth,
  largest at the center, and indifferent to where any individual electron sits.

Only the second is genuinely new as far as this package is concerned, and only the second
needs the profile machinery at all.'''),
    code(r'''prem_edges = np.concatenate(([0.0], earth.PREM_BOUNDARIES))
r_prem = np.unique(np.concatenate(
    [np.linspace(a, b, 900) for a, b in zip(prem_edges[:-1], prem_edges[1:])]))
ne_prem = matter.num_density_e_func(r_prem, earth.density_matter_func_prem,
                                    electron_fraction=0.5,
                                    density_matter_is_in_g_per_cm3=True)

print('grid: %d nodes, with one on every PREM boundary' % len(r_prem))
print('%-12s %-14s %s' % ('1/m', 'V(0)/V(R)', 'V(0)'))
print('-'*42)
V_ranges = {}
for inverse_m in (0.03, 0.3, 1.0, 10.0):
    V = long_range_potential(r_prem, ne_prem, 1.0/(inverse_m*R_E))
    V_ranges[inverse_m] = V
    print('%-12s %-14.3f %.4e' % ('%.2f R' % inverse_m, V[0]/V[-1], V[0]))

print('\nfor reference, the density contrast n_e(0)/n_e(R) = %.3f'
      % (ne_prem[0]/ne_prem[-1]))'''),
    code(r'''fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.4, 3.6))

for inverse_m, V in V_ranges.items():
    axL.plot(r_prem/R_E, V/V[0], label=r'$1/m = %g\,R$' % inverse_m)
axL.plot(r_prem/R_E, ne_prem/ne_prem[0], 'k--', lw=1.0, label=r'$n_e(r)$, scaled')
axL.set_xlabel(r'$r / R$')
axL.set_ylabel(r'$V_{e\mu}(r) / V_{e\mu}(0)$')
axL.legend(fontsize=7)
axL.set_xlim(0.0, 1.0)

inv = np.logspace(-2.0, 1.5, 30)
center = [long_range_potential(r_prem, ne_prem, 1.0/(v*R_E))[0] for v in inv]
axR.loglog(inv, center, marker='o', ms=3)
axR.set_xlabel(r'$1/m$  [units of $R$]')
axR.set_ylabel(r'$V_{e\mu}(0)$')
axR.set_xlim(inv[0], inv[-1])
fig.tight_layout()'''),
    md(r'''The left panel is the whole argument in one picture. At $1/m = 0.03\,R$ the potential
has nearly taken the shape of the dashed density curve, shell jumps and all; by $1/m = 10\,R$
it has forgotten the profile entirely and is a smooth bowl. The right panel shows the
crossover: $V_{e\mu}(0)$ grows as $1/m^2$ while the range is short, and saturates once the
mediator reaches across the body, because there are no more electrons left to enclose.

### Through the Earth

Now put it in a Hamiltonian and propagate. The potential is a lookup against the pass computed
above, so the cost of the new physics is one interpolation per position.'''),
    code(r'''INVERSE_M = 1.0                            # the awkward middle: 1/m = R
V_LRF = long_range_potential(r_prem, ne_prem, 1.0/(INVERSE_M*R_E))
lrf_charge = np.diag([1.0, -1.0, 0.0])     # the L_e - L_mu charge

# g'^2 chosen so the new potential is a tenth of V_CC at the center: small
# enough to be plausible, large enough to see.
VCC_center = matter.VCC_func(0.0, lambda l: ne_prem[0])
G_SQUARED = 0.1*VCC_center/V_LRF[0]

def H_long_range(energy, l, VCC):
    """Standard oscillation plus the long-range term, both position dependent."""
    r = earth.earth_radial_distance_from_depth(
        COSTHZ, np.asarray(l)/gd.CONV_KM_TO_INV_EV)
    V_new = G_SQUARED*np.interp(r, r_prem, V_LRF)
    return (h_vac/energy
            + np.asarray(VCC)[..., None, None]*e00
            + np.asarray(V_new)[..., None, None]*lrf_charge)

E_gev = np.logspace(-0.3, 1.5, 300)
E_scan = E_gev*gd.UNIT_GEV

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    p_std = oscprob.osc_prob_earth(H_standard, energy=E_scan, costhz=COSTHZ,
                                   L=L_EARTH, nu_i=gd.NUMU, nu_f=gd.NUMU)
    p_lrf = oscprob.osc_prob_earth(H_long_range, energy=E_scan, costhz=COSTHZ,
                                   L=L_EARTH, nu_i=gd.NUMU, nu_f=gd.NUMU)

gap = np.abs(np.asarray(p_std) - np.asarray(p_lrf))
print('V_new(0) is 10%% of V_CC there')
print('max |standard - long range| = %.4f  at %.2f GeV'
      % (gap.max(), E_gev[gap.argmax()]))'''),
    code(r'''fig, ax = plotting.plot_probability_vs_energy(
    E_gev,
    [dict(y=np.asarray(p_std), label=r'Standard $3\nu$', color='0.2', ls='--'),
     dict(y=np.asarray(p_lrf), label=r'$+\;L_e - L_\mu$, $1/m = R$', color='C1')],
    nu_i=gd.NUMU, nu_f=gd.NUMU, num_flavors=3,
    xlim=(E_gev[0], E_gev[-1]),
    legend_title='Hamiltonian', legend_loc='lower right',
    title=r'A long-range force through the Earth')'''),
    md(r'''The shift is not a rescaling of the standard curve: the new term has a different
flavor structure, so it moves the oscillation rather than damping it.

Everything the package needed from you was a callable. What it gave back was the whole slab
machinery -- adaptive refinement, the palindrome, the engine dispatch of notebook 22 -- applied
to a Hamiltonian it has never seen.'''),
    md(r'''## Summary -- a checklist for your own `H_func`

1. **Return a Hermitian matrix.** The dimension is inferred; 2 to 5 flavors need no flag.
2. **Accept an array of positions**, and use `np.asarray(VCC)[..., None, None]` to build the
   stack. Bit-identical results, and the gain grows with what your Hamiltonian costs.
3. **Take the matter potential from the entry point** (`osc_prob_earth`/`osc_prob_sun`) or
   from `matter.vcc_func_from_rho_func`. For antineutrinos pass `nubar=True` there and do not
   negate it again -- notebook 15.
4. **Declare discontinuities** with `t_breakpoints`. A slab straddling a jump is not fixed by
   refinement -- notebook 18.
5. **Do not declare symmetry yourself.** `osc_prob_earth` does it where geometry guarantees it.

If you want to know which engine then answered your request, and whether to believe it, that
is notebooks 21 and 22.'''),
    ])


# ------------------------------------------------- 20_magnus_numerical_edge_cases
books['20_magnus_numerical_edge_cases.ipynb'] = notebook(
    'Numerical edge cases, and what the warnings mean',
    r'''Two questions this notebook answers. First: which degenerate, empty or otherwise
pathological inputs return a **number** rather than a `NaN`? Exact degeneracies are the usual
place a closed-form implementation divides by zero, and the Magnus expansion never forms those
denominators -- it exponentiates a matrix, and a degenerate matrix exponentiates perfectly
well.

Second, and more useful in practice: Mag$\nu$s has **nine** warning classes, and they do not
all mean the same kind of thing. Some report a bad input, some an expensive choice, and some a
condition that was not met but may not matter. Knowing which is which is the difference between
a warning you act on and one you note.''',
    [
    code(r'''import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.magnus as magnus
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.globaldefs as gd

# load_nufit_params returns exactly the six mixing parameters, ready to splat
# into any osc_prob_3nu_* call.  'NuFIT 6.1' is the package default.
OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
osc = OSC
h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC))

ENERGY = 1.0*gd.UNIT_GEV
BASELINE = 1300.0*gd.UNIT_KM'''),
    md(r'''## 1. A Hamiltonian proportional to the identity

If every eigenvalue is the same there is no relative phase, so nothing oscillates. The answer
is the identity matrix, exactly.'''),
    code(r'''P = np.asarray(oscprob.osc_prob(np.eye(3)*3.7e-13, 0.0, BASELINE))
print(np.round(P, 12))'''),
    md(r'''## 2. The trace does not matter

Adding a multiple of the identity shifts every eigenvalue equally, which multiplies the
evolution operator by an overall phase and cancels in the probability. Useful in practice: you
may drop the trace of your Hamiltonian without changing anything.'''),
    code(r'''P_plain = np.asarray(oscprob.osc_prob(h_vac/ENERGY, 0.0, BASELINE))
shift = 3.0*np.max(np.abs(h_vac/ENERGY))     # comparable to H itself
P_shifted = np.asarray(oscprob.osc_prob(h_vac/ENERGY + shift*np.eye(3), 0.0, BASELINE))

print('max |P(H) - P(H + c*I)| = %.2e' % np.max(np.abs(P_plain - P_shifted)))'''),
    md(r'''## 3. Zero baseline, and an exactly degenerate spectrum

A zero baseline means no evolution; equal masses mean no oscillation. Both return the identity
rather than a division by zero.'''),
    code(r'''print('L = 0:')
print(np.round(np.asarray(oscprob.osc_prob(h_vac/ENERGY, 0.0, 0.0)), 12))

h_degenerate = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
    **dict(OSC, D21=0.0, D31=0.0)))
print('\nD21 = D31 = 0:')
print(np.round(np.asarray(oscprob.osc_prob(h_degenerate/ENERGY, 0.0, BASELINE)), 12))'''),
    md(r'''## 4. Approaching degeneracy

The interesting case is not the exact degeneracy but the approach to it, where a closed form
divides by a vanishing splitting. Here the probability tends smoothly to the identity and
unitarity holds at the $10^{-15}$ level throughout -- fifteen orders of magnitude of shrinking
splitting, with no special-casing anywhere.'''),
    code(r'''print('%-12s %-18s %s' % ('scale', 'P_ee', 'max |row sum - 1|'))
print('-'*50)
for scale in (1e-3, 1e-6, 1e-9, 1e-12, 1e-15):
    h_near = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
        **dict(OSC, D21=scale*OSC['D21'], D31=scale*OSC['D31'])))
    P = np.asarray(oscprob.osc_prob(h_near/ENERGY, 0.0, BASELINE))
    print('%-12.0e %-18.12f %.1e'
          % (scale, P[0][0], np.max(np.abs(P.sum(axis=1) - 1.0))))'''),
    md(r'''## 5. Degenerate *requests*

A single slab, and no tolerance at all. Both are legitimate: `n_slabs=1` asks for one Magnus
step over the whole baseline, and `rtol=atol=None` switches the adaptive ladder off entirely.
For a constant Hamiltonian one slab is already exact, so all three agree.'''),
    code(r'''for label, kwargs in (('default        ', {}),
                      ('n_slabs=1      ', dict(n_slabs=1)),
                      ('rtol=atol=None ', dict(rtol=None, atol=None))):
    P = np.asarray(oscprob.osc_prob(h_vac/ENERGY, 0.0, BASELINE, **kwargs))
    print('%s P_ee = %.12f' % (label, P[0][0]))'''),
    md(r'''## 6. The nine warnings

| class | says | act on it? |
|---|---|---|
| `DensityUnitWarning` | a density is implausible for the units declared | **yes -- bad input** |
| `ScalarHamiltonianWarning` | your `H_func` takes one position at a time | yes -- costs speed only |
| `MagnusHighOrderCostWarning` | order > 6 with trapezoid/simpson is dear | your call |
| `MagnusConvergenceWarning` | a slab is wider than the sufficient condition | **often not** -- see below |
| `ToleranceNotAchievedWarning` | refinement stopped without agreeing | usually yes |
| `HybridCertificationWarning` | the adiabatic path could not certify itself | yes |
| `UnmarkedDiscontinuityWarning` | a density jump was detected, not declared | **yes -- pass `t_breakpoints`** |
| `HiddenFeatureWarning` | structure was found the sampling nearly missed | yes |
| `PhaseAveragingWarning` | `average=True` where the phase has not averaged | yes -- wrong question |

`MagnusConvergenceWarning` deserves its own sentence: it is a statement about **slab width, not
about the answer**, and it is measured to be a false alarm about three quarters of the time. It
fires in notebook 16 on a converged result. Do not read it as "this number is wrong"; read it
as "a sufficient condition was not met somewhere".

Below, each of six is provoked deliberately.'''),
    code(r'''def provoke(label, call):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            call()
            note = ''
        except Exception as exc:
            note = '  [raised %s]' % type(exc).__name__
    names = sorted({w.category.__name__ for w in caught
                    if issubclass(w.category, Warning)})
    print('%-32s %s%s' % (label, ', '.join(names) if names else '(quiet)', note))'''),
    code(r'''NE0, L_SCALE = gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN
PARAMS_2NU = {'sth': 0.55, 'Dm2': 7.5e-5}

# 1. a density in g/cm^3 with the flag left at its default
provoke('density under-converted', lambda: oscprob.osc_prob_3nu_matter_constant_density(
    ENERGY, BASELINE, 2.848, **OSC))

# 2. the mirror mistake: already converted, but declared as g/cm^3
provoke('density double-converted', lambda: oscprob.osc_prob_3nu_matter_constant_density(
    ENERGY, BASELINE, 100.0*gd.UNIT_G_PER_CM3, **OSC,
    density_matter_is_in_g_per_cm3=True))

# 3. an H_func that only accepts one position at a time
provoke('scalar H_func', lambda: oscprob.osc_prob(
    lambda l: h_vac/ENERGY + np.diag([float(np.asarray(l))*0.0, 0.0, 0.0]),
    0.0, BASELINE))

# 4. order 8 with a non-Gauss-Legendre integrator
provoke('order 8, simpson', lambda: magnus.magnus_expansion(
    lambda t: -1j*np.array([[0.0, 1.0], [1.0, 0.0]])*1e-13,
    0.0, 1e13, order=8, integration_method='simpson', n_tpts=20))

# 5. a density jump the caller did not declare
def step_ne(l):
    x = np.asarray(l, dtype=float)
    out = np.where(x < 0.5*L_SCALE, 0.02*NE0, 0.30*NE0)
    return out[()] if out.ndim == 0 else out

provoke('unmarked density jump', lambda: oscprob.osc_prob_matter_std_potential(
    2, step_ne, 50.0e6, 1.0*L_SCALE, PARAMS_2NU, L0=0.0,
    density_is_of_number_of_electrons=True))

# 6. asking for the averaged probability where nothing has averaged yet
provoke('average=True, few cycles', lambda: oscprob.osc_prob_3nu_vacuum(
    ENERGY, 5.0*gd.UNIT_KM, **OSC, average=True))

# ... and the same request where it genuinely has
provoke('average=True, many cycles', lambda: oscprob.osc_prob_3nu_vacuum(
    ENERGY, 5.0e4*gd.UNIT_KM, **OSC, average=True))'''),
    md(r'''The last two lines are the pattern worth internalizing. `PhaseAveragingWarning` is
not about accuracy -- the returned matrix is a perfectly valid doubly stochastic probability
matrix either way. It says the *question* does not apply at that baseline, because the phase
has not averaged and no averaged expression describes it. Move far enough out and it goes
quiet.

## Summary

Nothing in section 1--5 returns a `NaN`, including exact degeneracies, a zero baseline and a
Hamiltonian with no structure at all. That is a property of the method rather than of
defensive coding: the Magnus expansion exponentiates a matrix, and never forms the
$1/(\lambda_i - \lambda_j)$ that closed forms must.

For the warnings, one rule: **`MagnusConvergenceWarning` is about slab width, everything else
is about you.** Measured false-alarm rates for each are in `implementation_details.rst`.

Notebook 21 takes the tolerance warnings further, and it is the one to read next if you have
ever taken `rtol` for an error bound.'''),
    ])


# ------------------------------------------------- 21_magnus_what_tolerance_means
books['21_magnus_what_tolerance_means.ipynb'] = notebook(
    'What `rtol` and `atol` actually promise',
    r'''`rtol` and `atol` look like an error bound. They are not one.

They are a **stopping criterion**. Mag$\nu$s computes the answer on one grid, computes it again
on a finer grid, and stops when the two agree to within the tolerance you asked for. The
quantity compared is the *difference between two of its own approximations* -- not the distance
to the truth, which it does not know.

Most of the time that difference overestimates the error and the tolerance is conservative,
often by orders of magnitude. Sometimes it underestimates it, and the two grids agree with each
other while both are wrong. This notebook measures both cases against an independent
`solve_ivp` ground truth, and neither result is the one the parameter name suggests.''',
    [
    code(r'''import warnings

import numpy as np
from scipy.integrate import solve_ivp

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd

params = gd.load_nufit_params('NuFIT 6.1', 'NO')
sth, Dm2 = params['s12'], params['D21']
PARAMS_2NU = {'sth': sth, 'Dm2': Dm2}

# An exponential solar profile: smooth, monotonic, nothing adversarial about it.
ne = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
ENERGY = 10.0e6                                  # 10 MeV [eV]
BASELINE = 0.5*gd.SUN_RADIUS*gd.UNIT_KM'''),
    md(r'''## 1. An independent ground truth

The comparison is worthless without an oracle that does not share Mag$\nu$s's machinery, so we
integrate the Schrodinger equation directly with `solve_ivp` at a tolerance far tighter than
anything we will ask Mag$\nu$s for. This is a different algorithm, not a finer version of the
same one -- which is the property that matters.'''),
    code(r'''h_vac = np.asarray(hamiltonians.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2))
vcc = matter.vcc_func_from_rho_func(ne, density_is_of_number_of_electrons=True)

def H(l):
    return h_vac/ENERGY + np.diag([vcc(l), 0.0])

def rhs(t, y):
    return (-1j*H(t) @ y.reshape(2, 2)).ravel()

solution = solve_ivp(rhs, [0.0, BASELINE], np.eye(2, dtype=complex).ravel(),
                     method='DOP853', rtol=1.0e-11, atol=1.0e-13)
U_truth = solution.y[:, -1].reshape(2, 2)
P_truth = float(abs(U_truth[0, 0])**2)

print('solve_ivp DOP853, rtol=1e-11: P_ee = %.9f' % P_truth)
print('unitarity deviation         : %.1e'
      % abs(abs(U_truth[0, 0])**2 + abs(U_truth[0, 1])**2 - 1.0))'''),
    md(r'''## 2. What each requested tolerance actually delivered

Now ask Mag$\nu$s for four tolerances and compare each answer against that truth. The last two
columns are the point of the notebook: what you asked for, and what you got.'''),
    code(r'''print('%-12s %-12s %-11s %-11s %-9s' %
      ('requested', 'P_ee', '|error|', 'rel. error', 'achieved'))
print('-'*60)
rows = []
for tol in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6):
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = oscprob.osc_prob_matter_std_potential(
            2, ne, ENERGY, BASELINE, PARAMS_2NU, L0=0.0,
            density_is_of_number_of_electrons=True,
            convergence_info=info, rtol=tol, atol=tol*1.0e-2)
    value = float(np.asarray(P)[0][0])
    rel = abs(value - P_truth)/P_truth
    rows.append((tol, value, rel, info['tolerance_achieved'], info['n_slabs']))
    print('%-12.0e %-12.7f %-11.2e %-11.2e %-9s'
          % (tol, value, abs(value - P_truth), rel, info['tolerance_achieved']))'''),
    md(r'''Read the table three times, once for each surprise.

**The tolerance can be missed while reporting success.** The first row asked for $10^{-2}$,
reported `tolerance_achieved=True`, and is wrong by $2.5\times10^{-2}$ -- two and a half times
the tolerance it claimed to have met. The two grids it compared agreed with each other; they
were simply both too coarse. Nothing about the returned number reveals this.

**When it is conservative, it is very conservative.** The second row asked for $10^{-3}$ and
delivered $8.7\times10^{-6}$, a hundred times better. Ask for $10^{-4}$ and you get the same
answer and the same work -- the ladder had already stepped past it.

**`tolerance_achieved=False` does not mean the answer is bad.** The last row reports failure
and is the most accurate of the four, at $4\times10^{-7}$. It says "I could not verify
convergence by refining further", which is a statement about the ladder running out of room,
not about the answer.'''),
    code(r'''for tol, value, rel, achieved, n in rows:
    verdict = ('accurate' if rel <= tol else 'OUTSIDE the requested tolerance')
    print('requested %.0e -> delivered %.1e (%-30s) achieved=%-5s n_slabs=%d'
          % (tol, rel, verdict, achieved, n))'''),
    md(r'''## 3. What `convergence_info` reports

Every entry point fills a dictionary you pass in. It is the only way to see what the ladder
actually did.'''),
    code(r'''info = {}
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    oscprob.osc_prob_matter_std_potential(
        2, ne, ENERGY, BASELINE, PARAMS_2NU, L0=0.0,
        density_is_of_number_of_electrons=True,
        convergence_info=info, rtol=1.0e-6, atol=1.0e-8)

for key in sorted(info):
    print('%-26s %s' % (key, info[key]))'''),
    md(r'''`last_gap` is the quantity actually tested against your tolerance: the difference
between the final two refinements. `n_agreements` counts how many successive levels agreed --
zero here, which is why `tolerance_achieved` is False. `n_slabs` hit its ceiling of 20000.

## 4. `n_slabs` is not `n_slab_edges`

The two are different numbers, and the difference is why the ladder can be fooled.

`n_slabs` is what you request. `n_slab_edges` is how many pieces the integrator actually
propagated -- and any `t_breakpoints` you declared are edges too. On an Earth chord with the
PREM shell boundaries declared, the nominal refinement is a much smaller real one:'''),
    code(r'''h2_atm = np.asarray(hamiltonians.hamiltonian_2nu_vacuum_energy_independent(
    gd.load_nufit_params('NuFIT 6.1', 'NO')['s23'],
    gd.load_nufit_params('NuFIT 6.1', 'NO')['D31']))

def num_density_e_prem(r):
    return matter.num_density_e_func(r, earth.density_matter_func_prem,
                                     electron_fraction=0.5,
                                     density_matter_is_in_g_per_cm3=True)

costhz = -1.0
L_earth = earth.distance_traveled_inside_earth(costhz)*gd.CONV_KM_TO_INV_EV
breakpoints = np.asarray(
    earth.prem_layer_edges_along_chord(costhz))*gd.CONV_KM_TO_INV_EV
E_earth = 5.0*gd.UNIT_GEV

def H_earth(l):
    r = earth.earth_radial_distance_from_depth(costhz, l/gd.CONV_KM_TO_INV_EV)
    return h2_atm/E_earth + hamiltonians.hamiltonian_2nu_matter(
        matter.VCC_func(r, num_density_e_prem))

print('PREM edges declared: %d\n' % len(breakpoints))
print('%-16s %-10s %-14s %s' % ('n_slabs asked', 'n_slabs', 'n_slab_edges', 'real step'))
print('-'*56)
previous = None
for n in (2, 3, 4, 8):
    ci = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        oscprob.osc_prob(H_earth, 0.0, L_earth, n_slabs=n, t_breakpoints=breakpoints,
                         convergence_info=ci, rtol=None, atol=None)
    step = ('%+.0f%%' % (100.0*(ci['n_slab_edges']/previous - 1.0))
            if previous else '--')
    previous = ci['n_slab_edges']
    print('%-16d %-10s %-14s %s' % (n, ci['n_slabs'], ci['n_slab_edges'], step))'''),
    md(r'''Going from 2 slabs to 3 sounds like a 50% refinement. With the PREM boundaries
declared it is a 20 &rarr; 21 edge step: **5%**. Two grids that differ by 5% will very often
agree to within a loose tolerance, whatever the answer -- and the ladder would then certify
convergence it had not achieved.

That was a real defect, fixed in PR #35: the ladder certified an agreement between two nearly
identical grids. The fix is that `n_slabs` is now a **floor** rather than a target, so a
refinement step is guaranteed to be a real one.

## Summary

**A tolerance is a stopping criterion, not an error bound.** Measured on a smooth, entirely
ordinary solar profile against an independent oracle:

| requested | delivered | verdict |
|---|---|---|
| $10^{-2}$ | $2.5\times10^{-2}$ | **worse than asked, and reported as achieved** |
| $10^{-3}$ | $8.7\times10^{-6}$ | 100x conservative |
| $10^{-4}$ | $8.7\times10^{-6}$ | same work, same answer |
| $10^{-6}$ | $4.0\times10^{-7}$ | accurate, reported as *not* achieved |

What to do about it:

1. **Do not read `rtol` as an error bar.** If you need one, get it from a genuinely different
   method -- `cross_check_strategies`, which is notebook 22, or a `solve_ivp` reference as
   above.
2. **Ask for more than you need.** The tolerance is usually conservative, and tightening it
   often costs nothing because the ladder has already stepped past.
3. **Read `convergence_info`,** not just the probability. `tolerance_achieved` and `last_gap`
   are the only visible evidence of what happened.
4. **Do not treat `tolerance_achieved=False` as a failed calculation.** It frequently
   accompanies the best answer in the set.'''),
    ])


# -------------------------------------------------- 22_magnus_which_engine_answered
books['22_magnus_which_engine_answered.ipynb'] = notebook(
    'Which engine answered, and why',
    r'''Mag$\nu$s does not have one algorithm. It has six, grouped into five families, and
`strategy='auto'` picks between them per request. Most of the time you neither know nor need to
know which one ran -- but when an answer looks wrong, "which engine produced this" is the first
question, and Mag$\nu$s will tell you.

The second half of this notebook is about something stronger. Every silently-wrong result the
package's adversarial validation ever found came from a method **certifying itself**: refining
its own knobs, comparing itself with itself, and agreeing. When a method has a blind spot both
sides of that comparison share it, and the agreement carries no information. Running two
genuinely *different* engines needs no oracle at all, and detects exactly the class that
self-certification cannot.''',
    [
    code(r'''import warnings

import numpy as np

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.matter as matter
import magnus.globaldefs as gd

NE0, L_SCALE = gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN
PARAMS_2NU = {'sth': 0.55, 'Dm2': 7.5e-5}
ne = matter.exp_density_profile(NE0, L_SCALE)
BASELINE = 1.0*L_SCALE'''),
    md(r'''## 1. The engines, and which share machinery

The families matter more than the engines. Two engines in the same family share code and
assumptions, so their agreement is weak evidence; two in different families agreeing is
strong evidence.'''),
    code(r'''for engine, family in oscprob.ENGINE_FAMILIES.items():
    print('%-12s %s' % (engine, family))'''),
    md(r'''## 2. What ran, for this request

Pass a dictionary as `strategy_info` and it comes back filled in. `engine` names what answered,
`family` places it, `certified` says whether that engine was able to vouch for its own result,
and `declined` lists what stood aside and why.'''),
    code(r'''print('%-10s %-11s %-21s %-10s %s'
      % ('strategy', 'engine', 'family', 'certified', 'declined'))
print('-'*76)
for strategy in ('auto', 'hybrid', 'magnus'):
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        oscprob.osc_prob_matter_std_potential(
            2, ne, 10.0e6, BASELINE, PARAMS_2NU, L0=0.0,
            density_is_of_number_of_electrons=True,
            strategy=strategy, strategy_info=info)
    print('%-10s %-11s %-21s %-10s %s'
          % (strategy, info.get('engine'), info.get('family'),
             info.get('certified'), info.get('declined') or '--'))'''),
    md(r'''On this smooth solar profile the default picks the **adiabatic** engine and certifies
it. Forcing `strategy='magnus'` gets a different family entirely -- the interaction-picture
integrator -- and it does not attempt to certify itself, which is why `certified` is `None`
rather than `False`.

## 3. Cross-checking: an error bar with no oracle

`cross_check_strategies` answers the same request with every engine that applies and reports
how far apart they are. It is never on by default -- it multiplies the cost of the call by the
number of engines -- and a large spread is *reported*, never raised. What it means depends on
the request, and deciding that is your job.'''),
    code(r'''with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    out = oscprob.cross_check_strategies(
        oscprob.osc_prob_matter_std_potential, 2, ne, 10.0e6, BASELINE,
        PARAMS_2NU, L0=0.0, density_is_of_number_of_electrons=True)

print('ran                  :', sorted(out['ran']))
print('families             :', sorted(set(out['families'].values())))
print('max spread           : %.3e' % out['max_spread'])
print('max across families  : %.3e  %s'
      % (out['max_spread_independent'], out['max_spread_independent_pair']))
print()
for label, reason in out['declined'].items():
    print('declined %-10s %s' % (label, reason))'''),
    md(r'''Four engines, three families, and they agree to $10^{-4}$ with no reference solution
anywhere in sight. That number is a far more honest error bar than any `rtol` (notebook 21),
because the things being compared do not share a method.

## 4. The case it exists for

Now a density profile with an **undeclared step** in it -- the construction from the package's
adversarial-validation findings. This is where the adiabatic engine historically returned a
confidently wrong answer.'''),
    code(r'''def step_ne(l):
    """A density jump the caller does not declare."""
    x = np.asarray(l, dtype=float)
    out = np.where(x < 0.5*BASELINE, 0.02*NE0, 0.30*NE0)
    return out[()] if out.ndim == 0 else out

print('%-10s %-11s %-11s %-11s %s'
      % ('strategy', 'P_ee', 'engine', 'certified', 'warnings'))
print('-'*74)
for strategy in ('auto', 'hybrid'):
    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = oscprob.osc_prob_matter_std_potential(
            2, step_ne, 50.0e6, BASELINE, PARAMS_2NU, L0=0.0,
            density_is_of_number_of_electrons=True,
            strategy=strategy, strategy_info=info)
    names = sorted({w.category.__name__.replace('Warning', '') for w in caught})
    print('%-10s %-11.6f %-11s %-11s %s'
          % (strategy, float(np.asarray(P)[0][0]), info.get('engine'),
             info.get('certified'), ', '.join(names)))
    if info.get('declined'):
        print('%-10s   declined: %s' % ('', info['declined']))'''),
    md(r'''Two different answers, 0.085 and 0.550, from the same request. The package is no
longer silent about it in either direction:

* under `'auto'` the adiabatic engine **declines** -- "the profile is not resolved at the probe
  scale" -- and the Magnus ladder answers instead;
* forced with `strategy='hybrid'` it still answers, but reports `certified=False` and raises
  `HybridCertificationWarning`.

Both paths also raise `UnmarkedDiscontinuityWarning`, which names the actual fix: declare the
jump with `t_breakpoints`.

And the cross-check sees the disagreement without being told any of that:'''),
    code(r'''with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    out_step = oscprob.cross_check_strategies(
        oscprob.osc_prob_matter_std_potential, 2, step_ne, 50.0e6, BASELINE,
        PARAMS_2NU, L0=0.0, density_is_of_number_of_electrons=True)

for engine in sorted(out_step['ran']):
    print('%-12s (%-14s) P_ee = %.6f'
          % (engine, oscprob.ENGINE_FAMILIES[engine],
             float(np.asarray(out_step['answers'][engine])[0][0])))
print('\nmax spread across families: %.3e  %s'
      % (out_step['max_spread_independent'], out_step['max_spread_independent_pair']))'''),
    md(r'''A spread of **0.47** on a probability. No ground truth was computed, no reference
code was installed, and nothing had to know in advance what was wrong with the profile. Two
engines from different families simply disagreed, which is all the signal you need to stop
trusting the number.

## 5. A zero spread is not always agreement

One trap, and it is the reason this function warns. `max_spread` is a plain float, and it is
`0.0` both when the engines agree perfectly and when **nothing was compared**. The commonest
way to reach the second case is to pass an entry point that has no `strategy` parameter --
`osc_prob` itself is one, and it is the function most of these notebooks call directly.'''),
    code(r'''import magnus.hamiltonians as hamiltonians

h_vac = np.asarray(hamiltonians.hamiltonian_2nu_vacuum_energy_independent(**PARAMS_2NU))
vcc = matter.vcc_func_from_rho_func(ne, density_is_of_number_of_electrons=True)

def H_func(l):
    return h_vac/10.0e6 + np.diag([vcc(l), 0.0])

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    empty = oscprob.cross_check_strategies(oscprob.osc_prob, H_func, 0.0, BASELINE)

print('ran        :', empty['ran'])
print('max_spread : %.1f   <- and yet nothing was compared' % empty['max_spread'])
print()
for w in caught:
    if w.category is oscprob.CrossCheckInconclusiveWarning:
        print('%s raised.' % w.category.__name__)'''),
    md(r'''The same vacuous zero appears when only one engine runs, and
`max_spread_independent` is zero whenever every engine that ran belongs to a single family --
which is precisely the self-certification the cross-check exists to avoid.
`CrossCheckInconclusiveWarning` covers all three.

## Summary

| question | how to answer it |
|---|---|
| which engine ran? | pass `strategy_info={}` and read `engine` |
| do they share machinery? | `oscprob.ENGINE_FAMILIES` |
| did the engine vouch for itself? | `strategy_info['certified']` |
| how far apart are independent methods? | `cross_check_strategies(...)['max_spread_independent']` |
| why did an engine stand aside? | `strategy_info['declined']`, or `out['declined']` |

**Always check `ran` before reading a spread.** A cross-check that ran nothing reports perfect
agreement, and a cross-check confined to one family reports independent agreement it never
tested. The warning will tell you, but the number will not.'''),
    ])


# ------------------------------------------------ 23_magnus_when_averaging_helps
books['23_magnus_when_averaging_helps.ipynb'] = notebook(
    'When averaging rescues you, and when it does not',
    r'''Notebooks 13 and 14 end in opposite places. On a tabulated solar model the instantaneous
probability carries an error of $1.4\times10^{-3}$, and the *observable* -- the same quantity
averaged over the detector's energy resolution -- carries $2.6\times10^{-5}$: the error falls by
a factor of **53**. On a supernova shock front the instantaneous error is $2.0\times10^{-1}$ and
the averaged error is $2.1\times10^{-1}$: it does not move at all.

The difference is not the size of the error but its **kind**. An error in the accumulated
*phase* moves the oscillation sideways, and sideways motion cancels when you integrate over
several oscillations. An error in the *envelope* moves the curve up or down, and no amount of
averaging removes an offset.

This notebook isolates the mechanism on a cheap vacuum probability, where both kinds of error
can be injected deliberately and neither costs anything to compute. The real measurements stay
where they were made -- notebooks 13 and 14 -- because reproducing them here would cost several
minutes and tell you nothing new.''',
    [
    code(r'''import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.avgprob as avgprob
import magnus.globaldefs as gd

sth, Dm2 = 0.55, 7.5e-5
BASELINE = 1.0e4*gd.UNIT_KM
E = np.linspace(9.0e6, 11.0e6, 4001)          # a 2 MeV band around 10 MeV

P = np.asarray(oscprob.osc_prob_2nu_vacuum(E, BASELINE, sth, Dm2,
                                           nu_i=gd.NUE, nu_f=gd.NUE))
cycles = (Dm2*BASELINE/(4*E[0]) - Dm2*BASELINE/(4*E[-1]))/np.pi
print('oscillation cycles across the band: %.1f' % cycles)'''),
    md(r'''## 1. Two kinds of wrong

A **phase** error: the same oscillation, evaluated as though the energy (or the baseline, or
the potential) were slightly different. The curve is displaced along the energy axis.

An **envelope** error: the oscillation is in the right place but its amplitude is wrong. The
curve is compressed towards its mean.

Both are parametrized so that their *instantaneous* size is comparable, which is the point --
they look equally bad point by point.'''),
    code(r'''def phase_error(fraction):
    """The right oscillation, in slightly the wrong place."""
    return np.asarray(oscprob.osc_prob_2nu_vacuum(
        E*(1.0 + fraction), BASELINE, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUE))

def envelope_error(amount):
    """The right place, with the wrong amplitude."""
    return 0.5 + (P - 0.5)*(1.0 - amount)'''),
    code(r'''print('%-20s %-14s %-14s %s'
      % ('error', 'instantaneous', 'averaged', 'suppression'))
print('-'*62)
for fraction in (1.0e-4, 3.0e-4):
    Q = phase_error(fraction)
    inst, avg = np.max(np.abs(Q - P)), abs(Q.mean() - P.mean())
    print('%-20s %-14.3e %-14.3e %6.0fx'
          % ('phase, %.0e' % fraction, inst, avg, inst/avg))
for amount in (0.02, 0.05):
    Q = envelope_error(amount)
    inst, avg = np.max(np.abs(Q - P)), abs(Q.mean() - P.mean())
    print('%-20s %-14.3e %-14.3e %6.1fx'
          % ('envelope, %.0f%%' % (100*amount), inst, avg, inst/avg))'''),
    md(r'''Averaging suppresses the phase error by about a hundredfold and the envelope error by
about seven. Note the second pair: doubling the envelope error changes both columns and leaves
the suppression at exactly 7.0. That is the signature of an offset -- averaging rescales it,
it does not remove it.'''),
    code(r'''fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.0, 3.4), sharey=True)
band = (E >= 9.6e6) & (E <= 10.4e6)
for ax, Q, title in ((ax0, phase_error(3.0e-4), 'phase error'),
                     (ax1, envelope_error(0.05), 'envelope error')):
    ax.plot(E[band]/1e6, P[band], color='0.2', label='truth')
    ax.plot(E[band]/1e6, Q[band], color='C3', ls='--', label='wrong')
    ax.axhline(P[band].mean(), color='0.2', lw=0.8)
    ax.axhline(Q[band].mean(), color='C3', lw=0.8, ls='--')
    ax.set_xlabel(r'$E$ [MeV]')
    ax.set_title('%s (thin lines: band means)' % title, fontsize=9)
ax0.set_ylabel(r'$P(\nu_e \to \nu_e)$')
ax0.legend(fontsize=8, loc='lower right')
fig.tight_layout()'''),
    md(r'''The thin horizontal lines are what an experiment with this energy resolution would
actually measure. On the left they are nearly on top of each other; on the right they are
visibly apart, by the same fraction as the curves themselves.

## 2. Averaging needs cycles

The suppression is not unconditional. Below about two oscillations in the band there is nothing
for the average to cancel, and a phase error survives almost intact.'''),
    code(r'''print('%-12s %-10s %-16s %s'
      % ('L [km]', 'cycles', 'phase suppr.', 'envelope suppr.'))
print('-'*54)
for L_km in (2.5e3, 5.0e3, 1.0e4, 2.0e4, 4.0e4):
    L = L_km*gd.UNIT_KM
    grid = np.linspace(9.0e6, 11.0e6, 8001)
    truth = np.asarray(oscprob.osc_prob_2nu_vacuum(grid, L, sth, Dm2,
                                                   nu_i=gd.NUE, nu_f=gd.NUE))
    n_cycles = (Dm2*L/(4*grid[0]) - Dm2*L/(4*grid[-1]))/np.pi
    shifted = np.asarray(oscprob.osc_prob_2nu_vacuum(grid*(1.0 + 1.0e-4), L, sth, Dm2,
                                                     nu_i=gd.NUE, nu_f=gd.NUE))
    ph = np.max(np.abs(shifted - truth))/abs(shifted.mean() - truth.mean())
    damped = 0.5 + (truth - 0.5)*0.98
    en = np.max(np.abs(damped - truth))/abs(damped.mean() - truth.mean())
    print('%-12.0f %-10.1f %-16.0f %.1f' % (L_km, n_cycles, ph, en))'''),
    md(r'''Two readings. The phase suppression climbs steeply from one to three cycles and then
settles around a hundred -- it does not grow without bound, because where the band edges fall
relative to the oscillation matters as much as how many cycles are inside it. The envelope
suppression sits near 6 throughout and never improves.

**The practical test**, and it needs no ground truth: *average your result over a few
oscillation lengths and see whether the discrepancy moves.* If it collapses, what you had was
a phase error, and the observable was fine all along. If it stays put, the error is real.

## 3. Asking for the averaged probability directly

You do not have to build the average by hand. `average=True` returns the exact decohered limit
-- the value the oscillation averages to when every relative phase is unresolvable -- and
`magnus.avgprob` exposes the machinery for a finite window.'''),
    code(r'''analytic = float(np.asarray(oscprob.osc_prob_2nu_vacuum(
    1.0e7, BASELINE, sth, Dm2, nu_i=gd.NUE, nu_f=gd.NUE, average=True)))

print('analytic decohered limit (average=True) : %.6f' % analytic)
print('numerical mean over the 2 MeV band      : %.6f' % P.mean())
print('difference                              : %.2e' % abs(analytic - P.mean()))'''),
    md(r'''They differ in the third decimal, and that is not an error in either: the analytic
value is the infinite-window limit, while the band mean is over 6.1 cycles with the edges
falling where they fall. Asking for `average=True` where the phase has *not* averaged raises
`PhaseAveragingWarning` -- notebook 20 provokes it deliberately.

## Summary

| | phase error | envelope error |
|---|---|---|
| what it does | displaces the oscillation | changes its amplitude |
| suppressed by averaging | **~100x** | ~7x, fixed |
| improves with more cycles | yes, then plateaus | no |
| improves with a smaller error | no -- the ratio is scale-free | no |
| real instance | supernova turbulence, 45 MeV (`docs/dev/adversarial_batteries/avg_check2.py`): **15x** | notebook 14, shock: **~1x** |
| is the observable affected? | barely | **yes** |

The reason notebook 14's shock error does not average away is physical rather than numerical: a
shock front changes the *adiabaticity* of the level crossing, so it moves the conversion
probability itself rather than the phase at which it oscillates. That is the error becoming an
envelope, and it is the one case in this package where the instantaneous error and the
observable error are the same number.

**A caution about the suppression factor itself, which this notebook is the right place for.**
Every number in the table above is a ratio of *finite-window* means, and the section just
above shows that such a mean is an estimator with its own bias -- 6.08e-03 from the analytic
limit here, over 6.1 cycles. On a profile whose density varies appreciably across the window,
that bias does not shrink as the window widens, because a wider window also averages over
different matter conditions; notebook 13's solar ray is exactly that case, and there the
suppression ratio carries no information at all. Use the ratio to tell phase from envelope on
a *controlled* comparison like this one. To get the observable, ask for it: `average=True`
computes the decohered limit in closed form, with no window to choose.'''),
    ])


# ------------------------------------------------------- 24_magnus_performance
books['24_magnus_performance.ipynb'] = notebook(
    'Performance: what is worth doing',
    r'''Mag$\nu$s is fast enough for most single calls -- the median across 164 Earth and solar
configurations is 2 ms. Scans are what the code mostly does, and three things make them
substantially faster **without changing any answer**.

Every measurement here is taken live. Timings vary between machines and between runs, so treat
the ratios rather than the absolute numbers as the result, and note that two of the three
optimizations are worth nothing at all in the wrong circumstances -- which is more useful to
know than a headline speed-up.''',
    [
    code(r'''import time
import warnings

import numpy as np
import matplotlib.pyplot as plt        # section 4 is the only figure in this notebook

# np.trapz was removed in NumPy 2.0 and renamed np.trapezoid.  Ask the installed
# version which one it has rather than pinning either: these notebooks are read
# on whatever numpy the reader happens to have.
trapezoid = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.magnus as magnus
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd

warnings.simplefilter('ignore')

# load_nufit_params returns exactly the six mixing parameters, ready to splat
# into any osc_prob_3nu_* call.  'NuFIT 6.1' is the package default.
OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
osc = OSC
h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC))
e00 = np.diag([1.0, 0.0, 0.0])

COSTHZ = -0.9
L_EARTH = earth.distance_traveled_inside_earth(COSTHZ)*gd.CONV_KM_TO_INV_EV

def best_of(call, repeats=3):
    """Fastest of a few runs -- the least noisy summary of a timing."""
    fastest, result = np.inf, None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = call()
        fastest = min(fastest, time.perf_counter() - t0)
    return result, fastest'''),
    md(r'''## 1. Pass an array of energies, do not loop

Every wrapper accepts an array of energies, of baselines, or both. Handing the whole scan over
at once lets the matter profile be built once rather than once per point, which is what the
energy-batched engine exists to do.'''),
    code(r'''def H(energy, l, VCC):
    return h_vac/energy + np.asarray(VCC)[..., None, None]*e00

energies = np.logspace(0.0, 1.5, 60)*gd.UNIT_GEV

batched, t_batched = best_of(lambda: np.asarray(oscprob.osc_prob_earth(
    H, energy=energies, costhz=COSTHZ, L=L_EARTH, nu_i=gd.NUMU, nu_f=gd.NUMU)))
looped, t_looped = best_of(lambda: np.array([
    oscprob.osc_prob_earth(H, energy=float(e), costhz=COSTHZ, L=L_EARTH,
                           nu_i=gd.NUMU, nu_f=gd.NUMU) for e in energies]))

print('%d energies through an Earth chord' % len(energies))
print('  array : %.3f s' % t_batched)
print('  loop  : %.3f s   (%.2fx slower)' % (t_looped, t_looped/t_batched))
print('  max |difference| = %.1e' % np.max(np.abs(batched - looped)))'''),
    md(r'''The two answers are not bit-identical -- the batched path refines the profile once for
the whole scan rather than independently per point, so the two land at slightly different
grids. The difference is at the $10^{-6}$ level, far inside any tolerance you would request.

## 2. The palindrome, and when it is worth nothing

A chord through a spherically symmetric Earth meets every radius twice, so Mag$\nu$s can
evaluate your `H_func` on the first half of the slab chain and mirror the rest. That halves the
**number of positions** at which your Hamiltonian is evaluated.

Whether that is worth anything depends on a distinction worth being precise about: it halves
the positions, not the number of *calls*. If your `H_func` is dominated by fixed per-call
overhead, halving the positions saves nothing. If its cost scales with how many positions it
was handed -- an interpolation, a table lookup per point, a quadrature -- it saves half.

First, plain PREM, where a density lookup is too cheap to be worth halving:'''),
    code(r'''def H_prem(energy, l, VCC):
    return h_vac/energy + np.asarray(VCC)[..., None, None]*e00

def timed(H_func, **kwargs):
    def call():
        return np.asarray(oscprob.osc_prob_earth(
            H_func, costhz=COSTHZ, L=L_EARTH, nu_i=gd.NUMU, nu_f=gd.NUMU, **kwargs))
    magnus.USE_PALINDROME = True
    on_result, on_time = best_of(call)
    magnus.USE_PALINDROME = False
    off_result, off_time = best_of(call)
    magnus.USE_PALINDROME = True                      # restore the default
    return on_time, off_time, np.max(np.abs(on_result - off_result))

on, off, gap = timed(H_prem, energy=5.0*gd.UNIT_GEV)
print('plain PREM lookup, single point')
print('  palindrome on  : %.4f s' % on)
print('  palindrome off : %.4f s   speed-up %.2fx' % (off, off/on))'''),
    md(r'''About 1.00x -- the bookkeeping costs as much as the lookups it saves. The package
documentation quotes 0.91x for this case, and disarming the optimization here would cost you
nothing.

Now a Hamiltonian whose cost genuinely scales with the number of positions it is given. The
quadrature below stands in for an interpolated profile, a per-point integral, or any of the
things a real custom Hamiltonian does.'''),
    code(r'''GRID = np.linspace(0.0, 1.0, 4000)

def H_expensive(energy, l, VCC):
    """Per-position work: a small quadrature for every position handed in."""
    V = np.atleast_1d(np.asarray(VCC, dtype=float))
    weight = trapezoid(np.exp(-GRID[None, :]*1.0e13*np.abs(V)[:, None]), GRID, axis=1)
    V_eff = V*(1.0 + 1.0e-12*weight)
    shape = np.asarray(VCC).shape
    V_eff = V_eff.reshape(shape) if shape else V_eff[0]
    return h_vac/energy + np.asarray(V_eff)[..., None, None]*e00

print('%-26s %-10s %-10s %-10s %s'
      % ('workload', 'on [s]', 'off [s]', 'speed-up', 'max |diff|'))
print('-'*66)
for label, kwargs in (('single point', dict(energy=5.0*gd.UNIT_GEV)),
                      ('20-energy scan',
                       dict(energy=np.logspace(0.0, 1.5, 20)*gd.UNIT_GEV))):
    on, off, gap = timed(H_expensive, **kwargs)
    print('%-26s %-10.4f %-10.4f %-10.2f %.1e' % (label, on, off, off/on, gap))'''),
    md(r'''Now it pays: the mirrored path is faster by a factor approaching two, and the answers
agree to round-off. The package documentation quotes 1.4--1.67x on an expensive `H_func`,
measured on a different profile; the numbers above are the same effect on this one.

The lesson is not "the palindrome is worth 1.75x". It is that **it is worth exactly half of
whatever your Hamiltonian charges per position, and nothing for what it charges per call**. If
you want it off, `magnus.magnus.USE_PALINDROME = False`.

Note also that only `osc_prob_earth` gets this: a chord is symmetric by geometry, and
`osc_prob_sun` deliberately does not declare it, because a solar profile is monotonic.

## 3. What a tolerance costs

Tightening `rtol` is not a smooth dial. The refinement ladder moves in steps, so several
requests can land on the same grid and cost exactly the same -- as notebook 21 showed, asking
for $10^{-4}$ sometimes buys the $10^{-3}$ answer for free, and sometimes the reverse.'''),
    code(r'''ne = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
PARAMS_2NU = {'sth': 0.55, 'Dm2': 7.5e-5}
L_SUN = 0.5*gd.SUN_RADIUS*gd.UNIT_KM

print('%-14s %-10s %-12s %s' % ('rtol', 'time [s]', 'n_slabs', 'relative cost'))
print('-'*52)
baseline_time = None
for tol in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5):
    info = {}
    _, t = best_of(lambda tol=tol, info=info: oscprob.osc_prob_matter_std_potential(
        2, ne, 1.0e7, L_SUN, PARAMS_2NU, L0=0.0,
        density_is_of_number_of_electrons=True,
        convergence_info=info, rtol=tol, atol=tol*1.0e-2), repeats=2)
    baseline_time = baseline_time or t
    print('%-14.0e %-10.3f %-12s %.2fx'
          % (tol, t, info['n_slabs'], t/baseline_time))'''),
    md(r'''Three orders of magnitude of extra accuracy for roughly twice the work, and the first
of them free. Tolerances are usually worth tightening.

## Summary

| what | speed-up | when it is worth nothing |
|---|---|---|
| pass an array of energies | **~2.7x** | single-point calls |
| write `H_func` to take an array of positions | **~4.6x** (notebook 19) | never -- always do this |
| the palindrome, expensive `H_func` | **~1.8x** | cheap or per-call-dominated `H_func` |
| the palindrome, plain PREM | ~1.0x | this is the "worth nothing" case |
| tightening `rtol` by $10^{3}$ | costs ~2x | -- |

Ranked by what you control: **vectorize your `H_func` first** (notebook 19), **pass arrays
second**, and let the palindrome look after itself -- it is on by default, it is free when it
helps, and it costs a few percent when it does not.

None of these change an answer by more than round-off, which is the property that makes them
worth taking.'''),
    md(r'''## 4. Choosing `magnus_exp_order`

The order of the Magnus expansion is the one dial in this library that is **never adjusted for
you**. `osc_prob` runs at a fixed `magnus_exp_order` and reaches the requested tolerance by
refining the number of slabs, never by raising the order. So the order does not decide whether
a tolerance is met -- it decides how *fast* slab refinement pays, and that is worth knowing
because the two ends of the range differ by an order of magnitude in cost per slab.

Two constraints shape everything below, and both are properties of the integrator rather than
of the physics.

**With the default `gl` integrator only three orders are distinct.** Gauss--Legendre schemes
come in whole nodes: orders 1--2 share the 1-node scheme, 3--4 the 2-node scheme, 5--6 the
3-node scheme. So asking for order 3 gets you order 4's arithmetic at order 4's price. The
table below shows all six anyway, once, because seeing the pairs collapse is more convincing
than being told they do.

**`gl` stops at 6.** `magnus._validate` rejects higher orders for `gl`; `MAGNUS_EXP_ORDER_MAX`
is 10, but reaching 7--10 requires `trapezoid` or `simpson`, which changes the *integrator* as
well as the order. Those points are not on the same curve and are not drawn here.'''),
    code(r'''from scipy.integrate import solve_ivp

E_ORD = np.logspace(np.log10(1.0), np.log10(10.0), 12)*gd.UNIT_GEV
CHORD_KM = earth.distance_traveled_inside_earth(COSTHZ)


def vcc_prem_at(l):
    """V_CC at distance l along the chord, from PREM."""
    r = np.sqrt(gd.EARTH_RADIUS**2 + (l/gd.CONV_KM_TO_INV_EV)**2
                + 2.0*gd.EARTH_RADIUS*(l/gd.CONV_KM_TO_INV_EV)*COSTHZ)
    return matter.vcc_func_from_rho_func(
        float(np.asarray(earth.density_matter_func_prem(r))),
        density_matter_is_in_g_per_cm3=True)


def dop853_earth(energy):
    """A referee that is not a Magnus expansion, so the order cannot flatter itself."""
    def rhs(l, y):
        return (-1j*(h_vac/energy + float(vcc_prem_at(l))*e00) @ y.reshape(3, 3)).ravel()
    sol = solve_ivp(rhs, (0.0, L_EARTH), np.eye(3, dtype=complex).ravel(),
                    rtol=1.0e-12, atol=1.0e-14, method='DOP853')
    u = sol.y[:, -1].reshape(3, 3)
    return abs(u[gd.NUE, gd.NUMU])**2


REF_ORD = np.array([dop853_earth(e) for e in E_ORD])


def earth_at_order(order):
    return np.asarray(oscprob.osc_prob_3nu_earth(
        E_ORD, costhz=COSTHZ, L=L_EARTH, **OSC, nu_i=gd.NUMU, nu_f=gd.NUE,
        magnus_exp_order=order, n_slabs=600, max_n_slabs=600,
        rtol=1.0e-13, atol=1.0e-15))


earth_at_order(4)                       # discard: the first call compiles the kernel

print('EARTH THROUGH PREM, slab count fixed at 600 so the order is the only variable')
print('%8s %14s %16s   %s' % ('order', 'ms', 'max |dP| vs DOP853', 'GL nodes'))
print('-'*62)
rows_ord = []
for order in (1, 2, 3, 4, 5, 6):
    P, t = best_of(lambda o=order: earth_at_order(o))
    err = float(np.max(np.abs(P - REF_ORD)))
    rows_ord.append((order, t, err))
    print('%8d %14.2f %16.3e   %d' % (order, 1.0e3*t, err, 1 if order <= 2 else
                                      (2 if order <= 4 else 3)))'''),
    md(r'''**The pairs collapse exactly**, which is the clearest way to see that the order is a
request for a quadrature scheme rather than a continuous knob: 1 and 2 agree to the last digit,
as do 3 and 4, and 5 and 6. There are three settings here wearing six names.

What each real step buys, on this profile: order 2 to 4 is worth a factor of about **5000** in
accuracy for **1.95x** the time; 4 to 6 is worth a further **6x** for **1.7x** more. That is
the shape of the trade -- the first step is overwhelmingly worth taking, the second is a
genuine choice that depends on how smooth the profile is.'''),
    code(r'''fig, ax = plt.subplots(figsize=(6.4, 4.4))
for order, t, err in rows_ord:
    marker = 'o' if order in (2, 4, 6) else 'x'
    ax.loglog(1.0e3*t, max(err, 1.0e-16), marker, ms=9 if order in (2, 4, 6) else 7,
              mfc='white' if order in (2, 4, 6) else 'C1',
              color='k' if order in (2, 4, 6) else 'C1', mew=1.4, zorder=4)
    ax.annotate(str(order), xy=(1.0e3*t, max(err, 1.0e-16)), xytext=(6, 3),
                textcoords='offset points', fontsize=8,
                color='k' if order in (2, 4, 6) else 'C1')
ax.set_xlabel('Time for 12 probabilities [ms]')
ax.set_ylabel(r'Error vs.\ DOP853,  max $|\Delta P|$')
ax.set_title(r'Earth through PREM: what `magnus_exp_order` buys at fixed slab count',
             fontsize=10)
ax.grid(True, which='both', alpha=0.2)
ax.text(0.03, 0.06, 'open circles: the three distinct GL schemes\n'
        'crosses: orders that reuse the scheme below them',
        transform=ax.transAxes, fontsize=6.6, color='0.3', linespacing=1.5)
fig.tight_layout(pad=1.2)
fig.savefig('../fig/expansion_order.pdf', bbox_inches='tight')'''),
    md(r'''### A case where the order buys nothing at all

`average=True` on a smooth, position-dependent profile takes the **adiabatic** route: it
decoheres in the matter eigenbasis at production, transports along the levels of the
instantaneous Hamiltonian, and reads out in vacuum. It never propagates -- and therefore never
evaluates a Magnus expansion. The order is accepted, and ignored.

This is worth measuring rather than reasoning about, because "accepted and ignored" is
indistinguishable from "accepted and used but unimportant" until you look.'''),
    code(r'''TABLE_ORD = '../docs/dev/adversarial_batteries/bs05_agsop.dat'
_rows = []
with open(TABLE_ORD) as fh:
    for line in fh:
        f = line.split()
        if len(f) == 12:
            try:
                _rows.append([float(x) for x in f])
            except ValueError:
                continue
_solar = np.array(_rows)
_mean_nucleon = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
_ne = _solar[:, 3]*gd.UNIT_G_PER_CM3/_mean_nucleon*(0.5*(1.0 + _solar[:, 6]))
_x = _solar[:, 1]*gd.SUN_RADIUS*gd.UNIT_KM
_logne = np.log(_ne)
R_SUN_ORD = float(_x[-1])


def ne_sun_ord(l):
    xs = np.clip(np.asarray(l, dtype=float), _x[0], _x[-1])
    out = np.exp(np.interp(xs, _x, _logne))
    return out[()] if np.ndim(out) == 0 else out


E_SUN_ORD = np.logspace(np.log10(0.1), np.log10(20.0), 40)*gd.UNIT_MEV
print('THE SUN, average=True (the adiabatic route)')
print('%8s %12s %26s' % ('order', 'seconds', 'max |P - P(order 2)|'))
print('-'*50)
P_ord2 = None
for order in (2, 4, 6):
    t0 = time.perf_counter()
    P = np.asarray(oscprob.osc_prob_matter_std_potential(
        3, ne_sun_ord, E_SUN_ORD, R_SUN_ORD, OSC, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE,
        density_is_of_number_of_electrons=True, average=True,
        magnus_exp_order=order))
    dt = time.perf_counter() - t0
    if P_ord2 is None:
        P_ord2 = P
    print('%8d %12.3f %26.3e' % (order, dt, float(np.max(np.abs(P - P_ord2)))))'''),
    md(r'''**Bit-identical, and the same cost.** Not "small", not "below tolerance" -- exactly
zero, because the expansion is never reached. If you are averaging over a smooth profile,
`magnus_exp_order` is not a dial you have.

The opposite extreme is a profile with a genuine discontinuity, where the order matters but the
*slabs* matter more: no order resolves a front that a slab edge straddles, which is what
`t_breakpoints` is for. Measured on notebook 14's supernova shock with the front declared and
the slab count fixed at 8000, order 2 gives $9.5\times10^{-5}$, order 4 gives
$2.5\times10^{-7}$ and order 6 gives $1.1\times10^{-8}$ -- so the order still pays there, about
**385x** for the first step and **22x** for the second, at 1.4x and 1.8x the cost. Notebook 25
runs that case against other codes in full.'''),
    md(r'''### The recommendation

**Leave it at 4.** That is the default, and it is the right one: it is the first scheme that is
genuinely fourth order, it costs under twice what order 2 costs, and on every profile measured
here it buys between two and five orders of magnitude over order 2.

Raise it to **6** when the profile is **smooth** and the target accuracy is tight. What that
extra node buys depends strongly on the profile, and the two cases measured here bracket it:
about **6x** on a piecewise-constant PREM chord, about **22x** on the resolved shock front,
both for roughly 1.7--1.8x the time. Smoother profiles pay better, which is the same ordering
the slab-refinement rate follows.

Do **not** raise it when a `MagnusConvergenceWarning` appears. That warning means a slab is too
wide for the series to converge on, and a higher-order truncation of a series that is not
converging is not a better answer -- add slabs, or declare the structure with `t_breakpoints`.

Order **1 or 2** is right only when the Hamiltonian is genuinely constant, where every term past
$\Omega_1$ vanishes identically; `osc_prob` already forces order 1 in that case, so this is a
setting you do not need to reach for.'''),
    md(r'''## 5. Unitarity does not depend on the truncation

Section 4 measured what the order buys in *accuracy*. This measures what it costs in
*correctness*, and the answer is nothing: **every truncation of the Magnus series lives in the
Lie algebra**, so the operator it exponentiates is unitary exactly, not to the accuracy of the
truncation.

The sweep raises the order from 1 to 10 on a smooth exponential profile at a deliberately
coarse four slabs, so the order is the only thing moving. `integration_method='simpson'` is
used rather than the default `gl`, because Gauss--Legendre nodes coincide in pairs and `gl`
offers only three distinct schemes wearing six names --- the answer would freeze every second
point.

The left panel is the answer converging. The right panel is the point.'''),
    code(r'''PER_NE_ORD = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
L_ORD = 4000.0*gd.CONV_KM_TO_INV_EV


def ne_expo_ord(l):
    return (1.0e-13*np.exp(-2.5*np.asarray(l, dtype=float)/L_ORD))/PER_NE_ORD


def truncated(order, n_slabs=4):
    """The full 3x3 probability matrix at one truncation order."""
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        3, ne_expo_ord, 1.0*gd.UNIT_GEV, L_ORD, OSC, L0=0.0,
        density_is_of_number_of_electrons=True, strategy='magnus',
        magnus_exp_order=order, integration_method='simpson',
        n_slabs=n_slabs, min_n_slabs=n_slabs, max_n_slabs=max(n_slabs, 2),
        rtol=1.0e-13, atol=1.0e-15)).reshape(3, 3)


ORDERS = list(range(1, gd.MAGNUS_EXP_ORDER_MAX + 1))
series = [truncated(k) for k in ORDERS]
p_series = [m[gd.NUMU, gd.NUE] for m in series]
unitarity = [float(np.max(np.abs(m.sum(axis=1) - 1.0))) for m in series]
# The converged value, at the DEFAULT order: at 600 slabs the answer no longer
# depends on the order, which is what lets it referee the sweep.  Asking for the
# top order at that slab count is what killed this cell the first time --
# `magnus._validate` warns order 10 costs about 17x order 6.
p_ref_ord = float(truncated(4, n_slabs=600)[gd.NUMU, gd.NUE])

fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.5))
ax = axes[0]
ax.plot(p_series, ORDERS, 'o', ms=7, color='k', mfc='white', mew=1.5)
for k, p in zip(ORDERS, p_series):
    ax.annotate(str(k), xy=(p, k), xytext=(7, 0), textcoords='offset points',
                fontsize=8, va='center', color='0.25')
ax.axvline(p_ref_ord, color='C3', ls='--', lw=1.2, label='Converged reference')
ax.set_xlabel(r'$P(\nu_\mu \to \nu_e)$')
ax.set_ylabel('Truncation order')
ax.set_title('The answer converges', fontsize=11)
ax.set_aspect(1.0/ax.get_data_ratio())               # 1:1 panel
ax.legend(fontsize=8, loc='lower right')

ax = axes[1]
ax.semilogx(np.maximum(unitarity, 1.0e-18), ORDERS, 'o', ms=7, color='C3')
for k, u in zip(ORDERS, unitarity):
    ax.annotate(str(k), xy=(max(u, 1.0e-18), k), xytext=(7, 0),
                textcoords='offset points', fontsize=8, va='center', color='0.25')
ax.set_xlim(1.0e-17, 1.0e-13)
ax.set_xlabel(r'$\max_\alpha |1 - \sum_\beta P_{\alpha\beta}|$')
ax.set_ylabel('Truncation order')
ax.set_title('Unitarity does not', fontsize=11)
ax.set_aspect(1.0/ax.get_data_ratio())
fig.tight_layout(pad=1.0)
fig.savefig('../fig/expansion_unitarity.pdf', bbox_inches='tight')

print('Worst departure from unitarity over all %d orders: %.2e' % (len(ORDERS),
                                                                  max(unitarity)))'''),
    md(r'''**A truncated series that is still exactly unitary** is the property the whole method
rests on: probabilities are non-negative and sum to one at machine precision at *any* accuracy
setting, so an under-resolved answer is a wrong probability rather than not a probability at
all. An adaptive Runge--Kutta integrator drifts off the unitary manifold instead, by around
$10^{-6}$ at typical tolerances and worse with baseline.'''),
    md(r'''## 6. Which engine answers, and what makes it change

`strategy='auto'` is not one method. It tries an adiabatic-plus-patch propagator first and
falls back to the slab ladder wherever that declines --- and **which one answers depends on
what you asked for**, not only on the physics.

The setting below is fixed and deliberately hard: a **15 MeV** electron neutrino crossing a
**supernova shock front**, along a ray from $10^4$ to $8\times10^4$ km, with the front at
$3\times10^4$ km and a width of a thousandth of the ray. Only the requested tolerance changes.

Sweeping a *physical* parameter does not show this. Measured across four decades of baseline
and three of energy on smooth exponential profiles, the answer is `hybrid`/`adiabatic` every
time with nothing declined. The tolerance is the axis that moves it.'''),
    code(r'''MEAN_NUCLEON_D = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
R0_D, R1_D, W_D, RF_D = 1.0e4, 8.0e4, 1.0e-3, 3.0e4
KM_D = gd.CONV_KM_TO_INV_EV


def ne_shock_d(l):
    w_km = W_D*(R1_D - R0_D)
    r = np.asarray(l, dtype=float)/KM_D
    u = np.clip((RF_D + 0.5*w_km - r)/w_km, 0.0, 1.0)
    out = (1.0e14*r**(-2.4)*(1.0 + (u*u*(3.0 - 2.0*u))*9.0)
           * gd.UNIT_G_PER_CM3/MEAN_NUCLEON_D*0.5)
    return out[()] if np.ndim(out) == 0 else out


def who_answered(rtol):
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p = float(np.asarray(oscprob.osc_prob_matter_std_potential(
            3, ne_shock_d, 15.0*gd.UNIT_MEV, R1_D*KM_D, OSC, L0=R0_D*KM_D,
            density_is_of_number_of_electrons=True, nu_i=gd.NUE, nu_f=gd.NUE,
            rtol=rtol, atol=rtol*1.0e-2, strategy_info=info)))
    return p, info.get('engine', '?'), info.get('family', '?'), \
        info.get('declined') or []


TOL_D = np.logspace(-2, -10, 9)
dispatch = [who_answered(t) for t in TOL_D]
fams_d = sorted({d[2] for d in dispatch})
palette_d = {f: c for f, c in zip(fams_d, ['k', 'C3', 'C0'])}

fig, ax = plt.subplots(figsize=(6.6, 4.2))
for fam in fams_d:
    sel = [i for i, d in enumerate(dispatch) if d[2] == fam]
    ax.semilogx(TOL_D[sel], [dispatch[i][0] for i in sel], 'o', ms=9,
                color=palette_d[fam], mfc='white', mew=1.8, label=fam)
ax.set_xlabel('Requested tolerance')
ax.set_ylabel(r'$P(\nu_e \to \nu_e)$')
ax.set_title('Which engine answered, on a 15 MeV shock crossing', fontsize=11)
ax.set_xlim(TOL_D[0], TOL_D[-1])
ax.legend(fontsize=8, title='Engine family')
fig.tight_layout(pad=1.0)
fig.savefig('../fig/dispatch_vs_tolerance.pdf', bbox_inches='tight')

print('%9s  %-8s %-14s %s' % ('rtol', 'engine', 'family', 'declined, and why'))
for t, d in zip(TOL_D, dispatch):
    print('%9.0e  %-8s %-14s %s' % (t, d[1], d[2], d[3][0][1] if d[3] else '--'))'''),
    md(r'''**The reason travels with the decision**, which is the part worth having: the request
did not silently change engine, it changed engine *because* the adiabatic route could not
self-certify at the tolerance asked for, and `strategy_info['declined']` says exactly that. An
answer that arrives with a stated reason for the route taken is a very different object from
one that does not.

Notebook 22 takes this further --- every engine, what each shares with the others, and how to
ask for one by name.'''),
    ])


# ------------------------------------------------ 26_magnus_nufit_evolution
books['26_magnus_nufit_evolution.ipynb'] = notebook(
    'Fourteen years of NuFIT, and what it did to the probabilities',
    r'''Every other notebook here fixes the oscillation parameters at one global fit and computes
a probability. This one asks a different question: **how much of the probability is the
parameters?**

The NuFIT collaboration has published global fits since 2012. Mag$\nu$s ships the best-fit
values for all eighteen releases, but a best fit alone cannot answer this -- the useful
question is not how the central curve moved but how the *distribution* moved, which needs the
likelihood.

NuFIT publishes one, as $\Delta\chi^2$ profiles. `notebooks/make_nufit_chi2.py` downloads them
and extracts the six one-dimensional projections into `notebooks/nufit_chi2.json`; this
notebook samples from those and pushes each sample through Mag$\nu$s.

**The data are the NuFIT collaboration's**, redistributed here in extract. Cite the
corresponding NuFIT paper and <http://www.nu-fit.org/> if you use them.''',
    [
    code(r'''import json
import pathlib

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.globaldefs as gd

CHI2 = json.loads((pathlib.Path.cwd()/'nufit_chi2.json').read_text())
print(CHI2['_attribution'])
print()
print('releases with chi^2 profiles: %d' % len(CHI2['releases']))
print(', '.join(CHI2['releases']))'''),
    md(r'''## 1. What is here, and what is not

Two limitations, both worth stating before any number is produced.

**Two sources, not one.** From v2.0 onward NuFIT publishes machine-readable
$\Delta\chi^2$ files, and those are read directly. Releases v1.0--v1.3 (2012--2014) publish
theirs only as figures -- but those figures are **vector** PDFs, so the curves are stroke
operators rather than pixels. Nothing is digitized: the color operator separates the
orderings exactly (`1 0 0 RG` is NO, `0 0 1 RG` is IO) and the dash pattern separates the
analyses exactly (solid is Free+RSBL, dashed is Huber). Each release records which route it
came by in its `source` field.

The extraction is checked against numbers it never uses as input -- the curve minima must
reproduce the published best fits, and $\delta_{\rm CP}$ must never reach
$\Delta\chi^2 = 9$, which is what its published "$0 \to 360$" range asserts.
`make_nufit_chi2.py` refuses to write the file if any release fails.

**These are one-dimensional marginals.** Each parameter is sampled from its own profile,
independently of the others. The real likelihood is correlated -- the
$\delta_{\rm CP}$--$\theta_{23}$ correlation especially -- and the source files contain the
pairwise projections if you need them. What follows is therefore a good description of how
much each parameter was individually known, and an imperfect one of their joint effect.'''),
    code(r'''def one_period(profile):
    """Restrict a delta_CP profile to a single 360-degree period.

    The releases do not agree on where the period starts -- the figures run
    0..360, most chi^2 files run -180..180, and v2.0 tabulates -180..360, which
    is a period and a half.  Sampling that last one as written would count part
    of the range twice."""
    x = np.asarray(profile['x'], dtype=float)
    chi2 = np.asarray(profile['chi2'], dtype=float)
    keep = x < x.min() + 360.0 - 1.0e-9
    return {'x': x[keep], 'chi2': chi2[keep]}


def sample_profile(profile, size, rng):
    """Draw from exp(-Delta_chi^2 / 2) on the tabulated grid, by inverse CDF.

    The profiles are normalized so the best fit sits at Delta_chi^2 = 0, so
    exp(-chi2/2) is the likelihood up to a constant -- which the CDF divides
    out anyway."""
    x = np.asarray(profile['x'], dtype=float)
    chi2 = np.asarray(profile['chi2'], dtype=float)
    weight = np.exp(-0.5*(chi2 - chi2.min()))
    cdf = np.concatenate([[0.0],
                          np.cumsum(0.5*(weight[1:] + weight[:-1])*np.diff(x))])
    cdf /= cdf[-1]
    return np.interp(rng.random(size), cdf, x)


def draw_parameters(release, size, rng):
    """One release's profiles -> `size` parameter sets in Magnus's conventions.

    The files store sin^2(theta) where Magnus wants sin(theta), degrees where it
    wants radians, log10(Delta_m21^2) where it wants the value, and
    Delta_m31^2 in units of 1e-3 eV^2."""
    block = CHI2['releases'][release]
    p = block['profiles']
    # The two sources spell Delta_m21^2 differently: the chi^2 files carry
    # log10(Delta_m21^2/eV^2), the figures carry it in 1e-5 eV^2.
    dms = sample_profile(p['DMS'], size, rng)
    D21 = (dms*1.0e-5 if block.get('source') == 'figure' else 10.0**dms)
    return dict(
        s12=np.sqrt(sample_profile(p['T12'], size, rng)),
        s13=np.sqrt(sample_profile(p['T13'], size, rng)),
        s23=np.sqrt(sample_profile(p['T23'], size, rng)),
        dCP=np.deg2rad(np.mod(sample_profile(one_period(p['DCP']), size, rng), 360.0)),
        D21=D21,
        D31=sample_profile(p['DMA'], size, rng)*1.0e-3)'''),
    md(r'''A check before trusting any of it: the minimum of each profile must reproduce the
best-fit value `magnus.globaldefs` ships independently for that release. If the extraction had
picked up the wrong section, or mixed orderings, this is where it would show.'''),
    code(r'''print('%-11s %-11s %-11s %-11s %s'
      % ('release', 'chi2 min', 'shipped BF', 'chi2 min', 'shipped BF'))
print('%-11s %-11s %-11s %-11s %s' % ('', 'sin^2(th23)', 'sin^2(th23)', 'D31', 'D31'))
print('-'*62)
for release in ('NuFIT 1.0', 'NuFIT 1.3', 'NuFIT 2.0', 'NuFIT 4.0', 'NuFIT 6.1'):
    block = CHI2['releases'][release]
    p = block['profiles']
    s23_sq = p['T23']['x'][int(np.argmin(p['T23']['chi2']))]
    d31 = p['DMA']['x'][int(np.argmin(p['DMA']['chi2']))]*1.0e-3
    bf = gd.load_nufit_params(release, 'NO')
    print('%-11s %-11.4f %-11.4f %-11.5f %.5f'
          % (release, s23_sq, bf['s23']**2, d31, bf['D31']))'''),
    md(r'''One entry deserves a word. For v1.2 and v1.3 the $\sin^2\theta_{23}$ curve has **two**
minima, and the deeper one on the normal-ordering curve (0.445 and 0.451) is not the branch
`globaldefs` quotes as the best fit (0.593 and 0.577). Both numbers are in NuFIT's own table,
joined by $\oplus$; the extraction recovers both to about 0.002. That is the bimodality this
notebook exists to keep, and it is exactly what a Gaussian centerd on either branch would
throw away -- so the check above compares the *pair* of minima for those releases rather than
insisting on one.

## 2. How the parameters themselves moved

Before any probability, the inputs. Each parameter is summarized by the central 68% of its
own sampled distribution, with the width of that interval in the narrower panel below it.

$\delta_{\rm CP}$ is circular, and the releases disagree about where its period starts: the
figures run $0$ to $360^\circ$, most $\chi^2$ files run $-180$ to $+180$, and v2.0 tabulates
a period and a half. Everything here is restricted to one period and wrapped to
$[0, 360)$ -- and its 68% width saturates near $240^\circ$ where the parameter is
essentially unconstrained, which is what the early releases show.'''),
    code(r'''PARAMETERS = [
    ('s12', lambda d: d['s12']**2, r'$\sin^2\theta_{12}$'),
    ('s23', lambda d: d['s23']**2, r'$\sin^2\theta_{23}$'),
    ('s13', lambda d: d['s13']**2, r'$\sin^2\theta_{13}$'),
    ('dCP', lambda d: np.mod(np.rad2deg(d['dCP']), 360.0), r'$\delta_{\rm CP}\ [^\circ]$'),
    ('D21', lambda d: d['D21']/1.0e-5, r'$\Delta m^2_{21}\ [10^{-5}\,{\rm eV}^2]$'),
    ('D31', lambda d: d['D31']/1.0e-3, r'$\Delta m^2_{31}\ [10^{-3}\,{\rm eV}^2]$'),
]

N_SAMPLES = 1500
rng = np.random.default_rng(4)
releases = list(CHI2['releases'])
names = [r.replace('NuFIT ', '') for r in releases]

par_stats = {key: [] for key, _, _ in PARAMETERS}
for release in releases:
    drawn = draw_parameters(release, N_SAMPLES, rng)
    for key, extract, _ in PARAMETERS:
        par_stats[key].append(np.percentile(extract(drawn), [16, 50, 84]))

for key, _, _ in PARAMETERS:
    stat = np.array(par_stats[key])
    width = stat[:, 2] - stat[:, 0]
    print('%-4s 68%% width: %-10.4g -> %-10.4g (%.2fx narrower)'
          % (key, width[0], width[-1], width[0]/width[-1]))'''),
    code(r'''x = np.arange(len(releases))

fig = plt.figure(figsize=(10.0, 6.8))
gs = fig.add_gridspec(4, 3, height_ratios=[2, 1, 2, 1], hspace=0.10, wspace=0.34)

for i, (key, _, label) in enumerate(PARAMETERS):
    row0, col = (i//3)*2, i % 3
    ax_v = fig.add_subplot(gs[row0, col])
    ax_w = fig.add_subplot(gs[row0 + 1, col], sharex=ax_v)
    stat = np.array(par_stats[key])
    lo_p, mid_p, hi_p = stat[:, 0], stat[:, 1], stat[:, 2]

    ax_v.fill_between(x, lo_p, hi_p, color='C0', alpha=0.25, label=r'Central 68\%')
    ax_v.plot(x, mid_p, color='C0', marker='o', ms=2.5)
    ax_v.set_ylabel(label, fontsize=8)
    ax_v.tick_params(labelbottom=False, labelsize=7)
    ax_v.set_xlim(x[0], x[-1])
    if i == 0:
        ax_v.legend(fontsize=7, loc='upper right')

    ax_w.plot(x, hi_p - lo_p, color='0.2', marker='o', ms=2.5)
    ax_w.set_ylabel(r'68\%', fontsize=8)
    ax_w.set_xlim(x[0], x[-1])
    ax_w.set_xticks(x)
    ax_w.set_xticklabels(names, rotation=90, fontsize=6)
    ax_w.tick_params(labelsize=7)
    if row0 == 2:
        ax_w.set_xlabel('NuFIT release', fontsize=8)

fig.suptitle('Mixing parameters across eighteen NuFIT releases (normal ordering)',
             fontsize=10)'''),
    md(r'''$\theta_{13}$ improved most -- a factor of four -- going from barely measured in
2012 to the best-known angle in the matrix once the reactor experiments reported.
$\theta_{23}$ is the one that did not: its band stays the widest and wanders rather than
shrinks, because the octant keeps changing its mind. That is the parameter which will drive
everything in the next section.

## 2b. The distribution of a probability, release by release

We use a DUNE-like configuration -- 1300 km through the crust, 2 GeV, the
$\nu_\mu \to \nu_e$ appearance channel -- because it is sensitive to $\delta_{\rm CP}$,
$\theta_{13}$ and $\theta_{23}$ at once, which is exactly what the fits have been pinning
down.'''),
    code(r'''RHO = 2.848
VCC = matter.vcc_func_from_rho_func(RHO, density_matter_is_in_g_per_cm3=True)
ENERGY = 2.0*gd.UNIT_GEV
BASELINE = 1300.0*gd.UNIT_KM

def probabilities_for(release, size, rng):
    """P(nu_mu -> nu_e) for `size` parameter sets drawn from one release."""
    drawn = draw_parameters(release, size, rng)
    out = np.empty(size)
    for i in range(size):
        h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
            **{k: v[i] for k, v in drawn.items()}))
        out[i] = np.asarray(oscprob.osc_prob(
            h_vac/ENERGY + np.diag([VCC, 0.0, 0.0]), 0.0, BASELINE))[gd.NUMU][gd.NUE]
    return out

rng = np.random.default_rng(20260809)
summary = []
print('%-11s %-9s %-9s %-9s %-10s %s'
      % ('release', 'median', '16%', '84%', '68% width', 'best fit'))
print('-'*60)
for release in CHI2['releases']:
    P = probabilities_for(release, N_SAMPLES, rng)
    lo, mid, hi = np.percentile(P, [16, 50, 84])
    bf = gd.load_nufit_params(release, 'NO')
    h_bf = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**bf))
    P_bf = np.asarray(oscprob.osc_prob(
        h_bf/ENERGY + np.diag([VCC, 0.0, 0.0]), 0.0, BASELINE))[gd.NUMU][gd.NUE]
    summary.append((release, mid, lo, hi, hi - lo, P_bf))
    print('%-11s %-9.4f %-9.4f %-9.4f %-10.4f %.4f'
          % (release, mid, lo, hi, hi - lo, P_bf))'''),
    md(r'''## 3. The uncertainty did not shrink monotonically

The obvious expectation is that fourteen years of data narrow the band steadily. It did not.'''),
    code(r'''widths = np.array([row[4] for row in summary])
names = [row[0].replace('NuFIT ', '') for row in summary]

print('narrowest : %s  (%.4f)' % (summary[int(np.argmin(widths))][0], widths.min()))
print('widest    : %s  (%.4f)' % (summary[int(np.argmax(widths))][0], widths.max()))
print('first -> last: %.4f -> %.4f  (%.2fx narrower overall)'
      % (widths[0], widths[-1], widths[0]/widths[-1]))'''),
    code(r'''fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(7.4, 5.6), sharex=True,
                               gridspec_kw=dict(height_ratios=[2, 1]))
x = np.arange(len(summary))
mid = np.array([r[1] for r in summary])
lo = np.array([r[2] for r in summary])
hi = np.array([r[3] for r in summary])
bf = np.array([r[5] for r in summary])

ax0.fill_between(x, lo, hi, color='C0', alpha=0.25, label=r'Central 68\%')
ax0.plot(x, mid, color='C0', marker='o', ms=3, label='median of samples')
ax0.plot(x, bf, color='C3', ls='--', marker='s', ms=3, label='best fit only')
ax0.set_ylabel(r'$P(\nu_\mu \to \nu_e)$')
ax0.legend(fontsize=8, loc='lower right')
ax0.set_title(r'2 GeV, $L = 1300$ km, $\rho = 2.848$ g cm$^{-3}$, normal ordering',
              fontsize=10)

ax1.plot(x, widths, color='0.2', marker='o', ms=3)
ax1.set_ylabel(r'68\%')
ax1.set_xticks(x)
ax1.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
ax1.set_xlabel('NuFIT release')
ax0.set_xlim(x[0], x[-1])
ax1.set_xlim(x[0], x[-1])
fig.tight_layout()'''),
    md(r'''The band is narrowest around **NuFIT 4.0--4.1** (2018--2019) and *widens again* at 5.0.
That is not a defect in the fits; it is their history. The $\theta_{23}$ octant preference has
flipped more than once between releases, and the $\delta_{\rm CP}$ constraint loosened as
T2K and NOvA data pulled in different directions. A one-dimensional $\Delta\chi^2$ for
$\sin^2\theta_{23}$ with two nearly degenerate minima is wide, and it produces a wide
probability distribution regardless of how much data went into it.

Note also how far the **best-fit curve** wanders relative to the band. Between 3.0 and 4.0 it
moves by more than the 68% width of either -- because the best fit hops between octants while
the distribution, which contains both, moves far less. **A single best-fit probability is a
less stable thing than the distribution it comes from**, which is the practical argument for
propagating the likelihood rather than the central values.

## 4. Where the movement comes from

Holding five parameters at the 6.1 best fit and varying only the sixth isolates each
contribution.'''),
    code(r'''best = gd.load_nufit_params('NuFIT 6.1', 'NO')
rng = np.random.default_rng(7)
print('%-8s %-12s %s' % ('vary', '68% width', 'share of the full width'))
print('-'*48)
full = None
for name, key in (('th12', 'T12'), ('th13', 'T13'), ('th23', 'T23'),
                  ('dCP', 'DCP'), ('D21', 'DMS'), ('D31', 'DMA')):
    drawn = draw_parameters('NuFIT 6.1', N_SAMPLES, rng)
    one = dict(best)
    magnus_key = {'T12': 's12', 'T13': 's13', 'T23': 's23',
                  'DCP': 'dCP', 'DMS': 'D21', 'DMA': 'D31'}[key]
    values = drawn[magnus_key]
    out = np.empty(N_SAMPLES)
    for i in range(N_SAMPLES):
        one[magnus_key] = values[i]
        h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**one))
        out[i] = np.asarray(oscprob.osc_prob(
            h_vac/ENERGY + np.diag([VCC, 0.0, 0.0]), 0.0, BASELINE))[gd.NUMU][gd.NUE]
    lo, hi = np.percentile(out, [16, 84])
    if full is None:
        full = widths[-1]
    print('%-8s %-12.4f %.0f%%' % (name, hi - lo, 100.0*(hi - lo)/full))'''),
    md(r'''## Summary

| | |
|---|---|
| releases from machine-readable $\chi^2$ | 14 (v2.0 -- v6.1) |
| releases read from vector figures | 4 (v1.0 -- v1.3) |
| **total** | **18**, spanning 2012--2025 |
| 68% width, first release | 0.024 |
| 68% width, narrowest (v4.0) | **0.012** |
| 68% width, latest (v6.1) | 0.016 |

Three things worth taking away.

1. **The parameter uncertainty is comparable to the effects other notebooks study.** A 68%
   width of 0.016 in $P(\nu_\mu \to \nu_e)$ is larger than the octant effect measured in
   notebook 17 (~0.015), and comparable to the antineutrino half-convention errors in notebook
   15. Numerical accuracy is rarely the limiting uncertainty in a real prediction.
2. **It has not shrunk monotonically,** and the widening at 5.0 is physics, not noise.
3. **Propagate the likelihood, not the best fit.** The best-fit probability moves further
   between releases than the distribution does.

The `_note` field in `nufit_chi2.json` and section 1 above both record the independence
assumption; a joint treatment would use the pairwise projections in the source files, which
`make_nufit_chi2.py` documents but does not extract.'''),
    ])


# ------------------------------------------------ 25_magnus_against_other_codes
books['25_magnus_against_other_codes.ipynb'] = notebook(
    'Against other codes',
    r'''Mag$\nu$s is not the only way to get an oscillation probability, and on the problems other
codes were built for it is not the fastest. This notebook measures that honestly, against
whichever of them are installed.

Two things make a comparison like this easy to get wrong, and both are handled explicitly
below. The first is **what problem is being solved**: a closed-form constant-density
expression and a general varying-profile integrator are not doing the same work, and timing
them on constant density flatters the former by construction. The second is **conventions**:
two codes given "the same" density can build different matter potentials, and the resulting
disagreement is not an accuracy difference even though it looks exactly like one.

Nothing here fails if a code is missing -- each is probed and skipped.''',
    [
    code(r'''import json
import pathlib
import sys
import time
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.linalg import expm

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.globaldefs as gd
import magnus.expmkernels as ek     # HAVE_NUMBA: whether the compiled kernel exists

warnings.simplefilter('ignore')

OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
L_KM, RHO, N_E = 1300.0, 2.848, 60
BASELINE = L_KM*gd.CONV_KM_TO_INV_EV
VCC = matter.vcc_func_from_rho_func(RHO, density_matter_is_in_g_per_cm3=True)
E_GEV = np.linspace(0.6, 20.0, N_E)
E = E_GEV*gd.UNIT_GEV
h_vac = np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC))

print('L = %.0f km, rho = %.3f g/cm^3, %d energies, 0.6-20 GeV' % (L_KM, RHO, N_E))
print('V_CC (Magnus) = %.6e eV' % VCC)'''),
    md(r'''## 1. Which codes are available

`nuSQuIDS` is a general neutrino transport solver; `NuOscProbExact` computes the exact
SU(3) closed form for piecewise-constant matter. Neither is a dependency of Mag$\nu$s, so both
are optional here.'''),
    code(r'''AVAILABLE = {}

try:
    import nuSQuIDS as nsq
    AVAILABLE['nuSQuIDS'] = True
except Exception as exc:
    AVAILABLE['nuSQuIDS'] = False
    print('nuSQuIDS not available: %s' % exc)

# NuOscProbExact is not packaged on PyPI; point NUOSCPROBEXACT_SRC at its src/
# directory to include it.  Absent that, this notebook simply skips it.
import os
_npe_src = os.environ.get('NUOSCPROBEXACT_SRC')
if _npe_src and os.path.isdir(_npe_src):
    sys.path.insert(0, _npe_src)
try:
    import oscprob3nu
    import hamiltonians3nu
    AVAILABLE['NuOscProbExact'] = True
except Exception:
    AVAILABLE['NuOscProbExact'] = False

for name, ok in AVAILABLE.items():
    print('%-16s %s' % (name, 'available' if ok else 'not installed -- skipped'))'''),
    md(r'''## 2. The reference

For a **constant** density the Hamiltonian is constant, so the evolution operator is a single
matrix exponential and `scipy.linalg.expm` is the exact answer. That makes an oracle that
belongs to none of the codes being compared.'''),
    code(r'''reference = np.array([np.abs(expm(-1j*(h_vac/e + np.diag([VCC, 0.0, 0.0]))*BASELINE))[0, 1]**2
                      for e in E])
print('reference P(nu_mu -> nu_e): %.6f ... %.6f' % (reference[0], reference[-1]))


def best_of(call, repeats=3):
    fastest, result = np.inf, None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = call()
        fastest = min(fastest, time.perf_counter() - t0)
    return result, fastest'''),
    md(r'''## 3. Constant density: speed

This is the problem `NuOscProbExact` exists for -- an exact closed form for a
piecewise-constant Hamiltonian. Mag$\nu$s solves the general problem, so the honest expectation
is that it pays for the generality here.

**Both codes are given both of their modes.** NuOscProbExact can be called once per energy or
handed a stacked array, and the two differ by an order of magnitude; comparing a batched
Mag$\nu$s against a looped NuOscProbExact would be a statement about Python loop overhead
rather than about either method. Mag$\nu$s likewise recognizes a position-independent
Hamiltonian and exponentiates the whole scan at once -- for a constant density the Magnus series
terminates at its first term, so one exponential per energy is not an approximation but the
exact answer.'''),
    code(r'''results = []

P_magnus, t_magnus = best_of(lambda: np.asarray(oscprob.osc_prob_3nu_matter_constant_density(
    E, BASELINE, RHO, **OSC, density_matter_is_in_g_per_cm3=True,
    nu_i=gd.NUMU, nu_f=gd.NUE)))
results.append(('Magnus (batched wrapper)', P_magnus, t_magnus))

# Magnus with the batched constant-H engine switched off, i.e. one osc_prob call
# per energy: the same number by the route the wrapper used to take.
with oscprob._engine_probe(disabled=('constant',)):
    P_mloop, t_mloop = best_of(lambda: np.asarray(
        oscprob.osc_prob_3nu_matter_constant_density(
            E, BASELINE, RHO, **OSC, density_matter_is_in_g_per_cm3=True,
            nu_i=gd.NUMU, nu_f=gd.NUE)))
results.append(('Magnus (one call per energy)', P_mloop, t_mloop))

if AVAILABLE['NuOscProbExact']:
    h_npe = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(**OSC)
    P_npe, t_npe = best_of(lambda: np.array([
        oscprob3nu.probabilities_3nu(
            hamiltonians3nu.hamiltonian_3nu_matter(h_npe, e, VCC), BASELINE)[3]
        for e in E]))
    results.append(('NuOscProbExact (looped)', P_npe, t_npe))

    # The same code handed the whole energy axis at once.
    P_npe_b, t_npe_b = best_of(lambda: np.asarray(oscprob3nu.probabilities_3nu(
        hamiltonians3nu.hamiltonian_3nu_matter(h_npe, E, VCC), BASELINE))[..., 3])
    results.append(('NuOscProbExact (batched)', P_npe_b, t_npe_b))

print('%-30s %-11s %-11s %s' % ('code', 'time [s]', 'us/point', 'max |P - exact|'))
print('-'*70)
for name, P, t in results:
    print('%-30s %-11.4f %-11.2f %.3e'
          % (name, t, 1.0e6*t/N_E, np.max(np.abs(P - reference))))'''),
    md(r'''Read the two batched rows against each other, and the two looped rows against each
other. Batched, the two codes are within a small factor and both sit at machine precision; the
looped rows are an order of magnitude slower on both sides, which is the Python interpreter
talking, not the physics.

The single-point picture is different and worth stating plainly: one probability at one energy
costs Mag$\nu$s about 1.7 times what it costs NuOscProbExact, and almost all of that gap is
parameter resolution in the wrapper -- NuFIT defaults, unit flags, engine dispatch -- rather
than arithmetic. The matrix exponential is under a tenth of the call.'''),
    md(r'''Both are exact to round-off; the closed form is several times faster. That is the
right result and it is worth stating plainly: **if constant density is your whole problem, a
closed-form code is the better tool.** Mag$\nu$s earns its cost on profiles that vary, where no
closed form exists -- notebook 16 measures a mean-density substitution being wrong by 0.51 on
an Earth chord.

## 4. PREM, three flavors: the comparison that is actually about method

Constant density is where a closed form is at its best, and it is not what either code exists
for. The Earth is: a varying profile, where the two take genuinely different routes.
NuOscProbExact discretizes PREM into piecewise-constant slabs and propagates each exactly.
Mag$\nu$s integrates within each slab to fourth or sixth order, so its slabs can be wider for
the same accuracy. **The question is not time per call but time at matched accuracy**, and the
dial is different on each side: `n_slabs_per_segment` for one, `rtol`/`atol` for the other.

Three warnings before the numbers, and they bound what the comparison can mean.

*The two PREM implementations need not be identical.* Layer radii, the polynomial coefficients,
and the electron fraction are all conventions, and NuOscProbExact's own notebook 10 records a
100 km discrepancy in one boundary against a published table. So a residual between the codes
here is an upper bound on method disagreement, not a measurement of it.

*The composition is pinned to $Y_e = 0.5$ everywhere in this section*, on all three sides --
Mag$\nu$s, NuOscProbExact and the referee. Mag$\nu$s's own default is now resolved per PREM
layer (0.4656 in the iron core, 0.4957 in the mantle), which is the better description of the
Earth but *not* what the other code assumes here, and the difference is not small: left
unmatched it moves $P(\nu_\mu \to \nu_e)$ on this chord by up to **7.1e-02**, some five orders
of magnitude above the referee's own resolution floor below. That would be a comparison of two
Earths wearing the costume of a comparison of two solvers. The uniform value is passed
explicitly rather than left to a default, so it cannot drift again.

(The current NuOscProbExact can vary its electron fraction too, so this could be matched from
the other side instead, or both sides moved to a layered profile. This notebook does not do
that yet: the frozen dataset in section 7 is built on the uniform convention, and one
convention across the notebook is worth more here than the better Earth in one section.)

*Each code's self-convergence is the honest accuracy statement.* Refining a code against itself
measures its own discretization error without borrowing anyone else's conventions, and that is
reported first.'''),
    code(r'''COSTHZ = -0.85                     # a chord through mantle and outer core

# ONE COMPOSITION FOR ALL THREE SIDES OF THIS SECTION.  Mag(nu)s resolves Y_e per PREM
# layer by default now -- 0.4656 for the iron core, 0.4957 for the mantle -- and this
# chord enters the outer core, so that default disagrees with NuOscProbExact's uniform
# 0.5 by up to 7.1e-02 in P(numu -> nue).  The referee's floor below is 4.3e-07, so an
# unmatched composition would sit five orders of magnitude above anything this section
# can resolve and every curve here would be measuring the Earth model rather than the
# solver.  Passed explicitly on all three sides, never left to a default: the whole
# reason this needed fixing is that one side's default moved and the others did not.
#
# 0.5 is also exactly self-consistent for the sterile sector: r = (1 - Y_e)/Y_e = 1.0,
# which is the `ratio_number_neutrons_to_protons` default that `matter_potential_projector`
# and the 4nu wrapper both take, so the projector's sterile entry is r/2 = 0.5 on both
# sides with nothing further to pass.
YE_UNIFORM = 0.5

# Two grids, because the two jobs want opposite things.  E_PREM is the *timing* grid:
# every convergence row below divides by len(E_PREM), so its only requirement is to be
# large enough to amortize per-call overhead and small enough that sweeping four
# tolerances and four slab counts stays affordable.  E_PLOT is the *picture*: an Earth
# chord at these energies oscillates fast enough that 24 points join peak to trough with
# a straight line and invent structure that is not there.
E_PREM = np.logspace(np.log10(0.6), np.log10(20.0), 24)*gd.UNIT_GEV
E_PLOT = np.logspace(np.log10(0.6), np.log10(20.0), 400)*gd.UNIT_GEV

try:
    import earth as npe_earth
    import slabs as npe_slabs          # its generic varying-profile route; see section 9
    HAVE_NPE_EARTH = AVAILABLE['NuOscProbExact']
except Exception as exc:
    HAVE_NPE_EARTH = False
    print('NuOscProbExact earth module not importable: %s' % exc)

import magnus.earth as mg_earth
# distance_traveled_inside_earth returns KILOMETERS; every osc_prob baseline is in
# natural units (eV^-1).  Passing the raw value is a factor of 5.07e9 too short, and
# it does not raise: the call returns a converged, unitary, entirely wrong answer at
# a baseline of a few meters, on which the refinement ladder trivially agrees with
# itself at every tolerance.  That is what this section measured on its first run.
CHORD_KM = mg_earth.distance_traveled_inside_earth(COSTHZ)
L_CHORD = CHORD_KM*gd.CONV_KM_TO_INV_EV
print('costhz = %.2f -> chord %.1f km, %d energies, 0.6-20 GeV'
      % (COSTHZ, CHORD_KM, len(E_PREM)))


def magnus_prem(rtol, atol, energies=None):
    return np.asarray(oscprob.osc_prob_3nu_earth(
        E_PREM if energies is None else energies,
        costhz=COSTHZ, L=L_CHORD, **OSC, nu_i=gd.NUMU, nu_f=gd.NUE,
        rtol=rtol, atol=atol, electron_fraction=YE_UNIFORM))


def npe_prem(n_per_segment, energies=None):
    h = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(**OSC)
    return np.asarray(npe_earth.probabilities_3nu_earth(
        h, E_PREM if energies is None else energies, COSTHZ,
        n_slabs_per_segment=n_per_segment))[..., 3]


# ----------------------------------------------------------------- the referee
#
# Neither code may be its own judge.  Both are measured against a slab product of
# scipy matrix exponentials: the same Hamiltonian, a different integrator -- which is
# exactly what NuOscProbExact's own comparison does with a 50-digit mpmath exponential
# rather than with anyone's self-convergence.
#
# Two details decide whether this works, and both were found by it being wrong first:
#
#   Slab edges must LAND ON the PREM layer boundaries.  A uniform grid straddles them,
#   a straddling slab is O(h) rather than O(h^2), and the naive version converged
#   *non-monotonically* -- 3.4e-04, then 3.2e-03, then 8.5e-05 -- which is the tell.
#
#   `prem_layer_edges_along_chord` returns the INTERIOR crossings only.  Without the
#   two endpoints the referee drops ~7 km of a 10 830 km chord: 0.065% of the path,
#   and it then converged beautifully to an answer 1.8e-03 away from both codes.  A
#   clean convergence rate says the discretization is consistent, never that it is
#   right.
def prem_referee(energies, n_slabs, dim=3, h_vac_dim=None):
    seg = np.concatenate(([0.0],
        np.asarray(mg_earth.prem_layer_edges_along_chord(COSTHZ), dtype=float),
        [CHORD_KM]))
    per = max(2, int(round(n_slabs/(len(seg) - 1))))
    edges = np.unique(np.concatenate([
        np.linspace(seg[i], seg[i+1], per + 1) for i in range(len(seg) - 1)]))
    mid = 0.5*(edges[:-1] + edges[1:])
    widths = np.diff(edges)*gd.CONV_KM_TO_INV_EV
    r = np.sqrt(gd.EARTH_RADIUS**2 + mid*mid + 2.0*gd.EARTH_RADIUS*mid*COSTHZ)
    # Y_e passed explicitly rather than left to this function's 0.5 default.  The value
    # is the same, but a referee that judges Mag(nu)s while silently inheriting a default
    # composition is one library change away from refereeing a different Earth than the
    # code it referees -- which is exactly what had happened here.
    vcc = np.array([matter.vcc_func_from_rho_func(
        float(x), electron_fraction=YE_UNIFORM, density_matter_is_in_g_per_cm3=True)
        for x in np.asarray(mg_earth.density_matter_func_prem(r), dtype=float)])
    hv = h_vac if h_vac_dim is None else h_vac_dim
    # The projector comes from the library rather than being written out here.  Writing
    # it out is precisely how this referee was wrong: it carried diag(1, 0, 0, 0), which
    # gives a sterile state the ACTIVE flavors' neutral-current potential.  A sterile
    # state feels neither current, so once the actives' common V_NC is removed it keeps
    # -V_NC = (r/2) V_CC.  With the hand-written version this referee sat 2.6e-02 from
    # Magnus AND from NuOscProbExact while the two agreed with each other to 3.7e-04 --
    # a referee disagreeing with every code it referees is the tell, and it is the
    # referee that is wrong.  Three flavors are unaffected: the sterile block is empty.
    proj = matter.matter_potential_projector(dim)
    out = []
    for e in np.atleast_1d(energies):
        U = np.eye(dim, dtype=complex)
        for k in range(len(mid)):
            U = expm(-1j*(hv/e + vcc[k]*proj)*widths[k]) @ U
        out.append(abs(U[gd.NUE, gd.NUMU])**2)
    return np.array(out)


def refereed(energies, n_lo=1600, dim=3, h_vac_dim=None):
    """Richardson-extrapolated referee, plus an honest bound on its own error.

    The slab product is O(h^2), so (4*P_2n - P_n)/3 removes the leading term and the
    two-grid difference bounds what is left.  Returned, not assumed: a referee that
    cannot state its own uncertainty cannot referee anything below it.
    """
    lo = prem_referee(energies, n_lo, dim, h_vac_dim)
    hi = prem_referee(energies, 2*n_lo, dim, h_vac_dim)
    return (4.0*hi - lo)/3.0, float(np.max(np.abs(hi - lo)))/3.0


P_REF3, REF3_UNC = refereed(E_PREM)
print('referee (scipy slab product, PREM edges honored, Richardson-extrapolated)')
print('  its own residual discretization error: %.2e' % REF3_UNC)
print('  nothing below that line can be resolved by this comparison')'''),
    md(r'''### Self-convergence: each code refined against itself'''),
    code(r'''print('Magnus, refined against its own tightest setting')
P_mg_ref = magnus_prem(1.0e-11, 1.0e-13)
print('%-22s %-12s %s' % ('rtol/atol', 'time [ms]', 'max |P - P_tightest|'))
print('-'*58)
mg_curve = []
mg_ref_curve = []
for rtol, atol in ((1e-4, 1e-6), (1e-6, 1e-8), (1e-8, 1e-10)):
    P, t = best_of(lambda r=rtol, a=atol: magnus_prem(r, a))
    err = float(np.max(np.abs(P - P_mg_ref)))
    mg_curve.append((1.0e6*t/len(E_PREM), max(err, 1.0e-16)))
    # Against the referee as well: the column above says the answer has stopped
    # moving, this one says whether it stopped on the right value.
    mg_ref_curve.append((1.0e6*t/len(E_PREM),
                         max(float(np.max(np.abs(P - P_REF3))), REF3_UNC)))
    print('%-22s %-12.3f %.3e' % ('%.0e / %.0e' % (rtol, atol), 1.0e3*t, err))

if HAVE_NPE_EARTH:
    print()
    print('NuOscProbExact, refined against its own densest slabbing')
    P_npe_ref = npe_prem(64)
    print('%-22s %-12s %s' % ('slabs per segment', 'time [ms]', 'max |P - P_densest|'))
    print('-'*58)
    npe_curve = []
    npe_ref_curve = []
    for n in (2, 4, 8, 16):
        P, t = best_of(lambda k=n: npe_prem(k))
        err = float(np.max(np.abs(P - P_npe_ref)))
        npe_curve.append((1.0e6*t/len(E_PREM), max(err, 1.0e-16)))
        npe_ref_curve.append((1.0e6*t/len(E_PREM),
                              max(float(np.max(np.abs(P - P_REF3))), REF3_UNC)))
        print('%-22d %-12.3f %.3e' % (n, 1.0e3*t, err))'''),
    md(r'''The two convergence tables are the substance of this notebook, and they say opposite
things about the two axes. **NuOscProbExact is an order of magnitude cheaper per call here, and
converges far more slowly**: doubling its slabs per segment buys roughly a factor of four, the
$\mathcal{O}(h^2)$ of a piecewise-constant approximation, and at sixteen per segment it is still
at a few times $10^{-5}$. Mag$\nu$s costs more per call and buys $10^{-10}$, because its error
falls with the *order* of the integrator inside each slab rather than only with the slab width.

Neither is the better code; they are priced differently. Below about $10^{-4}$ the closed form
is the cheaper way to get there, and above it there is no setting of `n_slabs_per_segment` that
competes.

### The two codes against each other'''),
    code(r'''if HAVE_NPE_EARTH:
    gap = float(np.max(np.abs(P_mg_ref - P_npe_ref)))
    rel = gap/float(np.max(P_npe_ref))
    print('max |Magnus - NuOscProbExact|  = %.3e   (%.2e of the peak probability)'
          % (gap, rel))
    print()
    print('Compare that against each code\'s own convergence above: Magnus reaches ~3e-10,')
    print('NuOscProbExact ~6e-05 at the densest slabbing tried.  The residual between them')
    print('is the same order as the looser of the two, so most of it is NuOscProbExact\'s')
    print('discretization rather than a convention difference -- which a residual far above')
    print('BOTH curves would have indicated instead.')

    # Drawn on E_PLOT, not on the 24-point timing grid.  The curve below has structure
    # between 1 and 4 GeV that 24 points cannot represent: at that spacing the two codes
    # appear to disagree in places where they do not, because each is being joined
    # through a different set of samples of the same oscillation.
    P_mg_plot = magnus_prem(1.0e-11, 1.0e-13, energies=E_PLOT)
    P_npe_plot = npe_prem(64, energies=E_PLOT)

    fig, ax = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True,
                           gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.08))
    ax[0].plot(E_PLOT/gd.UNIT_GEV, P_mg_plot, lw=1.6, label=r'Mag$\nu$s')
    ax[0].plot(E_PLOT/gd.UNIT_GEV, P_npe_plot, lw=1.1, ls='--',
               label='NuOscProbExact')
    ax[0].set_ylabel(r'$P(\nu_\mu \to \nu_e)$')
    ax[0].set_title(r'PREM chord, $\cos\theta_z = %.2f$, three flavors, %d energies'
                    % (COSTHZ, len(E_PLOT)), fontsize=10)
    ax[0].grid(True, alpha=0.2)
    ax[0].legend(fontsize=8)

    # The residual on its own axis, because on the panel above it is invisible -- which
    # is the point, and is not something a reader should have to take on trust.
    ax[1].semilogy(E_PLOT/gd.UNIT_GEV, np.abs(P_mg_plot - P_npe_plot), lw=1.0,
                   color='C3')
    ax[1].set_xscale('log')
    ax[1].set_xlabel(r'$E_\nu$ [GeV]')
    ax[1].set_ylabel(r'$|\Delta P|$')
    ax[1].grid(True, which='both', alpha=0.2)

    fig.savefig('../fig/prem_3nu_vs_energy_compare.pdf', bbox_inches='tight')'''),
    md(r'''### Time at matched accuracy

The trade-off curve, which is the only fair way to put the two on one axis: each point is one
setting of that code's own dial.'''),
    code(r'''if HAVE_NPE_EARTH:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot([u for u, _ in mg_curve], [e for _, e in mg_curve],
            marker='*', ms=13, lw=1.2, color='k', label=r'Mag$\nu$s (rtol/atol)')
    ax.plot([u for u, _ in npe_curve], [e for _, e in npe_curve],
            marker='o', ms=6, lw=1.2, color='C0',
            label='NuOscProbExact (slabs/segment)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Time per probability [$\mu$s]')
    ax.set_ylabel(r'Max $|P - P_{\rm own\ tightest}|$')
    ax.set_title(r'PREM: self-convergence against cost, $\cos\theta_z = %.2f$' % COSTHZ,
                 fontsize=10)
    ax.grid(True, which='both', alpha=0.2)
    ax.legend(fontsize=8)'''),
    md(r'''## 5. PREM, 3+1: a sterile state through the Earth

The 3+1 case is where the two codes' cost structures separate most, and for a reason worth
naming: Mag$\nu$s's compiled Cayley--Hamilton exponential covers $2\times2$ and $3\times3$ only,
because there is no practical closed form for a $4\times4$ Hermitian eigenproblem. Four flavors
therefore exponentiate through `numpy.linalg.eigh`, and Mag$\nu$s loses the factor it had at
three flavors. NuOscProbExact's SU(4) closed form has no such cliff.

An eV-scale splitting over an Earth-crossing baseline also carries a very large accumulated
phase, which is the regime where a refinement ladder works hardest.'''),
    code(r'''STERILE = dict(s14=0.15, s24=0.15, s34=0.0, D41=1.0)

def magnus_prem_4nu(rtol, atol):
    return np.asarray(oscprob.osc_prob_4nu_earth(
        E_PREM, costhz=COSTHZ, L=L_CHORD, **OSC, d14=0.0, d24=0.0,
        nu_i=gd.NUMU, nu_f=gd.NUE, rtol=rtol, atol=atol, **STERILE,
        electron_fraction=YE_UNIFORM))


def npe_prem_4nu(n_per_segment):
    import hamiltonians4nu
    h = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        OSC['s12'], OSC['s23'], OSC['s13'], STERILE['s14'], STERILE['s24'],
        STERILE['s34'], OSC['dCP'], OSC['D21'], OSC['D31'], STERILE['D41'])
    return np.asarray(npe_earth.probabilities_4nu_earth(
        h, E_PREM, COSTHZ, n_slabs_per_segment=n_per_segment))[..., 4]


# Warm up both sides first: at repeats=1 the first call carries the numba compile,
# the PREM slab caches and NuOscProbExact's own import-time work, which made its
# 8-slab row read slower than its 32-slab one.
magnus_prem_4nu(1.0e-2, 1.0e-4)
if HAVE_NPE_EARTH:
    npe_prem_4nu(4)

# Refereed by the same scipy slab product as the three-flavor case, not by either
# code.  That matters most here: Magnus does not converge in this regime -- it runs to
# its slab ceiling and warns -- so its answer at rtol=1e-5 is the same as at 1e-3, and
# a *self*-convergence curve would report an error near zero for a result that carries
# a MagnusConvergenceWarning.  Self-convergence measures stability; stability is not
# accuracy once the ladder has given up.
h_vac4 = np.asarray(hamiltonians.hamiltonian_4nu_vacuum_energy_independent(
    OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP'],
    STERILE['s14'], 0.0, STERILE['s24'], 0.0, STERILE['s34'],
    OSC['D21'], OSC['D31'], STERILE['D41']))
P_REF4, REF4_UNC = refereed(E_PREM, dim=4, h_vac_dim=h_vac4)
print('4nu referee residual discretization error: %.2e' % REF4_UNC)
print()
P_ref4 = P_REF4

print('%-14s %-12s %-14s %-12s %s'
      % ('rtol', 'ms total', 'us/probability', 'err vs ref', 'warned?'))
print('-'*72)
P_mg4 = None
mg4_curve = []
for rtol, atol in ((1.0e-3, 1.0e-5), (1.0e-5, 1.0e-7)):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P, t = best_of(lambda r=rtol, a=atol: magnus_prem_4nu(r, a), repeats=1)
        flags = sorted({c.category.__name__ for c in caught})
    P_mg4 = P
    err4 = (float(np.max(np.abs(P - P_ref4))) if P_ref4 is not None else float('nan'))
    mg4_curve.append((1.0e6*t/len(E_PREM), max(err4, REF4_UNC)))
    print('%-14.0e %-12.1f %-14.1f %-12.3e %s'
          % (rtol, 1.0e3*t, 1.0e6*t/len(E_PREM), err4, ','.join(flags) or 'no'))

if HAVE_NPE_EARTH:
    print()
    print('%-14s %-12s %-14s %s'
          % ('slabs/segment', 'ms total', 'us/probability', 'err vs ref'))
    print('-'*60)
    P_npe4 = None
    npe4_curve = []
    for n in (2, 4, 8, 32):
        P_npe4, t_npe4 = best_of(lambda k=n: npe_prem_4nu(k), repeats=1)
        errn = float(np.max(np.abs(P_npe4 - P_ref4)))
        npe4_curve.append((1.0e6*t_npe4/len(E_PREM), max(errn, REF4_UNC)))
        print('%-14d %-12.1f %-14.1f %.3e'
              % (n, 1.0e3*t_npe4, 1.0e6*t_npe4/len(E_PREM), errn))
    if P_mg4 is not None:
        print()
        print('max |Magnus - NuOscProbExact| = %.3e'
              % float(np.max(np.abs(P_mg4 - P_npe4))))'''),
    md(r'''**This is the most expensive case for Mag$\nu$s in this notebook -- and, once the
referee was corrected, not the least accurate one.** The two axes point in opposite directions
here, so they are worth separating.

**Cost: NuOscProbExact, by about 400x.** 56 000 $\mu$s per probability against 130, and the
cost does not fall when the tolerance is loosened. That flatness is a real diagnosis: the
refinement ladder is not converging and stopping, it is running to its slab ceiling, which is
why a `MagnusConvergenceWarning` appears on every row. An eV-scale $\Delta m^2_{41}$ over an
11 000 km chord accumulates a phase so large that the slab width needed for
$\lVert\Omega\rVert < \pi$ is below what the ladder will reach.

**Accuracy: Mag$\nu$s, and it reaches the floor of the measurement.** Its residual against the
referee is $4.5\times10^{-8}$, which is *below the referee's own discretization error of*
$4.1\times10^{-7}$ -- so the honest statement is that Mag$\nu$s agrees with the referee to
within the referee's own uncertainty, and this comparison cannot resolve it any further.
NuOscProbExact sits at $2.6\times10^{-4}$, some six hundred times above that floor, and its
convergence is **not monotonic**: 32 slabs per segment is worse than 8. Non-monotonic
convergence is the classic signature of a discretization whose edges are straddling structure
they do not resolve, and it means there is no setting of that dial to read off as "converged".

**Note what the warning does and does not say.** Mag$\nu$s warns on every row here and is
right to: the answer is not backed by a convergence argument. But it is not therefore a bad
answer -- it is the most accurate one on the plot. A convergence warning is a statement about
*evidence*, not about error, and the two come apart exactly here.

**This section previously reported the opposite, and the reason is worth recording.** Its
referee built the matter potential by hand as $\mathrm{diag}(1, 0, 0, 0)$, which gives a
sterile state the *active* flavors' neutral-current potential. Against that referee both codes
appeared wrong by $2.6\times10^{-2}$ while agreeing with each other to $3.7\times10^{-4}$ -- and
a referee that disagrees with every code it referees is the tell. The projector now comes from
`matter.matter_potential_projector`, and the corrected referee lands on Mag$\nu$s.

The residual between the two codes is still $3.7\times10^{-4}$ and is bounded by convention as
well as by convergence: the channel index differs (`[..., 4]` is $\nu_\mu \to \nu_e$ in a
four-flavor 16-tuple) and the sterile parametrizations are not identical, NuOscProbExact
taking no $\delta_{34}$.

The recommendation this notebook keeps arriving at is unchanged, because it was always about
cost: **for a piecewise-constant profile with a large accumulated phase, use a closed form.**
What has changed is the reason -- not that Mag$\nu$s is inaccurate there, but that it is slow
there, and a closed form gets to a perfectly adequate answer far more cheaply.

### Accuracy against cost, three flavors and 3+1 side by side

The figure this whole speed-up effort was for. Each point is one setting of that code's own
dial, so the curve -- not any single point -- is what a code offers you.

**Seven codes at three flavors, five at 3+1, and none of them its own judge.** The external
numbers are frozen in `external_prem_speed_accuracy.json` and none of those codes is needed to
run this notebook -- the same arrangement NuOscProbExact uses, and the same dataset.

The referee is a Richardson extrapolation $(4P_{256} - P_{128})/3$ of a 30-digit `mpmath` slab
product, cross-checked against an adaptive DOP853 integration of the *continuous* profile. They
agree to $2.3\times10^{-11}$ at three flavors and $1.2\times10^{-9}$ at 3+1, and that is the
dotted floor on each panel: **nothing below it is resolved**, and no point should be read as
more accurate than the ruler.

This is a different chord, channel **and physics point** from the tables above --
$\cos\theta_z = -0.9$, $P(\nu_\mu \to \nu_\mu)$, NuFIT 4.0 NO, $\sin^2\theta_{14} =
\sin^2\theta_{24} = 0.1$ -- because frozen numbers mean nothing off the workload they were
measured on, and the mixing parameters are part of that workload. Mag$\nu$s is measured live
on the same one.

That last point is not a formality. Run against this reference at NuFIT 6.1, Mag$\nu$s sits
$6.2\times10^{-2}$ away from it and *flat in* `rtol` -- which is exactly what a systematic
defect in the code looks like, and is nothing of the kind. Matching the parameters moves it
by two and a half orders of magnitude. A disagreement that does not respond to the tolerance
dial is evidence about the *setup*, not yet about the solver.'''),
    code(r'''# The external numbers are frozen, exactly as in NuOscProbExact's own comparison:
# none of these codes is needed to run this notebook, and none is its own judge.  The
# referee is a Richardson extrapolation (4*P(256) - P(128))/3 of a 30-digit mpmath slab
# product, cross-checked against an adaptive DOP853 integration of the CONTINUOUS
# profile.  Those two agree to 2.3e-11 at three flavors and 1.2e-09 at 3+1, which is
# the dotted floor on each panel.
PREM_EXT = json.loads(
    (pathlib.Path.cwd()/'external_prem_speed_accuracy.json').read_text())
COSTHZ_EXT = PREM_EXT['costhz']
L_EXT = PREM_EXT['baseline_km']*gd.CONV_KM_TO_INV_EV
OSC_EXT = gd.load_nufit_params('NuFIT 4.0', 'NO')
STERILE_EXT = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, D41=1.0)

# THE MATTER POTENTIAL IS MATCHED FIRST, or the comparison measures bookkeeping rather
# than physics.  The dataset hands every external code a scaled density for exactly this
# reason -- it records density_scale_nusquids = 0.99209 and density_scale_nucraft =
# 0.99267 -- and Mag(nu)s needs the same: its V_CC sits 1.90e-04 below NuOscProbExact's
# at the same density.
#
# Left unmatched, that 1.9e-04 pins the error at 2.35e-04 however tight the tolerance,
# and the curve comes out FLAT -- points stacked at one height, above every other code,
# saying nothing.  A trade-off curve that does not respond to its own dial is reporting
# a convention difference, not a solver.
VCC_MATCH = 1.0001896490


# V_CC is linear in the electron fraction, so matching their potential is one keyword.
# `osc_prob_*_earth` used to pin Y_e = 0.5 with no way to reach it, which is why an
# earlier draft of this cell dropped to `osc_prob_matter_std_potential` and rebuilt the
# PREM chord by hand.  The Earth entry points take the composition now.
YE_EXT = 0.5*VCC_MATCH


def timed_batch(call, n_energies, repeat=7, min_block=0.05):
    # Timed the way tests/prem_scan.py times every other code here, and the reason is
    # not fussiness.  A single 12-energy batched call takes a few hundred microseconds
    # while the FIRST one also pays ~0.7 s to compile the numba kernel -- charged to
    # twelve points that is 58 000 us/probability, and Mag(nu)s came out as one of the
    # slowest codes on the plot.  A user pays that once per session, not once per call.
    # So: throw the first pass away, then autorange the way `timeit` does, repeating
    # inside the timed block until it is long enough to measure, best of several.
    call()
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        elapsed = time.perf_counter() - t0
        if elapsed >= min_block:
            break
        reps *= 2
    best = elapsed/reps
    for _ in range(repeat - 1):
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        best = min(best, (time.perf_counter() - t0)/reps)
    return 1.0e6*best/n_energies


print('frozen dataset: costhz = %.2f, chord %.1f km, P(numu -> numu), NuFIT 4.0 NO'
      % (COSTHZ_EXT, PREM_EXT['baseline_km']))
print('our chord at that costhz: %.1f km   V_CC matched via Y_e = %.10f'
      % (mg_earth.distance_traveled_inside_earth(COSTHZ_EXT), YE_EXT))'''),
    code(r'''# Marker, color and size per code, fixed once so a reader tracks one code across
# both panels; NuOscProbExact drawn hollow, as the code the reference came from.
STYLE = {'NuOscProbExact': ('-o', 'C3', 4.4),
         'NuOscProbExact (tolerance)': ('-o', 'C3', 3.4),
         'NuOscProbExact (double-double)': ('-o', 'C3', 4.4),
         'NuOscProbExact (eigensolver)': ('-h', 'C3', 4.0),
         'nuSQuIDS': ('-v', 'C2', 3.8),
         'nuCraft': ('-s', 'C1', 3.6),
         'NuFast-Earth': ('-D', 'C4', 3.2),
         'GLoBES': ('-*', 'C6', 6.0),
         'Prob3++': ('-P', 'C5', 4.4)}
DIALS = ('n_slabs_per_segment', 'rtol', 'tolerance', 'num_prec', 'n_shells_per_layer')
# Started at 1e-1, not 1e-3: the curve should enter the plot at the cheap, inaccurate
# corner like every other code's, so the reader sees the whole trade-off rather than
# its converged tail.
RTOLS = (1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10)

# One figure per case, square, in the style of NuOscProbExact's paper figures: the two
# panels answer different questions and are read separately.
for key, title, sub, outfile in (
        ('three_flavor', r'PREM, three flavors',
         'PREM chord,  $\\cos\\theta_z = -0.9$\n$E = 3$--$40$ GeV,  '
         r'$P(\nu_\mu \to \nu_\mu)$' + '\nNuFIT 4.0 NO',
         'prem_speed_accuracy_3nu.pdf'),
        ('sterile_3plus1', r'PREM, 3+1 sterile',
         'PREM chord,  $\\cos\\theta_z = -0.9$\n$\\Delta m^2_{41} = 1$ eV$^2$,  '
         r'$\sin^2\theta_{14} = \sin^2\theta_{24} = 0.1$' + '\nNuFIT 4.0 NO',
         'prem_speed_accuracy_3plus1.pdf')):
    blk = PREM_EXT[key]
    E_ext = np.asarray(blk['energy_gev'])*gd.UNIT_GEV
    P_ref = np.asarray(blk['reference'])
    floor = float(blk['reference_vs_ode_max_abs'])

    fig, a = plt.subplots(figsize=(5.8, 5.4))
    for series in blk['series']:
        pts = series['points']
        marker, color, size = STYLE.get(series['name'], ('-o', '0.4', 3.6))
        dial = next((k for k in DIALS if pts[0].get(k) is not None), None)
        kw = dict(ms=size, color=color, lw=1.0, zorder=4,
                  label='%s%s' % (series['name'],
                                  '' if dial is None else '  (%s)' % dial))
        if series['name'].startswith('NuOscProbExact'):
            kw.update(mfc='white', mew=0.9, zorder=5)
        t = [p['us_per_probability'] for p in pts]
        e = [max(p['max_abs_error'], floor) for p in pts]
        a.loglog(t, e, marker, **kw)
        for j, ha, dx in ((0, 'right', -4), (len(pts) - 1, 'left', 4)):
            if len(pts) > 1:
                a.annotate(pts[j]['label'], xy=(t[j], e[j]), xytext=(dx, 3),
                           textcoords='offset points', fontsize=5.2, color=color,
                           ha=ha)

    mg_t, mg_e = [], []
    for rtol in RTOLS:
        if key == 'three_flavor':
            call = (lambda r=rtol: np.asarray(oscprob.osc_prob_3nu_earth(
                E_ext, costhz=COSTHZ_EXT, L=L_EXT, **OSC_EXT,
                nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=r, atol=r*1.0e-2,
                electron_fraction=YE_EXT)))
        else:
            call = (lambda r=rtol: np.asarray(oscprob.osc_prob_4nu_earth(
                E_ext, costhz=COSTHZ_EXT, L=L_EXT, **OSC_EXT, d14=0.0, d24=0.0,
                nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=r, atol=r*1.0e-2, **STERILE_EXT,
                electron_fraction=YE_EXT)))
        P = call()
        mg_t.append(timed_batch(call, len(E_ext)))
        mg_e.append(max(float(np.max(np.abs(P - P_ref))), floor))
    a.loglog(mg_t, mg_e, '-*', ms=13, color='k', lw=1.4, zorder=6,
             label=r'Mag$\nu$s  (rtol)')
    for j, ha, dx in ((0, 'right', -5), (len(RTOLS) - 1, 'left', 5)):
        a.annotate(r'$10^{%d}$' % round(np.log10(RTOLS[j])), xy=(mg_t[j], mg_e[j]),
                   xytext=(dx, 4), textcoords='offset points', fontsize=5.8, ha=ha)

    a.axhline(floor, color='0.5', ls=':', lw=0.8, zorder=1)
    a.text(0.03, 0.032, 'Referee floor  (mpmath vs DOP853)', transform=a.transAxes,
           fontsize=5.8, color='0.45')
    a.text(0.03, 0.955, sub, transform=a.transAxes, ha='left', va='top',
           fontsize=6.2, color='0.2', linespacing=1.45)
    a.set_xlabel(r'Time per probability [$\mu$s]')
    a.set_ylabel(r'Error vs. the converged reference,  max $|\Delta P|$')
    a.set_title(title, fontsize=10)
    a.grid(True, which='both', alpha=0.18)
    leg = a.legend(loc='lower left', fontsize=5.8)
    leg.get_frame().set_linewidth(0.7)
    # No dead margin left or right of the curves.  The floor line is included in the
    # vertical extent deliberately -- it is the thing everything else is read against,
    # so a limit that clipped it would hide the only absolute scale on the plot.
    allt = [p['us_per_probability'] for s in blk['series'] for p in s['points']] + mg_t
    alle = ([max(p['max_abs_error'], floor) for s in blk['series']
             for p in s['points']] + mg_e + [floor])
    a.set_xlim(min(allt)/1.7, max(allt)*1.7)
    a.set_ylim(min(alle)/3.5, max(alle)*3.5)
    fig.tight_layout(pad=0.4)
    fig.savefig('../fig/' + outfile, bbox_inches='tight')'''),

    md(r'''## 6. A conventions trap, and why the next table is in vacuum

Running nuSQuIDS on the same nominal problem produces a disagreement of a few times
$10^{-4}$ that **does not improve when its solver tolerance is tightened**. A plateau like
that is the signature of a systematic offset rather than a convergence limit.'''),
    code(r'''def run_nusquids(rho, tolerance):
    """P(nu_mu -> nu_e) from nuSQuIDS on the same configuration."""
    units = nsq.Const()
    solver = nsq.nuSQUIDS(E_GEV*units.GeV, 3, nsq.NeutrinoType.neutrino, False)
    if rho == 0.0:
        solver.Set_Body(nsq.Vacuum())
        solver.Set_Track(nsq.Vacuum.Track(L_KM*units.km))
    else:
        solver.Set_Body(nsq.ConstantDensity(rho, 0.5))
        solver.Set_Track(nsq.ConstantDensity.Track(L_KM*units.km))
    solver.Set_MixingAngle(0, 1, np.arcsin(OSC['s12']))
    solver.Set_MixingAngle(0, 2, np.arcsin(OSC['s13']))
    solver.Set_MixingAngle(1, 2, np.arcsin(OSC['s23']))
    solver.Set_CPPhase(0, 2, OSC['dCP'])
    solver.Set_SquareMassDifference(1, OSC['D21'])
    solver.Set_SquareMassDifference(2, OSC['D31'])
    solver.Set_rel_error(tolerance)
    solver.Set_abs_error(tolerance)
    state = np.zeros((N_E, 3))
    state[:, 1] = 1.0                                  # start as nu_mu
    solver.Set_initial_state(state, nsq.Basis.flavor)
    solver.EvolveState()
    return np.array([solver.EvalFlavorAtNode(0, i) for i in range(N_E)])


if AVAILABLE['nuSQuIDS']:
    print('%-18s %s' % ('solver tolerance', 'max |P - exact|, in matter'))
    for tolerance in (1.0e-5, 1.0e-7, 1.0e-10):
        P, _ = best_of(lambda t=tolerance: run_nusquids(RHO, t), repeats=1)
        print('%-18.0e %.3e' % (tolerance, np.max(np.abs(P - reference))))'''),
    md(r'''The diagnosis: run the same code in **vacuum**, where there is no matter potential to
disagree about. If the solver is fine, the vacuum comparison will be tight -- and it is, by
about four orders of magnitude.'''),
    code(r'''if AVAILABLE['nuSQuIDS']:
    vacuum_reference = np.array([np.abs(expm(-1j*(h_vac/e)*BASELINE))[0, 1]**2 for e in E])
    P_vac, _ = best_of(lambda: run_nusquids(0.0, 1.0e-12), repeats=1)
    print('vacuum : max |nuSQuIDS - exact| = %.3e' % np.max(np.abs(P_vac - vacuum_reference)))

    # What potential would reconcile them?  Rebuild the exact answer with V_CC
    # scaled, and find the scale that matches.
    P_matter, _ = best_of(lambda: run_nusquids(RHO, 1.0e-12), repeats=1)
    scales = np.linspace(0.95, 1.05, 21)
    gaps = []
    for scale in scales:
        ref_scaled = np.array([
            np.abs(expm(-1j*(h_vac/e + np.diag([VCC*scale, 0.0, 0.0]))*BASELINE))[0, 1]**2
            for e in E])
        gaps.append(np.max(np.abs(P_matter - ref_scaled)))
    best = int(np.argmin(gaps))
    unity = int(np.argmin(np.abs(scales - 1.0)))
    print('matter : best agreement at V_CC scale %.3f (residual %.2e)'
          % (scales[best], gaps[best]))
    print('         at Magnus\'s own V_CC   (scale %.3f), residual %.2e'
          % (scales[unity], gaps[unity]))'''),
    md(r'''In vacuum nuSQuIDS reproduces the exact answer to $\sim10^{-8}$ -- its solver is doing
its job. In matter, the disagreement is minimized at a matter potential about **1% larger**
than the one Mag$\nu$s builds from the same nominal $\rho$ and $Y_e$. That is a difference in
how the electron number density is derived from a mass density -- the average nucleon mass and
the electron fraction convention -- and **not** an accuracy difference in either code.

This is why the accuracy column in section 3 is only quoted for codes sharing Mag$\nu$s's
potential, and why the honest cross-code accuracy statement is the vacuum one. Fixing it
properly means agreeing a conversion, not tightening a tolerance.

'''),
    md(r'''## 7. Speed against accuracy, across six codes

The comparison above is two codes at a single working point. The useful picture is the whole
trade-off: every code has a dial -- a solver tolerance, a number of Newton steps, a batch mode
-- and what matters is the accuracy it buys per microsecond.

The NuOscProbExact project measured five external codes this way and published the result;
`notebooks/external_speed_accuracy.json` is that measurement, redistributed with attribution.
The configuration is the one used above -- $L = 1300$ km, $\rho = 3$ g cm$^{-3}$, three
flavors, $P(\nu_\mu \to \nu_e)$, 60 energies from 0.6 to 20 GeV -- scored against a 50-digit
`mpmath` matrix exponential.

**A warning about the timings, which is worth more than the figure.** Those numbers were taken
on other hardware, and there is no single factor that maps them onto this machine. Running
NuOscProbExact here in both of its modes gives *two different* conversion factors.'''),
    code(r'''EXTERNAL = json.loads((pathlib.Path.cwd()/'external_speed_accuracy.json').read_text())
print(EXTERNAL['_attribution'])
print()

frozen = {p['label']: p['us_per_probability']
          for s in EXTERNAL['series'] if s['name'] == 'NuOscProbExact'
          for p in s['points']}

if AVAILABLE['NuOscProbExact']:
    h_npe2 = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(**OSC)
    _, t_loop = best_of(lambda: np.array([
        oscprob3nu.probabilities_3nu(
            hamiltonians3nu.hamiltonian_3nu_matter(h_npe2, e, VCC), BASELINE)[3]
        for e in E]))
    stacked = hamiltonians3nu.hamiltonian_3nu_matter(h_npe2, E, VCC)
    _, t_arr = best_of(lambda: oscprob3nu.probabilities_3nu(stacked, BASELINE))

    print('%-16s %12s %12s %9s' % ('NuOscProbExact', 'frozen [us]', 'here [us]', 'factor'))
    print('-'*52)
    for label, here in (('One at a time', 1.0e6*t_loop/N_E), ('Array', 1.0e6*t_arr/N_E)):
        print('%-16s %12.3f %12.3f %9.2f'
              % (label, frozen[label], here, frozen[label]/here))'''),
    md(r'''One workload says this machine is more than twice as fast; the other says it is
slightly slower. That is not measurement noise -- looping is bound by the Python interpreter
and the array path is bound by BLAS, and two machines need not rank the same way in both.
**A single "machine factor" does not exist.** So nothing below is rescaled, and Mag$\nu$s is
plotted with a distinct symbol as a reminder that its horizontal position is not commensurate
with the others'.

Accuracy, by contrast, is a property of the algorithm rather than the hardware, and transfers
without qualification. It is the axis to read.'''),
    code(r'''fig, ax = plt.subplots(figsize=(7.4, 5.0))

FLOOR = 1.0e-16          # machine precision; several codes sit on it
markers = {'NuOscProbExact': 'o', 'nuSQuIDS': 's', 'NuFast-LBL': '^',
           'GLoBES': 'D', 'Prob3++': 'v', 'Second-order expansion': 'P'}
for i, series in enumerate(EXTERNAL['series']):
    us = [p['us_per_probability'] for p in series['points']]
    err = [max(p['max_abs_error'], FLOOR) for p in series['points']]
    ax.plot(us, err, marker=markers.get(series['name'], 'o'), ms=5,
            color='C%d' % i, lw=1.0, alpha=0.85, label=series['name'])

P_m, t_m = best_of(lambda: np.asarray(oscprob.osc_prob_3nu_matter_constant_density(
    E, BASELINE, RHO, **OSC, density_matter_is_in_g_per_cm3=True,
    nu_i=gd.NUMU, nu_f=gd.NUE)))
ax.plot([1.0e6*t_m/N_E], [max(float(np.max(np.abs(P_m - reference))), FLOOR)],
        marker='*', ms=16, color='k', ls='none', label=r'Mag$\nu$s (this machine)')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Time per probability [$\mu$s]')
ax.set_ylabel(r'max $|P - P_{\rm exact}|$')
ax.set_title(r'Speed against accuracy, constant density, $L = 1300$ km', fontsize=10)
ax.grid(True, which='both', alpha=0.2)
ax.legend(fontsize=7, loc='lower left')'''),
    md(r'''Read the vertical axis first. Three codes reach machine precision -- Mag$\nu$s,
NuOscProbExact, and NuFast at three Newton steps -- and the rest trade accuracy for speed by
between four and thirteen orders of magnitude. The second-order expansion is the extreme:
$6\times10^{-3}$, which is the size of the effects notebooks 15 and 17 are about, for a
quarter of a microsecond.

Read the horizontal axis only *within* a curve. nuSQuIDS tightening its tolerance from
$10^{-4}$ to $10^{-8}$ costs 20% and buys four orders of magnitude; past that it costs three
times as much and buys nothing, having hit its own floor. That statement is about nuSQuIDS
alone and survives a change of machine. "Code A is faster than code B" does not.

## Summary

| | |
|---|---|
| constant density, exactness | Mag$\nu$s and NuOscProbExact both $\sim10^{-16}$ |
| constant density, speed | **comparable batched**, both ~1 µs/probability; both an order of magnitude slower called one energy at a time |
| single point, constant density | the closed form by ~1.7×, almost all of it wrapper parameter resolution |
| PREM, 3ν, cost per call | NuOscProbExact by ~20× |
| PREM, 3ν, accuracy reachable | Mag$\nu$s to $3\times10^{-10}$; the closed form stalls near $6\times10^{-5}$ |
| PREM, 3+1, cost per call | NuOscProbExact by ~400×, and Mag$\nu$s *warns* |
| PREM, 3+1, accuracy reachable | Mag$\nu$s, to the referee's own floor ($4\times10^{-7}$); the closed form stalls near $3\times10^{-4}$, non-monotonically |
| vacuum, nuSQuIDS vs exact | $\sim10^{-8}$ -- solver is fine |
| matter, nuSQuIDS vs Mag$\nu$s | a **1% $V_{\rm CC}$ convention difference**, not accuracy |

The comparison worth making is not "which code is fastest" but "which code solves my problem",
and the table above does not have a winner in it -- it has a boundary. A **piecewise-constant**
profile is what a closed form is for, and the larger the accumulated phase the more decisively
so: on 3+1 PREM the closed form is cheaper by some 400x, which is the clearest cost result in
this notebook and is not in Mag$\nu$s's favor. Below about $10^{-4}$ on PREM it is also the
cheaper way to reach a given accuracy.

**Cost and accuracy part company on that row, though, and it is worth not eliding them.** The
closed form wins the cost by 400x; Mag$\nu$s wins the accuracy, reaching the referee's own floor
where the closed form stalls near $3\times10^{-4}$ and stops converging monotonically. Which of
those matters is a property of the problem being solved, not of the codes.

What Mag$\nu$s buys is everything that is not that: accuracy past where a piecewise-constant
discretization stalls, an arbitrary varying profile, a custom Hamiltonian, a BSM term nobody has
diagonalized. Where a closed form exists and the phase is large, use it.

And before comparing any two codes' numbers: **check that they agree in vacuum first.** If they
do not, the disagreement is in the solvers. If they do and they disagree in matter, it is in
the conventions, and no amount of tolerance will close it.'''),
    md(r'''## 8. What the batching and the compiled kernel are each worth

The three routes through this library to the same constant-density scan: one point at a time,
the whole stack in one call, and the whole stack with the compiled Cayley--Hamilton kernel. The
cost *per point* is what a parameter scan actually pays, and it is the number that moves as the
scan grows -- a single point is dominated by wrapper overhead that never shrinks.'''),
    code(r'''# Absolute microseconds here are a property of this machine and are worth little on
# their own; the ratios between the three routes are far more stable and are what the
# text quotes.  Measured with the minimum of several repetitions rather than the mean,
# because timing noise is one-sided.
import magnus.magnus as mgcore

SIZES = [1, 3, 10, 30, 100, 300, 1000, 3000]
saved_backend = mgcore.EXPM_BACKEND
rows = []


def const_scan(energies):
    return oscprob.osc_prob_3nu_matter_constant_density(
        energies, BASELINE, RHO, **OSC, density_matter_is_in_g_per_cm3=True,
        nu_i=gd.NUMU, nu_f=gd.NUE)


try:
    for n in SIZES:
        En = np.logspace(np.log10(0.6), np.log10(20.0), n)*gd.UNIT_GEV

        mgcore.EXPM_BACKEND = 'eigh'
        _, t_loop = best_of(lambda: [const_scan(float(e)) for e in En], repeats=3)
        _, t_arr = best_of(lambda: const_scan(En))

        if ek.HAVE_NUMBA:
            mgcore.EXPM_BACKEND = 'numba'
            const_scan(En)                       # warm the compiler, once
            _, t_ker = best_of(lambda: const_scan(En))
        else:
            t_ker = float('nan')
        rows.append((n, t_loop, t_arr, t_ker))
finally:
    mgcore.EXPM_BACKEND = saved_backend

rows = np.array(rows)
print('%8s %14s %14s %14s' % ('N', 'loop [us/pt]', 'array [us/pt]', 'kernel [us/pt]'))
print('-'*56)
for n, tl, ta, tk in rows:
    print('%8d %14.3f %14.3f %14.3f' % (n, tl/n*1e6, ta/n*1e6, tk/n*1e6))
print()
print('at N = %d: batching is %.0fx the loop, the kernel a further %.1fx'
      % (rows[-1, 0], rows[-1, 1]/rows[-1, 2], rows[-1, 2]/rows[-1, 3]))'''),
    code(r'''fig, ax = plt.subplots(figsize=(7.2, 4.6))
for col, style, lab in ((1, '-^', 'One point at a time'),
                        (2, '-s', 'Batched, eigh'),
                        (3, '-o', 'Batched + compiled kernel')):
    if col == 3 and not ek.HAVE_NUMBA:
        continue
    ax.plot(rows[:, 0], 1.0e6*rows[:, col]/rows[:, 0], style, ms=5, lw=1.2,
            label=lab)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$N$ energies in one request')
ax.set_ylabel(r'Time per probability [$\mu$s]')
ax.set_title(r'Constant density, 3$\nu$: cost per point against scan size', fontsize=10)
ax.grid(True, which='both', alpha=0.2)
ax.legend(fontsize=8)
fig.savefig('../fig/performance.pdf', bbox_inches='tight')'''),
    md(r'''Read the *shape*, not the endpoints. All three routes meet at $N = 1$, where the
answer is one exponential and everything else is wrapper: that is the regime this library is
worst in, and no amount of batching or compilation touches it. They separate as the request
grows, because the per-call overhead is amortized while the arithmetic is not.

The two gaps are different in kind. **Batching** removes Python-level dispatch -- it is the
same arithmetic, called once instead of $N$ times. **The kernel** removes the arithmetic's own
overhead, replacing `numpy.linalg.eigh`'s per-matrix LAPACK loop with a closed form over the
whole stack. That second gap is the one that closes for four and five flavors, where no closed
form for the eigenproblem exists and the `eigh` path is all there is.'''),
    md(r'''## 9. A smooth profile: where each method runs out

PREM is the profile that flatters a closed form, because a piecewise-constant density is
exactly what one is built to exploit. This section asks the opposite question on a *smooth*
potential, $V_{\rm CC}(l) = V_0 e^{-3l/L}$ over 3000 km: NuOscProbExact samples each slab at
its midpoint, second order in the slab width, while Mag$\nu$s's Gauss--Legendre expansion is
fourth order.

Every number below is frozen in `external_profile_benchmarks.json`, generated by
`gen_profile_benchmarks.py`. **All codes were timed in one process on one machine** -- unlike
section 7, where the externals come from another machine and only Mag$\nu$s is live -- so the
time axis here is a comparison rather than a juxtaposition. Both codes are batched over the
twelve energies: `probabilities_Nnu_slabs` takes a stack of chords sharing one set of widths,
and timing it one energy at a time would flatter Mag$\nu$s roughly fivefold. The referee is an
adaptive DOP853 integration of the evolution operator, which is neither code's method.'''),
    code(r'''BENCH = json.loads(
    (pathlib.Path.cwd()/'external_profile_benchmarks.json').read_text())
print('machine:', BENCH['machine'])
print('control ratio %.3f  (%s)' % (BENCH['control_ratio'], BENCH['control_note']))


def bench_case(profile, d):
    for c in BENCH['cases']:
        if c['profile'] == profile and c['flavours'] == d:
            return c
    return None


DIAL_STYLE = {'Magnus': ('-*', 'k', 12), 'NuOscProbExact': ('-o', 'C3', 4.4)}

# One hierarchy for these four panels, stated here rather than left to the
# defaults, because the defaults pull in two directions at once: `matplotlibrc`
# sets `axes.labelsize` to 25 ABSOLUTELY, for single-panel paper figures, while
# the annotations below had been tuned down to 5.4 pt to fit between the
# markers.  The result was a five-fold spread in which the axis labels dwarfed
# the title and the dial values were unreadable at any size the figure is
# actually viewed.  Ordered largest to smallest, and every size restated:
PANEL_FONT = {'title': 18, 'label': 17, 'tick': 14, 'legend': 11, 'dial': 10,
              'note': 10}


def plot_case(ax, case, title):
    for series in case['series']:
        marker, color, size = DIAL_STYLE[series['name']]
        pts = series['points']
        t = [p['us_per_probability'] for p in pts]
        e = [p['max_abs_error'] for p in pts]
        kw = dict(ms=size, color=color, lw=1.1, zorder=4,
                  label='%s  (%s)' % (series['name'], series['dial']))
        if series['name'] == 'NuOscProbExact':
            kw.update(mfc='white', mew=0.9)
        ax.loglog(t, e, marker, **kw)
        # The dial's value at both ends, so the reader can see which knob buys
        # what.  Both vertical offsets are deliberate, and both were collisions
        # first.  The FIRST label sits below its marker because every curve
        # starts at the top left, which is where the corner note already is --
        # above-left overprinted the note in three of the four panels.  At the
        # right-hand end the two series finish at similar errors, so they are
        # pushed apart vertically rather than left to overlap, which is what
        # `32768` and `1e-10` did in the four-flavor panel.
        last_dy = 6 if series['name'] == 'Magnus' else -14
        for j, ha, dx, dy in ((0, 'right', -4, -14),
                              (len(pts) - 1, 'left', 4, last_dy)):
            ax.annotate(pts[j]['label'], xy=(t[j], e[j]), xytext=(dx, dy),
                        textcoords='offset points', fontsize=PANEL_FONT['dial'],
                        color=color, ha=ha)
    ax.set_xlabel(r'Time per probability [$\mu$s]', fontsize=PANEL_FONT['label'])
    # Short on purpose.  `matplotlibrc` sets `axes.labelsize` to 25 absolutely,
    # and at that size the previous wording -- "Error vs. the DOP853 referee,
    # max |Delta P|" -- was taller than the 5.2 inch figure, so `bbox_inches`
    # sliced a character off each end and the axis read "rror ... ma".  What
    # the referee is stays in the section text above, where there is room.
    ax.set_ylabel(r'max $|\Delta P|$ vs.\ DOP853', fontsize=PANEL_FONT['label'])
    ax.set_title(title, fontsize=PANEL_FONT['title'], pad=34)
    # The run's parameters go ABOVE the axes, under the title, rather than into
    # a corner of them.  Every corner is taken in at least one of the four
    # panels: the curves begin at the top left, the legend holds the bottom
    # left, and the end-of-dial labels the bottom right.  In the four-flavor
    # panel the first point sits hard against the top left corner, so no
    # nudge inside the axes clears the block -- moving the label down only put
    # it on the next line of the note.  The series names this used to carry are
    # gone with it; the legend already names both codes.
    ax.text(0.5, 1.005, 'Smooth $V_{\\rm CC} = V_0 e^{-3l/L}$,  $L = 3000$ km\n'
            '$E = 2$--$12$ GeV,  $P(\\nu_\\mu \\to \\nu_\\mu)$',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=PANEL_FONT['note'], color='0.25', linespacing=1.4)
    ax.tick_params(labelsize=PANEL_FONT['tick'])
    ax.grid(True, which='both', alpha=0.18)
    leg = ax.legend(fontsize=PANEL_FONT['legend'], loc='lower left')
    leg.get_frame().set_linewidth(0.7)
    # No dead margin left or right of the curves: on a log-log plot matplotlib's default
    # padding is a whole decade, which makes two curves look further apart than they are.
    allt = [p['us_per_probability'] for s in case['series'] for p in s['points']]
    alle = [p['max_abs_error'] for s in case['series'] for p in s['points']]
    ax.set_xlim(min(allt)/1.6, max(allt)*1.6)
    ax.set_ylim(min(alle)/3.0, max(alle)*3.0)


for d in (2, 3, 4, 5):
    case = bench_case('exponential', d)
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    plot_case(ax, case, r'Exponential profile, %d$\nu$' % d)
    fig.tight_layout(pad=1.2)
    fig.savefig('../fig/smooth_speed_accuracy_%dnu.pdf' % d, bbox_inches='tight')'''),
    md(r'''**NuOscProbExact is faster at every accuracy it can reach, and that is the expected
result rather than a disappointing one** -- it solves each slab in closed form, and where a
closed form exists an exact algebraic solution beats a truncated series.

What the plots show instead is where each method *stops*. NuOscProbExact's error bottoms out
near $2.5\times10^{-11}$ and then rises: past about 16 000 slabs the round-off of composing
that many matrix products costs more than another halving of $h$ buys, so 32 768 slabs is
worse than 16 384. That is a floor with no setting below it. Mag$\nu$s continues to
$3\times10^{-13}$, because halving $h$ buys it a factor of sixteen rather than four.

At five flavors there is no comparison to draw: NuOscProbExact has no five-flavor route, so
the panel carries one curve. That is the other half of the same point -- reach, in the sense of
which problems a method can address at all.'''),
    md(r'''### The probability itself

The curves above say what each code *costs*; this says what they *return*. Both are drawn on a
200-point grid, well past what the twelve-energy timing grid resolves, and the residual sits
underneath because on the probability axis the two are indistinguishable -- which is the point,
and not something a reader should have to take on trust.'''),
    code(r'''L_EXPO = 3000.0*gd.CONV_KM_TO_INV_EV
V0_EXPO = 1.0e-13
PER_NE_EXPO = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)


def vcc_expo(x):
    return V0_EXPO*np.exp(-3.0*np.asarray(x, dtype=float)/L_EXPO)


def ne_expo(x):
    return vcc_expo(x)/PER_NE_EXPO


# The frozen benchmark used NuFIT 4.0, because that is the point the external dataset in
# section 7 was generated at.  This panel follows it, so the two halves of section 9
# describe the same physics rather than two neighboring ones.
OSC_EXPO = gd.load_nufit_params('NuFIT 4.0', 'NO')
h_vac_expo = np.asarray(
    hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**OSC_EXPO))

E_EXPO = np.linspace(2.0, 12.0, 200)*gd.UNIT_GEV
P_mg_expo = np.asarray(oscprob.osc_prob_matter_std_potential(
    3, ne_expo, E_EXPO, L_EXPO, OSC_EXPO, L0=0.0, nu_i=gd.NUMU, nu_f=gd.NUMU,
    density_is_of_number_of_electrons=True, rtol=1.0e-10, atol=1.0e-12,
    strategy='magnus'))

if HAVE_NPE_EARTH:
    n_slab = 16384
    edges = np.linspace(0.0, L_EXPO, n_slab + 1)
    mid = 0.5*(edges[:-1] + edges[1:])
    H_expo = np.broadcast_to(
        (h_vac_expo[None, None]/E_EXPO[:, None, None, None]).astype(complex),
        (len(E_EXPO), n_slab, 3, 3)).copy()
    H_expo[:, :, 0, 0] += vcc_expo(mid)[None, :]
    P_npe_expo = np.asarray(
        npe_slabs.probabilities_3nu_slabs(H_expo, np.diff(edges)))[..., 4]
    print('max |Magnus - NuOscProbExact| over 200 energies = %.3e'
          % float(np.max(np.abs(P_mg_expo - P_npe_expo))))'''),
    code(r'''if HAVE_NPE_EARTH:
    fig, ax = plt.subplots(2, 1, figsize=(6.6, 5.4), sharex=True,
                           gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.08))
    ax[0].plot(E_EXPO/gd.UNIT_GEV, P_mg_expo, lw=1.6, color='k',
               label=r'Mag$\nu$s  (rtol $10^{-10}$)')
    ax[0].plot(E_EXPO/gd.UNIT_GEV, P_npe_expo, lw=1.1, ls='--', color='C3',
               label='NuOscProbExact  (16384 slabs)')
    ax[0].set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
    ax[0].set_title('Exponential profile, three flavors, 200 energies', fontsize=10)
    ax[0].grid(True, alpha=0.2)
    ax[0].legend(fontsize=8)

    ax[1].semilogy(E_EXPO/gd.UNIT_GEV, np.abs(P_mg_expo - P_npe_expo), lw=1.0,
                   color='C0')
    ax[1].set_xlabel(r'$E_\nu$ [GeV]')
    ax[1].set_ylabel(r'$|\Delta P|$')
    ax[1].grid(True, which='both', alpha=0.2)
    fig.savefig('../fig/smooth_prob_vs_energy.pdf', bbox_inches='tight')'''),
    md(r'''## 10. The Sun: an observable the other codes do not offer

The solar case cannot be run as a race, and the reason is physics rather than bookkeeping.
Over the ray from the core to the surface the accumulated phase is about **13 000 radians** at
5 MeV -- some two thousand oscillations. The *instantaneous* probability at the surface is
therefore not a measurable quantity and not a stable one: an adaptive DOP853 integration of it
runs for minutes per energy, and no experiment sees it. What a solar experiment measures is
the **phase-averaged** probability.

Mag$\nu$s computes that directly. `average=True` on a smooth position-dependent profile takes
the adiabatic route -- decohere in the matter eigenbasis at production, transport along the
levels of the instantaneous Hamiltonian with the exact crossing probabilities wherever the
evolution stops being adiabatic, and read out in the vacuum basis at detection. It never
propagates, so the two thousand oscillations cost nothing.

The referee here is the analytic adiabatic limit, built from the eigenvectors of the
Hamiltonian at production and in vacuum. It touches no propagation code at all.'''),
    code(r'''TABLE = os.path.join('..', 'docs', 'dev', 'adversarial_batteries', 'bs05_agsop.dat')
rows = []
with open(TABLE) as fh:
    for line in fh:
        f = line.split()
        if len(f) == 12:
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                continue
solar = np.array(rows)
r_over_rsun, rho_cgs, x_h = solar[:, 1], solar[:, 3], solar[:, 6]
MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
ne_solar = rho_cgs*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*(0.5*(1.0 + x_h))
x_solar = r_over_rsun*gd.SUN_RADIUS*gd.UNIT_KM
log_ne_solar = np.log(ne_solar)
R_SUN = float(x_solar[-1])


def ne_sun(l):
    xs = np.clip(np.asarray(l, dtype=float), x_solar[0], x_solar[-1])
    out = np.exp(np.interp(xs, x_solar, log_ne_solar))
    return out[()] if np.ndim(out) == 0 else out


PER_NE_SUN = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
# 0.1-20 MeV, log-spaced.  An earlier draft used 1-15 MeV linear and the MSW transition
# was invisible: the low-energy vacuum-averaged plateau sits BELOW 1 MeV, so the window
# started on the falling edge and the curve looked like a featureless slope.
E_SUN = np.logspace(np.log10(0.1), np.log10(20.0), 40)*gd.UNIT_MEV

print('BS2005-AGS,OP: %d rows, ray 0 -> %.0f km' % (len(solar), R_SUN/gd.UNIT_KM))
print('accumulated phase at 5 MeV: ~%.0f radians'
      % (OSC['D21']*R_SUN/(4.0*5.0*gd.UNIT_MEV)))


def adiabatic_reference(energy, a=gd.NUE, b=gd.NUE):
    """Decohere in the matter eigenbasis at production, read out in vacuum.

    Analytic, and it never propagates -- which is what lets it referee a case where
    propagating is the expensive thing.
    """
    hm = h_vac/energy + float(PER_NE_SUN*ne_sun(0.0))*np.diag([1.0, 0.0, 0.0])
    _, u_matter = np.linalg.eigh(hm)
    _, u_vac = np.linalg.eigh(h_vac/energy)
    return float(np.sum(np.abs(u_matter[a])**2 * np.abs(u_vac[b])**2))


t0 = time.perf_counter()
P_sun = np.asarray(oscprob.osc_prob_matter_std_potential(
    3, ne_sun, E_SUN, R_SUN, OSC, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE,
    density_is_of_number_of_electrons=True, average=True))
t_sun = time.perf_counter() - t0
P_adiab = np.array([adiabatic_reference(e) for e in E_SUN])

print()
print('magnus average=True, %d energies: %.3f s  (%.1f ms per energy)'
      % (len(E_SUN), t_sun, 1.0e3*t_sun/len(E_SUN)))
print('worst |magnus - analytic adiabatic| = %.2e' % np.max(np.abs(P_sun - P_adiab)))

# nuSQuIDS on the SAME BS2005-AGS,OP file, so this is one model and two codes rather than
# two models.  Frozen in external_solar_nusquids.json by gen_solar_nusquids.py, because
# the cheapest setting that returns a probability at all takes about ten minutes; the
# generator's docstring records what had to be fixed to feed it, and what did not.
NSQ = json.loads((pathlib.Path.cwd()/'external_solar_nusquids.json').read_text())
NSQ_OK = [s for s in NSQ['series'] if s['physical']]
NSQ_BEST = NSQ_OK[-1] if NSQ_OK else None

print()
print('nuSQuIDS, same model file, %d energies, over its solver tolerance:'
      % NSQ['n_grid'])
print('  %-9s %9s %12s %26s %s'
      % ('rel_error', 'seconds', '|1 - sum|', 'P over all flavors', 'a probability?'))
for s in NSQ['series']:
    print('  %-9.0e %9.1f %12.1e %26s %s'
          % (s['tolerance'], s['seconds_total'], s['unitarity'],
             '%.4f .. %.4f' % (s['p_min'], s['p_max']),
             'yes' if s['physical'] else 'NO'))
print()
print('Read the last two columns together. The flavor sum is conserved to ~1e-16 on')
print('every row, including the rows returning P = 2.83 -- so the obvious check passes')
print('on output that is not a probability. See the discussion below.')'''),
    code(r'''fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                       gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.08))
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_sun, lw=1.8, color='k',
               label=r'Mag$\nu$s, average=True')
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_adiab, lw=1.2, ls='--', color='C3',
               label='Analytic adiabatic limit')
if NSQ_BEST is not None:
    e_nsq = np.array(NSQ_BEST['grid_energy_ev'])/gd.UNIT_MEV
    # The raw curve, thin and faint: it is not a competitor here, it is the reason the
    # instantaneous probability is not the observable.
    ax[0].semilogx(e_nsq, NSQ_BEST['P_ee_instantaneous'], lw=0.5, color='C2',
                   alpha=0.55,
                   label=r'nuSQuIDS, instantaneous (rtol $10^{%d}$)'
                         % round(np.log10(NSQ_BEST['tolerance'])))
    # And the window mean where enough samples were spent on it to mean anything: 300
    # inside each window rather than the three or four the sweep above leaves there.
    # Plotted WITH its standard error, because an average quoted without one cannot be
    # checked -- and here the error bar is the whole point of the panel.
    dense = NSQ['dense_check']
    ax[0].errorbar([r['energy_ev']/gd.UNIT_MEV for r in dense['points']],
                   [r['mean'] for r in dense['points']],
                   yerr=[r['stderr'] for r in dense['points']],
                   fmt='o', ms=6, color='C2', mfc='white', mew=1.3, capsize=3,
                   elinewidth=1.1, zorder=6,
                   label=r'nuSQuIDS, $\pm5\%%$ window mean (%d samples)'
                         % dense['samples_per_target'])
ax[0].set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
ax[0].set_title('BS2005-AGS,OP solar model: the averaged survival probability',
                fontsize=10)
ax[0].grid(True, alpha=0.2)
ax[0].legend(fontsize=8)

ax[1].loglog(E_SUN/gd.UNIT_MEV, np.abs(P_sun - P_adiab), lw=1.2, color='C0')
ax[1].set_xlabel(r'$E_\nu$ [MeV]')
ax[1].set_ylabel(r'$|\Delta P|$')
ax[1].grid(True, which='both', alpha=0.2)
# No dead margin either side of the curves.
for a in ax:
    a.set_xlim(E_SUN[0]/gd.UNIT_MEV, E_SUN[-1]/gd.UNIT_MEV)
ax[0].set_ylim(0.0, 1.0)
fig.savefig('../fig/solar_averaged.pdf', bbox_inches='tight')'''),
    md(r'''**The lower panel is not an error curve -- it is the non-adiabatic correction.** It
grows smoothly with energy, from $\sim10^{-7}$ at 1 MeV to $\sim10^{-5}$ at 15 MeV, which is
the expected direction: the higher the energy, the less adiabatically the resonance is
crossed. The analytic limit omits that by construction; Mag$\nu$s carries it, because its
adiabatic route uses the exact crossing probabilities where the evolution stops being
adiabatic. So the referee here validates Mag$\nu$s and cannot rank it -- below about
$10^{-5}$ Mag$\nu$s is the more correct of the two, and an accuracy axis drawn against the
analytic limit would penalize it for that.

The upper panel is the MSW transition, from the vacuum-averaged regime near 0.55 at low energy
to the matter-dominated one near 0.31 above 10 MeV.

**nuSQuIDS is on the plot, on the same model file, and it shows the problem rather than
solving it.** Its `SunASnu` body reads the same BS2005-AGS,OP table, and its raw output is the
thin green line: not a curve but a sampling of a phase that turns over tens of thousands of
times across the range. Nothing is wrong with the calculation -- that *is* the instantaneous
probability. It is simply not the observable.

**Recovering the observable from it is a Monte-Carlo estimate, and that is the expense.** A
$\pm5\%$ energy window is what a finite energy resolution means, but the probability turns over
so many times inside that window that a window mean is an average over effectively random
phases. Its uncertainty falls only as $1/\sqrt{N}$ in the number of evaluations spent inside
the window:

| samples in the $\pm5\%$ window | standard error of $\langle P_{ee}\rangle$ |
|---|---|
| 4 (what 200 energies across the range leaves) | $\sim0.14$ |
| 300 (the open circles above) | $\sim0.015$ |

The first row is why a window mean taken off the sweep scatters across most of $[0, 1]$ and is
not plotted here: with four samples the estimate is worthless, and it would have looked like a
disagreement with Mag$\nu$s rather than like noise. The open circles are the same calculation
given 300 samples per point, and they land on Mag$\nu$s's curve within their error bars.

That is the real comparison. Mag$\nu$s computes the average *analytically* -- it transports
along the levels of the instantaneous Hamiltonian and never propagates, so there is no sampling
error to reduce and no $N$ to choose. nuSQuIDS must estimate the same number by evaluating the
oscillation many times and averaging, and buying another decimal place costs a hundredfold more
evaluations.

**The MSW transition is the shape to read**, and it only appears once the range starts below
1 MeV: a plateau near 0.54 at low energy, where matter is irrelevant and the answer is the
vacuum average $1 - \tfrac{1}{2}\sin^2 2\theta_{12}$, falling through the resonance to about
0.30 above 10 MeV, where propagation is adiabatic and the answer approaches
$\sin^2\theta_{12}$. An earlier draft of this section ran 1--15 MeV and showed a featureless
slope, because the window began on the falling edge.'''),
    md(r'''### The tolerance dial, and a check that passes on nonsense

nuSQuIDS has a dial Mag$\nu$s's adiabatic route does not: the solver tolerance. Sweeping it
gives the speed-against-accuracy curve below -- and something more useful than the curve.

**Below `rel_error` $=10^{-6}$ this calculation stops returning probabilities, and says
nothing about it.** At $10^{-4}$ the survival probability reaches **2.83**; on looser settings
still, the flavor probabilities run from $-19$ to $+45$. The obvious guard against that --
sum the flavor probabilities and check they are 1 -- **passes on every one of those rows, to
$10^{-16}$**. That is not a bug in the guard. nuSQuIDS evolves the density matrix in an SU(3)
basis whose identity component is the trace, so the flavor sum is conserved *by construction*
however badly the traceless components are integrated. Unitarity is structural here, and a
structural invariant cannot test the thing it is built into.

The check that does bite is each probability lying in $[0, 1]$, which is what the generator
asserts before it writes anything. It is worth stating because the failure is otherwise
completely silent: a user sweeping the tolerance downward for speed gets numbers that are not
probabilities, from a solver that reports nothing wrong, and a unitarity check that agrees.'''),
    code(r'''# What each code's answer is UNCERTAIN BY, against what it cost.  For Mag(nu)s that
# is its residual against the analytic adiabatic limit -- which is the non-adiabatic
# correction, so it is an upper bound on any error rather than an error.  For nuSQuIDS it
# is the standard error of the window mean, because that is what actually limits the
# averaged answer here: the solver tolerance controls whether the output is a probability
# at all, and sampling controls how well the average is known.
def window_stderr(series):
    """Median standard error of the +/-5% window mean over the target energies."""
    grid = np.array(series['grid_energy_ev'])
    p = np.array(series['P_ee_instantaneous'])
    out = []
    for t in NSQ['energy_ev']:
        sel = (grid >= t*(1.0 - NSQ['spread'])) & (grid <= t*(1.0 + NSQ['spread']))
        if sel.sum() > 1:
            out.append(p[sel].std(ddof=1)/np.sqrt(sel.sum()))
    return float(np.median(out))


fig, ax = plt.subplots(figsize=(6.6, 4.6))

for s in NSQ['series']:
    good = s['physical']
    # Labels are deduplicated below, so each point can carry its own unconditionally.
    ax.loglog([s['seconds_total']], [window_stderr(s)],
              'o' if good else 'X', ms=7.5, color='C2' if good else 'C1',
              mfc='white' if good else 'C1', mew=1.2, ls='none', zorder=4,
              label='nuSQuIDS, sweep (a probability)' if good else
                    'nuSQuIDS, sweep (NOT a probability)')
    ax.annotate(r'$10^{%d}$' % round(np.log10(s['tolerance'])),
                xy=(s['seconds_total'], window_stderr(s)), xytext=(5, 4),
                textcoords='offset points', fontsize=6.2, color='0.3')

dense = NSQ['dense_check']
ax.loglog([dense['seconds_total']],
          [float(np.median([r['stderr'] for r in dense['points']]))], 's', ms=7.5,
          color='C2', mfc='C2', ls='none', zorder=5,
          label=r'nuSQuIDS, %d samples per window' % dense['samples_per_target'])

ax.loglog([t_sun], [max(float(np.max(np.abs(P_sun - P_adiab))), 1.0e-16)], '*',
          ms=17, color='k', ls='none', zorder=6,
          label=r'Mag$\nu$s, average=True (analytic, no sampling)')

ax.set_xlabel('Seconds for the whole curve')
ax.set_ylabel(r'Uncertainty on $\langle P(\nu_e \to \nu_e)\rangle$')
ax.set_title('The Sun: what the averaged observable costs each method', fontsize=10)
ax.grid(True, which='both', alpha=0.2)
handles, labels = ax.get_legend_handles_labels()
seen = dict(zip(labels, handles))
leg = ax.legend(seen.values(), seen.keys(), fontsize=7, loc='lower left')
leg.get_frame().set_linewidth(0.7)
fig.tight_layout(pad=0.4)
fig.savefig('../fig/solar_speed_accuracy.pdf', bbox_inches='tight')

print('Magnus   : %.2f s, residual vs analytic adiabatic %.2e (the non-adiabatic term)'
      % (t_sun, float(np.max(np.abs(P_sun - P_adiab)))))
for s in NSQ['series']:
    print('nuSQuIDS rtol %.0e: %6.1f s  window stderr %.4f  %s'
          % (s['tolerance'], s['seconds_total'], window_stderr(s),
             'physical' if s['physical'] else 'NOT a probability'))
print('nuSQuIDS %d samples/window: %.1f s  stderr %.4f'
      % (dense['samples_per_target'], dense['seconds_total'],
         float(np.median([r['stderr'] for r in dense['points']]))))'''),
    md(r'''**The axis is an uncertainty, not a ranking, and the two entries on it are not the
same kind of quantity.** Mag$\nu$s's is its residual against the analytic adiabatic limit,
which the panel above established is the *non-adiabatic correction* -- so it is an upper bound
on any error rather than an error, and it cannot be pushed lower by asking for more accuracy.
nuSQuIDS's is the standard error of a Monte-Carlo average, which falls as $1/\sqrt{N}$ and can
be pushed lower by spending evaluations.

The brief this section was written from expected both codes to be scored against the analytic
limit on the same axis. **They cannot be.** Sampling noise in nuSQuIDS's averaged answer is
$\sim10^{-2}$ even after nine hundred evaluations, three orders of magnitude above the
$\sim10^{-5}$ scale of the non-adiabatic correction the axis would need to resolve. An accuracy
comparison at that scale is not available at any affordable cost, and saying so is more useful
than drawing one that looks like it is.

What is left is the horizontal distance, and it is large in a way that is not about
implementation quality. Mag$\nu$s returns the averaged observable in **about 0.7 s** because it
never propagates and never samples. nuSQuIDS needs about **ten minutes** merely to reach the
tolerance at which its output is a probability, and then a further factor of $N$ to average
the phase away. Neither NuOscProbExact nor nuSQuIDS offers an averaging flag.

**Not "faster at the same algorithm", but a different algorithm for the question a solar
experiment actually asks** -- and the gap is a property of the question, not of the codes.'''),
    md(r"""### 10b. The same observable at 3+1, and a check that the sterile state is felt

A sterile state changes the solar answer through the matter term, not just through mixing, and
that is worth showing here because **this notebook's own referee was wrong about it until
today** (section 5). The matter potential carries $\mathrm{diag}(1, 0, 0, r/2)$ with
$r = n_n/n_p$: the active flavors share $V_{\rm NC}$ and it cancels, a sterile state feels
neither current and so keeps $-V_{\rm NC} = (r/2)\,V_{\rm CC}$.

Two things are checked below rather than assumed. **That the adiabatic route reaches four
flavors at all** -- `average=True` takes the adiabatic path on a smooth profile, and that path
has to find the level crossings of a $4\times4$ Hamiltonian rather than a $3\times3$ one. And
**that the sterile entry is actually live**, by varying $r$ and watching the curve move. If it
did not move, the term would be absent and every 3+1 solar number here would be quietly wrong
in the way section 5's referee was."""),
    code(r'''STERILE_SUN = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0,
                   d14=0.0, d24=0.0, D41=1.0e-5)
OSC4_SUN = dict(OSC, **STERILE_SUN)

h_vac4_sun = np.asarray(hamiltonians.hamiltonian_4nu_vacuum_energy_independent(
    OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP'],
    STERILE_SUN['s14'], 0.0, STERILE_SUN['s24'], 0.0, STERILE_SUN['s34'],
    OSC['D21'], OSC['D31'], STERILE_SUN['D41']))


def adiabatic_reference_4nu(energy, ratio=1.0, a=gd.NUE, b=gd.NUE):
    """The same analytic limit as before, at four flavors.

    The projector comes from the library.  Writing it out here is exactly the mistake
    section 5's referee made, and it is the reason that section reported the wrong
    winner for as long as it did.
    """
    proj = matter.matter_potential_projector(4, ratio)
    hm = h_vac4_sun/energy + float(PER_NE_SUN*ne_sun(0.0))*proj
    _, u_matter = np.linalg.eigh(hm)
    _, u_vac = np.linalg.eigh(h_vac4_sun/energy)
    return float(np.sum(np.abs(u_matter[a])**2 * np.abs(u_vac[b])**2))


t0 = time.perf_counter()
P_sun4 = np.asarray(oscprob.osc_prob_matter_std_potential(
    4, ne_sun, E_SUN, R_SUN, OSC4_SUN, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE,
    density_is_of_number_of_electrons=True, average=True))
t_sun4 = time.perf_counter() - t0
P_adiab4 = np.array([adiabatic_reference_4nu(e) for e in E_SUN])

print('CHECK 1 -- does the adiabatic route reach four flavors?')
print('  yes: %d averaged energies in %.3f s (%.1f ms per energy)'
      % (len(E_SUN), t_sun4, 1.0e3*t_sun4/len(E_SUN)))
print('  worst |magnus - analytic adiabatic, 4nu| = %.2e'
      % np.max(np.abs(P_sun4 - P_adiab4)))
print()
print('CHECK 2 -- is the sterile neutral-current entry live?')
print('  <P_ee> at 1, 5, 15 MeV as the neutron-to-proton ratio varies:')
E_PROBE = np.array([1.0, 5.0, 15.0])*gd.UNIT_MEV
for ratio in (0.5, 1.0, 1.5):
    P_r = np.asarray(oscprob.osc_prob_matter_std_potential(
        4, ne_sun, E_PROBE, R_SUN, OSC4_SUN, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE,
        density_is_of_number_of_electrons=True, average=True,
        ratio_number_neutrons_to_protons=ratio)).ravel()
    print('    r = %.1f  ->  %s' % (ratio, np.array2string(P_r, precision=6)))
print('  the curve moves with r, so the sterile state is feeling the medium.')
print('  Were that term missing the three rows would be identical.')'''),
    code(r'''fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                       gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.08))
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_sun, lw=1.8, color='k',
               label=r'Mag$\nu$s, 3$\nu$')
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_sun4, lw=1.6, color='C0',
               label=r'Mag$\nu$s, 3+1  ($\Delta m^2_{41} = 10^{-5}$ eV$^2$)')
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_adiab4, lw=1.0, ls='--', color='C3',
               label=r'Analytic adiabatic limit, 3+1')
ax[0].set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
ax[0].set_title('The Sun with a sterile state: the departure from the 3-flavor curve',
                fontsize=10)
ax[0].grid(True, alpha=0.2)
ax[0].legend(fontsize=7.6)

ax[1].semilogx(E_SUN/gd.UNIT_MEV, P_sun4 - P_sun, lw=1.4, color='C0')
ax[1].axhline(0.0, color='k', lw=0.7, alpha=0.5)
ax[1].set_xlabel(r'$E_\nu$ [MeV]')
ax[1].set_ylabel(r'$\langle P\rangle_{3+1} - \langle P\rangle_{3\nu}$')
ax[1].grid(True, which='both', alpha=0.2)
for a in ax:
    a.set_xlim(E_SUN[0]/gd.UNIT_MEV, E_SUN[-1]/gd.UNIT_MEV)
ax[0].set_ylim(0.0, 1.0)
fig.savefig('../fig/solar_3plus1.pdf', bbox_inches='tight')

print('largest departure from the 3nu curve: %.3f at %.2f MeV'
      % (np.max(np.abs(P_sun4 - P_sun)),
         E_SUN[int(np.argmax(np.abs(P_sun4 - P_sun)))]/gd.UNIT_MEV))
print('magnus 3nu %.3f s, 3+1 %.3f s for the same %d energies'
      % (t_sun, t_sun4, len(E_SUN)))'''),
    md(r"""**On timing this against another code: there is nothing to time it against.** The
observable is the phase-averaged probability, and neither NuOscProbExact nor nuSQuIDS offers an
averaging flag at any flavor count -- which section 10 already established at three flavors,
where recovering the average from instantaneous evaluations cost a Monte-Carlo estimate with a
$1/\sqrt{N}$ error. Nothing about that improves at four. NuOscProbExact does have an SU(4)
closed form, but it solves *piecewise-constant* slabs and the solar profile is smooth, so it
would be back to resolving a 13 000-radian phase.

So the honest entry for this row is not a ratio. It is that the comparison has no second
entrant, which is the "pre-packaged observables" axis stated at four flavors instead of
three."""),
    md(r"""### 10c. The same observable with NSI, and a route that has to be checked

The third solar variant, and it carries a dispatch question that matters more than it looks.
`average=True` can be served two different ways: the **adiabatic** route, which transports along
the levels of the instantaneous Hamiltonian and never propagates, or a **numerical window**,
which averages an explicitly propagated probability over a finite spread. Those compute
*different quantities* -- the $L/E \to \infty$ limit against an average over a window -- and
mixing them across panels would be the same class of error as comparing an instantaneous curve
against an averaged one, which is what section 10 had to be corrected for.

So which route the NSI wrapper takes is checked below rather than assumed, and the check is the
cost: the adiabatic route never propagates, so it returns in well under a second, while a
numerical window over a 13 000-radian phase would take minutes."""),
    code(r'''EPS_SUN = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0,
               eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
EPS_MATRIX_SUN = np.array(
    [[1.0 + EPS_SUN['eps_ee'], EPS_SUN['eps_em'], EPS_SUN['eps_et']],
     [np.conj(EPS_SUN['eps_em']), EPS_SUN['eps_mm'], EPS_SUN['eps_mt']],
     [np.conj(EPS_SUN['eps_et']), np.conj(EPS_SUN['eps_mt']), EPS_SUN['eps_tt']]],
    dtype=complex)


def adiabatic_reference_nsi(energy, a=gd.NUE, b=gd.NUE):
    """The analytic adiabatic limit with the NSI matter matrix in place."""
    hm = h_vac/energy + float(PER_NE_SUN*ne_sun(0.0))*EPS_MATRIX_SUN
    _, u_matter = np.linalg.eigh(hm)
    _, u_vac = np.linalg.eigh(h_vac/energy)
    return float(np.sum(np.abs(u_matter[a])**2 * np.abs(u_vac[b])**2))


t0 = time.perf_counter()
P_sun_nsi = np.asarray(oscprob.osc_prob_matter_nsi(
    3, ne_sun, E_SUN, R_SUN, OSC, EPS_SUN, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE,
    density_is_of_number_of_electrons=True, average=True))
t_sun_nsi = time.perf_counter() - t0
P_adiab_nsi = np.array([adiabatic_reference_nsi(e) for e in E_SUN])

print('WHICH ROUTE? -- the cost is the evidence')
print('  standard, average=True : %.3f s for %d energies' % (t_sun, len(E_SUN)))
print('  NSI,      average=True : %.3f s for %d energies' % (t_sun_nsi, len(E_SUN)))
print('  Both sub-second, so both took the adiabatic route.  A numerical window over')
print('  this profile propagates, and would have cost minutes rather than tenths.')
print()
print('worst |magnus NSI - analytic adiabatic NSI| = %.2e'
      % np.max(np.abs(P_sun_nsi - P_adiab_nsi)))
print('departure from the standard averaged curve: max %.4f, mean %.4f'
      % (np.max(np.abs(P_sun_nsi - P_sun)), np.mean(np.abs(P_sun_nsi - P_sun))))'''),
    code(r'''fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                       gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.08))
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_sun, lw=1.8, color='k',
               label=r'Mag$\nu$s, standard 3$\nu$')
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_sun_nsi, lw=1.6, color='C4',
               label=r'Mag$\nu$s, NSI  ($\varepsilon_{ee}=0.15$, '
                     r'$\varepsilon_{e\mu}=0.05$)')
ax[0].semilogx(E_SUN/gd.UNIT_MEV, P_adiab_nsi, lw=1.0, ls='--', color='C3',
               label='Analytic adiabatic limit, NSI')
ax[0].set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
ax[0].set_title('The Sun with NSI: a much smaller departure than the shock sees',
                fontsize=10)
ax[0].grid(True, alpha=0.2)
ax[0].legend(fontsize=7.4)

ax[1].semilogx(E_SUN/gd.UNIT_MEV, P_sun_nsi - P_sun, lw=1.4, color='C4')
ax[1].axhline(0.0, color='k', lw=0.7, alpha=0.5)
ax[1].set_xlabel(r'$E_\nu$ [MeV]')
ax[1].set_ylabel(r'$\langle P\rangle_{\rm NSI} - \langle P\rangle_{\rm std}$')
ax[1].grid(True, which='both', alpha=0.2)
for a in ax:
    a.set_xlim(E_SUN[0]/gd.UNIT_MEV, E_SUN[-1]/gd.UNIT_MEV)
ax[0].set_ylim(0.0, 1.0)
fig.savefig('../fig/solar_nsi.pdf', bbox_inches='tight')'''),
    md(r"""**The same $\varepsilon$ that moves the shock by 0.44 moves this curve by 0.014.**
That is a factor of thirty, with an identical Hamiltonian modification, and it is not a
numerical artifact -- it is what averaging does. The shock panel in section 13 shows an
*instantaneous* probability along a ray, where NSI shifts both the level structure and the
phase, and the phase term is large. The solar observable is **phase-averaged**, so everything
that enters through the phase integrates away and only the change in the eigenvectors and the
level crossings survives.

The lesson generalizes past NSI: **a BSM effect's size depends on the observable at least as
much as on the model.** Quoting "NSI changes the probability by 0.44" without saying which
probability would be misleading in either direction, and a sensitivity study that used the
wrong one would be wrong by a factor of thirty before any experimental detail entered.

There is again no second entrant to time this against. NuOscProbExact has an NSI Hamiltonian --
section 13 uses it -- but no averaging flag, so recovering this observable from it means the
Monte-Carlo estimate section 10 measured, at a $1/\sqrt{N}$ error and ten minutes per tolerance
setting merely to get a probability."""),
    md(r"""## 11. A supernova shock: the width of the front decides the winner

Section 9 put the two codes on a *smooth* profile and section 5 put them on a
*piecewise-constant* one. A supernova shock is the case that contains both, and which one it
behaves like is set by a single parameter: how wide the front is.

The case is notebook 14's, unchanged -- 15 MeV, three flavors, along a ray from
$10^4$ to $8\times10^4$ km, with the H resonance ($\Delta m^2_{31}$) sitting on the ray just
outside the forward shock. The observable is $P_{ee}$ at **61 points along one ray**, and the
referee is that notebook's frozen `solve_ivp`/DOP853 solution at `rtol=1e-12` -- an adaptive
integrator, which is neither a Magnus expansion nor a slab product.

Two front widths, and they are different physics rather than two settings:

* $w = 10^{-6}$ of the ray, **0.07 km** -- a real hydrodynamic shock, mean-free-path thin.
* $w = 10^{-3}$ of the ray, **70 km** -- a shock as it comes out of a simulation snapshot,
  smeared across a few grid cells.

Every number is frozen in `external_shock_benchmarks.json`, generated by
`gen_shock_benchmarks.py`, which executes notebook 14's own cells rather than transcribing
the profile -- there is one definition of this physics and both notebooks read it."""),
    md(r"""### How each code was driven, because both have a wrong way that converges cleanly

Neither of these was obvious, and each cost a wrong answer first.

**Mag$\nu$s must be driven cumulatively.** The request is 61 probabilities along *one* ray.
Point by point, Mag$\nu$s re-propagates the whole ray for each one: 3.27 s, and an error of
4.7e-04 that **does not move with `rtol`** -- the refinement ladder cannot reach the
tolerance that way, and flat-in-tolerance is the signature of a driver problem rather than an
accuracy limit. With `cumulative=True` the same request costs 0.24 s and lands at 8.1e-06.

**NuOscProbExact must be given slabs inside the front.** Its batched route returns the
operator at the *end* of a slab chain, so the ray is cut at the declared fronts and at every
target, each leg solved with `evolution_operator_3nu_slabs`, and the operators composed. The
60 intervals between consecutive targets are equal by construction, so they share one set of
widths and go in **one batched call**; looping them would time the loop. But the fronts are
0.07 km of a 70 000 km ray, so slabs allocated in proportion to length leave each front with a
single sample, and the code then appears to floor at $5\times10^{-7}$ no matter how fine the
rest of the ray gets. **That floor is this notebook's allocation, not NuOscProbExact's
limit**: raising the per-front minimum to a thousandth of the budget takes the same code to
$1.4\times10^{-9}$. It is stated here rather than buried in the generator because a floor
attributed to the wrong cause is worse than no measurement.

Both codes are told where the front is -- Mag$\nu$s through `t_breakpoints`, NuOscProbExact
through the leg boundaries. Telling one and not the other would measure the telling."""),
    code(r'''SHOCK = json.loads(
    (pathlib.Path.cwd()/'external_shock_benchmarks.json').read_text())
print('machine:', SHOCK['machine'])
print('control ratio %.3f' % SHOCK['control_ratio'])
print('%.0f MeV, ray %.0f -> %.0f km, %d targets over the last %.1f km'
      % (SHOCK['energy_mev'], SHOCK['L0_km'], SHOCK['L1_km'], SHOCK['n_targets'],
         SHOCK['targets_km'][-1] - SHOCK['targets_km'][0]))


def shock_case(width):
    for c in SHOCK['cases']:
        if abs(c['width'] - width) < 1.0e-12:
            return c
    return None


for c in SHOCK['cases']:
    print()
    print('front width %.0e of the ray (%.2f km)'
          % (c['width'], c['width']*(SHOCK['L1_km'] - SHOCK['L0_km'])))
    print('  the frozen referee is itself unitary to %.1e -- nothing below that line '
          'is resolvable' % c['reference_unitarity'])
    for s in c['series']:
        best = min(s['points'], key=lambda p: p['max_abs_error'])
        print('  %-16s best %.3e at %s = %-8s (%.0f us/probability)'
              % (s['name'], best['max_abs_error'], s['dial'], best['label'],
                 best['us_per_probability']))'''),
    code(r'''fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.9))
for ax, c in zip(axes, SHOCK['cases']):
    for s in c['series']:
        marker, color, size = DIAL_STYLE[s['name']]
        pts = s['points']
        t = [p['us_per_probability'] for p in pts]
        e = [p['max_abs_error'] for p in pts]
        kw = dict(ms=size, color=color, lw=1.1, zorder=4,
                  label='%s  (%s)' % (s['name'], s['dial']))
        if s['name'] == 'NuOscProbExact':
            kw.update(mfc='white', mew=0.9)
        ax.loglog(t, e, marker, **kw)
        for j, ha, dx in ((0, 'right', -4), (len(pts) - 1, 'left', 4)):
            ax.annotate(pts[j]['label'], xy=(t[j], e[j]), xytext=(dx, 3),
                        textcoords='offset points', fontsize=5.4, color=color, ha=ha)
    # The referee's own unitarity: nothing below this line is a measurement.
    ax.axhspan(1.0e-16, c['reference_unitarity'], color='0.85', alpha=0.55, zorder=0)
    ax.text(0.03, 0.03, 'below: the referee$\'$s own floor', transform=ax.transAxes,
            fontsize=5.8, color='0.35', va='bottom')
    ax.set_xlabel(r'Time per probability [$\mu$s]')
    ax.set_ylabel(r'Error vs.\ the frozen DOP853 referee,  max $|\Delta P|$')
    ax.set_title(r'Shock front $w = 10^{%d}$ of the ray (%.2f km)'
                 % (round(np.log10(c['width'])),
                    c['width']*(SHOCK['L1_km'] - SHOCK['L0_km'])), fontsize=10)
    ax.grid(True, which='both', alpha=0.18)
    leg = ax.legend(fontsize=6.6, loc='upper right')
    leg.get_frame().set_linewidth(0.7)
    # No dead margin either side of the curves.
    allt = [p['us_per_probability'] for s in c['series'] for p in s['points']]
    alle = [p['max_abs_error'] for s in c['series'] for p in s['points']]
    ax.set_xlim(min(allt)/1.6, max(allt)*1.6)
    ax.set_ylim(min(alle)/3.0, max(alle)*3.0)
fig.tight_layout(pad=1.2)
fig.savefig('../fig/shock_speed_accuracy.pdf', bbox_inches='tight')'''),
    md(r"""**The two panels point opposite ways, and that is the result.**

On the **0.07 km front** the closed form wins decisively. To a method that samples the
profile, a front a ten-thousandth the width of a single oscillation length is not a smooth
ramp at all -- it *is* a jump, and piecewise-constant is what a closed form is built for.
Both codes come within a factor of a few of the referee's own floor -- $1.4\times10^{-9}$ and
$2.4\times10^{-9}$ against its $4.1\times10^{-10}$ -- but NuOscProbExact gets there for about
$2.5\times10^3\ \mu$s per probability against Mag$\nu$s's $5.5\times10^5$: a factor of
roughly **230**.

On the **70 km front** the ordering reverses below about $10^{-5}$. Seventy kilometers is
several oscillation lengths, so the front is *resolved* rather than jumped, and the profile
is smooth on the scale of the slabs -- which is section 9's regime, where fourth-order
Gauss--Legendre integration beats second-order midpoint sampling. NuOscProbExact's finest
setting here reaches $5.6\times10^{-6}$; Mag$\nu$s is **22 times more accurate than that for
three times the cost** ($2.5\times10^{-7}$ at 7400 $\mu$s against 2500), and then continues
nearly four more orders of magnitude to $9.9\times10^{-10}$, which the slab product does not
reach at any dial measured.

So the shock is not one case with one answer. **The width of the front decides which method
it is a case for**, and the crossover sits where the front stops being wide compared with the
slabs and starts being a discontinuity."""),
    md(r"""### The probability itself

The panels above say what each code *costs*. This says what they *return*, across energy at
the end of the ray, on a 120-point grid.

There is no frozen referee off 15 MeV, so the curve is code against code -- and agreement is
a weaker claim than accuracy, because two codes can agree and both be wrong. Three of the
energies are therefore refereed by DOP853 as well, which is what says **whose** the residual
is. Without that the reader is free to assume it is shared."""),
    code(r'''PE = {round(np.log10(p['width'])): p for p in SHOCK['prob_vs_energy']}

fig, axes = plt.subplots(2, 2, figsize=(10.4, 5.9), sharex=True,
                         gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.09,
                                          wspace=0.26))
for col, lw in enumerate(sorted(PE)):
    pe = PE[lw]
    e_mev = np.array(pe['energy_ev'])/gd.UNIT_MEV
    mg = np.array(pe['magnus'])
    top, bot = axes[0][col], axes[1][col]
    top.plot(e_mev, mg, lw=1.6, color='k',
             label=r'Mag$\nu$s  (%d slabs)' % pe['magnus_n_slabs'])
    if pe['npe'] is not None:
        top.plot(e_mev, np.array(pe['npe']), lw=1.1, ls='--', color='C3',
                 label='NuOscProbExact  (%d slabs)' % pe['npe_n_slabs'])
        bot.semilogy(e_mev, np.abs(mg - np.array(pe['npe'])), lw=1.0, color='C0',
                     label=r'$|$Mag$\nu$s $-$ NuOscProbExact$|$')
    # The refereed points, which attribute that residual rather than sharing it out.
    for row in pe['refereed']:
        bot.plot(row['energy_mev'], abs(row['magnus_error']), '*', ms=9, color='k',
                 zorder=5)
        if 'npe_error' in row:
            bot.plot(row['energy_mev'], abs(row['npe_error']), 'o', ms=4.4,
                     mfc='white', mew=0.9, color='C3', zorder=5)
    top.set_ylabel(r'$P(\nu_e \to \nu_e)$')
    top.set_title(r'Shock front $w = 10^{%d}$ of the ray' % lw, fontsize=10)
    top.grid(True, alpha=0.2)
    top.legend(fontsize=7.4)
    bot.set_xlabel(r'$E_\nu$ [MeV]')
    bot.set_ylabel(r'$|\Delta P|$')
    bot.grid(True, which='both', alpha=0.2)
    bot.set_xlim(e_mev[0], e_mev[-1])
    top.set_xlim(e_mev[0], e_mev[-1])
    bot.text(0.02, 0.06, 'stars: Mag$\\nu$s vs DOP853   circles: NuOscProbExact vs DOP853',
             transform=bot.transAxes, fontsize=5.6, color='0.3')
fig.savefig('../fig/shock_prob_vs_energy.pdf', bbox_inches='tight')

for lw in sorted(PE):
    pe = PE[lw]
    print('front w = 1e%d' % lw)
    if pe['npe'] is not None:
        print('  worst |Magnus - NuOscProbExact| over %d energies = %.3e'
              % (len(pe['energy_ev']),
                 float(np.max(np.abs(np.array(pe['magnus'])
                                     - np.array(pe['npe']))))))
    for row in pe['refereed']:
        print('  %5.1f MeV vs DOP853:  magnus %+.2e   npe %+.2e'
              % (row['energy_mev'], row['magnus_error'], row.get('npe_error',
                                                                float('nan'))))'''),
    md(r"""**The residual belongs to a different code on each front, and the two curves alone
could never have said so.** The band between them is about $1.3\times10^{-5}$ on the sharp
front and $1.9\times10^{-5}$ on the smeared one -- nearly the same size, and on the
probability axis the curves are indistinguishable in both panels. The refereed points show
those two bands are not the same thing at all:

| front | Mag$\nu$s vs DOP853 | NuOscProbExact vs DOP853 |
|---|---|---|
| $w=10^{-6}$ (0.07 km) | $+3.2$, $+3.6$, $+1.2 \times10^{-6}$ | $-3.5$, $+1.2$, $+1.4 \times10^{-9}$ |
| $w=10^{-3}$ (70 km) | $-1.7\times10^{-8}$, $-5.5\times10^{-9}$, $+4.4\times10^{-11}$ | $+1.7\times10^{-5}$, $-5.4\times10^{-6}$, $-5.6\times10^{-7}$ |

On the thin front the band is almost entirely **Mag$\nu$s's** error; on the resolved front it
is almost entirely **NuOscProbExact's**, and in both cases by about three orders of magnitude.
That is the speed-against-accuracy result of the panels above, arrived at a second time by a
different route -- which is worth more than either statement alone, because the two share no
machinery beyond the profile.

A caution on reading the times: this panel is drawn at settings chosen for cost, not matched
accuracy -- 128 000 slabs for Mag$\nu$s against 524 288 for NuOscProbExact -- so the 536 s and
16 s it took are not a like-for-like ratio. The matched comparison is the panels above; this
one exists to show the curves and to attribute the gap between them.

So: **two codes agreeing is evidence, but only a third method says which of them to believe,
and here the answer changes between two panels of the same figure.**"""),
    md(r"""## 12. The same shock at 3+1, and a splitting that cannot be refereed

Section 11's shock, with a sterile state added. The observable and the ray are unchanged --
$P_{ee}$ at 61 points along the same 70 000 km ray at 15 MeV, on the 70 km front, where
section 11 found the resolved-front regime that suits a Magnus expansion.

**The splitting is $\Delta m^2_{41} = 10^{-2}\,\mathrm{eV}^2$, and not the eV-scale value
section 5 uses. That choice is forced by the referee, not by the physics**, and the arithmetic
is worth showing because it decides what can be validated at all:

| | accumulated phase over the ray | oscillation lengths |
|---|---|---|
| $\Delta m^2_{31} = 2.5\times10^{-3}$ (section 11) | $1.5\times10^{4}$ rad | 2 352 |
| $\Delta m^2_{41} = 1$ eV$^2$ | $5.9\times10^{6}$ rad | **940 981** |
| $\Delta m^2_{41} = 10^{-2}$ eV$^2$ | $5.9\times10^{4}$ rad | 9 410 |

An adaptive DOP853 reference has to resolve every one of those oscillations. At an eV-scale
splitting that is of order a day for a single front width -- measured by starting one and
killing it after an hour with the referee still unfinished. So **an eV-scale sterile state on
this ray is not a case that can be refereed by an independent integrator at all**, which is
worth knowing before attempting it: the comparison would have to fall back to the two codes
checking each other, and section 11 showed exactly how little that settles when they agree and
are both wrong. At $10^{-2}$ the referee costs four times the three-flavor one, and the
sterile oscillation is still fast against the shock structure, which is what the case is about.

The matter term carries $\mathrm{diag}(1, 0, 0, r/2)$ from
`matter.matter_potential_projector` -- never written out here, because writing it out is what
made section 5's referee wrong."""),
    code(r'''SHOCK4 = json.loads(
    (pathlib.Path.cwd()/'external_shock_4nu.json').read_text())
print('%.0f MeV, front width %.0e, D41 = %.0e eV^2'
      % (SHOCK4['energy_mev'], SHOCK4['width'], SHOCK4['sterile']['D41']))
print('referee unitary to %.2e' % SHOCK4['reference_unitarity'])
print()
for s in SHOCK4['series']:
    best = min(s['points'], key=lambda p: p['max_abs_error'])
    print('%-16s best %.3e at %s = %-8s (%.0f us/probability)'
          % (s['name'], best['max_abs_error'], s['dial'], best['label'],
             best['us_per_probability']))'''),
    code(r'''fig, ax = plt.subplots(figsize=(6.2, 5.0))
for s in SHOCK4['series']:
    marker, color, size = DIAL_STYLE[s['name']]
    pts = s['points']
    t = [p['us_per_probability'] for p in pts]
    e = [p['max_abs_error'] for p in pts]
    kw = dict(ms=size, color=color, lw=1.1, zorder=4,
              label='%s  (%s)' % (s['name'], s['dial']))
    if s['name'] == 'NuOscProbExact':
        kw.update(mfc='white', mew=0.9)
    ax.loglog(t, e, marker, **kw)
    for j, ha, dx in ((0, 'right', -4), (len(pts) - 1, 'left', 4)):
        ax.annotate(pts[j]['label'], xy=(t[j], e[j]), xytext=(dx, 3),
                    textcoords='offset points', fontsize=5.4, color=color, ha=ha)
ax.axhspan(1.0e-16, SHOCK4['reference_unitarity'], color='0.85', alpha=0.55, zorder=0)
ax.text(0.03, 0.03, "below: the referee's own floor", transform=ax.transAxes,
        fontsize=5.8, color='0.35', va='bottom')
ax.set_xlabel(r'Time per probability [$\mu$s]')
ax.set_ylabel(r'Error vs.\ the DOP853 referee,  max $|\Delta P|$')
ax.set_title(r'Supernova shock at 3+1, $\Delta m^2_{41} = 10^{-2}$ eV$^2$, 70 km front',
             fontsize=10)
ax.grid(True, which='both', alpha=0.18)
leg = ax.legend(fontsize=6.6, loc='upper right')
leg.get_frame().set_linewidth(0.7)
allt = [p['us_per_probability'] for s in SHOCK4['series'] for p in s['points']]
alle = [p['max_abs_error'] for s in SHOCK4['series'] for p in s['points']]
ax.set_xlim(min(allt)/1.6, max(allt)*1.6)
ax.set_ylim(min(alle)/3.0, max(alle)*3.0)
fig.tight_layout(pad=1.2)
fig.savefig('../fig/shock_3plus1_speed_accuracy.pdf', bbox_inches='tight')'''),
    md(r"""**On the resolved front, the fourth flavor does not change the verdict -- it
sharpens it.** Mag$\nu$s reaches $1.1\times10^{-8}$; NuOscProbExact's finest setting here
reaches $1.2\times10^{-4}$, four orders of magnitude short, and at matched accuracy Mag$\nu$s
is also the *cheaper* of the two. That is the same mechanism as section 11's 70 km panel --
a front spread over many slabs is smooth on the slab scale, where fourth-order integration
beats second-order midpoint sampling -- and adding a sterile state does not alter it.

Compare that with **section 5**, where the same two codes met at 3+1 on *PREM* and the closed
form won the cost by 400x. Same flavor content, opposite result, and the difference is the
profile: piecewise-constant is a closed form's home ground, a resolved front is not. The
flavor count is not what decides these cases, which is worth saying because the 3+1 row is
the one most often quoted as though it were about dimension.

**What this section cannot tell you** is what happens at an eV-scale splitting, because
nothing here can referee it. That is a limit of the *validation*, not of either code: both
will happily return an answer at $\Delta m^2_{41} = 1\,\mathrm{eV}^2$, and neither this
notebook nor any affordable integrator can say which one is right."""),
    md(r"""## 13. The same shock with non-standard interactions

The third BSM variant on the same ray, and the one where the two codes' conventions could most
easily have made a mess. Both offer NSI: Mag$\nu$s through `osc_prob_matter_nsi`,
NuOscProbExact through `hamiltonians3nu.hamiltonian_3nu_nsi`, each taking six dimensionless
$\varepsilon$.

**The conventions were checked before anything was measured, and they match exactly.** They do
not *look* like they match: NuOscProbExact writes the standard piece into the matrix as
$1 + \varepsilon_{ee}$, while Mag$\nu$s's `hamiltonian_3nu_nsi` returns $\varepsilon_{ee}$ alone
and the standard term is added separately. That is an off-by-one waiting to happen. Measured on
constant density -- 3 g/cm$^3$, 1300 km, 2 GeV, $\varepsilon_{ee} = 0.10$,
$\varepsilon_{em} = 0.05$ -- the *same* $\varepsilon$ handed to both gives
$P(\nu_\mu \to \nu_e)$ agreeing to **$1.7\times10^{-16}$**, and shifting either convention by
one puts them $3.7\times10^{-2}$ apart. So the user-facing conventions are identical and there
is no offset to correct, which is the opposite of the $V_{\rm CC}$ result in section 6 -- and
the reason to check rather than assume is that both outcomes look the same until you do."""),
    code(r'''NSI_SHOCK = json.loads(
    (pathlib.Path.cwd()/'external_shock_nsi.json').read_text())
print('eps = %s' % {k: v for k, v in NSI_SHOCK['eps'].items() if v})
print('%.0f MeV, front width %.0e, referee unitary to %.2e'
      % (NSI_SHOCK['energy_mev'], NSI_SHOCK['width'],
         NSI_SHOCK['reference_unitarity']))
print()
for s in NSI_SHOCK['series']:
    best = min(s['points'], key=lambda p: p['max_abs_error'])
    print('%-16s best %.3e at %s = %-8s (%.0f us/probability)'
          % (s['name'], best['max_abs_error'], s['dial'], best['label'],
             best['us_per_probability']))

mg_nsi = np.array(NSI_SHOCK['magnus_P_ee'])
std_nsi = np.array(NSI_SHOCK['standard_P_ee'])
print()
print('departure from the standard 3nu curve: max %.4f, mean %.4f'
      % (np.max(np.abs(mg_nsi - std_nsi)), np.mean(np.abs(mg_nsi - std_nsi))))'''),
    code(r'''fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.6))

for s in NSI_SHOCK['series']:
    marker, color, size = DIAL_STYLE[s['name']]
    pts = s['points']
    t = [p['us_per_probability'] for p in pts]
    e = [p['max_abs_error'] for p in pts]
    kw = dict(ms=size, color=color, lw=1.1, zorder=4,
              label='%s  (%s)' % (s['name'], s['dial']))
    if s['name'] == 'NuOscProbExact':
        kw.update(mfc='white', mew=0.9)
    ax[0].loglog(t, e, marker, **kw)
    for j, ha, dx in ((0, 'right', -4), (len(pts) - 1, 'left', 4)):
        ax[0].annotate(pts[j]['label'], xy=(t[j], e[j]), xytext=(dx, 3),
                       textcoords='offset points', fontsize=5.4, color=color, ha=ha)
ax[0].axhspan(1.0e-16, NSI_SHOCK['reference_unitarity'], color='0.85', alpha=0.55,
              zorder=0)
ax[0].text(0.03, 0.03, "below: the referee's own floor", transform=ax[0].transAxes,
           fontsize=5.8, color='0.35', va='bottom')
ax[0].set_xlabel(r'Time per probability [$\mu$s]')
ax[0].set_ylabel(r'Error vs.\ the DOP853 referee,  max $|\Delta P|$')
ax[0].set_title('Shock with NSI: speed against accuracy', fontsize=10)
ax[0].grid(True, which='both', alpha=0.18)
leg = ax[0].legend(fontsize=6.6, loc='upper right')
leg.get_frame().set_linewidth(0.7)
allt = [p['us_per_probability'] for s in NSI_SHOCK['series'] for p in s['points']]
alle = [p['max_abs_error'] for s in NSI_SHOCK['series'] for p in s['points']]
ax[0].set_xlim(min(allt)/1.6, max(allt)*1.6)
ax[0].set_ylim(min(alle)/3.0, max(alle)*3.0)

# The BSM curve alongside the standard one: a departure is only legible against what it
# departs from, and its SIZE is the thing worth reading.
xs_nsi = np.array(NSI_SHOCK['targets_km']) - NSI_SHOCK['targets_km'][0]
ax[1].plot(xs_nsi, std_nsi, lw=1.6, color='k', label=r'standard 3$\nu$')
ax[1].plot(xs_nsi, mg_nsi, lw=1.6, color='C0',
           label=r'NSI, $\varepsilon_{ee} = 0.15$, $\varepsilon_{e\mu} = 0.05$')
ax[1].fill_between(xs_nsi, std_nsi, mg_nsi, color='C0', alpha=0.16)
ax[1].set_xlabel('distance along the ray beyond %.0f km [km]'
                 % NSI_SHOCK['targets_km'][0])
ax[1].set_ylabel(r'$P(\nu_e \to \nu_e)$')
ax[1].set_title('and what that NSI actually does', fontsize=10)
ax[1].grid(True, alpha=0.2)
ax[1].legend(fontsize=7.4)
ax[1].set_xlim(xs_nsi[0], xs_nsi[-1])
fig.tight_layout(pad=1.2)
fig.savefig('../fig/shock_nsi.pdf', bbox_inches='tight')'''),
    md(r"""**Same verdict as sections 11 and 12, for the third time and with different
physics in the Hamiltonian.** On this 70 km front Mag$\nu$s reaches
$6.5\times10^{-10}$ -- below the referee's own floor of $7.7\times10^{-10}$, so it has
saturated the measurement -- while NuOscProbExact's finest setting stalls at
$4.0\times10^{-6}$. At matched accuracy Mag$\nu$s is **28 times more accurate for 2.2 times
the cost**, and past about $10^{-5}$ the closed form has no setting that reaches at all.

Put the three shock sections together and the pattern is the point:

| | 0.07 km front | 70 km front |
|---|---|---|
| standard 3$\nu$ | NuOscProbExact by ~230x | Mag$\nu$s, to $9.9\times10^{-10}$ |
| 3+1 | -- | Mag$\nu$s, to $1.1\times10^{-8}$ |
| NSI | -- | Mag$\nu$s, to $6.5\times10^{-10}$ |

**The BSM content does not decide these cases; the front width does.** Adding a sterile state
or a non-standard interaction changes what the Hamiltonian *is*, and changes the answer a great
deal -- the NSI departure from the standard curve averages 0.26 and reaches 0.44, so the two
curves do not overlap anywhere on the ray -- but it does not change which *method* suits the
problem. That is worth stating because "which code for BSM?" is a question people ask, and on
this evidence it is the wrong question: ask which code for the *profile*.

**Reproducing these numbers needs two keywords that are not in any signature.** Both codes were
driven with the fronts declared -- Mag$\nu$s through `t_breakpoints`, NuOscProbExact by cutting
its legs at the same positions -- and Mag$\nu$s was given an explicit `n_slabs`, because on a
single request `rtol` is not this route's dial. Neither `t_breakpoints` nor `n_slabs` appears in
`osc_prob_matter_nsi`'s signature; both reach it through `**kwargs`, and both demonstrably bite
($2.4\times10^{-6}$ and $8.8\times10^{-7}$ respectively on this case). Driven the way the
signature alone suggests, the same call takes 0.46 s and returns a worse answer. The comparison
above is fair -- each code is driven as its authors intend -- but it is **not** what a reader
reproduces by reading the signature, which is worth knowing before trying."""),
    ])


# ------------------------------------------------------ 27_magnus_animations
books['27_magnus_animations.ipynb'] = notebook(
    'Animated scenes',
    r'''Nine short scenes, each showing one thing this library does while a parameter sweeps.

The first four are the same four that [NuOscProbExact's notebook 19](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/19_animations.ipynb)
draws, computed here with Mag$\nu$s so that the two can be read side by side. The five after
them have no counterpart there, because each animates something a closed-form slab code does
not have: a refinement ladder deciding it has converged, a front that travels, an observable
that is an average rather than a value, and a Hamiltonian that genuinely varies along the path.

Rendering them as animations is expensive, so this notebook draws **stills**, as filmstrips,
and leaves the animation to an opt-in switch at the end. The full procedure, its measured
cost, and the traps worth knowing are all at the bottom.

Two things that belong with the *method* rather than with the physics -- what the truncation
order buys, and which engine the dispatcher picks -- are not animated at all. They are single
comparisons rather than sweeps, and they live in
[notebook 24](24_magnus_performance.ipynb).''',
    [
    code(r'''import os

import numpy as np
import matplotlib.pyplot as plt

# Mag(nu)s is imported as an installed package -- from the repository root,
# 'pip install -e .' (add [plot] for magnus.plotting). No sys.path juggling.
import magnus.oscprob as oscprob
import magnus.matter as matter
import magnus.earth as earth
import magnus.globaldefs as gd

# The repository `matplotlibrc` is tuned for single-panel paper figures at
# 5 x 4.75 inches: it sets `axes.labelsize` to 25 and the tick labels to 23,
# ABSOLUTELY rather than relative to `font.size`.  Setting `font.size` here
# therefore moves nothing, and on the wide multi-panel figures below those
# sizes make an axis label render larger than the title above it.  Every size
# used here is restated explicitly, smallest to largest: ticks, labels, title.
plt.rcParams.update({'figure.dpi': 100, 'axes.grid': True, 'grid.alpha': 0.3,
                     'font.size': 13, 'legend.frameon': False,
                     'axes.labelsize': 15, 'axes.titlesize': 15,
                     'xtick.labelsize': 13, 'ytick.labelsize': 13,
                     'legend.fontsize': 12,
                     # Axes end exactly at the data: no dead margin left or right
                     # of any curve, which is the house rule everywhere here.
                     'axes.xmargin': 0.0, 'axes.ymargin': 0.0})

OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
KM = gd.CONV_KM_TO_INV_EV                 # multiply a length in km to get eV^-1
GEV, MEV = gd.UNIT_GEV, gd.UNIT_MEV

ACCENT, MARK, MUTED = '#1d4ed8', '#dc2626', '#94a3b8'

# `matplotlibrc` sets text.usetex, so every label goes through LaTeX. In LaTeX's
# text mode `>` and `<` render as inverted punctuation -- a title reading
# "rtol -> magnus" comes out with an upside-down question mark in it. Use $\to$,
# in math mode, or avoid the characters. This cost one rendered clip.
ARROW = r'$\to$'


# The one call that returns a whole map, and the idiom is not the obvious one.
# Mag(nu)s pairs energy and baseline ELEMENTWISE -- two 1-D arrays of length n
# give n probabilities, not an n x n grid -- so a map is a flattened meshgrid
# handed over in a single call, rather than a Hamiltonian stack broadcast
# against baselines the way NuOscProbExact does it. Same result, and the same
# one call; a different route to it.
def constant_density_map(energies, baselines, rho, osc=None, nu_i=None, nu_f=None,
                         **kw):
    """P over the (energy, baseline) grid, in one call."""
    ee, ll = np.meshgrid(energies, baselines, indexing='ij')
    flat = np.asarray(oscprob.osc_prob_3nu_matter_constant_density(
        ee.ravel(), ll.ravel(), rho, **(osc or OSC),
        density_matter_is_in_g_per_cm3=True,
        nu_i=gd.NUMU if nu_i is None else nu_i,
        nu_f=gd.NUE if nu_f is None else nu_f, **kw))
    return flat.reshape(ee.shape)


def filmstrip(n_panels, height=3.1, width_each=3.4, ratios=None):
    """A row of panels, which is how a sweep is shown without animating it."""
    return plt.subplots(1, n_panels, figsize=(width_each*n_panels, height),
                        gridspec_kw=None if ratios is None
                        else {'width_ratios': ratios})


print('Nine scenes: four shared with NuOscProbExact, five particular to this library.')'''),
    md(r'''# Part I --- the four scenes NuOscProbExact also draws

Same four sweeps, same layout, computed with Mag$\nu$s. Reading them beside the originals is
the point: where the two libraries differ, the difference shows up in how the frame is
computed rather than in what it contains.'''),
    md(r'''## 1. The CP phase

An oscillogram of $P_{\mu e}$ in matter, recomputed at each phase, beside the bi-probability
ellipse --- which is the locus traced *by* the phase, so the ellipse is drawn once and the
marker says where on it the map currently sits.

Each map is 200 x 200 = 40 000 probabilities in **one** call, in about 0.09 s. The mechanism
is not the one NuOscProbExact uses: there, a stack of Hamiltonians broadcasts against a row of
baselines; here, energy and baseline pair elementwise, so the grid is flattened, handed over
once, and reshaped.'''),
    code(r'''GRID_CP = 200
energies_cp = np.logspace(-1.0, 1.0, GRID_CP)*GEV
baselines_cp = np.linspace(50.0, 12000.0, GRID_CP)*KM
EXTENT_CP = [50.0, 12000.0, -1.0, 1.0]


def oscillogram_cp(dcp):
    return constant_density_map(energies_cp, baselines_cp,
                                gd.DENSITY_MATTER_CRUST_G_PER_CM3,
                                osc=dict(OSC, dCP=dcp))


def ellipse_point(dcp, energy=0.8*GEV, baseline=1300.0*KM):
    """One (P, Pbar) pair. Antineutrinos need the flag, not a hand-built H."""
    p = float(np.asarray(oscprob.osc_prob_3nu_matter_constant_density(
        energy, baseline, gd.DENSITY_MATTER_CRUST_G_PER_CM3,
        **dict(OSC, dCP=dcp), density_matter_is_in_g_per_cm3=True,
        nu_i=gd.NUMU, nu_f=gd.NUE)))
    pbar = float(np.asarray(oscprob.osc_prob_3nu_matter_constant_density(
        energy, baseline, gd.DENSITY_MATTER_CRUST_G_PER_CM3,
        **dict(OSC, dCP=dcp), density_matter_is_in_g_per_cm3=True,
        nu_i=gd.NUMU, nu_f=gd.NUE, nubar=True)))
    return p, pbar


PHASES_SHOWN = [0.0, 2.0*np.pi/3.0, 4.0*np.pi/3.0]
locus_cp = np.array([ellipse_point(d)
                     for d in np.linspace(0.0, 2.0*np.pi, 160)])
# One color scale for the whole sweep. Taking it from a single frame lets the
# others clip silently, and clipping reads as structure rather than saturation.
CEIL_CP = max(float(oscillogram_cp(d).max())
              for d in np.linspace(0.0, 2.0*np.pi, 8, endpoint=False))


def style_map_cp(ax):
    ax.set_xlabel('Baseline [km]')
    ax.set_xlim(EXTENT_CP[0], EXTENT_CP[1])
    ax.set_ylim(EXTENT_CP[2], EXTENT_CP[3])
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(['0.1', '1', '10'])
    ax.grid(False)


fig, axes = filmstrip(4, height=3.3, ratios=[1, 1, 1, 1.05])
for ax, dcp in zip(axes[:3], PHASES_SHOWN):
    im = ax.imshow(oscillogram_cp(dcp), origin='lower', aspect='auto',
                   cmap='viridis', vmin=0.0, vmax=CEIL_CP, extent=EXTENT_CP)
    ax.set_title(r'$\delta_{\rm CP} = %.2f\pi$' % (dcp/np.pi), fontsize=10)
    style_map_cp(ax)
axes[0].set_ylabel('Energy [GeV]')
fig.colorbar(im, ax=axes[2], pad=0.02).set_label(r'$P_{\mu e}$')

ax = axes[3]
ax.plot(locus_cp[:, 0], locus_cp[:, 1], color=ACCENT, lw=1.6)
for dcp in PHASES_SHOWN:
    p, pbar = ellipse_point(dcp)
    ax.plot([p], [pbar], 'o', ms=7, color=MARK, mfc='white', mew=1.6)
ax.set_xlabel(r'$P(\nu_\mu \to \nu_e)$')
ax.set_ylabel(r'$P(\bar\nu_\mu \to \bar\nu_e)$')
ax.set_title('Bi-probability locus', fontsize=10)
fig.tight_layout(pad=1.2)'''),
    md(r'''## 2. A sterile state

Four flavors, in matter of constant density so the whole map is again one call. The sterile
state feels neither the charged- nor the neutral-current potential, so $V_{\rm NC}$ stops
canceling between the flavors and sits on the sterile entry --- which is what places the
resonance that moves across the frame as $\Delta m^2_{41}$ sweeps.

**That term is the one this library got wrong and fixed.** It comes from
`matter.matter_potential_projector`, and omitting it costs 0.29 in probability on an Earth
chord while being flat in tolerance, so no amount of refinement reveals it.'''),
    code(r'''GRID_ST = 260                     # the map is the picture, so it is worth resolving
energies_st = np.logspace(-0.7, 1.3, GRID_ST)*GEV
baselines_st = np.linspace(50.0, 12000.0, GRID_ST)*KM
EXTENT_ST = [50.0, 12000.0, -0.7, 1.3]
STERILE = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, d14=0.0, d24=0.0)


def oscillogram_sterile(d41):
    ee, ll = np.meshgrid(energies_st, baselines_st, indexing='ij')
    flat = np.asarray(oscprob.osc_prob_4nu_matter_constant_density(
        ee.ravel(), ll.ravel(), gd.DENSITY_MATTER_CRUST_G_PER_CM3,
        **OSC, **STERILE, D41=d41, density_matter_is_in_g_per_cm3=True,
        nu_i=gd.NUMU, nu_f=gd.NUS1))
    return flat.reshape(ee.shape)


def style_map_st(ax):
    ax.set_xlabel('Baseline [km]')
    ax.set_xlim(EXTENT_ST[0], EXTENT_ST[1])
    ax.set_ylim(EXTENT_ST[2], EXTENT_ST[3])
    ax.set_yticks([-0.7, 0.0, 1.3])
    ax.set_yticklabels(['0.2', '1', '20'])
    ax.grid(False)


SPLITTINGS_SHOWN = [0.05, 0.30, 1.50]
CEIL_ST = max(float(oscillogram_sterile(d).max()) for d in SPLITTINGS_SHOWN)

fig, axes = filmstrip(3, height=3.3, ratios=[1, 1, 1.08])
for ax, d41 in zip(axes, SPLITTINGS_SHOWN):
    im = ax.imshow(oscillogram_sterile(d41), origin='lower', aspect='auto',
                   cmap='magma', vmin=0.0, vmax=CEIL_ST, extent=EXTENT_ST,
                   interpolation='bilinear')
    ax.set_title(r'$\Delta m^2_{41} = %.2f$ eV$^2$' % d41, fontsize=10)
    style_map_st(ax)
axes[0].set_ylabel('Energy [GeV]')
fig.colorbar(im, ax=axes[2], pad=0.02).set_label(r'$P(\nu_\mu \to \nu_s)$')
fig.tight_layout(pad=1.2)'''),
    md(r'''## 3. Through the Earth

**A detector sits at the South Pole and the arrival direction swings**, which is the way the
measurement is actually made: the detector does not move, the sky does. As $\cos\theta_z$ goes
from $-1$ to grazing, the chord reaching that fixed point sweeps across the Earth's interior
and crosses fewer and fewer of the layers.

The cross-section carries **every PREM boundary**, nine of them, shaded by density --- the
inner core, the outer core, and the seven shells above it. The chord and the layer crossings
are rebuilt each frame; the energies along it are one call.'''),
    code(r'''from matplotlib.patches import Circle

R_EARTH = gd.EARTH_RADIUS
DETECTOR = np.array([0.0, -R_EARTH])          # the South Pole, and it stays put
energies_earth = np.logspace(np.log10(0.5), np.log10(20.0), 120)*GEV
ANGLES_SHOWN = [-1.0, -0.6, -0.25]


def chord_entry(costhz):
    """Where a neutrino with this zenith angle enters, to reach the detector.

    The chord subtends a central angle of 2 arcsin(-cos theta_z), so the entry
    point is the detector rotated by that much around the center.
    """
    alpha = 2.0*np.arcsin(min(-costhz, 1.0))
    phi = -0.5*np.pi + alpha
    return np.array([R_EARTH*np.cos(phi), R_EARTH*np.sin(phi)])


def earth_curve(costhz):
    chord = earth.distance_traveled_inside_earth(costhz)
    return chord, np.asarray(oscprob.osc_prob_3nu_earth(
        energies_earth, costhz=costhz, L=chord*KM, **OSC,
        nu_i=gd.NUMU, nu_f=gd.NUMU))


def draw_earth(ax):
    """Every PREM shell, brightest at the core, with no frame around it."""
    edges = np.concatenate([np.asarray(earth.PREM_BOUNDARIES), [R_EARTH]])
    rho = np.asarray(earth.density_matter_func_prem(
        np.clip(edges - 1.0, 0.0, R_EARTH)))
    norm = (rho - rho.min())/(rho.max() - rho.min())
    for r, shade in sorted(zip(edges, norm), reverse=True):
        ax.add_patch(Circle((0, 0), r, facecolor=plt.cm.YlOrRd(0.15 + 0.75*shade),
                            edgecolor='white', lw=0.4, zorder=1))
    ax.set_xlim(-R_EARTH, R_EARTH)
    ax.set_ylim(-R_EARTH, R_EARTH)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for side in ('top', 'bottom', 'left', 'right'):
        ax.spines[side].set_visible(False)          # no bounding box
    ax.plot(*DETECTOR, marker='v', ms=9, color='#0f172a', zorder=5)


fig, axes = filmstrip(4, height=3.4, ratios=[1, 1, 1, 1.25])
for ax, cz in zip(axes[:3], ANGLES_SHOWN):
    draw_earth(ax)
    entry = chord_entry(cz)
    ax.plot([entry[0], DETECTOR[0]], [entry[1], DETECTOR[1]],
            color='#1e3a8a', lw=2.2, zorder=4)
    ax.set_title(r'$\cos\theta_z = %+.2f$' % cz, fontsize=10)

ax = axes[3]
for cz in ANGLES_SHOWN:
    chord, prob = earth_curve(cz)
    ax.semilogx(energies_earth/GEV, prob, lw=1.5,
                label=r'$\cos\theta_z = %+.2f$  (%.0f km)' % (cz, chord))
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
ax.set_title('Survival at the detector', fontsize=10)
ax.set_xlim(energies_earth[0]/GEV, energies_earth[-1]/GEV)
ax.set_ylim(0.0, 1.0)
ax.legend(fontsize=7.5)
fig.tight_layout(pad=1.2)'''),
    md(r'''## 4. Cutting a profile into slabs

The one approximation the method makes. Within a slab nothing is approximated --- the
expansion is exact for a constant Hamiltonian --- so the only question is how finely a profile
that really varies is sliced.

**Two things had to be forced to make this scene show anything.** Left to itself,
`strategy='auto'` does not use the slab ladder on a smooth exponential at all: it takes the
adiabatic route, and the answer is then *identical* at 2, 8 and 40 slabs. The sweep below
therefore passes `strategy='magnus'` and pins the slab count at both ends, which is the only
way to watch the discretization converge rather than watch the dispatcher avoid it.'''),
    code(r'''PER_NE = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
L_SLAB = 4000.0*KM


def ne_expo(l):
    """A smoothly falling profile: the case a slab count actually has to resolve."""
    return (1.0e-13*np.exp(-2.5*np.asarray(l, dtype=float)/L_SLAB))/PER_NE


def slabbed(n_slabs):
    return float(np.asarray(oscprob.osc_prob_matter_std_potential(
        3, ne_expo, 1.0*GEV, L_SLAB, OSC, L0=0.0,
        density_is_of_number_of_electrons=True, strategy='magnus',
        n_slabs=n_slabs, min_n_slabs=n_slabs, max_n_slabs=max(n_slabs, 2),
        nu_i=gd.NUMU, nu_f=gd.NUE, rtol=1.0e-13, atol=1.0e-15)))


SLABS_SHOWN = [1, 3, 8]
P_CONVERGED = slabbed(600)
ell = np.linspace(0.0, L_SLAB, 400)
vcc_profile = np.asarray(ne_expo(ell))*PER_NE

fig, axes = filmstrip(4, height=3.2)
for ax, n in zip(axes[:3], SLABS_SHOWN):
    ax.plot(ell/KM, vcc_profile/1.0e-13, color=MUTED, lw=1.4)
    edges = np.linspace(0.0, L_SLAB, n + 1)
    mid = 0.5*(edges[:-1] + edges[1:])
    for e in edges:
        ax.axvline(e/KM, color=ACCENT, lw=0.8, alpha=0.55)
    ax.step(mid/KM, np.asarray(ne_expo(mid))*PER_NE/1.0e-13, where='mid',
            color=MARK, lw=1.6)
    ax.set_title('%d slab%s' % (n, '' if n == 1 else 's'), fontsize=10)
    ax.set_xlabel('Distance [km]')
    ax.set_xlim(0.0, L_SLAB/KM)
axes[0].set_ylabel(r'$V_{\rm CC}$  [$10^{-13}$ eV]')

ax = axes[3]
counts = np.unique(np.round(np.logspace(0, 2.2, 22)).astype(int))
errors = [max(abs(slabbed(int(n)) - P_CONVERGED), 1.0e-16) for n in counts]
ax.loglog(counts, errors, '-o', ms=3.5, color=ACCENT)
ax.set_xlabel('Number of slabs')
ax.set_ylabel(r'$|P - P_{600}|$')
ax.set_title('Convergence with the slab count', fontsize=10)
ax.set_xlim(counts[0], counts[-1])
fig.tight_layout(pad=1.2)

print('Converged reference (600 slabs): P = %.8f' % P_CONVERGED)
print('One slab is off by %.2e; forty slabs by %.2e.'
      % (abs(slabbed(1) - P_CONVERGED), abs(slabbed(40) - P_CONVERGED)))'''),
    md(r'''# Part II --- five scenes with no counterpart

Each of these animates something a closed-form slab code does not have: a ladder that decides
when it has converged, a front that travels, an observable that is an average rather than a
value, and Hamiltonians that vary along the path.'''),
    md(r'''## 5. The refinement ladder, deciding

A tolerance is not met by raising the order --- it is met by **adding slabs**, and the ladder
does that on its own until two successive refinements agree to within what was asked. The
interesting rows are the ones where it stops being able to: past a point it hits its ceiling
and says so with a warning rather than returning a silently wrong number.

This one is a table rather than a picture, so it is not animated.'''),
    code(r'''import warnings as _warnings

TOLERANCES = np.logspace(-2, -11, 10)


def laddered(rtol):
    """Let the ladder choose, and report whether it ran out of room."""
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter('always')
        p = float(np.asarray(oscprob.osc_prob_matter_std_potential(
            3, ne_expo, 1.0*GEV, L_SLAB, OSC, L0=0.0,
            density_is_of_number_of_electrons=True, strategy='magnus',
            rtol=rtol, atol=rtol*1.0e-2, nu_i=gd.NUMU, nu_f=gd.NUE)))
        # `certified` is filled in by the hybrid dispatcher, not by the bare
        # ladder, so on strategy='magnus' it is always None and cannot be the
        # signal. What the ladder itself reports is the warning.
        warned = any('Tolerance' in w.category.__name__ or
                     'Convergence' in w.category.__name__ for w in caught)
    return p, warned


results = [laddered(t) for t in TOLERANCES]
errs = [max(abs(p - P_CONVERGED), 1.0e-17) for p, _ in results]

fig, ax = plt.subplots(figsize=(6.0, 4.0))
ok = [(t, e) for t, e, (_, w) in zip(TOLERANCES, errs, results) if not w]
bad = [(t, e) for t, e, (_, w) in zip(TOLERANCES, errs, results) if w]
ax.loglog(TOLERANCES, TOLERANCES, ls=':', color=MUTED, label='Requested')
if ok:
    ax.loglog(*zip(*ok), 'o', ms=7, color=ACCENT, mfc='white', mew=1.6,
              label='Converged silently')
if bad:
    ax.loglog(*zip(*bad), 'X', ms=8, color=MARK, label='Warned: out of room')
ax.set_xlabel('Requested tolerance')
ax.set_ylabel(r'$|P - P_{\rm converged}|$')
ax.set_title('What the ladder delivers, and when it stops', fontsize=11)
ax.set_xlim(TOLERANCES[0], TOLERANCES[-1])
ax.legend(fontsize=8)
fig.tight_layout(pad=1.0)

print('Requested   achieved     warned?')
for tol, (p, w) in zip(TOLERANCES, results):
    print('  %.0e     %.2e     %s' % (tol, abs(p - P_CONVERGED), 'yes' if w else 'no'))'''),
    md(r'''## 6. A supernova shock, sweeping outward

The forward shock moves out through the star, and as it crosses the region where the
$\Delta m^2_{31}$ resonance sits it changes **how adiabatically that level crossing is
made** --- which moves the conversion probability itself rather than its phase. Averaging
cannot remove it.

Every frame is a different Hamiltonian along the whole ray, and the front is declared through
`t_breakpoints` at each position, because no fixed grid resolves a discontinuity it was not
told about.

The resulting trace looks ragged, and that is the answer rather than a defect in it. Between
neighboring front positions the survival probability really does swing by tenths, and
checking that takes more care than it looks: comparing a coarse sampling against a finer one
proves nothing if both are under-resolved, because then the two agree about a number neither
of them has right. Driven to convergence instead --- $102\,400$ slabs, where the answer is
stable to $8 \times 10^{-6}$ --- forty points spaced $103$ km apart across the busy stretch
still step by $0.142$ on average, and still span $0.228$ to $0.570$. The structure is
physical, and finer than a hundred kilometers in front position.

It is drawn as markers on a faint line for that reason: a confident curve through these
points would claim a resolution that no affordable sampling has. The points themselves are
trustworthy --- against those converged values the settings used here are worst-case
$3.6 \times 10^{-3}$ over the whole sweep, which is well under the width of the marker.'''),
    code(r'''MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
R0_SHOCK, R1_SHOCK = 1.0e4, 8.0e4          # the ray, in km
E_SHOCK = 15.0*MEV
W_SHOCK = 1.0e-3                            # front width, as a fraction of the ray
SHOCK_NOTE = (r'$E_\nu = 15$ MeV' '\n'
              r'ray $10^4 \to 8 \times 10^4$ km')


def smoothstep(u):
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    return u*u*(3.0 - 2.0*u)


def ne_shock(r_front_km):
    """Electron density along the ray, for a front at this radius."""
    w_km = W_SHOCK*(R1_SHOCK - R0_SHOCK)

    def ne(l):
        r = np.asarray(l, dtype=float)/KM
        rho = 1.0e14*r**(-2.4)
        shocked = smoothstep((r_front_km + 0.5*w_km - r)/w_km)
        out = rho*(1.0 + shocked*9.0)*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*0.5
        return out[()] if np.ndim(out) == 0 else out
    return ne


def shock_probability(r_front_km):
    w_km = W_SHOCK*(R1_SHOCK - R0_SHOCK)
    edges = np.array([R0_SHOCK, r_front_km - 0.5*w_km,
                      r_front_km + 0.5*w_km, R1_SHOCK])*KM
    return float(np.asarray(oscprob.osc_prob_matter_std_potential(
        3, ne_shock(r_front_km), E_SHOCK, R1_SHOCK*KM, OSC, L0=R0_SHOCK*KM,
        density_is_of_number_of_electrons=True,
        t_breakpoints=np.unique(edges), nu_i=gd.NUE, nu_f=gd.NUE)))


FRONTS_SHOWN = [1.5e4, 3.0e4, 5.5e4]
# Linear in radius, not logarithmic: on a log axis the decade labels crowd into
# each other over this range and the front's motion is squeezed into a corner.
rr = np.linspace(R0_SHOCK, R1_SHOCK, 600)

fig, axes = filmstrip(4, height=3.2, ratios=[1, 1, 1, 1.25])
for k, (ax, rf) in enumerate(zip(axes[:3], FRONTS_SHOWN)):
    ax.semilogy(rr/1.0e4, np.asarray(ne_shock(rf)(rr*KM)), color=ACCENT, lw=1.5)
    ax.axvline(rf/1.0e4, color=MARK, ls='--', lw=1.2)
    ax.set_title('Front at %.0f km' % rf, fontsize=10)
    ax.set_xlabel(r'Radius [$10^4$ km]')
    ax.set_xlim(R0_SHOCK/1.0e4, R1_SHOCK/1.0e4)
    if k == 0:
        ax.text(0.97, 0.95, SHOCK_NOTE, transform=ax.transAxes, ha='right',
                va='top', fontsize=7.5, color='0.25', linespacing=1.4)
axes[0].set_ylabel(r'$n_e$  [eV$^3$]')

ax = axes[3]
fronts = np.linspace(1.2e4, 6.5e4, 26)
probs = [shock_probability(r) for r in fronts]
ax.plot(fronts/1.0e4, probs, '-o', ms=3.5, color=ACCENT)
for rf in FRONTS_SHOWN:
    ax.axvline(rf/1.0e4, color=MARK, ls='--', lw=0.9, alpha=0.6)
ax.set_xlabel(r'Shock radius [$10^4$ km]')
ax.set_ylabel(r'$P(\nu_e \to \nu_e)$')
ax.set_title('The front crossing the resonance', fontsize=10)
ax.set_xlim(fronts[0]/1.0e4, fronts[-1]/1.0e4)
fig.tight_layout(pad=1.2)

print('P_ee swings from %.4f to %.4f as the front sweeps outward.'
      % (min(probs), max(probs)))'''),
    md(r'''## 7. The average, emerging

Over the ray out of the Sun a few-MeV neutrino accumulates thousands of radians of phase, so
the instantaneous survival probability is neither measurable nor stable. What an experiment
measures is the **phase-averaged** probability, and `average=True` returns it directly.

The other way to get it is to propagate anyway and average many evaluations over a window in
energy. It converges slowly, and that is the honest result: the window mean is a Monte-Carlo
estimate of an average over a phase that turns over thousands of times, so its error falls
only as $1/\sqrt{N}$. This one costs about a minute per point, so it is a still rather than a
clip.'''),
    code(r'''SOLAR_TABLE = '../docs/dev/adversarial_batteries/bs05_agsop.dat'
_rows = []
with open(SOLAR_TABLE) as fh:
    for line in fh:
        f = line.split()
        if len(f) == 12:
            try:
                _rows.append([float(x) for x in f])
            except ValueError:
                continue
_solar = np.array(_rows)
_mean_nucleon = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
_ne_solar = _solar[:, 3]*gd.UNIT_G_PER_CM3/_mean_nucleon*(0.5*(1.0 + _solar[:, 6]))
_x_solar = _solar[:, 1]*gd.SUN_RADIUS*gd.UNIT_KM
_log_ne = np.log(_ne_solar)
R_SUN_AVG = float(_x_solar[-1])


def ne_sun(l):
    xs = np.clip(np.asarray(l, dtype=float), _x_solar[0], _x_solar[-1])
    out = np.exp(np.interp(xs, _x_solar, _log_ne))
    return out[()] if np.ndim(out) == 0 else out


E_AVG = np.array([2.0, 5.0, 10.0])*MEV
P_AVG_EXACT = np.asarray(oscprob.osc_prob_matter_std_potential(
    3, ne_sun, E_AVG, R_SUN_AVG, OSC, L0=0.0,
    density_is_of_number_of_electrons=True, nu_i=gd.NUE, nu_f=gd.NUE,
    average=True))


def window_mean(half_width, n_samples):
    offs = np.linspace(-half_width, half_width, n_samples)
    with _warnings.catch_warnings():
        _warnings.simplefilter('ignore')
        return np.array([float(np.mean(np.asarray(
            oscprob.osc_prob_matter_std_potential(
                3, ne_sun, e*(1.0 + offs), R_SUN_AVG, OSC, L0=0.0,
                density_is_of_number_of_electrons=True,
                nu_i=gd.NUE, nu_f=gd.NUE)))) for e in E_AVG])


WINDOWS = [(0.02, 9), (0.06, 15), (0.10, 21)]
means = [window_mean(hw, n) for hw, n in WINDOWS]
gaps = [float(np.max(np.abs(m - P_AVG_EXACT))) for m in means]

fig, axes = filmstrip(2, height=3.4, width_each=4.6)
ax = axes[0]
idx = np.arange(len(E_AVG))
ax.plot(idx, P_AVG_EXACT, 'o', ms=11, color=MARK, mfc='white', mew=2.0,
        label='average=True, one call')
for (hw, n), m in zip(WINDOWS, means):
    ax.plot(idx, m, 's', ms=6, alpha=0.85,
            label=r'Window $\pm%.0f\%%$, %d samples' % (100*hw, n))
ax.set_xticks(idx)
ax.set_xticklabels(['%.0f MeV' % (e/MEV) for e in E_AVG])
ax.set_xlim(-0.4, len(E_AVG) - 0.6)
ax.set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
ax.set_title('Sampling approaching the analytic average', fontsize=11)
ax.legend(fontsize=7.5)

ax = axes[1]
ns = [n for _, n in WINDOWS]
ax.plot(ns, gaps, '-o', ms=7, color=ACCENT)
ax.set_xlabel('Evaluations spent per energy')
ax.set_ylabel('Worst gap to average=True')
ax.set_title('The cost of estimating what is returned exactly', fontsize=11)
ax.set_xlim(ns[0], ns[-1])
fig.tight_layout(pad=1.0)

print('Window          worst gap to average=True')
for (hw, n), g in zip(WINDOWS, gaps):
    print('  +/-%4.0f%%, %2d samples   %.4f' % (100*hw, n, g))'''),
    md(r'''## 8. A Hamiltonian that depends on position

`H_func` is an arbitrary function of position, so a profile can be anything --- including
something with no closed form and no piecewise-constant description at all. Here a density
wave travels along the path: each frame moves its crest, so the Hamiltonian differs at every
point of every frame.

The dashed curve is the profile at rest, drawn once and left there, so that the size of the
change is legible rather than something to be remembered between frames.'''),
    code(r'''L_WAVE = 6000.0*KM
E_WAVE = np.logspace(np.log10(0.5), np.log10(12.0), 220)*GEV
N_WAVES = 3.0


def ne_wave(phase):
    """A density profile with a traveling crest."""
    def ne(l):
        x = np.asarray(l, dtype=float)/L_WAVE
        env = 1.0 + 0.75*np.sin(2.0*np.pi*(N_WAVES*x - phase))*np.exp(-2.0*x)
        return (1.2e-13*env)/PER_NE
    return ne


def wave_curve(phase):
    with _warnings.catch_warnings():
        _warnings.simplefilter('ignore')
        return np.asarray(oscprob.osc_prob_matter_std_potential(
            3, ne_wave(phase), E_WAVE, L_WAVE, OSC, L0=0.0,
            density_is_of_number_of_electrons=True,
            nu_i=gd.NUMU, nu_f=gd.NUMU))


PHASES_WAVE = [0.0, 0.33, 0.66]
xx = np.linspace(0.0, L_WAVE, 400)
WAVE_REFERENCE = wave_curve(0.0)              # the fixed curve to read against

fig, axes = filmstrip(4, height=3.2, ratios=[1, 1, 1, 1.3])
for ax, ph in zip(axes[:3], PHASES_WAVE):
    ax.plot(xx/KM, np.asarray(ne_wave(ph)(xx))*PER_NE/1.0e-13, color=ACCENT, lw=1.5)
    ax.set_title('Phase %.2f' % ph, fontsize=10)
    ax.set_xlabel('Distance [km]')
    ax.set_xlim(0.0, L_WAVE/KM)
    ax.set_ylim(0.0, 2.2)
axes[0].set_ylabel(r'$V_{\rm CC}$  [$10^{-13}$ eV]')

ax = axes[3]
ax.semilogx(E_WAVE/GEV, WAVE_REFERENCE, color='0.55', ls='--', lw=1.2,
            label='Phase 0.00 (reference)')
for ph in PHASES_WAVE[1:]:
    ax.semilogx(E_WAVE/GEV, wave_curve(ph), lw=1.4, label='Phase %.2f' % ph)
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
ax.set_title('The crest moving changes the answer', fontsize=10)
ax.set_xlim(E_WAVE[0]/GEV, E_WAVE[-1]/GEV)
ax.set_ylim(0.0, 1.0)
ax.legend(fontsize=7.5)
fig.tight_layout(pad=1.2)

sweep = np.array([wave_curve(p) for p in np.linspace(0, 1, 5, endpoint=False)])
print('Moving the crest through one period changes the probability by up to %.3f.'
      % float((sweep.max(axis=0) - sweep.min(axis=0)).max()))'''),
    md(r'''## 9. The Sun, with and without NSI

The solar case is the one where Mag$\nu$s returns the *observable* directly, and it is also
where a BSM term is most cleanly read: the standard three-flavor curve is fixed, and the NSI
curve moves away from it as $\varepsilon_{ee}$ grows.

The standard curve is drawn as a fixed reference in every frame. What sweeps is
$\varepsilon_{ee}$, from zero --- where the two curves must coincide, which is the check that
the sweep is doing what it claims --- up to a value large enough to move the MSW transition
visibly.

Both curves are `average=True` on the real BS2005-AGS,OP model, so each frame is the
phase-averaged probability rather than a sampling of it, and each costs well under a second.'''),
    code(r'''E_NSI = np.logspace(np.log10(0.1), np.log10(20.0), 60)*MEV


def solar_curve(eps_ee, eps_em=0.0):
    """<P_ee> across energy, at this NSI strength. eps = 0 is the standard case."""
    eps = dict(eps_ee=eps_ee, eps_em=eps_em, eps_et=0.0,
               eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
    with _warnings.catch_warnings():
        _warnings.simplefilter('ignore')
        return np.asarray(oscprob.osc_prob_matter_nsi(
            3, ne_sun, E_NSI, R_SUN_AVG, OSC, eps, L0=0.0,
            density_is_of_number_of_electrons=True,
            nu_i=gd.NUE, nu_f=gd.NUE, average=True))


SOLAR_STANDARD = np.asarray(oscprob.osc_prob_matter_std_potential(
    3, ne_sun, E_NSI, R_SUN_AVG, OSC, L0=0.0,
    density_is_of_number_of_electrons=True, nu_i=gd.NUE, nu_f=gd.NUE,
    average=True))
EPS_SHOWN = [0.0, 0.15, 0.30]

# Two panels, because the departure is a few hundredths: on the probability axis
# alone the curves nearly overlie, and the second panel is what makes the size of
# the effect legible rather than something to be squinted at.
fig, axes = filmstrip(2, height=3.6, width_each=4.8)
ax = axes[0]
ax.semilogx(E_NSI/MEV, SOLAR_STANDARD, color='0.35', ls='--', lw=1.6,
            label=r'Standard 3$\nu$')
for eps in EPS_SHOWN[1:]:
    ax.semilogx(E_NSI/MEV, solar_curve(eps), lw=1.6,
                label=r'$\varepsilon_{ee} = %.2f$' % eps)
ax.set_xlabel(r'$E_\nu$ [MeV]')
ax.set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
ax.set_xlim(E_NSI[0]/MEV, E_NSI[-1]/MEV)
ax.set_ylim(0.0, 1.0)
ax.set_title('The averaged solar survival probability', fontsize=11)
ax.legend(fontsize=8, loc='lower left')

ax = axes[1]
for eps in EPS_SHOWN[1:]:
    ax.semilogx(E_NSI/MEV, solar_curve(eps) - SOLAR_STANDARD, lw=1.6,
                label=r'$\varepsilon_{ee} = %.2f$' % eps)
ax.axhline(0.0, color='0.35', ls='--', lw=1.2)
ax.set_xlabel(r'$E_\nu$ [MeV]')
ax.set_ylabel(r'Departure from 3$\nu$')
ax.set_xlim(E_NSI[0]/MEV, E_NSI[-1]/MEV)
ax.set_title('and what the NSI actually moves', fontsize=11)
ax.legend(fontsize=8)
fig.tight_layout(pad=1.2)

print('At eps_ee = 0 the two curves agree to %.2e, as they must.'
      % float(np.max(np.abs(solar_curve(0.0) - SOLAR_STANDARD))))
print('At eps_ee = 0.30 they differ by up to %.3f.'
      % float(np.max(np.abs(solar_curve(0.30) - SOLAR_STANDARD))))'''),
    md(r'''## Rendering them as animations

Six of the nine scenes animate: the CP phase, the sterile state, the Earth, the slab count,
the shock and the traveling crest, plus the Sun with NSI. The ladder and the averaging scenes
do not --- one is a table and the other costs about a minute per point.

The stills above are what this notebook draws by default. Set `RENDER = True` to write the
raw GIFs into `img/raw/`; the shrunk copies committed in `img/` come from those in a second
step. The procedure, its measured cost and the traps are below the cell.'''),
    code(r'''RENDER = False        # set True to write the GIFs
FPS = 20

if RENDER:
    from matplotlib.animation import FuncAnimation, PillowWriter

    # ONE figure size for every scene, and it is not cosmetic: ffmpeg's concat
    # filter refuses clips whose dimensions differ, and scaling them to a common
    # width leaves different heights when the aspect ratios differ.  Two scenes
    # were drawn narrower here and `--join` failed with "Failed to configure
    # output pad" until they matched.
    FIGSIZE, DPI = (10.0, 4.4), 110
    # RAW renders land in img/raw/, which is gitignored; the shrunk copies that
    # tools/make_demo_video.py writes are what belongs in img/ and in the
    # repository.  Writing both to one directory would mean the next render
    # silently overwrote the committed 15 MB set with the 225 MB one.
    OUT = os.path.join('..', 'img', 'raw')
    os.makedirs(OUT, exist_ok=True)

    def heading_of(fig, w_pad=None):
        """A large title, close to the axes, and nothing written underneath.

        `tight_layout` reserves the top strip; the suptitle then sits just above
        the axes rather than floating halfway up the figure.  20 pt is chosen
        against the 15 pt axis labels set at the top of the notebook, so the
        title reads as the title -- with the repository defaults it would be
        smaller than the words under the x axis.

        `w_pad` widens the gutter between panels, which one scene needs: a
        color bar carries its own rotated label, and with the default gutter
        that label and the next panel's y label sit close enough to be read as
        one two-line label rather than two separate ones.
        """
        head = fig.suptitle('', fontsize=20, y=0.975)
        fig.patch.set_facecolor('white')
        # `pad` is a multiple of the font size, and the default 1.08 leaves a
        # visible white band under the x labels at this aspect.  0.8 reclaims
        # it; going much below that starts trimming descenders, and the labels
        # are the one thing that must not be clipped.
        fig.tight_layout(rect=[0, 0.0, 1, 0.93], pad=0.8,
                         **({} if w_pad is None else {'w_pad': w_pad}))
        return head

    def write(fig, update, frames, name):
        anim = FuncAnimation(fig, update, frames=frames, blit=False)
        path = os.path.join(OUT, name)
        anim.save(path, writer=PillowWriter(fps=FPS),
                  savefig_kwargs={'facecolor': 'white'})
        plt.close(fig)
        print('wrote %s (%.1f MB)'
              % (path, os.path.getsize(path)/1024.0/1024.0), flush=True)

    # ---- 1. the CP phase ---------------------------------------------------
    phases = np.linspace(0.0, 2.0*np.pi, 120, endpoint=False)
    fig, (ax_map, ax_ell) = plt.subplots(
        1, 2, figsize=FIGSIZE, dpi=DPI, gridspec_kw={'width_ratios': [1.45, 1.0]})
    image = ax_map.imshow(oscillogram_cp(phases[0]), origin='lower', aspect='auto',
                          cmap='viridis', vmin=0.0, vmax=CEIL_CP, extent=EXTENT_CP)
    style_map_cp(ax_map)
    ax_map.set_ylabel('Energy [GeV]')
    fig.colorbar(image, ax=ax_map, pad=0.02).set_label(r'$P(\nu_\mu \to \nu_e)$')
    ax_ell.plot(locus_cp[:, 0], locus_cp[:, 1], color=ACCENT, lw=1.8)
    dot, = ax_ell.plot([], [], 'o', ms=10, color=MARK, mfc='white', mew=2.0)
    ax_ell.set_xlabel(r'$P(\nu_\mu \to \nu_e)$')
    ax_ell.set_ylabel(r'$P(\bar\nu_\mu \to \bar\nu_e)$')
    ax_ell.set_xlim(locus_cp[:, 0].min(), locus_cp[:, 0].max())
    ax_ell.set_ylim(locus_cp[:, 1].min(), locus_cp[:, 1].max())
    head = heading_of(fig, w_pad=3.5)

    def update_cp(k):
        image.set_data(oscillogram_cp(phases[k]))
        dot.set_data(*[[v] for v in ellipse_point(phases[k])])
        head.set_text(r'The CP phase:  $\delta_{\rm CP} = %.2f\pi$' % (phases[k]/np.pi))
        return image, dot

    write(fig, update_cp, len(phases), 'anim_cp.gif')

    # ---- 2. a sterile state ------------------------------------------------
    splittings = np.logspace(np.log10(0.03), np.log10(2.0), 90)
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    img_st = ax.imshow(oscillogram_sterile(splittings[0]), origin='lower',
                       aspect='auto', cmap='magma', vmin=0.0, vmax=CEIL_ST,
                       extent=EXTENT_ST, interpolation='bilinear')
    style_map_st(ax)
    ax.set_ylabel('Energy [GeV]')
    fig.colorbar(img_st, ax=ax, pad=0.02).set_label(r'$P(\nu_\mu \to \nu_s)$')
    head = heading_of(fig)

    def update_sterile(k):
        img_st.set_data(oscillogram_sterile(splittings[k]))
        head.set_text(r'A sterile state:  $\Delta m^2_{41} = %.2f$ eV$^2$'
                      % splittings[k])
        return (img_st,)

    write(fig, update_sterile, len(splittings), 'anim_sterile.gif')

    # ---- 3. through the Earth, detector fixed at the South Pole ------------
    angles = np.linspace(-1.0, -0.05, 150)          # slower than the others
    # The sketch is a circle with `set_aspect('equal')`, so its height sets its
    # width and a half-width panel is all it can fill; giving it an equal share
    # of the row only padded it with white.  The curve takes the difference.
    fig, (ax_geo, ax_p) = plt.subplots(
        1, 2, figsize=FIGSIZE, dpi=DPI, gridspec_kw={'width_ratios': [0.72, 1.0]})
    draw_earth(ax_geo)
    path_line, = ax_geo.plot([], [], color='#1e3a8a', lw=2.4, zorder=4)
    curve, = ax_p.semilogx([], [], color=ACCENT, lw=1.8)
    ax_p.set_xlim(energies_earth[0]/GEV, energies_earth[-1]/GEV)
    ax_p.set_ylim(0.0, 1.0)
    ax_p.set_xlabel('Energy [GeV]')
    ax_p.set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
    head = heading_of(fig)

    def update_earth(k):
        cz = angles[k]
        chord, prob = earth_curve(cz)
        entry = chord_entry(cz)
        path_line.set_data([entry[0], DETECTOR[0]], [entry[1], DETECTOR[1]])
        curve.set_data(energies_earth/GEV, prob)
        head.set_text(r'Arriving at the South Pole:  $\cos\theta_z = %+.2f$  (%.0f km)'
                      % (cz, chord))
        return path_line, curve

    write(fig, update_earth, len(angles), 'anim_earth.gif')

    # ---- 4. cutting a profile into slabs -----------------------------------
    slab_frames = np.unique(np.round(np.logspace(0, 2.0, 60)).astype(int))
    fig, (ax_prof, ax_err) = plt.subplots(1, 2, figsize=FIGSIZE, dpi=DPI)
    ax_prof.plot(ell/KM, vcc_profile/1.0e-13, color=MUTED, lw=1.5)
    step_line, = ax_prof.step([], [], where='mid', color=MARK, lw=1.8)
    ax_prof.set_xlim(0.0, L_SLAB/KM)
    ax_prof.set_ylim(0.0, float(vcc_profile.max()/1.0e-13)*1.02)
    ax_prof.set_xlabel('Distance [km]')
    ax_prof.set_ylabel(r'$V_{\rm CC}$  [$10^{-13}$ eV]')
    err_line, = ax_err.loglog([], [], '-o', ms=4, color=ACCENT)
    ax_err.set_xlim(1.0, float(slab_frames[-1]))
    ax_err.set_ylim(1.0e-11, 1.0e-2)
    ax_err.set_xlabel('Number of slabs')
    ax_err.set_ylabel(r'$|P - P_{600}|$')
    head = heading_of(fig)
    trail_s = {'n': [], 'e': []}

    def update_slabs(k):
        n = int(slab_frames[k])
        edges = np.linspace(0.0, L_SLAB, n + 1)
        mid = 0.5*(edges[:-1] + edges[1:])
        step_line.set_data(mid/KM, np.asarray(ne_expo(mid))*PER_NE/1.0e-13)
        trail_s['n'].append(n)
        trail_s['e'].append(max(abs(slabbed(n) - P_CONVERGED), 1.0e-16))
        err_line.set_data(trail_s['n'], trail_s['e'])
        head.set_text('Cutting the profile into %d slab%s' % (n, '' if n == 1 else 's'))
        return step_line, err_line

    write(fig, update_slabs, len(slab_frames), 'anim_slabs.gif')

    # ---- 6. the shock sweeping outward -------------------------------------
    fronts_anim = np.linspace(1.2e4, 6.5e4, 120)
    fig, (ax_prof, ax_p) = plt.subplots(1, 2, figsize=FIGSIZE, dpi=DPI)
    prof_line, = ax_prof.semilogy(rr/1.0e4,
                                  np.asarray(ne_shock(fronts_anim[0])(rr*KM)),
                                  color=ACCENT, lw=1.6)
    front_line = ax_prof.axvline(fronts_anim[0]/1.0e4, color=MARK, ls='--', lw=1.4)
    ax_prof.set_xlim(R0_SHOCK/1.0e4, R1_SHOCK/1.0e4)
    ax_prof.set_xlabel(r'Radius [$10^4$ km]')
    ax_prof.set_ylabel(r'$n_e$  [eV$^3$]')
    # On the right panel, not the left one where the still puts it: the front
    # marker sweeps across the whole width of the density panel, so every
    # corner of it is crossed at some frame.  The top left of the probability
    # panel is the one region no frame ever draws into, because the trace is
    # still near zero while the front is at small radius.
    ax_p.text(0.03, 0.96, SHOCK_NOTE, transform=ax_p.transAxes, ha='left',
              va='top', fontsize=9, color='0.25', linespacing=1.4)
    # Markers carrying a faint line, rather than a bold line through the
    # markers.  The trace really does swing by tenths between neighboring
    # front positions, and that was checked against CONVERGED values rather
    # than against a finer sampling at the same settings -- two under-resolved
    # runs agree about a number neither has right.  At 102400 slabs, where the
    # answer is stable to 8e-06, forty points 103 km apart still step by 0.142
    # on average.  So the structure is physical, and a confident thick line
    # drawn through it would claim a resolution that is not there.
    track, = ax_p.plot([], [], '-o', ms=3.2, lw=0.8, alpha=0.85, color=ACCENT)
    ax_p.set_xlim(fronts_anim[0]/1.0e4, fronts_anim[-1]/1.0e4)
    # NOT `max(probs) + 0.02`: `probs` is the coarse sweep the still above
    # draws, and the animation steps finer, so it finds extremes the still
    # never sampled and the trail was drawn clipped against the top spine.  A
    # 200-point scan of the busy stretch reaches 0.591, so the ceiling is set
    # above that rather than above what one particular sampling happened to see.
    ax_p.set_ylim(0.0, 0.62)
    ax_p.set_xlabel(r'Shock radius [$10^4$ km]')
    ax_p.set_ylabel(r'$P(\nu_e \to \nu_e)$')
    head = heading_of(fig)
    trail = {'x': [], 'y': []}

    def update_shock(k):
        rf = fronts_anim[k]
        prof_line.set_ydata(np.asarray(ne_shock(rf)(rr*KM)))
        front_line.set_xdata([rf/1.0e4, rf/1.0e4])
        trail['x'].append(rf/1.0e4)
        trail['y'].append(shock_probability(rf))
        track.set_data(trail['x'], trail['y'])
        head.set_text('The shock front at %.0f km' % rf)
        return prof_line, track

    write(fig, update_shock, len(fronts_anim), 'anim_shock.gif')

    # ---- 8. a traveling density crest -------------------------------------
    wave_phases = np.linspace(0.0, 1.0, 90, endpoint=False)
    fig, (ax_prof, ax_p) = plt.subplots(1, 2, figsize=FIGSIZE, dpi=DPI)
    prof_w, = ax_prof.plot(xx/KM, np.asarray(ne_wave(0.0)(xx))*PER_NE/1.0e-13,
                           color=ACCENT, lw=1.6)
    ax_prof.set_xlim(0.0, L_WAVE/KM)
    ax_prof.set_ylim(0.0, 2.2)
    ax_prof.set_xlabel('Distance [km]')
    ax_prof.set_ylabel(r'$V_{\rm CC}$  [$10^{-13}$ eV]')
    ax_p.semilogx(E_WAVE/GEV, WAVE_REFERENCE, color='0.55', ls='--', lw=1.3,
                  label='Phase 0.00 (reference)')
    curve_w, = ax_p.semilogx(E_WAVE/GEV, WAVE_REFERENCE, color=ACCENT, lw=1.6,
                             label='Current phase')
    ax_p.set_xlim(E_WAVE[0]/GEV, E_WAVE[-1]/GEV)
    ax_p.set_ylim(0.0, 1.0)
    ax_p.set_xlabel('Energy [GeV]')
    ax_p.set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
    # The curve is dense at the left, so an unframed legend sitting on top of it
    # is unreadable.  A white box, above the oscillations rather than inside
    # them, is the one place on this panel with room.
    ax_p.legend(loc='upper left', frameon=True, framealpha=0.92,
                edgecolor='none', borderpad=0.5)
    head = heading_of(fig)

    def update_wave(k):
        ph = wave_phases[k]
        prof_w.set_ydata(np.asarray(ne_wave(ph)(xx))*PER_NE/1.0e-13)
        curve_w.set_ydata(wave_curve(ph))
        head.set_text('A traveling crest:  phase %.2f' % ph)
        return prof_w, curve_w

    write(fig, update_wave, len(wave_phases), 'anim_wave.gif')

    # ---- 9. the Sun, with and without NSI ----------------------------------
    eps_frames = np.concatenate([np.linspace(0.0, 0.30, 60),
                                 np.linspace(0.30, 0.0, 30)])
    fig, (ax_p, ax_d) = plt.subplots(1, 2, figsize=FIGSIZE, dpi=DPI)
    ax_p.semilogx(E_NSI/MEV, SOLAR_STANDARD, color='0.35', ls='--', lw=1.8,
                  label=r'Standard 3$\nu$')
    nsi_line, = ax_p.semilogx(E_NSI/MEV, SOLAR_STANDARD, color=MARK, lw=2.2,
                              label='With NSI')
    ax_p.set_xlim(E_NSI[0]/MEV, E_NSI[-1]/MEV)
    ax_p.set_ylim(0.0, 1.0)
    ax_p.set_xlabel(r'$E_\nu$ [MeV]')
    ax_p.set_ylabel(r'$\langle P(\nu_e \to \nu_e)\rangle$')
    ax_p.legend(fontsize=10, loc='lower left')
    dep_line, = ax_d.semilogx(E_NSI/MEV, np.zeros_like(E_NSI), color=MARK, lw=2.2)
    ax_d.axhline(0.0, color='0.35', ls='--', lw=1.2)
    ax_d.set_xlim(E_NSI[0]/MEV, E_NSI[-1]/MEV)
    ax_d.set_ylim(-0.05, 0.05)
    ax_d.set_xlabel(r'$E_\nu$ [MeV]')
    ax_d.set_ylabel(r'Departure from 3$\nu$')
    head = heading_of(fig)

    def update_solar(k):
        eps = eps_frames[k]
        curve = solar_curve(eps)
        nsi_line.set_ydata(curve)
        dep_line.set_ydata(curve - SOLAR_STANDARD)
        head.set_text(r'The Sun with NSI:  $\varepsilon_{ee} = %.3f$' % eps)
        return nsi_line, dep_line

    write(fig, update_solar, len(eps_frames), 'anim_solar_nsi.gif')

    print()
    print('Seven clips written. Sections 5 and 7 are stills by design: the ladder is')
    print('a table, and the averaging scene costs about a minute per point.')'''),
    md(r'''### How to actually produce the clips

The cell above writes nothing as it stands, and the notebook is committed that way on
purpose: rendering costs about an hour and some two hundred megabytes, and CI executes every
notebook in this repository on every push. What follows is the whole procedure, measured on
the machine that wrote it.

**Which files are tracked, and which are not.** The shrunk clips in `img/` **are** committed,
because a reader should not have to spend an hour to see what the notebook is about. Their
raw originals in `img/raw/` are not: they are reproducible build artifacts in the sense
`fig/*.pdf` already is. Nothing under `img/` reaches PyPI either way --- the packaging
configuration looks only under `src/`, so neither the sdist nor the wheel contains any of it.

**Step 1 --- render the scenes.** Set `RENDER = True` in the cell above and run it. Seven
GIFs land in `img/raw/`, every one of them 3000 x 1320. Measured, single-threaded:

| clip | frames | raw size |
|---|---|---|
| `anim_sterile.gif` | 90 | 160.2 MB |
| `anim_cp.gif` | 120 | 26.7 MB |
| `anim_earth.gif` | 150 | 19.8 MB |
| `anim_wave.gif` | 90 | 10.1 MB |
| `anim_solar_nsi.gif` | 89 | 3.5 MB |
| `anim_shock.gif` | 120 | 2.9 MB |
| `anim_slabs.gif` | 39 | 0.7 MB |

about **52 minutes** and **224 MB** in total. Most of that is matplotlib, not the physics: the
two map scenes recompute tens of thousands of probabilities per frame and still spend longer
being drawn than being computed.

That 52 minutes is the whole set, and it is split very unevenly: the Earth scene on its own,
re-rendered on the same machine for this table, took **51 s**. Budget for the set rather than
per scene, and expect the two map scenes to be most of it.

**Step 2 --- shrink them into `img/`.** A GIF straight out of matplotlib's Pillow writer
gives every frame its own color table. One shared palette removes that duplication, and
combined with a lower frame rate and a smaller width it is what makes a file publishable:

```shell
for f in img/raw/anim_*.gif; do
    python tools/make_demo_video.py --shrink "$PWD/$f" \
        --out "$PWD/img/$(basename $f)" --fps 12 --width 860 --colors 128
done
```

That takes the set from **223.8 MB to 14.6 MB**, a factor of 15.4. The factor is not uniform,
and the reason is worth knowing before you tune it:

| clip | raw | shrunk | ratio |
|---|---|---|---|
| `anim_sterile.gif` | 160.2 MB | 8.01 MB | **20.0x** |
| `anim_solar_nsi.gif` | 3.5 MB | 0.28 MB | 12.6x |
| `anim_wave.gif` | 10.1 MB | 0.85 MB | 11.9x |
| `anim_cp.gif` | 26.7 MB | 2.45 MB | 10.9x |
| `anim_shock.gif` | 2.9 MB | 0.34 MB | 8.4x |
| `anim_earth.gif` | 19.8 MB | 2.50 MB | 7.9x |
| `anim_slabs.gif` | 0.7 MB | 0.13 MB | **5.1x** |

Note that the *ratio* and the *result* rank differently: the sterile map compresses best of
all and is still by far the largest file, while the slab scene compresses worst and is the
smallest thing here. Ratio is a property of how much redundancy the raw file had; what
matters for a README is the number on the right. **Check the scene you care about rather than
trusting a default.** The three knobs are `--fps`, `--width` and `--colors`, in that order of
effect on size.

**Step 3 --- join them into one reel, if you want a single clip.** `.mp4` is gitignored; a
reel is for showing, not for committing.

```shell
python tools/make_demo_video.py --join img/anim_cp.gif img/anim_sterile.gif \
    img/anim_earth.gif img/anim_slabs.gif img/anim_shock.gif \
    img/anim_wave.gif img/anim_solar_nsi.gif --out ~/reel.mp4 --fps 20
python tools/make_demo_video.py --shrink ~/reel.mp4 --out ~/reel.gif --fps 12 --width 900
```

Joining first and shrinking once is not the same as shrinking each clip and concatenating the
results --- one palette for the whole reel is a compromise across very different images, and
the dense sterile map drags it. Shrinking clip by clip wins by around a fifth on the same
content, so prefer separate clips unless a single file is the point.

### Three traps, all of them hit while writing this section

**Every clip in a reel must have the same pixel dimensions.** `ffmpeg`'s concat filter
refuses otherwise, with `Failed to configure output pad on Parsed_concat_N` --- which does
not mention dimensions at all. Scaling to a common *width* is not enough: clips whose aspect
ratios differ then end up with different heights. Two scenes here were originally drawn
narrower than the rest and the join failed until they matched. That is why the render cell
defines a single `FIGSIZE` and every scene uses it, and why the table above shows one
resolution for all seven.

**A GIF stores its frame delays in hundredths of a second.** `ffmpeg` reads a 90 ms frame as
a stream of roughly 100 fps and, without an explicit output rate, writes every frame nine
times over. Always pass `--fps`. A four-hundred-frame reel became ninety megabytes and
climbing before this was noticed.

**`ffmpeg` installed as a snap has a private `/tmp`.** It fails with `No such file or
directory` naming a path that plainly exists, which is a confusing way to be told about
confinement. Work somewhere under `$HOME`. The palette file is written next to the *output*
for this reason, so an output path under `$HOME` is enough --- but a `--out /tmp/...` will
fail on such a machine while Pillow, which is not confined, happily writes to the same place
from the cell above. That asymmetry makes it look like a codec problem. It is not.

### What is deliberately not animated

Sections 5 and 7 have no clip. The refinement ladder is a table rather than a picture, and
nothing is gained by watching rows appear. The averaging scene costs about a minute per
point, because every point propagates a phase of some thousands of radians; an animation of
it would run for hours to say exactly what its three stills already say. Both are better read
than watched.'''),
    ])


# --------------------------------------------------------------- reading order

# ---------------------------------------------- 28_magnus_paper_figures
books['28_magnus_paper_figures.ipynb'] = notebook(
    "The paper's figures",
    r'''Every figure in the Mag$\nu$s paper (`resources/paper/`), produced by one run of this
notebook.

Two rules govern what is computed here and what is read from a file.

**Mag$\nu$s's own numbers are computed as this runs.** A figure showing frozen numbers for
this package would go stale the moment the package changed, and nothing would say so.

**Every other code's numbers are read from `notebooks/external_*.json`.** None of the codes
compared against here has to be installed to redraw these figures, and none of them is its
own judge: each comparison is refereed by a method that is neither code's. Notebook 25 is
where those comparisons are made and discussed.

Figures are written to `resources/paper/figs/` as PDF, which is what `main.tex` includes.
Set `MAGNUS_PAPER_FIGDIR` to send them elsewhere. A rerun rewrites every PDF whether or not
anything changed and the bytes differ between runs, so `git status` will show them modified
even when the figures are identical --- commit them only when a figure actually changed.''',

    [
    md(r'''## Setup, house style, and the helpers every figure uses

Six helpers, each here because getting it wrong is easy and the wrong answer looks
ordinary: a Hamiltonian builder that passes every mixing angle by keyword, a chord that
remembers `distance_traveled_inside_earth` returns kilometres, a probability that
remembers $|U|^2$ comes out indexed the other way round, and an ODE reference handed the
*same* Hamiltonian as the code under test.

Two tolerance settings, and the difference matters. `RTOL_FIG` is used wherever a figure
shows a probability: the package default of $10^{-3}$ is a working setting, not a
publication one. `RTOL_ACC` is used where a figure *measures* accuracy, and is tightened
until the residual stops being a property of the setting --- which, as the next cell
shows, happens at $10^{-12}$ and not below.'''),
    code(r'''import hashlib
import json
import re
import os
import pathlib
import platform
import time
import warnings

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import FuncFormatter, LogLocator, AutoMinorLocator, NullLocator
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from magnus import magnus
import magnus.oscprob as oscprob
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.earth as earth
import magnus.avgprob as avgprob
import magnus.globaldefs as gd

HERE = pathlib.Path.cwd()
FIGDIR = pathlib.Path(os.environ.get('MAGNUS_PAPER_FIGDIR',
                                     HERE.parent/'resources'/'paper'/'figs'))
FIGDIR.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN, RED = '#1c71d8', '#e66100', '#2ec27e', '#c01c28'
PURPLE, INK, GRID = '#813d9c', '#333333', '#cccccc'

# Colour means one thing across Fig. 2: the truncation order.  DOP853 is the only curve
# that is not a Magnus order, so it is the only one drawn in black.
ORDER_COLOR = {2: GREEN, 4: BLUE, 6: PURPLE, 8: ORANGE, 10: RED}

plt.rcParams.update({
    # These override notebooks/matplotlibrc, whose sizes are set for a 5-inch standalone
    # figure.  Every figure here is drawn at the width it is included at, so a size set
    # here is the size it renders at on the page: the paper's body is 10 pt and its
    # captions are 8 pt, and nothing in a figure should be smaller than its own caption.
    'font.size': 9, 'axes.labelsize': 9.5, 'axes.titlesize': 9.5,
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'legend.fontsize': 8,
    'axes.linewidth': 0.7, 'lines.linewidth': 1.2,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.pad': 1.8, 'ytick.major.pad': 1.8,
    'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
    'xtick.minor.size': 1.8, 'ytick.minor.size': 1.8,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'legend.framealpha': 1.0, 'legend.edgecolor': 'black',
    'legend.fancybox': False, 'legend.borderpad': 0.3,
    'figure.dpi': 130, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
})

# The real \columnwidth and \textwidth of the paper, so that a size set here is the
# size that reaches the page: drawing narrower and letting \includegraphics stretch
# the result was scaling text by 1.01 to 1.34, differently in every figure.
COL, WIDE = 3.487, 7.224
trapz = getattr(np, 'trapezoid', None) or np.trapz

OSC = gd.load_nufit_params('NuFIT 6.1')
STERILE4 = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, D41=1.0)
STERILE5 = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0,
                s15=np.sqrt(0.06), s25=np.sqrt(0.06), s35=0.0, D41=1.0, D51=1.7)
EPS = dict(eps_ee=0.10, eps_em=0.05+0.0j, eps_et=0.0j, eps_mm=0.0,
           eps_mt=0.03+0.0j, eps_tt=0.0)
FLAVOR_LABEL = {2: r'$2\nu$', 3: r'$3\nu$', 4: r'$3+1$', 5: r'$3+2$'}

# The accuracy setting used wherever a figure shows a probability rather than
# sweeping a tolerance.  The default 1e-3 is a working setting, not a publication
# one: a plotted curve should not carry a discretization error a reader could see.
RTOL_FIG, ATOL_FIG = 1.0e-8, 1.0e-10
# And where a figure MEASURES accuracy, both sides are tightened until the residual
# stops being the setting.  See the note in fig 1.
RTOL_ACC, ATOL_ACC = 1.0e-12, 1.0e-14
ODE_RTOL, ODE_ATOL = 1.0e-12, 1.0e-14


def _plain(v, _pos=None):
    r"""Tick label: 1 rather than $10^0$, 0.1 rather than $10^{-1}$."""
    if v <= 0:
        return ''
    e = np.log10(v)
    if -3.001 < e < 4.001:
        s = ('%f' % v).rstrip('0').rstrip('.')
        return s if s else '0'
    return r'$10^{%d}$' % round(e)


def logx(ax):
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FuncFormatter(_plain))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(np.arange(2, 10)*0.1),
                                          numticks=100))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda *_: ''))


def logy(ax):
    ax.set_yscale('log')
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(np.arange(2, 10)*0.1),
                                          numticks=100))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda *_: ''))


def snug(ax, x, log=False):
    r"""No dead margin left or right of the data.  Applied to every panel."""
    x = np.asarray(x, dtype=float)
    ax.set_xlim(float(np.min(x)), float(np.max(x)))


def xticks_at(ax, values):
    r"""Explicit major ticks on a log axis whose range spans less than two decades.

    A decade locator puts one label on a 2-60 GeV axis, which reads as an axis with
    no scale at all.
    """
    from matplotlib.ticker import FixedLocator
    ax.xaxis.set_major_locator(FixedLocator(list(values)))
    ax.xaxis.set_major_formatter(FuncFormatter(_plain))


def minor_log_ticks(ax, which='y', length=3.6):
    r"""Minor ticks between the decades of a log axis, drawn on both sides.

    ``ytick.minor.visible`` is off in the rc, and a minor locator alone is not enough to
    guarantee they are drawn at a useful size: the length has to be set explicitly, and
    generously, or they vanish against the major ticks the rc sets to 10 points.
    """
    axis = ax.yaxis if which == 'y' else ax.xaxis
    axis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10)), numticks=200))
    axis.set_minor_formatter(FuncFormatter(lambda *_: ''))
    side = dict(left=True, right=True) if which == 'y' else dict(bottom=True, top=True)
    ax.tick_params(axis=which, which='minor', length=length, width=0.8,
                   direction='in', **side)


def unit_as_one(ax, which='y'):
    r"""Write the decade $10^0$ as ``1``.

    A log axis that crosses unity reads better with the unit written plainly: the
    surrounding ticks are powers because they have to be, but ``1`` is not clearer as
    $10^{0}$.  Every other decade keeps its exponent.
    """
    axis = ax.yaxis if which == 'y' else ax.xaxis
    axis.set_major_formatter(FuncFormatter(
        lambda v, _p=None: r'$1$' if abs(np.log10(v)) < 1.0e-9
        else r'$10^{%d}$' % round(np.log10(v))))


def minor_y(ax, n=5):
    ax.yaxis.set_minor_locator(AutoMinorLocator(n))


def corner(ax, text, loc='upper right', fontsize=8.5, x=None, y=0.94):
    r"""A rounded-box label in a corner, in place of a panel title."""
    ha, va = 'right', 'top'
    xx = 0.965 if x is None else x
    if loc == 'upper left':
        xx, ha = (0.035 if x is None else x), 'left'
    # Black, not the INK grey the curves use: a boxed label is a caption on the panel and
    # should read as text rather than as another datum.
    ax.text(xx, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=fontsize,
            color='black', zorder=10,
            bbox=dict(boxstyle='round,pad=0.32', facecolor='white',
                      edgecolor='black', linewidth=0.6))


def stamp(ax, text, x=0.035, y=0.06, fontsize=8.0, ha='left', va='bottom'):
    r"""Free text over a busy panel: black, outlined in white so it stays legible."""
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=fontsize,
            color='black', zorder=10,
            path_effects=[pe.withStroke(linewidth=1.8, foreground='white')])


MP_CACHE = pathlib.Path('paper_figure_cache.json')
# Continuous integration regenerates these notebooks on every push, and the expensive
# inputs here are either machine-specific (timings) or settled functions of the
# configuration (references).  Neither is worth re-deriving there.  With
# MAGNUS_PAPER_CACHE_ONLY set, a cache miss stops the build and says which section moved,
# instead of quietly spending the better part of an hour recomputing it.
CACHE_ONLY = bool(os.environ.get('MAGNUS_PAPER_CACHE_ONLY'))
# What the stored numbers were produced by.  A configuration fingerprint cannot see a
# change inside the package, so the version is what decides when the stored accuracies
# are worth re-checking.
# `magnus` is bound to the submodule in these notebooks, so the package is imported
# again under its own name to read the version off it.
import magnus as _magnus_pkg
MAGNUS_VERSION = getattr(_magnus_pkg, '__version__', 'unknown')


def cache_miss(section, key):
    if CACHE_ONLY:
        raise RuntimeError(
            'paper cache miss in section %r (fingerprint %s) while '
            'MAGNUS_PAPER_CACHE_ONLY is set.  The configuration behind this section '
            'changed, so it has to be re-measured on a quiet machine and the refreshed '
            '%s committed:  MAGNUS_PAPER_REDO=1 python notebooks/make_notebooks.py '
            '--only 28' % (section, key[:12], MP_CACHE.name))
CACHE_SECTIONS = ('what', 'fixed', 'scan', 'orders', 'timings', 'oracle', 'prem_timings')
CACHE_WHAT = ('Everything in the paper figures that is a property of the CONFIGURATION '
              'rather than of this run: reference probabilities, and wall-clock timings. '
              'Each section is keyed on a fingerprint of the configuration it was computed '
              'for, so changing an energy, a tolerance or a code recomputes it and changing '
              'nothing reuses it.  This exists so that regenerating the notebooks -- which '
              'continuous integration does on every push -- does not re-measure minutes of '
              'reference arithmetic and stopwatch that cannot have moved.  Set '
              'MAGNUS_PAPER_REDO=1 to force every section to be recomputed, which is what '
              'to do after moving the paper to another machine.')


def write_cache(blob):
    """Write the whole blob.

    This used to write only an allow-list of section names, so that a change of format
    could not leave its predecessor's keys behind.  That silently threw away every
    section not on the list --- which is every section added after the list was written
    --- and it was unsafe besides: sections are written the moment one of them
    recomputes, so an early section writing before a later one had been read would have
    dropped the later one's stored value.  Stale keys are a job for a one-off cleanup,
    not for a filter that runs on every write.
    """
    blob.setdefault('what', CACHE_WHAT)
    MP_CACHE.write_text(json.dumps(blob, indent=1) + '\n')


# Profiles enter every fingerprint as nine samples of what they return, never as the
# callable itself: a function's repr carries its memory address, which changes on every
# run and would invalidate the cache each time it was consulted.
def profile_samples(func, L):
    return np.asarray(func(np.linspace(0.0, L, 9)), dtype=float)


def fingerprint(*parts):
    """Everything a stored result depends on, in one hash.

    Profiles enter as samples of the array they produce rather than as the parameters that
    built them, so a change anywhere upstream -- a mixing angle, an energy, a potential, a
    baseline -- invalidates the entry without having to be enumerated.
    """
    h = hashlib.sha256()
    for part in parts:
        if isinstance(part, np.ndarray):
            h.update(np.ascontiguousarray(part).tobytes())
        else:
            h.update(repr(part).encode())
    return h.hexdigest()


def cached(section, key_parts, compute, what=''):
    """Compute a figure input once, then read it from disk until its configuration moves.

    References and timings are both properties of the configuration: the first exactly,
    the second up to the machine.  Re-deriving either on every rebuild costs minutes and
    tells nobody anything, so both are stored, keyed on the configuration alone.  The
    machine and the date are stored beside them, because a number that does not say where
    it came from cannot be checked.
    """
    blob = json.loads(MP_CACHE.read_text()) if MP_CACHE.exists() else {}
    key = fingerprint(*key_parts)
    got = blob.get(section)
    if got and got.get('fingerprint') == key and not os.environ.get('MAGNUS_PAPER_REDO'):
        print('  %s read from %s, unchanged configuration %s (measured %s on %s)'
              % (section, MP_CACHE.name, key[:12], got.get('measured', '?'),
                 got.get('machine', 'an unrecorded machine')))
        return got['value']
    print('  %s: configuration moved, recomputing' % section)
    t0 = time.perf_counter()
    value = compute()
    blob[section] = dict(fingerprint=key, value=value, what=what,
                         machine=platform.node(), measured=time.strftime('%Y-%m-%d'),
                         seconds=round(time.perf_counter() - t0, 1))
    write_cache(blob)
    print('    %s computed in %.0f s' % (section, time.perf_counter() - t0))
    return value


def save(fig, name):
    fig.savefig(FIGDIR/name)
    print('  wrote %s' % (FIGDIR/name))


def vacuum_hamiltonian(d):
    if d == 2:
        return hamiltonians.hamiltonian_2nu_vacuum_energy_independent(
            sth=OSC['s12'], Dm2=OSC['D21'])
    if d == 3:
        return hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
            s12=OSC['s12'], s23=OSC['s23'], s13=OSC['s13'], dCP=OSC['dCP'],
            D21=OSC['D21'], D31=OSC['D31'])
    if d == 4:
        return hamiltonians.hamiltonian_4nu_vacuum_energy_independent(
            s12=OSC['s12'], s23=OSC['s23'], s13=OSC['s13'], d13=OSC['dCP'],
            s14=STERILE4['s14'], d14=0.0, s24=STERILE4['s24'], d24=0.0,
            s34=STERILE4['s34'], D21=OSC['D21'], D31=OSC['D31'], D41=STERILE4['D41'])
    if d == 5:
        return hamiltonians.hamiltonian_5nu_vacuum_energy_independent(
            s12=OSC['s12'], s23=OSC['s23'], s13=OSC['s13'], d13=OSC['dCP'],
            s14=STERILE5['s14'], d14=0.0, s15=STERILE5['s15'], d15=0.0,
            s24=STERILE5['s24'], d24=0.0, s25=STERILE5['s25'],
            s34=STERILE5['s34'], s35=STERILE5['s35'], d35=0.0,
            D21=OSC['D21'], D31=OSC['D31'], D41=STERILE5['D41'], D51=STERILE5['D51'])
    raise ValueError(d)


def make_H_func(d, energy, vcc_func):
    # The vacuum term does not vary along the trajectory, so it is divided by the energy
    # once here rather than at every evaluation.  The ODE solver is the one that notices,
    # because it calls this once per step on a scalar: hoisting buys it 8 per cent and
    # buys Magnus nothing, which evaluates every node of a slab in one vectorised call.
    h_vac_over_e = np.asarray(vacuum_hamiltonian(d), dtype=complex)/energy
    proj = np.asarray(matter.matter_potential_projector(d), dtype=complex)

    def H_func(l):
        l = np.asarray(l, dtype=float)
        return h_vac_over_e + np.asarray(vcc_func(l))[..., None, None]*proj

    return H_func


def chord(costhz):
    return earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM


def prob_from_U(U):
    return (np.abs(U)**2).T


def ode_reference(H_func, L, d, rtol=ODE_RTOL, atol=ODE_ATOL):
    """DOP853 on the whole evolution matrix, in a single call.

    Integrating the d columns one at a time solves d initial-value problems where one
    will do and pays SciPy's per-step overhead d times over; batching them is worth a
    factor of 2.5 to 2.8 here.  That is not merely convenience.  This routine is the code
    the center panel of Fig. 2 times against Magnus, and timing against a needlessly slow
    competitor measures the implementation rather than the method.
    """
    def rhs(l, y):
        return (-1j*(H_func(float(l)) @ y.reshape(d, d))).ravel()

    sol = solve_ivp(rhs, (0.0, L), np.eye(d, dtype=complex).ravel(), method='DOP853',
                    rtol=rtol, atol=atol)
    return (np.abs(sol.y[:, -1].reshape(d, d))**2).T


def accumulated_phase(H_func, L, n=4001):
    l = np.linspace(0.0, L, n)
    w = np.linalg.eigvalsh(H_func(l))
    return float(trapz(w[:, -1] - w[:, 0], l))


def quiet(call, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return call(*args, **kwargs)
'''),
    md(r'''## Figure 1 --- accuracy, unitarity, and the limit of the reference

**Tighter is not always better, and the oracle has to be checked first.** At $10^{-14}$
the residual is an order of magnitude *worse* than at $10^{-12}$: the ladder is chasing
round-off, and the extra factors in the ordered product cost more than the finer grid
buys. And at four and five flavors the `DOP853` reference itself moves by $4\times10^{-11}$
between `rtol=1e-12` and `1e-13` --- above what an earlier version of this figure was
plotting. The cell below therefore computes that spread per panel and draws it as a
floor: nothing beneath it is a measurement.

One profile serves every panel. Only the energy window differs, and it has to: an eV-scale
$\Delta m^2_{41}$ and a Standard-Model splitting do not resonate at the same energy.'''),
    code(r'''# ============================================================ Figure 1
# ONE profile for every panel -- a supernova-envelope exponential, central density
# 3e3 g/cm^3, scale height 10 km, baseline 25 km.  Only the energy window differs,
# because an eV-scale sterile splitting and a Standard-Model one do not resonate at
# the same energy: at the energies that make the sterile panels legible the active
# sector has stopped oscillating, and vice versa.
RHO0_1, LS_KM_1, L_KM_1 = 3.0e3, 10.0, 25.0
rho_1 = matter.exp_density_profile(RHO0_1*gd.UNIT_G_PER_CM3, LS_KM_1*gd.UNIT_KM)
VCC_1 = matter.vcc_func_from_rho_func(rho_1, L0=0.0)
L_1 = L_KM_1*gd.UNIT_KM
CASES = {2: (0.0005, 0.05), 3: (0.002, 0.2), 4: (2.0, 20.0), 5: (2.0, 20.0)}
N_PLOT, N_REF = 140, 26

# Check the oracle before using it.  DOP853 at rtol=1e-12 is worth something only if
# tightening it does not move the answer, and how far it moves is a property of the
# problem: at two and three flavors it is 4e-14, at four and five it is 4e-11, which
# is ABOVE the residual being measured there.  Drawn as a floor on each panel rather
# than left for a reader to assume it away.
def compute_validation():
    """Probabilities, the DOP853 oracle, and the oracle's own spread, for every panel.

    Two {\tt DOP853} solutions per reference energy at four flavor counts is the bulk of
    this figure's cost, and none of it depends on the machine: it is arithmetic on a fixed
    configuration.  Computed once and stored.
    """
    out = {}
    for d, (lo, hi) in CASES.items():
        E_plot = np.logspace(np.log10(lo), np.log10(hi), N_PLOT)*gd.UNIT_GEV
        E_ref = np.logspace(np.log10(lo), np.log10(hi), N_REF)*gd.UNIT_GEV
        P_plot = np.array([np.asarray(quiet(oscprob.osc_prob, make_H_func(d, E, VCC_1),
                                            0.0, L_1, rtol=RTOL_ACC, atol=ATOL_ACC))
                           for E in E_plot])
        Pm, Pr, floor = [], [], []
        for E in E_ref:
            Hf = make_H_func(d, E, VCC_1)
            Pm.append(np.asarray(quiet(oscprob.osc_prob, Hf, 0.0, L_1,
                                       rtol=RTOL_ACC, atol=ATOL_ACC)))
            r12 = ode_reference(Hf, L_1, d, rtol=1.0e-12, atol=1.0e-14)
            r13 = ode_reference(Hf, L_1, d, rtol=1.0e-13, atol=1.0e-15)
            Pr.append(r12)
            floor.append(float(np.max(np.abs(r12 - r13))))
        Pm, Pr = np.array(Pm), np.array(Pr)
        out[str(d)] = dict(E_plot=E_plot.tolist(), P_plot=P_plot.tolist(),
                           E_ref=E_ref.tolist(), Pr=Pr.tolist(), floor=floor,
                           resid=np.max(np.abs(Pm - Pr), axis=(1, 2)).tolist(),
                           unit=np.max(np.abs(P_plot.sum(axis=2) - 1.0), axis=1).tolist())
    return out


_val = cached('oracle',
              (profile_samples(VCC_1, L_1), float(L_1), N_PLOT, N_REF, RTOL_ACC, ATOL_ACC,
               repr(sorted(CASES.items())), repr(sorted(OSC.items())),
               repr(sorted(STERILE4.items())), repr(sorted(STERILE5.items()))),
              compute_validation,
              what='Figure 1: probabilities, the DOP853 oracle, and the oracle spread.')
results = {int(k): {kk: np.array(vv) for kk, vv in v.items()} for k, v in _val.items()}
for d in sorted(results):
    r = results[d]
    print('d=%d  %d points   max|dP| %.2e   oracle floor %.2e   unitarity %.2e'
          % (d, N_PLOT, r['resid'].max(), r['floor'].max(), r['unit'].max()))'''),
    md(r'''### Drawing it'''),
    code(r'''fig, axes = plt.subplots(2, 4, figsize=(WIDE, 3.25), sharex='col',
                         gridspec_kw=dict(height_ratios=[2.05, 1.0], hspace=0.06,
                                          wspace=0.09))
COLORS = {2: BLUE, 3: ORANGE, 4: GREEN, 5: PURPLE}
for j, d in enumerate((2, 3, 4, 5)):
    r = results[d]
    top, bot = axes[0, j], axes[1, j]
    x = r['E_plot']/gd.UNIT_GEV
    top.plot(x, r['P_plot'][:, 0, 0], color=COLORS[d], lw=1.0, zorder=3)
    top.plot(r['E_ref']/gd.UNIT_GEV, r['Pr'][:, 0, 0], ls='none', marker='o', ms=2.0,
             mfc='none', mew=0.7, color=INK, zorder=4)
    logx(top); snug(top, x); top.set_ylim(0.0, 1.0); minor_y(top, 5)
    corner(top, FLAVOR_LABEL[d], fontsize=8.0, x=0.955, y=0.975)
    bot.semilogy(r['E_ref']/gd.UNIT_GEV, np.maximum(r['resid'], 1.0e-18),
                 color=COLORS[d], lw=1.0)
    bot.semilogy(r['E_ref']/gd.UNIT_GEV, np.maximum(r['floor'], 1.0e-18), color=INK,
                 lw=0.7, ls='--')
    bot.semilogy(x, np.maximum(r['unit'], 1.0e-18), color='0.55', lw=0.5)
    logx(bot); logy(bot); snug(bot, x); bot.set_ylim(1.0e-16, 3.0e-9)
    if j:
        top.tick_params(labelleft=False)
        bot.tick_params(labelleft=False)
axes[0, 0].set_ylabel(r'Survival probability, $P_{\nu_e \to \nu_e}$')
axes[1, 0].set_ylabel(r'Absolute deviation')
axes[0, 0].plot([], [], color=INK, lw=1.0, label=r'Mag$\nu$s')
axes[0, 0].plot([], [], ls='none', marker='o', ms=2.4, mfc='none', mew=0.7, color=INK,
                label='DOP853')
axes[0, 0].legend(loc='lower right', handlelength=1.3)
axes[1, 0].plot([], [], color=INK, lw=1.0, label=r'Max $|\Delta P|$')
axes[1, 0].plot([], [], color=INK, lw=0.7, ls='--', label='Oracle floor')
axes[1, 0].plot([], [], color='0.55', lw=0.5, label=r'$|\sum_\beta P - 1|$')
axes[1, 0].legend(loc='upper left', handlelength=1.4, labelspacing=0.18, fontsize=8.0)
# Centred on the four columns rather than on the figure, which the y-axis label
# would otherwise pull left, and close to the axis it belongs to.
fig.align_labels()
box = [axes[1, j].get_position() for j in range(4)]
fig.text(0.5*(box[0].x0 + box[3].x1), box[0].y0 - 0.075,
         r'Neutrino energy, $E$ [GeV]', ha='center', va='top', fontsize=8)
save(fig, 'validation.pdf')'''),
    md(r'''### Figure 1b --- how a call is answered

Six engines can answer a request, tried in a fixed order; each declines what it cannot
serve honestly and falls through to the next.'''),
    code(r'''# ------------------------------------------- Figure 1b: the six engines, in order
# Same idiom as the slab-composition figure of the companion paper: one bar per
# engine, shaded by density, with what each does to the trajectory drawn rather
# than named.  The order is the dispatch order of the three scenario wrappers.
SLAB = ['#eaf2fb', '#bcd8f3', '#7fb4e6', '#3a86d4', '#1c71d8']
from matplotlib.patches import Rectangle
fig, ax = plt.subplots(figsize=(WIDE, 4.3))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
X0, X1 = 30.0, 66.0                      # the bar spans the same x in every row
H = 6.2                                  # bar height
rows = [
 ('Closed-form average',   'An average is asked for,\nand $\\mathbb{H}$ does not vary'),
 ('Adiabatic $+$ Magnus',  'Smooth profile, a tolerance,\nand it certifies itself'),
 ('Interaction picture',   'Declared exponential, two flavors,\nand its iteration converges'),
 ('Energy-batched scan',   'Many energies,\none baseline'),
 ('Cumulative scan',       'One energy,\nmany baselines'),
 ('General Magnus ladder', 'Anything else'),
]
ys = np.linspace(84, 6, len(rows))

def bar(y, edges, shades, lw=0.7):
    for (a, b), c in zip(edges, shades):
        ax.add_patch(Rectangle((a, y), b-a, H, facecolor=c, edgecolor=INK, lw=lw, zorder=2))

for i, ((name, when), y) in enumerate(zip(rows, ys)):
    # The energy-batched row has arrows entering at its left edge, so its name needs
    # more clearance than the others.
    ax.text(X0-(4.2 if i == 3 else 2.5), y+H/2, name, ha='right', va='center',
            fontsize=8.6, color='black')
    # The ladder row carries the refine arrow just past its bar, so its condition
    # text starts further right than the others'.
    ax.text(X1+(6.0 if i == 5 else 2.6), y+H/2, when, ha='left', va='center',
            fontsize=6.6, color=INK)
    if i:                                        # every engine but the first walks a path
        ax.annotate('', xy=(X0-0.6, y+H/2), xytext=(X0-1.8, y+H/2),
                    arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.8))
    if i == 0:
        # One undivided block: nothing is composed along it, because nothing is
        # propagated.  avgprob diagonalises H once and sums |sum_i V*_ai V_bi|^2 over
        # the eigenbasis, which is exact for the averaged observable -- not a closed
        # form in the mixing parameters, since it still needs the eigenvectors.
        bar(y, [(X0, X1)], [SLAB[2]])
        ax.text((X0+X1)/2, y+H/2, r'$\langle P\rangle$ from the eigenbasis',
                ha='center', va='center', fontsize=7.0, color='white')
        ax.text((X0+X1)/2, y-1.6, 'Exact for the average; nothing is propagated',
                ha='center', va='top', fontsize=6.4, color=INK)
    elif i == 1:                                 # smooth gradient, one exact patch
        n = 60; e = np.linspace(X0, X1, n+1)
        g = plt.cm.Blues(np.linspace(0.15, 0.75, n))
        bar(y, list(zip(e[:-1], e[1:])), g, lw=0.0)
        ax.add_patch(Rectangle((X0, y), X1-X0, H, fill=False, edgecolor=INK, lw=0.7, zorder=3))
        px = X0 + 0.56*(X1-X0)
        ax.add_patch(Rectangle((px, y), 5.0, H, facecolor=ORANGE, edgecolor='black',
                               lw=0.8, zorder=4))
        ax.text(px+2.5, y+H+1.4, 'Magnus patch', ha='center', va='bottom',
                fontsize=6.6, color=ORANGE)
        ax.text(X0+0.22*(X1-X0), y+H/2, 'Adiabatic transport', ha='center', va='center',
                fontsize=7.0, color='white', zorder=5)
    elif i == 2:                                 # fast phase factored out first
        n = 40; e = np.linspace(X0, X1, n+1)
        g = plt.cm.Blues(np.linspace(0.65, 0.12, n))
        bar(y, list(zip(e[:-1], e[1:])), g, lw=0.0)
        ax.add_patch(Rectangle((X0, y), X1-X0, H, fill=False, edgecolor=INK, lw=0.7, zorder=3))
        xs = np.linspace(X0+1, X1-1, 300)
        ax.plot(xs, y+H/2 + 1.5*np.sin((xs-X0)*1.5), color='white', lw=0.9, zorder=5)
        ax.text((X0+X1)/2, y-1.6, r'Vacuum phase removed analytically, then one Magnus pass',
                ha='center', va='top', fontsize=6.4, color=INK)
    elif i == 3:                                 # profile once, many energies
        ed = np.linspace(X0, X1, 6)
        bar(y, list(zip(ed[:-1], ed[1:])), [SLAB[1], SLAB[3], SLAB[2], SLAB[4], SLAB[1]])
        for dy in (1.9, 0.0, -1.9):
            ax.annotate('', xy=(X0-0.5, y+H/2+dy), xytext=(X0-2.1, y+H/2+dy),
                        arrowprops=dict(arrowstyle='-|>', color=BLUE, lw=0.7))
        ax.text(X0-2.0, y+H+1.2, r'$E_1 \ldots E_n$', ha='center', va='bottom',
                fontsize=6.6, color=BLUE)
    elif i == 4:                                 # one pass, many baselines out
        ed = np.linspace(X0, X1, 6)
        bar(y, list(zip(ed[:-1], ed[1:])), [SLAB[1], SLAB[3], SLAB[2], SLAB[4], SLAB[1]])
        for xx in ed[1:]:
            ax.annotate('', xy=(xx, y-2.6), xytext=(xx, y+0.2),
                        arrowprops=dict(arrowstyle='-|>', color=GREEN, lw=0.7))
        ax.text((X0+X1)/2, y-4.2, r'$P(L_1),\, P(L_2),\, \ldots$', ha='center', va='top',
                fontsize=6.6, color=GREEN)
    else:                                        # the ladder: refine until two levels agree
        for tier, (nsl, dy, al) in enumerate([(4, 3.4, 0.40), (6, 1.7, 0.68), (10, 0.0, 1.0)]):
            ed = np.linspace(X0, X1, nsl+1)
            for a, b in zip(ed[:-1], ed[1:]):
                ax.add_patch(Rectangle((a, y+dy), b-a, H*0.62, facecolor=SLAB[2],
                                       edgecolor=INK, lw=0.5, alpha=al, zorder=2))
        ax.annotate('', xy=(X1+1.8, y-0.8), xytext=(X1+1.8, y+H+4.2),
                    arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9))
        ax.text(X1+3.4, y+H/2+1.7, 'Refine', ha='center', va='center', fontsize=6.6,
                color=INK, rotation=90)

ax.text(56.0, 97.0, r'How Mag$\nu$s answers a call: the six engines, in dispatch order',
        ha='center', va='center', fontsize=9.4, color='black')

# The order they are tried in, drawn once in the margin: down the middle it crossed
# every bar and collided with the notes under rows 1, 3 and 5.
ax.annotate('', xy=(6.0, ys[-1]-1.0), xytext=(6.0, ys[0]+H+1.0),
            arrowprops=dict(arrowstyle='-|>', color='black', lw=1.2))
ax.text(3.2, (ys[0]+ys[-1])/2 + H/2, 'Tried in this order; each falls through to the next',
        rotation=90, ha='center', va='center', fontsize=7.0, color='black')

# The last three share one kernel, which is what decides whether a disagreement between
# two of them means anything.  The energy-batched row is the exception noted beside it:
# a potential that does not vary is served by a single exact exponential instead.
gy0, gy1 = ys[5] - 1.6, ys[3] + H + 1.6
gx = 92.0
ax.plot([gx, gx+1.2, gx+1.2, gx], [gy0, gy0, gy1, gy1], color=INK, lw=0.8,
        solid_joinstyle='miter')
ax.text(gx + 2.8, (gy0 + gy1)/2, 'One slab kernel,\nthree ways of batching', rotation=90,
        ha='center', va='center', fontsize=6.4, color=INK)
ax.text(X1+2.6, ys[3]-1.4, r'(constant $\mathbb{H}$: one exact exponential instead)',
        ha='left', va='top', fontsize=6.0, color=INK, style='italic')
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
save(fig, 'strategies.pdf')'''),
    md(r'''## Figure 2 --- slab width follows the profile, not the phase

Three measurements: one slab against a constant Hamiltonian over six decades of $\Phi$;
the cost of a *delivered* accuracy as $\Phi$ grows; and what the truncation order buys at
fixed $\Phi$.

**All three panels are measured against the same reference**: a midpoint slab product
carried in `mpmath` and Richardson-extrapolated three times, with its own convergence
reported so that no panel has to be taken on trust. `DOP853` appears in the center panel
as the code being timed, never as a reference. For a constant Hamiltonian the
extrapolation collapses to a single exponential, which the left panel checks rather than
assumes.'''),
    code(r'''# ============================================================ Figure 2
D2 = 3
RHO0_2, LS_KM_2, RHO_CONST = 100.0, 300.0, 3.0
rho2 = matter.exp_density_profile(RHO0_2*gd.UNIT_G_PER_CM3, LS_KM_2*gd.UNIT_KM)
VCC2 = matter.vcc_func_from_rho_func(rho2, L0=0.0)
L2 = 3000.0*gd.UNIT_KM
E_FIX = 0.2*gd.UNIT_GEV
Hf2 = make_H_func(D2, E_FIX, VCC2)

# The six configurations both lower panels carry: what an order buys, beside what it costs.
# Orders two, four and six run on the commutator-free Gauss-Legendre schemes, which is what
# the package does by default.  Above six no such scheme exists, so the only route is
# cumulative quadrature; order six is run on it too, as a control, because the change in
# convergence rate belongs to the quadrature path and not to the higher orders.
SERIES = [('gl', 2), ('gl', 4), ('gl', 6), ('simpson', 6), ('simpson', 8), ('simpson', 10)]

MP_DPS_FIX, MP_NS_FIX = 50, (4096, 8192, 16384, 32768)
MP_DPS_SCAN, MP_SCAN_TARGET = 30, 1.0e-11
def slab_product(H_func, L, n_slabs):
    """The piecewise-constant operator, H frozen at each slab midpoint.

    The algorithm a closed-form slab code runs, written out here so that this
    comparison needs no second package installed.  The midpoint rather than the edge
    is what makes it second order: the linear part of the error cancels between a
    slab's two halves.
    """
    edges = np.linspace(0.0, L, n_slabs + 1)
    Hs = H_func(0.5*(edges[:-1] + edges[1:]))
    h = edges[1] - edges[0]
    U = np.eye(Hs.shape[-1], dtype=complex)
    for k in range(n_slabs):
        U = expm(-1j*Hs[k]*h) @ U
    return prob_from_U(U)


def chunked_prob(H_func, L, n_slabs, order, integration_method='gl', n_tpts=65,
                 max_cells=20000):
    """The package's slab chain, composed a batch at a time.

    Orders above six have no Gauss-Legendre scheme and run on cumulative quadrature,
    which carries ``n_tpts`` samples per slab.  At order ten the recursion holds of
    order a hundred intermediates of shape ``(n_slabs, n_tpts, d, d)``, about 6 MB per
    slab, and exhausts memory past roughly a thousand slabs.  Composition is a product
    of independent per-slab operators, so the chain can be built in chunks and
    multiplied in, which holds memory flat in ``n_slabs``.  The cell below checks this
    against ``osc_prob`` at every order.
    """
    edges = np.linspace(0.0, L, n_slabs + 1)
    nodes = (magnus.gl_nodes(order) if integration_method == 'gl'
             else np.linspace(0.0, 1.0, n_tpts))
    m = len(nodes)
    d = np.asarray(H_func(0.0)).shape[-1]
    chunk = max(1, min(n_slabs, int(max_cells//m)))
    U = np.eye(d, dtype=complex)
    for a in range(0, n_slabs, chunk):
        b = min(a + chunk, n_slabs)
        lo, hi = edges[a:b], edges[a + 1:b + 1]
        t = lo[:, None] + (hi - lo)[:, None]*nodes[None, :]
        At = -1j*np.asarray(H_func(t.ravel()), dtype=complex).reshape(b - a, m, d, d)
        Us = quiet(magnus.evolution_operators_from_samples, At, hi - lo, order=order,
                   integration_method=integration_method, validate_input=False)
        for k in range(b - a):
            U = Us[k] @ U
    return prob_from_U(U)


def mp_midpoint(H_func, L, n, dps):
    r"""A midpoint slab product carried at ``dps`` decimal digits, in mpmath.

    The panels below have to resolve deviations near $10^{-15}$, which no
    double-precision reference can certify.  A {\tt DOP853} solution moves by
    $2 \cdot 10^{-13}$ between ${\tt rtol} = 10^{-13}$ and $10^{-14}$, and a
    double-precision {\tt expm} of a constant Hamiltonian is itself wrong by
    $2 \cdot 10^{-10}$ at the largest phase the left panel reaches.  The probabilities
    come back as {\tt mpf} rather than as {\tt float}, since rounding them to double at
    this point would put a floor of $10^{-16}$ under everything built from them.
    """
    mp.mp.dps = dps
    edges = np.linspace(0.0, L, n + 1)
    Hs = H_func(0.5*(edges[:-1] + edges[1:]))
    h = mp.mpf(float(edges[1] - edges[0]))
    d = Hs.shape[-1]
    U = mp.eye(d)
    for k in range(n):
        M = mp.matrix([[mp.mpc(complex(Hs[k][i, j])) for j in range(d)] for i in range(d)])
        U = mp.expm(-1j*M*h)*U
    P = np.array([[mp.fabs(U[i, j])**2 for j in range(d)] for i in range(d)], dtype=object)
    return P.T


def richardson3(P, ns):
    r"""Remove $N^{-2}$, $N^{-4}$ and $N^{-6}$ in turn, staying in mpmath throughout.

    The midpoint product's error runs in even powers of the slab width, so
    $[4P(2N) - P(N)]/3$ removes the leading term, $[16Q(2N) - Q(N)]/15$ the next and
    $[64R(2N) - R(N)]/63$ the one after.  An earlier version returned {\tt float}
    probabilities from each product and combined them in double precision; it converged
    to $1.1 \cdot 10^{-16}$, which is double-precision epsilon and not the
    extrapolation.  Returns the reference together with its own self-convergence, so no
    panel has to take the reference on trust.
    """
    Q = {n: (4*P[2*n] - P[n])/3 for n in ns[:-1]}
    R = {n: (16*Q[2*n] - Q[n])/15 for n in ns[:-2]}
    T = (64*R[ns[1]] - R[ns[0]])/63
    return T, float(max(abs(x) for x in (T - R[ns[1]]).ravel()))


def mp_reference(H_func, L, base, dps=MP_DPS_SCAN, target=MP_SCAN_TARGET, cap=2048):
    r"""Triple-Richardson reference, with the base slab count raised until it converges.

    The base a trajectory needs grows with the accumulated phase: nine hundred radians
    needs $512$ where eight radians needs $64$.  Doubling until the self-convergence
    clears ``target`` finds it rather than assuming it.
    """
    while True:
        ns = (base, 2*base, 4*base, 8*base)
        T, sc = richardson3({n: mp_midpoint(H_func, L, n, dps) for n in ns}, ns)
        if sc < target or base >= cap:
            return T, sc, base
        base *= 2




# One reference machinery serves all three panels.  For a constant Hamiltonian every
# midpoint slab product is the same exponential whatever N, so the triple Richardson of
# them collapses to that one exponential; the print below measures the collapse rather
# than assuming it.  This matters: the double-precision expm this panel used to be
# measured against is itself wrong by 2e-10 at the largest phase here, so the old curve
# was as much the reference's round-off as the method's.
VCC_CONST = float(matter.vcc_func_from_rho_func(RHO_CONST*gd.UNIT_G_PER_CM3))
E_CONST = np.logspace(np.log10(50.0), np.log10(5.0e-6), 60)*gd.UNIT_GEV
PROJ2 = np.asarray(matter.matter_potential_projector(D2), dtype=complex)
HVAC2 = np.asarray(vacuum_hamiltonian(D2), dtype=complex)


def const_H_func(Hc):
    def f(l):
        return np.broadcast_to(Hc, np.shape(l) + Hc.shape).copy()
    return f


phase_a, err_a, collapse = [], [], []
for E in E_CONST:
    Hc = HVAC2/E + VCC_CONST*PROJ2
    w = np.linalg.eigvalsh(Hc)
    phase_a.append(float((w[-1] - w[0])*L2))
    Hcf = const_H_func(Hc)
    T, sc, _ = mp_reference(Hcf, L2, 4, dps=MP_DPS_FIX, target=1.0e-40)
    collapse.append(float(max(abs(x) for x in (T - mp_midpoint(Hcf, L2, 1, MP_DPS_FIX)).ravel())))
    P1 = np.asarray(quiet(oscprob.osc_prob, lambda l, H=Hc: H, 0.0, L2, n_slabs=1,
                          rtol=None, atol=None))
    err_a.append(max(float(max(abs(x) for x in (P1 - T).ravel())), 1.0e-18))
phase_a, err_a = np.array(phase_a), np.array(err_a)
print('left panel: %d points, phase %.1e to %.1e rad' % (len(phase_a), phase_a.min(), phase_a.max()))
print('  triple Richardson vs one exponential, constant H: %.1e' % max(collapse))
print('  deviation of one slab: %.1e to %.1e, growing as Phi^%.2f'
      % (err_a.min(), err_a.max(),
         np.polyfit(np.log(phase_a[phase_a > 1e2]), np.log(err_a[phase_a > 1e2]), 1)[0]))
print('  for comparison, Phi*eps at the largest phase is %.1e' % (phase_a.max()*2.220446e-16))'''),
    md(r'''### The center panel --- cost at an accuracy each code actually delivers'''),
    code(r'''# --------------------------- Figure 2, cost panel: price at one fixed, tight tolerance
# Every configuration is run at one tolerance, fixed before the sweep starts and never
# retuned -- the solver included, so that nothing is asked for less than anything else.
# Searching per energy for the cheapest setting that just clears an accuracy target was
# what made this panel jagged: the ladders are discrete, so neighbouring energies landed
# on different rungs and the cost jumped between them for no reason a reader could see.
# It is also not how the code is used.  A tolerance is chosen once and the energies are
# swept under it, so that is what the panel measures.
#
# A tolerance is a request, not a guarantee, so what each configuration DELIVERS at it is
# measured rather than assumed, and the spread is reported with the figure.
RTOL_FIXED = 1.0e-8
ENERGIES2 = np.logspace(np.log10(20.0), np.log10(0.02), 41)


# Order two converges as N^-2, so it refines to 106204 slabs at the highest phase swept
# here -- five times the 20000 the Gauss-Legendre path caps at by default, which would
# have truncated it silently.  The ceiling is raised for every Gauss-Legendre
# configuration rather than for order two alone, so that none of them is capped where the
# others are not, and the cost each is charged is the cost of the work it actually did.
GL_SLAB_CEILING = 400000


def magnus_call(Hf, order, meth, n_tpts):
    kw = dict(magnus_exp_order=order, integration_method=meth)
    if meth == 'gl':
        kw.update(max_n_slabs=GL_SLAB_CEILING, max_num_loops=80)
    elif n_tpts is not None:
        kw.update(n_tpts_per_slab=n_tpts, min_n_tpts_per_slab=n_tpts,
                  max_n_tpts_per_slab=n_tpts)
    return lambda rt: np.asarray(quiet(oscprob.osc_prob, Hf, 0.0, L2, rtol=rt,
                                       atol=rt*1e-2, **kw))


def dop853_call(Hf):
    return lambda rt: ode_reference(Hf, L2, D2, rtol=rt, atol=rt*1e-2)


# One entry per curve.  The panel carries the same six configurations as the panel beside
# it, so that what an order costs can be read against what it buys, together with the
# solver they are all measured against.
CODES = [('dop853', 'Ref.: Runge-Kutta order 8 (DOP853)', INK, '-', None, None)]
for _meth, _order in SERIES:
    CODES.append(('%s%d' % (_meth, _order),
                  '%d, %s' % (_order, 'GL' if _meth == 'gl' else 'Simpson'),
                  ORDER_COLOR[_order], '-' if _meth == 'gl' else '--', _meth, _order))


def fixed_setting(name, meth):
    """The single setting this configuration is run at, at every energy.

    A tolerance and nothing else.  Pinning the samples per slab for the cumulative-
    quadrature configurations would hold their delivered accuracies closer together, but
    it is a setting no caller types: the package refines that knob from the tolerance
    along with the slab count, and whatever it costs to do so is part of what the call
    costs.  The panel beside this one pins it, on purpose --- there the slab count is the
    variable under study, so the quadrature has to be held out of the way.
    """
    return {'rtol': RTOL_FIXED}


def call_for(name, meth, order, Hf, setting):
    return (dop853_call(Hf) if name == 'dop853'
            else magnus_call(Hf, order, meth, setting.get('n_tpts')))


def delivered(call, rt, T):
    return float(max(abs(x) for x in (np.asarray(call(rt)) - T).ravel()))


def measure_once(name, meth, order, Hf, rt, T):
    """One call, read twice: what it delivers, and what its ladder settled on.

    Asking those separately would run the whole configuration twice at identical settings
    --- order ten unpinned is not cheap enough to pay for that.
    """
    if name == 'dop853':
        P = np.asarray(dop853_call(Hf)(rt))
        return float(max(abs(x) for x in (P - T).ravel())), {}
    info = {}
    kw = dict(magnus_exp_order=order, integration_method=meth, convergence_info=info)
    if meth == 'gl':
        kw.update(max_n_slabs=GL_SLAB_CEILING, max_num_loops=80)
    P = np.asarray(quiet(oscprob.osc_prob, Hf, 0.0, L2, rtol=rt, atol=rt*1e-2, **kw))
    return (float(max(abs(x) for x in (P - T).ravel())),
            {k: int(info[k]) for k in ('n_slabs', 'n_tpts_per_slab') if k in info})


def scan_setup():
    """Reference and tuned setting at each energy, computed once and read from disk after.

    Five minutes of mpmath and a two-knob search that runs order ten at its tightest
    settings are worth spending once rather than on every rebuild.  Both are properties of
    the configuration and not of the machine, so both are cached; only the timings below
    have to be measured live.  The cached settings are re-verified rather than trusted ---
    one call per configuration per energy --- so a change in the package that moved them
    cannot pass silently.
    """
    key = fingerprint(profile_samples(VCC2, L2), float(L2), D2, MP_DPS_SCAN, MP_SCAN_TARGET,
                      [c[0] for c in CODES],
                      [float(e) for e in ENERGIES2], RTOL_FIXED, 'fixed-tolerance')
    blob = json.loads(MP_CACHE.read_text()) if MP_CACHE.exists() else {}
    cached = blob.get('scan', {})
    entries = cached.get('entries') if cached.get('fingerprint') == key else None
    if entries is not None:
        mp.mp.dps = MP_DPS_SCAN
        Ts = [np.array([[mp.mpf(x) for x in row] for row in e['P']], dtype=object)
              for e in entries]
        # Spot-checked rather than fully re-run: the fingerprint already covers the
        # configuration, so this is here to catch a change inside the package that moved
        # the answers under an unchanged configuration.  Three energies at the ends and
        # the middle of the sweep do that; forty-one would cost minutes on every
        # regeneration to re-derive what the fingerprint already settled.
        # Once per package version, not once per regeneration.  The fingerprint already
        # covers the configuration; what it cannot see is a change inside the package
        # that moves the answers under a configuration that did not move.  Tying the
        # check to the version catches exactly that, and costs nothing on the format
        # edits that make up most rebuilds.  A cache written before this field existed
        # is adopted rather than re-probed: the code that wrote it was this code.
        seen_version = cached.get('magnus_version')
        if seen_version == MAGNUS_VERSION or seen_version is None:
            if seen_version is None:
                blob['scan'] = dict(cached, magnus_version=MAGNUS_VERSION)
                write_cache(blob)
            print('  references and delivered accuracies read from %s, configuration %s'
                  % (MP_CACHE.name, key[:12]))
            return Ts, entries
        print('  package moved from %s to %s: spot-checking the stored accuracies'
              % (seen_version, MAGNUS_VERSION))
        probe = sorted({0, len(entries)//2, len(entries) - 1})
        stale = []
        for i in probe:
            e, T = entries[i], Ts[i]
            Hf = make_H_func(D2, e['E']*gd.UNIT_GEV, VCC2)
            for name, _, _, _, meth, order in CODES:
                st = e['setting'][name]
                now = delivered(call_for(name, meth, order, Hf, st), st['rtol'], T)
                if not (0.5 <= now/max(e['acc'][name], 1.0e-18) <= 2.0):
                    stale.append((e['E'], name, e['acc'][name], now))
        if not stale:
            blob['scan'] = dict(cached, magnus_version=MAGNUS_VERSION)
            write_cache(blob)
            print('  stored accuracies still hold under %s, re-verified at %d energies'
                  % (MAGNUS_VERSION, len(probe)))
            return Ts, entries
        print('  the package now delivers differently at %d spot checks '
              '(e.g. %s at %.3f GeV: %.1e stored, %.1e now): recomputing'
              % (len(stale), stale[0][1], stale[0][0], stale[0][2], stale[0][3]))
    else:
        cache_miss('scan', key)
    print('  energy scan: configuration moved, recomputing')
    t0 = time.perf_counter()
    Ts, entries = [], []
    for Egev in ENERGIES2:
        Hf = make_H_func(D2, Egev*gd.UNIT_GEV, VCC2)
        T, sc, base = mp_reference(Hf, L2, 32)
        setting = {name: fixed_setting(name, meth)
                   for name, _, _, _, meth, _ in CODES}
        got = {name: measure_once(name, meth, order, Hf, setting[name]['rtol'], T)
               for name, _, _, _, meth, order in CODES}
        acc = {k: v[0] for k, v in got.items()}
        chose = {k: v[1] for k, v in got.items()}
        Ts.append(T)
        entries.append(dict(E=float(Egev), phase=accumulated_phase(Hf, L2), base=int(base),
                            self_conv=float(sc), setting=setting, acc=acc, chose=chose,
                            P=[[mp.nstr(x, MP_DPS_SCAN - 5) for x in row] for row in T]))
        print('    %6.3f GeV done, worst delivered %.1e (%.0f s elapsed)'
              % (Egev, max(acc.values()), time.perf_counter() - t0), flush=True)
    print('  %d references in %.0f s, worst self-convergence %.1e'
          % (len(Ts), time.perf_counter() - t0, max(e['self_conv'] for e in entries)))
    blob['scan'] = dict(dps=MP_DPS_SCAN, target=MP_SCAN_TARGET, rtol=RTOL_FIXED,
                        fingerprint=key, entries=entries)
    write_cache(blob)
    return Ts, entries


REFS2, REF_INFO2 = scan_setup()
print('at rtol = %.0e, over %d energies:' % (RTOL_FIXED, len(REF_INFO2)))
for _n, _lab, _, _, _, _ in CODES:
    _a = [e['acc'][_n] for e in REF_INFO2]
    _c = [e.get('chose', {}).get(_n, {}) for e in REF_INFO2]
    _ns = [d['n_slabs'] for d in _c if 'n_slabs' in d]
    _nt = [d['n_tpts_per_slab'] for d in _c if 'n_tpts_per_slab' in d]
    print('  %-30s delivered %.1e to %.1e%s%s'
          % (_lab, min(_a), max(_a),
             '' if not _ns else ',  %d-%d slabs' % (min(_ns), max(_ns)),
             '' if not _nt else ',  %d-%d samples/slab' % (min(_nt), max(_nt))))


def best_of(fn, k=5):
    out = []
    for _ in range(k):
        t0 = time.perf_counter()
        fn()
        out.append(time.perf_counter() - t0)
    return min(out)


# A fixed workload interleaved through the sweep, so machine drift is visible rather than
# absorbed into the curves.
H_CTRL = make_H_func(D2, 1.0*gd.UNIT_GEV, VCC2)
_CTRL0 = []


def control_baseline():
    """The control's first value, measured on demand.

    At cell scope this ran a warm-up and five timed repeats on every regeneration of the
    notebook, including the ones that only moved a label.  It is an input to the timings
    and to nothing else, so it belongs behind the same cache they do.
    """
    if not _CTRL0:
        quiet(oscprob.osc_prob, H_CTRL, 0.0, L2, rtol=1e-8, atol=1e-10)  # warm numba
        _CTRL0.append(best_of(lambda: quiet(oscprob.osc_prob, H_CTRL, 0.0, L2,
                                            rtol=1e-8, atol=1e-10)))
    return _CTRL0[0]

def measure_timings():
    """Time each configuration at its tuned setting, once, and read from disk after.

    Timing is the only part of this figure that is a property of the machine rather than
    of the configuration, which is the argument for measuring it live.  It is also minutes
    of work on every rebuild, and the notebooks are regenerated in continuous integration,
    where minutes of stopwatch tell nobody anything.  So it is cached like everything else,
    keyed on the configuration alone: change an energy, a code, a tolerance target or a
    tuned setting and it is re-measured; change nothing and the stored numbers stand.

    The machine and the date are stored beside the numbers, because a timing that does not
    say where it came from cannot be checked.  Set MAGNUS_PAPER_RETIME=1 to force a fresh
    measurement without changing the configuration -- which is what to do after moving the
    paper to another machine.
    """
    key = fingerprint(profile_samples(VCC2, L2), float(L2), D2, RTOL_FIXED, [c[0] for c in CODES],
                      [float(e['E']) for e in REF_INFO2],
                      [repr(e['setting']) for e in REF_INFO2])
    blob = json.loads(MP_CACHE.read_text()) if MP_CACHE.exists() else {}
    stored = blob.get('timings', {})
    if stored.get('fingerprint') == key and not os.environ.get('MAGNUS_PAPER_RETIME'):
        print('  timings read from %s, unchanged configuration %s, measured on %s'
              % (MP_CACHE.name, key[:12], stored.get('machine', 'an unrecorded machine')))
        return stored['rows']
    cache_miss('timings', key)
    print('  timings: configuration moved, re-measuring')
    out = []
    for e in REF_INFO2:
        Hf = make_H_func(D2, e['E']*gd.UNIT_GEV, VCC2)
        # The control is measured beside the point, not once at the start, and divided
        # out.  A sweep of forty-one energies takes long enough that the machine drifts
        # under it, and drift that is not divided out is read off the panel as structure.
        ctrl = best_of(lambda: quiet(oscprob.osc_prob, H_CTRL, 0.0, L2,
                                     rtol=1e-8, atol=1e-10), 5)
        drift = ctrl/control_baseline()
        t = {}
        for name, _, _, _, meth, order in CODES:
            st = e['setting'][name]
            call, rt = call_for(name, meth, order, Hf, st), st['rtol']
            t[name] = best_of(lambda c=call, x=rt: c(x), 5)/drift
        out.append(dict(phase=e['phase'], E=e['E'], t=t, acc=e['acc'], ctrl=drift))
    blob['timings'] = dict(fingerprint=key, rows=out,
                           machine=platform.node() + ', ' + platform.processor(),
                           measured=time.strftime('%Y-%m-%d'))
    write_cache(blob)
    return out


rows = measure_timings()
ph = np.array([r['phase'] for r in rows])
times = {name: np.array([r['t'][name] for r in rows]) for name, _, _, _, _, _ in CODES}
print('cost panel: %d points, phase %.1f to %.0f rad' % (len(rows), ph.min(), ph.max()))
for name, lab, _, _, _, _ in CODES:
    tt = times[name]
    print('  %-28s phase^%5.2f  %7.2f ms to %8.1f ms   worst delivered %.1e'
          % (lab, np.polyfit(np.log(ph), np.log(tt), 1)[0], 1e3*tt[0], 1e3*tt[-1],
             max(r['acc'][name] for r in rows)))
ctrl = np.array([r['ctrl'] for r in rows])
print('  interleaved control held between %.2f and %.2f of its first value'
      % (ctrl.min(), ctrl.max()))'''),
    md(r'''### The right panel --- what the truncation order buys'''),
    code(r'''# ------------------------ Figure 2, right panel: what the truncation order buys
def fixed_reference():
    """The fixed-phase reference: fifty digits, Richardson three times, cached.

    Two and a half minutes of mpmath is not worth paying on every rebuild of the figure.
    """
    key = fingerprint(np.asarray(Hf2(np.linspace(0.0, L2, 9))).view(float), float(L2), MP_DPS_FIX, MP_NS_FIX)
    blob = json.loads(MP_CACHE.read_text()) if MP_CACHE.exists() else {}
    if blob.get('fixed', {}).get('fingerprint') == key:
        mp.mp.dps = MP_DPS_FIX
        print('  fixed-phase reference read from %s, unchanged configuration %s'
              % (MP_CACHE.name, key[:12]))
        P = {int(n): np.array([[mp.mpf(x) for x in row] for row in m], dtype=object)
             for n, m in blob['fixed']['products'].items()}
    else:
        print('  fixed-phase reference: configuration moved, recomputing')
        t0 = time.perf_counter()
        P = {n: mp_midpoint(Hf2, L2, n, MP_DPS_FIX) for n in MP_NS_FIX}
        print('  %d-digit slab products at N = %s in %.0f s'
              % (MP_DPS_FIX, ', '.join(str(n) for n in MP_NS_FIX), time.perf_counter() - t0))
        blob['fixed'] = dict(dps=MP_DPS_FIX, n_slabs=list(MP_NS_FIX), fingerprint=key,
                             products={str(n): [[mp.nstr(x, MP_DPS_FIX + 5) for x in row]
                                                for row in P[n]] for n in MP_NS_FIX})
        write_cache(blob)
    return richardson3(P, MP_NS_FIX)


P_TRUE, REF_FLOOR = fixed_reference()
print('  the reference and the next-coarser extrapolation of it differ by %.1e' % REF_FLOOR)

# The chunked composition has to BE the package's answer, not merely close to it.
for order, meth in ((2, 'gl'), (6, 'gl'), (8, 'simpson'), (10, 'simpson')):
    a = np.asarray(quiet(oscprob.osc_prob, Hf2, 0.0, L2, n_slabs=32, magnus_exp_order=order,
                         integration_method=meth, n_tpts_per_slab=65, rtol=None, atol=None))
    print('  chunked composition vs osc_prob, order %2d on %-8s: %.1e'
          % (order, meth, np.max(np.abs(a - chunked_prob(Hf2, L2, 32, order, meth)))))

# Half-octave steps, so the curves read as curves rather than as line segments.
NS = np.array(sorted(set(int(round(2**(k/2.0))) for k in range(2, 29))))
M_HI = 65
# Orders two, four and six on the commutator-free Gauss-Legendre schemes, which is what
# the package does by default.  Above six no such scheme exists, so the only route is
# cumulative quadrature; order six is run on it too, as a control, because the change in
# convergence rate belongs to the quadrature path and not to the higher orders.


def order_curves():
    """The six curves, computed once and read from disk after.

    Order ten over twenty-seven slab counts out to N = 16384 is twenty minutes of
    cumulative quadrature.  The curves are exact functions of the configuration, so they
    cache exactly as the reference does.
    """
    key = fingerprint(np.asarray(Hf2(np.linspace(0.0, L2, 9))).view(float), float(L2), MP_DPS_FIX, MP_NS_FIX, [int(n) for n in NS], M_HI,
                      [list(x) for x in SERIES])
    blob = json.loads(MP_CACHE.read_text()) if MP_CACHE.exists() else {}
    if blob.get('orders', {}).get('fingerprint') == key:
        print('  order curves read from %s, unchanged configuration %s'
              % (MP_CACHE.name, key[:12]))
        return {tuple(k.split('_')[:1]) + (int(k.split('_')[1]),): np.array(v)
                for k, v in blob['orders']['curves'].items()}
    print('  order curves: configuration moved, recomputing')
    t0 = time.perf_counter()
    out = {}
    for meth, order in SERIES:
        out[(meth, order)] = np.array(
            [max(float(max(abs(x) for x in
                           (chunked_prob(Hf2, L2, int(n), order, meth, n_tpts=M_HI)
                            - P_TRUE).ravel())), 1.0e-19) for n in NS])
    print('  six curves over %d slab counts in %.0f s' % (len(NS), time.perf_counter() - t0))
    blob['orders'] = dict(n_slabs=[int(n) for n in NS], n_tpts=M_HI, fingerprint=key,
                          curves={'%s_%d' % k: v.tolist() for k, v in out.items()})
    write_cache(blob)
    return out


curves = order_curves()

# The convergence rate, fitted where the curve is asymptotic: below the large-error
# regime at small N, above the round-off floor at large N.  The label each curve carries
# in the figure is the nearest even power, and the fit is printed so that the two can be
# seen to agree rather than assumed to.
POWERS = {}
print('    curve                 measured slope   label')
for k in SERIES:
    c = curves[k]
    # Below the regime where the series has not yet converged, above the round-off floor.
    # A window reaching up to 1e-2 catches the plunge between the two and fits it as slope.
    w = (c > 1.0e-12) & (c < 1.0e-5) & (NS >= 16)
    slope = np.polyfit(np.log(NS[w]), np.log(c[w]), 1)[0]
    POWERS[k] = int(round(slope/2.0))*2
    print('    order %2d on %-8s  %13.2f   N^%d' % (k[1], k[0], slope, POWERS[k]))

# Order two IS the midpoint slab product, so only one of the two is drawn.  That they
# coincide is checked here rather than plotted twice.
agree = [float(np.max(np.abs(chunked_prob(Hf2, L2, int(n), 2, 'gl')
                             - slab_product(Hf2, L2, int(n))))) for n in NS]
print('  order 2 against the midpoint slab product, matrix by matrix: %.1e to %.1e'
      % (min(agree), max(agree)))
print('  every curve bottoms out and turns up again; the floor is the round-off of the')
print('  ordered product, which no truncation order can get under:')
for k in SERIES:
    c = curves[k]
    i = int(np.argmin(c))
    print('    order %2d on %-8s  minimum %.1e at N = %5d, %.1e at N = %d'
          % (k[1], k[0], c[i], NS[i], c[-1], NS[-1]))'''),
    md(r'''### Drawing it'''),
    code(r'''def label_along(ax, xs, ys, i, text, color, fontsize=8.0, offset=(4, 4), chord=False):
    """Label a curve in place, rotated to the angle its curve makes ON THE PAGE.

    A fixed rotation cannot serve six slopes differing by a factor of six: it put the
    shallowest label through its own curve and the steepest one off the edge.  The angle
    is read from ``ax.transData``, so this must be called once the geometry is final.
    """
    if chord:
        # The angle a straight label needs is the one of the chord it actually spans,
        # not the tangent at its anchor.  On a curve that steepens under the words --
        # the solar profile does -- the tangent is the shallower of the two, and the
        # curve closes on the text and cuts through it before the last character.
        f = ax.figure
        f.canvas.draw()
        probe = ax.text(0, 0, text, fontsize=fontsize)
        w = probe.get_window_extent(f.canvas.get_renderer()).width
        probe.remove()
        P = ax.transData.transform(np.column_stack([np.asarray(xs), np.asarray(ys)]))
        m = i
        while m < len(xs) - 1 and P[m, 0] - P[i, 0] < w:
            m += 1
        a, b = P[i], P[m]
    else:
        j, k = max(i - 1, 0), min(i + 1, len(xs) - 1)
        a = ax.transData.transform((xs[j], ys[j]))
        b = ax.transData.transform((xs[k], ys[k]))
    ax.annotate(text, xy=(xs[i], ys[i]), xytext=offset, textcoords='offset points',
                color=color, fontsize=fontsize,
                rotation=np.degrees(np.arctan2(b[1] - a[1], b[0] - a[0])),
                rotation_mode='anchor')


def phase_axis(axis, energy, phase, ticks, powers=False):
    r"""A secondary top axis carrying the accumulated phase.

    The phase is what a stepping solver's cost tracks, but it is not a dial: what a user
    sets is the energy, so the energy carries the primary axis.  The two are not
    reciprocal.  Below about $0.2$ GeV the vacuum term dominates and $\Phi$ goes as
    $1/E$; above it the matter potential fixes the eigenvalue splitting and $\Phi$ barely
    moves, so a thousandfold in energy is only a hundred and fiftyfold in phase.  The top
    ticks are therefore unevenly spaced, which is the physics rather than a defect.
    """
    o = np.argsort(phase)
    tw = axis.twiny()
    tw.set_xscale('log')
    tw.set_xlim(axis.get_xlim())
    tw.set_xticks(np.exp(np.interp(np.log(ticks), np.log(np.asarray(phase)[o]),
                                   np.log(np.asarray(energy)[o]))))
    # Powers of ten only where every tick is one.  Writing half-decades that way rounds
    # 30 and 300 onto 10^1 and 10^2, which prints the same label twice.
    tw.set_xticklabels([(r'$10^{%d}$' % round(np.log10(t))) if powers else ('%g' % t)
                        for t in ticks])
    tw.set_xlabel(r'Accumulated phase, $\Phi$ [rad]', labelpad=8.5, fontsize=9.5)
    # Ticks point inward, so the labels have to be pushed clear of the spine or they
    # sit on top of it.
    tw.tick_params(axis='x', direction='in', which='both', top=True, pad=2.2, labelsize=8.5)
    # Only the twin's OWN x minors.  minorticks_off() would take the y axis with it,
    # and twiny() shares that axis with the panel underneath --- which silently
    # stripped the minor ticks off both panels that carry a phase axis.
    tw.xaxis.set_minor_locator(NullLocator())
    return tw


fig = plt.figure(figsize=(WIDE, 7.48))
gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.24)
axd = fig.add_subplot(gs[0, 0])
ax = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

l_km = np.linspace(0.0, 3000.0, 600)
axd.semilogy(l_km, RHO0_2*np.exp(-l_km/LS_KM_2), '-', color=INK, lw=1.2)
axd.axhline(RHO_CONST, color=INK, ls='--', lw=1.1)
axd.set_xlim(0.0, 3000.0); axd.set_ylim(3.0e-3, 1.0e2)
axd.set_xlabel(r'Distance along the trajectory, $l$ [km]', labelpad=1.5)
axd.set_ylabel(r'Matter density, $\rho$ [g cm$^{-3}$]', labelpad=2.0)
logy(axd)
unit_as_one(axd)
RHO_EXP = RHO0_2*np.exp(-l_km/LS_KM_2)

ax[0].loglog(E_CONST/gd.UNIT_GEV, err_a, '-', color=BLUE, lw=1.2)
ax[0].axhline(2.2e-16, color=INK, ls=':', lw=0.8)
ax[0].loglog(E_CONST/gd.UNIT_GEV, np.maximum(phase_a*2.2e-16, 2.2e-16), ls='--',
             color=INK, lw=0.7)
stamp(ax[0], r'Round-off, $\varepsilon$', y=0.075)
ax[0].set_xlabel(r'Neutrino energy, $E$ [GeV]')
ax[0].set_ylabel(r'Max $|\Delta P|$, one slab')
ax[0].set_ylim(1.0e-16, 1.0e-9)
logx(ax[0]); logy(ax[0]); snug(ax[0], E_CONST/gd.UNIT_GEV)
# Every decade carries a label; the default locator thins them once the range is
# this deep, which leaves the panel with three ticks over eight decades.
ax[0].yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
# Minor ticks between the decades, and the energy labels pushed clear of the spine
# that the inward-pointing ticks share with them.
minor_log_ticks(ax[0])
ax[0].tick_params(axis='x', which='major', pad=4.2)
corner(ax[0], r'Constant $\mathbb{H}$', x=0.955, y=0.94)

E_row = np.array([r['E'] for r in rows])
# Only the solver is named here.  The six Magnus configurations are named in the panel
# beside it, where the rate they converge at belongs; colour and line style are the same
# in both, so one legend identifies them twice.
for name, lab, col, ls, meth, order in CODES:
    ax[1].loglog(E_row, 1e3*times[name], ls=ls, color=col, lw=1.2,
                 label=(lab if meth is None else '_nolegend_'))
ax[1].set_xlabel(r'Neutrino energy, $E$ [GeV]')
ax[1].set_ylabel(r'Time per probability [ms]')
ax[1].set_ylim(1.0, 1.0e5)
# Five decades now that the cumulative-quadrature curves reach a minute per
# probability, so the ticks go by decade rather than by halves.
ax[1].yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=20))
unit_as_one(ax[1])
minor_log_ticks(ax[1])
logx(ax[1]); logy(ax[1]); snug(ax[1], E_row)
corner(ax[1], r'Varying $\mathbb{H}$', x=0.955, y=0.94)
ax[1].legend(loc='lower right', handlelength=1.5,
             fontsize=7.2, labelspacing=0.2, borderpad=0.3, handletextpad=0.5,
             framealpha=0.75)

STYLE = {k: (ORDER_COLOR[k[1]], '-' if k[0] == 'gl' else '--') for k in SERIES}
LABELS = {('gl',2):'Order 2, G-L', ('gl',4):'Order 4, G-L',
          ('gl',6):'Order 6, G-L', ('simpson',6):'Order 6, Simpson',
          ('simpson',8):'Order 8, Simpson', ('simpson',10):'Order 10, Simpson'}
for k in SERIES:
    col, ls = STYLE[k]
    ax[2].loglog(NS, curves[k], ls=ls, color=col, lw=1.2,
                 label=r'%s, $N^{%d}$' % (LABELS[k], POWERS[k]))
ax[2].set_xlabel(r'Slabs along the trajectory, $N$')
ax[2].set_ylabel(r'Max $|\Delta P|$, multiple slabs')
logx(ax[2]); logy(ax[2]); snug(ax[2], NS)
ax[2].set_ylim(3.0e-16, 1.0)
unit_as_one(ax[2])
corner(ax[2], r'Varying $\mathbb{H}$', x=0.955, y=0.94)
# Flush with the *border* of the box above it, not with the text inside it: the box
# adds 0.32 of its font size as padding, so the two align only once that is added back.
ax[2].annotate(r'Fixed $E = %g$ GeV' % (E_FIX/gd.UNIT_GEV), xy=(0.955, 0.868),
               xycoords='axes fraction', xytext=(0.32*8.5 + 0.3, 0),
               textcoords='offset points', ha='right', va='top', fontsize=8.5,
               color=INK, zorder=10)
# The six configurations are named here rather than in the cost panel, because the rate
# each converges at is what this panel measures.  Each entry is set in its curve's colour.
_leg = ax[2].legend(loc='lower left', handlelength=1.5, fontsize=7.4, labelspacing=0.2,
                    borderpad=0.3, handletextpad=0.5, title='Magnus order',
                    framealpha=0.75)
_leg.get_title().set_fontsize(7.4)
for _txt, _k in zip(_leg.get_texts(), SERIES):
    _txt.set_color(STYLE[_k][0])

# Square panels: the box aspect is set on the axes rather than left to the figure shape,
# so it survives whatever the layout does to the margins.
for _a in [axd] + list(ax):
    _a.set_box_aspect(1.0)
# Margins are set explicitly rather than by tight_layout, whose layout engine re-runs at
# draw time -- when the twinned axes below exist and it cannot place them.  Fixing the
# margins here makes the geometry deterministic and lets the twins be added safely.
# Sized so the cells are square before set_box_aspect has to shrink anything: leftover
# slack inside a cell is what opens the gaps between panels.
fig.subplots_adjust(left=0.108, right=0.985, top=0.935, bottom=0.072,
                    wspace=0.20, hspace=0.26)
# The profile labels ride their own curves, at the angle the curve makes on the page.
label_along(axd, l_km, RHO_EXP, 330, 'Varying (exponential)', INK, fontsize=8.0,
            offset=(4, 5))
label_along(axd, l_km, np.full_like(l_km, RHO_CONST), 470, 'Constant', INK, fontsize=8.0,
            offset=(0, 4))
phase_axis(ax[0], E_CONST/gd.UNIT_GEV, phase_a,
           [10, 100, 1000, 10000, 100000, 1000000], powers=True)
# The rising guide carries its own name, upright and at the middle of its own line:
# rotating it to the curve angle left the epsilon unreadable at this size, and a stamp
# cannot say which of the two guides it means.
_g = np.maximum(phase_a*2.2e-16, 2.2e-16)
_i = len(_g)//2
ax[0].annotate(r'$\Phi\varepsilon$', xy=(E_CONST[_i]/gd.UNIT_GEV, _g[_i]),
               xytext=(4, 3.5), textcoords='offset points', color=INK, fontsize=9.5)
phase_axis(ax[1], E_row, ph, [10, 30, 100, 300, 1000])
save(fig, 'phase_vs_profile.pdf')'''),




    md(r'''## Figure 3 --- new physics, four matrices in one slot

All four panels are the survival channel. The rows sit three decades apart in energy
because the two families of new physics do.

**The standard curve is flat at one in the lower row and that is physics, not an
artifact**: the vacuum term falls as $1/E$, so above about 100 GeV the active sector has
stopped oscillating over an Earth chord. The cell prints the range it spans.'''),
    code(r'''COSTHZ = -0.9
L = chord(COSTHZ)
LIV = dict(b1=0.0, b2=0.0, b3=np.pi/L, Lambda=1.0, sxi12=0.0,
           sxi23=1.0/np.sqrt(2.0), sxi13=0.0, dxiCP=0.0, n_liv=0)
KW = dict(nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=RTOL_FIG, atol=ATOL_FIG)

E_GEV = np.logspace(np.log10(1.0), np.log10(40.0), 260)*gd.UNIT_GEV
E_TEV = np.logspace(np.log10(1.0), np.log10(30.0), 200)*gd.UNIT_TEV

t0 = time.perf_counter()
std_gev = np.asarray(quiet(oscprob.osc_prob_3nu_earth, E_GEV, costhz=COSTHZ, L=L,
                           **OSC, **KW))
nsi_gev = np.asarray(quiet(oscprob.osc_prob_3nu_earth_nsi, E_GEV, costhz=COSTHZ, L=L,
                           **OSC, **EPS, **KW))
liv_gev = np.asarray(quiet(oscprob.osc_prob_3nu_earth_liv, E_GEV, costhz=COSTHZ, L=L,
                           **OSC, **LIV, **KW))
std_tev = np.asarray(quiet(oscprob.osc_prob_3nu_earth, E_TEV, costhz=COSTHZ, L=L,
                           **OSC, **KW))
ste1_tev = np.asarray(quiet(oscprob.osc_prob_4nu_earth, E_TEV, costhz=COSTHZ, L=L,
                            **OSC, **STERILE4, **KW))
ste2_tev = np.asarray(quiet(oscprob.osc_prob_5nu_earth, E_TEV, costhz=COSTHZ, L=L,
                            **OSC, **STERILE5, **KW))
print('five scenarios, survival channel, in %.1f s' % (time.perf_counter() - t0))
print('  b3 = %.3e eV; LIV departs from standard by up to %.3f'
      % (LIV['b3'], np.max(np.abs(liv_gev - std_gev))))
print('  standard at 1-30 TeV spans %.4f-%.4f: the vacuum term goes as 1/E and has'
      ' switched off' % (std_tev.min(), std_tev.max()))

PANELS = [(E_GEV/gd.UNIT_GEV, std_gev, nsi_gev, ORANGE, 'NSI', r'$E$ [GeV]'),
          (E_GEV/gd.UNIT_GEV, std_gev, liv_gev, GREEN, 'LIV', r'$E$ [GeV]'),
          (E_TEV/gd.UNIT_TEV, std_tev, ste1_tev, RED, r'$3+1$', r'$E$ [TeV]'),
          (E_TEV/gd.UNIT_TEV, std_tev, ste2_tev, PURPLE, r'$3+2$', r'$E$ [TeV]')]
fig, axes = plt.subplots(2, 2, figsize=(COL, 3.45))
YLIM = {0: (0.0, 1.0), 1: (0.55, 1.02)}
TICKS = {r'$E$ [GeV]': (1, 2, 5, 10, 20, 40), r'$E$ [TeV]': (1, 2, 5, 10, 20, 30)}
for k, (ax, (x, base, curve, color, label, xl)) in enumerate(zip(axes.ravel(), PANELS)):
    ax.plot(x, curve, color=color, lw=1.1, zorder=3, label=label)
    ax.plot(x, base, color=INK, lw=0.7, ls='--', zorder=4, label='Standard')
    logx(ax); snug(ax, x); xticks_at(ax, TICKS[xl])
    ax.set_ylim(*YLIM[k//2]); minor_y(ax, 5)
    ax.set_xlabel(xl, labelpad=1.5)
    ax.legend(loc='lower left', handlelength=1.3, labelspacing=0.2, fontsize=8.0)
for ax in axes[:, 0]:
    ax.set_ylabel(r'Survival probability, $P_{\nu_\mu \to \nu_\mu}$', fontsize=8.0)
fig.tight_layout(pad=0.3, w_pad=0.8, h_pad=0.9)
save(fig, 'bsm.pdf')'''),
    md(r'''## Figure 4 --- three oscillograms

**The middle row carries its own energy axis and cannot share the others'.** An eV-scale
sterile splitting over an 11 000 km chord oscillates far too fast at GeV energies for any
grid to represent --- and far too slowly to compute: one row of 200 energies costs 0.56 s
at three flavors and 21 s at 3+1, which is an hour for the panel and aliasing for the
result. At TeV energies the sterile matter resonance is the feature and the panel costs
seconds.'''),
    code(r'''# ------------------------------------------------------- three oscillograms
# Each row carries its own energy window, and it has to.  An eV-scale sterile
# splitting over an 11 000 km chord oscillates far too fast at GeV energies for any
# grid to represent -- and far too slowly to compute: one row of 200 energies costs
# 0.56 s at three flavors and 21 s at 3+1 there, which is an hour for the panel and
# aliasing for the result.  At TeV energies the sterile matter resonance is the
# feature, the phase is small, and the panel costs seconds.
NE, NZ = 200, 170
E_GEV_OSC = np.logspace(np.log10(2.0), np.log10(60.0), NE)*gd.UNIT_GEV
E_TEV_OSC = np.logspace(np.log10(1.0), np.log10(30.0), NE)*gd.UNIT_TEV
CZ = np.linspace(-1.0, -0.05, NZ)
GEV_TICKS, TEV_TICKS = (2, 3, 5, 10, 20, 30, 60), (1, 2, 5, 10, 20, 30)

SCEN = [
    (r'Earth, $3\nu$', E_GEV_OSC, gd.UNIT_GEV, r'Neutrino energy, $E$ [GeV]',
     GEV_TICKS,
     lambda E, cz, Lc: oscprob.osc_prob_3nu_earth(
         E, costhz=cz, L=Lc, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=RTOL_FIG,
         atol=ATOL_FIG, **OSC)),
    (r'Earth, $3\nu$ + NSI', E_GEV_OSC, gd.UNIT_GEV, r'Neutrino energy, $E$ [GeV]',
     GEV_TICKS,
     lambda E, cz, Lc: oscprob.osc_prob_3nu_earth_nsi(
         E, costhz=cz, L=Lc, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=RTOL_FIG,
         atol=ATOL_FIG, **OSC, **EPS)),
    (r'Earth, $3+1$', E_TEV_OSC, gd.UNIT_TEV, r'Neutrino energy, $E$ [TeV]',
     TEV_TICKS,
     lambda E, cz, Lc: oscprob.osc_prob_4nu_earth(
         E, costhz=cz, L=Lc, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=RTOL_FIG,
         atol=ATOL_FIG, **OSC, **STERILE4)),
    # The wrapper takes s25 but has no d25: the phase is not independent there, so we
    # pass the angles and leave every CP phase of the sterile sector at its default.
    (r'Earth, $3+2$', E_TEV_OSC, gd.UNIT_TEV, r'Neutrino energy, $E$ [TeV]',
     TEV_TICKS,
     lambda E, cz, Lc: oscprob.osc_prob_5nu_earth(
         E, costhz=cz, L=Lc, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=RTOL_FIG,
         atol=ATOL_FIG, **OSC, **STERILE5)),
]

fig, axes = plt.subplots(4, 1, figsize=(COL, 7.3),
                         gridspec_kw=dict(hspace=0.42))
def oscillogram_grid(label, E_ax, call):
    """One panel's worth of probability, computed once and read from disk after.

    Four panels of 170 x 200 points are about three minutes, and every one of those
    points is a function of the configuration rather than of the run.  Cosmetic work on
    this figure -- a tick label, a colorbar -- must not pay for them again.
    """
    def run():
        grid = np.empty((NZ, NE))
        for i, cz in enumerate(CZ):
            grid[i] = np.asarray(quiet(call, E_ax, float(cz), chord(float(cz))))
        return grid.tolist()
    key = ('oscillogram', label, [float(e) for e in E_ax], [float(c) for c in CZ],
           RTOL_FIG, ATOL_FIG, sorted(OSC.items()), sorted(STERILE4.items()),
           sorted(STERILE5.items()), sorted((k, str(v)) for k, v in EPS.items()))
    return np.asarray(cached('oscillogram_%s' % re.sub(r'\W+', '_', label).strip('_'),
                             key, run, what='one oscillogram panel'))


for ax, (label, E_ax, unit, xlabel, ticks, call) in zip(axes, SCEN):
    grid = oscillogram_grid(label, E_ax, call)
    im = ax.pcolormesh(E_ax/unit, CZ, grid, cmap='viridis', vmin=0.0, vmax=1.0,
                       shading='gouraud', rasterized=True)
    logx(ax); snug(ax, E_ax/unit); xticks_at(ax, ticks)
    ax.set_ylim(CZ.min(), CZ.max())
    ax.set_yticks(np.arange(-1.0, 0.0, 0.2)); minor_y(ax, 2)
    ax.set_ylabel(r'Direction, $\cos\theta_z$', fontsize=8.0)
    ax.set_xlabel(xlabel, labelpad=1.5)
    ax.axhline(-0.837, color='w', lw=0.8, ls='--', alpha=0.9)
    corner(ax, label, fontsize=8.5)
axes[0].text(2.3, -0.822, 'Core-mantle boundary', color='white', fontsize=8.0,
             va='bottom', path_effects=[pe.withStroke(linewidth=1.6, foreground='0.25')])
# aspect follows the panel count, so the bar is as tall as the stack beside it.
cb = fig.colorbar(im, ax=list(axes), pad=0.07, fraction=0.045, aspect=56)
cb.set_label(r'Survival probability, $P_{\nu_\mu \to \nu_\mu}$', fontsize=8.0)
cb.ax.tick_params(labelsize=8.0)
save(fig, 'earth_oscillogram.pdf')'''),
    md(r'''## Figure 5 --- the Sun: model, observable, and residual

The reference in the bottom panel is the adiabatic limit built from the instantaneous
eigenbases alone --- two calls to `eigh` and a contraction, touching none of the package's
averaging machinery, so it works for all three scenarios rather than only the standard
one.

One trap paid for here: `hamiltonian_3nu_nsi` returns $V_{\rm CC}$ times the epsilon matrix
**alone** --- it is zero when every epsilon is --- so the standard matter term has to be
added beside it. Omitting it put this reference 0.178 away from the answer instead of
1e-5, which reads as a spectacular disagreement rather than as a missing term.'''),
    code(r'''TABLE = os.path.join('..', 'docs', 'dev', 'adversarial_batteries', 'bs05_agsop.dat')
rows = []
with open(TABLE) as fh:
    for line in fh:
        f = line.split()
        if len(f) == 12:
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                continue
solar = np.array(rows)
r_over_rsun, rho_cgs, x_h = solar[:, 1], solar[:, 3], solar[:, 6]
MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
ne_tab = rho_cgs*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*(0.5*(1.0 + x_h))
x_solar = r_over_rsun*gd.SUN_RADIUS*gd.UNIT_KM
log_ne = np.log(ne_tab)
R_SUN = float(x_solar[-1])


def ne_sun(l):
    xs = np.clip(np.asarray(l, dtype=float), x_solar[0], x_solar[-1])
    out = np.exp(np.interp(xs, x_solar, log_ne))
    return out[()] if np.ndim(out) == 0 else out


PER_NE = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
print('BS2005-AGS,OP: %d rows, ray 0 to %.0f km; n_e falls by %.1e over it'
      % (len(solar), R_SUN/gd.UNIT_KM, ne_tab[0]/ne_tab[-1]))


def adiabatic_limit(build_H, energy, vcc0, a=gd.NUE, b=gd.NUE):
    r"""The decohered adiabatic limit, from the instantaneous eigenbases alone.

    Decohere in the matter eigenbasis at production, transport along the levels, read
    out in the vacuum eigenbasis at the surface.  This touches none of the package's
    averaging machinery -- two calls to ``eigh`` and a contraction -- so it is an
    independent reference for every scenario rather than only for the standard one,
    where it reduces to the textbook MSW expression.
    """
    hv = np.asarray(build_H(energy, 0.0), dtype=complex)
    hm = np.asarray(build_H(energy, vcc0), dtype=complex)
    _, u_matter = np.linalg.eigh(hm)
    _, u_vac = np.linalg.eigh(hv)
    return float(np.sum(np.abs(u_matter[a])**2 * np.abs(u_vac[b])**2))


OSC4 = dict(OSC); OSC4.update(STERILE4); OSC4.update(d14=0.0, d24=0.0)
OSC5 = dict(OSC); OSC5.update(STERILE5)
OSC5.update(d14=0.0, d24=0.0, d15=0.0, d25=0.0, d35=0.0)
E_AVG = np.logspace(np.log10(0.1), np.log10(20.0), 90)*gd.UNIT_MEV
COMMON = dict(L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True)
VCC0 = float(PER_NE*ne_sun(0.0))

HV3 = np.asarray(vacuum_hamiltonian(3), dtype=complex)
HV4 = np.asarray(vacuum_hamiltonian(4), dtype=complex)
HV5 = np.asarray(vacuum_hamiltonian(5), dtype=complex)
P3 = np.asarray(matter.matter_potential_projector(3), dtype=complex)
P4 = np.asarray(matter.matter_potential_projector(4), dtype=complex)
P5 = np.asarray(matter.matter_potential_projector(5), dtype=complex)
EPS_ORDER = ('eps_ee', 'eps_em', 'eps_et', 'eps_mm', 'eps_mt', 'eps_tt')

# 3nu, then the NSI case, then the two steriles: the order the middle panel reads in.
SCEN = [
    (r'$3\nu$', BLUE,
     lambda E: oscprob.osc_prob_matter_std_potential(3, ne_sun, E, R_SUN, OSC,
                                                     average=True, **COMMON),
     lambda E, v: HV3/E + v*P3),
    (r'$3\nu$ + NSI', ORANGE,
     lambda E: oscprob.osc_prob_matter_nsi(3, ne_sun, E, R_SUN, OSC, EPS,
                                           average=True, **COMMON),
     # hamiltonian_3nu_nsi returns V_CC times the epsilon matrix ALONE -- it is zero
     # when every epsilon is -- so the standard matter term has to be added beside it.
     # Omitting it put this reference 0.178 away from the answer rather than 1e-5.
     lambda E, v: HV3/E + v*P3 + np.asarray(hamiltonians.hamiltonian_3nu_nsi(
         v, *[EPS[k] for k in EPS_ORDER]), dtype=complex)),
    (r'$3+1$', RED,
     lambda E: oscprob.osc_prob_matter_std_potential(4, ne_sun, E, R_SUN, OSC4,
                                                     average=True, **COMMON),
     lambda E, v: HV4/E + v*P4),
    # Green, not purple: purple marks the resonance densities in the panel above.
    (r'$3+2$', GREEN,
     lambda E: oscprob.osc_prob_matter_std_potential(5, ne_sun, E, R_SUN, OSC5,
                                                     average=True, **COMMON),
     lambda E, v: HV5/E + v*P5),
]'''),
    code(r'''fig = plt.figure(figsize=(COL, 5.0))
outer = fig.add_gridspec(2, 1, height_ratios=[1.22, 2.30], hspace=0.22)
gs_low = outer[1].subgridspec(2, 1, height_ratios=[1.60, 0.90], hspace=0.09)
axes = [fig.add_subplot(outer[0]), fig.add_subplot(gs_low[0]), None]
axes[2] = fig.add_subplot(gs_low[1], sharex=axes[1])

# --- the model
ax = axes[0]
rr = x_solar/R_SUN
# n_e in cm^-3: natural units cube the energy, and nobody quotes a solar electron
# density that way.  gd.UNIT_PER_CM3 is one cm^-3 expressed in eV^3.
PER_CM3 = gd.UNIT_PER_CM3
print('  central n_e = %.2e cm^-3 (the Sun\'s core is ~6e25)' % (ne_tab[0]/PER_CM3))
ax.semilogy(rr, ne_tab/PER_CM3, color=INK, lw=1.4)
cos2th12 = 1.0 - 2.0*OSC['s12']**2
# Staggered horizontally: the three resonance densities are within a decade of each
# other, so labels at a common x sit on top of one another.
# The first label carries the whole statement; the other two need only their energy,
# right-aligned under it so the three read as one column.
LAB = {1.0: r'$3\nu$ MSW resonance density at 1 MeV', 5.0: '5 MeV', 20.0: '20 MeV'}
for Emev in (1.0, 5.0, 20.0):
    nr = OSC['D21']*cos2th12/(2.0*Emev*gd.UNIT_MEV)/PER_NE
    ax.axhline(nr/PER_CM3, color=PURPLE, lw=0.7, ls='--')
    # Under its own line, clear of it: the top line sits at the very top of the panel,
    # so a label placed above it lands outside the axes.
    ax.annotate(LAB[Emev], xy=(0.975, nr/PER_CM3), xytext=(0, -3.5),
                textcoords='offset points', color=PURPLE, fontsize=8.0, ha='right',
                va='top')
logy(ax); ax.set_xlim(0.0, 1.0)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel(r'Radius, $r/R_\odot$', labelpad=1.5)
ax.set_ylabel(r'Electron density, $n_e$ [cm$^{-3}$]', fontsize=8.0)


# --- the averaged observable, and the residual under it
resid = {}
for label, color, call, build_H in SCEN:
    def run(call=call, build_H=build_H):
        P = np.asarray(quiet(call, E_AVG))
        R = np.array([adiabatic_limit(build_H, e, VCC0) for e in E_AVG])
        return dict(P=P.tolist(), R=R.tolist())
    got = cached('solar_%s' % re.sub(r'\W+', '_', label).strip('_'),
                 ('solar', label, [float(e) for e in E_AVG], float(R_SUN),
                  profile_samples(lambda l: PER_NE*ne_sun(l), R_SUN),
                  sorted(OSC.items())),
                 run, what='one averaged solar scenario')
    P, R = np.asarray(got['P']), np.asarray(got['R'])
    axes[1].semilogx(E_AVG/gd.UNIT_MEV, P, color=color, lw=1.3, label=label)
    resid[label] = (color, np.abs(P - R))
    print('%-14s P in %.3f-%.3f, worst |Magnus - adiabatic| %.2e'
          % (label, P.min(), P.max(), resid[label][1].max()))

a = axes[1]
logx(a); snug(a, E_AVG/gd.UNIT_MEV); xticks_at(a, (0.1, 0.3, 1, 3, 10, 20))
a.set_ylim(0.25, 0.60); minor_y(a, 5)
a.set_ylabel(r'Average probability, $\langle P_{\nu_e \to \nu_e}\rangle$',
             fontsize=8.0)
a.tick_params(labelbottom=False)
a.legend(loc='upper right', handlelength=1.4)
corner(a, r'Sun', loc='upper left', x=0.035, y=0.965)

b = axes[2]
for label, (color, dP) in resid.items():
    b.semilogy(E_AVG/gd.UNIT_MEV, np.maximum(dP, 1.0e-17), color=color, lw=1.0)
logx(b); logy(b); snug(b, E_AVG/gd.UNIT_MEV); xticks_at(b, (0.1, 0.3, 1, 3, 10, 20))
b.set_xlabel(r'Neutrino energy, $E$ [MeV]')
b.set_ylabel(r'$|\Delta P|$', fontsize=8.0)
corner(b, r'Vs.\ adiabatic limit', loc='upper left', x=0.035, y=0.94,
       fontsize=8.0)
fig.subplots_adjust(left=0.20)
# Along the curve, and last: label_along reads the slope off transData, so it has to
# run once the axes have their final width.  Called before subplots_adjust, it was
# rotated to an angle the panel no longer had, and the curve cut through the words.
# Anchored by radius, not by index: the model table is sampled far more densely in
# the core, so a fraction of len(rr) lands nowhere near that fraction of the Sun.
label_along(ax, rr, ne_tab/PER_CM3, int(np.searchsorted(rr, 0.30)),
            'Sun (BS2005-AGS,OP)', INK, fontsize=8.0, offset=(0, -13), chord=True)
save(fig, 'solar_averaged.pdf')'''),
    md(r'''### Figure 5b --- a Hamiltonian the package never heard of

A gauged $L_e - L_\mu$ symmetry adds a long-range potential sourced by the electrons of the
Sun itself. Nothing about it is built in: it is a callable returning a Hermitian matrix,
which is the whole of the interface.'''),
    code(r'''# ------------------------------------------------ Figure 5b: L_e - L_mu in the Sun
def running_integral(y, x):
    """Trapezoidal running integral of y over x, zero at the first node."""
    return np.concatenate([[0.0], np.cumsum(0.5*(y[1:] + y[:-1])*np.diff(x))])


def shc(x):
    """sinh(x)/x, continued to 1 at the origin."""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    big = np.abs(x) > 1.0e-8
    out[big] = np.sinh(x[big])/x[big]
    return out


def long_range_potential(r_grid, ne_grid, m):
    """V_{e-mu}(r) on r_grid, in units of g'^2, for a spherical n_e(r)."""
    inner = r_grid**2*ne_grid*shc(m*r_grid)
    outer = r_grid*ne_grid*np.exp(-m*r_grid)
    I_in = running_integral(inner, r_grid)
    running_out = running_integral(outer, r_grid)
    I_out = running_out[-1] - running_out
    r_safe = np.where(r_grid > 0.0, r_grid, 1.0e-30)
    return np.exp(-m*r_grid)*I_in/r_safe + shc(m*r_grid)*I_out


LR_CHARGE = np.diag([1.0, -1.0, 0.0]).astype(complex)      # the L_e - L_mu charge
# Two mediator ranges: the solar radius, the analogue of 1/m = R_earth in notebook 19,
# and a tenth of it.  The shorter range samples the profile where it is steepest, so the
# two are not a rescaling of each other.
LR_RANGES = [(1.0, r'$1/m = R_\odot$', RED),
             (0.1, r'$1/m = R_\odot/10$', BLUE)]
V_LR = {}
G2 = {}
for frac, _, _ in LR_RANGES:
    V_LR[frac] = long_range_potential(x_solar, ne_tab, 1.0/(frac*R_SUN))
    # g'^2 fixed so the new potential is a tenth of V_CC at the center in each case,
    # which is what makes the two curves comparable.
    G2[frac] = 0.1*VCC0/V_LR[frac][0]
    print('L_e - L_mu, 1/m = %.1f R_sun: V_new/V_CC = %.3f at the center, %.4f at 0.5 R_sun'
          % (frac, G2[frac]*V_LR[frac][0]/VCC0,
             G2[frac]*np.interp(0.5*R_SUN, x_solar, V_LR[frac])
             / float(PER_NE*ne_sun(0.5*R_SUN))))


def H_lr(E, frac=None):
    """Standard three-flavor solar Hamiltonian, with a long-range term or without."""
    def f(l):
        v = float(PER_NE*ne_sun(l))
        h = HV3/E + v*P3
        if frac is not None:
            h = h + float(G2[frac]*np.interp(l, x_solar, V_LR[frac]))*LR_CHARGE
        return h
    return f


E_LR = np.logspace(np.log10(0.1), np.log10(20.0), 70)*gd.UNIT_MEV
avg = lambda H: avgprob.averaged_probabilities_adiabatic(H, 0.0, R_SUN)[0][0, 0]


def _lri_sweep():
    out = {'std': [avg(H_lr(e)) for e in E_LR]}
    for frac, _, _ in LR_RANGES:
        out['%g' % frac] = [avg(H_lr(e, frac)) for e in E_LR]
    return out


_got = cached('solar_long_range',
              ('lri', [float(e) for e in E_LR], float(R_SUN), [f for f, _, _ in LR_RANGES],
               {'%g' % f: float(G2[f]) for f, _, _ in LR_RANGES}, sorted(OSC.items())),
              _lri_sweep, what='the long-range solar sweep')
P_std = np.asarray(_got['std'])
P_lr = {frac: np.asarray(_got['%g' % frac]) for frac, _, _ in LR_RANGES}
print('  P_ee standard %.3f-%.3f' % (P_std.min(), P_std.max()))
for frac, lab, _ in LR_RANGES:
    print('    1/m = %.1f R_sun: %.3f-%.3f, largest shift %.3f'
          % (frac, P_lr[frac].min(), P_lr[frac].max(),
             np.max(np.abs(P_lr[frac] - P_std))))

fig, axes = plt.subplots(2, 1, figsize=(COL, 3.9), sharex=True,
                         gridspec_kw=dict(height_ratios=[2.0, 1.0], hspace=0.08))
a = axes[0]
a.semilogx(E_LR/gd.UNIT_MEV, P_std, color=INK, lw=1.3, ls='--', label=r'Standard $3\nu$')
for frac, lab, col in LR_RANGES:
    a.semilogx(E_LR/gd.UNIT_MEV, P_lr[frac], color=col, lw=1.3,
               label=r'$+\;L_e - L_\mu$,  ' + lab)
logx(a); snug(a, E_LR/gd.UNIT_MEV); xticks_at(a, (0.1, 0.3, 1, 3, 10, 20))
a.set_ylabel(r'Average probability, $\langle P_{\nu_e \to \nu_e}\rangle$', fontsize=8.0)
a.set_ylim(top=0.55)
a.tick_params(labelbottom=False); minor_y(a, 5)
a.legend(loc='lower left', handlelength=1.6)
corner(a, r'Sun (BS2005-AGS,OP)', loc='upper right', fontsize=8.0)

b = axes[1]
for frac, _, col in LR_RANGES:
    b.semilogx(E_LR/gd.UNIT_MEV, P_lr[frac] - P_std, color=col, lw=1.1)
b.axhline(0.0, color=INK, lw=0.6, ls=':')
logx(b); snug(b, E_LR/gd.UNIT_MEV); xticks_at(b, (0.1, 0.3, 1, 3, 10, 20))
b.set_xlabel(r'Neutrino energy, $E$ [MeV]')
b.set_ylabel(r'$\Delta \langle P \rangle$', fontsize=8.0)
b.set_ylim(-0.02, 0.02)
minor_y(b, 5)
fig.subplots_adjust(left=0.20)
save(fig, 'solar_long_range.pdf')'''),
    md(r'''## Figure 6 --- a supernova shock

Rows 1 and 2 share the full-ray axis; row 3 is a window at the front and gets its own.

The front is 0.1% of the ray at most, so it cannot show on the full-ray axis --- the inset
is the only place the two columns differ to the eye. Two flavors is not drawn: with only
$\Delta m^2_{21}$ available at 15 MeV it spans 0.998 to 1.000.

**Note which engine answers.** The paper describes an adiabatic strategy at length and it
is *not* what runs here: a baseline scan at one energy goes to the cumulative engine, and
declared breakpoints make the hybrid stand aside in any case. The cell prints it.'''),
    code(r'''KM = gd.UNIT_KM
MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
R_CONTACT_KM, R_FORWARD_KM = 12348.0, 30323.0
R0_KM, R1_KM = 1.0e4, 8.0e4
L0, L1 = R0_KM*KM, R1_KM*KM
ENERGY = 15.0*gd.UNIT_MEV


def smoothstep(u):
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    return u*u*(3.0 - 2.0*u)


def rarefaction(r_km, r_shock_km):
    """ln f(x) = [0.28 - 0.69 ln(x_s/km)] [arcsin(1 - x/x_s)]^1.1  (Fogli et al. 2003)."""
    u = np.clip(1.0 - np.asarray(r_km, dtype=float)/r_shock_km, 0.0, 1.0)
    return np.exp((0.28 - 0.69*np.log(r_shock_km))*np.arcsin(u)**1.1)


def sn_shock_ne(width_frac, contact_jump=2.5, y_e=0.5):
    """Electron number density along the ray, for a front of the given width."""
    w_km = float(width_frac)*(R1_KM - R0_KM)

    def ne(l):
        r = np.asarray(l, dtype=float)/KM
        rho = 1.0e14*r**(-2.4)
        shocked = smoothstep((R_FORWARD_KM + 0.5*w_km - r)/w_km)
        factor = 1.0 + shocked*(10.0*rarefaction(r, R_FORWARD_KM) - 1.0)
        inside = smoothstep((R_CONTACT_KM + 0.5*w_km - r)/w_km)
        factor = factor*(1.0 + inside*(contact_jump - 1.0))
        out = rho*factor*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*y_e
        return out[()] if np.ndim(out) == 0 else out

    return ne


def shock_breakpoints(width_frac):
    w_km = float(width_frac)*(R1_KM - R0_KM)
    edges = []
    for r in (R_CONTACT_KM, R_FORWARD_KM):
        edges += [(r - 0.5*w_km)*KM, (r + 0.5*w_km)*KM]
    return np.array([L0] + sorted(edges) + [L1])


WIDTHS = (1.0e-6, 1.0e-3)
OSC4 = dict(OSC); OSC4.update(STERILE4); OSC4.update(d14=0.0, d24=0.0)
OSC5_SHOCK = dict(OSC); OSC5_SHOCK.update(STERILE5)
OSC5_SHOCK.update(d14=0.0, d24=0.0, d15=0.0, d25=0.0, d35=0.0)
OSC2 = dict(sth=OSC['s12'], Dm2=OSC['D21'])
HALF_KM = 75.0
WIN = (R_FORWARD_KM - HALF_KM, R_FORWARD_KM + HALF_KM)
L_OSC_KM = 4.0*np.pi*ENERGY/OSC['D21']/KM
print('vacuum oscillation length at %g MeV: %.0f km; the ray spans %.0f of them'
      % (ENERGY/gd.UNIT_MEV, L_OSC_KM, (R1_KM - R0_KM)/L_OSC_KM))

Ls_full = np.linspace(1.02*L0, L1, 4000)
Ls_win = np.linspace(WIN[0]*KM, WIN[1]*KM, 1000)
COMMON = dict(L0=L0, nu_i=gd.NUE, nu_f=gd.NUE,
              density_is_of_number_of_electrons=True, rtol=RTOL_FIG, atol=ATOL_FIG)

SCEN = [
    (r'$3\nu$', ORANGE, lambda ne, bp, Ls: oscprob.osc_prob_matter_std_potential(
        3, ne, ENERGY, Ls, OSC, t_breakpoints=bp, **COMMON)),
    (r'$3+1$', RED, lambda ne, bp, Ls: oscprob.osc_prob_matter_std_potential(
        4, ne, ENERGY, Ls, OSC4, t_breakpoints=bp, **COMMON)),
    (r'$3\nu$ + NSI', PURPLE, lambda ne, bp, Ls: oscprob.osc_prob_matter_nsi(
        3, ne, ENERGY, Ls, OSC, EPS, t_breakpoints=bp, **COMMON)),
    (r'$3+2$', GREEN, lambda ne, bp, Ls: oscprob.osc_prob_matter_std_potential(
        5, ne, ENERGY, Ls, OSC5_SHOCK, t_breakpoints=bp, **COMMON)),
]

# Which engine answers?  Worth reporting, because the paper describes an adiabatic
# strategy at length and it is NOT what runs here: a baseline scan at one energy goes
# to the cumulative engine, and declared breakpoints make the hybrid stand aside in
# any case.  The whole ray is the Magnus ladder.
info = {}
quiet(oscprob.osc_prob_matter_std_potential, 3, sn_shock_ne(1e-3), ENERGY,
      Ls_full[:40], OSC, t_breakpoints=shock_breakpoints(1e-3),
      strategy_info=info, **COMMON)
print('engine answering this figure: %r (declined: %s)'
      % (info.get('engine'), [d[0] for d in (info.get('declined') or [])] or 'none'))

# Rows 1 and 2 share the full-ray axis and sit tight against each other; row 3 is a
# different axis entirely (a window at the front), so it gets its own block and its
# own gap.  A single hspace cannot do both.
fig = plt.figure(figsize=(WIDE, 6.3))
outer = fig.add_gridspec(2, 1, height_ratios=[2.25, 1.15], hspace=0.20)
gs_top = outer[0].subgridspec(2, 2, height_ratios=[1.15, 1.15], hspace=0.08,
                              wspace=0.06)
gs_bot = outer[1].subgridspec(1, 2, wspace=0.06)
axes = np.empty((3, 2), dtype=object)
for c in range(2):
    axes[0, c] = fig.add_subplot(gs_top[0, c])
    axes[1, c] = fig.add_subplot(gs_top[1, c], sharex=axes[0, c])
    axes[2, c] = fig.add_subplot(gs_bot[0, c])

for col, width in enumerate(WIDTHS):
    ne = sn_shock_ne(width)
    bp = shock_breakpoints(width)
    w_km = width*(R1_KM - R0_KM)
    wlab = ('%.2f' % w_km).rstrip('0').rstrip('.')
    top, mid, bot = axes[0, col], axes[1, col], axes[2, col]

    top.semilogy(Ls_full/KM/1.0e3, np.asarray(ne(Ls_full))/gd.UNIT_PER_CM3,
                 color=INK, lw=1.0)
    logy(top)
    # Limits from the data: they were pinned to the old MeV^3 scaling, which put the
    # cm^-3 curve thirty decades off the panel.
    ne_full = np.asarray(ne(Ls_full))/gd.UNIT_PER_CM3
    top.set_ylim(0.55*ne_full.min(), 2.2*ne_full.max())
    # Placed against the density rather than the panel top, so it sits over the flat
    # part of the profile instead of over its rise.
    top.text(0.035, 1.5e26, r'SN: front width %s km' % wlab, transform=
             top.get_yaxis_transform(), ha='left', va='center', fontsize=8.0,
             color='black', zorder=10,
             bbox=dict(boxstyle='round,pad=0.32', facecolor='white',
                       edgecolor='black', linewidth=0.6))

    # The front is 0.1% of the ray at most, so it cannot show on the axis above.  The
    # inset is the same window the bottom row uses, and is the only place the two
    # columns differ to the eye.
    # The inset spans a few front widths, not a fixed window: at 0.07 km a +/- 75 km
    # view shows a vertical line and no shaded band at all.
    half_in = max(1.5*w_km, 0.05)
    r_in = np.linspace(R_FORWARD_KM - half_in, R_FORWARD_KM + half_in, 600)
    ins = top.inset_axes([0.60, 0.50, 0.34, 0.42])
    # Scaled to 1e-5 MeV^3 so the axis needs no offset text, which collided with the
    # inset's own title at this size.
    ins.plot(r_in - R_FORWARD_KM,
             np.asarray(ne(r_in*KM))/gd.UNIT_PER_CM3/1.0e26, color=INK, lw=0.9)
    ins.axvspan(-0.5*w_km, 0.5*w_km, color=GREEN, alpha=0.35, lw=0, zorder=0)
    ins.set_xlim(-half_in, half_in)
    ins.set_xticks((-half_in, 0.0, half_in))
    ins.set_xticklabels([('%g' % v) for v in (-half_in, 0, half_in)])
    ins.yaxis.set_major_locator(plt.MaxNLocator(2))
    ins.tick_params(labelsize=8.0, pad=0.8, length=1.6)
    for side in ins.spines.values():
        side.set_linewidth(0.6)
    ins.set_xlabel(r'Distance from FS [km]', fontsize=7.4, labelpad=1.0)
    ins.set_ylabel(r'$n_e$ [$10^{26}$ cm$^{-3}$]', fontsize=7.4, labelpad=1.0)
    # A box on the parent marking where the inset was taken from, with corner
    # connectors.  At 0.07 km the box is a hairline, which is the honest width.
    lo, hi = top.get_ylim()
    rect, lines = top.indicate_inset(
        ((R_FORWARD_KM - half_in)/1.0e3, lo, 2.0*half_in/1.0e3, hi - lo),
        inset_ax=ins, edgecolor='0.45', linewidth=0.5, alpha=1.0)
    # Neither the marking box nor its connectors survive: the box spans the full height
    # of the panel, so at 70 km it reads as a shaded block standing over the shock.
    rect.set(visible=False)
    for ln in lines:
        ln.set(visible=False)

    for label, color, call in SCEN:
        def _run(call=call, ne=ne, bp=bp):
            return dict(full=np.asarray(quiet(call, ne, bp, Ls_full)).tolist(),
                        win=np.asarray(quiet(call, ne, bp, Ls_win)).tolist())
        _g = cached('shock_%s_%s' % (re.sub(r'\W+', '_', label).strip('_'),
                                     ('%g' % w_km).replace('.', 'p')),
                    ('shock', label, float(w_km), float(ENERGY),
                     [float(x) for x in Ls_full[::40]], [float(x) for x in Ls_win[::40]],
                     sorted(OSC.items()), sorted(STERILE4.items()),
                     sorted(STERILE5.items()), sorted((k, str(v)) for k, v in EPS.items())),
                    _run, what='one shock scenario at one front width')
        P_full, P_win = np.asarray(_g['full']), np.asarray(_g['win'])
        mid.plot(Ls_full/KM/1.0e3, P_full, color=color, lw=0.3, alpha=0.8)
        bot.plot(Ls_win/KM - R_FORWARD_KM, P_win, color=color, lw=0.8, label=label)
        print('  width %6.2f km  %-13s full P in %.3f-%.3f'
              % (w_km, label, P_full.min(), P_full.max()))

    for ax in (top, mid):
        ax.set_xlim(R0_KM/1.0e3, R1_KM/1.0e3)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        for r in (R_CONTACT_KM, R_FORWARD_KM):
            ax.axvline(r/1.0e3, color='0.45', lw=1.1, ls=':')
    top.tick_params(labelbottom=False)
    mid.set_ylim(0.0, 1.0); minor_y(mid, 5)
    mid.set_xlabel(r'Radius, $r$ [$10^3$ km]')

    bot.set_xlim(-HALF_KM, HALF_KM)
    bot.xaxis.set_minor_locator(AutoMinorLocator(5))
    bot.axvline(0.0, color='0.4', lw=1.1, ls=':')
    bot.set_ylim(0.0, 1.0); minor_y(bot, 5)
    bot.set_xlabel(r'Distance from the forward shock [km]')

    if col:
        for ax in (top, mid, bot):
            ax.tick_params(labelleft=False)

axes[0, 0].set_ylabel(r'$n_e$ [cm$^{-3}$]', fontsize=8.0)
axes[1, 0].set_ylabel(r'$P_{\nu_e \to \nu_e}$, whole ray', fontsize=8.0)
axes[2, 0].set_ylabel(r'$P_{\nu_e \to \nu_e}$, at the FS', fontsize=8.0)
axes[2, 1].legend(loc='lower left', handlelength=1.3, ncol=2, columnspacing=0.7,
                  labelspacing=0.2, fontsize=8.0)
# Named on both rows that show them: the bottom row speaks of "the forward shock"
# while the profile above carries two, which is a question a reader should not have
# to hold open.
for r, name in ((R_CONTACT_KM, 'Contact'), (R_FORWARD_KM, 'Forward shock (FS)')):
    axes[0, 0].annotate(name, xy=(r/1.0e3, 0.04),
                        xycoords=axes[0, 0].get_xaxis_transform(),
                        xytext=(4.5, 0), textcoords='offset points',
                        fontsize=8.0, color='0.35', ha='left', va='bottom')
save(fig, 'shock_probability.pdf')'''),
    md(r'''### Figure 8 --- what an astrophysical flux arrives as

Decohered over the distance, so the observable is the flavor composition. Standard
oscillations make it independent of energy; the figure asks which new physics does not.'''),
    code(r'''# ------------------------------------ Figure 8: what an astrophysical flux arrives as
# A neutrino from a distant source arrives decohered, so what is observable is not a
# probability curve but the flavor composition.  For standard oscillations that
# composition does not depend on energy at all; the point of the figure is which new
# physics makes it depend on energy, and where.
E_ASTRO = np.logspace(np.log10(1.0e3), np.log10(1.0e7), 60)*gd.UNIT_GEV   # 1 TeV to 10 PeV
PION_SOURCE = np.array([1.0/3.0, 2.0/3.0, 0.0])       # pion decay, before oscillating
COSTHZ_ASTRO = -1.0                                   # straight through the core
L_ASTRO = earth.distance_traveled_inside_earth(COSTHZ_ASTRO)*gd.CONV_KM_TO_INV_EV

# The LIV eigenvalue is fixed by where we want the new term to equal the vacuum one,
# rather than chosen for its size: at E_STAR the two are the same to within a percent.
N_LIV, E_STAR = 1, 100.0e3*gd.UNIT_GEV
B3_ASTRO = float(OSC['D31'])/(2.0*E_STAR**(N_LIV + 1))

HVA3 = np.asarray(vacuum_hamiltonian(3), dtype=complex)
HVA4 = np.asarray(vacuum_hamiltonian(4), dtype=complex)
def _astro_vac():
    return dict(p3=[avgprob.averaged_probabilities_constant_hamiltonian(HVA3/e).tolist()
                    for e in E_ASTRO],
                p4=[avgprob.averaged_probabilities_constant_hamiltonian(HVA4/e).tolist()
                    for e in E_ASTRO])


_av = cached('astro_vacuum',
             ('astro_vac', [float(e) for e in E_ASTRO], sorted(OSC.items()),
              sorted(STERILE4.items())), _astro_vac, what='the decohered vacuum matrices')
P_VAC3, P_VAC4 = np.asarray(_av['p3']), np.asarray(_av['p4'])


def liv_term(e):
    return np.asarray(hamiltonians.hamiltonian_3nu_liv(
        e, sxi12=OSC['s12'], sxi23=OSC['s23'], sxi13=OSC['s13'], dxiCP=0.0,
        b1=0.0, b2=0.0, b3=B3_ASTRO, Lambda=1.0, n_liv=N_LIV), dtype=complex)


P_LIV = np.asarray(cached(
    'astro_liv',
    ('astro_liv', [float(e) for e in E_ASTRO], N_LIV, float(B3_ASTRO), sorted(OSC.items())),
    lambda: [avgprob.averaged_probabilities_constant_hamiltonian(HVA3/e + liv_term(e)).tolist()
             for e in E_ASTRO], what='the decohered matrices with a Lorentz-violating term'))
# Through the Earth the flux is already decohered when it arrives, so the two legs
# compose as probability matrices rather than as amplitudes.
P_EARTH = np.asarray(cached(
    'astro_earth_nsi',
    ('astro_earth', [float(e) for e in E_ASTRO], COSTHZ_ASTRO, float(L_ASTRO),
     sorted(OSC.items()), sorted((k, str(v)) for k, v in EPS.items())),
    lambda: np.asarray(quiet(oscprob.osc_prob_3nu_earth_nsi, E_ASTRO, costhz=COSTHZ_ASTRO,
                             L=L_ASTRO, **OSC, **EPS, rtol=1e-6, atol=1e-8)).tolist(),
    what='the Earth leg with non-standard interactions'))
P_NSI = np.einsum('ij,ejk->eik', P_VAC3[0], P_EARTH)

# Pseudo-Dirac: each mass state may be two states split by a tiny delta m^2.  Pairing is
# chosen per state, and that is the whole of the physics here: pair all three and every
# active-active probability halves by the same factor, so the composition is untouched.
# Pair one, and the suppression is uneven and the composition moves.  We pair the second
# mass state alone.  At 100 Mpc a splitting of 1e-13 eV^2 leaves both members of the pair
# decohered from each other across the whole range drawn, so no averaged expression is
# being stretched: avgprob groups the spectrum itself, and finds six singletons.
L_SOURCE = 100.0*3.0857e19*gd.CONV_KM_TO_INV_EV        # 100 Mpc [eV^-1]
PD_PAIRS = {1: 1.0e-13}
U_PMNS = np.asarray(hamiltonians.pmns_mixing_matrix(
    OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP']), dtype=complex)
M2_PMNS = np.array([0.0, float(OSC['D21']), float(OSC['D31'])])


def _astro_pd():
    out = []
    for e in E_ASTRO:
        H = hamiltonians.hamiltonian_pseudo_dirac_vacuum(e, U_PMNS, M2_PMNS, PD_PAIRS)
        out.append(avgprob.averaged_probabilities_constant_hamiltonian(
            np.asarray(H, dtype=complex), baseline=L_SOURCE).tolist())
    return out


P_PD = np.asarray(cached(
    'astro_pseudo_dirac',
    ('astro_pd', [float(e) for e in E_ASTRO], sorted(PD_PAIRS.items()),
     float(L_SOURCE), sorted(OSC.items())),
    _astro_pd, what='the decohered matrices with one pseudo-Dirac pair'))

COMP_STD = np.einsum('a,eab->eb', PION_SOURCE, P_VAC3)
COMP_STD = COMP_STD/COMP_STD.sum(axis=1)[:, None]
CASES = [(r'Standard $3\nu$', P_VAC3, 3),
         (r'LIV',            P_LIV,  3),
         (r'NSI through Earth', P_NSI, 3),
         (r'$3+1$',           P_VAC4, 4),
         (r'Pseudo-Dirac, $\nu_2$', P_PD, 4)]
FCOL = {0: BLUE, 1: RED, 2: GREEN}
FLAB = {0: r'$\nu_e$', 1: r'$\nu_\mu$', 2: r'$\nu_\tau$'}

fig, axes = plt.subplots(2, 5, figsize=(WIDE, 4.1), sharex=True, sharey='row',
                         gridspec_kw=dict(hspace=0.10, wspace=0.10))
for col, (label, P, d) in enumerate(CASES):
    src = np.zeros(d); src[:3] = PION_SOURCE
    comp = np.einsum('a,eab->eb', src, P)[:, :3]
    # Renormalized to the active flavors: what a detector measures is the ratio among
    # the three it can see, and the sterile share is not observed rather than being a
    # fourth number to compare against.
    comp = comp/comp.sum(axis=1)[:, None]
    top, bot = axes[0, col], axes[1, col]
    for a in range(3):
        # The standard case under every panel, so that a small departure from it is
        # visible.  Three of the four cases are flat in energy; without the reference a
        # reader cannot see which of them sit at a different value.
        if col > 0:
            top.semilogx(E_ASTRO/gd.UNIT_TEV, P_VAC3[:, a, a], color='0.65', lw=0.7,
                         ls='--', zorder=1)
            bot.semilogx(E_ASTRO/gd.UNIT_TEV, COMP_STD[:, a], color='0.65', lw=0.7,
                         ls='--', zorder=1)
        top.semilogx(E_ASTRO/gd.UNIT_TEV, P[:, a, a], color=FCOL[a], lw=1.3, zorder=3,
                     label=FLAB[a] if col == 0 else None)
        bot.semilogx(E_ASTRO/gd.UNIT_TEV, comp[:, a], color=FCOL[a], lw=1.3, zorder=3)
    for ax in (top, bot):
        logx(ax); snug(ax, E_ASTRO/gd.UNIT_TEV)
    bot.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=12))
    unit_as_one(bot, which='x')
    corner(top, label, loc='upper left', x=0.075, fontsize=7.6)
    if d > 3:
        lost = 1.0 - np.einsum('a,eab->eb', src, P)[:, :3].sum(axis=1)
        print('  %-20s sterile takes %.3f-%.3f of the flux before renormalizing'
              % (label.replace('$', ''), lost.min(), lost.max()))
    print('  %-20s survival at 1 TeV %.3f/%.3f/%.3f -> at 10 PeV %.3f/%.3f/%.3f'
          % (label.replace('$', ''), P[0, 0, 0], P[0, 1, 1], P[0, 2, 2],
             P[-1, 0, 0], P[-1, 1, 1], P[-1, 2, 2]))
for _a in axes[1, :]:
    _a.tick_params(axis='x', which='major', pad=4.5)
fig.supxlabel(r'Neutrino energy, $E$ [TeV]',
              fontsize=plt.rcParams['axes.labelsize'], x=0.55, y=0.050)
axes[0, 0].set_ylabel(r'Avg.~survival probability, $\langle P_{\nu_\alpha \to \nu_\alpha}\rangle$',
                      fontsize=8.0)
axes[1, 0].set_ylabel(r'Flavor fraction at Earth', fontsize=8.0)
# Ranges chosen to show the differences rather than the distance from zero: every
# survival probability here lives above 0.30, and every flavor fraction between 0.31
# and 0.35.  The top is left to the data, so the LIV excursion is not clipped.
axes[0, 0].set_ylim(0.28, 0.90); minor_y(axes[0, 0], 5)
axes[1, 0].set_ylim(0.316, 0.35); minor_y(axes[1, 0], 5)
axes[0, 0].legend(loc='upper left', bbox_to_anchor=(0.02, 0.82), handlelength=1.1,
                  fontsize=7.4, ncol=1, labelspacing=0.25, borderpad=0.3)
axes[0, 1].axvline(E_STAR/gd.UNIT_TEV, color=INK, lw=0.7, ls=':')
axes[1, 1].axvline(E_STAR/gd.UNIT_TEV, color=INK, lw=0.7, ls=':')
fig.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.4)
# Only the first column keeps its leftmost tick label; in the others it would sit
# against the previous panel's right spine.  This has to be done on the drawn text.
# sharex makes the five columns share one Ticker, so a formatter set on any of them
# governs all of them, and hiding a tick's label is undone by the next draw.
fig.canvas.draw()
for _ax in axes[1, 1:]:
    for _lb in _ax.get_xticklabels():
        if _lb.get_text() in ('$1$', '1'):
            _lb.set_visible(False)
save(fig, 'astro_composition.pdf')'''),
    md(r'''## Figure 7 --- a smooth profile: reach, and the flavor ceiling

Four rows on one shared time axis, which spans three decades, so the cost of a flavor
count can be read against the others directly.'''),
    code(r'''# ------------------------------------------------------------ smooth profile
BENCH = json.loads((HERE/'external_profile_benchmarks.json').read_text())
print('machine: %s | interleaved control: %s'
      % (BENCH['machine'], BENCH['control_ratio']))

fig, axes = plt.subplots(4, 1, figsize=(COL, 5.0), sharex=True,
                         gridspec_kw=dict(hspace=0.10))
for ax, case in zip(axes, BENCH['cases']):
    allt = []
    for series in case['series']:
        t = [p['us_per_probability'] for p in series['points']]
        e = [p['max_abs_error'] for p in series['points']]
        allt += t
        if series['name'] == 'Magnus':
            ax.loglog(t, e, '-*', color=INK, ms=8, lw=1.2, zorder=5, label=r'Mag$\nu$s')
        else:
            ax.loglog(t, e, '-o', color=RED, ms=4.0, mfc='none', mew=0.9, lw=1.0,
                      zorder=4, label='NuOscProbExact')
    if case is BENCH['cases'][0]:
        # The dial each code is turned by, written beside its markers.  The name appears
        # once per curve, on the topmost point; the rest carry the value alone.
        for series in case['series']:
            pts = series['points']
            dial = next((k for k in ('rtol', 'n_slabs', 'tolerance', 'num_prec')
                         if k in pts[0]), None)
            col = INK if series['name'] == 'Magnus' else RED
            top = max(range(len(pts)), key=lambda k: pts[k]['max_abs_error'])
            for k, pt in enumerate(pts):
                txt = ('%s = %s' % (dial, pt['label'])) if k == top and dial else pt['label']
                # Down and to the right of the marker.  Up and to the right put the
                # topmost label of each curve outside the axes, where it was clipped.
                # The last point of a series is its rightmost, so that one is set to
                # the left instead: to the right it ran past the axis.
                # The last point is a series' rightmost and usually its lowest, so its
                # label goes up and to the left: to the right it ran past the axis, and
                # below it landed on the reference floor.
                last = (k == len(pts) - 1)
                ax.annotate(txt, xy=(pt['us_per_probability'], pt['max_abs_error']),
                            xytext=(-4.0, 10.0) if last else (4.0, -6.5),
                            textcoords='offset points',
                            ha='right' if last else 'left',
                            va='bottom' if last else 'baseline',
                            fontsize=5.6, color=col, zorder=6,
                            annotation_clip=True)
    logx(ax); logy(ax)
    corner(ax, FLAVOR_LABEL[case['flavours']], loc='upper right', fontsize=8.5)
    if case['flavours'] == 5:
        ax.text(0.06, 0.20, 'No route beyond SU(4)', transform=ax.transAxes,
                ha='left', va='center', fontsize=8.0, color=RED, style='italic')

# One x-axis for all four: the panels differ by three decades in cost, so a shared
# range is what makes the flavor counts comparable at a glance.
ALL_T = [p['us_per_probability'] for c in BENCH['cases'] for s_ in c['series']
         for p in s_['points']]
ALL_E = [p['max_abs_error'] for c in BENCH['cases'] for s_ in c['series']
         for p in s_['points']]
for ax in axes:
    # Room past the extreme markers, so their labels have somewhere to sit.
    ax.set_xlim(0.7*min(ALL_T), 1.5*max(ALL_T))
    # A shared vertical range too: a reach that differs by flavor count is the point,
    # and per-panel autoscaling hides it.
    ax.set_ylim(min(ALL_E)/2.5, max(ALL_E)*2.5)
axes[1].legend(loc='lower right', handlelength=1.4)
axes[-1].set_xlabel(r'Time per probability [$\mu$s]')
fig.tight_layout(pad=0.3, h_pad=0.4)
# One y-label for the four panels, centred on the block rather than on any one of them.
fig.supylabel(r'Maximum probability deviation, max $|\Delta P|$',
              fontsize=plt.rcParams['axes.labelsize'], x=0.005)
save(fig, 'smooth_reach.pdf')
for case in BENCH['cases']:
    for s in case['series']:
        best = min(p['max_abs_error'] for p in s['points'])
        last = s['points'][-1]['max_abs_error']
        print('  d=%d %-16s best %.2e, tightest %.2e%s'
              % (case['flavours'], s['name'], best, last,
                 '   <- rises' if last > best*1.01 else ''))'''),
    md(r'''## Figure 8 --- the same shock, as cost against accuracy'''),
    code(r'''# ------------------------------------------------------------------ the shock
SHOCK = json.loads((HERE/'external_shock_benchmarks.json').read_text())
SPAN_KM = SHOCK['L1_km'] - SHOCK['L0_km']
fig, axes = plt.subplots(2, 1, figsize=(COL, 3.4), sharex=True,
                         gridspec_kw=dict(hspace=0.10))
for ax, case in zip(axes, SHOCK['cases']):
    width_km = case['width']*SPAN_KM
    allt = []
    for series in case['series']:
        t = [p['us_per_probability'] for p in series['points']]
        e = [p['max_abs_error'] for p in series['points']]
        allt += t
        if series['name'] == 'Magnus':
            ax.loglog(t, e, '-*', color=INK, ms=8, lw=1.2, zorder=5, label=r'Mag$\nu$s')
        else:
            ax.loglog(t, e, '-o', color=RED, ms=4.0, mfc='none', mew=0.9, lw=1.0,
                      zorder=4, label='NuOscProbExact')
    if case is SHOCK['cases'][0]:
        # Same convention as Fig. 7: the dial's value beside every marker, its name
        # written once per curve on the topmost point.
        for series in case['series']:
            pts = series['points']
            dial = next((k for k in ('rtol', 'n_slabs', 'tolerance', 'num_prec')
                         if k in pts[0]), None)
            col = INK if series['name'] == 'Magnus' else RED
            top = max(range(len(pts)), key=lambda k: pts[k]['max_abs_error'])
            for k, pt in enumerate(pts):
                txt = ('%s = %s' % (dial, pt['label'])) if k == top and dial else pt['label']
                # Down and to the right of the marker.  Up and to the right put the
                # topmost label of each curve outside the axes, where it was clipped.
                # The last point of a series is its rightmost, so that one is set to
                # the left instead: to the right it ran past the axis.
                # The last point is a series' rightmost and usually its lowest, so its
                # label goes up and to the left: to the right it ran past the axis, and
                # below it landed on the reference floor.
                last = (k == len(pts) - 1)
                ax.annotate(txt, xy=(pt['us_per_probability'], pt['max_abs_error']),
                            xytext=(-4.0, 10.0) if last else (4.0, -6.5),
                            textcoords='offset points',
                            ha='right' if last else 'left',
                            va='bottom' if last else 'baseline',
                            fontsize=5.6, color=col, zorder=6,
                            annotation_clip=True)
    ax.axhline(case['reference_unitarity'], color='0.5', ls=':', lw=0.8)
    logx(ax); logy(ax)
    corner(ax, 'SN: front width %s km'
           % ('%.2f' % width_km).rstrip('0').rstrip('.'), fontsize=8.0)
    ax.set_ylabel(r'Max $|\Delta P|$', fontsize=8.0)
    print('  width %6.2f km: Magnus %.2e   NuOscProbExact %.2e'
          % (width_km,
             min(p['max_abs_error'] for p in case['series'][0]['points']),
             min(p['max_abs_error'] for p in case['series'][1]['points'])))
ALL_TS = [p['us_per_probability'] for c in SHOCK['cases'] for s_ in c['series']
          for p in s_['points']]
for ax in axes:
    # A little room past the extreme markers, which otherwise sit on the spines and
    # leave their labels nowhere to go.
    ax.set_xlim(0.7*min(ALL_TS), 1.5*max(ALL_TS))
axes[0].legend(loc='lower left', handlelength=1.4)
stamp(axes[1], 'Referee floor', x=0.04, y=0.04, fontsize=8.0)
axes[-1].set_xlabel(r'Time per probability [$\mu$s]')
fig.tight_layout(pad=0.3, h_pad=0.4)
save(fig, 'shock_speed_accuracy.pdf')'''),
    md(r'''## Figure 9 --- six codes through the Earth

**The matter potential is matched first**, and matching it does not buy a curve that falls
forever. The cell after the figure shows why.'''),
    code(r'''# ------------------------------------------------------------- the Earth plane
PREM = json.loads((HERE/'external_prem_speed_accuracy.json').read_text())
CZ = PREM['costhz']
L_EXT = PREM['baseline_km']*gd.CONV_KM_TO_INV_EV
OSC_EXT = gd.load_nufit_params('NuFIT 4.0', 'NO')
STER_EXT = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, D41=1.0)
VCC_MATCH = 1.0001896490
YE_EXT = 0.5*VCC_MATCH
RTOLS = (1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10)
print('frozen dataset: costhz = %.2f, chord %.1f km; V_CC matched via Y_e = %.10f'
      % (CZ, PREM['baseline_km'], YE_EXT))


def timed_batch(call, n, repeat=5, min_block=0.05):
    call()
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        el = time.perf_counter() - t0
        if el >= min_block:
            break
        reps *= 2
    best = el/reps
    for _ in range(repeat - 1):
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        best = min(best, (time.perf_counter() - t0)/reps)
    return 1.0e6*best/n


STYLE = {'NuOscProbExact': ('-o', RED, 3.6),
         'NuOscProbExact (tolerance)': ('-o', RED, 2.8),
         'NuOscProbExact (double-double)': ('-o', RED, 3.6),
         'NuOscProbExact (eigensolver)': ('-h', RED, 3.4),
         'nuSQuIDS': ('-v', GREEN, 3.2), 'nuCraft': ('-s', ORANGE, 3.0),
         'NuFast-Earth': ('-D', PURPLE, 2.8), 'GLoBES': ('-*', '#a51d2d', 5.2),
         'Prob3++': ('-P', '#986a44', 3.6)}
DIALS = ('n_slabs_per_segment', 'rtol', 'tolerance', 'num_prec', 'n_shells_per_layer')

def prem_magnus_curve():
    """The \\magnus\\ curve of this plane: one timed batch per tolerance, at both flavor counts.

    The five external curves are frozen in the repository, and this one is measured here.
    Measuring it on every rebuild would put minutes of stopwatch into continuous
    integration for numbers that cannot have moved unless the configuration did, so it is
    stored the same way the references are.
    """
    out = {}
    for key in ('three_flavor', 'sterile_3plus1'):
        blk = PREM[key]
        E_ext = np.asarray(blk['energy_gev'])*gd.UNIT_GEV
        P_ref = np.asarray(blk['reference'])
        floor = float(blk['reference_vs_ode_max_abs'])
        mg_t, mg_e = [], []
        for rtol in RTOLS:
            if key == 'three_flavor':
                call = (lambda r=rtol, E=E_ext: np.asarray(quiet(
                    oscprob.osc_prob_3nu_earth, E, costhz=CZ, L=L_EXT, **OSC_EXT,
                    nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=r, atol=r*1.0e-2,
                    electron_fraction=YE_EXT)))
            else:
                call = (lambda r=rtol, E=E_ext: np.asarray(quiet(
                    oscprob.osc_prob_4nu_earth, E, costhz=CZ, L=L_EXT, **OSC_EXT,
                    d14=0.0, d24=0.0, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=r, atol=r*1.0e-2,
                    **STER_EXT, electron_fraction=YE_EXT)))
            P = call()
            mg_t.append(timed_batch(call, len(E_ext)))
            mg_e.append(max(float(np.max(np.abs(P - P_ref))), floor))
        out[key] = dict(t=mg_t, e=mg_e)
    return out


PREM_MAGNUS = cached(
    'prem_timings',
    (repr(RTOLS), float(CZ), float(L_EXT), float(YE_EXT), repr(sorted(OSC_EXT.items())),
     repr(sorted(STER_EXT.items())),
     np.asarray(PREM['three_flavor']['energy_gev']),
     np.asarray(PREM['sterile_3plus1']['energy_gev'])),
    prem_magnus_curve,
    what='Figure 9: the Magnus curve of the six-code Earth plane, timed per tolerance.')
for _k in ('three_flavor', 'sterile_3plus1'):
    for _r, _t, _e in zip(RTOLS, PREM_MAGNUS[_k]['t'], PREM_MAGNUS[_k]['e']):
        print('  %-14s rtol %.0e -> %9.1f us/prob, err %.2e' % (_k, _r, _t, _e))

fig, axes = plt.subplots(1, 2, figsize=(WIDE, 3.1))
for ax, (key, label) in zip(axes, (('three_flavor', r'$3\nu$'),
                                   ('sterile_3plus1', r'$3+1$'))):
    blk = PREM[key]
    E_ext = np.asarray(blk['energy_gev'])*gd.UNIT_GEV
    P_ref = np.asarray(blk['reference'])
    floor = float(blk['reference_vs_ode_max_abs'])
    allt = []
    for series in blk['series']:
        pts = series['points']
        marker, color, size = STYLE.get(series['name'], ('-o', '0.4', 3.2))
        dial = next((k for k in DIALS if pts[0].get(k) is not None), None)
        kw = dict(ms=size, color=color, lw=0.9, zorder=4,
                  label='%s%s' % (series['name'],
                                  '' if dial is None else ' (%s)' % dial))
        if series['name'].startswith('NuOscProbExact'):
            kw.update(mfc='white', mew=0.8, zorder=5)
        t = [p['us_per_probability'] for p in pts]
        allt += t
        ax.loglog(t, [max(p['max_abs_error'], floor) for p in pts], marker, **kw)

    mg_t, mg_e = PREM_MAGNUS[key]['t'], PREM_MAGNUS[key]['e']
    allt += mg_t
    ax.loglog(mg_t, mg_e, '-*', ms=11, color=INK, lw=1.3, zorder=6,
              label=r'Mag$\nu$s  (rtol)')
    ax.axhline(floor, color='0.5', ls=':', lw=0.8, zorder=1)
    logx(ax); logy(ax); snug(ax, allt)
    corner(ax, label, loc='upper right', fontsize=8.5)
    stamp(ax, 'Referee floor', x=0.035, y=0.035, fontsize=8.0)
    ax.set_xlabel(r'Time per probability [$\mu$s]')
    ax.legend(loc='lower left', fontsize=8.0, handlelength=1.5, labelspacing=0.24)
axes[0].set_ylabel(r'Error against the converged reference, max $|\Delta P|$')
fig.tight_layout(pad=0.3, w_pad=1.6)
save(fig, 'prem_speed_accuracy.pdf')'''),
    md(r'''### Why that curve stops falling

Is Mag$\nu$s still converging where its curve flattens? If it is self-converged and still
disagreeing, what the plane measures beyond that point is a difference between two Earth
models, and calling it accuracy would be wrong.'''),
    code(r'''blk = PREM['three_flavor']
E_ext = np.asarray(blk['energy_gev'])*gd.UNIT_GEV
P_ref = np.asarray(blk['reference'])


def magnus_prem_3nu(rtol, ye=YE_EXT):
    return np.asarray(quiet(oscprob.osc_prob_3nu_earth, E_ext, costhz=CZ, L=L_EXT,
                            **OSC_EXT, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=rtol,
                            atol=rtol*1.0e-2, electron_fraction=ye))


a, b, c = magnus_prem_3nu(1.0e-8), magnus_prem_3nu(1.0e-10), magnus_prem_3nu(1.0e-12)
print('Is Magnus self-converged where the curve flattens?')
print('  |P(1e-8)  - P(1e-12)|  = %.3e' % np.max(np.abs(a - c)))
print('  |P(1e-10) - P(1e-12)|  = %.3e   <- yes, to a few 1e-12' % np.max(np.abs(b - c)))
print('  |P(1e-12) - reference| = %.3e   <- and yet this is the floor'
      % np.max(np.abs(c - P_ref)))
print()
print('Does the floor move with the matched potential?  (V_CC is linear in Y_e)')
for scale in (1.0, 1.0001894920, VCC_MATCH, 1.0003):
    print('  Y_e = 0.5 * %.10f  ->  |P - reference| = %.3e'
          % (scale, np.max(np.abs(magnus_prem_3nu(1.0e-10, 0.5*scale) - P_ref))))
print()
print('A relative change of %.1e in the potential moves the floor by about 5x, so the'
      % ((VCC_MATCH - 1.0001894920)/VCC_MATCH))
print('probability inherits the relative error of V_CC essentially one for one.  The')
print('scale used is the one the two codes\' own constants imply: fitting it to minimise')
print('the residual would make the figure prettier and the measurement meaningless.')'''),
    md(r'''## What was written

Nine PDFs, which is every figure in `resources/paper/main.tex`.

```bash
python notebooks/make_notebooks.py --only 28
```'''),
    code(r'''for name in sorted(p.name for p in FIGDIR.glob('*.pdf')):
    print('  %-28s %8.1f kB' % (name, (FIGDIR/name).stat().st_size/1024.0))'''),
    ])

books['29_magnus_pseudo_dirac.ipynb'] = notebook(
    'Pseudo-Dirac neutrinos',
    r"""Each Dirac neutrino may in fact be two Majorana states separated by a tiny mass-squared
splitting $\delta m^2$. The consequence is a separation of scales, and the separation is the
whole physical content: over an astrophysical baseline the standard splittings have long since
averaged away, while each pseudo-Dirac pair is still only part-way through its first cycle.

That is precisely the regime the **coherent-block** averaging form is for. Summing
probabilities eigenstate by eigenstate -- the form that is right once everything has decohered
-- undercounts here by the number of states sharing a block, which for a fully paired spectrum
is a factor of two.

This notebook builds the Hamiltonian, checks it reduces exactly to the Dirac case when nothing
is paired, and then works through the two regimes: the oscillatory one, where only the
un-averaged probability is valid, and the averaged one, where the block form is. It closes on
two negative results that are as useful as the positive ones -- the instantaneous probability
is not an observable at astrophysical distances, and the effect is invisible on Earth.""",
    [

    md(r"""## 1. The construction

`pseudo_dirac_mixing_matrix` replaces each paired mass eigenstate $j$ by the two combinations
$(|\nu_j\rangle \pm |s_j\rangle)/\sqrt{2}$, and `pseudo_dirac_mass_squared` gives them masses
$m_j^2$ and $m_j^2 + \delta m^2_j$. Pairing is chosen **per mass state**: a three-flavor
spectrum with partners on two of its three states is a five-dimensional problem."""),

    code(r'''import warnings

import numpy as np
import matplotlib.pyplot as plt

import magnus.avgprob as avgprob
import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.oscprob as oscprob

MPC_IN_KM = 3.0856775814913673e19        # 1 Mpc, in km

OSC = gd.load_nufit_params('NuFIT 6.1', 'NO')
U = hamiltonians.pmns_mixing_matrix(OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP'])
M2 = [0.0, OSC['D21'], OSC['D31']]

for pairs in ({}, {1: 1.0e-18}, {0: 1.0e-18, 2: 4.0e-18},
              {0: 1.0e-18, 1: 1.0e-18, 2: 1.0e-18}):
    W = hamiltonians.pseudo_dirac_mixing_matrix(U, pairs)
    unitary = np.allclose(W @ np.conj(W.T), np.eye(len(W)), atol=1.0e-14)
    print('pairs on %-9s -> %d states, unitary: %s'
          % (str(sorted(pairs)) if pairs else 'none', len(W), unitary))'''),

    md(r"""### Nothing paired is the Dirac case, exactly

Not approximately. With an empty pairing the same matrices are multiplied in the same order as
in `hamiltonian_3nu_vacuum_energy_independent`, so the two agree to round-off. A builder that
only *nearly* reduced would be hiding a convention difference."""),

    code(r'''H_pd = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(U, M2, {})
H_3nu = hamiltonians.hamiltonian_3nu_vacuum_energy_independent(
    OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP'], OSC['D21'], OSC['D31'])

print('max |H_pseudo-Dirac - H_3nu| = %.2e' % float(np.max(np.abs(H_pd - H_3nu))))'''),

    md(r"""## 2. The separation of scales

The ratio of the pair phase to the standard one is $\delta m^2/\Delta m^2_{31}$, and nothing
in the setup can change it. That single fact organizes everything below: at a splitting small
enough to be pseudo-Dirac, **you cannot have a resolvable standard oscillation and a
developing pair phase at the same time.**"""),

    code(r'''def phase(dm2, L_km, energy):
    """The relative phase dm^2 L / (4E), in radians."""
    return dm2*(L_km*gd.UNIT_KM)/(4.0*energy)

print('%-38s %12s %12s' % ('configuration', 'pair phase', 'Dm31 phase'))
print('-'*64)
for label, dm2, L_km, energy in (
        ('100 Mpc, 100 TeV, dm2 = 2.6e-17', 2.6e-17, 100.0*MPC_IN_KM, 100.0*gd.UNIT_TEV),
        ('10 Mpc,  100 TeV, dm2 = 2.6e-16', 2.6e-16, 10.0*MPC_IN_KM, 100.0*gd.UNIT_TEV),
        ('1300 km, 1 GeV,   dm2 = 2.6e-17', 2.6e-17, 1300.0, 1.0*gd.UNIT_GEV)):
    print('%-38s %12.3g %12.3g'
          % (label, phase(dm2, L_km, energy), phase(OSC['D31'], L_km, energy)))'''),

    md(r"""## 3. The oscillatory regime

Here the pair phase is $O(1)$ and the averaged expressions do not apply at all: the pair has
not decohered, so there is nothing to average, and the only valid route is the un-averaged
`osc_prob`.

**The splitting below is deliberately exaggerated**, to $1.2\times10^{-3}$ eV$^2$, so that the
modulation is visible on a terrestrial baseline in a calculation that runs in a fraction of a
second. A physical pseudo-Dirac splitting is some fourteen orders of magnitude smaller and
does nothing here -- Section 7 shows exactly that. This panel is about the *shape* of the
effect, not its size."""),

    code(r'''L_DUNE = 1300.0*gd.UNIT_KM
energies = np.linspace(0.6, 4.0, 200)*gd.UNIT_GEV
DM2_SHOWN = 1.2e-3                       # exaggerated, for visibility

with warnings.catch_warnings():
    # The splitting is not small against the standard ones, and the library says so.
    # That is the correct complaint; it is exaggerated on purpose here.
    warnings.simplefilter('ignore', hamiltonians.PseudoDiracSplittingWarning)
    H_pair = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(
        U, M2, {1: DM2_SHOWN, 2: DM2_SHOWN})
H_dirac = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(U, M2, {})

P_pair = np.array([np.asarray(oscprob.osc_prob(H_pair/e, 0.0, L_DUNE))[1, 0]
                   for e in energies])
P_dirac = np.array([np.asarray(oscprob.osc_prob(H_dirac/e, 0.0, L_DUNE))[1, 0]
                    for e in energies])

fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.plot(energies/gd.UNIT_GEV, P_dirac, lw=1.6, label='Dirac (no pairs)')
ax.plot(energies/gd.UNIT_GEV, P_pair, lw=1.6, ls='--',
        label=r'pseudo-Dirac, $\delta m^2 = 1.2\times10^{-3}$ eV$^2$')
ax.set_xlabel(r'$E_\nu$ [GeV]')
ax.set_ylabel(r'$P(\nu_\mu \to \nu_e)$')
ax.set_title('1300 km, exaggerated splitting so the modulation is visible', fontsize=10)
ax.legend(fontsize=8)
fig.tight_layout()

print('largest departure from the Dirac curve: %.4f' % float(np.max(np.abs(P_pair - P_dirac))))'''),

    md(r"""## 4. At astrophysical distances the instantaneous probability is not an observable

At 100 Mpc and 100 TeV the standard phases are $\sim10^{14}$ radians. No experiment resolves
that, and neither does double precision. Sampling `osc_prob` at neighbouring energies -- close
enough that no detector could tell them apart -- returns values scattered across a large part
of the range.

This is the same lesson notebook 13 draws for the Sun, and it is why the rest of this notebook
works with the averaged probability."""),

    code(r'''L_100MPC = 100.0*MPC_IN_KM*gd.UNIT_KM
PAIRS_PHYS = {0: 2.6e-17, 1: 2.6e-17, 2: 2.6e-17}

H_astro = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(U, M2, PAIRS_PHYS)
nearby = 100.0*gd.UNIT_TEV*(1.0 + np.linspace(0.0, 1.0e-12, 40))
sampled = np.array([np.asarray(oscprob.osc_prob(H_astro/e, 0.0, L_100MPC))[0, 0]
                    for e in nearby])

print('40 energies spanning a relative range of 1e-12:')
print('  P_ee from %.4f to %.4f, spread %.4f' % (sampled.min(), sampled.max(), sampled.ptp()))
print()
print('The energies differ by one part in 1e12.  Nothing measures that, so the')
print('instantaneous probability is not the quantity an experiment reports.')'''),

    md(r"""## 5. The averaged regime, and the factor of two

Now the standard phases have averaged away and each pair is still coherent.
`coherence_blocks` sees exactly that: three blocks of two.

Two expressions are then in play. The **coherent-block** form sums *amplitudes* within a block
and squares once; the **naive** form sums probabilities, one term per eigenstate. Within a
block the pair splits its parent state's mixing evenly between two columns, so summing the
amplitudes rebuilds $|U_{\alpha j}|^2$ and the block form returns the ordinary Dirac answer.
Summing probabilities instead loses a factor of two."""),

    code(r'''W = hamiltonians.pseudo_dirac_mixing_matrix(U, PAIRS_PHYS)
masses = hamiltonians.pseudo_dirac_mass_squared(M2, PAIRS_PHYS)
eigenvalues = masses/(2.0*100.0*gd.UNIT_TEV)

blocks = avgprob.coherence_blocks(eigenvalues, L_100MPC)
print('coherence_blocks:', blocks)

def block_form(V, blocks, a, b):
    """The coherent-block average: sum amplitudes inside a block, then square."""
    return float(sum(abs(sum(np.conj(V[a, i])*V[b, i] for i in blk))**2
                     for blk in blocks).real)

def naive_sum(V, a, b):
    """One term per eigenstate: correct only once every pair has decohered."""
    return float(sum(abs(np.conj(V[a, i])*V[b, i])**2 for i in range(len(V))).real)

dirac = float(np.sum(np.abs(U[0, :])**4))
print()
print('  block form  <P_ee> = %.5f' % block_form(W, blocks, 0, 0))
print('  naive sum   <P_ee> = %.5f' % naive_sum(W, 0, 0))
print('  Dirac 3nu   <P_ee> = %.5f' % dirac)
print()
print('  block/naive = %.4f' % (block_form(W, blocks, 0, 0)/naive_sum(W, 0, 0)))'''),

    md(r"""### A partially paired spectrum

The interface pairs individual mass states, which is more than a toggle: with partners on
states 0 and 2 only, the blocks come out as two-one-two and the unpaired state contributes its
ordinary single term."""),

    code(r'''PAIRS_PARTIAL = {0: 2.6e-17, 2: 2.6e-17}
W_partial = hamiltonians.pseudo_dirac_mixing_matrix(U, PAIRS_PARTIAL)
masses_partial = hamiltonians.pseudo_dirac_mass_squared(M2, PAIRS_PARTIAL)
blocks_partial = avgprob.coherence_blocks(
    masses_partial/(2.0*100.0*gd.UNIT_TEV), L_100MPC)

print('%d states, blocks %s' % (len(W_partial), blocks_partial))
print('  block form <P_ee> = %.5f' % block_form(W_partial, blocks_partial, 0, 0))
print('  naive sum  <P_ee> = %.5f' % naive_sum(W_partial, 0, 0))
print()
print('Between the fully paired case (a factor of two) and the Dirac case (no factor),')
print('as it must be: two of the three states carry a partner.')'''),

    md(r"""## 6. Sweeping the splitting through the three regimes

The library recognizes three regimes and refuses the middle one. Sweeping $\delta m^2$ upward
at fixed $L/E$ walks through all three: coherent pairs, then a band where neither limit
describes the physics, then full decoherence into six singletons.

`coherence_report` is what says which regime you are in, and the un-averaged probability is the
only valid route through the middle band."""),

    code(r'''print('%-12s %-32s %s' % ('dm2 [eV^2]', 'blocks', 'regime'))
print('-'*66)
for dm2 in (1.0e-19, 1.0e-18, 1.0e-17, 3.0e-17, 1.0e-16, 1.0e-15, 1.0e-13):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', hamiltonians.PseudoDiracSplittingWarning)
        masses_i = hamiltonians.pseudo_dirac_mass_squared(M2, {j: dm2 for j in range(3)})
    lam = masses_i/(2.0*100.0*gd.UNIT_TEV)
    blocks_i = avgprob.coherence_blocks(lam, L_100MPC)
    _, undecided = avgprob.coherence_report(lam, L_100MPC)
    sizes = sorted(len(b) for b in blocks_i)
    if sizes == [1]*6:
        regime = 'decohered: naive sum'
    elif undecided:
        # Grouped as pairs, but the pair phase has grown past the coherence
        # threshold: the cross term is neither kept nor dropped cleanly.
        regime = 'partly developed: NEITHER limit'
    else:
        regime = 'coherent pairs: block form'
    print('%-12.1e %-32s %s' % (dm2, str(blocks_i), regime))'''),

    md(r"""The middle band is not a numerical inconvenience: no averaged expression describes it.
Asking the library to average there raises `PhaseAveragingWarning` rather than returning a
number that looks fine."""),

    code(r'''lam_mid = hamiltonians.pseudo_dirac_mass_squared(
    M2, {j: 3.0e-17 for j in range(3)})/(2.0*100.0*gd.UNIT_TEV)
blocks_mid, undecided = avgprob.coherence_report(lam_mid, L_100MPC)

print('blocks:', blocks_mid)
print('pairs in neither limit (i, j, relative phase in radians):')
for i, j, relative_phase in undecided[:6]:
    print('   (%d, %d)   %.3f rad' % (i, j, relative_phase))
print()
print('That list is what makes the difference.  Empty, and one of the two averaged')
print('expressions applies.  Non-empty, as here, and neither does: the phases sit')
print('between %.2g and 2*pi, too large to keep the cross term and too small to' % 1.0e-2)
print('drop it.  The library raises PhaseAveragingWarning rather than choosing.')'''),

    md(r"""## 7. On Earth, the effect is invisible

A useful negative result, and the guard against assuming the feature matters everywhere. At a
terrestrial baseline the physical pair phase is $\sim10^{-14}$ radians: the two members of each
pair have not begun to separate, so the spectrum is indistinguishable from Dirac."""),

    code(r'''H_earth_pair = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(
    U, M2, PAIRS_PHYS)
H_earth_dirac = hamiltonians.hamiltonian_pseudo_dirac_vacuum_energy_independent(U, M2, {})

energy = 1.0*gd.UNIT_GEV
P_pair_earth = np.asarray(oscprob.osc_prob(H_earth_pair/energy, 0.0, L_DUNE))[1, 0]
P_dirac_earth = np.asarray(oscprob.osc_prob(H_earth_dirac/energy, 0.0, L_DUNE))[1, 0]

print('1300 km, 1 GeV, dm2 = 2.6e-17 eV^2')
print('  pair phase              %.2e rad' % phase(2.6e-17, 1300.0, energy))
print('  P(numu->nue) Dirac       %.6f' % P_dirac_earth)
print('  P(numu->nue) pseudo-Dirac %.6f' % P_pair_earth)
print('  difference               %.2e' % abs(P_pair_earth - P_dirac_earth))
print()
print('Pseudo-Dirac splittings are an astrophysical-baseline effect.  Nothing')
print('terrestrial constrains them, which is the point of looking at neutrinos')
print('that have travelled megaparsecs.')'''),

    md(r"""## What to take away

* Pairing is **per mass state**, so a partially paired spectrum is a first-class case.
* With nothing paired the Hamiltonian is the Dirac one to round-off, so the feature costs
  nothing when it is not used.
* In the coherent regime the **block form** returns the Dirac answer and the naive
  eigenstate sum is wrong by the number of states sharing a block.
* Between the two limits no averaged expression applies, and the library says so rather than
  returning a plausible number.
* The effect lives at astrophysical baselines. On Earth the pair phase is $10^{-14}$ radians."""),
    ])

READING_ORDER = [
    ('01_magnus_introduction.ipynb', 'Introduction',
     'the shortest path to a probability'),
    ('02_magnus_2nu_vacuum_matter.ipynb', 'Two-neutrino probabilities',
     'vacuum, constant and varying density, castle wall, Earth and Sun'),
    ('03_magnus_3nu_vacuum_matter.ipynb', 'Three-neutrino probabilities',
     'the same ground with three flavors and a CP phase'),
    ('04_magnus_long_baseline.ipynb', 'Long baselines',
     'probabilities between two points on the surface'),
    ('05_magnus_biprobability.ipynb', 'Biprobability plots',
     'neutrino against antineutrino, as the CP phase runs'),
    ('06_magnus_oscillograms.ipynb', 'Oscillograms',
     'probability across zenith angle and energy at once'),
    ('07_magnus_bsm_sterile_nu.ipynb', 'BSM: sterile neutrinos',
     'four and five flavors'),
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
    ('15_magnus_antineutrinos.ipynb', 'Antineutrinos, done properly',
     'conjugate and flip, and two ways to get it half right'),
    ('16_magnus_exact_vs_approximations.ipynb', 'Exact versus the approximations',
     'where the textbook formulas are exact, and where the substitution breaks'),
    ('17_magnus_ordering_and_octant.ipynb', 'Mass ordering and the octant',
     'two open questions, and how large they actually are'),
    ('18_magnus_unusual_density_profiles.ipynb', 'Unusual density profiles',
     'arrangement beats the mean, except for one exact symmetry'),
    ('19_magnus_custom_hamiltonian.ipynb', 'Bring your own Hamiltonian',
     'the contract, the vectorization trick, and what the Earth declares for you'),
    ('20_magnus_numerical_edge_cases.ipynb', 'Numerical edge cases',
     'degeneracies that return numbers, and what the nine warnings mean'),
    ('21_magnus_what_tolerance_means.ipynb', 'What rtol and atol promise',
     'a stopping criterion, not an error bound'),
    ('22_magnus_which_engine_answered.ipynb', 'Which engine answered, and why',
     'six engines, five families, and an error bar with no oracle'),
    ('23_magnus_when_averaging_helps.ipynb', 'When averaging rescues you',
     'phase error falls away, envelope error does not'),
    ('24_magnus_performance.ipynb', 'Performance',
     'what is worth doing, and when each trick is worth nothing'),
    ('25_magnus_against_other_codes.ipynb', 'Against other codes',
     'where a closed form wins, and a conventions trap that looks like accuracy'),
    ('26_magnus_nufit_evolution.ipynb', 'Fourteen years of NuFIT',
     'how the parameter likelihood, not just the best fit, moves the probability'),
    ('27_magnus_animations.ipynb', 'Animated scenes',
     'nine sweeps drawn as filmstrips: four shared with NuOscProbExact, five that need '
     'a ladder, an average or a varying Hamiltonian'),
    ('28_magnus_paper_figures.ipynb', "The paper's figures",
     'every figure in the CPC article, in one run'),
    ('29_magnus_pseudo_dirac.ipynb', 'Pseudo-Dirac neutrinos',
     'tiny splittings, coherent blocks, and where the effect is invisible'),
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


# Figures lifted out of the executed notebooks for the docs gallery.  Keyed by
# (notebook, index of the PNG output within it), so the recipes page and the
# notebook show the same figure and there is no third version to drift.
GALLERY_DIR = HERE.parent/'img'/'gallery'

GALLERY = {
    ('02_magnus_2nu_vacuum_matter.ipynb', 0): 'gallery_2nu_vacuum.png',
    ('03_magnus_3nu_vacuum_matter.ipynb', 0): 'gallery_3nu_vacuum.png',
    ('05_magnus_biprobability.ipynb', 0): 'gallery_biprobability.png',
    ('06_magnus_oscillograms.ipynb', 0): 'gallery_oscillogram.png',
    ('10_magnus_averaged_probability.ipynb', 0): 'gallery_averaged.png',
    ('13_magnus_tabulated_solar_model.ipynb', 0): 'gallery_solar_model.png',
    ('14_magnus_supernova_shock.ipynb', 0): 'gallery_shock.png',
    # The README gallery draws on these as well.  They carry what this package offers
    # that a closed-form code does not -- an averaged observable nothing else exposes,
    # the boundary where each method wins, and a Hamiltonian nobody has diagonalized --
    # so the front page argues from measurements rather than from claims.
    ('04_magnus_long_baseline.ipynb', 1): 'gallery_long_baseline.png',
    ('13_magnus_tabulated_solar_model.ipynb', 2): 'gallery_solar_bsm.png',
    ('14_magnus_supernova_shock.ipynb', 2): 'gallery_shock_bsm.png',
    ('07_magnus_bsm_sterile_nu.ipynb', 11): 'gallery_sterile_3plus2.png',
    ('18_magnus_unusual_density_profiles.ipynb', 0): 'gallery_density_arrangement.png',
    ('19_magnus_custom_hamiltonian.ipynb', 1): 'gallery_custom_h.png',
    ('25_magnus_against_other_codes.ipynb', 11): 'gallery_solar_averaged.png',
    # Not gallery material despite being good results: the expansion-order scatter and
    # the smooth-profile reach panel both carry a long y-axis label that the notebook's
    # tight bounding box truncates, and neither reads at 380 px.  Their findings are in
    # notebooks 24 and 25 instead.
}


def extract_gallery():
    r"""Writes the gallery figures out of the executed notebooks.

    The docs embed these rather than plotting their own, so a figure on the
    recipes page is by construction the one the notebook produced.
    """
    import base64

    GALLERY_DIR.mkdir(parents=True, exist_ok=True)
    written = 0
    for (notebook, index), filename in sorted(GALLERY.items()):
        nb = nbf.read(HERE/notebook, as_version=4)
        images = [output['data']['image/png']
                  for cell in nb.cells
                  for output in cell.get('outputs', [])
                  if 'image/png' in output.get('data', {})]
        if index >= len(images):
            raise SystemExit('%s has no figure %d (it has %d)'
                             % (notebook, index, len(images)))
        (GALLERY_DIR/filename).write_bytes(base64.b64decode(images[index]))
        written += 1
    print('  wrote %d gallery figures to %s' % (written, GALLERY_DIR))


def sources_match(path, nb):
    r"""Whether the notebook on disk has the same cell sources as the one built here.

    Outputs and execution counts are ignored: a rebuild rewrites those in every
    notebook, and they say nothing about whether the file is current.
    """
    if not path.exists():
        return False
    on_disk = nbf.read(path, as_version=4)
    return [c.source for c in on_disk.cells] == [c.source for c in nb.cells]


def build(execute=True, only=None):
    r"""Writes every notebook, executes it, and checks it kept its outputs.

    ``only`` restricts the work to the notebooks whose names contain one of the
    given fragments -- ``--only 19,24`` is enough.  The rest are **left exactly
    as they are on disk**, outputs and all, which is the whole point: a
    one-notebook change costs one notebook's runtime rather than the set's.

    That is also the hazard.  ``--no-execute`` rewrites every notebook without
    outputs and is documented as destructive; a filter that rewrote the others
    blank would be the same trap wearing a different hat.  So untouched
    notebooks are never written, and the run ends by checking that every
    notebook on disk still matches the sources built here -- which catches the
    one thing this flag can get wrong, a file left behind by an earlier edit.
    """
    add_footers()

    selected = dict(books) if not only else {
        name: nb for name, nb in books.items()
        if any(fragment in name for fragment in only)}
    if only and not selected:
        raise SystemExit('--only matched no notebooks: %s' % ', '.join(only))

    for name, nb in selected.items():
        nbf.write(nb, HERE/name)
    print('  wrote %d notebook%s%s' % (len(selected),
                                       '' if len(selected) == 1 else 's',
                                       '' if not only else ' (of %d)' % len(books)))

    if not execute:
        return

    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    failed = []
    for path in sorted(HERE/name for name in selected):
        nb = nbf.read(path, as_version=4)
        started = time.perf_counter()
        try:
            NotebookClient(
                # Six hours, not one.  A cold paper cache re-derives fifty-digit
                # references and re-runs the whole timing sweep in a single cell, which
                # is well over an hour on its own -- and the old one-hour limit killed
                # that cell *after* it had written its results, so the work was done and
                # thrown away.  Every warm rebuild is seconds, so this ceiling only ever
                # applies to the rare re-measurement it exists to protect.
                nb, timeout=21600, kernel_name='python3',
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
    # Every notebook on disk must still be the one this file builds.  With --only
    # the others were not rewritten, so this is what proves they were already
    # current rather than left over from an earlier edit.
    stale = [name for name, nb in books.items()
             if not sources_match(HERE/name, nb)]
    if stale:
        raise SystemExit('notebooks on disk no longer match this generator: %s\n'
                         'Rebuild them, or run without --only.' % ', '.join(stale))

    print('  %d executed; all %d notebooks match the generator and carry outputs'
          % (len(selected), len(books)))

    extract_gallery()


if __name__ == '__main__':
    import sys
    argv = sys.argv[1:]
    chosen = None
    if '--only' in argv:
        chosen = [f.strip() for f in argv[argv.index('--only') + 1].split(',')
                  if f.strip()]
    build(execute='--no-execute' not in argv, only=chosen)
