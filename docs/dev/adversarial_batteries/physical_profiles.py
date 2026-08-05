# -*- coding: utf-8 -*-
r"""Physically-motivated density profiles, for the reachability question.

Three tranches of robustness work produced findings ranging from *wrong by 0.54 in probability*
to *four orders of accuracy left on the table*.  Every catastrophic number came from either a
deliberate adversarial construction (an unmarked density step, a Gaussian narrower than the
probe grid, a sawtooth) or from the fuzzers' **random Fourier sum**, which is a mathematically
convenient way to make a smooth positive function and has no physics in it at all.  Against the
profiles the package actually ships for -- 164 Earth + solar configurations, 42 workloads on
solar and multi-resonance -- the record is clean.

So the robustness case rested on profiles nobody would compute.  This module builds a population
a referee would accept as physically motivated, so the existing instruments can be re-run against
it.  **A negative result is a real result here**: if the physical families show no silent miss,
that is the answer, and it should not be "fixed" by making the profiles more adversarial until
something breaks.

Every profile returns an **electron number density in eV^3**, matching ``harness.py``, so
``harness.vcc_of`` (which passes ``density_is_of_number_of_electrons=True``) is the only route to
a potential.  Positions are in natural units (eV^-1); ``gd.UNIT_KM`` converts.

Where each family sits on the "physically motivated" spectrum -- the distinction §7.3 of the
handover asks to be explicit about -- is recorded per family in ``PROVENANCE`` and repeated in
each builder's docstring:

===================  =====================================================================
family               provenance
===================  =====================================================================
``tabulated``        the package's own analytic profile, resampled and interpolated.  The
                     *shape* is the package's; what is physical is the **user behaviour**
                     being modelled (loading a table and interpolating it).
``bs05``             a real published solar model, used directly.  Highest credibility in
                     the population; the only family that is a measurement rather than a
                     parametrization.
``sn_shock``         an analytic parametrization **fitted to simulations in the published
                     literature**, with the shock radii taken from a named simulation.
``sn_turbulence``    the construction used in the published literature, with the paper's
                     spectrum, damping scale and amplitudes.
``earth_crust``      PREM (a real model) plus **invented** extra crustal layers.  The layers
                     are plausible, not measured; what is physical is that a user with their
                     own crust model would not pass ``t_breakpoints``.
===================  =====================================================================

References
----------
.. [Fogli2003] G. L. Fogli, E. Lisi, A. Mirizzi, D. Montanino, "Analysis of energy- and
   time-dependence of supernova shock effects on neutrino crossing probabilities",
   Phys. Rev. D 68, 033005 (2003) [hep-ph/0304056].  Static progenitor profile
   :math:`\rho_0(x) \simeq 10^{14}\,(x/\mathrm{km})^{-2.4}\ \mathrm{g\,cm^{-3}}`, forward-shock
   jump :math:`\xi = V_+/V_- \simeq 10`, and the rarefaction shape
   :math:`\ln f(x) = [0.28 - 0.69\ln(x_s/\mathrm{km})]\,[\arcsin(1-x/x_s)]^{1.1}`.
.. [Kneller2015] J. P. Kneller, N. V. Kabadi, "Sensitivity of neutrinos to the supernova
   turbulence power spectrum: point source statistics", Phys. Rev. D 92, 013009 (2015)
   [arXiv:1410.5698].  Turbulence as :math:`n_e \to n_e(1+F)` with
   :math:`F(r) = C_\star \tanh(\cdot)\tanh(\cdot)\sum_a \sqrt{V_a}\,[A_a\cos q_a r +
   B_a \sin q_a r]`, inverse-power-law spectrum
   :math:`E(q) = \frac{\alpha-1}{2q_\star}(q_\star/|q|)^\alpha\,\Theta(|q|-q_\star)` with
   :math:`\alpha = 5/3`, damping scale :math:`\lambda_\star = 100` km.  Its Fig. 1 gives the
   three discontinuities of the :math:`M = 10.8\,M_\odot` model of Fischer et al. at
   :math:`t = 3` s post-bounce: reverse shock 1734 km, contact discontinuity 12 348 km,
   forward shock 30 323 km.
.. [Bahcall2005] J. N. Bahcall, A. M. Serenelli, S. Basu, ApJ 621, L85 (2005)
   [astro-ph/0412440].  Model BS2005-AGS,OP, table ``bs05_agsop.dat``.
.. [PREM] A. M. Dziewonski, D. L. Anderson, Phys. Earth Planet. Inter. 25, 297 (1981), as
   already implemented in ``magnus.earth``.
"""

import os

import numpy as np

import harness as H
import magnus.earth as earth
import magnus.globaldefs as gd

# ----------------------------------------------------------------------
# Units.  The package converts g/cm^3 to eV^3 with an average nucleon mass rather than the
# atomic mass unit (matter.num_density_e_func with ratio_number_neutrons_to_protons = 1), so
# that convention is reproduced here exactly: a profile that disagreed with the package's own
# conversion by 0.1% would put every resonance in a slightly different place than the
# package's own solar profile does, for no reason.
# ----------------------------------------------------------------------

KM = gd.UNIT_KM
_MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)


def ne_from_rho_cgs(rho, y_e=0.5):
    r"""Electron number density [eV^3] from a mass density [g cm^-3]."""
    return np.asarray(rho, dtype=float)*gd.UNIT_G_PER_CM3/_MEAN_NUCLEON*y_e


def _smoothstep(u):
    """C^1 ramp from 0 to 1 over u in [0, 1], saturating outside.

    Used for shock fronts of finite width.  A tanh would never actually reach its asymptote,
    which makes "the feature is w wide" untrue by a factor that depends on w; this reaches it.
    """
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    return u*u*(3.0 - 2.0*u)


# ======================================================================
# 2.1  Tabulated density with interpolation kinks
# ======================================================================

def tabulated_ne(n_nodes, l0, l1, kind='linear', base=None):
    r"""The package's own solar exponential, sampled at ``n_nodes`` and interpolated.

    **This is the most common way a real user's profile stops being smooth, and the least
    glamorous.**  A user loads a density table -- a stellar model, a coarse PREM, a simulation
    snapshot -- and interpolates it.  Linear interpolation gives a :math:`C^0`-but-not-:math:`C^1`
    **kink at every node**; a cubic spline is :math:`C^2` but rings.  The previous session
    measured a hand-built kink at 1.448e-02, silently, before the fixes; this asks whether
    ordinary interpolation reproduces that.

    Provenance: the underlying shape is the package's own analytic profile, so nothing here is
    a claim about solar physics.  What is physical is the *user behaviour*.

    Parameters
    ----------
    n_nodes : int
        Table size.  Sweeping this moves the kink spacing across the probe grid: at
        ``n_nodes = 200`` the nodes coincide with ``n_probe0``, and at 5000 they are finer than
        ``max_n_probe = 6400``.
    kind : {'linear', 'cubic'}
        ``'linear'`` for the kinks, ``'cubic'`` for the ringing.
    """
    nodes = np.linspace(float(l0), float(l1), int(n_nodes))
    base = H.solar_ne() if base is None else base
    values = np.asarray([float(base(x)) for x in nodes])

    if kind == 'linear':
        def ne(l):
            return H.scalarize(np.interp(np.asarray(l, dtype=float), nodes, values))
        return ne
    if kind == 'cubic':
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(nodes, values, extrapolate=True)

        def ne(l):
            # A spline can undershoot to negative values on a steep table; the wrappers reject
            # a negative density, and that is correct validation, not a defect.  Clip at a tiny
            # positive floor so the family is testable at every node count.
            return H.scalarize(np.maximum(cs(np.asarray(l, dtype=float)), 1.0e-30*H.NE0))
        return ne
    raise ValueError(kind)


# ======================================================================
# 2.5  A real tabulated solar model
# ======================================================================

BS05_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bs05_agsop.dat')

# Trajectory for the tabulated-solar family: a radial ray from the centre.  Held to one solar
# scale height so the cost and the numbers are directly comparable with the 42-workload
# population, which uses exactly this span.
BS05_L0, BS05_L1 = 0.0, 1.0*H.L_SCALE


def load_bs05():
    r"""Radius [R_sun] and electron number density [eV^3] from the BS2005-AGS,OP table.

    ``n_e = \rho\,N_A\,(1+X)/2`` for fully ionized H + He, with the hydrogen mass fraction ``X``
    read from the table rather than assumed: :math:`Y_e` runs from 0.68 at the centre to 0.86 at
    the surface, so a fixed ``electron_fraction=0.5`` would be wrong by up to 70 %.  That is why
    this returns a number density directly instead of a mass density for
    ``vcc_func_from_rho_func`` to convert -- that function takes a *scalar* electron fraction.

    The file is the published table [Bahcall2005]_; the header and the two trailing
    ``Lsun=``/``Rsun=`` lines are skipped by requiring twelve float fields.
    """
    rows = []
    with open(BS05_PATH) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) != 12:
                continue
            try:
                rows.append([float(x) for x in fields])
            except ValueError:            # the column-heading line
                continue
    table = np.array(rows)
    r_over_rsun = table[:, 1]
    rho, x_hydrogen = table[:, 3], table[:, 6]
    return r_over_rsun, ne_from_rho_cgs(rho, y_e=0.5*(1.0 + x_hydrogen))


def bs05_ne(kind='linear'):
    r"""Electron number density along a radial ray, from the real solar model.

    Interpolated in :math:`\log n_e` -- which is what anyone reading a stellar model does, and
    what keeps the result positive -- so the kinks are in the log, exactly as a user's would be.

    **The node spacing of a real model is the interesting number.**  The table carries 563 nodes
    inside one solar scale height, a median spacing of 4.1e-3 of that span, against the
    ``n_probe0 = 200`` grid's 5.0e-3.  A published solar model is tabulated at almost exactly
    the density of the package's own first probe grid.

    Provenance: highest in this population.  A real published model, used directly.
    """
    r, ne = load_bs05()
    x = r*gd.SUN_RADIUS*KM                    # radius in natural units
    log_ne = np.log(ne)
    # The table starts at r = 0.00161 R_sun and ends at 0.983, so both ends need a rule.  Clamp
    # rather than extrapolate: the density is nearly flat in the core, and a cubic extrapolated
    # below the first node undershoots it by 30 %, which would be a finding about the
    # extrapolation rather than about the model.
    x_lo, x_hi = float(x[0]), float(x[-1])
    if kind == 'linear':
        def f(l):
            return H.scalarize(np.exp(np.interp(np.asarray(l, dtype=float), x, log_ne)))
        return f
    if kind == 'cubic':
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(x, log_ne, extrapolate=True)

        def f(l):
            xs = np.clip(np.asarray(l, dtype=float), x_lo, x_hi)
            return H.scalarize(np.exp(cs(xs)))
        return f
    raise ValueError(kind)


# ======================================================================
# 2.2  Supernova shock front
# ======================================================================

# [Kneller2015]_ Fig. 1: the three discontinuities of the M = 10.8 Msun model of Fischer et al.
# at t = 3 s post-bounce.
R_REVERSE_SHOCK_KM = 1734.0
R_CONTACT_KM = 12348.0
R_FORWARD_SHOCK_KM = 30323.0

# The trajectory.  Truncated at 1e4 km rather than started at the reverse shock because the
# oracle cost is set by the accumulated matter phase, which diverges as r^-1.4 inwards: the
# full [1734, 8e4] km ray costs ~10x what [1e4, 8e4] does, for one extra feature.  What is kept
# is the contact discontinuity and the forward shock -- **two narrow features on one
# trajectory**, which is what the window-merge logic and the hidden-feature scan need exercised.
# The span is 7e4 km = 1.06 solar scale heights, so costs are comparable to the solar batteries.
SN_R0_KM, SN_R1_KM = 1.0e4, 8.0e4
SN_L0, SN_L1 = SN_R0_KM*KM, SN_R1_KM*KM

# The H resonance for the package's default parameters sits at 4.1e4 km at 15 MeV, 5.5e4 km at
# 30 MeV and 6.4e4 km at 45 MeV -- inside this ray at every energy in the SN band, and just
# outside the forward shock, which is the configuration the shock-effect literature studies.
SN_ENERGIES = (15.0e6, 30.0e6, 45.0e6)

FOGLI_RHO0_G_PER_CM3 = 1.0e14      # rho_0(x) = 1e14 (x/km)^-2.4  [Fogli2003]_
FOGLI_EXPONENT = -2.4
FOGLI_SHOCK_JUMP = 10.0            # xi = V+/V- ~ 10  [Fogli2003]_
CONTACT_JUMP = 2.5
r"""float: Density ratio across the contact discontinuity.

[Fogli2003]_ parametrizes only the forward shock; Fig. 1 of [Kneller2015]_ shows a contact
discontinuity of order a factor of a few.  **2.5 is chosen, not fitted** -- what this family
measures is the *width* sweep, and the answer must not turn on this number.  ``sweep_contact``
in ``physical_battery.py`` varies it to confirm it does not.
"""


def _fogli_rho0(r_km):
    r"""Static progenitor density [g cm^-3], :math:`10^{14}(x/\mathrm{km})^{-2.4}` [Fogli2003]_."""
    return FOGLI_RHO0_G_PER_CM3*np.asarray(r_km, dtype=float)**FOGLI_EXPONENT


def _fogli_rarefaction(r_km, r_shock_km):
    r""":math:`f(x)`, the shape of the rarefaction zone behind the forward shock.

    :math:`\ln f(x) = [0.28 - 0.69\ln(x_s/\mathrm{km})]\,[\arcsin(1-x/x_s)]^{1.1}` [Fogli2003]_,
    which is 1 at the shock and falls to ~1e-5 deep behind it (the hot bubble).  Defined for
    :math:`x \le x_s`; the caller only evaluates it there.
    """
    x = np.asarray(r_km, dtype=float)
    u = np.clip(1.0 - x/float(r_shock_km), 0.0, 1.0)
    coeff = 0.28 - 0.69*np.log(float(r_shock_km))
    return np.exp(coeff*np.arcsin(u)**1.1)


def sn_shock_ne(width_frac, r_forward_km=R_FORWARD_SHOCK_KM, r_contact_km=R_CONTACT_KM,
                contact_jump=CONTACT_JUMP, y_e=0.5):
    r"""Post-bounce supernova envelope with a forward shock and a contact discontinuity.

    The textbook physical case where a real density profile is near-discontinuous, and a
    well-studied MSW problem.  Static progenitor power law, jump factor and rarefaction shape are
    all from [Fogli2003]_; the shock radii are the ones [Kneller2015]_ reads off a
    :math:`10.8\,M_\odot` simulation at :math:`t = 3` s.  The only invented number is
    :data:`CONTACT_JUMP`.

    The two fronts are given a finite width so the sweep has an axis: each is a :math:`C^1`
    ramp of width ``width_frac`` times the trajectory.  At ``width_frac = 1e-6`` (0.07 km) the
    profile is a discontinuity for every purpose; at 1e-2 (700 km) it is resolvable.

    Provenance: an analytic parametrization fitted to simulations in the literature, with radii
    from a named simulation.  Not a simulation snapshot.
    """
    w_km = float(width_frac)*(SN_R1_KM - SN_R0_KM)

    def ne(l):
        r_km = np.asarray(l, dtype=float)/KM
        rho = _fogli_rho0(r_km)
        # Forward shock: unshocked outside, shocked (xi * rarefaction) inside.
        shocked = _smoothstep((r_forward_km + 0.5*w_km - r_km)/w_km)
        factor = 1.0 + shocked*(FOGLI_SHOCK_JUMP*_fogli_rarefaction(r_km, r_forward_km) - 1.0)
        # Contact discontinuity: a further step deeper in.
        inside = _smoothstep((r_contact_km + 0.5*w_km - r_km)/w_km)
        factor = factor*(1.0 + inside*(contact_jump - 1.0))
        return H.scalarize(ne_from_rho_cgs(rho*factor, y_e=y_e))
    return ne


def sn_shock_breakpoints(width_frac, r_forward_km=R_FORWARD_SHOCK_KM,
                         r_contact_km=R_CONTACT_KM):
    """The two front locations as ``t_breakpoints``, to answer "does declaring them cure it?"."""
    w_km = float(width_frac)*(SN_R1_KM - SN_R0_KM)
    edges = []
    for r in (r_contact_km, r_forward_km):
        edges += [(r - 0.5*w_km)*KM, (r + 0.5*w_km)*KM]
    return np.array([SN_L0] + sorted(edges) + [SN_L1])


# ======================================================================
# 2.3  Density fluctuations with a Kolmogorov spectrum
# ======================================================================

TURBULENCE_DAMPING_KM = 100.0      # lambda_star  [Kneller2015]_
TURBULENCE_ALPHA = 5.0/3.0         # Kolmogorov
TURBULENCE_N_MODES = 45            # N_q; [Kneller2015]_ requires 40-50 for the dynamic range
TURBULENCE_DYNAMIC_RANGE = 1.0e5   # q_max/q_star, the 40-50 dB of [Kneller2015]_


def sn_turbulent_ne(rng, c_star=0.1, r_lo_km=R_CONTACT_KM, r_hi_km=R_FORWARD_SHOCK_KM,
                    alpha=TURBULENCE_ALPHA, n_modes=TURBULENCE_N_MODES,
                    dynamic_range=TURBULENCE_DYNAMIC_RANGE, base=None):
    r"""Kolmogorov density fluctuations between the contact discontinuity and the forward shock.

    :math:`n_e \to n_e\,(1 + F(r))` with, verbatim from [Kneller2015]_ eq. (2),

    .. math::
        F(r) = C_\star \tanh\!\Big(\frac{r-r_<}{\lambda_\star}\Big)
                       \tanh\!\Big(\frac{r_>-r}{\lambda_\star}\Big)
               \sum_{a=1}^{N_q}\sqrt{V_a}\,[A_a\cos(q_a r) + B_a\sin(q_a r)] ,

    the two ``tanh`` factors damping the fluctuations at the edges of the turbulence region over
    :math:`\lambda_\star = 100` km.  Wavenumbers are drawn by the randomization method's
    "variant C": wavenumber space is split into :math:`N_q` regions, one :math:`q_a` is drawn
    from each using the normalized spectrum
    :math:`E(q) \propto (q_\star/q)^{5/3}\,\Theta(|q|-q_\star)` as a pdf, and :math:`V_a` is the
    integral of :math:`E` over that region.  :math:`q_\star = 1/[2(r_> - r_<)]`, the driving
    scale, is twice the size of the turbulence domain.

    **Why this family is here and a random Fourier sum is not.**  A Fourier sum with a handful of
    modes has a shortest wavelength and nothing below it.  A power-law spectrum deliberately has
    power at **every** scale including below the probe grid: with the paper's 40-50 dB dynamic
    range the shortest mode is 0.36 km, or 2e-5 of the turbulence region, against 10.9 km for
    the 6400-point reference grid.  At :math:`C_\star = 10\ \%` the finest grid the package can
    reach therefore sees under a third of this profile's total variation.

    **It is nonetheless not detected**, by either instrument -- the resolution test declares it
    resolved and ``find_hidden_features`` returns a concentration two orders below its threshold.
    That is a property of the statistic rather than of its threshold: concentration asks whether
    hidden variation piles up in *one* reference interval, and a power-law spectrum spreads it
    evenly over all 6400 of them.  Recorded here so the docstring does not claim a coverage the
    measurement does not support.

    Provenance: the construction used in the published literature, with the paper's spectrum,
    damping scale and amplitude range (it uses :math:`C_\star` = 0.1 %, 1 % and 10 %).

    Parameters
    ----------
    c_star : float
        RMS amplitude of :math:`\delta n_e/n_e`.  The :math:`V_a` are normalized to sum to one
        here, so ``c_star`` *is* the rms rather than being proportional to it.
    base : Callable, optional
        Underlying profile.  Default: the shock profile at a resolvable width, so the turbulence
        sits on the physical envelope it belongs to rather than on a bare power law.
    """
    r_lo, r_hi = float(r_lo_km), float(r_hi_km)
    q_star = 1.0/(2.0*(r_hi - r_lo))                    # [km^-1]
    q_max = q_star*float(dynamic_range)

    # Split [q_star, q_max] into n_modes geometric regions; draw one q from each with E as pdf,
    # and integrate E over the region for V_a.  With E ~ q^-alpha,
    #   int_{q1}^{q2} E dq  =  (1/2) q_star^{alpha-1} (q1^{1-alpha} - q2^{1-alpha}).
    edges = np.geomspace(q_star, q_max, int(n_modes) + 1)
    q1, q2 = edges[:-1], edges[1:]
    v = 0.5*q_star**(alpha - 1.0)*(q1**(1.0 - alpha) - q2**(1.0 - alpha))
    v = v/v.sum()                                       # so c_star is exactly the rms
    # Inverse-transform sample one q per region from E restricted to that region.
    u = rng.random(int(n_modes))
    qa = (q1**(1.0 - alpha) - u*(q1**(1.0 - alpha) - q2**(1.0 - alpha)))**(1.0/(1.0 - alpha))
    amp_a = rng.normal(size=int(n_modes))
    amp_b = rng.normal(size=int(n_modes))
    sqrt_v = np.sqrt(v)

    base = sn_shock_ne(1.0e-3) if base is None else base

    def ne(l):
        r_km = np.asarray(l, dtype=float)/KM
        env = (np.tanh(np.clip((r_km - r_lo)/TURBULENCE_DAMPING_KM, 0.0, None))
               * np.tanh(np.clip((r_hi - r_km)/TURBULENCE_DAMPING_KM, 0.0, None)))
        ph = qa*r_km[..., None]
        f = c_star*env*(sqrt_v*(amp_a*np.cos(ph) + amp_b*np.sin(ph))).sum(axis=-1)
        # 1 + F can go negative for large c_star; the wrappers reject a negative density and
        # that is correct validation.  Floor it so the family is testable at every amplitude,
        # and note that the floor is never reached at the paper's amplitudes.
        return H.scalarize(np.asarray(base(l))*np.maximum(1.0 + f, 1.0e-6))
    return ne


# ======================================================================
# 2.4  Earth with a non-PREM crust
# ======================================================================

# PREM's own crust is 2.6 g/cm^3 (upper) and 2.9 (lower) over the outer 24.4 km, plus a 3 km
# ocean at 1.02.  A user with their own crust model, a 3-D tomographic slice, or an added
# sediment layer supplies something else -- and, crucially, would not pass t_breakpoints,
# because they are not using the package's osc_prob_earth path that fills them in from
# PREM_BOUNDARIES.  These are plausible continental-crust values, NOT a measured model.
EXTRA_CRUST_LAYERS = (
    (0.0, 4.0, 2.2),        # sediment: depth from surface [km], to [km], density [g/cm^3]
    (4.0, 18.0, 2.75),      # upper crust
    (18.0, 38.0, 3.05),     # lower crust, thicker than PREM's
)


def earth_crust_ne(costhz, extra_layers=EXTRA_CRUST_LAYERS, y_e=0.5):
    r"""PREM along a chord, with the outer shells replaced by a user's own crust model.

    ``osc_prob_earth`` passes PREM's layer edges as ``t_breakpoints``, which is why Earth has been
    clean across 164 configurations.  This is the same trajectory **without** that help, and with
    crustal layers PREM does not have -- so the edges fall where the package has no reason to
    expect them.  A chord crosses each layer **twice**, once on the way in and once on the way
    out, so three extra layers put six extra discontinuities on the path.

    Provenance: PREM is a real model; the extra layers are plausible continental-crust values,
    invented.  What is physical is the user behaviour, not the crust.

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle; must be negative for an upgoing trajectory.
    """
    r_earth = gd.EARTH_RADIUS

    def ne(l):
        # TRAP GUARD: density_matter_func_prem takes a radius from Earth's CENTRE in km, and its
        # second argument is `tol`, not costhz.  The route is always through
        # earth_radial_distance_from_depth(costhz, l_in_km).
        l_km = np.asarray(l, dtype=float)/KM
        r_km = np.asarray(earth.earth_radial_distance_from_depth(costhz, l_km))
        r_km = np.clip(r_km, 0.0, r_earth)
        rho = np.asarray(earth.density_matter_func_prem(r_km), dtype=float)
        depth = r_earth - r_km
        for (d_lo, d_hi, value) in extra_layers:
            rho = np.where((depth >= d_lo) & (depth < d_hi), value, rho)
        return H.scalarize(ne_from_rho_cgs(rho, y_e=y_e))
    return ne


def earth_chord_span(costhz):
    """(l0, l1) in natural units for a chord at this zenith angle."""
    return 0.0, earth.distance_traveled_inside_earth(costhz)*KM


# ======================================================================
# The population
# ======================================================================

def _fam(label, ne, l0, l1, kind, hides, energies, provenance, t_breakpoints=None):
    return dict(label=label, ne=ne, l0=float(l0), l1=float(l1), kind=kind,
                hides=bool(hides), energies=tuple(energies), provenance=provenance,
                t_breakpoints=t_breakpoints)


SOLAR_ENERGIES = (30.0e6, 100.0e6)
EARTH_ENERGIES = (30.0e6, 100.0e6)


def families(which='all', seed=20260804):
    """The physical population.

    Returns a list of dicts with keys ``label``, ``ne``, ``l0``, ``l1``, ``kind``, ``hides``,
    ``energies``, ``provenance`` and ``t_breakpoints``.

    ``kind`` is one of ``'tabulated'``, ``'bs05'``, ``'sn_shock'``, ``'sn_turbulence'``,
    ``'earth_crust'``.  ``hides`` says whether the family genuinely carries structure below the
    finest grid the package lays down -- the split the P2 (no new false positives) and P3
    (detection rate) criteria need, and it is set from the construction, not from what the
    detector says about it.
    """
    rng = np.random.default_rng(seed)
    out = []

    # ---- 2.1 interpolation kinks.  N crosses n_probe0 = 200 and max_n_probe = 6400.
    for n in (20, 50, 200, 1000, 5000):
        out.append(_fam(
            'tabulated linear N=%d' % n, tabulated_ne(n, 0.0, 1.0*H.L_SCALE, 'linear'),
            0.0, 1.0*H.L_SCALE, 'tabulated', n > 6400, SOLAR_ENERGIES,
            'package exponential resampled; models a user loading a table'))
    for n in (50, 1000):
        out.append(_fam(
            'tabulated cubic N=%d' % n, tabulated_ne(n, 0.0, 1.0*H.L_SCALE, 'cubic'),
            0.0, 1.0*H.L_SCALE, 'tabulated', False, SOLAR_ENERGIES,
            'as above, C2 spline instead of C0 linear'))

    # ---- 2.5 the real solar model.
    for k in ('linear', 'cubic'):
        out.append(_fam(
            'BS05(AGS,OP) %s' % k, bs05_ne(k), BS05_L0, BS05_L1, 'bs05', False,
            SOLAR_ENERGIES, 'published standard solar model, used directly'))

    # ---- 2.2 supernova shock, width sweep.
    #
    # hides = False at EVERY width, which is not what the brief expected and is worth stating
    # plainly.  A shock front is a **monotone step**: however narrow it is, the two reference
    # nodes that bracket it see its full height, so the total variation measured on the
    # 6400-point grid equals the total variation measured on a 32x finer one -- ratio 1.000 at
    # w = 1e-2 through 1e-6 alike (``physical_battery.py sub_grid`` measures this rather than
    # taking it on trust).  Nothing is hidden from the grid; what a narrow shock defeats is
    # *resolution*, which is ``adiabatic._profile_is_resolved``'s job, not
    # ``find_hidden_features``'.  The scan is a concentration test for variation the coarse grid
    # cannot see, and a step is not that.
    for w in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6):
        out.append(_fam(
            'SN shock w=%.0e' % w, sn_shock_ne(w), SN_L0, SN_L1, 'sn_shock',
            False, SN_ENERGIES,
            'Fogli et al. 2003 profile; shock radii from Kneller et al. 2015 Fig. 1'))

    # ---- 2.3 Kolmogorov turbulence, at the paper's amplitudes.
    for c in (0.01, 0.1):
        out.append(_fam(
            'SN turbulence C*=%g' % c, sn_turbulent_ne(rng, c_star=c), SN_L0, SN_L1,
            'sn_turbulence', True, SN_ENERGIES,
            'Kneller et al. 2015 eq. (2), alpha = 5/3, lambda* = 100 km'))

    # ---- 2.4 Earth with a non-PREM crust, no t_breakpoints.
    for costhz in (-1.0, -0.8, -0.4):
        l0, l1 = earth_chord_span(costhz)
        out.append(_fam(
            'Earth crust costhz=%.1f' % costhz, earth_crust_ne(costhz), l0, l1,
            'earth_crust', False, EARTH_ENERGIES,
            'PREM plus three invented crustal layers, breakpoints NOT declared'))

    if which != 'all':
        out = [f for f in out if f['kind'] == which]
    return out


SMOOTH_KINDS = ('tabulated', 'bs05', 'earth_crust')
r"""tuple: Families with no sub-grid structure, so a ``find_hidden_features`` firing is a
false positive.  Note that ``tabulated`` and ``earth_crust`` are *not* smooth in the
:math:`C^1` sense -- they have kinks and steps -- but every one of those features is wider than
the reference grid, which is the only thing the hidden-feature test claims to be about."""

FEATURED_KINDS = ('sn_shock', 'sn_turbulence')
r"""tuple: Families whose narrowest feature is below the probe grid.

Note that this is **not** the same split as ``hides``.  ``FEATURED_KINDS`` is about the
*resolution* test -- can a 200- or 6400-point probe grid see the feature at all -- and both
supernova families defeat it.  ``hides`` is about ``find_hidden_features``, which asks the
narrower question of whether variation exists that the reference grid cannot see even in
principle; a shock front is a monotone step and fails that test, so only the turbulence has
``hides = True``.  Keeping the two apart is the point: the earlier tranches conflated
"unresolved" with "invisible", and they are different failures with different cures."""


PROVENANCE = {
    'tabulated': 'package shape, physical user behaviour',
    'bs05': 'real published model, used directly',
    'sn_shock': 'literature parametrization, literature radii',
    'sn_turbulence': 'literature construction and spectrum',
    'earth_crust': 'real PREM plus invented crustal layers',
}


if __name__ == '__main__':
    fams = families()
    print('# %d physical families\n' % len(fams))
    print('%-26s %-14s %6s %10s %10s  %s'
          % ('family', 'kind', 'hides', 'span/km', 'ne(l0)', 'provenance'))
    for f in fams:
        span_km = (f['l1'] - f['l0'])/KM
        print('%-26s %-14s %6s %10.3e %10.3e  %s'
              % (f['label'], f['kind'], f['hides'], span_km, float(f['ne'](f['l0'])),
                 PROVENANCE[f['kind']]))
