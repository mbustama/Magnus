Standard solar models
=====================

The Sun entry points -- ``osc_prob_{2,3,4,5}nu_sun[_nsi|_liv]`` and
:func:`~magnus.oscprob.osc_prob_sun` -- describe the Sun by an exponential fit to its
electron density unless told otherwise.  ``density_profile`` puts a published standard
solar model in its place, from tables that ship with the package.

.. versionadded:: 1.1.1

.. list-table::
   :header-rows: 1
   :widths: 20 38 20 22

   * - ``density_profile``
     - Model
     - Tabulated range [:math:`R_\odot`]
     - Reference
   * - ``'BP2000'``
     - Bahcall, Pinsonneault & Basu
     - 0.0065 -- 0.949
     - :cite:`Bahcall:2000nu`
   * - ``'BP04'``
     - Bahcall & Pinsonneault
     - 0.0065 -- 0.947
     - :cite:`Bahcall:2004fg`
   * - ``'BS05-OP'``
     - Bahcall, Serenelli & Basu, GS98 composition
     - 0.0016 -- 0.983
     - :cite:`Bahcall:2004pz`
   * - ``'BS05-AGS-OP'``
     - Bahcall, Serenelli & Basu, AGS05 composition
     - 0.0016 -- 0.983
     - :cite:`Bahcall:2004pz`
   * - ``'B16-GS98'``
     - Vinyoles et al., GS98 composition
     - 0.0005 -- 1
     - :cite:`Vinyoles:2016djt`
   * - ``'B16-AGSS09met'``
     - Vinyoles et al., AGSS09met composition
     - 0.0005 -- 1
     - :cite:`Vinyoles:2016djt`
   * - ``'B23-GS98'``, ``'B23-AGSS09'``, ``'B23-C11'``, ``'B23-AAG21'``, ``'B23-MB22m'``,
       ``'B23-MB22p'``
     - Herrera & Serenelli, one per solar composition
     - 0 -- 1
     - :cite:`Herrera:2023b23`

Names are matched without regard to case, and ``'exp'`` is the default.

.. jupyter-execute::

    import magnus.globaldefs as gd
    from magnus import oscprob

    R_sun = gd.SUN_RADIUS*gd.UNIT_KM

    for profile in ('exp', 'BP04', 'B16-GS98', 'B23-AAG21'):
        P_ee = oscprob.osc_prob_3nu_sun(energy=10*gd.UNIT_MEV, L=R_sun, L0=0.0,
                                        nu_i=gd.NUE, nu_f=gd.NUE, average=True,
                                        density_profile=profile)
        print('%-10s <P_ee> = %.4f' % (profile, P_ee))
    # exp 0.2971, BP04 0.3185, B16-GS98 0.3189, B23-AAG21 0.3192

``L0`` is the radius where the neutrino is produced, 0 or more, and a scalar: ``L``
may be an array of end points, but ``L0`` may not, so for several production points call
once per point.

From the command line, ``magnus --environment sun`` takes the same names through
``--density-profile`` (see :doc:`cli`).  `Notebook 13
<https://github.com/mbustama/Magnus/blob/main/notebooks/13_magnus_tabulated_solar_model.ipynb>`_
compares all twelve on the averaged observable: they agree to :math:`2.4\times10^{-3}`,
and the exponential fit is off by up to 0.1, at 2.5 MeV (0.02 at 10 MeV, as above).

.. figure:: ../../img/paper/solar_models.svg
   :width: 90%
   :alt: Two tabulated solar models and the exponential fit, and their averaged survival probabilities

   Top: the tabulated BS2005-AGS,OP and B16-GS98 models and the exponential fit the
   Sun functions use by default, with the ratio of each to BS2005-AGS,OP.  Bottom:
   :math:`\langle P_{\nu_e \to \nu_e}\rangle` at three flavors on each profile, and
   its difference from the BS2005-AGS,OP result.  The two tables agree to
   :math:`1.3 \times 10^{-3}`; the fit is low by up to 0.1.  From the Magνs paper.


What is tabulated, and what is computed from it
-----------------------------------------------

Each table holds three columns of the authors' file, copied as written: the radius, the
mass density :math:`\rho` and the hydrogen mass fraction :math:`X`.  The electron density
is

.. math::

   n_e = \rho N_A \left(\frac{X}{m_{\rm H}} + \frac{2(1 - X)}{m_{\rm He}}\right),

with :math:`m_{\rm H} = 1.00783` and :math:`m_{\rm He} = 4.00260` the atomic masses of
:math:`^1`\ H and :math:`^4`\ He in atomic mass units: everything heavier than hydrogen is
counted as helium-4.  Against the full composition of BS05(OP) (:math:`^1`\ H,
:math:`^4`\ He, :math:`^3`\ He, :math:`^{12}`\ C, :math:`^{14}`\ N, :math:`^{16}`\ O) this is
within 0.002% at the centre and 0.07% at worst, near :math:`0.3\,R_\odot` where
:math:`^3`\ He peaks.  The textbook :math:`\rho N_A (1 + X)/2`, which Bahcall's own
electron-density tables follow, is 0.4 to 0.8% higher, and the mean free-nucleon mass
used until 1.2.0, :math:`\rho (1 + X)/(2 m_N)`, put :math:`n_e` 0.1 to 0.4% lower (issue
#168).  Notebooks 13 and 28 use the same formula; inside the table, ``'BS05-AGS-OP'`` gives
notebook 28's profile bit for bit.  :math:`n_e` is interpolated linearly in its logarithm.

For the wrappers with sterile states, the neutral-current term needs the
neutron-to-proton ratio, helium's neutrons over all the protons, from the same atomic
masses: :math:`n_n/n_p = [2(1 - X)/m_{\rm He}]/[X/m_{\rm H} + 2(1 - X)/m_{\rm He}]`.  With a solar model and
``ratio_number_neutrons_to_protons`` left at its default, ``None``, the wrappers take it
from the same table, radius by radius; a value or a callable passed explicitly is used
instead.  The exponential fit carries no composition, and keeps the 1.0 it always had.
:func:`~magnus.oscprob.osc_prob_sun` passes only the electron density to ``H_func``; a
Hamiltonian with sterile states that wants the composition can take it from
:func:`magnus.solarmodels.neutron_to_proton_ratio_profile`.


Outside the table
-----------------

Only the B23 tables start at the centre, and the older ones stop short of the surface.
Below the first row, the density holds its first value, since the core is flat.  Past the
last row, it continues along the logarithmic slope of the last tabulated interval, so that
it keeps falling rather than stopping at a constant or dropping to zero.  That is a
continuation, not a model of the outer Sun: continued from 0.95 :math:`R_\odot`, BP2000
and BP04 reach about :math:`1.6\times10^{-3}` g cm\ :sup:`-3` at the surface, where B16 and
B23, which are tabulated there, give about :math:`1.7\times10^{-7}`.

``stop_at_table_edge=True`` refuses the continuation.  A baseline that ends past the last
row returns NaN, with a :class:`~magnus.oscprob.SolarModelRangeWarning` naming the edge,
and the others are computed as usual:

.. jupyter-execute::

    import warnings
    import numpy as np

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P_ee = oscprob.osc_prob_3nu_sun(energy=10*gd.UNIT_MEV,
                                        L=np.array([0.5, 0.9, 1.0])*R_sun, L0=0.0,
                                        nu_i=gd.NUE, nu_f=gd.NUE, average=True,
                                        density_profile='BP04', stop_at_table_edge=True)
    print(P_ee)
    print(caught[0].category.__name__)

With ``stop_at_table_edge=True``, a path that starts past the last row has nothing
inside the table to compute, and is refused with a ``ValueError``; without it, the
profile is continued and a number is returned.  ``stop_at_table_edge`` with the
exponential fit, which has no last row, is refused as well.  :func:`magnus.solarmodels.table_edge` gives the edge in the units of
``L``.


Coherent probabilities through a table
--------------------------------------

The interpolated profile has a kink at every row, and the Bahcall tables, printed to four
significant figures, step through the flat core.  The phase-averaged probability
(``average=True``) does not notice: on all twelve models at 1, 5, 10 and 20 MeV, from the
centre to 0.9 :math:`R_\odot`, it raised no warning.

A coherent probability over most of the Sun can notice.  The hybrid engine may then fail
to certify it and hand it to the slab ladder, which runs out of slabs and says so with
:class:`~magnus.oscprob.ToleranceNotAchievedWarning`.  On the same grid, 4 of the 48
coherent calls did (BP04 and BS05-OP at 10 and 20 MeV).  The answers were still within the
default tolerance.  On eleven cases, those four among them, every default answer was within
:math:`10^{-3}` of a reference converged to :math:`10^{-9}`, and the ones that warned were
the closest, within :math:`3\times10^{-5}`.  The warning is the package saying it could not
verify that.

For a coherent probability that is verified, put ``t_breakpoints`` at the rows, so that
no slab straddles a kink, and start the slab ladder fine enough to resolve the
oscillation:

.. jupyter-execute::

    from magnus import solarmodels

    rows = solarmodels.load_solar_model('BS05-AGS-OP')['r_over_r_sun']*R_sun

    P_ee = oscprob.osc_prob_3nu_sun(energy=10*gd.UNIT_MEV, L=0.95*R_sun, L0=0.0,
                                    nu_i=gd.NUE, nu_f=gd.NUE,
                                    density_profile='BS05-AGS-OP', strategy='magnus',
                                    t_breakpoints=rows, n_slabs=200_000,
                                    max_n_slabs=10_000_000)
    print('P_ee = %.9f' % P_ee)    # 0.126766302

On the same eleven cases this came within :math:`4\times10^{-9}` of the reference, with no
:class:`~magnus.oscprob.ToleranceNotAchievedWarning`.  Two still raised
:class:`~magnus.magnus.MagnusConvergenceWarning`, which reports a slab width rather than an
error (see :doc:`diagnostics`).  The measurements are in
``docs/dev/adversarial_batteries/solar_model_engines.py`` and ``solar_model_coherent.py``, in a
source checkout (a PyPI install does not ship them).


Provenance
----------

:func:`magnus.solarmodels.solar_model_info` gives each model's reference, the source it was
taken from, the date, the SHA-256 of the original file, and its terms of use.  Each table,
in ``magnus/data/solar_models/``, carries the same in its header, above the authors'
original header.  The Bahcall tables come from the IAS archive, whose terms ask only for a
note on how they are used; the B16 tables from the Internet Archive's copy of the authors'
page, which is no longer online; and the B23 tables from their Zenodo release, under
CC-BY-4.0.  If you use a model, cite its paper.

``tools/build_solar_model_tables.py`` rebuilds the tables from the original files, and with
``--check`` verifies that the committed ones are exactly what those files give.  It refuses
any original whose hash has changed.
