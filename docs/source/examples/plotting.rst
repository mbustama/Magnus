.. _ex-sec-plotting:

Pre-packaged plotting routines
------------------------------

The ``plotting`` module of Magνs draws several of the standard figures of the field, each with a single call; Table :ref:`ex-tab-plotting` lists them, and the sections that follow use them. Drawn directly in ``matplotlib``, the same figure takes 25–40 lines of grid specification, tick locators, legend keywords, and ``savefig`` calls.

.. _ex-tab-plotting:

.. table:: The plotting functions

   +------------------------------+---------------------------------------------------+
   | Function (``plot_``)         | Draws                                             |
   +==============================+===================================================+
   | ``probability_vs_energy``    | Probability against energy                        |
   +------------------------------+---------------------------------------------------+
   | ``probability_vs_baseline``  | Probability against baseline                      |
   +------------------------------+---------------------------------------------------+
   | ``curves``                   | Any swept variable, optional error subpanel       |
   +------------------------------+---------------------------------------------------+
   | ``curves_stacked``           | One panel per case, over a shared abscissa        |
   +------------------------------+---------------------------------------------------+
   | ``oscillogram``              | Probability over energy and arrival direction     |
   +------------------------------+---------------------------------------------------+
   | ``biprobability``            | Antineutrino against neutrino probability         |
   +------------------------------+---------------------------------------------------+
   | ``probability_with_average`` | Oscillating probability, with its average         |
   +------------------------------+---------------------------------------------------+
   | ``probability_with_profile`` | Density profile, above the probabilities along it |
   +------------------------------+---------------------------------------------------+

The plotting functions of ``magnus.plotting``, each named with the prefix ``plot_``. Each returns the ``matplotlib`` figure it built and its axes. See :ref:`ex-sec-plotting` for details.

Most of these functions take the abscissa first, then a list of curves, each a dictionary with its data, ``y``, and its ``label``. For instance, the probability of :math:`\nu_\mu \to \nu_e` against energy, in vacuum and at constant density, with ``P_scan`` the 200-energy scan of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, is drawn by

.. code-block:: python

   import magnus.plotting as plotting

   P_vac = oscprob.osc_prob_3nu_vacuum(Es, L)
   fig, ax = plotting.\
       plot_probability_vs_energy(
       Es/gd.UNIT_GEV,
       [dict(y=P_vac[:, gd.NUMU, gd.NUE],
             label='Vacuum'),
        dict(y=P_scan[:, gd.NUMU, gd.NUE],
             label='Matter')],
       nu_i=gd.NUMU, nu_f=gd.NUE)


The flavor indices ``nu_i`` and ``nu_f`` only label the vertical axis, here :math:`P_{\nu_\mu \to \nu_e}`; the curves are the data passed. The API documentation :cite:p:`MagnusDocs` gives the arguments of every ``plotting`` function.

Each figure comes out in a preset style, so a set of figures is consistent without any styling passed. Since each function returns the ``matplotlib`` figure and its axes, the style can be changed afterwards. For instance, the figure above is resized, given a grid, and saved by

.. code-block:: python

   fig.set_size_inches(6, 3)
   ax.grid(True)
   fig.savefig('numu_to_nue.pdf')
