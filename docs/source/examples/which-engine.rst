.. _ex-sec-which-engine:

Strategy report after a call
----------------------------

Every wrapper and scenario function takes a dictionary as ``strategy_info`` and fills it in place with a report on how the probability was computed. Three functions do not take it, because they evaluate a closed-form expression and use no engine: ``osc_prob_2nu_vacuum_std``, ``osc_prob_3nu_vacuum_std``, and ``osc_prob_2nu_matter_std``. Nor does the direct call ``osc_prob``, which reports on its refinement ladder through ``convergence_info`` instead (:doc:`/methodology`).

The report has seven entries: ``'engine'``, the engine that answered, named as in :doc:`/engines`; ``'family'``, the group of engines that share its method, and so can err in the same way (:doc:`/engines`); ``'certified'``, whether the hybrid engine certified its answer, or ``None`` if another engine answered; ``'declined'``, the engines that attempted the request and gave up, each with its reason; ``'trace'``, every engine tried, in order, with its details; ``'hidden_feature'``, the result of a scan of the profile for a feature narrower than the grids on which the engines sample it, or ``None`` if the call declares where the profile changes, as on the Earth, or if the density is constant; and ``'sampling'``, described below. For instance,

.. code-block:: python

   L = earth.distance_traveled_inside_earth(
       -0.9)*gd.UNIT_KM
   info = {}
   P = oscprob.osc_prob_3nu_earth(
       10.0*gd.UNIT_GEV, costhz=-0.9, L=L,
       nu_i=gd.NUMU, nu_f=gd.NUMU,
       strategy_info=info)
   sorted(info)
   # ['certified', 'declined', 'engine',
   #  'family', 'hidden_feature', 'sampling',
   #  'trace']
   info['engine']                 # 'magnus'
   info['family']          # 'magnus-ladder'


This single energy went to the general Magnus ladder.

The ``'sampling'`` entry shows how many points a scan along the trajectory needs to resolve the oscillation. It gives the shortest oscillation length along the trajectory, how many oscillations of that length fit in the trajectory, and how many points a scan needs to sample them twice per cycle, as the Nyquist criterion requires. A scan with fewer points is aliased: each value is correct, but a curve drawn through them is not the probability. For this call,

.. code-block:: python

   s = info['sampling']
   s['oscillation_length']/gd.UNIT_KM
                                  # 3253.4
   s['cycles_over_trajectory']    # 3.5
   s['nyquist_points']            # 9


On the same chord, an array of energies, a phase average, and a constant density each go to a different engine:

.. code-block:: python

   kw = dict(costhz=-0.9, L=L, nu_i=gd.NUMU,
             nu_f=gd.NUMU, strategy_info=info)
   E = np.linspace(5.0, 15.0, 3)*gd.UNIT_GEV
   oscprob.osc_prob_3nu_earth(E, **kw)
   info['engine']              # 'separable'

   oscprob.osc_prob_3nu_earth(
       10.0*gd.UNIT_GEV, **kw, average=True)
   info['engine']                # 'average'

   oscprob.osc_prob_3nu_matter_constant_density(
       10.0*gd.UNIT_GEV, L, 4.0,
       density_matter_is_in_g_per_cm3=True,
       nu_i=gd.NUMU, nu_f=gd.NUMU,
       strategy_info=info)
   info['engine']               # 'constant'


When an engine declines a request, ``'declined'`` gives the reason. In this example, along an Earth diameter, the density falls as :math:`3\,e^{-l/(50~{\rm m})}` g cm\ :math:`^{-3}`, with :math:`l` the distance traveled. That is narrower than the spacing of the grid on which the hybrid engine probes the profile, so the hybrid engine declines, with a warning, and the interaction-picture engine answers:

.. code-block:: python

   info = {}
   P = oscprob.osc_prob_2nu_matter_exp_density(
       5.0*gd.UNIT_MEV, 12742.0*gd.UNIT_KM,
       0.0, 3.0, 0.05*gd.UNIT_KM, sth=0.5,
       Dm2=2.5e-3,
       density_matter_is_in_g_per_cm3=True,
       nu_i=0, nu_f=1, strategy_info=info)
   info['declined']
   # [('hybrid', 'the profile is not
   #   resolved at the probe scale')]
   info['engine']                # 'ip_exp'


With a fall-off length of 50 km instead of 50 m, the hybrid engine resolves the profile and certifies its own answer:

.. code-block:: python

   # 50.0*gd.UNIT_KM, not 0.05*gd.UNIT_KM
   info['engine']                # 'hybrid'
   info['certified']             # True
   info['declined']              # []
