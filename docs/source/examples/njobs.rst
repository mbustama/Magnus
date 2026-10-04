.. _ex-sec-njobs:

Running a scan in parallel
--------------------------

A scan of many probabilities can be shared among processes. Every wrapper accepts the argument ``n_jobs`` for that, the number of processes to use, with ``n_jobs = 1`` as default and ``n_jobs = -1`` for all cores:

.. code-block:: python

   import numpy as np
   import magnus.oscprob as oscprob
   import magnus.globaldefs as gd

   N = 5000
   E = np.logspace(-0.3, 1.3, N)*gd.UNIT_GEV
   L = np.linspace(2e3, 11467.8, N)*gd.UNIT_KM
   kw = dict(costhz=-0.9, rtol=1e-6, atol=1e-8,
             nu_i=gd.NUMU, nu_f=gd.NUMU)

   # Every point has its own baseline, so no
   # batched engine takes this scan
   P = oscprob.osc_prob_3nu_earth(
       E, L=L, n_jobs=4, **kw)


The first point of the scan runs by itself, in the main process, and fixes two numbers: how many slabs the trajectory needs, :math:`N_{\rm slabs}`, and how many sample points each slab needs. Those two numbers, set two rungs below where the first point landed, become the starting point of the refinement ladder for every later point, so each worker begins close to the answer. The rest of the scan goes to the workers, never more of them than there are points left or cores. Each worker takes a chunk of consecutive points and computes them one after another; ``joblib`` sizes the chunks to run for between :math:`0.2` and :math:`2` s, so cheap points are gathered in bulk and an expensive point is handed over on its own. If the first point shows that the rest of the scan would take less than about 1 s in a single process, Magνs finishes it there and starts no workers, since starting them would cost more.

Two conditions decide whether using ``n_jobs`` larger than 1 helps. The first is the shape of the scan. A batched engine answers only at ``n_jobs`` :math:`= 1`, so asking for more processes sends the scan down the per-point path instead, which for a batchable scan is the slower of the two (:doc:`/engines`). A scan whose points share one baseline should keep the default of ``n_jobs`` :math:`= 1` and let the batching of :ref:`ex-sec-batched-calls` do the work:

.. code-block:: python

   # One baseline for all energies: the batched
   # engine takes it, at the default n_jobs = 1
   P = oscprob.osc_prob_3nu_earth(
       E, L=11467.8*gd.UNIT_KM, **kw)


Parallelism is for scans that no batched engine accepts, such as the first one above, where every energy has a different baseline.

The second condition is length. Starting the workers costs time, paid once per call, so the scan has to be long enough to recover it. On the scan above, four workers more than double the speed at :math:`2\,000` points or more. At a few hundred points, they bring no gain, and can even slow the scan down. :doc:`/performance` shows how the gain grows with the size of the scan.

Parallelism changes more than the speed. A serial scan starts each point from where its neighbor converged; a parallel scan starts every point from where the first converged. A ladder that starts on a different rung can stop on a different rung, so two runs that differ in ``n_jobs`` agree to the requested tolerance, but not always to the last digit. On the chord above at ``rtol`` :math:`= 10^{-6}`, a scan of 40 points, short enough to stay in one process, agrees exactly, while 500 points differ by about :math:`10^{-7}`. The user should hold ``n_jobs`` fixed when two runs must match to the last digit.

``n_jobs`` affects only a scan of several probabilities. A call for a single probability accepts it and ignores it; its slabs go to the batched kernel in a single call, and no worker is started:

.. code-block:: python

   # One point: n_jobs is accepted and does
   # nothing; the call runs in one process
   P = oscprob.osc_prob_3nu_earth(
       10.0*gd.UNIT_GEV, L=11467.8*gd.UNIT_KM,
       n_jobs=4, **kw)


