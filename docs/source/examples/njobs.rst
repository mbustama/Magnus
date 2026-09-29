.. _ex-sec-njobs:

Running a scan in parallel
--------------------------

A scan of many probabilities can be shared among processes. Every wrapper accepts the argument ``n_jobs`` for that, the number of processes to use, with ``n_jobs = 1`` as default:

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


The first point in the scan runs by itself, in the main process. That run fixes the parameters. Its refinement ladder settles on two numbers: how many slabs the trajectory needs (:math:`N_{\rm slabs}`), and how many collocation points each slab needs to evaluate the Magnus expansion. Those two numbers become the starting floor for every later point, set two rungs below where the first point landed. Each worker then begins its ladder close to the answer. The rest of the scan goes to the workers. Each worker takes a chunk of consecutive points and computes them one after another. ``joblib`` picks the batch size, aiming for chunks that run between :math:`0.2` and :math:`2` s, so cheap points are gathered in bulk and an expensive point is handed over on its own.

Two conditions decide whether using ``n_jobs`` larger than 1 helps. The first is the shape of the scan. A batched engine answers only at ``n_jobs`` :math:`= 1`, so asking for more processes sends the request down the per-point path instead, which for a batchable scan is the slower of the two (:doc:`/engines`). A scan whose points share one baseline should keep the default of ``n_jobs`` :math:`= 1` and let the batching of :ref:`ex-sec-batched-calls` do the work:

.. code-block:: python

   # One baseline for all energies: the batched
   # engine takes it, at the default n_jobs = 1
   P = oscprob.osc_prob_3nu_earth(
       E, L=11467.8*gd.UNIT_KM, **kw)


Parallelism is for scans that no batched engine accepts, such as the one above, where every energy carries a different baseline.

The second condition is length. Starting the workers costs time, paid once per call, so the scan has to be long enough to absorb this overhead. Measured on the scan above, four workers return :math:`2.4` times the speed of one process at :math:`5\,000` points, :math:`2.3` at :math:`2\,000` and :math:`1.8` at :math:`500`. A scan of a few hundred points called once can finish slower than it would have in a single process. :doc:`/performance` and :doc:`/performance` measure how that gain grows with the size of the scan and where it stops paying.

Parallelism changes more than the speed. A serial scan warm-starts each point from its neighbor. A parallel scan warm-starts every point from the first one. A ladder that starts on a different rung can stop on a different rung, so two runs that differ only in ``n_jobs`` agree to the tolerance asked for, but no better. The gap between them grows with the length of the scan. Measured on the chord above at ``rtol`` :math:`= 10^{-6}`, two points agree exactly, three differ by :math:`8 \cdot 10^{-10}`, and forty by :math:`10^{-8}`. Hold ``n_jobs`` fixed when two runs have to match to the last digit.

All of the above concerns a scan of several probabilities, which is the only thing ``n_jobs`` affects. A call for a single probability accepts it and ignores it. The slabs of one calculation go to the batched kernel in a single call, and no worker is started:

.. code-block:: python

   # Many points: each worker receives points
   P = oscprob.osc_prob_3nu_earth(
       E, L=L, n_jobs=4, **kw)

   # One point: n_jobs is accepted and does
   # nothing; the call runs in one process
   P = oscprob.osc_prob_3nu_earth(
       10.0*gd.UNIT_GEV, L=11467.8*gd.UNIT_KM,
       n_jobs=4, **kw)


(Distributing slabs was tried and retired. A slab is far smaller than a probability, and evaluating all of them in one batched call beats handing them out one at a time.)
