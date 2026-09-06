How to cite
============

If Magνs contributed to work you are publishing, please cite it. A citation is
what makes the effort of maintaining a package visible, and it lets a reader
reproduce what you did.

Cite the software
------------------

Until the accompanying paper is available, cite the software itself, including
the version you used — results can depend on it, and Magνs records its own:

.. code-block:: python

   import magnus
   print(magnus.__version__)

.. code-block:: bibtex

   @software{Magnus,
     author  = {Bustamante, Mauricio},
     title   = {{Mag$\nu$s: neutrino oscillation probabilities via the
                 Magnus expansion}},
     url     = {https://github.com/mbustama/Magnus},
     version = {1.0.0rc1},
     year    = {2026}
   }

Replace ``version`` with the one you actually ran.

What to say in the text
------------------------

Enough for a reader to know what was computed and how precisely. In practice
that is three things:

#. **The version**, as above.
#. **The tolerance you asked for** (``rtol``/``atol``, or the fixed ``n_slabs``
   and ``n_tpts_per_slab`` if you disabled the refinement). Note that these are
   a stopping criterion rather than a bound on the error — see
   :ref:`what-rtol-atol-control` — so quoting them describes the *request*, not
   the achieved accuracy.
#. **The strategy**, if you did not use the default. ``strategy='auto'`` and
   ``strategy='magnus'`` can differ by far more than the tolerance on solar
   configurations, so it is worth stating which one produced the numbers.

If accuracy is central to your result, ``convergence_info`` reports what the
refinement ladder actually did, including whether it converged or hit a cap.

Citing the method
------------------

The Magnus expansion itself, and the Gauss–Legendre collocation integrators
Magνs uses by default, are due to others. The :doc:`references` page has the
full bibliography; the two worth citing alongside the software are the review by
Blanes, Casas, Oteo and Ros, and the high-order integrators of Blanes, Casas and
Ros.

Related software
-----------------

If your Hamiltonian is constant or piecewise constant, `NuOscProbExact
<https://github.com/mbustama/NuOscProbExact>`_ solves that case in closed form
and has its own citation; see :ref:`use-nuoscprobexact-instead`.
