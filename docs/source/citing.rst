How to cite
===========

If Magνs contributed to work you are publishing, please cite it. A citation is
what makes the effort of maintaining a package visible, and it lets a reader
reproduce what you did. Please cite both the paper that describes Magνs and the
version of the software you used.

Cite the paper
--------------

The paper, `arXiv:2610.07159 <https://arxiv.org/abs/2610.07159>`_, describes the
method, its implementation, and its validation. Its entry on `INSPIRE
<https://inspirehep.net/literature/3212165>`_ is

.. code-block:: bibtex

   @article{Bustamante:2026yce,
       author = "Bustamante, Mauricio",
       title = "{Magnus: neutrino oscillation probabilities for any Hermitian Hamiltonian, any number of flavors, and any matter profile}",
       eprint = "2610.07159",
       archivePrefix = "arXiv",
       primaryClass = "hep-ph",
       month = "10",
       year = "2026"
   }

After the paper appears in a journal, take the updated entry from INSPIRE.

Cite the software
-----------------

Cite also the version of the software you used, since results can depend on it.
Magνs records its own version:

.. code-block:: python

   import magnus
   print(magnus.__version__)

.. The version below is |release|, which conf.py reads from pyproject.toml -- the one place
   the number is written.  A parsed-literal rather than a bibtex code-block, because a
   code-block does not expand substitutions and this page would carry a second copy of it.

.. parsed-literal::

   @software{Magnus,
     author  = {Bustamante, Mauricio},
     title   = {{Mag$\\nu$s: neutrino oscillation probabilities for any
                 Hermitian Hamiltonian, any number of flavors, and any
                 matter profile}},
     doi     = {10.5281/zenodo.23160508},
     url     = {https://github.com/mbustama/Magnus},
     version = {|release|},
     year    = {2026}
   }

Replace ``version`` with the one you used.

The DOI above, `10.5281/zenodo.23160508 <https://doi.org/10.5281/zenodo.23160508>`_, is the
Zenodo concept DOI: it resolves to the latest version.  Each version also has a DOI of its
own, listed on that page; version 1.2.0 is
`10.5281/zenodo.23160509 <https://doi.org/10.5281/zenodo.23160509>`_.

What to say in the text
-----------------------

State enough for a reader to know what was computed and how precisely. In
practice, that is three things:

#. **The version**, as above.
#. **The tolerance you asked for** (``rtol``/``atol``, or the fixed ``n_slabs``
   and ``n_tpts_per_slab`` if you disabled the refinement). Note that these are
   a stopping rule rather than a bound on the error — see
   :ref:`what-rtol-atol-control` — so quoting them describes the *request*, not
   the achieved accuracy.
#. **The strategy**, if you did not use the default. The strategies ``'auto'``
   and ``'magnus'`` can differ by far more than the tolerance on solar
   configurations, so it is worth stating which one produced the numbers.

If accuracy is central to your result, pass ``convergence_info={}`` to the call
and read the dict afterwards: it reports what the refinement ladder did,
including whether it converged or hit a cap.

Citing the method
-----------------

The Magnus expansion itself, and the Gauss–Legendre collocation integrators
Magνs uses by default, are due to others. The :doc:`references` page has the
full bibliography; the two worth citing alongside the software are the review by
Blanes, Casas, Oteo and Ros :cite:p:`Blanes2009`, and the high-order integrators
of Blanes, Casas and Ros :cite:p:`Blanes2000`.

Related software
----------------

If your Hamiltonian is constant or piecewise constant, `NuOscProbExact
<https://github.com/mbustama/NuOscProbExact>`_ solves that case in closed form
and has its own citation; see :doc:`comparison` for when to use it.
