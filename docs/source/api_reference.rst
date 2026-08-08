API Reference
==============

Generated from the docstrings, module by module. Every public function,
class and module-level constant Magνs ships appears here with its full
signature, its parameters and what it returns.

This is the exhaustive view. Three other pages are usually the faster way in:

* :doc:`recipes` — what the package can compute, with the code that computes it.
* :doc:`functions` — the whole ``osc_prob_*`` family laid out by environment and
  flavour count, for when you know roughly what you want but not its name.
* :doc:`implementation_details` — how the engines choose between themselves, and
  the population every tuned constant was measured on.

Constants are documented where they live rather than collected into a table of
their own, because each one's docstring carries the measurement that set it —
what was swept, on which workloads, and what the alternatives cost.

.. The module pages are listed individually rather than through
   ``api/magnus/index``: that page is titled after the package, which put a
   redundant "magnus" level between this page and the modules in the sidebar.
   sphinx-autoapi does not keep its generated .rst files, so a :glob: cannot
   reach them either.  A module added without being listed here is reported by
   the -W build as a document that is not in any toctree.

.. toctree::
   :maxdepth: 2

   api/magnus/adiabatic/index
   api/magnus/avgprob/index
   api/magnus/cli/index
   api/magnus/earth/index
   api/magnus/expansionterms/index
   api/magnus/globaldefs/index
   api/magnus/hamiltonians/index
   api/magnus/magnus/index
   api/magnus/matter/index
   api/magnus/oscprob/index
   api/magnus/oscprobstd/index
   api/magnus/plotting/index
