.. _ex-sec-cli:

Magνs from the command line
---------------------------

Magνs installs a console script, ``magnus``, that computes one probability, or one probability matrix, without writing Python. It chooses one of the wrappers of :doc:`/functions` through three options, ``--flavors``, ``--environment``, and ``--scenario``. The physical inputs are options named after the arguments of the wrapper, e.g., ``--energy``, ``--baseline``, and ``--rho``, and the oscillation parameters default to the NuFIT 6.1 values. Listing :ref:`The command-line interface <ex-lst-cli>` shows two calls: the probability matrix of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, in matter of constant density, and the :math:`\nu_\mu \to \nu_e` probability along the chord from Fermilab to Homestake. ``magnus --help`` and the documentation site  :cite:p:`MagnusDocs` list every option.

.. _ex-lst-cli:

**The command-line interface.** Two probabilities computed from the command line. The first call returns the probability matrix of Listing :ref:`Vacuum and constant density, three ways <ex-lst-constant>`, with a row for each initial flavor and a column for each final flavor. The second returns one channel along the chord from Fermilab to Homestake, whose length is computed from the two named locations. The oscillation parameters take their NuFIT 6.1 defaults.

.. code-block:: bash

   $ magnus --flavors 3 --environment matter \
       --density-profile constant --rho 3.0 \
       --energy 1 --baseline 1000
   Magνs 1.2.0 -- osc_prob_3nu_matter_constant_density
   E = 1 GeV, L = 1000 km

               nu_e   nu_mu  nu_tau
   nu_e      0.9853  0.0135  0.0012
   nu_mu     0.0136  0.9863  0.0001
   nu_tau    0.0011  0.0002  0.9987

   $ magnus --flavors 3 --environment earth \
       --loc-ini fermilab --loc-fin homestake \
       --energy 2.5 --nu-i mu --nu-f e
   Magνs 1.2.0 -- osc_prob_3nu_earth
   E = 2.5 GeV

   P = 0.0724


Every environment except the Earth needs ``--baseline``. Through the Earth, the pair ``--loc-ini`` and ``--loc-fin`` fixes both the chord and its length, as in the second call of Listing :ref:`The command-line interface <ex-lst-cli>`. The option ``--costhz`` fixes only the direction of the chord, so it needs ``--baseline`` beside it, or ``--detector-depth`` or ``--source-depth``, from which the length is computed.

The numerical settings are options too: ``--rtol`` and ``--atol`` set the tolerances of :doc:`/methodology`, ``--magnus-exp-order`` and ``--integration-method`` the order and the quadrature rule of :doc:`/methodology`, and ``--strategy`` the ``strategy`` keyword of :doc:`/engines`. The tolerances default to :math:`10^{-3}`, as in the Python functions, whereas most figures of this paper use :math:```rtol`` = 10^{-8}` and :math:```atol`` = 10^{-10}`. For the second call of Listing :ref:`The command-line interface <ex-lst-cli>`, tightening them to those values changes the probability by :math:`2 \cdot 10^{-10}`.

Adding ``--json`` prints the result as a JSON object instead, for use in a pipeline (Listing :ref:`Command-line output as JSON <ex-lst-cli-json>`). The object names the function that computed the result and the options that selected it, and gives the energy and the baseline in natural units.

.. _ex-lst-cli-json:

**Command-line output as JSON.** The :math:`\nu_\mu \to \nu_e` entry of the matrix of Listing :ref:`The command-line interface <ex-lst-cli>`, printed as JSON: the function that computed it, the options that selected that function, the energy and the baseline in natural units, and the probability.

.. code-block:: bash

   $ magnus --flavors 3 --environment matter \
       --density-profile constant --rho 3.0 \
       --energy 1 --baseline 1000 \
       --nu-i mu --nu-f e --json
   {
     "function": "osc_prob_3nu_matter_constant_density",
     "flavors": 3,
     "environment": "matter",
     "scenario": "std",
     "nubar": false,
     "energy_eV": 1000000000.0,
     "baseline_eV-1": 5067730000000.0,
     "probability": 0.013589091782265706
   }


The script computes one probability at a time, at one energy and one baseline. It offers no scans, no averaged probability, no density profile read from a file, and no Hamiltonian of one’s own; those are done in Python, with the calls of the preceding sections. The script is meant for a quick number, a check against another code, or a probability inside a shell pipeline. (Where the console script is not on the path, ``python -m magnus`` runs the same program.)
