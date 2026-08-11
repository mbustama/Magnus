# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""__init__.py

Subpackage initializer for magnus.hamiltonians. Explicitly imports and
re-exports the public names of hamiltonians2nu.py, hamiltonians3nu.py,
hamiltonians4nu.py, and hamiltonians5nu.py: the mixing matrices and
vacuum/matter/NSI/LIV Hamiltonians for 2, 3, 4, and 5 neutrino flavors.

Routine listings
----------------

    (none; only re-exports the four hamiltonians{2,3,4,5}nu.py modules'
    public names -- see __all__ below for the exact list)
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from .hamiltonians2nu import (
    mixing_matrix_2nu,
    hamiltonian_2nu_vacuum_energy_independent,
    hamiltonian_2nu_vacuum_energy_independent_td,
    hamiltonian_2nu_vacuum,
    hamiltonian_2nu_vacuum_td,
    hamiltonian_2nu_matter,
    hamiltonian_2nu_matter_td,
    hamiltonian_2nu_nsi,
    hamiltonian_2nu_nsi_td,
    hamiltonian_2nu_liv,
    hamiltonian_2nu_liv_energy_independent,
)
from .hamiltonians3nu import (
    pmns_mixing_matrix,
    mixing_matrix_3x3,
    hamiltonian_3nu_vacuum_energy_independent,
    hamiltonian_3nu_vacuum_energy_independent_td,
    hamiltonian_3nu_vacuum,
    hamiltonian_3nu_vacuum_td,
    hamiltonian_3nu_matter,
    hamiltonian_3nu_matter_td,
    hamiltonian_3nu_nsi,
    hamiltonian_3nu_nsi_td,
    hamiltonian_3nu_liv,
    hamiltonian_3nu_liv_energy_independent,
)
from .hamiltonians4nu import (
    mixing_matrix_4x4,
    hamiltonian_4nu_vacuum_energy_independent,
    hamiltonian_4nu_vacuum_energy_independent_td,
    hamiltonian_4nu_vacuum,
    hamiltonian_4nu_vacuum_td,
    hamiltonian_4nu_matter,
    hamiltonian_4nu_matter_td,
    hamiltonian_4nu_nsi,
    hamiltonian_4nu_nsi_td,
    hamiltonian_4nu_liv,
    hamiltonian_4nu_liv_energy_independent,
)
from .hamiltonians5nu import (
    mixing_matrix_5x5,
    hamiltonian_5nu_vacuum_energy_independent,
    hamiltonian_5nu_vacuum_energy_independent_td,
    hamiltonian_5nu_vacuum,
    hamiltonian_5nu_vacuum_td,
    hamiltonian_5nu_matter,
    hamiltonian_5nu_matter_td,
    hamiltonian_5nu_nsi,
    hamiltonian_5nu_liv,
    hamiltonian_5nu_liv_energy_independent,
)

__all__ = [
    # hamiltonians2nu
    'mixing_matrix_2nu',
    'hamiltonian_2nu_vacuum_energy_independent',
    'hamiltonian_2nu_vacuum_energy_independent_td',
    'hamiltonian_2nu_vacuum',
    'hamiltonian_2nu_vacuum_td',
    'hamiltonian_2nu_matter',
    'hamiltonian_2nu_matter_td',
    'hamiltonian_2nu_nsi',
    'hamiltonian_2nu_nsi_td',
    'hamiltonian_2nu_liv',
    'hamiltonian_2nu_liv_energy_independent',
    # hamiltonians3nu
    'pmns_mixing_matrix',
    'mixing_matrix_3x3',
    'hamiltonian_3nu_vacuum_energy_independent',
    'hamiltonian_3nu_vacuum_energy_independent_td',
    'hamiltonian_3nu_vacuum',
    'hamiltonian_3nu_vacuum_td',
    'hamiltonian_3nu_matter',
    'hamiltonian_3nu_matter_td',
    'hamiltonian_3nu_nsi',
    'hamiltonian_3nu_nsi_td',
    'hamiltonian_3nu_liv',
    'hamiltonian_3nu_liv_energy_independent',
    # hamiltonians4nu
    'mixing_matrix_4x4',
    'hamiltonian_4nu_vacuum_energy_independent',
    'hamiltonian_4nu_vacuum_energy_independent_td',
    'hamiltonian_4nu_vacuum',
    'hamiltonian_4nu_vacuum_td',
    'hamiltonian_4nu_matter',
    'hamiltonian_4nu_matter_td',
    'hamiltonian_4nu_nsi',
    'hamiltonian_4nu_nsi_td',
    'hamiltonian_4nu_liv',
    'hamiltonian_4nu_liv_energy_independent',
    # hamiltonians5nu
    'mixing_matrix_5x5',
    'hamiltonian_5nu_vacuum_energy_independent',
    'hamiltonian_5nu_vacuum_energy_independent_td',
    'hamiltonian_5nu_vacuum',
    'hamiltonian_5nu_vacuum_td',
    'hamiltonian_5nu_matter',
    'hamiltonian_5nu_matter_td',
    'hamiltonian_5nu_nsi',
    'hamiltonian_5nu_liv',
    'hamiltonian_5nu_liv_energy_independent',
]
