"""
fmlip_relay.backends._ase_base
------------------------------
Generic mixin that implements ``compute()`` on top of any ASE calculator
stored in ``self._calc``.  Subclasses only need to construct ``self._calc``
in ``__init__``; everything else is handled here.

Not part of the public API.
"""

from __future__ import annotations

import numpy as np
from .base import BackendBase


class _ASEComputeMixin(BackendBase):
    """
    Mixin that implements ``compute()`` on top of a stored ASE calculator
    ``self._calc``.  Subclasses only need to construct ``self._calc`` in
    ``__init__``.
    """

    _calc = None   # must be set by subclass __init__

    def compute(self,
                atomic_numbers: np.ndarray,
                positions:      np.ndarray,
                cell:           np.ndarray,
                pbc:            np.ndarray,
                compute_stress: bool,
                charge:         int,
                spin:           int,
                ) -> tuple[float, np.ndarray, np.ndarray]:
        from ase import Atoms

        atoms = Atoms(
            numbers=atomic_numbers,
            positions=positions,
            cell=cell,
            pbc=pbc,
        )
        # Spin is the multiplicity (2S+1), which is >= 1; a stray 0 (e.g. an
        # uninitialised caller default) collapses to the ground-state singlet.
        atoms.info["charge"] = charge
        atoms.info["spin"] = spin if spin >= 1 else 1
        atoms.calc = self._calc

        energy = float(atoms.get_potential_energy())
        forces = atoms.get_forces().astype(np.float64)
        stress = (
            atoms.get_stress(voigt=False).astype(np.float64)
            if compute_stress
            else np.zeros((3, 3), dtype=np.float64)
        )
        return energy, forces, stress
