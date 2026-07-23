"""
fmlip_relay.backends.dummy
-------------------------
Random-number backend for integration testing.
Requires no external dependencies beyond NumPy.
"""

import numpy as np
from .base import BackendBase


class DummyBackend(BackendBase):
    """
    Returns plausible-looking but meaningless random energies and forces.
    Useful for validating the Fortran ↔ socket ↔ Python plumbing without
    needing any real ML potential installed.
    """

    def __init__(self, seed: int = 42):
        self._rng = np.random.default_rng(seed)

    def compute(self, atomic_numbers, positions, cell, pbc, compute_stress, charge=0, spin=1):
        natoms = positions.shape[0]
        energy = float(self._rng.uniform(-10.0, -1.0))
        forces = self._rng.uniform(-0.5, 0.5, (natoms, 3))
        stress = np.zeros((3, 3))
        return energy, forces, stress

    @property
    def name(self) -> str:
        return "dummy"
