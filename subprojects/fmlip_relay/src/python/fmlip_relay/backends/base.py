"""
fmlip_relay.backends.base
------------------------
Abstract base class that every backend must implement.

A backend is responsible for one thing only: given an atomic configuration,
return energy, forces, and (optionally) stress.  All I/O and process
management is handled by the server layer.
"""

from abc import ABC, abstractmethod
import numpy as np


class BackendBase(ABC):
    """
    Minimal interface every fmlip_relay backend must satisfy.

    Subclasses should load their model in ``__init__`` so that loading
    happens exactly once at server startup.
    """

    @abstractmethod
    def compute(self,
                atomic_numbers: np.ndarray,
                positions:      np.ndarray,
                cell:           np.ndarray,
                pbc:            np.ndarray,
                compute_stress: bool,
                charge: int,
                spin: int,
                ) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Evaluate the potential.

        Parameters
        ----------
        atomic_numbers : (N,)   int32  – atomic numbers (Z)
        positions      : (N, 3) float64 – Cartesian coordinates, Angstrom
        cell           : (3, 3) float64 – lattice row vectors, Angstrom
        pbc            : (3,)   bool    – periodic boundary flags
        compute_stress : bool           – whether to compute the stress tensor
        charge         : int32          - total charge (neutral = 0)
        spin           : int32          - spin multiplicity 2S+1
                                          (ground-state singlet = 1)

        Returns
        -------
        energy : float         – potential energy, eV
        forces : (N, 3) float64 – forces,  eV / Angstrom
        stress : (3, 3) float64 – stress tensor, eV / Angstrom^3
                                  (zero matrix when compute_stress is False)
        """

    @property
    def name(self) -> str:
        """Human-readable backend name for log messages."""
        return type(self).__name__

    @property
    def device(self) -> str | None:
        """Compute device (e.g. ``"cpu"``, ``"cuda:0"``), or None if N/A."""
        return getattr(self, "_device", None)

    @property
    def precision(self) -> str | None:
        """Floating-point precision (e.g. ``"float32"``), or None if N/A."""
        return getattr(self, "_dtype", None)
