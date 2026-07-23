"""
tests/conftest.py
-----------------
Shared pytest fixtures used across all test modules.
"""

import numpy as np
import pytest


# ── small geometry helpers ────────────────────────────────────────────────────

def _fcc_cell_and_positions(n_cells: int, alat: float):
    """Return (positions, cell) for an FCC supercell."""
    basis = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ])
    pos = []
    for ix in range(n_cells):
        for iy in range(n_cells):
            for iz in range(n_cells):
                for b in basis:
                    pos.append((b + [ix, iy, iz]) * alat)
    cell = np.eye(3) * n_cells * alat
    return np.array(pos, dtype=np.float64), cell


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def two_atoms_far():
    """Two atoms well within the default LJ cutoff (r = 4.0 Å ≈ σ*1.18)."""
    positions = np.array([[0.0, 0.0, 0.0],
                          [4.0, 0.0, 0.0]], dtype=np.float64)
    cell      = np.eye(3) * 100.0          # large box → effectively open BC
    pbc       = np.array([False, False, False])
    atomic_numbers = np.array([18, 18], dtype=np.int32)   # Ar
    return atomic_numbers, positions, cell, pbc


@pytest.fixture
def two_atoms_at_sigma():
    """Two atoms placed exactly at σ (E_LJ = 0 before shift)."""
    sigma = 3.40
    positions = np.array([[0.0, 0.0, 0.0],
                          [sigma, 0.0, 0.0]], dtype=np.float64)
    cell      = np.eye(3) * 100.0
    pbc       = np.array([False, False, False])
    atomic_numbers = np.array([18, 18], dtype=np.int32)
    return atomic_numbers, positions, cell, pbc


@pytest.fixture
def argon_fcc():
    """64-atom Argon FCC supercell (2×2×2 unit cells, a = 5.26 Å)."""
    pos, cell = _fcc_cell_and_positions(n_cells=2, alat=5.26)
    atomic_numbers = np.full(len(pos), 18, dtype=np.int32)
    pbc = np.array([True, True, True])
    return atomic_numbers, pos, cell, pbc


@pytest.fixture
def argon_cluster():
    """Same geometry but with pbc=False (molecular cluster)."""
    pos, cell = _fcc_cell_and_positions(n_cells=2, alat=5.26)
    atomic_numbers = np.full(len(pos), 18, dtype=np.int32)
    pbc = np.array([False, False, False])
    return atomic_numbers, pos, cell, pbc
