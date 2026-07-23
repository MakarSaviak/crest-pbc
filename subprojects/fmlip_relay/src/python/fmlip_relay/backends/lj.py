"""
fmlip_relay.backends.lj
----------------------
Lennard-Jones pair potential with periodic boundary conditions.

    E = 4ε Σ_{i<j} [ (σ/r)^12 − (σ/r)^6 ]

PBC is handled via the minimum-image convention over all 27 (or fewer,
for non-periodic axes) cell images.  A smooth cutoff is applied:

    E_cut(r) = E(r) − E(r_cut) − (r − r_cut) · dE/dr|_{r_cut}

so that both energy and forces go to zero at r = r_cut continuously.

No external dependencies beyond NumPy.

Default parameters approximate Argon (Lennard-Jones 1924):
  ε = 0.0104 eV   (≈ 120 K in units of k_B)
  σ = 3.40  Å
  r_cut = 2.5 σ

Usage (via CLI):
    fmlip-relay-server --port 54321 --backend lj
    fmlip-relay-server --port 54321 --backend lj \\
        --lj-epsilon 0.0104 --lj-sigma 3.40 --lj-cutoff 8.5
"""

from __future__ import annotations

import numpy as np
from .base import BackendBase

# Default Argon-like parameters
_DEFAULT_EPSILON = 0.0104   # eV
_DEFAULT_SIGMA   = 3.40     # Å
_DEFAULT_RCUT    = 8.50     # Å  (≈ 2.5 σ)


class LJBackend(BackendBase):
    """
    Lennard-Jones potential with PBC and a smooth force-shifted cutoff.

    Parameters
    ----------
    epsilon : float
        Well depth ε in eV.
    sigma : float
        Zero-crossing distance σ in Å.
    r_cut : float
        Cutoff radius in Å.  Interactions beyond r_cut are ignored.
    """

    def __init__(self,
                 epsilon: float = _DEFAULT_EPSILON,
                 sigma:   float = _DEFAULT_SIGMA,
                 r_cut:   float = _DEFAULT_RCUT):
        if r_cut <= 0:
            raise ValueError("r_cut must be positive")
        if sigma <= 0 or epsilon <= 0:
            raise ValueError("sigma and epsilon must be positive")

        self._eps   = float(epsilon)
        self._sig   = float(sigma)
        self._rcut  = float(r_cut)

        # Pre-compute shift terms so we can apply them cheaply per pair
        inv_rc  = sigma / r_cut
        inv_rc6 = inv_rc ** 6
        inv_rc12 = inv_rc6 ** 2
        self._e_shift = 4.0 * epsilon * (inv_rc12 - inv_rc6)

        # dE/dr at r_cut (negative for attractive region)
        self._f_shift = 4.0 * epsilon * (-12.0 * inv_rc12 + 6.0 * inv_rc6) / r_cut

    # ── public ────────────────────────────────────────────────────────────────

    def compute(self,
                atomic_numbers: np.ndarray,
                positions:      np.ndarray,
                cell:           np.ndarray,
                pbc:            np.ndarray,
                compute_stress: bool,
                charge:         int=0,
                spin:           int=1,
                ) -> tuple[float, np.ndarray, np.ndarray]:

        natoms = positions.shape[0]
        forces = np.zeros((natoms, 3), dtype=np.float64)
        stress = np.zeros((3, 3),      dtype=np.float64)
        energy = 0.0

        # Build list of image translation vectors to consider
        images = _image_vectors(cell, pbc)

        rcut2 = self._rcut ** 2

        for i in range(natoms - 1):
            # Displacement vectors from atom i to all j > i, all images
            # shape: (natoms-i-1, 1, 3) + (1, nimages, 3) → (natoms-i-1, nimages, 3)
            dR0 = positions[i+1:] - positions[i]          # (nj, 3)
            dR  = dR0[:, None, :] + images[None, :, :]    # (nj, nimages, 3)

            r2  = np.einsum('...k,...k->...', dR, dR)     # (nj, nimages)

            # Mask: within cutoff and (for zero image) i≠j is already guaranteed
            mask = r2 < rcut2                              # (nj, nimages)

            if not np.any(mask):
                continue

            r2m   = r2[mask]
            dRm   = dR[mask]                               # (npairs, 3)
            js    = np.where(mask)[0] + i + 1              # atom indices

            # LJ energy and scalar force magnitude / r
            inv_r2  = self._sig**2 / r2m
            inv_r6  = inv_r2 ** 3
            inv_r12 = inv_r6 ** 2

            e_pair = 4.0 * self._eps * (inv_r12 - inv_r6)
            # Force-shift correction: E → E - E(rc) - (r-rc)·F(rc)
            r_dist = np.sqrt(r2m)
            e_pair -= self._e_shift + self._f_shift * (r_dist - self._rcut)
            energy += e_pair.sum()

            # f_scalar = (dU_shifted/dr) / r = (dU/dr - f_shift) / r
            # dU/dr = 4ε(-12σ¹²/r¹³ + 6σ⁶/r⁷) = 4ε(-12·inv_r12 + 6·inv_r6) / r
            f_scalar = (4.0 * self._eps * (-12.0 * inv_r12 + 6.0 * inv_r6) / r2m
                        - self._f_shift / r_dist)         # (npairs,)

            f_vec = f_scalar[:, None] * dRm               # (npairs, 3)

            # Accumulate forces (Newton's third law)
            np.add.at(forces, i,   f_vec.sum(axis=0))
            np.add.at(forces, js, -f_vec)

            # Stress tensor: σ_ab = -1/V Σ r_a F_b  (virial, outer product)
            if compute_stress:
                # outer product summed over pairs: shape (3,3)
                stress -= np.einsum('pi,pj->ij', dRm, f_vec)

        if compute_stress:
            vol = abs(np.linalg.det(cell))
            if vol > 1e-12:
                stress /= vol

        return float(energy), forces, stress

    @property
    def name(self) -> str:
        return f"lj(eps={self._eps:.4f} eV, sig={self._sig:.3f} Å, rc={self._rcut:.3f} Å)"


# ── helpers ───────────────────────────────────────────────────────────────────

def _image_vectors(cell: np.ndarray, pbc: np.ndarray) -> np.ndarray:
    """
    Return an (nimages, 3) array of lattice translation vectors
    for the minimum-image search.

    For each axis: if periodic, consider images -1, 0, +1;
                   if not,      consider only 0.
    """
    ranges = [[-1, 0, 1] if pbc[k] else [0] for k in range(3)]
    offsets = np.array(
        [[i, j, k] for i in ranges[0] for j in ranges[1] for k in ranges[2]],
        dtype=np.float64,
    )                                    # (nimages, 3) in fractional coords
    return offsets @ cell                # (nimages, 3) in Cartesian Å
