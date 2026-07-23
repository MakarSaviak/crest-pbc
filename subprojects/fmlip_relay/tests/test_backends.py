"""
tests/test_backends.py
----------------------
Tests for:
  - BackendBase contract (cannot instantiate abstract class)
  - DummyBackend output shapes, dtypes, and reproducibility
  - LJBackend constructor validation, physics correctness
"""

import numpy as np
import pytest

from fmlip_relay.backends.base  import BackendBase
from fmlip_relay.backends.dummy import DummyBackend
from fmlip_relay.backends.lj    import LJBackend, _image_vectors


# ── BackendBase ───────────────────────────────────────────────────────────────

class TestBackendBase:
    def test_cannot_instantiate(self):
        with pytest.raises(TypeError):
            BackendBase()

    def test_name_default(self):
        """Concrete subclass without name override uses class name."""
        class Minimal(BackendBase):
            def compute(self, *a, **kw):
                return 0.0, np.zeros((1,3)), np.zeros((3,3))
        assert Minimal().name == "Minimal"


# ── DummyBackend ──────────────────────────────────────────────────────────────

class TestDummyBackend:
    @pytest.fixture
    def backend(self):
        return DummyBackend(seed=0)

    def _call(self, backend, natoms=4, compute_stress=False):
        Z   = np.ones(natoms, dtype=np.int32)
        pos = np.zeros((natoms, 3))
        cell = np.eye(3) * 10.0
        pbc  = np.array([False, False, False])
        return backend.compute(Z, pos, cell, pbc, compute_stress)

    def test_output_shapes(self, backend):
        N = 7
        energy, forces, stress = self._call(backend, natoms=N)
        assert isinstance(energy, float)
        assert forces.shape == (N, 3)
        assert stress.shape == (3, 3)

    def test_forces_dtype(self, backend):
        _, forces, stress = self._call(backend, natoms=3)
        assert forces.dtype == np.float64
        assert stress.dtype == np.float64

    def test_stress_zeros(self, backend):
        _, _, stress = self._call(backend, compute_stress=False)
        np.testing.assert_array_equal(stress, np.zeros((3, 3)))

    def test_reproducibility_with_same_seed(self):
        b1 = DummyBackend(seed=99)
        b2 = DummyBackend(seed=99)
        e1, f1, _ = self._call(b1)
        e2, f2, _ = self._call(b2)
        assert e1 == e2
        np.testing.assert_array_equal(f1, f2)

    def test_different_seeds_give_different_results(self):
        b1 = DummyBackend(seed=1)
        b2 = DummyBackend(seed=2)
        e1, _, _ = self._call(b1)
        e2, _, _ = self._call(b2)
        assert e1 != e2

    def test_name(self, backend):
        assert backend.name == "dummy"


# ── LJBackend – constructor ───────────────────────────────────────────────────

class TestLJBackendConstructor:
    def test_defaults_are_argon(self):
        b = LJBackend()
        assert b._eps  == pytest.approx(0.0104)
        assert b._sig  == pytest.approx(3.40)
        assert b._rcut == pytest.approx(8.50)

    def test_negative_epsilon_raises(self):
        with pytest.raises(ValueError):
            LJBackend(epsilon=-1.0)

    def test_zero_sigma_raises(self):
        with pytest.raises(ValueError):
            LJBackend(sigma=0.0)

    def test_negative_cutoff_raises(self):
        with pytest.raises(ValueError):
            LJBackend(r_cut=-1.0)

    def test_name_contains_params(self):
        b = LJBackend(epsilon=0.01, sigma=3.0, r_cut=7.5)
        assert "lj(" in b.name
        assert "eV"  in b.name
        assert "Å"   in b.name


# ── LJBackend – output contract ───────────────────────────────────────────────

class TestLJBackendOutputContract:
    @pytest.fixture
    def lj(self):
        return LJBackend()

    def test_output_shapes(self, lj, argon_fcc):
        Z, pos, cell, pbc = argon_fcc
        N = len(pos)
        energy, forces, stress = lj.compute(Z, pos, cell, pbc, True)
        assert isinstance(energy, float)
        assert forces.shape == (N, 3)
        assert stress.shape == (3, 3)

    def test_forces_dtype_float64(self, lj, argon_fcc):
        Z, pos, cell, pbc = argon_fcc
        _, forces, stress = lj.compute(Z, pos, cell, pbc, False)
        assert forces.dtype == np.float64
        assert stress.dtype == np.float64

    def test_stress_zero_when_not_requested(self, lj, argon_fcc):
        Z, pos, cell, pbc = argon_fcc
        _, _, stress = lj.compute(Z, pos, cell, pbc, compute_stress=False)
        np.testing.assert_array_equal(stress, np.zeros((3, 3)))


# ── LJBackend – physics ───────────────────────────────────────────────────────

class TestLJBackendPhysics:
    """
    Physical correctness tests.  We use a small epsilon / sigma so that
    the numbers are easy to reason about.
    """

    @pytest.fixture
    def lj(self):
        return LJBackend(epsilon=1.0, sigma=1.0, r_cut=4.0)

    # ── Two-atom analytical comparison ───────────────────────────────────────

    def test_two_atoms_energy_analytic(self, lj):
        """
        For two atoms at separation r, the force-shifted energy is:
            E = 4ε[(σ/r)^12 - (σ/r)^6] - E(rc) - (r-rc)·dE/dr|rc
        Compare against manually computed value at r = 1.5σ.
        """
        eps, sig, rc = 1.0, 1.0, 4.0
        r = 1.5

        # analytic unshifted
        inv_r  = sig / r
        e_raw  = 4 * eps * (inv_r**12 - inv_r**6)

        # shift terms
        inv_rc  = sig / rc
        e_shift = 4 * eps * (inv_rc**12 - inv_rc**6)
        f_shift = 4 * eps * (-12*inv_rc**12 + 6*inv_rc**6) / rc
        e_analytic = e_raw - e_shift - f_shift * (r - rc)

        pos = np.array([[0.0, 0.0, 0.0],
                        [r,   0.0, 0.0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)

        energy, _, _ = lj.compute(Z, pos, cell, pbc, False)
        assert energy == pytest.approx(e_analytic, rel=1e-10)

    def test_energy_at_cutoff_is_zero(self, lj):
        """Force-shifted potential must give E≈0 at r = r_cut."""
        rc = 4.0
        pos = np.array([[0.0, 0.0, 0.0],
                        [rc,  0.0, 0.0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)
        energy, _, _ = lj.compute(Z, pos, cell, pbc, False)
        assert energy == pytest.approx(0.0, abs=1e-12)

    def test_forces_at_cutoff_are_zero(self, lj):
        """Force-shifted potential must give F≈0 at r = r_cut."""
        rc = 4.0
        pos = np.array([[0.0, 0.0, 0.0],
                        [rc,  0.0, 0.0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)
        _, forces, _ = lj.compute(Z, pos, cell, pbc, False)
        np.testing.assert_allclose(forces, np.zeros((2, 3)), atol=1e-12)

    # ── Newton's third law ────────────────────────────────────────────────────

    def test_newton_third_law_open(self, lj, two_atoms_far):
        """Sum of forces must vanish for any configuration."""
        Z, pos, cell, pbc = two_atoms_far
        _, forces, _ = lj.compute(Z, pos, cell, pbc, False)
        np.testing.assert_allclose(forces.sum(axis=0), np.zeros(3), atol=1e-14)

    def test_newton_third_law_periodic(self, lj, argon_fcc):
        Z, pos, cell, pbc = argon_fcc
        _, forces, _ = lj.compute(Z, pos, cell, pbc, False)
        np.testing.assert_allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)

    def test_newton_third_law_cluster(self, lj, argon_cluster):
        Z, pos, cell, pbc = argon_cluster
        _, forces, _ = lj.compute(Z, pos, cell, pbc, False)
        np.testing.assert_allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)

    # ── Force–energy consistency (finite differences) ─────────────────────────

    def test_forces_match_finite_difference(self, lj):
        """F = -∇E verified by central finite differences on each Cartesian DOF.
        Uses a geometry well within r_cut=4.0 for the custom lj fixture (sigma=1.0).
        """
        # Two atoms at r=1.5σ, comfortably inside r_cut=4.0
        pos  = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)
        dx   = 1e-5

        _, forces, _ = lj.compute(Z, pos.copy(), cell, pbc, False)

        for atom in range(2):
            for coord in range(3):
                p_fwd = pos.copy(); p_fwd[atom, coord] += dx
                p_bwd = pos.copy(); p_bwd[atom, coord] -= dx
                e_fwd, _, _ = lj.compute(Z, p_fwd, cell, pbc, False)
                e_bwd, _, _ = lj.compute(Z, p_bwd, cell, pbc, False)
                fd = -(e_fwd - e_bwd) / (2 * dx)
                assert forces[atom, coord] == pytest.approx(fd, rel=1e-5, abs=1e-10)

    def test_forces_match_finite_difference_periodic(self, argon_fcc):
        """Same finite-difference check for a periodic system with Argon params."""
        lj = LJBackend()   # default Argon: sigma=3.40, r_cut=8.50
        Z, pos, cell, pbc = argon_fcc
        dx = 1e-5

        _, forces, _ = lj.compute(Z, pos.copy(), cell, pbc, False)
        for atom in [0, 5, 10]:
            for coord in range(3):
                p_fwd = pos.copy(); p_fwd[atom, coord] += dx
                p_bwd = pos.copy(); p_bwd[atom, coord] -= dx
                e_fwd, _, _ = lj.compute(Z, p_fwd, cell, pbc, False)
                e_bwd, _, _ = lj.compute(Z, p_bwd, cell, pbc, False)
                fd = -(e_fwd - e_bwd) / (2 * dx)
                assert forces[atom, coord] == pytest.approx(fd, rel=1e-5, abs=1e-10)

    # ── PBC vs open boundary ──────────────────────────────────────────────────

    def test_pbc_differs_from_open_boundary(self, lj, argon_fcc):
        """A periodic system should give a different energy than the same
        cluster with pbc=False (images add interactions)."""
        Z, pos, cell, _ = argon_fcc
        pbc_on  = np.array([True,  True,  True])
        pbc_off = np.array([False, False, False])
        e_pbc,  _, _ = lj.compute(Z, pos, cell, pbc_on,  False)
        e_open, _, _ = lj.compute(Z, pos, cell, pbc_off, False)
        assert e_pbc != pytest.approx(e_open, rel=1e-6)

    # ── Stress tensor properties ──────────────────────────────────────────────

    def test_stress_is_symmetric(self, lj, argon_fcc):
        Z, pos, cell, pbc = argon_fcc
        _, _, stress = lj.compute(Z, pos, cell, pbc, compute_stress=True)
        np.testing.assert_allclose(stress, stress.T, atol=1e-12)

    def test_stress_nonzero_for_periodic(self, lj, argon_fcc):
        Z, pos, cell, pbc = argon_fcc
        _, _, stress = lj.compute(Z, pos, cell, pbc, compute_stress=True)
        assert not np.allclose(stress, np.zeros((3, 3)))

    def test_stress_zero_for_single_pair_open(self, lj, two_atoms_far):
        """Stress is not requested → must be exactly zero."""
        Z, pos, cell, pbc = two_atoms_far
        _, _, stress = lj.compute(Z, pos, cell, pbc, compute_stress=False)
        np.testing.assert_array_equal(stress, np.zeros((3, 3)))

    # ── Energy monotonicity near equilibrium ─────────────────────────────────

    def test_energy_minimum_near_r_min(self, lj):
        """LJ minimum is at r = 2^(1/6) σ.  Check E(r_min) < E(r_min ± δ)."""
        r_min = 2**(1/6)   # in units of σ=1
        delta = 0.05
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)

        def E(r):
            pos = np.array([[0.0, 0.0, 0.0], [r, 0.0, 0.0]])
            e, _, _ = lj.compute(Z, pos, cell, pbc, False)
            return e

        assert E(r_min) < E(r_min - delta)
        assert E(r_min) < E(r_min + delta)

    # ── Two-atom force direction ──────────────────────────────────────────────

    def test_force_direction_repulsive(self, lj):
        """At r < r_min = 2^(1/6)σ ≈ 1.122 atoms repel:
        force on atom 0 (at origin) points away from atom 1 (at +x), i.e. F[0,0] < 0."""
        r = 0.9   # < r_min for sigma=1.0
        pos = np.array([[0.0, 0.0, 0.0], [r, 0.0, 0.0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)
        _, forces, _ = lj.compute(Z, pos, cell, pbc, False)
        assert forces[0, 0] < 0.0   # atom 0 pushed in -x (away from atom 1)

    def test_force_direction_attractive(self, lj):
        """At r > r_min = 2^(1/6)σ ≈ 1.122 atoms attract:
        force on atom 0 (at origin) points toward atom 1 (at +x), i.e. F[0,0] > 0."""
        r = 1.5   # > r_min for sigma=1.0
        pos = np.array([[0.0, 0.0, 0.0], [r, 0.0, 0.0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        Z    = np.array([1, 1], dtype=np.int32)
        _, forces, _ = lj.compute(Z, pos, cell, pbc, False)
        assert forces[0, 0] > 0.0   # atom 0 pulled in +x (toward atom 1)


# ── _image_vectors ────────────────────────────────────────────────────────────

class TestImageVectors:
    def test_all_pbc_gives_27_images(self):
        cell = np.eye(3) * 5.0
        pbc  = np.array([True, True, True])
        imgs = _image_vectors(cell, pbc)
        assert imgs.shape == (27, 3)

    def test_no_pbc_gives_1_image(self):
        cell = np.eye(3) * 5.0
        pbc  = np.array([False, False, False])
        imgs = _image_vectors(cell, pbc)
        assert imgs.shape == (1, 3)
        np.testing.assert_array_equal(imgs[0], [0.0, 0.0, 0.0])

    def test_one_periodic_axis_gives_3_images(self):
        cell = np.eye(3) * 5.0
        pbc  = np.array([True, False, False])
        imgs = _image_vectors(cell, pbc)
        assert imgs.shape == (3, 3)

    def test_two_periodic_axes_gives_9_images(self):
        cell = np.eye(3) * 5.0
        pbc  = np.array([True, True, False])
        imgs = _image_vectors(cell, pbc)
        assert imgs.shape == (9, 3)

    def test_image_vectors_scale_with_cell(self):
        """The ±1 image should be exactly one cell vector away."""
        L = 7.3
        cell = np.eye(3) * L
        pbc  = np.array([True, False, False])
        imgs = _image_vectors(cell, pbc)
        xs = sorted(set(imgs[:, 0]))
        assert xs == pytest.approx([-L, 0.0, L])
