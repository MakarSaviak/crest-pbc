"""
fmlip_relay.test_mode
---------------------
Single-point test evaluation for the ``--test`` CLI flag.

Two public entry points:

    read_geometry(path)          → (atomic_numbers, positions, cell, pbc)
    run_test(backend, path)      → int (exit code: 0 = ok, 1 = error)

Geometry reading tries ASE first (XYZ, extended-XYZ, CIF, POSCAR, …) and
falls back to a self-contained plain-XYZ parser so that ``--test`` works even
without ASE installed (e.g. with the ``dummy`` or ``lj`` backend).
"""

from __future__ import annotations

import logging
import os
import sys
import time

import numpy as np


# ── element table for the plain-XYZ fallback ─────────────────────────────────

_SYMBOL_TO_Z: dict[str, int] = {
    "H": 1,  "He": 2, "Li": 3, "Be": 4, "B": 5,  "C": 6,  "N": 7,
    "O": 8,  "F": 9,  "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14,
    "P":15,  "S":16,  "Cl":17, "Ar":18, "K": 19, "Ca":20, "Sc":21,
    "Ti":22, "V": 23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28,
    "Cu":29, "Zn":30, "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35,
    "Kr":36, "Rb":37, "Sr":38, "Y": 39, "Zr":40, "Nb":41, "Mo":42,
    "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49,
    "Sn":50, "Sb":51, "Te":52, "I": 53, "Xe":54, "Cs":55, "Ba":56,
    "Au":79, "Hg":80, "Pb":82, "Bi":83,
}


# ── public API ────────────────────────────────────────────────────────────────

def read_geometry(
    path: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int | None, int | None]:
    """
    Read a molecular geometry file and return arrays suitable for
    ``BackendBase.compute()``.

    Tries ``ase.io.read`` first (broad format support); falls back to a
    minimal plain-XYZ parser when ASE is not installed.

    Parameters
    ----------
    path : str
        Path to the geometry file.

    Returns
    -------
    atomic_numbers : np.ndarray  shape (N,)   int32
    positions      : np.ndarray  shape (N, 3) float64  [Angstrom]
    cell           : np.ndarray  shape (3, 3) float64  [Angstrom]
    pbc            : np.ndarray  shape (3,)   bool
    charge         : int | None  total charge if present in the file, else None
    spin           : int | None  spin multiplicity (2S+1) if present, else None
    """
    # ── ASE path (preferred) ──────────────────────────────────────────────────
    try:
        import ase.io
        atoms = ase.io.read(path)
        # extended-XYZ and similar formats may carry these in atoms.info
        charge = atoms.info.get("charge")
        spin   = atoms.info.get("spin", atoms.info.get("multiplicity"))
        return (
            atoms.get_atomic_numbers().astype(np.int32),
            atoms.get_positions().astype(np.float64),
            np.array(atoms.get_cell(), dtype=np.float64),
            np.array(atoms.get_pbc(), dtype=bool),
            None if charge is None else int(charge),
            None if spin   is None else int(spin),
        )
    except ImportError:
        pass  # ASE not available — fall through to plain-XYZ parser

    # ── Minimal plain-XYZ fallback ────────────────────────────────────────────
    # Format: line 1 = natoms, line 2 = comment, lines 3… = symbol x y z
    with open(path) as fh:
        lines = [ln.rstrip() for ln in fh if ln.strip()]

    try:
        natoms = int(lines[0])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Cannot parse '{path}' as XYZ: {exc}") from exc

    atom_lines = lines[2 : 2 + natoms]
    if len(atom_lines) < natoms:
        raise ValueError(
            f"XYZ file claims {natoms} atoms but only "
            f"{len(atom_lines)} coordinate lines found."
        )

    symbols, coords = [], []
    for ln in atom_lines:
        parts = ln.split()
        sym = parts[0].capitalize()
        if sym not in _SYMBOL_TO_Z:
            raise ValueError(
                f"Unknown element symbol '{sym}' in plain-XYZ fallback parser. "
                "Install ASE (pip install ase) for broader format support."
            )
        symbols.append(_SYMBOL_TO_Z[sym])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return (
        np.array(symbols, dtype=np.int32),
        np.array(coords,  dtype=np.float64),
        np.zeros((3, 3),  dtype=np.float64),
        np.zeros(3,       dtype=bool),
        None,   # plain-XYZ carries no charge
        None,   # plain-XYZ carries no spin
    )


def run_test(backend, geometry_path: str,
             charge: int | None = None, spin: int | None = None) -> int:
    """
    Load *geometry_path*, call ``backend.compute()`` once, print a summary,
    and return an exit code (0 = success, 1 = error).

    Parameters
    ----------
    backend       : BackendBase instance (already initialised)
    geometry_path : str  path to the input geometry file
    charge        : int | None  total charge override; if None, the value from
                    the geometry file is used, falling back to 0 (neutral).
    spin          : int | None  spin multiplicity (2S+1) override; if None, the
                    value from the geometry file is used, falling back to 1
                    (ground-state singlet).
    """
    log = logging.getLogger(__name__)

    # ── read geometry ─────────────────────────────────────────────────────────
    log.info("Reading geometry from '%s' …", geometry_path)
    try:
        atomic_numbers, positions, cell, pbc, file_charge, file_spin = \
            read_geometry(geometry_path)
    except Exception as exc:
        print(f"\n[ERROR] Could not read geometry: {exc}", file=sys.stderr)
        return 1

    natoms = len(atomic_numbers)
    log.info("Geometry loaded: %d atoms, pbc=%s", natoms, pbc.tolist())

    # ── resolve charge / spin: CLI override > geometry file > default ──────────
    used_charge, charge_src = _resolve(charge, file_charge, 0)
    used_spin,   spin_src   = _resolve(spin,   file_spin,   1)
    if used_spin < 1:
        used_spin = 1   # multiplicity 2S+1 is >= 1; clamp stray 0
    log.info("Using charge=%d (%s), spin=%d (%s)",
             used_charge, charge_src, used_spin, spin_src)

    # ── single-point evaluation ───────────────────────────────────────────────
    compute_stress = bool(np.any(pbc))
    log.info("Running single-point calculation (compute_stress=%s) …", compute_stress)

    t0 = time.perf_counter()
    try:
        energy, forces, stress = backend.compute(
            atomic_numbers, positions, cell, pbc,
            compute_stress=compute_stress,
            charge=used_charge, spin=used_spin,
        )
    except Exception as exc:
        print(f"\n[ERROR] Backend computation failed: {exc}", file=sys.stderr)
        log.exception("Backend raised an exception during test evaluation.")
        return 1
    elapsed = time.perf_counter() - t0

    # ── print summary ─────────────────────────────────────────────────────────
    _print_summary(
        geometry_path, backend, natoms, pbc, cell,
        energy, forces, stress, compute_stress, elapsed,
        atomic_numbers,
        used_charge, charge_src, used_spin, spin_src,
    )
    return 0


def _resolve(override, from_file, default) -> tuple[int, str]:
    """Pick a value by precedence (CLI override > file > default) and label its source."""
    if override is not None:
        return override, "cli"
    if from_file is not None:
        return from_file, "file"
    return default, "default"


# ── internal helpers ──────────────────────────────────────────────────────────

def _print_summary(
    geometry_path: str,
    backend,
    natoms: int,
    pbc: np.ndarray,
    cell: np.ndarray,
    energy: float,
    forces: np.ndarray,
    stress: np.ndarray,
    compute_stress: bool,
    elapsed: float,
    atomic_numbers: np.ndarray,
    charge: int,
    charge_src: str,
    spin: int,
    spin_src: str,
) -> None:
    # Element symbols — ASE if available, raw Z otherwise
    try:
        from ase.data import chemical_symbols
        symbols = [chemical_symbols[z] for z in atomic_numbers]
    except ImportError:
        symbols = [str(z) for z in atomic_numbers]

    fmax  = float(np.linalg.norm(forces, axis=1).max())
    fmean = float(np.linalg.norm(forces, axis=1).mean())
    frms  = float(np.sqrt(np.mean(forces ** 2)))

    sep = "─" * 56
    print(f"\n{sep}")
    print(f"  fmlip-relay  single-point test")
    print(sep)
    print(f"  File         : {os.path.basename(geometry_path)}")
    print(f"  Backend      : {backend.name}")
    device    = backend.device
    precision = backend.precision
    if device is not None:
        print(f"  Device       : {device}")
    if precision is not None:
        print(f"  Precision    : {precision}")
    print(f"  Atoms        : {natoms}")
    print(f"  Charge       : {charge:+d}  ({charge_src})")
    print(f"  Multiplicity : {spin}  (2S+1, {spin_src})")
    print(f"  PBC          : {pbc.tolist()}")
    if compute_stress:
        cell_abc = np.linalg.norm(cell, axis=1)
        print(f"  Cell (|a|)   : {cell_abc[0]:.4f}  {cell_abc[1]:.4f}  {cell_abc[2]:.4f}  Å")
    print(sep)
    print(f"  Energy       : {energy:+.6f}  eV        ({energy/natoms:+.6f} eV/atom)")
    print(f"  |F| max      : {fmax:.6f}  eV/Å")
    print(f"  |F| mean     : {fmean:.6f}  eV/Å")
    print(f"  |F| rms      : {frms:.6f}  eV/Å")
    if compute_stress:
        print(f"  Stress (eV/Å³) :")
        for row in stress:
            print(f"    {row[0]:+.4e}  {row[1]:+.4e}  {row[2]:+.4e}")
    print(sep)
    print(f"  Wall time    : {elapsed*1e3:.1f} ms")
    print(sep)

    # Per-atom force table for small structures or when DEBUG logging is active
    if natoms <= 20 or logging.getLogger("fmlip_relay").isEnabledFor(logging.DEBUG):
        print(f"\n  {'#':>4}  {'sym':>3}  {'Fx':>12}  {'Fy':>12}  {'Fz':>12}  |F|")
        for i, (sym, f) in enumerate(zip(symbols, forces)):
            fnorm = float(np.linalg.norm(f))
            print(
                f"  {i+1:>4}  {sym:>3}  "
                f"{f[0]:+12.6f}  {f[1]:+12.6f}  {f[2]:+12.6f}  {fnorm:.6f}"
            )
        print()
