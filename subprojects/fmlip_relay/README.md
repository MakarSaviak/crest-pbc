# fmlip-relay

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Persistent ML potential server for Fortran clients.**  
Zero model-loading overhead on repeated calls. Supports multiple parallel instances via OpenMP or MPI.

---

## Architecture

```
┌──────────────────────────────────────────────────────────────┐
│  Fortran Driver (MD, MC, geometry optimiser, …)              │
│  use fmlip_relay_client                                      │
│    mlip_init      → spawn fmlip-relay-server, wait READY     │
│    mlip_compute   → send coords, recv energy/forces/stress   │
│    mlip_ping      → liveness check                           │
│    mlip_finalize  → send QUIT, kill process                  │
└────────────────────────────┬─────────────────────────────────┘
                             │  127.0.0.1:PORT  (binary TCP)
                  ┌──────────┴──────────┐
                  │  fmlip-relay-server │
                  │  server.py          │  ← protocol + event loop
                  │  backends/*.py      │  ← model loaded once at startup
                  └─────────────────────┘
```

The Fortran side spawns one Python subprocess per instance, each on its own port. Because the model is loaded once at startup, repeated `mlip_compute` calls incur only socket round-trip latency (~1 ms on loopback).

---

## Installation

### Python package

```bash
# Core (dummy + LJ backends, no ML deps)
pip install .

# With MACE support
pip install ".[mace]"

# With UMA support (Meta FAIR / fairchem)
pip install ".[uma]"
# UMA checkpoints are gated on the Hugging Face Hub — authenticate once:
huggingface-cli login

# Development install (includes pytest)
pip install -e ".[dev]"
```

> **Note:** The `[mace]` and `[uma]` extras pull in the upstream packages
> (`mace-torch` / `fairchem-core`) with default settings, which may not match
> your platform (CUDA version, PyTorch build, etc.). For a working ML setup,
> especially GPU support, follow the official installation instructions on the
> respective project repositories:
> [MACE](https://github.com/ACEsuit/mace) and
> [fairchem/UMA](https://github.com/facebookresearch/fairchem).

### Fortran/C library (CMake)

```bash
cmake -S . -B build -DWITH_OpenMP=ON -DBUILD_EXAMPLE=ON
cmake --build build
cmake --install build --prefix ~/.local
```

Consuming in another CMake project after installation:

```cmake
find_package(fmlip_relay REQUIRED)
target_link_libraries(my_code PRIVATE fmlip_relay::fmlip_relay)
```

---

## Quick Start

### 1. Start a server manually (optional — Fortran does this automatically)

```bash
# Dummy backend
fmlip-relay-server --port 54321 --backend dummy

# Lennard-Jones backend (Argon defaults)
fmlip-relay-server --port 54321 --backend lj
fmlip-relay-server --port 54321 --backend lj --lj-epsilon 0.0104 --lj-sigma 3.40 --lj-cutoff 8.50

# MACE backend
fmlip-relay-server --port 54321 --backend mace --model /path/to/model.model \
    --device cpu --dtype float64

# UMA backend (Meta FAIR / fairchem) — choose a task head
fmlip-relay-server --port 54321 --backend uma --uma-task omol --device cuda
fmlip-relay-server --port 54321 --backend uma --uma-model uma-m-1 --uma-task omat
```

### 2. Build and run the Fortran example

```bash
# CMake
cmake -S . -B build -DBUILD_EXAMPLE=ON && cmake --build build
./build/run_example dummy
./build/run_example lj
./build/run_example lj 0.0104 3.40 8.50
./build/run_example mace /path/to/model.model

# Or with the plain Makefile
cd example && make
./build/run_example dummy
```

### 3. Minimal Fortran usage

```fortran
use fmlip_relay_client

integer :: ierr
integer :: pbc(3), atomic_numbers(64)
real(8) :: pos(3,64), cell(3,3), energy, forces(3,64), stress(3,3)
integer :: charge, spin

! Spawn server — blocks until Python prints "READY"
call mlip_init(1, 54321, &
    "fmlip-relay-server --port 54321 --backend mace --model /path/model.model", &
    60, ierr)

! Evaluate
pbc = [1, 1, 1]
atomic_numbers = 13   ! aluminium
charge = 0 
spin = 1
call mlip_compute(1, 64, atomic_numbers, pos, cell, pbc, 1, charge, spin, &
                  energy, forces, stress, ierr)

! Teardown
call mlip_finalize(1, ierr)
```

For molecular clusters (no PBC), pass `pbc = [0, 0, 0]` and `compute_stress = 0`. The cell is unused in that case but should be a non-degenerate matrix to avoid accidental division by zero.

---

## Backends

| Backend | CLI name | Extra deps | Notes |
|---------|----------|------------|-------|
| Dummy   | `dummy`  | —          | Random numbers; validates plumbing without a model |
| Lennard-Jones | `lj` | — | Force-shifted cutoff; PBC via minimum-image; Argon defaults |
| MACE    | `mace`   | `mace-torch`, `ase` | Custom-trained `.model` file; requires `pip install ".[mace]"` |
| MACE-MP | `mace_mp`| `mace-torch`, `ase` | Materials Project foundation model (89 elements) |
| MACE-OFF| `mace_off`| `mace-torch`, `ase` | Organic foundation force field (10 elements) |
| UMA     | `uma`    | `fairchem-core`, `ase` | Meta FAIR universal model; multi-task (`--uma-task omol\|omat\|omc\|oc20\|odac`); gated HF checkpoints; `pip install ".[uma]"` |

### UMA parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--uma-model` | `uma-s-1p2` | Checkpoint: `uma-s-1`, `uma-s-1p1`, `uma-s-1p2`, `uma-m-1` |
| `--uma-task`  | `omol`     | Domain head: `omol`, `omat`, `omc`, `oc20`, `odac` |
| `--device`    | `cpu`      | PyTorch device; UMA runs best on `cuda` |

The `omol` task uses the total charge and spin multiplicity forwarded through
the relay protocol (`atoms.info["charge"]` / `atoms.info["spin"]`).

### Lennard-Jones parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--lj-epsilon` | `0.0104` eV | Well depth ε (Argon-like) |
| `--lj-sigma`   | `3.40` Å   | Zero-crossing distance σ |
| `--lj-cutoff`  | `8.50` Å   | Cutoff radius (≈ 2.5 σ) |

---

## Wire Protocol

All values are little-endian.

### Request (Fortran → Python)

```
Offset      Type           Field
────────────────────────────────────────────────────────────
 0          int32          msg_type   (1=COMPUTE, 2=QUIT, 3=PING)
 4          int32          natoms
 8          int32  × N     atomic_numbers  (Z per atom)
 8+4N       float64 × 3N   positions       [atom-major, Angstrom]
 8+28N      float64 × 9    cell            [row-major,  Angstrom]
 8+28N+72   int32  × 3     pbc             (1=periodic per axis)
 8+28N+84   int32          compute_stress  (0/1)
 8+28N+88   int32          charge
 8+28N+92   int32          spin 
```

### Response (Python → Fortran)

```
Offset   Type           Field
─────────────────────────────────────────────
 0       int32          status    (0=OK, 1=ERROR)
 4       float64        energy    [eV]
12       float64 × 3N   forces    [eV/Å, atom-major]
12+24N   float64 × 9    stress    [eV/Å³, row-major]
```

---

## Parallel Instances (OpenMP)

Each instance owns an independent socket fd, so OpenMP threads can call different instances simultaneously without locking.

```fortran
! Start all servers serially before the parallel region
! (avoids races when loading from a shared filesystem)
do i = 1, N
  write(cmd, '("fmlip-relay-server --port ",I0," --backend mace --model ",A)') &
        BASE_PORT+i, trim(model_path)
  call mlip_init(i, BASE_PORT+i, trim(cmd), 120, ierr)
end do

!$omp parallel do num_threads(N) private(iid, ...) schedule(static,1)
do i = 1, N
  iid = omp_get_thread_num() + 1
  call mlip_compute(iid, natoms, atomic_numbers, pos, cell, pbc, 0, E, F, S, ierr)
end do
!$omp end parallel do

call mlip_finalize_all(ierr)
```

For MPI, give each rank its own port range: `port = base_port + mpi_rank * 100 + local_instance`.

---

## Adding a New Backend

1. Create `src/python/fmlip_relay/backends/mybackend.py`:

```python
from .base import BackendBase
import numpy as np

class MyBackend(BackendBase):
    def __init__(self, **kwargs):
        pass  # load model here — called once at server startup

    def compute(self, atomic_numbers, positions, cell, pbc, compute_stress):
        # atomic_numbers : (N,)   int32
        # positions      : (N, 3) float64, Angstrom
        # cell           : (3, 3) float64, Angstrom
        # pbc            : (3,)   bool
        # returns        : (energy: float, forces: (N,3), stress: (3,3))
        ...
```

2. Register it in `src/python/fmlip_relay/backends/__init__.py`:

```python
from .mybackend import MyBackend
_REGISTRY["mybackend"] = MyBackend
```

3. Add a CLI argument group in `src/python/fmlip_relay/__main__.py` if the backend needs parameters, and handle construction in `_build_backend`.

---

## Testing

```bash
pip install -e ".[dev]"
pytest
```

Tests are organised into three modules:

- **`test_protocol.py`** — wire encoding/decoding, `recv_exactly`/`send_all` edge cases, request round-trips, response layout
- **`test_backends.py`** — `DummyBackend` output contracts; `LJBackend` constructor validation, analytical energy, force-shift continuity at cutoff, Newton's 3rd law, finite-difference force consistency, stress symmetry, PBC vs open-boundary behaviour
- **`test_server.py`** — full TCP integration: PING, sequential COMPUTEs, QUIT, backend exceptions → `STATUS_ERROR`, unknown message type, abrupt client disconnect

---

## Troubleshooting

**`mlip_init` times out waiting for READY**  
Verify `fmlip-relay-server` is on `PATH` (`which fmlip-relay-server`). Pass `--loglevel DEBUG --log /tmp/fb.log` in the server command to capture Python-side startup errors.

**Connection refused after READY**  
The Fortran side retries the connect up to 10 times with 100 ms gaps. If the model is large and startup is slow, increase `timeout_sec` in the `mlip_init` call.

**Wrong energies / garbled data**  
Confirm `positions` is declared `(3, natoms)` in Fortran (xyz-major). Both sides assume little-endian, which is correct on all x86/ARM hardware.

**Stress is NaN for a cluster**  
Pass `compute_stress = 0` for non-periodic systems. The stress is formally undefined without a well-defined volume.

---

## Project Structure

```
fmlip_relay/
├── CMakeLists.txt               # CMake build (library + example)
├── pyproject.toml               # Python package (pip install)
├── pytest.ini
├── README.md
├── cmake/
│   └── fmlip_relay-config.cmake.in   # find_package() support
├── config/                      # CMake build configuration
│   ├── CMakeLists.txt
│   └── modules/
│       └── fmlip_relay-utils.cmake
├── src/
│   ├── c/
│   │   └── socket_utils.c       # POSIX socket/process wrappers
│   ├── fortran/
│   │   └── mlip_client.f90      # Fortran client module
│   └── python/
│       └── fmlip_relay/
│           ├── __init__.py
│           ├── __main__.py      # CLI entry point (fmlip-relay-server)
│           ├── protocol.py      # Wire constants + pack/unpack helpers
│           ├── server.py        # Backend-agnostic event loop
│           └── backends/
│               ├── base.py      # Abstract base class
│               ├── dummy.py     # Random backend (no deps, for testing)
│               ├── lj.py        # Lennard-Jones with PBC + force-shifted cutoff
│               ├── mace.py      # MACE backend (mace-torch + ase)
│               └── ...          # additional implementations at demand
├── example/
│   ├── CMakeLists.txt
│   ├── Makefile                 # Simple alternative to CMake
│   └── example_usage.f90
└── tests/
    ├── conftest.py
    ├── test_protocol.py         # Wire encoding/decoding
    ├── test_backends.py         # Backend contracts + LJ physics
    └── test_server.py           # Full TCP integration tests
```

---

## License

[MIT](LICENSE) © Philipp Pracht

This project is open-source and free to use for both personal and commercial projects. For more details, see the [LICENSE](LICENSE) file.
