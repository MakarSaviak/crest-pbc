"""
fmlip_relay.protocol
-------------------
Wire protocol constants and low-level I/O helpers shared by the server
and any pure-Python test clients.

All values are little-endian.

Request layout (Fortran → Python)
──────────────────────────────────
  int32          msg_type
  int32          natoms
  int32  × N     atomic_numbers
  float64 × 3N   positions   [atom-major, Angstrom]
  float64 × 9    cell        [row-major,  Angstrom]
  int32  × 3     pbc         (1 = periodic)
  int32          compute_stress (0/1)
  int32          total charge          (neutral default 0)
  int32          spin multiplicity 2S+1 (ground-state singlet default 1)

Response layout (Python → Fortran)
────────────────────────────────────
  int32          status
  float64        energy      [eV]
  float64 × 3N   forces      [eV/Å, atom-major]
  float64 × 9    stress      [eV/Å³, row-major]
"""

import struct
import socket
import numpy as np

# ── message types ─────────────────────────────────────────────────────────────
MSG_COMPUTE = 1
MSG_QUIT    = 2
MSG_PING    = 3

# ── status codes ──────────────────────────────────────────────────────────────
STATUS_OK    = 0
STATUS_ERROR = 1

# ── I/O primitives ────────────────────────────────────────────────────────────

def recv_exactly(conn: socket.socket, n: int) -> bytes:
    """Block until exactly *n* bytes have been read from *conn*."""
    buf = bytearray()
    while len(buf) < n:
        chunk = conn.recv(n - len(buf))
        if not chunk:
            raise ConnectionError("Client disconnected unexpectedly")
        buf.extend(chunk)
    return bytes(buf)


def send_all(conn: socket.socket, data: bytes) -> None:
    """Send all of *data*, retrying on short writes."""
    total = 0
    mv = memoryview(data)
    while total < len(data):
        sent = conn.send(mv[total:])
        if sent == 0:
            raise ConnectionError("Socket broken during send")
        total += sent


# ── request reader ────────────────────────────────────────────────────────────

class ComputeRequest:
    """Decoded COMPUTE request."""
    __slots__ = ("natoms", "atomic_numbers", "positions", "cell", "pbc", "compute_stress", "charge", "spin")

    def __init__(self, natoms, atomic_numbers, positions, cell, pbc, compute_stress, charge, spin):
        self.natoms         = natoms
        self.atomic_numbers = atomic_numbers   # (N,)   int32
        self.positions      = positions        # (N, 3) float64, Angstrom
        self.cell           = cell             # (3, 3) float64, Angstrom
        self.pbc            = pbc              # (3,)   bool
        self.compute_stress = compute_stress   # bool
        self.charge         = charge           # int32
        self.spin           = spin             # int32


def read_compute_request(conn: socket.socket) -> ComputeRequest:
    """Read and decode one COMPUTE payload (msg_type already consumed)."""
    natoms, = struct.unpack('<i', recv_exactly(conn, 4))
    if natoms <= 0:
        raise ValueError(f"Invalid natoms={natoms}")

    atomic_numbers = np.frombuffer(recv_exactly(conn, natoms * 4),  dtype='<i4').copy()
    positions      = np.frombuffer(recv_exactly(conn, natoms * 3*8),dtype='<f8').reshape(natoms, 3).copy()
    cell           = np.frombuffer(recv_exactly(conn, 9 * 8),       dtype='<f8').reshape(3, 3).copy()
    pbc_raw        = np.array(struct.unpack('<3i', recv_exactly(conn, 12)), dtype=bool)
    compute_stress = bool(struct.unpack('<i', recv_exactly(conn, 4))[0])
    charge = struct.unpack('<i', recv_exactly(conn, 4))[0]
    spin   = struct.unpack('<i', recv_exactly(conn, 4))[0]

    return ComputeRequest(natoms, atomic_numbers, positions, cell, pbc_raw, compute_stress, charge, spin)


# ── response writer ───────────────────────────────────────────────────────────

def write_ok_response(conn: socket.socket,
                      energy: float,
                      forces: np.ndarray,
                      stress: np.ndarray) -> None:
    resp  = struct.pack('<i', STATUS_OK)
    resp += struct.pack('<d', energy)
    resp += forces.astype('<f8').tobytes()
    resp += stress.astype('<f8').tobytes()
    send_all(conn, resp)


def write_error_response(conn: socket.socket, natoms: int) -> None:
    resp  = struct.pack('<i', STATUS_ERROR)
    resp += struct.pack('<d', 0.0)
    resp += np.zeros((natoms, 3), dtype='<f8').tobytes()
    resp += np.zeros((3, 3),      dtype='<f8').tobytes()
    send_all(conn, resp)
