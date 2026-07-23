"""
tests/test_server.py
--------------------
Integration tests for fmlip_relay.server.

Spins up a real TCP server (on a random loopback port) in a background
thread, then drives it with a raw socket client.  Tests cover:
  - PING round-trip
  - COMPUTE with DummyBackend
  - QUIT shuts the server down cleanly
  - Backend exception → STATUS_ERROR response (server stays alive)
  - Unknown message type drops the connection (server stays alive)
  - Multiple sequential COMPUTE requests on one connection
"""

import socket
import struct
import threading
import time

import numpy as np
import pytest

from fmlip_relay.backends.dummy import DummyBackend
from fmlip_relay.backends.lj    import LJBackend
from fmlip_relay.protocol import (
    MSG_COMPUTE, MSG_QUIT, MSG_PING,
    STATUS_OK, STATUS_ERROR,
    recv_exactly,
)
from fmlip_relay import server


# ── helpers ───────────────────────────────────────────────────────────────────

def _free_port() -> int:
    """Ask the OS for a free loopback port."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _encode_compute(natoms, atomic_numbers, positions, cell, pbc, compute_stress,
                    charge=0, spin=1):
    buf  = struct.pack('<i', MSG_COMPUTE)
    buf += struct.pack('<i', natoms)
    buf += atomic_numbers.astype('<i4').tobytes()
    buf += positions.astype('<f8').tobytes()
    buf += cell.astype('<f8').tobytes()
    buf += struct.pack('<3i', *[int(p) for p in pbc])
    buf += struct.pack('<i',  int(compute_stress))
    buf += struct.pack('<i',  int(charge))
    buf += struct.pack('<i',  int(spin))
    return buf


def _read_response(conn, natoms):
    status, = struct.unpack('<i', recv_exactly(conn, 4))
    energy, = struct.unpack('<d', recv_exactly(conn, 8))
    forces  = np.frombuffer(recv_exactly(conn, natoms*3*8), dtype='<f8').reshape(natoms, 3)
    stress  = np.frombuffer(recv_exactly(conn, 9*8),        dtype='<f8').reshape(3, 3)
    return status, energy, forces, stress


def _simple_config(natoms=4):
    """Return (Z, pos, cell, pbc) for a minimal non-interacting test config."""
    Z    = np.ones(natoms, dtype=np.int32)
    pos  = np.zeros((natoms, 3), dtype=np.float64)
    cell = np.eye(3, dtype=np.float64) * 100.0
    pbc  = np.array([False, False, False])
    return Z, pos, cell, pbc


# ── fixture: running server ───────────────────────────────────────────────────

class _ServerHandle:
    """Wraps a server thread and a connected client socket."""
    def __init__(self, port, backend):
        self.port    = port
        self.backend = backend
        self._thread = None

    def start(self):
        self._thread = threading.Thread(
            target=server.run,
            args=(self.port, self.backend),
            daemon=True,
        )
        self._thread.start()
        # Wait until the port is actually listening
        for _ in range(50):
            try:
                s = socket.socket()
                s.connect(("127.0.0.1", self.port))
                s.close()
                break
            except ConnectionRefusedError:
                time.sleep(0.05)
        else:
            raise RuntimeError("Server did not start in time")

    def connect(self) -> socket.socket:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect(("127.0.0.1", self.port))
        s.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        return s

    def join(self, timeout=2.0):
        if self._thread:
            self._thread.join(timeout=timeout)


@pytest.fixture
def dummy_server():
    port = _free_port()
    handle = _ServerHandle(port, DummyBackend(seed=0))
    handle.start()
    yield handle
    # Teardown: send QUIT if the server is still running
    try:
        c = handle.connect()
        c.sendall(struct.pack('<i', MSG_QUIT))
        c.close()
    except OSError:
        pass
    handle.join()


@pytest.fixture
def lj_server():
    port = _free_port()
    handle = _ServerHandle(port, LJBackend(epsilon=1.0, sigma=1.0, r_cut=4.0))
    handle.start()
    yield handle
    try:
        c = handle.connect()
        c.sendall(struct.pack('<i', MSG_QUIT))
        c.close()
    except OSError:
        pass
    handle.join()


# ── PING ──────────────────────────────────────────────────────────────────────

class TestPing:
    def test_ping_returns_ok(self, dummy_server):
        c = dummy_server.connect()
        c.sendall(struct.pack('<i', MSG_PING))
        status, = struct.unpack('<i', recv_exactly(c, 4))
        c.close()
        assert status == STATUS_OK

    def test_multiple_pings_on_same_connection(self, dummy_server):
        c = dummy_server.connect()
        for _ in range(5):
            c.sendall(struct.pack('<i', MSG_PING))
            status, = struct.unpack('<i', recv_exactly(c, 4))
            assert status == STATUS_OK
        c.close()


# ── COMPUTE ───────────────────────────────────────────────────────────────────

class TestCompute:
    def test_compute_returns_ok_status(self, dummy_server):
        N = 4
        Z, pos, cell, pbc = _simple_config(N)
        c = dummy_server.connect()
        c.sendall(_encode_compute(N, Z, pos, cell, pbc, False))
        status, energy, forces, stress = _read_response(c, N)
        c.close()
        assert status == STATUS_OK

    def test_compute_returns_correct_shapes(self, dummy_server):
        N = 6
        Z, pos, cell, pbc = _simple_config(N)
        c = dummy_server.connect()
        c.sendall(_encode_compute(N, Z, pos, cell, pbc, True))
        status, energy, forces, stress = _read_response(c, N)
        c.close()
        assert isinstance(energy, float)
        assert forces.shape == (N, 3)
        assert stress.shape == (3, 3)

    def test_compute_lj_newton_third_law(self, lj_server):
        """Forces from the server must sum to zero."""
        N = 4
        Z   = np.full(N, 18, dtype=np.int32)
        pos = np.array([[0,0,0],[3,0,0],[0,3,0],[3,3,0]], dtype=np.float64)
        cell = np.eye(3) * 100.0
        pbc  = np.array([False, False, False])
        c = lj_server.connect()
        c.sendall(_encode_compute(N, Z, pos, cell, pbc, False))
        status, _, forces, _ = _read_response(c, N)
        c.close()
        assert status == STATUS_OK
        np.testing.assert_allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)

    def test_multiple_sequential_computes(self, dummy_server):
        """Server must handle multiple COMPUTE requests on one connection."""
        N = 2
        Z, pos, cell, pbc = _simple_config(N)
        c = dummy_server.connect()
        results = []
        for _ in range(3):
            c.sendall(_encode_compute(N, Z, pos, cell, pbc, False))
            status, energy, forces, stress = _read_response(c, N)
            results.append((status, energy))
        c.close()
        assert all(s == STATUS_OK for s, _ in results)

    def test_compute_then_ping_on_same_connection(self, dummy_server):
        N = 2
        Z, pos, cell, pbc = _simple_config(N)
        c = dummy_server.connect()

        c.sendall(_encode_compute(N, Z, pos, cell, pbc, False))
        status, _, _, _ = _read_response(c, N)
        assert status == STATUS_OK

        c.sendall(struct.pack('<i', MSG_PING))
        ping_status, = struct.unpack('<i', recv_exactly(c, 4))
        assert ping_status == STATUS_OK
        c.close()

    def test_charge_and_spin_accepted(self, dummy_server):
        """Non-default charge and spin must not cause an error response."""
        N = 2
        Z, pos, cell, pbc = _simple_config(N)
        c = dummy_server.connect()
        c.sendall(_encode_compute(N, Z, pos, cell, pbc, False, charge=-1, spin=3))
        status, _, _, _ = _read_response(c, N)
        c.close()
        assert status == STATUS_OK


# ── QUIT ──────────────────────────────────────────────────────────────────────

class TestQuit:
    def test_quit_stops_server(self, dummy_server):
        c = dummy_server.connect()
        c.sendall(struct.pack('<i', MSG_QUIT))
        c.close()
        dummy_server.join(timeout=3.0)
        assert not dummy_server._thread.is_alive(), \
            "Server thread did not stop after QUIT"

    def test_server_unreachable_after_quit(self, dummy_server):
        c = dummy_server.connect()
        c.sendall(struct.pack('<i', MSG_QUIT))
        c.close()
        dummy_server.join(timeout=3.0)

        # New connection should be refused
        time.sleep(0.1)
        with pytest.raises((ConnectionRefusedError, OSError)):
            s = socket.socket()
            s.settimeout(1.0)
            s.connect(("127.0.0.1", dummy_server.port))
            s.close()


# ── Error handling ────────────────────────────────────────────────────────────

class TestErrorHandling:
    def test_backend_exception_returns_error_status(self):
        """A backend that always raises must yield STATUS_ERROR, server stays alive."""
        class BrokenBackend(DummyBackend):
            def compute(self, *a, **kw):
                raise RuntimeError("intentional failure")

        port = _free_port()
        handle = _ServerHandle(port, BrokenBackend())
        handle.start()

        N = 2
        Z, pos, cell, pbc = _simple_config(N)
        c = handle.connect()
        c.sendall(_encode_compute(N, Z, pos, cell, pbc, False))
        status, energy, forces, stress = _read_response(c, N)
        c.close()

        assert status == STATUS_ERROR
        assert energy == pytest.approx(0.0)
        np.testing.assert_array_equal(forces, np.zeros((N, 3)))

        # Server should still be responsive
        c2 = handle.connect()
        c2.sendall(struct.pack('<i', MSG_PING))
        ping_status, = struct.unpack('<i', recv_exactly(c2, 4))
        c2.sendall(struct.pack('<i', MSG_QUIT))
        c2.close()
        assert ping_status == STATUS_OK
        handle.join()

    def test_unknown_message_drops_connection(self, dummy_server):
        """An unknown message type must drop the client without killing the server."""
        c = dummy_server.connect()
        c.sendall(struct.pack('<i', 255))   # unknown type
        c.close()

        # Server should still accept a new connection
        time.sleep(0.1)
        c2 = dummy_server.connect()
        c2.sendall(struct.pack('<i', MSG_PING))
        status, = struct.unpack('<i', recv_exactly(c2, 4))
        c2.close()
        assert status == STATUS_OK

    def test_client_disconnect_mid_request_does_not_crash_server(self, dummy_server):
        """Abrupt client disconnect during a COMPUTE must not kill the server."""
        c = dummy_server.connect()
        # Send only the msg_type and natoms – then vanish
        c.sendall(struct.pack('<i', MSG_COMPUTE))
        c.sendall(struct.pack('<i', 4))
        c.close()

        time.sleep(0.2)

        # Server should still be alive
        c2 = dummy_server.connect()
        c2.sendall(struct.pack('<i', MSG_PING))
        status, = struct.unpack('<i', recv_exactly(c2, 4))
        c2.close()
        assert status == STATUS_OK
