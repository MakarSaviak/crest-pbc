"""
tests/test_protocol.py
----------------------
Tests for fmlip_relay.protocol: constants, wire encoding/decoding,
recv_exactly / send_all edge-cases.

All I/O uses socket.socketpair() so no real network is required.
"""

import socket
import struct
import threading

import numpy as np
import pytest

from fmlip_relay.protocol import (
    MSG_COMPUTE, MSG_QUIT, MSG_PING,
    STATUS_OK, STATUS_ERROR,
    recv_exactly, send_all,
    ComputeRequest, read_compute_request,
    write_ok_response, write_error_response,
)


# ── helpers ───────────────────────────────────────────────────────────────────

def socketpair():
    """Return a connected (writer, reader) socket pair."""
    return socket.socketpair(socket.AF_UNIX, socket.SOCK_STREAM)


def _encode_compute_payload(natoms, atomic_numbers, positions, cell, pbc, compute_stress, charge, spin):
    """Encode a COMPUTE request body (without the leading msg_type int32)."""
    buf  = struct.pack('<i', natoms)
    buf += atomic_numbers.astype('<i4').tobytes()
    buf += positions.astype('<f8').tobytes()
    buf += cell.astype('<f8').tobytes()
    buf += struct.pack('<3i', *[int(p) for p in pbc])
    buf += struct.pack('<i',  int(compute_stress))
    buf += struct.pack('<i', int(charge))
    buf += struct.pack('<i', int(spin))
    return buf


# ── protocol constants ────────────────────────────────────────────────────────

class TestConstants:
    def test_message_type_values(self):
        assert MSG_COMPUTE == 1
        assert MSG_QUIT    == 2
        assert MSG_PING    == 3

    def test_status_values(self):
        assert STATUS_OK    == 0
        assert STATUS_ERROR == 1

    def test_values_are_distinct(self):
        assert len({MSG_COMPUTE, MSG_QUIT, MSG_PING}) == 3


# ── recv_exactly ──────────────────────────────────────────────────────────────

class TestRecvExactly:
    def test_receives_full_payload(self):
        w, r = socketpair()
        data = b"\x01\x02\x03\x04\x05\x06\x07\x08"
        w.sendall(data)
        w.close()
        assert recv_exactly(r, 8) == data
        r.close()

    def test_receives_in_multiple_chunks(self):
        """Send in two halves; recv_exactly must stitch them together."""
        w, r = socketpair()
        payload = bytes(range(16))

        def _send():
            w.sendall(payload[:8])
            w.sendall(payload[8:])
            w.close()

        t = threading.Thread(target=_send)
        t.start()
        result = recv_exactly(r, 16)
        t.join()
        r.close()
        assert result == payload

    def test_raises_on_eof(self):
        w, r = socketpair()
        w.close()
        with pytest.raises(ConnectionError):
            recv_exactly(r, 4)
        r.close()

    def test_zero_bytes(self):
        w, r = socketpair()
        w.close()
        assert recv_exactly(r, 0) == b""
        r.close()


# ── send_all ──────────────────────────────────────────────────────────────────

class TestSendAll:
    def test_sends_full_payload(self):
        w, r = socketpair()
        data = bytes(range(64))
        send_all(w, data)
        w.close()
        received = b""
        while chunk := r.recv(1024):
            received += chunk
        r.close()
        assert received == data

    def test_raises_on_broken_socket(self):
        w, r = socketpair()
        r.close()
        with pytest.raises((ConnectionError, BrokenPipeError, OSError)):
            # keep sending until the OS surfaces the broken pipe
            for _ in range(1000):
                send_all(w, b"\x00" * 65536)
        w.close()


# ── ComputeRequest ────────────────────────────────────────────────────────────

class TestComputeRequest:
    def _make(self, natoms=2):
        return ComputeRequest(
            natoms         = natoms,
            atomic_numbers = np.array([18]*natoms, dtype=np.int32),
            positions      = np.zeros((natoms, 3), dtype=np.float64),
            cell           = np.eye(3, dtype=np.float64),
            pbc            = np.array([False, False, False]),
            compute_stress = False,
            charge         = 0,
            spin           = 1,
        )

    def test_attributes_accessible(self):
        req = self._make(4)
        assert req.natoms == 4
        assert req.atomic_numbers.shape == (4,)
        assert req.positions.shape      == (4, 3)
        assert req.cell.shape           == (3, 3)
        assert req.pbc.shape            == (3,)
        assert req.compute_stress is False
        assert req.charge == 0
        assert req.spin   == 1


# ── read_compute_request round-trip ──────────────────────────────────────────

class TestReadComputeRequest:
    def _write_and_read(self, natoms, atomic_numbers, positions, cell, pbc, compute_stress,
                        charge=0, spin=1):
        w, r = socketpair()
        payload = _encode_compute_payload(
            natoms, atomic_numbers, positions, cell, pbc, compute_stress, charge, spin
        )
        w.sendall(payload)
        w.close()
        req = read_compute_request(r)
        r.close()
        return req

    def test_roundtrip_open_boundary(self):
        N = 3
        Z   = np.array([1, 6, 8],  dtype=np.int32)
        pos = np.random.default_rng(0).uniform(0, 5, (N, 3))
        cell = np.eye(3) * 10.0
        pbc  = np.array([False, False, False])

        req = self._write_and_read(N, Z, pos, cell, pbc, False, charge=0, spin=1)

        assert req.natoms == N
        np.testing.assert_array_equal(req.atomic_numbers, Z)
        np.testing.assert_allclose(req.positions, pos)
        np.testing.assert_allclose(req.cell,      cell)
        np.testing.assert_array_equal(req.pbc,    pbc)
        assert req.compute_stress is False
        assert req.charge == 0
        assert req.spin   == 1

    def test_roundtrip_periodic(self):
        N = 5
        Z   = np.full(N, 13, dtype=np.int32)
        pos = np.random.default_rng(1).uniform(0, 4, (N, 3))
        cell = np.diag([4.05, 4.05, 4.05])
        pbc  = np.array([True, True, True])

        req = self._write_and_read(N, Z, pos, cell, pbc, True, charge=-1, spin=2)

        assert req.compute_stress is True
        np.testing.assert_array_equal(req.pbc, pbc)
        assert req.charge == -1
        assert req.spin   == 2

    def test_invalid_natoms_raises(self):
        w, r = socketpair()
        w.sendall(struct.pack('<i', 0))   # natoms = 0
        w.close()
        with pytest.raises(ValueError, match="natoms"):
            read_compute_request(r)
        r.close()

    def test_negative_natoms_raises(self):
        w, r = socketpair()
        w.sendall(struct.pack('<i', -1))
        w.close()
        with pytest.raises(ValueError, match="natoms"):
            read_compute_request(r)
        r.close()


# ── response writers ──────────────────────────────────────────────────────────

class TestResponseWriters:
    def _read_response(self, sock, natoms):
        status,  = struct.unpack('<i', recv_exactly(sock, 4))
        energy,  = struct.unpack('<d', recv_exactly(sock, 8))
        forces   = np.frombuffer(recv_exactly(sock, natoms*3*8), dtype='<f8').reshape(natoms, 3)
        stress   = np.frombuffer(recv_exactly(sock, 9*8),        dtype='<f8').reshape(3, 3)
        return status, energy, forces, stress

    def test_ok_response_layout(self):
        N = 4
        w, r = socketpair()
        E = -3.14
        F = np.ones((N, 3), dtype=np.float64) * 0.5
        S = np.eye(3, dtype=np.float64) * 0.01

        write_ok_response(w, E, F, S)
        w.close()

        status, energy, forces, stress = self._read_response(r, N)
        r.close()

        assert status == STATUS_OK
        assert energy == pytest.approx(E)
        np.testing.assert_allclose(forces, F)
        np.testing.assert_allclose(stress, S)

    def test_error_response_layout(self):
        N = 3
        w, r = socketpair()
        write_error_response(w, N)
        w.close()

        status, energy, forces, stress = self._read_response(r, N)
        r.close()

        assert status == STATUS_ERROR
        assert energy == pytest.approx(0.0)
        np.testing.assert_array_equal(forces, np.zeros((N, 3)))
        np.testing.assert_array_equal(stress, np.zeros((3, 3)))

    def test_ok_response_endianness(self):
        """Verify the status field is a 4-byte little-endian int."""
        w, r = socketpair()
        write_ok_response(w, 0.0, np.zeros((1, 3)), np.zeros((3, 3)))
        w.close()
        raw = recv_exactly(r, 4)
        r.close()
        assert struct.unpack('<i', raw)[0] == STATUS_OK
