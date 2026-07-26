# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""decode() on the shipped HAWK-256 public keys: the bit-exact Golomb–Rice
port must produce a Hermitian, positive-definite, det-1 Gram matrix Q."""

from __future__ import annotations

from pathlib import Path

import pytest

from hawk_recover import ring
from hawk_recover.decode import DecodeError, decode
from hawk_recover.params import HAWK256

PKS = sorted(Path(__file__).resolve().parents[2].glob("test_pk_*.hex"))


@pytest.mark.parametrize("pk_path", PKS, ids=lambda p: p.name)
def test_decode_shipped_pk(pk_path: Path) -> None:
    assert PKS, "test_pk_*.hex fixtures missing"
    raw = bytes.fromhex(pk_path.read_text().strip())
    pk = decode(raw)
    assert pk.params == HAWK256
    n = pk.params.n

    q00, q01, q10, q11 = pk.Q_list
    # q00 is c-fixed and q00[n/2] = 0 (structure decode reconstructs).
    assert q00 == ring.nc_conj(q00, n)
    assert q00[n // 2] == 0
    # q10 = c(q01); det Q = q00·q11 − q01·q10 = 1.
    assert q10 == ring.nc_conj(q01, n)
    detQ = ring.nc_mul(q00, q11, n)
    q01q10 = ring.nc_mul(q01, q10, n)
    detQ = [a - b for a, b in zip(detQ, q01q10, strict=True)]
    assert detQ == [1] + [0] * (n - 1)

    # Q(L) over the cyclotomic field is Hermitian, and q00 = |f|²+|g|² is
    # totally positive (positive under every complex embedding of K^+).
    L, _z, c, _tau, _ = ring.field_and_homs(pk.params.ell)
    Q = pk.Q(L)
    assert ring.star(Q, c) == Q
    assert all(e.real() > 0 for e in Q[0, 0].complex_embeddings())


def test_decode_rejects_wrong_length() -> None:
    from hawk_recover.types import VerificationError

    with pytest.raises((ValueError, VerificationError)):
        decode(b"\x00" * 449)


def test_decode_rejects_corrupt_pk() -> None:
    raw = bytearray(bytes.fromhex(PKS[0].read_text().strip()))
    # Trailing padding must be zero; setting a padding byte trips DecodeError.
    raw[-1] = 0xFF
    with pytest.raises(DecodeError):
        decode(bytes(raw))
