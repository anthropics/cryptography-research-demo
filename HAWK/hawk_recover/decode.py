# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
#
# Derived from decode_public_key() in hawk_vrfy.c of the HAWK development
# repository (https://github.com/hawk-sign/dev), MIT-licensed:
# Copyright (c) 2022 Hawk Project. See ../../NOTICE for the notice.

"""HAWK public-key decoding: bytes → Q = B*·B over O_L.

Exact port of `decode_public_key` (hawk_vrfy.c): segregated Golomb–Rice coding;
q00[0] carries `16 − lim_bits` extra low bits in the q00 block's trailing slack.
"""

from __future__ import annotations

import numpy as np

from hawk_recover import ring
from hawk_recover.params import HawkParams
from hawk_recover.types import PublicKey, require


class DecodeError(ValueError):
    pass


def _ntz8(x: int) -> int:
    """Trailing-zero count of the low 8 bits (ntz(0)=8)."""
    for k in range(8):
        if (x >> k) & 1:
            return k
    return 8


def _decode_gr_seg(
    buf: bytes, logn: int, low: int, lim_bits: int
) -> tuple[list[int], int, int]:
    """Segregated GR decode of n=2^logn values (hawk_vrfy.c:decode_gr).

    Layout: [n/8 sign bytes][low·n/8 low bytes][variable unary part].
    Returns (values, bytes_consumed, ignored_bits_in_last_byte)."""
    n = 1 << logn
    sign_bytes = n >> 3
    low_bytes = (low * n) >> 3
    voff = sign_bytes + low_bytes
    if len(buf) < voff:
        raise DecodeError(f"buffer too short for sign+low ({len(buf)} < {voff})")

    # variable (unary high) part
    lim_hi = 1 << (lim_bits - low)
    d = [0] * n
    acc, acc_off = 0, 0
    for u in range(n):
        while acc == 0:
            if acc_off >= lim_hi:
                raise DecodeError(f"high part at index {u} exceeds limit")
            if voff >= len(buf):
                raise DecodeError("buffer exhausted in variable part")
            acc |= buf[voff] << acc_off
            voff += 1
            acc_off += 8
        k = _ntz8(acc & 0xFF)
        if k == 8:
            k += _ntz8((acc >> 8) & 0xFF)
        if k >= lim_hi:
            raise DecodeError(f"high part {k} exceeds limit {lim_hi}")
        d[u] = k << low
        acc >>= k + 1
        acc_off -= k + 1
    ni = acc_off

    # sign + fixed-width low part
    lmask = (1 << low) - 1
    for u in range(0, n, 8):
        sbb = buf[u >> 3]
        base = sign_bytes + (u * low >> 3)
        lpp = int.from_bytes(buf[base : base + low], "little")
        for i in range(8):
            lp = (lpp >> (i * low)) & lmask
            v = d[u + i] | lp
            d[u + i] = ~v if (sbb >> i) & 1 else v  # ~v = -(v+1) as signed

    return d, voff, ni


def _decode_q00(
    buf: bytes, logn: int, low: int, lim_bits: int
) -> tuple[list[int], int]:
    """Decode n/2 values of q00, then patch q00[0] with its extra bits
    (hawk_vrfy.c:decode_q00)."""
    d, voff, ni = _decode_gr_seg(buf, logn - 1, low, lim_bits)
    eb_len = 16 - lim_bits
    if eb_len <= ni:
        last = buf[voff - 1] >> (8 - ni)
        eb = last & ((1 << eb_len) - 1)
        rest = last >> eb_len
    else:
        if voff >= len(buf):
            raise DecodeError("buffer exhausted reading q00[0] extra bits")
        eb = buf[voff - 1] >> (8 - ni)
        last = buf[voff]
        voff += 1
        eb |= (last << ni) & ((1 << eb_len) - 1)
        rest = last >> (eb_len - ni)
    if rest != 0:
        raise DecodeError("nonzero bits after q00[0] extra field")
    d[0] = (d[0] << eb_len) | eb
    return d, voff


def _nc_conj(a: list[int], n: int) -> list[int]:
    """Image under c: ζ ↦ ζ⁻¹."""
    return [a[0]] + [-a[n - k] for k in range(1, n)]


def _nc_divexact(num: list[int], den: list[int], n: int) -> list[int]:
    """Exact division num/den in Z[x]/(x^n+1). Solves the well-conditioned
    n×n negacyclic system in float64, rounds, and verifies exactly."""
    M = np.zeros((n, n), dtype=np.float64)
    for j in range(n):
        for i, di in enumerate(den):
            k = i + j
            if k < n:
                M[k, j] += di
            else:
                M[k - n, j] -= di
    q_f = np.linalg.solve(M, np.asarray(num, dtype=np.float64))
    q = [int(round(x)) for x in q_f]
    require(
        max(abs(x - y) for x, y in zip(q_f, q, strict=True)) < 1e-3,
        "float/exact division mismatch",
    )
    if ring.nc_mul(q, den, n) != num:
        raise DecodeError("q11 not integral (q00 ∤ 1+|q01|²)")
    return q


def decode(pk: bytes, params: HawkParams | None = None) -> PublicKey:
    """Decode a HAWK public key (spec-format; no header byte)."""
    p = params or HawkParams.from_pk_len(len(pk))
    n, logn = p.n, p.ell - 1
    require(
        len(pk) == p.pk_bytes, f"pk length {len(pk)} ≠ {p.pk_bytes} for ell={p.ell}"
    )

    q00_half, off = _decode_q00(pk, logn, p.k_gr_q00, p.lim_bits_q00)
    q01, off2, ni = _decode_gr_seg(pk[off:], logn, p.k_gr_q01, p.lim_bits_q01)
    if ni and (pk[off + off2 - 1] >> (8 - ni)) != 0:
        raise DecodeError("nonzero trailing bits in q01 block")
    if any(pk[off + off2 :]):
        raise DecodeError("nonzero padding bytes")

    # q00 is c-fixed: q00[0]; q00[1..n/2-1]; q00[n/2]=0; q00[n-k]=-q00[k].
    q00 = list(q00_half) + [0] + [-q00_half[n // 2 - k] for k in range(1, n // 2)]
    require(len(q00) == n, f"q00 has {len(q00)} coefficients, expected {n}")

    q10 = _nc_conj(q01, n)
    num = ring.nc_mul(q01, q10, n)
    num[0] += 1
    q11 = _nc_divexact(num, q00, n)

    return PublicKey(params=p, Q_list=(q00, q01, q10, q11), raw=pk)
