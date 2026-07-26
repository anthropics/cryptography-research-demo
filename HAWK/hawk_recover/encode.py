# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Encode a recovered basis B_rec to HAWK sk bytes (spec §3.4).

B_rec is any GL₂(O_L) matrix with B_rec*·B_rec = Q. HAWK's sk convention is
B = [[f,F],[g,G]] with (f,g) the short column, which col 0 of any such B_rec
already is (q00 = |f|²+|g|² is the small diagonal of Q). We scale row 0 by
(det B_rec)⁻¹ to force det 1, apply a c-unitary left variant U ∈ μ₄²⋊S₂ that
keeps det(U·B) = 1, and pack f,g as int8 and F,G mod 2 (spec §3.4).
"""

from __future__ import annotations

import struct
from typing import Any

from hawk_recover.infra import shim
from hawk_recover.types import PublicKey, RecoveredBasis, SecretKey, require


def _coeffs(x, n: int) -> list[int]:
    from sage.all import ZZ

    co = list(x)
    return [int(ZZ(co[k])) if k < len(co) else 0 for k in range(n)]


def _maxabs(*vs: list[int]) -> int:
    return max(abs(x) for v in vs for x in v)


def _pack_mod2(v: list[int]) -> bytes:
    n = len(v)
    require(n % 8 == 0, f"mod-2 pack needs n divisible by 8, got {n}")
    return bytes(
        sum(((v[8 * i + j] & 1) << j) for j in range(8)) for i in range(n // 8)
    )


def encode(rb: RecoveredBasis, pk: PublicKey) -> SecretKey:
    from hawk_recover import ring

    p = rb.params
    n = p.n
    L, z, c, _, _ = ring.field_and_homs(p.ell)
    I2 = ring.gidentity(L, 2)

    B = rb.B(L)
    Q = pk.Q(L)

    # Normalize det to 1 by scaling row 0: left-multiply by diag(d⁻¹, 1),
    # which preserves B*·B = Q since |d|² = det Q = 1.
    d = B.det()
    require(abs(int((d * c(d)).norm())) == 1, "det B_rec not a unit: N(d·c(d))≠1")
    B = ring.gmat(L, [[d**-1, 0], [0, 1]]) * B
    require(B.det() == 1, "det ≠ 1 after normalization")
    require(ring.star(B, c) * B == Q, "B*·B ≠ Q after det-normalization")

    # Search c-unitary U ∈ μ₄²⋊S₂ for det(U·B)=1 with minimal max|col0|.
    # Safety net only: since q00 is the small diagonal of Q = B*·B, col0 of any
    # such B is already the short column, and ζ^a rotation and row-swap both
    # preserve max|col0|; the loop reduces to "take the first det-1 variant".
    best: tuple[int, Any] | None = None
    powers = (0, n // 2, n, 3 * n // 2)
    for swap in (False, True):
        for a in powers:
            for b in powers:
                U = ring.gmat(L, [[z**a, 0], [0, z**b]])
                if swap:
                    U = ring.gmat(L, [[0, 1], [1, 0]]) * U
                if ring.star(U, c) * U != I2:
                    continue
                Bt = U * B
                if Bt.det() != 1:
                    continue
                f, g = _coeffs(Bt[0, 0], n), _coeffs(Bt[1, 0], n)
                mx = _maxabs(f, g)
                if best is None or mx < best[0]:
                    best = (mx, Bt)
    require(best is not None, "no c-unitary variant found (unreachable)")
    _, Bf = best
    require(ring.star(Bf, c) * Bf == Q, "B*·B ≠ Q for the selected variant")

    f = _coeffs(Bf[0, 0], n)
    g = _coeffs(Bf[1, 0], n)
    F = _coeffs(Bf[0, 1], n)
    G = _coeffs(Bf[1, 1], n)
    require(_maxabs(f, g) < 128, f"(f,g) exceed int8 range: max={_maxabs(f, g)}")

    core = (
        struct.pack(f"{n}b", *f)
        + struct.pack(f"{n}b", *g)
        + _pack_mod2(F)
        + _pack_mod2(G)
    )
    hp = shim.hpub(pk.raw, p)
    sk = SecretKey(params=p, f=f, g=g, F=F, G=G, core=core, hpub=hp)
    sk.check()
    return sk
