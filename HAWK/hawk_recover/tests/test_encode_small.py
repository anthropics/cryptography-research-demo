# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""The pre-packing head of encode() at small ℓ: det-normalization to a valid
HAWK basis and the key equation |f|²+|g|² = q00 that evidence/verify checks.

The final int8 packing of (f,g) is only defined for genuine HAWK keys (whose
short column is tiny by construction); a synthetic shear-product key has large
(f,g), so packing is exercised end-to-end by the real HAWK-256 run, not here.
"""

from __future__ import annotations

import random

import pytest

from hawk_recover import ring
from hawk_recover.params import HawkParams
from hawk_recover.tower import recover_basis
from hawk_recover.types import Intertwiner, PublicKey


def _instance(ell: int, seed: int):
    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, tau, _ = ring.field_and_homs(ell)
    rng = random.Random(seed)

    def small():
        return sum(rng.randint(-1, 1) * z**k for k in range(n))

    B = ring.gidentity(L, 2)
    for _ in range(3):
        B = B * ring.gmat(L, [[1, small()], [0, 1]]) * ring.gmat(L, [[1, 0], [small(), 1]])
    Q = ring.star(B, c) * B
    Q_list = tuple(
        [int(co) for co in list(Q[r, cc])] + [0] * (n - len(list(Q[r, cc])))
        for r in range(2)
        for cc in range(2)
    )
    pk = PublicKey(params=p, Q_list=Q_list, raw=b"\x00")  # raw only feeds hpub
    return pk, B, B.inverse() * ring.mat_apply(B, tau)


@pytest.mark.parametrize("ell", [5, 6])
def test_recovered_basis_normalizes_and_binds_to_pk(ell: int) -> None:
    from sage.all import ZZ

    from hawk_recover.encode import _coeffs

    pk, _, Vtau = _instance(ell, seed=700 + ell)
    p = pk.params
    n = p.n
    L, _z, c, _tau, _ = ring.field_and_homs(ell)
    itw = Intertwiner(params=p, wt=None, V=Vtau, la_rank=n // 2)
    rb = recover_basis(itw, pk)

    # encode's det-normalization step: scale row 0 by (det B)⁻¹ → det 1 and
    # B*·B = Q still holds (|det B|² = det Q = 1).
    B = rb.B(L)
    d = B.det()
    assert abs(ZZ((d * c(d)).norm())) == 1
    B = ring.gmat(L, [[d**-1, 0], [0, 1]]) * B
    assert B.det() == 1
    assert ring.star(B, c) * B == pk.Q(L)

    # Col 0 = (f,g) satisfies the key equation |f|²+|g|² = q00. (For a real
    # HAWK key col 0 is also the *short* column — q00 is the small diagonal —
    # but a synthetic shear-product Q has no small diagonal, so shortness and
    # the int8 packing are exercised only by the real HAWK-256 run.)
    f, g = _coeffs(B[0, 0], n), _coeffs(B[1, 0], n)
    ff = ring.nc_mul(f, ring.nc_conj(f, n), n)
    gg = ring.nc_mul(g, ring.nc_conj(g, n), n)
    assert [a + b for a, b in zip(ff, gg, strict=True)] == list(pk.Q_list[0])
