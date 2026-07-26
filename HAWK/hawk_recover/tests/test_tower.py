# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Tower-descent test: synthetic key at ℓ=5–6 with the secret-derived V_τ
preset — exercises module_basis, the recursion, sub-level SVP, and the base
case in seconds."""

from __future__ import annotations

import random

import pytest

from hawk_recover import ring
from hawk_recover.params import HawkParams
from hawk_recover.tower import recover_basis
from hawk_recover.types import Intertwiner, PublicKey


def _synthetic_key(ell: int, seed: int):
    """Random B ∈ GL₂(O_L) (product of elementary shears), its Gram, and V_τ.

    All matrix algebra uses ring.gmat/mat_apply (generic implementation) — see
    the Matrix_cyclo_dense warning in ring.py."""
    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, tau, _ = ring.field_and_homs(ell)
    rng = random.Random(seed)

    def ro():
        return sum(rng.randint(-1, 1) * z**k for k in range(n))

    B = ring.gidentity(L, 2)
    for _ in range(3):
        B = B * ring.gmat(L, [[1, ro()], [0, 1]]) * ring.gmat(L, [[1, 0], [ro(), 1]])
    Q = ring.star(B, c) * B
    Vtau = B.inverse() * ring.mat_apply(B, tau)
    Q_list = tuple(
        [int(co) for co in list(Q[r, cc])] + [0] * (n - len(list(Q[r, cc])))
        for r in range(2)
        for cc in range(2)
    )
    pk = PublicKey(params=p, Q_list=Q_list, raw=b"")
    return pk, B, Vtau


@pytest.mark.parametrize("ell", [5, 6])
def test_tower_smoke_small(ell: int) -> None:
    pk, _, Vtau = _synthetic_key(ell, seed=100 + ell)
    p = pk.params
    L, _, c, tauH, _ = ring.field_and_homs(ell)
    assert Vtau * ring.mat_apply(Vtau, tauH) == ring.gidentity(L, 2)

    itw = Intertwiner(params=p, wt=None, V=Vtau, la_rank=p.n // 2)
    rb = recover_basis(itw, pk, verbose=True)
    rb.check()  # eq_G, in_OL, det_unit all hold

    # And explicitly: the recovered basis reproduces Q.
    Brec = rb.B(L)
    assert ring.star(Brec, c) * Brec == pk.Q(L)
