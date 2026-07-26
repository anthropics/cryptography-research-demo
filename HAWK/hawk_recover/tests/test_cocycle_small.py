# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Small-ℓ smoke test: generate a random key at ell∈{5,6}, build Λ_τ, and
check that the secret cocycle V_τ = B⁻¹·τ(B) lies in it with ‖V_τ‖² = λ₁²."""

from __future__ import annotations

import random

import pytest

from hawk_recover import ring
from hawk_recover.cocycle import cocycle_lattice
from hawk_recover.params import HawkParams
from hawk_recover.types import PublicKey


@pytest.mark.parametrize("ell", [5, 6])
def test_vtau_is_shortest(ell: int) -> None:
    from sage.all import ZZ, vector

    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, tau, _ = ring.field_and_homs(ell)
    rng = random.Random(1000 + ell)

    def ro():
        return sum(rng.randint(-1, 1) * z**k for k in range(n))

    B = ring.gidentity(L, 2)
    for _ in range(3):
        B = B * ring.gmat(L, [[1, ro()], [0, 1]]) * ring.gmat(L, [[1, 0], [ro(), 1]])
    Q = ring.star(B, c) * B
    Q_list = tuple(
        [int(co) for co in list(Q[r, cc])] + [0] * (n - len(list(Q[r, cc])))
        for r in range(2)
        for cc in range(2)
    )
    pk = PublicKey(params=p, Q_list=Q_list, raw=b"")
    lat = cocycle_lattice(pk)

    Vtau = B.inverse() * ring.mat_apply(B, tau)
    amb: list[int] = []
    for r in range(2):
        for cc in range(2):
            co = list(Vtau[r, cc])
            amb += [int(ZZ(co[k] if k < len(co) else 0)) for k in range(n)]
    w = lat.LB.solve_left(vector(ZZ, amb))
    assert all(x.is_integer() for x in w), "V_τ not integral in Λ_τ basis"
    w = vector(ZZ, w)
    assert w * lat.Gram_tau * w == p.lam1sq
