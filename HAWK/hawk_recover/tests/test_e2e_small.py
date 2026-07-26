# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Full-pipeline test at small ℓ: cocycle → svp → intertwiner → tower.

Synthetic keys at ℓ=5–6 (n=16, 32) run the entire attack in seconds.
decode/encode need the real HAWK encoding (defined only for n ≥ 256), so they
are exercised by the CLI run on a real key, not here.
"""

from __future__ import annotations

import random

import pytest

from hawk_recover import ring
from hawk_recover.cocycle import cocycle_lattice
from hawk_recover.intertwiner import intertwiner
from hawk_recover.params import HawkParams
from hawk_recover.svp import find_shell_vectors
from hawk_recover.tower import recover_basis
from hawk_recover.types import PublicKey


def _synthetic_pk(ell: int, seed: int):
    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, _, _ = ring.field_and_homs(ell)
    rng = random.Random(seed)

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
    return PublicKey(params=p, Q_list=Q_list, raw=b"")


def test_intertwiner_when_cocycle_is_first_shell_vector(tmp_path) -> None:
    """Seed placing ±V_τ as shell vector 0: every pair through it fails
    (Remark 6.2's hypothesis), so V_τ must be taken directly."""
    from sage.all import ZZ, vector

    from hawk_recover.tower import _is_aside

    pk = _synthetic_pk(5, seed=1039)
    n = pk.params.n
    lat = cocycle_lattice(pk)
    shell = find_shell_vectors(lat, workdir=tmp_path)
    # Premise: this seed puts the a-side cocycle at shell vector 0.
    assert _is_aside(vector(ZZ, shell.W[0] * lat.LB), n)

    itw = intertwiner(shell, pk, lat)
    assert list(itw.wt) == list(shell.W[0])  # taken directly, not via a pair
    L, _, _, _, _ = ring.field_and_homs(5)
    itw.check(L)


@pytest.mark.parametrize("ell", [5, 6])
def test_e2e_synthetic(ell: int, tmp_path) -> None:
    # --- synthetic key ---
    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, _, _ = ring.field_and_homs(ell)
    rng = random.Random(4242 + ell)

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

    # --- the four attack steps (paper §2.1) ---
    lat = cocycle_lattice(pk, verbose=True)  # step 1
    shell = find_shell_vectors(lat, workdir=tmp_path, verbose=True)  # step 2
    assert len(shell.W) >= 2
    itw = intertwiner(shell, pk, lat, verbose=True)  # step 3
    rb = recover_basis(itw, pk, verbose=True)  # step 4

    # --- the recovered basis solves the LIP instance ---
    rb.check()
    Brec = rb.B(L)
    assert ring.star(Brec, c) * Brec == Q
