#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Numerically verify the paper's §4 statements on a random key at small ℓ:
Lemma 4.1 (V_τ ∈ Λ_B^(τ), det 1, on the shell), Lemma 4.2 (φ_τ carries
Λ_B^(τ) onto Λ_I^(τ) preserving det), Lemma 4.3 (rank n and the W(a,b)
parametrization), Lemma 4.4 (λ₁² and the count 2(n/2+1)), and Prop 4.5 (the
explicit orthogonal basis and its Gram, i.e. the isometry class).

Norms use this package's Gram convention Gram_τ = 8·Q^(τ) (types.CocycleLattice),
so the paper's minimum n/4 reads λ₁² = 2n and Prop 4.5's summands read
2n·Z^{n/2+1} ⊕ 4n·Z^{n/2−1}.

Usage: nix-shell nix/shell.nix --run 'sage -python scripts/verify_lemmas.py [ell ...]'
       (default ells: 4 5 6)
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from typing import Any  # noqa: E402

from hawk_recover import ring  # noqa: E402
from hawk_recover.cocycle import cocycle_lattice  # noqa: E402
from hawk_recover.params import HawkParams  # noqa: E402
from hawk_recover.tower import _amb_to_matrix  # noqa: E402
from hawk_recover.types import (  # noqa: E402
    CocycleLattice,
    PublicKey,
    SageMatrix,
    SageVector,
    require,
)


def to_amb(W: SageMatrix, n: int) -> list[int]:
    """2×2 matrix over L → ambient integer vector, index (2r+c)·n + k."""
    from sage.all import ZZ

    return [int(ZZ(list(W[r, c])[k])) for r in range(2) for c in range(2) for k in range(n)]


def coords(lat: CocycleLattice, amb: list[int]) -> SageVector:
    """Integer coordinates of an ambient vector in the LB basis; requires membership."""
    from sage.all import ZZ, vector

    w = lat.LB.solve_left(vector(ZZ, amb))
    require(all(x.is_integer() for x in w), "vector not in Λ_τ")
    return vector(ZZ, w)


def random_key(ell: int, seed: int) -> tuple[HawkParams, tuple[Any, ...], SageMatrix, PublicKey]:
    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, tau, sigma = ring.field_and_homs(ell)
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
    return p, (L, z, c, tau, sigma), B, PublicKey(params=p, Q_list=Q_list, raw=b"")


def identity_pk(p: HawkParams) -> PublicKey:
    one = [1] + [0] * (p.n - 1)
    zero = [0] * p.n
    return PublicKey(params=p, Q_list=(one, zero, zero, one), raw=b"")


def adj(Y: SageMatrix, L: Any) -> SageMatrix:
    return ring.gmat(L, [[Y[1, 1], -Y[0, 1]], [-Y[1, 0], Y[0, 0]]])


def count_shortest(lat: CocycleLattice) -> tuple[int, int]:
    """(λ₁², #vectors of norm λ₁²) by exact enumeration under Gram_τ."""
    from sage.all import ZZ, QuadraticForm

    U = lat.Gram_tau.LLL_gram()  # unimodular, Uᵀ·G·U is the LLL-reduced Gram
    Gt = U.transpose() * lat.Gram_tau * U
    sv = QuadraticForm(ZZ, 2 * Gt).short_vector_list_up_to_length(lat.lam1sq + 1)
    nonzero = [k for k in range(1, lat.lam1sq + 1) if sv[k]]
    require(bool(nonzero), f"no nonzero vector of norm ≤ {lat.lam1sq}: λ₁² exceeds it")
    return nonzero[0], len(sv[nonzero[0]])


def verify(ell: int, seed: int = 0) -> None:
    p, (L, z, c, tau, sigma), B, pk = random_key(ell, seed)
    n = p.n
    I2 = ring.gidentity(L, 2)
    Q = pk.Q(L)
    lat_B = cocycle_lattice(pk)
    lat_I = cocycle_lattice(identity_pk(p))
    print(f"ℓ={ell}  n={n}  (seed {seed})")

    # Lemma 4.1: V_τ = B⁻¹τ(B) has det 1, satisfies (C1), (C2), lies in Λ_B^(τ)
    # and attains the minimum λ₁² = 2n.
    V = B.inverse() * ring.mat_apply(B, tau)
    require(V.det() == 1, "det V_τ ≠ 1")
    require(ring.mat_apply(V, tau) == adj(V, L), "(C1) fails")
    require(Q * V == ring.mat_apply(V, sigma).transpose() * ring.mat_apply(Q, tau), "(C2) fails")
    w = coords(lat_B, to_amb(V, n))
    require(w * lat_B.Gram_tau * w == p.lam1sq, "V_τ not on the shell")
    print("  Lemma 4.1  ✓  V_τ ∈ Λ_B^(τ), det 1, (C1)&(C2), norm λ₁²")

    # Lemma 4.2: φ_τ(Y) = B·Y·τ(B)⁻¹ carries each basis vector of Λ_B^(τ) to
    # a point of Λ_I^(τ) with the same det, and the images' coordinate matrix
    # in Λ_I's basis is unimodular — so φ_τ(Λ_B) = Λ_I exactly (onto).
    from sage.all import matrix

    tB_inv = ring.mat_apply(B, tau).inverse()
    img_coords = []
    for row in lat_B.LB.rows():
        Y = _amb_to_matrix(row, L, z, n)
        W = B * Y * tB_inv
        require(W.det() == Y.det(), "φ_τ changes det")
        img_coords.append(coords(lat_I, to_amb(W, n)))  # image lies in Λ_I^(τ)
    C_img = matrix(img_coords)
    require(abs(C_img.det()) == 1, "φ_τ(Λ_B^(τ)) ≠ Λ_I^(τ) (not onto)")
    print("  Lemma 4.2  ✓  φ_τ(Λ_B^(τ)) = Λ_I^(τ) exactly, det-preserving")

    # Lemma 4.3 + Prop 4.5: the explicit orthogonal basis of Λ_I^(τ) —
    # I; W(0,ζ^j) for odd j; W(a,0) for a = ζ^j + (−1)^{j+1}ζ^{n−j} —
    # is a Z-basis (unimodular change from LB) with Gram diag(2n,…,4n,…).
    def Wab(a: Any, b: Any) -> SageMatrix:
        return ring.gmat(L, [[a, b], [sigma(b), tau(a)]])

    basis = [I2] + [Wab(0, z**j) for j in range(1, n, 2)]
    basis += [Wab(z**j + (-1) ** (j + 1) * z**(n - j), 0) for j in range(1, n // 2)]
    require(len(basis) == n, "rank ≠ n")
    C = matrix([coords(lat_I, to_amb(Wb, n)) for Wb in basis])
    require(abs(C.det()) == 1, "explicit family is not a basis of Λ_I^(τ)")
    G = C * lat_I.Gram_tau * C.transpose()
    p_, q_ = n // 2 + 1, n // 2 - 1
    s2, s4 = 2 * n, 4 * n
    expect = [s2] * p_ + [s4] * q_
    require(G.is_diagonal() and sorted(G.diagonal()) == sorted(expect), "Gram ≠ diag(2n,…,4n,…)")
    print("  Lemma 4.3  ✓  rank n; the W(a,b) family is a Z-basis (unimodular over LB)")
    print("  Prop 4.5   ✓  Gram_τ = diag(%d × %d, %d × %d)" % (p_, s2, q_, s4))
    print("                i.e. Λ ≅ %d·Z^%d ⊕ %d·Z^%d" % (s2, p_, s4, q_))

    # Lemma 4.4: λ₁² = 2n, attained by exactly 2(n/2+1) vectors — checked on
    # the keyed lattice Λ_B^(τ) by exact enumeration.
    lam1sq, count = count_shortest(lat_B)
    require(lam1sq == p.lam1sq, f"λ₁²={lam1sq} ≠ {p.lam1sq}")
    require(count == 2 * (n // 2 + 1), f"#shortest={count} ≠ {2 * (n // 2 + 1)}")
    print(f"  Lemma 4.4  ✓  λ₁² = {lam1sq} = 2n, attained by {count} = 2(n/2+1) vectors")


if __name__ == "__main__":
    ells = [int(a) for a in sys.argv[1:]] or [4, 5, 6]
    for ell in ells:
        verify(ell)
        print()
    print("all lemma checks passed")
